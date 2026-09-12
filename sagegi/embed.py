r"""Tile embedding extraction from a frozen genome foundation model.

Why tiles?
----------
SSG-LUGIA scans a 10 kb window with a 100 bp step, i.e. every base is seen by ~100
windows.  Embedding each window independently would cost ~100 x (genome length) of
transformer compute — about 48,500 CPU-hours for the 485 Mbp IslandPick benchmark on a
commodity laptop, and still ~200 GPU-hours.  That is why no genome-scale study has used
a foundation model this way before.

SAGE-GI instead embeds **non-overlapping tiles once** and reconstructs any window
embedding as the overlap-weighted mean of the tiles it covers (:func:`pool_windows`).
Compute becomes *linear* in genome length — a ~100x reduction — while the window grid,
and therefore the resolution of the downstream anomaly signal, is unchanged.  Pooling with
fractional weights keeps the signal continuous at the 100 bp window step even though tiles
are 1 kb apart.

The default tile length of 1 kb aggregates to the 10 kb sequence length DNABERT-S was
contrastively trained on, which is also SSG-LUGIA's window length.
"""
from __future__ import annotations
import os
from dataclasses import dataclass
from pathlib import Path

import numpy as np


@dataclass
class EmbedConfig:
    model_path: str
    tile: int = 1000
    stride: int | None = None          # defaults to `tile` (non-overlapping)
    batch: int = 16
    max_tokens: int = 4096
    device: str = "auto"
    dtype: str = "auto"                # fp16 on CUDA, fp32 on CPU

    @property
    def step(self) -> int:
        return self.stride or self.tile


def resolve_device(cfg: EmbedConfig):
    import torch
    if cfg.device != "auto":
        return torch.device(cfg.device)
    return torch.device("cuda" if torch.cuda.is_available() else "cpu")


def force_pytorch_attention() -> None:
    """Disable the vendored Triton flash-attention kernel.

    DNABERT-2/S bundle a `flash_attn_triton.py` written against an old Triton release; it
    breaks with current versions. At 1 kb tiles (~220 tokens) attention is a negligible
    share of the FLOPs, so the PyTorch path costs almost nothing and is reproducible.
    """
    import sys
    for name, mod in list(sys.modules.items()):
        if name.endswith("bert_layers") and hasattr(mod, "flash_attn_qkvpacked_func"):
            mod.flash_attn_qkvpacked_func = None


def make_half_safe(model) -> int:
    """Make the PyTorch attention path work in half precision.

    Two things in the vendored model code are hard-wired to fp32: the extended attention
    mask is built with ``.to(dtype=torch.float32)`` (its "fp16 compatibility" comment holds
    only for the Triton kernel, which casts internally), and the ALiBi tensor is a plain
    attribute rather than a registered buffer, so ``.half()`` never converts it. Both are
    summed into ``bias``; in fp16 that promotes the attention scores back to fp32 and the
    following ``matmul(attention_probs, v)`` raises "expected scalar type Float but found
    Half". Casting ``bias`` to the module dtype on the way in fixes both at once without
    editing the downloaded file.
    """
    import torch

    def pre_hook(mod, args):
        if len(args) >= 6 and torch.is_tensor(args[5]):
            dt = mod.Wqkv.weight.dtype
            if args[5].dtype != dt:
                args = list(args)
                args[5] = args[5].to(dt)
                args = tuple(args)
        return args

    n = 0
    for mod in model.modules():
        if type(mod).__name__ == "BertUnpadSelfAttention":
            mod.register_forward_pre_hook(pre_hook)
            n += 1
    return n


def load_model(cfg: EmbedConfig):
    """Load a frozen DNABERT checkpoint.

    ``AutoModel.from_pretrained`` cannot load DNABERT-2 on transformers >= ~4.45: it calls
    ``AutoModel.register``, which rejects the checkpoint because its model classes inherit
    ``config_class`` from the stock ``BertPreTrainedModel`` while its own
    ``configuration_bert.BertConfig`` subclasses ``PretrainedConfig`` directly. Loading the
    remote class explicitly skips that registration and works for both checkpoints.
    """
    import torch
    from transformers import AutoConfig, AutoTokenizer
    from transformers.dynamic_module_utils import get_class_from_dynamic_module

    dev = resolve_device(cfg)
    tok = AutoTokenizer.from_pretrained(cfg.model_path, trust_remote_code=True)
    conf = AutoConfig.from_pretrained(cfg.model_path, trust_remote_code=True)
    klass = get_class_from_dynamic_module("bert_layers.BertModel", cfg.model_path)
    model = klass.from_pretrained(cfg.model_path, config=conf)
    force_pytorch_attention()

    if cfg.dtype == "fp16" or (cfg.dtype == "auto" and dev.type == "cuda"):
        model = model.half()
        make_half_safe(model)
    model.eval().to(dev)
    if dev.type == "cpu":
        torch.set_num_threads(int(os.environ.get("SAGEGI_THREADS", os.cpu_count() or 4)))
    return tok, model, dev


def tile_bounds(genome_len: int, tile: int, step: int) -> np.ndarray:
    """Start offsets of the tiles covering the genome (last partial tile is dropped)."""
    n = max(0, (genome_len - tile) // step + 1)
    return np.arange(n, dtype=np.int64) * step


def embed_genome(seq: str, cfg: EmbedConfig, tok=None, model=None, dev=None,
                 progress=None) -> tuple[np.ndarray, np.ndarray]:
    """Return ``(embeddings (n_tiles, H) float16, tile_starts)`` for one genome."""
    import torch
    if model is None:
        tok, model, dev = load_model(cfg)
    starts = tile_bounds(len(seq), cfg.tile, cfg.step)
    out = None
    with torch.no_grad():
        for b0 in range(0, len(starts), cfg.batch):
            chunk = starts[b0:b0 + cfg.batch]
            batch = [seq[s:s + cfg.tile] for s in chunk]
            enc = tok(batch, return_tensors="pt", padding=True, truncation=True,
                      max_length=cfg.max_tokens)
            ids = enc["input_ids"].to(dev)
            mask = enc["attention_mask"].to(dev)
            hidden = model(input_ids=ids, attention_mask=mask)[0]
            m = mask.unsqueeze(-1).to(hidden.dtype)
            pooled = (hidden * m).sum(1) / m.sum(1).clamp(min=1)   # mean over real tokens
            pooled = pooled.float().cpu().numpy().astype(np.float16)
            if out is None:
                out = np.empty((len(starts), pooled.shape[1]), dtype=np.float16)
            out[b0:b0 + len(chunk)] = pooled
            if progress is not None:
                progress(min(b0 + cfg.batch, len(starts)), len(starts))
    return out, starts


def save(path: str | Path, emb: np.ndarray, starts: np.ndarray, cfg: EmbedConfig,
         genome_len: int) -> None:
    np.savez_compressed(path, emb=emb, starts=starts, tile=cfg.tile, step=cfg.step,
                        genome_len=genome_len, model=os.path.basename(cfg.model_path))


def pool_windows(emb: np.ndarray, tile: int, step: int, window_starts: np.ndarray,
                 w: int) -> np.ndarray:
    r"""Overlap-weighted mean of tile embeddings over each window.

    For a window :math:`[s, s+w)` the pooled embedding is
    :math:`\sum_i o_i e_i / \sum_i o_i` where :math:`o_i` is the number of bases the window
    shares with tile *i*.  With non-overlapping tiles this is computed in O(1) per window
    from a prefix sum, with fractional weights for the two partial end tiles, so the result
    varies smoothly as the window slides by 100 bp.
    """
    if step != tile:
        raise NotImplementedError("closed-form pooling assumes non-overlapping tiles")
    e = emb.astype(np.float32)
    n_tiles, H = e.shape
    pref = np.zeros((n_tiles + 1, H), dtype=np.float32)
    np.cumsum(e, axis=0, out=pref[1:])

    s = np.asarray(window_starts, dtype=np.int64)
    t = s + w
    a = np.clip(s // tile, 0, n_tiles - 1)                 # first overlapping tile
    b = np.clip((t - 1) // tile, 0, n_tiles - 1)           # last overlapping tile

    # full interior tiles a+1 .. b-1
    lo, hi = np.minimum(a + 1, b), b
    acc = (pref[hi] - pref[lo]) * float(tile)      # each interior tile contributes `tile` bases
    wt = (hi - lo).astype(np.float32) * tile

    # partial contribution of the first tile
    ov_a = (np.minimum(t, (a + 1) * tile) - np.maximum(s, a * tile)).astype(np.float32)
    ov_a = np.clip(ov_a, 0, tile)
    acc += e[a] * ov_a[:, None]
    wt += ov_a

    # partial contribution of the last tile (only when it differs from the first)
    diff = b > a
    ov_b = (np.minimum(t, (b + 1) * tile) - np.maximum(s, b * tile)).astype(np.float32)
    ov_b = np.clip(ov_b, 0, tile) * diff
    acc += e[b] * ov_b[:, None]
    wt += ov_b

    return acc / np.maximum(wt, 1e-6)[:, None]
