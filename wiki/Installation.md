# Installation

## Requirements

- Python **3.10+** (tested on 3.13)
- ~2 GB disk for the two model checkpoints, ~500 MB for the benchmark genomes
- A GPU is optional. It is needed only for embedding extraction, and the provided notebook
  does that step on free Colab.

## Install

```bash
git clone https://github.com/saimoon-oman/SAGE-GI.git
cd SAGE-GI
pip install -r requirements.txt
```

## Why `transformers` is pinned

`requirements.txt` pins `transformers>=4.40,<5.0` and `tokenizers<0.21`. This is not
conservatism — it is required:

> DNABERT-2 and DNABERT-S ship **custom modelling code** (a MosaicBERT variant with ALiBi)
> loaded via `trust_remote_code=True`. That code was written against the transformers 4.x
> loader. Transformers 5 initialises modules on the `meta` device before loading weights, and
> the bundled code creates real tensors during `__init__`, so loading fails with
> `RuntimeError: Tensor on device meta is not on the expected device cpu!`

If you see that error, you are on transformers 5. Downgrade:

```bash
pip install "transformers==4.46.3" "tokenizers<0.21"
```

## Where large files go

Heavy intermediates never touch the repository. They live under `$SAGEGI_WORK`
(default `C:/sagegi_work`, or set it yourself):

```bash
export SAGEGI_WORK=/scratch/sagegi     # Linux/macOS
setx SAGEGI_WORK C:\sagegi_work        # Windows
```

```
$SAGEGI_WORK/
  genomes/   FASTA, one per accession
  models/    DNABERT-S, DNABERT-2 checkpoints
  emb/       tile embeddings, one .npz per (model, tile size, genome)
  comp/      cached compositional features and anomaly scores
  sweep/     per-genome hyper-parameter sweep tables
```

Only small tidy CSVs land in `results/`.

## Verify the install

```bash
python tests/test_sagegi.py
```

Expect `21/21 tests passed`. These cover k-mer counting, the window-pooling identity, the
SSG-LUGIA feature definitions, detrending, the change-point search, chimera construction and
the evaluation metrics — i.e. every place a silent error would corrupt downstream numbers.

## Optional: check the models load

```bash
python -c "
from transformers import AutoTokenizer, AutoModel
p='zhihan1996/DNABERT-S'
t=AutoTokenizer.from_pretrained(p, trust_remote_code=True)
m=AutoModel.from_pretrained(p, trust_remote_code=True)
print('ok', sum(q.numel() for q in m.parameters())/1e6, 'M params')"
```

A warning about Triton (`defaulting MosaicBERT attention implementation to pytorch`) is
expected and harmless — at 1 kb tiles attention is a negligible share of the compute.

**Next:** **[[Quick Start]]**
