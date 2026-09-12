#!/usr/bin/env python
"""Download every chromosome used in this study from NCBI Entrez.

Sources
-------
* 118 IslandPick chromosomes  (Langille et al. 2008, BMC Bioinformatics 9:329)
* Salmonella enterica serovar Typhi CT18   NC_003198.1   (case study)
* Corynebacterium diphtheriae NCTC13129    NC_002935.2   (HGT-gene evaluation)
* Pseudomonas aeruginosa LESB58            NC_011770.1   (HGT-gene evaluation)

FASTA files are written to ``$SAGEGI_WORK/genomes/<accession>.fna`` and are skipped when
already present, so the script is safe to re-run.
"""
from __future__ import annotations
import json, sys, threading, time, urllib.request, urllib.error
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from sagegi.paths import GENOMES, BENCH

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
EMAIL = "1024052020@grad.cse.buet.ac.bd"
TOOL = "SAGE-GI"

EXTRA = ["NC_003198.1", "NC_002935.2", "NC_011770.1"]


def fetch(acc: str, retries: int = 5) -> str:
    url = (f"{EUTILS}/efetch.fcgi?db=nuccore&id={acc}&rettype=fasta&retmode=text"
           f"&tool={TOOL}&email={EMAIL}")
    for attempt in range(retries):
        try:
            with urllib.request.urlopen(url, timeout=180) as r:
                txt = r.read().decode("utf-8", "replace")
            if txt.startswith(">") and len(txt) > 1000:
                return txt
            raise ValueError(f"short/invalid payload ({len(txt)} chars)")
        except Exception as exc:                       # noqa: BLE001
            wait = 3 * (attempt + 1)
            print(f"  ! {acc} attempt {attempt+1}/{retries} failed: {exc}; retry in {wait}s",
                  flush=True)
            time.sleep(wait)
    raise RuntimeError(f"could not download {acc}")


_GATE = threading.Semaphore(1)
_LAST = [0.0]


def _throttled_fetch(acc: str) -> str:
    """Serialise request *starts* to stay under the 3 req/s Entrez limit."""
    with _GATE:
        wait = 0.4 - (time.time() - _LAST[0])
        if wait > 0:
            time.sleep(wait)
        _LAST[0] = time.time()
    return fetch(acc)


def main() -> None:
    accs = json.loads((BENCH / "islandpick_accessions.json").read_text())
    accs = [a for a in accs if a != "Total"] + EXTRA
    todo = [a for a in accs
            if not (GENOMES / f"{a}.fna").exists()
            or (GENOMES / f"{a}.fna").stat().st_size <= 1000]
    print(f"{len(accs)} accessions, {len(todo)} to download", flush=True)

    done = 0
    with ThreadPoolExecutor(3) as ex:
        futs = {ex.submit(_throttled_fetch, a): a for a in todo}
        for f in as_completed(futs):
            acc = futs[f]
            done += 1
            try:
                (GENOMES / f"{acc}.fna").write_text(f.result(), encoding="utf-8")
                print(f"[{done}/{len(todo)}] {acc} ok", flush=True)
            except Exception as exc:                    # noqa: BLE001
                print(f"[{done}/{len(todo)}] {acc} FAILED {exc}", flush=True)

    have = sorted(p.stem for p in GENOMES.glob("*.fna"))
    print(f"done: {len(have)} FASTA files in {GENOMES}", flush=True)


if __name__ == "__main__":
    main()
