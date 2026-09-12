#!/usr/bin/env python
"""Count the manuscript body words the way a journal editor would.

Briefings in Bioinformatics caps a Problem Solving Protocol at 5,000 words. That count is of
the body: not the abstract, not table content, not figure captions, and not the reference
list. This strips LaTeX markup and those elements so the number is meaningful.

    python scripts/92_wordcount.py
"""
from __future__ import annotations
import re
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
LIMIT = 5000
TARGET = 4700          # leave headroom for revision


def strip_latex(s: str) -> str:
    s = re.sub(r"(?m)^\s*%.*$", "", s)                                  # comments
    s = re.sub(r"\\caption\{(?:[^{}]|\{[^{}]*\})*\}", " ", s)           # captions don't count
    s = re.sub(r"\\(?:begin|end)\{[^}]*\}", " ", s)                     # environment delimiters
    s = re.sub(r"\\(?:input|includegraphics|label|ref|cite)\s*\{[^}]*\}", " ", s)
    s = re.sub(r"\\[a-zA-Z]+\*?(?:\[[^\]]*\])?", " ", s)                # remaining commands
    s = re.sub(r"[{}$&~^_\\]", " ", s)
    return s


def count(text: str) -> int:
    return len([w for w in strip_latex(text).split() if any(c.isalnum() for c in w)])


def main() -> None:
    total = 0
    for f in sorted((REPO / "paper" / "sections").glob("*.tex")):
        n = count(f.read_text(encoding="utf-8"))
        total += n
        print(f"  {f.name:28s} {n:5d}")
    print(f"  {'BODY TOTAL':28s} {total:5d}   (limit {LIMIT})")
    over = total - TARGET
    if over > 0:
        print(f"  {'over target by':28s} {over:5d}   (target {TARGET} for revision headroom)")
    else:
        print(f"  {'headroom':28s} {-over:5d}")
    sys.exit(0)


if __name__ == "__main__":
    main()
