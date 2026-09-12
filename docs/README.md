# Project website

Static site served at <https://saimoon-oman.github.io/SAGE-GI/>. No build step, no external
requests — everything is served from this folder.

It lives in `docs/` because GitHub Pages' branch deployment offers only two choices of
folder: the repository root or `/docs`. An arbitrary folder name would not work.

## Publishing

1. Push the repository to GitHub.
2. **Settings → Pages → Build and deployment → Source: Deploy from a branch.**
3. Branch `main`, folder **`/docs`**. Save.

The site is live within a minute or two. `.nojekyll` is present so GitHub serves the files
as-is instead of running them through Jekyll.

## Keeping it current

`scripts/91_make_tables.py` writes `results/table1_main.csv`; copy it into `assets/` and the
results table on the page fills itself in from that file.

```bash
cp results/table1_main.csv docs/assets/
cp figures/fig*.png        docs/assets/
cp paper/main.pdf          docs/paper.pdf
```

## Local preview

```bash
python -m http.server 8000 --directory docs
```
