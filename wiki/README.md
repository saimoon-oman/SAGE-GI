# Wiki pages

These Markdown files are the content of the GitHub wiki at
<https://github.com/saimoon-oman/SAGE-GI/wiki>. They live in the main repository so they are
versioned with the code; GitHub keeps wikis in a *separate* git repository, so they have to be
copied across once.

## Publishing them

**Option A — clone the wiki repo (recommended).**

```bash
# create the wiki first by adding any page through the GitHub web UI, then:
git clone https://github.com/saimoon-oman/SAGE-GI.wiki.git
cp wiki/*.md SAGE-GI.wiki/
cd SAGE-GI.wiki
git add . && git commit -m "Add SAGE-GI wiki" && git push
```

**Option B — paste through the web UI.** Repository → **Wiki** → **New Page**. Use the file
name without `.md` as the page title, exactly as spelled below, or the `[[links]]` will break.

## Pages and required titles

| File | Page title to use |
|---|---|
| `Home.md` | `Home` |
| `Installation.md` | `Installation` |
| `Quick-Start.md` | `Quick Start` |
| `Method-Overview.md` | `Method Overview` |
| `Benchmarks-and-Data.md` | `Benchmarks and Data` |
| `Pipeline-Reference.md` | `Pipeline Reference` |
| `API-Reference.md` | `API Reference` |
| `Reproducing-the-Paper.md` | `Reproducing the Paper` |
| `Troubleshooting.md` | `Troubleshooting` |
| `_Sidebar.md` | `_Sidebar` (renders as the navigation sidebar) |
| `_Footer.md` | `_Footer` (renders as the footer) |

GitHub maps a page titled `Quick Start` to the file `Quick-Start.md`, which is why the file
names use hyphens and the `[[Quick Start]]` links use spaces. Keep both as they are.
