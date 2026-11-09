# Levels agree, changes do not: built-up growth in Bujumbura, 2016 to 2025, from three open products

How much has Bujumbura grown, and where? Four free products (Dynamic World, ESA WorldCover, GHSL, WorldPop) are compared for one city, with
terrain and distance analysis of new growth. No new classifier is trained, because no building labels exist for the city.
*Status: preprint draft, not peer reviewed.* Paper: [`paper/paper.pdf`](paper/paper.pdf).

## Main findings
* Built-up **levels agree** across products (correlation 0.87 to 0.88 on a 100 m grid), but totals for about 2020 range from 47 to 82 km2 (GHSL 47.3; Dynamic World 73.7; WorldCover 67.4 and 81.5).
* **Changes agree poorly**: Dynamic World vs GHSL change correlation is 0.17 at 100 m and 0.36 at 1 km; WorldPop change is unrelated to either. WorldCover's two releases differ by 14 km2 in one year, so it cannot be used for change.
* Dynamic World shows 25.7 km2 of new built land (2016 to 2018 vs 2023 to 2025; spatial bootstrap 21 to 30 km2) against 5.0 km2 in GHSL for 2015 to 2020.
* Growth is mostly infill and edge expansion (52% within 100 m of earlier built-up land, 0.6% beyond 1 km) and is over-represented on slopes of 15 degrees or more (7.3% of new area vs 2.5% of earlier built-up land), not on the low plain below 790 m.
* There are no ground labels, so the paper cannot say which product is closest to the truth.

## Reproduce
1. `pip install -r requirements.txt`.
2. The four rasters used in the paper are in `data/`. To regenerate them, run `bujumbura_urban_export.js` in the Earth Engine Code Editor.
3. From the repository root: `python analysis.py` (writes `results/results.json` and `results/masks.npz`), then `python make_figs.py` (writes `paper/figures/`). `explore.py` prints the data-quality checks.
4. `cd paper && pdflatex paper && bibtex paper && pdflatex paper && pdflatex paper`.

Data terms and attribution: `DATA.md`. Code MIT, paper CC BY 4.0, cite via `CITATION.cff`. Contact: Kana Amédée, kanaamedee00@gmail.com.
