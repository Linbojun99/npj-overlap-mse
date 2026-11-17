# npj-overlap-mse

R code for *Overlap-area implementation between Pacific tuna commissions alters the performance of a South Pacific albacore management procedure*.

Requires R ≥ 4.3 with `dplyr`, `tidyr`, `readr` and `jsonlite`. Download the WCPFC and IATTC public-domain longline data ([WCPFC](https://www.wcpfc.int/sustainability/scientific-data/wcpfc-public-domain-aggregated-catcheffort-data-download-page), [IATTC](https://www.iattc.org/en-US/Data/Public-domain)) into `data/raw/`.

Run the scripts in `R/` in numerical order from the repository root. `07_experiment.R` is run once for each of the 12 configurations listed in `00_settings.R`, and `10_gear_structure.R` once for each of its six structures; the Article uses `N_REP=300`.

```bash
CONFIG=reference N_REP=300 Rscript R/07_experiment.R
```
