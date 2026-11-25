# npj-overlap-mse

R code for *Allocation rules determine how errors in reported overlap catch affect South Pacific albacore management procedures*.

Requires R ≥ 4.3 with `dplyr`, `tidyr`, `readr` and `jsonlite`, and `ggplot2`, `patchwork`, `ggtext` and `ragg` for the figures. Download the WCPFC and IATTC public-domain longline data ([WCPFC](https://www.wcpfc.int/sustainability/scientific-data/wcpfc-public-domain-aggregated-catcheffort-data-download-page), [IATTC](https://www.iattc.org/en-US/Data/Public-domain)) into `data/raw/`.

Run the scripts from the repository root in this order. `07_experiment.R` runs steps 03–06 for one of the 12 configurations listed in `00_settings.R`; the Article uses `N_REP=300`.

```bash
Rscript R/01_public_data.R
Rscript R/02_overlap_share.R
for c in reference adult_M_0.24 adult_M_0.36 adult_M_0.54 steepness_0.65 steepness_0.95 \
         connectivity_0.5 connectivity_2 selectivity_2.5 selectivity_4.5 errors movement_information; do
  CONFIG=$c N_REP=300 Rscript R/07_experiment.R > output/07_$c.log &
done; wait
Rscript R/08_effects.R                    # effects, Supplementary Tables 9 and 14
N_REP=300 Rscript R/09_fixed_shares.R     # fixed reported shares under MP-B
N_REP=300 Rscript R/11_external_catch.R   # catch outside the scope of the MP
N_REP=300 Rscript R/12_unbiased_error.R   # unbiased reporting error
Rscript R/13_figures.R                    # Figs 1–6 and Supplementary Figs 1–4
```

Fig. 1a shows EEZ boundaries when `EEZ_SHP` points to the Marine Regions EEZ v12 shapefile with 0–360° longitudes (requires `sf`). `10_gear_structure.R` is an additional analysis with separate longline and troll fleets; no figure or table of the Article uses it.
