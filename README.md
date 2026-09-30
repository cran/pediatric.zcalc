[![R-CMD-check](https://github.com/bips-hb/pediatric.zcalc/actions/workflows/r.yml/badge.svg)](https://github.com/bips-hb/pediatric.zcalc/actions/workflows/r.yml)
# 📦 pediatric.zcalc

**pediatric.zcalc** provides functions to compute individual percentile ranks and z-scores for clinical biomarkers in children, adolescents and young adults, 
based on age-, sex-, and height-specific reference data from the IDEFICS study and the Biomarkers4Pediatrics collaboration.
It also supports the computation of a composite Metabolic Syndrome (MetS) score and associated monitoring/action levels for health monitoring.

## 🔧 Installation

You can install the development version of `pediatric.zcalc` from GitHub using:

```r
# Install devtools if not already installed
install.packages("devtools")

# Install directly from GitHub
devtools::install_github("bips-hb/pediatric.zcalc")
```

## 🚀 Example

```r
library(pediatric.zcalc)

# Input: data frame with raw values
df <- data.frame(
  sex = c("f", "m"),
  age = c(8, 15),          # years
  height = c(120, 125),    # cm
  waist = c(55, 60),       # cm
  homa = c(1.2, 1.4),      # HOMA-IR, dimensionless
  sbp = c(100, 105),       # mmHg
  dbp = c(65, 70),         # mmHg
  trg = c(80, 90),         # mg/dL
  hdl = c(50, 45),         # mg/dL
  crp = c(6, 2)            # mg/L
)

# Calculate z-scores and MetS
results <- ScoreCalc(df, return_values = c("z.score", "MetS"))

# View output
print(results)
```

## 📚 Functions

- `get_scores()` — Get percentiles or z-scores for a single variable.
- `ScoreCalc()` — Calculate scores for multiple variables and optionally MetS.
- `MetSScore()` — Calculate MetS from z-scores.
- `action_levels()` — Assign action levels based on percentile cutoffs.

## 📖 Reference

This package uses internal parameter tables based on the **IDEFICS study**, a European cohort focused on childhood obesity and metabolic health, and **Biomarkers4Pediatrics**,  an international multicohort pediatric biomarker collaboration:

-  Ahrens, W., Moreno, L. A., Mårild, S., Molnár, D., Siani, A., De Henauw, S., Böhmann, J., Günther, K., Hadjigeorgiou, C., Iacoviello, L., Lissner, L., Veidebaum, T., Pohlabeln, H. & Pigeot, I. on behalf of the IDEFICS consortium (2014). [Metabolic syndrome in young children: definitions and results of the IDEFICS study.](https://doi.org/10.1038/ijo.2014.130) *International Journal of Obesity*, 38 (Suppl 2), S4–S14.

- [http://www.biomarkers4pediatrics.eu/](https://www.bips-institut.de/en/biomarkers4pediatrics.html) 

The reference data used for all variables other than CRP in this package were originally published in:
- [bips-hb/IDEFICS-Score_Calculator](https://github.com/bips-hb/IDEFICS-Score_Calculator/) — a repository containing the IDEFICS tables and an R script version of the scoring logic.

This package improves upon and packages that logic cleanly for programmatic and research use.

