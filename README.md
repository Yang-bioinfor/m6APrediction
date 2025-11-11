m6APrediction: A Machine Learning Tool for Predicting m6A Modification
Sites in RNA
================

- [Quick Start](#quick-start)
- [Batch Prediction (Multiple
  Samples)](#batch-prediction-multiple-samples)
- [Single-Sample Prediction](#single-sample-prediction)
- [Model Performance](#model-performance)
- [Features](#features)
- [Contributing](#contributing)
- [License](#license)
- [Citation](#citation)




    # m6APrediction <small>(Version 1.0.0)</small>

    [![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
    [![R build status](https://github.com/Yang-bioinfor/m6APrediction/workflows/R-CMD-check/badge.svg)](https://github.com/Yang-bioinfor/m6APrediction/actions)

    ---

    **m6APrediction** is an R package for predicting **N6-methyladenosine (m6A) modification sites** in RNA sequences. m6A is a common internal RNA modification involved in gene expression regulation, RNA stability, and disease mechanisms such as cancer.

    The package leverages a **random forest machine learning model** trained on sequence and structural features, achieving high predictive accuracy validated by ROC and Precision-Recall Curve analyses.

    ---

    ## Key Functions

    * **`prediction_multiple()`**: Batch prediction for multiple samples from a feature data frame.
    * **`prediction_single()`**: Single-sample prediction using individual feature inputs.

    **Features used by the model include:**
    GC content, RNA type/region, exon length, distance to splice junctions, evolutionary conservation, and one-hot encoded 5-mer DNA sequences.
    Output includes predicted probabilities and binary statuses ("Positive" or "Negative").

    ---

    ## Installation

    **Recommended method:**

    ```r
    # Install remotes if needed
    if (!require("remotes")) install.packages("remotes")

    # Install m6APrediction from GitHub
    remotes::install_github("Yang-bioinfor/m6APrediction")

**Alternative using devtools:**

``` r
if (!require("devtools")) install.packages("devtools")
devtools::install_github("Yang-bioinfor/m6APrediction")
```

**Dependencies:** `randomForest`, `stats` (installed automatically). **R
version requirement:** R \>= 3.5.0

------------------------------------------------------------------------

## Quick Start

``` r
library(m6APrediction)

# Load pre-trained model
model <- readRDS(system.file("extdata", "rf_fit.rds", package = "m6APrediction"))

# Load example input data (multiple samples)
example_df <- read.csv(system.file("extdata", "m6A_input_example.csv", 
                                   package = "m6APrediction"))
head(example_df)
```

------------------------------------------------------------------------

## Batch Prediction (Multiple Samples)

``` r
# Predict m6A status for multiple samples
results <- prediction_multiple(model, example_df, positive_threshold = 0.5)

# View DNA 5-mer, predicted probability, and status
head(results[, c("DNA_5mer", "predicted_m6A_prob", "predicted_m6A_status")])
```

> **Note:** The output augments the input with `predicted_m6A_prob`
> (probability of “Positive”) and `predicted_m6A_status`
> (“Negative”/“Positive”).

------------------------------------------------------------------------

## Single-Sample Prediction

``` r
single_result <- prediction_single(
  ml_fit = model,
  gc_content = 0.45,
  RNA_type = "mRNA",
  RNA_region = "CDS",
  exon_length = 120,
  distance_to_junction = 30,
  evolutionary_conservation = 0.8,
  DNA_5mer = "ATGCC",
  positive_threshold = 0.5
)

print(single_result)
```

> Returns a named vector: `predicted_m6A_prob` and
> `predicted_m6A_status`.

------------------------------------------------------------------------

## Model Performance

The random forest model was trained and validated on RNA sequence
datasets. Cross-validation demonstrated high accuracy:

- **ROC AUC \> 0.90**
- **PRC AUPRC \> 0.85**

> ⚠️ **Tip:** To generate plots locally, use `pROC` and `PRROC` packages
> with your validation data. Commit PNG files to display on GitHub.

------------------------------------------------------------------------

## Features

| Feature | Description |
|----|----|
| **One-Hot Encoding** | Built-in `dna_encoding()` for 5-mer sequences |
| **Flexible Threshold** | Adjustable `positive_threshold` for classification |
| **Reproducible** | Includes pre-trained model and example data |
| **Lightweight** | Minimal dependencies; runs on standard R setups (R \>= 3.5.0) |

------------------------------------------------------------------------

## Contributing

Contributions are welcome!

- Fork the repository
- Make changes and submit a pull request
- Open a GitHub issue for questions or bugs

> **Tip:** Follow R coding style (`styler`/`lintr`) and test with
> example datasets.

------------------------------------------------------------------------

## License

MIT License. See
[LICENSE](https://github.com/Yang-bioinfor/m6APrediction/blob/main/LICENSE)
for details.

------------------------------------------------------------------------

## Citation

Hongyao Yang. (2025). *m6APrediction: A Machine Learning Tool for
Predicting m6A Modification Sites in RNA*. R package version 1.0.0.
<https://github.com/Yang-bioinfor/m6APrediction>

**BibTeX example:**

``` bibtex
@software{Yang2025m6APrediction,
  author = {Hongyao Yang},
  title = {m6APrediction: A Machine Learning Tool for Predicting m6A Modification Sites in RNA},
  year = {2025},
  url = {https://github.com/Yang-bioinfor/m6APrediction}
}
```
