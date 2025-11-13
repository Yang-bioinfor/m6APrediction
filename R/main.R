#' @import randomForest
NULL  # 防止 roxygen2 把 import 放在函数里

#' One-hot encode DNA 5-mers
#'
#' Convert a character vector of DNA 5-mers (e.g., "ATGCC") into a one-hot
#' encoded data frame where each nucleotide position becomes a factor column
#' with levels "A", "T", "C", "G".
#'
#' @param dna_strings Character vector of equal-length DNA strings (all must
#' have the same length).
#'
#' @return A data frame with `nchar(dna_strings[1])` columns named
#' `nt_pos1`, `nt_pos2`, …, each column being a factor with levels
#' `c("A", "T", "C", "G")`.
#'
#' @examples
#' dna <- c("ATGCC", "GCTAA")
#' dna_encoding(dna)
#'
#' @export
dna_encoding <- function(dna_strings) {
  nn <- nchar(dna_strings[1])
  seq_m <- matrix(unlist(strsplit(dna_strings, "")), ncol = nn, byrow = TRUE)
  colnames(seq_m) <- paste0("nt_pos", 1:nn)
  seq_df <- as.data.frame(seq_m)
  seq_df[] <- lapply(seq_df, factor, levels = c("A", "T", "C", "G"))
  return(seq_df)
}

#' Predict m6A status for multiple samples
#'
#' Given a data frame of features and a fitted random forest model,
#' return the original data frame augmented with two new columns:
#' `predicted_m6A_prob` (probability of class "Positive") and
#' `predicted_m6A_status` ("Positive" or "Negative").
#'
#' @param ml_fit A fitted `randomForest` model (must contain class "Positive").
#' @param feature_df Data frame that **must** contain the columns:
#' `gc_content`, `RNA_type`, `RNA_region`, `exon_length`,
#' `distance_to_junction`, `evolutionary_conservation`, and `DNA_5mer`.
#' @param positive_threshold Numeric threshold (default 0.5) above which a
#' sample is classified as "Positive".
#'
#' @return The input `feature_df` with two additional columns:
#' `predicted_m6A_prob` (numeric) and `predicted_m6A_status`
#' (factor with levels `c("Negative", "Positive")`).
#'
#' @examples
#' # Load example model and dataset included in the package
#' model <- readRDS(system.file("extdata", "rf_fit.rds",
#'                              package = "m6APrediction"))
#' df <- read.csv(system.file("extdata", "m6A_input_example.csv",
#'                            package = "m6APrediction"))
#' head(prediction_multiple(model, df))
#'
#' @importFrom stats predict
#' @export
prediction_multiple <- function(ml_fit, feature_df,
                                positive_threshold = 0.5) {
  needed <- c("gc_content", "RNA_type", "RNA_region", "exon_length",
              "distance_to_junction", "evolutionary_conservation", "DNA_5mer")
  stopifnot(all(needed %in% colnames(feature_df)))

  # Add one-hot encoded DNA features
  feature_df <- cbind(feature_df, dna_encoding(feature_df$DNA_5mer))

  # Convert categorical variables to factors with consistent levels
  feature_df$RNA_type <- factor(feature_df$RNA_type,
                                levels = c("mRNA", "lincRNA", "lncRNA", "pseudogene"))
  feature_df$RNA_region <- factor(feature_df$RNA_region,
                                  levels = c("CDS", "intron", "3'UTR", "5'UTR"))

  # Predict probabilities
  probs <- predict(ml_fit, newdata = feature_df, type = "prob")[, "Positive"]

  # Convert probabilities to binary labels
  status <- ifelse(probs > positive_threshold, "Positive", "Negative")

  feature_df$predicted_m6A_prob <- probs
  feature_df$predicted_m6A_status <- factor(status,
                                            levels = c("Negative", "Positive"))

  return(feature_df)
}

#' Predict m6A status for a single sample
#'
#' Wrapper around `prediction_multiple()` that accepts individual feature
#' values instead of a data frame. Returns a named vector with the predicted
#' probability and status.
#'
#' @param ml_fit Fitted `randomForest` model (same as in `prediction_multiple()`).
#' @param gc_content Numeric GC content.
#' @param RNA_type Character; one of `"mRNA"`, `"lincRNA"`, `"lncRNA"`,
#' `"pseudogene"`.
#' @param RNA_region Character; one of `"CDS"`, `"intron"`, `"3'UTR"`,
#' `"5'UTR"`.
#' @param exon_length Numeric length of the exon.
#' @param distance_to_junction Numeric distance to nearest splice junction.
#' @param evolutionary_conservation Numeric conservation score.
#' @param DNA_5mer Character 5-mer DNA sequence (e.g., `"ATGCC"`).
#' @param positive_threshold Numeric threshold (default 0.5).
#'
#' @return Named vector of length 2:
#' `predicted_m6A_prob` (numeric) and `predicted_m6A_status`
#' (character `"Positive"` or `"Negative"`).
#'
#' @examples
#' # Load example model included in the package
#' model <- readRDS(system.file("extdata", "rf_fit.rds",
#'                              package = "m6APrediction"))
#' prediction_single(
#'   ml_fit = model,
#'   gc_content = 0.45,
#'   RNA_type = "mRNA",
#'   RNA_region = "CDS",
#'   exon_length = 120,
#'   distance_to_junction = 30,
#'   evolutionary_conservation = 0.8,
#'   DNA_5mer = "ATGCC"
#' )
#'
#' @export
prediction_single <- function(ml_fit,
                              gc_content,
                              RNA_type,
                              RNA_region,
                              exon_length,
                              distance_to_junction,
                              evolutionary_conservation,
                              DNA_5mer,
                              positive_threshold = 0.5) {

  # Create a single-row data frame
  single_df <- data.frame(
    gc_content = gc_content,
    RNA_type = factor(RNA_type,
                      levels = c("mRNA", "lincRNA", "lncRNA", "pseudogene")),
    RNA_region = factor(RNA_region,
                        levels = c("CDS", "intron", "3'UTR", "5'UTR")),
    exon_length = exon_length,
    distance_to_junction = distance_to_junction,
    evolutionary_conservation = evolutionary_conservation,
    DNA_5mer = DNA_5mer,
    stringsAsFactors = FALSE
  )

  # Add one-hot encoding
  single_df <- cbind(single_df, dna_encoding(single_df$DNA_5mer))

  # Use the multiple-sample prediction function
  pred_df <- prediction_multiple(ml_fit, single_df, positive_threshold)

  # Return a clean named vector
  returned_vector <- c(
    predicted_m6A_prob = pred_df$predicted_m6A_prob[1],
    predicted_m6A_status = as.character(pred_df$predicted_m6A_status[1])
  )
  return(returned_vector)
}
