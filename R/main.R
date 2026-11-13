
#' Encode DNA Sequences
#'
#' @description
#' Converts DNA sequences into categorical nucleotide features
#' for machine learning prediction.
#'
#' @param dna_strings A character vector containing DNA sequences
#' of equal length. The m6A model uses 5-mer sequences.
#'
#' @return A data frame containing one factor column for each
#' nucleotide position, named nt_pos1, nt_pos2, and so on.
#'
#' @examples
#' m6APrediction:::dna_encoding(c("GGACA", "GAACA"))
dna_encoding <- function(dna_strings){
  nn <- nchar(dna_strings[1])
  seq_m <- matrix(unlist(strsplit(dna_strings, "")),
                  ncol = nn, byrow = TRUE)
  colnames(seq_m) <- paste0("nt_pos", 1:nn)
  seq_df <- as.data.frame(seq_m)
  seq_df[] <- lapply(seq_df, factor,
                     levels = c("A", "T", "C", "G"))
  return(seq_df)
}


#' Predict m6A Sites for Multiple Samples
#'
#' @description
#' Predicts m6A modification probabilities and classification
#' labels for multiple RNA sites using a trained random forest
#' model. DNA sequences are encoded into nucleotide features
#' before prediction.
#'
#' @param ml_fit A trained random forest classification model.
#' @param feature_df A data frame containing gc_content,
#' RNA_type, RNA_region, exon_length, distance_to_junction,
#' evolutionary_conservation, and DNA_5mer.
#' @param positive_threshold Probability threshold used to
#' classify positive m6A sites. Default is 0.5.
#'
#' @return A data frame containing the original features and
#' two additional columns: predicted_m6A_prob (the positive
#' class probability) and predicted_m6A_status
#' (Positive or Negative).
#'
#' @import randomForest
#' @importFrom stats predict
#' @export
#'
#' @examples
#' ml_fit <- readRDS(system.file(
#'   "extdata", "rf_fit.rds",
#'   package = "m6APrediction", mustWork = TRUE))
#'
#' example_df <- read.csv(system.file(
#'   "extdata", "m6A_input_example.csv",
#'   package = "m6APrediction", mustWork = TRUE))
#'
#' prediction_multiple(ml_fit, example_df[1:3, ])

prediction_multiple <- function(ml_fit, feature_df,
                                positive_threshold = 0.5){

  stopifnot(all(c(
    "gc_content", "RNA_type", "RNA_region",
    "exon_length", "distance_to_junction",
    "evolutionary_conservation", "DNA_5mer"
  ) %in% colnames(feature_df)))

  feature_df$RNA_type <- factor(
    feature_df$RNA_type,
    levels = c("mRNA", "lincRNA", "lncRNA", "pseudogene")
  )

  feature_df$RNA_region <- factor(
    feature_df$RNA_region,
    levels = c("CDS", "intron", "3'UTR", "5'UTR")
  )

  infered_prob <- predict(
    ml_fit,
    newdata = cbind(
      feature_df,
      dna_encoding(feature_df$DNA_5mer)
    ),
    type = "prob"
  )

  feature_df$predicted_m6A_prob <- infered_prob[, "Positive"]

  feature_df$predicted_m6A_status <- "Negative"

  feature_df$predicted_m6A_status[
    feature_df$predicted_m6A_prob > positive_threshold
  ] <- "Positive"

  return(feature_df)
}


#' Predict m6A Status for a Single RNA Site
#'
#' @description
#' Predicts the m6A modification probability and classification
#' label for a single RNA site using a trained random forest
#' model.
#'
#' @param ml_fit A trained random forest classification model.
#' @param gc_content Numeric GC-content feature.
#' @param RNA_type RNA biotype, such as mRNA or lincRNA.
#' @param RNA_region RNA region, such as CDS, intron, 3'UTR,
#' or 5'UTR.
#' @param exon_length Numeric exon-length feature.
#' @param distance_to_junction Numeric distance-to-junction feature.
#' @param evolutionary_conservation Numeric conservation feature.
#' @param DNA_5mer A character string containing a 5-mer DNA sequence.
#' @param positive_threshold Probability threshold for positive
#' classification. Default is 0.5.
#'
#' @return A named character vector containing
#' predicted_m6A_prob and predicted_m6A_status.
#' The probability is converted to character because
#' the vector also contains a classification label.
#'
#' @export
#'
#' @examples
#' ml_fit <- readRDS(system.file(
#'   "extdata", "rf_fit.rds",
#'   package = "m6APrediction", mustWork = TRUE))
#'
#' prediction_single(
#'   ml_fit = ml_fit,
#'   gc_content = 0.6,
#'   RNA_type = "mRNA",
#'   RNA_region = "CDS",
#'   exon_length = 12,
#'   distance_to_junction = 5,
#'   evolutionary_conservation = 0.8,
#'   DNA_5mer = "GGACA"
#' )

prediction_single <- function(
    ml_fit, gc_content, RNA_type, RNA_region,
    exon_length, distance_to_junction,
    evolutionary_conservation, DNA_5mer,
    positive_threshold = 0.5){

  pred_outcome <- prediction_multiple(
    ml_fit = ml_fit,
    data.frame(
      gc_content = gc_content,
      RNA_type = factor(
        RNA_type,
        levels = c("mRNA", "lincRNA", "lncRNA", "pseudogene")
      ),
      RNA_region = factor(
        RNA_region,
        levels = c("CDS", "intron", "3'UTR", "5'UTR")
      ),
      exon_length = exon_length,
      distance_to_junction = distance_to_junction,
      evolutionary_conservation = evolutionary_conservation,
      DNA_5mer = DNA_5mer
    ),
    positive_threshold = positive_threshold
  )

  returned_vector <- c(
    pred_outcome$predicted_m6A_prob,
    pred_outcome$predicted_m6A_status
  )

  names(returned_vector) <- c(
    "predicted_m6A_prob",
    "predicted_m6A_status"
  )

  return(returned_vector)
}
