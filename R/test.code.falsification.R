.test.code.falsification <- function(barcodes=character(), d=1, metric=c("hamming","seqlev","levenshtein","phaseshift"), cores=detectCores()/2) {
  metric <- match.arg(metric)

  # Calculated mutated closure of barcodes
  cf <- .code_falsification(barcodes, d, metric, cores)

  # Demultiplex them
  dm <- demultiplex(cf$mutation, barcodes, metric=metric)

  if (sum(cf$barcode != dm) > 0) {
    show("The barcode set could not be correctly demultiplexed.")
  } else {
    show("Correct.")
  }
  return(sum(cf$barcode != dm) / nrow(cf))
}

