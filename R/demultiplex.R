# TODO: return data frame of barcode, read, distance
demultiplex <- function(reads, barcodes, metric = c("hamming", "seqlev", "levenshtein", 
    "phaseshift"), cost_sub = 1, cost_indel = 1) {
    metric <- match.arg(metric)
    
    barcode_length <- range(nchar(barcodes))
    
    if (!is.character(barcodes)) {
        stop("barcodes must be a vector of characters")
    }
    
    if (!is.character(reads)) {
        stop("reads must be a vector of characters")
    }
    
    if (!(length(barcodes) > 1)) {
        stop("At least 2 barcodes need to be supplied.")
    }
    
    if (!(length(reads) > 0)) {
        stop("At least one read need to be supplied.")
    }
    
    if (barcode_length[1] != barcode_length[2]) {
        stop("Barcodes need to be of the same length. (At the moment!)")
    }
    
    if (metric == "hamming") {
        if (min(nchar(reads)) < barcode_length[1]) {
            stop("All reads must at least be as long as the length of the barcodes.")
        }
        # We do not need anything after first couple of bases for hamming
        tmp_reads <- substr(reads, 1, barcode_length[1])
        return(data.frame(read = reads, .demultiplex(barcodes, tmp_reads, metric, 
            cost_sub, cost_indel), stringsAsFactors = FALSE))
    } else {
        return(data.frame(read = reads, .demultiplex(barcodes, reads, metric, cost_sub, 
            cost_indel), stringsAsFactors = FALSE))
    }
    
}
 
