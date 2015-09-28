
analyse.barcodes <- function(barcodes, metric = c("hamming", "seqlev", "levenshtein"), 
    cores = detectCores()/2, cost_sub = 1, cost_indel = 1) {
    metric <- match.arg(metric, several.ok = TRUE)
    
    if (!is.character(barcodes)) {
        stop("barcodes must be a vector of characters")
    }
    
    if (!(length(barcodes) > 1)) {
        stop("At least 2 barcodes need to be supplied.")
    }
    
    range_barcodes <- range(nchar(barcodes))
    
    if (range_barcodes[1] != range_barcodes[2]) {
        stop("All barcodes must be of equal length.")
    }
    
    # Test that only C,G,A, and T were used in the barcodes
    
    barcodes <- toupper(barcodes)
    
    barcode_chars <- unique(unlist(strsplit(barcodes, split = "")))
    
    if (sum(!(barcode_chars %in% c("C", "A", "G", "T"))) > 0) {
        stop("Only characters C,A,G, and T are allowed in the barcodes. (Lower case characters are accepted but make no difference)")
    }
    
    metric_vectors <- sapply(metric, .analyse_barcodes, barcodes = barcodes, cores = cores, 
        cost_sub = cost_sub, cost_indel = cost_indel, simplify = FALSE)
    
    metric_results <- lapply(metric_vectors, function(metric_res) {
        res <- c(mean_distance = mean(metric_res), median_distance = median(metric_res), 
            min_distance = min(metric_res), max_distance = max(metric_res), error_correction = floor((min(metric_res) - 
                1)/2), error_detection = floor((min(metric_res) - 1)))
        return(res)
    })
    metric_results <- data.frame(Description = c("Mean Distance", "Median Distance", 
        "Minimum Distance", "Maximum Distance", "Guaranteed Error Correction", "Guaranteed Error Detection"), 
        do.call(cbind.data.frame, metric_results))
    rownames(metric_results) <- NULL
    
    return(metric_results)
}
 
