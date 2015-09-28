
barcode.set.distances <- function(barcodes, metric = c("hamming", "seqlev", "levenshtein"), 
    cores = detectCores()/2, cost_sub = 1, cost_indel = 1) {
    metric <- match.arg(metric)
    
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
    
    x <- .barcode_set_distances(metric, barcodes, cores, cost_sub, cost_indel)
    
    return(pack(Matrix(x, dimnames = list(barcodes, barcodes))))
    
} 
