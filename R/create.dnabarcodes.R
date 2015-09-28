
create.dnabarcodes.clique <- function(n, heuristic = "clique", ...) {
    create.dnabarcodes(n, heuristic = heuristic, ...)
}

create.dnabarcodes.conway <- function(n, heuristic = "conway", ...) {
    create.dnabarcodes(n, heuristic = heuristic, ...)
}

create.dnabarcodes.sampling <- function(n, heuristic = "sampling", iterations = 20000, 
    ...) {
    create.dnabarcodes(n, heuristic = heuristic, iterations = iterations, ...)
}

create.dnabarcodes.ashlock <- function(n, heuristic = "ashlock", iterations = 100, 
    population = 200, ...) {
    create.dnabarcodes(n, heuristic = heuristic, iterations = iterations, population = population, 
        ...)
}

create.dnabarcodes <- function(n, dist = 3, metric = c("hamming", "seqlev", "levenshtein", 
    "phaseshift"), heuristic = c("conway", "clique", "sampling", "ashlock"), filter.triplets = TRUE, 
    filter.gc = TRUE, filter.self_complementary = TRUE, pool = character(), iterations = 100, 
    population = 200, cores = detectCores()/2, use_cache = FALSE, cost_sub = 1, cost_indel = 1) {
    
    # Match argument to accepted list
    metric <- match.arg(metric)
    heuristic <- match.arg(heuristic)
    
    # Make sure that n and dist are integers
    n <- round(n)
    dist <- round(dist)
    
    if (!is.character(pool)) {
        stop("pool must be a vector of characters")
    }
    
    if (length(pool) > 0) {
        # Test the correct length of pool characters
        if (sum(nchar(pool) != n) > 0) {
            stop("All sequences in pool need to have the same length as parameter n")
        }
        
        # Test that only C,G,A, and T were used in the pool
        pool <- toupper(pool)
        pool_chars <- unique(unlist(strsplit(pool, split = "")))
        if (sum(!(pool_chars %in% c("C", "A", "G", "T"))) > 0) {
            stop("Only characters C,A,G, and T are allowed in the pool. (Lower case characters are accepted but make no difference)")
        }
    }
    
    if (dist > n) {
        stop("Distance d is greater than barcode length n. That will not work with any distance metric.")
    }
    
    if (dist < 0 | n <= 0) {
        stop("Non-positive barcode lengths or negative distances do not work. Try positive numbers.")
    }
    
    if (n > 20) {
        warning("n is greater than 20. I hope you know what you are doing. (Try n=7 or n=8 to get working results)")
    }
    
    # Do not filter gc content of n is 3 or smaller. That won't work anyway!
    filter.gc <- filter.gc & (n > 3)
    
    .create_dnabarcodes(n, dist, metric, heuristic, filter.triplets, filter.gc, filter.self_complementary, 
        pool, iterations, population, cores, use_cache, cost_sub, cost_indel)
}
 
