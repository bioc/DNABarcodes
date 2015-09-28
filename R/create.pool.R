
create.pool <- function(n, filter.triplets = TRUE, filter.gc = TRUE, filter.self_complementary = TRUE, 
    cores = detectCores()/2) {
    # Make sure that n and dist are integers
    n <- round(n)
    
    if (n < 1) {
        stop("Non-positive barcode lengths do not work. Try positive numbers.")
    }
    
    if (n > 20) {
        warning("n is greater than 20. I hope you know what you are doing. (Try n=7 or n=8 to get working results)")
    }
    
    # Do not filter gc content of n is 3 or smaller. That won't work anyway!
    filter.gc <- filter.gc & (n > 3)
    
    .create_pool(n, filter.triplets, filter.gc, filter.self_complementary, cores)
}
 
