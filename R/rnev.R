#' Draw nested extreme value random variables
#'
#' `rnev` draws nested extreme value random variables.
#'
#' @param sigma should be a scalar.
#' @param nests should be a length-num_products vector listing each product's nest.
#' @param N the number of vectors;
#' @param sigma the parameter that measures within-nest correlation;
#' @param nests an integer-valued vector indicating the nest of each alternative;
#' @param tol the tolerance on the requirement that the approximate pdf be real-valued, non-negative, and equal to zero at its boundary;
#' @param lower_int an argument passed to fourierin (default depends on sigma);
#' @param upper_int an argument passed to fourierin (default depends on sigma);
#' @param lower_eval an argument passed to fourierin (default = -10);
#' @param upper_eval an argument passed to fourierin (default = 30);
#' @param resolution an argument passed to fourierin (default = 2^14).
#' @return An N-by-length(nests) matrix, with each row being a draw from the nested extreme value distribution.
#' @export
rnev = function(N, sigma, nests,
                tol=1e-3,
                lower_int = NULL, upper_int = NULL,
                lower_eval = -10, upper_eval = 30,
                resolution=2^14
                ) {
    # If lower_int and upper_int are not provided, set equal to \pm 100*sigma
    if (is.null(lower_int)) lower_int = -400*min(sigma, 0.1)
    if (is.null(upper_int)) upper_int =  400*min(sigma, 0.1)
    # If sigma < 0.01, throw a warning.
    if (sigma < 0.01) warning("sigma < 0.01 may result in a pdf with value > `tol' at the start of its domain, because the limit of integration required to avoid this behaviour would result in a NaN characteristic function")
    # Check nests is a vector of integers.
    if (!is.vector(nests) | !all(as.integer(nests) == nests)) {
        stop("`nests` must be a vector of integers")
    }
    # Assign each nest to a sequential integer.
    num_nests = length(unique(nests))
    xwalk = cbind(unique(nests), 1:num_nests)
    new_nests = sapply(nests, function(nest) xwalk[which(xwalk[,1] == nest), 2])
    # Check that...
    stopifnot(all(new_nests >= 1))
    stopifnot(num_nests == max(new_nests))
    # Count products
    num_products = length(new_nests)
    # Define the characteristic function.
    characteristic = function(z) exp(log(pracma::gammaz(1 - 1i*z/sigma)) - log(pracma::gammaz(1 - 1i*(1-sigma)*z/sigma)))
    # Integrate to find the PDF.
    out = fourierin::fourierin(f = characteristic, lower_int = lower_int, upper_int = upper_int,
                    lower_eval = lower_eval, upper_eval = upper_eval,
                    const_adj = -1, freq_adj = -1, resolution = resolution)
    # Check that...
     # pdf is defined
    if (any(is.nan(out$values))) stop("pdf not defined -- check whether characteristic function can be calculated at lower_int and upper_int")
     # imaginary values are small
    if (max(abs(Im(out$values))) > tol) stop("Imaginary values of the pdf have magnitude greater than `tol`")
     # pdf is almost non-negative
    if (min(Re(out$values) < tol)) stop("The pdf has values < -`tol'")
     # pdf end points are zero
    if (Re(out$values[1]) > tol) stop("The pdf has value > `tol' on the start of its domain")
    if (Re(out$values[resolution]) > tol) stop("The pdf has value > `tol' on the end of its domain")
    # Censor negative and imaginary parts out of the pdf.
    out$values_censored = sapply(out$values, function(v) max(0, Re(v)))
    # Integrate over the pdf to find the CDF.
    out$cdf = cumsum(out$values_censored)
    out$cdf = out$cdf/out$cdf[resolution]
    # Draw using inverse transform sampling, interpolating the inverse CDF.
    u = stats::runif(N*num_nests)
    index_before = sapply(u, function(ui) which(out$cdf < ui)[sum(out$cdf < ui)] )
    index_after  = index_before + 1
    cdf_before = out$cdf[index_before]
    cdf_after  = out$cdf[index_after]
    weight_on_after = (u - cdf_before)/(cdf_after - cdf_before)
    xi = weight_on_after*out$w[index_after] + (1-weight_on_after)*out$w[index_before]
    zeta = xi*sigma
    # Allocate zeta over draws.
    zeta_reshape = matrix(zeta, nrow=N, ncol=num_nests)
    zeta_allocated = zeta_reshape[, new_nests]
    # Now just add the standard GEV.
    gumbels = matrix(extraDistr::rgumbel(num_products*N), nrow=N)
    return((1-sigma)*gumbels + zeta_allocated)
}

