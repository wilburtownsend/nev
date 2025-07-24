#' Draw nested extreme value random variables
#'
#' `rnev` draws nested extreme value random variables.
#'
#' @param sigma should be a scalar.
#' @param nests should be a length-num_products vector listing each product's nest.
#' @param N the number of vectors;
#' @param sigma the parameter that measures within-nest correlation, expressed either as a length-num_nests vector or as a scalar (in which case sigma is assumed constant across nests);
#' @param nests a vector of positive integers, indicating the nest of each alternative;
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
                resolution=2^15
                ) {
    # Check nests is a vector of integers.
    if (!is.vector(nests) | !all(as.integer(nests) == nests) | any(nests <= 0)) {
        stop("`nests` must be a vector of positive integers")
    }
    # Count products and nests
    num_products = length(nests)
    num_nests = max(nests)
    # Check sigma is a vector of either length-1 or length-num_nests.
    if (!is.vector(nests) | !is.numeric(sigma) | !(length(sigma) %in% c(1, num_nests))) {
        stop("`sigma` must be a scalar or a length-num_nests vector")
    }
    # Check all sigma is in (0,1].
    if ((min(sigma) <= 0) | (max(sigma) > 1)) {
        stop("`sigma` must be in (0, 1]")
    }
    # If any sigma < 0.01, throw a warning.
    if (min(sigma) < 0.01) warning("sigma < 0.01 may result in a pdf with value > `tol' at the start of its domain, because the limit of integration required to avoid this behaviour would result in a NaN characteristic function")

    # This function constructs the CDF, given a particular sigma value.
    construct_cdf = function(sigma) {
        # If lower_int and upper_int are not provided, set equal to \pm 100*sigma.
        if (is.null(lower_int)) lower_int = -400*min(sigma, 0.1)
        if (is.null(upper_int)) upper_int =  400*min(sigma, 0.1)
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
        # Return the pdf.
        return(list(x=out$w, P=out$cdf))
    }

    # This function does inverse transform sampling to draw from a CDF.
    # The argument CDF should be a list, with a vector x and a vector P.
    inverse_transform_sampling = function(N, cdf) {
        u = stats::runif(N)
        index_before = findInterval(u, cdf$P, all.inside=TRUE)
        index_after  = index_before + 1
        cdf_before = cdf$P[index_before]
        cdf_after  = cdf$P[index_after]
        weight_on_after = (u - cdf_before)/(cdf_after - cdf_before)
        Y = weight_on_after*cdf$x[index_after] + (1-weight_on_after)*cdf$x[index_before]
        return(Y)
    }

    # If only one value of sigma is provided, draw xi given that....
    if (length(sigma) == 1) {
        xi_cdf = construct_cdf(sigma)
        xi = inverse_transform_sampling(N*num_nests, xi_cdf)
        zeta = matrix(xi, nrow=N, ncol=num_nests)*sigma
    } else {
    # Othwerwise, draw the xi for each sigma
        zeta = matrix(nrow=N, ncol=num_nests)
        for (nest in 1:num_nests) {
            xi_cdf = construct_cdf(sigma[nest])
            xi = inverse_transform_sampling(N, xi_cdf)
            zeta[, nest] = xi*sigma[nest]
        }
    }
    # Allocate the zeta over draws.
    zeta_allocated = zeta[, nests]
    # Now just add the standard GEV.
    gumbels = matrix(extraDistr::rgumbel(num_products*N), nrow=N)
    return((1-sigma)*gumbels + zeta_allocated)
}

