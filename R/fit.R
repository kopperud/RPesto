#' fits the bds model
#'
#'
#'
#'
#'
#'
#' @export
fit_bds <- function(phy, rho, num_classes = 6, sd = 0.587, tol = 1e-6){
    newick_string <- ape::write.tree(phy)

    phylogeny <- Phylogeny$new(newick_string)

    ## fit constant-rate model
    bd_opt <- stats::optim(
        par = c(0.1, 0.05),
    #fn bd_likelihood(&self, lambda: f64, mu: f64, rho: f64, tol: f64) -> f64{
    #fn bd_likelihood(&self, lambda: f64, mu: f64, rho: f64, tol: f64) -> f64{
        fn = function(x) phylogeny$bd_likelihood(x[1], x[2], rho, tol),
        lower = c(0.00001, 0.00001),
        upper = c(2.0, 2.0),
        method = "L-BFGS-B",
        control = list(fnscale = -1),
        hessian = FALSE,
     )

    lambdaml <- bd_opt$par[1]
    muml <- bd_opt$par[2]


    ## fit BDS model
    bds_opt <- stats::optim(
        par = c(0.005),
    #pub fn bds_likelihood(&self, lambda_hat: f64, mu_hat: f64, eta: f64, rho: f64, sd: f64, n: usize, tol: f64) -> f64{
        fn = function(x) phylogeny$bds_likelihood(lambdaml, muml, x[1], rho, sd, num_classes, tol),
        lower = c(1e-8),
        upper = c(1.0),
        method = "Brent",
        control = list(fnscale = -1),
        hessian = FALSE,
     )

    res <- list(
            "lambda_hat" = bd_opt$par[1],
            "mu_hat" = bd_opt$par[2],
            "eta" = bds_opt$par[1]
            )

    return(res)
}
