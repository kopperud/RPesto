#' fits the bds model
#'
#' @param phy an object of type phylo
#' @param sampling_fraction the probability that each species was sampled in the tree
#' @param lambda_hat the overall scale of the log-normal base distribution for the speciation rates. If not specified, the function will estimate it using ML
#' @param mu_hat the overall scale of the log-normal base distribution for the extinction rates. If not specified, the function will estimate it using ML
#' @param eta the shift rate parameter. If not specified, the function will estimate it using ML
#' @param num_classes the number of rate class discretizations (n), such that rate categories is k = n^2
#' @param sd the spread parameter for the log-normal base distribution
#' @param tol the local error threshold in the numerical ODE solver (per delta_t time step)
#'
#'
#' @export
fit_bds <- function(phy, sampling_fraction, lambda_hat, mu_hat, eta, num_classes = 6, sd = 0.587, tol = 1e-6){
    newick_string <- ape::write.tree(phy)
    phylogeny <- Phylogeny$new(newick_string)

    if (missing(lambda_hat) & missing(mu_hat)){
        ## fit constant-rate model
        bd_opt <- stats::optim(
            par = c(0.1, 0.05),
            fn = function(x) phylogeny$bd_likelihood(x[1], x[2], sampling_fraction, tol, FALSE),
            lower = c(0.00001, 0.00001),
            upper = c(2.0, 2.0),
            method = "L-BFGS-B",
            control = list(fnscale = -1),
            hessian = FALSE,
         )

        lambda_hat <- bd_opt$par[1]
        mu_hat <- bd_opt$par[2]
    }else{
        if (missing(lambda_hat) || missing(mu_hat)){
            stop("either specify lambda_hat and mu_hat, or neither") 
        }
    }

    if (missing(eta)){
        ## fit BDS model
        bds_opt <- stats::optim(
            par = c(0.005),
            fn = function(x) phylogeny$bds_likelihood(lambda_hat, mu_hat, x[1], sampling_fraction, sd, num_classes, tol, FALSE),
            lower = c(1e-8),
            upper = c(0.1),
            method = "Brent",
            control = list(fnscale = -1),
            hessian = FALSE,
         )

        eta <- bds_opt$par[1]
    }

    # store D(t) as a cubic spline 
    phylogeny$bds_likelihood(lambda_hat, mu_hat, eta, sampling_fraction, sd, num_classes, tol, TRUE)

    # calculate F(t)
    phylogeny$bds_preorder(lambda_hat, mu_hat, eta, sampling_fraction, sd, num_classes, tol);

    ## calculate netdiv per branch
    phylogeny$branch_rates(lambda_hat, mu_hat, eta, sampling_fraction, sd, num_classes);

    ## calculate no. shifts per branch
    phylogeny$number_of_shifts(lambda_hat, mu_hat, eta, sampling_fraction, sd, num_classes, tol);

    ## calculate Bayes factors per branch 
    phylogeny$bayes_factors(lambda_hat, mu_hat, eta, sampling_fraction, sd, num_classes, tol);

    ## write newick string
    s <- phylogeny$write_newick()

    ## parse newick string
    txt <- textConnection(s)
    tree <- treeio::read.beast.newick(txt)

    res <- list(
                "model" = c(
                            "lambda_hat" = lambda_hat,
                            "mu_hat" = mu_hat,
                            "eta" = eta
                            ),
                "td" = tree
                )

    return(res)
}
