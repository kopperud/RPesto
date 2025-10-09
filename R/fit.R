#' fits the bds model
#'
#' @param phy an object of type phylo
#' @param sampling_fraction the probability that each species was sampled in the tree
#' @param lambda_hat the overall scale of the log-normal base distribution for the speciation rates. If not specified, the function will estimate it using ML
#' @param mu_hat the overall scale of the log-normal base distribution for the extinction rates. If not specified, the function will estimate it using ML
#' @param eta the shift rate parameter. If not specified, the function will estimate it using ML
#' @param num_speciation_classes the number of speciation rate class discretizations, such that rate categories is k = n_sp * n_mu
#' @param num_extinction_classes the number of extinction rate class discretizations, such that rate categories is k = n_sp * n_mu
#' @param sd the spread parameter for the log-normal base distribution
#' @param tol the local error threshold in the numerical ODE solver (per delta_t time step)
#' @param condition_survival whether or not to condition on the survival of the left and right lineages descending from the root (default TRUE)
#' @param condition_root_speciation whether or not to condition on that there was a speciation event at the root node (default TRUE)
#' @param condition_marginal whether or condition using the marginal or per-category approach (default FALSE, i.e., to condition per rate category)
#' @param extinction_approximation whether or not to approximate the extinction probability calculations, by assuming that rate shift events are not allowed on extinct lineages (default FALSE)
#'
#'
#' @export
fit_bds <- function(
        phy, 
        sampling_fraction,
        lambda_hat,
        mu_hat,
        eta,
        num_speciation_classes = 6,
        num_extinction_classes = 6,
        sd = 0.587,
        tol = 1e-6,
        condition_survival = TRUE,
        condition_root_speciation = TRUE,
        condition_marginal = FALSE,
        extinction_approximation = FALSE,
        verbose = FALSE
    ){
    phy$node.label <- NULL
    newick_string <- ape::write.tree(phy)
    phylogeny <- Phylogeny$new(newick_string)

    if (missing(lambda_hat) & missing(mu_hat)){
        ## fit constant-rate model
        if (verbose) print("fitting constant-rate model")
        bd_opt <- stats::optim(
            par = c(0.1, 0.05),
            fn = function(x) phylogeny$bd_likelihood(x[1], x[2], sampling_fraction, tol, FALSE, condition_survival, condition_root_speciation),
            lower = c(0.00001, 0.00001),
            upper = c(10.0, 10.0),
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
        if (verbose) print("fitting shift rate parameter")
        bds_opt <- stats::optim(
            par = c(0.005),
            fn = function(x) phylogeny$bds_likelihood(lambda_hat, mu_hat, x[1], sampling_fraction, sd, num_speciation_classes, num_extinction_classes, tol, FALSE, condition_survival, condition_root_speciation, condition_marginal, extinction_approximation),
            lower = c(1e-8),
            upper = c(0.1),
            method = "Brent",
            control = list(fnscale = -1),
            hessian = FALSE,
         )

        eta <- bds_opt$par[1]
    }

    # store D(t) as a cubic spline 
    if (verbose) print("calculating dense postorder")
    phylogeny$bds_likelihood(lambda_hat, mu_hat, eta, sampling_fraction, sd, num_speciation_classes, num_extinction_classes, tol, TRUE, condition_survival, condition_root_speciation, condition_marginal, extinction_approximation)

    # calculate F(t)
    if (verbose) print("calculating F(t)")
    phylogeny$bds_preorder(lambda_hat, mu_hat, eta, sampling_fraction, sd, num_speciation_classes, num_extinction_classes, tol);

    ## calculate netdiv per branch
    if (verbose) print("calculating netdiv per branch")
    phylogeny$branch_rates(lambda_hat, mu_hat, eta, sampling_fraction, sd, num_speciation_classes, num_extinction_classes);

    ## calculate no. shifts per branch
    if (verbose) print("calculating number of shifts")
    phylogeny$number_of_shifts(lambda_hat, mu_hat, eta, sampling_fraction, sd, num_speciation_classes, num_extinction_classes, tol);

    ## calculate Bayes factors per branch 
    if (verbose) print("calculating Bayes factors")
    phylogeny$bayes_factors(lambda_hat, mu_hat, eta, sampling_fraction, sd, num_speciation_classes, num_extinction_classes, tol);

    ## calculate tip rates
    if (verbose) print("calculating tip rates")
    tip_rates <- phylogeny$tip_rates(lambda_hat, mu_hat, eta, sampling_fraction, sd, num_speciation_classes, num_extinction_classes);

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
                "td" = tree,
                "tip_rates" = tip_rates
                )

    return(res)
}


#' fits the constant-rate birth-death model
#' 
#' @param phy an object of type phylo
#' @param sampling_fraction the probability that each species was sampled in the tree
#' @param tol the local error threshold in the numerical ODE solver (per delta_t time step)
#' @param condition_survival whether or not to condition on the survival of the left and right lineages descending from the root (default TRUE)
#' @param condition_root_speciation whether or not to condition on that there was a speciation event at the root node (default TRUE)
#'
#'
#' @export
fit_cbd <- function(
        phy, 
        sampling_fraction,
        tol = 1e-6,
        condition_survival = TRUE,
        condition_root_speciation = TRUE,
        extinction_approximation = FALSE,
        verbose = FALSE
    ){
    phy$node.label <- NULL
    newick_string <- ape::write.tree(phy)
    phylogeny <- Phylogeny$new(newick_string)

    ## fit constant-rate model
    if (verbose) print("fitting constant-rate model")
    bd_opt <- stats::optim(
        par = c(0.1, 0.05),
        fn = function(x) phylogeny$bd_likelihood(x[1], x[2], sampling_fraction, tol, FALSE, condition_survival, condition_root_speciation),
        lower = c(0.00001, 0.00001),
        upper = c(10.0, 10.0),
        method = "L-BFGS-B",
        control = list(fnscale = -1),
        hessian = FALSE,
     )

    lambda_hat <- bd_opt$par[1]
    mu_hat <- bd_opt$par[2]

    res <- list(
            "lambda" = lambda_hat,
            "mu" = mu_hat
            )

    return(res)
}
