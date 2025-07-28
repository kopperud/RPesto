use crate::categories::rate_categories;

#[derive(Debug, Default)]
pub struct ConstantBD {
    pub lambda: f64,
    pub mu: f64,
    pub rho: f64,
}

#[derive(Debug, Default)]
pub struct ShiftBD {
    pub lambda_hat: f64,
    pub mu_hat: f64,
    pub eta: f64,
    pub rho: f64,
    pub k: usize,
    pub lambda: Vec<f64>,
    pub mu: Vec<f64>,
    pub extinction_approximation: bool,
}


        //let model = ShiftBD::new(lambda_hat, mu_hat, eta, rho, sd, n);
impl ShiftBD{
    pub fn new(lambda_hat: f64, mu_hat: f64, eta: f64, rho: f64, sd: f64, n_lambda: usize, n_mu: usize, extinction_approximation: bool) -> ShiftBD{

        let k = n_lambda * n_mu;

        let (lambda, mu) = rate_categories(lambda_hat, mu_hat, sd, n_lambda, n_mu);

        let model = ShiftBD{
            lambda_hat, 
            mu_hat,
            eta,
            rho,
            k,
            lambda,
            mu,
            extinction_approximation,
        };

        return model;
    }
}




