
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
    pub k: usize,
    pub lambdas: f64,
    pub mus: f64,
    pub rho: f64,
}


