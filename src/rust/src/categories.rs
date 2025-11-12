use core::f64;

//use statrs::distribution::{LogNormal,ContinuousCDF};
use crate::erf::erf_inv;


pub fn lognormal_probability(x: f64, location: f64, sigma: f64) -> f64{
    let a = 1.0 / (x * sigma * (2.0 * f64::consts::PI).sqrt());

    let b = (-0.5 * (x.ln() - location).powi(2) / (sigma * sigma)).exp();

    let p = a * b;

    p
}


pub fn lognormal_quantile(x: f64, location: f64, sigma: f64) -> f64{
    let p = lognormal_probability(x, location, sigma);

    let a = (2.0 * sigma * sigma).sqrt();
    let b = erf_inv(2.0 * p - 1.0);

    let q = (location + a * b).exp();

    q
}


/*
mod tests {
    macro_rules! assert_delta {
        ($x:expr, $y:expr, $d:expr) => {
            if !($x - $y < $d || $y - $x < $d) { panic!(); }
        }
    }

    use std::f64;
    use crate::categories::*;

    #[test]
    fn test_lognormal() {
        let x: f64 = 2.0_f64.ln();
        let sigma = 0.8;
        let location = 1.5_f64.ln();
        
        let a = lognormal_quantile(x, location, sigma);
        let b = lognormal_quantile2(x, location, sigma);

        assert_delta!(a, b, 1e-32);
    }
}
*/


/*
pub fn lognormal_quantile(x: f64, location: f64, sigma: f64) -> f64{
    let d = LogNormal::new(location, sigma).unwrap();
    let q = d.inverse_cdf(x);
    return q;
}
*/


pub fn make_quantiles(location: f64, sigma: f64, n: usize) -> Vec<f64>{
    let mut quantiles = Vec::new();

    //let d = LogNormal::new(location, sigma).expect("expected to create lognormal distribution");

    let step = 0.5;

    for i in 1..(n+1){
        let p = ((i as f64) - step)/(n as f64);
        //let q = d.inverse_cdf(p);
        let q = lognormal_quantile(p, location, sigma);
        quantiles.push(q);
    }

    return quantiles;
}

pub fn rate_categories(lambda: f64, mu: f64, sd: f64, n_lambda: usize, n_mu: usize) -> (Vec<f64>, Vec<f64>){
    let lambda_quantiles = make_quantiles(lambda.ln(), sd, n_lambda);
    let mu_quantiles = make_quantiles(mu.ln(), sd, n_mu);

    let mut lambdas = Vec::new();
    let mut mus = Vec::new();

    for i in 0..n_lambda{
        for j in 0..n_mu{
            lambdas.push(lambda_quantiles[i]);
            mus.push(mu_quantiles[j]);
        }
    }

    let res = (lambdas, mus);

    return res;
}

