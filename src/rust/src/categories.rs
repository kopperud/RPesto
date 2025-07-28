use statrs::distribution::{LogNormal,ContinuousCDF};

pub fn lognormal_quantile(x: f64, location: f64, sigma: f64) -> f64{
    let d = LogNormal::new(location, sigma).unwrap();
    let q = d.inverse_cdf(x);
    return q;
}


pub fn make_quantiles(location: f64, sigma: f64, n: usize) -> Vec<f64>{
    let mut quantiles = Vec::new();

    let d = LogNormal::new(location, sigma).expect("expected to create lognormal distribution");

    let step = 0.5;

    for i in 1..(n+1){
        let p = ((i as f64) - step)/(n as f64);
        let q = d.inverse_cdf(p);
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

