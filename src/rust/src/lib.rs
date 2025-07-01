use extendr_api::prelude::*;
use microbench::{self, Options};

use crate::tokenizer::tokenize;
use crate::parser::*;
use crate::extinction::*;
use crate::branch_probability::*;
use crate::utils::*;
use crate::odesolver::*;
use crate::likelihood::*;
use crate::preorder::*;
use crate::models::*;
use crate::categories::*;
use crate::marginal_probability::*;
use crate::tree::Node;
use crate::branchrates::*;

pub mod spline;
pub mod marginal_probability;
pub mod categories;
pub mod tokenizer;
pub mod tree;
pub mod parser;
pub mod extinction;
pub mod likelihood;
pub mod preorder;
pub mod branch_probability;
pub mod utils;
pub mod odesolver;
pub mod branchrates;
pub mod height;
pub mod models;

/// Return string `"Hello world!"` to R.
/// @export
#[extendr]
fn hello_world() -> &'static str {
    "Hello world!"
}

/// @export
#[extendr]
fn print_me_a_tree(s: String) -> (){
    let tree = parse_tree(s);

    println!("{:?}", tree);
}





fn sequence(from: f64, to: f64, num: usize) -> Vec<f64> {
    let mut v = Vec::new();
    
    let delta = (to - from) / (num as f64);

    let mut val = from;
    for _ in 0..num{
        v.push(val);
        val += delta;
    }
    v.push(val);

    return v;
}

/// @export
#[extendr]
fn extinction_probability(lambda: f64, mu: f64, t: f64, tol: f64) -> extendr_api::List{
    let ode = Extinction{mu, lambda};

    let u0 = vec![0.0];
    let t0 = 0.0;

    //let options = Options::default();
    //microbench::bench(&options, "solving ODE", || ode.solve_dopri45(u0.clone(), t0, t, false, 10, tol) );
    let (times, probs) = ode.solve_dopri45(u0, t0, t, true, 10, tol);
    //let n_rows = probs.len();
    //let m = RMatrix::new_matrix(n_rows, 1, |r, c| probs[r][c]);
    
    //let spline = MonotonicCubicSpline::new(times, probs, 1);


    let mut m = Vec::new();

    for p in probs{
        m.push(p[0]);
    }

    let res = list!(t = &times, probs = &m);
    //microbench::bench(&options, "making an R list", || {list!(t = &times, probs = &m)});

    return res;
}

/// @export
#[extendr]
fn branch_probability(lambda: f64, mu: f64, t: f64, tol : f64) -> extendr_api::List{

    //let height = 60.0;
    let ode = BranchProbability::new(lambda, mu);
    let u0 = vec![1.0];
    let t0 = 0.0;

    let (times, probs) = ode.solve_dopri45(u0, t0, t, true, 10, tol);

    let mut m = Vec::new();
    for p in probs{
        m.push(p[0]);
    }

    let res = list!(t = &times, probs = &m);
    //microbench::bench(&options, "making an R list", || {list!(t = &times, probs = &m)});

    return res;
}

/// @export
#[extendr]
fn branch_probability2(lambda_hat: f64, mu_hat: f64, eta: f64, sd: f64, n: usize, t: f64, tol : f64) -> extendr_api::List{

    let height = 100.0;
    let (lambda, mu) = rate_categories(lambda_hat, mu_hat, sd, n);
    let rho = 1.0;
    let k = n*n;

    //return list!(lambda = lambda, mu = mu);

    let ode = BranchProbabilityMultiState::new(lambda.clone(), mu.clone(), eta);

    let mut u0 = vec![1.0 - rho; k];
    u0.extend(vec![rho; k]);

    let t0 = 0.0;

    let (times, probs) = ode.solve_dopri45(u0, t0, t, true, 4, tol);

    let nrows = probs.len();
    let ncols = probs[0].len();

    let m = extendr_api::matrix::RMatrix::new_matrix(nrows, ncols, |r,c| probs[r][c]);

    let res = list!(t = &times, probs = m);

    return res;
}

/*
/// @export
#[extendr]
fn foo(lambda_hat: f64, mu_hat: f64, eta: f64, rho: f64, sd: f64, n: usize, height: f64, tol: f64) -> (){

    let (lambda, mu) = rate_categories(lambda_hat, mu_hat, sd, n);

    let ode = ExtinctionMultiState{
        lambda: lambda.clone(),
        mu: mu.clone(),
        eta
    };

    let k = n*n;
    let u0 = vec![1.0 - rho; k];
    let t0 = 0.0;
    let t1 = height;
    let dense = true;
    let n_steps_init = 10;

    let options = Options::default();

    let (times, sol) = ode.solve_dopri45(u0.clone(), t0, t1, dense, n_steps_init, tol);
    microbench::bench(&options, "solving ODE", || ode.solve_dopri45(u0.clone(), t0, t1, false, 10, tol) );
    
    let extinction_probability = MonotonicCubicSpline::new(times, sol, k);

    microbench::bench(&options, "interpolate at t=1.0", || extinction_probability.interpolate(1.0) );
    microbench::bench(&options, "interpolate at t=60.0", || extinction_probability.interpolate(60.0) );

    let mut du = Vec::new();
    let mut u = Vec::new();
    for _ in 0..k{
        du.push(0.0);
        u.push(0.0);
    }
    microbench::bench(&options, "gradient for extinction SSE", || ode.gradient(&mut du, &u, &0.0) );
}
*/


#[derive(Debug)]
#[extendr]
/// @export
struct Phylogeny {
    pub tree: Box<Node>,
}

#[extendr]
impl Phylogeny {
    fn new(newick: String) -> Self {
        Self {
            tree: parse_tree(newick)
        }
    }

    fn print(&self) -> () {
        println!("{:?}", self.tree); 
    }

    fn bd_likelihood(&mut self, lambda: f64, mu: f64, rho: f64, tol: f64, store: bool) -> f64{
        let model = ConstantBD{lambda, mu, rho};
        let lnl = model.likelihood(&mut self.tree, tol, store);
        return lnl;
    }

    pub fn bds_likelihood(&mut self, lambda_hat: f64, mu_hat: f64, eta: f64, rho: f64, sd: f64, n: usize, tol: f64, store: bool) -> f64{
        let model = ShiftBD::new(lambda_hat, mu_hat, eta, rho, sd, n);
        let lnl = model.likelihood(&mut self.tree, tol, store);

        //let options = Options::default();
        //microbench::bench(&options, "calculating likelihood", || model.likelihood(&tree, tol) );
        return lnl;
    }

    pub fn bds_preorder(&mut self, lambda_hat: f64, mu_hat: f64, eta: f64, rho: f64, sd: f64, n: usize, tol: f64) -> (){
        let model = ShiftBD::new(lambda_hat, mu_hat, eta, rho, sd, n);
        model.preorder(&mut self.tree, tol);
    }

     pub fn branch_rates(&mut self, lambda_hat: f64, mu_hat: f64, eta: f64, rho: f64, sd: f64, n: usize) -> (){
        let model = ShiftBD::new(lambda_hat, mu_hat, eta, rho, sd, n);
        model.net_diversification(&mut self.tree );
    }

    //pub fn marginal_probabilities(&mut self) -> (){
        //model.marginal_probability(&mut self.tree);
    //}
}


// Macro to generate exports.
// This ensures exported functions are registered with R.
// See corresponding C code in `entrypoint.c`.
extendr_module! {
    mod RPesto;
    fn hello_world;
    fn print_me_a_tree;
    fn extinction_probability; 
    fn branch_probability; 
    fn branch_probability2;
    impl Phylogeny;
}

//extendr_mo!(






















