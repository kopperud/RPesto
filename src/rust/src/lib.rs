#![allow(non_snake_case)]


use extendr_api::prelude::*;

use crate::parser::*;
use crate::extinction::*;
use crate::odesolver::*;
use crate::likelihood::*;
use crate::preorder::*;
use crate::models::*;
use crate::tree::Node;
use crate::branchrates::*;
use crate::writenewick::*;
use crate::number_of_shifts::*;
use crate::bayes_factor::*;
use crate::conditioning::*;

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
pub mod writenewick;
pub mod number_of_shifts;
pub mod bayes_factor;
pub mod conditioning;

/// Return string `"Hello world!"` to R.
/// @export
#[extendr]
fn hello_world() -> &'static str {
    "Hello world!"
}


#[extendr]
fn extinction_probability(lambda: f64, mu: f64, t: f64, tol: f64) -> extendr_api::List{
    let ode = Extinction{mu, lambda};

    let u0 = vec![0.0];
    let t0 = 0.0;

    let equation = EquationType::Probability;
    let (times, probs) = ode.solve_dopri45(u0, t0, t, true, 10, tol, equation).expect("asd");


    let mut m = Vec::new();

    for p in probs{
        m.push(p[0]);
    }

    let res = list!(t = &times, probs = &m);

    return res;
}

use crate::branch_probability::*;

/// @export
#[extendr]
fn branch_probability_bds(lambda_hat: f64, mu_hat: f64, eta: f64, rho: f64, sd: f64, n_lambda: usize, n_mu: usize, t: f64, tol: f64, extinction_approximation: bool) -> extendr_api::List{
//fn extinction_probability_bds(lambda_hat: f64, mu_hat: f64, eta: f64, sd: f64, sampling_probability: f64, t: f64, tol: f64, extinction_approximation: bool) -> extendr_api::List{

    let model = ShiftBD::new(lambda_hat, mu_hat, eta, rho, sd, n_lambda, n_mu, extinction_approximation);

    let k = model.mu.len();

    //let ode = ExtinctionMultiState{mu: model.mu, lambda: model.lambda, eta, extinction_approximation};
    let ode = BranchProbabilityMultiState::new(model.lambda.clone(), model.mu.clone(), eta, extinction_approximation);

    let mut u0 = vec![0.0; k*2];

    for i in 0..k{
        u0[i] = 1.0 - rho;
        u0[k+i] = rho;
    }
        
    let dense = true; // solve ODE with dense output
    let n_steps_init = 4;

    let t0 = 0.0;
    let t1 = t;

    let equation = EquationType::ProbabilityDensity;
    let (times, sol) = ode.solve_dopri45(u0, t0, t1, dense, n_steps_init, tol, equation).expect("could not calculate branch probabiltiy E(t) and D(t) in likelihood function");


    let n_times = times.len();

    let m = RMatrix::new_matrix(n_times, 2*k, |r, c| sol[r][c]);

    let res = list!(t = &times, probs = m);

    return res;
}

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

    fn bd_likelihood(&mut self, lambda: f64, mu: f64, rho: f64, tol: f64, store: bool, condition_survival: bool, condition_root_speciation: bool) -> f64{
        let model = ConstantBD{lambda, mu, rho};

        let mut conditions: Vec<Condition> = Vec::new();
        if condition_survival{
            conditions.push(Condition::Survival);
        }
        if condition_root_speciation{
            conditions.push(Condition::RootSpeciation);
        }

        // this does not matter for cbdp
        let condition_marginal = false; 

        let lnl = model.likelihood(&mut self.tree, conditions, tol, condition_marginal, store);
        return lnl;
    }

    pub fn bds_likelihood(&mut self, lambda_hat: f64, mu_hat: f64, eta: f64, rho: f64, sd: f64, n_lambda: usize, n_mu: usize, tol: f64, store: bool, condition_survival: bool, condition_root_speciation: bool, condition_marginal: bool, extinction_approximation: bool) -> f64{
        let model = ShiftBD::new(lambda_hat, mu_hat, eta, rho, sd, n_lambda, n_mu, extinction_approximation);

        let mut conditions: Vec<Condition> = Vec::new();

        if condition_survival{
            conditions.push(Condition::Survival);
        }
        if condition_root_speciation{
            conditions.push(Condition::RootSpeciation);
        }

        let lnl = model.likelihood(&mut self.tree, conditions, tol, condition_marginal, store);

        return lnl;
    }

    pub fn bds_preorder(&mut self, lambda_hat: f64, mu_hat: f64, eta: f64, rho: f64, sd: f64, n_lambda: usize, n_mu: usize, tol: f64 ) -> (){
        let model = ShiftBD::new(lambda_hat, mu_hat, eta, rho, sd, n_lambda, n_mu, false);
        model.preorder(&mut self.tree, tol);
    }

    pub fn branch_rates(&mut self, lambda_hat: f64, mu_hat: f64, eta: f64, rho: f64, sd: f64, n_lambda: usize, n_mu: usize) -> (){
        let model = ShiftBD::new(lambda_hat, mu_hat, eta, rho, sd, n_lambda, n_mu, false);
        model.net_diversification(&mut self.tree );
        model.speciation(&mut self.tree );
        model.extinction(&mut self.tree );
    }

    pub fn number_of_shifts(&mut self, lambda_hat: f64, mu_hat: f64, eta: f64, rho: f64, sd: f64, n_lambda: usize, n_mu: usize, tol: f64) -> (){
        let model = ShiftBD::new(lambda_hat, mu_hat, eta, rho, sd, n_lambda, n_mu, false);
        model.number_of_shifts(&mut self.tree, tol);
    }

    pub fn bayes_factors(&mut self, lambda_hat: f64, mu_hat: f64, eta: f64, rho: f64, sd: f64, n_lambda: usize, n_mu: usize, tol: f64) -> (){
        let model = ShiftBD::new(lambda_hat, mu_hat, eta, rho, sd, n_lambda, n_mu, false);
        model.bayes_factors(&mut self.tree, tol);
    }

    pub fn write_newick(&mut self) -> String{
        let newick = self.tree.writenewick();
        return newick;
    }

}


// Macro to generate exports.
// This ensures exported functions are registered with R.
// See corresponding C code in `entrypoint.c`.
extendr_module! {
    mod RPesto;
    fn hello_world;
    fn extinction_probability; 
    fn branch_probability_bds; 
    impl Phylogeny;
}























