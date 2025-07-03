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
    let (times, probs) = ode.solve_dopri45(u0, t0, t, true, 10, tol, equation);


    let mut m = Vec::new();

    for p in probs{
        m.push(p[0]);
    }

    let res = list!(t = &times, probs = &m);

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

    fn bd_likelihood(&mut self, lambda: f64, mu: f64, rho: f64, tol: f64, store: bool) -> f64{
        let model = ConstantBD{lambda, mu, rho};
        let lnl = model.likelihood(&mut self.tree, tol, store);
        return lnl;
    }

    pub fn bds_likelihood(&mut self, lambda_hat: f64, mu_hat: f64, eta: f64, rho: f64, sd: f64, n: usize, tol: f64, store: bool) -> f64{
        let model = ShiftBD::new(lambda_hat, mu_hat, eta, rho, sd, n);
        let lnl = model.likelihood(&mut self.tree, tol, store);

        return lnl;
    }

    pub fn bds_preorder(&mut self, lambda_hat: f64, mu_hat: f64, eta: f64, rho: f64, sd: f64, n: usize, tol: f64) -> (){
        let model = ShiftBD::new(lambda_hat, mu_hat, eta, rho, sd, n);
        model.preorder(&mut self.tree, tol);
    }

    pub fn branch_rates(&mut self, lambda_hat: f64, mu_hat: f64, eta: f64, rho: f64, sd: f64, n: usize) -> (){
        let model = ShiftBD::new(lambda_hat, mu_hat, eta, rho, sd, n);
        model.net_diversification(&mut self.tree );
        model.speciation(&mut self.tree );
        model.extinction(&mut self.tree );
    }

    pub fn number_of_shifts(&mut self, lambda_hat: f64, mu_hat: f64, eta: f64, rho: f64, sd: f64, n: usize, tol: f64) -> (){
        let model = ShiftBD::new(lambda_hat, mu_hat, eta, rho, sd, n);
        model.number_of_shifts(&mut self.tree, tol);
    }

    pub fn bayes_factors(&mut self, lambda_hat: f64, mu_hat: f64, eta: f64, rho: f64, sd: f64, n: usize, tol: f64) -> (){
        let model = ShiftBD::new(lambda_hat, mu_hat, eta, rho, sd, n);
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
    impl Phylogeny;
}























