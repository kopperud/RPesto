use extendr_api::prelude::*;
use microbench::{self, Options};


use crate::spline::*;
use crate::tokenizer::tokenize;
use crate::parser::*;
use crate::extinction::*;
use crate::branch_probability::*;
use crate::utils::*;
use crate::odesolver::*;
use crate::likelihood::*;

pub mod spline;
pub mod tokenizer;
pub mod tree;
pub mod parser;
pub mod extinction;
pub mod likelihood;
pub mod branch_probability;
pub mod utils;
pub mod odesolver;
pub mod height;

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

/// @export
#[extendr]
fn likelihood(lambda: f64, mu: f64, phy: String, tol: f64) -> f64{
    let lnl = likelihood_cbdp(lambda, mu, phy, tol);

    return lnl;
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

    //let n_times = 100;
    //let times2 = sequence(0.0, t, n_times);
    for p in probs{
        m.push(p[0]);
    //for t in times2.iter(){
        //m.push(spline.interpolate(*t)[0]);
    }

    let res = list!(t = &times, probs = &m);
    //microbench::bench(&options, "making an R list", || {list!(t = &times, probs = &m)});

    return res;
}

/// @export
#[extendr]
fn branch_probability(lambda: f64, mu: f64, t: f64, tol : f64) -> extendr_api::List{

    let height = 60.0;
    let ode = BranchProbability::new(lambda, mu, height, tol);
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

const P: f64 = 0.47047;
const A1: f64 = 0.3480242;
const A2: f64 = -0.0958798;
const A3: f64 = 0.7478556;


fn error_function(x: f64) -> f64{
    let t = 1.0 / (1.0 + P*x);

    let res = 1.0 - (A1*t + A2*t*t + A3*t*t*t) * f64::exp(-(x*x));

    return res;
}

/// @export
#[extendr]
pub fn lognormal_quantile(x: f64, scale: f64, sigma: f64) -> f64{
    let res = 0.5 * (1.0 + error_function((x.ln() - scale)/(sigma * f64::sqrt(2.0))));
    return res;
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
    fn likelihood; 
    fn lognormal_quantile; 
}
