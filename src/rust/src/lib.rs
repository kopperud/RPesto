use extendr_api::prelude::*;


use crate::spline::*;
use crate::tokenizer::tokenize;
use crate::parser::*;
use crate::extinction::*;
use crate::utils::*;
use crate::odesolver::*;

pub mod spline;
pub mod tokenizer;
pub mod tree;
pub mod parser;
pub mod extinction;
pub mod utils;
pub mod odesolver;

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

fn likelihood(lambda: f64, mu: f64, phy: String) -> f64{
    let tree = parse_tree(phy);

    let lnl = 0.0;

    return lnl;
}

// Macro to generate exports.
// This ensures exported functions are registered with R.
// See corresponding C code in `entrypoint.c`.
extendr_module! {
    mod RPesto;
    fn hello_world;
    fn print_me_a_tree;
}
