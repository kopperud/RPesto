use crate::tree::*;
use crate::height::*;
use crate::branch_probability::*;
use crate::models::*;
use crate::utils::*;
use core::f64;
use std::fmt;

use extendr_api::prelude::*;

#[derive(Debug,IntoDataFrameRow)]
pub struct TipRateRow {
  species: String,
  netdiv: f64,
  lambda: f64,
  mu: f64,
  relext: f64,
}

#[derive(Debug, Clone)]
pub struct TipRatesError;

impl fmt::Display for TipRatesError {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "could not solve eq")
    }
}

fn reassign_nan(x: &mut f64) -> (){
    if x.is_nan(){
        *x = 999.0;
    }
}

pub trait TipRates<T>{
    fn tip_rates( &self, tree: &mut Box<Node>) -> Vec<TipRateRow>;
    fn tip_rates_pre( &self, node: &mut Box<Node>, items: &mut Vec<TipRateRow>, time: f64) -> ();
}

impl TipRates<BranchProbabilityMultiState> for ShiftBD{
    fn tip_rates( &self, tree: &mut Box<Node>) -> Vec<TipRateRow>{
        let height = treeheight(&tree);
        let time = height;

        let mut items: Vec<TipRateRow> = Vec::new();
        
        self.tip_rates_pre(tree, &mut items, time);

        return items;
    }

    fn tip_rates_pre(&self, node: &mut Box<Node>, items: &mut Vec<TipRateRow>, time: f64) -> (){
        let t1 = time - node.length;

        let n_children = node.children.len();
        let is_tip = n_children == 0;

        if is_tip{
            let d_t1 = node.subtree_probability.as_ref().unwrap().interpolate(t1);
            let f_t1 = node.forward_probability.as_ref().unwrap().interpolate(t1);

            let mut marginal_probability_t1 = Vec::new();
            for i in 0..self.k{
                let p = d_t1[i] * f_t1[i];
                marginal_probability_t1.push(p);
            }
            normalize(&mut marginal_probability_t1);

            let mut tip_netdiv = 0.0;
            let mut tip_lambda = 0.0;
            let mut tip_mu = 0.0;
            let mut tip_relext = 0.0;

            for i in 0..self.k{
                let p = marginal_probability_t1[i];

                let ri = self.lambda[i] - self.mu[i];
                tip_netdiv += p * ri;

                tip_lambda += p * self.lambda[i];

                tip_mu += p * self.mu[i];
                
                let relext_i = self.mu[i] / self.lambda[i];
                tip_relext += p * relext_i;
            }

            // not sure if this is a good idea ..
            // if the value is NaN then R parses it as a string
            reassign_nan(&mut tip_netdiv);
            reassign_nan(&mut tip_lambda);
            reassign_nan(&mut tip_mu);
            reassign_nan(&mut tip_relext);

            let item = TipRateRow{
                    species: node.label.clone(),
                    netdiv: tip_netdiv,
                    lambda: tip_lambda,
                    mu: tip_mu,
                    relext: tip_relext,
                };


            items.push(item);

        }

        for child_node in node.children
            .iter_mut()
            {
            self.tip_rates_pre(child_node, items, time - node.length);
        }
    }
}

 
