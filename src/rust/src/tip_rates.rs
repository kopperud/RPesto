use crate::tree::*;
use crate::height::*;
use crate::branch_probability::*;
use crate::models::*;
use crate::utils::*;
use std::collections::HashMap;

use core::f64;
use std::fmt;

#[derive(Debug, Clone)]
pub struct TipRatesError;

impl fmt::Display for TipRatesError {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "could not solve eq")
    }
}

pub trait TipRates<T>{
    fn tip_rates( &self, tree: &mut Box<Node>) -> HashMap<String, f64>;
    fn tip_rates_pre( &self, node: &mut Box<Node>, map: &mut HashMap<String, f64>, time: f64) -> ();
}

impl TipRates<BranchProbabilityMultiState> for ShiftBD{
    fn tip_rates( &self, tree: &mut Box<Node>) -> HashMap<String,f64>{
        let height = treeheight(&tree);
        let time = height;

        let mut map: HashMap<String, f64> = HashMap::new();
        //map.entry("homo sapiens".to_string()).or_insert(3.0);
        
        self.tip_rates_pre(tree, &mut map, time);

        return map;
    }

    fn tip_rates_pre(&self, node: &mut Box<Node>, map: &mut HashMap<String, f64>, time: f64) -> (){
        //let t0 = time;
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

            for i in 0..self.k{
                let p = marginal_probability_t1[i];
                let ri = self.lambda[i] - self.mu[i];
                tip_netdiv += p * ri;
            }

            // not sure if this is a good idea ..
            // if the value is NaN then R parses it as a string
            if tip_netdiv.is_nan(){
                tip_netdiv = -999.0;
            }

            //node.tip_netdiv = Some(tip_netdiv);  
            map.entry(node.label.clone()).or_insert(tip_netdiv);
        }else{
            //node.tip_netdiv = Some(0.0);
        }

        for child_node in node.children
            .iter_mut()
            {
            self.tip_rates_pre(child_node, map, time - node.length);
        }
    }
}

 
