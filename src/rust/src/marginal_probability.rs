use crate::odesolver::Solve;
use crate::tree::*;
use crate::height::*;
use crate::branch_probability::*;
use crate::models::*;
use crate::spline::*;
use crate::utils::*;

use rayon::iter::ParallelIterator;
use rayon::prelude::*;

// the likelihood trait 
pub trait Marginal<T>{
    fn preorder( &self, tree: &mut Box<Node> ) -> ();
    fn preorder_po( &self, node: &mut Box<Node>, time: f64 ) -> ();
}

impl Marginal<BranchProbabilityMultiState> for ShiftBD{
    fn preorder( &self, tree: &mut Box<Node> ) -> (){
        let height = treeheight(&tree);
        let time = height;

        let k = self.lambda.len();
        
        let root_prior= vec![1.0 / (k as f64); k];
        let u_dense = tree.u_dense.clone().unwrap().clone();

        let last_index = u_dense.len();
        let root_d = &u_dense[last_index-1];

        let mut marginal_probability= Vec::new();
        for i in 0..k{
            let p = root_d[i] * root_prior[i];
            marginal_probability.push(p);
        }

        // normalize
        let s = marginal_probability.iter().fold(0.0, |acc, x| acc + x);
        for i in 0..k{
            marginal_probability[i] = marginal_probability[i] / s;
        }


        self.preorder_po(tree, time);
    }

    fn preorder_po(&self, node: &mut Box<Node>, time: f64) -> (){
        // should do some checks here if there is None in t and u

        let k = self.k;

        let x = ForwardProbability::new(self.lambda.clone(), self.mu.clone(), self.eta, node.extinction_probability.clone().unwrap());

        println!("asd6");
        /*
        let t0 = time;
        let t1 = time - node.length;

        let u_dense = node.u_dense.clone().unwrap().clone();
        let last_index = u_dense.len();
        let d_old = &u_dense[last_index-1];
        */
        let t0 = time;
        let t1 = time - node.length;
        let n_steps_init = 5;

        let d = node.subtree_probability
            .as_ref()
            .unwrap()
            .interpolate(t0);

        let f = node.subtree_probability
            .as_ref()
            .unwrap()
            .interpolate(t0);

        let mut si = Vec::new();
        for (di, fi) in d.iter().zip(f){
            si.push(di * fi);
        }

        // normalize such that f sums to 1
        normalize(&mut marginal_probability_young);

        println!("asd3");

        for child_node in node.children.iter_mut(){
            println!("asd4");
            self.preorder_po(child_node, time - node.length);
        }

        
    }
}


