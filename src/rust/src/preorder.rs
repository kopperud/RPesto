use crate::odesolver::{EquationType, Solve};
use crate::tree::*;
use crate::height::*;
use crate::branch_probability::*;
use crate::models::*;
use crate::spline::*;
use crate::utils::*;

use rayon::iter::ParallelIterator;
use rayon::prelude::*;

// the likelihood trait 
pub trait Preorder<T>{
    fn preorder( &self, tree: &mut Box<Node>, tol: f64) -> ();
    fn preorder_po( &self, node: &mut Box<Node>, time: f64, parent_marginal_probability: Vec<f64>, tol: f64) -> ();
}

impl Preorder<BranchProbabilityMultiState> for ShiftBD{
    fn preorder( &self, tree: &mut Box<Node>, tol: f64) -> (){
        let height = treeheight(&tree);
        let time = height;

        let k = self.lambda.len();
        
        let root_prior= vec![1.0 / (k as f64); k];

        //let last_index = u_dense.len();
        //let root_d = &u_dense[last_index-1];
        let root_d = tree
            .subtree_probability
            .as_ref()
            .unwrap()
            .interpolate(time);

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


        self.preorder_po(tree, time, marginal_probability, tol);
    }

    fn preorder_po(&self, node: &mut Box<Node>, time: f64, parent_marginal_probability: Vec<f64>, tol: f64) -> (){
        let k = self.k;

        let x = ForwardProbability::new(self.lambda.clone(), self.mu.clone(), self.eta, node.extinction_probability.clone().unwrap());

        let t0 = time;
        let t1 = time - node.length;
        let n_steps_init = 5;

        let d_old = node.subtree_probability
            .as_ref()
            .unwrap()
            .interpolate(t0);

        let mut marginal_probability= Vec::new();
        for i in 0..k{
            let p = parent_marginal_probability[i] / d_old[i];
            marginal_probability.push(p);
        }
        normalize(&mut marginal_probability);

        let u0 = marginal_probability;

        let (times, sol) = x.solve_dopri45(u0, t0, t1, true, n_steps_init, tol, EquationType::ProbabilityDensity);


        let sol_last = sol.last().cloned().unwrap();
        let forward_probability = MonotonicCubicSpline::new(times.clone(), sol, self.k, false);

        node.forward_probability = Some(forward_probability);
        // marginal probability at t1, youngest on branch

        let d1 = node.subtree_probability.as_ref().unwrap().interpolate(t1);


        let mut marginal_probability_young = Vec::new();
        for i in 0..self.k{
            let p = d1[i] * sol_last[i];
            marginal_probability_young.push(p);
        }

        // normalize such that f sums to 1
        normalize(&mut marginal_probability_young);

        for child_node in node.children.iter_mut(){
            self.preorder_po(child_node, time - node.length, marginal_probability_young.clone(), tol);
        }
    }
}


