use crate::odesolver::Solve;
use crate::tree::*;
use crate::height::*;
use crate::branch_probability::*;
use crate::models::*;
use crate::spline::*;

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


        self.preorder_po(tree, time, marginal_probability, tol);
    }

    fn preorder_po(&self, node: &mut Box<Node>, time: f64, parent_marginal_probability: Vec<f64>, tol: f64) -> (){
        // should do some checks here if there is None in t and u
        //let t = node.t_dense.clone().unwrap();
        //let u = node.u_dense.clone().unwrap();

        //let num = u[0].len();
        //let k = num / 2;
        let k = self.k;

        let x = ForwardProbability::new(self.lambda.clone(), self.mu.clone(), self.eta, node.extinction_probability.clone().unwrap());
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

        let d_old = node.subtree_probability
            .as_ref()
            .unwrap()
            .interpolate(t0);

        let mut marginal_probability= Vec::new();
        for i in 0..k{
            let p = parent_marginal_probability[i] / d_old[i];
            marginal_probability.push(p);
        }

        let u0 = marginal_probability;
        //let u0 = vec![1.0; k];

        let (times, f) = x.solve_dopri45(u0, t0, t1, true, n_steps_init, tol);

        println!("t0 = {}", t0);
        println!("t1 = {}", t1);

        let z = t0 == t1;
        println!("equal(t0, t1) = {}", z);

        //let z = spline.interpolate(0.0);
        println!("times = {:?}", times);
        println!("f = {:?}", f);
    }
}


