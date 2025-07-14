use crate::odesolver::*;
use crate::tree::*;
use crate::height::*;
use crate::branch_probability::*;
use crate::models::*;
use crate::spline::*;
use crate::utils::*;
use core::f64;
use std::fmt;

#[derive(Debug, Clone)]
pub struct BayesFactorError;

impl fmt::Display for BayesFactorError {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "could not solve ODE")
    }
}

type Result<T> = std::result::Result<T, BayesFactorError>;



pub struct BayesFactorProblem{
    pub lambda: Vec<f64>,
    pub mu: Vec<f64>,
    pub eta: f64,
    pub k: usize,
    pub subtree_probability: MonotonicCubicSpline,
    pub forward_probability: MonotonicCubicSpline,
}

impl BayesFactorProblem{
    pub fn new(
        lambda: Vec<f64>,
        mu: Vec<f64>,
        eta: f64,
        subtree_probability: MonotonicCubicSpline,
        forward_probability: MonotonicCubicSpline,
    ) -> BayesFactorProblem{
        let k = lambda.len();
        let res = BayesFactorProblem{
            lambda,
            mu,
            eta,
            k,
            subtree_probability,
            forward_probability};
        return res;
    }
}

#[allow(nonstandard_style)]
impl Gradient for BayesFactorProblem{
    fn gradient(&self, dlogX: &mut Vec<f64>, _logX: &Vec<f64>, t: &f64 ) -> (){
        let r = self.eta / (self.k as f64 - 1.0);

        let Dt = self.subtree_probability.interpolate(*t);

        let sum_Dt = Dt.iter().fold(0.0, |acc, x| acc + x);

        for j in 0..self.k{
            dlogX[j] = r * (sum_Dt - Dt[j]) / Dt[j];
        }
    }
}


pub trait BayesFactor<T>{
    fn bayes_factors( &self, tree: &mut Box<Node>, tol: f64) -> ();
    fn bayes_factors_pre( &self, node: &mut Box<Node>, time: f64, tol: f64) -> ();

    fn prob_no_shifts( &self, node: &mut Box<Node>, t0: f64, t1: f64, n_slices: usize) -> Vec<f64>;
    fn prob_no_shifts_iterative( &self, node: &mut Box<Node>, t0: f64, t1: f64, tol: f64) -> Result<Vec<f64>>;
}

impl BayesFactor<BranchProbabilityMultiState> for ShiftBD{
    fn bayes_factors( &self, tree: &mut Box<Node>, tol: f64) -> (){
        let height = treeheight(&tree);
        let time = height;
        
        self.bayes_factors_pre(tree, time, tol);
    }

    fn bayes_factors_pre(&self, node: &mut Box<Node>, time: f64, tol: f64) -> (){
        let t0 = time;
        let t1 = time - node.length;

        let x = self.prob_no_shifts_iterative(node, t0, t1, tol);

        match x{
            Ok(log_prob_no_shifts) => {

                let d_t0 = node.subtree_probability.as_ref().unwrap().interpolate(t0);
                let f_t0 = node.forward_probability.as_ref().unwrap().interpolate(t0);

                let mut marginal_probability_t0 = Vec::new();
                for i in 0..self.k{
                    let p = d_t0[i] * f_t0[i];
                    marginal_probability_t0.push(p);
                }
                normalize(&mut marginal_probability_t0);

                let mut probability_no_shifts = 0.0;
                for i in 0..self.k{
                    let p = marginal_probability_t0[i] * log_prob_no_shifts[i].exp();
                    probability_no_shifts += p;
                }

                let probability_minimum_one_shift = 1.0 - probability_no_shifts;

                let prior_prob_no_shifts = (-(t0 - t1) * self.eta).exp();
                let prior_prob_minimum_one_shift = 1.0 - prior_prob_no_shifts;

                let mut bayes_factor = (probability_minimum_one_shift / prior_prob_minimum_one_shift) / (probability_no_shifts / prior_prob_no_shifts);

                // not sure if this is a good idea ..
                // if bf is NaN then R parses it as a string
                if bayes_factor.is_nan(){
                    bayes_factor = 0.0;
                }

                node.bayes_factor = Some(bayes_factor);  
            },
            Err(_) => {
                println!("could not calculate bayes factor for branch");
                node.bayes_factor = Some(0.0);  
            }
        }

        for child_node in node.children
            .iter_mut()
            {
            self.bayes_factors_pre(child_node, time - node.length, tol);
        }
    }


    fn prob_no_shifts( &self, node: &mut Box<Node>, t0: f64, t1: f64, n_slices: usize) -> Vec<f64>{
        let k = self.k;
        let r = self.eta / (k as f64 - 1.0);

        let delta_t = (t1 - t0) / (n_slices as f64);

        let mut times = Vec::new();
        for i in 0..n_slices{
            //let t = t0 + 0.5*delta_t + delta_t * (i as f64);
            let t = t0 + delta_t * (i as f64);
            times.push(t);
        }

        let mut logX = vec![0.0; k];

        for t in times{
            let d = node.subtree_probability.clone().expect("asd").interpolate(t);
            let sum_d = d.iter().fold(0.0, |acc, x| acc + x);

            for j in 0..k{
                //dlogX[j] = r * (sum_Dt - Dt[j]) / Dt[j];
                logX[j] += delta_t * r * (sum_d / d[j] - 1.0);
            }
        }

        return logX;
    }
    fn prob_no_shifts_iterative( &self, node: &mut Box<Node>, t0: f64, t1: f64, tol: f64) -> Result<Vec<f64>>{

        let bayes_factor_problem = BayesFactorProblem::new(self.lambda.clone(), self.mu.clone(), self.eta, node.subtree_probability.clone().unwrap(), node.forward_probability.clone().unwrap());

        let mut n= 8;
        let mut res = vec![0.0; self.k];
        let mut mean_err = f64::INFINITY;
        let tol1 = tol * 100.0;

        while mean_err > tol1{
            n *= 2;

            let u0 = vec![0.0; self.k];
            let (_, logX1) = bayes_factor_problem.solve_rk4(u0, t0, t1, false, n-1).expect("asd");

            let u0 = vec![0.0; self.k];
            let (_, logX2) = bayes_factor_problem.solve_rk4(u0, t0, t1, false, n).expect("asd");

            let low_order_estimate = &logX1[0];
            let high_order_estimate = &logX2[0];

            let mut sq_err = 0.0;
            for (x, y) in high_order_estimate.iter().zip(low_order_estimate){
                sq_err += (x.exp() - y.exp()).powi(2);
            }
            mean_err = sq_err.sqrt();
            //println!("n = {}, mean_err = {}", n, mean_err);

            if mean_err < tol1{
                res = high_order_estimate.clone();
            }

            if n > 2049{
                return Err(BayesFactorError);
            }
        }

        return Ok(res);
    }
}
 
