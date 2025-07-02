use crate::odesolver::*;
use crate::tree::*;
use crate::height::*;
use crate::branch_probability::*;
use crate::models::*;
use crate::spline::*;
use crate::utils::*;

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
}

impl BayesFactor<BranchProbabilityMultiState> for ShiftBD{
    fn bayes_factors( &self, tree: &mut Box<Node>, tol: f64) -> (){
        let height = treeheight(&tree);
        let time = height;
        
        self.bayes_factors_pre(tree, time, tol);
    }

    fn bayes_factors_pre(&self, node: &mut Box<Node>, time: f64, tol: f64) -> (){

        let bayes_factor_problem = BayesFactorProblem::new(self.lambda.clone(), self.mu.clone(), self.eta, node.subtree_probability.clone().unwrap(), node.forward_probability.clone().unwrap());

        let n_steps_init = 5;

        let u0 = vec![0.0; self.k];

        let t0 = time;
        let t1 = time - node.length;

        let (_, log_prob_no_shifts) = bayes_factor_problem.solve_dopri45(u0, t0, t1, false, n_steps_init, tol, EquationType::LogProbability);


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
            let p = marginal_probability_t0[i] * log_prob_no_shifts[0][i].exp();
            probability_no_shifts += p;
        }

        let probability_minimum_one_shift = 1.0 - probability_no_shifts;

        let prior_prob_no_shifts = (-(t0 - t1) * self.eta).exp();
        let prior_prob_minimum_one_shift = 1.0 - prior_prob_no_shifts;

        let mut bayes_factor = (probability_minimum_one_shift / prior_prob_minimum_one_shift) / (probability_no_shifts / prior_prob_no_shifts);

        // not sure if this is a good idea ..
        // if bf is NaN then R parses it as a string
        if bayes_factor.is_nan(){
            bayes_factor = 1.0;
        }

        node.bayes_factor = Some(bayes_factor);  


        for child_node in node.children.iter_mut(){
            self.bayes_factors_pre(child_node, time - node.length, tol);
        }
    }
}
 
