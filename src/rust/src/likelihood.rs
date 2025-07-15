use crate::extinction::*;
use crate::odesolver::{EquationType, Solve};
use crate::tree::*;
use crate::height::*;
use crate::branch_probability::*;
use crate::models::*;
use crate::spline::*;
use crate::conditioning::*;

use rayon::iter::ParallelIterator;
use rayon::prelude::*;

// the likelihood trait 
pub trait Likelihood<T>{
    fn likelihood( &self, tree: &mut Box<Node>, conditions: Vec<Condition>, tol: f64, store: bool) -> f64;
    fn likelihood_po( &self, node: &mut Box<Node>, ode: &T, time: f64, tol: f64, store: bool) -> (Vec<f64>, f64);
}

// the likelihood implementation for constant-rate birth death model
impl Likelihood<BranchProbability> for ConstantBD{
    fn likelihood(&self, tree: &mut Box<Node>, conditions: Vec<Condition>, tol: f64, store: bool) -> f64{
        let height = treeheight(&tree);

        let time = height;

        let ode = BranchProbability::new(self.lambda, self.mu);

        let (u, sf) = self.likelihood_po(tree, &ode, time, tol, store);
        let mut p = u[0];

        if conditions.contains(&Condition::Survival){
            let ode = Extinction{lambda: self.lambda, mu: self.mu};
            let u0 = vec![1.0 - self.rho];
            let (_, w) = ode.solve_dopri45(u0, 0.0, time, false, 5, tol, EquationType::Probability).expect("could not calculate extinction probability in cbdp");
            let e = w[0][0];

            p = p / ((1.0 - e) * (1.0 - e));
        }

        if conditions.contains(&Condition::RootSpeciation){
            p = p / self.lambda;
        }

        let lnl = p.ln() + sf;
        return lnl;
    }

    fn likelihood_po(&self, node: &mut Box<Node>, ode: &BranchProbability, time: f64, tol: f64, store: bool) -> (Vec<f64>, f64){

        let mut u = vec![0.0, 1.0];
        let mut log_sf = 0.0;

        let child_time = time - node.length;

        let r: Vec<(Vec<f64>, f64)> = node
            .children
            .iter_mut()
            .par_bridge()
            //.children.par_iter()
            .map(|child| {
                let x = self.likelihood_po(child, ode, child_time, tol, store);
                return x;
            })
        .collect();

        for (child_u, child_sf) in r.iter(){
            u[0] = child_u[0];
            u[1] *= child_u[1];
            log_sf += child_sf;
        }

        let n_children = node.children.len();
        if n_children > 1{
            for _ in 0..(n_children-1){
                u[1] *= self.lambda;
            }
        }else{
            u[0] = 1.0 - self.rho;
            u[1] = self.rho;
        }
            
        let u0 = u;
        let dense = false;
        let n_steps_init = 4;

        let t0 = child_time;
        let t1 = time;

        let (_times, sol) = ode.solve_dopri45(u0, t0, t1, dense, n_steps_init, tol, EquationType::ProbabilityDensity).expect("could not calculate branch probability in likelihood for cbdp");
    
        let mut p = sol[0].clone();


        log_sf += p[1].ln();
        p[1] = 1.0;
            
        return (p, log_sf);
    }
}



// the likelihood implementation for birth-death-shift model
impl Likelihood<BranchProbabilityMultiState> for ShiftBD{
    fn likelihood(&self, tree: &mut Box<Node>, conditions: Vec<Condition>, tol: f64, store: bool) -> f64{
        let height = treeheight(&tree);

        let time = height;

        let ode = BranchProbabilityMultiState::new(self.lambda.clone(), self.mu.clone(), self.eta, self.extinction_approximation);

        let (p, sf) = self.likelihood_po(tree, &ode, time, tol, store);

        let mut e: Vec<f64> = Vec::new();

        if conditions.contains(&Condition::Survival) || conditions.contains(&Condition::MarginalSurvival){
            let ode = ExtinctionMultiState{
                lambda: self.lambda.clone(),
                mu: self.mu.clone(),
                eta: self.eta, 
                extinction_approximation: self.extinction_approximation,
            };
            let u0 = vec![1.0 - self.rho; self.k];
            let (_, w) = ode.solve_dopri45(u0, 0.0, time, false, 5, tol, EquationType::Probability).expect("could not calculate extinction probability");
            e.extend(w.last().unwrap());
        }


        let root_prior = vec![1.0 / (self.k as f64); self.k];

        let mut lnl = 0.0;
        let mut pr = 0.0;
        for i in 0..self.k{
            let mut x = root_prior[i] * p[self.k+i];

            if conditions.contains(&Condition::Survival){
                x = x / ((1.0 - e[i]) * (1.0 - e[i]));
            }


            if conditions.contains(&Condition::RootSpeciation){
                x = x / self.lambda[i];
            }
            pr += x;
        }

        lnl += pr.ln() + sf;

        if conditions.contains(&Condition::MarginalSurvival){
            let mut marginal_extinction_prob = 0.0;
            for i in 0..self.k{
                marginal_extinction_prob += root_prior[i] * e[i];
            }

            let marginal_survival_prob = 1.0 - marginal_extinction_prob;

            let log_msp = marginal_survival_prob.ln();
            lnl -= 2.0 * log_msp;
        }

        return lnl;
    }

    fn likelihood_po(&self, node: &mut Box<Node>, ode: &BranchProbabilityMultiState, time: f64, tol: f64, store: bool) -> (Vec<f64>, f64){

        let mut u = vec![0.0; self.k];
        let b = vec![1.0; self.k];
        u.extend(b);

        let mut log_sf = 0.0;

        let child_time = time - node.length;

        let r: Vec<(Vec<f64>, f64)> = node
            .children
            .iter_mut()
            .par_bridge()
            .map(|child| {
                let x = self.likelihood_po(child, ode, child_time, tol, store);
                return x;
            })
        .collect();

        for (child_u, child_sf) in r.iter(){

            for i in 0..self.k{
                u[i] = child_u[i];
                u[self.k+i] *= child_u[self.k+i];
            }
            log_sf += child_sf;
        }


        let n_children = node.children.len();
        if n_children > 1{
            for _ in 0..(n_children-1){
                for i in 0..self.k{
                    u[self.k+i] *= self.lambda[i];
                }
            }
        }else{
            for i in 0..self.k{
                u[i] = 1.0 - self.rho;
                u[self.k+i] = self.rho;
            }
        }
            
        let u0 = u;
        let dense = store; // if "store", then also solve ODE with dense output
        let n_steps_init = 4;

        let t0 = child_time;
        let t1 = time;

        let (times, sol) = ode.solve_dopri45(u0, t0, t1, dense, n_steps_init, tol, EquationType::ProbabilityDensity).expect("could not calculate branch probabiltiy E(t) and D(t) in likelihood function");

        let n = sol.len()-1;

        let alpha: f64 = sol[n][(self.k+1)..(2*self.k)].iter().sum();

        let mut p = Vec::new();
        //let p = sol.last();
        for i in 0..self.k{
            p.push(sol[n][i]);
        }

        for i in 0..self.k{
            p.push(sol[n][self.k+i] / alpha);
        }

        if store{
            let (e, d) = split_e_and_d(sol, self.k);

            let extinction_probability = MonotonicCubicSpline::new(times.clone(), e, self.k, true);
            let branch_probability = MonotonicCubicSpline::new(times.clone(), d, self.k, true);

            node.extinction_probability = Some(extinction_probability);
            node.subtree_probability = Some(branch_probability);
        }

        log_sf += alpha.ln();
            
        return (p, log_sf);
    }
}


fn split_e_and_d(sol: Vec<Vec<f64>>, k: usize) -> (Vec<Vec<f64>>, Vec<Vec<f64>>){
    let mut e = Vec::new();
    let mut d = Vec::new();
    for sol_i in sol{
        let mut ei = Vec::new();
        let mut di = Vec::new();
        for j in 0..k{
            ei.push(sol_i[j]);
            di.push(sol_i[j+k]);
        }
        e.push(ei);
        d.push(di);
    }

    return (e, d);
}

