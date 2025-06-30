use crate::odesolver::Solve;
use crate::tree::*;
use crate::height::*;
use crate::branch_probability::*;
use crate::models::*;
use crate::spline::*;

use rayon::iter::ParallelIterator;
use rayon::prelude::*;

// the likelihood trait 
pub trait Likelihood<T>{
    fn likelihood( &self, tree: &mut Box<Node>, tol: f64, store: bool) -> f64;
    fn likelihood_po( &self, node: &mut Box<Node>, ode: &T, time: f64, tol: f64, store: bool) -> (Vec<f64>, f64);
}

// the likelihood implementation for constant-rate birth death model
impl Likelihood<BranchProbability> for ConstantBD{
    fn likelihood(&self, tree: &mut Box<Node>, tol: f64, store: bool) -> f64{
        let height = treeheight(&tree);

        let time = height;

        let ode = BranchProbability::new(self.lambda, self.mu);

        let (p, sf) = self.likelihood_po(tree, &ode, time, tol, store);


        let lnl = p[0].ln() + sf;
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

        let (times, sol) = ode.solve_dopri45(u0, t0, t1, dense, n_steps_init, tol);
    
        let mut p = sol[0].clone();

        if store{
            node.u_dense = Some(sol.clone());
            node.t_dense = Some(times);

        }

        log_sf += p[1].ln();
        p[1] = 1.0;
            
        return (p, log_sf);
    }
}



// the likelihood implementation for birth-death-shift model
impl Likelihood<BranchProbabilityMultiState> for ShiftBD{
    fn likelihood(&self, tree: &mut Box<Node>, tol: f64, store: bool) -> f64{
        let height = treeheight(&tree);

        let time = height;

        let ode = BranchProbabilityMultiState::new(self.lambda.clone(), self.mu.clone(), self.eta);

        let (p, sf) = self.likelihood_po(tree, &ode, time, tol, store);

        let root_prior = vec![1.0 / (self.k as f64); self.k];

        let mut lnl = 0.0;
        let mut pr = 0.0;
        for i in 0..self.k{
            pr += root_prior[i] * p[self.k+i];
        }
        lnl += pr.ln() + sf;

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

        let (times, sol) = ode.solve_dopri45(u0, t0, t1, dense, n_steps_init, tol);

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
            //node.u_dense = Some(p.clone());
            node.u_dense = Some(sol.clone());
            node.t_dense = Some(times.clone());

            let (e, d) = split_e_and_d(sol, self.k);

            let extinction_probability = MonotonicCubicSpline::new(times.clone(), e, self.k);
            let branch_probability = MonotonicCubicSpline::new(times.clone(), d, self.k);

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

