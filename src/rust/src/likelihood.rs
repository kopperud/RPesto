use crate::odesolver::Solve;
use crate::parser::*;
use crate::tree::*;
use crate::height::*;
use crate::branch_probability::*;
use crate::models::*;


// the likelihood trait 
pub trait Likelihood<T>{
    fn likelihood( &self, tree: &Box<Node>, tol: f64 ) -> f64;
    fn likelihood_po( &self, node: &Box<Node>, ode: &T, time: f64, tol: f64) -> (Vec<f64>, f64);
}

// the likelihood implementation
impl Likelihood<BranchProbability> for ConstantBD{
    fn likelihood(&self, tree: &Box<Node>, tol: f64) -> f64{
        let height = treeheight(&tree);

        let time = height;

        let ode = BranchProbability::new(self.lambda, self.mu, height, tol);

        let (p, sf) = self.likelihood_po(&tree, &ode, time, tol);

        let lnl = p[0].ln() + sf;
        return lnl;
    }

    fn likelihood_po(&self, node: &Box<Node>, ode: &BranchProbability, time: f64, tol: f64) -> (Vec<f64>, f64){

        let mut u = 1.0;
        let mut log_sf = 0.0;

        let child_time = time - node.length;

        for child in node.children.iter(){
            let (child_u, child_sf) = self.likelihood_po(child, ode, child_time, tol);

            u *= child_u[0];
            log_sf += child_sf;
        }

        let n_children = node.children.len();
        if n_children > 1{
            for _ in 0..(n_children-1){
                u *= self.lambda;
            }
        }else{
            u = self.rho;
        }
            
        let u0 = vec![u];
        let dense = false;
        let n_steps_init = 4;

        let t0 = child_time;
        let t1 = time;

        //println!("t0: {}, t1: {}", t0, t1);

        let (_, sol) = ode.solve_dopri45(u0, t0, t1, dense, n_steps_init, tol);

        //println!("u0: {:?}", u);
        //println!("sol: {:?}", sol[0][0]);

        log_sf += sol[0][0].ln();

        let p = vec![1.0];
            
        return (p, log_sf);
    }
}



// the likelihood implementation for SSE
impl Likelihood<BranchProbabilityMultiState> for ShiftBD{
    fn likelihood(&self, tree: &Box<Node>, tol: f64) -> f64{
        let height = treeheight(&tree);

        let time = height;

        let ode = BranchProbabilityMultiState::new(self.lambda.clone(), self.mu.clone(), self.eta, self.rho, height, self.k, tol);

        let (p, sf) = self.likelihood_po(&tree, &ode, time, tol);

        let mut lnl = 0.0;
        let root_prior = vec![1.0 / (self.k as f64); self.k];

        for i in 0..self.k{
            lnl += (root_prior[i] * p[i]).ln();
            //println!("lnl = {}", lnl);
        }
        lnl += sf;
        //println!("lnl = {}", lnl);

        return lnl;
    }

    fn likelihood_po(&self, node: &Box<Node>, ode: &BranchProbabilityMultiState, time: f64, tol: f64) -> (Vec<f64>, f64){

        let mut u = vec![1.0; self.k];
        let mut log_sf = 0.0;

        let child_time = time - node.length;

        for child in node.children.iter(){
            let (child_u, child_sf) = self.likelihood_po(child, ode, child_time, tol);

            for i in 0..self.k{
                u[i] *= child_u[i];
            }
            log_sf += child_sf;
        }

        let n_children = node.children.len();
        if n_children > 1{
            for _ in 0..(n_children-1){
                for i in 0..self.k{
                    u[i] *= self.lambda[i];
                }
            }
        }else{
            for i in 0..self.k{
                u[i] = self.rho;
            }
        }
            
        let u0 = u;
        let dense = false;
        let n_steps_init = 4;

        let t0 = child_time;
        let t1 = time;

        println!("t0: {}, t1: {}", t0, t1);

        let (_, sol) = ode.solve_dopri45(u0, t0, t1, dense, n_steps_init, tol);

        //println!("u0: {:?}", u);
        println!("sol: {:?}", sol[0]);

        let alpha: f64 = sol[0].iter().sum();

        let z = alpha < 0.0;
        if z {
            println!("log alpha = {}", alpha.ln());
            println!("p(before normalizing) = {:?}", sol[0]);
        }

        //let p = sol[0].clone();
        let mut p = Vec::new();
        for i in 0..self.k{
            p.push(sol[0][i] / alpha);
        }
        //log_sf += sol[0].iter();
        log_sf += alpha.ln();
            
        return (p, log_sf);
    }
}



