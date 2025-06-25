use crate::odesolver::Solve;
use crate::parser::*;
use crate::tree::*;
use crate::height::*;
use crate::branch_probability::*;
use crate::models::*;


// the likelihood trait 
pub trait Likelihood{
    fn likelihood( &self, tree: &Box<Node>, tol: f64 ) -> f64;
    fn likelihood_po( &self, node: &Box<Node>, ode: &BranchProbability, time: f64, tol: f64) -> (Vec<f64>, f64);
}

// the likelihood implementation
impl Likelihood for ConstantBD{
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




