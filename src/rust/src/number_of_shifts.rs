use crate::odesolver::*;
use crate::tree::*;
use crate::height::*;
use crate::branch_probability::*;
use crate::models::*;
use crate::spline::*;
use crate::utils::*;

pub struct ShiftProblem{
    pub lambda: Vec<f64>,
    pub mu: Vec<f64>,
    pub eta: f64,
    pub k: usize,
    pub subtree_probability: MonotonicCubicSpline,
    pub forward_probability: MonotonicCubicSpline,
}

impl ShiftProblem{
    pub fn new(
        lambda: Vec<f64>,
        mu: Vec<f64>,
        eta: f64,
        subtree_probability: MonotonicCubicSpline,
        forward_probability: MonotonicCubicSpline,
    ) -> ShiftProblem{
        let k = lambda.len();
        let res = ShiftProblem{
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
impl Gradient for ShiftProblem{
    fn gradient(&self, dN: &mut Vec<f64>, _N: &Vec<f64>, t: &f64 ) -> (){
        let r = self.eta / (self.k as f64 - 1.0);

        let Dt = self.subtree_probability.interpolate(*t);
        let Ft = self.forward_probability.interpolate(*t);

        let mut St = Vec::new();
        for (Dti, Fti) in Dt.iter().zip(Ft){
            St.push(Dti * Fti);
        }
        normalize(&mut St);

        let s1 = Dt.iter().fold(0.0, |acc, x| acc + x);
        let mut s2 = 0.0;
        for j in 0..self.k{
            s2 += St[j] / Dt[j];
        }

        dN[0] = - r * (s1 * s2 - 1.0);
    }
}

// the likelihood trait 
pub trait NumberOfShifts<T>{
    fn number_of_shifts( &self, tree: &mut Box<Node>, tol: f64) -> ();
    fn number_of_shifts_pre( &self, node: &mut Box<Node>, time: f64, tol: f64) -> ();
}

impl NumberOfShifts<BranchProbabilityMultiState> for ShiftBD{
    fn number_of_shifts( &self, tree: &mut Box<Node>, tol: f64) -> (){
        let height = treeheight(&tree);
        let time = height;
        
        self.number_of_shifts_pre(tree, time, tol);
    }

    fn number_of_shifts_pre(&self, node: &mut Box<Node>, time: f64, tol: f64) -> (){

        let shift_problem = ShiftProblem::new(self.lambda.clone(), self.mu.clone(), self.eta, node.subtree_probability.clone().unwrap(), node.forward_probability.clone().unwrap());

        let n_steps_init = 5;

        let u0 = vec![0.0];
        let t0 = time;
        let t1 = time - node.length;

        let (_, no_shifts_v) = shift_problem.solve_dopri45(u0, t0, t1, false, n_steps_init, tol);

        let no_shifts = no_shifts_v[0][0];

        //node.forward_probability = Some(forward_probability);
        node.number_of_shifts = Some(no_shifts);

        for child_node in node.children.iter_mut(){
            self.number_of_shifts_pre(child_node, time - node.length, tol);
        }
    }
}
   
