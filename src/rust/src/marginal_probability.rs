use crate::odesolver::Solve;
use crate::tree::*;
use crate::branch_probability::*;
use crate::models::*;
use crate::utils::*;

// the marginal probability trait 
pub trait Marginal<T>{
    fn marginal_probability( &self, node: & Box<Node>, time: f64) -> Vec<f64>;
}

impl Marginal<BranchProbabilityMultiState> for ShiftBD{
    fn marginal_probability( &self, node: & Box<Node>, time: f64) -> Vec<f64> {
        let mut s = Vec::new();

        // the probability of the subtree given
        // that the category is j at time t
        let d = node.subtree_probability
            .as_ref()
            .unwrap()
            .interpolate(time);

        // the probability of the category
        // being j, given the part of the
        // tree that is not descended from t
        let f = node.forward_probability
            .as_ref()
            .unwrap()
            .interpolate(time);

        // un-normalized marginal probability of rate category j
        for (di, fi) in d.iter().zip(f){
            s.push(di * fi);
        }

        // normalize such that s(t) sums to one
        // over the rate categories
        normalize(&mut s);

        return s;
    }
}

