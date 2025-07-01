use crate::tree::*;
use crate::utils::*;

// the marginal probability trait 
pub trait Marginal{
    fn marginal_probability( &self, time: f64) -> Vec<f64>;
}

impl Marginal for Node{
    fn marginal_probability( &self, time: f64) -> Vec<f64> {
        let mut s = Vec::new();

        // the probability of the subtree given
        // that the category is j at time t
        let d = self.subtree_probability
            .as_ref()
            .expect("expected to be able to unwrap d")
            .interpolate(time);

        // the probability of the category
        // being j, given the part of the
        // tree that is not descended from t
        let f = self.forward_probability
            .as_ref()
            .expect("expected to be able to unwrap f")
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

