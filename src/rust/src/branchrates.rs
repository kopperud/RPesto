use crate::marginal_probability::Marginal;
use crate::height::*;
use crate::models::*;
use crate::tree::*;

// the marginal probability trait 
pub trait BranchRates{
    fn net_diversification( &self, node: &mut Box<Node> ) -> ();
    fn net_diversification_preorder( &self, node: &mut Box<Node>, time: f64) -> ();
}

impl BranchRates for ShiftBD{
    fn net_diversification( &self, tree: &mut Box<Node>) -> (){
        let time = treeheight(tree);

        self.net_diversification_preorder( tree, time );
    }

    fn net_diversification_preorder( &self, node: &mut Box<Node>, time: f64) -> () {
        let t0 = time;
        let t1 = t0 - node.length;

        let n = 6;

        let mut r = 0.0;

        // t0 is older time
        // t1 is younger time
        // meaning that, t0 >= t1
        // and delta_t =< 0.0
        let delta_t = (t1 - t0) / ((n-1) as f64);

        let mut time = t0;

        for _ in 0..n{
            let s = node.marginal_probability(time);
            let k = s.len();

            let mut ri = 0.0;
            for i in 0..k{
                ri += s[i] * (self.lambda[i] - self.mu[i]);
            }
            r += ri;
            time += delta_t;
        }

        let net_div = r / (n as f64);

        //assign to node
        node.r = Some(net_div);

        for child_node in node.children.iter_mut(){
            self.net_diversification_preorder(child_node, t1);
        }
    }
}

