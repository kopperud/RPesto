use crate::marginal_probability::Marginal;
use crate::height::*;
use crate::models::*;
use crate::tree::*;

use itertools::izip;

// the marginal probability trait 
pub trait BranchRates{
    fn net_diversification( &self, node: &mut Box<Node> ) -> ();
    fn net_diversification_preorder( &self, node: &mut Box<Node>, time: f64) -> ();

    fn speciation( &self, node: &mut Box<Node> ) -> ();
    fn speciation_preorder( &self, node: &mut Box<Node>, time: f64) -> ();

    fn extinction( &self, node: &mut Box<Node> ) -> ();
    fn extinction_preorder( &self, node: &mut Box<Node>, time: f64) -> ();
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

        // delta_netdiv: the change in netdiv across the branch
        let netdiv_old: f64 = izip!(
            node.marginal_probability(t0),
            &self.lambda,
            &self.mu,
        )
            .map(|(p,l,m)| p*(l-m))
            .sum();

        let netdiv_young: f64 = izip!(
            node.marginal_probability(t1),
            &self.lambda,
            &self.mu,
        )
            .map(|(p,l,m)| p*(l-m))
            .sum();

        let delta_netdiv = netdiv_young - netdiv_old;

        node.delta_netdiv = Some(delta_netdiv);

        for child_node in node.children.iter_mut(){
            self.net_diversification_preorder(child_node, t1);
        }
    }

    fn speciation( &self, tree: &mut Box<Node>) -> (){
        let time = treeheight(tree);
            self.speciation_preorder( tree, time );
        }

    fn speciation_preorder( &self, node: &mut Box<Node>, time: f64) -> () {
        let t0 = time;
        let t1 = t0 - node.length;

        let n = 6;
        let mut lambda_acc = 0.0;

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
                ri += s[i] * (self.lambda[i]);
            }
            lambda_acc += ri;
            time += delta_t;
        }

        let lambda_mean = lambda_acc / (n as f64);

        //assign to node
        node.lambda = Some(lambda_mean);

        for child_node in node.children.iter_mut(){
            self.speciation_preorder(child_node, t1);
        }
    }

    fn extinction( &self, tree: &mut Box<Node>) -> (){
        let time = treeheight(tree);
            self.extinction_preorder( tree, time );
        }

    fn extinction_preorder( &self, node: &mut Box<Node>, time: f64) -> () {
        let t0 = time;
        let t1 = t0 - node.length;

        let n = 6;
        let mut mu_acc = 0.0;

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
                ri += s[i] * (self.mu[i]);
            }
            mu_acc += ri;
            time += delta_t;
        }

        let mu_mean = mu_acc / (n as f64);

        //assign to node
        node.mu = Some(mu_mean);

        for child_node in node.children.iter_mut(){
            self.extinction_preorder(child_node, t1);
        }
    }
}


