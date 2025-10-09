/*
#[derive(Debug, Default)]
pub struct Branch {
    pub index: i32,
    pub time: f64,
    pub outbounds: RefCell<Rc<Node>>,
}
*/

use extendr_api::prelude::*;

use crate::spline::MonotonicCubicSpline;

#[derive(Debug, Default)]
#[extendr]
/// @export
pub struct Node {
    pub label: String,
    pub length: f64,
    pub children: Vec<Box<Node>>,

    // these are functions that return probabilities
    pub extinction_probability: Option<MonotonicCubicSpline>,
    pub subtree_probability: Option<MonotonicCubicSpline>,
    pub forward_probability: Option<MonotonicCubicSpline>,

    // these are summary statistics that we want in the newick output
    pub r: Option<f64>,
    pub delta_netdiv: Option<f64>,
    pub lambda: Option<f64>,
    pub delta_lambda: Option<f64>,
    pub mu: Option<f64>,
    pub delta_mu: Option<f64>,
    pub number_of_shifts: Option<f64>,
    pub bayes_factor: Option<f64>,
    pub tip_netdiv: Option<f64>,
}


#[extendr]
impl Node {
    pub fn new(label: String, length: f64) -> Self {
        let children: Vec<Box<Node>> = Vec::new(); 
        Self {
            label,
            length,
            children,
            extinction_probability: None,
            subtree_probability: None,
            forward_probability: None,
            r: None,
            delta_netdiv: None,
            lambda: None,
            delta_lambda: None,
            mu: None,
            delta_mu: None,
            number_of_shifts: None,
            bayes_factor: None,
            tip_netdiv: None,
        }
    }

    /*
    fn reindex(&self) -> (){
        let node_index: i32 = 1;
        let tip_index: i32 = 1;
        
    }

    fn reindex_postorder(&mut self, tip_index: &mut i32) -> (){
        let n_children = self.children.len();
        if n_children == 0{
            self.index = tip_index.clone();
            (*tip_index) += 1;
        }else{
           for child in self.children.iter(){
               child.reindex_postorder(tip_index);
           }
        }
    }

    fn reindex_preorder(&self, node_index: &mut i32) -> (){
    }
    */

}

extendr_module! {
    mod tree; 
    impl Node;
}
