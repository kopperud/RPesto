/*
#[derive(Debug, Default)]
pub struct Branch {
    pub index: i32,
    pub time: f64,
    pub outbounds: RefCell<Rc<Node>>,
}
*/

use extendr_api::prelude::*;

#[derive(Debug, Default)]
#[extendr]
/// @export
pub struct Node {
    pub label: String,
    pub length: f64,
    //pub index: i32,
    pub children: Vec<Box<Node>>,
    pub u_old: Option<Vec<f64>>,
}


#[extendr]
impl Node {
    pub fn new(label: String, length: f64) -> Self {
        let children: Vec<Box<Node>> = Vec::new(); 
        //let index = 0;
        //Self {label, length, index, children }
        Self {label, length, children, u_old: None}
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
    //mod RPesto;
    mod tree; 
    impl Node;
}
