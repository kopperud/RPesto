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
    pub children: Vec<Box<Node>>,
}


#[extendr]
impl Node {
    pub fn new(label: String, length: f64) -> Self {
        let children: Vec<Box<Node>> = Vec::new(); 

        Self {label, length, children }
    }
}

extendr_module! {
    //mod RPesto;
    mod tree; 
    impl Node;
}
