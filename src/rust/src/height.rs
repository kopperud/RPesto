use crate::tree::*;

pub fn treeheight(tree: &Box<Node>) -> f64{
    let nd = node_depths(tree);

    let height = nd
        .iter()
        .cloned()
        .fold(-1./0. /* -inf */, f64::max);

    return height;
}

fn node_depths(tree: &Box<Node>) -> Vec<f64> {
    let mut nd: Vec<f64> = Vec::new();

    let time = 0.0;
    depths(&mut nd, tree, time);

    return nd;
}

fn depths(nd: &mut Vec<f64>, node: &Box<Node>, time: f64) -> (){
    nd.push(time);

    for child in node.children.iter(){
        depths(nd, child, time + child.length);
    }
}

