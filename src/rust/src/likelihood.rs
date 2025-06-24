use crate::odesolver::Solve;
use crate::parser::*;
use crate::tree::*;
use crate::height::*;
use crate::branch_probability::*;

pub fn likelihood_cbdp(lambda: f64, mu: f64, s: String, tol: f64) -> f64{
    let tree = parse_tree(s);
    let height = treeheight(&tree);

    let time = height;

    let ode = BranchProbability::new(lambda, mu, height, tol);

    let (p, sf) = likelihood_cbdp_po(lambda, mu, &ode, &tree, time, tol);

    let lnl = p.ln() + sf;
    return lnl;
}

fn likelihood_cbdp_po(
            lambda: f64,
            mu: f64,
            ode: &BranchProbability,
            node: &Box<Node>,
            time: f64,
            tol: f64
) -> (f64, f64){

    let mut u = 1.0;
    let mut log_sf = 0.0;

    let child_time = time - node.length;

    for child in node.children.iter(){
        let (child_u, child_sf) = likelihood_cbdp_po(lambda, mu, ode, child, child_time, tol);

        u *= child_u;
        //u *= lambda;
        log_sf += child_sf;
    }

    let n_children = node.children.len();
    if n_children > 1{
        for _ in 0..(n_children-1){
            u *= lambda;
        }
    }

    //u /= lambda;
    

    let u0 = vec![u];
    let dense = false;
    let n_steps_init = 4;

    let t0 = child_time;
    let t1 = time;

    //println!("t0: {}, t1: {}", t0, t1);

    let (_, sol) = ode.solve_dopri45(u0, t0, t1, dense, n_steps_init, tol);

    //println!("u0: {:?}", u);
    //println!("sol: {:?}", sol[0][0]);

    let mut p = sol[0][0];

    log_sf += p.ln();

    p = 1.0;
        
    return (p, log_sf);
}
