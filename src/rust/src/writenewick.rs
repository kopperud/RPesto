use crate::tree::*;

pub trait WriteNewick{
    fn writenewick( &self ) -> String;
    fn writenewick_preorder( &self, s: &mut String) -> ();
}

impl WriteNewick for Node{
    fn writenewick( &self ) -> String{
        let mut s = "".to_string();

        self.writenewick_preorder(&mut s);

        s.push(';');

        return s;
    }

    fn writenewick_preorder( &self, s: &mut String) -> () {
        let n_children = self.children.len();

        let is_tip = self.children.is_empty();        

        if is_tip{
            s.push_str(self.label.as_str());
        }else{
            s.push('(');
            for (i, child_node) in self.children.iter().enumerate(){
                child_node.writenewick_preorder(s);

                if i < (n_children-1){
                    s.push(',');
                }
            }
            s.push(')');
        }

        let mut items = Vec::new();

        add_variable(self.r, &mut items, "mean_netdiv");
        add_variable(self.delta_netdiv, &mut items, "delta_netdiv");
        add_variable(self.lambda, &mut items, "mean_lambda");
        add_variable(self.mu, &mut items, "mean_mu");
        add_variable(self.number_of_shifts, &mut items, "nshift");
        add_variable(self.bayes_factor, &mut items, "shift_bf");

        if !items.is_empty(){
            let joined = items.join(",");
            let x = format!("[&{}]", joined);
            s.push_str(x.as_str());
        }
        
        s.push(':');
        let brlen = self.length.to_string();
        s.push_str(brlen.as_str());
    }
}


fn add_variable(variable: Option<f64>, items: &mut Vec<String>, label: &str) -> (){
    match variable{
       Some(x) => items.push(format!("{}={}", label, x)),
       _ => (),
    }
}


