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

        //s.push('[');
        let x = format!("[&netdiv={},nshifts={}]", self.r.unwrap(), self.number_of_shifts.unwrap());
        s.push_str(x.as_str());
        //s.push(']');
        s.push(':');
        let brlen = self.length.to_string();
        s.push_str(brlen.as_str());
    }
}



