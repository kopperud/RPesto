use regex::Regex;
use once_cell::sync::Lazy;

pub fn stripcomments(contents: &str) -> String {
    //let re = Regex::new(r"\[.*?\]").unwrap();
    static RE: Lazy<Regex> = Lazy::new(|| Regex::new(r"\[.*?\]").unwrap());
    let stripped_contents = RE.replace_all(contents, "");

    stripped_contents.to_string()
}

pub fn find_newick_string(contents: String) -> String {
    let lparen = contents.find('(').expect("expected to find opening parenthesis. are you sure the file has tree?");

    let semicolon = contents.rfind(';').expect("expected to find closing semicolon (;), are you sure your file has newick trees?");

    let res = contents
        .get(lparen..(semicolon+1))
        .unwrap()
        .to_string();

    res
}

// normalize
// meaning sum to 1
pub fn normalize(x: &mut Vec<f64>) -> (){
    let mut sum = 0.0;
    for xi in x.iter(){
        sum += xi;
    }

    for xi in x.iter_mut(){
        (*xi) = (*xi) / sum;
    }
}
