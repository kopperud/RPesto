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
