use std::collections::VecDeque;

pub fn tokenize(s: &str) -> VecDeque<&str> {
    //let mut tokens: Vec<&str> = Vec::new();
    let mut tokens: VecDeque<&str> = VecDeque::new();
    const SPECIAL_TOKENS: &[char] = &['(', ')', ',', ';'];
    let chars: Vec<char> = s.chars().collect();

    let mut start = 0;
    let mut end = 0;

    while end < chars.len(){
        let c = chars[end];
        if SPECIAL_TOKENS.contains(&c){
            if start != end {
                tokens.push_back(&s[start..end]);
            }
            if c != 0xA as char{
                tokens.push_back(&s[end..end+1]);
            }
            end += 1;
            start = end;
        }else{
            end += 1;
        }
    }

    if start != end{
        tokens.push_back(&s[start..end]);
    }

    // remove the root branch length, if applicable
    let n_tokens = tokens.len();
    let second_last_token = *tokens.get(n_tokens-2)
        .expect("should have been able to get second-to-last token");
    if second_last_token.contains(':'){
        tokens.remove(n_tokens-2);
    }

    tokens
}
