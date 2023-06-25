use rug::{Complete, Integer};
use std::io;

mod lanczos;
mod linalg;
mod mod_sqrt;
mod qs;

fn main() {
    let stdin = io::stdin();
    let mut buffer = String::new();
    let _ = stdin.read_line(&mut buffer);
    let n = Integer::parse(buffer).unwrap().complete();
    let (p, q) = qs::factorize(&n);
    println!(
        "Factored n with\np = {}\nq = {}",
        p.to_string(),
        q.to_string()
    );
}
