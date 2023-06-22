use rug::{self, Integer};
use std::io;

pub mod linalg;

fn main() {
    let stdin = io::stdin();
    let mut buffer = String::new();
    let _ = stdin.read_line(&mut buffer);
    let n = Integer::parse(buffer);
}
