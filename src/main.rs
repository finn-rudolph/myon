mod gfpolynomial;
mod lanczos;
mod linalg;
mod nfs;
mod nt;
mod params;
mod polynomial;
mod sqrt;

use std::env;
use std::io;
use std::io::Write;

use rug::Complete;
use rug::Integer;

fn main() {
    env::set_var("RUST_LOG", "debug");
    env_logger::init();

    print!("Enter number to be factored: ");
    let _ = io::stdout().flush();
    let mut buf = String::new();
    io::stdin()
        .read_line(&mut buf)
        .expect("Failed to read integer.");
    let n = Integer::parse(buf).unwrap().complete();
    let factors = nfs::factorize(&n);

    println!("Found the following factorizations:\n");
    for a in factors {
        println!("{} * {}", &a, (&n / &a).complete());
    }
}
