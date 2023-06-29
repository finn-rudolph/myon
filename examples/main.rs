use std::env;
use std::io;

use log::error;

fn main() {
    env::set_var("RUST_LOG", "info");
    env_logger::init();

    let stdin = io::stdin();
    let mut buffer = String::new();
    let _ = stdin.read_line(&mut buffer);
    let factorization = match myon::factorize(&buffer) {
        Ok(v) => v,
        Err(e) => {
            error!("{}", e);
            return;
        }
    };

    print!("Number factored as ");
    for (i, (p, e)) in factorization.iter().enumerate() {
        print!("{}^{}", p, e);
        if i + 1 < factorization.len() {
            print!(" * ");
        }
    }
    print!("\n");
}
