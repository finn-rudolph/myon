mod lanczos;
mod linalg;
mod nfs;
mod nt;
mod params;

use std::env;
use std::io;

use rug::{ops::Pow, Integer};

fn read_int<T: std::str::FromStr>() -> T
where
    <T as std::str::FromStr>::Err: std::fmt::Debug,
{
    let mut buf = String::new();
    io::stdin()
        .read_line(&mut buf)
        .expect("Failed to read integer.");
    buf.trim().parse().expect("Input is not a valid integer.")
}

fn main() {
    env::set_var("RUST_LOG", "info");
    env_logger::init();

    println!("You can factor numbers of the form r^e - s, where r and |s| are preferably small.");
    print!("r: ");
    let r = read_int::<u32>();
    print!("e: ");
    let e = read_int::<u32>();
    print!("s: ");
    let s = read_int::<i32>();

    let n: Integer = Integer::from(r).pow(e) - s;
    let a = nfs::factorize(r, e, s);

    println!("Number factored as {} * {}", &a, n / &a);
}
