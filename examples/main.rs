use std::io;

fn main() {
    let stdin = io::stdin();
    let mut buffer = String::new();
    let _ = stdin.read_line(&mut buffer);
    let factorization = nfs::factorize(&buffer);

    print!("Number factored as ");
    for (i, (p, e)) in factorization.iter().enumerate() {
        print!("{}^{}", p, e);
        if i + 1 < factorization.len() {
            print!(" * ");
        }
    }
    print!("\n");
}
