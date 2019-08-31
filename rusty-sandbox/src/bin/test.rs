// Todo convert sam to paf.
extern crate bio_utils;
extern crate rusty_sandbox;
fn main() -> std::io::Result<()> {
    fn get_input() -> std::io::Result<String> {
        let mut input = String::new();
        std::io::stdin().read_line(&mut input)?;
        Ok(input.trim().to_string())
    }
    let input = get_input()?;
    print!("{}",input);
    // let args: Vec<_> = std::env::args().collect();
    // let sams = rusty_sandbox::open_sam_into_hashmap(&args[1])?;
    // eprintln!("hai. {}", sams.len());
    // for sam in &sams["16766"]{
    //     let (start, end) = sam.get_range();
    //     let qname = sam.q_name();
    //     let rname = sam.r_name();
    //     println!("{}\t{}\t{}\t{}", qname, rname, start, end);
    // }
    Ok(())
}
