extern crate bio;
extern crate serde;
#[macro_use]
extern crate serde_derive;
extern crate serde_json;
use bio::io::fasta;
#[derive(Serialize,Deserialize)]
struct Contig{
    name:String,
    length:usize,
}

fn main()->std::io::Result<()>{
    let args:Vec<_> = std::env::args().collect();
    let result:Vec<_> = fasta::Reader::from_file(&std::path::Path::new(&args[1]))?
    .records()
        .filter_map(|e|e.ok())
        .map(|rec|
             Contig{name:rec.id().to_string(),
                    length:rec.seq().len()})
        .collect();
    println!("{}",serde_json::ser::to_string_pretty(&result).unwrap());
    Ok(())
}
