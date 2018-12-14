extern crate bio;
use std::path::Path;
fn main()->std::io::Result<()>{
    let args:Vec<_> = std::env::args().collect();
    let mut assembly:Vec<_> = bio::io::fasta::Reader::from_file(&Path::new(&args[1]))?
    .records()
        .filter_map(|e|e.ok())
        .collect();
    assembly.sort_by_key(|e|e.seq().len());
    let name = args[2].clone();
    let total:usize = assembly.iter().map(|e|e.seq().len()).sum();
    let nfif = {
        let (mut sofar,mut index) = (0,0);
        while sofar < total/2{
            sofar += assembly[index].seq().len();
            index += 1;
        }
        assembly[index-1].seq().len()
    };
    let num_of_contig = assembly.len();
    let mean = total/num_of_contig;
    let num_of_n = assembly.iter().flat_map(|e|e.seq().iter())
        .filter(|e| match **e as char {
            'a' | 'A' | 't' | 'T' | 'g' | 'G' | 'c' | 'C' => false,
            _ => true
        })
        .count();
    let num_of_gap:usize = assembly.iter().map(|e|num_of_gap(e.seq())+1).sum::<usize>() - 1;
    println!("{}\t{}\t{}\t{}\t{}\t{}\t{}",
             name,
             total,
             nfif,
             mean,
             num_of_contig,
             num_of_n,
             num_of_gap
    );
    Ok(())
}

fn num_of_gap(seq:&[u8])->usize{
    let (gaps,_) = seq.iter().map(|e|match *e as char {
        'a' | 'A' | 't' | 'T' | 'g' | 'G' | 'c' | 'C' => true,
        _ => false
    })
        .fold((0,false),|(gaps,previous),current|
              if previous  && !current {
                  (gaps+1,current)
              }else{
                  (gaps,current)
              });
    gaps
}

#[test]
fn gaps(){
    let test:Vec<u8> = vec!['a','a','n','n','t','n','n','c','n']
        .into_iter().map(|e|e as u8).collect();
    assert_eq!(num_of_gap(&test),3);
}
