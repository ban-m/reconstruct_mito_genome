// Peak call from start and stop TSV.
// input:tsv file containing [contig ID]/t[position]\t[# of start read]\t[# of stop read]
use std::fs::File;
use std::io::Read;
use std::path::Path;
fn main() -> std::io::Result<()> {
    let args: Vec<_> = std::env::args().collect();
    let input = open_tsv_file(&args[1])?;
    let mut current_contig_id = String::new();
    let mut chunk: Vec<(usize, usize, usize)> = vec![];
    println!("refname\tposition\tnumber_of_start_read\tnumber_of_stop_read");
    for line in input.lines().skip(1) {
        let contents: Vec<&str> = line.split('\t').collect();
        if contents[0] != current_contig_id {
            if !chunk.is_empty(){
                for (position, num_of_start, num_of_stop) in peak_call(&chunk) {
                    println!(
                        "{}\t{}\t{}\t{}",
                        current_contig_id, position, num_of_start, num_of_stop
                    );
                }
            }
            current_contig_id = contents[0].to_string();
            chunk.clear();
        } else {
            let position:usize = contents[1].parse().unwrap();
            let start_num:usize = contents[3].parse().unwrap();
            let stop_num:usize = contents[4].parse().unwrap();
            chunk.push((position, start_num, stop_num));
        }
    }
    for (position, num_of_start, num_of_stop) in peak_call(&chunk) {
        println!(
            "{}\t{}\t{}\t{}",
            current_contig_id, position, num_of_start, num_of_stop
        );
    }
    Ok(())
}

fn peak_call(input: &Vec<(usize, usize, usize)>) -> Vec<(usize, usize, usize)> {
    // peak calling. First,filtering out positions with 'coverage' less than 30,
    // then, aggregate remainig peaks.
    let mut peaks = vec![];
    let mut current_position = 0;
    let (mut max_position, mut max_start, mut max_stop) = (0, 0, 0);
    for (position, start, stop) in input
        .into_iter()
        .filter(|(_, start_read, stop_read)| start_read > &30 || stop_read > &30)
    {
        if current_position + 100 < *position {
            peaks.push((max_position, max_start, max_stop));
            max_position = *position;
            max_start = *start;
            max_stop = *stop;
            current_position = *position;
        } else {
            if start > &max_start && stop > &max_stop {
                max_position = *position;
                max_start = *start;
                max_stop = *stop;
            }
            current_position = *position;
        }
    }
    peaks.push((max_position, max_start, max_stop));
    peaks
}

fn open_tsv_file(file: &str) -> std::io::Result<String> {
    let mut file = File::open(&Path::new(file))?;
    let mut input = String::new();
    file.read_to_string(&mut input)?;
    Ok(input)
}
