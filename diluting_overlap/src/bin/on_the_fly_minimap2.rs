extern crate bio;
extern crate libc;
use bio::io::fastq;
use std::ffi::CString;
use std::io::{BufWriter, Write};
fn main() {
    let fifo1 = unsafe {
        let file = CString::new("./test1").unwrap();
        libc::mkfifo(file.as_ptr(), 0o644);
        "./test1"
    };
    let fifo2 = unsafe {
        let file = CString::new("./test2").unwrap();
        libc::mkfifo(file.as_ptr(), 0o644);
        "./test2"
    };
    eprintln!("Setup done.");
    let test = "/grid/ban-m/arabidopsis_thaliana/sequel/dilution/0.fq";
    let test: Vec<_> = fastq::Reader::from_file(test)
        .unwrap()
        .records()
        .filter_map(|e| e.ok())
        .collect();
    let child = std::process::Command::new("minimap2")
        .args(&["-x", "ava-pb"])
        .args(&["-t", "12"])
        .arg(fifo1)
        .arg(fifo2)
        .stdout(std::process::Stdio::piped())
        .spawn()
        .expect("failed to exec");
    {
        eprintln!("Open...");
        let mut input1 = BufWriter::new(
            std::fs::OpenOptions::new()
                .write(true)
                .read(false)
                .create(false)
                .open(fifo1)
                .expect("open fail."),
        );
        for record in &test{
            write_retry(record,&mut input1);
        }
    }
    {
        let mut input2 = BufWriter::new(
            std::fs::OpenOptions::new()
                .write(true)
                .read(false)
                .create(false)
                .open(fifo2)
                .unwrap(),
        );
        for record in &test {
            write_retry(record, &mut input2);
        }
    }
    let out = child.wait_with_output().expect("failed to get out");
    println!("{}", String::from_utf8_lossy(&out.stdout));
    eprintln!("Removing...");
    std::process::Command::new("rm")
        .arg(fifo1)
        .arg(fifo2)
        .output()
        .expect("");
}

fn write_retry<W:std::io::Write>(record:&fastq::Record, wtr:&mut BufWriter<W>){
    let mut output = vec![b'>'];
    output.extend(record.id().as_bytes());
    output.push(b'\n');
    output.extend(record.seq());
    output.push(b'\n');
    let mut is_cond = 0;
    while let Err(_) = wtr.write_all(&output){
        is_cond += 1;
    }
    if is_cond != 0 {
        eprintln!("Write {} sucess in {} trial", record.id(), is_cond);
    }
    
}
