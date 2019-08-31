extern crate bio;
extern crate libc;
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
    let child = std::process::Command::new("cat")
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
                .unwrap(),
        );
        eprintln!("Yes");
        writeln!(&mut input1, "this is a test").unwrap();
        writeln!(&mut input1, "this is a test 2nd").unwrap();
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
        writeln!(&mut input2, "this is a test").unwrap();
    }
    let out = child.wait_with_output().expect("failed to get out");
    println!("{:?}", String::from_utf8_lossy(&out.stdout));
    eprintln!("Removing...");
    std::process::Command::new("rm")
        .arg(fifo1)
        .arg(fifo2)
        .output()
        .expect("");
}
//    let test = "/grid/ban-m/arabidopsis_thaliana/sequel/dilution/0.fq";
    // let test = fastq::Reader::from_file(test)
    //     .unwrap()
    //     .records()
    //     .filter_map(|e| e.ok());

    // let child = std::process::Command::new("minimap2")
    //     .args(&["-x", "ava-pb"])
    //     .arg(fifo1)
    //     .arg(fifo2)
    //     .stdout(std::process::Stdio::piped())
    //     .spawn()
    //     .expect("failed to exec");
        // for (idx,rec) in test.enumerate() {
        //     if idx % 100 == 0{
        //         eprintln!("Yes.{}",idx);
        //     }
