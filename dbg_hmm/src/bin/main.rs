extern crate dbg_hmm;
extern crate edlib_sys;
fn main() {
    let target = b"ATATATAT";
    let query = b"ATATAGTAT";
    let ops = edlib_sys::global(target, query);
    println!("{:?}", ops);
}
