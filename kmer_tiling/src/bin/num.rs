fn main(){
    let args:Vec<_> = std::env::args().collect();
    let k:usize = args[1].parse().unwrap();
    let dist:usize = args[2].parse().unwrap();
    let res = num_of_dist_k(k,dist);
    println!("{}",res);
}


fn log_fact(n:usize)->f64{
    (1..=n).map(|e| (e as f64).ln()).sum()
}

fn num_of_dist_k(k:usize,dist:usize)->u64{
    (0..=dist).map(|indel|{
        let (ins,del,sub) = (indel, indel, dist - indel);
        let mat = k - ins - sub;
        let pat = log_fact(ins + del + sub + mat) - log_fact(ins) - log_fact(del) - log_fact(sub) - log_fact(mat);
        let tot = pat + (ins + sub) as f64* 4f64.ln();
        tot.exp().floor() as u64
    }).sum()
}
