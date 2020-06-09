use last_decompose::create_windows;
fn main() {
    env_logger::Builder::from_env(env_logger::Env::default().default_filter_or("debug")).init();
    let len = 4000;
    let covs: Vec<_> = (0..len)
        .map(|idx| match idx {
            x if 100 < x && x < 500 => 10,
            _ => 100,
        })
        .collect();
    for (c, s, e) in create_windows(0, len, &covs) {
        println!("{}:{}-{}", c, s, e);
    }
}
