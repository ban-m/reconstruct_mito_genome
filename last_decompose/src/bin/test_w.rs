use last_decompose::create_windows;
fn main() {
    let covs = vec![10; 5000];
    for range in create_windows(0, 5000, &covs) {
        eprintln!("{:?}", range);
    }
}
