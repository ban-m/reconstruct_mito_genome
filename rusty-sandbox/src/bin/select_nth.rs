extern crate rand;
fn main() {
    use rand::{thread_rng, Rng};
    let mut rng = thread_rng();
    for t in 0..100 {
        let len = rng.gen_range(0, 100_000);
        let vector: Vec<usize> = (0..len).map(|_| rng.gen()).collect();
        let mut sorted = vector.clone();
        sorted.sort();
        for i in 0..len {
            let nth = select_nth_by(&vector, i, |x, y| x.partial_cmp(y).unwrap());
            assert_eq!(nth, sorted[i]);
        }
        eprintln!("{},{},OK", t, len);
    }
}

fn select_nth_by<T: Clone, F: Fn(&T, &T) -> Ordering>(xs: &[T], n: usize, cmp: F) -> T {
    use std::cmp::Ordering::*;
    assert!(xs.len() >= n);
    let pivot = &xs[xs.len() / 2];
    let (small, same) = xs
        .iter()
        .fold((0, 0), |(small, same), x| match cmp(pivot, x) {
            Equal => (small, same + 1),
            Greater => (small + 1, same),
            Less => (small, same),
        });
    // Recursive call.
    if n < small {
        // We can remove elements more than `pivot` from `xs`.
        let xs: Vec<_> = xs
            .iter()
            .filter(|x| cmp(pivot, x) == Greater)
            .cloned()
            .collect();
        select_nth_by(&xs, n, cmp)
    } else if small + same <= n {
        let xs: Vec<_> = xs
            .iter()
            .filter(|x| cmp(pivot, x) == Less)
            .cloned()
            .collect();
        select_nth_by(&xs, n - small - same, cmp)
    } else {
        pivot.clone()
    }
}
