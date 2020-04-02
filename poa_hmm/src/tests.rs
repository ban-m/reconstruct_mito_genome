use super::gen_sample::*;
use super::rayon::prelude::*;
use super::*;
use rand::seq::SliceRandom;
#[test]
fn it_works() {
    assert_eq!(2 + 2, 4);
}

#[test]
fn init() {
    let _poa = POA::default();
    let seed = b"ACGTAGCTGATCGTAC";
    let _poa = POA::new(seed, 1.);
    let poa = POA::new(seed, 0.4);
    assert!(poa.nodes.iter().any(|n| n.is_head));
}
#[test]
fn add() {
    let seed = b"ACGTAGCTGATCGTAC";
    POA::new(seed, 1.).add_default(seed, 1.);
}
#[test]
fn create() {
    let seed = b"ACGTAGCTGATCGTAC";
    let res = POA::new(seed, 1.).add_default(seed, 1.);
    assert_eq!(res.nodes.len(), seed.len());
}
#[test]
fn insertion_test() {
    let seed = b"ACGTAGCTGATCGTAC";
    let seq2 = b"ACGTAGCTGATTTCGTAC";
    let res = POA::new(seed, 1.).add_default(seq2, 1.);
    assert_eq!(res.nodes.len(), seed.len() + 2);
    let seed = b"ACGTAGCTGATCGTAC";
    let seq2 = b"ACGTAGCTGATTTCGTAC";
    let res = POA::new(seed, 1.).add_default(seq2, 1.);
    assert_eq!(res.nodes.len(), seed.len() + 2);
    let seed = b"ACGTAGCTGATCGTAC";
    let seq2 = b"ACGTAGCTGATCGTAC";
    let res = POA::new(seed, 1.).add_default(seq2, 1.);
    assert_eq!(res.nodes.len(), seed.len());
    let seed = b"ACGTAGCTGATCGTAC";
    let seq2 = b"ACGTAGCTGATCGTACG";
    let res = POA::new(seed, 1.).add_default(seq2, 1.);
    assert_eq!(res.nodes.len(), seed.len() + 1);
    let seq2 = b"ACGTAGCTGATCGTAC";
    let seed = b"AACGTAGCTGATCGTAC";
    let res = POA::new(seed, 1.).add_default(seq2, 1.);
    assert_eq!(res.nodes.len(), seed.len());
    let seed = b"ACGTAGCTGATCGTAC";
    let seq2 = b"AAACGTAGCTGATCGTAC";
    let res = POA::new(seed, 1.).add_default(seq2, 1.);
    assert_eq!(res.nodes.len(), seed.len() + 2);
}

#[test]
fn mismatch_test() {
    let seed = b"ACGTAGCTGATCGTAC";
    let seq2 = b"ACGTAGCTGATCGGAC";
    let res = POA::new(seed, 1.).add_default(seq2, 1.);
    assert_eq!(res.nodes.len(), seed.len() + 1);
    let seed = b"ACGTAGCTGATCGTAC";
    let seq2 = b"ACGTAGCTGATCGTCC";
    let res = POA::new(seed, 1.).add_default(seq2, 1.);
    assert_eq!(res.nodes.len(), seed.len() + 1);
    let seed = b"ACGTAGCTGAGTAC";
    let seq2 = b"ACGAGCTGATCTAC";
    let res = POA::new(seed, 1.).add_default(seq2, 1.);
    assert_eq!(res.nodes.len(), seed.len() + 2);
    let seed = b"ACGTAGCTGATCGTAC";
    let seq2 = b"CCGTAGCTGATCGTAC";
    let res = POA::new(seed, 1.).add_default(seq2, 1.);
    assert_eq!(res.nodes.len(), seed.len() + 1);
    let seed = b"ACGTAGCTGATCGTAC";
    let seq2 = b"ACGTAGCTGATCGTAG";
    let res = POA::new(seed, 1.).add_default(seq2, 1.);
    assert_eq!(res.nodes.len(), seed.len() + 1);
}
#[test]
fn deletion_test() {
    let seed = b"ACGTAGCTGATCGTAC";
    let seq2 = b"ACGTAGCTGTCGTAC";
    let res = POA::new(seed, 1.).add_default(seq2, 1.);
    assert_eq!(res.nodes.len(), seed.len());
    assert_eq!(res.num_edges(), seed.len());
    let seed = b"ACGTAGCTGATCGTAC";
    let seq2 = b"ACGAGCTGATCGTAC";
    let res = POA::new(seed, 1.).add_default(seq2, 1.);
    assert_eq!(res.nodes.len(), seed.len());
    assert_eq!(res.num_edges(), seed.len());
    let seed = b"ACGTAGCTGATCGTAC";
    let seq2 = b"ACGAGCTGATCTAC";
    let res = POA::new(seed, 1.).add_default(seq2, 1.);
    assert_eq!(res.nodes.len(), seed.len());
    assert_eq!(res.num_edges(), seed.len() + 1);
    let seed = b"ACGTAGCTGATCGTAC";
    let seq2 = b"ACGTAGCTGATCGTA";
    let res = POA::new(seed, 1.).add_default(seq2, 1.);
    assert_eq!(res.nodes.len(), seed.len());
    assert_eq!(res.num_edges(), seed.len() - 1);
    let seed = b"ACGTAGCTGATCGTAC";
    let seq2 = b"CGTAGCTGATCGTAC";
    let res = POA::new(seed, 1.).add_default(seq2, 1.);
    assert_eq!(res.nodes.len(), seed.len());
    assert_eq!(res.num_edges(), seed.len() - 1);
}

#[test]
fn merge_test() {
    let seq1 = b"ACAAT";
    let seq2 = b"ATGCT";
    let seq3 = b"ACACT";
    let res = POA::new(seq1, 1.)
        .add_default(seq2, 1.)
        .add_default(seq3, 1.);
    assert_eq!(res.nodes.len(), 8);
    let seq1 = b"GATC";
    let seq2 = b"GCTC";
    let seq3 = b"GAGTC";
    let seq4 = b"GCGTC";
    let res = POA::new(seq1, 1.)
        .add_default(seq2, 1.)
        .add_default(seq3, 1.)
        .add_default(seq4, 1.);
    assert_eq!(res.nodes.len(), 6);
    let seq1 = b"GATC";
    let seq2 = b"GCTC";
    let seq3 = b"GAGTC";
    let seq4 = b"GCGTC";
    let res = POA::new(seq1, 1.)
        .add_default(seq2, 1.)
        .add_default(seq3, 1.)
        .add_default(seq4, 1.);
    assert_eq!(res.nodes.len(), 6);
}

#[test]
fn initialize() {
    let test = [
        b"CAGTGCTAGTCGATGTCA".to_vec(),
        b"CAGTGCTAGTCGATGTCA".to_vec(),
        b"CAGTGCTAGTCGATGTCA".to_vec(),
        b"CAGTGCTAGTCGATGTCA".to_vec(),
        b"CAGTGCTAGTCGATGTCA".to_vec(),
        b"CAGTGCTAGTCGATGTCA".to_vec(),
        b"CAGTGCTAGTCGATGTCA".to_vec(),
        b"CAGTGCTAGTCGATGTCA".to_vec(),
        b"CA".to_vec(),
        b"TTTTTTGTGTGACTGTACGTGACG".to_vec(),
        b"CACACACACGTGTACGTGTTGGGGGGCTAAA".to_vec(),
    ];
    test.iter()
        .fold(POA::default(), |x, y| x.add_default(y, 1.));
}
#[test]
fn forward() {
    let test = [
        b"CAGTGCTAGTCGATGTCA".to_vec(),
        b"CA".to_vec(),
        b"TTTTTTGTGTGACTGTACGTGACG".to_vec(),
        b"CACACACACGTGTACGTGTTGGGGGGCTAAA".to_vec(),
        b"CACACACACGTGTACGTGTTGGGGGGCTAAA".to_vec(),
        b"CACACACACGTGTACGTGTTGGGGGGCTAAA".to_vec(),
    ];
    let m = POA::generate_vec(&test);
    eprintln!("{}", m);
    let lk = m.forward(b"CACACAGCAGTCAGTGCA", &DEFAULT_CONFIG);
    eprintln!("{}", lk);
    assert!(lk < 0.);
}


#[test]
fn connectivity_check() {
    let bases = b"ACTG";
    let mut rng: Xoshiro256StarStar = SeedableRng::seed_from_u64(1212132);
    let template: Vec<_> = (0..30)
        .filter_map(|_| bases.choose(&mut rng))
        .copied()
        .collect();
    let model1: Vec<Vec<_>> = (0..20)
        .map(|_| introduce_randomness(&template, &mut rng, &PROFILE))
        .collect();
    let mut m = POA::default();
    for seq in model1 {
        m = m.add_default(&seq, 1.);
        assert!(is_connected(&m));
    }
}
fn is_connected(m: &POA) -> bool {
    let mut edges = vec![vec![]; m.nodes.len()];
    for (from, n) in m.nodes.iter().enumerate() {
        for &to in n.edges.iter() {
            edges[from].push(to);
            edges[to].push(from);
        }
    }
    let mut dfs = vec![false; edges.len()];
    let mut stack = vec![0];
    'dfs: while !stack.is_empty() {
        let n = *stack.last().unwrap();
        if !dfs[n] {
            dfs[n] = true;
        }
        for &to in edges[n].iter() {
            if !dfs[to] {
                stack.push(to);
                continue 'dfs;
            }
        }
        stack.pop();
    }
    dfs.iter().all(|&e| e)
}

fn alignment<F>(xs: &[u8], ys: &[u8], del: f64, ins: f64, score: F) -> f64
where
    F: Fn(u8, u8) -> f64,
{
    let mut dp = vec![vec![0.; ys.len() + 1]; xs.len() + 1];
    for i in 0..xs.len() {
        dp[i + 1][0] = ins * (i + 1) as f64;
    }
    // for j in 0..ys.len() {
    //     dp[0][j + 1] = del * (j + 1) as f64;
    // }
    for i in 0..xs.len() {
        for j in 0..ys.len() {
            dp[i + 1][j + 1] = (dp[i][j] + score(xs[i], ys[j]))
                .max(dp[i + 1][j] + del)
                .max(dp[i][j + 1] + ins);
        }
    }
    // for line in &dp {
    //     let line: Vec<_> = line.iter().map(|x| format!("{:4.0}", x)).collect();
    // eprintln!("{}", line.join(" "));
    //}
    *dp[xs.len()]
        .iter()
        .max_by(|a, b| a.partial_cmp(&b).unwrap())
        .unwrap()
}

#[test]
fn alignment_check() {
    let bases = b"ACTG";
    let mut rng: Xoshiro256StarStar = SeedableRng::seed_from_u64(1212132);
    use rand::Rng;
    for _ in 0..100 {
        let x: Vec<_> = (0..(rng.gen_range(100, 150)))
            .filter_map(|_| bases.choose(&mut rng))
            .copied()
            .collect();
        let y: Vec<_> = (0..(rng.gen_range(100, 150)))
            .filter_map(|_| bases.choose(&mut rng))
            .copied()
            .collect();
        let score = |x, y| if x == y { -1. } else { -2. };
        let opt = alignment(&y, &x, -2., -2., score.clone());
        eprintln!("OPT:{}", opt);
        let score = |x, y| if x == y { -1 } else { -2 };
        let (poa_score, _) = POA::new(&x, 1.).align(&y, (-2, -2, score));
        let poa_score = poa_score as f64;
        eprintln!("POA:{}", poa_score);
        assert!((opt - poa_score).abs() < 0.001, "{},{}", opt, poa_score);
    }
}

#[test]
fn forward_check() {
    let bases = b"ACTG";
    let mut rng: Xoshiro256StarStar = SeedableRng::seed_from_u64(1212132);
    let template: Vec<_> = (0..30)
        .filter_map(|_| bases.choose(&mut rng))
        .copied()
        .collect();
    let model1: Vec<Vec<_>> = (0..20)
        .map(|_| introduce_randomness(&template, &mut rng, &PROFILE))
        .collect();
    let m = POA::generate_vec(&model1);
    let lk = m.forward(&template, &DEFAULT_CONFIG);
    eprintln!("{:?}", m);
    assert!(lk < 0., "{}", lk);
}

#[test]
fn random_check() {
    let bases = b"ACTG";
    let mut rng: Xoshiro256StarStar = SeedableRng::seed_from_u64(1212132);
    let template: Vec<_> = (0..150)
        .filter_map(|_| bases.choose(&mut rng))
        .copied()
        .collect();
    let model1: Vec<Vec<_>> = (0..50)
        .map(|_| introduce_randomness(&template, &mut rng, &PROFILE))
        .collect();
    eprintln!("OK");
    let model2: Vec<Vec<_>> = (0..50)
        .map(|_| {
            (0..150)
                .filter_map(|_| bases.choose(&mut rng))
                .copied()
                .collect()
        })
        .collect();
    let model1 = POA::generate_vec(&model1);
    let model2 = POA::generate_vec(&model2);
    eprintln!("{}/{}", model1, model2);
    let likelihood1 = model1.forward(&template, &DEFAULT_CONFIG);
    let likelihood2 = model2.forward(&template, &DEFAULT_CONFIG);
    assert!(likelihood1 > likelihood2, "{},{}", likelihood1, likelihood2);
}

#[test]
fn random_check_forward() {
    let bases = b"ACTG";
    for i in 0..10 {
        let mut rng: Xoshiro256StarStar = SeedableRng::seed_from_u64(i as u64);
        let template: Vec<_> = (0..150)
            .filter_map(|_| bases.choose(&mut rng))
            .copied()
            .collect();
        let model1: Vec<Vec<_>> = (0..50)
            .map(|_| introduce_randomness(&template, &mut rng, &PROFILE))
            .collect();
        let model1 = POA::generate_vec(&model1);
        eprintln!("{:?}", model1);
        for _ in 0..100 {
            let query = introduce_randomness(&template, &mut rng, &PROFILE);
            model1.forward(&query, &DEFAULT_CONFIG);
        }
    }
}

#[test]
fn hard_test() {
    let bases = b"ACTG";
    let mut rng: Xoshiro256StarStar = SeedableRng::seed_from_u64(1212132);
    let template1: Vec<_> = (0..150)
        .filter_map(|_| bases.choose(&mut rng))
        .copied()
        .collect();
    let p = Profile {
        sub: 0.005,
        ins: 0.005,
        del: 0.005,
    };
    let template2 = introduce_randomness(&template1, &mut rng, &p);
    let model1: Vec<Vec<_>> = (0..50)
        .map(|_| introduce_randomness(&template1, &mut rng, &PROFILE))
        .collect();
    let model2: Vec<Vec<_>> = (0..50)
        .map(|_| introduce_randomness(&template2, &mut rng, &PROFILE))
        .collect();
    let model1 = POA::generate_vec(&model1);
    let model2 = POA::generate_vec(&model2);
    eprintln!("{}", String::from_utf8_lossy(&template1));
    eprintln!("{}", String::from_utf8_lossy(&template2));
    let template1 = introduce_randomness(&template1, &mut rng, &PROFILE);
    let likelihood1 = model1.forward(&template1, &DEFAULT_CONFIG);
    let likelihood2 = model2.forward(&template1, &DEFAULT_CONFIG);
    assert!(
        likelihood1 > likelihood2,
        "{}\n{},{},{},{}",
        String::from_utf8_lossy(&template1),
        likelihood1,
        likelihood2,
        model1,
        model2
    );
    let template2 = introduce_randomness(&template2, &mut rng, &PROFILE);
    let likelihood1 = model1.forward(&template2, &DEFAULT_CONFIG);
    let likelihood2 = model2.forward(&template2, &DEFAULT_CONFIG);
    assert!(
        likelihood1 - 3. < likelihood2,
        "{},{},{},{}",
        likelihood1,
        likelihood2,
        model1,
        model2
    );
}
#[test]
fn mix_test_prior() {
    let bases = b"ACTG";
    let mut rng: Xoshiro256StarStar = SeedableRng::seed_from_u64(1212132);
    let len = 150;
    let template1: Vec<_> = (0..len)
        .filter_map(|_| bases.choose(&mut rng))
        .copied()
        .collect();
    let p = Profile {
        sub: 0.03,
        ins: 0.03,
        del: 0.03,
    };
    let cov = 30;
    let score = |x, y| if x == y { 3 } else { -4 };
    let parameters = (-6, -6, &score);
    for _ in 0..5 {
        let template2 = introduce_randomness(&template1, &mut rng, &p);
        let model1: Vec<Vec<_>> = (0..cov)
            .map(|_| introduce_randomness(&template1, &mut rng, &PROFILE))
            .collect();
        let model2: Vec<Vec<_>> = (0..cov)
            .map(|_| introduce_randomness(&template2, &mut rng, &PROFILE))
            .collect();
        let dataset: Vec<_> = model1
            .iter()
            .map(|e| e.as_slice())
            .chain(model2.iter().map(|e| e.as_slice()))
            .collect();
        let weight1 = vec![vec![0.8; cov], vec![0.2; cov]].concat();
        let weight2 = vec![vec![0.2; cov], vec![0.8; cov]].concat();
        let model1 = POA::generate(&dataset, &weight1, parameters);
        let model2 = POA::generate(&dataset, &weight2, parameters);
        eprintln!("{}\n{}", model1, model2);
        let num = 50;
        let correct = (0..num)
            .filter(|_| {
                let q = introduce_randomness(&template1, &mut rng, &PROFILE);
                let lk1 = model1.forward(&q, &DEFAULT_CONFIG);
                let lk2 = model2.forward(&q, &DEFAULT_CONFIG);
                eprintln!("1\t{:.3}\t{:.3}", lk1, lk2);
                lk1 > lk2
            })
            .count();
        assert!(correct >= num * 4 / 5, "{}", correct);
        let correct = (0..num)
            .filter(|_| {
                let q = introduce_randomness(&template2, &mut rng, &PROFILE);
                let lk1 = model1.forward(&q, &DEFAULT_CONFIG);
                let lk2 = model2.forward(&q, &DEFAULT_CONFIG);
                eprintln!("2\t{:.3}\t{:.3}", lk1, lk2);
                lk1 < lk2
            })
            .count();
        assert!(correct >= num * 4 / 5, "{}", correct);
    }
}


#[ignore]
#[test]
fn abundance_test_prior() {
    let bases = b"ACTG";
    let mut rng: Xoshiro256StarStar = SeedableRng::seed_from_u64(1_219);
    let p = Profile {
        sub: 0.03,
        ins: 0.03,
        del: 0.03,
    };
    let len = 150;
    let cov = 20;
    let ratio = 5;
    let errors = PROFILE;
    let score = |x, y| if x == y { 3 } else { -4 };
    let parameters = (-6, -6, &score);
    for _ in 0..3 {
        let template1: Vec<_> = (0..len)
            .filter_map(|_| bases.choose(&mut rng))
            .copied()
            .collect();
        let template2 = introduce_randomness(&template1, &mut rng, &p);
        let data1: Vec<Vec<_>> = (0..(ratio + 1) * cov)
            .map(|_| introduce_randomness(&template1, &mut rng, &PROFILE))
            .collect();
        let mut data2: Vec<Vec<_>> = (0..ratio * cov)
            .map(|_| introduce_randomness(&template1, &mut rng, &PROFILE))
            .collect();
        eprintln!("{}/{}", template1.len(), template2.len());
        data2.extend((0..cov).map(|_| introduce_randomness(&template2, &mut rng, &PROFILE)));
        let data1: Vec<_> = data1.iter().map(|e| e.as_slice()).collect();
        let data2: Vec<_> = data2.iter().map(|e| e.as_slice()).collect();
        let total = (ratio + 1) * cov;
        let weight = vec![1.; total];
        let model1 = POA::generate(&data1, &weight, parameters);
        let model2 = POA::generate(&data2, &weight, parameters);
        eprintln!("{}", model1);
        eprintln!("{}", String::from_utf8_lossy(&template1));
        eprintln!("{}", String::from_utf8_lossy(&model1.consensus()));
        eprintln!("----------------------");
        eprintln!("{}", model2);
        eprintln!("{}", String::from_utf8_lossy(&template2));
        eprintln!("{}", String::from_utf8_lossy(&model2.consensus()));
        let num = 50;
        let correct = (0..num)
            .filter(|_| {
                let q = introduce_randomness(&template1, &mut rng, &errors);
                let lk1 = model1.forward(&q, &DEFAULT_CONFIG);
                let lk2 = model2.forward(&q, &DEFAULT_CONFIG);
                //eprintln!("1\t{:.3}\t{:.3}", lk1, lk2);
                lk1 > lk2
            })
            .count();
        //eprintln!("1:{}", correct);
        assert!(correct >= num * 6 / 10, "1:{}", correct);
        let correct = (0..num)
            .filter(|_| {
                let q = introduce_randomness(&template2, &mut rng, &errors);
                let lk1 = model1.forward(&q, &DEFAULT_CONFIG);
                let lk2 = model2.forward(&q, &DEFAULT_CONFIG);
                //eprintln!("2\t{:.3}\t{:.3}", lk1, lk2);
                lk1 < lk2
            })
            .count();
        // eprintln!("2:{}", correct);
        assert!(correct >= num * 6 / 10, "2:{}", correct);
    }
}


use rand::Rng;
#[ignore]
#[test]
fn single_error_test() {
    let bases = b"ACTG";
    let coverage = 150;
    let start = 20;
    let step = 3;
    let len = 150;
    let results: Vec<_> = (start..coverage)
        .step_by(step)
        .map(|cov| {
            let mut rng: Xoshiro256StarStar = SeedableRng::seed_from_u64(1_234_567);
            let seed = rng.gen_range(0, 100_000);
            let mut rng: Xoshiro256StarStar = SeedableRng::seed_from_u64(seed);
            let template1: Vec<_> = (0..len)
                .filter_map(|_| bases.choose(&mut rng))
                .copied()
                .collect();
            let template2 = introduce_errors(&template1, &mut rng, 1, 0, 0);
            let sub = check(&template1, &template2, &mut rng, cov);
            let template2 = introduce_errors(&template1, &mut rng, 0, 1, 0);
            let del = check(&template1, &template2, &mut rng, cov);
            let template2 = introduce_errors(&template1, &mut rng, 0, 0, 1);
            let ins = check(&template1, &template2, &mut rng, cov);
            (cov, (sub, del, ins))
        })
        .collect();
    let (sub, del, ins) = results
        .iter()
        .fold((0, 0, 0), |(x, y, z), &(_, (a, b, c))| {
            (x + a, y + b, z + c)
        });
    for (cov, res) in results {
        eprintln!("Cov:{},Sub:{},Del:{},Ins:{}", cov, res.0, res.1, res.2);
    }
    eprintln!("Tot:{}", (start..coverage).step_by(step).count() * 100);
    eprintln!("Sub:{},Del:{},Ins:{}", sub, del, ins);
    assert!(false);
}

fn check<R: rand::Rng>(t1: &[u8], t2: &[u8], rng: &mut R, cov: usize) -> usize {
    let model1: Vec<_> = (0..cov)
        .map(|_| introduce_randomness(&t1, rng, &PROFILE))
        .collect();
    let model2: Vec<_> = (0..cov)
        .map(|_| introduce_randomness(&t2, rng, &PROFILE))
        .collect();
    let seqs: Vec<_> = model1
        .iter()
        .chain(model2.iter())
        .map(|e| e.as_slice())
        .collect();
    let score = |x, y| if x == y { 3 } else { -4 };
    let parameters = (-6, -6, &score);
    let weight1 = vec![vec![1.; cov], vec![0.; cov]].concat();
    let weight2 = vec![vec![0.; cov], vec![1.; cov]].concat();
    let m1 = POA::generate(&seqs, &weight1, parameters);
    let m2 = POA::generate(&seqs, &weight2, parameters);
    eprintln!("{}\t{}\t{}", cov, m1, m2);
    let tests: Vec<_> = (0..100)
        .map(|e| {
            if e % 2 == 0 {
                (e, introduce_randomness(&t1, rng, &PROFILE))
            } else {
                (e, introduce_randomness(&t2, rng, &PROFILE))
            }
        })
        .collect();
    let correct = tests
        .par_iter()
        .filter(|(e, q)| {
            if e % 2 == 0 {
                m1.forward(&q, &DEFAULT_CONFIG) > m2.forward(&q, &DEFAULT_CONFIG)
            } else {
                m1.forward(&q, &DEFAULT_CONFIG) < m2.forward(&q, &DEFAULT_CONFIG)
            }
        })
        .count();
    correct
}

#[test]
fn low_coverage_test() {
    let mut rng: Xoshiro256StarStar = SeedableRng::seed_from_u64(121212);
    let template1: Vec<_> = b"CAGTGTCAGTGCTAGCT".to_vec();
    let template2: Vec<_> = b"CAGTGTCTGTGCTAGCT".to_vec();
    let model1: Vec<Vec<_>> = (0..10)
        .map(|_| introduce_randomness(&template1, &mut rng, &PROFILE))
        .collect();
    let model2: Vec<Vec<_>> = (0..10)
        .map(|_| introduce_randomness(&template2, &mut rng, &PROFILE))
        .collect();
    for m in &model1 {
        eprintln!("1:{}", String::from_utf8_lossy(m));
    }
    for m in &model2 {
        eprintln!("2:{}", String::from_utf8_lossy(m));
    }
    let test1 = introduce_randomness(&template1, &mut rng, &PROFILE);
    let test2 = introduce_randomness(&template2, &mut rng, &PROFILE);
    eprintln!("1:{}", String::from_utf8_lossy(&test1));
    eprintln!("2:{}", String::from_utf8_lossy(&test2));
    let model1 = POA::generate_vec(&model1);
    eprintln!("Model1:\n{:?}", model1);
    let model2 = POA::generate_vec(&model2);
    eprintln!("Model2:\n{:?}", model2);
    {
        let likelihood1 = model1.forward(&test1, &DEFAULT_CONFIG);
        let likelihood2 = model2.forward(&test1, &DEFAULT_CONFIG);
        assert!(likelihood1 > likelihood2, "{},{}", likelihood1, likelihood2);

    }
    {
        let likelihood1 = model1.forward(&test2, &DEFAULT_CONFIG);
        let likelihood2 = model2.forward(&test2, &DEFAULT_CONFIG);
        assert!(likelihood1 < likelihood2, "{},{}", likelihood1, likelihood2);
    }
}

#[test]
fn high_coverage_test() {
    let mut rng: Xoshiro256StarStar = SeedableRng::seed_from_u64(1212123324);
    let template1 = generate_seq(&mut rng, 100);
    let template2 = introduce_errors(&template1, &mut rng, 1, 1, 1);
    // let template1: Vec<_> = b"CAGTGTCAGTGCTAGCT".to_vec();
    // let template2: Vec<_> = b"CAGTGTCTGTGCTAGCT".to_vec();
    let model1: Vec<Vec<_>> = (0..200)
        .map(|_| introduce_randomness(&template1, &mut rng, &PROFILE))
        .collect();
    let model2: Vec<Vec<_>> = (0..200)
        .map(|_| introduce_randomness(&template2, &mut rng, &PROFILE))
        .collect();
    let test1 = introduce_randomness(&template1, &mut rng, &PROFILE);
    let test2 = introduce_randomness(&template2, &mut rng, &PROFILE);
    eprintln!("1:{}", String::from_utf8_lossy(&template1));
    eprintln!("2:{}", String::from_utf8_lossy(&template2));
    let dataset: Vec<_> = model1
        .iter()
        .chain(model2.iter())
        .map(|e| e.as_slice())
        .collect();
    let weight1 = vec![vec![1.; 200], vec![0.; 200]].concat();
    let weight2 = vec![vec![0.; 200], vec![1.; 200]].concat();
    let model1 = POA::generate_by(&dataset, &weight1, &DEFAULT_CONFIG);
    eprintln!("Model1:{:?}", model1);
    let model2 = POA::generate_by(&dataset, &weight2, &DEFAULT_CONFIG);
    eprintln!("Model2:{:?}", model2);
    {
        let likelihood1 = model1.forward(&test1, &DEFAULT_CONFIG);
        let likelihood2 = model2.forward(&test1, &DEFAULT_CONFIG);
        assert!(likelihood1 > likelihood2, "{},{}", likelihood1, likelihood2);
    }
    {
        let likelihood1 = model1.forward(&test2, &DEFAULT_CONFIG);
        let likelihood2 = model2.forward(&test2, &DEFAULT_CONFIG);
        assert!(likelihood1 < likelihood2, "{},{}", likelihood1, likelihood2);
        eprintln!("2:{:.4}\t{:.4}", likelihood1, likelihood2);
    }
}
