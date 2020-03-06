use super::construct_initial_weights;
//use super::INIT_PICK_PROB;
use super::entropy;
use super::get_max_coverage;
use super::{serialize, to_pos};
use super::{utils, ERead, Read};
use super::{FACTOR, INIT_BETA};
use poa_hmm::*;
use rand::{thread_rng, Rng, SeedableRng};
use rand_xoshiro::Xoshiro256StarStar;
use rayon::prelude::*;
const BETA_STEP: f64 = 1.1;
const DEFAULT_WEIGHT: f64 = 0.;
// const BETA_OFFSET: f64 = 0.1;
const ENTROPY_THR: f64 = 0.2;
const PICK_PROB: f64 = 0.02;

pub struct AlnParam<F>
where
    F: Fn(u8, u8) -> i32,
{
    ins: i32,
    del: i32,
    score: F,
}

fn score(x: u8, y: u8) -> i32 {
    if x == y {
        3
    } else {
        -4
    }
}
pub const DEFAULT_ALN: AlnParam<fn(u8, u8) -> i32> = AlnParam {
    ins: -6,
    del: -6,
    score: score,
};

struct ModelFactory<'a> {
    chunks: Vec<Vec<&'a [u8]>>,
    weights: Vec<Vec<f64>>,
}

impl<'a> ModelFactory<'a> {
    fn new(chain_len: usize, data: &'a [Read]) -> Self {
        let mut chunks: Vec<Vec<&[u8]>> = vec![vec![]; chain_len];
        let weights: Vec<Vec<f64>> = vec![vec![]; chain_len];
        for read in data.iter() {
            for x in read.iter() {
                chunks[x.0].push(&x.1);
            }
        }
        Self { chunks, weights }
    }
    fn generate_model<F>(
        &mut self,
        ws: &[Vec<f64>],
        reads: &[Read],
        cl: usize,
        &AlnParam {
            ins,
            del,
            ref score,
        }: &AlnParam<F>,
    ) -> Vec<POA>
    where
        F: Fn(u8, u8) -> i32 + std::marker::Sync,
    {
        assert!(self.weights.iter().all(|ws| ws.is_empty()));
        for (read, w) in reads.iter().zip(ws) {
            for &(pos, _) in read.iter() {
                self.weights[pos].push(w[cl]);
            }
        }
        assert_eq!(self.weights.len(), self.chunks.len());
        let res: Vec<_> = self
            .chunks
            .par_iter()
            .zip(self.weights.par_iter())
            .map(|(chunks, ws)| POA::generate_w_param(chunks, ws, ins, del, score))
            .collect();
        self.weights.iter_mut().for_each(|ws| ws.clear());
        res
    }
    fn update_model<F>(
        &mut self,
        updates: &[bool],
        ws: &[Vec<f64>],
        reads: &[Read],
        cl: usize,
        &AlnParam {
            ins,
            del,
            ref score,
        }: &AlnParam<F>,
        _class: usize,
    ) -> Vec<POA>
    where
        F: Fn(u8, u8) -> i32 + std::marker::Sync,
    {
        let default_weight = DEFAULT_WEIGHT; // / class as f64;
        assert!(self.weights.iter().all(|ws| ws.is_empty()));
        for ((read, w), &b) in reads.iter().zip(ws).zip(updates) {
            let w = if b { 0. } else { w[cl] + default_weight };
            for &(pos, _) in read.iter() {
                self.weights[pos].push(w);
            }
        }
        assert_eq!(self.weights.len(), self.chunks.len());
        let res: Vec<_> = self
            .chunks
            .par_iter()
            .zip(self.weights.par_iter())
            .map(|(chunks, ws)| POA::generate_w_param(chunks, ws, ins, del, score))
            .collect();
        self.weights.iter_mut().for_each(|ws| ws.clear());
        res
    }
}

fn convert(c: &dbg_hmm::Config) -> poa_hmm::Config {
    poa_hmm::Config {
        mismatch: c.mismatch,
        base_freq: c.base_freq,
        p_match: c.p_match,
        p_ins: c.p_ins,
        p_del: c.p_del,
        p_extend_ins: c.p_extend_ins,
        p_extend_del: c.p_extend_del,
        p_del_to_ins: c.p_del_to_ins,
    }
}

pub fn soft_clustering_poa<F>(
    data: &[ERead],
    label: &[u8],
    forbidden: &[Vec<u8>],
    cluster_num: usize,
    answer: &[u8],
    config: &dbg_hmm::Config,
    aln: &AlnParam<F>,
) -> Vec<Vec<f64>>
where
    F: Fn(u8, u8) -> i32 + std::marker::Sync,
{
    assert!(cluster_num > 1);
    let config = &convert(config);
    let (matrix_pos, chain_len) = to_pos(data);
    let data = serialize(data, &matrix_pos);
    let id: u64 = thread_rng().gen::<u64>() % 100_000;
    let seed = data.len() as u64;
    let mut weights_of_reads =
        crate::construct_initial_weights(label, forbidden, cluster_num, data.len(), seed);
    //let pick_up_len = 2 * PICK_PROB.recip().floor() as usize;
    // let max_coverage = get_max_coverage(&data, chain_len);
    // let beta_step = 1. + (BETA_STEP - 1.) / (max_coverage as f64).log10();
    debug!("Chain length is {}", chain_len);
    let mut mf = ModelFactory::new(chain_len, &data);
    let mut models: Vec<Vec<POA>> = (0..cluster_num)
        .map(|cl| mf.generate_model(&weights_of_reads, &data, cl, aln))
        .collect();
    let (border, datasize) = (label.len(), data.len() as f64);
    let mut rng: Xoshiro256StarStar = SeedableRng::seed_from_u64(data.len() as u64 * 23);
    let pick_num = (datasize * PICK_PROB).floor() as usize;
    let mut picks: Vec<_> = (0..data.len()).skip(border).collect();
    picks.shuffle(&mut rng);
    use rand::seq::SliceRandom;
    let max_entropy = (cluster_num as f64).ln() * datasize;
    let mut ws: Vec<f64> = (0..cluster_num)
        .map(|i| weights_of_reads.iter().map(|g| g[i]).sum::<f64>() / datasize)
        .collect();
    // super::updates_flags(&mut updates, &weights_of_reads, &mut rng, PICK_PROB, border);
    // let mut beta =
    //     search_initial_beta(&data, label, forbidden, cluster_num, config, chain_len, aln);
    let weights = super::variant_calling::variant_call_poa(&models, &data, config, true);
    let mut variants: Vec<_> = weights.iter().map(|&b| b * b).collect();
    let mut loop_num = 1;
    let mut previous_lk = calc_lk(&models, &ws, &data, config);
    let mut saved = weights_of_reads.clone();
    debug!("LK\t0\t{}\t{}", id, previous_lk);
    loop {
        while !picks.is_empty() {
            // select reads to be updated.
            let mut updates = vec![false; data.len()];
            for _ in 0..pick_num {
                if let Some(pos) = picks.pop() {
                    updates[pos] = true;
                }
            }
            let max = variants.iter().fold(0., |x, &y| if x < y { y } else { x });
            let betas: Vec<_> = variants.iter().map(|v| v / max).collect();
            let w = &mut weights_of_reads;
            update_weights(
                w, &mut ws, border, &data, &models, &updates, &betas, config, None,
            );
            let wor = &weights_of_reads;
            //super::updates_flags(&mut updates, wor, &mut rng, PICK_PROB, border);
            models = (0..cluster_num)
                .map(|cl| mf.update_model(&updates, wor, &data, cl, aln, cluster_num))
                .collect();
        }
        picks = (0..data.len()).skip(border).collect();
        picks.shuffle(&mut rng);
        let wor = &weights_of_reads;
        report(id, wor, border, answer, &ws, 1.);
        let lk = calc_lk(&models, &ws, &data, config);
        debug!("LK\t{}\t{}\t{}", loop_num, id, lk);
        if lk <= previous_lk {
            break;
        }
        previous_lk = lk;
        saved = weights_of_reads.clone();
        loop_num += 1;
        // if wor.iter().map(|e| entropy(e)).sum::<f64>() / max_entropy < ENTROPY_THR {
        //     break 'outer;
        // }
        let weights = super::variant_calling::variant_call_poa(&models, &data, config, true);
        variants.iter_mut().zip(weights).for_each(|(v, w)| {
            *v = (*v + w * w) / 2.;
        });
    }
    let weights_of_reads = saved;
    report(id, &weights_of_reads, border, answer, &ws, 1.);
    weights_of_reads
}

fn calc_lk(models: &[Vec<POA>], ws: &[f64], data: &[Read], config: &Config) -> f64 {
    fn lk(r: &Read, m: &[POA], c: &Config) -> f64 {
        r.iter().map(|&(p, ref u)| m[p].forward(u, c)).sum::<f64>()
    }
    fn lk_of_read(read: &Read, ws: &[f64], models: &[Vec<POA>], config: &Config) -> f64 {
        let lks: Vec<_> = models
            .iter()
            .zip(ws.iter())
            .map(|(m, w)| w.ln() + lk(read, m, config))
            .collect();
        super::utils::logsumexp(&lks)
    }
    data.par_iter()
        .map(|read| lk_of_read(read, ws, models, config))
        .sum::<f64>()
}

// Update the K-dimensional vectors at each position towarding the axis(1,1,1,1..,1) from the center.
fn update_centrize_vectors<R: rand::Rng>(
    centrize_vectors: &mut [Vec<f64>],
    models: &[Vec<POA>],
    data: &[Read],
    config: &Config,
    chain_len: usize,
    rng: &mut R,
) {
    //Clear vectors.
    centrize_vectors
        .iter_mut()
        .for_each(|vec| vec.iter_mut().for_each(|x| *x = 0.));
    // Let's calcurate the average lk on each position.
    // To this end, we need aux. vector to store how many units
    // we use at each position. Let it be counts.
    let mut counts = vec![0.; chain_len];
    let pick_num = 25.min(data.len());
    use rand::seq::SliceRandom;
    for read in data.choose_multiple(rng, pick_num) {
        let likelihoods: Vec<_> = models
            .par_iter()
            .enumerate()
            .map(|(cluster, ms)| {
                let lks = read
                    .par_iter()
                    .map(|&(pos, ref unit)| (pos, ms[pos].forward(unit, config)))
                    .collect::<Vec<_>>();
                (cluster, lks)
            })
            .collect();
        for (cluster, lks) in likelihoods {
            for (pos, lk) in lks {
                counts[pos] += 1.;
                centrize_vectors[cluster][pos] += lk;
            }
        }
    }
    for vector in centrize_vectors.iter_mut() {
        vector
            .iter_mut()
            .zip(counts.iter())
            .for_each(|(sum, count)| *sum /= count);
    }
    // Here we've compute not c, but c/K.
    // Thus we can calculate <c,1>/K 1by just adding the likelihood of each culaster
    // This is because we've already add 1 for K times for
    // each position of the read. Thus the /K is not needed.
    // Compute <c,1>/k.
    let projection: Vec<_> = centrize_vectors
        .iter()
        .fold(vec![0.; chain_len], |mut acc, xs| {
            acc.iter_mut().zip(xs.iter()).for_each(|(a, &x)| *a += x);
            acc
        });
    // Then, compuate <c,1>/K 1 - c = <c/K, 1> 1 - c/K * K.
    let cluster_num = models.len() as f64;
    for vector in centrize_vectors.iter_mut() {
        vector
            .iter_mut()
            .zip(projection.iter())
            .for_each(|(cv, &p)| {
                *cv = p - *cv * cluster_num;
            });
    }
}

fn search_initial_beta<F>(
    data: &[Read],
    label: &[u8],
    forbidden: &[Vec<u8>],
    cluster: usize,
    c: &Config,
    chain_len: usize,
    aln: &AlnParam<F>,
) -> f64
where
    F: Fn(u8, u8) -> i32 + std::marker::Sync,
{
    let seed = data.iter().map(|e| e.len()).sum::<usize>() as u64;
    let mut rng: Xoshiro256StarStar = SeedableRng::seed_from_u64(seed);
    let weight_of_read = construct_initial_weights(label, forbidden, cluster, data.len(), seed);
    let border = label.len();
    let mut mf = ModelFactory::new(chain_len, data);
    let mut beta = INIT_BETA;
    let soe = weight_of_read.iter().map(|e| entropy(e)).sum::<f64>();
    let wor = &weight_of_read;
    loop {
        let betas = vec![beta; chain_len];
        let next = soe_after_sampling(
            &betas, data, wor, border, &mut rng, cluster, &mut mf, c, aln,
        );
        debug!("SEARCH\t{:.3}\t{:.3}\t{:.3}", beta, next, soe - next);
        if soe > next {
            break;
        } else {
            beta *= FACTOR;
        }
    }
    beta /= FACTOR;
    loop {
        let betas = vec![beta; chain_len];
        let next = soe_after_sampling(
            &betas, data, wor, border, &mut rng, cluster, &mut mf, c, aln,
        );
        debug!("SEARCH\t{:.3}\t{:.3}\t{:.3}", beta, next, soe - next);
        if soe < next {
            break;
        } else {
            beta /= FACTOR;
        }
    }
    beta
}

fn soe_after_sampling<R: Rng, F>(
    betas: &[f64],
    data: &[Read],
    wor: &[Vec<f64>],
    border: usize,
    rng: &mut R,
    cluster_num: usize,
    mf: &mut ModelFactory,
    c: &Config,
    aln: &AlnParam<F>,
) -> f64
where
    F: Fn(u8, u8) -> i32 + std::marker::Sync,
{
    let datasize = data.len() as f64;
    let mut updates = vec![false; data.len()];
    let mut ws: Vec<f64> = (0..cluster_num)
        .map(|i| wor.iter().map(|g| g[i]).sum::<f64>() / datasize)
        .collect();
    let mut models: Vec<Vec<POA>> = (0..cluster_num)
        .map(|cl| mf.generate_model(&wor, data, cl, aln))
        .collect();
    let mut wor = wor.to_vec();
    //for _ in 0..SAMPLING {
    for _ in 0..10 {
        use super::INIT_PICK_PROB;
        super::updates_flags(&mut updates, &wor, rng, INIT_PICK_PROB, border);
        update_weights(
            &mut wor, &mut ws, border, data, &models, &updates, &betas, c, None,
        );
        models = (0..cluster_num)
            .map(|cl| mf.generate_model(&wor, &data, cl, aln))
            .collect();
    }
    wor.iter().map(|e| entropy(e)).sum::<f64>()
}

fn report(
    id: u64,
    weight_of_read: &[Vec<f64>],
    border: usize,
    answer: &[u8],
    ws: &[f64],
    beta: f64,
) -> f64 {
    let correct = weight_of_read
        .iter()
        .skip(border)
        .zip(answer.iter())
        .filter(|&(weights, &ans)| {
            weights
                .iter()
                .filter(|&&g| g > weights[ans as usize])
                .count()
                == 0
        })
        .count();
    let count: Vec<_> = (0..ws.len())
        .map(|cl| {
            weight_of_read
                .iter()
                .filter(|r| r.iter().all(|&w| w <= r[cl]))
                .count()
        })
        .map(|e| format!("{}", e))
        .collect();
    let count = count.join("\t");
    let acc = correct as f64 / answer.len() as f64;
    let ws: Vec<_> = (0..ws.len())
        .map(|cl| weight_of_read.iter().map(|r| r[cl]).sum::<f64>())
        .collect();
    let pi: Vec<_> = ws.iter().map(|e| format!("{:.2}", *e)).collect();
    let pi = pi.join("\t");
    let soe = weight_of_read.iter().map(|e| entropy(e)).sum::<f64>();
    info!(
        "Summary\t{}\t{:.3}\t{}\t{}\t{:.3}\t{}\t{:.2}",
        id, soe, pi, count, beta, correct, acc
    );
    0.
}

fn update_weights(
    weight_of_read: &mut [Vec<f64>],
    ws: &mut [f64],
    border: usize,
    data: &[Read],
    models: &[Vec<POA>],
    updates: &[bool],
    betas: &[f64],
    c: &Config,
    centrize_vectors: Option<&[Vec<f64>]>,
) {
    data.par_iter()
        .zip(weight_of_read.par_iter_mut())
        .zip(updates.par_iter())
        .skip(border)
        .filter(|&(_, &b)| b)
        .for_each(|((read, weights), _)| {
            compute_log_probs(models, &ws, read, weights, c, &betas, &centrize_vectors);
            let tot = utils::logsumexp(weights);
            weights.iter_mut().for_each(|w| *w = (*w - tot).exp());
        });
    let datasize = data.len() as f64;
    ws.iter_mut().enumerate().for_each(|(cl, w)| {
        *w = weight_of_read.iter().map(|e| e[cl]).sum::<f64>() / datasize;
    });
    assert!((1. - ws.iter().sum::<f64>()).abs() < 0.001);
}

fn compute_log_probs(
    models: &[Vec<POA>],
    ws: &[f64],
    read: &Read,
    weights: &mut Vec<f64>,
    c: &Config,
    betas: &[f64],
    centrize_vectors: &Option<&[Vec<f64>]>,
) {
    assert_eq!(models.len(), ws.len());
    assert_eq!(models.len(), weights.len());
    assert!((1. - ws.iter().sum::<f64>()).abs() < 0.001);
    if let Some(centrize_vectors) = centrize_vectors {
        models
            .par_iter()
            .zip(weights.par_iter_mut())
            .zip(ws.par_iter())
            .zip(centrize_vectors.par_iter())
            .for_each(|(((ms, g), w), cv)| {
                *g = read
                    .par_iter()
                    .map(|&(pos, ref u)| {
                        let model = &ms[pos];
                        let beta = betas[pos];
                        let centrize = cv[pos];
                        beta * (model.forward(u, c) + centrize)
                    })
                    .sum::<f64>()
                    + w.ln();
            });
    } else {
        models
            .par_iter()
            .zip(weights.par_iter_mut())
            .zip(ws.par_iter())
            .for_each(|((ms, g), w)| {
                *g = read
                    .par_iter()
                    .map(|&(pos, ref u)| {
                        let model = &ms[pos];
                        let beta = betas[pos];
                        beta * model.forward(u, c)
                    })
                    .sum::<f64>()
                    + w.ln();
            });
    }
}
