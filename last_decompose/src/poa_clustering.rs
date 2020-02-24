use super::construct_initial_weights;
//use super::INIT_PICK_PROB;
use super::get_max_coverage;
use super::{entropy, ENTROPY_THR};
use super::{serialize, to_pos};
use super::{utils, ERead, Read, LOOP_NUM};
use super::{FACTOR, INIT_BETA};
use poa_hmm::*;
use rand::{thread_rng, Rng, SeedableRng};
use rand_xoshiro::Xoshiro256StarStar;
use rayon::prelude::*;
const BETA_STEP: f64 = 1.09;
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
    fn generate_model(
        &mut self,
        ws: &[Vec<f64>],
        reads: &[Read],
        cl: usize,
        config: &Config,
    ) -> Vec<POA> {
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
            .map(|(chunks, ws)| poa_hmm::PartialOrderAlignment::generate(chunks, ws, &config))
            .collect();
        self.weights.iter_mut().for_each(|ws| ws.clear());
        res
    }
    fn update_model(
        &mut self,
        ws: &[Vec<f64>],
        mask: &[bool],
        reads: &[Read],
        cl: usize,
        models: &mut Vec<POA>,
        config: &Config,
    ) {
        for ((read, w), _) in reads.iter().zip(ws).zip(mask) {
            //let w = if b { 0. } else { w[cl] };
            let w = w[cl];
            for &(pos, _) in read.iter() {
                self.weights[pos].push(w);
            }
        }
        assert_eq!(self.weights.len(), self.chunks.len());
        assert_eq!(self.weights.len(), models.len());
        self.chunks
            .par_iter()
            .zip(self.weights.par_iter())
            .zip(models.par_iter_mut())
            .for_each(|((chunks, ws), m)| {
                *m = poa_hmm::PartialOrderAlignment::generate(chunks, ws, &config);
            });
        self.weights.iter_mut().for_each(|ws| ws.clear());
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

pub fn soft_clustering_poa(
    data: &[ERead],
    label: &[u8],
    forbidden: &[Vec<u8>],
    cluster_num: usize,
    answer: &[u8],
    config: &dbg_hmm::Config,
) -> Vec<Vec<f64>> {
    assert!(cluster_num > 1);
    let config = &convert(config);
    let (matrix_pos, chain_len) = to_pos(data);
    let data = serialize(data, &matrix_pos);
    let id: u64 = thread_rng().gen::<u64>() % 100_000;
    let seed = data.len() as u64;
    let mut weights_of_reads =
        crate::construct_initial_weights(label, forbidden, cluster_num, data.len(), seed);
    // let pick_prob = INIT_PICK_PROB;
    // let pick_up_len = INIT_PICK_PROB.recip().floor() as usize * LOOP_NUM;
    let pick_up_len = LOOP_NUM;
    let max_coverage = get_max_coverage(&data, chain_len);
    let beta_step = 1. + (BETA_STEP - 1.) / (max_coverage as f64).log10();
    let mut beta = search_initial_beta(&data, label, forbidden, cluster_num, config, chain_len);
    debug!("Chain length is {}", chain_len);
    let mut mf = ModelFactory::new(chain_len, &data);
    let mut models: Vec<Vec<POA>> = (0..cluster_num)
        .map(|cl| mf.generate_model(&weights_of_reads, &data, cl, config))
        .collect();
    let updates = vec![false; data.len()];
    // let mut rng: Xoshiro256StarStar = SeedableRng::seed_from_u64(data.len() as u64 * 23);
    let (border, datasize) = (label.len(), data.len() as f64);
    let max_entropy = (cluster_num as f64).ln() * datasize;
    let mut ws: Vec<f64> = (0..cluster_num)
        .map(|i| weights_of_reads.iter().map(|g| g[i]).sum::<f64>() / datasize)
        .collect();
    // updates_flags(&mut updates, &weights_of_reads, &mut rng, pick_prob, border);
    'outer: loop {
        for _ in 0..pick_up_len {
            let before_soe = weights_of_reads.iter().map(|e| entropy(e)).sum::<f64>();
            let wor = &mut weights_of_reads;
            update_weights(wor, &mut ws, border, &data, &models, &updates, beta, config);
            // updates_flags(&mut updates, &weights_of_reads, &mut rng, pick_prob, border);
            models.iter_mut().enumerate().for_each(|(cluster, model)| {
                mf.update_model(&weights_of_reads, &updates, &data, cluster, model, config);
            });
            let soe = weights_of_reads.iter().map(|e| entropy(e)).sum::<f64>();
            let rate = soe / max_entropy;
            let diff = (before_soe - soe) / max_entropy;
            if rate < ENTROPY_THR || (diff < ENTROPY_THR && beta >= 1.0) {
                break 'outer;
            }
            let wr = &weights_of_reads;
            report(id, wr, border, answer, &ws, beta);
        }
        for ms in models.iter() {
            for (idx, m) in ms.iter().enumerate() {
                debug!("{}\t{}", idx, m);
            }
        }
        beta = (beta * beta_step).min(1.);
    }
    let wr = &weights_of_reads;
    report(id, wr, border, answer, &ws, beta);
    weights_of_reads
}

fn search_initial_beta(
    data: &[Read],
    label: &[u8],
    forbidden: &[Vec<u8>],
    cluster_num: usize,
    c: &Config,
    chain_len: usize,
) -> f64 {
    let seed = data.iter().map(|e| e.len()).sum::<usize>() as u64;
    let mut rng: Xoshiro256StarStar = SeedableRng::seed_from_u64(seed);
    let weight_of_read = construct_initial_weights(label, forbidden, cluster_num, data.len(), seed);
    let border = label.len();
    let mut mf = ModelFactory::new(chain_len, data);
    let mut beta = INIT_BETA;
    let soe = weight_of_read.iter().map(|e| entropy(e)).sum::<f64>();
    let wor = &weight_of_read;
    loop {
        let next = soe_after_sampling(beta, data, wor, border, &mut rng, cluster_num, &mut mf, c);
        debug!("SEARCH\t{:.3}\t{:.3}\t{:.3}", beta, next, soe - next);
        if soe > next {
            break;
        } else {
            beta *= FACTOR;
        }
    }
    beta /= FACTOR;
    loop {
        let next = soe_after_sampling(beta, data, wor, border, &mut rng, cluster_num, &mut mf, c);
        debug!("SEARCH\t{:.3}\t{:.3}\t{:.3}", beta, next, soe - next);
        if soe < next {
            break;
        } else {
            beta /= FACTOR;
        }
    }
    beta
}

fn soe_after_sampling<R: Rng>(
    beta: f64,
    data: &[Read],
    wor: &[Vec<f64>],
    border: usize,
    _rng: &mut R,
    cluster_num: usize,
    mf: &mut ModelFactory,
    c: &Config,
) -> f64 {
    let datasize = data.len() as f64;
    let updates = vec![false; data.len()];
    let mut ws: Vec<f64> = (0..cluster_num)
        .map(|i| wor.iter().map(|g| g[i]).sum::<f64>() / datasize)
        .collect();
    let mut models: Vec<Vec<POA>> = (0..cluster_num)
        .map(|cl| mf.generate_model(&wor, data, cl, c))
        .collect();
    let mut wor = wor.to_vec();
    //for _ in 0..SAMPLING {
    for _ in 0..10 {
        // updates_flags(&mut updates, &wor, rng, INIT_PICK_PROB, border);
        update_weights(&mut wor, &mut ws, border, data, &models, &updates, beta, c);
        models.iter_mut().enumerate().for_each(|(cluster, model)| {
            mf.update_model(&wor, &updates, data, cluster, model, c);
        });
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

// fn updates_flags<R: Rng>(
//     updates: &mut [bool],
//     weight_of_read: &[Vec<f64>],
//     rng: &mut R,
//     pick_prob: f64,
//     border: usize,
// ) {
//     updates.iter_mut().for_each(|e| *e = false);
//     let total = ((weight_of_read.len() - border) as f64 * pick_prob).floor() as usize;
//     let mut current = 0;
//     let max = weight_of_read.len();
//     while current < total {
//         let i = rng.gen_range(border, max);
//         current += (!updates[i]) as usize;
//         updates[i] = true;
//     }
// }

fn update_weights(
    weight_of_read: &mut [Vec<f64>],
    ws: &mut [f64],
    border: usize,
    data: &[Read],
    models: &[Vec<POA>],
    updates: &[bool],
    beta: f64,
    c: &Config,
) {
    data.par_iter()
        .zip(weight_of_read.par_iter_mut())
        .zip(updates.par_iter())
        .skip(border)
        // .filter(|&(_, &b)| b)
        .for_each(|((read, weights), _)| {
            compute_log_probs(models, &ws, read, weights, c, beta);
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
    beta: f64,
) {
    assert_eq!(models.len(), ws.len());
    assert_eq!(models.len(), weights.len());
    assert!((1. - ws.iter().sum::<f64>()).abs() < 0.001);
    models
        .par_iter()
        .zip(weights.par_iter_mut())
        .zip(ws.par_iter())
        .for_each(|((ms, g), w)| {
            *g = read
                .par_iter()
                .map(|&(pos, ref u)| {
                    let model = &ms[pos];
                    beta * model.forward(u, c)
                })
                .sum::<f64>()
                + w.ln();
        });
}
