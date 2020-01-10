extern crate dbg_hmm;
extern crate libc;
use dbg_hmm::{Config, DBGHMM};

#[no_mangle]
pub extern "C" fn dbg_hmm(
    seqs: *const *const libc::c_char,
    lengths: *const libc::size_t,
    weights: *const libc::c_double,
    num_seqs: libc::size_t,
    k: libc::size_t,
) -> *mut DBGHMM {
    let seqs: Vec<&[u8]> = unsafe {
        let lengths: &[libc::size_t] = std::slice::from_raw_parts(lengths, num_seqs);
        let seqs: &[*const libc::c_char] = std::slice::from_raw_parts(seqs, num_seqs);
        assert_eq!(seqs.len(), lengths.len());
        seqs.into_iter()
            .zip(lengths)
            .map(|(&s, &l)| std::slice::from_raw_parts(s as *const u8, l))
            .collect()
    };
    let weights: &[f64] = unsafe { std::slice::from_raw_parts(weights, num_seqs) };
    let res = DBGHMM::new_with_weight_prior(&seqs, weights, k);
    let res = Box::into_raw(Box::new(res));
    res
}

#[no_mangle]
pub extern "C" fn default_config() -> *mut Config {
    let p = dbg_hmm::DEFAULT_CONFIG;
    Box::into_raw(Box::new(p))
}

use libc::{c_char, size_t};
#[no_mangle]
pub extern "C" fn forward(
    ptr: *const DBGHMM,
    query: *const c_char,
    length: size_t,
    config: *const Config,
) -> f64 {
    let query = unsafe { std::slice::from_raw_parts(query as *const u8, length) };
    let model = if ptr.is_null() {
        eprintln!("ERROR!!!!");
        return 0.;
    } else {
        unsafe { &*ptr }
    };
    let config = if config.is_null() {
        eprintln!("ERROR!!!!");
        return 0.;
    } else {
        unsafe { &*config }
    };
    model.forward(query, config)
}

#[no_mangle]
pub extern "C" fn viterbi(
    ptr: *const DBGHMM,
    query: *const c_char,
    length: size_t,
    config: *const Config,
    seqs: *mut *const libc::c_char,
    states: *mut *const libc::c_char,
    result_length: *mut libc::size_t,
) -> libc::c_double {
    let query = unsafe { std::slice::from_raw_parts(query as *const u8, length) };
    let model = if ptr.is_null() {
        eprintln!("ERROR!!!!");
        return 0.;
    } else {
        unsafe { &*ptr }
    };
    let config = if config.is_null() {
        eprintln!("ERROR!!!!");
        return 0.;
    } else {
        unsafe { &*config }
    };
    let (score, seq) = model.viterbi(query, config);
    unsafe {
        *result_length = seq.len();
        let viterbi_path: Vec<_> = seq.iter().map(|e| e.0 as i8).collect();
        *seqs = viterbi_path.as_ptr();
        std::mem::forget(viterbi_path);
        let viterbi_states: Vec<_> = seq.iter().map(|e| e.1 as i8).collect();
        *states = viterbi_states.as_ptr();
        std::mem::forget(viterbi_states);
    }
    score
}

#[no_mangle]
pub extern "C" fn free(ptr: *mut DBGHMM) {
    let _t = if !ptr.is_null() {
        unsafe { Box::from_raw(ptr) }
    } else {
        return;
    };
}

#[cfg(test)]
mod tests {
    #[test]
    fn it_works() {
        assert_eq!(2 + 2, 4);
    }
}
