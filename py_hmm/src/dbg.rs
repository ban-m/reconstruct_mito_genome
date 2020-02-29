use dbg_hmm::DBGHMM;
/// Creates de Bruijn HMM from the given sequence and returns the result.
/// # Safety
/// This function is unsafe if invalid pointers are given.
#[no_mangle]
pub unsafe extern "C" fn dbg_hmm(
    seqs: *const *const libc::c_char,
    lengths: *const libc::size_t,
    weights: *const libc::c_double,
    num_seqs: libc::size_t,
    k: libc::size_t,
) -> *mut DBGHMM {
    let seqs: Vec<&[u8]> = {
        let lengths: &[libc::size_t] = std::slice::from_raw_parts(lengths, num_seqs);
        let seqs: &[*const libc::c_char] = std::slice::from_raw_parts(seqs, num_seqs);
        assert_eq!(seqs.len(), lengths.len());
        seqs.iter()
            .zip(lengths)
            .map(|(&s, &l)| std::slice::from_raw_parts(s as *const u8, l))
            .collect()
    };
    let weights: &[f64] = std::slice::from_raw_parts(weights, num_seqs);
    let res = DBGHMM::new_with_weight(&seqs, weights, k);
    Box::into_raw(Box::new(res))
}

/// Forward algorithm. Return Likelihood of the input.
#[no_mangle]
pub extern "C" fn forward_dbg(
    ptr: *const DBGHMM,
    query: *const libc::c_char,
    length: libc::size_t,
    config: *const crate::Config,
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
        dbg_hmm::Config::from(unsafe { *config })
    };
    model.forward(query, &config)
}

/// Free the pointer.
#[no_mangle]
pub extern "C" fn free_dbg(ptr: *mut DBGHMM) {
    let _t = if !ptr.is_null() {
        unsafe { Box::from_raw(ptr) }
    } else {
        return;
    };
}
