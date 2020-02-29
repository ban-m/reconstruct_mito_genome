use poa_hmm::POA;
/// Creates Partial Order alignment from the given sequence and returns the result.
#[no_mangle]
pub extern "C" fn poa_hmm(
    seqs: *const *const libc::c_char,
    lengths: *const libc::size_t,
    weights: *const libc::c_double,
    num_seqs: libc::size_t,
    config: *const crate::Config,
) -> *mut POA {
    let seqs: Vec<&[u8]> = unsafe {
        let lengths: &[libc::size_t] = std::slice::from_raw_parts(lengths, num_seqs);
        let seqs: &[*const libc::c_char] = std::slice::from_raw_parts(seqs, num_seqs);
        assert_eq!(seqs.len(), lengths.len());
        seqs.iter()
            .zip(lengths)
            .map(|(&s, &l)| std::slice::from_raw_parts(s as *const u8, l))
            .collect()
    };
    let weights: &[f64] = unsafe { std::slice::from_raw_parts(weights, num_seqs) };
    let config = if config.is_null() {
        eprintln!("ERROR!!!!");
        let res = POA::default();
        return Box::into_raw(Box::new(res));
    } else {
        poa_hmm::Config::from(unsafe { *config })
    };
    let res = POA::generate(&seqs, weights, &config);
    Box::into_raw(Box::new(res))
}

/// Forward algorithm. Return Likelihood of the input.
#[no_mangle]
pub extern "C" fn forward_poa(
    ptr: *const POA,
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
        poa_hmm::Config::from(unsafe { *config })
    };
    model.forward(query, &config)
}

/// Forward algorithm. Return Likelihood of the input.
#[no_mangle]
pub extern "C" fn consensus_poa(ptr: *const POA, result: *mut libc::c_char, len: libc::size_t) {
    let model = if ptr.is_null() {
        eprintln!("ERROR!!!!");
        return;
    } else {
        unsafe { &*ptr }
    };
    let result = unsafe { std::slice::from_raw_parts_mut(result, len) };
    let consensus = model.consensus().unwrap_or_else(|| vec![]);
    for (x, &y) in result.iter_mut().zip(consensus.iter()) {
        *x = y as i8;
    }
}

#[no_mangle]
pub extern "C" fn consensus_poa_len(ptr: *const POA) -> libc::size_t {
    let model = if ptr.is_null() {
        eprintln!("ERROR!!!!");
        return 0;
    } else {
        unsafe { &*ptr }
    };
    model.consensus().unwrap_or_else(|| vec![]).len()
}

/// Free the pointer.
#[no_mangle]
pub extern "C" fn free_poa(ptr: *mut POA) {
    let _t = if !ptr.is_null() {
        unsafe { Box::from_raw(ptr) }
    } else {
        return;
    };
}
