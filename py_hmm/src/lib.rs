extern crate dbg_hmm;
extern crate libc;
extern crate poa_hmm;
mod dbg;
pub use dbg::*;
mod poa;
pub use poa::*;

#[derive(Debug, Clone, Default, Copy)]
pub struct Config {
    /// Mismatch probability at given position. # mismatch/(#mism + #match)
    pub mismatch: f64,
    /// The base composition of the read. A,C,G, and T.
    pub base_freq: [f64; 4],
    /// probability of matching at match states.
    pub p_match: f64,
    /// probability of starting insertion at match states.
    pub p_ins: f64,
    /// probability of starting deletion at match states.
    pub p_del: f64,
    /// probability of extending insertion at insertion state.
    pub p_extend_ins: f64,
    /// same for deletion
    pub p_extend_del: f64,
    /// probability of jump deletion state to insertion state.
    pub p_del_to_ins: f64,
}

pub const DEFAULT_CONFIG: Config = Config {
    mismatch: 0.03,
    base_freq: [0.25, 0.25, 0.25, 0.25],
    p_match: 0.89,
    p_ins: 0.06,
    p_del: 0.05,
    p_extend_ins: 0.06,
    p_extend_del: 0.05,
    p_del_to_ins: 0.06,
};

impl std::convert::From<Config> for poa_hmm::Config {
    fn from(c: Config) -> poa_hmm::Config {
        poa_hmm::Config {
            mismatch: c.mismatch,
            base_freq: c.base_freq,
            p_match: c.p_match,
            p_ins: c.p_ins,
            p_del: c.p_del,
            p_extend_ins: c.p_extend_ins,
            p_extend_del: c.p_del_to_ins,
            p_del_to_ins: c.p_del_to_ins,
        }
    }
}

impl std::convert::From<Config> for dbg_hmm::Config {
    fn from(c: Config) -> dbg_hmm::Config {
        dbg_hmm::Config {
            mismatch: c.mismatch,
            base_freq: c.base_freq,
            p_match: c.p_match,
            p_ins: c.p_ins,
            p_del: c.p_del,
            p_extend_ins: c.p_extend_ins,
            p_extend_del: c.p_del_to_ins,
            p_del_to_ins: c.p_del_to_ins,
        }
    }
}

#[no_mangle]
pub extern "C" fn default_config() -> *mut Config {
    let p = DEFAULT_CONFIG;
    Box::into_raw(Box::new(p))
}

/// # Safety
/// This function should not be called with invalid pointer.
/// Usually, all you want to do is to free a configuration
/// allocated by `default_config()` defined above.
#[no_mangle]
pub unsafe extern "C" fn free_config(ptr: *mut Config) {
    if !ptr.is_null() {
        Box::from_raw(ptr);
    }
}

// #[no_mangle]
// pub extern "C" fn viterbi(
//     ptr: *const DBGHMM,
//     query: *const c_char,
//     length: size_t,
//     config: *const dbg_hmm::Config,
//     seqs: *mut *const libc::c_char,
//     states: *mut *const libc::c_char,
//     result_length: *mut libc::size_t,
// ) -> libc::c_double {
//     let query = unsafe { std::slice::from_raw_parts(query as *const u8, length) };
//     let model = if ptr.is_null() {
//         eprintln!("ERROR!!!!");
//         return 0.;
//     } else {
//         unsafe { &*ptr }
//     };
//     let config = if config.is_null() {
//         eprintln!("ERROR!!!!");
//         return 0.;
//     } else {
//         unsafe { &*config }
//     };
//     let (score, seq) = model.viterbi(query, config);
//     unsafe {
//         *result_length = seq.len();
//         let viterbi_path: Vec<_> = seq.iter().map(|e| e.0 as i8).collect();
//         *seqs = viterbi_path.as_ptr();
//         std::mem::forget(viterbi_path);
//         let viterbi_states: Vec<_> = seq.iter().map(|e| e.1 as i8).collect();
//         *states = viterbi_states.as_ptr();
//         std::mem::forget(viterbi_states);
//     }
//     score
// }

#[cfg(test)]
mod tests {
    #[test]
    fn it_works() {
        assert_eq!(2 + 2, 4);
    }
}
