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

pub const PACBIO_CONFIG: Config = Config {
    mismatch: 0.0341,
    base_freq: [0.28, 0.22, 0.22, 0.28],
    p_match: 0.9124,
    p_ins: 0.0606,
    p_del: 0.0270,
    p_extend_ins: 0.3014,
    p_extend_del: 0.0695,
    p_del_to_ins: 0.0086,
};


/// A configure struct.
#[derive(Debug, Clone, Default)]
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

impl Config {
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        mismatch: f64,
        p_match: f64,
        p_ins: f64,
        p_del: f64,
        p_extend_ins: f64,
        p_extend_del: f64,
        p_del_to_ins: f64,
        base_freq: [f64; 4],
    ) -> Self {
        Self {
            mismatch,
            base_freq,
            p_match,
            p_ins,
            p_del,
            p_extend_ins,
            p_extend_del,
            p_del_to_ins,
        }
    }
}
