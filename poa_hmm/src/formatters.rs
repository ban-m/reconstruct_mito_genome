use super::POA;

use std::fmt;
impl fmt::Debug for POA {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let result: Vec<_> = self
            .nodes
            .iter()
            .enumerate()
            .map(|(idx, e)| format!("{}\t{:?}", idx, e))
            .collect();
        write!(f, "{}", result.join("\n\n"))
    }
}

impl fmt::Display for POA {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let edges = self
            .nodes
            .iter()
            .map(|n| n.edges.iter().count())
            .sum::<usize>();
        write!(f, "{}\t{}\t{:.3}", self.nodes.len(), edges, self.weight)
    }
}
