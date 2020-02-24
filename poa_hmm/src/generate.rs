use crate::PartialOrderAlignment;
impl PartialOrderAlignment {
    pub fn consensus(&self) -> Option<Vec<u8>> {
        let mut node = self
            .nodes
            .iter()
            .filter(|e| e.is_head)
            .max_by(|a, b| a.head_weight().partial_cmp(&b.head_weight()).unwrap())?;
        let mut res = vec![node.base()];
        while node.has_edge() {
            let (idx, _) = node
                .edges()
                .iter()
                .zip(node.weights().iter())
                .max_by(|a, b| (a.1).partial_cmp(&(b.1)).unwrap())?;
            node = &self.nodes[*idx];
            res.push(node.base());
        }
        Some(res)
    }
}
