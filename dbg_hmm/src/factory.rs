use rand::{Rng, SeedableRng};
use rand_xoshiro::Xoshiro256StarStar;
use std::collections::HashMap;

use super::{find_union, Kmer, DBGHMM, SCALE, THR, THR_ON};

#[derive(Default)]
pub struct Factory {
    inner: HashMap<Vec<u8>, (f64, usize)>,
    is_safe: Vec<bool>,
    temp_index: Vec<usize>,
    edges: Vec<Vec<usize>>,
    fu: find_union::FindUnion,
}

impl Factory {
    fn len(&self) -> usize {
        self.inner.len()
    }
    fn clear(&mut self) {
        self.inner.clear();
        self.temp_index.clear();
        self.is_safe.clear();
        self.edges.clear();
        self.fu.clear();
    }
    fn is_empty(&self) -> bool {
        self.inner.is_empty()
            && self.is_safe.is_empty()
            && self.temp_index.is_empty()
            && self.edges.is_empty()
            && self.fu.is_empty()
    }
    pub fn generate(&mut self, dataset: &[Vec<u8>], k: usize) -> DBGHMM {
        // let counter = &mut self.inner;
        dataset.iter().for_each(|seq| {
            seq.windows(k)
                .for_each(|kmer| match self.inner.get_mut(kmer) {
                    Some(x) => x.0 += 1.,
                    None => {
                        self.inner.insert(kmer.to_vec(), (1., std::usize::MAX));
                    }
                })
        });
        let thr = self.calc_thr();
        let mut nodes = Vec::with_capacity(self.len());
        for seq in dataset.iter().filter(|seq| seq.len() > k) {
            for x in seq.windows(k + 1) {
                if let Some((from, to)) = self.get_indices(x, &mut nodes, thr, k) {
                    nodes[from].push_edge_with(x[k], to);
                }
            }
            if let Some(x) = self.inner.get(&seq[seq.len() - k..]) {
                if x.0 > thr && x.1 != std::usize::MAX {
                    nodes[x.1].is_tail = true;
                }
            }
            if let Some(x) = self.inner.get(&seq[..k]) {
                if x.0 > thr && x.1 != std::usize::MAX {
                    nodes[x.1].is_head = true;
                }
            }
        }
        let nodes = self.clean_up_nodes(nodes);
        let weight = dataset.len() as f64;
        self.clear();
        DBGHMM { nodes, k, weight }
    }
    pub fn generate_from_ref(&mut self, dataset: &[&[u8]], k: usize) -> DBGHMM {
        assert!(self.is_empty());
        dataset.iter().for_each(|seq| {
            seq.windows(k)
                .for_each(|kmer| match self.inner.get_mut(kmer) {
                    Some(x) => x.0 += 1.,
                    None => {
                        self.inner.insert(kmer.to_vec(), (1., std::usize::MAX));
                    }
                })
        });
        let thr = self.calc_thr();
        let mut nodes = Vec::with_capacity(self.len());
        for seq in dataset {
            for x in seq.windows(k + 1) {
                if let Some((from, to)) = self.get_indices(x, &mut nodes, thr, k) {
                    nodes[from].push_edge_with(x[k], to);
                }
            }
            if let Some(x) = self.inner.get(&seq[seq.len() - k..]) {
                if x.0 > thr && x.1 != std::usize::MAX {
                    nodes[x.1].is_tail = true;
                }
            }
            if let Some(x) = self.inner.get(&seq[..k]) {
                if x.0 > thr && x.1 != std::usize::MAX {
                    nodes[x.1].is_head = true;
                }
            }
        }
        let nodes = self.clean_up_nodes(nodes);
        let weight = dataset.len() as f64;
        self.clear();
        DBGHMM { nodes, k, weight }
    }
    pub fn generate_with_weight(&mut self, dataset: &[&[u8]], ws: &[f64], k: usize) -> DBGHMM {
        assert!(self.is_empty());
        let mut rng: Xoshiro256StarStar = SeedableRng::seed_from_u64(dataset.len() as u64);
        assert_eq!(dataset.len(), ws.len());
        let ep = 0.0001;
        dataset
            .iter()
            .zip(ws.iter())
            .filter(|&(_, &w)| w > ep)
            .for_each(|(seq, w)| {
                seq.windows(k)
                    .for_each(|kmer| match self.inner.get_mut(kmer) {
                        Some(x) => x.0 += w,
                        None => {
                            self.inner.insert(kmer.to_vec(), (*w, std::usize::MAX));
                        }
                    })
            });
        let thr = self.calc_thr_weight();
        let mut nodes = Vec::with_capacity(1_000);
        for (seq, &w) in dataset.iter().zip(ws.iter()).filter(|&(_, &w)| w > ep) {
            for x in seq.windows(k + 1) {
                if let Some((from, to)) = self.get_indices_exp(x, &mut nodes, thr, k, &mut rng) {
                    nodes[from].push_edge_with_weight(x[k], to, w);
                }
            }
            if let Some(x) = self.inner.get(&seq[seq.len() - k..]) {
                if x.1 != std::usize::MAX {
                    nodes[x.1].is_tail = true;
                }
            }
            if let Some(x) = self.inner.get(&seq[..k]) {
                if x.1 != std::usize::MAX {
                    nodes[x.1].is_head = true;
                }
            }
        }
        let weight = ws.iter().sum::<f64>();
        if weight < 1.0001 {
            self.clear();
            let nodes = vec![];
            return DBGHMM { nodes, k, weight };
        }
        let thr = nodes.iter().map(|e| e.tot).sum::<f64>() / nodes.len() as f64 - THR;
        let nodes = match self.clean_up_nodes_exp(nodes, thr, &mut rng) {
            Some(res) => res,
            None => panic!(
                "thr:{},weight:{},raw_kmers:{},wsdump:{:#?}",
                thr,
                weight,
                self.inner.len(),
                ws
            ),
        };
        self.clear();
        DBGHMM { nodes, k, weight }
    }
    fn clean_up_nodes(&mut self, nodes: Vec<Kmer>) -> Vec<Kmer> {
        let mut nodes = self.renaming_nodes(nodes);
        nodes.iter_mut().for_each(Kmer::finalize);
        nodes
    }
    fn clean_up_nodes_exp<R: Rng>(
        &mut self,
        nodes: Vec<Kmer>,
        thr: f64,
        rng: &mut R,
    ) -> Option<Vec<Kmer>> {
        let mut buffer: Vec<_> = Vec::with_capacity(nodes.len());
        let nodes = self.cut_lightweight_loop(nodes, thr, rng, &mut buffer);
        assert!(buffer.is_empty());
        let nodes = self.cut_tip(nodes, 2, &mut buffer);
        assert!(buffer.is_empty());
        let nodes = self.pick_largest_components(nodes, &mut buffer)?;
        let mut nodes = self.renaming_nodes(nodes);
        nodes.iter_mut().for_each(Kmer::finalize);
        Some(nodes)
    }
    fn pick_largest_components(
        &mut self,
        mut nodes: Vec<Kmer>,
        buffer: &mut Vec<Kmer>,
    ) -> Option<Vec<Kmer>> {
        assert!(self.temp_index.is_empty());
        assert!(self.is_safe.is_empty());
        assert!(self.fu.is_empty());
        self.fu.refresh(nodes.len());
        for (from, n) in nodes.iter().enumerate() {
            for &to in n.edges.iter().filter_map(|e| e.as_ref()) {
                self.fu.unite(from, to);
            }
        }
        let max_group_and_size = (0..nodes.len())
            .filter_map(|e| {
                let parent = self.fu.find(e).unwrap();
                if parent == e {
                    self.fu.size(e).map(|s| (e, s))
                } else {
                    None
                }
            })
            .max_by_key(|&(_, size)| size);
        let (max_group, _) = match max_group_and_size {
            Some(res) => res,
            None => {
                error!("No components:{:?}", nodes);
                return None;
            }
        };
        let mut index = 0;
        for i in 0..nodes.len() {
            self.temp_index.push(index);
            index += (self.fu.find(i).unwrap() == max_group) as usize;
        }
        let mut index = nodes.len();
        while let Some(mut node) = nodes.pop() {
            index -= 1;
            if self.fu.find(index).unwrap() == max_group {
                node.remove_if_not(&mut self.fu, max_group);
                node.rename_by(&self.temp_index);
                buffer.push(node);
            }
        }
        buffer.reverse();
        std::mem::swap(buffer, &mut nodes);
        self.temp_index.clear();
        self.fu.clear();
        Some(nodes)
    }
    #[allow(dead_code)]
    fn select_supported_node(edges: &[Vec<usize>], nodes: usize, is_safe: &[bool]) -> Vec<bool> {
        let mut is_supported = vec![false; nodes];
        for (from, edges) in edges.iter().enumerate() {
            for &to in edges.iter() {
                is_supported[to] |= is_safe[from];
                is_supported[from] |= is_safe[to];
            }
        }
        is_supported
    }
    fn cut_lightweight_loop<R: Rng>(
        &mut self,
        mut nodes: Vec<Kmer>,
        thr: f64,
        rng: &mut R,
        buffer: &mut Vec<Kmer>,
    ) -> Vec<Kmer> {
        self.is_safe.clear();
        self.temp_index.clear();
        for node in &nodes {
            self.temp_index
                .push(!rng.gen_bool(Self::prob(node.kmer_weight - thr)) as usize)
        }
        (0..nodes.len()).for_each(|_| self.is_safe.push(false));
        for (from, ref node) in nodes.iter().enumerate() {
            for &to in node.edges.iter().filter_map(|e| e.as_ref()) {
                self.is_safe[to] |= self.temp_index[from] == 1;
                self.is_safe[from] |= self.temp_index[to] == 1;
            }
        }
        self.temp_index.clear();
        let mut index = 0;
        for &b in &self.is_safe {
            self.temp_index.push(index);
            index += b as usize;
        }
        assert!(buffer.is_empty());
        let mut idx = nodes.len();
        while let Some(mut node) = nodes.pop() {
            idx -= 1;
            if self.is_safe[idx] {
                node.remove_if_not_supported(&self.is_safe);
                node.rename_by(&self.temp_index);
                buffer.push(node);
            }
        }
        buffer.reverse();
        self.is_safe.clear();
        self.temp_index.clear();
        std::mem::swap(&mut nodes, buffer);
        nodes
    }
    // Cut tips. We can assume that the graph is a
    // connected graph or a tree.
    fn cut_tip(&mut self, mut nodes: Vec<Kmer>, t: usize, buffer: &mut Vec<Kmer>) -> Vec<Kmer> {
        assert!(self.is_safe.is_empty());
        assert!(self.edges.is_empty());
        let pseudo_node = nodes.len();
        (0..nodes.len() + 1).for_each(|_| self.edges.push(vec![]));
        nodes.iter().for_each(|e| self.is_safe.push(e.is_tail));
        for (idx, node) in nodes.iter().enumerate() {
            for &to in node.edges.iter().filter_map(|e| e.as_ref()) {
                self.edges[idx].push(to);
                self.is_safe[idx] |= nodes[to].is_tail;
            }
            if node.is_tail {
                self.edges[idx].push(pseudo_node);
            }
        }
        let is_supported_forward = self.cut_tip_inner(nodes.len() + 1, t);
        self.edges.iter_mut().for_each(|eds| eds.clear());
        nodes
            .iter()
            .zip(self.is_safe.iter_mut())
            .for_each(|(n, b)| *b |= n.is_head);
        for (idx, node) in nodes.iter().enumerate() {
            for &to in node.edges.iter().filter_map(|e| e.as_ref()) {
                self.edges[to].push(idx);
                self.is_safe[to] |= node.is_head;
            }
            if node.is_head {
                self.edges[idx].push(pseudo_node);
            }
        }
        let is_supported_backward = self.cut_tip_inner(nodes.len() + 1, t);
        let ave = nodes.iter().map(|e| e.kmer_weight).sum::<f64>() / nodes.len() as f64;
        nodes
            .iter()
            .zip(self.is_safe.iter_mut())
            .for_each(|(n, b)| *b &= n.kmer_weight > ave / 2.0);
        is_supported_forward
            .into_iter()
            .zip(is_supported_backward)
            .zip(self.is_safe.iter_mut())
            .for_each(|((forward, backward), is_safe)| *is_safe = (forward & backward) | *is_safe);
        assert!(self.temp_index.is_empty());
        let mut idx = 0;
        for &b in &self.is_safe {
            self.temp_index.push(idx);
            idx += b as usize;
        }
        assert!(buffer.is_empty());
        let mut idx = nodes.len();
        while let Some(mut node) = nodes.pop() {
            idx -= 1;
            if self.is_safe[idx] {
                node.remove_if_not_supported(&self.is_safe);
                node.rename_by(&self.temp_index);
                buffer.push(node)
            }
        }
        buffer.reverse();
        std::mem::swap(buffer, &mut nodes);
        self.is_safe.clear();
        self.temp_index.clear();
        self.edges.clear();
        buffer.clear();
        nodes
    }
    // Cut tips. We can assume that the graph is a
    // connected graph or a tree.
    fn cut_tip_inner(&mut self, nodes: usize, t: usize) -> Vec<bool> {
        let edges = &self.edges;
        let t = t as i32;
        let dist_to_root = Self::dist_to_root(edges, nodes);
        let bad_list: Vec<_> = (0..nodes)
            .filter(|&e| edges[e].len() > 1)
            .flat_map(|e| edges[e].iter().filter(|&&to| dist_to_root[to] <= t))
            .copied()
            .collect();
        // If arrived, true.
        let mut is_arrived: Vec<bool> = vec![false; nodes];
        for i in bad_list {
            if is_arrived[i] {
                continue;
            }
            let mut stack = vec![i];
            'dfs: while !stack.is_empty() {
                let node = *stack.last().unwrap();
                is_arrived[node] |= true;
                for &to in &edges[node] {
                    if !is_arrived[to] {
                        stack.push(to);
                        continue 'dfs;
                    }
                }
                stack.pop().unwrap();
            }
        }
        is_arrived.iter_mut().for_each(|e| *e = !*e);
        is_arrived
    }
    fn dist_to_root(edges: &[Vec<usize>], nodes: usize) -> Vec<i32> {
        // 0 -> never arrived
        // 1 -> active(arrived, being traversed currently)
        // 2 -> inactive(arrived, have been traversed)
        let mut status = vec![0; nodes];
        let mut dist_to_root: Vec<i32> = vec![-1; nodes];
        for i in 0..nodes {
            if status[i] != 0 {
                continue;
            }
            let mut stack = vec![i];
            'dfs: while !stack.is_empty() {
                let node = *stack.last().unwrap();
                if status[node] == 0 {
                    status[node] = 1;
                }
                for &to in &edges[node] {
                    if status[to] == 0 {
                        stack.push(to);
                        continue 'dfs;
                    } else if status[to] == 1 {
                        // Loop detected.
                        dist_to_root[node] = 100000;
                    }
                }
                let last = stack.pop().unwrap();
                dist_to_root[last] = edges[node]
                    .iter()
                    .map(|&to| dist_to_root[to])
                    .max()
                    .unwrap_or(-1)
                    .max(dist_to_root[last])
                    + 1;
                // Deactivate
                status[last] = 2;
            }
        }
        dist_to_root
    }
    fn calc_thr(&self) -> f64 {
        let ave = self.inner.values().map(|e| e.0).sum::<f64>() / self.inner.len() as f64;
        if THR_ON {
            ave - THR
        } else {
            0.
        }
    }
    fn calc_thr_weight(&self) -> f64 {
        let (sum, denom) =
            self.inner.values().fold(
                (0., 0.),
                |(x, y), &(w, _)| if w >= 1. { (x + w, y + 1.) } else { (x, y) },
            );
        let ave = sum / denom;
        if THR_ON {
            // y = x - 1 - epsilon
            ave - THR
        } else {
            0.
        }
    }
    fn get_idx(&mut self, kmer: &[u8], nodes: &mut Vec<Kmer>, thr: f64) -> Option<usize> {
        match self.inner.get_mut(kmer) {
            Some(x) if x.0 <= thr => None,
            Some(x) if x.1 != std::usize::MAX => Some(x.1),
            Some(x) => {
                x.1 = nodes.len();
                nodes.push(Kmer::new(kmer, x.0));
                Some(nodes.len() - 1)
            }
            _ => unreachable!(),
        }
    }
    fn prob(x: f64) -> f64 {
        // if x.is_sign_positive(){
        //     0.
        // }else{
        ((SCALE * x).exp() + 1.).recip()
        // }
    }
    fn get_idx_exp<R: Rng>(
        &mut self,
        kmer: &[u8],
        nodes: &mut Vec<Kmer>,
        thr: f64,
        rng: &mut R,
    ) -> Option<usize> {
        match self.inner.get_mut(kmer) {
            Some(x) if rng.gen_bool(Self::prob(x.0 - thr)) => None,
            Some(x) if x.1 != std::usize::MAX => Some(x.1),
            Some(x) => {
                x.1 = nodes.len();
                nodes.push(Kmer::new(kmer, x.0));
                Some(nodes.len() - 1)
            }
            _ => unreachable!(),
        }
    }

    fn get_indices(
        &mut self,
        x: &[u8],
        nodes: &mut Vec<Kmer>,
        thr: f64,
        k: usize,
    ) -> Option<(usize, usize)> {
        let (prev, after) = (&x[..k], &x[1..]);
        let from = self.get_idx(prev, nodes, thr)?;
        let to = self.get_idx(after, nodes, thr)?;
        Some((from, to))
    }
    fn get_indices_exp<R: Rng>(
        &mut self,
        x: &[u8],
        nodes: &mut Vec<Kmer>,
        thr: f64,
        k: usize,
        rng: &mut R,
    ) -> Option<(usize, usize)> {
        let (prev, next) = (&x[..k], &x[1..]);
        let from = self.get_idx_exp(prev, nodes, thr, rng)?;
        let next = self.get_idx_exp(next, nodes, thr, rng)?;
        Some((from, next))
    }

    // Return topological sorted order. If there's a cycle,
    // it neglect the back-edge, and proceed with raise error messages to stderr(DEBUG MODE).
    fn topological_sort(nodes: &[Kmer]) -> Vec<usize> {
        let mut edges: Vec<Vec<_>> = vec![vec![]; nodes.len()];
        for (idx, node) in nodes.iter().enumerate() {
            for to in node.edges.iter().filter_map(|e| e.as_ref()) {
                edges[idx].push(*to);
            }
        }
        match Self::topological_sort_inner(&edges, nodes.len()) {
            Ok(res) => res,
            Err(res) => res,
        }
    }
    // Topological sorting.
    fn topological_sort_inner(
        edges: &[Vec<usize>],
        nodes: usize,
    ) -> Result<Vec<usize>, Vec<usize>> {
        // 0 -> never arrived
        // 1 -> active(arrived, being traversed currently)
        // 2 -> inactive(arrived, have been traversed)
        let mut status = vec![0; nodes];
        let mut order = vec![];
        let mut is_dag = true;
        for i in 0..nodes {
            if status[i] != 0 {
                continue;
            }
            let mut stack = vec![i];
            'dfs: while !stack.is_empty() {
                let node = *stack.last().unwrap();
                if status[node] == 0 {
                    // preorder
                    status[node] = 1;
                }
                for &to in &edges[node] {
                    if status[to] == 0 {
                        stack.push(to);
                        continue 'dfs;
                    } else if status[to] == 1 {
                        is_dag = false;
                    }
                }
                // No-op
                let last = stack.pop().unwrap();
                order.push(last);
                // Deactivate
                status[last] = 2;
            }
        }
        order.reverse();
        if is_dag {
            Ok(order)
        } else {
            Err(order)
        }
    }
    fn renaming_nodes(&mut self, mut nodes: Vec<Kmer>) -> Vec<Kmer> {
        assert!(self.temp_index.is_empty());
        (0..nodes.len()).for_each(|_| self.temp_index.push(0));
        let order = Self::topological_sort(&nodes);
        for (order, v) in order.into_iter().enumerate() {
            self.temp_index[v] = order;
        }
        let mut result = vec![None; nodes.len()];
        let mut idx = nodes.len();
        while let Some(mut node) = nodes.pop() {
            idx -= 1;
            node.rename_by(&self.temp_index);
            result[self.temp_index[idx]] = Some(node);
        }
        assert!(result.iter().all(|e| e.is_some()));
        assert!(nodes.is_empty());
        self.temp_index.clear();
        nodes.extend(result.into_iter().filter_map(|e| e));
        nodes
    }
    pub fn new() -> Self {
        let inner = std::collections::HashMap::new();
        let is_safe = vec![];
        let temp_index = vec![];
        let edges = vec![];
        let fu = find_union::FindUnion::new(0);
        Self {
            inner,
            is_safe,
            temp_index,
            edges,
            fu,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn tip_cut() {
        let edges = vec![
            vec![1, 2],
            vec![],
            vec![],
            vec![4],
            vec![0, 5],
            vec![9],
            vec![7],
            vec![8],
            vec![11],
            vec![6, 7, 10],
            vec![],
            vec![],
        ];
        let nodes = 12;
        assert_eq!(edges.len(), 12);
        let dist_to_root = Factory::dist_to_root(&edges, nodes);
        let answer = vec![1, 0, 0, 7, 6, 5, 3, 2, 1, 4, 0, 0];
        assert_eq!(dist_to_root, answer);
        // let is_supported = Factory::cut_tip_inner(&edges, nodes, 1 );
        // let answer = vec![
        //     false, false, false, true, true, true, true, true, true, true, false, true,
        // ];
        // assert_eq!(is_supported, answer);
        let edges = vec![
            vec![1],
            vec![2, 3],
            vec![],
            vec![4],
            vec![],
            vec![6],
            vec![7],
            vec![0, 8],
            vec![5],
        ];
        let nodes = 9;
        assert_eq!(edges.len(), nodes);
        let dist_to_root = Factory::dist_to_root(&edges, nodes);
        let answer = vec![3, 2, 0, 1, 0];
        eprintln!("{:?}", dist_to_root);
        for i in 0..5 {
            assert_eq!(dist_to_root[i], answer[i]);
        }
        for i in 5..nodes {
            assert!(dist_to_root[i] > 100);
        }
        // let is_supported = Factory::cut_tip_inner(&edges, nodes, 1);
        // let answer = vec![true, true, false, false, false, true, true, true, true];
        // assert_eq!(answer, is_supported);
    }
    #[test]
    fn select_supported_test() {
        let edges = vec![
            vec![1],
            vec![2],
            vec![3, 7],
            vec![],
            vec![5],
            vec![],
            vec![7],
            vec![4],
            vec![6],
            vec![8, 0],
            vec![9],
        ];
        let nodes = 11;
        let is_safe = vec![
            false, true, true, false, true, true, false, false, false, true, true,
        ];
        let supported = Factory::select_supported_node(&edges, nodes, &is_safe);
        let mut answer = vec![true; nodes];
        answer[6] = false;
        assert_eq!(supported, answer);
    }
    #[test]
    fn topological_sorting() {
        let edges = vec![
            vec![1],
            vec![2],
            vec![3],
            vec![4],
            vec![],
            vec![6],
            vec![2],
            vec![8],
            vec![9],
            vec![3],
        ];
        let order = Factory::topological_sort_inner(&edges, edges.len()).unwrap();
        assert_eq!(order, vec![7, 8, 9, 5, 6, 0, 1, 2, 3, 4]);
        let edges = vec![
            vec![1],
            vec![2],
            vec![3],
            vec![4, 8],
            vec![],
            vec![6],
            vec![2],
            vec![8],
            vec![9],
            vec![3],
        ];
        let order = Factory::topological_sort_inner(&edges, edges.len()).unwrap_err();
        assert_eq!(order, vec![7, 5, 6, 0, 1, 2, 3, 8, 9, 4]);
        let edges = vec![
            vec![1],
            vec![2, 7],
            vec![3],
            vec![4],
            vec![5],
            vec![6],
            vec![],
            vec![8],
            vec![5],
        ];
        let order = Factory::topological_sort_inner(&edges, edges.len()).unwrap();
        assert_eq!(order, vec![0, 1, 7, 8, 2, 3, 4, 5, 6]);
    }
}
