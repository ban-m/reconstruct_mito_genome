use crate::Base;
impl crate::PartialOrderAlignment {
    pub fn remove_node(mut self, thr_rank: usize) -> Self {
        if self.nodes.len() <= thr_rank {
            return self;
        }
        // eprintln!("Before:{},{}", self, thr_rank);
        let select = |node: &Base| node.weight();
        let thr = self.nodes.len().max(thr_rank) - thr_rank;
        let thr = select_nth_by(&self.nodes, thr, select).unwrap();
        let saved = self.clone();
        self = self.remove_lightweight(thr);
        self.trim_unreachable_nodes();
        if self.nodes.len() < 85 * thr_rank / 100 {
            // eprintln!("After:{},{}", self, thr_rank);
            return saved;
        }
        let saved = self.clone();
        let thr = {
            let median = select_nth_by(&self.nodes, self.nodes.len() / 2, select).unwrap();
            let deviation = |node: &Base| (median - node.weight()).abs();
            let mad = select_nth_by(&self.nodes, self.nodes.len() / 2, deviation).unwrap();
            median - 20. * mad.max(1.)
        };
        self = self.remove_lightweight(thr);
        self.trim_unreachable_nodes();
        if self.nodes.len() < 85 * thr_rank / 100 {
            // eprintln!("After:{},{}", self, thr_rank);
            return saved;
        }
        // eprintln!("Here:{},{}", self, thr_rank);
        self.nodes
            .iter_mut()
            .for_each(|node| node.remove_lightweight_edges());
        self.trim_unreachable_nodes();
        // eprintln!("After:{}", self);
        self
    }
    fn remove_lightweight(mut self, thr: f64) -> Self {
        let (_, mapping) = self
            .nodes
            .iter()
            .fold((0, vec![]), |(index, mut map), node| {
                let is_ok = node.weight() > thr;
                map.push((index, is_ok));
                (index + is_ok as usize, map)
            });
        let mut buffer = vec![];
        while let Some(mut node) = self.nodes.pop() {
            if node.weight() > thr {
                node.remove_if(&mapping);
                buffer.push(node);
            }
        }
        buffer.reverse();
        std::mem::swap(&mut buffer, &mut self.nodes);
        assert!(buffer.is_empty() && !self.nodes.is_empty());
        self
    }
    fn trim_unreachable_nodes(&mut self) {
        self.trim_unreachable_nodes_forward();
        self.trim_unreachable_nodes_reverse();
    }
    fn trim_unreachable_nodes_forward(&mut self) {
        let mut dfs_flag = vec![0; self.nodes.len()];
        let mut dfs_stack: Vec<_> = self
            .nodes
            .iter()
            .enumerate()
            .filter_map(|(idx, n)| if n.is_head { Some(idx) } else { None })
            .collect();
        assert!(!dfs_stack.is_empty());
        'dfs: while !dfs_stack.is_empty() {
            let node = *dfs_stack.last().unwrap();
            dfs_flag[node] = 1;
            for &to in self.nodes[node].edges.iter() {
                if dfs_flag[to] == 0 {
                    dfs_stack.push(to);
                    continue 'dfs;
                }
            }
            dfs_stack.pop().unwrap();
        }
        let (mapping, _) = dfs_flag.iter().fold((vec![], 0), |(mut map, index), &b| {
            map.push((index, b == 1));
            (map, index + b as usize)
        });
        let mut buffer = vec![];
        let mut idx = self.nodes.len();
        while let Some(mut node) = self.nodes.pop() {
            idx -= 1;
            if dfs_flag[idx] == 1 {
                node.remove_if(&mapping);
                buffer.push(node);
            }
        }
        buffer.reverse();
        assert!(self.nodes.is_empty());
        std::mem::swap(&mut self.nodes, &mut buffer);
    }
    fn trim_unreachable_nodes_reverse(&mut self) {
        let mut dfs_flag = vec![0; self.nodes.len()];
        let mut dfs_stack: Vec<_> = self
            .nodes
            .iter()
            .enumerate()
            .filter_map(|(idx, n)| if n.is_tail { Some(idx) } else { None })
            .collect();
        let edges = self.reverse_edges();
        'dfs: while !dfs_stack.is_empty() {
            let node = *dfs_stack.last().unwrap();
            dfs_flag[node] = 1;
            for &to in &edges[node] {
                if dfs_flag[to] == 0 {
                    dfs_stack.push(to);
                    continue 'dfs;
                }
            }
            dfs_stack.pop().unwrap();
        }
        let (mapping, _) = dfs_flag.iter().fold((vec![], 0), |(mut map, index), &b| {
            map.push((index, b == 1));
            (map, index + b as usize)
        });
        let mut buffer = vec![];
        let mut idx = self.nodes.len();
        while let Some(mut node) = self.nodes.pop() {
            idx -= 1;
            if dfs_flag[idx] == 1 {
                node.remove_if(&mapping);
                buffer.push(node);
            }
        }
        assert!(self.nodes.is_empty());
        buffer.reverse();
        std::mem::swap(&mut self.nodes, &mut buffer);
    }
}

// Return the n-th smallest element. O(log xs.len()) usually.
// If you force it to run in O(log xs.len()) almost always, use `rand` crate
// and randomize the pivot selection like `let pivot = xs.choose(&mut rng).unwrap();`.
use std::cmp::{PartialEq, PartialOrd};
fn select_nth_by<T: Clone, F: Fn(&T) -> K, K>(xs: &[T], n: usize, f: F) -> Option<K>
where
    K: PartialOrd + PartialEq,
{
    // Sainity check
    if xs.len() <= n {
        return None;
    }
    let pivot = f(&xs[xs.len() / 2]);
    let small = xs.iter().filter(|x| f(x) < pivot).count();
    let same = xs.iter().filter(|x| f(x) == pivot).count();
    // Recursive call.
    if n < small {
        // We can remove elements more than `pivot` from `xs`.
        let xs: Vec<_> = xs.iter().filter(|x| f(&x) < pivot).cloned().collect();
        select_nth_by(&xs, n, f)
    } else if small + same <= n {
        let xs: Vec<_> = xs.iter().filter(|x| f(&x) > pivot).cloned().collect();
        select_nth_by(&xs, n - small - same, f)
    } else {
        assert!(small <= n && n < small + same);
        return Some(pivot);
    }
}
