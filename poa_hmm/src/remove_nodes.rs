#![allow(dead_code)]
impl crate::PartialOrderAlignment {
    pub fn remove_node(mut self, thr: f64) -> Self {
        if self.nodes.len() < 2 {
            panic!("Invalid input:{}", self);
        }
        let saved = self.clone();
        let len = self.nodes.len();
        self = self
            .nodewise_remove(thr)
            .unwrap_or_else(|| panic!("N:{}", saved));
        self = self
            .edgewise_remove(thr)
            .unwrap_or_else(|| panic!("E:{}", saved));
        self.trim_unreachable_nodes();
        if self.nodes.len() < len / 10 {
            saved
        } else {
            self
        }
    }
    fn nodewise_remove(mut self, thr: f64) -> Option<Self> {
        let (_start, arrived, _) = self.traverse(thr)?;
        let to_remove: Vec<_> = arrived
            .into_iter()
            .zip(self.nodes.iter())
            .map(|(b, _)| !b)
            .collect();
        Some(self.remove(&to_remove))
    }
    fn edgewise_remove(mut self, thr: f64) -> Option<Self> {
        let (_, _, used_edges) = self.traverse(thr)?;
        self.nodes
            .iter_mut()
            .zip(used_edges)
            .for_each(|(n, e)| n.remove_edges(&e));
        Some(self)
    }
    fn traverse(&mut self, thr: f64) -> Option<(usize, Vec<bool>, Vec<Vec<bool>>)> {
        let mut arrived = vec![false; self.nodes.len()];
        let mut used_edges: Vec<_> = self
            .nodes
            .iter()
            .map(|n| vec![false; n.edges.len()])
            .collect();
        let (start, _) = self
            .nodes
            .iter()
            .enumerate()
            .filter(|(_, e)| e.is_head)
            .max_by(|&(_, a), &(_, b)| a.head_weight().partial_cmp(&b.head_weight()).unwrap())?;
        let mut queue = std::collections::VecDeque::new();
        queue.push_back(start);
        arrived[start] = true;
        let weight_thr = {
            let sum = self
                .nodes
                .iter()
                .flat_map(|n| n.weights.iter())
                .sum::<f64>();
            let ave = sum / self.nodes.iter().map(|n| n.edges.len()).sum::<usize>() as f64;
            ave * thr
        };
        while !queue.is_empty() {
            let idx = queue.pop_front().unwrap();
            let node = &self.nodes[idx];
            let max = *node
                .weights()
                .iter()
                .max_by(|a, b| a.partial_cmp(b).unwrap())
                .unwrap_or(&0.);
            for (edg, (&next, &w)) in node.edges().iter().zip(node.weights.iter()).enumerate() {
                let is_heavy = if node.is_tail {
                    weight_thr <= w
                } else {
                    weight_thr <= w || ((w - max) / w.max(max)).abs() < 0.000_1
                };
                if is_heavy {
                    if !arrived[next] {
                        queue.push_back(next);
                    }
                    arrived[next] = true;
                    used_edges[idx][edg] = true;
                }
            }
        }
        Some((start, arrived, used_edges))
    }
    fn remove(mut self, to_remove: &[bool]) -> Self {
        let (num, mapping) = to_remove
            .iter()
            .fold((0, vec![]), |(index, mut map), remove| {
                map.push((index, !remove));
                (index + (!remove) as usize, map)
            });
        self.nodes = self
            .nodes
            .into_iter()
            .zip(to_remove)
            .filter(|(_, &b)| !b)
            .map(|(mut node, _)| {
                node.remove_if(&mapping);
                node
            })
            .collect();
        assert_eq!(self.nodes.len(), num);
        self
    }
    fn trim_unreachable_nodes(&mut self) {
        self.trim_unreachable_nodes_forward();
        self.trim_unreachable_nodes_reverse();
    }
    fn trim_unreachable_nodes_forward(&mut self) {
        let mut arrived = vec![false; self.nodes.len()];
        let mut bfs_queue: std::collections::VecDeque<_> = self
            .nodes
            .iter()
            .enumerate()
            .filter_map(|(idx, n)| if n.is_head { Some(idx) } else { None })
            .collect();
        while !bfs_queue.is_empty() {
            let node = bfs_queue.pop_front().unwrap();
            if !arrived[node] {
                arrived[node] = true;
                for &to in self.nodes[node].edges.iter() {
                    bfs_queue.push_back(to);
                }
            }
        }
        let (mapping, _) = arrived.iter().fold((vec![], 0), |(mut map, index), &b| {
            map.push((index, b));
            (map, index + b as usize)
        });
        let mut buffer = vec![];
        let mut idx = self.nodes.len();
        while let Some(mut node) = self.nodes.pop() {
            idx -= 1;
            if arrived[idx] {
                node.remove_if(&mapping);
                buffer.push(node);
            }
        }
        buffer.reverse();
        assert!(self.nodes.is_empty());
        std::mem::swap(&mut self.nodes, &mut buffer);
    }
    fn trim_unreachable_nodes_reverse(&mut self) {
        let mut arrived = vec![false; self.nodes.len()];
        let mut bfs_queue: std::collections::VecDeque<_> = self
            .nodes
            .iter()
            .enumerate()
            .filter_map(|(idx, n)| if n.is_tail { Some(idx) } else { None })
            .collect();
        let edges = self.reverse_edges();
        if bfs_queue.is_empty() {
            return;
        }
        while !bfs_queue.is_empty() {
            let node = bfs_queue.pop_front().unwrap();
            if !arrived[node] {
                arrived[node] = true;
                for &to in &edges[node] {
                    bfs_queue.push_back(to);
                }
            }
        }
        let (mapping, _) = arrived.iter().fold((vec![], 0), |(mut map, index), &b| {
            map.push((index, b));
            (map, index + b as usize)
        });
        let mut buffer = vec![];
        let mut idx = self.nodes.len();
        while let Some(mut node) = self.nodes.pop() {
            idx -= 1;
            if arrived[idx] {
                node.remove_if(&mapping);
                buffer.push(node);
            }
        }
        assert!(self.nodes.is_empty());
        buffer.reverse();
        std::mem::swap(&mut self.nodes, &mut buffer);
    }
    pub fn remove_weight(mut self, f: f64) -> Self {
        let to_remove: Vec<_> = self.nodes.iter().map(|n| n.weight() <= f).collect();
        self = self.remove(&to_remove);
        self.trim_unreachable_nodes();
        self
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
        Some(pivot)
    }
}
