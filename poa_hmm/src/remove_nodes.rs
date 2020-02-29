#![allow(dead_code)]
use super::{FRAC, THR};
impl crate::PartialOrderAlignment {
    pub fn remove_node(mut self) -> Self {
        let saved = self.clone();
        self = self.nodewise_remove();
        self = self.edgewise_remove();
        // self = self.merge_nodes();
        self.trim_unreachable_nodes();
        if self.nodes.len() < 10 {
            saved
        } else {
            assert!(self.nodes.len() > 10);
            self
        }
    }
    fn nodewise_remove(mut self) -> Self {
        // Removing with node
        let (start, arrived, _) = self.traverse();
        let arrived_len = arrived.iter().filter(|&&b| b).count();
        let remaining_nodes: Vec<_> = self
            .nodes
            .iter()
            .zip(arrived.iter())
            .filter_map(|(n, &b)| if b { None } else { Some(n.weight()) })
            .collect();
        let thr_rank = remaining_nodes.len().max(arrived_len / FRAC) - arrived_len / FRAC;
        let thr = select_nth_by(&remaining_nodes, thr_rank, |&x| x).unwrap_or(1.);
        let median = select_nth_by(&self.nodes, self.nodes.len() / 2, |n| n.weight()).unwrap();
        let deviation = |n: &super::Base| (n.weight() - median).abs();
        let mad = select_nth_by(&self.nodes, self.nodes.len() / 2, deviation).unwrap();
        let thr = thr.max(median - mad * 10.);
        let to_remove: Vec<_> = self
            .nodes
            .iter()
            .zip(arrived)
            .enumerate()
            .map(|(idx, (n, a))| {
                if n.is_head || n.is_tail {
                    n.weight() <= thr && idx != start
                } else {
                    n.weight() <= thr && !a
                }
            })
            .collect();
        self.remove(&to_remove)
    }
    fn edgewise_remove(mut self) -> Self {
        // Removing with edges.
        let (_, used_nodes, used_edges) = self.traverse();
        let remaining_edges: Vec<_> = self
            .nodes
            .iter()
            .zip(used_edges.iter())
            .flat_map(|(n, e)| n.weights_except(e))
            .copied()
            .collect();
        let arrived_len = used_nodes.iter().filter(|&&b| b).count();
        let thr_rank = remaining_edges.len() - arrived_len / 10;
        let edge_thr = select_nth_by(&remaining_edges, thr_rank, |&x| x).unwrap_or(1.);
        self.nodes
            .iter_mut()
            .zip(used_edges)
            .for_each(|(n, e)| n.remove_edges(edge_thr, &e, 0.1));
        self
    }
    fn traverse(&mut self) -> (usize, Vec<bool>, Vec<Vec<bool>>) {
        let mut arrived = vec![false; self.nodes.len()];
        //let mut used_edges = vec![0; self.nodes.len()];
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
            .max_by(|&(_, a), &(_, b)| a.head_weight().partial_cmp(&b.head_weight()).unwrap())
            .expect(&format!("{}:{}", line!(), self));
        let mut queue = std::collections::VecDeque::new();
        queue.push_back(start);
        arrived[start] = true;
        while !queue.is_empty() {
            let idx = queue.pop_front().unwrap();
            let node = &self.nodes[idx];
            if !node.has_edge() {
                continue;
            }
            let max = node
                .weights
                .iter()
                .max_by(|a, b| a.partial_cmp(b).unwrap())
                .unwrap();
            for (edg, (&next, &w)) in node.edges().iter().zip(node.weights.iter()).enumerate() {
                if max * THR <= w && !arrived[next] {
                    arrived[next] = true;
                    //used_edges[idx] = next;
                    used_edges[idx][edg] = true;
                    queue.push_back(next);
                }
            }
        }
        (start, arrived, used_edges)
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
    // fn select_head(&mut self) {
    //     let mut is_head = vec![true; self.nodes.len()];
    //     for n in self.nodes.iter() {
    //         for &to in n.edges.iter() {
    //             is_head[to] = false;
    //         }
    //     }
    //     if self.nodes.iter().all(|e| !e.is_head) {
    //         self.nodes.iter_mut().zip(is_head).for_each(|(n, is_head)| {
    //             n.is_head = is_head;
    //         });
    //     }
    //     assert!(self.nodes.iter().any(|n| n.is_head));
    // }

    // fn select_head_tail(mut self) -> Self {
    //     let mut is_head = vec![true; self.nodes.len()];
    //     for n in self.nodes.iter() {
    //         for &to in n.edges.iter() {
    //             is_head[to] = false;
    //         }
    //     }
    //     assert!(is_head.iter().any(|&e| e), "{:?}", self);
    //     if self.nodes.iter().all(|e| !e.is_head) {
    //         self.nodes.iter_mut().zip(is_head).for_each(|(n, is_head)| {
    //             n.is_head = is_head;
    //         });
    //     }
    //     if self.nodes.iter().all(|e| !e.is_tail) {
    //         self.nodes.iter_mut().for_each(|n| {
    //             n.is_tail = n.edges.is_empty();
    //         });
    //     }
    //     assert!(self.nodes.iter().any(|n| n.is_head));
    //     assert!(self.nodes.iter().any(|n| n.is_tail));
    //     self
    // }
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
        Some(pivot)
    }
}
