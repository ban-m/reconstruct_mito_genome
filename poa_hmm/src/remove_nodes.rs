#![allow(dead_code)]
use super::THR;
impl crate::PartialOrderAlignment {
    pub fn remove_weight(mut self, f: f64) -> Self {
        let to_remove: Vec<_> = self.nodes.iter().map(|n| n.weight() <= f).collect();
        self = self.remove(&to_remove);
        self.trim_unreachable_nodes();
        self
    }
    pub fn remove_node(mut self) -> Self {
        let saved = self.clone();
        self = self.nodewise_remove();
        self = self.edgewise_remove();
        // self = self.merge_nodes();
        self.trim_unreachable_nodes();
        if self.nodes.len() < saved.nodes.len() / 10 {
            eprintln!("!!!!");
            saved
        } else {
            self
        }
    }
    fn nodewise_remove(mut self) -> Self {
        // Removing with node
        let (_start, arrived, _) = self.traverse();
        let to_remove: Vec<_> = arrived
            .into_iter()
            .zip(self.nodes.iter())
            //.map(|(b, n)| (n.is_tail && n.weight() <= 2.) || !b)
            .map(|(b, _)| !b)
            .collect();
        self.remove(&to_remove)
    }
    fn edgewise_remove(mut self) -> Self {
        // Removing with edges.
        let (_, _used_nodes, used_edges) = self.traverse();
        // let remaining_edges: Vec<_> = self
        //     .nodes
        //     .iter()
        //     .zip(used_edges.iter())
        //     .flat_map(|(n, e)| n.weights_except(e))
        //     .copied()
        //     .collect();
        // let arrived_len = used_nodes.iter().filter(|&&b| b).count();
        // let thr_rank = remaining_edges.len().max(arrived_len / FRAC) - arrived_len / FRAC;
        // let edge_thr = select_nth_by(&remaining_edges, thr_rank, |&x| x).unwrap_or(1.);
        self.nodes
            .iter_mut()
            .zip(used_edges)
            .for_each(|(n, e)| n.remove_edges(100., &e));
        self
    }
    fn traverse(&mut self) -> (usize, Vec<bool>, Vec<Vec<bool>>) {
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
            .max_by(|&(_, a), &(_, b)| a.head_weight().partial_cmp(&b.head_weight()).unwrap())
            .expect(&format!("{}:{}", line!(), self));
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
            ave * THR
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
                    weight_thr <= w || w == max
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
