impl crate::PartialOrderAlignment {
    pub fn remove_node(mut self) -> Self {
        let saved = self.clone();
        self = self.nodewise_remove();
        self = self.edgewise_remove();
        // self = self.remove_triangles();
        self = self.merge_nodes();
        self.trim_unreachable_nodes();
        if self.nodes.len() < 10 {
            return saved;
        } else {
            //self.adjust_weight(super::DEFAULT_WEIGHT);
            assert!(self.nodes.len() > 10);
            self
        }
    }
    fn adjust_weight(&mut self, weight: f64) {
        // The unwrap is safe.
        let max = self
            .nodes
            .iter()
            .map(|n| n.weight())
            .max_by(|a, b| a.partial_cmp(b).unwrap())
            .unwrap();
        let factor = weight / max;
        self.nodes
            .iter_mut()
            .for_each(|n| n.weight = n.weight.mul_add(factor, 1.));
        self.weight = self.weight.mul_add(factor, 1.);
    }
    fn merge_nodes(mut self) -> Self {
        let margeable_nodes = self.get_margeable_nodes();
        for (i, (_, b)) in margeable_nodes.into_iter().enumerate().filter(|x| (x.1).0) {
            let merge_idx: Vec<_> = self
                .nodes
                .iter()
                .enumerate()
                .filter(|&(_, n)| n.base() == b && n.edges.contains(&i))
                .map(|(idx, _)| idx)
                .collect();
            let (u, v) = match merge_idx.as_slice() {
                &[u, v] if u < v => (u, v),
                &[_, _] => panic!("u >= v, violating topological order"),
                _ => panic!("merge index's length is not 2."),
            };
            if self.not_reachable(u, v) {
                // Merge u and v.
                self.merge(u, v);
            }
        }
        self
    }
    // Return true if v is reachable from u.
    fn not_reachable(&self, u: usize, v: usize) -> bool {
        let mut stack = vec![u];
        while let Some(node) = stack.pop() {
            for &to in self.nodes[node].edges().iter() {
                if to == v {
                    return true;
                } else if to > v {
                    continue;
                } else {
                    stack.push(to);
                }
            }
        }
        false
    }
    // merge the v-th node into the u-th node.
    fn merge(&mut self, u: usize, v: usize) {
        let weight_and_edges: Vec<_> = self.nodes[v]
            .weights()
            .iter()
            .zip(self.nodes[v].edges().iter())
            .map(|(&w, &e)| (w, e))
            .collect();
        for (w, to) in weight_and_edges {
            match self.nodes[u].edges.iter().position(|&e| e == to) {
                Some(x) => self.nodes[u].weights[x] += w,
                None => {
                    self.nodes[u].edges.push(to);
                    self.nodes[u].weights.push(w);
                }
            }
        }
        self.nodes[v].weights.clear();
        self.nodes[v].edges.clear();
        self.nodes.iter_mut().for_each(|n| n.remove(v));
    }
    fn get_margeable_nodes(&self) -> Vec<(bool, u8)> {
        let mut base_counts = vec![(0, 0, 0, 0); self.nodes.len()];
        for n in self.nodes.iter() {
            for &to in n.edges().iter() {
                match n.base() {
                    b'A' => base_counts[to].0 += 1,
                    b'C' => base_counts[to].1 += 1,
                    b'G' => base_counts[to].2 += 1,
                    b'T' => base_counts[to].3 += 1,
                    _ => {}
                }
            }
        }
        fn is_margeable((a, c, g, t): (u8, u8, u8, u8)) -> (bool, u8) {
            if a == 2 {
                (true, b'A')
            } else if c == 2 {
                (true, b'C')
            } else if g == 2 {
                (true, b'G')
            } else if t == 2 {
                (true, b'T')
            } else {
                (false, 0)
            }
        }
        base_counts.into_iter().map(is_margeable).collect()
    }
    // fn remove_triangles(mut self) -> Self {
    //     let (used_nodes, used_edges) = self.traverse();
    //     for i in 0..self.nodes.len() {
    //         let target = {
    //             let node = &self.nodes[i];
    //             if !used_nodes[i]
    //                 || node.edges().len() != 2
    //                 || !node.edges.iter().all(|&to| used_nodes[to])
    //             {
    //                 continue;
    //             }
    //             let base = self.nodes[node.edges[0]].base();
    //             let is_all = node
    //                 .edges()
    //                 .iter()
    //                 .map(|&j| self.nodes[j].base())
    //                 .all(|e| e == base);
    //             let passed = used_edges[used_edges[i]];
    //             if is_all && node.edges.contains(&passed) {
    //                 self.nodes[i].edges().iter().position(|&e| e == passed)
    //             } else {
    //                 None
    //             }
    //         };
    //         if let Some(target) = target {
    //             self.nodes[i].edges.remove(target);
    //             self.nodes[i].weights.remove(target);
    //         }
    //     }
    //     self
    // }
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
        let thr_rank = remaining_nodes.len().max(arrived_len / 10) - arrived_len / 10;
        let thr = select_nth_by(&remaining_nodes, thr_rank, |&x| x).unwrap_or(1.);
        let to_remove: Vec<_> = self
            .nodes
            .iter()
            .zip(arrived)
            .enumerate()
            .map(|(idx, (n, a))| {
                if n.is_head || n.is_tail {
                    n.weight() < thr && idx != start
                } else {
                    n.weight() < thr && !a
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
            .flat_map(|(n, &e)| n.weights_except(e))
            .copied()
            .collect();
        let arrived_len = used_nodes.iter().filter(|&&b| b).count();
        let thr_rank = remaining_edges.len() - arrived_len / 10;
        // let ave = remaining_edges.iter().sum::<f64>() / remaining_edges.len() as f64;
        let edge_thr = select_nth_by(&remaining_edges, thr_rank, |&x| x).unwrap_or(1.);
        self.nodes
            .iter_mut()
            .zip(used_edges)
            .for_each(|(n, e)| n.remove_edges(edge_thr, e, 0.1));
        self
    }
    fn traverse(&mut self) -> (usize, Vec<bool>, Vec<usize>) {
        // if self.nodes.iter().all(|n| !n.is_head) {
        //     self.select_head();
        // }
        let mut arrived = vec![false; self.nodes.len()];
        let mut used_edges = vec![0; self.nodes.len()];
        let (start, mut node) = self
            .nodes
            .iter()
            .enumerate()
            .filter(|(_, e)| e.is_head)
            .max_by(|&(_, a), &(_, b)| a.head_weight().partial_cmp(&b.head_weight()).unwrap())
            .expect(&format!("{}:{}", line!(), self));
        arrived[start] = true;
        let mut idx = start;
        while node.has_edge() {
            let (&next, _) = node
                .edges()
                .iter()
                .zip(node.weights().iter())
                .max_by(|a, b| (a.1).partial_cmp(&(b.1)).unwrap())
                .expect(&format!("{}", line!()));
            arrived[next] = true;
            used_edges[idx] = next;
            // Update
            idx = next;
            node = &self.nodes[next];
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
