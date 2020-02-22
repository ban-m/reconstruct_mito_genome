impl crate::PartialOrderAlignment {
    pub fn remove_node(mut self) -> Self {
        let mut arrived = vec![false; self.nodes.len()];
        let (idx, mut node) = self
            .nodes
            .iter()
            .enumerate()
            .filter(|(_, e)| e.is_head)
            .max_by(|&(_, a), &(_, b)| a.head_weight().partial_cmp(&b.head_weight()).unwrap())
            .unwrap();
        arrived[idx] = true;
        let mut arrived_len = 1;
        while node.has_edge() {
            let (&idx, _) = node
                .edges()
                .iter()
                .zip(node.weights().iter())
                .max_by(|a, b| (a.1).partial_cmp(&(b.1)).unwrap())
                .unwrap();
            node = &self.nodes[idx];
            arrived[idx] = true;
            arrived_len += 1;
        }
        let remainings: Vec<_> = self
            .nodes
            .iter()
            .zip(arrived.iter())
            .filter_map(|(n, &b)| if b { None } else { Some(n.weight()) })
            .collect();
        let thr_rank = remainings.len() - arrived_len / 10;
        let thr = select_nth_by(&remainings, thr_rank, |&x| x).unwrap_or(1.);
        let to_remove: Vec<_> = self
            .nodes
            .iter()
            .zip(arrived)
            .map(|(n, a)| !(a | (n.weight() > thr)))
            .collect();
        self = self.remove(&to_remove);
        self.trim_unreachable_nodes();
        self
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
