use crate::PartialOrderAlignment;
impl PartialOrderAlignment {
    pub fn topological_sort(mut self) -> Self {
        let edges = self.edges();
        let mapping = Self::topological_sort_inner(&edges);
        let mut result = vec![None; self.nodes.len()];
        let mut idx = self.nodes.len();
        while let Some(mut node) = self.nodes.pop() {
            idx -= 1;
            node.rename_by(&mapping);
            result[mapping[idx]] = Some(node);
        }
        assert!(result.iter().all(|e| e.is_some()));
        assert!(self.nodes.is_empty());
        self.nodes.extend(result.into_iter().filter_map(|e| e));
        self
    }
    // Return topological sorted order. The graph never contains cycles.
    fn topological_sort_inner(edges: &[Vec<usize>]) -> Vec<usize> {
        // 0 -> never arrived
        // 1 -> active(arrived, being traversed currently)
        // 2 -> inactive(arrived, have been traversed)
        let len = edges.len();
        let mut dfs_flag = vec![0; len];
        let mut dfs_stack = vec![];
        let mut mapping = vec![];
        for i in 0..len {
            if dfs_flag[i] != 0 {
                continue;
            }
            assert!(dfs_stack.is_empty());
            dfs_stack.push(i);
            'dfs: while !dfs_stack.is_empty() {
                let node = *dfs_stack.last().unwrap();
                if dfs_flag[node] == 0 {
                    // preorder
                    dfs_flag[node] = 1;
                }
                for &to in &edges[node] {
                    if dfs_flag[to] == 0 {
                        dfs_stack.push(to);
                        continue 'dfs;
                    }
                }
                // No-op
                let last = dfs_stack.pop().unwrap();
                mapping.push(last);
                // Deactivate
                dfs_flag[last] = 2;
            }
        }
        mapping.reverse();
        assert!(dfs_stack.is_empty());
        // Use DFS stack as a temporary buffer.
        std::mem::swap(&mut mapping, &mut dfs_stack);
        // Filling zero so that the boundary would be violated.
        mapping.extend(std::iter::repeat(0).take(len));
        for (order, &v) in dfs_stack.iter().enumerate() {
            mapping[v] = order;
        }
        mapping
    }
}
