use crate::PartialOrderAlignment;
impl PartialOrderAlignment {
    pub fn topological_sort(mut self) -> Self {
        let edges = self.edges();
        self.topological_sort_inner(&edges);
        let mapping = &self.mapping;
        self.nodes.iter_mut().for_each(|n| n.rename_by(&mapping));
        let len = self.nodes.len();
        self.dfs_flag.clear();
        self.dfs_stack.clear();
        let current_position = &mut self.dfs_flag;
        let original_position = &mut self.dfs_stack;
        current_position.extend(0..len);
        original_position.extend(0..len);
        for (from, &order) in mapping.iter().enumerate() {
            let current_pos = current_position[from];
            self.nodes.swap(current_pos, order);
            let swapped = original_position[order];
            original_position[order] = from;
            original_position[current_pos] = swapped;
            current_position[from] = order;
            current_position[swapped] = current_pos;
        }
        //assert_eq!(current_position, mapping);
        self
    }
    // Return topological sorted order. The graph never contains cycles.
    fn topological_sort_inner(&mut self, edges: &[Vec<usize>]) {
        // 0 -> never arrived
        // 1 -> active(arrived, being traversed currently)
        // 2 -> inactive(arrived, have been traversed)
        let len = edges.len();
        self.dfs_flag.clear();
        self.dfs_stack.clear();
        self.mapping.clear();
        self.dfs_flag.extend(std::iter::repeat(0).take(len));
        let dfs_flag = &mut self.dfs_flag;
        let dfs_stack = &mut self.dfs_stack;
        let mapping = &mut self.mapping;
        for i in 0..len {
            if dfs_flag[i] != 0 {
                continue;
            }
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
                let last = dfs_stack.pop().unwrap();
                mapping.push(last);
                // Deactivate
                dfs_flag[last] = 2;
            }
        }
        mapping.reverse();
        // Use DFS stack as a temporary buffer.
        std::mem::swap(mapping, dfs_stack);
        // Filling zero so that the boundary would not be violated.
        mapping.extend(std::iter::repeat(0).take(len));
        for (order, &v) in dfs_stack.iter().enumerate() {
            // The v-th node should be in the order-position.
            mapping[v] = order;
        }
    }
}
