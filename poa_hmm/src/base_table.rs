pub const BASE_TABLE: [usize; 128] = [
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 1, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 1, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
];

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn test() {
        assert!(BASE_TABLE[b'a' as usize] == 0, "{}", b'a');
        assert!(BASE_TABLE[b'A' as usize] == 0, "{}", b'A');
        assert!(BASE_TABLE[b'c' as usize] == 1);
        assert!(BASE_TABLE[b'C' as usize] == 1);
        assert!(BASE_TABLE[b'g' as usize] == 2);
        assert!(BASE_TABLE[b'G' as usize] == 2);
        assert!(BASE_TABLE[b't' as usize] == 3);
        assert!(BASE_TABLE[b'T' as usize] == 3);
    }
}
