extern crate rand;
pub mod encoded_read;
use encoded_read::*;
#[derive(Debug, PartialEq, Eq)]
pub enum Ed {
    Match(usize),
    UpGap(usize),
    DownGap(usize),
}

#[derive(Debug)]
pub struct Alignment<'a> {
    pub read1: &'a EncodedRead,
    pub read2: &'a EncodedRead,
    pub score: i64,
    pub path: Vec<Ed>,
    pub start_1: usize,
    pub start_2: usize,
}

const CHUNKSIZE: usize = 30;
impl<'a> std::fmt::Display for Alignment<'a> {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        writeln!(
            f,
            "Read1:{}({} length)\nRead2:{}({} length)\nScore:{}",
            &self.read1.id,
            self.read1.len(),
            &self.read2.id,
            self.read2.len(),
            self.score
        )?;
        let (mut upper, mut middle, mut lower) = (vec![], vec![], vec![]);
        let (mut pos1, mut pos2) = (self.start_1, self.start_2);
        for ed in self.path.iter() {
            match ed {
                &Ed::Match(l) => {
                    for i in 0..l {
                        upper.push(self.read1[pos1 + i].format_digit());
                        lower.push(self.read2[pos2 + i].format_digit());
                        if self.read1[pos1 + i] == self.read2[pos2 + i] {
                            middle.push("| unt |".to_string());
                        } else {
                            middle.push("X unt X".to_string());
                        }
                    }
                    pos1 += l;
                    pos2 += l;
                }
                &Ed::UpGap(l) => {
                    for i in 0..l {
                        upper.push("| gap |".to_string());
                        middle.push("       ".to_string());
                        lower.push(self.read2[pos2 + i].format_digit());
                    }
                    pos2 += l;
                }
                &Ed::DownGap(l) => {
                    for i in 0..l {
                        upper.push(self.read1[pos1 + i].format_digit());
                        middle.push("       ".to_string());
                        lower.push("| gap |".to_string());
                    }
                    pos1 += l;
                }
            }
        }
        for ((uc, mc), lc) in upper
            .chunks(CHUNKSIZE)
            .zip(middle.chunks(CHUNKSIZE))
            .zip(lower.chunks(CHUNKSIZE))
        {
            writeln!(f, "{}\n{}\n{}", uc.join(" "), mc.join(" "), lc.join(" "))?;
            writeln!(f, "")?;
        }
        Ok(())
    }
}

impl<'a> Alignment<'a> {
    pub fn new<F>(a: &'a EncodedRead, b: &'a EncodedRead, gap: i64, matchfun: F) -> Self
    where
        F: Fn(&Unit, &Unit) -> i64,
    {
        let mut dp = vec![vec![0; b.len() + 1]; a.len() + 1];
        let mut direction: Vec<Vec<u8>> = vec![vec![0; b.len() + 1]; a.len() + 1];
        let mut max = 0;
        let mut argmax = (0, 0);
        for i in 0..a.len() {
            for j in 0..b.len() {
                let (i, j) = (i + 1, j + 1);
                let gap_in_1 = dp[i][j - 1] + gap;
                let gap_in_2 = dp[i - 1][j] + gap;
                let match_12 = dp[i - 1][j - 1] + matchfun(&a[i - 1], &b[j - 1]);
                if gap_in_1 < 0 && gap_in_2 < 0 && match_12 < 0 {
                    dp[i][j] = 0;
                    direction[i][j] = 4; // Meaning "stop"
                } else if gap_in_1 > gap_in_2 && gap_in_1 > match_12 {
                    dp[i][j] = gap_in_1;
                    direction[i][j] = 1; // Meaning Gap on 1
                } else if gap_in_2 > gap_in_1 && gap_in_2 > match_12 {
                    dp[i][j] = gap_in_2;
                    direction[i][j] = 2; // Meaning Gap on 2
                } else {
                    dp[i][j] = match_12;
                    direction[i][j] = 3; // Meaning Match
                }
                if dp[i][j] > max {
                    max = dp[i][j];
                    argmax = (i, j);
                }
            }
        }
        let score = max;
        let (path, start_1, start_2) = Self::recover_an_optimal_path(&direction, argmax);
        let (read1, read2) = (a, b);
        Alignment {
            read1,
            read2,
            score,
            path,
            start_1,
            start_2,
        }
    }
    fn recover_an_optimal_path(
        direction: &Vec<Vec<u8>>,
        (starti, startj): (usize, usize),
    ) -> (Vec<Ed>, usize, usize) {
        let (mut current_i, mut current_j) = (starti, startj);
        let mut current_op = direction[current_i][current_j];
        let mut operation_num = 0;
        let mut path = vec![];
        while current_i != 0 && current_j != 0 && direction[current_i][current_j] != 4 {
            if direction[current_i][current_j] == current_op {
                operation_num += 1;
                if current_op == 3 {
                    current_i = current_i - 1;
                    current_j = current_j - 1;
                } else if current_op == 1 {
                    current_j = current_j - 1;
                } else if current_op == 2 {
                    current_i = current_i - 1;
                } else {
                    unreachable!()
                }
            } else {
                match current_op {
                    1 => path.push(Ed::UpGap(operation_num)),
                    2 => path.push(Ed::DownGap(operation_num)),
                    3 => path.push(Ed::Match(operation_num)),
                    _ => unreachable!(),
                };
                operation_num = 0;
                current_op = direction[current_i][current_j];
            }
        }
        match current_op {
            0 => {}
            1 => path.push(Ed::UpGap(operation_num)),
            2 => path.push(Ed::DownGap(operation_num)),
            3 => path.push(Ed::Match(operation_num)),
            _ => unreachable!(),
        };
        path.reverse();
        (path, current_i, current_j)
    }
}

#[cfg(test)]
mod tests {
    fn gen_read(rng: &mut rand::rngs::ThreadRng) -> EncodedRead {
        let id: String = std::iter::repeat(())
            .map(|_| rng.sample(Alphanumeric))
            .take(8)
            .collect();
        let len = rng.gen_range(10, 100);
        let read: Vec<_> = std::iter::repeat(())
            .map(|_| generate(rng))
            .take(len)
            .collect();
        EncodedRead { id, read }
    }
    #[test]
    fn alignment() {
        let mut rng = thread_rng();
        let read: Vec<_> = std::iter::repeat(())
            .map(|_| gen_read(&mut rng))
            .take(10)
            .collect();
        for i in 0..read.len() {
            for j in 0..read.len() {
                Alignment::new(&read[i], &read[j], -5, |a, b| if a == b { 1 } else { -1 });
            }
        }
    }
    #[test]
    fn alignment_2() {
        let read1 = "m54113_160913_184949/11076505/0_25281:G-1245 6-12 6-13 6-14 6-15 6-16 G-2952 6-3 6-4 6-5 6-6 6-7 6-8 6-9 6-10 6-11 G-3735 6-0 6-1 6-2 G-4186 5-0 G-4700 4-19 4-20 G-5622 4-16 4-17 4-18 G-5622 4-2 4-3 4-4 4-5 4-6 4-7 4-8 4-9 4-10 4-11 4-12 4-13 4-14 G-8154 2-91 G-8702 4-0 4-1 G-9262 3-109 3-110 3-111 G-14658 3-83 3-84 3-85 3-86 3-87 3-88 3-89 3-90 3-91 3-92 3-93 3-94 3-95 3-96 3-97 3-98 3-99 3-100 3-101 3-102 3-103 3-104 3-105 3-106 3-107 3-108";
        let read2 = "m54113_160913_184949/21299573/0_16964:G-3547 6-23 6-24 6-25 6-26 6-27 6-28 6-29 6-30 6-31 6-32 6-33 6-34 6-35 6-36 6-37 6-38 G-6153 6-10 6-11 6-12 6-13 6-14 6-15 6-16 6-17 6-18 6-19 6-20 6-21 6-22 G-7093 6-6 6-7 6-8 6-9 G-8326 6-1 6-2 6-3 6-4 G-8722 5-0 G-9821 4-1 4-2 4-3 4-4 4-5 4-6 4-7 4-8 4-9 4-10 4-11 4-12 4-13 4-14 4-15 4-16 4-17 4-18 4-19 4-20 4-21 G-13318 2-91 2-92 G-13906 3-109 3-110 3-111 G-14892 3-105 3-106 3-107 3-108";
        let read1 = EncodedRead::new(read1).unwrap();
        let read2 = EncodedRead::new(read2).unwrap();
        Alignment::new(&read1, &read2, -5, |a, b| if a == b { 1 } else { -1 });
    }
    #[test]
    fn alignment_3() {
        let read1 = EncodedRead {
            id: "1".to_string(),
            read: vec![Unit::Encoded(10, 10); 10],
        };
        let read2 = EncodedRead {
            id: "2".to_string(),
            read: vec![Unit::Encoded(10, 10); 10],
        };
        let alignment = Alignment::new(&read1, &read2, -5, |a, b| if a == b { 1 } else { -1 });
        assert_eq!(alignment.path, vec![Ed::Match(10)]);
        let read1 = EncodedRead {
            id: "1".to_string(),
            read: vec![
                Unit::Encoded(10, 10),
                Unit::Encoded(10, 10),
                Unit::Encoded(10, 10),
                Unit::Encoded(10, 10),
                Unit::Encoded(10, 10),
                Unit::Encoded(9, 121),
                Unit::Encoded(10, 10),
                Unit::Encoded(10, 10),
                Unit::Encoded(10, 10),
                Unit::Encoded(10, 10),
                Unit::Encoded(10, 10),
                Unit::Encoded(10, 10),
            ],
        };
        let read2 = EncodedRead {
            id: "2".to_string(),
            read: vec![
                Unit::Encoded(10, 10),
                Unit::Encoded(10, 10),
                Unit::Encoded(10, 10),
                Unit::Encoded(10, 10),
                Unit::Encoded(10, 10),
                Unit::Encoded(8, 121),
                Unit::Encoded(10, 10),
                Unit::Encoded(10, 10),
                Unit::Encoded(10, 10),
                Unit::Encoded(10, 10),
                Unit::Encoded(10, 10),
                Unit::Encoded(10, 10),
            ],
        };
        let alignment = Alignment::new(&read1, &read2, -5, |a, b| if a == b { 1 } else { -1 });
        assert_eq!(alignment.path, vec![Ed::Match(12)]);
        let read1 = EncodedRead {
            id: "1".to_string(),
            read: vec![
                Unit::Encoded(10, 10),
                Unit::Encoded(10, 11),
                Unit::Encoded(10, 12),
                Unit::Encoded(10, 13),
                Unit::Encoded(10, 14),
                Unit::Encoded(10, 15),
                Unit::Encoded(10, 16),
                Unit::Encoded(10, 17),
                Unit::Encoded(10, 18),
                Unit::Encoded(10, 19),
                Unit::Encoded(10, 20),
            ],
        };
        let read2 = EncodedRead {
            id: "2".to_string(),
            read: vec![
                Unit::Encoded(10, 10),
                Unit::Encoded(10, 11),
                Unit::Encoded(10, 12),
                Unit::Encoded(10, 13),
                Unit::Encoded(10, 14),
                Unit::Encoded(8, 121),
                Unit::Encoded(10, 15),
                Unit::Encoded(10, 16),
                Unit::Encoded(10, 17),
                Unit::Encoded(10, 18),
                Unit::Encoded(10, 19),
                Unit::Encoded(10, 20),
            ],
        };
        let alignment = Alignment::new(&read1, &read2, -4, |a, b| if a == b { 1 } else { -1 });
        assert_eq!(
            alignment.path,
            vec![Ed::Match(5), Ed::UpGap(1), Ed::Match(6)]
        );
    }
}
