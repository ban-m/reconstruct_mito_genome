use std::collections::HashSet;
#[derive(Debug, Eq, PartialEq, Clone, Hash)]
pub enum Unit {
    Gap(usize),
    // (Unit, Subunit)
    Encoded(u8, u8, u16),
}

impl std::cmp::Ord for Unit {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        use std::cmp::Ordering::*;
        match (self, other) {
            (_, Unit::Gap(_)) => std::cmp::Ordering::Greater,
            (Unit::Gap(_), _) => std::cmp::Ordering::Less,
            (Unit::Encoded(c1, u1, s1), Unit::Encoded(c2, u2, s2)) => match c1.cmp(c2) {
                Equal => match u1.cmp(u2) {
                    Equal => s1.cmp(s2),
                    x => x,
                },
                x => x,
            },
        }
    }
}

impl std::cmp::PartialOrd for Unit {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}

impl Unit {
    fn new(unit: &str) -> Result<Unit, std::num::ParseIntError> {
        let mut contents = unit.split('-');
        let contig_or_gap = contents.next().unwrap();
        let is_gap = contig_or_gap.starts_with('G');
        if is_gap {
            let size = contents.next().unwrap().parse()?;
            Ok(Unit::Gap(size))
        } else {
            let contig = contig_or_gap.parse()?;
            let unit = contents.next().unwrap().parse()?;
            let subunit = contents.next().unwrap().parse()?;
            Ok(Unit::Encoded(contig, unit, subunit))
        }
    }
    pub fn get_unit(&self) -> Option<u8> {
        match self {
            &Unit::Gap(_) => None,
            &Unit::Encoded(_, u, _) => Some(u),
        }
    }
    pub fn to_output(&self) -> String {
        match self {
            &Unit::Gap(size) => format!("G-{}", size),
            &Unit::Encoded(contig, unit, subunit) => format!("{}-{}-{}", contig, unit, subunit),
        }
    }
}

impl std::fmt::Display for Unit {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        match self {
            &Unit::Gap(size) => write!(f, "G-{:08}", size),
            &Unit::Encoded(contig, unit, subunit) => {
                write!(f, "{:02}-{:03}-{:03}", contig, unit, subunit)
            }
        }
    }
}

#[derive(Eq, PartialEq, Debug)]
pub struct EncodedRead {
    pub id: String,
    pub read: Vec<Unit>,
}

impl std::fmt::Display for EncodedRead {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        let units: Vec<_> = self.read.iter().map(|e| format!("{}", e)).collect();
        write!(f, "{}:{}", self.id, units.join(" "))
    }
}

impl EncodedRead {
    pub fn is_empty(&self) -> bool {
        self.read.is_empty()
    }
    pub fn id(&self) -> &str {
        &self.id
    }
    pub fn to_output(&self) -> String {
        let units: Vec<_> = self.read.iter().map(|e| e.to_output()).collect();
        format!("{}:{}", self.id, units.join(" "))
    }
    fn sanity_check(line: &str) -> Option<()> {
        let m = line.matches(':').count();
        if m == 1 {
            Some(())
        } else if m == 0 {
            eprintln!("The line {} is invalid format.", line);
            eprintln!("No id-field delimiter ':' was found. Please check it.");
            None
        } else if m > 1 {
            eprintln!("The line {} is invalid format.", line);
            eprintln!("Multiple id-filed demilter ':' was found. Please check it.");
            None
        } else {
            unreachable!()
        }
    }
    pub fn new(line: &str) -> Option<EncodedRead> {
        Self::sanity_check(line)?;
        let mut line = line.split(':');
        let id = line.next()?.to_string();
        let mut read = vec![];
        for unit in line.next()?.split_whitespace() {
            read.push(Unit::new(unit).ok()?);
        }
        Some(EncodedRead { id, read })
    }
    pub fn len(&self) -> usize {
        self.read.len()
    }
    pub fn iter<'a>(&'a self) -> Iter<'a> {
        Iter::new(self)
    }
    pub fn contains(&self, unit: u8) -> bool {
        self.read
            .iter()
            .filter_map(|e| e.get_unit())
            .any(|u| u == unit)
    }
    pub fn color(&self) -> Option<u8> {
        use std::collections::HashMap;
        let mut counts = HashMap::new();
        for unit in self.read.iter().filter_map(|e| e.get_unit()) {
            let count = counts.entry(unit).or_insert(0);
            *count += 1;
        }
        counts.into_iter().max_by_key(|e| e.1).map(|e| e.0)
    }
}

pub struct Iter<'a> {
    inner: &'a EncodedRead,
    pointer: usize,
}

impl<'a> Iter<'a> {
    fn new(a: &'a EncodedRead) -> Self {
        Self {
            inner: a,
            pointer: 0,
        }
    }
}

impl<'a> std::iter::Iterator for Iter<'a> {
    type Item = &'a Unit;
    fn next(&mut self) -> Option<Self::Item> {
        if self.pointer >= self.inner.len() {
            None
        } else {
            self.pointer += 1;
            Some(&self.inner[self.pointer - 1])
        }
    }
}

impl std::ops::Index<usize> for EncodedRead {
    type Output = Unit;
    fn index(&self, index: usize) -> &Unit {
        &self.read[index]
    }
}

pub fn collect_units<'a>(reads: &'a Vec<EncodedRead>) -> Vec<(u8, Vec<&'a EncodedRead>)> {
    let units: HashSet<_> = reads
        .iter()
        .flat_map(|e| e.iter().filter_map(|e| e.get_unit()))
        .collect();
    units
        .iter()
        .map(|&unit| {
            (
                unit,
                reads.iter().filter(|read| read.contains(unit)).collect(),
            )
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use rand::distributions::Alphanumeric;
    use rand::{thread_rng, Rng};
    #[test]
    fn unit_test() {
        let mut rng = thread_rng();
        for gapsize in (0..1000).map(|_| rng.gen()) {
            assert_eq!(
                Unit::new(&format!("G-{}", gapsize)),
                Some(Unit::Gap(gapsize))
            );
        }
        for _ in 0..1000 {
            let unit = rng.gen();
            let subunit = rng.gen();
            let res = Unit::new(&format!("{}-{}", unit, subunit)).unwrap();
            assert_eq!(res, Unit::Encoded(unit, subunit));
        }
        assert_eq!(Unit::new("L-32"), None);
        assert_eq!(Unit::new("-12334"), None);
        assert_eq!(Unit::new(" sdf2303= - 43"), None);
    }
    fn generate(rng: &mut rand::rngs::ThreadRng) -> Unit {
        if rng.gen_bool(0.001) {
            Unit::Gap(rng.gen())
        } else {
            Unit::Encoded(rng.gen(), rng.gen())
        }
    }

    #[test]
    fn encoded_read_check() {
        let mut rng = thread_rng();
        for _ in 0..100 {
            let id: String = std::iter::repeat(())
                .map(|_| rng.sample(Alphanumeric))
                .take(8)
                .collect();
            let read: Vec<_> = std::iter::repeat(())
                .map(|_| generate(&mut rng))
                .take(1000)
                .collect();
            let format: Vec<_> = read.iter().map(|e| format!("{}", e)).collect();
            let format = id.clone() + ":" + &format.join(" ");
            let encoded_read = EncodedRead::new(&format).unwrap();
            for i in 0..encoded_read.len() {
                assert_eq!(&encoded_read[i], &read[i]);
            }
            assert_eq!(encoded_read, EncodedRead { id, read });
        }
    }
    #[test]
    fn print_check() {
        let mut rng = thread_rng();
        for _ in 0..100 {
            let read = gen_read(&mut rng);
            let as_string = format!("{}", read);
            let read_recover = EncodedRead::new(&as_string);
            assert_eq!(Some(read), read_recover);
        }
    }
}
