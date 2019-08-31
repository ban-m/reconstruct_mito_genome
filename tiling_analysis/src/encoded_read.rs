use std::collections::HashSet;
#[derive(Debug, Eq, PartialEq, Clone)]
pub enum Unit {
    Gap(usize),
    // (Unit, Subunit)
    Encoded(usize, usize),
}

impl Unit {
    fn new(unit: &str) -> Option<Unit> {
        let mut unit = unit.split('-');
        let unit_or_gap = unit.next()?;
        let is_gap = unit_or_gap.starts_with('G');
        let size = unit.next().and_then(|e| e.parse().ok())?;
        if is_gap {
            Some(Unit::Gap(size))
        } else {
            unit_or_gap.parse().ok().map(|res| Unit::Encoded(res, size))
        }
    }
    pub fn format_digit(&self) -> String {
        match self {
            &Unit::Gap(size) => format!("G-{:05}", size),
            &Unit::Encoded(unit, subunit) => format!("{:03}-{:03}", unit, subunit),
        }
    }
    pub fn get_unit(&self) -> Option<usize> {
        match self {
            &Unit::Gap(_) => None,
            &Unit::Encoded(u, _) => Some(u),
        }
    }
}

impl std::fmt::Display for Unit {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        match self {
            &Unit::Gap(size) => write!(f, "G-{}", size),
            &Unit::Encoded(unit, subunit) => write!(f, "{}-{}", unit, subunit),
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
    pub fn id(&self)->&str{
        &self.id
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
            read.push(Unit::new(unit)?);
        }
        Some(EncodedRead { id, read })
    }
    pub fn len(&self) -> usize {
        self.read.len()
    }
    pub fn iter<'a>(&'a self) -> Iter<'a> {
        Iter::new(self)
    }
    pub fn contains(&self, unit: usize) -> bool {
        self.read
            .iter()
            .filter_map(|e| e.get_unit())
            .any(|u| u == unit)
    }
    pub fn color(&self)->Option<usize>{
        use std::collections::HashMap;
        let mut counts = HashMap::new();
        for unit in self.read.iter().filter_map(|e|e.get_unit()){
            let count = counts.entry(unit).or_insert(0);
            *count += 1;
        }
        counts.into_iter().max_by_key(|e|e.1).map(|e|e.0)
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

pub fn collect_units<'a>(reads: &'a Vec<EncodedRead>) -> Vec<(usize, Vec<&'a EncodedRead>)> {
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
