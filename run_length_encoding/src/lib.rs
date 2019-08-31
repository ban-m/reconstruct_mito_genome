//! RunLength Encoder.

/// Runlength sequence. It is capable to switch the implementation or type for
/// inner struct, i.e., the way to represent each block of bases.
/// This struct allows you to access each index by usize. Thus,
/// You can do like `seq[i]` to retrieve the i-th runs of base.
/// To take the base, `seq[i].base()`, and to see the length `seq[i].length()` is
/// sufficient.
/// Note that currently there's no auxularity structure, thus, to indexing the
/// i-th sequence in the "decoded" sequence would take O(l) time,
/// where l is the length of the original sequence.
#[derive(Debug)]
pub struct RunLengthSeq<T: BaseWithLength + std::fmt::Debug> {
    id: usize,
    seq: Vec<T>,
}

/// The default implementation of run-length sequence.
/// # Example
///```
/// use runlength_sequence::DefaultSeq;
/// let seq = b"AATTGGCCATGC";
/// let id = 10;
/// let encoded_seq = DefaultSeq::new(id,seq);
///```
pub type DefaultSeq = RunLengthSeq<U8Base>;

impl<T> RunLengthSeq<T>
where
    T: BaseWithLength + std::fmt::Debug,
{
    /// Create new instance of run length encoded reads.
    /// It takes id of the read and the raw sequcence, which should be
    /// encoded by usual ASCII encoding. In other words,
    /// It should be a string over b"ATGCatgc". No other characters is allowed.
    /// # Example
    /// ```
    /// use runlength_sequence::{RunLengthSeq, U8Base, BaseWithLength};
    /// let seq = b"AATTGGCCATGC";
    /// let id = 10;
    /// let encoded_seq = RunLengthSeq::<U8Base>::new(id,seq);
    /// ```
    pub fn new(id: usize, raw_seq: &[u8]) -> Self {
        let mut seq = vec![];
        if raw_seq.is_empty() {
            Self { id, seq }
        } else {
            let mut base = raw_seq[0];
            let mut length = 1;
            for &x in &raw_seq[1..] {
                if base == x {
                    length += 1;
                } else {
                    seq.push(T::new(base, length));
                    base = x;
                    length = 1;
                }
            }
            seq.push(T::new(base, length));
            Self { id, seq }
        }
    }
    /// Get the id of the record.
    pub fn id(&self) -> usize {
        self.id
    }
    /// Get the length of the record. It is the same as "decoded" one's.
    pub fn len(&self) -> usize {
        self.seq.iter().map(|e| e.length() as usize).sum()
    }
    /// Get the length of the record. It is the same as "encoded" one's.
    /// If you want to iterate over the vector, please use `.seq()` instead.
    pub fn inner_len(&self) -> usize {
        self.seq.len()
    }
    /// Decoding to usual ASCII encoding(Maybe slow).
    pub fn decode(&self) -> Vec<u8> {
        self.seq
            .iter()
            .map(|base| (base.base(), base.length()))
            .fold(vec![], |mut seq, (base, length)| {
                (0..length).for_each(|_| seq.push(base));
                seq
            })
    }
    /// Get the n-th character in the original sequence.
    /// It would take O(l) time, where l is the length of the sequence.
    pub fn get_nth(&self, n: usize) -> u8 {
        let mut idx = 0;
        for x in &self.seq {
            if idx + x.length() as usize > n {
                return x.base();
            } else {
                idx += x.length() as usize;
            }
        }
        panic!(
            "Index out of range. The length is {} but the index is {}",
            self.inner_len(),
            n
        );
    }
    /// Get the reference to inner seqeunce from this struct.
    pub fn seq(&self) -> &Vec<T> {
        &self.seq
    }
}

impl<T> std::fmt::Display for RunLengthSeq<T>
where
    T: BaseWithLength + std::fmt::Debug,
{
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        let digit = 100;
        let len = self.inner_len();
        for i in 0..len / digit {
            let seq: String = self.seq[i * digit..(i + 1) * digit]
                .iter()
                .map(|e| e.base() as char)
                .collect();
            let num: String = self.seq[i * digit..(i + 1) * digit]
                .iter()
                .map(|e| (e.length().min(9) + b'0') as char)
                .collect();
            writeln!(f, "{}", seq)?;
            writeln!(f, "{}", num)?;
            writeln!(f)?;
        }
        let seq: String = self.seq[(len / digit) * digit..]
            .iter()
            .map(|e| e.base() as char)
            .collect();
        let num: String = self.seq[(len / digit) * digit..]
            .iter()
            .map(|e| (e.length().min(9) + b'0') as char)
            .collect();
        writeln!(f, "{}", seq)?;
        write!(f, "{}", num)
    }
}

impl<T> std::ops::Index<usize> for RunLengthSeq<T>
where
    T: BaseWithLength + std::fmt::Debug,
{
    type Output = T;
    fn index(&self, idx: usize) -> &Self::Output {
        &self.seq[idx]
    }
}

impl<T> std::ops::IndexMut<usize> for RunLengthSeq<T>
where
    T: BaseWithLength + std::fmt::Debug,
{
    fn index_mut<'a>(&'a mut self, idx: usize) -> &'a mut Self::Output {
        std::ops::IndexMut::index_mut(&mut self.seq, idx)
    }
}

/// This is the trait to represent each run of the base.
/// Each run consists of the combination of char(u8) and length(u8).
/// Thus, currently, no base with the length above 255 is represented correctly.
/// Also, currently alphabet exept b"ATCGacgt" is not supported.
pub trait BaseWithLength {
    /// Create a struct to represent base `base` with the length of `length`.
    fn new(base: u8, length: u8) -> Self;
    /// Retrieve the base insede.
    fn base(&self) -> u8;
    /// Retrieve the length inside, the upper bound is 255.
    /// It is not the maximum value. If you want to get the maxmum length can be represented,
    /// call .max() instead.
    fn length(&self) -> u8;
    /// Get the maxmum length of this struct. For example,
    /// `U8Base` uses the least 2 bits for representing base, while the rest for
    /// length counting. Thus, on calling `.max()` for U8Base, you get 0b111111.
    fn max(&self) -> u8;
    /// Get the complement of the base. It should be one of the b"ATCG".
    fn cmpl(&self) -> u8;
    /// Get the complement of the base. The return value is the same type as the original value.
    fn cmpl_with_len(&self) -> Self;
    /// Determine the equality regardless of sequence number.
    fn eq(&self,other:&Self)->bool;
}

pub type U8Base = u8;

impl BaseWithLength for U8Base {
    fn eq(&self, other:&Self)->bool{
        (self & 0b11) == (other & 0b11)
    }
    fn new(base: u8, length: u8) -> Self {
        let b = match base {
            b'A' | b'a' => 0b00,
            b'C' | b'c' => 0b01,
            b'G' | b'g' => 0b10,
            b'T' | b't' => 0b11,
            _ => 0,
        };
        b | length << 2
    }
    fn base(&self) -> u8 {
        match self & 0b11 {
            0b00 => b'A',
            0b01 => b'C',
            0b10 => b'G',
            0b11 => b'T',
            _ => b'-',
        }
    }
    fn length(&self) -> u8 {
        self >> 2
    }
    fn cmpl(&self) -> u8 {
       match self & 0b11 {
            0b00 => b'T',
            0b01 => b'G',
            0b10 => b'C',
            0b11 => b'A',
            _ => b'-',
        }
    }
    fn max(&self) -> u8 {
        0b111111
    }
    fn cmpl_with_len(&self) -> u8 {
        (self & 0b11111100) |  ((!self << 6) >> 6)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn it_works() {
        assert_eq!(2 + 2, 4);
    }
    #[test]
    fn compl() {
        for x in b"ACGT" {
            let base = U8Base::new(*x, 3);
            let cmpl = match base & 0b11 {
                0b00 => b'T',
                0b01 => b'G',
                0b10 => b'C',
                0b11 => b'A',
                _ => b'-',
            };
            assert_eq!(cmpl, base.cmpl_with_len().base())
        }
    }
}
