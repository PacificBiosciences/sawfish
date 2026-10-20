use std::ffi::CString;
use std::io::{BufRead, Result};

pub struct FastqRecord {
    pub name: CString,
    pub seq: CString,
    pub qual: CString,
}

pub struct FastqReader<R: BufRead> {
    reader: R,
}

impl<R: BufRead> FastqReader<R> {
    pub fn new(reader: R) -> Self {
        FastqReader { reader }
    }
}

impl<R: BufRead> Iterator for FastqReader<R> {
    type Item = Result<FastqRecord>;

    fn next(&mut self) -> Option<Self::Item> {
        let mut f = || {
            {
                let mut buf = [0u8; 1];
                if self.reader.read(&mut buf)? == 0 {
                    return Ok(None);
                }
                assert!(buf[0] == b'@');
            }
            let name = {
                let mut buf = vec![];
                let _ = self.reader.read_until(b'\n', &mut buf)?;
                let _ = buf.pop(); // pop the '\n'
                if buf.last().copied() == Some(b'\r') {
                    *buf.last_mut().unwrap() = 0;
                } else {
                    buf.push(0);
                }
                CString::from_vec_with_nul(buf).unwrap()
            };
            let seq = {
                let mut buf = vec![];
                let _ = self.reader.read_until(b'\n', &mut buf)?;
                let _ = buf.pop(); // pop the '\n'
                if buf.last().copied() == Some(b'\r') {
                    *buf.last_mut().unwrap() = 0;
                } else {
                    buf.push(0);
                }
                CString::from_vec_with_nul(buf).unwrap()
            };
            {
                let mut buf = vec![];
                let _ = self.reader.read_until(b'\n', &mut buf)?;
                let _ = buf.pop(); // pop the '\n'
                if buf.last().copied() == Some(b'\r') {
                    let _ = buf.pop();
                }
                assert!(buf.last().copied().unwrap() == b'+');
            }
            let qual = {
                let mut buf = vec![];
                let _ = self.reader.read_until(b'\n', &mut buf)?;
                let _ = buf.pop(); // pop the '\n'
                if buf.last().copied() == Some(b'\r') {
                    *buf.last_mut().unwrap() = 0;
                } else {
                    buf.push(0);
                }
                CString::from_vec_with_nul(buf).unwrap()
            };
            Ok(Some(FastqRecord { name, seq, qual }))
        };
        f().transpose()
    }
}
