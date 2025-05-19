#![allow(dead_code)]

use std::collections::BTreeMap;

/// Helper function to visualize sequence segments on each side of a breakend pos for debug output
///
/// # Arguments
/// * `pos`` - 0-indexed position immediatly after the breakend
/// * `window` - Number of bases flanking the breakend on each side to represent in the output
///
pub fn get_seq_pos_flanks(seq: &[u8], pos: i64, window: i64) -> (String, String) {
    (
        {
            let start = std::cmp::max(pos - window, 0) as usize;
            let end = std::cmp::max(pos, 0) as usize;
            std::str::from_utf8(&seq[start..end]).unwrap().to_string()
        },
        {
            let start = std::cmp::max(pos, 0) as usize;
            let end = std::cmp::max(pos + window, 0) as usize;
            std::str::from_utf8(&seq[start..end]).unwrap().to_string()
        },
    )
}

/// Track the total content within a fixed window_size, assuming observations at
/// monotonically increasing positions
///
pub struct SparseWindowSum {
    window_size: usize,
    sum: i32,
    pos_map: BTreeMap<i64, i32>,
}

impl SparseWindowSum {
    pub fn new(window_size: usize) -> Self {
        assert!(window_size > 1);
        Self {
            window_size,
            sum: 0,
            pos_map: BTreeMap::new(),
        }
    }

    pub fn sum(&self) -> i32 {
        self.sum
    }

    pub fn min_pos(&self) -> Option<i64> {
        let (pos, _) = self.pos_map.first_key_value()?;
        Some(*pos)
    }

    pub fn max_pos(&self) -> Option<i64> {
        let (pos, _) = self.pos_map.last_key_value()?;
        Some(*pos)
    }

    pub fn span(&self) -> Option<i64> {
        let min_pos = self.min_pos()?;
        let max_pos = self.max_pos()?;
        Some(1 + max_pos - min_pos)
    }

    pub fn is_empty_at_pos(&self, pos: i64) -> bool {
        if let Some((last_pos, _last_count)) = self.pos_map.last_key_value() {
            assert!(pos >= *last_pos);
            if (pos - *last_pos) >= self.window_size as i64 {
                return true;
            }
        }
        false
    }

    pub fn clear(&mut self) {
        self.sum = 0;
        self.pos_map.clear();
    }

    /// pos must be equal or higher than any previously pushed pos
    pub fn push(&mut self, pos: i64, count: i32) {
        if self.is_empty_at_pos(pos) {
            self.clear();
        } else {
            // Iterate through current members and trim every key below min_pos:
            let min_pos = 1 + pos - self.window_size as i64;
            let mut remove_pos = Vec::new();
            for (map_pos, map_count) in self.pos_map.iter() {
                if *map_pos >= min_pos {
                    break;
                }
                self.sum -= *map_count;
                remove_pos.push(*map_pos);
            }

            for map_pos in remove_pos {
                self.pos_map.remove(&map_pos);
            }
        }

        *self.pos_map.entry(pos).or_insert(0) += count;
        self.sum += count;
    }
}

/// Remove sorted indices from vector
///
/// Based on slightly diff design here
/// https://users.rust-lang.org/t/removing-multiple-indices-from-a-vector/65599/3
pub fn remove_sorted_indices<T>(v: &mut Vec<T>, indices: impl IntoIterator<Item = usize>) {
    let mut indices = indices.into_iter();
    let mut i = match indices.next() {
        None => return,
        Some(i) => i,
    };

    let mut tmp = Vec::new();
    std::mem::swap(v, &mut tmp);

    for (j, x) in tmp.into_iter().enumerate() {
        if j == i {
            if let Some(idx) = indices.next() {
                i = idx;
            }
        } else {
            v.push(x);
        }
    }
}

enum MatchState {
    Match,
    Mismatch,
    Gap,
}

/// Create line to highlight mismatches and gaps between alignment sequences
fn get_pairwise_alignment_mismatch_line(lines: &[&[u8]], start: usize, end: usize) -> String {
    let mut mismatches = String::new();
    for pos in start..end {
        use MatchState::*;
        let mut match_state = Match;
        let mut cons = None;
        for line in lines.iter().take(2) {
            if let Some(cons) = cons {
                if cons == b'-' || line[pos] == b'-' {
                    match_state = Gap;
                    break;
                } else if line[pos] != cons {
                    match_state = Mismatch;
                    break;
                }
            } else {
                cons = Some(line[pos]);
            }
        }
        mismatches.push(match match_state {
            Match => '|',
            Mismatch => 'X',
            Gap => ' ',
        });
    }
    mismatches
}

/// Handles wrapping a pairwise alignment over multiple lines and adding visualization help
///
/// Input is two lines for the top and bottom alignments, these must be the same length
///
/// If more than 2 lines of input are provided, the extra lines will be ignored
///
pub fn pairwise_alignment_printer(width: usize, lines: &[&[u8]]) {
    let ruler_tik_interval = 10;

    assert!(lines.len() > 1);
    let len = lines.first().unwrap().len();
    for line in lines.iter().take(2) {
        assert_eq!(line.len(), len);
    }

    let ruler = {
        let mut ruler = String::new();
        for i in 0..width {
            if (i + 1) % ruler_tik_interval == 0 {
                ruler.push('.');
            } else {
                ruler.push(' ');
            }
        }
        ruler
    };

    let rows = len.div_ceil(width);
    for row_index in 0..rows {
        let start = width * row_index;
        let end = std::cmp::min(start + width, len);
        let print_line = |index: usize| {
            eprintln!(
                "{}",
                std::str::from_utf8(&lines[index][start..end]).unwrap()
            );
        };

        eprintln!("{}", ruler);
        print_line(0);
        let mismatches = get_pairwise_alignment_mismatch_line(lines, start, end);
        eprintln!("{mismatches}");
        print_line(1);
        eprintln!();
    }
}

/// Handles wrapping long sequences into multiple rows
///
pub fn print_fasta(width: usize, lines: &[&[u8]]) {
    if lines.is_empty() {
        return;
    }

    for (seq_index, seq) in lines.iter().enumerate() {
        let len = seq.len();
        let rows = len.div_ceil(width);

        eprintln!("> {}", seq_index);
        for row_index in 0..rows {
            let start = width * row_index;
            let end = std::cmp::min(start + width, len);
            eprintln!("{}", std::str::from_utf8(&seq[start..end]).unwrap());
        }
    }
}

/// Mutate a to select the maximum value of a and b, where Some(x) is greater than None,
/// and Some(x) > Some(y) if x > y
pub fn option_f64_max(a: &mut Option<f64>, b: &Option<f64>) {
    *a = match b {
        Some(bval) => match *a {
            Some(aval) => Some(aval.max(*bval)),
            None => *b,
        },
        None => *a,
    };
}

/// Drop all true entries from vector
pub fn drop_true<T>(vec: &mut Vec<T>, drop_list: &[bool]) {
    assert_eq!(vec.len(), drop_list.len());
    let mut drop = drop_list.iter();
    vec.retain(|_| !*drop.next().unwrap())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_get_seq_pos_flanks() {
        let seq1 = b"ABCDEFGH";

        let (before, after) = get_seq_pos_flanks(seq1, 4, 2);
        assert_eq!(before, "CD");
        assert_eq!(after, "EF");
    }

    #[test]
    fn test_sparse_window_sum() {
        let mut sws = SparseWindowSum::new(3);

        assert_eq!(0, sws.sum());

        sws.push(100, 2);
        assert_eq!(2, sws.sum());
        sws.push(101, 2);
        assert_eq!(4, sws.sum());
        sws.push(102, 2);
        assert_eq!(6, sws.sum());
        sws.push(103, 2);
        assert_eq!(6, sws.sum());
        sws.push(103, 2);
        assert_eq!(8, sws.sum());
        assert_eq!(Some(101), sws.min_pos());

        sws.push(200, 2);
        assert_eq!(2, sws.sum());
    }

    #[test]
    fn test_remove_sorted_indices() {
        let mut v = vec![10, 20, 30, 40, 50];

        remove_sorted_indices(&mut v, [2, 4]);
        assert_eq!(v, vec![10, 20, 40]);
    }

    #[test]
    fn test_option_f64_max() {
        let mut a = None;
        let b = None;
        option_f64_max(&mut a, &b);
        assert_eq!(a, None);

        let mut a = None;
        let b = Some(1.0);
        option_f64_max(&mut a, &b);
        assert_eq!(a, Some(1.0));

        let mut a = Some(1.0);
        let b = None;
        option_f64_max(&mut a, &b);
        assert_eq!(a, Some(1.0));

        let mut a = Some(1.0);
        let b = Some(2.0);
        option_f64_max(&mut a, &b);
        assert_eq!(a, Some(2.0));
    }

    #[test]
    fn test_drop_true() {
        let mut v = vec![1, 2, 3, 4, 5];
        let p = vec![true, false, false, true, false];

        drop_true(&mut v, &p);
        assert_eq!(v, vec![2, 3, 5]);
    }
}
