pub fn comp_base(x: &u8) -> u8 {
    match *x {
        b'A' => b'T',
        b'T' => b'A',
        b'C' => b'G',
        b'G' => b'C',
        b'N' => b'N',
        b'a' => b't',
        b't' => b'a',
        b'c' => b'g',
        b'g' => b'c',
        b'n' => b'n',
        _ => b'N',
    }
}

/// Return reverse complement sequence
///
/// Complements [ACGTacgt] and preserves case, any other base character is complemented as 'N'
///
pub fn rev_comp(dna: &[u8]) -> Vec<u8> {
    dna.iter().rev().map(comp_base).collect::<Vec<_>>()
}

/// Reverse complement sequence in-place
///
/// Complements [ACGTacgt] and preserves case, any other base character is complemented as 'N'
///
pub fn rev_comp_in_place(dna: &mut [u8]) {
    let len = dna.len();
    let halflen = len - len / 2;
    for i in 0..halflen {
        dna[i] = comp_base(&dna[i]);
        let rev_i = len - 1 - i;
        if i != rev_i {
            dna[rev_i] = comp_base(&dna[rev_i]);
            dna.swap(i, rev_i);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_rev_comp() {
        let input = b"NNATGCG".to_vec();
        let expected_output = b"CGCATNN".to_vec();
        let output = rev_comp(&input);
        assert_eq!(output, expected_output);
    }

    #[test]
    fn test_rev_comp_in_place() {
        let mut input = b"NNATGCG".to_vec();
        let expected_output = b"CGCATNN".to_vec();
        rev_comp_in_place(&mut input);
        assert_eq!(input, expected_output);
    }
}
