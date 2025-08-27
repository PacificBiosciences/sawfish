use std::collections::HashMap;
use std::fs::File;

use bio::io::fasta;
use log::info;
use unwrap::unwrap;

#[derive(Default)]
pub struct GenomeRef {
    /// A map from chrom name to chrom sequence
    pub chroms: HashMap<String, Vec<u8>>,
}

impl GenomeRef {
    /// Scan all GenomeRef seqeunces and convert any characters not in the allowed list to the unknown character
    ///
    pub fn convert_disallowed_characters(&mut self, allowed: &[u8], unknown: u8) {
        let allowed_lut = {
            // Build a lookup table for all possible u8 values
            const TYPE_WIDTH: usize = (u8::MAX as usize) + 1;
            let mut x = [false; TYPE_WIDTH];
            for &c in allowed.iter() {
                x[c as usize] = true;
            }
            x
        };
        for seq in self.chroms.values_mut() {
            for c in seq.iter_mut().filter(|x| !allowed_lut[**x as usize]) {
                *c = unknown;
            }
        }
    }

    /// Convert all bases besides "ACGTN" to "N"
    pub fn simplify_ambiguous_dna_bases(&mut self) {
        self.convert_disallowed_characters(b"ACGTN", b'N');
    }
}

/// Read fasta file pointer into GenomeRef data structure
///
/// This method converts all input characters to upper-case
///
pub fn get_genome_ref_from_fasta_fp(file: File) -> GenomeRef {
    let reader = fasta::Reader::new(file);

    let mut genome_ref = GenomeRef::default();

    for result in reader.records() {
        let record = result.expect("Error during fasta record parsing");

        genome_ref
            .chroms
            .insert(record.id().to_string(), record.seq().to_ascii_uppercase());

        /*
        let chrom = record.id().to_string();
        let seq = genome_ref.chroms.get(chrom.as_str()).unwrap();
        let mini = seq.iter().take(10).copied().collect::<Vec<_>>();
        println!("Read chrom {} length {}", chrom, seq.len());
        println!("Read mini {}", std::str::from_utf8(&mini).unwrap());
         */
    }
    genome_ref
}

/// Read fasta file into GenomeRef data structure
///
/// This method converts all input characters to upper-case
///
pub fn get_genome_ref_from_fasta(filename: &str) -> GenomeRef {
    info!("Reading reference genome from file '{filename}'");

    let file = unwrap!(
        File::open(filename),
        "Unable to open reference fasta file: '{}'",
        filename,
    );

    get_genome_ref_from_fasta_fp(file)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::{Seek, SeekFrom, Write};

    #[test]
    fn test_get_genome_ref_from_fasta_fp() {
        let mut file = tempfile::tempfile().unwrap();

        let cname = "foo";
        let seq = "ACGTACGT";
        writeln!(file, ">{cname}").unwrap();
        writeln!(file, "{seq}").unwrap();
        file.seek(SeekFrom::Start(0)).unwrap();
        let result = get_genome_ref_from_fasta_fp(file);

        assert_eq!(result.chroms.len(), 1);
        assert_eq!(result.chroms.keys().next().unwrap(), cname);
        assert_eq!(
            std::str::from_utf8(result.chroms.values().next().unwrap()).unwrap(),
            seq
        );
    }

    #[test]
    fn test_simplify_ambiguous_dna_bases() {
        let mut chroms = HashMap::default();
        chroms.insert(String::from("foo"), b"ACGT1234acgtNNNNMMMM".to_vec());
        let mut genome_ref = GenomeRef { chroms };

        genome_ref.simplify_ambiguous_dna_bases();

        assert_eq!(genome_ref.chroms["foo"], b"ACGTNNNNNNNNNNNNNNNN".to_vec());
    }
}
