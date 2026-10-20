use std::cell::RefCell;
use std::fs::File;
use std::io::{BufRead, BufReader, Read};
use std::os::unix::fs::FileExt;
use std::sync::atomic::{AtomicU64, Ordering};

use log::{error, info, warn};
use rust_vc_utils::{ChromList, GenomeSegment, get_genome_ref_from_fasta};
use unwrap::unwrap;

/// Location of one chromosome sequence in the fasta file, as described by a samtools faidx index
struct FastaIndexEntry {
    offset: u64,
    line_bases: u64,
    line_width: u64,
}

/// Minimum size of the reference window read for `GenomeRef::with_seq`
const REF_WINDOW_SIZE: usize = 1 << 20;

/// Positional read access to the uncompressed content of a fasta file, which is either plain or bgzip compressed
enum FastaSource {
    Plain(File),
    Bgzf(BgzfFile),
}

impl FastaSource {
    /// Open the fasta, detecting bgzip compression from the file content rather than the file name
    fn open(filename: &str) -> std::io::Result<Self> {
        let file = File::open(filename)?;
        let mut magic = [0u8; 2];
        if file.read_exact_at(&mut magic, 0).is_ok() && magic == GZIP_MAGIC {
            Ok(Self::Bgzf(BgzfFile::new(file)?))
        } else {
            Ok(Self::Plain(file))
        }
    }

    fn read_exact_at(&self, buf: &mut [u8], offset: u64) -> std::io::Result<()> {
        match self {
            Self::Plain(file) => file.read_exact_at(buf, offset),
            Self::Bgzf(file) => file.read_exact_at(buf, offset),
        }
    }
}

const GZIP_MAGIC: [u8; 2] = [0x1f, 0x8b];

/// Fixed part of a bgzf block header, followed by the extra field of XLEN bytes
const BGZF_HEADER_SIZE: usize = 12;

/// Size of the crc32 and uncompressed size which end each bgzf block
const BGZF_TRAILER_SIZE: usize = 8;

fn invalid_bgzf(filename_hint: &str) -> std::io::Error {
    std::io::Error::new(
        std::io::ErrorKind::InvalidData,
        format!("{filename_hint}: file is gzip compressed but not in bgzip (bgzf) block format"),
    )
}

/// A bgzip file with random access to its uncompressed content
///
/// Each bgzf block is a complete gzip member, so the file is read by locating the blocks containing
/// the requested uncompressed range and inflating only those. The block layout is found by walking
/// the block headers once on open, so no .gzi index is needed.
///
struct BgzfFile {
    file: File,
    /// (compressed file offset, uncompressed offset) of each non-empty block's start
    blocks: Vec<(u64, u64)>,
    /// Total uncompressed size
    len: u64,
}

impl BgzfFile {
    /// Total size of the block starting at `offset`, read from its header
    fn read_block_size(file: &File, offset: u64) -> std::io::Result<usize> {
        let mut header = [0u8; BGZF_HEADER_SIZE];
        file.read_exact_at(&mut header, offset)?;
        // gzip magic, deflate method, and the extra field flag
        if header[..2] != GZIP_MAGIC || header[2] != 8 || header[3] & 4 == 0 {
            return Err(invalid_bgzf("bgzf block header"));
        }
        let extra_len = u16::from_le_bytes([header[10], header[11]]) as usize;
        let mut extra = vec![0u8; extra_len];
        file.read_exact_at(&mut extra, offset + BGZF_HEADER_SIZE as u64)?;
        let mut i = 0;
        while i + 4 <= extra.len() {
            let subfield_len = u16::from_le_bytes([extra[i + 2], extra[i + 3]]) as usize;
            if extra[i] == b'B' && extra[i + 1] == b'C' && subfield_len == 2 && i + 6 <= extra.len()
            {
                let size = u16::from_le_bytes([extra[i + 4], extra[i + 5]]) as usize + 1;
                if size < BGZF_HEADER_SIZE + extra_len + BGZF_TRAILER_SIZE {
                    break;
                }
                return Ok(size);
            }
            i += 4 + subfield_len;
        }
        Err(invalid_bgzf("bgzf block header"))
    }

    fn new(file: File) -> std::io::Result<Self> {
        let file_len = file.metadata()?.len();
        let mut blocks = Vec::new();
        let mut offset = 0u64;
        let mut uncompressed = 0u64;
        while offset < file_len {
            let size = Self::read_block_size(&file, offset)? as u64;
            if offset + size > file_len {
                return Err(invalid_bgzf("truncated bgzf block"));
            }
            let mut isize_bytes = [0u8; 4];
            file.read_exact_at(&mut isize_bytes, offset + size - 4)?;
            let block_len = u32::from_le_bytes(isize_bytes) as u64;
            if block_len > 0 {
                blocks.push((offset, uncompressed));
            }
            uncompressed += block_len;
            offset += size;
        }
        Ok(Self {
            file,
            blocks,
            len: uncompressed,
        })
    }

    /// Inflate the block starting at the compressed file offset
    fn read_block(&self, offset: u64) -> std::io::Result<Vec<u8>> {
        let size = Self::read_block_size(&self.file, offset)?;
        let mut block = vec![0u8; size];
        self.file.read_exact_at(&mut block, offset)?;
        let extra_len = u16::from_le_bytes([block[10], block[11]]) as usize;
        let payload = &block[(BGZF_HEADER_SIZE + extra_len)..(size - BGZF_TRAILER_SIZE)];
        let expected_len = u32::from_le_bytes(block[(size - 4)..].try_into().unwrap()) as usize;
        let mut data = Vec::with_capacity(expected_len);
        flate2::read::DeflateDecoder::new(payload).read_to_end(&mut data)?;
        if data.len() != expected_len {
            return Err(invalid_bgzf("bgzf block with inconsistent size"));
        }
        Ok(data)
    }

    fn read_exact_at(&self, buf: &mut [u8], offset: u64) -> std::io::Result<()> {
        if buf.is_empty() {
            return Ok(());
        }
        if offset + buf.len() as u64 > self.len {
            return Err(std::io::ErrorKind::UnexpectedEof.into());
        }
        let mut block_index = self.blocks.partition_point(|x| x.1 <= offset) - 1;
        let mut filled = 0;
        let mut pos = offset;
        while filled < buf.len() {
            let (block_offset, block_start) = self.blocks[block_index];
            let data = self.read_block(block_offset)?;
            let skip = (pos - block_start) as usize;
            let n = std::cmp::min(data.len() - skip, buf.len() - filled);
            buf[filled..(filled + n)].copy_from_slice(&data[skip..(skip + n)]);
            filled += n;
            pos += n as u64;
            block_index += 1;
        }
        Ok(())
    }
}

enum ChromSeqs {
    /// Sequence is read from the fasta file as it is requested, using its faidx index
    Indexed {
        filename: String,
        /// The fasta file, opened once and read with positional reads from any thread
        file: FastaSource,
        index: Vec<FastaIndexEntry>,
        lens: Vec<usize>,
    },

    /// All sequence is held in memory
    Loaded(Vec<Vec<u8>>),
}

/// Reference sequence for each chromosome in ref_chrom_list order
///
/// When the reference fasta has a faidx index, sequence is read from the fasta only as it is
/// requested, and is not kept in memory, so that memory use doesn't grow with the fraction of the
/// genome touched by the run.
///
pub struct GenomeRef {
    /// Unique id of this reference, to identify it in each thread's reference window
    id: u64,
    seqs: ChromSeqs,
}

static NEXT_GENOME_REF_ID: AtomicU64 = AtomicU64::new(0);

/// The most recent reference window read on this thread
struct RefWindow {
    genome_ref_id: u64,
    chrom_index: usize,
    start: usize,
    seq: Vec<u8>,
}

thread_local! {
    static REF_WINDOW: RefCell<Option<RefWindow>> = const { RefCell::new(None) };
}

/// Convert sequence to upper case and replace any characters other than ACGTN with N
///
/// This matches the sequence transformations applied when reading the full fasta through rust_vc_utils.
///
fn simplify_ref_base(c: u8) -> u8 {
    match c.to_ascii_uppercase() {
        x @ (b'A' | b'C' | b'G' | b'T' | b'N') => x,
        _ => b'N',
    }
}

/// Read chromosome sequence over [start, end) from the fasta file into `seq`
fn read_chrom_range(
    filename: &str,
    file: &FastaSource,
    entry: &FastaIndexEntry,
    start: usize,
    end: usize,
    seq: &mut [u8],
) {
    assert_eq!(seq.len(), end - start);
    if start == end {
        return;
    }
    let file_offset = |pos: u64| {
        entry.offset + (pos / entry.line_bases) * entry.line_width + pos % entry.line_bases
    };
    let file_start = file_offset(start as u64);
    let file_end = file_offset(end as u64 - 1) + 1;

    let mut buffer = vec![0u8; (file_end - file_start) as usize];
    unwrap!(
        file.read_exact_at(&mut buffer, file_start),
        "Unable to read chromosome sequence from reference fasta file: '{}'",
        filename
    );

    let mut seq_len = 0;
    for c in buffer.into_iter().filter(|&c| c != b'\n' && c != b'\r') {
        assert!(
            seq_len < seq.len(),
            "Unexpected chromosome sequence length read from reference fasta file '{filename}'"
        );
        seq[seq_len] = simplify_ref_base(c);
        seq_len += 1;
    }
    assert_eq!(
        seq_len,
        seq.len(),
        "Unexpected chromosome sequence length read from reference fasta file '{filename}'"
    );
}

impl GenomeRef {
    fn new(seqs: ChromSeqs) -> Self {
        Self {
            id: NEXT_GENOME_REF_ID.fetch_add(1, Ordering::Relaxed),
            seqs,
        }
    }

    /// Create from sequences already loaded in memory, in chrom_list order
    pub fn from_chrom_seqs(seqs: Vec<Vec<u8>>) -> Self {
        Self::new(ChromSeqs::Loaded(seqs))
    }

    pub fn chrom_count(&self) -> usize {
        match &self.seqs {
            ChromSeqs::Indexed { lens, .. } => lens.len(),
            ChromSeqs::Loaded(seqs) => seqs.len(),
        }
    }

    pub fn chrom_len(&self, chrom_index: usize) -> usize {
        match &self.seqs {
            ChromSeqs::Indexed { lens, .. } => lens[chrom_index],
            ChromSeqs::Loaded(seqs) => seqs[chrom_index].len(),
        }
    }

    /// Get chromosome sequence over [start, end)
    pub fn get_seq(&self, chrom_index: usize, start: usize, end: usize) -> Vec<u8> {
        assert!(start <= end && end <= self.chrom_len(chrom_index));
        match &self.seqs {
            ChromSeqs::Indexed {
                filename,
                file,
                index,
                ..
            } => {
                let mut seq = vec![0u8; end - start];
                read_chrom_range(filename, file, &index[chrom_index], start, end, &mut seq);
                seq
            }
            ChromSeqs::Loaded(seqs) => seqs[chrom_index][start..end].to_vec(),
        }
    }

    /// Run f on chromosome sequence over [start, end)
    ///
    /// This is intended for many requests in nearby regions on each thread, such as the reference
    /// for each read in a region. The reference is read in windows of at least REF_WINDOW_SIZE,
    /// and each thread keeps its most recent window.
    ///
    /// f must not call with_seq.
    ///
    pub fn with_seq<T>(
        &self,
        chrom_index: usize,
        start: usize,
        end: usize,
        f: impl FnOnce(&[u8]) -> T,
    ) -> T {
        let chrom_len = self.chrom_len(chrom_index);
        assert!(start <= end && end <= chrom_len);
        let (filename, file, index) = match &self.seqs {
            ChromSeqs::Indexed {
                filename,
                file,
                index,
                ..
            } => (filename, file, index),
            ChromSeqs::Loaded(seqs) => return f(&seqs[chrom_index][start..end]),
        };
        REF_WINDOW.with(|window| {
            let mut window = window.borrow_mut();
            let is_in_window = window.as_ref().is_some_and(|x| {
                x.genome_ref_id == self.id
                    && x.chrom_index == chrom_index
                    && x.start <= start
                    && end <= x.start + x.seq.len()
            });
            if !is_in_window {
                let window_end = std::cmp::min(
                    std::cmp::max(end, start.saturating_add(REF_WINDOW_SIZE)),
                    chrom_len,
                );
                let mut seq = window.take().map(|x| x.seq).unwrap_or_default();
                seq.resize(window_end - start, 0);
                read_chrom_range(
                    filename,
                    file,
                    &index[chrom_index],
                    start,
                    window_end,
                    &mut seq,
                );
                *window = Some(RefWindow {
                    genome_ref_id: self.id,
                    chrom_index,
                    start,
                    seq,
                });
            }
            let window = window.as_ref().unwrap();
            f(&window.seq[(start - window.start)..(end - window.start)])
        })
    }

    /// Run f on each consecutive chunk of up to chunk_size bases of chrom_index, in order
    ///
    /// For an indexed fasta only one chunk is held in memory at a time.
    ///
    pub fn for_each_chrom_chunk(
        &self,
        chrom_index: usize,
        chunk_size: usize,
        mut f: impl FnMut(&[u8]),
    ) {
        assert!(chunk_size > 0);
        let chrom_len = self.chrom_len(chrom_index);
        match &self.seqs {
            ChromSeqs::Indexed {
                filename,
                file,
                index,
                ..
            } => {
                let mut seq = Vec::new();
                for start in (0..chrom_len).step_by(chunk_size) {
                    let end = std::cmp::min(start + chunk_size, chrom_len);
                    seq.resize(end - start, 0);
                    read_chrom_range(filename, file, &index[chrom_index], start, end, &mut seq);
                    f(&seq);
                }
            }
            ChromSeqs::Loaded(seqs) => {
                for chunk in seqs[chrom_index].chunks(chunk_size) {
                    f(chunk);
                }
            }
        }
    }
}

/// Parse faidx index, returning None if it can't be read or isn't valid
///
/// Logs a warning if the index file exists but can't be used.
///
fn read_fasta_index(fai_filename: &str) -> Option<Vec<(String, u64, FastaIndexEntry)>> {
    let file = File::open(fai_filename).ok()?;
    let entries = parse_fasta_index(BufReader::new(file));
    if entries.is_none() {
        warn!(
            "Ignoring fasta index '{fai_filename}' because it is malformed, the whole reference will be loaded instead"
        );
    }
    entries
}

fn parse_fasta_index(reader: impl BufRead) -> Option<Vec<(String, u64, FastaIndexEntry)>> {
    let mut entries = Vec::new();
    for line in reader.lines() {
        let line = line.ok()?;
        let words = line.split('\t').collect::<Vec<_>>();
        if words.len() < 5 {
            return None;
        }
        let length = words[1].parse::<u64>().ok()?;
        let entry = FastaIndexEntry {
            offset: words[2].parse().ok()?,
            line_bases: words[3].parse().ok()?,
            line_width: words[4].parse().ok()?,
        };
        if entry.line_bases == 0 || entry.line_width < entry.line_bases {
            return None;
        }
        entries.push((words[0].to_string(), length, entry));
    }
    Some(entries)
}

/// True if the file starts with the gzip magic bytes
fn is_gzip_file(filename: &str) -> bool {
    let mut magic = [0u8; 2];
    File::open(filename)
        .and_then(|f| f.read_exact_at(&mut magic, 0))
        .is_ok()
        && magic == GZIP_MAGIC
}

/// Error message if the fasta is compressed but has no usable index
///
/// A compressed fasta can only be read through its index, because the whole-fasta loader doesn't decompress.
///
fn get_missing_index_error(
    ref_filename: &str,
    fai_filename: &str,
    fai_entries: &Option<Vec<(String, u64, FastaIndexEntry)>>,
) -> Option<String> {
    (fai_entries.is_none() && is_gzip_file(ref_filename)).then(|| {
        format!(
            "Reference fasta '{ref_filename}' is compressed, which requires a 'samtools faidx' index at '{fai_filename}'."
        )
    })
}

/// True if the index was modified before the fasta it describes
fn is_index_older_than_fasta(fai_filename: &str, ref_filename: &str) -> bool {
    let modified = |x: &str| std::fs::metadata(x).and_then(|m| m.modified()).ok();
    match (modified(fai_filename), modified(ref_filename)) {
        (Some(fai_time), Some(ref_time)) => fai_time < ref_time,
        _ => false,
    }
}

fn report_missing_chrom(chrom_label: &str) {
    error!(
        "Chromosome \"{chrom_label}\" specified in the read alignment file, but not in the reference fasta"
    );
}

fn report_inconsistent_chrom_length(chrom_label: &str, chrom_len: u64, fasta_len: usize) {
    error!(
        "Chromosome \"{chrom_label}\" specified with inconsistent length: {chrom_len} in the read alignment file, and {fasta_len} in the reference fasta"
    );
}

fn exit_on_reference_consistency_error() -> ! {
    error!("Exiting due to one or more reference consistency issues");
    std::process::exit(exitcode::DATAERR);
}

/// Setup reference sequence access for the chromosomes in ref_chrom_list, using the fasta faidx index
///
fn get_indexed_genome_ref(
    ref_filename: &str,
    ref_chrom_list: &ChromList,
    fai_entries: Vec<(String, u64, FastaIndexEntry)>,
) -> GenomeRef {
    let mut fai_entries = fai_entries
        .into_iter()
        .map(|(label, length, entry)| (label, (length, entry)))
        .collect::<std::collections::HashMap<_, _>>();

    let mut index = Vec::new();
    let mut lens = Vec::new();
    let mut reference_consistency_error = false;
    for chrom_info in ref_chrom_list.data.iter() {
        let chrom_label = &chrom_info.label;
        let chrom_len = chrom_info.length;
        match fai_entries.remove(chrom_label) {
            Some((fasta_len, entry)) => {
                if fasta_len != chrom_len {
                    report_inconsistent_chrom_length(chrom_label, chrom_len, fasta_len as usize);
                    reference_consistency_error = true;
                } else {
                    index.push(entry);
                    lens.push(chrom_len as usize);
                }
            }
            None => {
                report_missing_chrom(chrom_label);
                reference_consistency_error = true;
            }
        }
    }

    if reference_consistency_error {
        exit_on_reference_consistency_error();
    }

    let file = unwrap!(
        FastaSource::open(ref_filename),
        "Unable to open reference fasta file: '{}'",
        ref_filename
    );
    GenomeRef::new(ChromSeqs::Indexed {
        filename: ref_filename.to_string(),
        file,
        index,
        lens,
    })
}

/// Load reference sequence from fasta, stored as an array in ref_chrom_list order
///
/// If the fasta has a faidx index, chromosome sequence is read from the fasta as it is requested.
///
pub fn get_ref_chrom_array(ref_filename: &str, ref_chrom_list: &ChromList) -> GenomeRef {
    let fai_filename = format!("{ref_filename}.fai");
    let fai_entries = read_fasta_index(&fai_filename);
    if let Some(message) = get_missing_index_error(ref_filename, &fai_filename, &fai_entries) {
        error!("{message}");
        std::process::exit(exitcode::DATAERR);
    }
    if let Some(fai_entries) = fai_entries {
        if is_index_older_than_fasta(&fai_filename, ref_filename) {
            warn!(
                "Fasta index '{fai_filename}' is older than '{ref_filename}', so it may not match the fasta. Regenerate it with 'samtools faidx'."
            );
        }
        info!("Reading reference genome on demand from indexed fasta file '{ref_filename}'");
        return get_indexed_genome_ref(ref_filename, ref_chrom_list, fai_entries);
    }

    let mut genome_ref = get_genome_ref_from_fasta(ref_filename);
    genome_ref.simplify_ambiguous_dna_bases();

    // For each chrom listed in the asssembly-to-ref bam index, check that it is found in genome_ref,
    // and transfer the chrom sequence from genome_ref to a chrom array
    //
    let mut r = Vec::new();
    let mut reference_consistency_error = false;
    for chrom_info in ref_chrom_list.data.iter() {
        let chrom_label = &chrom_info.label;
        let chrom_len = chrom_info.length;
        match genome_ref.chroms.remove(chrom_label) {
            Some(x) => {
                if x.len() != chrom_len as usize {
                    report_inconsistent_chrom_length(chrom_label, chrom_len, x.len());
                    reference_consistency_error = true;
                } else {
                    r.push(x);
                }
            }
            None => {
                report_missing_chrom(chrom_label);
                reference_consistency_error = true;
            }
        };
    }

    if reference_consistency_error {
        exit_on_reference_consistency_error();
    }

    GenomeRef::from_chrom_seqs(r)
}

/// Get gap-compressed sequence identity from the bam record, aligned to chrom_index
///
/// If the record's CIGAR uses only =/X (no M) the identity is computed from the CIGAR alone, so
/// neither the read sequence nor the reference is accessed. Otherwise only the reference sequence
/// covered by the alignment is read.
///
pub fn get_gap_compressed_identity(
    record: &rust_htslib::bam::Record,
    reference: &GenomeRef,
    chrom_index: usize,
) -> f64 {
    if let Some(gci) = crate::bam_utils::get_gap_compressed_identity_from_raw_cigar_no_align_match(
        record.raw_cigar(),
    ) {
        return gci;
    }

    let read_seq = crate::bam_utils::get_simplified_dna_seq(record);
    // The scan caches the CIGAR of each record, and a read doesn't clear it
    let cigar = crate::bam_utils::get_record_cigar(record);
    let ref_start = record.pos();
    let ref_end = std::cmp::max(cigar.end_pos(), ref_start);

    // Alignment positions are given relative to the start of ref_seq
    let ignore_hard_clip = true;
    reference.with_seq(
        chrom_index,
        ref_start as usize,
        ref_end as usize,
        |ref_seq| {
            rust_vc_utils::bam_utils::cigar::get_gap_compressed_identity_from_alignment(
                0,
                cigar.as_slice(),
                &read_seq,
                ref_seq,
                ignore_hard_clip,
            )
        },
    )
}

/// Get the reference sequence corresponding to genome_segment
///
pub fn get_ref_segment_seq(reference: &GenomeRef, genome_segment: &GenomeSegment) -> Vec<u8> {
    let range = &genome_segment.range;
    assert!(range.start >= 0);
    assert!(range.end > range.start);

    reference.get_seq(
        genome_segment.chrom_index,
        range.start as usize,
        range.end as usize,
    )
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_fasta_index() {
        let good = "chr1\t100\t6\t50\t51\nchr2\t20\t120\t20\t21\n";
        let entries = parse_fasta_index(good.as_bytes()).unwrap();
        assert_eq!(entries.len(), 2);
        assert_eq!(entries[1].0, "chr2");
        assert_eq!(entries[1].1, 20);

        // Too few columns, non-numeric length, and line width below line bases
        assert!(parse_fasta_index("chr1\t100\t6\t50\n".as_bytes()).is_none());
        assert!(parse_fasta_index("chr1\tx\t6\t50\t51\n".as_bytes()).is_none());
        assert!(parse_fasta_index("chr1\t100\t6\t50\t49\n".as_bytes()).is_none());
    }

    #[test]
    fn test_is_gzip_file() {
        let dir = std::env::temp_dir().join(format!("sawfish_ref_test_{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let gz = dir.join("a.fa.gz");
        let plain = dir.join("a.fa");
        std::fs::write(&gz, [0x1f, 0x8b, 8, 4]).unwrap();
        std::fs::write(&plain, b">chr1\nACGT\n").unwrap();
        assert!(is_gzip_file(gz.to_str().unwrap()));
        assert!(!is_gzip_file(plain.to_str().unwrap()));
        assert!(!is_gzip_file(dir.join("missing").to_str().unwrap()));

        // The index is older than the fasta only if it was written before it
        let fai = dir.join("a.fa.fai");
        std::fs::write(&fai, b"chr1\t4\t6\t4\t5\n").unwrap();
        let old = std::time::SystemTime::now() - std::time::Duration::from_secs(3600);
        File::options()
            .write(true)
            .open(&fai)
            .unwrap()
            .set_modified(old)
            .unwrap();
        assert!(is_index_older_than_fasta(
            fai.to_str().unwrap(),
            plain.to_str().unwrap()
        ));
        assert!(!is_index_older_than_fasta(
            plain.to_str().unwrap(),
            fai.to_str().unwrap()
        ));
        std::fs::remove_dir_all(&dir).unwrap();
    }

    /// Check that on-demand sequence reads match the full fasta load
    fn check_indexed_ref_matches_full_load(ref_filename: &str) {
        check_indexed_ref_matches_plain_fasta(ref_filename, ref_filename);
    }

    /// Check that on-demand reads from the (possibly compressed) indexed fasta match the full load of
    /// the equivalent plain fasta
    fn check_indexed_ref_matches_plain_fasta(ref_filename: &str, plain_filename: &str) {
        let fai_entries = read_fasta_index(&format!("{ref_filename}.fai")).unwrap();
        let mut genome_ref = get_genome_ref_from_fasta(plain_filename);
        genome_ref.simplify_ambiguous_dna_bases();
        let mut chrom_list = ChromList::default();
        for (label, length, _) in fai_entries.iter() {
            chrom_list.add_chrom(label, *length);
        }
        let indexed_ref = get_indexed_genome_ref(ref_filename, &chrom_list, fai_entries);
        for (chrom_index, chrom_info) in chrom_list.data.iter().enumerate() {
            let expected = genome_ref.chroms.remove(&chrom_info.label).unwrap();
            assert_eq!(indexed_ref.chrom_len(chrom_index), expected.len());
            // Use an odd chunk size so that the last chunk is partial
            let mut chunked = Vec::new();
            indexed_ref
                .for_each_chrom_chunk(chrom_index, 1_000_003, |x| chunked.extend_from_slice(x));
            assert!(chunked == expected, "{}", chrom_info.label);

            // Check partial sequence, both read directly and from the thread's reference window,
            // including windows that start near the end of the previous one
            let len = expected.len();
            let mut ranges = vec![
                (len / 3, std::cmp::min(len / 3 + 1000, len)),
                (len / 3 + 500, std::cmp::min(len / 3 + 700, len)),
                (len.saturating_sub(10), len),
                (0, std::cmp::min(len, 1)),
                (0, len),
                (len / 2, len / 2),
            ];
            for start in [
                REF_WINDOW_SIZE - 5,
                REF_WINDOW_SIZE,
                2 * REF_WINDOW_SIZE - 3,
            ] {
                if start + 50 <= len {
                    ranges.push((start, start + 50));
                }
            }
            for (start, end) in ranges {
                let start = std::cmp::min(start, end);
                assert_eq!(
                    indexed_ref.get_seq(chrom_index, start, end),
                    &expected[start..end],
                    "{}",
                    chrom_info.label
                );
                assert!(
                    indexed_ref.with_seq(chrom_index, start, end, |x| x == &expected[start..end]),
                    "{}",
                    chrom_info.label
                );
            }
        }
    }

    /// Run with: SAWFISH_TEST_REF=ref.fasta cargo test --release -- --ignored compare_indexed_ref
    ///
    #[test]
    #[ignore]
    fn compare_indexed_ref() {
        check_indexed_ref_matches_full_load(&std::env::var("SAWFISH_TEST_REF").unwrap());
    }

    /// Build a fasta and its faidx index with the given chromosomes, each with a line width and
    /// line ending, and with sequence in mixed case and IUPAC codes
    fn write_test_fasta(path: &std::path::Path, chroms: &[(&str, usize, usize, &str)]) {
        let mut state = 0x4a7c_15f3_9e37_79b9u64;
        let mut rand = |n: usize| -> usize {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            (state % n as u64) as usize
        };
        let alphabet = b"ACGTACGTACGTACGTacgtacgtnNRYKMSWBDHV";
        let mut fasta = Vec::new();
        let mut fai = String::new();
        for (label, len, line_bases, eol) in chroms {
            fasta.extend_from_slice(format!(">{label} some description{eol}").as_bytes());
            let offset = fasta.len();
            let seq = (0..*len)
                .map(|_| alphabet[rand(alphabet.len())])
                .collect::<Vec<_>>();
            for line in seq.chunks(*line_bases) {
                fasta.extend_from_slice(line);
                fasta.extend_from_slice(eol.as_bytes());
            }
            fai.push_str(&format!(
                "{label}\t{len}\t{offset}\t{line_bases}\t{}\n",
                line_bases + eol.len()
            ));
        }
        std::fs::write(path, fasta).unwrap();
        std::fs::write(format!("{}.fai", path.display()), fai).unwrap();
    }

    #[test]
    fn test_indexed_ref_matches_full_load() {
        let dir = std::env::temp_dir().join(format!("sawfish_indexed_ref_{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let path = dir.join("ref.fa");
        write_test_fasta(
            &path,
            &[
                ("line60", 1000, 60, "\n"),
                ("crlf7", 100, 7, "\r\n"),
                // The last line is as long as the others, and a single line
                ("exact70", 140, 70, "\n"),
                ("oneline", 50, 50, "\n"),
                ("onebase", 1, 80, "\n"),
                ("onebase_per_line", 40, 1, "\n"),
                // Long enough for several reference windows, with and without CRLF
                ("long80", 2 * REF_WINDOW_SIZE + 12_345, 80, "\n"),
                ("long60crlf", REF_WINDOW_SIZE + 999, 60, "\r\n"),
            ],
        );
        check_indexed_ref_matches_full_load(path.to_str().unwrap());
        std::fs::remove_dir_all(&dir).unwrap();
    }

    /// Compress data as bgzf blocks of block_size uncompressed bytes, followed by the empty EOF block
    fn bgzf_compress(data: &[u8], block_size: usize) -> Vec<u8> {
        use std::io::Write;
        let mut out = Vec::new();
        let chunks = data.chunks(block_size).chain(std::iter::once(&[][..]));
        for chunk in chunks {
            let mut encoder =
                flate2::write::DeflateEncoder::new(Vec::new(), flate2::Compression::default());
            encoder.write_all(chunk).unwrap();
            let payload = encoder.finish().unwrap();
            let mut crc = flate2::Crc::new();
            crc.update(chunk);
            let block_len = BGZF_HEADER_SIZE + 6 + payload.len() + BGZF_TRAILER_SIZE;
            out.extend_from_slice(&[0x1f, 0x8b, 8, 4, 0, 0, 0, 0, 0, 0xff, 6, 0]);
            out.extend_from_slice(&[b'B', b'C', 2, 0]);
            out.extend_from_slice(&((block_len - 1) as u16).to_le_bytes());
            out.extend_from_slice(&payload);
            out.extend_from_slice(&crc.sum().to_le_bytes());
            out.extend_from_slice(&(chunk.len() as u32).to_le_bytes());
        }
        out
    }

    #[test]
    fn test_indexed_bgzip_ref_matches_full_load() {
        let dir = std::env::temp_dir().join(format!("sawfish_bgzf_ref_{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let plain = dir.join("ref.fa");
        write_test_fasta(
            &plain,
            &[
                ("line60", 1000, 60, "\n"),
                ("crlf7", 100, 7, "\r\n"),
                ("onebase", 1, 80, "\n"),
                ("long80", 2 * REF_WINDOW_SIZE + 12_345, 80, "\n"),
            ],
        );

        // The compressed file has a name without any gz extension, and block boundaries which fall
        // inside lines and records
        let compressed = dir.join("ref_compressed");
        let plain_bytes = std::fs::read(&plain).unwrap();
        std::fs::write(&compressed, bgzf_compress(&plain_bytes, 4093)).unwrap();
        std::fs::copy(
            format!("{}.fai", plain.display()),
            format!("{}.fai", compressed.display()),
        )
        .unwrap();

        let compressed_name = compressed.to_str().unwrap();
        assert!(is_gzip_file(compressed_name));
        check_indexed_ref_matches_plain_fasta(compressed_name, plain.to_str().unwrap());

        // Plain gzip without bgzf block structure can't be read at random
        let gzip = dir.join("plain.gz");
        std::fs::write(
            &gzip,
            [
                0x1f, 0x8b, 8, 0, 0, 0, 0, 0, 0, 3, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            ],
        )
        .unwrap();
        assert!(FastaSource::open(gzip.to_str().unwrap()).is_err());

        // A compressed fasta without a usable index is an error, a plain fasta is not, and neither
        // is a compressed fasta with an index
        let compressed_fai = format!("{compressed_name}.fai");
        let fai_entries = read_fasta_index(&compressed_fai);
        assert!(get_missing_index_error(compressed_name, &compressed_fai, &fai_entries).is_none());
        let no_fai_entries = read_fasta_index(&format!("{compressed_name}.missing.fai"));
        assert!(no_fai_entries.is_none());
        let message =
            get_missing_index_error(compressed_name, &compressed_fai, &no_fai_entries).unwrap();
        assert!(message.contains("samtools faidx"), "{message}");
        assert!(
            get_missing_index_error(plain.to_str().unwrap(), "x.fai", &no_fai_entries).is_none()
        );
        std::fs::remove_dir_all(&dir).unwrap();
    }
}
