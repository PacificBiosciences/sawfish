#!/usr/bin/env python3

import os, sys

def get_args():
    import argparse

    parser = argparse.ArgumentParser(
        description="Process GA4GH single-sample output vcfs from truvari into summary performance metrics")

    parser.add_argument("--query-vcf",
                        required=True,
                        metavar="FILE",
                        help="Input query GA4GH VCF (can be gzipped)")
    parser.add_argument("--truth-vcf",
                        required=True,
                        metavar="FILE",
                        help="Input truth GA4GH VCF (can be gzipped)")

    args = parser.parse_args()

    def error_exit(msg):
        raise Exception(msg)

    def check_required_file(filename, label):
        if not os.path.isfile(filename):
            error_exit(f"Can't find {label} file '{filename}'")

    check_required_file(args.query_vcf, "input query vcf")
    check_required_file(args.truth_vcf, "input truth vcf")

    return args

class VCFID :
    CHROM = 0
    POS = 1
    REF = 3
    ALT = 4
    QUAL = 5
    FILTER = 6
    INFO = 7
    FORMAT = 8
    SAMPLE = 9

def get_info_key_val(info_string,key) :
    """
    Parse the info field to find the given key. Return its value if found, or None if not found
    """
    import re
    match=re.search(b"%s=([^;\t]*);?" % (key) ,info_string)
    if match is None : return None
    return match.group(1);

def get_sample_key_val(word, key, sample_index = 0) :
    """
    Parse the format and sample fields to find the given key. Return its value if found, or None if not found
    """

    def get_format_tags(format_string):
        format_tags=format_string.split(b":")
        assert(len(format_tags))
        if format_tags[0] == "." : format_tags = []
        return format_tags

    sample_col_index=VCFID.SAMPLE+sample_index
    if len(word) <= sample_col_index:
        return None
    format_tags = get_format_tags(word[VCFID.FORMAT])
    key_index = format_tags.index(key)
    if key_index < 0:
        return None
    return word[sample_col_index].split(b":")[key_index]

def open_maybe_gz_file(filename, mode='r'):
    """
    Open an optionally gzipped file
    """
    import gzip
    if filename.endswith('.gz'):
        return gzip.open(filename, mode)
    else:
        return open(filename, mode)

def safe_frac(a,b):
    if b==0:
        return None
    else:
        return a/b

def get_sv_size(word):
    """
    Return absolute SV size
    """
    svlen = get_info_key_val(word[VCFID.INFO],b"SVLEN")
    if svlen is not None:
        return abs(int(svlen))

    return abs(len(word[VCFID.ALT]) - len(word[VCFID.REF]))

class SVType :
    DEL = 0
    INS = 1

def get_sv_type(word):
    svt = get_info_key_val(word[VCFID.INFO],b"SVTYPE")
    if svt is not None:
        if svt == b"DEL":
            return SVType.DEL
        elif svt == b"INS":
            return SVType.INS
        else:
            return None
    else:
        alt_len=len(word[VCFID.ALT])
        ref_len=len(word[VCFID.REF])

        if alt_len > ref_len:
            return SVType.INS
        else:
            return SVType.DEL

def process_vcf(vcf_filename, stats, is_truth):
    for line in open_maybe_gz_file(vcf_filename):
        # Skip header:
        if line.startswith(b"#"):
            continue

        word = line.strip().split(b"\t")
        bd = get_sample_key_val(word, b"BD")
        svlen = get_sv_size(word)
        svtype = get_sv_type(word)
        stats.update(svlen, svtype, bd, is_truth)

class TruthStats:
    def __init__(self):
        self.tp_base = 0
        self.tp_comp = 0
        self.fp = 0
        self.fn = 0
        self.tp_len = 0
        self.fp_len = 0
        self.fn_len = 0

    def update(self, bd, is_truth, size):
        if bd is None:
            return

        if is_truth:
            if bd == b"TP":
                self.tp_base += 1
            elif bd == b"FN":
                self.fn += 1
                self.fn_len += size
        else:
            if bd == b"TP":
                self.tp_comp += 1
                self.tp_len += size
            elif bd == b"FP":
                self.fp += 1
                self.fp_len += size

    def report(self, fp):
        fp.write("TP_base:  {} TP_comp: {} FP: {} FN: {}\n".format(self.tp_base, self.tp_comp, self.fp, self.fn))
        fp.write("TP_comp_bp: {} FP_bp: {} FN_bp: {}\n".format(self.tp_len, self.fp_len, self.fn_len))
        r=safe_frac(self.tp_base, self.tp_base+self.fn)
        p=safe_frac(self.tp_comp, self.tp_comp+self.fp)
        if p is None or r is None :
            f1 = None
        else:
            f1=safe_frac(2*p*r, p+r)

        def format_metric(m):
            if m is None:
                return "N/A   ";
            else:
                return f"{m:.4f}"

        fp.write("recall: {} prec: {} f1: {}\n".format(format_metric(r),format_metric(p),format_metric(f1)))
        fp.write("\n")

class SizeStratifiedTruthStats:
    def __init__(self, depth_cuts):
        assert(len(depth_cuts) > 1)
        self.depth_cuts = depth_cuts
        self.ins_depth_stats = [TruthStats() for _ in range(len(depth_cuts))]
        self.del_depth_stats = [TruthStats() for _ in range(len(depth_cuts))]
        self.indel_depth_stats = [TruthStats() for _ in range(len(depth_cuts))]
        self.total_stats = TruthStats()

    def update(self, size, type, bd, is_truth):
        self.total_stats.update(bd, is_truth, size)

        def update_indel_bins(i, size, bd, is_truth):
            self.indel_depth_stats[i].update(bd, is_truth, size)
            if type == SVType.DEL:
                self.del_depth_stats[i].update(bd, is_truth, size)
            elif type == SVType.INS:
                self.ins_depth_stats[i].update(bd, is_truth, size)

        for i in range(len(self.depth_cuts)-1):
            cut0=self.depth_cuts[i]
            cut1=self.depth_cuts[i+1]
            if size >= cut0 and size < cut1:
                update_indel_bins(i, size, bd, is_truth)
                return

        update_indel_bins(-1, size, bd, is_truth)

    def report(self, fp):
        fp.write("Total Stats:\n")
        self.total_stats.report(fp)

        def write_size_results(i, size_label):
            fp.write(f"Size {size_label} INDEL Stats:\n")
            self.indel_depth_stats[i].report(fp)
            fp.write(f"Size {size_label} DEL Stats:\n")
            self.del_depth_stats[i].report(fp)
            fp.write(f"Size {size_label} INS Stats:\n")
            self.ins_depth_stats[i].report(fp)

        for i in range(len(self.depth_cuts)-1):
            cut0=self.depth_cuts[i]
            cut1= self.depth_cuts[i+1]
            size_label=f"[{cut0},{cut1})"
            write_size_results(i, size_label)

        final_cut=self.depth_cuts[-1]
        size_label=f"{final_cut}+"
        write_size_results(-1, size_label)

def main():
    args = get_args()

    # Set first cut to zero because truvari has already handled the minimum size while producing the labeled VCF output:
    stats = SizeStratifiedTruthStats([0,500,5000])
    process_vcf(args.truth_vcf, stats, True)
    process_vcf(args.query_vcf, stats, False)
    stats.report(sys.stdout)

main()
