import argparse
import pandas as pd
import subprocess
import os
import re
import sys

# ================= CLI SETUP =================
parser = argparse.ArgumentParser(description="Remap manifest probes to a new reference genome.")
parser.add_argument("-i", "--manifest_file", type=str, required=True, help="Path to input manifest CSV")
parser.add_argument("-r", "--reference_genome", type=str, required=True, help="Path to reference genome FASTA")
parser.add_argument("-o", "--output_file", type=str, required=True, help="Path to write the updated manifest CSV")
args = parser.parse_args()

MANIFEST_FILE = args.manifest_file
REF_FASTA = args.reference_genome
OUTPUT_FILE = args.output_file
# =============================================

def get_header_and_footer_info(filename, start_marker='[Assay]', end_marker='[Controls]'):
    """Locates the data section within the Illumina manifest."""
    header_line = 0
    footer_line = None
    with open(filename, 'r') as f:
        for i, line in enumerate(f):
            stripped = line.strip()
            if stripped == start_marker:
                header_line = i + 1
            if stripped == end_marker:
                footer_line = i
                break
    nrows = None
    if footer_line is not None:
        nrows = footer_line - header_line
    return header_line, nrows

def parse_cigar_to_ref_pos(start_pos, cigar, query_index):
    """
    Maps a 0-based index in the Query sequence to a 1-based Coordinate in the Reference
    using the CIGAR string. Used for TopGenomicSeq fallback.
    """
    ops = re.findall(r'(\d+)([MIDNSHP=X])', cigar)
    curr_q = 0 
    curr_r = start_pos 
    
    for length, op in ops:
        length = int(length)
        if op in ['M', '=', 'X']:
            if query_index < curr_q + length:
                return curr_r + (query_index - curr_q)
            curr_q += length
            curr_r += length
        elif op in ['I', 'S']:
            if query_index < curr_q + length:
                # Target is inside insertion/clip; return junction
                return curr_r 
            curr_q += length
        elif op in ['D', 'N']:
            curr_r += length
    return curr_r

def get_probe_coordinate(pos_start, cigar_str, strand, assay_type):
    """
    Calculates the Variant Start Coordinate based on Probe Alignment.
    Standard: Start of Variant = Probe End + 1 (for Inf II).
    """
    ops = re.findall(r'(\d+)([MIDNSHP=X])', cigar_str)
    parsed_ops = [(int(n), op) for n, op in ops]
    
    if strand == '+':
        # Probe matches Ref 5'->3'. 
        # Physical 3' end = Start + Span - 1.
        ref_span = sum([n for n, op in parsed_ops if op in ['M', 'D', 'N', '=', 'X', 'S']])
        probe_end = pos_start + ref_span - 1
        
        # For Inf II, variant is the base AFTER the probe.
        # For Inf I, variant is the LAST base of the probe.
        return probe_end + 1 if assay_type == 'II' else probe_end
    else:
        # Probe matches RC(Ref). 
        # Physical 3' end corresponds to Alignment Start (POS).
        # Adjust for leading soft clips if any.
        leading_S = parsed_ops[0][0] if parsed_ops and parsed_ops[0][1] == 'S' else 0
        probe_end = pos_start - leading_S
        
        # For Inf II on Minus, "After" the probe is the lower coordinate.
        return probe_end - 1 if assay_type == 'II' else probe_end

def calculate_overlap(r1_start, r1_end, r2_start, r2_end):
    """Calculates overlap length between two genomic intervals."""
    return max(0, min(r1_end, r2_end) - max(r1_start, r2_start))

def get_alignment_span(pos, cigar):
    """Returns the genomic end position of an alignment."""
    ops = re.findall(r'(\d+)([MIDNSHP=X])', cigar)
    span = sum([int(n) for n, op in ops if op in 'MDN=X'])
    return pos + span - 1

def extract_candidates(top_seq):
    """Extracts Pre, AlleleA, AlleleB, Post from TopGenomicSeq."""
    match = re.search(r'(.*?)\[(.*?)/(.*?)\](.*)', top_seq)
    if match:
        return match.group(1), match.group(2), match.group(3), match.group(4)
    return None, None, None, None

def run_remapping():
    print(f"Reading manifest: {MANIFEST_FILE}...")
    skip_rows, nrows = get_header_and_footer_info(MANIFEST_FILE)
    dtype_dict = { "AddressA_ID": "string", "AddressB_ID": "string", "GenomeBuild": "string", "Chr": "string", "MapInfo": "Int64", "SourceVersion": "string", "BeadSetID": "string"}
    df = pd.read_csv(MANIFEST_FILE, dtype=dtype_dict, skiprows=skip_rows, nrows=nrows, low_memory=False)
    # Filter out empty rows if any
    df = df.dropna(subset=['Name', 'AlleleA_ProbeSeq'])


    print("Generating sequence files...")
    topseq_fasta = 'temp_topseq.fasta'
    probe_fasta = 'temp_probes.fasta'
    
    candidates_info = {} 
    
    with open(topseq_fasta, 'w') as ft, open(probe_fasta, 'w') as fp:
        for idx, row in df.iterrows():
            name = row['Name']
            # 1. Probe Fasta
            fp.write(f">{name}\n{row['AlleleA_ProbeSeq']}\n")
            
            # 2. TopSeq Fasta (Allele A and B candidates)
            top_seq = row.get('TopGenomicSeq', '')
            pre, a, b, post = extract_candidates(top_seq)
            if pre is not None:
                # Construct candidates. Remove '-' for alignment.
                seq_a = pre + (a if a != '-' else '') + post
                seq_b = pre + (b if b != '-' else '') + post
                
                ft.write(f">{name}_A\n{seq_a}\n")
                ft.write(f">{name}_B\n{seq_b}\n")
                
                candidates_info[name] = {
                    'PreLen': len(pre), 
                    'PostLen': len(post), 
                    'AlleleA': a, 
                    'AlleleB': b
                }

    print("Aligning sequences with minimap2...")
    topseq_sam = 'temp_topseq.sam'
    probe_sam = 'temp_probe.sam'
    
    # Align TopSeq (Context)
    subprocess.check_call(f"minimap2 -ax sr {REF_FASTA} {topseq_fasta} > {topseq_sam} 2> /dev/null", shell=True)
    # Align Probes (Coordinates) - allow secondary alignments (-N 5) for overlap checking
    subprocess.check_call(f"minimap2 -ax sr -N 5 {REF_FASTA} {probe_fasta} > {probe_sam} 2> /dev/null", shell=True)

    print("Parsing TopGenomicSeq Alignments (Ref/Alt Authority)...")
    topseq_results = {}
    
    def get_nm(cols):
        for tag in cols[11:]:
            if tag.startswith('NM:i:'): return int(tag.split(':')[2])
        return 999

    temp_top_results = {}
    with open(topseq_sam, 'r') as f:
        for line in f:
            if line.startswith('@') or line.startswith('[M'): continue
            cols = line.split('\t')
            if int(cols[1]) & 4: continue 
            
            qname_full = cols[0]
            qname = qname_full[:-2]
            which = qname_full[-1]
            flag = int(cols[1])
            
            if qname not in temp_top_results: temp_top_results[qname] = {}
            temp_top_results[qname][which] = {
                'NM': get_nm(cols), 'Chr': cols[2], 'Pos': int(cols[3]), 
                'Cigar': cols[5], 'MAPQ': int(cols[4]), 'Strand': '-' if flag & 16 else '+'
            }

    # Determine Winner (Ref/Alt)
    for nm, res in temp_top_results.items():
        if nm not in candidates_info: continue
        info = candidates_info[nm]
        
        res_a = res.get('A')
        res_b = res.get('B')
        nm_a = res_a['NM'] if res_a else 999
        nm_b = res_b['NM'] if res_b else 999
        
        # Winner is the one with lower Edit Distance (NM)
        winner = res_b if nm_b < nm_a else (res_a if res_a else res_b)
        
        # Set Ref/Alt based on winner
        final_ref = info['AlleleB'] if winner == res_b else info['AlleleA']
        final_alt = info['AlleleA'] if winner == res_b else info['AlleleB']
            
        if winner:
            span_end = get_alignment_span(winner['Pos'], winner['Cigar'])
            topseq_results[nm] = {
                'Ref': final_ref, 'Alt': final_alt,
                'Chr': winner['Chr'], 'Start': winner['Pos'], 'End': span_end,
                'MAPQ': winner['MAPQ'], 'Cigar': winner['Cigar'],
                'Strand': winner['Strand'], 'PreLen': info['PreLen'], 'PostLen': info['PostLen']
            }

    print("Parsing Probe Alignments (Coordinate Authority)...")
    probe_candidates = {}
    assay_types = {}
    
    # Pre-load Assay Types
    for idx, row in df.iterrows():
        is_ii = pd.isna(row.get('AlleleB_ProbeSeq', float('nan')))
        assay_types[row['Name']] = 'II' if is_ii else 'I'

    with open(probe_sam, 'r') as f:
        for line in f:
            if line.startswith('@') or line.startswith('[M'): continue
            cols = line.split('\t')
            flag = int(cols[1])
            if flag & 4: continue 
            pos = int(cols[3])
            cigar = cols[5]
            probe_candidates.setdefault(cols[0], []).append({
                'Chr': cols[2], 'Pos': pos, 'Cigar': cigar,
                'Strand': '-' if flag & 16 else '+', 
                'MAPQ': int(cols[4]), 'End': get_alignment_span(pos, cigar)
            })

    print("Resolving final coordinates...")
    new_rows = {}
    
    for idx, row in df.iterrows():
        nm = row['Name']
        ts_res = topseq_results.get(nm)
        
        c_chr, c_pos, c_strand, c_ref, c_alt = '0', 0, 'N/A', 'N', 'N'
        c_mapq_ts, c_mapq_pb = 0, 0
        
        if ts_res:
            c_ref, c_alt = ts_res['Ref'], ts_res['Alt']
            c_mapq_ts = ts_res['MAPQ']
            c_chr = ts_res['Chr']
            c_strand = ts_res['Strand'] 
            
            # 1. Try to find a Probe that matches TopSeq Strand and Overlaps
            selected_pb = None
            max_overlap = 0
            
            for pb in probe_candidates.get(nm, []):
                if pb['Chr'] != ts_res['Chr']: continue
                ov = calculate_overlap(ts_res['Start'], ts_res['End'], pb['Pos'], pb['End'])
                if ov > max_overlap:
                    max_overlap = ov
                    selected_pb = pb
            
            if selected_pb and selected_pb['Strand'] == c_strand:
                # Case A: Probe is consistent. Use it for high-precision coordinate.
                c_mapq_pb = selected_pb['MAPQ']
                assay = assay_types.get(nm, 'II')
                raw_pos = get_probe_coordinate(selected_pb['Pos'], selected_pb['Cigar'], c_strand, assay)
                
                # Deletion Correction (Minus Strand Only)
                # If Deletion (Ref > Alt) and Strand is -, we must subtract the deletion length
                # because 'raw_pos' points to the High Coordinate (Start of Event in Reverse).
                is_deletion = len(c_ref) > len(c_alt)
                if is_deletion and c_strand == '-':
                    del_len = len(c_ref) - len(c_alt)
                    raw_pos -= del_len
                
                c_pos = raw_pos
            else:
                # Case B: Probe missing or disagrees on Strand. Fallback to TopSeq.
                # If Strand is +, Variant starts after 'Pre' (index = PreLen).
                # If Strand is -, Variant starts after 'Post' (index = PostLen).
                if c_strand == '+':
                    target_idx = ts_res['PreLen']
                else:
                    target_idx = ts_res['PostLen']
                
                c_pos = parse_cigar_to_ref_pos(ts_res['Start'], ts_res['Cigar'], target_idx)
                
        new_rows[idx] = {
            'Chr_EquCab3': c_chr, 'MapInfo_EquCab3': c_pos, 'Strand_EquCab3': c_strand,
            'Ref_EquCab3': c_ref, 'Alt_EquCab3': c_alt,
            'MAPQ_TopGenomicSeq': c_mapq_ts, 'MAPQ_Probe': c_mapq_pb
        }

    # Write output
    cols = ['Chr_EquCab3', 'MapInfo_EquCab3', 'Strand_EquCab3', 'Ref_EquCab3', 'Alt_EquCab3', 'MAPQ_TopGenomicSeq', 'MAPQ_Probe']
    for c in cols: df[c] = pd.NA
    for idx, data in new_rows.items():
        for k,v in data.items(): df.at[idx, k] = v
        
    print(f"Saving to {OUTPUT_FILE}...")
    df.to_csv(OUTPUT_FILE, index=False)
    
    # Cleanup
    #for f in [topseq_fasta, probe_fasta, topseq_sam, probe_sam]:
    #    if os.path.exists(f): os.remove(f)
    print("Done.")

if __name__ == "__main__":
    run_remapping()
