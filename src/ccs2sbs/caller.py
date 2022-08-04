import time
import math
import pysam
import bisect
import natsort
import numpy as np
import basic.util
import basic.cslib
import basic.haplib
import basic.bamlib
import basic.vcflib
import multiprocessing as mp
from dataclasses import dataclass
from typing import Dict, List, Tuple
from collections import defaultdict, Counter


def update_allelecounts(
    read,
    tpos2allelecounts: Dict[int, np.ndarray],
    tpos2qbase2bq_lst: Dict[int, Dict[int, List[int]]],
    tpos2qbase2read_lst: Dict[int, Dict[int, List[str]]]
) -> None:

    tpos2qbase = {}
    tpos = read.tstart
    qpos = read.qstart
    basic.cslib.cs2tuple(read)
    if read.mapq >= 1 and read.is_primary:
        for cstuple in read.cstuple_lst:
            state, ref, alt, ref_len, alt_len, = cstuple
            if state == 1:  # match
                for i, alt_base in enumerate(alt):
                    tpos2allelecounts[tpos + i + 1][basic.util.base2idx[alt_base]] += 1
                    tpos2qbase[tpos + i + 1] = (alt_base, read.bq_int_lst[qpos + i])
                    tpos2qbase2read_lst[tpos + i + 1][basic.util.base2idx[alt_base]].append(read.qname)
                    tpos2qbase2bq_lst[tpos + i + 1][basic.util.base2idx[alt_base]].append(read.bq_int_lst[qpos + i])
            elif state == 2:  # sub
                tpos2allelecounts[tpos + 1][basic.util.base2idx[alt]] += 1
                tpos2qbase[tpos + 1] = (alt, read.bq_int_lst[qpos])
                tpos2qbase2read_lst[tpos + 1][basic.util.base2idx[alt]].append(read.qname)
                tpos2qbase2bq_lst[tpos + 1][basic.util.base2idx[alt]].append(read.bq_int_lst[qpos])
            elif state == 3:  # insertion
                tpos2allelecounts[tpos + 1][4] += 1
                tpos2qbase2read_lst[tpos + 1][4].append(read.qname)
                tpos2qbase2bq_lst[tpos + 1][4].append(read.bq_int_lst[qpos])
            elif state == 4:  # deletion
                for j in range(len(ref[1:])):
                    tpos2allelecounts[tpos + j + 1][5] += 1
                    tpos2qbase[tpos + j + 1] = ("-", 0)
                    tpos2qbase2read_lst[tpos + j + 1][5].append(read.qname)
                    tpos2qbase2bq_lst[tpos + j + 1][5].append(0)
            tpos += ref_len
            qpos += alt_len
    else:
        tpos2qbase = basic.cslib.cs2tpos2qbase(read)
    return tpos2qbase


def get_sbs_allelecounts(
    tpos: int,
    ref: str,
    alt: str,
    tpos2allelecounts: Dict[int, np.ndarray],
    tpos2qbase2bq_lst: Dict[int, Dict[int, List[int]]]
) -> Tuple[int, int, str, float, int, int, int]:
    
    ins_count = tpos2allelecounts[tpos][4]
    del_count = tpos2allelecounts[tpos][5]
    total_count = sum(tpos2allelecounts[tpos]) - ins_count
    ref_count = tpos2allelecounts[tpos][basic.util.base2idx[ref]]
    alt_count = tpos2allelecounts[tpos][basic.util.base2idx[alt]]  
    bq = "{:.1f}".format(sum(tpos2qbase2bq_lst[tpos][basic.util.base2idx[alt]]) / float(alt_count))
    vaf = alt_count / float(total_count)
    return bq, vaf, ref_count, alt_count, ins_count, del_count, total_count


def get_somatic_substitutions(
    chrom: str,
    chrom_len: int,
    bam_file: str,
    vcf_file: str,
    phased_vcf_file: str,
    common_snps: str,
    panel_of_normals: str,
    loci_lst: List[Tuple[str, int, int]],
    min_mapq: int,
    min_trim: float,
    qlen_mean: int,
    qlen_lower_limit: int,
    qlen_upper_limit: int,
    min_sequence_identity: float,
    min_hq_base_proportion: float,
    min_alignment_proportion: float,
    min_bq: int,
    mismatch_window: int,
    max_mismatch_count: int,
    min_ref_count: int,
    min_alt_count: int,
    md_threshold: int,
    phase: bool,
    ploidy: str,
    min_hap_count: int,
    tree_of_life_sample: bool,
    create_panel_of_normals: bool,
    chrom2tsbs_lst: Dict[str, List[List[Tuple[str, int, str, str, int, int, int, int, str]]]],
    chrom2tsbs_statistics: Dict[str, List[int]],
) -> List[Tuple[str, int, str, str, int, int, int, float, float]]:

    tsbs_lst = []
    if vcf_file.endswith(".vcf"):
        sample_snp_set = basic.vcflib.load_snp(chrom, vcf_file)
        
    common_snp_set = set()
    if not tree_of_life_sample and common_snps.endswith(".vcf"):
        common_snp_set = basic.vcflib.load_common_snp(chrom, common_snps)
       
    pon_sbs_set = set()
    if not create_panel_of_normals and not tree_of_life_sample and panel_of_normals.endswith(".vcf"):
        pon_sbs_set, _, = basic.vcflib.load_pon(chrom, panel_of_normals)
         
    if phase and ploidy == "diploid":
        (
            hpos_lst,
            hblock_lst,
            hetsnp_lst,
            hidx2hetsnp,
            hidx2hstate,
            hetsnp2bidx,
            hetsnp2hidx,
        ) = basic.vcflib.get_phased_hetsnps(phased_vcf_file, chrom, chrom_len)
  
    alignments = pysam.AlignmentFile(bam_file, "rb")
    for loci in loci_lst:
        chunkloci_lst = basic.util.chunkloci(loci)
        for chunkloci in chunkloci_lst:
            tsbs_candidate_lst = []
            read2tpos2qbase = defaultdict(dict)
            chunk_start, chunk_end = chunkloci[1:]
            tpos2allelecounts = defaultdict(lambda: np.zeros(6)) 
            tpos2qbase2bq_lst = defaultdict(lambda: {0: [], 1:[], 2:[], 3:[], 4:[], 5:[]})
            tpos2qbase2read_lst = defaultdict(lambda: {0: [], 1:[], 2:[], 3:[], 4:[], 5:[]})
            if vcf_file.endswith(".bgz"):
                sample_snp_set = basic.vcflib.load_bgz_snp((chrom, chunk_start - qlen_mean, chunk_end + qlen_mean), vcf_file)
            
            if not tree_of_life_sample and common_snps.endswith(".bgz"):
                common_snp_set = basic.vcflib.load_bgz_common_snp((chrom, chunk_start - qlen_mean, chunk_end + qlen_mean), common_snps)
                
            if not create_panel_of_normals and not tree_of_life_sample and panel_of_normals.endswith(".bgz"):
                pon_sbs_set, _ = basic.vcflib.load_bgz_pon((chrom, chunk_start - qlen_mean, chunk_end + qlen_mean), panel_of_normals) 
                
            for line in alignments.fetch(*chunkloci):
                read = basic.bamlib.BAM(line)
                read2tpos2qbase[read.qname] = update_allelecounts(read, tpos2allelecounts, tpos2qbase2bq_lst, tpos2qbase2read_lst)

                if not read.is_primary:
                    continue

                if read.mapq < min_mapq:
                    continue

                if read.qlen < qlen_lower_limit or read.qlen > qlen_upper_limit:
                    continue

                if read.query_alignment_proportion < min_alignment_proportion:
                    continue 

                if basic.bamlib.get_hq_base_proportion(read) < min_hq_base_proportion:
                    continue

                if basic.util.get_blast_sequence_identity(read) < min_sequence_identity:
                    continue

                basic.cslib.cs2mut(read)
                contamination_state, read_tsbs_lst, read_qsbs_lst, read_sbs_bq_lst = basic.util.get_sbs_candidates(read, sample_snp_set, common_snp_set)
                if contamination_state:
                    continue
                    
                if len(read_tsbs_lst) == 0 and len(read.tdbs_lst) == 0:
                    continue
                
                trimmed_qstart = math.floor(min_trim * read.qlen)
                trimmed_qend = math.ceil((1 - min_trim) * read.qlen)
                mismatch_lst = natsort.natsorted(read_tsbs_lst + read.mismatch_lst)
                mismatch_set = set(mismatch_lst)
                mpos_lst = [mismatch[1] for mismatch in mismatch_lst]
                for i, (tsbs, qsbs) in enumerate(zip(read_tsbs_lst, read_qsbs_lst)):
                    bq = read_sbs_bq_lst[i]
                    if bq < min_bq:
                        continue
                    
                    if not tree_of_life_sample and not create_panel_of_normals:
                        if tsbs in pon_sbs_set:
                            continue

                    tpos = tsbs[1]
                    qpos = qsbs[1]
                    if qpos < trimmed_qstart:
                        continue
                    elif qpos > trimmed_qend:
                        continue

                    mismatch_start, mismatch_end = basic.util.get_mismatch_range(tpos, qpos, read.qlen, mismatch_window)
                    idx = bisect.bisect_left(mpos_lst, mismatch_start)
                    jdx = bisect.bisect_right(mpos_lst, mismatch_end)
                    if tsbs in mismatch_set:
                        mismatch_count = (jdx - idx) - 1 
                    else:
                        mismatch_count = jdx - idx 

                    if mismatch_count > max_mismatch_count:
                        continue
                    tsbs_candidate_lst.append(tsbs)
                    
            tsbs2count = Counter(tsbs_candidate_lst)     
            for tsbs, count in tsbs2count.items():
                _, tpos, ref, alt = tsbs 
                if tpos <= chunk_start or tpos >= chunk_end:
                    continue
                
                bq, vaf, ref_count, alt_count, ins_count, del_count, total_count = get_sbs_allelecounts(tpos, ref, alt, tpos2allelecounts, tpos2qbase2bq_lst)
                if del_count != 0 or ins_count != 0:
                    continue

                if total_count > md_threshold:
                    continue

                if ref_count >= min_ref_count and alt_count >= min_alt_count:
                    if not phase: # unphased single base substitutions
                        tsbs_lst.append(
                            (chrom, tpos, ref, alt, "PASS", bq, total_count, ref_count, alt_count, vaf, ".")
                        )
                    else:
                        bidx2count = defaultdict(lambda: 0)
                        alt_hap2count = defaultdict(lambda: 0)
                        ref_read_lst = tpos2qbase2read_lst[tpos][basic.util.base2idx[ref]]
                        alt_read_lst = tpos2qbase2read_lst[tpos][basic.util.base2idx[alt]]
                        for alt_read in alt_read_lst:
                            bidx, alt_hap = basic.haplib.get_read_haplotype(
                                read2tpos2qbase[alt_read],
                                hpos_lst,
                                hetsnp_lst,
                                hidx2hstate,
                                hetsnp2bidx,
                                hetsnp2hidx,
                            )
                            bidx2count[bidx] += 1
                            alt_hap2count[alt_hap] += 1

                        if len(bidx2count.keys()) == 1 and len(alt_hap2count.keys()) == 1:
                            bidx = list(bidx2count.keys())[0]
                            if bidx == ".":
                                continue
                            
                            alt_hap = list(alt_hap2count.keys())[0]
                            if alt_hap == ".":
                                continue

                        ref_hap = "1" if alt_hap == "0" else "0"
                        hap2count = basic.haplib.get_region_hap2count( 
                            ref_read_lst, read2tpos2qbase, hblock_lst[bidx], hidx2hetsnp,
                        )
                        ref_hap_count = hap2count[ref_hap]
                        alt_hap_count = hap2count[alt_hap] 
                        if ref_hap_count >= min_hap_count and alt_hap_count >= min_hap_count:
                            phase_set = hidx2hetsnp[hblock_lst[bidx][0][0]][1]
                            tsbs_lst.append((chrom, tpos, ref, alt, "PASS", bq, total_count, ref_count, alt_count, vaf, phase_set))
                            
    chrom2tsbs_lst[chrom] = natsort.natsorted(tsbs_lst)
    alignments.close()


def call_somatic_substitutions(
    bam_file: str,
    vcf_file: str,
    phased_vcf_file: str,
    region: str,
    region_list: str,
    ploidy: str,
    min_mapq: int,
    min_trim: float,
    min_sequence_identity: float,
    min_hq_base_proportion: float,
    min_alignment_proportion: float,
    common_snps: str,
    panel_of_normals: str,
    min_bq: int,
    mismatch_window: int,
    max_mismatch_count: int,
    min_ref_count: int,
    min_alt_count: int,
    min_hap_count: int,
    threads: int,
    create_panel_of_normals: bool,
    tree_of_life_sample: bool,
    phase: bool,
    version: str,
    out_file: str,
) -> None:

    cpu_start = time.time() / 60
    basic.util.check_caller_input_exists(
        bam_file,
        vcf_file,
        phased_vcf_file,
        region,
        region_list,
        common_snps,
        panel_of_normals,
        phase,
        out_file,
    )

    _, tname2tsize = basic.bamlib.get_tname2tsize(bam_file)
    chrom_lst, chrom2loci_lst = basic.util.load_loci(region, region_list, tname2tsize)
    qlen_mean, qlen_lower_limit, qlen_upper_limit, md_threshold = basic.bamlib.get_thresholds(
        bam_file, chrom_lst, tname2tsize
    )
    vcf_header = basic.vcflib.get_basic_vcf_header(
        bam_file,
        vcf_file,
        phased_vcf_file,
        region,
        region_list,
        ploidy,         
        tname2tsize,
        common_snps,
        panel_of_normals,
        min_mapq,
        min_trim,
        min_sequence_identity,
        min_hq_base_proportion,
        min_alignment_proportion,
        min_bq,
        mismatch_window,
        max_mismatch_count,
        min_ref_count,
        min_alt_count,
        md_threshold,
        min_hap_count, 
        threads,
        create_panel_of_normals,
        tree_of_life_sample,
        phase,
        version,
        out_file,
    )
    print("calling substitutions with {} threads".format(threads))
    p = mp.Pool(threads)
    manager = mp.Manager()
    chrom2tsbs_lst = manager.dict()
    chrom2tsbs_statistics = manager.dict()
    get_somatic_substitutions_arg_lst = [
        (
            chrom,
            tname2tsize[chrom],
            bam_file,
            vcf_file,
            phased_vcf_file,
            common_snps,
            panel_of_normals,
            chrom2loci_lst[chrom],
            min_mapq,
            min_trim,
            qlen_mean,
            qlen_lower_limit,
            qlen_upper_limit,
            min_sequence_identity,
            min_hq_base_proportion,
            min_alignment_proportion,
            min_bq,
            mismatch_window,
            max_mismatch_count,
            min_ref_count,
            min_alt_count,
            md_threshold,
            phase,
            ploidy,
            min_hap_count,
            chrom2tsbs_lst,
            chrom2tsbs_statistics,
        )
        for chrom in chrom_lst
    ]
    p.starmap(
        get_somatic_substitutions, get_somatic_substitutions_arg_lst,
    )
    p.close()
    p.join()
    basic.vcflib.dump_basic_sbs(chrom_lst, chrom2tsbs_lst, phase, vcf_header, out_file)  
    basic.vcflib.dump_basic_statistics(chrom_lst, chrom2tsbs_statistics, out_file.replace(".vcf", ".stats"))
    print("finished calling and returning substitutions with {} threads".format(threads))
    cpu_end = time.time() / 60
    duration = cpu_end - cpu_start
    print("single molecule somatic mutation detection took {} minutes".format(duration))
