import os 
import numpy as np


def read_vcf(vcf_filepath, key='AD'):
    """
    Reads a VCF file and extracts the Allelic Depth (AD) for each sample at each site.

    Args:
        vcf_filepath (str): The path to the VCF file.

    Returns:
        tuple[list[str], list[list[tuple[int | None, int | None]]]]:
            A tuple containing:
            - A list of sample names.
            - A matrix (list of lists) where each inner list represents a variant site,
              and each element in the inner list is a tuple (ref_reads, alt_reads)
              for a given sample. If AD is not available or malformed for a sample
              at a site, (None, None) will be used.
            Lines where the 'AD' field is not present in the FORMAT column will be skipped.
    """
    sample_names = []
    ad_matrix = []

    assert os.path.exists(vcf_filepath), (f"Error: VCF file not found at '{vcf_filepath}'")

    with open(vcf_filepath, 'r') as f:
        site_names = []
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith('#CHROM'):
                parts = line.split('\t')
                sample_names = parts[9:]
                continue
            if line.startswith('#'):
                continue
            parts = line.split('\t')
            if len(parts) < 9:
                print(f"Warning: Skipping malformed line (too few columns): {line}")
                continue
            format_str = parts[8]
            format_fields = format_str.split(':')
            try:
                ad_index = format_fields.index(key)
            except ValueError:
                print(f"Info: Skipping line because {key} not found in FORMAT: CHR[{parts[0]}], POS[{parts[1]}]")
                continue
            chrom = parts[0]
            pos = parts[1]
            current_site_name = f"{chrom}:{pos}"
            site_names.append(current_site_name)
            site_ad_data = []
            for sample_gt_str in parts[9:]:
                ref_reads, alt_reads = 0, 0
                gt_fields = sample_gt_str.split(':')
                if ad_index < len(gt_fields):
                    ad_value_str = gt_fields[ad_index]
                    if ad_value_str != '.' and ',' in ad_value_str:
                        try:
                            ad_parts = ad_value_str.split(',')
                            ref_reads = int(ad_parts[0])
                            alt_reads = int(ad_parts[1])
                        except (ValueError, IndexError):
                            pass
                site_ad_data.append((ref_reads, alt_reads))
            ad_matrix.append(site_ad_data)
    return np.array(ad_matrix), sample_names, site_names

