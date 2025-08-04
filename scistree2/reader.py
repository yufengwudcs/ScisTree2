import os 
import numpy as np
import pandas as pd 

def read_vcf(vcf_filepath, key='AD'):
    """
    Reads a VCF file and extracts the Allelic Depth (AD) for each sample at each site.

    Args:
        vcf_filepath (str): The path to the VCF file.

    Returns:
        tuple[list[str], list[list[tuple[int | None, int | None]]]]:
            A tuple containing:
            - A matrix (list of lists) where each inner list represents a variant site,
              and each element in the inner list is a tuple (ref_reads, alt_reads)
              for a given sample. If AD is not available or malformed for a sample
              at a site, (None, None) will be used.
            Lines where the 'AD' field is not present in the FORMAT column will be skipped.
            - A list of sample names.
            - A list of site names.
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


def convert_2d_string_array_to_3d(input_2d_array):
    """
    Converts a 2D array (list of lists) of strings formatted as 'x|y'
    back into a 3D NumPy array of shape (rows, cols, 2).

    Args:
        input_2d_array (list): A 2D list of strings, where each string is 'x|y'.

    Returns:
        np.ndarray: A 3D NumPy array with the last dimension being of size 2.
                    Returns an empty array if the input list is empty.
    """
    num_rows = len(input_2d_array)
    num_cols = len(input_2d_array[0]) if num_rows > 0 else 0
    result_3d_array = np.empty((num_rows, num_cols, 2), dtype=object)

    for i in range(num_rows):
        for j in range(num_cols):
            parts = input_2d_array[i][j].split('|')
            if len(parts) == 2:
                try:
                    x_val = int(parts[0])
                    y_val = int(parts[1])
                    result_3d_array[i, j] = [x_val, y_val]
                except ValueError:
                    print(f"Warning: Could not convert '{input_2d_array[i][j]}' to integers. Storing as original string parts.")
                    result_3d_array[i, j] = parts
            else:
                print(f"Warning: Unexpected format for string '{input_2d_array[i][j]}'. Expected 'x|y'.")
                result_3d_array[i, j] = [None, None]
    return result_3d_array


def read_csv(csv_filepath, reads=False):
    """
    Reads a CSV file.
    Args:
        csv_filepath (str): The path to the CSV file.

    Returns:
        tuple[list[str], list[list[tuple[int | None, int | None]]]]:
            A tuple containing:
            - A matrix (list of lists) where each inner list represents a variant site,
              and each element in the inner list is a tuple (ref_reads, alt_reads).
            - A list of sample names.
            - A list of site names.
    """
    assert os.path.exists(csv_filepath), (f"Error: CSV file not found at '{csv_filepath}'")
    df = pd.read_csv(csv_filepath, index_col=0)
    sample_names = df.columns.to_list()
    site_names = df.index.to_list()
    return convert_2d_string_array_to_3d(df.values) if reads else df.values, sample_names, site_names



if __name__ == "__main__":
    a = read_csv('../tutorials/data/toy_probs.csv', reads=False)
    print(a)