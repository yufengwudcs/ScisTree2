import numpy as np
import pandas as pd 
from .reader import read_vcf, read_csv

class GenotypeProbability():
    def __init__(self, probs, cell_names=None, site_names=None):
        self.probs = probs
        self.nsite, self.ncell = probs.shape
        self.shape = (self.nsite, self.ncell)
        if cell_names:
            self.cell_names = cell_names
        else:
            self.cell_names = [f'c{i}' for i in range(1, self.ncell+1)]
        if site_names:
            self.site_names = site_names
        else:
            self.site_names = [f's{i}' for i in range(1, self.nsite+1)]  

    def __str__(self):
        df = pd.DataFrame(data=np.round(self.probs, 2), columns=self.cell_names, index=self.site_names)
        # return f'{self.ncell} cells, {self.nsite} sites'
        return df.__str__()
    
    def save(self, out_name='probs.csv'):
        df = pd.DataFrame(data=np.round(self.probs, 4), columns=self.cell_names, index=self.site_names)
        df.to_csv(out_name)


def from_probs(probs, cell_names=None, site_names=None, margin=1e-5):
    probs = np.clip(probs, a_min=margin, a_max=1-margin)
    return GenotypeProbability(probs, cell_names=cell_names, site_names=site_names)


def from_reads(reads, ado=0.2, seqerr=0.01, posterior=True, af=None, cell_names=None, site_names=None):
    probs = genotype_probability(reads, ado=ado, seqerr=seqerr, posterior=posterior, af=af)
    return from_probs(probs, cell_names=cell_names, site_names=site_names)


def from_vcf(vcf_path, ado=0.2, seqerr=0.01, posterior=True, af=None, key='AD'):
    # TODO: AF may exist!
    reads, cell_names, site_names = read_vcf(vcf_path, key=key)
    return from_reads(reads, ado=ado, seqerr=seqerr, posterior=posterior, af=af, cell_names=cell_names, site_names=site_names)


def from_csv(csv_path, source='probability', ado=0.2, seqerr=0.01, posterior=True, af=None):
    assert source in ['probability', 'read'], "source should be either 'probability' or 'read'."
    if source == 'probability':
        probs, cell_names, site_names = read_csv(csv_path, reads=False)
        return from_probs(probs, cell_names, site_names)
    else:
        reads, cell_names, site_names = read_csv(csv_path, reads=True)
        return from_reads(reads, ado=ado, seqerr=seqerr, posterior=posterior, af=af, cell_names=cell_names, site_names=site_names)


def concatenate_strings(val1, val2):
    return str(val1) + "|" + str(val2)


def get_mask(ref_counts, alt_counts):
    total_counts = ref_counts + alt_counts
    mask = total_counts == 0
    return mask


def genotype_probability(reads, ado=0.2, seqerr=0.01, posterior=True, af=None):
    ref_counts = reads[:, :, 0]
    alt_counts = reads[:, : ,1]
    mask = get_mask(ref_counts, alt_counts)    
    l00, l01, l11 = likelihood_GATK(ref_counts, alt_counts, ado, seqerr)
    if posterior:
        if af is None:
            af = allele_frequency(l00, l01, l11)
            prob = posterior_probability_GATK(l00, l01, l11, prior_ref=af)
        else:
            prob = posterior_probability_GATK(l00, l01, l11, prior_ref=af)
    else:
        prob = posterior_probability_GATK(l00, l01, l11)
    return np.ma.array(prob, mask=mask).filled(fill_value=0.5)


def allele_frequency(l00, l01, l11):
    """
    Calculate the allele frequency for each site.
    """
    ml_gt = np.argmax(np.concatenate([l00[:, :, np.newaxis], l01[:, :, np.newaxis], l11[:, :, np.newaxis]], axis=-1), axis=-1)
    af = np.mean(ml_gt, axis=-1) / 2
    af = af[:, np.newaxis]
    return 1 - af
    

def posterior_probability_GATK(l00, l01, l11, prior_ref=0.5, margin=1e-5):
    """
    Calculate the posterior probability of each SNV being a true positive.

    We use GATK likelihood with ADO included as described in CellCoal manual (https://dapogon.github.io/cellcoal/cellcoal.manual.v1.1.html#537_genotype_likelihoods)
        P(D|G={g_1, g_2}) = (1-ado)\prod_{i=1}^{r}P(b_i|G={g_1, g_2}) + 0.5*ado[\prod_{i=1}^{r}P(b_i|G={g_1}) + \prod_{i=1}^{r}P(b_i|G={g_2})]
    """
    g00 = prior_ref**2 * l00
    g01 = 2*(1 - prior_ref) * prior_ref * l01
    g11 = (1 - prior_ref)**2 * l11
    g = g00 / (g00 + g01 + g11)
    g = np.clip(np.round(g, 5), a_min=margin, a_max=1-margin)
    return g


def likelihood_GATK(ref_counts, alt_counts, ado=0.2, seqerr=0.01):
    # Q-phred score = -10 * log10(p)
    ref_counts = ref_counts.astype(float)
    alt_counts = alt_counts.astype(float)
    p00, p01, p10, p11 = np.log(1-seqerr), np.log(seqerr), np.log(seqerr), np.log(1-seqerr)
    l00 = np.exp(ref_counts * p00 + alt_counts * p01)
    l01 = (1 - ado) * np.exp(ref_counts * np.log(0.5 * np.exp(p00) + 0.5 * np.exp(p10)) + alt_counts * np.log(0.5 * np.exp(p01) + 0.5 * np.exp(p11))) \
            + (0.5 * ado) * (np.exp(ref_counts * p00 + alt_counts * p10) + np.exp(ref_counts * p10 + alt_counts * p11))
    l11 = np.exp(ref_counts*p10 + alt_counts*p11)
    return l00, l01, l11

