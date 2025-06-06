import numpy as np


def genotype_probability(reads, ado=0.2, seqerr=0.01, posterior=True):
    reads = reads.astype(float)
    ref_counts = reads[:, :, 0]
    alt_counts = reads[:, :, 1]
    l00, l01, l11 = likelihood_GATK(ref_counts, alt_counts, ado, seqerr)
    if posterior:
        af = allele_frequency(l00, l01, l11)
        prob = posterior_probability_GATK(l00, l01, l11, prior_ref=af)
    else:
        prob = posterior_probability_GATK(l00, l01, l11)
    return prob

"""
 Calculate the allele frequency for each site.
"""
def allele_frequency(l00, l01, l11):
    ml_gt = np.argmax(np.concatenate([l00[:, :, np.newaxis], l01[:, :, np.newaxis], l11[:, :, np.newaxis]], axis=-1), axis=-1)
    af = np.mean(ml_gt, axis=-1) / 2
    af = af[:, np.newaxis]
    return 1 - af
    

"""
Calculate the posterior probability of each SNV being a true positive.

We use GATK likelihood with ADO included as described in CellCoal manual (https://dapogon.github.io/cellcoal/cellcoal.manual.v1.1.html#537_genotype_likelihoods)
    P(D|G={g_1, g_2}) = (1-ado)\prod_{i=1}^{r}P(b_i|G={g_1, g_2}) + 0.5*ado[\prod_{i=1}^{r}P(b_i|G={g_1}) + \prod_{i=1}^{r}P(b_i|G={g_2})]
"""
def posterior_probability_GATK(l00, l01, l11, prior_ref=0.5, margin=1e-5):
    g00 = prior_ref**2 * l00
    g01 = 2*(1 - prior_ref) * prior_ref * l01
    g11 = (1 - prior_ref)**2 * l11
    g = g00 / (g00 + g01 + g11)
    g = np.clip(np.round(g, 5), a_min=margin, a_max=1-margin)
    return g

"""
xxx
"""
def likelihood_GATK(ref_counts, alt_counts, ado=0.2, seqerr=0.01):
    # Q-phred score = -10 * log10(p)
    p00, p01, p10, p11 = np.log(1-seqerr), np.log(seqerr), np.log(seqerr), np.log(1-seqerr)
    z = ref_counts * p00 + alt_counts * p01
    # print(type(z))
    l00 = np.exp(ref_counts * p00 + alt_counts * p01)
    l01 = (1 - ado) * np.exp(ref_counts * np.log(0.5 * np.exp(p00) + 0.5 * np.exp(p10)) + alt_counts * np.log(0.5 * np.exp(p01) + 0.5 * np.exp(p11))) \
            + (0.5 * ado) * (np.exp(ref_counts * p00 + alt_counts * p10) + np.exp(ref_counts * p10 + alt_counts * p11))
    l11 = np.exp(ref_counts*p10 + alt_counts*p11)
    return l00, l01, l11

