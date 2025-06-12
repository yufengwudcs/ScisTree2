import numpy as np

def tree_accuracy(tree1, tree2):
    count = 0
    splits1 = tree1.get_splits(return_label=True, contains_leaf='1')
    splits2 = tree2.get_splits(return_label=True, contains_leaf='1')
    for split in splits1:
        count += int(split in splits2)
    return count / len(splits1)

def split_accuracy(splits1, splits2):
    count = 0
    for split in splits1:
        count += int(split in splits2)
    return count / len(splits1)

def genotype_accuarcy(geno1, geno2):
    return np.mean(geno1 == geno2)


def is_covered(clade1, clade2):
    diff = clade1 - clade2
    if 0 in diff and 1 in diff and -1 in diff:
        return 0 # no relationship
    if 1 in diff:
        return 1 # clade1(snp1) happens before clade2(snp2)
    else:
        return -1 # clade2(snp2) happens before clade1(snp1)

def get_ancestor_descendant_pairs(geno):
    num_snp, _ = geno.shape
    mutations = {f's{_}': [] for _ in range(num_snp)}
    for i in range(num_snp):
        for j in range(i+1, num_snp):
            code = is_covered(geno[i], geno[j])
            if code == 1:
                mutations[f's{i}'].append(f's{j}')
            elif code == -1:
                mutations[f's{j}'].append(f's{i}')
    return mutations

def ancestor_descendant_error(mutation1, mutation2):
    count = 0
    total = 0
    for snp1 in mutation1:
        for snp2 in mutation1[snp1]:
            count += int(snp2 not in mutation2[snp1])
            total += 1
    # total = len(mutation1) * (len(mutation1) - 1) / 2
    return count / total

def different_lineage_error(mutation1, mutation2):
    count = 0
    total = 0
    muts = list(mutation1.keys())
    num_mut = len(muts)
    for i in range(num_mut):
        for j in range(i+1, num_mut):
            if muts[i] not in mutation1[muts[j]] and muts[j] not in mutation1[muts[i]]:
                count += int(muts[j] in mutation2[muts[i]]) + int(muts[i] in mutation2[muts[j]])
                total += 1
    return count / total