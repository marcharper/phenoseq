import glob
import sys
from collections import defaultdict
from itertools import combinations
from phenoseq.analyze import parse_args, read_genbank_annots, read_tag_files, map_snps, filter_nonsyn

ordered_tags = "ATCACG CGATGT TTAGGC TGACCA ACAGTG GCCAAT CAGATC ACTTGA GATCAG TAGCTT GGCTAC CTTGTA".split(' ')

sample_pool_dict = dict(zip([(1,), (2,), (3,), (4,), (5,), (6,), (1,2), (1,3), (1,4), (1,5), (1,6), (2,3), (2,4), (2,5), (2,6), (3, 4), (3, 5), (3,6), (4,5), (4,6), (5,6)], range(1,22)))

def load_data(exp_name, include_syn=False):
    if exp_name == "annika":
        tags = ordered_tags[6:]
        geneQualifier = "protein_id"
        gbfile = "ACNL01.gbk"
        fastafile = "ACNL01.fasta"
    elif exp_name == "sam":
        tags = ordered_tags[0:6]
        #geneQualifier = "locus_tag"
        geneQualifier = "protein_id"
        gbfile = "Staphylococcus_aureus_subsp_aureus_NCTC_8325.gbk"
        fastafile = "Staphylococcus_aureus_subsp_aureus_NCTC_8325.fasta"
    else:
        print "Valid experiment names: annika, sam"
        exit()

    tag_pool_dict = dict(zip(tags, range(1, 7)))
    tagFiles = glob.glob("vcf_files/%s/*.vcf" % exp_name)
    print 'reading gene annotations from', gbfile
    annodb, al, dna = read_genbank_annots(gbfile, fastafile=fastafile, geneQualifier=geneQualifier)
    print 'reading tag files:', tagFiles
    snps = read_tag_files(tagFiles, filterFunc=None)
    gsd = map_snps(snps, al, dna)
    if not include_syn:
        gsd = filter_nonsyn(gsd)

    print len(snps), sum([len(v) for v in gsd.values()])

    # snp_dict contains the pools for which each mutation replicates
    snp_dict = defaultdict(list)
    for k, v in gsd.items():
        for s in v:
            snp_dict[(k, s.pos, s.ref, s.alt)].append(tag_pool_dict[s.tag.split('_')[-1]])
    return (tags, tag_pool_dict, annodb, al, dna, snps, gsd, snp_dict)

def pool_histogram_counts(snp_dict):
    items = snp_dict.items()
    items.sort()
    hist = defaultdict(int)
    for k, v in items:
        hist[len(v)] += 1
    keys = hist.keys()
    keys.sort()
    for k in keys:
        print k, hist[k]

def print_sample_pool_map(sample_pool_dict):
    "print sample --> pool mapping"
    for k, v in sorted(sample_pool_dict.items()):
        if len(k) == 1:
            print ",".join(map(str,k)), v
    for k, v in sorted(sample_pool_dict.items()):
        if len(k) == 2:
            print ",".join(map(str,k)), v

if __name__ == '__main__':
    exp_name = sys.argv[1]
    #(tags, tag_pool_dict, annodb, al, dna, snps, gsd, snp_dict) = load_data(exp_name, include_syn=True)
    (tags, tag_pool_dict, annodb, al, dna, snps, gsd, snp_dict) = load_data(exp_name, include_syn=False)
    pool_histogram_counts(snp_dict)

    # Collect SNPs by ordered tuple of pools in which they appear
    pool_snps = defaultdict(list)
    for k, v in sorted(snp_dict.items()):
        pool_snps[tuple(sorted(v))].append(k)
    keys = pool_snps.keys()
    keys.sort()
    # Print "proper SNPs" that replicate in only 1 or 2 lanes (as expected)
    for l in [1,2]:
        for k in keys:
            if len(k) == l:
                strain_id = sample_pool_dict[k]
                for v in pool_snps[k]:
                    print strain_id, " ".join(map(str, v))
    # Print improper SNPs
    for l in [3,4,5]:
        for k in keys:
            if len(k) == l:
                # Figure out the possible strains
                possibles = []
                possibles.extend(k)
                for pair in combinations(k, 2):
                    possibles.append(sample_pool_dict[pair])
                possible_ids = ",".join(map(str,sorted(possibles)))
                strain_id = ",".join(map(str,sorted(k)))
                for v in pool_snps[k]:
                    print strain_id, possible_ids, " ".join(map(str, v))

    # Print reference SNPs occuring in all pools
    for k in keys:
        if len(k) == 6:
            strain_id = "ALL"
            for v in pool_snps[k]:
                print strain_id, " ".join(map(str, v))
    
    # Print list of strain ids in which unambiguous mutations were found
    print list(sorted([sample_pool_dict[x] for x in keys if len(x) <= 2]))
    
    # Print the number of mutations per coding region, counting a mutation from a particular strain just once
    exclude = [2]
    #exclude = []
    gene_dict = defaultdict(set)
    for l in [1,2]:
        for k in keys:
            if len(k) == l:
                strain_id = sample_pool_dict[k]
                if strain_id in exclude:
                    continue
                for v in pool_snps[k]:
                    (gene, _, _, _) = v
                    gene_dict[gene].add(strain_id)

    for k, v in sorted(gene_dict.items()):
        print k, len(v), v
    

