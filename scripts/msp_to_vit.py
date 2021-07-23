## Converts RFmix 2.3 output file .msp.tsv to RFmix 1.5 output format "Viterbi.txt"
## Usage: python msp_to_vit.py chr21.msp.tsv americans_subset_remove.chr21.bim chr21.viterbi.tsv

import sys
import pandas as pd


# read in input files to pandas dataframe
msp = pd.read_csv(sys.argv[1], sep='\t', header=[1,1])
bim = pd.read_csv(sys.argv[2], sep='\t', header = None)
outfile = sys.argv[3]

# make list of bim file snp positions
bim_snps = list(bim[3])

# this function accepts a snp position, and iterates through each
# row of the msp file to check if that position is inside the CRF range of that row
# when it finds the correct range, it returns the corresponding index in the msp file
def check_range(snp):
    nomatch= "NA"
    for x in range(len(msp)):
        if (int(snp) >= int(msp.iloc[x,1])) & (int(snp) <= int(msp.iloc[x,2])):
            return x
        else:
            continue
    return nomatch

# build a list of rows indeces from the msp that correspond to the location of the snp in the bim file
# if there is a snp from the bim file missing from the msp, exit
frames = [check_range(f) for f in bim_snps ]
if "NA" in frames:
    raise Exception("a snp is missing from msp file")

## next: for each value in frames, append [value,6:len(msp) -1] from msp file
def get_msp_values(index):
    list = msp.iloc[index,6:(len(msp.columns)-1)]
    int_list = [int(i) + 1 for i in list]
    return int_list

vit_list = [get_msp_values(i) for i in frames ]

vit_df = pd.DataFrame(vit_list)

vit_df.to_csv(outfile, sep = " ", index=False, header =False)
