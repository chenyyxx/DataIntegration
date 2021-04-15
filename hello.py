import numpy as np
import pandas as pd
from pyliftover import LiftOver

def read_data(input_path, output_path):
    df = pd.read_csv(input_path, compression='gzip', header=0, sep='\t',quotechar='"')
    result = df.loc[:,["#chrom","pos", "alt", "ref", "maf", "beta", "sebeta", "pval"]]
    res = result.rename({"#chrom":"Chr","pos":"BP", "alt":"A1", "ref":"A2", "maf":"EAF", "beta":"Beta", "sebeta":"Se", "pval":"P"},axis="columns")
    return res[0:100]

df = read_data("finngen_R4_AB1_ARTHROPOD.gz", "")

def modify_data(df, input_format, output_format):
    lo = LiftOver(input_format, output_format)
    modified_output = pd.DataFrame()
    unconvertible_output = pd.DataFrame()
    print(df)
    for i in range(df.shape[0]):
        # print(df.iloc[i])
        chrom = "chr" + str(df.iloc[i]["Chr"])
        # print(chrom)
        pos = df.iloc[i]["BP"]
        # print(pos)
        modified =lo.convert_coordinate(chrom, pos)
        print(modified)
        # if not modified: # mark as NA, and output a summary table noting down the uncessful convert, number of row containing rsid, not containing rsid
            # print(chrom)
            # print(pos)
            # print(i)
        # else:  # there might be cases when after conversion, its empty => cannot be converted
        new_chrom = modified[0][0]
        new_pos = modified[0][1]

        
        df["Chr"][i] = new_chrom
        df["BP"][i] = new_pos
        # print(df[i])
        # modified_output.append(df.iloc[i,:])
            
            # print(new_chrom)
            # print(new_pos)
    # print(modified_output)


modify_data(df, "hg38", "hg19")

# file_path_1 = ""
# file_path_finngen = "finngen_R4_AB1_ARTHROPOD.gz"
#in_file_1 = pd.read_csv(file_path_1)
# df = pd.read_csv(file_path_finngen, compression='gzip', header=0, sep='\t',quotechar='"')
# maybe utilize tbi

# result = df.loc[:,["#chrom","pos", "alt", "ref", "maf", "beta", "sebeta", "pval"]]
# res = result.rename({"#chrom":"Chr","pos":"BP", "alt":"A1", "ref":"A2", "maf":"EAF", "beta":"Beta", "sebeta":"Se", "pval":"P"},axis="columns")

# lo = LiftOver('hg19', 'hg38')
#print(lo.convert_coordinate('chr1', 1000000))


# save file 


# create new pd df

# add column into new df: concat (rearrange columns order)

# merge with file 1

# return new df

#print(res.head(10))