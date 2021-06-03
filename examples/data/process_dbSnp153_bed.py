import numpy as np
import pandas as pd
import pickle

def read_dbsnp153(input_path, output_path):
    df = pd.read_csv(input_path, sep="\t", header=None)
    header=["chrom","chromStart",	"chromEnd",	"name",	"ref",	"altCount",	"alts",	"shiftBases",	"freqSourceCount",	"minorAlleleFreq",	"majorAllele",	"minorAllele",	"maxFuncImpact",	"class",	"ucscNotes",	"_dataOffset",	"_dataLen"]
    # print(len(header))
    df.columns = header[:len(df.columns)]
    return df[["chrom", "chromStart" ,"chromEnd", "name", "ref", "alts"]]

def build_dict(df):
    result = {}
    for row in df.itertuples():
        chrom = row.chrom
        end_pos = row.chromEnd
        start_pos = row.chromStart
        rsid = row.name
        A1 = row.ref
        A2 = row.alts
        if start_pos == end_pos - 1:
            key = (chrom, end_pos)
            result[key] = [rsid, A1, A2]
    return result

"""Function to save python data structure on disk
    Args:
        obj (obj): the data structure/object to be saved on disk.
        name (str): the name for the obj to be saved as.
    Returns:
        return nothing
"""
def save_obj(obj, name ):
    with open(name + '.pkl', 'wb') as f:
        pickle.dump(obj, f, pickle.HIGHEST_PROTOCOL)




"""Function to load saved python data structure from disk

    Args:
        name (str): the name of the saved obj on disk to be loaded
    Returns:
        return the loaded object/ data structure
"""
def load_obj(obj_path ):
    with open(obj_path, 'rb') as f:
        return pickle.load(f)

if __name__ == "__main__":
    df = read_dbsnp153("dbSnp153.bed", "")
    print(df)
    result = build_dict(df)
    # print(result)
    print(len(result))
    save_obj(result, "dbSnp153")