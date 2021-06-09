import pandas as pd
import pickle

def read_dbsnp153(input_path):
    dtype = dict(chrom="string", chromStart="Int64", chromEnd="Int64", name="string", ref="string", alts="string")
    df = pd.read_csv(input_path, header=None, sep="\t", usecols=[0, 1,2,3,4,6], names=['chrom', 'chromStart', 'chromEnd', 'name','ref','alts'], dtype=dtype)
    return df

def read_index_dbsnp153(input_path):
    dtype = dict(chrom="string", chromStart="Int64", chromEnd="Int64", name="string", ref="string", alts="string")
    df = pd.read_csv(input_path, header=None, sep="\t", usecols=[0, 1,2,3,4,6], names=['chrom', 'chromStart', 'chromEnd', 'name','ref','alts'], index_col=['chrom', 'chromStart', 'chromEnd'] ,dtype=dtype)
    return df

def to_dic(input_path):
    bed = read_dbsnp153(input_path)
    db = bed.drop_duplicates(subset=['chrom', 'chromStart', 'chromEnd'], keep=False) 
    result = {}
    for row in db.itertuples():
        chrom = row.chrom
        end_pos = row.chromEnd
        start_pos = row.chromStart
        rsid = row.name
        A1 = row.ref
        A2 = row.alts
        if start_pos == end_pos-1:
            key = (chrom, end_pos)
            result[key] = [rsid, A1, A2]
    return result

# def to_dic_1(input_path):
#     bed = read_index_dbsnp153(input_path)
#     db = bed[~bed.index.duplicated(keep=False)] 
#     result = db.to_dict("index")
#     return result


def save_obj(obj, name ):
    with open(name + '.pkl', 'wb') as f:
        pickle.dump(obj, f, pickle.HIGHEST_PROTOCOL)



def load_obj(obj_path ):
    with open(obj_path, 'rb') as f:
        return pickle.load(f)

# input_path = "dbSnp153_X.bed"

# %time res1 = to_dic(input_path) 
# %time res2 = to_dic_1(input_path) 

from os import listdir
from os.path import isfile, join

mypath="dbSnp153"
onlyfiles = [f for f in listdir(mypath) if isfile(join(mypath, f))]

for f in onlyfiles:
    dic = to_dic(f)
    save_obj(dic, f[:-4])
    