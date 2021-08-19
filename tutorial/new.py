import pandas as pd
import pickle
import time
from dataintegrator import DataIntegrator as di
from os import listdir
from os.path import isfile, join



#code snippet
def main(df,):
    cur_chr = "1"
    cur_dict_name = "dbSnp153_" + cur_chr + ".pkl"
    cur_dict = load_obj(cur_dict_name)
    for row in df.itertuples():
        chrom = row.Chr
        if chrom != cur_chr: # reach a new chr
            cur_chr = chrom
            cur_dict_name = "dbSnp153_" + cur_chr + ".pkl"
            cur_dict = None
            cur_dict = load_obj(cur_dict_name)
        pos = row.BP
        rs_id = row.SNP
        key = ("chr"+chrom, pos)
        res = cur_dict[key]
        if res != None:
            #do sth
            pass
        else:
            #mark sth
            pass

def read_dbsnp153(input_path):
    dtype = dict(chrom="string", chromStart="Int64", chromEnd="Int64", name="string", ref="string", alts="string")
    df = pd.read_csv(input_path, header=None, sep="\t", usecols=[0, 1,2,3,4,6], names=['chrom', 'chromStart', 'chromEnd', 'name','ref','alts'], dtype=dtype)
    return df

def merge_add_rsid(df):
    new_df = df
    new_df["start"] = new_df.BP - 1
    new_df["Chr"] = new_df["Chr"].astype("string")
    new_df["chrom"] = "chr"+new_df.Chr
    print(new_df)
    mypath="data/dbSnp153"
    onlyfiles = [f for f in listdir(mypath) if isfile(join(mypath, f))]
    # print(onlyfiles)
    bed = read_dbsnp153(mypath + "/" + onlyfiles[0])
    print(bed)
    result = bed.merge(new_df, how="inner", left_on=["chrom", "chromStart", "chromEnd"] ,right_on=["chrom", "start", "BP"])
    print(result)
    for f in onlyfiles[1:]:
        bed = read_dbsnp153(mypath + "/" +f)
        # print(bed)
        chr_merge = bed.merge(new_df, how="inner", left_on=["chrom", "chromStart", "chromEnd"] ,right_on=["chrom", "start", "BP"])
        # print(chr_merge)
        result = pd.concat([result, chr_merge]).reset_index(drop=True)
        print(result)


    return result

def new_add_rsid(df,keep_all=False, inplace=False, show_comment=False, show_errors=False): 
    added_rsid = []
    comment = []
    cur_chr = "1"
    cur_dict_name = "data/dbSnp153"+"dbSnp153_" + cur_chr + ".pkl"
    cur_dict = load_obj(cur_dict_name)
    for row in df.itertuples():
        chrom = row.Chr
        if chrom != cur_chr: # reach a new chr
            cur_chr = chrom
            cur_dict_name = "data/dbSnp153"+"dbSnp153_" + cur_chr + ".pkl"
            cur_dict = None
            cur_dict = load_obj(cur_dict_name)
        pos = row.BP
        rs_id = row.SNP
        key = ("chr"+chrom, pos)
        if key in cur_dict: # the row in the data can be found in dbSnp153
            result_list = cur_dict[key]
            data_rs_id = result_list[0] 
            if pd.isna(rs_id): # rs_id is absence in original dataset
                added_rsid.append(data_rs_id)
                comment.append("added")
            elif rs_id == data_rs_id: # if rs_id in original dataset is the same as dnSnp153
                added_rsid.append(rs_id)
                comment.append("same")
            else: # find different rsid in dbSnp153, update with new
                added_rsid.append(data_rs_id)
                comment.append("different")
        else:
            added_rsid.append(pd.NA)
            comment.append("key not found")
    
    result = df.assign(added_rsid = added_rsid)

    if show_errors:
        if inplace or show_comment or keep_all:
            print("The `show_errors` flag cannot be used with the `inplace`, `comment`, and `keep_all` flags together.")
            return
        result = result.assign(comment=comment)
        mask = pd.isna(result["added_rsid"])
        result = result[mask]
        return result

    if inplace:
        result = result[["Chr", "BP" ,"added_rsid", "A1", "A2", "EAF", "Beta", "Se", "P"]].rename({"added_rsid": "SNP"},axis="columns")

    if show_comment:
        result = result.assign(comment=comment)


    if not keep_all:
        if inplace:
            result  = result.dropna(subset=["SNP"]).reset_index(drop=True)
        else:
            result  = result.dropna(subset=["added_rsid"]).reset_index(drop=True)

    

    return result
        



def load_obj(obj_path ):
    with open(obj_path, 'rb') as f:
        return pickle.load(f)



if __name__ == "__main__":
    input_path = "data/finngen_R4_AB1_ARTHROPOD.gz"
    # 1. Read Data
    # ----------------------------------------------------------------------------------------------------------------------------------------------------
    print("start reading data")
    df = di.read_data(input_path, "#chrom","pos", "rsids" ,"alt", "ref", "maf", "beta", "sebeta", "pval")
    print(df.shape)
    print("data successfully read")

    # 2. Filter Biallelic
    # ----------------------------------------------------------------------------------------------------------------------------------------------------
    print("start filtering bi-allelic cases")
    bi_allelic = di.filter_bi_allelic(df)

    # 3. drop duplicates
    # ----------------------------------------------------------------------------------------------------------------------------------------------------
    print("start deduplicating")
    dedup_df = di.deduplicate(bi_allelic)

    # 4. sort
    # ----------------------------------------------------------------------------------------------------------------------------------------------------
    print("start sorting")
    sorted_df = di.sort_by_chr_bp(dedup_df)

    # 6. query dbSnp153 for required information
    # ----------------------------------------------------------------------------------------------------------------------------------------------------
    print("start getting required info from dbSnp153")
    C = time.time()
    # dbSnp153 = di.query_data(sorted_df, "data/dbSnp153.bb")
    # di.save_obj(dbSnp153, "obj/dbSnp153_finngen")
    dbSnp153 = di.load_obj("obj/dbSnp153_finngen.pkl") # once save you can load it
    print("end querying")
    E = time.time()

     # 7 data processing (e.g.): Add rsID
    # ----------------------------------------------------------------------------------------------------------------------------------------------------
    print("start processing data: add rsID")
    added_rsid = di.add_rsid(sorted_df, dbSnp153)
    D = time.time()
    print(added_rsid)

    added_rsid = new_add_rsid(sorted_df)

    F = time.time()
    print(added_rsid)


    # added_rsid = di.add_rsid(sorted_df, dbSnp153, show_comment=True)
    # print("comment true")
    # print(added_rsid)

    # added_rsid = di.add_rsid(sorted_df, dbSnp153, inplace=True)
    # print("inplace true")
    # print(added_rsid)

    # added_rsid = di.add_rsid(sorted_df, dbSnp153, inplace=True, show_comment=True)
    # print("comment inplace true")
    # print(added_rsid)

    # added_rsid = di.add_rsid(sorted_df, dbSnp153, show_errors=True)
    # print("show_erros true")
    # print(added_rsid)

    print("Time used:" , D-C)
    print("Time used (query):" , E-C)
    print("Time used (process):" , D-E)
    print("Time used (process):" , F-D)