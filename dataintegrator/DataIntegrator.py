# libraries
import numpy as np
import pandas as pd
from pyliftover import LiftOver
import io
import os
import pyBigWig
import pickle
import time





"""Function to read and format data for process

    Args:
        input_path (str): The path of the input data.
        Chr_col_name (str): The name of the column in the input data representing chromosome name.
        BP_col_name (str): The name of the column in the input data representing base pair position.
        SNP_col_name (str): The name of the column in the input data representing rs ID.
        A1_col_name (str): The name of the column in the input data representing effect allele.
        A2_col_name (str): The name of the column in the input data representing non-effect allele.
        EAF_col_name (str): The name of the column in the input data representing effect size.
        Beta_col_name (str): The name of the column in the input data representing beta.
        Se_col_name (str): The name of the column in the input data representing standard deviation.
        P_col_name (str): The name of the column in the input data representing p-value.
        separate_by (str): How the input data is separated. Default to "\t" (tab separated).

    Returns:
        pandas.Data.Frame: return formatted data in the form of pandas Data.Frame

"""

def read_data( input_path, Chr_col_name, BP_col_name, SNP_col_name, A1_col_name, A2_col_name, EAF_col_name, Beta_col_name, Se_col_name, P_col_name, separate_by="\t"):
    raw_df = pd.read_csv(input_path, compression='gzip', header=0, sep=separate_by,quotechar='"')
    # print(raw_df)
    result = raw_df.loc[:,[Chr_col_name, BP_col_name, SNP_col_name, A1_col_name, A2_col_name, EAF_col_name, Beta_col_name, Se_col_name, P_col_name]]
    res = result.rename(
        {
            Chr_col_name:"Chr",
            BP_col_name:"BP", 
            SNP_col_name:"SNP", 
            A1_col_name:"A1", 
            A2_col_name:"A2", 
            EAF_col_name:"EAF", 
            Beta_col_name:"Beta", 
            Se_col_name:"Se", 
            P_col_name:"P"
        },axis="columns")
    dtype = dict(Chr="string", BP='Int64', SNP="string", A1="string", A2="string", EAF=float, Beta=float, Se=float, P=float)
    res = res.astype(dtype)
    res["Chr"] = res["Chr"].str.upper()
    res["Chr"] = res["Chr"].apply(lambda y: "X" if y=="23" else("Y" if y=="24" else y))
    res["A1"] = res["A1"].str.upper()
    res["A2"] = res["A2"].str.upper()
    res["SNP"] = res["SNP"].str.lower()
    return res





"""Function to filter only bi-allelic cases in the data

    Args:
        df (pandas.Data.Frame): The data frame to be filtered.
        rest (boolean): value indicating wether or not to keep (mark only) the non-bi-allelic cases. Default to False.

    Returns:
        pandas.Data.Frame: return filtered data in the form of pandas Data.Frame.

"""

def filter_bi_allelic(df, rest=False):
    len_mask = (df['A1'].str.len() == 1) & (df['A2'].str.len() == 1)
    val_mask = (df['A1'] != "I") & (df['A1'] != "D") & (df['A1'] != "R") & (df['A2'] != "I") & (df['A2'] != "D") & (df['A2'] != "R")
    mask = len_mask & val_mask
    if not rest:
        result = df[mask].reset_index(drop=True)
        return result
    else:
        result = df[~mask].reset_index(drop=True)
        return result






"""Function to drop rows in data containing dduplicate keys (Chr + BP)

    Args:
        df (pandas.Data.Frame): The data frame to be deduplicated.

    Returns:
        pandas.Data.Frame: return filtered data in the form of pandas Data.Frame.

"""

def deduplicate(df):
    result = df.drop_duplicates(subset=['Chr', 'BP'], keep=False)
    return result


"""Function to sort the data based on Chr and BP

    Args:
        df (pandas.Data.Frame): the data to be sorted
    Returns:
        pandas.Data.Frame: return the sorted data
"""
def sort_by_chr_bp(df):
    def mixs(v):
        try:
            return int(v)
        except ValueError:
            return v
    df = df.assign(chr_numeric = lambda x: x['Chr'].apply(lambda y: 22 if y=="X" else(23 if y=="Y" else int(y))))
    result = df.sort_values(by=["chr_numeric", "BP"]).drop(['chr_numeric'], axis=1).reset_index(drop=True)
    return result




"""Function to query required data from dbSnp153

    Args:
        df (pandas.Data.Frame): the data we want more info
        link (str): path or link of the '.bb' file of dbSnp153 
    Returns:
        pandas.Data.Frame: return complete information from dbSnp153 as a python dictionary
"""
# link = "http://hgdownload.soe.ucsc.edu/gbdb/hg38/snp/dbSnp153.bb"

def query_data(df, link="http://hgdownload.soe.ucsc.edu/gbdb/hg38/snp/dbSnp153.bb"):
    bb = pyBigWig.open(link)
    result = {}
    set_list = []
    for row in df.itertuples():
        chrom = "chr" + str(row.Chr)
        end_pos = row.BP
        start_pos =end_pos - 1
        # print(chrom, start_pos, end_pos)
        dat = bb.entries(chrom, start_pos, end_pos)
        if dat != None:  
            for i in dat:
                reference_start = i[0]
                reference_end = i[1]
                raw_string = i[2]
                if reference_start == start_pos and reference_end == end_pos:
                    key = (chrom[3:], reference_end)
                    result[key] = raw_string
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












"""Function to create query data for lift over

    Args:
        input_version (str): the genome build of the original data.
        output_version (str): the desired genome build you want to lift over to.
    Returns:
        python dictionary: return the input version, output version and the liftover chain file.

"""

def create_lo(input_version, output_version):
    lo = LiftOver(input_version, output_version)
    return {"input_version": input_version, "output_version":output_version, "lo":lo}





"""Function to lift over genome build

    Args:
        df (pandas.Data.Frame): the data to be lifted over
        lo_dict (python dictionary): the lift over dictionary return from the create_lo function
        keep_unconvertible (boolean): if true, the function will keep and mark the rows that are not convertible. Default to False.
        keep_original_version (boolean): if true, the function will keep the Chr + BP of original genome build. Default to False.
    Returns:
        pandas.Data.Frame: return the data being lifted over to the desired genome build
"""
def lift_over(df, lo_dict, keep_unconvertible=False, keep_original_version= False):
    reference_table = _lift_over_basic(df, lo_dict)
    result = _lift_over_merge(df, reference_table)
    if not keep_unconvertible:
        new_chr_name = reference_table.columns[2]
        new_pos_name = reference_table.columns[3]
        result = result.dropna(subset=[new_chr_name]).reset_index(drop=True)
    if not keep_original_version:
        new_chr_col_name = lo_dict['output_version']+"_chr"
        new_pos_col_name = lo_dict['output_version']+"_pos"
        result = result[[new_chr_col_name, new_pos_col_name, "SNP", "A1", "A2", "EAF", "Beta", "Se", "P"]].rename({new_chr_col_name:"Chr", new_pos_col_name:"BP"}, axis="columns")
    return result





"""Function to query and add rs ID for rows missing rsIDs.

    Args:
        df (pandas.Data.Frame): the data to be added rs_ids
        data (python dictionary): the dictionary containing required info from dbSnp153
    Returns:
        pandas.Data.Frame: return the data being added rs_ids.
"""
def add_rsid(df, data):
    added_rs_id = []
    for row in df.itertuples():
        chrom = row.Chr
        pos = row.BP
        rs_id = row.SNP
        key = (chrom, pos)
        if key in data: # the row in the data can be found in dbSnp153
            raw_string = data[key]
            parsed_string = raw_string.split('\t')
            data_rs_id = parsed_string[0] 
            if pd.isna(rs_id): # rs_id is absence in original dataset
                added_rs_id.append(data_rs_id)
                # print("find none")
            elif rs_id == data_rs_id: # if rs_id in original dataset is the same as dnSnp153
                added_rs_id.append(rs_id)
                # print("same")
            else: # find different rsid in dbSnp153, update with new
                added_rs_id.append(data_rs_id)
                # print("different")
        else:
            added_rs_id.append("key not found")
            # print("key not found")
    result = df.assign(added_rs_id = added_rs_id)
    # print(result)
    return result





"""Function to flip the input data to forward strand

    Args:
        df (pandas.Data.Frame): the data to be flipped to forward strand
        data (python dictionary): the dictionary containing required info from dbSnp153
        keep_unconvertible (boolean): if true, the function will keep and mark the rows that are not flipped. Default to False.

    Returns:
        pandas.Data.Frame: return the data being flipped to forward strand
"""
def flip_strand( df, data, keep_all=False):
    flipped_A1 = []
    flipped_A2 =  []
    comment = []
    for row in df.itertuples():
        chrom = row.Chr
        pos = row.BP
        A1 = row.A1
        A2 = row.A2
        key = (chrom, pos)
        if key in data: # check if key in dnSnp153
            cur_set = {A1, A2}
            raw_string = data[key]
            parsed_string = raw_string.split("\t")
            data_a1 = [i for i in parsed_string[1] if i != ","]
            data_a2 = [i for i in parsed_string[3] if i != ","]

            if len(data_a1) == 1 and len(data_a2) == 1:
                cur_set.add(data_a1[0])
                cur_set.add(data_a2[0])
                # print(cur_set)
                if len(cur_set) == 4: # flip
                    new_a1 = _flip(A1)
                    new_a2 = _flip(A2)
                    flipped_A1.append(new_a1)
                    flipped_A2.append(new_a2)
                    comment.append("flipped")
                elif len(cur_set) == 2: # do not flip
                    # print(i)
                    flipped_A1.append(A1)
                    flipped_A2.append(A2)
                    comment.append("keep original")
                else: # mark: what is this case? => original data T/C, dbsnp153 C/A: 10  94958283  rs111998500
                    flipped_A1.append("1")
                    flipped_A2.append("1")
                    print(data_a1)
                    print(data_a2)
                    # print(parsed_string[3])
                    # print(parsed_string[3].split(","))
                    comment.append("mark")
            else: # tri-alleic snps -> mark
                flipped_A1.append("2")
                flipped_A2.append("2")
                comment.append("dbSnp153: Indel")
        else: # key not found
            flipped_A1.append("3")
            flipped_A2.append("3")
            comment.append("Key not found") 

    result = df.assign(new_A1 = flipped_A1)
    result = result.assign(new_A2 = flipped_A2)
    result = result.assign(comment = comment)
    # print(result)
    if keep_all:
        return result
    else:
        return result.query('new_A1 != "3" & new_A1 != "2" & new_A1 != "1"').reset_index(drop=True)




"""Function to align effect allele
    this function will align the effect allele of input data based on a reference data

    Args:
        reference (pandas.Data.Frame): the reference table
        df (pandas.Data.Frame): the data to be aligned
        check_error_rows (boolean): if true, the function will output the rows that cannot be aligned. Default to False.
    Returns:
        pandas.Data.Frame: return the data with its effect allele being aligned with the reference table.
"""

def align_effect_allele( reference, df, check_error_rows=False):
    reference = reference[["Chr", "BP", "A1", "A2"]].rename({"A1":"reference_A1", "A2":"reference_A2"}, axis="columns")
    process = df[["Chr", "BP", "A1", "A2"]].rename({"A1":"process_A1", "A2":"process_A2"}, axis="columns")
    merge_table = pd.merge(process, reference, on=["Chr", "BP"], how="inner")
    if len(merge_table) == 0:
        print("reference data and process data have no records in common. Please check data source.")
        return
    nochange_mask = (merge_table["process_A1"] == merge_table["reference_A1"]) & (merge_table["process_A2"] == merge_table["reference_A2"])
    align_mask = (merge_table["process_A1"] == merge_table["reference_A2"]) & (merge_table["process_A2"] == merge_table["reference_A1"])
    error_mask = ~nochange_mask & ~align_mask

    key_to_nochange = merge_table[nochange_mask][["Chr", "BP"]]
    key_to_align = merge_table[align_mask][["Chr", "BP"]]
    key_to_error = merge_table[error_mask][["Chr", "BP"]]
    # print(key_to_error)

    nochange = pd.merge(df, key_to_nochange, on=["Chr", "BP"], how="inner")
    align = pd.merge(df, key_to_align, on=["Chr", "BP"], how="inner")
    aligned = _swap_effect_allele(align)
    # print("aligned")
    # print(aligned)
    error = pd.merge(df, key_to_error, on=["Chr", "BP"], how="inner")
    # print(error)
    # print(aligned)
    result = nochange.append(aligned).reset_index(drop=True)
    # print(result)
    sorted_result = sort_by_Chr(result)
    if check_error_rows:
        return pd.merge(error, merge_table[error_mask], on=["Chr", "BP"], how="inner")
    print(str(nochange.shape[0]) + " rows were left unchanged (already aligned)")
    print(str(aligned.shape[0]) + " rows were aligned successfully")
    print(str(error.shape[0]) + " rows failed to align, dropped from result! Set the check_error_rows flag to True to view them.")
    return sorted_result
        



        
"""Function to save the processed data in gz or csv

    Args:
        output_path (str): the path you want the data to be saved.
        df (pandas.Data.Frame): the processed data to be saved.
        name (str): the output name of the data.
        save_format (str): the saving format. Choose between 'gzip' or 'csv'. Default to gz.
    Returns:
        pandas.Data.Frame: return filtered data in the form of pandas Data.Frame

"""

def save_data(output_path, df, name, save_format="gzip"):
    if save_format == "gzip":
        df_out = output_path + "/" + name +".gz"
        try:
            df.to_csv(df_out, compression='gzip')
            return "successfully save"
        except:
            return "fail to save data"
    elif save_format == "csv": # csv
        df_out = output_path + "/" + name + ".csv"
        try:
            df.to_csv(df_out)
            return "successfully save"
        except:
            return "fail to save data"
    else:
        print("format not accepted, use 'gzip' or 'csv' for the `save_format` argument")
        return







# ---------------------------------------------------------------------------------------------
# Helper Functions

# helper function to flip strand for one row
def _flip( allele):
    new_allele = ""
    if allele == "A":
        new_allele = "T"
    if allele == "T":
        new_allele = "A"
    if allele == "C":
        new_allele = "G"
    if allele == "G":
        new_allele = "C"
    return new_allele

# helper function to swap effect allele and align effect size

def _swap_effect_allele( df):
    col_list = list(df)
    col_list[3], col_list[4] = col_list[4], col_list[3]
    df.columns = col_list
    df["Beta"] = -1 * df["Beta"]
    df["EAF"] = 1 - df["EAF"]
    df = df[["Chr", "BP", "SNP", "A1", "A2", "EAF", "Beta", "Se", "P"]]
    return df

# helper function to lift over
def _lift_over_basic( df, lo_dict):
    cols = list(df.columns)
    temp = []
    lo = lo_dict["lo"]
    input_version = lo_dict["input_version"]
    output_version = lo_dict["output_version"]
    for row in df.itertuples():
        chrom = "chr" + str(row.Chr)
        pos =row.BP
        modified = lo.convert_coordinate(chrom, pos)
        # print(i)
        if modified:
            new_chrom = modified[0][0][3:]
            new_pos = modified[0][1]
            temp.append([chrom[3:], pos, new_chrom, new_pos])
    new_chr_name = output_version + '_' + 'chr'
    new_pos_name = output_version + '_' + 'pos'
    temp_df = pd.DataFrame(temp, columns=["Chr", "BP", new_chr_name, new_pos_name])
    dtype= {new_chr_name: str, new_pos_name:'Int64'}
    temp_df.drop_duplicates(subset = ['Chr','BP'], inplace=True)
    temp_df = temp_df.astype(dtype)
    return temp_df


# helper function for liftover
def _lift_over_merge( df, reference_table):
    result = pd.merge(reference_table, df, on=["Chr", "BP"], how="right")
    return result
    



# ---------------------------------------------------------------------------------------------
# functions to be implemented

def select( col, value):
    pass
# these can be done with df.query(), no need to write new function.
# can be put into doc directly

def insert():
    pass

def delete():
    pass


def create_tbi_index( df):
    pass

def query_db():
    pass

# ---------------------------------------------------------------------------------------------




    
        
        


if __name__ == "__main__":
    # -------------------------------------------------------
    """
        STEPS:
            1. read the data to be processed and formatted in uniform form 
            2. filter only the bi-allelic case and leave out other cases
            3. deduplicate data to make sure unique key (chr+bp) exist. Remove rows that contains duplicate keys.
            4. sort
            5. query required info from dbSnp153
            6. lift over to the correct genome build
            7. process the data: add rsid/ align effect allele with reference/ lift over/ flip strand
    """
    
    # -------------------------------------------------------
    # TEST CALLS ON Stroke2018NG
    # setting parameters: examples
    # input_path = "29531354-GCST006910-EFO_1001976-build37.f.tsv.gz"
    # input_path = "29531354-GCST006910-EFO_1001976.h.tsv.gz"
    # output_path = "result"
    # input_format = "hg19"
    # output_format = "hg38"

    # create class instance
    # converter = DataConverter(input_path, output_path, input_format,output_format)
    # converter = DataConverter()
    # df = converter.read_data()

    # test call for read_data()
    # df = converter.read_data(input_path, "chromosome","base_pair_location", "variant_id" ,"effect_allele", "other_allele", "effect_allele_frequency", "beta", "standard_error", "p_value")

    # df = converter.read_data(input_path, "chromosome","base_pair_location", "hm_variant_id" ,"hm_effect_allele", "hm_other_allele", "hm_effect_allele_frequency", "hm_beta", "standard_error", "p_value")
    
    # print(new_df)
    # print(df)
    # test call for filter_by_allelic()
    # bi_allelic = converter.filter_bi_allelic(df[:100000])
    # new_bi_allelic = converter.new_bi_allelic(df[:100000])
    # print(bi_allelic)
    # print(new_bi_allelic)

    # dedup_bi_allelic = converter.deduplicate(bi_allelic)

    # ut = Utility()
    # E = time.time()
    # data = ut.query_data(bi_allelic, "dbSnp153.bb")
    # F = time.time()
    # print(F-E)

    # A = time.time()
    # print(converter.flip_strand(dedup_bi_allelic, data))
    # B = time.time()
    # print(B-A)

    # print(dedup_bi_allelic)

    # reference_path = "finngen_R4_AB1_ARTHROPOD.gz"
    # reference_df = converter.read_data(reference_path, "#chrom","pos", "rsids" ,"alt", "ref", "maf", "beta", "sebeta", "pval")
    # reference_subset = reference_df.query('Chr == "10"').reset_index(drop=True)
    
    # A = time.time()
    # reference_bi_allelic = converter.filter_bi_allelic(reference_subset)
    # B = time.time()
    # print(reference_bi_allelic)
    # print(B-A)


    # print(converter.filter_bi_allelic(reference_subset, rest=True))
    # dedup_finngen = converter.deduplicate(reference_bi_allelic)
    # print(dedup_finngen)
    # other = converter.filter_bi_allelic(df, rest=True)
    # print(bi_allelic)
    # # print(other)

    # # test call for create_lo()
    # lo_dict = converter.create_lo(input_format, output_format)
    # # print(lo_dict)

    # # test call for lift_over()
    # lo_result = converter.lift_over(bi_allelic, lo_dict, "Chr", "BP", keep_original_version=True)
    # lift_over_result = converter.lift_over(bi_allelic, lo_dict, "Chr", "BP", keep_original_version=False)
    # print(lift_over_result)

    # test call for query_data()
    
    # save the query from dbSnp153 as python dictionary
    # ut = Utility()
    # dbSnp153 = ut.query_data(bi_allelic, "dbSnp153.bb")
    # converter.save_obj(data, "dbSnp153")
    # ut = Utility()
    # dbSnp153 = ut.load_obj("dbSnp153")
    # print(dbSnp153)
    # key = ("10", 42272967)
    # if key in data:
    #     print(data[key])
    # else:
    #     print("key not found")
    # print(data)
    
    # test ccall for flip_strand()
    # flipped = converter.flip_strand(bi_allelic, dbSnp153, keep_all=True)
    # print(flipped)
    # # check cases
    # # print(flipped.query('new_A1 == "4"'))
    # print(flipped.query('new_A1 == "3"'))
    # print(flipped.query('new_A1 == "2"'))
    # print(flipped.query('new_A1 == "1"'))
    # flipped_not_keep = converter.flip_strand(bi_allelic, dbSnp153)
    # print(flipped_not_keep)
    
    
    # # test call for add_rsid()
    # print(converter.add_rsid(bi_allelic, dbSnp153))

    # # test call for save_data()
    # res = converter.add_rsid(df, data)
    # converter.save_data(input_path, output_path, res, "add_rsid")

    # # test call for swap effect allele
    # print(df)
    # print(converter.swap_effect_allele(df))

    # reference_path = "finngen_R4_AB1_ARTHROPOD.gz"
    # reference_df = converter.read_data(reference_path, "#chrom","pos", "rsids" ,"alt", "ref", "maf", "beta", "sebeta", "pval")
    # reference_subset = reference_df.query('Chr == "10"').reset_index(drop=True)
    
    # reference_bi_allelic = converter.filter_bi_allelic(reference_subset)
    # dedup_reference_bi_allelic = converter.deduplicate(reference_bi_allelic)
    # print(bi_allelic)
    # print(reference_bi_allelic)
    # print(dedup_bi_allelic)
    # print(dedup_reference_bi_allelic)
    # # print(lift_over_result)
    # print("aligned result")
    # aligned = converter.align_allele_effect_size(dedup_bi_allelic, dedup_reference_bi_allelic)
    # print(aligned)
    # error_rows = converter.align_effect_allele(dedup_reference_bi_allelic,dedup_bi_allelic, check_error_rows=True)
    # aligned = converter.align_effect_allele(dedup_reference_bi_allelic, dedup_bi_allelic)
    # print(error_rows)
    # print(aligned)


    # print("check")
    # print(lo_result[lo_result['hg38_pos'] == 279194])
    # print(bi_allelic[bi_allelic['BP'] == 279194])
    # print(reference_bi_allelic[reference_bi_allelic['BP'] == 279194])

    # print(0)
    # print(dbSnp153[("10", 279194)])
    # print(1)
    # print(dbSnp153[("10", 428660)])
    # print(2)
    # print(dbSnp153[("10", 498538)])
    # print(3)
    # print(dbSnp153[("10", 606191)])
    # print(4)
    # print(dbSnp153[("10", 761226)])
    # print(561)
    # print(dbSnp153[("10", 131345386)])
    # print(562)
    # print(dbSnp153[("10", 131368196)])
    # print(563)
    # print(dbSnp153[("10", 131379938)])
    # print(564)
    # print(dbSnp153[("10", 131433359)])
    # print(565)
    # print(dbSnp153[("10", 131734323)])
    # print(bi_allelic)
    # print(bi_allelic[bi_allelic['SNP'] == "rs11814979"])
    # print(lift_over_result[lift_over_result['BP'] == 28984072])


   



    # -------------------------------------------------------
    # TEST CALLS ON Finngen


    # setting parameters: examples
    # input_path = "finngen_R4_AB1_ARTHROPOD.gz"
    # output_path = "result"
    # input_format = "hg38"
    # output_format = "hg19"

    # # create class instance
    # converter = DataConverter(input_path, output_path, input_format,output_format)
    # # df = converter.read_data()

    # # test call for read_data()
    # df = converter.read_data(input_path, '\t', "#chrom","pos", "rsids" ,"alt", "ref", "maf", "beta", "sebeta", "pval")
    # print(df)

    # # test call for filter_by_allelic()
    # bi_allelic = converter.filter_bi_allelic(df)
    # other = converter.filter_bi_allelic(df, rest=True)
    # print(bi_allelic)
    # print(other)

    # # test call for query_data()
    # data = converter.query_data(df)
    
    # # test ccall for flip_strand()
    # print(converter.flip_strand(df, data))
    
    # # test call for add_rsid()
    # print(converter.add_rsid(df, data))

    # # test call for save_data()
    # res = converter.add_rsid(df, data)
    # converter.save_data(res, "add_rsid")

    # # test call for swap effect allele
    # print(df)
    # print(converter.swap_effect_allele(df))



    # 问题一：会否存在统一各data一部分forward strand一部分backward trand
    # 问题二：会否存在同一个data
    # 问题三：lift over的时候，知识chr + pos flip to new version， should I also change other such as a1 and a2