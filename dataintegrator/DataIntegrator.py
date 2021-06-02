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
        pandas.DataFrame: return formatted data in the form of pandas DataFrame

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
        df (pandas.DataFrame): The data frame to be filtered.
        rest (boolean): value indicating wether or not to keep (mark only) the non-bi-allelic cases. Default to False.

    Returns:
        pandas.DataFrame: return filtered data in the form of pandas DataFrame.

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
        df (pandas.DataFrame): The data frame to be deduplicated.

    Returns:
        pandas.DataFrame: return filtered data in the form of pandas DataFrame.

"""

def deduplicate(df):
    result = df.drop_duplicates(subset=['Chr', 'BP'], keep=False)
    return result


"""Function to sort the data based on Chr and BP

    Args:
        df (pandas.DataFrame): the data to be sorted
    Returns:
        pandas.DataFrame: return the sorted data
"""
def sort_by_chr_bp(df):
    def mixs(v):
        try:
            return int(v)
        except ValueError:
            return v
    df = df.assign(chr_numeric = lambda x: x['Chr'].apply(lambda y: 23 if y=="X" else(24 if y=="Y" else int(y))))
    result = df.sort_values(by=["chr_numeric", "BP"]).drop(['chr_numeric'], axis=1).reset_index(drop=True)
    return result




"""Function to query required data from dbSnp153

    Args:
        df (pandas.DataFrame): the data we want more info
        link (str): path or link of the '.bb' file of dbSnp153 
    Returns:
        pandas.DataFrame: return complete information from dbSnp153 as a python dictionary
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
        df (pandas.DataFrame): the data to be lifted over
        lo_dict (python dictionary): the lift over dictionary return from the create_lo function
        keep_all (boolean): if true, the function will keep and mark the rows that are not convertible. Default to False.
        inplace (boolean): if true, the function will keep the Chr + BP of original genome build. Default to False.
    Returns:
        pandas.DataFrame: return the data being lifted over to the desired genome build
"""
def lift_over(df, lo_dict, keep_all=False, inplace= False, comment=False):
    reference_table = _lift_over_basic(df, lo_dict)
    result = _lift_over_merge(df, reference_table)
    if not keep_all:
        new_chr_name = reference_table.columns[2]
        new_pos_name = reference_table.columns[3]
        result = result.dropna(subset=[new_chr_name]).reset_index(drop=True)
    if not inplace:
        new_chr_col_name = lo_dict['output_version']+"_chr"
        new_pos_col_name = lo_dict['output_version']+"_pos"
        result = result[[new_chr_col_name, new_pos_col_name, "SNP", "A1", "A2", "EAF", "Beta", "Se", "P"]].rename({new_chr_col_name:"Chr", new_pos_col_name:"BP"}, axis="columns")
    return result





"""Function to query and add rs ID for rows missing rsIDs.

    Args:
        df (pandas.DataFrame): the data to be added rs_ids
        data (python dictionary): the dictionary containing required info from dbSnp153
    Returns:
        pandas.DataFrame: return the data being added rs_ids.
"""
def add_rsid(df, data, keep_all=False, inplace=False, show_comment=False, show_errors=False):
    added_rsid = []
    comment = []
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






"""Function to flip the input data to forward strand

    Args:
        df (pandas.DataFrame): the data to be flipped to forward strand
        data (python dictionary): the dictionary containing required info from dbSnp153
        keep_unconvertible (boolean): if true, the function will keep and mark the rows that are not flipped. Default to False.

    Returns:
        pandas.DataFrame: return the data being flipped to forward strand
"""
def flip_strand( df, data, keep_all=False, inplace = False, show_comment=False, show_errors=False):
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
                    comment.append("same")
                else: # mark: what is this case? => original data T/C, dbsnp153 C/A: 10  94958283  rs111998500
                    flipped_A1.append(data_a1[0])
                    flipped_A2.append(data_a2[0])
                    # print(cur_set)
                    # print(data_a1)
                    # print(data_a2)
                    # print(A1)
                    # print(A2)
                    comment.append("different")
            else: # tri-alleic snps in dbSnp153 -> mark
                flipped_A1.append(pd.NA)
                flipped_A2.append(pd.NA)
                comment.append("dbSnp153: Indel")
        else: # key not found
            flipped_A1.append(pd.NA)
            flipped_A2.append(pd.NA)
            comment.append("Key not found") 

    result = df.assign(new_A1 = flipped_A1)
    result = result.assign(new_A2 = flipped_A2)


    if show_errors:
        if inplace or show_comment or keep_all:
            print("The `show_errors` flag cannot be used with the `inplace`, `comment`, and `keep_all` flags together.")
            return
        result = result.assign(comment=comment)
        mask = (result["comment"] != "flipped") & (result["comment"] != "same")
        result = result[mask]
        return result

    # print(result)
    if inplace:
        result = result[["Chr", "BP" , "new_A1", "new_A2", "EAF", "Beta", "Se", "P"]].rename({"new_A1": "A1", "new_A2":"A2"},axis="columns")
    if show_comment:
        result = result.assign(comment=comment)
    if not keep_all:
        if inplace:
            result  = result.dropna(subset=["A1", "A2"]).reset_index(drop=True)
        else:
            result  = result.dropna(subset=["new_A1", "new_A2"]).reset_index(drop=True)

    return result




"""Function to align effect allele
    this function will align the effect allele of input data based on a reference data

    Args:
        reference (pandas.DataFrame): the reference table
        df (pandas.DataFrame): the data to be aligned
        check_error_rows (boolean): if true, the function will output the rows that cannot be aligned. Default to False.
    Returns:
        pandas.DataFrame: return the data with its effect allele being aligned with the reference table.
"""

def align_effect_allele( reference, df, show_errors=False):
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
    sorted_result = sort_by_chr_bp(result)
    if show_errors:
        return pd.merge(error, merge_table[error_mask], on=["Chr", "BP"], how="inner")
    print(str(nochange.shape[0]) + " rows were left unchanged (already aligned)")
    print(str(aligned.shape[0]) + " rows were aligned successfully")
    print(str(error.shape[0]) + " rows failed to align, dropped from result! Set the check_error_rows flag to True to view them.")
    return sorted_result
        



        
"""Function to save the processed data in gz or csv

    Args:
        output_path (str): the path you want the data to be saved.
        df (pandas.DataFrame): the processed data to be saved.
        name (str): the output name of the data.
        save_format (str): the saving format. Choose between 'gzip' or 'csv'. Default to gz.
    Returns:
        pandas.DataFrame: return filtered data in the form of pandas DataFrame

"""

def save_data(output_path, df, name, save_format="gzip"):
    # TODO: add support for other compression format/ txt
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


