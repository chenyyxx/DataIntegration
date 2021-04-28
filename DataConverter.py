# libraries
import numpy as np
import pandas as pd
from pyliftover import LiftOver
import io
import os
import mysql.connector
import pyBigWig

# class
class DataConverter:

    # Class Constructor

    def __init__(self, input_path, output_path, input_version, output_version):
        self.input_path = input_path
        self.output_path = output_path
        self.input_version = input_version
        self.output_verion = output_version
        self.lo = LiftOver(input_version, output_version)

    # function to read data for process

    def read_data(self, input_path, separate_by, Chr_col_name, BP_col_name, SNP_col_name, A1_col_name, A2_col_name, EAF_col_name, Beta_col_name, Se_col_name, P_col_name):
        raw_df = pd.read_csv(self.input_path, compression='gzip', header=0, sep=separate_by,quotechar='"')
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
        return res[:1000]

    # function to output data and save as gzip

    def save_data(self, df, prefix):
        name = self.input_path.split(".")[0]
        print(name)

        df_out = self.output_path + "/" + prefix + "_" + name +".gz"
        df.to_csv(df_out, compression='gzip')

    # helper function for liftover

    def __lift_over_basic(self, df, chr_col_name, pos_col_name):
        cols = list(df.columns)
        temp = []
        for i in range(df.shape[0]):
            chrom = "chr" + str(df.iloc[i][chr_col_name])
            pos =df.iloc[i][pos_col_name]
            modified = self.lo.convert_coordinate(chrom, pos)
            print(i)
            if modified:
                new_chrom = modified[0][0][3:]
                new_pos = modified[0][1]
                temp.append([chrom[3:], pos, new_chrom, new_pos])
        new_chr_name = self.output_verion + '_' + 'chr'
        new_pos_name = self.output_verion + '_' + 'pos'
        temp_df = pd.DataFrame(temp, columns=[chr_col_name, pos_col_name, new_chr_name, new_pos_name])
        dtype= {new_chr_name: str, new_pos_name:'Int64'}
        temp_df.drop_duplicates(subset = ['Chr','BP'], inplace=True)
        temp_df = temp_df.astype(dtype)
        return temp_df
    
    # helper function for liftover

    def __lift_over_merge(self, df, reference_table):
        result = pd.merge(reference_table, df, on=["Chr", "BP"], how="right")
        return result

    # function for liftover

    def lift_over(self, df, chr_col_name, pos_col_name, keep_unconvertible=False):
        reference_table = self.__lift_over_basic(df, chr_col_name, pos_col_name)
        result = self.__lift_over_merge(df, reference_table)
        if not keep_unconvertible:
            new_chr_name = reference_table.columns[2]
            new_pos_name = reference_table.columns[3]
            result = result.dropna(subset=[new_chr_name]).reset_index(drop=True)
        return result

    # function to fill in missing rsid

    def add_rsid(self, df, data):
        added_rs_id = []
        for i in range(df.shape[0]):
            chrom = df.iloc[i]["Chr"]
            pos = df.iloc[i]["BP"]
            rs_id = df.iloc[i]["SNP"]
            key = (chrom, pos)
            if key in data:
                raw_string = data[key]
                parsed_string = raw_string.split('\t')
                data_rs_id = parsed_string[0] 
                if pd.isna(rs_id):
                    added_rs_id.append(data_rs_id)
                    # print("find none")
                elif rs_id == data_rs_id:
                    added_rs_id.append(rs_id)
                    # print("same")
                else: # find different rsid in dnpsnp153, update with new
                    added_rs_id.append(data_rs_id)
                    # print("different")
            else:
                added_rs_id.append("key not found")
                # print("key not found")
        result = df.assign(added_rs_id = added_rs_id)
        # print(result)
        return result


    # helper function for query dbsnp153 from USCS genome broswer
    
    def query_data(self, df):
        bb = pyBigWig.open("http://hgdownload.soe.ucsc.edu/gbdb/hg38/snp/dbSnp153.bb")
        result = {}
        set_list = []
        for i in range(df.shape[0]):
            chrom = "chr" + str(df.iloc[i]["Chr"])
            end_pos = df.iloc[i]["BP"]
            start_pos =end_pos - 1
            dat = bb.entries(chrom, start_pos, end_pos)
            if dat != None:
                reference_start = dat[-1][0]
                reference_end = dat[-1][1]
                raw_string = dat[-1][2]
                key = (chrom[3], reference_end)
                result[key] = raw_string
        # print(result)
        return result

    # helper function for flipping strand

    def flip(self, allele):
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

    # function to flip strand

    def flip_strand(self, df, data):
        flipped_A1 = []
        flipped_A2 =  []
        for i in range(df.shape[0]):
            chrom = df.iloc[i]["Chr"]
            pos = df.iloc[i]["BP"]
            A1 = df.iloc[i]["A1"]
            A2 = df.iloc[i]["A2"]
            if len(A1) == 1 and len(A2) == 1:
                key = (chrom, pos)
                if key in data:
                    cur_set = {A1, A2}
                    raw_string = data[key]
                    parsed_string = raw_string.split("\t")
                    data_a1 = parsed_string[1].split(",")
                    if "" in data_a1:
                        data_a1.remove("")
                    data_a2 = parsed_string[3].split(",")
                    if "" in data_a2:
                        data_a2.remove("")
                    if len(data_a1) == 1 and len(data_a2) == 1:
                        cur_set.add(data_a1[0])
                        cur_set.add(data_a2[0])
                        # print(cur_set)
                        if len(cur_set) == 4: # flip
                            new_a1 = self.flip(A1)
                            new_a2 = self.flip(A2)
                            flipped_A1.append(new_a1)
                            flipped_A2.append(new_a2)
                        elif len(cur_set) == 2: # do not flip
                            # print(i)
                            flipped_A1.append(A1)
                            flipped_A2.append(A2)
                        else: # mark
                            flipped_A1.append("1")
                            flipped_A2.append("1")
                    else: # tri-alleic snps -> mark
                        flipped_A1.append("3")
                        flipped_A2.append("3")
                else: # key not found
                    flipped_A1.append("4")
                    flipped_A2.append("4")
            else: # tri-alleic snps -> mark
                flipped_A1.append("3")
                flipped_A2.append("3")
        # print(flipped_A1)
        # print(len(flipped_A1))
        # print(flipped_A2)
        # print(len(flipped_A2))
        result = df.assign(new_A1 = flipped_A1)
        result = result.assign(new_A2 = flipped_A2)
        # print(result)
        return result


    # functions to align effect allele and effect size between two datasets

    ## assuming data set has no problem with consistent effect allele
    ## assuming both data set belongs to the same genome version

    def align_allele_effect_size(self, reference_data, process_data):
        reference = reference_data[["Chr", "BP", "A1", "Beta"]].rename({"A1":"reference_A1", "Beta":"reference_Beta"}, axis="columns")
        process = process_data[["Chr", "BP", "A1", "Beta"]].rename({"A1":"process_A1", "Beta":"process_Beta"}, axis="columns")
        merge_table = pd.merge(process, reference, by=["Chr, Pos"], how="inner")
        first_ref_A1 = merge_table.iloc[0]["reference_A1"]
        first_proc_A1 = merge_table.iloc[0]["process_A1"]
        if first_ref_A1 == first_proc_A1: # check the rest to make sure all equal
            for i in range(1, merge_table.shape[0]):
                ref_A1 = merge_table.iloc[i]["reference_A1"]
                proc_A1 = merge_table.iloc[i]["process_A1"]
                if ref_A1 != proc_A1:
                    print("data effect allele is not consistent")
                    return
            # all consitent: no need to change
            print("two data sets have same effect allele, no need to change")
            return

        else: # check the rest to make sure all not equal
            for i in range(1, merge_table.shape[0]):
                ref_A1 = merge_table.iloc[i]["reference_A1"]
                proc_A1 = merge_table.iloc[i]["process_A1"]
                if ref_A1 == proc_A1:
                    print("data effect allele is not consistent")
                    return
            # all consitent: need to change size 
            result = self.swap_effect_allele(process_data)
            return result
    
    # helper function to swap effect allele and align effect size

    def swap_effect_allele(self, df):
        col_list = list(df)
        col_list[3], col_list[4] = col_list[4], col_list[3]
        df.columns = col_list
        df["Beta"] -= 1.0
        df = df[["Chr", "BP", "SNP", "A1", "A2", "EAF", "Beta", "Se", "P"]]
        return df


        
    # ---------------------------------------------------------------------------------------------
    # functions to be implemented

    def select(self):
        pass

    def insert(self):
        pass

    def delete(self):
        pass
    
    def deduplicate(self):
        pass

    def create_tbi_index(self, df):
        pass

    def query_db(self):
        pass

    # ---------------------------------------------------------------------------------------------
    
    # unsued functions 
    # helper function to connect to gene browser database
    def connect_db(self):
        try:
            db = mysql.connector.connect(
                user="genome",
                host="genome-mysql.soe.ucsc.edu",
                database="hg38",
                port="3306")
            return db
        except:
            print("connection failed, try another port")

    
    
    # ---------------------------------------------------------------------------------------------


        
        
        


if __name__ == "__main__":
    # setting parameters: examples
    input_path = "finngen_R4_AB1_ARTHROPOD.gz"
    output_path = "result"
    input_format = "hg38"
    output_format = "hg19"

    # create class instance
    converter = DataConverter(input_path, output_path, input_format,output_format)
    # df = converter.read_data()

    # test call for read_data()
    df = converter.read_data(input_path, '\t', "#chrom","pos", "rsids" ,"alt", "ref", "maf", "beta", "sebeta", "pval")
    print(df)

    # test call for query_data()
    data = converter.query_data(df)
    
    # test ccall for flip_strand()
    print(converter.flip_strand(df, data))
    
    # test call for add_rsid()
    print(converter.add_rsid(df, data))

    # test call for save_data()
    res = converter.add_rsid(df, data)
    converter.save_data(res, "add_rsid")

    # test call for swap effect allele
    print(df)
    print(converter.swap_effect_allele(df))