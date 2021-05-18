# libraries
import numpy as np
import pandas as pd
from pyliftover import LiftOver
import io
import os
import mysql.connector
import pyBigWig
import pickle
from Utility import Utility

# class
class DataConverter:

    # Class Constructor

    # def __init__(self, input_path, output_path, input_version, output_version):
    #     self.input_path = input_path
    #     self.output_path = output_path
    #     self.input_version = input_version
    #     self.output_version = output_version
    #     # self.lo = LiftOver(input_version, output_version) # TODO: do not create instance of lift over here

    # function to filter only bi-allelic case
    def filter_bi_allelic(self, df, rest=False):
        list_df = list(df)
        result = []
        other_case = []
        columns = df.columns
        for i in range(df.shape[0]):  
            # TODO: instead of using range(), maybe try using itertuples() / zip() / zip + to_dict('list') to increase time efficiency
            # These could be done for other functions that contain for loops as well.
            print(i)
            A1 = df.iloc[i]["A1"]
            A2 = df.iloc[i]["A2"]
            row = list(df.iloc[i])
            if len(A1) == 1 and len(A2) == 1: # case where current snp is a bi-allelic case
                result.append(row)
            else: # otherwise
                other_case.append(row)
        bi = pd.DataFrame(result, columns=columns)
        dtype = dict(Chr="string", BP='Int64', SNP="string", A1="string", A2="string", EAF=float, Beta=float, Se=float, P=float)
        bi = bi.astype(dtype)
        other = pd.DataFrame(other_case, columns=columns)
        dtype = dict(Chr="string", BP='Int64', SNP="string", A1="string", A2="string", EAF=float, Beta=float, Se=float, P=float)
        other = other.astype(dtype)
        if not rest:
            return bi
        else:
            return other

    # function to read data for process

    def read_data(self, input_path, Chr_col_name, BP_col_name, SNP_col_name, A1_col_name, A2_col_name, EAF_col_name, Beta_col_name, Se_col_name, P_col_name, separate_by="\t"):
        """Example function with types documented in the docstring.

        `PEP 484`_ type annotations are supported. If attribute, parameter, and
        return types are annotated according to `PEP 484`_, they do not need to be
        included in the docstring:

        Args:
            param1 (int): The first parameter.
            param2 (str): The second parameter.

        Returns:
            bool: The return value. True for success, False otherwise.

        .. _PEP 484:
            https://www.python.org/dev/peps/pep-0484/

        """
        raw_df = pd.read_csv(input_path, compression='gzip', header=0, sep=separate_by,quotechar='"')
        print(raw_df)
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
        return res

    # function to output data and save as gzip

    def save_data(self, input_path, output_path, df, prefix):
        name = input_path.split(".")[0]
        print(name)

        df_out = output_path + "/" + prefix + "_" + name +".gz"
        df.to_csv(df_out, compression='gzip')

    # helper function for liftover

    def create_lo(self, input_version, output_version):
        lo = LiftOver(input_version, output_version)
        return {"input_version": input_version, "output_version":output_version, "lo":lo}

    def __lift_over_basic(self, df, lo_dict ,chr_col_name, pos_col_name):
        cols = list(df.columns)
        temp = []
        lo = lo_dict["lo"]
        input_version = lo_dict["input_version"]
        output_version = lo_dict["output_version"]
        for i in range(df.shape[0]):
            chrom = "chr" + str(df.iloc[i][chr_col_name])
            pos =df.iloc[i][pos_col_name]
            modified = lo.convert_coordinate(chrom, pos)
            # print(i)
            if modified:
                new_chrom = modified[0][0][3:]
                new_pos = modified[0][1]
                temp.append([chrom[3:], pos, new_chrom, new_pos])
        new_chr_name = output_version + '_' + 'chr'
        new_pos_name = output_version + '_' + 'pos'
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
    # the current lift over tool only get the chr + bp in new build version but leave the SNP and A1 A2 unchanged. 
    # Should I also change the SNP and A1 A2 to match the new build version?
    def lift_over(self, df, lo_dict, chr_col_name, pos_col_name, keep_unconvertible=False):
        reference_table = self.__lift_over_basic(df, lo_dict ,chr_col_name, pos_col_name)
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

    def flip_strand(self, df, data, keep_all=False):
        flipped_A1 = []
        flipped_A2 =  []
        comment = []
        for i in range(df.shape[0]):
            chrom = df.iloc[i]["Chr"]
            pos = df.iloc[i]["BP"]
            A1 = df.iloc[i]["A1"].upper()
            A2 = df.iloc[i]["A2"].upper()

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
                        new_a1 = self.flip(A1)
                        new_a2 = self.flip(A2)
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



    # functions to align effect allele and effect size between two datasets

    ## assuming data set has no problem with consistent effect allele
    ## assuming both data set belongs to the same genome version

    def align_allele_effect_size(self, reference_data, process_data):
        reference = reference_data[["Chr", "BP", "A1", "Beta"]].rename({"A1":"reference_A1", "Beta":"reference_Beta"}, axis="columns")
        process = process_data[["Chr", "BP", "A1", "Beta"]].rename({"A1":"process_A1", "Beta":"process_Beta"}, axis="columns")
        merge_table = pd.merge(process, reference, on=["Chr, BP"], how="inner")
        first_ref_A1 = merge_table.iloc[0]["reference_A1"]
        first_proc_A1 = merge_table.iloc[0]["process_A1"]
        if first_ref_A1 == first_proc_A1: # check the rest to make sure all equal
            for i in range(1, merge_table.shape[0]):
                print(i)
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
                print(i)
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
        df["Beta"] = -1 * df["Beta"]
        df["EAF"] = 1 - df["EAF"]
        df = df[["Chr", "BP", "SNP", "A1", "A2", "EAF", "Beta", "Se", "P"]]
        return df


        
    # ---------------------------------------------------------------------------------------------
    # functions to be implemented

    def select(self, col, value):
        pass
    # these can be done with df.query(), no need to write new function.
    # can be put into doc directly

    def sort_by_Chr(self, df):
        def mixs(num):
            try:
                ele = int(num)
                return (0, ele, '')
            except ValueError:
                return (1, num, '')
        df.sort_values(by=["Chr", "BP"], key = mixs)
        return df


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
    # -------------------------------------------------------
    """
        STEPS:
            1. read the data to be processed and formatted in uniform form 
            2. filter only the bi-allelic case and leave out
            3. process the data: add rsid/ align effect allele with reference/ lift over/ flip strand
    """
    
    # -------------------------------------------------------
    # TEST CALLS ON Stroke2018NG
    # setting parameters: examples
    # input_path = "29531354-GCST006910-EFO_1001976-build37.f.tsv.gz"
    input_path = "29531354-GCST006910-EFO_1001976.h.tsv.gz"
    output_path = "result"
    input_format = "hg19"
    output_format = "hg38"

    # create class instance
    # converter = DataConverter(input_path, output_path, input_format,output_format)
    converter = DataConverter()
    # df = converter.read_data()

    # test call for read_data()
    df = converter.read_data(input_path, "chromosome","base_pair_location", "variant_id" ,"effect_allele", "other_allele", "effect_allele_frequency", "beta", "standard_error", "p_value")
    # print(df)

    # test call for filter_by_allelic()
    bi_allelic = converter.filter_bi_allelic(df)
    # other = converter.filter_bi_allelic(df, rest=True)
    print(bi_allelic)
    # print(other)

    # test call for create_lo()
    # lo_dict = converter.create_lo(input_format, output_format)
    # print(lo_dict)

    # test call for lift_over()
    # lift_over_result = converter.lift_over(bi_allelic, lo_dict, "Chr", "BP")
    # print(lift_over_result)

    # test call for query_data()
    
    # save the query from dbSnp153 as python dictionary

    # data = converter.query_data(bi_allelic)
    # converter.save_obj(data, "dbSnp153")
    ut = Utility()
    dbSnp153 = ut.load_obj("dbSnp153")
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

    reference_path = "finngen_R4_AB1_ARTHROPOD.gz"
    reference_df = converter.read_data(reference_path, "#chrom","pos", "rsids" ,"alt", "ref", "maf", "beta", "sebeta", "pval")

    referecence_bi_allelic = converter.filter_bi_allelic(reference_df)
    # print(bi_allelic)
    print(referecence_bi_allelic)
    print("aligned result")
    aligned = converter.align_allele_effect_size(referecence_bi_allelic, bi_allelic)
    print(aligned)


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