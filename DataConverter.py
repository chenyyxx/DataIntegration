# libraries
import numpy as np
import pandas as pd
from pyliftover import LiftOver
import io
import os

# class
class DataConverter:
    def __init__(self, input_path, output_path, lo):
        self.input_path = input_path
        self.output_path = output_path
        self.lo = lo

    def read_data(self):
        df = pd.read_csv(self.input_path, compression='gzip', header=0, sep='\t',quotechar='"')
        result = df.loc[:,["#chrom","pos", "rsids" ,"alt", "ref", "maf", "beta", "sebeta", "pval"]]
        res = result.rename({"#chrom":"Chr","pos":"BP", "rsids":"SNP", "alt":"A1", "ref":"A2", "maf":"EAF", "beta":"Beta", "sebeta":"Se", "pval":"P"},axis="columns")
        return res
        # return df

    def process_data(self, df):
        # print(df)
        cols = list(df.columns)
        modified_output = pd.DataFrame(columns=cols)
        unconvertible_output = pd.DataFrame(columns=cols)
        # num_containing_rsids_modified = 0
        # num_containing_rsids_unconvertible = 0
        for i in range(df.shape[0]):
            chrom = "chr" + str(df.iloc[i]["Chr"])
            # print(chrom)
            pos = df.iloc[i]["BP"]
            # row = df[i:i+1]
            modified = lo.convert_coordinate(chrom, pos)
            row = df.iloc[i]
            if modified:
                # if pd.isna(row["SNP"]):
                #     num_containing_rsids_modified += 1
                new_chrom = modified[0][0]
                new_pos = modified[0][1]
                row["Chr"] = new_chrom[3:]
                row["BP"] = new_pos
                modified_output = modified_output.append(row, ignore_index=True)
            else:
                # print(pd.isna(row["SNP"]))
                # if pd.isna(row["SNP"]):
                #     num_containing_rsids_unconvertible += 1
                modified_output = modified_output.append(pd.Series(), ignore_index=True)
                unconvertible_output = unconvertible_output.append(row, ignore_index=True)
        # print(num_containing_rsids)
        return [modified_output, unconvertible_output]
    
    def save_data(self, df, prefix):
        name = self.input_path.split(".")[0]

        df_out = self.output_path + "/" + prefix + "_" + name +".csv"
        df.to_csv(df_out)

        meta_file_name = self.output_path + "/" + prefix + "_"+ name +"_meta.txt"
        buf = io.StringIO()
        df.info(buf = buf)
        s= buf.getvalue()
        with open(meta_file_name, "w") as meta:
            meta.write(s)

    def lift_over(self):
        pass

    def add_rsid(self):
        pass

    def select(self):
        pass

    def insert(self):
        pass

    def delete(self):
        pass

    def flip_strand(self):
        pass

    def deduplicate(self):
        pass

    def zip_data(self):
        pass

    def align_allele_effect_size(self):
        pass



        
        
        


if __name__ == "__main__":
    input_path = "finngen_R4_AB1_ARTHROPOD.gz"
    output_path = "result"
    input_format = "hg38"
    output_format = "hg19"
    lo = LiftOver(input_format, output_format)
    converter = DataConverter(input_path, output_path, lo)
    df = converter.read_data()
    # print(df)
    result = converter.process_data(df)
    converter.save_data(result[0], "modified")
    converter.save_data(result[1], "unconvertible")


    # output_path = "result"
    # input_format = "hg38"
    # output_format = "hg19"
    # lo = LiftOver(input_format, output_format)
    # all_file = os.listdir()
    # for f in all_file:
    #     print(f)
    #     input_path = f
    #     converter = DataConverter(input_path, output_path, lo)
    #     df = converter.read_data()
    #     # print(df)
    #     result = converter.process_data(df)
    #     converter.save_data(result[0], "modified")
    #     converter.save_data(result[1], "unconvertible")