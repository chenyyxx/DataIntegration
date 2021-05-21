# libraries
import numpy as np
import pandas as pd
import pyBigWig
import pickle

class Utility:
    def __init__(self):
        super().__init__()


    # link = "http://hgdownload.soe.ucsc.edu/gbdb/hg38/snp/dbSnp153.bb"

    def query_data(self, df, link):
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

    def save_obj(self, obj, name ):
        with open('obj/'+ name + '.pkl', 'wb') as f:
            pickle.dump(obj, f, pickle.HIGHEST_PROTOCOL)

    def load_obj(self, name ):
        with open('obj/' + name + '.pkl', 'rb') as f:
            return pickle.load(f)

    def read_data(self, path):
        pass

    def save_data(self, df):
        pass


if __name__ == "__main__":
    # content = []
    result = {}
    with open("dbSnp153.bed")as f:
        for line in f:
            splitted = line.strip().split()
            chrom = splitted[0]
            start = splitted[1]
            end = splitted[2]
            rest = splitted[3:]
            key = (chrom, start, end)
            result[key] = rest
            print(key)
            # content.append(line.strip().split())
            # print(line.strip().split())
    print(result)
    save_obj(result, dbSnp153)

