import sys
sys.path.append("..")

from data_converter.DataConverter import DataConverter
from data_converter.Utility import Utility
import time

def main():
    """
        STEPS:
            1. read the data to be processed and formatted in uniform form 
            2. filter only the bi-allelic case and leave out other cases
            3. deduplicate data to make sure unique key (chr+bp) exist. Remove rows that contains duplicate keys.
            4. sort
            5. lift over to the correct genome build
            6. query required info from dbSnp153
            7. process the data: add rsid/ align effect allele with reference/ lift over/ flip strand
            8. save data
    """
    # create DataConverter() instance
    # ----------------------------------------------------------------------------------------------------------------------------------------------------
    dc = DataConverter()

    # setting up global variables such paths and build version
    input_path = "data/29531354-GCST006910-EFO_1001976.h.tsv.gz"
    output_path = "result"

    print("start reading data")
    df = dc.read_data(input_path, "chromosome","base_pair_location", "variant_id" ,"effect_allele", "other_allele", "effect_allele_frequency", "beta", "standard_error", "p_value")
    print("data successfully read")

    # 2. filter biallelic
    # ----------------------------------------------------------------------------------------------------------------------------------------------------
    print("start filtering bi-allelic cases")
    bi_allelic = dc.filter_bi_allelic(df)

    # 3. drop duplicates
    # ----------------------------------------------------------------------------------------------------------------------------------------------------
    print("start deduplicating")
    dedup_df = dc.deduplicate(bi_allelic)

    # 4. sort
    # ----------------------------------------------------------------------------------------------------------------------------------------------------
    print("start sorting")
    sorted_df = dc.sort_by_Chr(dedup_df)
    # print(sorted_df)



    # 6. query dbSnp153 for required information
    # ----------------------------------------------------------------------------------------------------------------------------------------------------
    print("start getting required info from dbSnp153")
    C = time.time()
    ut = Utility()
    dbSnp153 = ut.query_data(sorted_df, "../dbSnp153.bb")
    ut.save_obj("stroke_dbSNp153")
    print("end querying")
    E = time.time()

    # 7.3 data processing (e.g.): Add rsID
    # ----------------------------------------------------------------------------------------------------------------------------------------------------
    print("start processing data: add rsID")
    flipped = dc.flip_strand(df, dbSnp153)
    D = time.time()
    print(flipped)
    print("Time used:" , D-C)
    print("Time used (query):" , E-C)
    print("Time used (process):" , D-E)

    

if __name__ == "__main__":
    main()