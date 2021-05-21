from DataConverter import DataConverter
from Utility import Utility
import time


def main():
    # create DataConverter() instance
    # ----------------------------------------------------------------------------------------------------------------------------------------------------
    dc = DataConverter()

    # setting up global variables such paths and build version
    # input_path = "29531354-GCST006910-EFO_1001976.h.tsv.gz"
    reference_path = "finngen_R4_AB1_ARTHROPOD.gz"

    print("start reading data")
    # df = dc.read_data(input_path, "chromosome","base_pair_location", "variant_id" ,"effect_allele", "other_allele", "effect_allele_frequency", "beta", "standard_error", "p_value")
    # df = dc.read_data(input_path, "chromosome","base_pair_location", "variant_id" ,"hm_effect_allele", "hm_other_allele", "hm_effect_allele_frequency", "hm_beta", "standard_error", "p_value")
    reference_df = dc.read_data(reference_path, "#chrom","pos", "rsids" ,"alt", "ref", "maf", "beta", "sebeta", "pval")
    print("data successfully read")

    # 2. filter biallelic
    # ----------------------------------------------------------------------------------------------------------------------------------------------------
    print("start filtering bi-allelic cases")
    # bi_allelic = dc.filter_bi_allelic(df)
    # print(bi_allelic)
    # print(dc.filter_bi_allelic(df, rest=True))
    reference_bi_allelic = dc.filter_bi_allelic(reference_df)

    # 3. drop duplicates
    # ----------------------------------------------------------------------------------------------------------------------------------------------------
    print("start deduplicating")
    # dedup_df = dc.deduplicate(bi_allelic)
    dedup_reference = dc.deduplicate(reference_bi_allelic)

    # 4. sort
    # ----------------------------------------------------------------------------------------------------------------------------------------------------
    print("start sorting")
    # sorted_df = dc.sort_by_Chr(dedup_df)
    sorted_reference = dc.sort_by_Chr(dedup_reference)[:1000000]
    # print(sorted_df)




    # 7.3 data processing (e.g.): Add rsID
    # ----------------------------------------------------------------------------------------------------------------------------------------------------
    print("start processing data: add rsID")
    A = time.time()
    new_added_rsid = dc.new_add_rsid(sorted_reference)
    B = time.time()
    print(new_added_rsid)
    print("Time used (new): ", B-A)



    # 6. query dbSnp153 for required information
    # ----------------------------------------------------------------------------------------------------------------------------------------------------
    print("start getting required info from dbSnp153")
    C = time.time()
    ut = Utility()
    dbSnp153 = ut.query_data(sorted_reference, "dbSnp153.bb")
    print("end querying")
    E = time.time()

    # 7.3 data processing (e.g.): Add rsID
    # ----------------------------------------------------------------------------------------------------------------------------------------------------
    print("start processing data: add rsID")
    added_rsid = dc.add_rsid(sorted_reference, dbSnp153)
    D = time.time()
    print(added_rsid)
    print("Time used (old):" , D-C)
    print("Time used (query):" , E-C)
    print("Time used (process):" , D-E)

    

if __name__ == "__main__":
    main()