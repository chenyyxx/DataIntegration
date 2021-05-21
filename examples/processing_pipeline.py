from DataConverter import DataConverter
from Utility import Utility

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
    input_path = "29531354-GCST006910-EFO_1001976.h.tsv.gz"
    reference_path = "finngen_R4_AB1_ARTHROPOD.gz"
    output_path = "result"
    input_format = "hg19"
    output_format = "hg38"


    # 1. read data
    # ----------------------------------------------------------------------------------------------------------------------------------------------------
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
    sorted_reference = dc.sort_by_Chr(dedup_reference)
    # print(sorted_df)


    # 5. lift over to the desired genome build
    # ----------------------------------------------------------------------------------------------------------------------------------------------------
    # print("start converting genome build")
    # lo_dict = dc.create_lo(input_format, output_format)
    # lift_over_result = dc.lift_over(sorted_df, lo_dict)
    # print(lift_over_result)

    # 6. query dbSnp153 for required information
    # ----------------------------------------------------------------------------------------------------------------------------------------------------
    print("start getting required info from dbSnp153")
    ut = Utility()
    dbSnp153 = ut.query_data(sorted_reference, "dbSnp153.bb")
    print("end querying")

    

    # 7.1 data processing (e.g.): Flip Strand
    # ----------------------------------------------------------------------------------------------------------------------------------------------------
    # print("start processing data: flip strand")
    # result=dc.flip_strand(dedup_df, dbSnp153)
    # print(result)

    # 7.2 data processing (e.g.): Align Effect Allele
    # ----------------------------------------------------------------------------------------------------------------------------------------------------
    # print("start processing data: align effect allele")
    # result = dc.align_effect_allele(sorted_reference, sorted_df)
    # errors = dc.align_effect_allele(sorted_reference, sorted_df, check_error_rows=True)
    # print(result)
    # print(errors)

    # 7.3 data processing (e.g.): Add rsID
    # ----------------------------------------------------------------------------------------------------------------------------------------------------
    print("start processing data: add rsID")
    added_rsid = dc.add_rsid(sorted_reference, dbSnp153)
    print(added_rsid)


    # 8. save data
    # ----------------------------------------------------------------------------------------------------------------------------------------------------
    # dc.save_data(input_path, output_path, result, "flipped_strand")


if __name__ == "__main__":
    main()