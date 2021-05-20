from DataConverter import DataConverter
from Utility import Utility

def main():
    """
        STEPS:
            1. read the data to be processed and formatted in uniform form 
            2. filter only the bi-allelic case and leave out other cases
            3. deduplicate data to make sure unique key (chr+bp) exist. Remove rows that contains duplicate keys.
            4. sort
            5. query required info from dbSnp153
            6. lift over to the correct genome build
            7. process the data: add rsid/ align effect allele with reference/ lift over/ flip strand
            8. save data
    """
    # create DataConverter() instance
    dc = DataConverter()

    # setting up global variables such paths and build version
    input_path = "29531354-GCST006910-EFO_1001976.h.tsv.gz"
    reference_path = "finngen_R4_AB1_ARTHROPOD.gz"
    output_path = "result"
    input_format = "hg19"
    output_format = "hg38"


    # read data
    print("start reading data")
    df = dc.read_data(input_path, "chromosome","base_pair_location", "variant_id" ,"effect_allele", "other_allele", "effect_allele_frequency", "beta", "standard_error", "p_value")
    # df = dc.read_data(input_path, "chromosome","base_pair_location", "variant_id" ,"hm_effect_allele", "hm_other_allele", "hm_effect_allele_frequency", "hm_beta", "standard_error", "p_value")
    # reference_df = dc.read_data(reference_path, "#chrom","pos", "rsids" ,"alt", "ref", "maf", "beta", "sebeta", "pval")
    print("data successfully read")

    # filter biallelic
    print("start filtering bi-allelic cases")
    bi_allelic = dc.filter_bi_allelic(df[:100000])
    # reference_bi_allelic = dc.filter_bi_allelic(reference_df)

    # drop duplicates
    print("start deduplicating")
    dedup_df = dc.deduplicate(bi_allelic)
    # dedup_reference = dc.deduplicate(reference_bi_allelic)

    # sort
    print("start sorting")

    # query dbSnp153 for required information
    print("start getting required info from dbSnp153")
    ut = Utility()
    dbSnp153 = ut.query_data(bi_allelic, "dbSnp153.bb")
    print("end querying")

    # data processing (e.g.): Flip Strand
    print("start processing data: flip strand")
    result=dc.flip_strand(dedup_df, dbSnp153)
    print(result)

    # save data
    dc.save_data(input_path, output_path, result, "flipped_strand")


if __name__ == "__main__":
    main()