

from dataintegrator import DataIntegrator as di
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
    # setting up global variables such paths and build version
    # ----------------------------------------------------------------------------------------------------------------------------------------------------
    input_path = "data/29531354-GCST006910-EFO_1001976.h.tsv.gz"
    output_path = "result"

    # 1. Read Data
    # ----------------------------------------------------------------------------------------------------------------------------------------------------
    print("start reading data")
    df = di.read_data(input_path, "chromosome","base_pair_location", "variant_id" ,"effect_allele", "other_allele", "effect_allele_frequency", "beta", "standard_error", "p_value")
    print(df.shape)
    print("data successfully read")

    # 2. filter biallelic
    # ----------------------------------------------------------------------------------------------------------------------------------------------------
    print("start filtering bi-allelic cases")
    bi_allelic = di.filter_bi_allelic(df)

    # 3. drop duplicates
    # ----------------------------------------------------------------------------------------------------------------------------------------------------
    print("start deduplicating")
    dedup_df = di.deduplicate(bi_allelic)

    # 4. sort
    # ----------------------------------------------------------------------------------------------------------------------------------------------------
    print("start sorting")
    sorted_df = di.sort_by_chr_bp(dedup_df)
    # print(sorted_df)



    # 6. query dbSnp153 for required information
    # ----------------------------------------------------------------------------------------------------------------------------------------------------
    print("start getting required info from dbSnp153")
    C = time.time()
    dbSnp153 = di.query_data(sorted_df, "data/dbSnp153.bb")
    # di.save_obj(dbSnp153, "obj/dbSnp153_stroke")
    # dbSnp153 = di.load_obj("obj/dbSnp153_stroke.pkl") # once save you can load it
    print("end querying")
    E = time.time()

    # 7 data processing (e.g.): flip strand
    # ----------------------------------------------------------------------------------------------------------------------------------------------------
    print("start processing data: flip_strand")
    # flipped = di.flip_strand(sorted_df, dbSnp153)
    D = time.time()
    # print(flipped)

    # flipped = di.flip_strand(sorted_df, dbSnp153, keep_all=True)
    # print("keep_all true")
    # print(flipped)

    # flipped = di.flip_strand(sorted_df, dbSnp153, show_comment=True)
    # print("comment true")
    # print(flipped)

    # flipped = di.flip_strand(sorted_df, dbSnp153, inplace=True,)
    # print("inplace true")
    # print(flipped)

    # flipped = di.flip_strand(sorted_df, dbSnp153, inplace=True, show_comment=True)
    # print("comment inplace true")
    # print(flipped)
    flipped = di.flip_strand(sorted_df, dbSnp153, show_errors=True)
    print("show errors true")
    print(flipped)


    print("Time used:" , D-C)
    print("Time used (query):" , E-C)
    print("Time used (process):" , D-E)

    

if __name__ == "__main__":
    main()