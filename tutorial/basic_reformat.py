from dataintegrator import DataIntegrator as di
import time

def main():
    """
        1. read the data to be processed and formatted in uniform form
        2. filter only the bi-allelic case and leave out other cases
        3. deduplicate data to make sure on unique key (chr+bp) exist. Remove the rows that contain duplicate keys
        4. sort
        5. lift over to the correct genome build (hg19)
    """

    # For example: performing reformatting on finngen data
    input_path = "data/finngen_R4_AB1_ARTHROPOD.gz"
    output_path = "result"
    input_format = "hg38"
    output_format = "hg19"

    # 1. read data
    # ----------------------------------------------------------------------------------------------------------------------------------------------------
    print("start reading data")
    df = di.read_data(input_path, "#chrom","pos", "rsids" ,"alt", "ref", "maf", "beta", "sebeta", "pval")
    print("data successfully read")

    # TODO:
    # add a function to only keep chr 1-23/XY
    # in case sth like 1_gl000192_random break the program.
    # can do this by doing filter bi-allelic again
    # optimize the read_data() input list: read a txt/JSON file instead of all the parameters
    # 

    # 2. filter biallelic
    # ----------------------------------------------------------------------------------------------------------------------------------------------------
    print("start filtering bi-allelic cases")
    bi_allelic = di.filter_bi_allelic(df)

    # 3. drop duplicates
    # ----------------------------------------------------------------------------------------------------------------------------------------------------
    print("start deduplicating")
    dedup_df = di.deduplicate(bi_allelic)

    # 4. getting required information for liftover
    # ----------------------------------------------------------------------------------------------------------------------------------------------------
    print("start getting required info for liftover")
    C = time.time()
    lo_dict = di.create_lo(input_format, output_format)
    E = time.time()

    # 5. data processing (e.g.): lift over
    # ----------------------------------------------------------------------------------------------------------------------------------------------------
    print("start processing data: lift over")
    result = di.lift_over(dedup_df, lo_dict)
    
    D = time.time()
    

    # 6. sort
    # ----------------------------------------------------------------------------------------------------------------------------------------------------
    print("start post processing")
    

    result_bi = di.filter_bi_allelic(result)
    print(result_bi)

    sorted_result = di.sort_by_chr_bp(result_bi)
    print(sorted_result)

    dedup_result = di.deduplicate(sorted_result)
    print(dedup_result)


    print("Time used:" , D-C)
    print("Time used (query):" , E-C)
    print("Time used (lift over):" , D-E)

    # 8. save data
    # ----------------------------------------------------------------------------------------------------------------------------------------------------
    di.save_data(output_path, result, "basic_reformat")


if __name__ == "__main__":
    main()