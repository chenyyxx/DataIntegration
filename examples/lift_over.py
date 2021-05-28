

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
    input_path = "data/finngen_R4_AB1_ARTHROPOD.gz"
    output_path = "result"
    input_format = "hg38"
    output_format = "hg19"

    # 1. Read Data
    # ----------------------------------------------------------------------------------------------------------------------------------------------------
    print("start reading data")
    df = di.read_data(input_path, "#chrom","pos", "rsids" ,"alt", "ref", "maf", "beta", "sebeta", "pval")
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
    print("start getting required info for liftover")
    C = time.time()
    lo_dict = di.create_lo(input_format, output_format)
    E = time.time()

    # 7 data processing (e.g.): lift over
    # ----------------------------------------------------------------------------------------------------------------------------------------------------
    print("start processing data: add rsID")
    result = di.lift_over(dedup_df, lo_dict)
    D = time.time()
    print(result)


    print("Time used:" , D-C)
    print("Time used (query):" , E-C)
    print("Time used (lift over):" , D-E)

    

if __name__ == "__main__":
    main()