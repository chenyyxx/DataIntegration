
import time
from dataintegrator import DataIntegrator as di

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

    # 1. Read Data
    # ----------------------------------------------------------------------------------------------------------------------------------------------------
    print("start reading data")
    df = di.read_data(input_path, "#chrom","pos", "rsids" ,"alt", "ref", "maf", "beta", "sebeta", "pval")
    print("data successfully read")

    # 2. Filter Biallelic
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
    di.save_obj(dbSnp153, "obj/dbSnp153_finngen")
    # dbSnp153 = di.load_obj("obj/dbSnp153_finngen.pkl") # once save you can load it
    print("end querying")
    E = time.time()


    # 7 data processing (e.g.): Add rsID
    # ----------------------------------------------------------------------------------------------------------------------------------------------------
    print("start processing data: add rsID")
    added_rsid = di.add_rsid(sorted_df, dbSnp153)
    D = time.time()
    print(added_rsid)
    print("Time used:" , D-C)
    print("Time used (query):" , E-C)
    print("Time used (process):" , D-E)

    

if __name__ == "__main__":
    main()