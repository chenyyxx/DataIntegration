# user input:
# 1. a list of parameters including the column names to be matched
# 2. a list of actions required
# 

from dataintegrator import DataIntegrator as di
import sys
import json
import time

def main():
    # print ('Number of arguments:', len(sys.argv), 'arguments.')
    # print ('Argument List:', str(sys.argv))
    
    # read command line arguments (read in as JSON file)
    with open(sys.argv[1], "r") as read_file:
        data = json.load(read_file)

    # read data
    print("start reading data")
    A = time.time()
    df = di.read_data(data["input_path"], data["Chr_col_name"],data["BP_col_name"], data["SNP_col_name"] ,data["A1_col_name"], data["A2_col_name"], data["EAF_col_name"], data["Beta_col_name"], data["Se_col_name"], data["P_col_name"])
    B = time.time()

    # liftover to desired build
    print("start liftover")
    lo_dict = di.create_lo(data["input_format"], data["output_format"])
    result = di.lift_over(df, lo_dict)
    C = time.time()

    # basic format: filter bi-allelic, deduplicate, and sort
    print("start post processing")
    result_bi = di.filter_bi_allelic(result)
    dedup_result = di.deduplicate(result_bi)
    sorted_result = di.sort_by_chr_bp(dedup_result)
    D = time.time()

    # save result
    di.save_data(data["output_path"], sorted_result, data["output_name"])

    print(sorted_result)
    print("Time used (total):" , D-A)
    print("Time used (read data):" , B-A)
    print("Time used (liftover):" , C-B)
    print("Time used (post processing):" , D-C)



if __name__ == "__main__":
    main()