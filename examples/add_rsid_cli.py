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
    df = di.read_formatted_data(data["input_path"])
    B = time.time()

    # liftover to desired build
    print("query data")
    dbSnp153 = di.query_data(df, data["dbSnp153_path"])
    C = time.time()

    # basic format: filter bi-allelic, deduplicate, and sort
    print("start adding rsid")
    added_rsid = di.add_rsid(df, dbSnp153,
        select_cols=data["select_cols"], 
        filter_rows=data["filter_rows"])
    D = time.time()

    # save result
    di.save_data(data["output_path"], added_rsid, data["output_name"])

    print(added_rsid)
    print("Time used (total):" , D-A)
    print("Time used (read data):" , B-A)
    print("Time used (query data):" , C-B)
    print("Time used (adding rsids):" , D-C)



if __name__ == "__main__":
    main()