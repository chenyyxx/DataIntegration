import pandas as pd
import pyBigWig
import pickle
import time
from statistics import mean

bb = pyBigWig.open("data/dbSnp153.bb")

def load_obj(obj_path ):
    with open(obj_path, 'rb') as f:
        return pickle.load(f)

pkl = load_obj("obj/dbSnp153_finngen.pkl")
# print(pkl)

chrom = "19"
start = 11397500
end = 11397501


pbw_time = []
pkl_time = []
for i in range(100000):
    A = time.time()
    resa = bb.entries("chr"+chrom, start, end)
    # print(resa)
    B = time.time()
    pbw_time.append(B-A)

    resb = pkl[(chrom, end)]
    # print(resb)
    C = time.time()
    pkl_time.append(C-B)

print("pbw: ", mean(pbw_time))
print("pkl: ", mean(pkl_time))
print(mean(pbw_time)/mean(pkl_time))