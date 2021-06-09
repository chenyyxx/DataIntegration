#!/bin/bash

for i in {21..22};
do  
    echo $i
    ./bigBedToBed dbSnp153.bb -chrom=chr$i dbSnp153_$i.bed
done

./bigBedToBed dbSnp153.bb -chrom=chrX dbSnp153_X.bed