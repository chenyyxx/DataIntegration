# DataIntegration
Data Integration tool set for genetic data


# Main functions
1. liftover genome build
Add chr and pos columns according to hg18, hg19 or GRCh38 in the file and then output the other version chr and pos into the new gwas summary statistics. Additional we should check and flag failed chr:pos;
2. add rsID
Add rsID according to chr and pos position in the file
Add rsID and chr and pos columns according to marker name like: chr:pos:a1:a2
3. flip strand to forward
According to chr and pos, you can get the two alleles of the forward strand, then compare with the input GWAS summary file, if all ATCG appears, we can flip the GWAS summary file using the A-T and C-G rule.
4. process alleles and effect size (if need to make all summary data using the same effect alleles for each SNP?);
Input two gwas summary files, then use the first one as reference and make sure the second file has the same effect allele. If it needs to switch then effect size need to multiply by -1.
5. filter options,
such as remove SNPs based on hwe pvalue or imputation quality or MAF, delete duplicates, non-biallelic SNPs, insertions, and deletions
6. sort by chr and BP in the output file
7. write out as txt.gz
use Python to output gz directly and index file if possible


# Dependencies
1. pyliftover
    * `pip install pyliftover`
2. numpy
    * `pip install numpy`
3. pandas
    * `pip install pandas`
4. mysql
    * `pip install mysql-connector-python`
5. pyBigWig
    * `pip install pyBigWig`





# Usage
To use the pacakge, please follow the steps below:
1. install dependencies mentioned above
2. Please make sure the system have mysql isntalled
3. Git clone/ download this repository
4. cd into the directory
5. In python, get started with the following steps
    - ```python 
        from . import DataConverter
      ```
    - Setting inital parameters: 
        * input_path
        * output_path
        * input_format (e.g. "hg19")
        * output_format (e.g. "hg38")
    - Create instance of the class with the following code chunk
        - ```python 
            converter = DataConverter(input_path, output_path, input_format,output_format)
        ```
    - Start using the provided functions, e.g.: 
        - ```python 
            df = converter.read_data(input_path, '\t', "#chrom","pos", "rsids" ,"alt", "ref", "maf",    "beta", "sebeta", "pval")
            print(df)
        ```


# Functions Provided

### Read Data

### Save Data

### Flip Over

### Add Rsid's

### Flip Strand

### Align Effect Allele and Effect Size between Two Datasets

**Functions to be Implemented**
### Insert/ Filter/ Delete

### Deduplicate

### Create Tbi Index

### Query UCSC Database (for previous dnSNP version)

