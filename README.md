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
        from DataConverter import DataConverter
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
Description: This function reads the data to be processed and output it in a formatted way
    
Args:
    - input_path (str): the complete or relative path of the input data
    - Chr_col_name (str): the column name in the original data representing chromosome
    - BP_col_name (str): the column name in the original data representing base pair position
    - SNP_col_name (str): the column name in the original data representing rsID
    - A1_col_name (str): the column name in the original data representing effect allele
    - A2_col_name (str): the column name in the original data representing non-effect allele
    - EAF_col_name (str): the column name in the original data represneting allele frequency for effect allele
    - Beta_col_name (str): the column name in the original data represneting effect size for effect allele
    - Se_col_name (str): the column name in the original data represneting standard error for effect size
    - P_col_name (str): the column name in the original data represneting p-value
    - separate_by (str): the delimiter of the original data, '\t' by default (tab separated)

Returns:
    - A pandas data frame formmated in the following ways:

    | Chr    | BP     | SNP    | A1     | A2     | EAF    | Beta   | Se     | P      |
    | ------ | ------ | ------ | ------ | ------ | ------ | ------ | ------ | ------ |
    | 1      | 438956 | rs4596 | G      | A      | 0.0021 | -0.538 | 0.5802 | 0.3533 |
    | X      | 704956 | rs1234 | T      | C      | 0.0242 | 0.1685 | 0.2469 | 0.0843 |

Example Usage:
```python
    df = converter.read_data(input_path, '\t', "#chrom","pos", "rsids" ,"alt", "ref", "maf", "beta", "sebeta", "pval")
    print(df)
```
### Save Data
Example:
```python
    data = converter.query_data(df)
    res = converter.add_rsid(df, data)
    converter.save_data(res, "add_rsid")
```

### Lift Over
Example:
```python
    converter.liftover(df, "Chr", "BP", )
```

### Add Rsid's
Example:
```python
    data = converter.query_data(df)
    converter.add_rsid(df, data)
```

### Flip Strand
Example:
```python
    data = converter.query_data(df)
    converter.flip_strand(df, data)
```

### Align Effect Allele and Effect Size between Two Datasets
Example:
```python
    print(df)
    print(converter.swap_effect_allele(df))
```

**Functions to be Implemented**
### Insert/ Filter/ Delete

### Deduplicate

### Create Tbi Index

### Query UCSC Database (for previous dnSNP version)

