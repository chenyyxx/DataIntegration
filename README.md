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
            converter = DataConverter()
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
- pandas.Data.Frame: return formatted data in the form of pandas Data.Frame in the following ways:

| Chr    | BP     | SNP    | A1     | A2     | EAF    | Beta   | Se     | P      |
| ------ | ------ | ------ | ------ | ------ | ------ | ------ | ------ | ------ |
| 1      | 438956 | rs4596 | G      | A      | 0.0021 | -0.538 | 0.5802 | 0.3533 |
| X      | 704956 | rs1234 | T      | C      | 0.0242 | 0.1685 | 0.2469 | 0.0843 |

Example Usage:
```python
    df = converter.read_data(input_path, "#chrom","pos", "rsids" ,"alt", "ref", "maf", "beta", "sebeta", "pval")
```
### Save Data
Description: function to save the processed data in the tsv form as a gz file

Args:
- output_path (str): the path you want the data to be saved.
- df (pandas.Data.Frame): the processed data to be saved.
- name (str): the output name of the data.
- save_format (str): the saving format. Choose between 'gzip' or 'csv'. Default to gz.

Returns:
pandas.Data.Frame: return filtered data in the form of pandas Data.Frame


Example:
```python
    # process data
    ut = Utility()
    dbSnp153 = ut.load_obj("dbSnp153")
    res = converter.add_rsid(df, dbSnp153)
    # save data to a gz file with file name prefix being add_rsid
    converter.save_data(res, "add_rsid")
```

### Lift Over
Description: Function to lift over genome build

Args:
- df (pandas.Data.Frame): the data to be lifted over
- lo_dict (python dictionary): the lift over dictionary return from the create_lo function
- keep_unconvertible (boolean): if true, the function will keep and mark the rows that are not convertible. Default to False.
- keep_original_version (boolean): if true, the function will keep the Chr + BP of original genome build. Default to False.

Returns:
pandas.Data.Frame: return the data being lifted over to the desired genome build

Example:
```python
    converter.liftover(df, "Chr", "BP", )
```

### Add Rsid's
Descriptions: Function to query and add rs ID for rows missing rsIDs.

Args:
- df (pandas.Data.Frame): the data to be added rs_ids
- data (python dictionary): the dictionary containing required info from dbSnp153

Returns:
- pandas.Data.Frame: return the data being added rs_ids.

Example:
```python
    data = converter.query_data(df)
    converter.add_rsid(df, data)
```

### Flip Strand
Descriptions: Function to flip the input data to forward strand

Args:
- df (pandas.Data.Frame): the data to be flipped to forward strand
- data (python dictionary): the dictionary containing required info from dbSnp153
- keep_unconvertible (boolean): if true, the function will keep and mark the rows that are not flipped. Default to False.

Returns:
- pandas.Data.Frame: return the data being flipped to forward strand

Example:
```python
    data = converter.query_data(df)
    converter.flip_strand(df, data)
```

### Align Effect Allele and Effect Size between Two Datasets

Descriptions: this function will align the effect allele of input data based on a reference data

Args:
- reference (pandas.Data.Frame): the reference table
- df (pandas.Data.Frame): the data to be aligned
- check_error_rows (boolean): if true, the function will output the rows that cannot be aligned. Default to False.

Returns:
- pandas.Data.Frame: return the data with its effect allele being aligned with the reference table.

Example:
```python
    print(df)
    print(converter.swap_effect_allele(df))
```

### Sort by Chr and BP
Description: Function to sort the data based on Chr and BP

Args:
- df (pandas.Data.Frame): the data to be sorted

Returns:
- pandas.Data.Frame: return the sorted data

Example:
```python
    print(df)
    print(converter.swap_effect_allele(df))
```

### Filter bi-allelic
Description: Function to filter only bi-allelic cases in the data

Args:
- df (pandas.Data.Frame): The data frame to be filtered.
- rest (boolean): value indicating wether or not to keep (mark only) the non-bi-allelic cases. Default to False.

Returns:
- pandas.Data.Frame: return filtered data in the form of pandas Data.Frame.

Example:
```python
    print(df)
    print(converter.swap_effect_allele(df))
```

### Deduplicate
Description: Function to drop rows in data containing dduplicate keys (Chr + BP)

Args:
- df (pandas.Data.Frame): The data frame to be deduplicated.

Returns:
- pandas.Data.Frame: return filtered data in the form of pandas Data.Frame.

Example:
```python
    print(df)
    print(converter.swap_effect_allele(df))
```

### Query UCSC Database for dbSNP153 info
Description: Function to query required data from dbSnp153

Args:
- df (pandas.Data.Frame): the data we want more info
- link (str): path or link of the '.bb' file of dbSnp153 

Returns:
- pandas.Data.Frame: return complete information from dbSnp153 as a python dictionary

Example:
```python
    print(df)
    print(converter.swap_effect_allele(df))
```

### Save Object
Description: Function to save python data structure on disk

Args:
- obj (obj): the data structure/object to be saved on disk.
- name (str): the name for the obj to be saved as.

Returns:
- return nothing

Example:
```python
    print(df)
    print(converter.swap_effect_allele(df))
```

### Load Object
Description: Function to load saved python data structure from disk

Args:
- name (str): the name of the saved obj on disk to be loaded

Returns:
- pandas.Data.Frame: return complete information from dbSnp153 as a python dictionary

Example:
```python
    print(df)
    print(converter.swap_effect_allele(df))
```


**Functions to be Implemented**
### Insert/ Filter/ Delete

### Create Tbi Index


