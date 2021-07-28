# DataIntegration
Data Integration tool set for genetic data





# Main functions (create a index and link)
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
4. pyBigWig
    * `pip install pyBigWig`





# Usage
To use the pacakge, please follow the steps below:
1. Install dependencies mentioned above
2. Install the package
    * `pip install dataintegrator`
3. To use the `query_data()` function, you ahve to provide the link/path to the local dbSnp153.bb downloaded from the UCSC website; otherwise, the function will query the data online instead of querying data from local file (which will significantly reduce the run time).
4. In python, get started with the following steps
    - ```python 
        from dataintegrator import DataIntegrator as di
      ```
    - Setting inital parameters: 
        * input_path (the path of the data to be processed)
        * output_path (the path you want the processed result to be saved to)
        * input_format (e.g. "hg19")
        * output_format (e.g. "hg38")

    - Start using the provided functions, e.g.: 
        - ```python 
            input_path = "<path_to_your_data_file>.gz"
            df = di.read_data(input_path, '\t', "#chrom","pos", "rsids" ,"alt", "ref", "maf",    "beta", "sebeta", "pval")
            #print(df)
        ```
5. To view example calls of the main functions, clone this repository and see the `.py` files under the `/examples` directory.

6. Data Integrating work flow:
    - The general data processing pipe line includes the following five steps:

    ```sequence
    Andrew->China: Says Hello
    Note right of China: China thinks \n about it
    China-->Andrew: How are you?
    Andrew->>China: I am good thanks!
    ```

    `read_data >>> clean data >>> (query data from dbSnp153 if needed) >>> process data >>> dave data`
    - Clean data with the following pipe lines:
    `filter bi-allelic cases >>> deduplicate data >>> sort data (recommended)`












# Functions Provided



## 1. Function to read data in formatted ways

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
- pandas.DataFrame: return formatted data in the form of pandas DataFrame in the following ways:

| Chr    | BP     | SNP    | A1     | A2     | EAF    | Beta   | Se     | P      |
| ------ | ------ | ------ | ------ | ------ | ------ | ------ | ------ | ------ |
| 1      | 438956 | rs4596 | G      | A      | 0.0021 | -0.538 | 0.5802 | 0.3533 |
| X      | 704956 | rs1234 | T      | C      | 0.0242 | 0.1685 | 0.2469 | 0.0843 |

Example Usage:
```python
    input_path = "<path_to_your_data_file>.gz"
    df = di.read_data(input_path, '\t', "#chrom","pos", "rsids" ,"alt", "ref", "maf",    "beta", "sebeta", "pval")
    #print(df)
```







## 2. Functions to clean data for further processing

### Filter bi-allelic
Description: Function to filter only bi-allelic cases in the data

Args:
- df (pandas.DataFrame): The data frame to be filtered.
- rest (boolean): value indicating wether or not to keep (mark only) the non-bi-allelic cases. Default to False.

Returns:
- pandas.DataFrame: return filtered data in the form of pandas DataFrame.

Example:
```python 
    bi_allelic = di.filter_bi_allelic(df)
    
    # if you want to check the non-bi-allelic cases, use the following command
    bi_allelic = di.filter_bi_allelic(df, rest=True)
```

### Deduplicate
Description: Function to drop rows in data containing dduplicate keys (Chr + BP)

Args:
- df (pandas.DataFrame): The data frame to be deduplicated.

Returns:
- pandas.DataFrame: return filtered data in the form of pandas DataFrame.

Example:
```python
    deduplicated = di.deduplicate(df)
```

### Sort by Chr and BP
Description: Function to sort the data based on Chr and BP

Args:
- df (pandas.DataFrame): the data to be sorted

Returns:
- pandas.DataFrame: return the sorted data

Example:
```python
    sorted_data = di.sort_by_chr_bp(df)
```






## 3. Functions to query data from dbSnp153 for further processing

### Query UCSC Database for dbSNP153 info
Description: Function to query required data from dbSnp153

Args:
- df (pandas.DataFrame): the data we want more info
- link (str): path or link of the '.bb' file of dbSnp153 

Returns:
- pandas.DataFrame: return complete information from dbSnp153 as a python dictionary

Example:
```python
    link = "<path_to_your_dbSnp153_path>.bb"
    dbSnp153 = di.query_data(df, link) # This will usually take longer time
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
    di.save_obj(dbSnp153, "obj/dbSnp153")
```

### Load Object
Description: Function to load saved python data structure from disk

Args:
- name (str): the name of the saved obj on disk to be loaded

Returns:
- pandas.DataFrame: return complete information from dbSnp153 as a python dictionary

Example:
```python
    dbSnp153 = di.load_obj("obj/dbSnp153.pkl")
```








## 4. Functions to process data

### Lift Over
Description: Function to lift over genome build

Args:
- df (pandas.DataFrame): the data to be lifted over
- lo_dict (python dictionary): the lift over dictionary return from the create_lo function
- keep_unconvertible (boolean): if true, the function will keep and mark the rows that are not convertible. Default to False.
- keep_original_version (boolean): if true, the function will keep the Chr + BP of original genome build. Default to False.

Returns:
pandas.DataFrame: return the data being lifted over to the desired genome build

Use two tables to illustrate output

Example:
```python
    input_format = "hg<**>"
    output_format = "hg<**>"
    lo = create_lo(input_format, output_format) # create chain file as reference for genome-build-lift-over.

    # drop unconvertible rows and keep only the result after lift over.
    lift_over = di.lift_over(df, lo)
    print(lift_over)

    # keep and mark rows that are not convertible
    lift_over = di.lift_over(df, lo, keep_unconvertible=True)

    # keep the original genome build version as separate columns in the data set
    lift_over = di.lift_over(df, lo, keep_orginal_version=True)
```

### Add Rsid's
Descriptions: Function to query and add rs ID for rows missing rsIDs.

Args:
- df (pandas.DataFrame): the data to be added rs_ids
- data (python dictionary): the dictionary containing required info from dbSnp153
- keep_all (boolean): value indicating whether the function should keep all rows in the original dataset. Default to False.
- inplace (boolean): value indicating whether the function should replace the original rsID column with the new added_rsid column. Default to True.
- show_comment (boolean): value indicating whether the function should add a column indicating the status of adding rsID. Default to False.
- show_errors (boolean): value indicating whether the function will output a table containing rows that cannot be properly added rsIDs. Default to False.

    1. "added" : missing rsID in orginal dataset. The Chr+BP key can be found in dbSNP153 and the rsID is successfully added
    2. "same": the original dataset have the same rsID as dbSnp153. No need to modify or add.
    3. "different": the original dataset have different rsID as compared to dbSnp153. Use dbSnp153 153 as reference to repalce the original.
    4. "key not found" : The Chr+BP key in original dataset cannot be found in dbSnp153. Fill in NA value. Mark in the comment column.

Returns:
- pandas.DataFrame: return the data being added rs_ids.


Example Call:
```python
    added_rsid = di.add_rsid(df, dbSnp153)
```

A few example output:

`inplace = False, show_comment=False, keep_all=False`
| Chr    | BP     | SNP    | A1     | A2     | EAF    | Beta   | Se     | P      | added_rsid|
| ------ | ------ | ------ | ------ | ------ | ------ | ------ | ------ | ------ | --------- |
| 1      | 438956 | <NA>   | G      | A      | 0.0021 | -0.538 | 0.5802 | 0.3533 | rs12445   |
| X      | 704956 | rs1234 | T      | C      | 0.0242 | 0.1685 | 0.2469 | 0.0843 | rs1234    |

`inplace = True, show_comment=False, keep_all=False`
| Chr    | BP     | SNP    | A1     | A2     | EAF    | Beta   | Se     | P      |
| ------ | ------ | ------ | ------ | ------ | ------ | ------ | ------ | ------ |
| 1      | 438956 | rs12445| G      | A      | 0.0021 | -0.538 | 0.5802 | 0.3533 |
| X      | 704956 | rs1234 | T      | C      | 0.0242 | 0.1685 | 0.2469 | 0.0843 |


`inplace = False, show_comment=True, keep_all=False`
| Chr    | BP     | SNP    | A1     | A2     | EAF    | Beta   | Se     | P      | added_rsid| comment|
| ------ | ------ | ------ | ------ | ------ | ------ | ------ | ------ | ------ | --------- | ------ |
| 1      | 438956 | <NA>   | G      | A      | 0.0021 | -0.538 | 0.5802 | 0.3533 | rs12445   | "added"|
| X      | 704956 | rs1234 | T      | C      | 0.0242 | 0.1685 | 0.2469 | 0.0843 | rs1234    | "same" |

`check_errors=True`
| Chr    | BP     | SNP    | A1     | A2     | EAF    | Beta   | Se     | P      | added_rsid| comment|
| ------ | ------ | ------ | ------ | ------ | ------ | ------ | ------ | ------ | --------- | ------ |
| 1      | 438956 | <NA>   | G      | A      | 0.0021 | -0.538 | 0.5802 | 0.3533 | <NA>      | "key not found"|
| X      | 704956 | <NA>   | T      | C      | 0.0242 | 0.1685 | 0.2469 | 0.0843 | <NA>      | "key not found" |

### Flip Strand
Descriptions: Function to flip the input data to forward strand

Args:
- df (pandas.DataFrame): the data to be flipped to forward strand
- data (python dictionary): the dictionary containing required info from dbSnp153
- keep_all (boolean): value indicating whether the function should keep all rows in the original dataset. Default to False.
- inplace (boolean): value indicating whether the function should replace the original A1 and A2 columns with the new_A1 and new_A2 columns. Default to False.
- show_comment (boolean): value indicating whether the function should add a column indicating the status of flipping strand. Default to False.
- show_errors (boolean): value indicating whether the function will output a table containing rows where strand cannot be properly flipped. Default to False.

    1. "flipped" : The Chr+BP key can be found in dbSNP153 and the strand is successfully flipped.
    2. "same": the original dataset uses the same strand as dbSnp153. No need to modify or add.
    3. "different": the original data set and its correspondence in dbSnp153 show completely different strand patter that cannot be flipped (e.g. T/C vs. C/A)
    4. "dbSnp153: Indel" : the A1 and A2 in dbSnp153 correponds to the Chr+BP in the processed data set contain Indel
    5. "key not found" : The Chr+BP key in original dataset cannot be found in dbSnp153. Fill in NA value. Mark in the comment column.

Returns:
- pandas.DataFrame: return the data being flipped to forward strand

Example:
```python
    flipped = di.flip_strand(df, dbSnp153)
```

A few examples:

`inplace = False, show_comment=False, keep_all=False`
| Chr    | BP     | SNP    | A1     | A2     | EAF    | Beta   | Se     | P      | new_A1 | new_A2 |
| ------ | ------ | ------ | ------ | ------ | ------ | ------ | ------ | ------ | ------ | ------ |
| 1      | 438956 | <NA>   | G      | A      | 0.0021 | -0.538 | 0.5802 | 0.3533 | A      | G      |
| X      | 704956 | rs1234 | T      | C      | 0.0242 | 0.1685 | 0.2469 | 0.0843 | C      | T      |

`inplace = True, show_comment=False, keep_all=False`
| Chr    | BP     | SNP    | A1     | A2     | EAF    | Beta   | Se     | P      |
| ------ | ------ | ------ | ------ | ------ | ------ | ------ | ------ | ------ |
| 1      | 438956 | <NA>   | A      | G      | 0.0021 | -0.538 | 0.5802 | 0.3533 |
| X      | 704956 | rs1234 | C      | T      | 0.0242 | 0.1685 | 0.2469 | 0.0843 |


`inplace = False, show_comment=True, keep_all=False`
| Chr    | BP     | SNP    | A1     | A2     | EAF    | Beta   | Se     | P      | new_A1 | new_A2 | comment|
| ------ | ------ | ------ | ------ | ------ | ------ | ------ | ------ | ------ | ------ | ------ | ------ |
| 1      | 438956 | <NA>   | G      | A      | 0.0021 | -0.538 | 0.5802 | 0.3533 | G      | A      | "kept original"|
| 13     | 704956 | <NA>   | T      | C      | 0.0242 | 0.1235 | 0.2469 | 0.0673 | C      | T      | "flipped" |
| 22     | 568952 | <NA>   | A      | G      | 0.0267 | 0.7485 | 0.7869 | 0.0843 | <NA>   | <NA>   | "key not found" |
| X      | 274586 | <NA>   | C      | T      | 0.0243 | 0.1357 | 0.2435 | 0.1243 | G      | T      | "different" |


`check_errors=True`
| Chr    | BP     | SNP    | A1     | A2     | EAF    | Beta   | Se     | P      | new_A1 | new_A2 | comment|
| ------ | ------ | ------ | ------ | ------ | ------ | ------ | ------ | ------ | ------ | ------ | ------ |
| 22     | 568952 | <NA>   | A      | G      | 0.0267 | 0.7485 | 0.7869 | 0.0843 | <NA>   | <NA>   | "key not found" |
| X      | 274586 | <NA>   | C      | T      | 0.0243 | 0.1357 | 0.2435 | 0.1243 | G      | T      | "different" |

### Align Effect Allele and Effect Size between Two Datasets

Descriptions: this function will align the effect allele of input data based on a reference data

Args:
- reference (pandas.DataFrame): the reference table
- df (pandas.DataFrame): the data to be aligned
- check_error_rows (boolean): if true, the function will output the rows that cannot be aligned. Default to False.

Returns:
- pandas.DataFrame: return the data with its effect allele being aligned with the reference table.

Example:
```python
    reference_path = "<path_to_your_reference_data_file>.gz"
    reference_df = di.read_data(input_path, "chromosome","base_pair_location", "variant_id" ,"effect_allele", "other_allele", "effect_allele_frequency", "beta", "standard_error", "p_value") # for example
    aligned = align_effect_allele(reference_df, df) # be sure df and reference_df are using the same genome build, and both data are properlly cleaned!!

    # if you want to see rows that cannot be aligned
    aligned = align_effect_allele(reference_df, df, check_error_rows=True)
```

## 5. Functions to save result

### Save Data
Description: function to save the processed data in the tsv form as a gz file

Args:
- output_path (str): the path you want the data to be saved.
- df (pandas.DataFrame): the processed data to be saved.
- name (str): the output name of the data.
- save_format (str): the saving format. Choose between 'gzip' or 'csv'. Default to gz.

Returns:
pandas.DataFrame: return filtered data in the form of pandas DataFrame


Example:
```python
    output_path = "<path_to_your_output_directory>" # "result" for example
    di.save_data(output_path, aligned, "aligned")

    # if you want to save the file as csv
    di.save_data(output_path, aligned, "csv")
```




**Functions to be Implemented**
### Insert/ Filter/ Delete

### Create Tbi Index


