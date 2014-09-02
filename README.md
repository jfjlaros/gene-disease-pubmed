# PubMed abstract query post processor
This program queries the PubMed database for a list of search terms and
extracts gene symbols from the abstract text and the title.

## Input format
The input consists of one search term per line, for example:

    mental retardation

## Output format
The output is a tab delimited file consisting of three columns:
- Gene name.
- Number of abstracts containing the gene name.
- Space separated list of PubMed IDs.

Example:

    SMARCB1	3	25169878 25169651 25168959
