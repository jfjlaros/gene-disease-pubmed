# PubMed abstract query post processor
This program queries the PubMed database for a list of search terms and
extracts gene symbols from the abstract text and the title.


## Installation
Via [pypi](https://pypi.python.org/pypi/gene-disease-pubmed):

    pip install gene-disease-pubmed

From source:

    git clone https://github.com/jfjlaros/gene-disease-pubmed.git
    cd gene-disease-pubmed
    pip install .


## Usage
The program requires an input file that contains one search term per line. An
example can be found [here](data/query.txt).

Because the NCBI requires anyone that uses their services to be able to be
contacted, a valid e-mail address must be supplied.

For an input file named `input.txt`, and the e-mail address being
`j.doe@podunk.us`, the following command can be used:

    gdp input.txt output.tsv j.doe@podunk.us

This will create a tab delimited file named `output.tsv` containing three
columns:

- Gene name.
- Number of abstracts containing the gene name.
- Space separated list of PubMed IDs.

For a full list of options, see the built-in help function:

    gdp -h
