#!/usr/bin/env python

import sys

from Bio import Entrez

Entrez.email = "J.F.J.Laros@lumc.nl"

search_terms = [
    "mental retardation",
    "mental retardation profound",
    "mental retardation range",
    "mental retardation ranging",
    "mental retardation syndrome",
    "mentally retarded",
    "learning disabilities",
    "learning problems",
    "mental deficiency",
    "mental delay",
    "mental detoriation",
    "mental impairment",
    "mental lethargy",
    "mental regression",
    "mental retardation congenital",
    "mental retardation gene",
    "development abnormalities",
    "developmental delay",
    "developmental disorder",
    "developmental retardation",
    "developmental stagnation",
    "developmental decline",
    "intellectual disability",
    "intellectual impairment",
    "intellectual regression",
    "development deficit",
    "development delayed",
    "development delay",
    "development failure",
    "developmental stagnation"
]

def _get_key(record, keys):
    """
    """
    if record.has_key(keys[0]):
        if len(keys) == 1:
            return record[keys[0]]
        return _get_key(record[keys[0]], keys[1:])
    #if
    return None
#_get_key

def gene_disease_pubmed(query):
    """
    """
    batch_size 10000
    search_results = Entrez.read(Entrez.esearch(db="pubmed",
        term=test_query, reldate=730, datetype="pdat", usehistory="y"))
    count = int(search_results["Count"])

    sys.stdout.write("Found {} results\n".format(count))
    for start in range(0, count, batch_size):
        sys.stdout.write(".")
        sys.stdout.flush()
        fetch_handle = Entrez.efetch(db="pubmed", rettype="medline",
            retmode="xml", retstart=start, retmax=batch_size,
            webenv=search_results["WebEnv"],
            query_key=search_results["QueryKey"])

        records = Entrez.parse(fetch_handle)
        for record in records:
            body = _get_key(record, ('MedlineCitation', 'Article'))

            if body:
                title = _get_key(body, ('ArticleTitle', ))
                text = _get_key(body, ('Abstract', 'AbstractText'))

                # Now search title and text for gene names.
            #if
        #for
        fetch_handle.close()
    #for
#gene_disease_pubmed

def main():
    """
    """
    query = " OR ".join(search_terms)

    gene_disease_pubmed(query)
#main

if __name__ == "__main__":
    main()
