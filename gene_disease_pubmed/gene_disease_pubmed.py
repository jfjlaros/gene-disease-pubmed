#!/usr/bin/env python

import collections
import json
import sys

import restkit

from Bio import Entrez

Entrez.email = "J.F.J.Laros@lumc.nl"

hgnc_resource = "http://rest.genenames.org/search/"
hgnc_query = "status:Approved"

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

def abstract_walker(query):
    """
    """
    batch_size = 100
    search_results = Entrez.read(Entrez.esearch(db="pubmed",
        term=query, reldate=730, datetype="pdat", usehistory="y"))
    count = int(search_results["Count"])

    for start in range(0, count, batch_size):
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
                if text:
                    text = text[0]
                else:
                    text = ""

                yield title, text
            #if
        #for
        fetch_handle.close()
    #for
#abstract_walker

def hgnc_genes():
    """
    """
    resource = restkit.Resource(hgnc_resource)
    response = resource.get(hgnc_query, headers={"Accept": "application/json"})
    records = json.load(response.body_stream())

    return map(lambda record: record['symbol'], records['response']['docs'])
#hgnc_genes

def gene_disease_pubmed(handle):
    """
    """
    hits = collections.defaultdict(int)

    sys.stdout.write("Retrieving gene list.\n")
    genes = hgnc_genes()

    sys.stdout.write("Starting analysis.\n")
    query = " OR ".join(search_terms)
    walker = abstract_walker(query)
    indicator = 0
    for abstract in walker:
        if not indicator % 100:
            sys.stdout.write(".")
            sys.stdout.flush()
        #if
        indicator += 1
        title = abstract[0].split()
        text = abstract[1].split()

        for gene in genes:
            if gene in title + text:
                hits[gene] += 1
    #for

    for hit in hits:
        handle.write("{}\t{}\n".format(hit, hits[hit]))
#gene_disease_pubmed

def main():
    """
    """
    gene_disease_pubmed(open("results.txt", "w"))
#main

if __name__ == "__main__":
    main()
