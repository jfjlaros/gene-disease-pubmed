#!/usr/bin/env python

import argparse
import collections
import json
import re
import sys

import restkit

from Bio import Entrez

from . import ProtectedFileType, usage, version

def _get_key(record, keys):
    """
    :arg record:
    :type record:
    :arg keys:
    :type keys: tuple(str)

    :returns:
    :rtype:
    """
    if record.has_key(keys[0]):
        if len(keys) == 1:
            return record[keys[0]]
        return _get_key(record[keys[0]], keys[1:])
    #if
    return None
#_get_key

def hgnc_genes(query="status:Approved"):
    """
    :arg query:
    :type query: str

    :returns:
    :rtype:
    """
    hgnc_resource = "http://rest.genenames.org/search/"

    resource = restkit.Resource(hgnc_resource)
    response = resource.get(query, headers={"Accept": "application/json"})

    return json.load(response.body_stream())
#hgnc_genes

def abstract_walker(search_results, batch_size=1000):
    """
    :arg search_results:
    :type search_results:
    :arg batch_size:
    :type batch_size: int

    :returns:
    :rtype:
    """
    for start in range(0, int(search_results["Count"]), batch_size):
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

                pmid = "0"
                for i in _get_key(record, ('PubmedData', 'ArticleIdList')):
                    if i.attributes['IdType'] == 'pubmed':
                        pmid = str(i)

                if text:
                    yield title, text[0], pmid
                else:
                    yield title, "", pmid
            #if
        #for
        fetch_handle.close()
    #for
#abstract_walker

def search_abstracts(search_terms, reldate=730):
    """
    :arg search_terms:
    :type search_terms: list(str)
    :arg reldate:
    :type reldate: int

    :returns:
    :rtype:
    """
    query = "\"{}\"".format("\" OR \"".join(search_terms))

    return Entrez.read(Entrez.esearch(db="pubmed", term=query, reldate=reldate,
        datetype="pdat", usehistory="y"))
#search_abstracts

def gene_disease_pubmed(input_handle, output_handle, log_handle, email,
        progress_indicator=100):
    """
    :arg input_handle:
    :type input_handle: handle
    :arg output_handle:
    :type output_handle: handle
    :arg log_handle:
    :type log_handle: handle
    :arg email:
    :type email: str
    :arg progress_indicator:
    :type progress_indicator: int
    """
    Entrez.email = email
    trim_pattern = "[\W_]+"

    log_handle.write("Retrieving gene list ... ")
    log_handle.flush()
    genes = set(map(lambda x: x['symbol'], hgnc_genes()['response']['docs']))
    log_handle.write("found {}.\n".format(len(genes)))

    log_handle.write("Searching abstracts ... ")
    log_handle.flush()
    search_terms = map(lambda x: x.strip(), input_handle.readlines())
    search_results = search_abstracts(search_terms)
    log_handle.write("found {}.\n".format(search_results["Count"]))

    walker = abstract_walker(search_results)
    log_handle.write("Starting analysis ({} abstracts per dot) ".format(
        progress_indicator))
    log_handle.flush()

    indicator = 0
    hits = collections.defaultdict(lambda: [0, []])
    trim_prefix = re.compile('^{}'.format(trim_pattern))
    trim_suffix = re.compile('{}$'.format(trim_pattern))
    for abstract in walker:
        if not indicator % progress_indicator:
            log_handle.write(".")
            log_handle.flush()
        #if
        indicator += 1
        title = abstract[0].split()
        text = abstract[1].split()
        pmid = abstract[2]
        total_text = set(map(lambda x: trim_prefix.sub('', trim_suffix.sub('',
            x)), set(title + text)))

        for gene in genes & total_text:
            hits[gene][0] += 1
            hits[gene][1].append(pmid)
    #for

    for hit in hits:
        output_handle.write("{}\t{}\t{}\n".format(hit, hits[hit][0],
            ' '.join(hits[hit][1])))
    log_handle.write(" done.\n")
#gene_disease_pubmed

def main():
    """
    """
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=usage[0], epilog=usage[1])
    parser.add_argument('-v', action="version", version=version(parser.prog))
    parser.add_argument('input_handle', metavar='INPUT',
        type=argparse.FileType('r'), help='file containing a list of queries')
    parser.add_argument('output_handle', metavar='OUTPUT',
        type=ProtectedFileType('w'), help='output file')
    parser.add_argument('email', metavar='EMAIL',
        type=str, help='email address (%(type)s)')
    parser.add_argument('-o', dest='log_handle', default=sys.stdout,
        type=argparse.FileType('w'), help='log file (default=<stdout>)')
    parser.add_argument('-p', dest='progress_indicator', default=100,
        type=int, help='report progress after this many abstracts '
        '(%(type)s default=%(default)s)')

    try:
        arguments = parser.parse_args()
    except IOError, error:
        parser.error(error)

    try:
        gene_disease_pubmed(arguments.input_handle, arguments.output_handle,
            arguments.log_handle, arguments.email,
            arguments.progress_indicator)
    except ValueError, error:
        parser.error(error)
#main

if __name__ == "__main__":
    main()
