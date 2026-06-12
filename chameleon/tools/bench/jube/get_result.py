#!/usr/bin/env python3
#
# @file get_result.py
#
# @copyright 2020-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
#                      Univ. Bordeaux. All rights reserved.
#
# @version 1.4.0
# @author Florent Pruvost
# @date 2025-12-19
#
import pandas
import click
from elasticsearch import Elasticsearch

@click.command()
@click.option("-e", "--elastic-url", default="http://localhost:9200", help="elasticsearch instance url")
@click.option("-t", "--team", required=True, help="team name")
@click.option("-p", "--project", required=True, help="project name")
@click.option("-c", "--commit", required=True, help="chameleon commit")
def main(
    elastic_url: str,
    team: str,
    project: str,
    commit: str
):
    """Get a result from an elasticsearch database, e.g.
    https://elasticsearch.bordeaux.inria.fr."""
    es = Elasticsearch(elastic_url)
    es_index = team + "-" + project + "_" + "perf"

    search_param = {
      "query": {
        "bool": {
          "must": [
            {"term": {"Commit_sha_chameleon": {"value": commit}}}
          ]
        }
      },
      "size": 1000
    }
    response = es.search(index=es_index, **search_param)
    elastic_docs = response["hits"]["hits"]

    docs = pandas.DataFrame()
    for num, doc in enumerate(elastic_docs):

        # get _source data dict from document
        source_data = doc["_source"]

        # get _id from document
        _id = doc["_id"]

        # create a Series object from doc dict object
        doc_data = pandas.Series(source_data, name = _id)
        doc_data = doc_data.drop(labels=['Commit_date_chameleon', 'Commit_sha_chameleon', 'Commit_sha_guix', 'Commit_sha_guix_hpc', 'Commit_sha_guix_hpcnonfree'])

        # append the Series object to the DataFrame object
        docs = pandas.concat([docs, doc_data.to_frame().T])

    docs = docs.astype({"Nmpi": int, "Nthread": int, "Ngpu": int, "P": int, "Q": int, "M": int, "N": int, "K": int})
    docs = docs.rename(columns=str.lower)
    docs.to_csv(path_or_buf="chameleon.csv", sep=",", index=False)

if __name__ == "__main__":
    main()
