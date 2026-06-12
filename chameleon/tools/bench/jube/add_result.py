#!/usr/bin/env python3
#
# @file add_result.py
#
# @copyright 2020-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
#                      Univ. Bordeaux. All rights reserved.
#
# @version 1.4.0
# @author Florent Pruvost
# @date 2025-12-19
#
import click
import csv
import json
import math
import sys
import time
from elasticsearch import Elasticsearch
from git import Repo
from typing import Any, Dict, List, Union

Row = Dict[str, Union[str, float]]


def open_csv(filename: str) -> List[Dict[str, str]]:
    """
    Open a csv file a return it as dictionary.
    First row is titles.
    """
    csv_rows = []
    with open(filename) as csv_data:
        reader = csv.DictReader(csv_data)
        titles = reader.fieldnames
        for row in reader:
            csv_rows.append(
                {
                    title: row[title]
                    for title in titles
                }
            )
    return csv_rows


def format_entry(row: Row, mpivendor: str, commit_chameleon: Repo, guix_commits: Dict) -> Dict[str, Any]:
    """"format a result"""
    commit_date_chameleon = str(time.strftime("%Y-%m-%d %H:%M:%S", time.gmtime(commit_chameleon.committed_date)))
    commit_sha_chameleon  = str(commit_chameleon.hexsha)
    hostname    = str(row.pop('hostname'))
    algorithm   = str(row.pop('algorithm'))
    precision   = str(row.pop('precision'))
    nmpi        = int(row.pop('nmpi'))
    p           = int(row.pop('p'))
    q           = int(row.pop('q'))
    nthread     = int(row.pop('nthr'))
    ngpu        = int(row.pop('ngpu'))
    m           = int(row.pop('m'))
    n           = int(row.pop('n'))
    k           = int(row.pop('k'))
    cputime     = float(row.pop('cputime'))
    gflops      = float(row.pop('gflops'))
    result = {
        "Commit_date_chameleon": commit_date_chameleon,
        "Commit_sha_chameleon": commit_sha_chameleon,
        "Commit_sha_guix": guix_commits["guix"],
        "Commit_sha_guix_hpc": guix_commits["guix-hpc"],
        "Commit_sha_guix_hpcnonfree": guix_commits["guix-hpc-non-free"],
        "Hostname": hostname,
        "MPIvendor": mpivendor,
        "Algorithm": algorithm,
        "Precision": precision,
        "Nmpi": nmpi,
        "P": p,
        "Q": q,
        "Nthread": nthread,
        "Ngpu": ngpu,
        "M": m,
        "N": n,
        "K": k,
        "Cputime": cputime,
        "Gflops": gflops
    }
    return result

def format_entry_stats(row: Row, mpivendor: str, commit_chameleon: Repo, commit_stats: str,
                       guix_commits: Dict, es: Elasticsearch, es_index: str):
    """"format a result and compute stats: mean and stdev of gflop"""
    err = 0

    # format measures entry
    result = format_entry(row, mpivendor, commit_chameleon, guix_commits)

    # prepare default gflops stats entry if not existing
    result_stats = {
        "Commit_date_chameleon": result['Commit_date_chameleon'],
        "Commit_sha_chameleon": result['Commit_sha_chameleon'],
        "Commit_sha_guix": result['Commit_sha_guix'],
        "Commit_sha_guix_hpc": result['Commit_sha_guix_hpc'],
        "Commit_sha_guix_hpcnonfree": result['Commit_sha_guix_hpcnonfree'],
        "Hostname": result['Hostname'],
        "MPIvendor": result['MPIvendor'],
        "Algorithm": result['Algorithm'],
        "Precision": result['Precision'],
        "Nmpi": result['Nmpi'],
        "P": result['P'],
        "Q": result['Q'],
        "Nthread": result['Nthread'],
        "Ngpu": result['Ngpu'],
        "M": result['M'],
        "N": result['N'],
        "K": result['K'],
        "mean": format(result['Gflops'], '.1f'),
        "stdev": format(result['Gflops']*0.1, '.1f')
    }

    if commit_stats != 'null':
        # search stats data for this commit and given input parameters
        search_param = {
        "query": {
          "bool": {
            "must": [
                { "match": { "Commit_sha_chameleon": commit_stats }},
                { "match": { "Hostname" : result['Hostname'] }},
                { "match": { "MPIvendor" : result['MPIvendor'] }},
                { "match": { "Algorithm" : result['Algorithm'] }},
                { "match": { "Precision" : result['Precision'] }},
                { "match": { "Nmpi" : result['Nmpi'] }},
                { "match": { "P" : result['P'] }},
                { "match": { "Q" : result['Q'] }},
                { "match": { "Nthread" : result['Nthread'] }},
                { "match": { "Ngpu" : result['Ngpu'] }},
                { "match": { "M" : result['M'] }},
                { "match": { "N" : result['N'] }},
                { "match": { "K" : result['K'] }}
            ]
          }
        },
        "size": 1,
        "_source": ["mean", "stdev"],
        }
        response = es.search(index=es_index, **search_param)
        elastic_docs2 = response["hits"]["hits"]

        if len(elastic_docs2) > 0:
            last_stats_data = elastic_docs2[0]["_source"]
            #print("last_stats_data ", last_stats_data)

            # compute formula from https://public.kitware.com/Wiki/CDash:Design#Test_Timing
            alpha = 0.3
            multiplier = 3

            previousMean = float(last_stats_data['mean'])
            previousSD = float(last_stats_data['stdev'])
            #print("previousMean ", previousMean)

            currentV = result['Gflops']
            # just to test: apply a perturbation
            #pert = random.uniform(-previousSD, previousSD)
            #currentV = result['Gflops'] + pert
            #print("currentV ", currentV)

            newMean = (1-alpha)*previousMean + alpha*currentV
            newSD = math.sqrt((1-alpha)*previousSD*previousSD + alpha*(currentV-newMean)*(currentV-newMean))

            # prepare stats data to put in database newMean and newSD
            result_stats = {
                "Commit_date_chameleon": result['Commit_date_chameleon'],
                "Commit_sha_chameleon": result['Commit_sha_chameleon'],
                "Commit_sha_guix": result['Commit_sha_guix'],
                "Commit_sha_guix_hpc": result['Commit_sha_guix_hpc'],
                "Commit_sha_guix_hpcnonfree": result['Commit_sha_guix_hpcnonfree'],
                "Hostname": result['Hostname'],
                "MPIvendor": result['MPIvendor'],
                "Algorithm": result['Algorithm'],
                "Precision": result['Precision'],
                "Nmpi": result['Nmpi'],
                "P": result['P'],
                "Q": result['Q'],
                "Nthread": result['Nthread'],
                "Ngpu": result['Ngpu'],
                "M": result['M'],
                "N": result['N'],
                "K": result['K'],
                "mean": format(newMean, '.1f'),
                "stdev": format(newSD, '.1f')
            }

            # check for regression
            thresholdSD = 0.1*previousMean
            if previousSD < thresholdSD:
                previousSD = thresholdSD
            maxAcceptableDiff = multiplier*previousSD
            diff = abs(currentV-previousMean)
            if diff > maxAcceptableDiff and currentV < previousMean:
                print("Regression: inputs %(Hostname)s, %(MPIvendor)s, %(Algorithm)s, %(Precision)s, %(Nmpi)s, %(P)s, %(Q)s, %(Nthread)s, %(Ngpu)s, %(M)s, %(N)s, %(K)s " % result)
                print("Regression: outputs Gflops={0}, previousMean={1}, diff={2}, maxAcceptableDiff={3}".format(currentV, previousMean, diff, maxAcceptableDiff))
                # make the script fail only if regression for not too small matrices
                if (result['M'] > 10000) and (result['N'] > 10000):
                    err = 1

    return [result_stats, err]

@click.command()
@click.option("-d", "--directory", default=".", help="git working directory")
@click.option("-e", "--elastic-url", default="http://localhost:9200", help="elasticsearch instance url")
@click.option("-t", "--team", required=True, help="team name")
@click.option("-p", "--project", required=True, help="project name")
@click.option("-m", "--mpi", required=True, help="MPI vendor (openmpi, nmad)")
@click.argument("csv-files", nargs=-1)
def main(
    directory: str,
    elastic_url: str,
    team: str,
    project: str,
    mpi: str,
    csv_files: str,
):
    """Add a result to an elasticsearch database."""

    repo = Repo(directory, search_parent_directories=True)
    commit_chameleon = repo.head.commit

    guix_commits = {"guix" : "",
                    "guix-hpc": "",
                    "guix-hpc-non-free" : ""}
    # collect guix commits info
    with open('guix.json') as f:
        guix_describe = json.load(f)
    for index_guix in guix_describe:
        if index_guix["name"] in guix_commits.keys() :
            guix_commits[ index_guix["name"] ] = index_guix["commit"]

    # measures in term of cputimes and gflops
    es = Elasticsearch(elastic_url)
    es_index = team + "-" + project + "_" + "perf"
    if not es.indices.exists(index=es_index):
        es.indices.create(es_index)

    # call this if mapping must be changed (e.g. add a new field)
    mapping_input = {
        "properties": {
            "Commit_date_chameleon": {"type": "date", "format": "yyyy-MM-dd HH:mm:ss"},
            "Commit_sha_chameleon": {"type": "keyword"},
            "Commit_sha_guix": {"type": "keyword"},
            "Commit_sha_guix_hpc": {"type": "keyword"},
            "Commit_sha_guix_hpcnonfree": {"type": "keyword"},
            "Hostname": {"type": "keyword"},
            "MPIvendor": {"type": "keyword"},
            "Algorithm": {"type": "keyword"},
            "Precision": {"type": "keyword"},
            "Nmpi": {"type": "integer"},
            "P": {"type": "integer"},
            "Q": {"type": "integer"},
            "Nthread": {"type": "integer"},
            "Ngpu": {"type": "integer"},
            "M": {"type": "integer"},
            "N": {"type": "integer"},
            "K": {"type": "integer"},
            "Cputime": {"type": "float"},
            "Gflops": {"type": "float"}
        }
    }
    es.indices.put_mapping(index=es_index, body=mapping_input)

    # load data from csv file
    requests = [
        request
        for file in csv_files
            for request in map(
                lambda row: format_entry(row, mpi, commit_chameleon, guix_commits),
                open_csv(file)
            )
    ]

    # insert measures in database
    for request in requests:
        es.index(index=es_index.lower(), document=request)

    # compute stats: mean and stdev of gflops measured
    # database for stats
    es_index_stats = team + "-" + project + "_" + "stats"
    if not es.indices.exists(index=es_index_stats):
        es.indices.create(es_index_stats)

    # call this if mapping must be changed (e.g. add a new field)
    mapping_input_stats = {
        "properties": {
            "Commit_date_chameleon": {"type": "date", "format": "yyyy-MM-dd HH:mm:ss"},
            "Commit_sha_chameleon": {"type": "keyword"},
            "Commit_sha_guix": {"type": "keyword"},
            "Commit_sha_guix_hpc": {"type": "keyword"},
            "Commit_sha_guix_hpcnonfree": {"type": "keyword"},
            "Hostname": {"type": "keyword"},
            "MPIvendor": {"type": "keyword"},
            "Algorithm": {"type": "keyword"},
            "Precision": {"type": "keyword"},
            "Nmpi": {"type": "integer"},
            "P": {"type": "integer"},
            "Q": {"type": "integer"},
            "Nthread": {"type": "integer"},
            "Ngpu": {"type": "integer"},
            "M": {"type": "integer"},
            "N": {"type": "integer"},
            "K": {"type": "integer"},
            "mean": {"type": "float"},
            "stdev": {"type": "float"}
        }
    }
    es.indices.put_mapping(index=es_index_stats, body=mapping_input_stats)

    # search last commit of database stats
    search_param = {
      "sort": [{"Commit_date_chameleon": {"order": "desc"}}],
      "size": 1,
      "_source": ["Commit_sha_chameleon"],
    }
    response = es.search(index=es_index_stats, **search_param)
    elastic_docs = response["hits"]["hits"]
    last_stats_commit = 'null'
    if len(elastic_docs) > 0:
        # search last records commit
        last_stats_commit = elastic_docs[0]["_source"]['Commit_sha_chameleon']
        print("Regression: mean and stdev taken from commit ", last_stats_commit)

    for file in csv_files:
        csvfile = open_csv(file)
        err = 0
        for row in csvfile:
            #print("row ", row)
            [entry, err2] = format_entry_stats(row, mpi, commit_chameleon, last_stats_commit,
                                               guix_commits, es, es_index_stats)
            err = max(err, err2)
            # insert updated mean and stdev in database
            #print("entry ", entry)
            es.index(index=es_index_stats.lower(), document=entry)
        if err == 1:
            sys.exit(1)

if __name__ == "__main__":
    main()
