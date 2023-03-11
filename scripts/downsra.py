# single thread, aspera version
# TODO: add multi-threading
import sys
import os
import json
import re
from time import sleep
from subprocess import run
import requests
import pandas as pd
info = sys.argv[1] # infomation csv file, must have ["Run", "name"] columns
outdir = sys.argv[2] # output directory
info = pd.read_csv(info)

# --- get fastq url from accession number ---
# generate request url
def newname(url):
    match = re.search('/(SRR[\d]+)_?\d?.fastq.gz', url)
    ext = url[match.end(1):]
    Run = match.group(1)
    name = info.loc[info["Run"] == Run, "name"].values[0]
    return name + ext

json_urls=[
    "https://www.ebi.ac.uk/ena/portal/api/filereport?accession="+accession+"&result=read_run&fields=fastq_aspera&format=json&download=true"
    for accession in info["Run"]
]
fastq_ftp_urls = []
for json_url, accession in zip(json_urls, info["Run"]):
    fastq_ftp_urls.append(
        json.loads(requests.get(url=json_url).text)[0]["fastq_aspera"].split(";")
    )
    sleep(0.3) # avoid too many requests
    print("Search url %s done." % accession)
# --- download fastq ---
# generate readable file names
fastq_ftp_urls_flatten = [j for i in fastq_ftp_urls for j in i]
outfiles = [os.path.join(outdir, newname(url)) for url in fastq_ftp_urls_flatten]
# generate download commands
keyp = "~/miniconda3/envs/fastq-downloader/etc/asperaweb_id_dsa.openssh"
commands = [
    "ascp  -vQT -l 500m -P33001 -k 1 -i {keyp} era-fasp@{url} {outfile}".format(
        keyp = keyp,
        url = url,
        outfile = outfile
    )
    for url, outfile in zip(fastq_ftp_urls_flatten, outfiles)
]
# download
for command, outfile in zip(commands, outfiles):
    if os.path.exists(outfile):
        print("File %s exists, skip." % outfile)
        continue
    print("Downloading %s" % outfile)
    run(command, shell=True)
    sleep(0.3) # avoid too many requests
print("All done.")