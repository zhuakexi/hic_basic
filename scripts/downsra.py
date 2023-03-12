# single thread, aspera version
# TODO: add multi-threading
import sys
import os
import json
import re
from time import sleep
from subprocess import run
from hashlib import md5
import requests
import pandas as pd
info = sys.argv[1] # infomation csv file, must have ["Run", "name"] columns
outdir = sys.argv[2] # output directory
if len(sys.argv) > 3:
    backend = sys.argv[3]
else:
    backend = "aspera"
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
    "https://www.ebi.ac.uk/ena/portal/api/filereport?accession="+accession+"&result=read_run&fields=fastq_aspera,fastq_ftp,fastq_md5&format=json&download=true"
    for accession in info["Run"]
]
fastq_aspera_urls = []
fastq_ftp_urls = []
fastq_md5s = []
print("Start searching urls...")
# check cache
if os.path.exists(os.path.join(outdir, "fastq_urls.json")):
    print("Cache found, load %s" % os.path.join(outdir, "fastq_urls.json"))
    cached = json.load(open(os.path.join(outdir, "fastq_urls.json")))
else:
    cached = {}
for json_url, accession in zip(json_urls, info["Run"]):
    if accession in cached:
        print("Cache found for %s" % accession)
        fastq_aspera_urls.append(cached[accession]["fastq_aspera_urls"])
        fastq_ftp_urls.append(cached[accession]["fastq_ftp_urls"])
        fastq_md5s.append(cached[accession]["fastq_md5s"])
        continue
    # search if not cached
    dat = json.loads(requests.get(url=json_url).text)
    fastq_aspera_urls.append(
        dat[0]["fastq_aspera"].split(";")
    )
    fastq_ftp_urls.append(
        dat[0]["fastq_ftp"].split(";")
    )
    fastq_md5s.append(
        dat[0]["fastq_md5"].split(";")
    )
    cached[accession] = {
        "fastq_aspera_urls": fastq_aspera_urls[-1],
        "fastq_ftp_urls": fastq_ftp_urls[-1],
        "fastq_md5s": fastq_md5s[-1]
    }
    sleep(0.3) # avoid too many requests
    print("Search url %s done." % accession)
print("Save cache...")
json.dump(cached, open(os.path.join(outdir, "fastq_urls.json"), "w"))

# --- download fastq ---

# generate readable file names
fastq_aspera_urls_flatten = [j for i in fastq_aspera_urls for j in i]
fastq_ftp_urls_flatten = [j for i in fastq_ftp_urls for j in i]
fastq_md5s_flatten = [j for i in fastq_md5s for j in i]
outfiles = [os.path.join(outdir, newname(url)) for url in fastq_aspera_urls_flatten]
# generate download commands
if backend == "aspera":
    keyp = "~/miniconda3/envs/fastq-downloader/etc/asperaweb_id_dsa.openssh"
    commands = [
        "ascp  -vQT -l 500m -P33001 -k 1 -i {keyp} era-fasp@{url} {outfile}".format(
            keyp = keyp,
            url = url,
            outfile = outfile
        )
        for url, outfile in zip(fastq_aspera_urls_flatten, outfiles)
    ]
elif backend == "wget":
    commands = [
        "wget -O {outfile} {url}".format(
            url = url,
            outfile = outfile
        )
        for url, outfile in zip(fastq_ftp_urls_flatten, outfiles)
    ]
elif backend == "curl":
    commands = [
        "curl -o {outfile} {url}".format(
            url = url,
            outfile = outfile
        )
        for url, outfile in zip(fastq_ftp_urls_flatten, outfiles)
    ]
else:
    raise ValueError("Unknown backend.")
# download
for command, outfile in zip(commands, outfiles):
    if os.path.exists(outfile):
        print("File %s exists, skip." % outfile)
        continue
    print("Downloading %s" % outfile)
    run(command, shell=True)
    sleep(0.3) # avoid too many requests
print("All done.")
# --- check md5 ---
damaged = []
for md5value, outfile in zip(fastq_md5s_flatten, outfiles):
    if os.path.exists(outfile):
        print("Checking %s" % outfile)
        # if have newer md5 res cache, use it
        if os.path.exists(outfile + ".md5") and os.path.getmtime(outfile + ".md5") > os.path.getmtime(outfile):
            print("Use cached md5 result.")
            with open(outfile + ".md5", "r") as f:
                md5res = f.read().strip()
            if md5res != md5value:
                print("MD5 mismatch for %s" % outfile)
                print("Expected: %s" % md5value)
                print("Got: %s" % md5res)
                damaged.append(outfile)
            with open(outfile + ".md5", "w") as f:
                f.write(md5res)
            continue
        with open(outfile, "rb") as f:
            md5res = md5(f.read()).hexdigest()
        if md5res != md5value:
            print("MD5 mismatch for %s" % outfile)
            print("Expected: %s" % md5value)
            print("Got: %s" % md5res)
            damaged.append(outfile)
        with open(outfile + ".md5", "w") as f:
            f.write(md5res)
if len(damaged) > 0:
    print("Damaged files found:")
    print("\t" + " ".join(damaged))
    exit(1)
else:
    print("All files are fine.")
    exit(0)