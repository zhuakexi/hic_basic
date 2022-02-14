# lousy helper functions to work with cooler, cooltools stuff
# you need a conda env that names `cooler`, with cooler and cooltools installed
import subprocess
def pairs2cool(filei,fileo,sizef,binsize):
    """
    Generate .cool file from 4DN .pairs file
    Input:
        filei: input pairs file
        fileo: output .cool file
        sizef: chrom size file, 2col(chrom:str, size:int) tsv
        binsize: binsize as you wish
    """
    subprocess.check_output(
        "conda run -n cooler cooler cload pairs -c1 2 -p1 3 -c2 4 -p2 5 %s:%d %s %s " % (sizef,binsize,filei,fileo),
        shell=True)
def mergecool(incools,outcool):
    """
    Merge cool files with same indices to get a consensus heatmap.
    Input:
        incools: list of .cool file path
        outcool: output .cool file path
    """
    subprocess.check_output(
        "conda run -n cooler cooler merge %s %s" %(outcool, incools),
        shell=True
    )
def cl_Balance(filei,threads=8):
    """
    cooler(cl) balance
    """
    subprocess.check_output(
        "conda run -n cooler cooler balance -p %d --force %s" % (threads,filei),
        shell=True
    )
def ct_callTAD(filei,fileo):
    """
    cooltools(ct) call TAD
    """
    subprocess.check_output(
        "conda run -n cooler cooltools diamond-insulation -o %s --append-raw-scores %s 100000" %(fileo,filei),
        shell=True
    )