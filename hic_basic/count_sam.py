import sys
import re
import json
# unique mapped and multi mapped count for every cell in the bam file as STAR definition
# ID must be umi_tools flavor, xxx:xxx:xxxx_cell_umibarcode
# time samtools view -F 256 Aligned.out.bam | python ~/dev/hic_basic/count_sam.py > ../count_mapped.tsv

cells = {}
idr = re.compile('\d+_(\w+)_[N,G,A,T,C]+')
for line in sys.stdin:
    # ID:0, NH:11
    eles = line.split()
    try:
        ID = idr.search(eles[0]).groups()[0]
    except AttributeError:
        print(eles[0])
    if ID not in cells:
        cells[ID] = [0,0]
    NH = int(eles[11].split(":")[2])
    if NH == 1:
        cells[ID][0] += 1
    else:
        cells[ID][1] += 1
for cell in cells:
    # result: cell_id, unique_mapped, multimapped
    print(cell, cells[cell][0], cells[cell][1],sep="\t")