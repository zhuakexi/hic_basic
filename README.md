# hic_basic
# pseudo bulk analysis
pairs --> .scool --> .cool(s)  
Assume you have a meta file `xx.meta.csv.gz` that stores all pairs file path and sample name in `pairs_c12` col.  
To generate a `xx.scool` file that stores all samples, you can use the following command:  

```
from hic_basic.coolstuff import pairs2cool, hic_pileup
from hic_basic.io import read_meta
filesp = read_meta("xx.meta.csv.gz")
pairs2scool(filesp["pairs_c12"].to_dict(), "xx.scool", "mm10.len.tsv", 20e3)
hic_pileup("xx.scool", {sample: cell_type, ...}, "{}.pileup.cool")
```
you need a chromosome length file, two columns, no header.  
you also need a dict that maps sample name to cell types.
# plot hic matrix
vmax and binsize are super important here.
当只能看到一条对角线时，用大binsize和小vmax。
```
from hic_basic.plot.hic import plot_cool
#This command plot whole genome heatmap with all intra and inter chromosomal interactions.
#Make sure you use a low resolution cooler file here like 1Mb binsize for mouse genome.
plot_cool("xx.cool","title", [slice(0,-1),slice(0,-1)],vmax=500)
# This works for mcool file with multiple resolutions.
plot_cool("xx.mcool::resolutions/1000000","title",[slice(0,-1),slice(0,-1)],vmax=500)
# This plot inter chromosomal interaction of chr1 and chr2.
plot_cool("xx.cool","title", ["chr1","chr2"],vmax=500)
# This plot intra chromosome interaction of chr19
plot_cool("xx.cool","title", "chr19",vmax=500)
```
# plot compartment