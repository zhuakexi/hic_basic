---
# mm10_TSS.csv.gz

calculate TSS sites for mm10

```
gtf = parse_gtf("/share/Data/ychi/genome/GRCm38/raw/gencode.vM25.annotation.gtf.gz")
gtf = gtf.query('feature == "exon"')
gtf = gtf.assign(exon_number = gtf["attributes"].str.extract('exon_number (\d+)')[0])
gtf = gtf.query('exon_number == "1"')
tss = [row["start"] if row["score"]=="+" else row["end"] for _, row in gtf.iterrows()]
gtf = gtf.assign(txStart = tss)
gtf = gtf.assign(transcript_id = gtf["attributes"].str.extract("transcript_id \"(\w+).")[0])
TSS = gtf[["gene_id","gene_name","transcript_id","seqname","txStart","source"]]
TSS = TSS.copy()
TSS["gene_id"] = TSS["gene_id"].str.extract('(\w+).')[0]
TSS.reset_index(drop=True,inplace=True)
TSS.to_csv("/shareb/ychi/notebook_data/Project/embryo_integrate/test_nb/mm10_TSS.csv.gz",index=False)
```
---