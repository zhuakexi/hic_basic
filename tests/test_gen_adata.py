import os
from numpy.testing import assert_allclose
from hic_basic.wet.adtools import gen_adata
from hic_basic.hicio import read_meta
from hic_basic.wet import basic_filter
# TODO: remove basic filter from here
def test_gen_adata(request, tmp_path):
    cache_dir = tmp_path / "test_gen_adata"
    qc = read_meta(os.path.join(request.fspath.dirname, "data", "testqc.csv.gz"))
    expr = os.path.join(request.fspath.dirname, "data", "counts.gene.tsv.gz")
    adata = gen_adata(qc, cache_dir, expr = ([[expr], dict(zip(qc.index.values, qc.index.values))],{"samplelist":qc.index.values}))
    assert adata.obs.shape[0] == 3
def test_gen_adata_fullinput(request, tmp_path):
    till220515meta = read_meta("/share/Data/ychi/notebook/Project/embryo_integrate/meta/till220515.csv.gz")
    cache_dir = "/shareb/ychi/repo/embryo_integrate/anndatas/till220515_cache/"
    qc = basic_filter(till220515meta)
    velo_files = [
    'anndatas/attached/velocyto/20210823_embryo_9IYWY.loom',
    'anndatas/attached/velocyto/20210912_embryo_0U368.loom',
    'anndatas/attached/velocyto/20211126_embryo_S57FO.loom',
    'anndatas/attached/velocyto/20210712_embryo_re_BP89K.loom',
    "anndatas/attached/velocyto/E0903.loom",
    "anndatas/attached/velocyto/E1022.loom",
    "anndatas/attached/velocyto/1225_e_late2_EIB9Q.loom",
    "anndatas/attached/velocyto/1228_e_late2_T0LER.loom",
    "anndatas/attached/velocyto/20211209_embryo_UBOV7.loom",
    "anndatas/attached/velocyto/220110_embryo_8EZD5.loom",
    "anndatas/attached/velocyto/220111_embryo_69EXU.loom",
    "anndatas/attached/velocyto/220128_embryo_CLONB.loom",
    "anndatas/attached/velocyto/220217_embryo_JUEJP.loom",
    "anndatas/attached/velocyto/220308_embryo_L9X35.loom",
    "anndatas/attached/velocyto/220309_embryo_V3X1M.loom",
    "/shareb/ychi/repo/embryo_integrate/velocyto/Velo/220330_embryo/220330_embryo_G0QHE.loom",
    "/shareb/ychi/repo/embryo_integrate/velocyto/Velo/220406_embryo/220406_embryo_SPN2K.loom",
    "/shareb/ychi/repo/embryo_integrate/velocyto/Velo/220419_embryo/220419_embryo_ZN6HN.loom",
    "/shareb/ychi/repo/embryo_integrate/velocyto/Velo/220427_embryo/220427_embryo_O4IQ9.loom",
    "/shareb/ychi/repo/embryo_integrate/velocyto/Velo/220509_embryo/220509_embryo_82CAL.loom",
    "/shareb/ychi/repo/embryo_integrate/velocyto/Velo/220515_embryo/220515_embryo_E51RI.loom"
    ]
    g1_files = [
    #"/share/home/ychi/data/repo/embryo_integrate/matrix_genome1/counts.gene.tsv.gz",
     "/shareb/ychi/repo/embryo_integrate/SNPsplit/result/block1.g1.gene_expr.tsv.gz",
    #"/share/home/ychi/data/repo/embryo_integrate/SNP_split/count_res/genom1/count_matrix/counts.gene.tsv.gz",
    "/shareb/ychi/repo/embryo_integrate/SNPsplit/result/block2.g1.gene_expr.tsv.gz",
    #"/share/home/ychi/data/repo/embryo_integrate/SNP_split/pipeline_out/result/220110_embryo.g1.gene_expr.tsv.gz",
    #"/share/home/ychi/data/repo/embryo_integrate/SNP_split/pipeline_out/result/220111_embryo.g1.gene_expr.tsv.gz",
    "/shareb/ychi/repo/embryo_integrate/SNPsplit/result/220110_embryo.g1.gene_expr.tsv.gz",
     "/shareb/ychi/repo/embryo_integrate/SNPsplit/result/220111_embryo.g1.gene_expr.tsv.gz",
    "/shareb/ychi/repo/embryo_integrate/SNPsplit/result/220128_embryo.g1.gene_expr.tsv.gz",
    "/shareb/ychi/repo/embryo_integrate/SNPsplit/result/220217_embryo.g1.gene_expr.tsv.gz",
    "/shareb/ychi/repo/embryo_integrate/SNPsplit/result/220309_embryo.g1.gene_expr.tsv.gz",
    "/shareb/ychi/repo/embryo_integrate/SNPsplit/result/220308_embryo.g1.gene_expr.tsv.gz",
    "/shareb/ychi/repo/embryo_integrate/SNPsplit/result/220330_embryo.g1.gene_expr.tsv.gz",
    "/shareb/ychi/repo/embryo_integrate/SNPsplit/result/220406_embryo.g1.gene_expr.tsv.gz",
    "/shareb/ychi/repo/embryo_integrate/SNPsplit/result/220419_embryo.g1.gene_expr.tsv.gz",
    "/shareb/ychi/repo/embryo_integrate/SNPsplit/result/220427_embryo.g1.gene_expr.tsv.gz",
    "/shareb/ychi/repo/embryo_integrate/SNPsplit/result/220509_embryo.g1.gene_expr.tsv.gz",
    "/shareb/ychi/repo/embryo_integrate/SNPsplit/result/220515_embryo.g1.gene_expr.tsv.gz"
    ]
    g2_files = [
    #"/share/home/ychi/data/repo/embryo_integrate/matrix_genome2/counts.gene.tsv.gz",
    "/shareb/ychi/repo/embryo_integrate/SNPsplit/result/block1.g2.gene_expr.tsv.gz",
    #"/share/home/ychi/data/repo/embryo_integrate/SNP_split/count_res/genom2/count_matrix/counts.gene.tsv.gz",
    "/shareb/ychi/repo/embryo_integrate/SNPsplit/result/block2.g2.gene_expr.tsv.gz",
    #"/share/home/ychi/data/repo/embryo_integrate/SNP_split/pipeline_out/result/220110_embryo.g2.gene_expr.tsv.gz",
    #"/share/home/ychi/data/repo/embryo_integrate/SNP_split/pipeline_out/result/220111_embryo.g2.gene_expr.tsv.gz",
    "/shareb/ychi/repo/embryo_integrate/SNPsplit/result/220110_embryo.g2.gene_expr.tsv.gz",
     "/shareb/ychi/repo/embryo_integrate/SNPsplit/result/220111_embryo.g2.gene_expr.tsv.gz",
    "/shareb/ychi/repo/embryo_integrate/SNPsplit/result/220128_embryo.g2.gene_expr.tsv.gz",
    "/shareb/ychi/repo/embryo_integrate/SNPsplit/result/220217_embryo.g2.gene_expr.tsv.gz",
    "/shareb/ychi/repo/embryo_integrate/SNPsplit/result/220309_embryo.g2.gene_expr.tsv.gz",
    "/shareb/ychi/repo/embryo_integrate/SNPsplit/result/220308_embryo.g2.gene_expr.tsv.gz",
    "/shareb/ychi/repo/embryo_integrate/SNPsplit/result/220330_embryo.g2.gene_expr.tsv.gz",
    "/shareb/ychi/repo/embryo_integrate/SNPsplit/result/220406_embryo.g2.gene_expr.tsv.gz",
    "/shareb/ychi/repo/embryo_integrate/SNPsplit/result/220419_embryo.g2.gene_expr.tsv.gz",
    "/shareb/ychi/repo/embryo_integrate/SNPsplit/result/220427_embryo.g2.gene_expr.tsv.gz",
    "/shareb/ychi/repo/embryo_integrate/SNPsplit/result/220509_embryo.g2.gene_expr.tsv.gz",
    "/shareb/ychi/repo/embryo_integrate/SNPsplit/result/220515_embryo.g2.gene_expr.tsv.gz"
    ]
    #adata = gen_adata(qc, cache_dir, expr = None, velo = velo_files, g1 = g1_files, g2 = g2_files)
    #adata = gen_adata(qc, cache_dir, expr = None, velo = velo_files, g1 = g1_files, g2 = g2_files, cdps=None, rs="/shareb/ychi/repo/embryo_integrate/anndatas/till220515_cache/repliscore_cache/", pm=None, g1cs = "/shareb/ychi/repo/embryo_integrate/CompartmentStrength/cs/", g2cs = "/shareb/ychi/repo/embryo_integrate/CompartmentStrength/cs/")
    adata = gen_adata(qc, cache_dir, expr = None, velo = velo_files, 
        g1 = g1_files, g2 = g2_files, cdps=None, 
        rs="/shareb/ychi/repo/embryo_integrate/anndatas/till220515_cache/repliscore_cache/", 
        pm=None, g1cs = "/shareb/ychi/repo/embryo_integrate/CompartmentStrength/cs/", 
        g2cs = "/shareb/ychi/repo/embryo_integrate/CompartmentStrength/cs/",
        annote = ["group","partition","cell_type"],
        g1_UMIs=None, g2_UMIs=None)
    assert adata.obs.shape[0] == qc.shape[0]
    assert len(adata.uns.keys()) == 3
def test_gen_adata_sperm(request, tmp_path):
    qc = read_meta("/share/Data/ychi/notebook/Project/sperm_integrate/meta/till14/till14.mm.csv.gz")
    cache_dir = tmp_path
    adata = gen_adata(qc, cache_dir, [], expr=None, chrom_hap_score = ["/shareb/ychi/repo/sperm12_check_sample_other_snpfile/rd/dump/","/shareb/ychi/repo/sperm13/rd/dump/","/shareb/ychi/repo/sperm14/rd/dump/"])
    assert adata.obs.shape[0] == qc.shape[0] - 1 # one sample has no expression
def test_gen_adata_fullinput2(request, tmp_path):
    g1_files = [
        #"/share/home/ychi/data/repo/embryo_integrate/matrix_genome1/counts.gene.tsv.gz",
        "/shareb/ychi/repo/embryo_integrate/SNPsplit/result/block1.g1.gene_expr.tsv.gz",
        #"/share/home/ychi/data/repo/embryo_integrate/SNP_split/count_res/genom1/count_matrix/counts.gene.tsv.gz",
        "/shareb/ychi/repo/embryo_integrate/SNPsplit/result/block2.g1.gene_expr.tsv.gz",
        #"/share/home/ychi/data/repo/embryo_integrate/SNP_split/pipeline_out/result/220110_embryo.g1.gene_expr.tsv.gz",
        #"/share/home/ychi/data/repo/embryo_integrate/SNP_split/pipeline_out/result/220111_embryo.g1.gene_expr.tsv.gz",
        "/shareb/ychi/repo/embryo_integrate/SNPsplit/result/220110_embryo.g1.gene_expr.tsv.gz",
        "/shareb/ychi/repo/embryo_integrate/SNPsplit/result/220111_embryo.g1.gene_expr.tsv.gz",
        "/shareb/ychi/repo/embryo_integrate/SNPsplit/result/220128_embryo.g1.gene_expr.tsv.gz",
        "/shareb/ychi/repo/embryo_integrate/SNPsplit/result/220217_embryo.g1.gene_expr.tsv.gz",
        "/shareb/ychi/repo/embryo_integrate/SNPsplit/result/220309_embryo.g1.gene_expr.tsv.gz",
        "/shareb/ychi/repo/embryo_integrate/SNPsplit/result/220308_embryo.g1.gene_expr.tsv.gz",
        "/shareb/ychi/repo/embryo_integrate/SNPsplit/result/220330_embryo.g1.gene_expr.tsv.gz",
        "/shareb/ychi/repo/embryo_integrate/SNPsplit/result/220406_embryo.g1.gene_expr.tsv.gz",
        "/shareb/ychi/repo/embryo_integrate/SNPsplit/result/220419_embryo.g1.gene_expr.tsv.gz",
        "/shareb/ychi/repo/embryo_integrate/SNPsplit/result/220427_embryo.g1.gene_expr.tsv.gz",
        "/shareb/ychi/repo/embryo_integrate/SNPsplit/result/220509_embryo.g1.gene_expr.tsv.gz",
        
        "/shareb/ychi/repo/embryo_integrate/SNPsplit/result/220515_embryo.g1.gene_expr.tsv.gz",
        "/shareb/ychi/repo/embryo_integrate/SNPsplit/result/220520_embryo.g1.gene_expr.tsv.gz",
        "/shareb/ychi/repo/embryo_integrate/SNPsplit/result/220530_embryo.g1.gene_expr.tsv.gz"
    ]
    g2_files = [
        #"/share/home/ychi/data/repo/embryo_integrate/matrix_genome2/counts.gene.tsv.gz",
        "/shareb/ychi/repo/embryo_integrate/SNPsplit/result/block1.g2.gene_expr.tsv.gz",
        #"/share/home/ychi/data/repo/embryo_integrate/SNP_split/count_res/genom2/count_matrix/counts.gene.tsv.gz",
        "/shareb/ychi/repo/embryo_integrate/SNPsplit/result/block2.g2.gene_expr.tsv.gz",
        #"/share/home/ychi/data/repo/embryo_integrate/SNP_split/pipeline_out/result/220110_embryo.g2.gene_expr.tsv.gz",
        #"/share/home/ychi/data/repo/embryo_integrate/SNP_split/pipeline_out/result/220111_embryo.g2.gene_expr.tsv.gz",
        "/shareb/ychi/repo/embryo_integrate/SNPsplit/result/220110_embryo.g2.gene_expr.tsv.gz",
        "/shareb/ychi/repo/embryo_integrate/SNPsplit/result/220111_embryo.g2.gene_expr.tsv.gz",
        "/shareb/ychi/repo/embryo_integrate/SNPsplit/result/220128_embryo.g2.gene_expr.tsv.gz",
        "/shareb/ychi/repo/embryo_integrate/SNPsplit/result/220217_embryo.g2.gene_expr.tsv.gz",
        "/shareb/ychi/repo/embryo_integrate/SNPsplit/result/220309_embryo.g2.gene_expr.tsv.gz",
        "/shareb/ychi/repo/embryo_integrate/SNPsplit/result/220308_embryo.g2.gene_expr.tsv.gz",
        "/shareb/ychi/repo/embryo_integrate/SNPsplit/result/220330_embryo.g2.gene_expr.tsv.gz",
        "/shareb/ychi/repo/embryo_integrate/SNPsplit/result/220406_embryo.g2.gene_expr.tsv.gz",
        "/shareb/ychi/repo/embryo_integrate/SNPsplit/result/220419_embryo.g2.gene_expr.tsv.gz",
        "/shareb/ychi/repo/embryo_integrate/SNPsplit/result/220427_embryo.g2.gene_expr.tsv.gz",
        "/shareb/ychi/repo/embryo_integrate/SNPsplit/result/220509_embryo.g2.gene_expr.tsv.gz",
        "/shareb/ychi/repo/embryo_integrate/SNPsplit/result/220515_embryo.g2.gene_expr.tsv.gz",
        "/shareb/ychi/repo/embryo_integrate/SNPsplit/result/220520_embryo.g2.gene_expr.tsv.gz",
        "/shareb/ychi/repo/embryo_integrate/SNPsplit/result/220530_embryo.g2.gene_expr.tsv.gz"
    ]
    velo_files = [
        'anndatas/attached/velocyto/20210823_embryo_9IYWY.loom',
        'anndatas/attached/velocyto/20210912_embryo_0U368.loom',
        'anndatas/attached/velocyto/20211126_embryo_S57FO.loom',
        'anndatas/attached/velocyto/20210712_embryo_re_BP89K.loom',
        "anndatas/attached/velocyto/E0903.loom",
        "anndatas/attached/velocyto/E1022.loom",
        "anndatas/attached/velocyto/1225_e_late2_EIB9Q.loom",
        "anndatas/attached/velocyto/1228_e_late2_T0LER.loom",
        "anndatas/attached/velocyto/20211209_embryo_UBOV7.loom",
        "anndatas/attached/velocyto/220110_embryo_8EZD5.loom",
        "anndatas/attached/velocyto/220111_embryo_69EXU.loom",
        "anndatas/attached/velocyto/220128_embryo_CLONB.loom",
        "anndatas/attached/velocyto/220217_embryo_JUEJP.loom",
        "anndatas/attached/velocyto/220308_embryo_L9X35.loom",
        "anndatas/attached/velocyto/220309_embryo_V3X1M.loom",
        "/shareb/ychi/repo/embryo_integrate/velocyto/Velo/220330_embryo/220330_embryo_G0QHE.loom",
        "/shareb/ychi/repo/embryo_integrate/velocyto/Velo/220406_embryo/220406_embryo_SPN2K.loom",
        "/shareb/ychi/repo/embryo_integrate/velocyto/Velo/220419_embryo/220419_embryo_ZN6HN.loom",
        "/shareb/ychi/repo/embryo_integrate/velocyto/Velo/220427_embryo/220427_embryo_O4IQ9.loom",
        "/shareb/ychi/repo/embryo_integrate/velocyto/Velo/220509_embryo/220509_embryo_82CAL.loom",
        "/shareb/ychi/repo/embryo_integrate/velocyto/Velo/220515_embryo/220515_embryo_E51RI.loom",
        "/shareb/ychi/repo/embryo_integrate/velocyto/Velo/220520_embryo/220520_embryo_TVC4B.loom",
        "/shareb/ychi/repo/embryo_integrate/velocyto/Velo/220530_embryo/220530_embryo_0XBL2.loom"
        ]
    till220530meta = read_meta("/share/Data/ychi/notebook/Project/embryo_integrate/meta/till220530.csv.gz")
    cache_dir = "/shareb/ychi/repo/embryo_integrate/anndatas/till220530_cache/"
    qc = basic_filter(till220530meta)
    adata = gen_adata(qc, cache_dir, expr = None, velo = velo_files, 
        g1 = g1_files, g2 = g2_files, cdps=None, 
        rs="/shareb/ychi/repo/embryo_integrate/anndatas/till220530_cache/repliscore_cache/", 
        pm=None, collect_hour = None, cell_type = None,
        annote = ["group","partition","cell_type"],
        g1_UMIs=None, g2_UMIs=None, g1cs = "/shareb/ychi/repo/embryo_integrate/CompartmentStrength/cs/", 
        g2cs = "/shareb/ychi/repo/embryo_integrate/CompartmentStrength/cs/")
    assert adata.obs.shape[0] == qc.shape[0]
    assert ("collect_hour" in adata.obs.columns) and ("cell_type" in adata.obs.columns)
    #assert len(adata.uns.keys()) == 3