import re
from abc import ABC, abstractmethod
from pathlib import Path

import pandas as pd
from .hicio import parse_bed, parse_gtf
# from .genome import GenomeIdeograph
from .data import ref_dir
from .sequence import count_CpG

class Feature(ABC):
    """
    Feature object to hold genome features.
    Attributes:
    TODO:1.Add different fn, like feature.fn, feature.bgr, ...
    """
    SUFFIX = {
        "bed" : "bed",
        "chrom.sizes" : "chrom.sizes",
    }
    def __init__(self, genome_name=None, feature_name=None, file_type="bed", db_dir=None, **kwargs):
        """
        Initialize Feature object with genome name and database directory.
        Input:
            genome_name: name of the genome; string
            feature_name: name of the feature; string
            db_dir: directory of the database; string or None
            **kwargs: additional parameters for the feature, will be used in compile method and fn names.
        """
        self.genome_name = genome_name
        self._data = None
        self._is_compiled = False
        self._is_temp = False
        # Set database directory
        if db_dir is None:
            # Shut down storage mode and use class as a well-defined parser
            self.db_dir = None
            self._is_temp = True
        else:
            self.db_dir = Path(db_dir)
            self.db_dir.mkdir(parents=True, exist_ok=True)
        self.feature_name = feature_name
        self.file_type = file_type

        # Set additional parameters from kwargs
        for key, value in kwargs.items():
            setattr(self, key, value)

        param_values = []
        for value in kwargs.values(): # from python 3.6, dict order is guaranteed
            # change values to safe string representation
            str_val = str(value)
            safe_str = re.sub(r'[^a-zA-Z0-9]+', '_', str_val)
            # TODO: for 20e3 and 20000, give same str representation
            param_values.append(safe_str)

        file_suffix = self.SUFFIX.get(file_type, ".csv.gz")
        base_name = f"{self.genome_name}.{self.feature_name}"
        if param_values:
            param_str = ".".join(param_values)
            base_name += f".{param_str}"
        if not self._is_temp:
            self.fn = self.db_dir / (f"{base_name}.{file_suffix}")
            # set to compiled state if file exists
            if self.fn.exists():
                self._is_compiled = True
                print(f"Feature {self.feature_name} for genome {self.genome_name} already compiled in {self.fn}.")
            else:
                self._is_compiled = False
                print(f"Feature {self.feature_name} for genome {self.genome_name} not compiled yet. Will compile to {self.fn}.")
        else:
            self.fn = None
            self._is_compiled = True

    @property
    def data(self, fn=None):
        """
        Lazy-loaded data property.
        """
        if self._data is None:
            self.load(fn=fn)
        return self._data
    @abstractmethod
    def compile(self, **kwargs):
        """
        Process standardized genomic data formats (e.g., GTF, VCF) and convert them into simplified intermediate formats like BED or CSV.

        This abstract method serves two primary purposes:
        1. Data Transformation: Extracts core features from complex genomic annotations and restructures them into standardized intermediate formats for downstream analysis.
        2. Caching Mechanism: Generates reference files (stored in `self.fn`) that can be directly consumed by external tools like pybedtools, BEDTools, or other bioinformatics pipelines.

        The implementation must be customized by subclasses to handle specific genome assemblies and data formats.
        When processing new genomes, this method pre-processes the raw data to create optimized intermediate representations that improve computational efficiency for subsequent operations.
        
        Set the `self._is_compiled` attribute to True after successful compilation.
        """
        pass
    @abstractmethod
    def load(self, fn=None):
        """
        Load the feature data to inner data structure from the compiled file.
        
        Sets the `self.data` attribute to the loaded data.
        """
        pass
class Chromosomes(Feature):
    """
    Chromosome feature object to hold chromosome sizes and natural order.
    Attributes:
        genome_name: name of the genome; string
        db_dir: directory of the database; string
        feature_name: name of the feature; string
        file_type: type of the file; string
        fn: path to the file; Path
    """
    def __init__(self, genome_name=None, db_dir=None):
        """
        Initialize Chromosome feature object.
        Input:
            genome_name: name of the genome; string
            db_dir: directory of the database; string or None
        """
        super().__init__(genome_name, "chromosomes", "chrom.sizes", db_dir)
    def compile(self, size_file=None, force=False):
        chromosomes = pd.read_table(
            size_file,
            names = ["chrom", "size"],
        )
        chromosomes.to_csv(self.fn, sep="\t", index=False, header=False)
        self._is_compiled = True
        return
    def load(self, fn=None):
        """
        Load chromosome sizes from compiled file.
        """
        if self._is_temp:
            assert fn is not None, "For temporary Chromosome object, fn must be provided."
        else:
            if not self._is_compiled:
                raise RuntimeError("Chromosome sizes not compiled yet. Please call compile() first.")
            else:
                fn = self.fn
        print(f"Loading chromosome sizes from {fn}...")
        data = pd.read_table(
            fn,
            names = ["chrom", "length"],
        )
        data["chrom"] = data["chrom"].astype(
            pd.CategoricalDtype(
                data["chrom"],
                ordered=True
            )
        )
        data = data.set_index("chrom").sort_index()
        self._data = data
        return self._data
class TSS(Feature):
    """
    TSS feature object to hold TSS regions.
    Attributes:
        genome_name: name of the genome; string
        db_dir: directory of the database; string
        feature_name: name of the feature; string
        file_type: type of the file; string
        fn: path to the file; Path
    """
    def __init__(self, genome_name=None, db_dir=None):
        super().__init__(genome_name, "TSS", "bed", db_dir)
    @staticmethod
    def get_tss_region_from_gtf(gtf_path, cache_file=None):
        """
        Get TSS region from GTF file.
        Input:
            gtf_path: path to the GTF file; string
            cache_file: path to the cache file; string or None
        Output:
            TSS table; pd.DataFrame with columns:
                gene_id, gene_name, transcript_id, seqname, txStart, source
            Same gene will have multiple rows if it has multiple transcripts.
        """
        print("Parsing GTF file for TSS regions...")
        gtf = parse_gtf(gtf_path)
        print("Extracting TSS regions...")
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
        if cache_file is not None:
            # cache TSS table
            print(f"Caching TSS table to {cache_file}...")
            TSS.to_csv(cache_file, index=False)
        print("TSS regions extracted.")
        return TSS
    def compile(self, gtf_path, force=False):
        """Compile TSS regions from GTF file."""
        # Implement TSS compilation logic here
        if self._is_compiled and not force:
            print(f"TSS regions already compiled in {self.fn}. Use force=True to recompile.")
            return
        elif not force and self.fn.exists():
            print(f"Skipping compilation, file {self.fn} already exists. Use force=True to recompile.")
            self._is_compiled = True
            return
        else:
            print(f"Compiling TSS regions from {gtf_path}...")
        data = self.get_tss_region_from_gtf(gtf_path, cache_file=None)
        data = data[["seqname", "txStart"]]
        data = data.rename(columns={"seqname": "chrom", "txStart": "end"})
        data = data.assign(
            start=data["end"] - 1,  # Convert to 0-based start
        )
        data["end"] = data["end"].clip(lower=0)  # Ensure end is not negative
        # save to standardized BED format
        data = data[["chrom", "start", "end"]]
        data.to_csv(self.fn, sep="\t", index=False, header=False)
        # set the _is_compiled attribute to True
        self._is_compiled = True
        print(f"TSS regions compiled and saved to {self.fn}.")
        return 
    
    def load(self):
        """Load TSS regions from compiled file."""
        # Implement TSS loading logic here
        if not self._is_compiled:
            raise RuntimeError("TSS regions not compiled yet. Please call compile() first.")
        print(f"Loading TSS regions from {self.fn}...")
        self._data = parse_bed(self.fn)
class CpG(Feature):
    """
    GpG ratio of each genomic bin.
    """
    def __init__(self, genome_name=None, db_dir=None, binsize=1e6, flavor="hickit"):
        """
        Initialize CpG feature object.
        Input:
            genome_name, db_dir: parameters for Feature class
            binsize: length of each genomic bin; int or float
            flavor: see GenomeIdeograph
        """
        super().__init__(genome_name, "CpG", "bed", db_dir, binsize=binsize, flavor=flavor)
    def compile(self, fasta, force=False):
        if self._is_compiled and not force:
            print(f"{self.feature_name} already compiled in {self.fn}. Use force=True to recompile.")
            return
        elif not force and self.fn.exists():
            print(f"Skipping compilation, file {self.fn} already exists. Use force=True to recompile.")
            self._is_compiled = True
            return
        else:
            print(f"Compiling {self.feature_name} from {fasta}...")
        from .genome import GenomeIdeograph
        genome = GenomeIdeograph(self.genome_name)
        # TODO: integrate GenomeIdeograph to Feature frames
        bins = genome.bins(
            self.binsize, bed=True, order=True, flavor=self.flavor
            )
        CpG = count_CpG(
            bins,
            fasta
        )
        CpG = CpG.assign(
            name = CpG["chrom"].astype("string") + ":" + CpG["start"].astype("string") + "-" + CpG["end"].astype("string")
        )
        # chrom, start, end, name, score
        CpG[["chrom", "start", "end", "name", "CpG"]].to_csv(
            self.fn, sep="\t", index=False, header=False)
        # set the _is_compiled attribute to True
        self._is_compiled = True
        return
    def load(self):
        """
        Load CpG regions from compiled file.
        """
        if not self._is_compiled:
            raise RuntimeError(f"{self.feature_name} not compiled yet. Please call compile() first.")
        print(f"Loading CpG regions from {self.fn}...")
        self._data = parse_bed(self.fn)
        return self._data
# class FeatureFactory:
#     REGISTRY = {
#         "tss": TSSFeature,
#         "enhancer": EnhancerFeature,
#         "gene_body": GeneBodyFeature
#     }
    
#     @classmethod
#     def create_feature(cls, feature_type: str, source_path: str, compiled_dir: str) -> GenomicFeature:
#         """创建特定类型的特征对象"""
#         if feature_type not in cls.REGISTRY:
#             raise ValueError(f"未知的特征类型: {feature_type}. 可用类型: {list(cls.REGISTRY.keys())}")
#         return cls.REGISTRY[feature_type](source_path, compiled_dir)
    
#     @classmethod
#     def register_feature(cls, name: str, feature_class):
#         """注册新的特征类型"""
#         if not issubclass(feature_class, GenomicFeature):
#             raise TypeError("特征类必须继承自GenomicFeature")
#         cls.REGISTRY[name] = feature_class