"""
GTF/GFF3 parser for gene track visualization.

This module provides parsers for GTF and GFF3 annotation files,
normalizing them into a unified internal data model.
"""
import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Union


@dataclass
class Interval:
    """Genomic interval with optional phase."""
    start: int
    end: int
    phase: Optional[str] = None
    attrs: Dict = field(default_factory=dict)


@dataclass
class Transcript:
    """Transcript structure with exons, CDS, and UTRs."""
    transcript_id: str
    transcript_name: Optional[str]
    parent_gene_id: str
    seqid: str
    strand: str
    start: int
    end: int
    biotype: Optional[str] = None
    exons: List[Interval] = field(default_factory=list)
    cds: List[Interval] = field(default_factory=list)
    utr5: List[Interval] = field(default_factory=list)
    utr3: List[Interval] = field(default_factory=list)
    utr: List[Interval] = field(default_factory=list)
    attrs: Dict = field(default_factory=dict)

    def get_label(self) -> str:
        """Get transcript label for display."""
        if self.transcript_name:
            return self.transcript_name
        return self.transcript_id


@dataclass
class Gene:
    """Gene structure with transcripts."""
    gene_id: str
    gene_name: Optional[str]
    seqid: str
    strand: str
    start: int
    end: int
    biotype: Optional[str] = None
    transcripts: List[Transcript] = field(default_factory=list)
    attrs: Dict = field(default_factory=dict)

    def get_label(self) -> str:
        """Get gene label for display."""
        if self.gene_name:
            return self.gene_name
        return self.gene_id


class GTFParser:
    """Parser for GTF format annotation files."""

    def __init__(self, filepath: Union[str, Path]):
        self.filepath = Path(filepath)
        self.genes: Dict[str, Gene] = {}
        self.transcripts: Dict[str, Transcript] = {}

    def parse(self) -> Dict[str, Gene]:
        """Parse GTF file and return genes dictionary."""
        with open(self.filepath, 'r') as f:
            for line in f:
                if line.startswith('#'):
                    continue
                fields = line.strip().split('\t')
                if len(fields) < 9:
                    continue

                seqid, source, feature_type, start, end, score, strand, phase, attributes = fields
                start, end = int(start), int(end)
                attrs = self._parse_attributes(attributes)

                self._process_feature(
                    seqid=seqid,
                    feature_type=feature_type,
                    start=start,
                    end=end,
                    strand=strand,
                    phase=phase,
                    attrs=attrs
                )

        # Build gene-transcript hierarchy
        self._build_hierarchy()
        return self.genes

    def _parse_attributes(self, attributes: str) -> Dict[str, str]:
        """Parse GTF attributes string into dictionary."""
        attrs = {}
        # Handle both key "value"; and key "value"; formats
        pattern = r'(\S+)\s+"([^"]*)"'
        for match in re.finditer(pattern, attributes):
            key, value = match.groups()
            attrs[key] = value
        return attrs

    def _process_feature(
        self,
        seqid: str,
        feature_type: str,
        start: int,
        end: int,
        strand: str,
        phase: str,
        attrs: Dict
    ):
        """Process a single GTF feature line."""
        if feature_type == 'gene':
            gene_id = attrs.get('gene_id', attrs.get('ID', ''))
            if not gene_id:
                return
            self.genes[gene_id] = Gene(
                gene_id=gene_id,
                gene_name=attrs.get('gene_name', attrs.get('Name')),
                seqid=seqid,
                strand=strand,
                start=start,
                end=end,
                biotype=attrs.get('gene_biotype', attrs.get('gene_type')),
                attrs=attrs
            )

        elif feature_type == 'transcript':
            transcript_id = attrs.get('transcript_id', attrs.get('ID', ''))
            gene_id = attrs.get('gene_id', attrs.get('Parent', ''))
            if not transcript_id:
                return
            self.transcripts[transcript_id] = Transcript(
                transcript_id=transcript_id,
                transcript_name=attrs.get('transcript_name'),
                parent_gene_id=gene_id,
                seqid=seqid,
                strand=strand,
                start=start,
                end=end,
                biotype=attrs.get('transcript_biotype', attrs.get('transcript_type')),
                attrs=attrs
            )

        elif feature_type in ('exon', 'CDS', 'UTR', 'five_prime_UTR', 'three_prime_UTR',
                               'start_codon', 'stop_codon'):
            transcript_id = attrs.get('transcript_id', attrs.get('Parent', ''))
            if not transcript_id:
                return

            interval = Interval(start=start, end=end, phase=phase, attrs=attrs)

            if transcript_id in self.transcripts:
                transcript = self.transcripts[transcript_id]
                if feature_type == 'exon':
                    transcript.exons.append(interval)
                elif feature_type == 'CDS':
                    transcript.cds.append(interval)
                elif feature_type == 'five_prime_UTR':
                    transcript.utr5.append(interval)
                elif feature_type == 'three_prime_UTR':
                    transcript.utr3.append(interval)
                elif feature_type == 'UTR':
                    transcript.utr.append(interval)

    def _build_hierarchy(self):
        """Build gene-transcript hierarchy and infer missing data."""
        # Associate transcripts with genes
        for transcript in self.transcripts.values():
            if transcript.parent_gene_id in self.genes:
                self.genes[transcript.parent_gene_id].transcripts.append(transcript)
            else:
                # Create gene if not exists
                gene = Gene(
                    gene_id=transcript.parent_gene_id,
                    gene_name=None,
                    seqid=transcript.seqid,
                    strand=transcript.strand,
                    start=transcript.start,
                    end=transcript.end,
                    transcripts=[transcript]
                )
                self.genes[transcript.parent_gene_id] = gene

        # Infer transcript spans from exons if needed
        for transcript in self.transcripts.values():
            if transcript.exons and (transcript.start is None or transcript.end is None):
                transcript.start = min(e.start for e in transcript.exons)
                transcript.end = max(e.end for e in transcript.exons)

        # Infer gene spans from transcripts
        for gene in self.genes.values():
            if gene.transcripts:
                gene.start = min(t.start for t in gene.transcripts)
                gene.end = max(t.end for t in gene.transcripts)


class GFF3Parser:
    """Parser for GFF3 format annotation files."""

    def __init__(self, filepath: Union[str, Path]):
        self.filepath = Path(filepath)
        self.genes: Dict[str, Gene] = {}
        self.transcripts: Dict[str, Transcript] = {}
        self._id_to_feature: Dict[str, Union[Gene, Transcript]] = {}

    def parse(self) -> Dict[str, Gene]:
        """Parse GFF3 file and return genes dictionary."""
        # First pass: collect all features by ID
        features_by_id: Dict[str, dict] = {}
        children_by_parent: Dict[str, List[dict]] = {}

        with open(self.filepath, 'r') as f:
            for line in f:
                if line.startswith('#'):
                    continue
                fields = line.strip().split('\t')
                if len(fields) < 9:
                    continue

                seqid, source, feature_type, start, end, score, strand, phase, attributes = fields
                start, end = int(start), int(end)
                attrs = self._parse_attributes(attributes)

                feature_id = attrs.get('ID')
                parent_id = attrs.get('Parent')

                feature_data = {
                    'seqid': seqid,
                    'source': source,
                    'type': feature_type,
                    'start': start,
                    'end': end,
                    'strand': strand,
                    'phase': phase,
                    'attrs': attrs
                }

                if feature_id:
                    features_by_id[feature_id] = feature_data
                    self._id_to_feature[feature_id] = None  # Placeholder

                if parent_id:
                    if parent_id not in children_by_parent:
                        children_by_parent[parent_id] = []
                    children_by_parent[parent_id].append(feature_data)

        # Second pass: build hierarchy
        for feature_id, feature_data in features_by_id.items():
            feature_type = feature_data['type']
            attrs = feature_data['attrs']

            if feature_type in ('gene', 'pseudogene'):
                gene = self._create_gene(feature_data)
                self.genes[gene.gene_id] = gene
                self._id_to_feature[feature_id] = gene

        # Process transcripts and subfeatures
        for parent_id, children in children_by_parent.items():
            if parent_id in self._id_to_feature:
                parent = self._id_to_feature[parent_id]
                if isinstance(parent, Gene):
                    for child in children:
                        if child['type'] in ('mRNA', 'transcript', 'ncRNA', 'lnc_RNA',
                                              'miRNA', 'rRNA', 'tRNA', 'snRNA', 'snoRNA'):
                            transcript = self._create_transcript(child, parent.gene_id)
                            parent.transcripts.append(transcript)
                            self.transcripts[transcript.transcript_id] = transcript
                            self._id_to_feature[child['attrs'].get('ID')] = transcript

        # Process subfeatures (exons, CDS, UTRs)
        for parent_id, children in children_by_parent.items():
            if parent_id in self.transcripts:
                transcript = self.transcripts[parent_id]
                for child in children:
                    self._process_subfeature(child, transcript)

        return self.genes

    def _parse_attributes(self, attributes: str) -> Dict[str, str]:
        """Parse GFF3 attributes string into dictionary."""
        attrs = {}
        for part in attributes.split(';'):
            part = part.strip()
            if '=' in part:
                key, value = part.split('=', 1)
                attrs[key.strip()] = value.strip()
        return attrs

    def _create_gene(self, feature_data: dict) -> Gene:
        """Create Gene from feature data."""
        attrs = feature_data['attrs']
        return Gene(
            gene_id=attrs.get('ID', ''),
            gene_name=attrs.get('Name', attrs.get('gene_name')),
            seqid=feature_data['seqid'],
            strand=feature_data['strand'],
            start=feature_data['start'],
            end=feature_data['end'],
            biotype=attrs.get('biotype', attrs.get('gene_biotype')),
            attrs=attrs
        )

    def _create_transcript(self, feature_data: dict, parent_gene_id: str) -> Transcript:
        """Create Transcript from feature data."""
        attrs = feature_data['attrs']
        return Transcript(
            transcript_id=attrs.get('ID', ''),
            transcript_name=attrs.get('Name', attrs.get('transcript_name')),
            parent_gene_id=parent_gene_id,
            seqid=feature_data['seqid'],
            strand=feature_data['strand'],
            start=feature_data['start'],
            end=feature_data['end'],
            biotype=attrs.get('biotype', attrs.get('transcript_biotype')),
            attrs=attrs
        )

    def _process_subfeature(self, feature_data: dict, transcript: Transcript):
        """Process subfeature (exon, CDS, UTR) and add to transcript."""
        feature_type = feature_data['type']
        attrs = feature_data['attrs']
        interval = Interval(
            start=feature_data['start'],
            end=feature_data['end'],
            phase=feature_data['phase'],
            attrs=attrs
        )

        if feature_type == 'exon':
            transcript.exons.append(interval)
        elif feature_type == 'CDS':
            transcript.cds.append(interval)
        elif feature_type == 'five_prime_UTR':
            transcript.utr5.append(interval)
        elif feature_type == 'three_prime_UTR':
            transcript.utr3.append(interval)
        elif feature_type == 'UTR':
            transcript.utr.append(interval)


def parse_annotation_file(filepath: Union[str, Path]) -> Dict[str, Gene]:
    """
    Parse GTF or GFF3 annotation file.

    Parameters
    ----------
    filepath : str or Path
        Path to annotation file.

    Returns
    -------
    Dict[str, Gene]
        Dictionary of gene_id -> Gene objects.
    """
    filepath = Path(filepath)
    suffix = filepath.suffix.lower()

    if suffix == '.gtf':
        parser = GTFParser(filepath)
        return parser.parse()
    elif suffix in ('.gff', '.gff3'):
        parser = GFF3Parser(filepath)
        return parser.parse()
    else:
        # Auto-detect based on content (skip comments)
        with open(filepath, 'r') as f:
            for line in f:
                if line.startswith('#'):
                    continue
                # Check first non-comment line
                if 'ID=' in line and 'Parent=' in line:
                    parser = GFF3Parser(filepath)
                else:
                    parser = GTFParser(filepath)
                return parser.parse()
        # Empty file
        return {}


def filter_genes_in_region(
    genes: Dict[str, Gene],
    seqid: str,
    start: int,
    end: int,
    padding: int = 0
) -> Dict[str, Gene]:
    """
    Filter genes that overlap with a genomic region.

    Parameters
    ----------
    genes : Dict[str, Gene]
        Dictionary of genes.
    seqid : str
        Chromosome/sequence ID.
    start : int
        Region start position.
    end : int
        Region end position.
    padding : int, optional
        Padding to add to region boundaries. Default is 0.

    Returns
    -------
    Dict[str, Gene]
        Filtered dictionary of genes.
    """
    filtered = {}
    for gene_id, gene in genes.items():
        if gene.seqid != seqid:
            continue
        # Check overlap
        if gene.end + padding >= start and gene.start - padding <= end:
            filtered[gene_id] = gene
    return filtered
