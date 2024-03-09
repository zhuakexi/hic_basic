import os
import random
import shutil
import string
import subprocess
from io import StringIO
from pathlib import Path
from string import Template
from tempfile import NamedTemporaryFile

import pandas as pd
from hic_basic.data import fetch_cent_chromlen
from hires_utils.hires_io import parse_3dg
from hires_utils.mmcif import threedg_to_cif, chrom_rm_suffix
class PyMolRender:
    def __init__(self, template_pml_file, png_file_path, tmpdir=None, conda="pymol"):
        """
        Generate a and run an intermediate pymol script to render pngs.
        Delete the intermediate script after rendering.
        Script name is randomly generated.
        Input:
            template_pml_file: template pymol script file path
            png_file_path: output png file path
            tmpdir: directory to save intermediate files
            conda: conda environment name, None for no conda
        """
        self.template_pml_file = template_pml_file
        self.png_file_path = png_file_path
        self.tmpdir = tmpdir
        self.conda = conda

    def __enter__(self):
        letters = string.ascii_lowercase
        self.script_file_path = self._generate_random_file_path(letters, ".pml")
        self.cif_file_path = self._generate_random_file_path(letters, ".cif")
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        if self.script_file_path.exists():
            self.script_file_path.unlink()
        if self.cif_file_path.exists():
            self.cif_file_path.unlink()

    def _generate_random_file_path(self, letters, extension):
        file_path = Path("".join(random.choice(letters) for _ in range(10)) + extension)
        if self.tmpdir is not None:
            return self.tmpdir / file_path
        else:
            return Path.cwd() / file_path
    def gen_cif(self, _3dg, *args, **kwargs):
        """
        Rewrite the intermediate cif file with the given kwargs.
        """
        threedg_to_cif(_3dg, self.cif_file_path, *args, **kwargs)
        return self.cif_file_path
    def gen_script(self, **kwargs):
        """
        Rewrite the intermediate pymol script with the given kwargs.
        """
        with open(self.template_pml_file, "r") as f:
            template = Template(f.read())
        script = template.substitute(
            cif = self.cif_file_path,
            png = self.png_file_path,
            **kwargs
            )
        with open(self.script_file_path, "w") as f:
            f.write(script)
    def render(self):
        if self.conda is not None:
            code = f"conda run -n {self.conda} pymol -cq {self.script_file_path}"
        else:
            code = f"pymol -cq {self.script_file_path}"
        subprocess.run(code, shell=True)
# --- basic modes --- #
def surface_territory_pymol(_3dg, png, tmpdir=None, conda="pymol", dcmap="dip_c"):
    """
    Render surface, color by each chromosome.
    Input:
        _3dg: _3dg file path
        png: output png file path
        dcmap: discrete color map.
            dip_c: each chromosome has a different rainbow color
        tmpdir: directory to save intermediate files
        conda: conda environment name, None for no conda
    """
    if dcmap == "dip_c":
        with PyMolRender(
            # using b_factor to achieve discrete color map
            Path(__file__).parent / "pml" / "surface_b_factor.pml",
            png, tmpdir, conda
            ) as render:
            _3dg = parse_3dg(_3dg)
            chroms = chrom_rm_suffix(_3dg.index.get_level_values(0)).unique()
            b_factor = pd.merge(
                _3dg.reset_index().assign(
                    new_chr = chrom_rm_suffix(_3dg.index.get_level_values(0))
                ),
                pd.DataFrame(
                    {"new_chr":chroms, "b_factor":range(1, len(chroms)+1)}
                ),
                on="new_chr"
            )[["new_chr","pos","b_factor"]]
            b_factor = b_factor.drop_duplicates(
                subset=["new_chr","pos"],
                keep="first"
            )
            b_factor = b_factor.rename(columns={"new_chr":"chr"})
            render.gen_cif(
                StringIO(_3dg.to_csv(sep="\t", index=True, header=False)),
                StringIO(b_factor.to_csv(sep="\t", index=False, header=False)),
                dupref = True
                )
            vmin, vmax = b_factor["b_factor"].min(), b_factor["b_factor"].max()
            render.gen_script(
                cmap=f"rainbow, all, {vmin}, {vmax}"
                )
            render.render()
    else:
        with PyMolRender(
            Path(__file__).parent / "pml" / "surface_territory.pml",
            png, tmpdir,conda
            ) as render:
            render.gen_cif(_3dg)
            render.gen_script()
            render.render()
    return png
def clip_territory_pymol(_3dg, png, clip=0, slab=2, tmpdir=None, conda="pymol", **args):
    """
    Render clip view, color by each chromosome.
    Input:
        _3dg: _3dg file path
        png: output png file path
        clip: clip position 0 for middle, negative for more back, positive for more front
        tmpdir: directory to save intermediate files
        conda: conda environment name, None for no conda
        **args: transparent to threedg_to_cif
    """
    with PyMolRender(
        Path(__file__).parent / "pml" / "clip_territory.pml",
        png, tmpdir,conda
        ) as render:
        render.gen_cif(_3dg, **args)
        render.gen_script(
            clip=clip,
            slab=slab
            )
        render.render()
def surface_b_pymol(_3dg, b_factor, png, cmap="magenta green, all, 0.005, 0.02", tmpdir=None, conda="pymol", **args):
    """
    Render surface, color by b factor.
    Input:
        _3dg: _3dg file path
        b_factor: a 3 column tsv (chrom pos b_factor) without header
        png: output png file path
        cmap: pymol color map
        tmpdir: directory to save intermediate files
        conda: conda environment name, None for no conda
        **args: transparent to threedg_to_cif
    """
    with PyMolRender(
        Path(__file__).parent / "pml" / "surface_b_factor.pml",
        png, tmpdir, conda
        ) as render:
        render.gen_cif(_3dg, b_factor, **args)
        render.gen_script(cmap=cmap)
        render.render()
def clip_b_pymol(_3dg, b_factor, png, clip=0, slab=2, cmap="magenta green, all, 0.005, 0.02", tmpdir=None, conda="pymol", **args):
    """
    Render clip view, color by b factor.
    Input:
        _3dg: _3dg file path
        b_factor: a 3 column tsv (chrom pos b_factor) without header
        png: output png file path
        clip: clip position 0 for middle, negative for more back, positive for more front
        cmap: pymol color map
        tmpdir: directory to save intermediate files
        conda: conda environment name, None for no conda
        **args: transparent to threedg_to_cif
    """
    with PyMolRender(
        Path(__file__).parent / "pml" / "clip_b_factor.pml",
        png, tmpdir, conda
        ) as render:
        render.gen_cif(_3dg, b_factor, **args)
        render.gen_script(cmap=cmap, clip=clip, slab=slab)
        render.render()
def highlight_surface_b_pymol(_3dg, b_factor, chain, png, cmap="magenta green, chain {}, 0.005, 0.02", tmpdir=None, conda="pymol", **args):
    """
    Render surface, color by b factor, highlight only one chain, other chains are transparent.
    Input:
        _3dg: _3dg file path
        b_factor: a 3 column tsv (chrom pos b_factor) without header
        chain: chain to highlight
        png: output png file path
        cmap: pymol color map
        tmpdir: directory to save intermediate files
        conda: conda environment name, None for no conda
        **args: transparent to threedg_to_cif
    """
    with PyMolRender(
        Path(__file__).parent / "pml" / "highlight_surface_b_factor.pml",
        png, tmpdir, conda
        ) as render:
        render.gen_cif(_3dg, b_factor, **args)
        render.gen_script(
            chain=chain,
            cmap=cmap.format(chain)
            )
        render.render()
# --- useful modes --- #
def clip_single_territory_pymol(_3dg_file, png, target_chroms=["chrX","chrY"], clip=0, slab=2, tmpdir=None, conda=None, **args):
    """
    Render clip view, only color target chromosomes.
    Input:
        _3dg_file: _3dg file path
        png: output png file path
        target_chroms: chromosomes to color; list
        clip: clip position 0 for middle, negative for more back, positive for more front
        tmpdir: directory to save intermediate files
        conda: conda environment name, None for no conda
        **args: transparent to threedg_to_cif
    """
    _3dg = parse_3dg(_3dg_file)
    b_factor = []
    for chrom, _ in _3dg.index:
        if chrom not in target_chroms:
            b = 0
        else:
            b = 1
        b_factor.append(b)
    b_factor = pd.DataFrame(
        {"chrom":_3dg.index.get_level_values(0), "pos":_3dg.index.get_level_values(1), "b_factor":b_factor}
        )
    with PyMolRender(
        Path(__file__).parent / "pml" / "clip_b_factor.pml",
        png, tmpdir, conda
        ) as render:
        render.gen_cif(
            StringIO(_3dg.to_csv(sep="\t", index=True, header=False)),
            StringIO(b_factor.to_csv(sep="\t", index=False, header=False)),
            **args
            )
        render.gen_script(
            cmap="white_magenta, all, 0, 1",
            clip=clip,
            slab=slab
            )
        render.render()
def centelo_relpos(positions_2col, genome, dupref=False):
    centromeres = fetch_cent_chromlen(genome)
    b_factor = []
    for chrom, pos in positions_2col:
        if dupref:
            chrom = chrom_rm_suffix(chrom)
        else:
            pass
        if chrom not in centromeres.index:
            # without reference
            print(f"Warning: {chrom} not in centromeres, set b factor to 0.5")
            relpos = 0.5
            b_factor.append(relpos)
            continue
        if pos < centromeres.loc[chrom, "start"]:
            # left arm
            relpos = (centromeres.loc[chrom, "start"] - pos) / (centromeres.loc[chrom, "start"] - 0 )
        elif pos > centromeres.loc[chrom, "end"]:
            # right arm
            relpos = (pos - centromeres.loc[chrom, "end"]) / (centromeres.loc[chrom, "chrom_length"] - centromeres.loc[chrom, "end"])
        else:
            # centromere
            relpos = 0
        b_factor.append(relpos)
    b_factor = pd.DataFrame(
        list(zip(*zip(*positions_2col), b_factor)),
        columns=["chrom","pos","b_factor"]
    )
    return b_factor
def surface_centelo_pymol(_3dg_file, png, genome="mm10", tmpdir=None, 
                          cif_name=None, dupref=False, conda="pymol", **args):
    """
    Render clip view, color centromere-telomere.
    Input:
        _3dg_file: _3dg file path
        png: output png file path
        genome: genome name, used to fetch centromere position
        tmpdir: directory to save intermediate files
        cif_name: cif file path to dump, if None, a random name will be generated and removed after rendering
        dupref: whether to annote diploid genome with haploid reference
        **args: transparent to threedg_to_cif
    """
    _3dg = parse_3dg(_3dg_file)
    b_factor = centelo_relpos(_3dg.index, genome, dupref)
    with PyMolRender(
        Path(__file__).parent / "pml" / "surface_b_factor.pml",
        png, tmpdir, conda
        ) as render:
        tmp_cif = render.gen_cif(
            StringIO(_3dg.to_csv(sep="\t", index=True, header=False)),
            StringIO(b_factor.to_csv(sep="\t", index=False, header=False)),
            **args
            )
        if cif_name is not None:
            shutil.copy(tmp_cif, cif_name)
        render.gen_script(
            cmap="blue_white_red, all, 0, 1"
            )
        render.render()
        return png
def clip_centelo_pymol(_3dg_file, png, genome="mm10", clip=0, slab=2, tmpdir=None, cif_name=None, dupref=False, conda="pymol", **args):
    """
    Render clip view, color centromere-telomere.
    Input:
        _3dg_file: _3dg file path
        png: output png file path
        genome: genome name, used to fetch centromere position
        tmpdir: directory to save intermediate files
        cif_name: cif file path to dump, if None, a random name will be generated and removed after rendering
        dupref: whether to annote diploid genome with haploid reference
        conda: conda environment name, None for no conda
        **args: transparent to threedg_to_cif
    """
    _3dg = parse_3dg(_3dg_file)
    b_factor = centelo_relpos(_3dg.index, genome, dupref)
    with PyMolRender(
        Path(__file__).parent / "pml" / "clip_b_factor.pml",
        png, tmpdir, conda
        ) as render:
        tmp_cif = render.gen_cif(
            StringIO(_3dg.to_csv(sep="\t", index=True, header=False)),
            StringIO(b_factor.to_csv(sep="\t", index=False, header=False)),
            **args
            )
        if cif_name is not None:
            shutil.copy(tmp_cif, cif_name)
        render.gen_script(
            cmap="blue_white_red, all, 0, 1",
            clip=clip,
            slab=slab
            )
        render.render()
        return png
def single_centelo_relpos(positions_2col, genome, target_chroms, dupref=False):
    centromeres = fetch_cent_chromlen(genome)
    b_factor = []
    for chrom, pos in positions_2col:
        if dupref:
            chrom = chrom_rm_suffix(chrom)
        else:
            pass

        if chrom not in centromeres.index:
            # without reference
            print(f"Warning: {chrom} not in centromeres, set b factor to 0.5")
            relpos = 0.5
            b_factor.append(relpos)
            continue

        if chrom not in target_chroms:
            # set all other chromosomes to 0.5
            relpos = 0.5
            b_factor.append(relpos)
            continue

        if pos < centromeres.loc[chrom, "start"]:
            # left arm
            relpos = (centromeres.loc[chrom, "start"] - pos) / (centromeres.loc[chrom, "start"] - 0 )
        elif pos > centromeres.loc[chrom, "end"]:
            # right arm
            relpos = (pos - centromeres.loc[chrom, "end"]) / (centromeres.loc[chrom, "chrom_length"] - centromeres.loc[chrom, "end"])
        else:
            # centromere
            relpos = 0
        b_factor.append(relpos)
    b_factor = pd.DataFrame(
        list(zip(*zip(*positions_2col), b_factor)),
        columns=["chrom","pos","b_factor"]
    )
    return b_factor
def clip_single_centelo_pymol(_3dg_file, png, target_chroms=["chrX","chrY"], genome="mm10",
                              clip=0, slab=2, tmpdir=None, cif_name=None, dupref=False, conda="pymol", **args):
    """
    Render clip view, color centromere-telomere of target chromosomes. Other chromosomes are set to 0.5/white.
    Input:
        _3dg_file: _3dg file path
        png: output png file path
        target_chroms: chromosomes to color; list
        genome: genome name, used to fetch centromere position
        clip: clip position 0 for middle, negative for more back, positive for more front
        tmpdir: directory to save intermediate files
        cif_name: cif file path to dump, if None, a random name will be generated and removed after rendering
        dupref: whether to annote diploid genome with haploid reference
        **args: transparent to threedg_to_cif
    """
    _3dg = parse_3dg(_3dg_file)
    b_factor = single_centelo_relpos(_3dg.index, genome, target_chroms, dupref)
    with PyMolRender(
        Path(__file__).parent / "pml" / "clip_b_factor.pml",
        png, tmpdir, conda
        ) as render:
        tmp_cif = render.gen_cif(
            StringIO(_3dg.to_csv(sep="\t", index=True, header=False)),
            StringIO(b_factor.to_csv(sep="\t", index=False, header=False)),
            **args
            )
        if cif_name is not None:
            shutil.copy(tmp_cif, cif_name)
        render.gen_script(
            cmap="blue_white_red, all, 0, 1",
            clip=clip,
            slab=slab
            )
        render.render()
        return png
if __name__ == "__main__":
    from io import StringIO

    import pandas as pd
    # clip territory
    from hires_utils.hires_io import parse_3dg
    # intermingling ratio
    intermingling_score = pd.read_csv(
        "/shareb/ychi/repo/sperm_struct/ds_pipeline/smk/intermingle/GMO1001.R3.scores.csv.gz",
        index_col=0
    )[["chrom","start","multi_chrom_intermingling"]]
    clip_b_pymol(
        "/shareb/ychi/repo/sperm40_GM/3dg_c/GMO1001.clean.20k.4.3dg",
        StringIO(intermingling_score.to_csv(sep="\t", index=False, header=False)),
        "/share/home/ychi/dev/hic_basic/tests/output/GMO1001.clean.20k.4.intermingling.ratio.png",
        cmap="white red, all, 0, 1.6",
        tmpdir="/share/home/ychi/dev/hic_basic/tests/output",
        dupref = False
        )