import os
import random
import string
import subprocess
from pathlib import Path
from string import Template
from tempfile import NamedTemporaryFile

from hires_utils.mmcif import threedg_to_cif
def surface_pymol(_3dg, png, tmpdir=None):
    """
    Generate a and run an intermediate pymol script to render surface pngs.
    Delete the intermediate script after rendering.
    Script name is randomly generated.
    Input:
        _3dg: _3dg file path
        png: output png file path
    """
    # Generate a random string as the intermediate pymol script name
    letters = string.ascii_lowercase
    script_file_path = Path("".join(random.choice(letters) for i in range(10)) + ".pml")
    if tmpdir is not None:
        script_file_path = tmpdir / script_file_path
    else:
        script_file_path = os.getcwd() / script_file_path
    # Generate a cif file from the 3dg file
    if tmpdir is not None:
        cif_file_path = tmpdir / _3dg.name.with_suffix(".cif")
    else:
        cif_file_path = _3dg.parent / _3dg.name.with_suffix(".cif")
    threedg_to_cif(_3dg, cif_file_path)
    # Generate the intermediate pymol script
    template_file_path = Path(__file__).parent / "surface.pml"
    with open(template_file_path, "r") as f:
        template = Template(f.read())
    script = template.substitute(
        cif=cif_file_path,
        png=png
        )

    # Write the intermediate pymol script
    with open(script_file_path, "w") as f:
        f.write(script)

    # Run the intermediate pymol script
    #code = f"conda run --live-stream -n pymol pymol -cq {script_name}"
    code = f"conda run -n pymol pymol -cq {script_file_path}"
    subprocess.run(code, shell=True)

    # Delete the intermediate pymol script
    os.remove(script_file_path)
def clip_territory_pymol(_3dg, png, tmpdir=None, **args):
    """
    Generate a and run an intermediate pymol script to render slicing pngs.
    Each chromosome has a different color, but currently only support 18 colors.
    Delete the intermediate script and cif file after rendering.
    Tmp file names are randomly generated.
    Input:
        _3dg: _3dg file path
        png: output png file path
        **args: transparent to threedg_to_cif
    """
    # Generate a random string as the intermediate pymol script name
    letters = string.ascii_lowercase
    script_file_path = Path("".join(random.choice(letters) for i in range(10)) + ".pml")
    cif_file_path = Path("".join(random.choice(letters) for i in range(10)) + ".cif")
    if tmpdir is not None:
        script_file_path = tmpdir / script_file_path
        cif_file_path = tmpdir / cif_file_path
    else:
        script_file_path = Path.cwd() / script_file_path
        cif_file_path = Path.cwd() / cif_file_path
    threedg_to_cif(_3dg, cif_file_path, **args)
    # Generate the intermediate pymol script
    template_file_path = Path(__file__).parent / "clip_territory.pml"
    with open(template_file_path, "r") as f:
        template = Template(f.read())
    script = template.substitute(
        cif=cif_file_path,
        png=png
        )

    # Write the intermediate pymol script
    with open(script_file_path, "w") as f:
        f.write(script)

    # Run the intermediate pymol script
    #code = f"conda run --live-stream -n pymol pymol -cq {script_name}"
    code = f"conda run -n pymol pymol -cq {script_file_path}"
    subprocess.run(code, shell=True)

    # Delete the intermediate pymol script
    os.remove(cif_file_path)
    os.remove(script_file_path)
def clip_b_pymol(_3dg, b_factor, png, cmap="magenta green, all, 0.005, 0.02", tmpdir=None, **args):
    """
    Generate a and run an intermediate pymol script to render surface pngs.
    Delete the intermediate script after rendering.
    Script name is randomly generated.
    Input:
        _3dg: _3dg file path
        b_factor: a 3 column tsv (chrom pos b_factor) without header
        png: output png file path
        **args: transparent to threedg_to_cif
    """
    # Generate a random string as the intermediate pymol script name
    letters = string.ascii_lowercase
    script_file_path = Path("".join(random.choice(letters) for i in range(10)) + ".pml")
    cif_file_path = Path("".join(random.choice(letters) for i in range(10)) + ".cif")
    if tmpdir is not None:
        script_file_path = tmpdir / script_file_path
        cif_file_path = tmpdir / cif_file_path
    else:
        script_file_path = Path.cwd() / script_file_path
        cif_file_path = Path.cwd() / cif_file_path
    threedg_to_cif(_3dg, cif_file_path, b_factor, **args)
    # Generate the intermediate pymol script
    template_file_path = Path(__file__).parent / "b_factor.pml"
    with open(template_file_path, "r") as f:
        template = Template(f.read())
    script = template.substitute(
        cif=cif_file_path,
        png=png,
        cmap=cmap
        )

    # Write the intermediate pymol script
    with open(script_file_path, "w") as f:
        f.write(script)

    # Run the intermediate pymol script
    #code = f"conda run --live-stream -n pymol pymol -cq {script_name}"
    code = f"conda run -n pymol pymol -cq {script_file_path}"
    subprocess.run(code, shell=True)

    # Delete the intermediate pymol script
    os.remove(cif_file_path)
    os.remove(script_file_path)
if __name__ == "__main__":
    from io import StringIO

    import pandas as pd
    # clip territory
    from hires_utils.hires_io import parse_3dg
    # replace () in chromosome name with _
    _3dg = parse_3dg("/shareb/ychi/repo/sperm40_GM/3dg_c/GMO1001.clean.20k.4.3dg").reset_index()
    _3dg["chr"] = _3dg["chr"].str.replace("[()]", "_", regex=True)
    _3dg = _3dg.set_index(["chr","pos"], drop=True)
    clip_territory_pymol(
        StringIO(_3dg.to_csv(sep="\t", index=True, header=False)),
        "/share/home/ychi/dev/hic_basic/tests/output/GMO1001.clean.20k.4.territory.png",
        tmpdir="/share/home/ychi/dev/hic_basic/tests/output"
        )
    # cpg
    clip_b_pymol(
        "/shareb/ychi/repo/sperm40_GM/3dg_c/GMO1001.clean.20k.4.3dg",
        "/share/home/ychi/software/dip-c/color/hg19.cpg.20k.txt",
        "/share/home/ychi/dev/hic_basic/tests/output/GMO1001.clean.20k.4.cpg.png",
        tmpdir="/share/home/ychi/dev/hic_basic/tests/output"
        )
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