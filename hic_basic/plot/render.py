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
def clip_b_pymol(_3dg, b_factor, png, cmap="magenta green, all, 0.005, 0.02", tmpdir=None):
    """
    Generate a and run an intermediate pymol script to render surface pngs.
    Delete the intermediate script after rendering.
    Script name is randomly generated.
    Input:
        _3dg: _3dg file path
        b_factor: a 3 column tsv (chrom pos b_factor) without header
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
        cif_file_path = (tmpdir / _3dg).name.with_suffix(".cif")
    else:
        cif_file_path = Path(_3dg).with_suffix(".cif")
    threedg_to_cif(_3dg, cif_file_path, b_factor)
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
    cpg = pd.read_table(
        "/share/home/ychi/software/dip-c/color/hg19.cpg.20k.txt",
        header=None,
        names=["chrom", "pos", "cpg"]
    )
    cpg["chrom"] = cpg["chrom"].apply(lambda x: "chr"+x)
    cpg = StringIO(cpg.to_csv(sep="\t", index=False, header=False))
    # wrap cpg to file-like object
    clip_b_pymol(
        "/shareb/ychi/repo/sperm40_GM/3dg_c/GMO1001.clean.20k.4.3dg",
        cpg,
        "/share/home/ychi/dev/hic_basic/tests/output/GMO1001.clean.20k.4.cpg.png",
        )