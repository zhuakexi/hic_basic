import os
import random
import shutil
import string
import subprocess
from io import StringIO
from pathlib import Path
from string import Template
from tempfile import NamedTemporaryFile
from xml.etree import ElementTree as ET

import numpy as np
import pandas as pd
from hic_basic.data import fetch_cent_chromlen
from hires_utils.hires_io import parse_3dg
from hires_utils.mmcif import threedg_to_cif, chrom_rm_suffix


### --- chimera render --- ###


def search_point(p, points, atol=1e-3):
    """
    Iterate over points and check if p is close to any point in points.
    Input:
        p: flattened np.ndarray
        points: dict of flattened np.ndarray
        atol: float, tolerance for checking if two points are close
    Output:
        If close, return True and index of point.
        If not, return False and None.
    NOTE: this is not very efficient. Use KDTree for large number of points.
    """
    for idx, point in points.items():
        if np.isclose(p, point, atol=atol).all():
            return True, idx
    return False, None
def dump_chimera_links(links, outfile, rgb:str="255,255,255", radius=0.1, atol=1e-3):
    """
    Dump links to Chimera marker file (.cmm file).
    Input:
        links: list of tuple of np.ndarray of shape (3,)
            [(corner1, corner2), ...]
        outfile: str, output file name
        rgb: str, color in rgb format
        radius: float, radius of marker
        atol: float, tolerance for checking if two points are close
    Output:
        None
    """
    r,g,b = rgb.split(",")
    root = ET.Element("marker_set", {"name": "marker set 1"})
    point_eles, link_eles = {}, []
    points = {}
    for link in links:
        indices = []
        for corner in link:
            marked, idx = search_point(corner, points, atol=atol)
            if not marked:
                # add new point, and get index
                idx = max(points.keys()) + 1 if points else 0
                points[idx] = corner
                point_eles[idx] = ET.SubElement(
                        root,
                        "marker",
                        {
                            "id": str(idx),
                            "x" : str(corner[0]),
                            "y" : str(corner[1]),
                            "z" : str(corner[2]),
                            "r" : r, "g" : g, "b" : b,
                            "radius" : str(radius)
                        }
                    )
            indices.append(idx)
        link_eles.append(
            ET.SubElement(
                root,
                "link",
                {
                    "id1": str(indices[0]),
                    "id2": str(indices[1]),
                    "r" : r, "g" : g, "b" : b,
                    "radius" : str(radius)
                }
            )
        )
    tree = ET.ElementTree(root)
    with open(outfile, "wb") as file:
        tree.write(file, encoding="utf-8", xml_declaration=True)
    print(f"Dumped to {outfile}")


### --- pymol render --- ###


import numpy as np

def get_rotation_commands(target_vector, retvec=False):
    """
    Calculate PyMOL turn commands to rotate the camera to face a given direction.

    Input:
        target_vector: A list or numpy array representing the unit vector towards which the camera should be oriented.
        retvec: If True, return the rotation degrees along x, y, z axis as well.
    Output:
        A list of string containing a sequence of 'turn' commands for PyMOL to adjust the camera's orientation accordingly.
    """
    # Assume initial camera direction is along the -z axis (towards the positive z-axis of the model space)
    initial_direction = np.array([0, 0, -1])
    
    # Target vector represents the desired viewing direction in model space
    target_direction = np.array(target_vector) / np.linalg.norm(target_vector)  # Normalize target vector
    
    # Compute the rotation matrix that aligns initial_direction with target_direction
    # Find the rotation axis and angle
    rotation_axis = np.cross(initial_direction, target_direction)
    rotation_axis_norm = np.linalg.norm(rotation_axis)
    if rotation_axis_norm == 0:
        # If target_direction is already aligned with initial_direction, no rotation is needed
        x_degrees, y_degrees, z_degrees = 0, 0, 0
        turn_commands = ["turn x, 0", "turn y, 0", "turn z, 0"]
        if retvec:
            return turn_commands, [x_degrees, y_degrees, z_degrees]
        else:
            return turn_commands
    rotation_axis /= rotation_axis_norm
    cos_angle = np.dot(initial_direction, target_direction)
    # trimming cos_angle to prevent floating point error
    angle = np.arccos(np.clip(cos_angle, -1.0, 1.0))
    
    # Convert rotation axis and angle into a rotation matrix
    c = np.cos(angle)
    s = np.sin(angle)
    t = 1 - c
    R = np.array([[t*rotation_axis[0]*rotation_axis[0] + c, t*rotation_axis[0]*rotation_axis[1] - s*rotation_axis[2], t*rotation_axis[0]*rotation_axis[2] + s*rotation_axis[1]],
                  [t*rotation_axis[0]*rotation_axis[1] + s*rotation_axis[2], t*rotation_axis[1]*rotation_axis[1] + c, t*rotation_axis[1]*rotation_axis[2] - s*rotation_axis[0]],
                  [t*rotation_axis[0]*rotation_axis[2] - s*rotation_axis[1], t*rotation_axis[1]*rotation_axis[2] + s*rotation_axis[0], t*rotation_axis[2]*rotation_axis[2] + c]])
    
    # Extract Euler angles from rotation matrix
    sy = np.sqrt(R[0,0] * R[0,0] +  R[1,0] * R[1,0])
    singular = sy < 1e-6
    if not singular:
        x = np.arctan2(R[2,1], R[2,2])
        y = np.arctan2(-R[2,0], sy)
        z = np.arctan2(R[1,0], R[0,0])
    else:  # Gimbal lock
        x = np.arctan2(-R[1,2], R[1,1])
        y = np.arctan2(-R[2,0], sy)
        z = 0

    # Convert radians to degrees
    x_degrees = np.degrees(x)
    y_degrees = np.degrees(y)
    z_degrees = np.degrees(z)

    # Construct the 'turn' command sequence
    turn_commands = [
        f"turn x, {x_degrees:.3f}",
        f"turn y, {y_degrees:.3f}",
        f"turn z, {z_degrees:.3f}"
    ]
    if retvec:
        return turn_commands, [x_degrees, y_degrees, z_degrees]
    return turn_commands
def apply_rotation_commands(df, commands_str):
    def parse_commands(commands_str):
        commands = []
        for cmd in commands_str.split(';'):
            cmd = cmd.strip()
            if not cmd:
                continue
            parts = cmd.split(',', 1)
            if len(parts) != 2:
                raise ValueError("Invalid command format")
            axis_part = parts[0].split()
            if len(axis_part) != 2 or axis_part[0].lower() != 'turn':
                raise ValueError("Command must start with 'turn'")
            axis = axis_part[1].strip().lower()
            if axis not in ['x', 'y', 'z']:
                raise ValueError("Axis must be x, y, or z")
            try:
                angle = float(parts[1].strip())
            except ValueError:
                raise ValueError("Angle must be a number")
            commands.append((axis, angle))
        return commands

    def get_rotation_matrix(axis, theta):
        c = np.cos(theta)
        s = np.sin(theta)
        if axis == 'x':
            return np.array([[1, 0, 0],
                            [0, c, -s],
                            [0, s, c]])
        elif axis == 'y':
            return np.array([[c, 0, s],
                            [0, 1, 0],
                            [-s, 0, c]])
        elif axis == 'z':
            return np.array([[c, -s, 0],
                            [s, c, 0],
                            [0, 0, 1]])
        else:
            raise ValueError("Invalid axis")

    def apply_rotation(df, R):
        index, columns = df.index, df.columns
        points = df[['x', 'y', 'z']].values.T  # Convert to 3xN array
        rotated = R @ points  # Apply rotation matrix
        new_df = pd.DataFrame(
            rotated.T, columns=columns, index=index)
        return new_df

    rotations = parse_commands(commands_str)
    current_df = df.copy()
    for axis, angle in rotations:
        theta = np.radians(angle)
        R = get_rotation_matrix(axis, theta)
        current_df = apply_rotation(current_df, R)
    return current_df

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
        default_kwargs = {
            "turn_cmd" : ""
        }
        kwargs = {**default_kwargs, **kwargs}
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
def surface_territory_pymol(_3dg, png, tmpdir=None, conda="pymol", dcmap="dip_c", turn_cmd="", cif_name=None):
    """
    Render surface, color by each chromosome.
    Input:
        _3dg: _3dg file path
        png: output png file path
        dcmap: discrete color map.
            dip_c: each chromosome has a different rainbow color
        tmpdir: directory to save intermediate files
        conda: conda environment name, None for no conda
        turn_cmd: PyMOL turn command to adjust camera orientation
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
            tmp_cif = render.gen_cif(
                StringIO(_3dg.to_csv(sep="\t", index=True, header=False)),
                StringIO(b_factor.to_csv(sep="\t", index=False, header=False)),
                dupref = True
                )
            if cif_name is not None:
                shutil.copy(tmp_cif, cif_name)
            vmin, vmax = b_factor["b_factor"].min(), b_factor["b_factor"].max()
            render.gen_script(
                cmap=f"rainbow, all, {vmin}, {vmax}",
                turn_cmd=turn_cmd
                )
            render.render()
    else:
        with PyMolRender(
            Path(__file__).parent / "pml" / "surface_territory.pml",
            png, tmpdir,conda
            ) as render:
            tmp_cif = render.gen_cif(_3dg)
            if cif_name is not None:
                shutil.copy(tmp_cif, cif_name)
            render.gen_script(turn_cmd=turn_cmd)
            render.render()
    return png
def clip_territory_pymol(_3dg, png, clip=0, slab=2, tmpdir=None, conda="pymol", turn_cmd="", **args):
    """
    Render clip view, color by each chromosome.
    Input:
        _3dg: _3dg file path
        png: output png file path
        clip: clip position 0 for middle, negative for more back, positive for more front
        tmpdir: directory to save intermediate files
        conda: conda environment name, None for no conda
        turn_cmd: PyMOL turn command to adjust camera orientation
        **args: transparent to threedg_to_cif
    """
    with PyMolRender(
        Path(__file__).parent / "pml" / "clip_territory.pml",
        png, tmpdir,conda
        ) as render:
        render.gen_cif(_3dg, **args)
        render.gen_script(
            clip=clip,
            slab=slab,
            turn_cmd=turn_cmd
            )
        render.render()
def surface_b_pymol(_3dg, b_factor, png, cmap="magenta green, all, 0.005, 0.02", tmpdir=None, conda="pymol", turn_cmd="", **args):
    """
    Render surface, color by b factor.
    Input:
        _3dg: _3dg file path
        b_factor: a 3 column tsv (chrom pos b_factor) without header
        png: output png file path
        cmap: pymol color map
        tmpdir: directory to save intermediate files
        conda: conda environment name, None for no conda
        turn_cmd: PyMOL turn command to adjust camera orientation
        **args: transparent to threedg_to_cif
    """
    with PyMolRender(
        Path(__file__).parent / "pml" / "surface_b_factor.pml",
        png, tmpdir, conda
        ) as render:
        render.gen_cif(_3dg, b_factor, **args)
        render.gen_script(cmap=cmap, turn_cmd=turn_cmd)
        render.render()
def clip_b_pymol(_3dg, b_factor, png, clip=0, slab=2, cmap="magenta green, all, 0.005, 0.02", turn_cmd="", tmpdir=None, conda="pymol", **args):
    """
    Render clip view, color by b factor.
    Input:
        _3dg: _3dg file path
        b_factor: a 3 column tsv (chrom pos b_factor) without header
        png: output png file path
        clip: clip position 0 for middle, negative for more back, positive for more front
        cmap: pymol color map
        turn_cmd: PyMOL turn command to adjust camera orientation
        tmpdir: directory to save intermediate files
        conda: conda environment name, None for no conda
        **args: transparent to threedg_to_cif
    """
    with PyMolRender(
        Path(__file__).parent / "pml" / "clip_b_factor.pml",
        png, tmpdir, conda
        ) as render:
        render.gen_cif(_3dg, b_factor, **args)
        render.gen_script(cmap=cmap, clip=clip, slab=slab, turn_cmd=turn_cmd)
        render.render()
    return png
def highlight_surface_b_pymol(_3dg, b_factor, chain, png, cmap="magenta green, chain {}, 0.005, 0.02", turn_cmd="", tmpdir=None, conda="pymol", **args):
    """
    Render surface, color by b factor, highlight only one chain, other chains are transparent.
    Input:
        _3dg: _3dg file path
        b_factor: a 3 column tsv (chrom pos b_factor) without header
        chain: chain to highlight
        png: output png file path
        cmap: pymol color map
        turn_cmd: PyMOL turn command to adjust camera orientation
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
            cmap=cmap.format(chain),
            turn_cmd=turn_cmd
            )
        render.render()
# --- useful modes --- #
def clip_single_territory_pymol(_3dg_file, png, target_chroms=["chrX","chrY"], clip=0, slab=2, tmpdir=None, conda=None, turn_cmd="", **args):
    """
    Render clip view, only color target chromosomes.
    Input:
        _3dg_file: _3dg file path
        png: output png file path
        target_chroms: chromosomes to color; list
        clip: clip position 0 for middle, negative for more back, positive for more front
        tmpdir: directory to save intermediate files
        conda: conda environment name, None for no conda
        turn_cmd: PyMOL turn command to adjust camera orientation
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
            slab=slab,
            turn_cmd=turn_cmd
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
                          cif_name=None, dupref=False, conda="pymol", turn_cmd="", **args):
    """
    Render clip view, color centromere-telomere.
    Input:
        _3dg_file: _3dg file path
        png: output png file path
        genome: genome name, used to fetch centromere position
        tmpdir: directory to save intermediate files
        cif_name: cif file path to dump, if None, a random name will be generated and removed after rendering
        dupref: whether to annote diploid genome with haploid reference
        turn_cmd: PyMOL turn command to adjust camera orientation
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
            dupref=False, # if dupref in relpos, don't do it again in mmcif
            **args
            )
        if cif_name is not None:
            shutil.copy(tmp_cif, cif_name)
        render.gen_script(
            cmap="blue_white_red, all, 0, 1",
            turn_cmd=turn_cmd
            )
        render.render()
        return png
def clip_centelo_pymol(_3dg_file, png, genome="mm10", clip=0, slab=2, tmpdir=None, 
                       cif_name=None, dupref=False, conda="pymol", turn_cmd="", **args):
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
        turn_cmd: PyMOL turn command to adjust camera orientation
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
            dupref=False, # if dupref in relpos, don't do it again in mmcif
            **args
            )
        if cif_name is not None:
            shutil.copy(tmp_cif, cif_name)
        render.gen_script(
            cmap="blue_white_red, all, 0, 1",
            clip=clip,
            slab=slab,
            turn_cmd=turn_cmd
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
                              clip=0, slab=2, tmpdir=None, cif_name=None, dupref=False, conda="pymol", turn_cmd="", **args):
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
        turn_cmd: PyMOL turn command to adjust camera orientation
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
            slab=slab,
            turn_cmd=turn_cmd
            )
        render.render()
        return png

from io import StringIO

import numpy as np
from hires_utils.hires_io import parse_3dg

from ..structure.measure import primary_views
from ..structure.pileup import sig_primary_coords

def pymol_primary_views(_3dg_file, view="lr", targets=None):
    """
    Rotate a 3dg structure relative to primary axis.
    Input:
        _3dg_file: 3dg file of a single-cell structure
        view: view to look at in pymol, ["lr", "ht", "dv"], means looking towards \
            left-right, head-tail or dorsal-ventral axis.
        targets: flip according to this.
    Return:
        _3dg: a rotated 3dg structure suitable for pymol rendering
    """
    # rotate structure on the fly
    _3dg_primary_views = primary_views(parse_3dg(_3dg_file))
    _3dg =  sig_primary_coords(
            _3dg_primary_views,
            targets if targets is not None else [1, 1, 1],
            parse_3dg(_3dg_file)
            )
    # now x, y, z = ht, dv, lr

# remove () in chrom names, because pymol does not like it
    _3dg = _3dg.reset_index()
    _3dg["chr"] = _3dg["chr"].str.replace("[()]", "_", regex=True)
    # remove bad particles if any
    if len(_3dg_primary_views["bad_particles"]) > 0:
        print("Bad particles found in", _3dg_file)
        _3dg = _3dg.loc[_3dg.index.difference(_3dg_primary_views["bad_particles"])]
    _3dg.set_index(["chr", "pos"], drop=True, inplace=True)

    if view == "lr":
        # pymol sees towards -z/-lr direction by default, do nothing
        _3dg = _3dg
    elif view == "dv":
        # rotate 90 degree around x/ht axis
        theta = np.radians(90)
        R = np.array([
            [1, 0, 0],
            [0, np.cos(theta), -np.sin(theta)],
            [0, np.sin(theta), np.cos(theta)]
        ])
        _3dg = _3dg.dot(R.T)
    elif view == "ht":
        # rotate 90 degree around y/dv axis
        theta = np.radians(90)
        R = np.array([
            [np.cos(theta), 0, np.sin(theta)],
            [0, 1, 0],
            [-np.sin(theta), 0, np.cos(theta)]
        ])
        _3dg = _3dg.dot(R.T)
    else:
        raise ValueError("Unsupported view")
    return _3dg
# --- single-cell rendering --- #
def render_surface_primary_view(_3dg_file, outpng, tmpdir, view="lr", targets=None, **kwargs):
    """
    Rendering primary views of a single-cell structure with pymol.
    Designed for a none-glomerulus nucleus.
    Input:
        _3dg_file: 3dg file of a single-cell structure
        outpng: output png file
        tmpdir: temporary directory
        view: view to render, ["lr", "ht", "dv"], means looking towards \
            left-right, head-tail or dorsal-ventral axis.
    Output:
        outpng path; and wrote outpng file
    """
    _3dg = pymol_primary_views(_3dg_file, view, targets)
    surface_territory_pymol(
        StringIO(_3dg.to_csv(sep="\t", index=True, header=False)),
        outpng,
        tmpdir,
        **kwargs
    )
    return str(outpng)
def render_surface_centelo_primary_view(_3dg_file, outpng, genome, tmpdir, view="lr", targets=None, **kwargs):
    """
    Color primary views of a single-cell structure with relative dist to centromeres and telomeres.
    Designed for a none-glomerulus nucleus.
    Input:
        _3dg_file: 3dg file of a single-cell structure
        outpng: output png file
        genome: genome name of the structure
        tmpdir: temporary directory
        view: view to render, ["lr", "ht", "dv"], means looking towards \
            left-right, head-tail or dorsal-ventral axis.
        targets: flip according to this.
    Output:
        outpng path; and wrote outpng file
    """
    _3dg = pymol_primary_views(_3dg_file, view, targets)
    surface_centelo_pymol(
        StringIO(_3dg.to_csv(sep="\t", index=True, header=False)),
        outpng,
        genome,
        tmpdir,
        **kwargs
    )
    return str(outpng)
# def render_clip_primary_view(_3dg_file, outpng, view="lr", **kwargs):
#     """
#     Rendering clip primary views of a single-cell structure with pymol, colored by territory.
#     Designed for a none-glomerulus nucleus.
#     Input:
#         _3dg_file: 3dg file of a single-cell structure
#         outpng: output png file
#         view: view to render, ["lr", "ht", "dv"], means looking towards \
#             left-right, head-tail or dorsal-ventral axis.
#     Output:
#         outpng path; and wrote outpng file
#     """
#     _3dg = pymol_primary_views(_3dg_file, view)
#     clip_territory_pymol(
#         StringIO(_3dg.to_csv(sep="\t", index=True, header=False)),
#         outpng,
#         **kwargs
#     )
#     return str(outpng)
def render_clip_primary_view(_3dg_file, outpng, view="lr", targets=None, **kwargs):
    """
    Rendering clip primary views of a single-cell structure with pymol, colored by territory.
    Designed for a none-glomerulus nucleus.
    Input:
        _3dg_file: 3dg file of a single-cell structure
        outpng: output png file
        view: view to render, ["lr", "ht", "dv"], means looking towards \
            left-right, head-tail or dorsal-ventral axis.
    Output:
        outpng path; and wrote outpng file
    """
    _3dg = pymol_primary_views(_3dg_file, view, targets)
    clip_territory_pymol(
        StringIO(_3dg.to_csv(sep="\t", index=True, header=False)),
        outpng,
        **kwargs
    )
    return str(outpng)
def render_clip_b_primary_view(_3dg_file, b_factor, outpng, view="lr", targets=None, **kwargs):
    """
    Rendering clip primary views of a single-cell structure with pymol, colored by territory.
    Designed for a none-glomerulus nucleus.
    Input:
        _3dg_file: 3dg file of a single-cell structure
        b_factor: reference b factor
        outpng: output png file
        view: view to render, ["lr", "ht", "dv"], means looking towards \
            left-right, head-tail or dorsal-ventral axis.
    Output:
        outpng path; and wrote outpng file
    """
    _3dg = pymol_primary_views(_3dg_file, view, targets)
    clip_b_pymol(
        StringIO(_3dg.to_csv(sep="\t", index=True, header=False)),
        b_factor,
        outpng,
        **kwargs
    )
    return str(outpng)
def render_clip_centelo_primary_view(_3dg_file, outpng, genome, clip, slab, tmpdir, view="lr", targets=None, **kwargs):
    """
    Color primary views of a single-cell structure with relative dist to centromeres and telomeres.
    Designed for a none-glomerulus nucleus.
    Input:
        _3dg_file: 3dg file of a single-cell structure
        outpng: output png file
        genome: genome name of the structure
        tmpdir: temporary directory
        view: view to render, ["lr", "ht", "dv"], means looking towards \
            left-right, head-tail or dorsal-ventral axis.
        targets: flip according to this.
    Output:
        outpng path; and wrote outpng file
    """
    _3dg = pymol_primary_views(_3dg_file, view, targets)
    clip_centelo_pymol(
        StringIO(_3dg.to_csv(sep="\t", index=True, header=False)),
        outpng,
        genome=genome,
        clip=clip,
        slab=slab,
        tmpdir=tmpdir,
        **kwargs
    )
    return str(outpng)

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