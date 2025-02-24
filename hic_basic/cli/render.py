from concurrent.futures import ThreadPoolExecutor
from itertools import repeat
from pathlib import Path
import pandas as pd

name = "render"
description = """Use this command to render multiple structures in parallel and output a png file for each structure.
This script is used to render clip of primary view (colored by CpG) of a none-glomerulus single-cell structure.
Also use this to render all kinds of clip-CpG views, for example round GM12878 cells."""

from ..hicio import read_meta
from ..plot.render import render_clip_b_primary_view, clip_b_pymol, render_surface_primary_view, surface_territory_pymol

def render_task(args):
    """
    Render a single structure.
    """
    sample_name, _3dg_file, b_factor, dupref, cmap, clip, slab, view, target, outdir = args
    outpng = outdir / f"{sample_name}.{view}c{clip}s{slab}.png"
    tmpdir = outdir
    if outpng.exists():
        print(f"Skip {sample_name} {clip} {view}")
        return outpng
    print(f"Processing {sample_name} {clip} {view}")
    #print(b_factor)
    if clip == "surf":
        if view == "rand":
            res = surface_territory_pymol(
                _3dg_file, outpng, outdir, conda=None
            )
        else:
            res = render_surface_primary_view(
                _3dg_file, outpng, outdir, view=view, targets=target, conda=None)
    else:
        if view == "rand":
            res =  clip_b_pymol(
            _3dg_file, b_factor, outpng, dupref=dupref, cmap=cmap, clip=clip, slab=slab, tmpdir=tmpdir, conda=None
            )
        else:
            res = render_clip_b_primary_view(
            _3dg_file, b_factor, outpng, view, target, dupref=dupref, cmap=cmap, clip=clip, slab=slab, tmpdir=tmpdir, conda=None)
    return res

def add_arguments(subparser):
    subparser.add_argument(
        '-i','--sample_table_file',
        type=str,
        required=True,
        help='A .csv sample table, must have `20k_g_struct1` or `gs` column')
    subparser.add_argument('-b','--b_factor', type=str, help='Reference b factor. If multiple B factors are used, separate them by comma and set --multiB')
    subparser.add_argument('--multiB', action='store_true', help='Whether to use multiple B factors. See --b_factor')
    subparser.add_argument('--dupref', action='store_true', help='Whether to duplicate the reference for diploid genome')
    subparser.add_argument('--cmap', type=str, default='', help='Colormap to use. Use -- to separate words. eg: viridis')
    subparser.add_argument(
        '-o','--outdir', 
        type=str,
        required=True,
        help='Directory to store the output png files')
    subparser.add_argument('--filesp_path', type=str, help='Path of the output filesp file')

    # view and clip selection
    subparser.add_argument(
        '--view', 
        choices=['lr', 'dv', 'ht', 'rand'],
        default='rand',
        help='View to render'
    )
    subparser.add_argument(
        '--target',
        choices=["1", "-1", "0"],
        type=str,
        default=None,
        help='Flip according to this after selecting the view.'
    )
    subparser.add_argument('--clips', nargs="*", help='Clips to render. 0 for middle, negative for further, positive for closer. eg: -1 0 1. `surf` for non-clip surface')
    subparser.add_argument('--view_clips', type=str, help='View and clip combinations to render, using - to separate. Overwrite the view and clip arguments. eg: lr:0,1-dv:0')
    subparser.add_argument('--slab', type=int, help='Thickness of the slab')
    subparser.add_argument('--nproc', type=int, default=4, help='Number of processes to use')

    # debug instructions
    subparser.add_argument('--samples', type=str, default=None, help='Do these samples only. Useful for debugging.')
    subparser.add_argument('--targets', type=str, default=None, help='Flip according to this after selecting the view. Require `samples` argument, must be same lengths with samples. Overwrite the target argument. Eg: 1:-1:-1,1:1:1')


def run(args):
    # --- parse arguments --- #
    # get the sample file list
    sample_table = read_meta(args.sample_table_file)
    if "20k_g_struct1" in sample_table.columns:
        gs = sample_table["20k_g_struct1"].dropna()
    elif "gs" in sample_table.columns:
        gs = sample_table["gs"].dropna()
    else:
        raise ValueError("No gs column found in sample table")

    # get the color
    b_factor = args.b_factor
    if args.multiB:
        b_factor = b_factor.split(",")
    dupref = args.dupref
    cmap = args.cmap.replace("--", " ")

    # get the view and clip
    view = args.view
    clips = args.clips
    slab = args.slab

    # prepare the output directory
    outdir = Path(args.outdir)
    if not outdir.exists():
        outdir.mkdir(parents=True)
    filesp_path = args.filesp_path

    # check debug instructions
    # whether to render specific view-clip combinations for each sample
    if args.view_clips:
        view_clips = {
            view_clip.split(":")[0]: [int(i) for i in view_clip.split(":")[1].split(",")]
            for view_clip in args.view_clips.split("--")
        }
    else:
        view_clips = {view: clips}
    if args.samples:
        samples = args.samples.split(",")
        gs = gs.loc[samples]
        if args.multiB:
            assert len(b_factor) == gs.shape[0]
    if args.targets:
        targets = args.targets.split(",")
        targets = [list(map(int, row.split(":"))) for row in targets]
        assert len(targets) == gs.shape[0]
    else:
        targets = repeat(None)

    nproc = args.nproc

    # 使用ThreadPoolExecutor来并行处理样本
    with ThreadPoolExecutor(max_workers=nproc) as executor:
        outpngs = {
            f"{view}c{clip}s{slab}" : list(executor.map(
                render_task,
                zip(gs.index, gs, b_factor if args.multiB else repeat(b_factor),
                    repeat(dupref), repeat(cmap), repeat(clip), repeat(slab), repeat(view), targets, repeat(outdir))
            ))
            for view, clips in view_clips.items()
            for clip in clips
        }