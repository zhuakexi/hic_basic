from concurrent.futures import ThreadPoolExecutor
from itertools import repeat
from pathlib import Path
import pandas as pd

name = "render"
description = """Use this command to render multiple structures in parallel and output a png file for each structure."""

from ..hicio import read_meta
from ..plot.render import (
    render_clip_b_primary_view,
    clip_b_pymol,
    render_surface_primary_view,
    surface_territory_pymol,
    surface_b_pymol,
    surface_centelo_pymol,
    render_clip_centelo_primary_view,
    clip_centelo_pymol
)
def render_task(args):
    """
    Render a single structure.
    """
    sample_name, _3dg_file, ref, dupref, cmap, clip, slab, view, target, \
        outdir, mode, force, turn_cmd = args
    outpng = outdir / f"{sample_name}.{view}c{clip}s{slab}.png"
    tmpdir = outdir
    if outpng.exists() and not force:
        print(f"Skip {sample_name} {clip} {view}")
        return outpng
    print(f"Processing {sample_name} {clip} {view}")
    if clip == "surf":
        # surface
        if view in ["rand", "custom"]:
            # random/custom-surface view
            if mode == "territory":
                res = surface_territory_pymol(
                    _3dg_file, outpng, outdir, turn_cmd=turn_cmd, conda=None
                )
            elif mode == "b_factor":
                res = surface_b_pymol(
                    _3dg_file, ref, outpng, dupref=dupref, cmap=cmap, tmpdir=tmpdir, turn_cmd=turn_cmd, conda=None
                )
            elif mode == "centelo":
                res = surface_centelo_pymol(
                    _3dg_file, outpng, genome=ref, dupref=dupref, tmpdir=tmpdir, turn_cmd=turn_cmd, conda=None
                )
            else:
                raise ValueError("Invalid mode")
        else:
            # primary-surface views
            # will ignore the turn_cmd
            if mode == "territory":
                res = render_surface_primary_view(
                    _3dg_file, outpng, outdir, view=view, targets=target, conda=None)
            elif mode == "b_factor":
                raise NotImplementedError("Specific-outside view not supported for b_factor mode")
            elif mode == "centelo":
                raise NotImplementedError("Specific-outside view not supported for centelo mode")
            else:
                raise ValueError("Invalid mode")
    else:
        # clip
        if view in ["rand","custom"]:
            # random/custom-clip view
            if mode == "territory":
                raise NotImplementedError("Random-inside view not supported for territory mode")
            elif mode == "b_factor":
                res = clip_b_pymol(
                _3dg_file, ref, outpng, dupref=dupref, cmap=cmap, clip=clip,
                slab=slab, tmpdir=tmpdir, turn_cmd=turn_cmd, conda=None
                )
            elif mode == "centelo":
                res = clip_centelo_pymol(
                    _3dg_file, outpng, genome=ref, dupref=dupref, clip=clip,
                    slab=slab, tmpdir=tmpdir, turn_cmd=turn_cmd, conda=None
                )
            else:
                raise ValueError("Invalid mode")
        else:
            # primary-clip views
            # will ignore the turn_cmd
            if mode == "territory":
                raise NotImplementedError("Specific-inside view not supported for territory mode")
            elif mode == "b_factor":
                res = render_clip_b_primary_view(
                    _3dg_file, ref, outpng, view, target, dupref=dupref, cmap=cmap, clip=clip, slab=slab, tmpdir=tmpdir, conda=None)
            elif mode == "centelo":
                res = render_clip_centelo_primary_view(
                    _3dg_file, outpng, ref, clip, slab, outdir, view=view, targets=target, conda=None, dupref=dupref
                )
            else:
                raise ValueError("Invalid mode")
    return res

def add_arguments(subparser):
    subparser.add_argument(
        '-i','--sample_table_file',
        type=str,
        required=True,
        help='A .csv sample table, must have `20k_g_struct1` or `gs` column')
    subparser.add_argument("--colname", type=str, default='20k_g_struct1', help='Column name of the structure file')
    subparser.add_argument('--cmap', type=str, default='', help='Colormap to use. Use -- to separate words. eg: viridis')
    subparser.add_argument(
        '-o','--outdir', 
        type=str,
        required=True,
        help='Directory to store the output png files')
    subparser.add_argument('--filesp_path', type=str, help='Path of the output filesp file')
    subparser.add_argument('--force', action='store_true', help='Whether to overwrite existing files')

    # mode_switcher
    subparser.add_argument(
        '--mode',
        choices=['b_factor', 'centelo', 'territory'],
        default='b_factor',
        help='Mode to render:\nb_factor: color according to a reference b factor file. Need to provide --b_factor.\ncentelo: render the structure with centrometer and telomere. Need to provide --genome'
    )

    # b_factor mode
    subparser.add_argument('-b','--b_factor', type=str, help='Reference b factor. If multiple B factors are used, separate them by comma and set --multiB')
    subparser.add_argument('--multiB', action='store_true', help='Whether to use multiple B factors. See --b_factor')
    subparser.add_argument('--dupref', action='store_true', help="""Whether to duplicate the reference for diploid genome.
    For haploid sample, never use this argument.
    For diploid sample, if you use haploid reference(b_factor that does not have chromosome name like chr1(mat)..., or hapliod genome (mm10 rather than mm10_dip) in centelo mode), you should set this argument.
    """)
    
    # centelo mode
    subparser.add_argument('--genome', type=str, help='Genome name. Use this to decide the centromere and telomere positions')
    # view and clip selection
    subparser.add_argument(
        '--view', 
        type=str,
        default='rand',
        help="""View to render. lr, dv, ht are for non-spherical samples like sperms. rand will use default pymol view. Otherwise, view will be treated as column name in the sample table.
        Under this condition, the column should have pymol command strings like `turn x, 90; turn y, 90; zoom 0.5`. If the column is not found, the default pymol view will be used."""
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
    if args.colname in sample_table.columns:
        gs = sample_table[args.colname].dropna()
    else:
        raise ValueError("Column name not found in sample table")

    # get the color
    if args.mode == "centelo":
        genome = args.genome
    elif args.mode == "b_factor":
        b_factor = args.b_factor
        if args.multiB:
            b_factor = b_factor.split(",")
    elif args.mode == "territory":
        pass
    else:
        raise ValueError("Invalid mode")

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
    elif view in ["lr", "dv", "ht", "rand"]:
        view_clips = {view: clips}
    else:
        # view is not given here, use value in the sample table
        view_clips = {"custom": clips}
    if args.samples:
        samples = args.samples.split(",")
        gs = gs.loc[samples]
        if args.multiB and args.mode == "b_factor":
            assert len(b_factor) == gs.shape[0]
    if args.targets:
        targets = args.targets.split(",")
        targets = [list(map(int, row.split(":"))) for row in targets]
        assert len(targets) == gs.shape[0]
    else:
        targets = repeat(None)

    # variables for parallel processing
    nproc = args.nproc
    if args.mode == "centelo":
        ref_list = repeat(genome)
    elif args.mode == "b_factor":
        ref_list = b_factor if args.multiB else repeat(b_factor)
    elif args.mode == "territory":
        ref_list = repeat(None)
    else:
        raise ValueError("Invalid mode")
    if args.view not in ["lr", "dv", "ht", "rand"]:
        # view is given in the sample table
        # use the view column name to infer the identity of the view in output png file name
        view_str_dict = dict.fromkeys(view_clips.keys(), args.view)
        turn_cmds = sample_table.loc[gs.index, args.view].tolist()
    else:
        # echo the view argument in the output png file name
        view_str_dict = {view: view for view in view_clips.keys()}
        turn_cmds = repeat("")

    # 使用ThreadPoolExecutor来并行处理样本
    with ThreadPoolExecutor(max_workers=nproc) as executor:
        outpngs = {
            f"{view_str_dict[view]}c{clip}s{slab}" : list(executor.map(
                render_task,
                zip(gs.index, gs, ref_list,
                    repeat(dupref), repeat(cmap), repeat(clip), repeat(slab), repeat(view),
                    targets, repeat(outdir), repeat(args.mode), repeat(args.force), turn_cmds)
            ))
            for view, clips in view_clips.items()
            for clip in clips
        }