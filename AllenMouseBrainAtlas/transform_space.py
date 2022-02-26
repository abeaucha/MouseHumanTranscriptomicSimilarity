#!/usr/bin/env python3

# %% Imports

import argparse
import atexit
import numpy as np
import os
import sys
import tempfile
from pyminc.volumes.factory import *
from subprocess import call

# %% Argument parsing

parser = argparse.ArgumentParser(description='Orient Allen Institute CCFv3 imaging data to standard spaces')
parser.add_argument('infile', type=str, help='Input volume')
parser.add_argument('outfile', type=str, help='Output volume')
parser.add_argument('--tmpdir', '-t', default='/tmp', type=str, help='Temporary directory')
parser.add_argument('--voxel_orientation', '-v', default='RAS', type=str, help='Voxel orientation; "RAS" (default) or "PIR"')
parser.add_argument('--world_space', '-w', default='MICe', type=str, help='World space; "MICe" (default), "CCFv3", or "stereotaxic"')
parser.add_argument('--expansion_factor', '-x', default=1.0, type=float, help='Factor to artificially expand volume by (e.g. to represent world coordinates in um rather than mm, set this to 1000')
parser.add_argument('--volume_type', default=None, type=str, help='Volume type (default is from the input file)')
parser.add_argument('--data_type', default=None, type=str, help='Data type  (default is from the input file)')
parser.add_argument('--labels', default=False, action='store_true', help='Labels? (default is from the input file)')
parser.add_argument('--clobber', default=False, action='store_true', help='Overwrite output file? (Default: false)')

args = parser.parse_args()

# %% Testing arguments
"""
class Args:
    infile = "/projects/yyee/tools/Allen_brain/data/mouse_ccf/average_template/average_template_50.nrrd"
    outfile = "/tmp/ts_test.mnc"
    tmpdir = "/tmp" 
    voxel_orientation = "RAS"
    world_space = "CCFv3" 
    expansion_factor = 10
    volume_type = None
    data_type = None
    labels = None


args = Args()
"""

# %% Test if outfile exists

if os.path.exists(args.outfile) and not args.clobber:
    print("Output file exists! Use -clobber to overwrite")
    sys.exit(1)

# %% Function definitions


def is_minc(infile):
    if os.path.splitext(infile)[1] == ".mnc":
        return(True)
    else:
        return(False)


def itk_convert(infile, outfile):
    call(["itk_convert", "--clobber", infile, outfile])


def remove_temp_files(file_list):
    for f in file_list:
        if os.path.exists(f.name):
            f.close()


def reorient_to_standard(dat):
    dat = np.rot90(dat, k=1, axes=(0, 2))
    dat = np.rot90(dat, k=1, axes=(0, 1))

    shape = dat.shape
    dat = np.ravel(dat)
    dat = np.reshape(dat, shape)

    return(dat)


def do_nothing(dat):
    return(dat)

# %% Coordinate definitions

# Centers are listed as x,y,z; reverse these when writing out
# Centers listed in um in CCFv3 coordinates
centers_RAS = {"MICe"   :   [5700, 7900, 2700],
               "CCFv3"  :   [0, 13200, 8000]}
centers_PIR = {"MICe"   :   [5300, 5300, 5700],
               "CCFv3"  :   [0, 0, 0]}

# Direction cosines
direction_cosines_RAS = {"MICe"     :   [[1, 0, 0], [0, 1, 0], [0, 0, 1]],
                         "CCFv3"    :   [[0, 0, 1], [-1, 0, 0], [0, -1, 0]]}
    
direction_cosines_PIR = {"MICe"     :   [[0, -1, 0], [0, 0, -1], [1, 0, 0]],
                         "CCFv3"    :   [[1, 0, 0], [0, 1, 0], [0, 0, 1]]}

# Map arguments to functions/dicts/values
map_voxel_orientations = {"RAS" :   reorient_to_standard,
                          "PIR":    do_nothing}

map_centers = {"RAS"    :   centers_RAS,
               "PIR"    :   centers_PIR}

map_dir_cosines = {"RAS"    :   direction_cosines_RAS,
                   "PIR"    :   direction_cosines_PIR}

size_10 = 1320*800*1140
size_25 = 528*320*456
size_50 = 264*160*228
size_100 = 132*80*114
size_200 = 58*41*67
map_resolutions = {size_10: 10,
                   size_25: 25,
                   size_50: 50,
                   size_100: 100,
                   size_200: 200}
# %% Preprocessing

temporary_files = []

if not is_minc(args.infile):
    tmp = tempfile.NamedTemporaryFile(dir=args.tmpdir,
                                      prefix='mousetools-tmpconv-',
                                      suffix='.mnc')
    infile = tmp.name
    itk_convert(args.infile, infile)
    temporary_files.append(tmp)
else:
    infile = args.infile


# %% Main script

# Read volume
vol = volumeFromFile(infile)

# Detect resolution
res = map_resolutions[vol.data.size]

# Voxel orientation
if args.voxel_orientation in map_voxel_orientations:
    new_data = map_voxel_orientations[args.voxel_orientation](vol.data)
else:
    print("Invalid voxel orientation")
    sys.exit(1)

# World coordinate system
centers = [args.expansion_factor*c/(1000)
           for c in map_centers[args.voxel_orientation][args.world_space]]
steps = [args.expansion_factor*res/1000] * 3
xdc = map_dir_cosines[args.voxel_orientation][args.world_space][0]
ydc = map_dir_cosines[args.voxel_orientation][args.world_space][1]
zdc = map_dir_cosines[args.voxel_orientation][args.world_space][2]

# Types
vtype = vol.volumeType if args.volume_type is None else args.volume_type
dtype = vol.dtype if args.data_type is None else args.data_type
labels = vol.labels if args.labels is None else args.labels
# %% Data output

if not is_minc(args.outfile):
    tmp = tempfile.NamedTemporaryFile(dir=args.tmpdir,
                                      prefix='mousetools-tmpconv-',
                                      suffix='.mnc')
    outfile = tmp.name
    temporary_files.append(tmp)
else:
    outfile = args.outfile

outvol = volumeFromDescription(outputFilename=outfile,
                               dimnames=["zspace", "yspace", "xspace"],
                               sizes=new_data.shape,
                               starts=[-c for c in reversed(centers)],
                               steps=[s for s in reversed(steps)],
                               x_dir_cosines=xdc,
                               y_dir_cosines=ydc,
                               z_dir_cosines=zdc,
                               volumeType=vtype,
                               dtype=dtype,
                               labels=labels)
outvol.data = new_data
outvol.writeFile()
outvol.closeVolume()

# %% Postprocessing

if not is_minc(args.outfile):
    itk_convert(outfile, args.outfile)

# %% On exit

atexit.register(remove_temp_files, file_list=temporary_files)

# %% TODO

# Set attribute in minc file

