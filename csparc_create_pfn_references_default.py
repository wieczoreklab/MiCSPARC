#!/usr/bin/env python3

import multiprocessing as mp
from multiprocessing import Pool
import gc
import mrcfile
import numpy as np
from scipy.ndimage import shift, rotate
from tqdm import tqdm
import click
from cryosparc.dataset import Dataset

def gen_helix(params):
    new_pf, twist, rise = params
    rotated = rotate(new_pf, twist, reshape=False, axes=(2, 1))
    rotated_shifted = shift(rotated, (rise, 0, 0))
    return rotated_shifted

def process_pfn(args):
    pfn, new_box, recenter, initial_pfn, map_path, pixel_size, parameters = args
    twist = 360 / pfn
    rise = parameters[pfn]
    shifted = shift(new_box, np.array(-recenter[::-1]) * pfn / initial_pfn, order=1)
    positions = list(range(-int(2*pfn), int(2*pfn)))
    summed = np.zeros_like(shifted)
    summed += shifted
    for i in positions:
        rotated = gen_helix((shifted, twist * i, rise * i / pfn))
        summed = np.maximum(summed, rotated)
        del rotated
    with mrcfile.new(f"{map_path[:-4]}_{pfn}pf.mrc", overwrite=True) as mrc:
        mrc.set_data(summed)
        mrc.voxel_size = pixel_size
    del summed
    gc.collect()
    return 

@click.command()
@click.option(
    "--i", "input_protofilament", required=True, help="Input volume .cs file."
)
@click.option(
    "--recenter",
    required=True,
    help="New center reported by CryoSPARC when recentering on a single protofilament.",
    type=str,
)
@click.option(
    "--apix",
    "original_pix",
    required=True,
    help="Pixel size during recentering on a single protofilament.",
    type=float,
)
@click.option(
    "--pfn",
    "initial_pfn",
    required=True,
    help="Protofilament number of the original MT.",
    type=float,
)


def main(input_protofilament,
    recenter,
    original_pix,
    initial_pfn,
    ):
    recenter = [float(x) for x in recenter.split(", ")]
    volume = Dataset.load(input_protofilament)
    tqdm.write("Reading...")
    map_path = volume["map/path"][0][1:]
    pf = mrcfile.read(map_path)
    pixel_size = volume["map/psize_A"][0]
    box_size = volume["map/shape"][0]
   
    
    # Load data
    try:
        with mrcfile.mmap(map_path, mode='r') as mrc:
            pf = mrc.data.copy()
    except:
        pf = mrcfile.read(map_path)
    
    box_size = pf.shape[0]
    pad_size = int(box_size // 2)


    # Calculate parameters
    recenter = (box_size - np.array(recenter)) * original_pix / pixel_size
    
    parameters = {
        11: 123 / pixel_size,
        12: 123 / pixel_size,
        13: 123 / pixel_size,
        14: 123 / pixel_size,
        15: 164 / pixel_size,
        16: 164 / pixel_size,
    }
    
    # Padding
    
    new_box = np.pad(pf, pad_size, mode='constant', constant_values=0)
    
    # Run multiprocessing
    with Pool(6) as pool:
        args = [(pfn, new_box, recenter, initial_pfn, map_path, pixel_size, parameters) 
                for pfn in parameters]
        for _ in tqdm(pool.imap_unordered(process_pfn, args), 
                      total=len(parameters),
                      desc="Processing protofilaments"):
            pass

if __name__ == '__main__':
    main()
