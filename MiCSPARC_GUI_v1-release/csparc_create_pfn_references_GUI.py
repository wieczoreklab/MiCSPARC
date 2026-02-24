import multiprocessing as mp

import click
import mrcfile
import numpy as np
import os
from cryosparc.dataset import Dataset
from scipy.ndimage import rotate, shift
from tqdm import tqdm


def gen_helix(params):
    new_pf, twist, rise = params

    rotated = rotate(new_pf, twist, reshape=False, axes=(2, 1))
    rotated_shifted = shift(rotated, (rise, 0, 0))

    return rotated_shifted


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
@click.option(
    "--j",
    "num_cpus",
    help="Maximum number of CPU cores to use for multiprocessing.",
    default=mp.cpu_count() // 2,
    show_default=True,
)
def main(
    input_protofilament,
    num_cpus,
    recenter,
    original_pix,
    initial_pfn,
):
    recenter = [float(x) for x in recenter.split(", ")]

    volume = Dataset.load(input_protofilament)

    tqdm.write("Reading...")
    # Clean and join the relative path with the location of the .cs file
    cs_dir = os.path.dirname(input_protofilament)
    relative_map_path = volume["map/path"][0].lstrip(">/ ").strip()
    map_path = os.path.abspath(os.path.join(cs_dir, relative_map_path))

    # Optional: fail early with a clearer error if the file is still not found
    if not os.path.isfile(map_path):
        raise FileNotFoundError(f"Could not find volume file: {map_path}")
    pf = mrcfile.read(map_path)
    pixel_size = volume["map/psize_A"][0]
    box_size = volume["map/shape"][0]

    recenter = abs(box_size - np.array(recenter)) * original_pix / pixel_size

    parameters = {
        11: 123 / pixel_size,
        12: 123 / pixel_size,
        13: 123 / pixel_size,
        14: 123 / pixel_size,
        15: 164 / pixel_size,
        16: 164 / pixel_size,
    }

    tqdm.write("Padding...")
    new_box = np.pad(pf, int(box_size[0] // 2))

    for pfn in tqdm(parameters):
        twist = 360 / pfn
        rise = parameters[pfn]

        tqdm.write("Recentering...")
        new_pf = shift(new_box, recenter[::-1] * pfn / initial_pfn)

        tqdm.write("Applying twists and rises...")
        positions = list(range(-int(pfn * 2), int(pfn * 2)))
        with mp.Pool(min(len(positions), num_cpus)) as pool:
            rotations = list(
                tqdm(
                    pool.imap_unordered(
                        gen_helix,
                        [(new_pf, twist * i, rise * i / pfn) for i in positions],
                    ),
                    total=len(positions),
                )
            )
            new_pf = shift(new_pf, (82 / pixel_size, 0, 0))
            rotations += list(
                tqdm(
                    pool.imap_unordered(
                        gen_helix,
                        [(new_pf, twist * i, rise * i / pfn) for i in positions],
                    ),
                    total=len(positions),
                )
            )

        tqdm.write("Writing...")
        with mrcfile.new(f"{map_path[:-4]}_{pfn}pf.mrc", overwrite=True) as mrc:
            mrc.set_data(np.max(rotations, axis=0))
            mrc.voxel_size = pixel_size


if __name__ == "__main__":
    main()
