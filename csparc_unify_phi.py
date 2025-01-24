import multiprocessing as mp

import click
import numpy as np
from cryosparc.dataset import Dataset
from matplotlib import pyplot as plt
from scipy.optimize import differential_evolution
from scipy.spatial.transform import Rotation as R
from tqdm import tqdm


def wrap_to_180(n):
    return ((n + 180) % 360) - 180


def find_modal_rot(coeffs, rots, pfn):
    x = np.arange(len(rots))
    y = rots
    twist = 360 / pfn

    close = []
    for i in range(pfn):
        m = coeffs[0]
        c = wrap_to_180(coeffs[1] + i * twist)

        close.append(
            np.count_nonzero(
                np.isclose(
                    y,
                    wrap_to_180(
                        np.polynomial.polynomial.polyval(
                            list(range(min(x), max(x) + 1)), [c, m]
                        )
                    ),
                    rtol=0,
                    atol=twist / 4,
                )
            )
        )

    return [m, wrap_to_180(coeffs[1] + np.argmax(close) * twist)]


def fit_rots(coeffs, rots, pfn):
    x = np.arange(len(rots))
    y = rots
    twist = 360 / pfn

    dists = np.empty((pfn, len(x)))
    for i in range(pfn):
        m = coeffs[0]
        c = wrap_to_180(coeffs[1] + i * twist)

        d = np.abs(
            y
            - wrap_to_180(
                np.polynomial.polynomial.polyval(
                    list(range(min(x), max(x) + 1)), [c, m]
                )
            )
        )

        dists[i] = d

    return dists.min(axis=0).sum()


def plot_rots(coeffs, rots, pfn):
    plt.scatter(np.arange(len(rots)), rots)
    plt.ylim(-180, 180)
    twist = 360 / pfn

    p = np.arange(len(rots))

    for i in range(pfn):
        m = coeffs[0]
        c = wrap_to_180(coeffs[1] + (i * twist))
        plt.plot(p, wrap_to_180(np.polynomial.polynomial.polyval(p, [c, m])))

    plt.show()


def process_mic(args):
    mic, pfn = args

    # Split filaments if their aligned psi is not smooth
    new_tubes = []
    new_tube_id = 1
    tubes = mic.split_by("filament/filament_uid")
    for filament_uid in tubes:
        tube = tubes[filament_uid]

        eulers = R.from_rotvec(tube["alignments3D/pose"]).as_euler("ZYZ", degrees=True)
        psi = eulers[:, 2]
        psi_diff = np.abs(np.diff(psi))

        new_tube_groups = (
            [0]
            + list(
                np.nonzero(np.min([psi_diff, np.abs(180 - psi_diff)], axis=0) > 10)[0]
            )
            + [len(tube)]
        )
        for n in range(len(new_tube_groups) - 1):
            new_tube = tube.take(range(new_tube_groups[n], new_tube_groups[n + 1]))
            new_tube["filament/filament_uid"] = new_tube_id
            new_tubes.append(new_tube)
            new_tube_id += 1

    mic = Dataset.append_many(*new_tubes)
    tubes = mic.split_by("filament/filament_uid")

    unified_tubes = []

    for filament_uid in tubes:
        tube = tubes[filament_uid]

        if len(tube) < 5:
            # Probably too short to fit confidently
            continue

        eulers = R.from_rotvec(tube["alignments3D/pose"]).as_euler("ZYZ", degrees=True)
        rots = eulers[:, 0]

        res = differential_evolution(
            fit_rots,
            x0=[0, rots[0] % (360 / pfn)],
            args=(rots, pfn),
            bounds=([-3, 3], [-(360 / pfn), (360 / pfn)]),
        )
        result = res.x
        result[1] = wrap_to_180(result[1])

        result = find_modal_rot(result, rots, pfn)

        corrected_rot = wrap_to_180(
            np.polynomial.polynomial.polyval(np.arange(len(rots)), result[::-1])
        )

        eulers[:, 0] = corrected_rot
        poses = R.from_euler("ZYZ", eulers, degrees=True).as_rotvec()
        tube["alignments3D/pose"] = poses

        unified_tubes.append(tube)

    return Dataset.append_many(*unified_tubes)


@click.command()
@click.option(
    "--i",
    "input_particles",
    required=True,
    help="Input .cs file.",
    type=click.Path(exists=True, dir_okay=False),
)
@click.option("--pfn", required=True, help="Protofilament number of the MT.", type=int)
@click.option(
    "--o",
    "output_name",
    help="Override output filename for new .cs and .csg files.",
    type=str,
)
@click.option(
    "--j",
    "num_cpus",
    help="Maximum number of CPU cores to use for multiprocessing.",
    default=mp.cpu_count() // 2,
    show_default=True,
)
def main(input_particles, output_name, pfn, num_cpus):
    particles = Dataset.load(input_particles)

    mics = particles.split_by("location/micrograph_path")

    unified_tubes = []

    with mp.Pool(min(len(mics), num_cpus)) as pool:
        unified_tubes = list(
            tqdm(
                pool.imap_unordered(
                    process_mic, [(mic, pfn) for mic in mics.values()], chunksize=1
                ),
                total=len(mics),
            )
        )

    unified_tubes_cs = Dataset.append_many(*unified_tubes)

    if not output_name:
        output_name = f"{input_particles[:-3]}_uphi"

    unified_tubes_cs.save(f"{output_name}.cs")

    with open(input_particles + "g") as f:
        csg = f.read()
    csg = csg.replace(input_particles, f"{output_name}.cs")
    with open(f"{output_name}.csg", "w") as f:
        f.write(csg)

    tqdm.write(f"Written {output_name}.csg")


if __name__ == "__main__":
    main()
