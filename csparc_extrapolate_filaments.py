import math
import multiprocessing as mp
import time
import warnings
import click
import numba as nb
import numpy as np
import yaml
from cryosparc.dataset import Dataset
from scipy.spatial import distance
from tqdm import tqdm

warnings.simplefilter("ignore", np.exceptions.RankWarning)


@nb.njit(fastmath=True)
def taylor(u: float, L: float, lu: float) -> float:
    return u + (L / lu)


@nb.njit(fastmath=True)
def newton(v: float, b: float, a: float, L: float, du: float, lu: float) -> float:
    dv = b + 2 * a * v
    lv = math.sqrt(1 + dv**2)
    return v + (
        (lv * (4 * a * L + du * lu - dv * lv + math.asinh(du) - math.asinh(dv)))
        / (2 * a * (2 + 2 * dv**2))
    )


# Recursively run newton until it converges on an answer
@nb.njit(fastmath=True)
def do_approximation(
    v: float, b: float, a: float, L: float, du: float, lu: float, its: int
) -> float:
    sig_figs = 100000

    last_newton = newton(v, b, a, L, du, lu)
    if abs(v - last_newton) <= 1 / sig_figs or its > 10:
        return math.floor(round(v * sig_figs)) / sig_figs
    else:
        return do_approximation(last_newton, b, a, L, du, lu, its + 1)


@nb.njit(fastmath=True)
def parabolic_arc_solver(a: float, b: float, u: float, L: int) -> float:
    # Adapted from https://gist.github.com/sikanrong/bd7b05b800a5086c1502e2c7033127ed

    du = b + 2 * a * u
    lu = math.sqrt(1 + du**2)

    # get a good first guess for newton from taylor
    first_guess = taylor(u, L, lu)

    return do_approximation(first_guess, b, a, L, du, lu, 0)


@nb.njit(fastmath=True)
def psi_prior(a: float, b: float, x: float):
    return -math.degrees(math.atan(b + (2 * a * x)))


@nb.njit(fastmath=True)
def generate_coords(fit, start_x, xlim, ylim, step):
    c, b, a = fit
    y = np.polynomial.polynomial.polyval(start_x, fit)

    corners = (
        int(math.sqrt((start_x**2) + (y**2)) // step),
        int(math.sqrt(((xlim - start_x) ** 2) + (y**2)) // step),
        int(math.sqrt((start_x**2) + ((ylim - y) ** 2)) // step),
        int(math.sqrt(((xlim - start_x) ** 2) + ((ylim - y) ** 2)) // step),
    )
    max_n = max(corners)

    extrapolated = []
    extrapolated_reverse = []

    for i in range(max_n):
        new_x = parabolic_arc_solver(a, b, start_x, i * step)
        if new_x > xlim or new_x < 0:
            break
        new_y = np.polynomial.polynomial.polyval(new_x, fit)
        if new_y > ylim or new_y < 0:
            break
        extrapolated.append([new_x, new_y])

    for i in range(-1, -max_n, -1):
        new_x = parabolic_arc_solver(a, b, start_x, i * step)
        if new_x > xlim or new_x < 0:
            break
        new_y = np.polynomial.polynomial.polyval(new_x, fit)
        if new_y > ylim or new_y < 0:
            break
        extrapolated_reverse.append([new_x, new_y])

    return np.array(extrapolated_reverse[::-1] + extrapolated)


def create_dataset(mic, new_coords, fil_num, xlim, ylim):
    num_coords = len(new_coords)

    dataset_copy = {}
    for key in [
        "location/micrograph_uid",
        "location/exp_group_id",
        "location/micrograph_path",
        "location/micrograph_shape",
        "location/micrograph_psize_A",
        "filament/inter_box_dist_A",
        "filament/filament_uid",
        "location/center_x_frac",
        "location/center_y_frac",
        "filament/position_A",
    ]:
        dataset_copy[key] = [mic[key][0]] * num_coords

    new_dataset = Dataset(dataset_copy)
    new_dataset["filament/filament_uid"] = [fil_num] * num_coords
    new_dataset["location/center_x_frac"] = new_coords[:, 0] / xlim
    new_dataset["location/center_y_frac"] = new_coords[:, 1] / ylim
    new_dataset["filament/position_A"] = (
        np.arange(num_coords) * mic["filament/inter_box_dist_A"][0]
    )

    return new_dataset


def calc_coords(dataset, pix_size, xlim, ylim):
    return np.column_stack(
        (
            (dataset["location/center_x_frac"] * xlim)
            + (dataset["alignments2D/shift"][:, 0] / pix_size),
            (dataset["location/center_y_frac"] * ylim)
            + (dataset["alignments2D/shift"][:, 1] / pix_size),
        )
    )


def extrapolate_filaments(mic):
    pix_size = mic["location/micrograph_psize_A"][0]
    box_step = mic["filament/inter_box_dist_A"][0] / pix_size
    xlim = mic["location/micrograph_shape"][0][1]
    ylim = mic["location/micrograph_shape"][0][0]

    extrapolated_tubes = []

    # Split filaments if their aligned psi is not smooth
    new_tubes = []
    new_tube_id = 1
    tubes = mic.split_by("filament/filament_uid")
    for fil_uid in tubes:
        tube = tubes[fil_uid]

        psi = np.degrees(tube["alignments2D/pose"])
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
    fil_nums = sorted(tubes, key=lambda x: len(tubes[x]), reverse=True)
    mic_coords = calc_coords(mic, pix_size, xlim, ylim)

    i = 0
    for fil_num in fil_nums:
        i += 1
        tube = mic.query({"filament/filament_uid": fil_num})

        if len(tube) < 3:
            # Probably too short to fit confidently
            continue

        coords = calc_coords(tube, pix_size, xlim, ylim)
        dist_matrix = distance.cdist(coords, coords)
        dist_order = dist_matrix.argsort(axis=None)[::-1]

        for points_index in dist_order:
            furthest_1, furthest_2 = sorted(
                np.unravel_index(points_index, dist_matrix.shape)
            )
            segment = coords[furthest_1 : furthest_2 + 1]
            fit = np.polynomial.polynomial.polyfit(
                segment[:, 0],
                segment[:, 1],
                1,
            )
            fit = np.append(fit, 0.0000000001)

            if np.count_nonzero(
                np.isclose(
                    segment[:, 1],
                    np.polynomial.polynomial.polyval(segment[:, 0], fit),
                    rtol=0,
                    atol=50 * pix_size,
                )
            ) > 0.9 * len(segment):
                new_coords = generate_coords(
                    fit, np.median(segment[:, 0]), xlim, ylim, box_step
                )
                new_dataset = create_dataset(mic, new_coords, fil_num, xlim, ylim)
                extrapolated_tubes.append(new_dataset)

                break

            else:
                fit = np.polynomial.polynomial.polyfit(segment[:, 0], segment[:, 1], 2)
                new_coords = generate_coords(
                    fit, np.median(segment[:, 0]), xlim, ylim, box_step
                )
                x_comp = segment[:, 1]
                y_comp = np.polynomial.polynomial.polyval(segment[:, 0], fit)

                if abs(fit[1]) > 10:
                    fit = np.polynomial.polynomial.polyfit(
                        segment[:, 1], segment[:, 0], 2
                    )
                    new_coords = generate_coords(
                        fit, np.median(segment[:, 1]), ylim, xlim, box_step
                    )
                    x_comp = segment[:, 0]
                    y_comp = np.polynomial.polynomial.polyval(segment[:, 1], fit)
                    new_coords = new_coords[:, ::-1]

                if np.count_nonzero(
                    np.isclose(x_comp, y_comp, rtol=0, atol=50 * pix_size)
                ) > 0.9 * len(segment):
                    new_dataset = create_dataset(mic, new_coords, fil_num, xlim, ylim)
                    extrapolated_tubes.append(new_dataset)

                    break

        else:
            continue

        dist_matrix = distance.cdist(mic_coords, new_coords)
        close_particles = np.any(dist_matrix < 150 * pix_size, axis=1) & ~np.isin(
            mic["filament/filament_uid"], fil_nums[:i]
        )
        far_picks = np.all(dist_matrix > 150 * pix_size, axis=1) & (
            mic["filament/filament_uid"] == fil_num
        )

        fil_uids = mic["filament/filament_uid"]
        fil_uids[close_particles] = fil_num
        mic["filament/filament_uid"] = fil_uids

        if np.any(far_picks):
            overflow_fil_num = max(fil_nums) + 1

            fil_uids = mic["filament/filament_uid"]
            fil_uids[far_picks] = overflow_fil_num
            mic["filament/filament_uid"] = fil_uids

            fil_nums.append(overflow_fil_num)

    return Dataset.append_many(*extrapolated_tubes)


@click.command()
@click.option(
    "--i",
    "input_particles",
    required=True,
    help="Input .cs file.",
    type=click.Path(exists=True, dir_okay=False),
)
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
def main(input_particles, output_name, num_cpus):
    particles = Dataset.load(input_particles)
    mics = particles.split_by("location/micrograph_path")

    # Initial JIT compilation before multiprocessing
    start = time.time()
    psi_prior(1, 1, 1)
    small_mic = max(mics.values(), key=len).slice(10)
    extrapolate_filaments(small_mic)
    tqdm.write(f"Finished numba pre-compiling ({time.time() - start:.4f}s)")

    tqdm.write("Processing remaining micrographs")

    with mp.Pool(min(len(mics), num_cpus)) as pool:
        extrapolated = list(
            tqdm(pool.imap_unordered(extrapolate_filaments, mics.values()), total=len(mics))
        )

    tqdm.write("Writing final results group")

    result = Dataset.append_many(*extrapolated)

    if not output_name:
        output_name = f"{input_particles[:-3]}_extrapolated"

    result.save(f"{output_name}.cs")
    with open(input_particles + "g") as f:
        csg = yaml.safe_load(f.read())
    new_results_group = {
        "location": csg["results"]["location"],
        "filament": csg["results"]["filament"],
    }
    new_results_group["location"]["metafile"] = new_results_group["location"][
        "metafile"
    ].replace(input_particles, f"{output_name}.cs")
    new_results_group["filament"]["metafile"] = new_results_group["filament"][
        "metafile"
    ].replace(input_particles, f"{output_name}.cs")
    new_csg = {
        "created": csg["created"],
        "group": csg["group"],
        "results": new_results_group,
        "version": csg["version"],
    }
    with open(f"{output_name}.csg", "w") as f:
        yaml.dump(new_csg, f)

    tqdm.write(f"Written {output_name}.csg")


if __name__ == "__main__":
    main()
