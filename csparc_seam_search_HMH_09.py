import math
import multiprocessing as mp

import click
import matplotlib
import numpy as np
from cryosparc.dataset import Dataset
from matplotlib import pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.figure import Figure  # Adde this in here
from scipy.spatial.transform import Rotation as R
from tqdm import tqdm




matplotlib.use("Agg")
font = {"size": 10}
matplotlib.rc("font", **font)
plt.rcParams['figure.max_open_warning'] = 5000

def calc_sphere_coordinates(
    x: float, y: float, z: float, z_rot: float, y_rot: float, z2_rot: float
):
    vec = [x, y, z]

    rotation = R.from_euler("zyz", [z_rot, y_rot, z2_rot], degrees=True)
    vec = rotation.apply(vec)

    return vec


def class_posterior(particles):
    return particles["alignments3D_multi/class"][
        np.arange(len(particles)),
        np.argmax(particles["alignments3D_multi/class_posterior"], axis=1),
    ]


def count_values(arr):
    unique, counts = np.unique(arr, return_counts=True)
    class_counts = dict(zip(unique, counts))

    return class_counts


def count_max(counts):
    return max(counts, key=counts.get)


def process_mic(args):
    mic, conf_threshold = args
    new_tubes = {}
    new_pfs = []
    figs = []

    tubes = mic.split_by("filament/filament_uid")

    for tube_id in tubes:
        tube = tubes[tube_id]

        pfn = tube["sym_expand/helix_num_rises"][0]
        if pfn not in new_tubes:
            new_tubes[pfn] = []

        tube_probs = []

        src_uids = np.unique(tube["sym_expand/src_uid"])

        # n = len(src_uids)
        # n_rows = math.floor(math.sqrt(n))
        # n_cols = math.ceil(n / n_rows)

        # fig, ax = plt.subplots(n_rows, n_cols, sharey=True, sharex=True)
        # fig.supxlabel("Protofilament index")
        # fig.supylabel("Register class / seam score")

        # ax_i = 0
        for src_uid in src_uids:
            expanded_group = tube.query({"sym_expand/src_uid": src_uid})
            if len(expanded_group) != pfn:
                continue

            probs = []

            order = np.argsort(expanded_group["sym_expand/idx"])
            classes = np.array(class_posterior(expanded_group))[order]

            # Score where seam is at index 0 (no seam break)
            class_counts = count_values(classes)
            tube_class = count_max(class_counts)
            probs.append(class_counts[tube_class] / pfn)

            for i in range(1, pfn):
                group1 = classes[:i]
                group2 = classes[i:]

                group1_counts = count_values(group1)
                group2_counts = count_values(group2)

                group1_class = count_max(group1_counts)
                group2_class = count_max(group2_counts)

                group1_class_prob = group1_counts.get(group1_class, 0) / len(group1)
                group2_class_prob = group2_counts.get(group2_class, 0) / len(group2)

                if group1_class_prob > group2_class_prob:
                    group2_class = int(not group1_class)
                else:
                    group1_class = int(not group2_class)

                seam_prob = (
                    group1_counts.get(group1_class, 0)
                    + group2_counts.get(group2_class, 0)
                ) / pfn

                probs.append(seam_prob)

            # ax[np.unravel_index(ax_i, (n_rows, n_cols))].bar(range(pfn), probs)
            # ax[np.unravel_index(ax_i, (n_rows, n_cols))].scatter(
            #     range(pfn), classes / 2 + 0.5
            # )

            tube_probs.append(probs)

            # ax_i += 1

        # figs.append(fig)

        fig = plt.figure()
        plt.bar(range(pfn), np.mean(tube_probs, axis=0))
        plt.text(0.5, 0.95, f"Filament UID: {src_uid}", ha='center', va='center', transform=plt.gca().transAxes, fontsize=10)  # Add the UID to the plot title
        figs.append(fig)

        seam_pos = np.argmax(np.mean(tube_probs, axis=0))

        if np.mean(tube_probs, axis=0)[seam_pos] < conf_threshold:
            continue

        new_tubes[pfn].append(tube.query({"sym_expand/idx": seam_pos}))

        src_uids = tube.split_by("sym_expand/src_uid")

        for src_uid in src_uids:
            expanded_group = src_uids[src_uid]
            if len(expanded_group) != pfn:
                continue

            if seam_pos == 0:
                counts = count_values(class_posterior(expanded_group))
                modal_class = count_max(counts)
                expanded_group["alignments3D_multi/class_posterior"] = (
                    expanded_group["alignments3D_multi/class"] == modal_class
                ).astype(int)
                new_pfs.append(expanded_group)
                continue

            group1 = expanded_group.query({"sym_expand/idx": list(range(seam_pos))})
            group2 = expanded_group.query(
                {"sym_expand/idx": list(range(seam_pos, pfn + 1))}
            )

            group1_counts = count_values(class_posterior(group1))
            group2_counts = count_values(class_posterior(group2))

            group1_class = count_max(group1_counts)
            group2_class = count_max(group2_counts)

            group1_class_prob = group1_counts.get(group1_class, 0) / len(group1)
            group2_class_prob = group2_counts.get(group2_class, 0) / len(group2)

            if group1_class_prob > group2_class_prob:
                group2_class = int(not group1_class)
            else:
                group1_class = int(not group2_class)

            group1["alignments3D_multi/class_posterior"] = (
                group1["alignments3D_multi/class"] == group1_class
            ).astype(int)
            group2["alignments3D_multi/class_posterior"] = (
                group2["alignments3D_multi/class"] == group2_class
            ).astype(int)

            new_pfs.append(group1)
            new_pfs.append(group2)

    return new_tubes, new_pfs, figs


@click.command()
@click.option("--i", "input_particles", required=True, help="Input .cs file.")
@click.option(
    "--recenter",
    required=True,
    help="New center reported by CryoSPARC in Volume Aligment Tools when recentering on a single protofilament, check pixel size. Separate centers of different pfns with ;",
    type=str,
)
@click.option(
    "--recenter_init_pxsize",
    required=True,
    help="Pixel size during New center reported by CryoSPARC in Volume Aligment Tools when recentering on a single protofilament. A/px ;",
    type=float,  # Changed from str to float
)
@click.option(
    "--recenter_final_pxsize",
    required=True,
    help="Final Pixel size for recentering. New center will be located at voxel coordinates. A/px ;",
    type=float,  # Changed from str to float
)
@click.option(
    "--o", "output_name", help="Override output filename for new .cs and .csg files."
)
@click.option(
    "--conf",
    "conf_threshold",
    help="Confidence threshold for seam assignment.",
    default=0.5,
    show_default=True,
)
@click.option(
    "--j",
    "num_cpus",
    help="Maximum number of CPU cores to use for multiprocessing.",
    default=mp.cpu_count() // 2,
    show_default=True,
)
def main(input_particles, output_name, conf_threshold, num_cpus, recenter, recenter_init_pxsize, recenter_final_pxsize ):
    tqdm.write("Reading particles...")
    particles = Dataset.load(input_particles)

    mics = particles.split_by("location/micrograph_path")

    if not output_name:
        output_name = f"{input_particles[:-3]}_seamed_{conf_threshold * 100:.0f}"

    pdf = PdfPages(f"{output_name}.pdf")

    new_tubes = {}
    pfs = []
    tqdm.write("Starting seam search...")
    with mp.Pool(min(len(mics), num_cpus)) as pool:
        for new_tube, new_pfs, figs in tqdm(
            pool.imap_unordered(
                process_mic, [(mic, conf_threshold) for mic in mics.values()]
            ),
            total=len(mics),
        ):
            for pfn in new_tube:
                if pfn not in new_tubes:
                    new_tubes[pfn] = []
                new_tubes[pfn].extend(new_tube[pfn])
            pfs.extend(new_pfs)

            for fig in figs:
                pdf.savefig(fig)
                plt.close(fig)

    pdf.close()

    tqdm.write("Writing reassigned particles...")
    pfs_dataset = Dataset.append_many(*pfs)
    pfs_dataset.save(f"{input_particles[:-3]}_assigned.cs")
    with open(input_particles + "g") as f:
        csg = f.read()
    csg = csg.replace(input_particles, f"{input_particles[:-3]}_assigned.cs")
    with open(f"{input_particles[:-3]}_assigned.csg", "w") as f:
        f.write(csg)

    for pfn in new_tubes:
        # tqdm.write(f"Writing {pfn}pf MT...")
        result = Dataset.append_many(*new_tubes[pfn])
        result.save(f"{output_name}_{pfn}pf.cs")
        with open(input_particles + "g") as f:
            csg = f.read()
        csg = csg.replace(input_particles, f"{output_name}_{pfn}pf.cs")
        with open(f"{output_name}_{pfn}pf.csg", "w") as f:
            f.write(csg)
        tqdm.write(f"Written {output_name}_{pfn}pf.csg")

    recenters = recenter.split(";")

    for coords in recenters:
        coords = list(map(float, coords.split(",")))  # Parse coordinates
        coords_final = [coord * recenter_init_pxsize / recenter_final_pxsize for coord in coords]  # Scale
        box_size = particles["blob/shape"][0][0]  # Retrieve box size
        coords_final = (box_size / 2) + (box_size - np.array(coords_final))  # Recenter
        tqdm.write(f"Recenter to {', '.join(map(str, coords_final))} px")  # Display updated coordinates

if __name__ == "__main__":
    main()
