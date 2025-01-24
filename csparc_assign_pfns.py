import click
import numpy as np
from cryosparc.dataset import Dataset
from scipy import stats
from tqdm import tqdm


def class_posterior(particles):
    return particles["alignments3D_multi/class"][
        np.arange(len(particles)),
        np.argmax(particles["alignments3D_multi/class_posterior"], axis=1),
    ]


@click.command()
@click.option("--i", "input_particles", required=True, help="Input .cs file.")
@click.option(
    "--conf",
    "conf_threshold",
    help="Confidence threhold for pf number assignment.",
    default=0.7,
    show_default=True,
)
@click.option(
    "--o", "output_name", help="Override output filename for new .cs and .csg files."
)
def main(input_particles, output_name, conf_threshold):
    particles = Dataset.load(input_particles)
    mics = particles.split_by("location/micrograph_path")

    new_tubes = []

    for mic_name in tqdm(mics):
        mic = mics[mic_name]

        fil_uids = np.unique(mic["filament/filament_uid"])
        new_fil_uids = fil_uids.tolist()

        for fil_uid in fil_uids:
            tube = mic.query({"filament/filament_uid": fil_uid})
            classes = class_posterior(tube)

            smoothened = []
            for i in range(0, len(tube)):
                window = classes[
                    np.clip(i - 3, 0, None) : np.clip(i + 4, None, len(tube))
                ]
                mode, count = stats.mode(window)
                smoothened.append(mode)

            changes = np.diff(smoothened)
            split_pos = (
                [0] + (np.argwhere(changes) + 1).flatten().tolist() + [len(tube)]
            )

            for i in range(len(split_pos) - 1):
                split = tube.take(range(split_pos[i], split_pos[i + 1]))
                modal_class = smoothened[split_pos[i]]

                if (
                    np.count_nonzero(
                        classes[split_pos[i] : split_pos[i + 1]] == modal_class
                    )
                    / len(split)
                    < conf_threshold
                ):
                    continue

                split["alignments3D_multi/class_posterior"] = (
                    split["alignments3D_multi/class"] == modal_class
                ).astype(int)

                new_fil_uid = max(new_fil_uids) + 1
                split["filament/filament_uid"] = new_fil_uid

                new_tubes.append(split)

    new_tubes_cs = Dataset.append_many(*new_tubes)

    if not output_name:
        output_name = f"{input_particles[:-3]}_assigned"

    new_tubes_cs.save(f"{output_name}.cs")

    with open(input_particles + "g") as f:
        csg = f.read()
    csg = csg.replace(input_particles, f"{output_name}.cs")
    with open(f"{output_name}.csg", "w") as f:
        f.write(csg)

    tqdm.write(f"Written {output_name}.csg")


if __name__ == "__main__":
    main()
