import click
import numpy as np
from cryosparc.dataset import Dataset
from tqdm import tqdm


@click.command()
@click.option(
    "--i",
    "input_particles",
    required=True,
    help="Input .cs file.",
    type=click.Path(exists=True, dir_okay=False),
)
@click.option(
    "--keep",
    "kept_classes",
    required=True,
    help="Comma-delimited classes to keep.",
    type=str,
)
@click.option(
    "--o",
    "output_name",
    help="Override output filename for new .cs and .csg files.",
    type=str,
)
def main(input_particles, output_name, kept_classes):
    kept_classes = set(int(x) for x in kept_classes.split(","))
    particles = Dataset.load(input_particles)

    src_uids = np.unique(particles["sym_expand/src_uid"])

    selected_tubes = []

    for src_uid in tqdm(src_uids):
        expanded_group = particles.query({"sym_expand/src_uid": src_uid})
        pfn = expanded_group["sym_expand/helix_num_rises"][0]

        if len(expanded_group) != pfn:
            continue

        classes = set(
            np.argmax(expanded_group["alignments3D_multi/class_posterior"], axis=1)
        )
        if not classes - kept_classes:
            selected_tubes.append(src_uid)

    selected_particles = particles.query({"sym_expand/src_uid": selected_tubes})

    if not output_name:
        output_name = f"{input_particles[:-3]}_selected"

    selected_particles.save(f"{output_name}.cs")

    with open(input_particles + "g") as f:
        csg = f.read()
    csg = csg.replace(input_particles, f"{output_name}.cs")
    with open(f"{output_name}.csg", "w") as f:
        f.write(csg)

    tqdm.write(f"Written {output_name}.csg")


if __name__ == "__main__":
    main()
