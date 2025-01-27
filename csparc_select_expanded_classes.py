import sys

import numpy as np
from cryosparc.dataset import Dataset
from tqdm import tqdm

input_particles_star = sys.argv[1]
pfn = int(sys.argv[2])
kept_classes = set(int(x) for x in sys.argv[3].split(","))
particles = Dataset.load(input_particles_star)

src_uids = np.unique(particles["sym_expand/src_uid"])

selected_tubes = []

for src_uid in tqdm(src_uids):
    expanded_group = particles.query({"sym_expand/src_uid": src_uid})
    if len(expanded_group) != pfn:
        continue
    classes = set(
        np.argmax(expanded_group["alignments3D_multi/class_posterior"], axis=1)
    )
    if not classes - kept_classes:
        selected_tubes.append(src_uid)

selected_particles = particles.query({"sym_expand/src_uid": selected_tubes})

selected_particles.save(f"{input_particles_star[:-3]}_selected.cs")

with open(input_particles_star + "g") as f:
    csg = f.read()
csg = csg.replace(input_particles_star, f"{input_particles_star[:-3]}_selected.cs")
with open(f"{input_particles_star[:-3]}_selected.csg", "w") as f:
    f.write(csg)

tqdm.write(f"Written {output_name}.csg")
