import click
import numpy as np
from cryosparc.dataset import Dataset
from matplotlib import pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from scipy.spatial.transform import Rotation as R
from tqdm import tqdm


@click.command()
@click.option("--i", "input_particles", required=True, help="Input .cs file.")
@click.option(
    "--o", "output_name", help="Override output filename for new .cs and .csg files."
)
def main(input_particles, output_name):
    particles = Dataset.load(input_particles)

    mics = particles.split_by("location/micrograph_path")

    if not output_name:
        output_name = f"{input_particles[:-3]}_upsi"

    pdf = PdfPages(f"{output_name}.pdf")

    unified_tubes = []

    for mic_name in tqdm(mics):
        tubes = mics[mic_name].split_by("filament/filament_uid")

        for filament_uid in tubes:
            tube = tubes[filament_uid]

            if len(tube) < 3:
                continue

            eulers = R.from_rotvec(tube["alignments3D/pose"]).as_euler(
                "ZYZ", degrees=True
            )
            psi = eulers[:, 2]

            if not (abs(np.diff(psi)) > 1).any():
                unified_tubes.append(tube)
                continue

            vals, bins = np.histogram(psi, bins="doane", range=(-180, 180))
            mode_ind = vals.argmax()

            modal_mask = (psi > bins[mode_ind] - 10) & (psi < bins[mode_ind + 1] + 10)
            modal_points = psi[modal_mask]

            if len(modal_points) < 3:
                continue

            fig, ax = plt.subplots(1, 2, sharey=True)
            fig.supxlabel("Particle index along tube")
            fig.supylabel("Psi angle (degrees)")
            plt.ylim(-180, 180)

            ax[0].scatter(range(len(psi)), psi)

            fit = (
                np.polynomial.polynomial.Polynomial.fit(
                    np.argwhere(modal_mask).flatten(),
                    modal_points,
                    2,
                )
                .convert()
                .coef
            )
            corrected_psi = np.polynomial.polynomial.polyval(np.arange(len(tube)), fit)
            eulers[:, 2] = corrected_psi
            poses = R.from_euler("ZYZ", eulers, degrees=True).as_rotvec()
            tube["alignments3D/pose"] = poses

            ax[1].scatter(range(len(corrected_psi)), corrected_psi)

            pdf.savefig(fig)
            plt.close(fig)

            unified_tubes.append(tube)

    pdf.close()

    unified_tubes_cs = Dataset.append_many(*unified_tubes)

    unified_tubes_cs.save(f"{output_name}.cs")

    with open(input_particles + "g") as f:
        csg = f.read()
    csg = csg.replace(input_particles, f"{output_name}.cs")
    with open(f"{output_name}.csg", "w") as f:
        f.write(csg)

    tqdm.write(f"Written {output_name}.csg")


if __name__ == "__main__":
    main()
