🧊 MiCSPARC: A GUI-Driven Microtubule Processing Pipeline
MiCSPARC is a microtubule (MT) processing pipeline built around the CryoSPARC platform to facilitate high-resolution structure determination of both decorated and undecorated microtubules.

Microtubules are polar filaments formed by αβ-tubulin heterodimers. Their structural arrangement often disrupts helical symmetry due to a discontinuity known as the seam, complicating their analysis via single-particle analysis (SPA). Moreover, cryo-EM images of microtubules tend to be noisy and pseudo-helical, which makes it difficult to distinguish between nearly identical tubulin subunits.

MiCSPARC introduces a user-friendly, GUI-based pipeline that integrates with CryoSPARC to automate key tasks including filament tracing, protofilament number assignment, register correction, and seam localization. This greatly accelerates cryo-EM processing workflows and enhances the reproducibility of structural studies. MiCSPARC has been validated on diverse datasets and offers a reliable, accessible solution for structural biologists.

Install Instructions

1. Set up a Conda environment using the provided 'environment.yml' file.                                                                                                                                                                     $ conda env create -f environment_GUI.yml
2. Before running MiCSPARC, ensure all dependencies are installed and that 'cryosparc-tools' package matches your CryoSPARC version (major and minor release). For example, if using CryoSPARC v4.7.1, install the latest v4.7.x version of cryosparc-tools. You can edit environment_GUI.yml to match your version.                                                                                                                                                                                                               - cryosparc-tools==4.7.1
3. Activate the environment and launch the MiCSPARC GUI:                                                                                                                                                                                             $ conda activate micsparc-GUI
4. Edit 'cs_config.yml'to ensure that contains your CryoSPARC credentials: port, email, host, license, password
5. Some may use test_connection.py to check taht everything works fine                                                                                                                                                                              $ python test_connection.py
6. Launch the GUI                                                                                                                                                                                                                                    $ python micsparc_gui_launcher_v1-release.py 

PROTOCOL
1.	Motion correction, CTF estimation, curate exposures
Generate initial templates for filament tracer
2.	Filament tracer — ∼300 Å*, 82 Å separation, 300–400 Å template-free diameter
If the filament tracer is not picking any microtubules at this stage, increasing the sd of gaussian blur (e.g. to 2.0) and/or increasing both hysteresis values in increments of 1 can help.
3.	Inspect picks, extract — ∼550–600 Å box*, bin to 4–5 Å/px
4.	2D classification
5.	Select classes with single clear tube
•	If tube edges are close to circular mask (esp. with large decorators), expand box size in later extractions

Generate initial picks
6.	Filament tracer — Templates from step 5, ∼300 Å*, 82 Å separation
7.	Inspect picks, extract
8.	2D classification — 50–100 classes depending on number of picks, 2 final iterations
9.	2D classification 
50 classes, 2 final iterations, disable sigma annealing (start annealing sigma at iteration 200)
10.	2D classification — 1 class → export particles group
Extrapolate initial picks
The filament tracer will generally not pick along the entire microtubule, and may introduce errors in filament grouping. This can be rectified using the filament extrapolation script.
11.	Run MiCSPARC GUI – Extrapolate Filaments
Enter the exported directory (it is easiest to run the scripts directly within the exported directory in order to be able to reimport the results group)
$ conda activate micsparc
python /path/to/python micsparc_gui_launcher_v1-beta.py 
 
12.	Import MiCSPARC Result into CryoSPARC
Ensure that the script directory contains a valid cs_config.yml file with the correct CryoSPARC instance credentials. Also, make sure the installed version of cryosparc-tools matches your CryoSPARC major and minor version. For example, if your CryoSPARC version is v4.7.1, use the latest v4.7.x release of cryosparc-tools. Then, browse for your MiCSPARC result file (*.csg) and specify the CryoSPARC project and workspace where you'd like to import the result group.Cs_conig.yml
base_port: 61000
email: hugo.munoz@mol.biol.ethz.ch
host: localhost
license: asdhkdjjdskjsjdgjkljgljdkljbkfdjlk
password: "*********"

 
13.	Extract — bin to ∼2 Å/px
14.	2D classification — 50–100 classes depending on number of picks, 2 final iterations
15.	2D classification — 50 classes, 2 final iterations, disable sigma annealing

Protofilament number sorting
Generate references (potentially optional if suitable references already exist, but highly recommended)
References are required for the classification in 3D into groups of different protofilament numbers. This classification tends to be very dependent on references with the correct factors such as decoration state, positioning, and heterogeneity, and potentially even expansion or compression of the tube. Thus, it is recommended to create synthetic references directly from the dataset, as follows:
16.	Helical refinements with theoretical helical parameters corresponding to potential protofilament numbers in the dataset
Theoretical helical parameters are calculated in the form:
rise= n-start × tubulin spacing / protofilament no. 
twist= 360◦ + supertwist / protofilament no.
Supertwist can generally be taken as 0 at this stage. Common theoretical microtubule parameters:
Protofilament (pfn)	Rise (Å)	Twist (°)
11-3	11.2	−32.7
12-3	10.3	−30.0
13-3	9.46	−27.7
14-3	8.79	−25.7
15-4	10.9	−24.0
16-4	10.3	−22.5

17.	Symmetry expand the best helical refinement by refined helical parameters
The “best” refinement is not necessarily that with the highest resolution, due to the high symmetry applied, but rather the one with the visually best-resolved protofilaments (no smearing around the tube, ideally distinct tubulin dimers)
18.	Recenter on a protofilament (volume alignment tools with recentering to mask CoM) Creating a protofilament mask with ChimeraX segmentation:
•	Volume tools lowpass 15 Å, download map, open in ChimeraX, set an appropriate masking threshold
•	Tools > Volume Data > Segment Map
•	Select a single protofilament, group
•	File (in the Segment window) > Save selected regions to .mrc file (this is temporary, it does not matter where they are saved)
•	cmd volume resample [new volume] onGrid #1
•	Save the resampled volume, import into CryoSPARC
•	Volume tools lowpass 15 Å, output as mask, initial threshold of full map
19.	Downsample — 1/2 box size
20.	Reconstruct, local refinement — default sd rotations and shifts, allow recentering
21.	Inverse mask around single protofilament (volume tools 1/2 box and invert mask from step 18)
22.	Particle subtraction
23.	Reconstruct, local refinement — default sd rotations and shifts, allow recentering, enforce non-negativity → export volume group
24.	Run Create References in MiCSPARC GUI 
$ python /path/to/micsparc_GUI
The --recenter coords are those reported by CryoSPARC in the volume alignment step 18 “New center will be located at voxel coordinates:”. 
 
Protofilament number assignment
25.	Import the microtubule reference volumes of 11, 12, 13, 14, 15, and 16 protofilaments into cryoSPARC.
26.	Heterogeneous refinement — force hard classification → export particles group
27.	Run Assign PFNs in MiCSPARC GUI 
$ python /path/to/micsparc_GUI.py 
 

Analys results:

 
 
28.	Import MiCSPARC Result into CryoSPARC
 
29.	Heterogeneous reconstruction
30.	Split volume groups

Rough MT alignment
Having each particle along the microtubule in the same orientation is important for later symmetry expansion and seam searching steps, as it allows us to assign each group of expansions as the same protofilament along the microtubule.
For each good class:
31.	Helical refinement of particles with theoretical helical parameters — minimise per- particle scale, use NU refinement
31.1.	Particles may need to be unbinned if they hit binned nyquist
32.	CTF refinement
33.	Helical refinement w/ no helical symmetry — minimise per-particle scale, use NU refinement
→ export particles group
34.	Run psi unification script
$ python /path/to/micsparc_GUI.py  
Note that the Results of the Psi Unification can be analised in the Analyse Result Tab  
 
35.	Import MiCSPARC Result into CryoSPARC
 
36.	Local refinement — C1 helical volume and mask, 5sd rotations, 3sd shifts, allow rotation recentering → export particles group
37.	Run phi unification script
$ python /path/to/micsparc_GUI.py
 
38.	Import MiCSPARC Result into CryoSPARC
 
Check the Results 
 
 
39.	Local refinement — C1 helical volume and mask, 5sd rotations, 3sd shifts, allow rotation recentering, refine_gs_resplit→ export particles group
40.	Run phi unification again (Phi 2 times in total)
$  python /path/to/micsparc_GUI.py
41.	Import new results group
42.	Local refinement — C1 helical volume and mask, 3sd rotations, 2sd shifts, don’t allow recentering
43.	Symmetry expansion — refined helical parameters, order = pfn
44.	Local refinement — 3sd rotations, 2sd shifts, don’t allow recentering
•	Undecorated: If this does not reach <4 Å resolution, later steps may prove difficult/im- possible. Decorated particles will be more lenient depending on the size of the decorator.

Protofilament alignment
45.	Create mask around single protofilament — if refining multiple pfns, it can be helpful to select the protofilament that coincides best across all pfns
With large decorators, select the protofilament with the most mixed population of decorator registers to ensure no gaps in the mask (or use a high dilation). Selecting the “worst” protofilament will also give the best split for the register classification, which will generally also allow seam search to work better.
46.	Volume alignment tools — recenter to mask CoM
47.	Volume tools — crop realigned mask to 1/2 box size, invert
48.	Extract realigned particles — 1/2 box size, recenter using aligned shifts
49.	Reconstruct, local refinement — 3sd rotations, 2sd shifts, don’t allow recentering
50.	Particle subtraction — inverted mask
51.	Reconstruct, local refinement — 3sd rotations, 2sd shifts, allow recentering, enforce non- negativity
If refining multiple classes:
52.	Pick a base refinement, align other pfns with align 3D maps (may need to be roughly aligned manually with volume alignment and ChimeraX)
53.	Local refinement — combine all pfns, 3sd rotations, 2sd shifts, allow rotation recentering, enforce non-negativity

Register correction
54.	3D classification — 10 classes, no per-particle scale, filter 4 Å for undecorated, up to 12 Å for large decorators
Lower resolution limits (closer to 12 than 4) will allow the classification to work better for decorators.
55.	Identify good references
For undecorated: at least one class with distinct α-tubulin S9-S10 loop, ideally at least two which exhibit a difference in register
For decorated: at least one class with clear and convincing spacing of decorator, ideally at least two which exhibit a difference in register (if applicable)
56.	Reconstruct, local refinement for each register, combining particles from appropriate classes — 3sd rotations, 2sd shifts, allow rotation recentering, enforce non-negativity
57.	If only one register was observed, create the other register using volume alignment
— recenter in Z by 41 Å
58.	3D classification w/ both refined registers — 2 classes, input initialization mode, filter 4 Å for undecorated, up to 12 Å for large decorators, force hard classification, no per-particle scale This step should result in a ∼50/50 split of particles
59.	Local refinement both classes separately — 3sd rotations, 2sd shifts, allow recentering, enforce non-negativity
60.	Volume alignment — Shift worse class by 41 Å in Z
61.	Align 3D both maps w/ particles from shifted class
62.	Local refine both classes together — 3sd rotations, 2sd shifts, allow recentering, enforce non-negativity
63.	CTF refinement
64.	Local refinement — 3sd rotations, 2sd shifts, allow recentering, enforce non-negativity For a final highest-resolution average of a protofilament:
65.	Duplicate removal
66.	Reference-based motion refinement
67.	Local refinement — 3sd rotations, 2sd shifts, allow recentering, enforce non-negativity

MT reconstruction
68.	Heterogeneous reconstruction using pre-duplicate removal particles (step 64) and alignments3D_multi field of the 2-class 3D classification → export particles group
69.	Run seam search script
$ python /path/to/csparc_seam_search.py ...
Using the recenter coords from step 46 *Check pixel size of this step and final pixel size of 48.
Take note of the recenter coordinates
 
70.	Import results group/s (if multiple pfns)
(e.g. /path/to/CS-project/exports/groups/JX_particles/JX_particles_seamed_50_13pf.csg)
71.	Volume alignment — recenter on center coords provided by step 69
72.	Extract — return to original MT box size
73.	Local refinement — volume and mask of symmetrised helical refinement, 3sd rotations, 2sd shifts, don’t allow recentering
74.	CTF refinement
75.	Local refinement — 3sd rotations, 2sd shifts, don’t allow recentering
76.	Duplicate removal
77.	3D classification has being proved beneficial for unveiling a clear decoration of every PT.
78.	Reference-based motion refinement
79.	Local refinement — 3sd rotations, 2sd shifts, don’t allow recentering


Tips

* increase size as appropriate for large decorations  
sd - standard deviation
pf -protofilament
pfn - protofilament number w/ - with

All scripts can be used with help function to show required and optional parameters

python /path/to/script.py --help    

For all reconstructions that are followed by local refinements, use the particles, volume, and mask generated by the reconstruction as input for the refinement.

All local refinements may benefit from having “FSC noise substitution” enabled, but this has not been vigorously tested.

