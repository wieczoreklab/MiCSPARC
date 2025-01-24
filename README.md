# MiCSPARC
Microtubule Image Processing in CSPARC Pipeline

* increase size for large decorations
sd - standard deviation
pfn - protofilament number
w/ - with
All scripts can be used with python /path/to/script.py --help to show required and optional
parameters
Create a new conda environment with the provided environment.yml

Picking
1. Motion correction, CTF estimation, curate exposures
Generate initial templates for filament tracer
2. Filament tracer — ∼300Å*, 82Å separation, 300-400Å template-free diameter
3. Inspect picks, extract — ∼550-600Å box*, bin to 4-5Å/px
4. 2D classification
5. Pick classes with single clear tube
• If tube edges are close to circular mask (esp. with large decorators), expand box size in later
extractions
Generate seeds for filament extrapolation
6. Filament tracer — Templates from step 5, ∼300Å*, 82Å separation
7. Inspect picks, extract
8. 2D classification — 50-100 classes depending on number of picks, 2 final iterations
9. 2D classification — 50 classes, 2 final iterations, disable sigma annealing (start annealing sigma
at iteration 200)
10. 2D classification — 1 class → export particles
Final particle picks
11. Run filament extrapolation script
12. Import new results group
13. Extract — bin to ∼2Å/px
14. 2D classification — 50-100 classes depending on number of picks, 2 final iterations
15. 2D classification — 50 classes, 2 final iterations, disable sigma annealing

Protofilament number sorting
Generate references (potentially optional, but recommended)
16. Helical refinements with theoretical helical parameters
17. Symmetry expand best helical refinement with refined helical parameters
18. Recenter on top-most protofilament (volume alignment tools with recentering to mask CoM)
19. Downsample — 1/2 box size
20. Reconstruct, local refinement w/ recentering — default sd rotations and shifts, allow recentering
21. Inverse mask around single protofilament (volume tools 1/2 box and invert on mask from step18)
22. Particle subtraction
23. Reconstruct, local refinement w/ recentering — default sd rotations and shifts, allow recentering,
enforce non-negativity → export
24. Run reference generation script
Protofilament number assignment
25. Import reference pf volumes—if using external references: it is important to match the decoration
state of the sample, reference with different decorator is untested
26. Heterogeneous refinement — force hard classification → export
27. Run pf number assignment script
28. Import results group
29. Heterogeneous reconstruction, split volume groups

MT rough alignment
For each good class:
30. Helical refinement of particles with theoretical helical parameters — minimise per-particle scale
30.1. Particles may need to be unbinned if they hit binned nyquist
31. CTF refinement
32. Helical refinement w/ no helical symmetry → export
33. Run psi unification script
34. Import results group
35. Local refinement — C1 helical volume and mask, 5sd rotations, 3sd shifts, allow recentering →
export
36. Run phi unification script
37. Import results group
38. Local refinement — C1 helical volume and mask, 5sd rotations, 3sd shifts, allow recentering →
export
39. Run phi unification script
40. Import results group
41. Local refinement — C1 helical volume and mask, 5sd rotations, 3sd shifts, don’t allow recentering
42. Symmetry expansion — refined helical parameters, order = pfn
43. Local refinement — 5sd rotations, 3sd shifts, don’t allow recentering
• Undecorated: If this does not reach <4Å resolution, later steps may prove difficult/impossible.
Decorated is more lenient.

Protofilament alignment
44. Create mask around single protofilament — lowpass 15Å; if refining multiple pfns, can be helpful
to select the protofilament that coincides best across all pfns
45. Volume alignment tools — recenter to mask CoM
46. Volume tools — crop realigned mask to ∼1/2 box size, invert
47. Extract realigned particles — ∼1/2 box size, recenter using aligned shifts
48. Reconstruct, local refinement — 5sd rotations, 3sd shifts, don’t allow recentering
49. Particle subtraction — inverted mask
50. Reconstruct, local refinement — 5sd rotations, 3sd shifts, allow recentering
If refining multiple classes:
51. Pick a base refinement, align other pfns with align 3D maps (may need to be done manually with
volume alignment and ChimeraX)
52. Local refinement — all pfns’ particles, 5sd rotations, 3sd shifts, allow recentering
Register correction
53. 3D classification — 10 classes, no per-particle scale, filter 4Å
53.5. If applicable, remove unwanted classes using class selection script
54. Identify at least one class with distinct α-tubulin S9-S10 loop, ideally at least two which exhibit
a difference in register
55. Reconstruct, local refinement for each register — 5sd rotations, 3sd shifts, allow recentering,
enforce non-negativity
56. If only one register observed, create the other register using ChimeraX and volume alignment —
recenter in Z by 41Å
57. 3D classification w/ both refined registers — 2 classes, input initialization mode, filter 4Å (or up
to 12Å depending on decorator size), force hard classification, no per-particle scale
58. Local refinement both classes—5sd rotations, 3sd shifts, allow recentering, enforce non-negativity
59. Volume alignment — Shift worse class by 41Å in Z
60. Align 3D both maps w/ particles from shifted class
61. Local refine both classes together — 5sd rotations, 3sd shifts, allow recentering, enforce nonnegativity
62. CTF refinement
63. Polishing
64. Local refinement — 5sd rotations, 3sd shifts, allow recentering, enforce non-negativity

MT reconstruction
65. Heterogeneous reconstruction using particles of last local refinement and alignments3D_multi
field of 2-class 3D classification → export particles
66. Run seam search script
67. Import results group(s if multiple pfns)
68. Volume alignment — recenter on given center coords
69. Extract — original MT box size
70. Local refinement — volume and mask of symmetrised helical refinement, 5sd rotations, 3sd shifts,
don’t allow recentering
71. CTF refinement
72. Polishing
73. Local refinement — 5sd rotations, 3sd shifts, don’t allow recentering
