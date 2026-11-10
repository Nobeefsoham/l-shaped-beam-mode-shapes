# L-Shaped Fixed-Free Beam Mode Shapes

MATLAB / GNU Octave script using the transfer matrix method to estimate natural frequencies and plot mode shapes of a two-segment L-shaped fixed-free steel beam.

## Run

Open `l_shaped_beam_mode_shapes.m` in MATLAB or GNU Octave, set the current folder to the script directory, and run:

```matlab
l_shaped_beam_mode_shapes
```

The script prints frequencies in rad/s and Hz and opens a separate figure for each of up to three detected modes.

## Model

Two 0.5 m segments joined at a 90-degree corner. Each segment has Young's modulus 200 GPa, density 7850 kg/m^3, width 20 mm, and thickness 2 mm. The first end is fixed and the final end is free. Cross-sectional properties are A = Width * Thk and I = Width * Thk^3 / 12.

## Numerical interpretation

The frequency sweep is 1:1:1000 rad/s. The script detects determinant sign changes and reports the grid value immediately before each detected change; it does not refine roots. Frequencies and the corresponding shapes are therefore approximate, despite the original comment referring to an exact natural frequency. The plotting routine explicitly handles two segments. Deflections are normalized for visibility and do not represent a forced-response amplitude.

The original script is reproduced as supplied by the repository owner. It was not executed or independently validated as part of publication.

## Citation

Nobeefsoham. (2026). *L-Shaped Fixed-Free Beam Mode Shapes* [Source code]. GitHub. https://github.com/Nobeefsoham/l-shaped-beam-mode-shapes

Record the specific Git commit or cite a permanent file link to identify the version used. `CITATION.cff` provides citation metadata with the GitHub username as the author identifier.
