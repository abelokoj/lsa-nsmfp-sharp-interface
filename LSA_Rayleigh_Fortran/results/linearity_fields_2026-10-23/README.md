# Field-direction linearity test (2026-10-23)

At k R₀ = 0.7, the twin-difference fields of three amplitudes, ε = 2.5e-3, 5e-3 and 1e-2, were compared after scaling each by its ε. The comparison was made at matched times, up to the time at which the largest amplitude reaches a/R₀ = 0.1. Three measures were used, each against the smallest amplitude:

- the relative difference of the scaled fields;
- the angle between them;
- the ratio of their norms.

The cases were run on three meshes, each with its own twin: 8 cells (`amp08_*`), 16 cells (`sw070_16`, `lin050_16`, `lin100_16`) and 32 cells (`amp32_*`). The 16-cell runs `lin050_16` and `lin100_16` were re-run on 22–23 October to regenerate their plotfiles. Their amplitude records are identical to the October originals.

Produced on the cluster by `scripts/linearity_fields.sh`, which calls `plotting/linearity_fields.py`. The band-2 results use `--band 2`.

Largest values over the window, given as ε = 5e-3 / ε = 1e-2:

| Mesh | Observable | Relative difference | Angle | Norm ratio |
|---|---|---|---|---|
| 8 | velocity | 3.96 / 8.52 % | 2.23 / 4.84° | 0.88 / 1.80 % |
| 8 | level set, band 2Δx | 16.04 / 17.01 % | 9.11 / 9.65° | 1.28 / 1.44 % |
| 8 | level set, whole domain | 51.6 / 75.3 % | 22.3 / 32.1° | 49.2 / 73.3 % |
| 16 | velocity | 4.13 / 8.60 % | 2.33 / 4.90° | 0.85 / 1.67 % |
| 16 | level set, band 2Δx | 5.73 / 6.10 % | 3.28 / 3.49° | 0.16 / 0.19 % |
| 16 | level set, whole domain | 58.3 / 77.5 % | 32.2 / 46.1° | 48.7 / 73.1 % |
| 32 | velocity | 3.58 / 7.22 % | 2.03 / 4.13° | 0.61 / 1.04 % |
| 32 | level set, band 2Δx | 2.03 / 2.61 % | 1.16 / 1.50° | 0.15 / 0.20 % |
| 32 | level set, whole domain | 57.4 / 74.6 % | 34.8 / 46.2° | 48.6 / 72.7 % |

Reading:

- **Velocity and narrow-band level set.** The scaled fields agree in size within 1.8 % (velocity) and 1.4 % (band level set; 0.2 % at 16 and 32 cells), and in direction within 4.9° (velocity) and 9.7°, 3.5° and 1.5° (band level set at 8, 16 and 32 cells). The largest relative differences that remain are 17.0, 6.1 and 2.6 % for the band level set, and 7.2–8.6 % for the velocity at ε = 1e-2.
- **How the residuals vary.** The level-set residual decreases with mesh refinement and changes little between ε = 5e-3 and 1e-2 (16.0 against 17.0 % at 8 cells). This is consistent with discretization error of the interface reconstruction, but that has not been tested separately. The velocity residual roughly doubles when ε doubles (3.6–4.1 % against 7.2–8.6 %) and changes little with the mesh, which is consistent with a contribution quadratic in ε near the end of the window. The residuals have not been decomposed further.
- **Whole-domain level set: not proportional to ε, already at t = 0.** The largest norm mismatches are 48.6–49.2 % (ε = 5e-3) and 72.7–73.3 % (ε = 1e-2); at 16 cells they are already 47 % and 72 % at t = 0. The solver stores the level set only within about 0.3 R₀ of the interface and holds it at a fixed value beyond. In the `meas16` level-set export at t = 0, the cells on the edge of that region differ by up to 0.14 between the perturbed and the base run, about thirty times ε, while inside the region the difference is ε sin(kz). A part of the difference that does not scale with ε gives a norm mismatch of 1 − ε_ref/ε, 50 % and 75 %, close to the measured values; the edge cells were inspected only in the `meas16` export. The narrow band used by the level-set decomposition lies inside the stored region, which is why the observable is restricted to it.

Files: `lin<mesh>_export.*` (velocity), `lin<mesh>_export_ls.*` (level set, whole domain) and `lin<mesh>_export_ls_band2.*` (level set, band). The `.txt` files hold the per-time tables and the `.csv` files the same tables.
