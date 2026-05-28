# Marmara Tomography Review — Figure Source Code

Reproducible figure scripts for the review paper:

> **Ahmed, S.M.S. & Güneyli, H.** (2026). Seismic Tomography of the Marmara–Istanbul Region: Advances in Crustal Imaging, Fault Behaviour, and Hazard Implications.

---

## Figures

| Script | Description |
|--------|-------------|
| `fig1_tectonic_map.py` | Tectonic setting map with fault coupling segments, historical earthquakes, repeating-earthquake clusters, and tectonic zone boundaries |
| `fig2_3d_structure.py` | 3-D perspective structural view of bathymetry, fault plane, microseismicity, and Moho |
| `fig3_depth_slices.py` | Synthetic P-wave velocity perturbation (δVp/Vp %) maps at 3, 7, 12, and 18 km depth |
| `fig4_cross_sections.py` | Along-strike (W→E) and across-strike (N→S) Vp cross-sections with Moho and fault-zone damage |
| `fig5_workflow_checker.py` | LET inversion workflow flowchart and conceptual checkerboard resolution test |
| `fig6_coupling_vpvs.py` | Three-panel along-strike profile: Vp, Vp/Vs ratio, and fault coupling with repeater density |
| `fig7_migration.py` | Timeline of progressive eastward migration of M > 5 events along the MMF (2011–2025) |

## Quick Start

```bash
# Clone the repository
git clone https://github.com/<your-username>/marmara-tomography-review.git
cd marmara-tomography-review

# Install dependencies
pip install matplotlib numpy scipy

# Generate all figures at once
python generate_all_figures.py

# Or generate a single figure
python fig3_depth_slices.py

# Custom DPI and output directory
python generate_all_figures.py --dpi 600 --outdir high_res_figs/
```

Output PNGs are saved to `output/` by default (created automatically).

## Requirements

| Package | Version |
|---------|---------|
| Python | ≥ 3.8 |
| matplotlib | ≥ 3.5 |
| numpy | ≥ 1.21 |
| scipy | ≥ 1.7 |

## Repository Structure

```
marmara-tomography-review/
├── README.md
├── generate_all_figures.py      # Master script
├── fig1_tectonic_map.py
├── fig2_3d_structure.py
├── fig3_depth_slices.py
├── fig4_cross_sections.py
├── fig5_workflow_checker.py
├── fig6_coupling_vpvs.py
├── fig7_migration.py
└── output/                      # Generated PNGs (gitignored)
```

## Key References

These figures synthesise results from the following published studies:

- Bécel, A. et al. (2009). *Tectonophysics*, 467, 1–21. [doi:10.1016/j.tecto.2008.10.022](https://doi.org/10.1016/j.tecto.2008.10.022)
- Bayrakçı, G. et al. (2013). *Geophysical Journal International*, 194, 1335–1357. [doi:10.1093/gji/ggt211](https://doi.org/10.1093/gji/ggt211)
- Laigle, M. et al. (2008). *Geophysical Journal International*, 174, 1037–1051. [doi:10.1111/j.1365-246X.2008.03835.x](https://doi.org/10.1111/j.1365-246X.2008.03835.x)
- Polat, G. et al. (2016). *Earth, Planets and Space*, 68, 132. [doi:10.1186/s40623-016-0503-4](https://doi.org/10.1186/s40623-016-0503-4)
- Yamamoto, Y. et al. (2017). *Journal of Geophysical Research*, 122, 2069–2084. [doi:10.1002/2016JB013608](https://doi.org/10.1002/2016JB013608)
- Karabulut, H. (2024). *Geophysical Journal International*, 237, 1208–1221. [doi:10.1093/gji/ggae109](https://doi.org/10.1093/gji/ggae109)
- Karabulut, H. (2025). *Geophysical Journal International*, 241, 986–1008. [doi:10.1093/gji/ggaf024](https://doi.org/10.1093/gji/ggaf024)
- Becker, D. et al. (2023). *Geophysical Research Letters*, 50, e2022GL101471. [doi:10.1029/2022GL101471](https://doi.org/10.1029/2022GL101471)
- Martínez-Garzón, P. et al. (2025). *Science*. [doi:10.1126/science.adz0072](https://doi.org/10.1126/science.adz0072)
- Chen, X. et al. (2025). *Geophysical Research Letters*, 52, e2024GL111460. [doi:10.1029/2024GL111460](https://doi.org/10.1029/2024GL111460)

## Notes

- All figures are **synthetic/schematic** representations constructed from published velocity models, fault parameters, and seismicity patterns. They are not reproductions of original data.
- Each script is self-contained and can be run independently.
- Each script exposes a `create_figure(save_path, dpi)` function for programmatic use.

## License

MIT — see [LICENSE](LICENSE) for details.
