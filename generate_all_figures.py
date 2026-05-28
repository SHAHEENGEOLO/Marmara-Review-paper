#!/usr/bin/env python3
"""
generate_all_figures.py
========================
Master script to generate all seven figures for:

  Ahmed & Güneyli — "Seismic Tomography of the Marmara–Istanbul Region:
  Advances in Crustal Imaging, Fault Behaviour, and Hazard Implications"

Requirements:
    pip install matplotlib numpy scipy

Usage:
    python generate_all_figures.py          # default 300 dpi, saves to ./output/
    python generate_all_figures.py --dpi 600 --outdir figures/
"""

import argparse
import os
import sys

# ── import each figure module ──────────────────────────────────────────
import fig1_tectonic_map
import fig2_3d_structure
import fig3_depth_slices
import fig4_cross_sections
import fig5_workflow_checker
import fig6_coupling_vpvs
import fig7_migration


FIGURES = [
    ('fig1_tectonic_map',    fig1_tectonic_map),
    ('fig2_3d_structure',    fig2_3d_structure),
    ('fig3_depth_slices',    fig3_depth_slices),
    ('fig4_cross_sections',  fig4_cross_sections),
    ('fig5_workflow_checker', fig5_workflow_checker),
    ('fig6_coupling_vpvs',   fig6_coupling_vpvs),
    ('fig7_migration',       fig7_migration),
]


def main():
    parser = argparse.ArgumentParser(
        description='Generate all figures for the Marmara Tomography review.')
    parser.add_argument('--dpi', type=int, default=300,
                        help='Output resolution (default: 300)')
    parser.add_argument('--outdir', type=str, default='output',
                        help='Directory for saved PNGs (default: output/)')
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    print(f'\n{"=" * 60}')
    print(f'  Generating 7 figures  |  DPI = {args.dpi}  |  → {args.outdir}/')
    print(f'{"=" * 60}\n')

    for name, module in FIGURES:
        path = os.path.join(args.outdir, f'{name}.png')
        try:
            module.create_figure(save_path=path, dpi=args.dpi)
        except Exception as exc:
            print(f'  ✗ {name}: {exc}', file=sys.stderr)

    print(f'\n  Done — {len(FIGURES)} figures saved to {args.outdir}/\n')


if __name__ == '__main__':
    main()
