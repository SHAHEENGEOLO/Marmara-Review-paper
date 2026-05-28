"""
Figure 7: Progressive Eastward Migration of M > 5 Events Along the MMF
========================================================================
Timeline plot showing the decade-long (2011–2025) eastward migration of
moderate-to-large earthquakes along the Main Marmara Fault, culminating
in the April 2025 MW 6.2 event near the locked segment south of Istanbul.

References:
    - Martínez-Garzón et al. (2025), Science, doi:10.1126/science.adz0072
    - Chen et al. (2025), GRL, doi:10.1029/2024GL111460
    - Becker et al. (2023), GRL, doi:10.1029/2022GL101471

Authors: S.M.S. Ahmed & H. Güneyli
"""

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np


def create_figure(save_path='fig7_migration.png', dpi=300):
    fig, ax = plt.subplots(1, 1, figsize=(7.5, 3.5), dpi=dpi)

    # ---- Event data (approximate, after Martínez-Garzón et al. 2025) ----
    events = [
        # (year, longitude_along_MMF, magnitude, coupling_segment)
        (2011, 27.8, 5.1, 'creeping'),
        (2014, 28.1, 5.1, 'transition'),
        (2016, 28.3, 4.9, 'transition'),
        (2019, 28.6, 5.8, 'transition'),
        (2025, 28.9, 6.2, 'locked'),
    ]

    colours = {
        'creeping':   '#2196F3',
        'transition': '#FF9800',
        'locked':     '#F44336',
    }

    # ---- Background segment bands ----
    ax.axhspan(27.5, 28.0, color='#2196F3', alpha=0.08)
    ax.axhspan(28.0, 28.5, color='#FF9800', alpha=0.08)
    ax.axhspan(28.5, 29.2, color='#F44336', alpha=0.08)

    ax.text(2010, 27.75, 'Creeping',   fontsize=6, color='#1565C0',
            fontstyle='italic')
    ax.text(2010, 28.25, 'Transition', fontsize=6, color='#E65100',
            fontstyle='italic')
    ax.text(2010, 28.85, 'Locked',     fontsize=6, color='#C62828',
            fontstyle='italic')

    # ---- Plot events ----
    for year, lon, mag, seg in events:
        c = colours[seg]
        size = mag**2.5 * 8
        ax.scatter(year, lon, s=size, c=c, edgecolors='k',
                   linewidths=0.5, zorder=5)
        ax.annotate(f'M{mag}', xy=(year, lon),
                    xytext=(year + 0.5, lon + 0.08),
                    fontsize=6, fontweight='bold', color=c)

    # ---- Migration arrow ----
    ax.annotate(
        '', xy=(2026, 29.1), xytext=(2010.5, 27.7),
        arrowprops=dict(arrowstyle='->', color='#333', lw=1.5,
                        linestyle='--',
                        connectionstyle='arc3,rad=0.2'))
    ax.text(2018, 28.6, 'Eastward\nmigration', fontsize=7,
            color='#333', fontstyle='italic', rotation=25, ha='center')

    # ---- Istanbul reference line ----
    ax.axhline(y=29.0, color='#F44336', linestyle=':', linewidth=1,
               alpha=0.5)
    ax.text(2026.5, 29.0, '← Istanbul', fontsize=6.5, color='#C62828',
            va='center', fontweight='bold')

    # ---- Axes ----
    ax.set_xlabel('Year', fontsize=9)
    ax.set_ylabel('Longitude along MMF (°E)', fontsize=9)
    ax.set_xlim(2009.5, 2027)
    ax.set_ylim(27.4, 29.3)
    ax.tick_params(labelsize=7)
    ax.set_title(
        'Progressive Eastward Migration of M > 5 Events Along the MMF',
        fontsize=9, fontweight='bold')
    ax.grid(True, alpha=0.2)

    fig.tight_layout()
    fig.savefig(save_path, dpi=dpi, bbox_inches='tight', facecolor='white')
    plt.close(fig)
    print(f'Saved → {save_path}')


if __name__ == '__main__':
    create_figure()
