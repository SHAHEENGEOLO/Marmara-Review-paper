"""
Figure 1: Tectonic Setting Map of the Marmara Sea Region
=========================================================
Displays the Main Marmara Fault (MMF) colour-coded by coupling status
(creeping/transition/locked), historical earthquakes, repeating earthquake
clusters, basin locations, and tectonic zone boundaries.

References:
    - Becker et al. (2023), GRL, doi:10.1029/2022GL101471
    - Martínez-Garzón et al. (2025), Science, doi:10.1126/science.adz0072
    - Yamamoto et al. (2017), JGR, doi:10.1002/2016JB013608

Authors: S.M.S. Ahmed & H. Güneyli
"""

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np


def create_figure(save_path='fig1_tectonic_map.png', dpi=300):
    fig, ax = plt.subplots(1, 1, figsize=(7.5, 5.5), dpi=dpi)

    # ----- Sea of Marmara outline (simplified polygon) -----
    sea_x = [26.6, 27.0, 27.4, 28.0, 28.6, 29.2, 29.4,
             29.3, 29.0, 28.5, 28.0, 27.5, 27.0, 26.6]
    sea_y = [40.55, 40.45, 40.35, 40.35, 40.38, 40.5, 40.6,
             40.72, 40.75, 40.78, 40.75, 40.72, 40.65, 40.55]
    ax.fill(sea_x, sea_y, color='#B8D4E3', alpha=0.6, zorder=1)
    ax.text(28.0, 40.55, 'Sea of Marmara', fontsize=8,
            fontstyle='italic', color='#2C5F7C', ha='center', zorder=5)

    # Land background colour
    ax.set_facecolor('#F5F0E8')

    # ----- Sub-basin markers -----
    basins = [
        ('Tekirdağ\nBasin', 27.1, 40.48),
        ('Central\nBasin',  27.85, 40.45),
        ('Çınarcık\nBasin', 29.0, 40.55),
    ]
    for name, bx, by in basins:
        ax.plot(bx, by, 'v', color='#4A7FA5', markersize=6, zorder=4)
        ax.text(bx, by - 0.06, name, fontsize=5.5,
                ha='center', color='#4A7FA5', zorder=5)

    # ----- Main Marmara Fault segments -----
    # Western segment (creeping)
    wx = np.array([26.7, 27.0, 27.4, 27.7])
    wy = np.array([40.58, 40.52, 40.47, 40.44])
    ax.plot(wx, wy, '-', color='#2196F3', linewidth=3.5,
            zorder=6, label='Creeping segment')

    # Central segment (transition)
    cx = np.array([27.7, 28.0, 28.3, 28.6])
    cy = np.array([40.44, 40.42, 40.42, 40.44])
    ax.plot(cx, cy, '-', color='#FF9800', linewidth=3.5,
            zorder=6, label='Transition zone')

    # Eastern segment (locked)
    ex = np.array([28.6, 28.9, 29.15, 29.3])
    ey = np.array([40.44, 40.48, 40.55, 40.62])
    ax.plot(ex, ey, '-', color='#F44336', linewidth=3.5,
            zorder=6, label='Locked segment')

    # Northern and southern secondary branches
    ax.plot([28.8, 29.0, 29.1], [40.65, 40.7, 40.72],
            '--', color='#888', linewidth=1.5, zorder=3)
    ax.plot([28.6, 29.0, 29.3], [40.38, 40.35, 40.38],
            '--', color='#888', linewidth=1.5, zorder=3)

    # Ganos Fault extension
    ax.plot([26.5, 26.7], [40.65, 40.58], '-', color='#888',
            linewidth=2, zorder=3)
    ax.text(26.5, 40.68, 'Ganos\nFault', fontsize=5,
            color='#666', ha='center')

    # ----- Cities -----
    cities = [('Istanbul', 28.97, 41.01),
              ('Bursa', 29.0, 40.18),
              ('Tekirdağ', 27.5, 40.98)]
    for name, cx_, cy_ in cities:
        ax.plot(cx_, cy_, 's', color='#333', markersize=5, zorder=7)
        ax.text(cx_ + 0.05, cy_ + 0.03, name, fontsize=7,
                fontweight='bold', zorder=7)

    # ----- Historical / recent earthquakes -----
    hist_eq = [
        (27.2, 40.6,  '1912\nM7.3'),
        (29.4, 40.7,  '1999\nM7.4'),
        (28.9, 40.45, '2025\nM6.2'),
    ]
    for hx, hy, label in hist_eq:
        ax.plot(hx, hy, '*', color='#FFD700', markersize=12,
                markeredgecolor='#333', markeredgewidth=0.5, zorder=8)
        ax.text(hx + 0.08, hy + 0.02, label, fontsize=5,
                color='#333', zorder=8)

    # ----- Repeating earthquake clusters -----
    rep_x = [27.2, 27.5, 27.8, 28.0, 28.2]
    rep_y = [40.54, 40.48, 40.44, 40.43, 40.42]
    ax.scatter(rep_x, rep_y, s=15, c='#E91E63', marker='o',
               edgecolors='#333', linewidths=0.3, zorder=7,
               label='Repeating earthquakes')

    # ----- Plate motion arrows -----
    ax.annotate('', xy=(27.0, 40.2), xytext=(27.7, 40.2),
                arrowprops=dict(arrowstyle='->', color='#555', lw=1.5))
    ax.text(27.35, 40.15, 'Anatolian plate\n~24 mm/yr W',
            fontsize=5.5, ha='center', color='#555')

    ax.annotate('', xy=(28.5, 41.1), xytext=(28.5, 41.0),
                arrowprops=dict(arrowstyle='->', color='#555', lw=0.8))
    ax.text(28.5, 41.12, 'Eurasian plate (fixed)',
            fontsize=5.5, ha='center', color='#555')

    # ----- Fault label -----
    ax.text(27.8, 40.36, 'Main Marmara Fault (NAF)', fontsize=6,
            fontstyle='italic', rotation=-3, color='#333', zorder=5)

    # ----- Tectonic zone labels -----
    zone_box = dict(boxstyle='round,pad=0.2', facecolor='#F5F0E8',
                    alpha=0.8, edgecolor='#ccc')
    ax.text(28.3, 40.9, 'Istanbul Zone', fontsize=6.5, color='#666',
            fontstyle='italic', ha='center', bbox=zone_box)
    ax.text(28.5, 40.25, 'Sakarya Zone', fontsize=6.5, color='#666',
            fontstyle='italic', ha='center', bbox=zone_box)

    # ----- Axes and legend -----
    ax.set_xlim(26.3, 29.6)
    ax.set_ylim(40.05, 41.2)
    ax.set_xlabel('Longitude (°E)', fontsize=9)
    ax.set_ylabel('Latitude (°N)', fontsize=9)
    ax.legend(loc='lower right', fontsize=6, framealpha=0.9,
              edgecolor='#ccc')
    ax.set_aspect(1.3)
    ax.tick_params(labelsize=7)

    # Scale bar
    ax.plot([29.0, 29.4], [40.1, 40.1], 'k-', linewidth=2)
    ax.text(29.2, 40.07, '~30 km', fontsize=5.5, ha='center')

    fig.tight_layout()
    fig.savefig(save_path, dpi=dpi, bbox_inches='tight', facecolor='white')
    plt.close(fig)
    print(f'Saved → {save_path}')


if __name__ == '__main__':
    create_figure()
