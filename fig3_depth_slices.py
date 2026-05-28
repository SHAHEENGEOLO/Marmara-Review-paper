"""
Figure 3: Synthetic P-wave Velocity Perturbation Maps at Multiple Depths
=========================================================================
Four depth slices (3, 7, 12, 18 km) showing basin low-Vp anomalies,
structural-high positive anomalies, and fault-parallel velocity contrasts.

References:
    - Laigle et al. (2008), GJI, doi:10.1111/j.1365-246X.2008.03835.x
    - Polat et al. (2016), EPS, doi:10.1186/s40623-016-0503-4
    - Karabulut (2025), GJI, doi:10.1093/gji/ggaf024

Authors: S.M.S. Ahmed & H. Güneyli
"""

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LinearSegmentedColormap
from scipy.ndimage import gaussian_filter


def create_figure(save_path='fig3_depth_slices.png', dpi=300):
    fig, axes = plt.subplots(2, 2, figsize=(7.5, 7), dpi=dpi)

    cmap_vel = LinearSegmentedColormap.from_list(
        'vel', ['#2196F3', '#90CAF9', '#FFFFFF', '#EF9A9A', '#F44336'])

    # Sea of Marmara outline
    sea_x = [26.6, 27.0, 27.4, 28.0, 28.6, 29.2, 29.4,
             29.3, 29.0, 28.5, 28.0, 27.5, 27.0, 26.6]
    sea_y = [40.55, 40.45, 40.35, 40.35, 40.38, 40.5, 40.6,
             40.72, 40.75, 40.78, 40.75, 40.72, 40.65, 40.55]

    # MMF trace
    fx = np.array([26.7, 27.0, 27.4, 27.7, 28.0,
                   28.3, 28.6, 28.9, 29.15, 29.3])
    fy = np.array([40.58, 40.52, 40.47, 40.44, 40.42,
                   40.42, 40.44, 40.48, 40.55, 40.62])

    # Model grid
    lon_g = np.linspace(26.5, 29.5, 100)
    lat_g = np.linspace(40.2, 41.0, 80)
    LG, LT = np.meshgrid(lon_g, lat_g)

    depths = [3, 7, 12, 18]
    titles = ['z = 3 km', 'z = 7 km', 'z = 12 km', 'z = 18 km']

    for idx, (ax, depth, title) in enumerate(zip(axes.flat, depths, titles)):

        # Random background
        rng = np.random.default_rng(idx + 10)
        dv = rng.standard_normal((80, 100)) * 0.3
        dv = gaussian_filter(dv, sigma=4)

        # Scale factor: anomalies weaken with depth
        scale = max(0.2, 1.0 - depth / 25.0)

        # Basin low-Vp anomalies
        dv -= 5 * scale * np.exp(-((LG - 27.1)**2  / 0.15 +
                                    (LT - 40.48)**2 / 0.015))   # Tekirdağ
        dv -= 7 * scale * np.exp(-((LG - 27.85)**2 / 0.10 +
                                    (LT - 40.45)**2 / 0.012))   # Central
        dv -= 6 * scale * np.exp(-((LG - 29.0)**2  / 0.20 +
                                    (LT - 40.55)**2 / 0.015))   # Çınarcık

        # Structural-high positive anomalies
        dv += 4   * scale * np.exp(-((LG - 27.5)**2 / 0.06 +
                                      (LT - 40.52)**2 / 0.008))
        dv += 3.5 * scale * np.exp(-((LG - 28.4)**2 / 0.06 +
                                      (LT - 40.45)**2 / 0.008))

        # Deep high-V lower crust
        if depth > 10:
            dv += 3 * np.exp(-((LG - 28.0)**2 / 0.8 +
                                (LT - 40.6)**2 / 0.05))

        # Fault-parallel velocity contrast
        fault_y = 40.44 + 0.07 * np.sin((LG - 27.5) * 2)
        dv += 2 * scale * np.tanh((LT - fault_y) * 8)

        im = ax.pcolormesh(LG, LT, dv, cmap=cmap_vel,
                           vmin=-8, vmax=8, shading='gouraud')

        # Overlays
        ax.plot(fx, fy, 'k-', linewidth=1.5)
        ax.plot(sea_x, sea_y, 'k--', linewidth=0.5, alpha=0.4)

        if idx == 0:
            ax.text(27.1, 40.42, 'TB', fontsize=6, color='#1565C0',
                    fontweight='bold', ha='center')
            ax.text(27.85, 40.38, 'CB', fontsize=6, color='#1565C0',
                    fontweight='bold', ha='center')
            ax.text(29.0, 40.48, 'ÇB', fontsize=6, color='#1565C0',
                    fontweight='bold', ha='center')

        ax.set_xlim(26.5, 29.5)
        ax.set_ylim(40.2, 41.0)
        ax.set_aspect(1.3)
        ax.set_title(f'({chr(97 + idx)}) {title}',
                     fontsize=8, fontweight='bold')
        ax.tick_params(labelsize=6)
        if idx >= 2:
            ax.set_xlabel('Longitude (°E)', fontsize=7)
        if idx % 2 == 0:
            ax.set_ylabel('Latitude (°N)', fontsize=7)

    # Global colour bar
    cbar_ax = fig.add_axes([0.92, 0.15, 0.02, 0.7])
    cbar = fig.colorbar(im, cax=cbar_ax)
    cbar.set_label('δVp/Vp (%)', fontsize=8)
    cbar.ax.tick_params(labelsize=6)

    fig.suptitle(
        'Synthetic P-wave Velocity Perturbation Maps at Multiple Depths',
        fontsize=10, fontweight='bold', y=0.98)
    fig.subplots_adjust(right=0.9, hspace=0.25, wspace=0.2)
    fig.savefig(save_path, dpi=dpi, bbox_inches='tight', facecolor='white')
    plt.close(fig)
    print(f'Saved → {save_path}')


if __name__ == '__main__':
    create_figure()
