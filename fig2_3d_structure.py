"""
Figure 2: 3-D Structural Perspective of the Marmara Sea Region
================================================================
Three-dimensional perspective view showing simplified bathymetry with
sub-basins, the MMF surface trace and fault plane at depth, schematic
microseismicity, and the Moho discontinuity.

References:
    - Bécel et al. (2009), Tectonophysics, doi:10.1016/j.tecto.2008.10.022
    - Yamamoto et al. (2017), JGR, doi:10.1002/2016JB013608
    - Karabulut (2025), GJI, doi:10.1093/gji/ggaf024

Authors: S.M.S. Ahmed & H. Güneyli
"""

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from scipy.ndimage import gaussian_filter
from matplotlib import cm
from mpl_toolkits.mplot3d.art3d import Poly3DCollection


def create_figure(save_path='fig2_3d_structure.png', dpi=300):
    fig = plt.figure(figsize=(7.5, 5.5), dpi=dpi)
    ax = fig.add_subplot(111, projection='3d')

    # ---- Grid for surface ----
    lon = np.linspace(26.6, 29.4, 120)
    lat = np.linspace(40.3, 40.8, 80)
    LON, LAT = np.meshgrid(lon, lat)

    # Approximate sea mask (elliptical)
    sea_mask = ((LON - 28.0)**2 / 2.5 +
                (LAT - 40.55)**2 / 0.12) < 1

    # ---- Bathymetry (basins as negative anomalies) ----
    bathy = np.zeros_like(LON)
    bathy -= 0.8  * np.exp(-((LON - 27.1)**2  / 0.08 +
                              (LAT - 40.48)**2 / 0.005))   # Tekirdağ
    bathy -= 1.2  * np.exp(-((LON - 27.85)**2 / 0.06 +
                              (LAT - 40.45)**2 / 0.004))   # Central
    bathy -= 1.0  * np.exp(-((LON - 29.0)**2  / 0.10 +
                              (LAT - 40.55)**2 / 0.006))   # Çınarcık
    # Structural highs
    bathy += 0.30 * np.exp(-((LON - 27.5)**2 / 0.04 +
                              (LAT - 40.5)**2 / 0.005))
    bathy += 0.25 * np.exp(-((LON - 28.4)**2 / 0.04 +
                              (LAT - 40.45)**2 / 0.005))

    # Land elevation (flat)
    land_elev = 0.15 * np.ones_like(LON)
    land_elev[sea_mask] = 0

    # Composite surface
    surface = np.where(sea_mask, bathy, land_elev)
    surface = gaussian_filter(surface, sigma=2)

    # ---- Plot surfaces ----
    land_surface = np.ma.masked_where(sea_mask, surface)
    sea_surface  = np.ma.masked_where(~sea_mask, surface)

    ax.plot_surface(LON, LAT, sea_surface, cmap='ocean_r',
                    alpha=0.75, edgecolor='none', antialiased=True)
    ax.plot_surface(LON, LAT, land_surface, color='#C4B99A',
                    alpha=0.6, edgecolor='none', antialiased=True)

    # ---- MMF surface trace ----
    fault_lon = np.array([26.7, 27.0, 27.4, 27.7, 28.0,
                          28.3, 28.6, 28.9, 29.15, 29.3])
    fault_lat = np.array([40.58, 40.52, 40.47, 40.44, 40.42,
                          40.42, 40.44, 40.48, 40.55, 40.62])
    fault_z = np.zeros_like(fault_lon) - 0.05
    ax.plot(fault_lon, fault_lat, fault_z, 'r-', linewidth=2.5,
            zorder=10, label='MMF trace')

    # ---- Fault plane extending to depth ----
    fp_lon = np.array([27.5, 27.5, 28.5, 28.5])
    fp_lat = np.array([40.47, 40.47, 40.43, 40.43])
    fp_z   = np.array([0, -1.5, -1.5, 0])
    verts  = [list(zip(fp_lon, fp_lat, fp_z))]
    poly   = Poly3DCollection(verts, alpha=0.25, facecolor='red',
                               edgecolor='red', linewidth=0.5)
    ax.add_collection3d(poly)

    # ---- Basin labels ----
    ax.text(27.1,  40.48, -0.9, 'TB', fontsize=7, color='white',
            fontweight='bold', ha='center')
    ax.text(27.85, 40.45, -1.3, 'CB', fontsize=7, color='white',
            fontweight='bold', ha='center')
    ax.text(29.0,  40.55, -1.1, 'ÇB', fontsize=7, color='white',
            fontweight='bold', ha='center')

    # ---- Moho interface (flat plane) ----
    moho_lon = np.array([26.6, 26.6, 29.4, 29.4])
    moho_lat = np.array([40.3, 40.8, 40.8, 40.3])
    moho_z   = np.array([-2.6, -2.6, -2.6, -2.6])
    moho_verts = [list(zip(moho_lon, moho_lat, moho_z))]
    moho_poly  = Poly3DCollection(moho_verts, alpha=0.15,
                                   facecolor='#8B4513',
                                   edgecolor='#8B4513', linewidth=0.5)
    ax.add_collection3d(moho_poly)
    ax.text(27.0, 40.35, -2.5, 'Moho (~26–30 km)',
            fontsize=6, color='#8B4513')

    # ---- Schematic microseismicity ----
    rng = np.random.default_rng(42)
    eq_lon = rng.uniform(27.0, 29.2, 60)
    eq_lat = rng.uniform(40.4, 40.6, 60)
    eq_z   = rng.uniform(-1.5, -0.2, 60)
    ax.scatter(eq_lon, eq_lat, eq_z, s=3, c='yellow', alpha=0.6,
               edgecolors='orange', linewidths=0.3,
               label='Microseismicity')

    # ---- Axes ----
    ax.set_xlabel('Longitude (°E)', fontsize=7, labelpad=5)
    ax.set_ylabel('Latitude (°N)', fontsize=7, labelpad=5)
    ax.set_zlabel('Depth (relative)', fontsize=7, labelpad=3)
    ax.view_init(elev=35, azim=-55)
    ax.set_zlim(-3, 0.5)
    ax.tick_params(labelsize=5)
    ax.legend(fontsize=5.5, loc='upper left')
    ax.set_title('3-D Structural Overview of the Marmara Sea Region',
                 fontsize=9, fontweight='bold', pad=10)

    fig.tight_layout()
    fig.savefig(save_path, dpi=dpi, bbox_inches='tight', facecolor='white')
    plt.close(fig)
    print(f'Saved → {save_path}')


if __name__ == '__main__':
    create_figure()
