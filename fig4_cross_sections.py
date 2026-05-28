"""
Figure 4: Synthetic Velocity Cross-Sections Through the Marmara Sea Crust
==========================================================================
(a) Along-strike W→E section showing basin-to-basin Vp variation with Moho.
(b) Across-strike N→S section through the Central Basin showing fault-zone
    damage and velocity contrast between the Istanbul and Sakarya zones.

References:
    - Laigle et al. (2008), GJI, doi:10.1111/j.1365-246X.2008.03835.x
    - Polat et al. (2016), EPS, doi:10.1186/s40623-016-0503-4
    - Bécel et al. (2009), Tectonophysics, doi:10.1016/j.tecto.2008.10.022
    - Karabulut (2025), GJI, doi:10.1093/gji/ggaf024

Authors: S.M.S. Ahmed & H. Güneyli
"""

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np


def create_figure(save_path='fig4_cross_sections.png', dpi=300):
    fig, (ax_a, ax_b) = plt.subplots(2, 1, figsize=(7.5, 6), dpi=dpi)

    # ================================================================
    # Panel (a): Along-strike cross-section (W → E, ~150 km)
    # ================================================================
    dist_h = np.linspace(0, 150, 200)
    depth  = np.linspace(0, 30, 100)
    DH, DP = np.meshgrid(dist_h, depth)

    # Background Vp increasing with depth
    vp = 5.0 + 0.1 * DP + 0.003 * DP**2 - 5e-5 * DP**3
    vp = np.clip(vp, 4.5, 8.0)

    # Basin low-Vp anomalies
    vp -= 1.2 * np.exp(-((DH - 20)**2  / (2 * 12**2) +
                          (DP - 3)**2   / (2 * 4**2)))    # Tekirdağ
    vp -= 1.8 * np.exp(-((DH - 65)**2  / (2 * 10**2) +
                          (DP - 4)**2   / (2 * 5**2)))    # Central
    vp -= 1.4 * np.exp(-((DH - 120)**2 / (2 * 15**2) +
                          (DP - 3)**2   / (2 * 4**2)))    # Çınarcık

    # Structural-high positive anomalies
    vp += 0.5 * np.exp(-((DH - 42)**2 / (2 * 8**2) +
                          (DP - 8)**2  / (2 * 5**2)))     # Western High
    vp += 0.4 * np.exp(-((DH - 95)**2 / (2 * 8**2) +
                          (DP - 8)**2  / (2 * 5**2)))     # Central High

    # Moho (undulating beneath basins)
    moho = 28 - 3 * np.exp(-((dist_h - 80)**2) / (2 * 40**2))

    im_a = ax_a.pcolormesh(DH, DP, vp, cmap='RdYlBu_r',
                           vmin=4.5, vmax=7.5, shading='gouraud')
    ax_a.plot(dist_h, moho, 'k--', linewidth=1.5, label='Moho')

    # Schematic seismicity
    rng = np.random.default_rng(123)
    eq_d = rng.uniform(5, 150, 80)
    eq_z = rng.uniform(3, 18, 80)
    ax_a.scatter(eq_d, eq_z, s=5, c='white', edgecolors='k',
                 linewidths=0.3, alpha=0.7, label='Seismicity')

    # Labels
    ax_a.text(20,  2, 'TB', fontsize=7, fontweight='bold',
              color='white', ha='center')
    ax_a.text(65,  2, 'CB', fontsize=7, fontweight='bold',
              color='white', ha='center')
    ax_a.text(120, 2, 'ÇB', fontsize=7, fontweight='bold',
              color='white', ha='center')
    ax_a.text(42, 6, 'WH', fontsize=6, color='#333', ha='center')
    ax_a.text(95, 6, 'CH', fontsize=6, color='#333', ha='center')

    ax_a.set_ylim(30, 0)
    ax_a.set_xlim(0, 150)
    ax_a.set_ylabel('Depth (km)', fontsize=8)
    ax_a.set_title('(a) Along-strike Vp cross-section (W → E)',
                   fontsize=9, fontweight='bold')
    ax_a.tick_params(labelsize=7)
    ax_a.legend(fontsize=6, loc='lower right')

    cb_a = plt.colorbar(im_a, ax=ax_a, pad=0.02, shrink=0.8)
    cb_a.set_label('Vp (km/s)', fontsize=7)
    cb_a.ax.tick_params(labelsize=6)

    # ================================================================
    # Panel (b): Across-strike cross-section (N → S, ~60 km)
    # ================================================================
    dist_ns = np.linspace(0, 60, 150)
    DN, DPn = np.meshgrid(dist_ns, depth)

    vp_ns = 5.5 + 0.08 * DPn + 0.003 * DPn**2

    # Central Basin anomaly
    vp_ns -= 1.5 * np.exp(-((DN - 30)**2 / (2 * 8**2) +
                              (DPn - 4)**2 / (2 * 5**2)))

    # Cross-fault velocity contrast
    vp_ns += 0.4 * np.tanh((DN - 30) * 0.3) * (1 - np.exp(-DPn / 5))

    # Narrow fault-zone damage (low-V)
    vp_ns -= 0.8 * np.exp(-((DN - 30)**2  / (2 * 2**2)) *
                            np.exp(-((DPn - 8)**2 / (2 * 10**2))))

    # Moho (N-S)
    moho_ns = 28 - 2 * np.exp(-((dist_ns - 30)**2) / (2 * 15**2))

    im_b = ax_b.pcolormesh(DN, DPn, vp_ns, cmap='RdYlBu_r',
                           vmin=4.5, vmax=7.5, shading='gouraud')
    ax_b.plot(dist_ns, moho_ns, 'k--', linewidth=1.5, label='Moho')
    ax_b.plot([30, 30], [0, 18], 'r-', linewidth=2, label='MMF')
    ax_b.annotate('NAF', xy=(30, 1), fontsize=7, color='red',
                  fontweight='bold', ha='center')

    # Zone labels
    lbl_box = dict(boxstyle='round', facecolor='#333', alpha=0.5)
    ax_b.text(10, 5, 'Istanbul Zone\n(N)', fontsize=6.5, color='white',
              ha='center', fontweight='bold', bbox=lbl_box)
    ax_b.text(50, 5, 'Sakarya Zone\n(S)', fontsize=6.5, color='white',
              ha='center', fontweight='bold', bbox=lbl_box)

    ax_b.set_ylim(30, 0)
    ax_b.set_xlim(0, 60)
    ax_b.set_xlabel('Distance N→S (km)', fontsize=8)
    ax_b.set_ylabel('Depth (km)', fontsize=8)
    ax_b.set_title(
        '(b) Across-strike Vp cross-section through Central Basin (N → S)',
        fontsize=9, fontweight='bold')
    ax_b.tick_params(labelsize=7)
    ax_b.legend(fontsize=6, loc='lower right')

    cb_b = plt.colorbar(im_b, ax=ax_b, pad=0.02, shrink=0.8)
    cb_b.set_label('Vp (km/s)', fontsize=7)
    cb_b.ax.tick_params(labelsize=6)

    fig.tight_layout()
    fig.savefig(save_path, dpi=dpi, bbox_inches='tight', facecolor='white')
    plt.close(fig)
    print(f'Saved → {save_path}')


if __name__ == '__main__':
    create_figure()
