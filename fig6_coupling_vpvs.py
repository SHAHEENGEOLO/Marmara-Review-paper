"""
Figure 6: Along-Strike Geophysical Parameter Profile
======================================================
Three-panel figure showing:
(a) P-wave velocity at ~5 km depth along the MMF (W→E).
(b) Vp/Vs ratio at comparable depth with elevated values flagged.
(c) Fault coupling status with repeating-earthquake density overlay.

References:
    - Polat et al. (2016), EPS, doi:10.1186/s40623-016-0503-4
    - Becker et al. (2023), GRL, doi:10.1029/2022GL101471
    - Yamamoto et al. (2017), JGR, doi:10.1002/2016JB013608

Authors: S.M.S. Ahmed & H. Güneyli
"""

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np


def create_figure(save_path='fig6_coupling_vpvs.png', dpi=300):
    fig, (ax_a, ax_b, ax_c) = plt.subplots(
        3, 1, figsize=(7.5, 6.5), dpi=dpi,
        gridspec_kw={'height_ratios': [1.2, 0.6, 1]})

    x = np.linspace(0, 150, 300)

    # ================================================================
    # Panel (a): Vp at 5 km depth
    # ================================================================
    vp_bg = 6.2 + 0.1 * np.sin(x / 30)
    vp_bg -= 0.4 * np.exp(-((x -  25)**2) / (2 * 12**2))   # Tekirdağ
    vp_bg -= 0.6 * np.exp(-((x -  65)**2) / (2 * 10**2))   # Central
    vp_bg -= 0.5 * np.exp(-((x - 120)**2) / (2 * 15**2))   # Çınarcık
    vp_bg += 0.3 * np.exp(-((x -  45)**2) / (2 *  8**2))   # W. High
    vp_bg += 0.25* np.exp(-((x -  95)**2) / (2 *  8**2))   # C. High

    ax_a.fill_between(x, 5.2, vp_bg, alpha=0.3, color='#42A5F5')
    ax_a.plot(x, vp_bg, 'k-', linewidth=1.5)
    ax_a.axhline(y=6.2, color='#999', linestyle='--', linewidth=0.8,
                 label='Regional average')

    for label, xp, yp, c in [
            ('Tekirdağ\nBasin', 25, 5.85, '#1565C0'),
            ('Central\nBasin',  65, 5.70, '#1565C0'),
            ('Çınarcık\nBasin',120, 5.75, '#1565C0'),
            ('W. High', 45, 6.55, '#C62828'),
            ('C. High', 95, 6.50, '#C62828')]:
        ax_a.annotate(label, xy=(xp, yp), fontsize=5 if 'High' in label else 6,
                      ha='center', color=c, fontweight='bold')

    ax_a.set_ylabel('Vp (km/s)\nat 5 km depth', fontsize=7)
    ax_a.set_ylim(5.2, 6.8)
    ax_a.set_xlim(0, 150)
    ax_a.set_xticklabels([])
    ax_a.tick_params(labelsize=6)
    ax_a.set_title('(a) Along-strike P-wave velocity', fontsize=8,
                   fontweight='bold')
    ax_a.legend(fontsize=5.5)

    # ================================================================
    # Panel (b): Vp/Vs ratio
    # ================================================================
    rng = np.random.default_rng(99)
    vpvs = (1.73
            + 0.05 * np.exp(-((x -  25)**2) / (2 * 12**2))
            + 0.08 * np.exp(-((x -  65)**2) / (2 * 10**2))
            + 0.06 * np.exp(-((x - 120)**2) / (2 * 15**2))
            - 0.03 * np.exp(-((x -  45)**2) / (2 *  8**2))
            - 0.02 * np.exp(-((x -  95)**2) / (2 *  8**2))
            + rng.standard_normal(300) * 0.005)

    ax_b.fill_between(x, 1.70, vpvs, where=vpvs > 1.76,
                      alpha=0.3, color='#FF5722', label='Elevated Vp/Vs')
    ax_b.fill_between(x, 1.70, vpvs, where=vpvs <= 1.76,
                      alpha=0.2, color='#4CAF50')
    ax_b.plot(x, vpvs, 'k-', linewidth=1)
    ax_b.axhline(y=1.73, color='#999', linestyle='--', linewidth=0.8)

    ax_b.set_ylabel('Vp/Vs', fontsize=7)
    ax_b.set_ylim(1.70, 1.85)
    ax_b.set_xlim(0, 150)
    ax_b.set_xticklabels([])
    ax_b.tick_params(labelsize=6)
    ax_b.set_title('(b) Along-strike Vp/Vs ratio', fontsize=8,
                   fontweight='bold')
    ax_b.text(65, 1.82, 'Fluid-rich?', fontsize=5.5, color='#BF360C',
              ha='center', fontstyle='italic')
    ax_b.legend(fontsize=5, loc='upper right')

    # ================================================================
    # Panel (c): Coupling status + repeater density
    # ================================================================
    # Colour-coded segments
    ax_c.axhspan(0, 1, xmin=0/150,  xmax=50/150,
                 color='#2196F3', alpha=0.7)
    ax_c.axhspan(0, 1, xmin=50/150, xmax=85/150,
                 color='#FF9800', alpha=0.7)
    ax_c.axhspan(0, 1, xmin=85/150, xmax=150/150,
                 color='#F44336', alpha=0.7)

    ax_c.text(25,   0.5, 'CREEPING\n(~10–15 mm/yr)',
              fontsize=7, ha='center', va='center',
              fontweight='bold', color='white')
    ax_c.text(67.5, 0.5, 'TRANSITION',
              fontsize=7, ha='center', va='center',
              fontweight='bold', color='white')
    ax_c.text(117,  0.5, 'LOCKED\n(strain accumulating)',
              fontsize=7, ha='center', va='center',
              fontweight='bold', color='white')

    # Repeater density curve (twin axis)
    rep_x = np.linspace(0, 150, 50)
    rep_y = (8 * np.exp(-(rep_x - 20)**2 / (2 * 20**2))
             + 5 * np.exp(-(rep_x - 55)**2 / (2 * 15**2))
             + 0.5)

    ax_tw = ax_c.twinx()
    ax_tw.plot(rep_x, rep_y, 'ko-', markersize=2, linewidth=1.2,
               label='Repeater density')
    ax_tw.set_ylabel('Repeater count', fontsize=6, color='#333')
    ax_tw.set_ylim(0, 12)
    ax_tw.tick_params(labelsize=5)
    ax_tw.legend(fontsize=5, loc='upper right')

    ax_c.annotate('→ Istanbul', xy=(148, 0.85), fontsize=7,
                  fontweight='bold', color='#333', ha='right')
    ax_c.set_xlim(0, 150)
    ax_c.set_ylim(0, 1)
    ax_c.set_xlabel('Distance along fault from W to E (km)', fontsize=8)
    ax_c.set_yticks([])
    ax_c.tick_params(labelsize=6)
    ax_c.set_title('(c) Fault coupling status and repeater density',
                   fontsize=8, fontweight='bold')

    fig.tight_layout()
    fig.savefig(save_path, dpi=dpi, bbox_inches='tight', facecolor='white')
    plt.close(fig)
    print(f'Saved → {save_path}')


if __name__ == '__main__':
    create_figure()
