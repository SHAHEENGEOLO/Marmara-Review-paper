"""
Figure 5: LET Inversion Workflow and Checkerboard Resolution Test
==================================================================
(a) Flowchart of the standard local earthquake tomography workflow as
    applied in Marmara studies.
(b) Conceptual checkerboard recovery test showing well-resolved anomalies
    in the upper crust and smeared recovery at depth/edges.

References:
    - Koulakov (2009), BSSA, doi:10.1785/0120080013
    - Polat et al. (2016), EPS, doi:10.1186/s40623-016-0503-4

Authors: S.M.S. Ahmed & H. Güneyli
"""

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
from matplotlib.colors import LinearSegmentedColormap
from scipy.ndimage import gaussian_filter


def create_figure(save_path='fig5_workflow_checker.png', dpi=300):
    fig, (ax_a, ax_b) = plt.subplots(1, 2, figsize=(7.5, 4.2), dpi=dpi)

    # ================================================================
    # Panel (a): LET Workflow Flowchart
    # ================================================================
    steps = [
        ('Earthquake\nCatalog',                   0.5, 0.92, '#E3F2FD', '#1565C0'),
        ('Phase Picking\n(P & S arrivals)',        0.5, 0.76, '#E3F2FD', '#1565C0'),
        ('1-D Velocity\nModel (initial)',          0.5, 0.60, '#FFF3E0', '#E65100'),
        ('3-D Iterative\nInversion\n(LOTOS/tomoDD)', 0.5, 0.42, '#E8F5E9', '#2E7D32'),
        ('Resolution\nAssessment\n(checkerboard, DWS)', 0.5, 0.22, '#FCE4EC', '#C62828'),
        ('Final Vp, Vs,\nVp/Vs Models',            0.5, 0.06, '#F3E5F5', '#6A1B9A'),
    ]

    for label, x, y, fc, ec in steps:
        box = mpatches.FancyBboxPatch(
            (x - 0.20, y - 0.06), 0.40, 0.12,
            boxstyle="round,pad=0.02",
            facecolor=fc, edgecolor=ec, linewidth=1.5)
        ax_a.add_patch(box)
        ax_a.text(x, y, label, ha='center', va='center',
                  fontsize=5.5, fontweight='bold', color=ec)

    # Downward arrows between steps
    for i in range(len(steps) - 1):
        ax_a.annotate(
            '', xy=(0.5, steps[i + 1][2] + 0.065),
            xytext=(0.5, steps[i][2] - 0.065),
            arrowprops=dict(arrowstyle='->', color='#555', lw=1.2))

    # Iteration feedback loop
    ax_a.annotate(
        '', xy=(0.78, 0.60), xytext=(0.78, 0.42),
        arrowprops=dict(arrowstyle='->', color='#999', lw=1,
                        connectionstyle='arc3,rad=0.3'))
    ax_a.text(0.88, 0.51, 'Iterate', fontsize=5, color='#999',
              rotation=90, ha='center', va='center')

    ax_a.set_xlim(0, 1)
    ax_a.set_ylim(-0.02, 1.0)
    ax_a.axis('off')
    ax_a.set_title('(a) LET Inversion Workflow',
                   fontsize=8, fontweight='bold', pad=5)

    # ================================================================
    # Panel (b): Checkerboard Recovery Test
    # ================================================================
    nx, ny = 8, 8
    checker = np.zeros((ny, nx))
    for i in range(ny):
        for j in range(nx):
            checker[i, j] = 1 if (i + j) % 2 == 0 else -1

    recovered = gaussian_filter(checker, sigma=0.6)

    # Degrade edges to simulate poor peripheral recovery
    mask = np.ones_like(recovered)
    mask[0, :]  = 0.3
    mask[-1, :] = 0.3
    mask[:, 0]  = 0.3
    mask[:, -1] = 0.3
    recovered_masked = recovered * mask

    cmap_ck = LinearSegmentedColormap.from_list(
        'vel', ['#2196F3', '#FFFFFF', '#F44336'])

    im = ax_b.imshow(recovered_masked, cmap=cmap_ck,
                     vmin=-1, vmax=1, aspect='equal',
                     interpolation='bilinear')

    ax_b.set_xlabel('Distance along fault (km)', fontsize=7)
    ax_b.set_ylabel('Depth (km)', fontsize=7)
    ax_b.set_xticks([0, 4, 7])
    ax_b.set_xticklabels(['0', '40', '70'], fontsize=6)
    ax_b.set_yticks([0, 3, 7])
    ax_b.set_yticklabels(['0', '15', '35'], fontsize=6)
    ax_b.set_title('(b) Checkerboard Recovery',
                   fontsize=8, fontweight='bold', pad=5)

    cbar = plt.colorbar(im, ax=ax_b, shrink=0.7, pad=0.05)
    cbar.set_label('dVp/Vp (%)', fontsize=6)
    cbar.ax.tick_params(labelsize=5)

    # Annotations
    ax_b.text(4, 3.5, 'Well\nresolved', fontsize=5.5, ha='center',
              va='center', fontweight='bold', color='#333',
              bbox=dict(boxstyle='round', facecolor='white',
                        alpha=0.7, edgecolor='#aaa'))
    ax_b.text(1, 7, 'Smeared', fontsize=5.5, ha='center',
              va='center', color='#333', fontstyle='italic',
              bbox=dict(boxstyle='round', facecolor='white',
                        alpha=0.7, edgecolor='#ccc'))

    fig.tight_layout()
    fig.savefig(save_path, dpi=dpi, bbox_inches='tight', facecolor='white')
    plt.close(fig)
    print(f'Saved → {save_path}')


if __name__ == '__main__':
    create_figure()
