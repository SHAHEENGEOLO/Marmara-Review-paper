"""
Figure 8: Medical CT vs Seismic Tomography — Imaging the Unseen
================================================================
Improved by Dr. Shaheen Mohammed Saleh Ahmed
Split-panel conceptual figure comparing medical CT scanning of the
human body to seismic tomography of the Marmara crust, highlighting
the analogy: both reveal hidden structure through wave transmission.
"""
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyArrowPatch, Arc, Circle, Ellipse, Wedge
import numpy as np
from matplotlib.colors import LinearSegmentedColormap
from scipy.ndimage import gaussian_filter

DPI = 300

fig = plt.figure(figsize=(10, 7), dpi=DPI)

# ══════════════════════════════════════════════
# LEFT PANEL — Medical CT Scan Analogy
# ══════════════════════════════════════════════
ax_left = fig.add_axes([0.03, 0.08, 0.44, 0.82])
ax_left.set_xlim(-1.5, 1.5)
ax_left.set_ylim(-1.6, 1.6)
ax_left.set_aspect('equal')
ax_left.axis('off')

# CT scanner ring
theta = np.linspace(0, 2*np.pi, 200)
r_outer = 1.35
r_inner = 1.15
ax_left.fill(r_outer*np.cos(theta), r_outer*np.sin(theta), fc='#37474F', ec='#263238', lw=2, zorder=1)
ax_left.fill(r_inner*np.cos(theta), r_inner*np.sin(theta), fc='#ECEFF1', ec='none', zorder=2)

# Body cross-section (ellipse with internal structures)
body = Ellipse((0, 0), 1.6, 1.8, fc='#FFCCBC', ec='#BF360C', lw=1.5, zorder=3, alpha=0.85)
ax_left.add_patch(body)

# Spine
spine = Ellipse((0, -0.25), 0.18, 0.22, fc='#F5F5F5', ec='#9E9E9E', lw=1, zorder=4)
ax_left.add_patch(spine)

# Lungs
lung_l = Ellipse((-0.35, 0.15), 0.35, 0.55, fc='#E1BEE7', ec='#7B1FA2', lw=0.8, zorder=4, alpha=0.7)
lung_r = Ellipse((0.35, 0.15), 0.35, 0.55, fc='#E1BEE7', ec='#7B1FA2', lw=0.8, zorder=4, alpha=0.7)
ax_left.add_patch(lung_l)
ax_left.add_patch(lung_r)

# Heart
heart = Ellipse((0.05, 0.1), 0.22, 0.28, fc='#EF5350', ec='#C62828', lw=1, zorder=5, alpha=0.8)
ax_left.add_patch(heart)

# Tumour (hidden danger)
tumour = Circle((-0.2, 0.3), 0.08, fc='#FFD600', ec='#F57F17', lw=1.5, zorder=6)
ax_left.add_patch(tumour)
ax_left.annotate('Hidden\nlesion', xy=(-0.2, 0.3), xytext=(-0.9, 0.9),
                fontsize=6, fontweight='bold', color='#F57F17',
                arrowprops=dict(arrowstyle='->', color='#F57F17', lw=1.2),
                ha='center', va='center', zorder=10)

# X-ray beams
for angle_deg in [30, 90, 150, 210, 270, 330]:
    a = np.radians(angle_deg)
    x0 = 1.25 * np.cos(a)
    y0 = 1.25 * np.sin(a)
    ax_left.annotate('', xy=(-x0*0.05, -y0*0.05), xytext=(x0, y0),
                    arrowprops=dict(arrowstyle='->', color='#43A047', lw=0.8, alpha=0.5))

# Source/detector labels
ax_left.plot(1.25, 0, '>', color='#43A047', markersize=6, zorder=8)
ax_left.text(1.25, -0.12, 'X-ray\nsource', fontsize=4.5, ha='center', color='#2E7D32')
ax_left.plot(-1.25, 0, 's', color='#1565C0', markersize=5, zorder=8)
ax_left.text(-1.25, -0.12, 'Detector', fontsize=4.5, ha='center', color='#1565C0')

ax_left.set_title('(a) Medical CT Scan', fontsize=11, fontweight='bold', pad=10, color='#263238')
ax_left.text(0, -1.45, 'X-rays traverse body → absorption map\n→ reveals hidden internal structures',
            fontsize=6.5, ha='center', color='#555', fontstyle='italic')


# ══════════════════════════════════════════════
# RIGHT PANEL — Seismic Tomography of the Earth
# ══════════════════════════════════════════════
ax_right = fig.add_axes([0.52, 0.08, 0.46, 0.82])
ax_right.set_xlim(0, 150)
ax_right.set_ylim(35, -2)  # depth increasing downward
ax_right.set_facecolor('#F5F0E8')

# Velocity model (cross-section)
x = np.linspace(0, 150, 200)
z = np.linspace(0, 35, 120)
X, Z = np.meshgrid(x, z)

# Background velocity
Vp = 5.0 + 0.08*Z + 0.002*Z**2
Vp = np.clip(Vp, 4.5, 7.8)

# Basin anomalies (the "hidden pathologies")
Vp -= 1.5 * np.exp(-((X-30)**2/(2*12**2) + (Z-4)**2/(2*4**2)))   # Tekirdağ
Vp -= 2.0 * np.exp(-((X-70)**2/(2*10**2) + (Z-5)**2/(2*5**2)))   # Central — deep
Vp -= 1.3 * np.exp(-((X-115)**2/(2*14**2) + (Z-3)**2/(2*4**2)))  # Çınarcık

# Locked asperity (high-V)
Vp += 0.8 * np.exp(-((X-115)**2/(2*8**2) + (Z-12)**2/(2*4**2)))

# Fluid pocket (very low V — the "hidden danger")
Vp -= 1.8 * np.exp(-((X-90)**2/(2*5**2) + (Z-8)**2/(2*3**2)))

# Fault damage zone
Vp -= 0.6 * np.exp(-((X-70)**2/(2*2**2))) * np.exp(-Z/15)

cmap = LinearSegmentedColormap.from_list('tomo',
    ['#1565C0', '#42A5F5', '#90CAF9', '#FFF9C4', '#FFCC80', '#EF6C00', '#BF360C'])

im = ax_right.pcolormesh(X, Z, Vp, cmap=cmap, vmin=4.0, vmax=7.5, shading='gouraud', zorder=1)

# Sea surface
ax_right.fill_between([0, 150], -2, 0, color='#B3D9F2', alpha=0.4, zorder=2)
ax_right.text(75, -1, 'Sea of Marmara', fontsize=7, ha='center', color='#1A5276',
             fontstyle='italic', zorder=3)

# Station triangles on surface
stn_x = [5, 20, 35, 50, 65, 80, 95, 110, 125, 140]
for sx in stn_x:
    ax_right.plot(sx, 0, 'v', color='#2E7D32', markersize=5, zorder=5)

# OBS on seafloor
obs_x = [25, 45, 65, 85, 105]
for ox in obs_x:
    ax_right.plot(ox, 0.5, 'D', color='#0D47A1', markersize=4, zorder=5)

# Seismic rays (like X-rays in CT)
np.random.seed(77)
for _ in range(12):
    eq_x = np.random.uniform(20, 130)
    eq_z = np.random.uniform(5, 18)
    stn = np.random.choice(stn_x)
    ax_right.plot([eq_x, stn], [eq_z, 0], '-', color='#43A047', lw=0.4, alpha=0.35, zorder=3)
    ax_right.plot(eq_x, eq_z, 'o', color='#FFD600', markersize=2, markeredgecolor='#F57F17',
                 markeredgewidth=0.3, zorder=4)

# Moho
moho = 28 - 3*np.exp(-((x-75)**2)/(2*35**2))
ax_right.plot(x, moho, 'k--', lw=1.5, zorder=4, label='Moho')
ax_right.text(15, 29, 'Moho', fontsize=6, color='#333')

# Fault line
fx = np.array([10, 30, 50, 70, 90, 110, 130, 145])
fz = np.array([1, 3, 6, 9, 11, 10, 7, 4])
ax_right.plot(fx, fz, 'r-', lw=2, zorder=4, label='MMF')

# Label the "hidden danger" — fluid pocket
ax_right.annotate('Fluid-rich\nzone?', xy=(90, 8), xytext=(130, 18),
                fontsize=7, fontweight='bold', color='#F57F17',
                arrowprops=dict(arrowstyle='->', color='#F57F17', lw=1.5),
                ha='center', zorder=8,
                bbox=dict(boxstyle='round,pad=0.2', fc='#FFF9C4', ec='#F57F17', alpha=0.9))

# Locked asperity label
ax_right.annotate('Locked\nasperity', xy=(115, 12), xytext=(140, 25),
                fontsize=6, fontweight='bold', color='#BF360C',
                arrowprops=dict(arrowstyle='->', color='#BF360C', lw=1.2),
                ha='center', zorder=8,
                bbox=dict(boxstyle='round,pad=0.2', fc='#FFCCBC', ec='#BF360C', alpha=0.9))

# Basin labels
ax_right.text(30, 2, 'TB', fontsize=7, fontweight='bold', color='white', ha='center', zorder=6)
ax_right.text(70, 2, 'CB', fontsize=7, fontweight='bold', color='white', ha='center', zorder=6)
ax_right.text(115, 2, 'ÇB', fontsize=7, fontweight='bold', color='white', ha='center', zorder=6)

# Legend
from matplotlib.lines import Line2D
legend_items = [
    Line2D([0], [0], marker='v', color='w', markerfacecolor='#2E7D32', markersize=6, label='Land station'),
    Line2D([0], [0], marker='D', color='w', markerfacecolor='#0D47A1', markersize=5, label='OBS'),
    Line2D([0], [0], marker='o', color='w', markerfacecolor='#FFD600', markeredgecolor='#F57F17', markersize=5, label='Earthquake'),
    Line2D([0], [0], color='#43A047', lw=1, alpha=0.5, label='Seismic ray'),
    Line2D([0], [0], color='red', lw=2, label='MMF'),
]
ax_right.legend(handles=legend_items, fontsize=5.5, loc='lower left', framealpha=0.9, edgecolor='#ccc')

cb = fig.colorbar(im, ax=ax_right, pad=0.02, shrink=0.75, aspect=25)
cb.set_label('Vp (km/s)', fontsize=8)
cb.ax.tick_params(labelsize=6)

ax_right.set_xlabel('Distance along fault (km)', fontsize=9)
ax_right.set_ylabel('Depth (km)', fontsize=9)
ax_right.tick_params(labelsize=7)
ax_right.set_title('(b) Seismic Tomography of the Marmara Crust', fontsize=11, fontweight='bold', pad=10, color='#263238')
ax_right.text(75, 34, 'Seismic waves traverse crust → velocity map\n→ reveals hidden faults, fluids, and locked zones',
            fontsize=6.5, ha='center', color='#555', fontstyle='italic')

# ══════════════════════════════════════════════
# CENTRAL ARROW connecting the two panels
# ══════════════════════════════════════════════
ax_arrow = fig.add_axes([0.44, 0.35, 0.1, 0.3])
ax_arrow.set_xlim(0, 1)
ax_arrow.set_ylim(0, 1)
ax_arrow.axis('off')
ax_arrow.annotate('', xy=(0.9, 0.5), xytext=(0.1, 0.5),
                 arrowprops=dict(arrowstyle='->', color='#333', lw=2.5,
                                connectionstyle='arc3,rad=0'))
ax_arrow.text(0.5, 0.72, 'Same\nPrinciple', fontsize=7.5, ha='center', va='center',
             fontweight='bold', color='#333')
ax_arrow.text(0.5, 0.28, 'Waves →\nImage', fontsize=6.5, ha='center', va='center',
             fontstyle='italic', color='#666')

fig.savefig('fig8_ct_vs_seismic.png', dpi=DPI, bbox_inches='tight', facecolor='white')
plt.close(fig)
print('Done — fig8_ct_vs_seismic.png')
