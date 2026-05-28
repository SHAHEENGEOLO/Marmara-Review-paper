"""
Figure 1: Tectonic Map of the Marmara Sea Region (from real shapefiles)
========================================================================
Uses actual fault segment geometries and sub-basin polygons from GIS data.
Faults colour-coded by Fault_Segm attribute mapped to coupling behaviour.
"""
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.lines import Line2D
from matplotlib.colors import Normalize
from matplotlib import cm
import numpy as np
import geopandas as gpd

# ── Shapefile paths (adjust to your local directory structure) ──
BASIN_SHP = 'shapefiles/Marmara_SubBasins_edges_aligned_with_attributes/'
FAULT_SHP = 'shapefiles/Marmara_Faults_GeoParams_VERIFIED/'

# ── Load & reproject to WGS84 ──
basins = gpd.read_file(BASIN_SHP).to_crs(epsg=4326)
faults = gpd.read_file(FAULT_SHP).to_crs(epsg=4326)

# ── Segment → colour & label mapping ──
seg_style = {
    'Cinarcik_Basin_segment':  {'color': '#D32F2F', 'label': 'Çınarcık Basin seg. (locked, C=0.97)', 'zorder': 7},
    'Kumburgaz_Basin_segment': {'color': '#C62828', 'label': 'Kumburgaz Basin seg. (locked, C=1.00)', 'zorder': 7},
    'Central_High_segment':    {'color': '#FF9800', 'label': 'Central High seg. (partial, C=0.53)', 'zorder': 6},
    'Westren_High_segment':    {'color': '#FFC107', 'label': 'Western High seg. (partial, C=0.30)', 'zorder': 6},
    'Central_Basin_segment':   {'color': '#1E88E5', 'label': 'Central Basin seg. (creeping, C=0.22)', 'zorder': 5},
    'Tekirdag_Basin_segment':  {'color': '#42A5F5', 'label': 'Tekirdağ Basin seg. (partial, C=0.43)', 'zorder': 5},
}

# ── Basin fill colours ──
basin_fill = {
    'Basin': {'fc': '#B3D9F2', 'ec': '#5A9AC6', 'alpha': 0.55},
    'High':  {'fc': '#E8D5B7', 'ec': '#A08060', 'alpha': 0.50},
}

# ── Figure ──
fig, ax = plt.subplots(figsize=(10, 6.5), dpi=300)
ax.set_facecolor('#F7F3ED')

# Plot basins
for _, row in basins.iterrows():
    style = basin_fill.get(row['Type'], basin_fill['Basin'])
    gpd.GeoSeries([row.geometry]).plot(
        ax=ax, facecolor=style['fc'], edgecolor=style['ec'],
        linewidth=1.2, alpha=style['alpha'], zorder=2)

# Basin labels at centroids
for _, row in basins.iterrows():
    c = row.geometry.centroid
    code = row['Code']
    name = row['Name']
    fs = 6.5 if row['Type'] == 'Basin' else 5.5
    clr = '#1A5276' if row['Type'] == 'Basin' else '#6D4C1A'
    ax.text(c.x, c.y, f"{name}\n({row['Area_km2']:.0f} km²)",
            fontsize=fs, ha='center', va='center', color=clr,
            fontweight='bold', fontstyle='italic', zorder=10,
            bbox=dict(boxstyle='round,pad=0.15', fc='white', alpha=0.65, ec='none'))

# Plot faults coloured by segment
for seg_name, style in seg_style.items():
    subset = faults[faults['Fault_Segm'] == seg_name]
    if len(subset) == 0:
        continue
    subset.plot(ax=ax, color=style['color'], linewidth=0.9,
                alpha=0.85, zorder=style['zorder'])

# ── Cities ──
cities = [
    ('Istanbul', 28.97, 41.015),
    ('Bursa', 29.06, 40.20),
    ('Tekirdağ', 27.51, 40.98),
    ('Yalova', 29.27, 40.655),
    ('Silivri', 28.25, 41.075),
]
for name, cx, cy in cities:
    ax.plot(cx, cy, 's', color='#222', markersize=4.5, zorder=11)
    ax.text(cx + 0.04, cy + 0.025, name, fontsize=6.5, fontweight='bold',
            color='#222', zorder=11)

# ── Historical earthquakes ──
hist_eq = [
    (27.2, 40.62, '1912 M7.3\n(Ganos)', '#FFD600'),
    (29.37, 40.72, '1999 M7.4\n(İzmit)', '#FFD600'),
    (28.90, 40.50, '2025 M6.2', '#FF6D00'),
]
for hx, hy, label, clr in hist_eq:
    ax.plot(hx, hy, '*', color=clr, markersize=14, markeredgecolor='#333',
            markeredgewidth=0.6, zorder=12)
    ax.text(hx + 0.06, hy - 0.05, label, fontsize=5.5, color='#333',
            fontweight='bold', zorder=12)

# ── Plate motion arrows ──
ax.annotate('', xy=(27.05, 40.28), xytext=(27.7, 40.28),
            arrowprops=dict(arrowstyle='->', color='#444', lw=1.8))
ax.text(27.37, 40.22, 'Anatolian Plate\n~24 mm/yr WSW', fontsize=6,
        ha='center', color='#444', fontstyle='italic')

ax.annotate('', xy=(28.9, 41.15), xytext=(28.9, 41.08),
            arrowprops=dict(arrowstyle='->', color='#444', lw=1.0))
ax.text(28.9, 41.17, 'Eurasian Plate (fixed)', fontsize=5.5,
        ha='center', color='#444', fontstyle='italic')

# ── Tectonic zone labels ──
zone_box = dict(boxstyle='round,pad=0.25', facecolor='#F7F3ED', alpha=0.8, edgecolor='#BBB')
ax.text(28.3, 41.08, 'ISTANBUL ZONE', fontsize=7, color='#555',
        fontstyle='italic', ha='center', fontweight='bold', bbox=zone_box)
ax.text(28.8, 40.28, 'SAKARYA ZONE', fontsize=7, color='#555',
        fontstyle='italic', ha='center', fontweight='bold', bbox=zone_box)

# ── Legend ──
legend_handles = []
for seg_name, style in seg_style.items():
    legend_handles.append(Line2D([0], [0], color=style['color'], lw=2.5, label=style['label']))
legend_handles.append(mpatches.Patch(facecolor='#B3D9F2', edgecolor='#5A9AC6', label='Sub-basin (marine)'))
legend_handles.append(mpatches.Patch(facecolor='#E8D5B7', edgecolor='#A08060', label='Structural high'))
legend_handles.append(Line2D([0], [0], marker='*', color='w', markerfacecolor='#FFD600',
                             markeredgecolor='#333', markersize=10, label='Historical earthquake'))

ax.legend(handles=legend_handles, loc='lower left', fontsize=5.5,
          framealpha=0.92, edgecolor='#ccc', ncol=1,
          title='Fault Segments (by coupling)', title_fontsize=6)

# ── Scale bar ──
ax.plot([29.0, 29.35], [40.18, 40.18], 'k-', linewidth=2.5, zorder=10)
ax.text(29.175, 40.155, '~25 km', fontsize=6, ha='center', fontweight='bold')

# ── Axes ──
ax.set_xlim(26.6, 29.65)
ax.set_ylim(40.10, 41.22)
ax.set_xlabel('Longitude (°E)', fontsize=10)
ax.set_ylabel('Latitude (°N)', fontsize=10)
ax.tick_params(labelsize=8)
ax.set_aspect(1.35)
ax.grid(True, alpha=0.15, linewidth=0.5)

fig.tight_layout()
fig.savefig('fig1_real_shapefile.png', dpi=300, bbox_inches='tight', facecolor='white')
plt.close(fig)
print('Done — fig1_real_shapefile.png')
