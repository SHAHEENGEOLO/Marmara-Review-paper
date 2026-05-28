"""
Figure 7 / Figure 2 style: Real 3-D GIS map of the Marmara Sea region
=======================================================================

This updated enlarged-core version replaces the synthetic bathymetry and hand-drawn
fault trace with the real input files:

1) F7_2024.tif bathymetry/elevation raster
2) Marmara_Faults_GeoParams.shp fault segments
3) Marmara_SubBasins_edges_aligned_with_attributes.shp sub-basin polygons

Main improvements:
- Larger publication format with tighter focus on the inner Marmara basin–fault core.
- Real TIFF values are used for the 3-D surface and colorbar.
- Faults are plotted from the real shapefile and colored by HAZ_SCOR.
- Sub-basin polygons/edges are plotted from the real shapefile.
- Stronger 3-D visual effects: hillshade, fault shadows, semi-transparent
  basin panels, optional down-dip fault curtains, contour floor, north arrow,
  and clear legends.

Author: Dr. Shaheen Mohammed Saleh Ahmed
"""

from __future__ import annotations

import os
from pathlib import Path
from typing import Iterable, Sequence

import geopandas as gpd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
import numpy as np
import rasterio
from rasterio.enums import Resampling
from rasterio.windows import Window, from_bounds, intersection
from scipy.ndimage import gaussian_filter
from shapely.geometry import LineString, MultiLineString, Polygon, MultiPolygon, box
from matplotlib import cm
from matplotlib.colors import LightSource, Normalize
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from mpl_toolkits.mplot3d.art3d import Poly3DCollection


# ---------------------------------------------------------------------
# 1. INPUT PATHS
# ---------------------------------------------------------------------
# Your Windows paths. Keep the r before the string to avoid backslash errors.
WINDOWS_RASTER_PATH = r"C:\Users\shahe\Desktop\manuscipts for publication\Marmara Review paper\code of paper\F7_2024.tif\F7_2024.tif"
WINDOWS_FAULTS_PATH = r"C:\Users\shahe\Desktop\manuscipts for publication\Marmara Review paper\code of paper\GIS FILES\Marmara Faults\Marmara_Faults_GeoParams_VERIFIED\Marmara_Faults_GeoParams.shp"
WINDOWS_BASINS_PATH = r"C:\Users\shahe\Desktop\multi modal\GIS FILES\Marmara_SubBasins_edges_aligned_with_attributes\Marmara_SubBasins_edges_aligned_with_attributes.shp"

# Sandbox fallback paths, useful only inside ChatGPT after upload.
SANDBOX_RASTER_PATH = "/mnt/data/F7_2024(1).tif"
SANDBOX_FAULTS_PATH = "/mnt/data/Marmara_Faults_GeoParams(1).shp"
SANDBOX_BASINS_PATH = "/mnt/data/Marmara_SubBasins_edges_aligned_with_attributes(2).shp"


# ---------------------------------------------------------------------
# 2. SETTINGS YOU CAN SAFELY CHANGE
# ---------------------------------------------------------------------
OUTPUT_PNG = "Marmara_F7_real_3D_structural_map_big_inner_core.png"
OUTPUT_PDF = None  # Set to "Marmara_F7_real_3D_structural_map.pdf" if you also need PDF

FIGSIZE = (18, 11.5)        # Larger canvas so the Marmara core fills the figure
DPI = 300                   # Publication-quality export
MAX_GRID_SIZE = 620         # Sharper surface for the enlarged inner Marmara view
VERTICAL_EXAGGERATION = 1.55  # Stronger 3-D relief for the inner Marmara basins
SMOOTH_SIGMA = 0.7          # Gentle smoothing for 3-D display only
CROP_PADDING_DEG = 0.015    # Very small margin so the Marmara basin/fault core appears larger
CENTER_ZOOM_FACTOR = 1.00  # Crop already focuses on basin/fault core; keep 1.00 to avoid cutting real features
CENTER_ZOOM_Y_SHIFT = 0.00 # Positive shifts crop north; negative shifts south
ROBUST_COLOR = True         # True = stronger color contrast using 2-98 percentiles
DRAW_FAULT_CURTAINS = True  # 3-D down-dip fault panels for the longest/highest-risk faults
MAX_CURTAIN_FAULTS = 25     # Number of fault curtains to draw, to avoid clutter
FAULT_DEPTH_SCALE = 0.10    # Visual scale for LOCK_DEP km in 3-D curtains
VIEW_ELEVATION = 38
VIEW_AZIMUTH = -64
MAP_BOX_ZOOM = 1.32       # Visual enlargement of the 3-D map inside the figure canvas

# Label/legend readability controls.
DECLUTTER_LABELS = True       # Move labels away from crowded geometries.
LABEL_Z_LIFT = 0.42           # Vertical lift for offset labels.
LABEL_LEADER_LIFT = 0.22      # Lift at feature anchor before leader line.
MAX_BASIN_LABELS = 12         # Keep all basin labels if <= this value.
SHOW_KEY_FAULT_LABELS = False # Default False: keep fault names out of the map to avoid overlap.
MAX_FAULT_LABELS = 5          # Used only if SHOW_KEY_FAULT_LABELS=True.
SHOW_INTERNAL_REGION_LABEL = False  # Keep the map clean; title is outside the 3-D axes.
FOCUS_ON_BASIN_FAULT_CORE = True  # Use real basin+fault bounds, not broad outer TIFF context.
INNER_CORE_PADDING_DEG = 0.015     # Extra padding around real basins and faults; increase to 0.03 for more context.


# ---------------------------------------------------------------------
# 3. PATH AND GEOMETRY HELPERS
# ---------------------------------------------------------------------
def resolve_existing_path(*candidates: str | os.PathLike) -> Path:
    """Return the first existing path from a list of possible paths."""
    for candidate in candidates:
        p = Path(candidate)
        if p.exists():
            return p
    checked = "\n".join(str(Path(c)) for c in candidates)
    raise FileNotFoundError(f"None of these paths exist:\n{checked}")


def iter_lines(geom) -> Iterable[LineString]:
    """Yield LineString objects from lines or polygon boundaries."""
    if geom is None or geom.is_empty:
        return
    if isinstance(geom, LineString):
        yield geom
    elif isinstance(geom, MultiLineString):
        yield from geom.geoms
    elif isinstance(geom, Polygon):
        yield LineString(geom.exterior.coords)
    elif isinstance(geom, MultiPolygon):
        for poly in geom.geoms:
            yield LineString(poly.exterior.coords)
    elif hasattr(geom, "geoms"):
        # Handles GeometryCollection objects that may be produced by map clipping.
        for part in geom.geoms:
            yield from iter_lines(part)


def iter_polygons(geom) -> Iterable[Polygon]:
    """Yield Polygon objects from Polygon or MultiPolygon geometries."""
    if geom is None or geom.is_empty:
        return
    if isinstance(geom, Polygon):
        yield geom
    elif isinstance(geom, MultiPolygon):
        yield from geom.geoms
    elif hasattr(geom, "geoms"):
        # Handles GeometryCollection objects that may be produced by map clipping.
        for part in geom.geoms:
            yield from iter_polygons(part)


def thin_xy(x: Sequence[float], y: Sequence[float], max_points: int = 350):
    """Reduce very dense geometry for faster 3-D plotting while preserving shape."""
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    if len(x) <= max_points:
        return x, y
    idx = np.linspace(0, len(x) - 1, max_points).astype(int)
    return x[idx], y[idx]


def add_3d_text(
    ax,
    x,
    y,
    z,
    text,
    fontsize=11,
    color="white",
    ha="center",
    va="center",
    bbox=True,
    box_edge="white",
):
    """Bold 3-D label with halo and optional dark rounded box."""
    bbox_kwargs = None
    if bbox:
        bbox_kwargs = dict(
            boxstyle="round,pad=0.28,rounding_size=0.12",
            facecolor=(0, 0, 0, 0.62),
            edgecolor=box_edge,
            linewidth=1.25,
        )
    t = ax.text(
        x, y, z, text,
        fontsize=fontsize,
        color=color,
        fontweight="bold",
        ha=ha,
        va=va,
        zorder=140,
        bbox=bbox_kwargs,
    )
    t.set_path_effects([pe.withStroke(linewidth=3.2, foreground="black")])
    return t


def add_offset_3d_label(
    ax,
    anchor_x,
    anchor_y,
    anchor_z,
    label_x,
    label_y,
    label_z,
    text,
    color="white",
    fontsize=10.5,
    ha="center",
):
    """Place a label away from its feature and connect it using a leader line."""
    # leader line shadow
    ax.plot(
        [anchor_x, label_x], [anchor_y, label_y], [anchor_z, label_z - 0.04],
        color="black", linewidth=2.7, alpha=0.60, zorder=118,
    )
    # bright leader line
    ax.plot(
        [anchor_x, label_x], [anchor_y, label_y], [anchor_z, label_z - 0.04],
        color=color, linewidth=1.25, alpha=0.96, zorder=121,
    )
    ax.scatter(
        [anchor_x], [anchor_y], [anchor_z],
        s=34, color=color, edgecolors="black", linewidths=0.8, alpha=0.95, zorder=122,
    )
    return add_3d_text(
        ax, label_x, label_y, label_z, text,
        fontsize=fontsize, color="white", ha=ha, va="center", bbox=True, box_edge=color,
    )


def short_label_name(name: str) -> str:
    """Make long geological names shorter so labels do not overlap."""
    replacements = {
        "Basin": "B.",
        "Central": "Cent.",
        "Western": "W.",
        "Eastern": "E.",
        "Çınarcık": "Çın.",
        "Tekirdağ": "Tek.",
        "Kumburgaz": "Kumb.",
        "High": "High",
    }
    out = str(name)
    for old, new in replacements.items():
        out = out.replace(old, new)
    return out


# ---------------------------------------------------------------------
# 4. MAIN FIGURE FUNCTION
# ---------------------------------------------------------------------
def create_real_3d_marmara_figure(
    raster_path: str | os.PathLike | None = None,
    faults_path: str | os.PathLike | None = None,
    basins_path: str | os.PathLike | None = None,
    output_png: str | os.PathLike = OUTPUT_PNG,
    output_pdf: str | os.PathLike = OUTPUT_PDF,
):
    """Create a large, readable 3-D Marmara map from real GIS files."""

    # Resolve paths. If you run this on your PC, it will use your Windows paths.
    raster_path = Path(raster_path) if raster_path else resolve_existing_path(WINDOWS_RASTER_PATH, SANDBOX_RASTER_PATH)
    faults_path = Path(faults_path) if faults_path else resolve_existing_path(WINDOWS_FAULTS_PATH, SANDBOX_FAULTS_PATH)
    basins_path = Path(basins_path) if basins_path else resolve_existing_path(WINDOWS_BASINS_PATH, SANDBOX_BASINS_PATH)

    print("Using raster:", raster_path)
    print("Using faults:", faults_path)
    print("Using basins:", basins_path)

    faults = gpd.read_file(faults_path)
    basins = gpd.read_file(basins_path)

    with rasterio.open(raster_path) as src:
        raster_crs = src.crs
        if raster_crs is None:
            raise ValueError("Raster has no CRS. Please define the CRS before plotting.")

        # Reproject vectors to the raster CRS. Your TIFF is EPSG:4326 and the shapefiles are UTM.
        if faults.crs is None:
            raise ValueError("Fault shapefile has no CRS. Define it before plotting.")
        if basins.crs is None:
            raise ValueError("Sub-basin shapefile has no CRS. Define it before plotting.")

        faults = faults.to_crs(raster_crs)
        basins = basins.to_crs(raster_crs)

        # Light simplification in geographic degrees for faster 3-D rendering;
        # it preserves topology and does not change the attribute values.
        faults["geometry"] = faults.geometry.simplify(0.00012, preserve_topology=True)
        basins["geometry"] = basins.geometry.simplify(0.00020, preserve_topology=True)

        # Crop the raster tightly around the real Marmara basin-fault core.
        # This makes the TIFF surface, sub-basins, and faults visually larger
        # instead of leaving broad empty map context around the central features.
        if FOCUS_ON_BASIN_FAULT_CORE:
            all_bounds = np.vstack([faults.total_bounds, basins.total_bounds])
            pad = INNER_CORE_PADDING_DEG
        else:
            all_bounds = np.vstack([faults.total_bounds, basins.total_bounds])
            pad = CROP_PADDING_DEG

        minx = float(np.min(all_bounds[:, 0]) - pad)
        miny = float(np.min(all_bounds[:, 1]) - pad)
        maxx = float(np.max(all_bounds[:, 2]) + pad)
        maxy = float(np.max(all_bounds[:, 3]) + pad)

        # Clamp to raster bounds.
        rb = src.bounds
        minx = max(minx, rb.left)
        maxx = min(maxx, rb.right)
        miny = max(miny, rb.bottom)
        maxy = min(maxy, rb.top)

        # -------------------------------------------------------------
        # Inner Marmara zoom
        # -------------------------------------------------------------
        # The first crop is based on all vector data. This second crop
        # gently zooms into the centre of the TIFF/vector overlap so the
        # Marmara Sea surface fills more of the figure. Use values between
        # 0.78 and 0.92 for a little zoom without losing too much context.
        if CENTER_ZOOM_FACTOR and CENTER_ZOOM_FACTOR < 1.0:
            cx = (minx + maxx) / 2.0
            cy = (miny + maxy) / 2.0 + CENTER_ZOOM_Y_SHIFT * (maxy - miny)
            new_dx = (maxx - minx) * float(CENTER_ZOOM_FACTOR)
            new_dy = (maxy - miny) * float(CENTER_ZOOM_FACTOR)
            minx = max(cx - new_dx / 2.0, rb.left)
            maxx = min(cx + new_dx / 2.0, rb.right)
            miny = max(cy - new_dy / 2.0, rb.bottom)
            maxy = min(cy + new_dy / 2.0, rb.top)

        # Clip vector layers to the zoomed map extent. This keeps off-frame
        # fault and basin labels from being calculated and prevents overlap.
        plot_extent = box(minx, miny, maxx, maxy)
        faults = faults[faults.geometry.intersects(plot_extent)].copy()
        basins = basins[basins.geometry.intersects(plot_extent)].copy()
        faults["geometry"] = faults.geometry.intersection(plot_extent)
        basins["geometry"] = basins.geometry.intersection(plot_extent)
        faults = faults[~faults.geometry.is_empty & faults.geometry.notnull()].copy()
        basins = basins[~basins.geometry.is_empty & basins.geometry.notnull()].copy()
        if faults.empty or basins.empty:
            raise ValueError(
                "The zoomed map extent removed all faults or basins. Increase CENTER_ZOOM_FACTOR, for example 0.90."
            )

        win = from_bounds(minx, miny, maxx, maxy, src.transform)
        win = win.round_offsets().round_lengths()
        full_win = Window(0, 0, src.width, src.height)
        win = intersection(win, full_win)

        scale = max(float(win.width), float(win.height)) / float(MAX_GRID_SIZE)
        scale = max(1.0, scale)
        out_width = max(2, int(round(float(win.width) / scale)))
        out_height = max(2, int(round(float(win.height) / scale)))

        elev_m = src.read(
            1,
            window=win,
            out_shape=(out_height, out_width),
            resampling=Resampling.bilinear,
            masked=True,
        ).astype("float64")

        # Correct transform after display resampling.
        transform = src.window_transform(win) * rasterio.Affine.scale(
            float(win.width) / out_width,
            float(win.height) / out_height,
        )

    elev_m = np.ma.filled(elev_m, np.nan)
    elev_m[~np.isfinite(elev_m)] = np.nan

    # Very gentle smoothing only for surface rendering. Original raster values still control the map.
    valid_median = float(np.nanmedian(elev_m))
    elev_for_smooth = np.where(np.isfinite(elev_m), elev_m, valid_median)
    if SMOOTH_SIGMA and SMOOTH_SIGMA > 0:
        elev_for_plot_m = gaussian_filter(elev_for_smooth, sigma=SMOOTH_SIGMA)
    else:
        elev_for_plot_m = elev_for_smooth

    z_surface = (elev_for_plot_m / 1000.0) * VERTICAL_EXAGGERATION

    # Grid coordinates at cell centers.
    cols = np.arange(out_width) + 0.5
    rows = np.arange(out_height) + 0.5
    lon = transform.c + transform.a * cols
    lat = transform.f + transform.e * rows
    LON, LAT = np.meshgrid(lon, lat)

    def sample_z(xs, ys, lift=0.045):
        """Nearest-neighbour z sampling from the plotted surface."""
        xs = np.asarray(xs, dtype=float)
        ys = np.asarray(ys, dtype=float)
        cc, rr = (~transform) * (xs, ys)
        ci = np.clip(np.rint(cc).astype(int), 0, out_width - 1)
        ri = np.clip(np.rint(rr).astype(int), 0, out_height - 1)
        zz = z_surface[ri, ci]
        zz = np.where(np.isfinite(zz), zz, (valid_median / 1000.0) * VERTICAL_EXAGGERATION)
        return zz + lift

    # Color normalization.
    if ROBUST_COLOR:
        vmin, vmax = np.nanpercentile(elev_m, [2, 98])
    else:
        vmin, vmax = np.nanmin(elev_m), np.nanmax(elev_m)
    norm = Normalize(vmin=float(vmin), vmax=float(vmax))
    surface_cmap = plt.get_cmap("turbo")

    # Hillshade for stronger 3-D relief.
    ls = LightSource(azdeg=315, altdeg=45)
    face_rgb = ls.shade(
        elev_for_plot_m,
        cmap=surface_cmap,
        norm=norm,
        vert_exag=0.85,
        blend_mode="soft",
    )

    # -----------------------------------------------------------------
    # Plot
    # -----------------------------------------------------------------
    fig = plt.figure(figsize=FIGSIZE, dpi=DPI)
    ax = fig.add_subplot(111, projection="3d")
    fig.patch.set_facecolor("white")

    # Main 3-D raster surface.
    ax.plot_surface(
        LON,
        LAT,
        z_surface,
        facecolors=face_rgb,
        rstride=1,
        cstride=1,
        linewidth=0,
        antialiased=True,
        shade=False,
        alpha=0.98,
    )

    # Contour floor for extra 3-D depth impression.
    z_floor = float(np.nanmin(z_surface)) - 0.55
    ax.contour(
        LON,
        LAT,
        z_surface,
        levels=13,
        zdir="z",
        offset=z_floor,
        cmap="Greys",
        linewidths=0.7,
        alpha=0.55,
    )

    # Sub-basin transparent polygons and bright edges.
    basin_colors = ["#00E5FF", "#FFD166", "#8AFF00", "#FF6B6B", "#B967FF", "#00FFA3"]
    basin_legend_items = []
    basin_label_records = []
    for i, row in basins.iterrows():
        name = str(row.get("Name", row.get("Code", f"Basin {i+1}")))
        code = str(row.get("Code", i + 1))
        color = basin_colors[i % len(basin_colors)]
        geom = row.geometry

        # Semi-transparent colored panels.
        for poly in iter_polygons(geom):
            x, y = poly.exterior.xy
            x, y = thin_xy(x, y, max_points=250)
            z = sample_z(x, y, lift=0.075)
            verts = [list(zip(x, y, z))]
            panel = Poly3DCollection(
                verts,
                facecolor=color,
                edgecolor=color,
                linewidth=1.3,
                alpha=0.20,
            )
            ax.add_collection3d(panel)
            ax.plot(x, y, z + 0.012, color="black", linewidth=3.2, alpha=0.45)
            ax.plot(x, y, z + 0.020, color=color, linewidth=2.2, alpha=0.98)

        # Store basin label anchors. Labels are placed later in separated slots
        # so they do not overlap the basin edges, faults, or each other.
        try:
            c = geom.representative_point()
            anchor_z = sample_z([c.x], [c.y], lift=LABEL_LEADER_LIFT)[0]
            label = f"{code} | {short_label_name(name)}"
            basin_label_records.append(
                dict(x=c.x, y=c.y, z=anchor_z, label=label, color=color, name=name, code=code)
            )
        except Exception:
            pass

        basin_legend_items.append(Patch(facecolor=color, edgecolor=color, alpha=0.40, label=f"{code}: {name}"))

    # Decluttered sub-basin labels: placed in two clean rows near map margins.
    # This avoids the common problem of overlapping labels inside narrow Marmara basins.
    dx = maxx - minx
    dy = maxy - miny
    if DECLUTTER_LABELS and basin_label_records:
        label_records = sorted(basin_label_records[:MAX_BASIN_LABELS], key=lambda r: r["x"])
        n = len(label_records)
        top_records = label_records[::2]
        bottom_records = label_records[1::2]

        def place_label_row(records, y_position, top=True):
            if not records:
                return
            for k, rec in enumerate(records):
                x_position = minx + (k + 0.70) / (len(records) + 0.40) * dx
                z_position = sample_z([x_position], [y_position], lift=LABEL_Z_LIFT + (0.04 if top else 0.00))[0]
                add_offset_3d_label(
                    ax,
                    rec["x"], rec["y"], rec["z"],
                    x_position, y_position, z_position,
                    rec["label"],
                    color=rec["color"],
                    fontsize=9.4,
                    ha="center",
                )

        place_label_row(top_records, maxy - 0.045 * dy, top=True)
        place_label_row(bottom_records, miny + 0.055 * dy, top=False)
    else:
        for rec in basin_label_records:
            add_3d_text(ax, rec["x"], rec["y"], rec["z"], rec["label"], fontsize=9.4, color="white", box_edge=rec["color"])

    # Fault segments colored by hazard score when available.
    hazard_col = "HAZ_SCOR" if "HAZ_SCOR" in faults.columns else None
    if hazard_col:
        haz_values = faults[hazard_col].astype(float).to_numpy()
        haz_norm = Normalize(vmin=float(np.nanmin(haz_values)), vmax=float(np.nanmax(haz_values)))
        fault_cmap = plt.get_cmap("plasma")
    else:
        haz_norm = Normalize(0, 1)
        fault_cmap = plt.get_cmap("plasma")

    # Draw fault curtains first so bright fault traces stay on top.
    if DRAW_FAULT_CURTAINS:
        sort_col = hazard_col if hazard_col else ("length_km" if "length_km" in faults.columns else None)
        if sort_col:
            curtain_faults = faults.sort_values(sort_col, ascending=False).head(MAX_CURTAIN_FAULTS)
        else:
            curtain_faults = faults.head(MAX_CURTAIN_FAULTS)

        for _, row in curtain_faults.iterrows():
            lock_depth = float(row.get("LOCK_DEP", 5.0)) if "LOCK_DEP" in faults.columns else 5.0
            geom = row.geometry
            color_value = float(row.get(hazard_col, 1.0)) if hazard_col else 1.0
            color = fault_cmap(haz_norm(color_value))
            for line in iter_lines(geom):
                x, y = line.xy
                x, y = thin_xy(x, y, max_points=70)
                if len(x) < 2:
                    continue
                z_top = sample_z(x, y, lift=0.115)
                z_bottom = z_top - lock_depth * FAULT_DEPTH_SCALE
                curtain_verts = list(zip(x, y, z_top)) + list(zip(x[::-1], y[::-1], z_bottom[::-1]))
                curtain = Poly3DCollection(
                    [curtain_verts],
                    facecolor=color,
                    edgecolor=color,
                    linewidth=0.45,
                    alpha=0.18,
                )
                ax.add_collection3d(curtain)

    fault_label_records = []
    for _, row in faults.iterrows():
        geom = row.geometry
        color_value = float(row.get(hazard_col, 1.0)) if hazard_col else 1.0
        color = fault_cmap(haz_norm(color_value))
        length_km = float(row.get("length_km", 1.0)) if "length_km" in faults.columns else 1.0
        lw = 1.4 + min(2.4, 0.12 * np.sqrt(max(length_km, 0.1)))
        for line in iter_lines(geom):
            x, y = line.xy
            x, y = thin_xy(x, y, max_points=80)
            if len(x) < 2:
                continue
            z = sample_z(x, y, lift=0.16)
            # Shadow + bright top line = stronger 3-D readability.
            ax.plot(x, y, z - 0.06, color="black", linewidth=lw + 2.5, alpha=0.42, zorder=50)
            ax.plot(x, y, z, color=color, linewidth=lw, alpha=1.0, zorder=80)

        # Record only the most important fault labels; do not label every segment.
        if SHOW_KEY_FAULT_LABELS:
            try:
                p = geom.representative_point()
                score = float(row.get(hazard_col, 0.0)) if hazard_col else length_km
                raw_name = str(row.get("fault_name", row.get("Fault_Segm", "Fault")))
                short_name = raw_name.replace("aciklari", "offshore").replace("marmara denizi", "Marmara").replace("marmara sea", "Marmara")
                if len(short_name) > 28:
                    short_name = short_name[:25] + "…"
                fault_label_records.append(
                    dict(x=p.x, y=p.y, z=sample_z([p.x], [p.y], lift=0.27)[0], label=short_name, score=score, color=color)
                )
            except Exception:
                pass

    # Highlight seismic-gap / historical-rupture categories if present.
    if "Seismic_Ga" in faults.columns:
        try:
            gap = faults[faults["Seismic_Ga"].astype(str).str.contains("gap", case=False, na=False)]
            for _, row in gap.head(120).iterrows():
                for line in iter_lines(row.geometry):
                    x, y = line.xy
                    x, y = thin_xy(x, y, max_points=50)
                    z = sample_z(x, y, lift=0.23)
                    ax.plot(x, y, z, color="#00FFFF", linewidth=0.75, alpha=0.82, zorder=95)
        except Exception:
            pass

    # A small number of important fault labels, separated along the side margin.
    # Too many fault labels would hide the real fault network.
    if SHOW_KEY_FAULT_LABELS and fault_label_records:
        key_faults = sorted(fault_label_records, key=lambda r: r["score"], reverse=True)[:MAX_FAULT_LABELS]
        for k, rec in enumerate(key_faults):
            x_position = maxx - 0.020 * dx
            y_position = miny + (k + 0.80) / (len(key_faults) + 0.75) * dy
            z_position = sample_z([x_position], [y_position], lift=LABEL_Z_LIFT + 0.12)[0]
            add_offset_3d_label(
                ax, rec["x"], rec["y"], rec["z"],
                x_position, y_position, z_position,
                rec["label"],
                color="#FFEA00",
                fontsize=8.7,
                ha="right",
            )

    # Optional internal region label. Disabled by default to avoid covering map features.
    if SHOW_INTERNAL_REGION_LABEL:
        midx = (minx + maxx) / 2
        midy = maxy - 0.045 * (maxy - miny)
        add_3d_text(
            ax, midx, midy, sample_z([midx], [midy], lift=0.42)[0],
            "MARMARA SEA STRUCTURAL SYSTEM", fontsize=16, color="#FFFFFF", bbox=True,
        )

    # Axes and view.
    # Compact axis labels prevent collision with the colorbars and annotations.
    ax.set_xlabel("Longitude (°E)", fontsize=11.5, fontweight="bold", labelpad=10)
    ax.set_ylabel("Latitude (°N)", fontsize=11.5, fontweight="bold", labelpad=10)
    ax.set_zlabel("Relief (km, VE)", fontsize=10.0, fontweight="bold", labelpad=6)
    ax.tick_params(axis="both", labelsize=9.0, pad=2)
    ax.tick_params(axis="z", labelsize=8.5, pad=2)

    ax.set_xlim(minx, maxx)
    ax.set_ylim(miny, maxy)
    ax.set_zlim(z_floor, float(np.nanmax(z_surface)) + 0.95)
    ax.view_init(elev=VIEW_ELEVATION, azim=VIEW_AZIMUTH)
    # The zoom argument enlarges the 3-D map inside the axes without changing
    # the real coordinates or attribute values. Older Matplotlib versions may
    # not support zoom, so a safe fallback is included.
    try:
        ax.set_box_aspect((maxx - minx, maxy - miny, 0.60), zoom=MAP_BOX_ZOOM)
    except TypeError:
        ax.set_box_aspect((maxx - minx, maxy - miny, 0.60))

    # Pane styling for a clean scientific 3-D frame.
    ax.xaxis.pane.set_facecolor((0.96, 0.96, 0.96, 0.35))
    ax.yaxis.pane.set_facecolor((0.96, 0.96, 0.96, 0.35))
    ax.zaxis.pane.set_facecolor((0.90, 0.90, 0.90, 0.22))
    ax.grid(True, alpha=0.35)

    # Title outside the 3-D axes so it never overlaps the map surface.
    fig.suptitle(
        "Real 3-D Bathymetric–Tectonic Map of the Marmara Sea Region\n"
        "Enlarged inner Marmara core | F7_2024 TIFF + Marmara faults + sub-basin boundaries",
        fontsize=19.5,
        fontweight="bold",
        y=0.978,
    )

    # Colorbars are placed in dedicated right-side axes, away from the 3-D map.
    mappable = cm.ScalarMappable(norm=norm, cmap=surface_cmap)
    mappable.set_array(elev_m)
    cax1 = fig.add_axes([0.855, 0.38, 0.016, 0.36])
    cbar = fig.colorbar(mappable, cax=cax1)
    cbar.set_label("Real TIFF elevation / bathymetry (m)", fontsize=10.8, fontweight="bold", labelpad=10)
    cbar.ax.tick_params(labelsize=9.2)

    # Fault hazard colorbar.
    if hazard_col:
        haz_map = cm.ScalarMappable(norm=haz_norm, cmap=fault_cmap)
        haz_map.set_array(haz_values)
        cax2 = fig.add_axes([0.935, 0.38, 0.016, 0.36])
        cbar2 = fig.colorbar(haz_map, cax=cax2)
        cbar2.set_label("Fault hazard score (HAZ_SCOR)", fontsize=10.8, fontweight="bold", labelpad=10)
        cbar2.ax.tick_params(labelsize=9.2)

    # Legends are placed outside the 3-D map panel to prevent overlap with labels/features.
    fault_handle = Line2D([0], [0], color=fault_cmap(0.95), lw=4, label="Real Marmara fault segments")
    gap_handle = Line2D([0], [0], color="#00FFFF", lw=2.2, label="Seismic-gap highlight")
    curtain_handle = Patch(facecolor=fault_cmap(0.70), edgecolor=fault_cmap(0.70), alpha=0.22, label="3-D fault-depth curtain")
    contour_handle = Line2D([0], [0], color="0.35", lw=1.2, label="Projected relief contours")
    leader_handle = Line2D([0], [0], color="white", lw=1.4, marker="o", markerfacecolor="white",
                           markeredgecolor="black", label="Offset labels + leader lines")

    fig.legend(
        handles=[fault_handle, gap_handle, curtain_handle, contour_handle, leader_handle],
        loc="lower center",
        bbox_to_anchor=(0.45, 0.065),
        ncol=5,
        fontsize=10.3,
        title="Map layers and annotation style",
        title_fontsize=11.5,
        frameon=True,
        framealpha=0.94,
        facecolor="white",
        edgecolor="black",
    )

    fig.legend(
        handles=basin_legend_items[:8],
        loc="lower center",
        bbox_to_anchor=(0.45, 0.013),
        ncol=min(4, max(1, len(basin_legend_items[:8]))),
        fontsize=8.7,
        title="Sub-basins / highs",
        title_fontsize=10.2,
        frameon=True,
        framealpha=0.94,
        facecolor="white",
        edgecolor="black",
    )

    # North arrow and caption note on the figure canvas.
    fig.text(0.885, 0.865, "N\n↑", ha="center", va="center", fontsize=22, fontweight="bold")
    fig.text(
        0.50,
        0.125,
        "Note: enlarged inner-core view; real coordinates and values are preserved, while labels use leader lines to avoid overlap.",
        ha="center",
        va="center",
        fontsize=11.5,
        color="black",
    )

    fig.subplots_adjust(left=0.005, right=0.835, top=0.90, bottom=0.155)
    fig.savefig(output_png, dpi=DPI, bbox_inches="tight", facecolor="white")
    if output_pdf:
        fig.savefig(output_pdf, dpi=DPI, bbox_inches="tight", facecolor="white")
    plt.close(fig)

    print(f"Saved PNG: {output_png}")
    if output_pdf:
        print(f"Saved PDF: {output_pdf}")


if __name__ == "__main__":
    create_real_3d_marmara_figure()
