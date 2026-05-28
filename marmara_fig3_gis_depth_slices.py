"""
Figure 3: GIS-constrained synthetic P-wave velocity perturbation maps at multiple depths
=====================================================================================

This script upgrades the original schematic/synthetic Figure 3 code by using the
provided GIS layers directly:

1. F7_2024.tif bathymetry/elevation raster
2. Marmara_Faults_GeoParams.shp verified Marmara fault segments
3. Marmara_SubBasins_edges_aligned_with_attributes.shp Marmara sub-basin polygons

The produced figure is still a synthetic conceptual velocity perturbation model
because no real 3-D seismic tomography volume is provided. However, the anomaly
locations, basin shapes, fault overlays, and map extent are constrained by the
provided GIS files rather than manually drawn outlines.

Author: S.M.S. Ahmed & H. Güneyli
Updated GIS-based plotting version
"""

from __future__ import annotations

import os
import warnings
from pathlib import Path
from typing import Dict, Iterable, Tuple

import geopandas as gpd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import rasterio
from matplotlib.colors import LinearSegmentedColormap, Normalize
from matplotlib.lines import Line2D
from matplotlib.patches import Patch, Polygon as MplPolygon
from matplotlib.collections import LineCollection
from pyproj import Transformer
from rasterio.windows import from_bounds
from scipy.ndimage import gaussian_filter
from scipy.spatial import cKDTree
from shapely.geometry import LineString, MultiLineString

try:
    # Shapely >= 2.0, fast vectorized point-in-polygon operation
    from shapely import contains_xy
except Exception:  # pragma: no cover
    contains_xy = None

warnings.filterwarnings("ignore", category=UserWarning)

# =============================================================================
# 1. USER PATHS AND FALLBACK SANDBOX PATHS
# =============================================================================

# Your original Windows paths. Keep these unchanged for running on your computer.
WINDOWS_TIF = r"C:\Users\shahe\Desktop\manuscipts for publication\Marmara Review paper\code of paper\F7_2024.tif\F7_2024.tif"
WINDOWS_FAULTS = r"C:\Users\shahe\Desktop\manuscipts for publication\Marmara Review paper\code of paper\GIS FILES\Marmara Faults\Marmara_Faults_GeoParams_VERIFIED\Marmara_Faults_GeoParams.shp"
WINDOWS_BASINS = r"C:\Users\shahe\Desktop\multi modal\GIS FILES\Marmara_SubBasins_edges_aligned_with_attributes\Marmara_SubBasins_edges_aligned_with_attributes.shp"

# Sandbox fallback paths used only inside ChatGPT/Python testing.
SANDBOX_TIF = "/mnt/data/F7_2024(1).tif"
SANDBOX_FAULTS = "/mnt/data/Marmara_Faults_GeoParams(1).shp"
SANDBOX_BASINS = "/mnt/data/Marmara_SubBasins_edges_aligned_with_attributes(2).shp"

# Output files
OUTPUT_PNG = "Marmara_Fig3_GIS_depth_slices_decluttered.png"
OUTPUT_PDF = None  # Example: "Marmara_Fig3_GIS_depth_slices_decluttered.pdf"


def first_existing_path(*paths: str) -> str:
    """Return the first existing path. If none exists, return the first path."""
    for path in paths:
        if path and Path(path).exists():
            return path
    return paths[0]


TIF_PATH = first_existing_path(WINDOWS_TIF, SANDBOX_TIF)
FAULTS_PATH = first_existing_path(WINDOWS_FAULTS, SANDBOX_FAULTS)
BASINS_PATH = first_existing_path(WINDOWS_BASINS, SANDBOX_BASINS)

# =============================================================================
# 2. STYLE AND MODEL SETTINGS
# =============================================================================

FIGSIZE = (18.0, 12.2)
DPI = 350

# Number of nodes in each synthetic depth-slice grid.
NX = 430
NY = 230

# Focus on inner Marmara basin/fault core. Increase padding for more regional context.
INNER_CORE_PADDING_DEG = 0.050

# Map readability options.
# Values > 1 visually stretch latitude so the narrow Marmara basin core fills the panels better.
MAP_VERTICAL_STRETCH = 2.35
SHOW_BASIN_LABELS_IN_ALL_PANELS = True
SHOW_FAULT_LABELS = False      # False is recommended to avoid overlap.
SHOW_BATHY_CONTOURS = True
SHOW_RASTER_HILLSHADE = True
SHOW_GRID = True

# Synthetic depth slices in km.
DEPTHS_KM = [3, 7, 12, 18]
DEPTH_TITLES = ["z = 3 km", "z = 7 km", "z = 12 km", "z = 18 km"]

# Color map: negative perturbations blue, positive perturbations red.
CMAP_VEL = LinearSegmentedColormap.from_list(
    "marmara_vp_energy",
    ["#08306B", "#08519C", "#2171B5", "#6BAED6", "#DDEEFF",
     "#FFFFFF", "#FFF1D6", "#FCAE91", "#FB6A4A", "#CB181D", "#67000D"],
    N=256,
)
VP_VMIN, VP_VMAX = -8.0, 8.0

# Basin/high color scheme for GIS polygons.
BASIN_EDGE_COLORS = {
    "TB": "#FFD400",   # Tekirdag Basin
    "CB": "#FF8C00",   # Central Basin
    "KB": "#00E5FF",   # Kumburgaz Basin
    "CiB": "#00E676",  # Cinarcik Basin
    "CH": "#FF00FF",   # Central High
    "WH": "#A3FF12",   # Western High
}
DEFAULT_BASIN_COLOR = "#FFFF66"

# Codes treated as sedimentary/basin low-velocity domains and structural highs.
LOW_VP_CODES = {"TB", "CB", "KB", "CiB"}
HIGH_VP_CODES = {"CH", "WH"}

# Label positions are intentionally offset from the polygons to avoid overlapping.
# The first two values are x/y offsets in degrees from the representative point.
BASIN_LABEL_OFFSETS = {
    "TB": (-0.03,  0.075),
    "WH": (-0.12, -0.065),
    "CB": (-0.02,  0.075),
    "KB": ( 0.02, -0.075),
    "CH": ( 0.10,  0.060),
    "CiB": (0.10, -0.070),
}

PANEL_EXPLANATIONS = {
    3:  "Shallow basin low-Vp dominates",
    7:  "Basin anomalies attenuate",
    12: "Mid-crustal contrast emerges",
    18: "Deep positive domain appears",
}

# =============================================================================
# 3. DATA LOADING AND PREPARATION
# =============================================================================


def clean_geometries(gdf: gpd.GeoDataFrame) -> gpd.GeoDataFrame:
    """Remove empty geometries and repair invalid geometries where possible."""
    gdf = gdf.copy()
    gdf = gdf[~gdf.geometry.is_empty & gdf.geometry.notna()]
    try:
        gdf["geometry"] = gdf.geometry.make_valid()
    except Exception:
        gdf["geometry"] = gdf.geometry.buffer(0)
    return gdf


def load_gis_layers() -> Tuple[gpd.GeoDataFrame, gpd.GeoDataFrame, str]:
    """Load raster CRS plus reprojected fault and basin layers."""
    if not Path(TIF_PATH).exists():
        raise FileNotFoundError(f"TIFF file was not found: {TIF_PATH}")
    if not Path(FAULTS_PATH).exists():
        raise FileNotFoundError(f"Fault shapefile was not found: {FAULTS_PATH}")
    if not Path(BASINS_PATH).exists():
        raise FileNotFoundError(f"Basin shapefile was not found: {BASINS_PATH}")

    with rasterio.open(TIF_PATH) as src:
        raster_crs = src.crs or "EPSG:4326"

    faults = clean_geometries(gpd.read_file(FAULTS_PATH))
    basins = clean_geometries(gpd.read_file(BASINS_PATH))

    # Reproject vector GIS files to raster CRS for direct plotting on raster coordinates.
    if faults.crs is None:
        raise ValueError("Fault shapefile has no CRS. Please define the CRS before plotting.")
    if basins.crs is None:
        raise ValueError("Sub-basin shapefile has no CRS. Please define the CRS before plotting.")

    faults = faults.to_crs(raster_crs)
    basins = basins.to_crs(raster_crs)

    return faults, basins, str(raster_crs)


def compute_inner_extent(faults: gpd.GeoDataFrame, basins: gpd.GeoDataFrame) -> Tuple[float, float, float, float]:
    """Compute a tight extent around the real basin and fault GIS layers."""
    fb = faults.total_bounds
    bb = basins.total_bounds
    xmin = min(fb[0], bb[0]) - INNER_CORE_PADDING_DEG
    ymin = min(fb[1], bb[1]) - INNER_CORE_PADDING_DEG
    xmax = max(fb[2], bb[2]) + INNER_CORE_PADDING_DEG
    ymax = max(fb[3], bb[3]) + INNER_CORE_PADDING_DEG
    return xmin, xmax, ymin, ymax


def read_raster_background(extent: Tuple[float, float, float, float], out_shape: Tuple[int, int]) -> np.ndarray:
    """Read and resample raster inside the figure extent."""
    xmin, xmax, ymin, ymax = extent
    height, width = out_shape
    with rasterio.open(TIF_PATH) as src:
        window = from_bounds(xmin, ymin, xmax, ymax, src.transform)
        arr = src.read(
            1,
            window=window,
            out_shape=(height, width),
            boundless=True,
            fill_value=np.nan,
            masked=True,
        ).astype("float32")
        arr = np.asarray(arr.filled(np.nan), dtype="float32")
    # Raster rows are north-to-south; flip to match ascending latitude grid.
    return arr[::-1, :]


def robust_normalize(arr: np.ndarray, pmin: float = 2.0, pmax: float = 98.0) -> np.ndarray:
    """Normalize an array by robust percentiles."""
    finite = np.isfinite(arr)
    if finite.sum() == 0:
        return np.zeros_like(arr, dtype=float)
    lo, hi = np.nanpercentile(arr[finite], [pmin, pmax])
    if not np.isfinite(lo) or not np.isfinite(hi) or hi <= lo:
        return np.zeros_like(arr, dtype=float)
    return np.clip((arr - lo) / (hi - lo), 0, 1)


def make_hillshade(z: np.ndarray, azimuth: float = 315.0, altitude: float = 45.0) -> np.ndarray:
    """Simple hillshade from raster values."""
    z_fill = np.nan_to_num(z, nan=np.nanmedian(z[np.isfinite(z)]) if np.isfinite(z).any() else 0.0)
    dy, dx = np.gradient(z_fill)
    slope = np.pi / 2.0 - np.arctan(np.hypot(dx, dy))
    aspect = np.arctan2(-dx, dy)
    az = np.deg2rad(azimuth)
    alt = np.deg2rad(altitude)
    shaded = np.sin(alt) * np.sin(slope) + np.cos(alt) * np.cos(slope) * np.cos(az - aspect)
    return robust_normalize(shaded, 1, 99)


# =============================================================================
# 4. SYNTHETIC GIS-CONSTRAINED VELOCITY MODEL
# =============================================================================


def polygon_mask(geom, lon_grid: np.ndarray, lat_grid: np.ndarray) -> np.ndarray:
    """Return boolean mask for points inside a polygon geometry."""
    if contains_xy is not None:
        return np.asarray(contains_xy(geom, lon_grid, lat_grid), dtype=bool)

    # Fallback for older shapely versions. Slower, but robust enough for this grid.
    from shapely.geometry import Point
    flat = [geom.contains(Point(x, y)) for x, y in zip(lon_grid.ravel(), lat_grid.ravel())]
    return np.asarray(flat, dtype=bool).reshape(lon_grid.shape)


def extract_fault_vertices(faults_original_crs: gpd.GeoDataFrame) -> np.ndarray:
    """Extract vertices from fault LineString/MultiLineString geometries."""
    coords = []
    for geom in faults_original_crs.geometry:
        if geom is None or geom.is_empty:
            continue
        if isinstance(geom, LineString):
            coords.extend(list(geom.coords))
        elif isinstance(geom, MultiLineString):
            for line in geom.geoms:
                coords.extend(list(line.coords))
        else:
            try:
                # In case a geometry collection appears.
                for part in geom.geoms:
                    if isinstance(part, LineString):
                        coords.extend(list(part.coords))
            except Exception:
                pass
    if not coords:
        raise ValueError("No fault vertices could be extracted from the fault shapefile.")
    return np.asarray(coords, dtype=float)


def build_fault_distance_fields(
    lon_grid: np.ndarray,
    lat_grid: np.ndarray,
    raster_crs: str,
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Build nearest-fault distance and north/south side fields.

    The shapefile is originally projected, so distances are computed in its own CRS
    where coordinates are metric. This is more stable than computing distances in degrees.
    """
    faults_metric = clean_geometries(gpd.read_file(FAULTS_PATH))
    if faults_metric.crs is None:
        raise ValueError("Fault shapefile has no CRS. Cannot compute metric distances.")

    fault_vertices = extract_fault_vertices(faults_metric)
    tree = cKDTree(fault_vertices[:, :2])

    transformer = Transformer.from_crs(raster_crs, faults_metric.crs, always_xy=True)
    gx, gy = transformer.transform(lon_grid.ravel(), lat_grid.ravel())
    points_metric = np.column_stack([gx, gy])

    dist, nearest_idx = tree.query(points_metric, k=1)
    nearest_y = fault_vertices[nearest_idx, 1]
    side = np.tanh((points_metric[:, 1] - nearest_y) / 9000.0)

    return dist.reshape(lon_grid.shape), side.reshape(lon_grid.shape)


def make_synthetic_depth_slice(
    depth_km: float,
    idx: int,
    lon_grid: np.ndarray,
    lat_grid: np.ndarray,
    basins: gpd.GeoDataFrame,
    bathy_grid: np.ndarray,
    fault_distance_m: np.ndarray,
    fault_side: np.ndarray,
) -> np.ndarray:
    """Create one GIS-constrained synthetic δVp/Vp depth slice."""
    rng = np.random.default_rng(1200 + idx)

    # Smooth deterministic texture avoids a flat artificial appearance.
    dv = gaussian_filter(rng.normal(0.0, 0.45, lon_grid.shape), sigma=7)

    # Depth attenuation: shallow basin anomalies strong; deep lower-crustal anomaly emerges.
    shallow_scale = np.exp(-depth_km / 13.0)
    deep_scale = 1.0 / (1.0 + np.exp(-(depth_km - 11.0) / 2.2))

    # Use bathymetry/TIFF as a weak structural control. Deep marine depressions strengthen low-Vp zones.
    depth_norm = robust_normalize(-bathy_grid, 5, 96)  # deeper water = larger value
    depth_norm = gaussian_filter(np.nan_to_num(depth_norm, nan=0.0), sigma=3)
    dv -= (2.1 * shallow_scale) * depth_norm

    # GIS basin/high polygons control anomaly locations.
    for _, row in basins.iterrows():
        code = str(row.get("Code", "")).strip()
        geom = row.geometry
        if geom is None or geom.is_empty:
            continue

        mask = polygon_mask(geom, lon_grid, lat_grid).astype(float)
        smooth_mask = gaussian_filter(mask, sigma=4)

        if code in LOW_VP_CODES:
            # Basin low-Vp anomaly caused by sedimentary fill and structural depressions.
            strength = {
                "TB": 5.1,
                "CB": 5.8,
                "KB": 4.4,
                "CiB": 5.3,
            }.get(code, 4.5)
            dv -= strength * shallow_scale * smooth_mask
        elif code in HIGH_VP_CODES:
            # Structural highs have relatively positive perturbation.
            strength = {"CH": 4.0, "WH": 3.7}.get(code, 3.0)
            dv += strength * shallow_scale * smooth_mask

    # Fault-parallel velocity contrast; strongest near mapped fault vertices.
    fault_zone = np.exp(-((fault_distance_m / 16000.0) ** 2))
    dv += (2.3 * shallow_scale) * fault_side * fault_zone

    # Additional deep positive domain that appears at 12–18 km.
    center_x = 0.5 * (lon_grid.min() + lon_grid.max())
    center_y = 0.5 * (lat_grid.min() + lat_grid.max())
    regional_deep_high = np.exp(-(((lon_grid - center_x) ** 2) / 0.85 + ((lat_grid - center_y) ** 2) / 0.025))
    dv += 3.4 * deep_scale * regional_deep_high

    # Light smoothing produces tomography-like continuity.
    dv = gaussian_filter(dv, sigma=2.3)
    return np.clip(dv, VP_VMIN, VP_VMAX)


# =============================================================================
# 5. DRAWING HELPERS
# =============================================================================


def _line_segments_from_geometries(gdf: gpd.GeoDataFrame):
    """Convert LineString/MultiLineString geometries to coordinate arrays for LineCollection."""
    segments = []
    for geom in gdf.geometry:
        if geom is None or geom.is_empty:
            continue
        if isinstance(geom, LineString):
            segments.append(np.asarray(geom.coords))
        elif isinstance(geom, MultiLineString):
            for line in geom.geoms:
                segments.append(np.asarray(line.coords))
        else:
            try:
                for part in geom.geoms:
                    if isinstance(part, LineString):
                        segments.append(np.asarray(part.coords))
            except Exception:
                pass
    return segments


def plot_faults_with_halo(ax, faults: gpd.GeoDataFrame) -> None:
    """Plot faults with a bright halo so they remain readable over color fields."""
    segments = _line_segments_from_geometries(faults)
    if not segments:
        return
    ax.add_collection(LineCollection(segments, colors="white", linewidths=3.7, alpha=0.95, zorder=8))
    ax.add_collection(LineCollection(segments, colors="black", linewidths=2.15, alpha=0.95, zorder=9))
    ax.add_collection(LineCollection(segments, colors="#FF0033", linewidths=0.95, alpha=0.95, zorder=10))


def _polygon_parts(geom):
    """Yield polygon geometries from Polygon or MultiPolygon-like objects."""
    if geom is None or geom.is_empty:
        return
    if geom.geom_type == "Polygon":
        yield geom
    elif geom.geom_type == "MultiPolygon":
        for part in geom.geoms:
            yield part
    else:
        try:
            for part in geom.geoms:
                if part.geom_type == "Polygon":
                    yield part
        except Exception:
            return


def plot_basins(ax, basins: gpd.GeoDataFrame) -> None:
    """Plot sub-basin polygon boundaries and transparent fills without slow GeoPandas plotting."""
    for _, row in basins.iterrows():
        code = str(row.get("Code", "")).strip()
        color = BASIN_EDGE_COLORS.get(code, DEFAULT_BASIN_COLOR)
        for poly in _polygon_parts(row.geometry):
            exterior = np.asarray(poly.exterior.coords)
            ax.add_patch(
                MplPolygon(
                    exterior,
                    closed=True,
                    facecolor=color,
                    edgecolor=color,
                    linewidth=2.0,
                    alpha=0.18,
                    zorder=7,
                )
            )
            # Draw holes if present, using a thin white outline.
            for interior in poly.interiors:
                hole = np.asarray(interior.coords)
                ax.plot(hole[:, 0], hole[:, 1], color="white", linewidth=0.8, alpha=0.6, zorder=8)

def add_decluttered_basin_labels(ax, basins: gpd.GeoDataFrame, panel_idx: int) -> None:
    """
    Place basin labels away from polygons using boxed annotations and leader lines.

    The offsets are staggered so labels do not sit on top of the basin polygons, fault
    lines, or color anomalies. This is the key readability improvement.
    """
    for _, row in basins.iterrows():
        code = str(row.get("Code", "")).strip()
        if not code:
            continue
        name = str(row.get("Name", code)).strip()
        point = row.geometry.representative_point()
        x, y = point.x, point.y
        dx, dy = BASIN_LABEL_OFFSETS.get(code, (0.04, 0.04))
        tx, ty = x + dx, y + dy
        color = BASIN_EDGE_COLORS.get(code, DEFAULT_BASIN_COLOR)

        # Keep labels just inside each panel extent.
        xlim = ax.get_xlim()
        ylim = ax.get_ylim()
        tx = min(max(tx, xlim[0] + 0.035), xlim[1] - 0.035)
        ty = min(max(ty, ylim[0] + 0.020), ylim[1] - 0.020)

        ax.annotate(
            code,
            xy=(x, y),
            xytext=(tx, ty),
            fontsize=10.5,
            fontweight="bold",
            color="black",
            ha="center",
            va="center",
            bbox=dict(boxstyle="round,pad=0.25", fc="white", ec=color, lw=1.5, alpha=0.94),
            arrowprops=dict(arrowstyle="-", lw=1.2, color=color, alpha=0.95, shrinkA=4, shrinkB=4),
            zorder=20,
        )

        # Optional very small full-name tooltip-style text for the first panel only.
        if panel_idx == 0:
            ax.text(
                tx,
                ty - 0.018,
                name.replace("Westren", "Western"),
                fontsize=5.8,
                color="black",
                ha="center",
                va="top",
                bbox=dict(boxstyle="round,pad=0.15", fc="white", ec="none", alpha=0.74),
                zorder=19,
            )


def add_north_arrow(ax, xy=(0.055, 0.83), size=0.085) -> None:
    """Add a simple north arrow in axes coordinates."""
    ax.annotate(
        "N",
        xy=(xy[0], xy[1] + size),
        xytext=(xy[0], xy[1]),
        xycoords="axes fraction",
        textcoords="axes fraction",
        ha="center",
        va="center",
        fontsize=11,
        fontweight="bold",
        arrowprops=dict(arrowstyle="-|>", lw=1.5, color="black"),
        zorder=30,
    )


def format_axes(ax, extent: Tuple[float, float, float, float], idx: int) -> None:
    """Format one map axis."""
    xmin, xmax, ymin, ymax = extent
    ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymin, ymax)
    ax.set_aspect(MAP_VERTICAL_STRETCH, adjustable="box")

    if SHOW_GRID:
        ax.grid(color="white", linestyle="--", linewidth=0.45, alpha=0.55, zorder=3)
        ax.grid(which="minor", color="white", linestyle=":", linewidth=0.25, alpha=0.30)
        ax.minorticks_on()

    ax.tick_params(axis="both", which="major", labelsize=9, width=1.0, length=4)
    ax.tick_params(axis="both", which="minor", labelsize=7, width=0.7, length=2)

    if idx >= 2:
        ax.set_xlabel("Longitude (°E)", fontsize=11, fontweight="bold")
    else:
        ax.set_xlabel("")
        ax.tick_params(labelbottom=False)
    if idx % 2 == 0:
        ax.set_ylabel("Latitude (°N)", fontsize=11, fontweight="bold")
    else:
        ax.set_ylabel("")
        ax.tick_params(labelleft=False)

    for spine in ax.spines.values():
        spine.set_linewidth(1.2)
        spine.set_color("black")


def add_panel_title(ax, letter: str, title: str, explanation: str) -> None:
    """Add separated panel title and depth explanation."""
    ax.text(
        0.012,
        1.030,
        f"({letter}) {title}",
        transform=ax.transAxes,
        fontsize=13,
        fontweight="bold",
        ha="left",
        va="bottom",
        color="black",
        bbox=dict(boxstyle="round,pad=0.22", fc="white", ec="black", lw=0.7, alpha=0.93),
        zorder=40,
    )
    ax.text(
        0.012,
        0.972,
        explanation,
        transform=ax.transAxes,
        fontsize=8.4,
        fontstyle="italic",
        ha="left",
        va="top",
        color="#333333",
        bbox=dict(boxstyle="round,pad=0.18", fc="white", ec="none", alpha=0.78),
        zorder=40,
    )


# =============================================================================
# 6. MAIN FIGURE FUNCTION
# =============================================================================


def create_figure(save_path: str = OUTPUT_PNG, dpi: int = DPI) -> None:
    faults, basins, raster_crs = load_gis_layers()
    extent = compute_inner_extent(faults, basins)
    xmin, xmax, ymin, ymax = extent

    # Model grid.
    lon_g = np.linspace(xmin, xmax, NX)
    lat_g = np.linspace(ymin, ymax, NY)
    LG, LT = np.meshgrid(lon_g, lat_g)

    # Raster background resampled to model grid.
    bathy_grid = read_raster_background(extent, out_shape=(NY, NX))
    hillshade = make_hillshade(bathy_grid)

    # Metric distance to nearest mapped fault.
    fault_distance_m, fault_side = build_fault_distance_fields(LG, LT, raster_crs)

    # Create the four synthetic GIS-constrained depth slices.
    slices = [
        make_synthetic_depth_slice(depth, i, LG, LT, basins, bathy_grid, fault_distance_m, fault_side)
        for i, depth in enumerate(DEPTHS_KM)
    ]

    fig, axes = plt.subplots(2, 2, figsize=FIGSIZE, dpi=dpi)
    axes_flat = axes.flat

    last_im = None
    for idx, (ax, depth, title, dv) in enumerate(zip(axes_flat, DEPTHS_KM, DEPTH_TITLES, slices)):
        # Hillshade behind the velocity field adds real TIFF morphology without dominating the colors.
        if SHOW_RASTER_HILLSHADE:
            ax.imshow(
                hillshade,
                extent=[xmin, xmax, ymin, ymax],
                origin="lower",
                cmap="Greys",
                alpha=0.32,
                interpolation="bilinear",
                zorder=1,
            )

        last_im = ax.imshow(
            dv,
            extent=[xmin, xmax, ymin, ymax],
            origin="lower",
            cmap=CMAP_VEL,
            vmin=VP_VMIN,
            vmax=VP_VMAX,
            alpha=0.92,
            interpolation="bilinear",
            zorder=2,
        )

        # Bathymetric contour context from the TIFF. Values are real raster values.
        if SHOW_BATHY_CONTOURS:
            finite = np.isfinite(bathy_grid)
            if finite.any():
                levels = np.nanpercentile(bathy_grid[finite], [10, 25, 40, 55, 70, 85])
                ax.contour(
                    LG,
                    LT,
                    bathy_grid,
                    levels=np.unique(np.round(levels, 1)),
                    colors="#3A3A3A",
                    linewidths=0.35,
                    alpha=0.42,
                    zorder=5,
                )

        plot_basins(ax, basins)
        plot_faults_with_halo(ax, faults)

        if SHOW_BASIN_LABELS_IN_ALL_PANELS or idx == 0:
            add_decluttered_basin_labels(ax, basins, idx)

        # Add north arrow only once to reduce panel clutter.
        if idx == 0:
            add_north_arrow(ax)

        add_panel_title(ax, chr(97 + idx), title, PANEL_EXPLANATIONS[depth])
        format_axes(ax, extent, idx)

    # Shared colorbar, separated from panels.
    cbar_ax = fig.add_axes([0.925, 0.172, 0.022, 0.665])
    cbar = fig.colorbar(last_im, cax=cbar_ax)
    cbar.set_label("δVp/Vp (%)", fontsize=12.5, fontweight="bold", labelpad=12)
    cbar.ax.tick_params(labelsize=10)
    cbar.outline.set_linewidth(1.1)

    # Separate legend outside map area to prevent overlap with the panels.
    legend_handles = [
        Line2D([0], [0], color="black", lw=3.0, label="Mapped Marmara fault segments"),
        Line2D([0], [0], color="#FF0033", lw=1.5, label="Fault highlight / surface trace"),
        Patch(facecolor="#08306B", edgecolor="black", alpha=0.8, label="Negative velocity perturbation"),
        Patch(facecolor="#CB181D", edgecolor="black", alpha=0.8, label="Positive velocity perturbation"),
    ]
    for code in ["TB", "CB", "KB", "CiB", "CH", "WH"]:
        legend_handles.append(
            Patch(
                facecolor=BASIN_EDGE_COLORS.get(code, DEFAULT_BASIN_COLOR),
                edgecolor=BASIN_EDGE_COLORS.get(code, DEFAULT_BASIN_COLOR),
                alpha=0.35,
                label=f"{code} sub-basin/high boundary",
            )
        )

    fig.legend(
        handles=legend_handles,
        loc="lower center",
        bbox_to_anchor=(0.50, 0.018),
        ncol=5,
        fontsize=9.2,
        frameon=True,
        fancybox=True,
        framealpha=0.96,
        edgecolor="black",
    )

    # Title and subtitle.
    fig.suptitle(
        "GIS-Constrained Synthetic P-wave Velocity Perturbation Maps at Multiple Depths",
        fontsize=18,
        fontweight="bold",
        y=0.973,
    )
    fig.text(
        0.5,
        0.943,
        "Real Marmara TIFF bathymetry/elevation, verified fault shapefile, and sub-basin polygon shapefile used for spatial control",
        ha="center",
        va="center",
        fontsize=11.2,
        color="#333333",
    )

    # Leave clear room for external legend and colorbar.
    fig.subplots_adjust(left=0.055, right=0.905, top=0.905, bottom=0.120, hspace=0.125, wspace=0.075)

    fig.savefig(save_path, dpi=dpi, bbox_inches="tight", facecolor="white")
    if OUTPUT_PDF:
        fig.savefig(OUTPUT_PDF, dpi=dpi, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"Saved → {save_path}")
    if OUTPUT_PDF:
        print(f"Saved → {OUTPUT_PDF}")


if __name__ == "__main__":
    create_figure()
