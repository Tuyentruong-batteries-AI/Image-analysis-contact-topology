"""
eds_topology.py — contact-topology analysis of SEM–EDS C/F maps
=================================================================

What it does
------------
PVDF contains carbon, so a C map alone cannot show where the CNTs are.
Plotted as C signal (x) versus F signal (y), the pixels of a CNT/PVDF map
fill a triangle with three corners:

    background  (bare Si, or dark pores/shadows in a film)
    CNT         (C without F)
    PVDF        (C with F)

Every pixel is described by where it sits in that triangle ("unmixing").
Its CNT fraction  x = CNT / (CNT + PVDF)  sorts it into
CNT-rich (x > 0.6), mixed, or PVDF-rich (x < 0.4). Topology metrics are
computed from those classes, and every metric is recomputed over ten
analysis settings to show how sensitive it is.

Input   : Oxford AZtec "Site Report – Layered Image and Maps" (.docx)
Outputs : Excel workbook + CSV (data you can re-plot in Origin/Excel),
          TIFF class maps (ImageJ), PNG/SVG/PDF figures (SVG/PDF are
          editable in PowerPoint, Illustrator or Inkscape).

Requires: numpy, scipy, scikit-image, matplotlib, pandas, openpyxl, pillow
(all included in Anaconda and Google Colab).

Command line (batch):
    python eds_topology.py samples.csv  -o results
where samples.csv has columns: file, sample, group, elements, background
(elements = map order in the report, e.g. "C;F;O"; background = low or si).
"""
from __future__ import annotations

import io
import os
import re
import sys
import zipfile

import numpy as np
from PIL import Image
from scipy import ndimage as ndi
from scipy import stats
from skimage.filters import apply_hysteresis_threshold, threshold_otsu
from skimage.graph import MCP_Geometric
from skimage.morphology import remove_small_objects, skeletonize

NATIVE_WIDTH = 512          # AZtec map width in pixels (report images are upsampled)

# Default analysis settings and the ten variants used for the sensitivity ranges
DEFAULTS = dict(sigma=1.0, t_lo=0.4, t_hi=0.6, lo_frac=0.55, endpct=(1, 99))
SETTINGS = [
    {},                                  # baseline
    dict(sigma=0.7), dict(sigma=1.5),    # smoothing
    dict(t_lo=0.35, t_hi=0.65), dict(t_lo=0.45, t_hi=0.55),   # class thresholds
    dict(lo_frac=0.4), dict(lo_frac=0.75),                    # material mask
    dict(endpct=(0.5, 99.5)), dict(endpct=(3, 97)), dict(endpct=(5, 95)),  # corners
]

def _drop_small(mask, n):
    """Remove connected objects of n pixels or fewer (same result on old and new scikit-image)."""
    try:
        return remove_small_objects(mask, max_size=n)        # scikit-image ≥ 0.26
    except TypeError:
        return remove_small_objects(mask, min_size=n + 1)    # older versions


CLASS_NAMES = {0: 'low signal', 1: 'PVDF-rich', 2: 'mixed', 3: 'CNT-rich'}
CLASS_COLORS = {0: '#3a3a38', 1: '#1baf7a', 2: '#2a78d6', 3: '#eb6834'}

METRIC_INFO = {
    # key: (label, unit, what it means)
    'corr_CF': ('C–F correlation within material', 'Pearson r',
                'Negative: CNT and PVDF occupy separate regions. Positive: polymer mixed through the CNTs.'),
    'F_in_CNT': ('F in CNT-rich domains', 'fraction of PVDF-rich F',
                 '≤0: CNT domains carry no F above background (polymer-free). >0: polymer present in CNT domains.'),
    'low_signal_frac': ('Low-signal (pore/shadow) area', 'fraction of field',
                        'Open pores or shadowed depressions.'),
    'PVDF_largest_network': ('PVDF-rich area in largest connected network', 'fraction',
                             'Continuity of the polymer network.'),
    'PVDF_within_5px': ('PVDF-rich area within 5 px of a CNT-rich domain', 'fraction',
                        'Pixel proximity of polymer to CNT; high values can also mean polymer blended inside CNT agglomerates.'),
    'PVDF_median_dist_px': ('Median path distance, PVDF-rich → CNT-rich', 'map pixels',
                            'Distance measured inside the material (never across pores or substrate).'),
    'area_CNT': ('CNT-rich area', 'fraction of material', 'Relative amount of CNT-dominated pixels.'),
    'area_mixed': ('Mixed area', 'fraction of material', ''),
    'area_PVDF': ('PVDF-rich area', 'fraction of material', ''),
    'CNT_perimeter_bordering_material': ('CNT-domain edge bordering material', 'fraction',
                                         'Particles on Si only: edge facing material rather than bare substrate.'),
    'free_CNT_area': ('CNT-rich area in objects with no PVDF', 'fraction',
                      'Particles on Si only: free CNT agglomerates.'),
}


# ----------------------------------------------------------------------------- reading
def read_aztec_report(path):
    """Return the images of an AZtec site report in document order.

    Typical order: [layered image, map 1, map 2, ..., map sum spectrum].
    """
    z = zipfile.ZipFile(path)
    rels = z.read('word/_rels/document.xml.rels').decode('utf8')
    doc = z.read('word/document.xml').decode('utf8')
    rid2target = {}
    for el in re.findall(r'<Relationship\b[^>]*>', rels):
        rid = re.search(r'\bId="([^"]+)"', el)
        tgt = re.search(r'\bTarget="([^"]+)"', el)
        if rid and tgt:
            rid2target[rid.group(1)] = tgt.group(1)
    images = []
    for rid in re.findall(r'r:embed="([^"]+)"', doc):
        t = rid2target[rid]
        name = t.lstrip('/') if t.startswith('/') else 'word/' + t
        images.append(Image.open(io.BytesIO(z.read(name))).convert('RGB'))
    return images


def load_maps(path, elements, native_width=NATIVE_WIDTH):
    """Load the element maps of one report as 2-D float arrays on the native map grid.

    elements: element of each map image in report order, e.g. ['C', 'F', 'O'].
              Use 'X' for maps to ignore (e.g. misidentified Re/Os/Pt).
    """
    imgs = read_aztec_report(path)
    map_imgs = imgs[1:-1]                     # drop layered image and spectrum
    if len(map_imgs) != len(elements):
        raise ValueError(f'{os.path.basename(path)}: report has {len(map_imgs)} maps but '
                         f'{len(elements)} elements were given ({elements}).')
    arrs = [np.asarray(im).astype(float).max(axis=-1) for im in map_imgs]   # colour → intensity
    union = np.max(np.stack(arrs), axis=0) > 0
    ys, xs = np.where(union)
    y0, y1, x0, x1 = ys.min(), ys.max() + 1, xs.min(), xs.max() + 1         # crop black frame
    h = int(round(native_width * (y1 - y0) / (x1 - x0)))
    out = {}
    for el, a in zip(elements, arrs):
        if el.upper() == 'X':
            continue
        crop = Image.fromarray(a[y0:y1, x0:x1].astype(np.float32), mode='F')
        out[el] = np.asarray(crop.resize((native_width, h), Image.BOX))     # box-average to native grid
    return out


def layered_image(path):
    """The AZtec layered (overlay) image, cropped to the data area."""
    im = np.asarray(read_aztec_report(path)[0])
    m = im.sum(-1) > 30
    ys, xs = np.where(m)
    return im[ys.min():ys.max() + 1, xs.min():xs.max() + 1]


# ----------------------------------------------------------------------------- unmixing
def find_corners(C, F, background='low', Si=None, bgpct=3, endpct=(1, 99)):
    """Locate the three corners (end-members) of the C–F triangle.

    background='low': background = median of the lowest `bgpct` % of total signal (films).
    background='si' : background = median of the brightest 30 % of the Si map (particles on Si).
    The CNT and PVDF corners lie along the two extreme directions of the scatter
    as seen from the background corner (lowest and highest angle), at the 97th
    percentile of distance.
    Returns (S, CNT, PVDF), each an array [C, F].
    """
    c, f = C.ravel(), F.ravel()
    if background == 'si':
        if Si is None:
            raise ValueError("background='si' needs a Si map")
        sub = Si.ravel() > np.percentile(Si, 70)
    else:
        cn, fn = c / np.percentile(c, 99.5), f / np.percentile(f, 99.5)
        sub = (cn + fn) < np.percentile(cn + fn, bgpct)
    S = np.array([np.median(c[sub]), np.median(f[sub])])
    v = np.stack([c - S[0], f - S[1]], 1)
    mag = np.hypot(v[:, 0], v[:, 1])
    ang = np.arctan2(v[:, 1], v[:, 0])
    big = mag > 0.35 * np.percentile(mag, 99.5)

    def corner(a0):
        sel = big & (np.abs(ang - a0) < 0.06)
        return S + np.percentile(mag[sel], 97) * np.array([np.cos(a0), np.sin(a0)])

    return S, corner(np.percentile(ang[big], endpct[0])), corner(np.percentile(ang[big], endpct[1]))


def unmix(C, F, S, CNT, PVDF):
    """Barycentric coordinates of every pixel in the triangle (background, CNT, PVDF).

    Returns (f_CNT, f_PVDF): abundance of each end-member, clipped to ≥ 0 and
    renormalised so that f_CNT + f_PVDF ≤ 1 (the rest is background).
    """
    T = np.array([[CNT[0] - S[0], PVDF[0] - S[0]],
                  [CNT[1] - S[1], PVDF[1] - S[1]]])
    w = np.linalg.inv(T) @ np.stack([C - S[0], F - S[1]]).reshape(2, -1)
    w = np.clip(w.reshape(2, *C.shape), 0, None)
    s = np.maximum(w.sum(0), 1)
    return w[0] / s, w[1] / s


def segment(maps, sigma=1.0, t_lo=0.4, t_hi=0.6, lo_frac=0.55, endpct=(1, 99),
            background='low', bgpct=3, mask_foreign=None):
    """Smooth, unmix and classify one sample. Returns a dict of arrays."""
    sm = {k: ndi.gaussian_filter(v, sigma) for k, v in maps.items()}
    C, F = sm['C'], sm['F']
    S, CNTc, PVDFc = find_corners(C, F, background, sm.get('Si'), bgpct, endpct)
    fc, fp = unmix(C, F, S, CNTc, PVDFc)
    mat = fc + fp
    th = threshold_otsu(mat)
    M = _drop_small(apply_hysteresis_threshold(mat, lo_frac * th, th), 12)
    foreign = np.zeros_like(M)
    if mask_foreign is None:
        mask_foreign = background == 'si'
    if mask_foreign and 'O' in sm:     # O-rich objects that are neither C nor F (dust, salt)
        O = sm['O']
        cand = _drop_small(ndi.binary_opening((O - np.median(O)) / O.std() > 3), 40)
        lab, n = ndi.label(cand)
        for i in range(1, n + 1):
            if np.median(mat[lab == i]) < 0.5:
                foreign |= lab == i
        foreign = ndi.binary_dilation(foreign, iterations=2)
    M &= ~foreign
    x = fc / np.maximum(mat, 1e-6)
    CNT = _drop_small(M & (x > t_hi), 4)
    PV = M & (x < t_lo)
    MIX = M & ~CNT & ~PV
    labels = np.zeros(M.shape, np.uint8)
    labels[PV], labels[MIX], labels[CNT] = 1, 2, 3
    return dict(C=C, F=F, x=x, M=M, CNT=CNT, PV=PV, MIX=MIX, labels=labels,
                foreign=foreign, corners=(S, CNTc, PVDFc))


# ----------------------------------------------------------------------------- metrics
def path_distance_to(target, allowed):
    """Distance (px) from every pixel to the nearest `target` pixel, moving only through `allowed`."""
    cost = np.where(allowed, 1.0, np.inf)
    d, _ = MCP_Geometric(cost).find_costs(list(zip(*np.where(target))))
    return d


def metrics(g, background='low'):
    """Topology metrics for one segmentation (dict from `segment`)."""
    M, CNT, PV, MIX = g['M'], g['CNT'], g['PV'], g['MIX']
    n = M.sum()
    r = dict(low_signal_frac=1 - M.mean(),
             area_CNT=CNT.sum() / n, area_mixed=MIX.sum() / n, area_PVDF=PV.sum() / n,
             corr_CF=float(np.corrcoef(g['C'][M], g['F'][M])[0, 1]))
    lab, nl = ndi.label(PV, np.ones((3, 3)))
    r['PVDF_largest_network'] = np.bincount(lab[lab > 0]).max() / PV.sum() if nl else np.nan
    d = path_distance_to(CNT, M)
    dp = d[PV]
    r['PVDF_within_5px'] = float(np.mean(dp <= 5))
    r['PVDF_median_dist_px'] = float(np.median(dp[np.isfinite(dp)]))
    if background == 'si':
        tot = mat_edge = 0
        for sh, edge in [((0, 1), (slice(None), 0)), ((0, -1), (slice(None), -1)),
                         ((1, 0), (0, slice(None))), ((-1, 0), (-1, slice(None)))]:
            e = CNT & ~np.roll(CNT, sh, (0, 1))
            e[edge] = False
            tot += e.sum()
            mat_edge += (e & np.roll(M, sh, (0, 1))).sum()
        r['CNT_perimeter_bordering_material'] = mat_edge / tot
        mlab, _ = ndi.label(M, np.ones((3, 3)))
        clab, nc = ndi.label(CNT, np.ones((3, 3)))
        sizes = ndi.sum(np.ones_like(clab), clab, range(1, nc + 1))
        comp = [np.bincount(mlab[clab == i]).argmax() for i in range(1, nc + 1)]
        has_pv = np.array([(PV & (mlab == c)).any() for c in comp])
        r['free_CNT_area'] = sizes[~has_pv].sum() / sizes.sum()
    return r, d


def F_in_CNT(maps, background='low'):
    """F signal at the CNT corner relative to the PVDF corner (both above background).

    Baseline = median over smoothing σ ∈ {0.7, 1, 1.5} × corner percentiles
    {0.5/99.5, 1/99, 3/97}; range = min–max over those nine combinations.
    """
    vals = []
    for sig in (0.7, 1.0, 1.5):
        C = ndi.gaussian_filter(maps['C'], sig)
        F = ndi.gaussian_filter(maps['F'], sig)
        Si = ndi.gaussian_filter(maps['Si'], sig) if 'Si' in maps else None
        for ep in ((0.5, 99.5), (1, 99), (3, 97)):
            S, cnt, pv = find_corners(C, F, background, Si, 3, ep)
            vals.append((cnt[1] - S[1]) / (pv[1] - S[1]))
    return float(np.median(vals)), float(min(vals)), float(max(vals))


def analyze_sample(path, elements, sample=None, group='', background='low', verbose=True):
    """Full analysis of one report: baseline segmentation + metrics with sensitivity ranges."""
    sample = sample or os.path.splitext(os.path.basename(path))[0]
    maps = load_maps(path, elements)
    runs, base, base_d = [], None, None
    for k, s in enumerate(SETTINGS):
        g = segment(maps, **{**DEFAULTS, **s}, background=background)
        r, d = metrics(g, background)
        runs.append(r)
        if k == 0:
            base, base_d = g, d
    keys = runs[0].keys()
    m = {k: (float(runs[0][k]), float(np.nanmin([r[k] for r in runs])),
             float(np.nanmax([r[k] for r in runs]))) for k in keys}
    m['F_in_CNT'] = F_in_CNT(maps, background)
    if verbose:
        print(f'{sample:>12s} [{group}]  r = {m["corr_CF"][0]:+.2f}   F in CNT = {m["F_in_CNT"][0]:+.2f}   '
              f'low-signal = {100 * m["low_signal_frac"][0]:.0f}%   PVDF network = {100 * m["PVDF_largest_network"][0]:.0f}%')
    return dict(sample=sample, group=group, file=os.path.basename(path), elements=list(elements),
                background=background, maps=maps, seg=base, dist=base_d, metrics=m,
                layered=layered_image(path))


def analyze_batch(table, folder='.', verbose=True):
    """table: list of dicts with keys file, sample, group, elements (list or 'C;F;O'), background."""
    out = []
    for row in table:
        els = row['elements']
        els = [e.strip() for e in re.split(r'[;,\s]+', els) if e.strip()] if isinstance(els, str) else els
        out.append(analyze_sample(os.path.join(folder, row['file']), els, row['sample'],
                                  row.get('group', ''), row.get('background', 'low') or 'low', verbose))
    return out


def group_stats(results, metric_keys=None):
    """Mean ± SD per group and Welch's t-test when there are exactly two groups (n ≥ 2 each)."""
    import pandas as pd
    groups = [g for g in dict.fromkeys(r['group'] for r in results) if g]
    metric_keys = metric_keys or [k for k in METRIC_INFO if k in results[0]['metrics']]
    rows = []
    for k in metric_keys:
        row = {'metric': k, 'label': METRIC_INFO[k][0], 'unit': METRIC_INFO[k][1]}
        vals = {}
        for g in groups:
            v = np.array([r['metrics'][k][0] for r in results if r['group'] == g])
            vals[g] = v
            row[f'{g} n'] = len(v)
            row[f'{g} mean'] = v.mean()
            row[f'{g} SD'] = v.std(ddof=1) if len(v) > 1 else np.nan
            row[f'{g} min'] = v.min()
            row[f'{g} max'] = v.max()
        if len(groups) == 2 and all(len(v) >= 2 for v in vals.values()):
            a, b = vals[groups[0]], vals[groups[1]]
            row['Welch t-test p'] = stats.ttest_ind(a, b, equal_var=False).pvalue
            row['Mann-Whitney p (exact)'] = stats.mannwhitneyu(a, b, alternative='two-sided', method='exact').pvalue
            lo_a = min(r['metrics'][k][1] for r in results if r['group'] == groups[0])
            hi_a = max(r['metrics'][k][2] for r in results if r['group'] == groups[0])
            lo_b = min(r['metrics'][k][1] for r in results if r['group'] == groups[1])
            hi_b = max(r['metrics'][k][2] for r in results if r['group'] == groups[1])
            row['groups overlap across settings?'] = 'yes' if not (hi_a < lo_b or hi_b < lo_a) else 'no'
        rows.append(row)
    return pd.DataFrame(rows)


# ----------------------------------------------------------------------------- export
def _style():
    import matplotlib
    matplotlib.rcParams.update({'svg.fonttype': 'none', 'pdf.fonttype': 42, 'ps.fonttype': 42,
                                'font.family': 'DejaVu Sans', 'font.size': 8,
                                'axes.spines.top': False, 'axes.spines.right': False})


def class_rgb(labels):
    from matplotlib.colors import to_rgb
    rgb = np.zeros(labels.shape + (3,))
    for k, c in CLASS_COLORS.items():
        rgb[labels == k] = to_rgb(c)
    return rgb


def cf_overlay(maps):
    n = lambda a: np.clip(ndi.gaussian_filter(a, 0.8) / np.percentile(ndi.gaussian_filter(a, 0.8), 99.5), 0, 1)
    return np.dstack([n(maps['C']), n(maps['F']), np.zeros_like(maps['C'])])


def plot_triangle(res, ax=None):
    """C–F scatter (2-D histogram) with the three corners marked."""
    import matplotlib.pyplot as plt
    from matplotlib.colors import LogNorm
    g = res['seg']
    ax = ax or plt.gca()
    h = ax.hist2d(g['C'].ravel(), g['F'].ravel(), bins=120, cmap='Greys', norm=LogNorm())
    h[3].set_rasterized(True)          # keeps SVG/PDF small; axes and text stay editable
    S, CNT, PV = g['corners']
    ax.plot([S[0], CNT[0], PV[0], S[0]], [S[1], CNT[1], PV[1], S[1]], '-', color='#d03b3b', lw=1)
    for p_, lab in [(S, 'background'), (CNT, 'CNT'), (PV, 'PVDF')]:
        ax.plot(*p_, 'o', color='#d03b3b', ms=5)
        ax.annotate(lab, p_, xytext=(4, 4), textcoords='offset points', fontsize=7, color='#d03b3b')
    ax.set_xlabel('C signal')
    ax.set_ylabel('F signal')
    ax.set_title(res['sample'], fontsize=8, loc='left')
    return ax


def export_all(results, outdir='EDS_topology_results', scatter_points=5000, seed=0):
    """Write every result to `outdir`:

    EDS_topology_results.xlsx   all tables (one sheet each) — open in Excel, copy into Origin
    csv/                        the same tables as CSV
    maps/<sample>_*.png|tif     class map (PNG for slides, TIFF for ImageJ), CNT-fraction map (32-bit TIFF)
    data/<sample>_scatter.csv   random pixels (C, F, class) for re-plotting the triangle
    raw/<sample>_*.csv          full pixel matrices: C and F maps (raw and smoothed), CNT fraction,
                                class labels (0 low signal, 1 PVDF-rich, 2 mixed, 3 CNT-rich)
    figures/*.svg|pdf|png       editable figures
    """
    import matplotlib
    matplotlib.use('Agg', force=False)
    import matplotlib.pyplot as plt
    import pandas as pd
    _style()
    for sub in ('csv', 'maps', 'data', 'raw', 'figures'):
        os.makedirs(os.path.join(outdir, sub), exist_ok=True)
    rng = np.random.default_rng(seed)
    keys = [k for k in METRIC_INFO if k in results[0]['metrics']]

    # ---- tables
    summary = pd.DataFrame([{'sample': r['sample'], 'group': r['group'], 'file': r['file'],
                             **{k: r['metrics'][k][0] for k in keys}} for r in results])
    ranges = pd.DataFrame([{'sample': r['sample'], 'group': r['group'], 'metric': k,
                            'label': METRIC_INFO[k][0], 'baseline': r['metrics'][k][0],
                            'min': r['metrics'][k][1], 'max': r['metrics'][k][2]}
                           for r in results for k in keys])
    gstats = group_stats(results, keys)
    bins = np.linspace(0, 1, 26)
    xdist = pd.DataFrame({'CNT fraction x (bin centre)': (bins[1:] + bins[:-1]) / 2})
    dgrid = np.arange(0, 61)
    cdf = pd.DataFrame({'path distance (px)': dgrid})
    corners = []
    for r in results:
        g = r['seg']
        h, _ = np.histogram(g['x'][g['M']], bins=bins, density=True)
        xdist[r['sample']] = h
        dp = r['dist'][g['PV']]
        cdf[r['sample']] = [np.mean(dp <= v) for v in dgrid]
        S, C_, P_ = g['corners']
        corners.append({'sample': r['sample'], 'background C': S[0], 'background F': S[1],
                        'CNT corner C': C_[0], 'CNT corner F': C_[1],
                        'PVDF corner C': P_[0], 'PVDF corner F': P_[1],
                        'elements (report order)': ';'.join(r['elements']), 'background mode': r['background']})
    corners = pd.DataFrame(corners)
    defs = pd.DataFrame([{'metric': k, 'label': v[0], 'unit': v[1], 'meaning': v[2]} for k, v in METRIC_INFO.items() if k in keys])
    settings = pd.DataFrame([{'setting': i, **{**DEFAULTS, **s}} for i, s in enumerate(SETTINGS)])
    settings['endpct'] = settings['endpct'].astype(str)
    sheets = {'Summary': summary, 'Ranges': ranges, 'Group_stats': gstats, 'x_distribution': xdist,
              'Distance_CDF': cdf, 'Corners': corners, 'Metric_definitions': defs, 'Settings': settings}
    with pd.ExcelWriter(os.path.join(outdir, 'EDS_topology_results.xlsx'), engine='openpyxl') as xw:
        for name, df in sheets.items():
            df.to_excel(xw, sheet_name=name, index=False)
            ws = xw.sheets[name]
            for col in ws.columns:
                ws.column_dimensions[col[0].column_letter].width = min(48, max(10, max(len(str(c.value or '')) for c in col) + 2))
    for name, df in sheets.items():
        df.to_csv(os.path.join(outdir, 'csv', f'{name}.csv'), index=False)

    # ---- per-sample maps and scatter data
    for r in results:
        g, s = r['seg'], r['sample'].replace(' ', '_')
        Image.fromarray(g['labels']).save(os.path.join(outdir, 'maps', f'{s}_classes_labels.tif'))
        Image.fromarray(np.where(g['M'], g['x'], np.nan).astype(np.float32), mode='F').save(
            os.path.join(outdir, 'maps', f'{s}_CNT_fraction.tif'))
        plt.imsave(os.path.join(outdir, 'maps', f'{s}_classes.png'), class_rgb(g['labels']))
        plt.imsave(os.path.join(outdir, 'maps', f'{s}_CF_overlay.png'), cf_overlay(r['maps']))
        # raw pixel matrices (rows = image rows, columns = image columns) — open directly in Excel
        for nm, arr in [('C_map', r['maps']['C']), ('F_map', r['maps']['F']),
                        ('C_smoothed', g['C']), ('F_smoothed', g['F']),
                        ('CNT_fraction', np.where(g['M'], g['x'], np.nan)), ('class_labels', g['labels'])]:
            np.savetxt(os.path.join(outdir, 'raw', f'{s}_{nm}.csv'), arr, delimiter=',',
                       fmt='%d' if nm == 'class_labels' else '%.4g')
        idx = rng.choice(g['C'].size, size=min(scatter_points, g['C'].size), replace=False)
        pd.DataFrame({'C': g['C'].ravel()[idx], 'F': g['F'].ravel()[idx],
                      'class': [CLASS_NAMES[v] for v in g['labels'].ravel()[idx]]}).to_csv(
            os.path.join(outdir, 'data', f'{s}_scatter.csv'), index=False)

    # ---- figures
    def save(fig, name):
        for ext in ('svg', 'pdf', 'png'):
            fig.savefig(os.path.join(outdir, 'figures', f'{name}.{ext}'), dpi=300, bbox_inches='tight')
        plt.close(fig)

    n = len(results)
    fig, axes = plt.subplots(2, n, figsize=(2.0 * n, 3.0), gridspec_kw=dict(hspace=0.08, wspace=0.05))
    axes = np.array(axes).reshape(2, n)
    for j, r in enumerate(results):
        axes[0, j].imshow(cf_overlay(r['maps'])); axes[0, j].set_axis_off()
        axes[0, j].set_title(f"{r['group']} · {r['sample']}", fontsize=7, loc='left')
        axes[1, j].imshow(class_rgb(r['seg']['labels']), interpolation='nearest'); axes[1, j].set_axis_off()
    from matplotlib.patches import Patch
    fig.legend(handles=[Patch(color=CLASS_COLORS[k], label=CLASS_NAMES[k]) for k in (3, 2, 1, 0)],
               loc='lower center', ncol=4, frameon=False, fontsize=7, bbox_to_anchor=(0.5, -0.02))
    fig.text(0.01, 0.98, 'Top: C (red) / F (green) maps. Bottom: unmixed classes.', fontsize=7, va='top')
    save(fig, 'class_maps')

    fig, axes = plt.subplots(1, n, figsize=(2.2 * n, 2.2))
    for ax, r in zip(np.atleast_1d(axes), results):
        plot_triangle(r, ax)
    fig.tight_layout()
    save(fig, 'CF_triangles')

    show = [k for k in ('corr_CF', 'F_in_CNT', 'low_signal_frac', 'PVDF_largest_network', 'PVDF_within_5px') if k in keys]
    groups = list(dict.fromkeys(r['group'] for r in results))
    fig, axes = plt.subplots(1, len(show), figsize=(2.1 * len(show), 0.35 * n + 1.2), sharey=True)
    markers = ['o', 's', '^', 'D']
    for ax, k in zip(np.atleast_1d(axes), show):
        for i, r in enumerate(results):
            b, lo, hi = r['metrics'][k]
            gi = groups.index(r['group'])
            ax.plot([lo, hi], [i, i], color='0.4', lw=1)
            ax.plot(b, i, markers[gi % 4], color='black' if gi == 0 else 'white', mec='black', ms=5)
        if k in ('corr_CF', 'F_in_CNT'):
            ax.axvline(0, color='0.6', lw=0.8)
        import textwrap
        ax.set_title(textwrap.fill(METRIC_INFO[k][0], 26), fontsize=7, loc='left')
        ax.set_xlabel(METRIC_INFO[k][1], fontsize=7)
    ax0 = np.atleast_1d(axes)[0]
    ax0.set_yticks(range(n))
    ax0.set_yticklabels([f"{r['group']} · {r['sample']}" for r in results], fontsize=7)
    ax0.invert_yaxis()
    fig.text(0.01, -0.02, 'Symbol: baseline analysis. Bar: range over the 10 analysis settings.', fontsize=6.5)
    fig.tight_layout()
    save(fig, 'metrics_dotplot')

    fig, ax = plt.subplots(figsize=(3.2, 2.4))
    for r in results:
        ax.plot(cdf['path distance (px)'], cdf[r['sample']], label=f"{r['group']} · {r['sample']}")
    ax.set_xlabel('Path distance to nearest CNT-rich domain (map px)')
    ax.set_ylabel('Cumulative fraction of PVDF-rich area')
    ax.legend(fontsize=6, frameon=False)
    save(fig, 'distance_CDF')
    print(f'Exported to {os.path.abspath(outdir)}')
    return sheets


# ----------------------------------------------------------------------------- CLI
def read_table(csv_path):
    import pandas as pd
    df = pd.read_csv(csv_path).fillna('')
    return df.to_dict('records')


if __name__ == '__main__':
    import argparse
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('table', help='CSV with columns file, sample, group, elements, background')
    ap.add_argument('-o', '--out', default='EDS_topology_results')
    ap.add_argument('--folder', default=None, help='folder containing the .docx reports (default: folder of the CSV)')
    a = ap.parse_args()
    folder = a.folder or os.path.dirname(os.path.abspath(a.table))
    res = analyze_batch(read_table(a.table), folder)
    export_all(res, a.out)
    print(group_stats(res).to_string())
