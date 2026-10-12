# Image-analysis-contact-topology
EDS Contact-Topology Analysis: Methodology
C/F unmixing of SEM–EDS maps of CNT/PVDF interlayers — principle, algorithm, parameters, outputs, and how to run the Python toolkit
Toolkit files: eds_topology.py (analysis functions and command-line batch tool), EDS_topology_walkthrough.ipynb (step-by-step notebook), samples.csv (sample table), README.txt. Compiled 7 October 2026 for Tuyen T.T. Truong.
1. Purpose
The analysis measures how CNT and PVDF are arranged in an interlayer from routine SEM–EDS element maps: whether the polymer forms a continuous network around exposed CNT domains, or is blended through compact CNT agglomerates. It turns the maps into numbers that can be compared between samples (for example shear-assembled FS-CNT/PVDF-SDC versus cast CNT/PVDF) and tested for robustness.
2. Principle: the C–F triangle
PVDF, (CH2–CF2)n, emits both C Kα and F Kα X-rays. CNTs emit only C Kα. Pores, shadowed regions, or a bare substrate emit neither. A C map alone therefore cannot locate the CNTs, because the polymer is also bright in carbon.
If every pixel is plotted by its C signal (x) against its F signal (y), the cloud of pixels fills a triangle whose three corners are the three pure end-members: background, CNT (high C, no F), and PVDF (C and F). A pixel’s position inside the triangle gives its composition. Pixels near the CNT corner are CNT-dominated; near the PVDF corner, polymer-dominated; in between, a mixture within the pixel (≈1 µm sampling depth and width).
 
Figure 1. C–F scatter of each map (grey, log density) with the three automatically located corners (red). SDC films (left three) fill the triangle and crowd toward the corners; cast films (right three) lie along a line from background toward a C+F composition, because the polymer is mixed through the CNTs.
3. Data requirements
•	Oxford AZtec “Site Report – Layered Image and Maps” exported as .docx (the toolkit reads the map images inside the report). Raw count maps exported from AZtec can be used in the same way after small changes to load_maps.
•	C Kα and F Kα maps are required. O is optional (used to mask O-rich foreign particles on Si). Si is required only for particles on a Si wafer.
•	All maps to be compared must be acquired with the same mapping mode, accelerating voltage, field width (magnification), pixel grid, dwell time and number of frames. A map acquired in a different mode is not comparable (Site 3 was excluded for this reason).
•	Record the field width of each map: µm per pixel = field width / 512.
•	Remove misidentified elements in AZtec before exporting (e.g. Re, Os, Pt on Si, which are Si-line overlaps), or mark those maps as X in the sample table.
•	Recommended sampling: ≥3 independent samples per group and ≥2 maps per sample, at preset positions rather than hand-picked.
4. Algorithm
1.	Read the report. Images are taken in document order: layered image, element maps, spectrum. The user assigns the element of each map (the report does not store element names with the pictures).
2.	Convert to intensity and crop. Each map picture is converted to intensity (maximum of the R, G, B channels), and the black frame is removed (bounding box of non-zero pixels across all maps).
3.	Resample to the native grid. The report pictures are up-sampled about 2.75×; they are box-averaged back to 512 pixels wide (512 × 384).
4.	Smooth. Gaussian filter, σ = 1 pixel, to reduce counting noise.
5.	Locate the corners. Background: median (C, F) of the 3% of pixels with the lowest normalized C + F (films), or of the 30% brightest Si pixels (particles on Si). CNT and PVDF corners: seen from the background corner, the cloud spans a fan of angles; the CNT corner lies along the 1st-percentile angle (most C, least F) and the PVDF corner along the 99th-percentile angle, each at the 97th percentile of distance from background among pixels in that direction.
6.	Unmix. Each pixel p is written as p = S + a(CNT − S) + b(PVDF − S) and solved for the abundances a (CNT) and b (PVDF). Negative values are set to 0; if a + b > 1 both are scaled so that a + b = 1. The CNT fraction of the material is x = a/(a + b).
7.	Material mask. Pixels with high total abundance a + b (Otsu threshold with hysteresis, low threshold = 0.55 × Otsu). Objects of ≤12 pixels are removed. Everything else is low signal (pore or shadow).
8.	Classes. CNT-rich: x > 0.6 (objects of ≤4 pixels removed); PVDF-rich: x < 0.4; mixed: in between.
9.	Metrics (Section 5) are computed from the classes and the smoothed maps.
10.	Sensitivity. Steps 4–9 are repeated under ten settings (Section 6); every metric is reported as baseline, minimum and maximum.
5. Metrics
Table 1. Metric definitions and interpretation
Metric (Excel column)	Definition	Interpretation
C–F correlation (corr_CF)	Pearson r between smoothed C and F over material pixels	Negative: CNT and PVDF occupy separate regions. Positive: polymer mixed through the CNTs
F in CNT-rich domains (F_in_CNT)	(FCNT corner − Fbackground) / (FPVDF corner − Fbackground). Baseline = median over σ ∈ {0.7, 1, 1.5} × corner percentiles {0.5/99.5, 1/99, 3/97}; range = min–max of these 9	≤ 0: CNT domains carry no polymer signal (exposed at EDS resolution). > 0: polymer present in CNT domains
Low-signal area (low_signal_frac)	Fraction of the field outside the material mask	Open pores or shadowed depressions
PVDF in largest network (PVDF_largest_network)	Fraction of PVDF-rich pixels in the largest 8-connected PVDF-rich region	Continuity of the polymer network
PVDF within 5 px of CNT (PVDF_within_5px)	Fraction of PVDF-rich pixels whose path distance to the nearest CNT-rich pixel, travelling only through material, is ≤ 5 px	Proximity. High values can also mean polymer blended inside CNT agglomerates; read together with F_in_CNT
Median PVDF→CNT distance (PVDF_median_dist_px)	Median of the same path distances (px)	Typical distance from polymer to the nearest CNT domain
Class areas (area_CNT, area_mixed, area_PVDF)	Fraction of material pixels in each class	Relative within a map (corners are set per map); use to check that compositions are similar
Particles only: CNT_perimeter_bordering_material, free_CNT_area	Edge of CNT domains facing material rather than substrate; CNT area in objects that contain no PVDF	Attachment of CNTs to the polymer scaffold
6. Sensitivity analysis
Table 2. The ten analysis settings
Setting	σ (px)	Class thresholds (x)	Mask low threshold (× Otsu)	Corner percentiles
0 (baseline)	1.0	0.4 / 0.6	0.55	1 / 99
1	0.7	0.4 / 0.6	0.55	1 / 99
2	1.5	0.4 / 0.6	0.55	1 / 99
3	1.0	0.35 / 0.65	0.55	1 / 99
4	1.0	0.45 / 0.55	0.55	1 / 99
5	1.0	0.4 / 0.6	0.40	1 / 99
6	1.0	0.4 / 0.6	0.75	1 / 99
7	1.0	0.4 / 0.6	0.55	0.5 / 99.5
8	1.0	0.4 / 0.6	0.55	3 / 97
9	1.0	0.4 / 0.6	0.55	5 / 95
A difference between two groups is called robust when no sample of one group overlaps with any sample of the other across all ten settings (column “groups overlap across settings?” in the Group_stats sheet).
7. Statistics
•	Group mean ± SD, minimum and maximum are reported for each metric.
•	Welch’s t-test (unequal variances) compares two groups. With n = 3 per group, the exact two-sided Mann–Whitney p cannot be lower than 0.10 whatever the data, so report mean ± SD and ranges, and name the test and n when quoting p.
•	Maps from the same sample are not independent; average them per sample before comparing groups.
8. Validation on the current data
Running the toolkit on the six reports (three FS-CNT/PVDF-SDC interlayers, three cast CNT/PVDF films) reproduces the values in the 6 October contact-topology report exactly.
Table 3. Baseline metrics per sample (from Summary sheet)
Sample	Group	C–F r	F in CNT	Low signal (%)	PVDF network (%)	PVDF ≤5 px (%)	CNT-rich (%)
film 2	SDC	−0.39	−0.26	7	48	7	29
film 3	SDC	−0.36	−0.10	8	49	12	58
sdc 50	SDC	−0.42	−0.23	7	55	27	27
cst 2	Cast	0.83	0.21	28	38	43	44
cast3	Cast	0.75	0.30	23	29	43	36
cast 50	Cast	0.77	0.19	33	16	52	28
Table 4. Group statistics (from Group_stats sheet)
Metric	SDC mean ± SD	Cast mean ± SD	Welch p	Overlap across settings
C–F correlation within material	−0.39 ± 0.03	0.78 ± 0.04	5 × 10−6	no
F in CNT-rich domains	−0.20 ± 0.08	0.24 ± 0.06	0.003	no
Low-signal (pore/shadow) area	7% ± 1%	28% ± 5%	0.020	no
PVDF-rich area in largest connected network	51% ± 4%	27% ± 11%	0.054	yes
PVDF-rich area within 5 px of a CNT-rich domain	16% ± 10%	46% ± 5%	0.022	yes
CNT-rich area	38% ± 17%	36% ± 8%	0.879	yes
 
Figure 2. C (red)/F (green) maps and unmixed classes for the six samples (orange CNT-rich, blue mixed, green PVDF-rich, dark grey low signal).

 9. How to run
9.1 Option A — Anaconda (on your computer)
1.	Install Anaconda (free) if it is not installed. It already contains every package the toolkit needs.
2.	Make a folder, e.g. EDS_analysis, and copy into it: eds_topology.py, EDS_topology_walkthrough.ipynb, samples.csv, and a sub-folder reports with the AZtec .docx files.
3.	Open Anaconda Navigator → Jupyter Notebook (or JupyterLab), browse to the folder, and open EDS_topology_walkthrough.ipynb.
4.	Run the cells from top to bottom with Shift + Enter. Each step explains what it does and what you should see.
5.	Results are written to the folder EDS_topology_results.
9.2 Option B — Google Colab (no installation)
1.	Go to colab.research.google.com → Upload → choose EDS_topology_walkthrough.ipynb.
2.	Run the first code cell (“Colab only”) and upload eds_topology.py, samples.csv and the .docx reports when asked.
3.	Run the remaining cells. The last export cell downloads the results as a .zip file.
9.3 Option C — one command (batch)
From the Anaconda Prompt, in the folder that contains eds_topology.py:
python eds_topology.py samples.csv --folder reports -o EDS_topology_results
 
9.4 The sample table (samples.csv)
Column	Meaning	Example
file	AZtec report file name (inside the reports folder)	CNt-PVDF_sdc film 2_2026-10-06_18-41-35.docx
sample	Short name used in tables and figures	film 2
group	Sample type	SDC, Cast, Post-mixed
elements	Element of each map in report order, separated by “;” (X = ignore)	C;F;O
background	low for films; si for particles on a Si wafer	low
10. Outputs and how to edit them
Output	Content	Use
EDS_topology_results.xlsx	Sheets: Summary (one row per sample), Ranges (baseline/min/max), Group_stats, x_distribution, Distance_CDF, Corners, Metric_definitions, Settings	Copy into your data workbook or Origin
csv/*.csv	The same tables as CSV	Origin: File → Import → Single ASCII
raw/<sample>_C_map.csv, _F_map.csv	Full pixel matrices (384 rows × 512 columns) of the raw maps	Open in Excel; re-plot as a contour/heat map in Origin
raw/<sample>_C_smoothed.csv, _F_smoothed.csv	Smoothed maps used for unmixing	Same
raw/<sample>_CNT_fraction.csv	x = CNT/(CNT + PVDF) per pixel (empty outside material)	Same
raw/<sample>_class_labels.csv	0 low signal, 1 PVDF-rich, 2 mixed, 3 CNT-rich	Same
data/<sample>_scatter.csv	5000 random pixels: C, F, class	Re-plot the C–F triangle in Origin
maps/.png, .tif	Class map and C/F overlay (PNG); class labels (8-bit TIFF) and CNT fraction (32-bit TIFF)	Slides; ImageJ/Fiji
figures/.svg, .pdf, *.png	Class maps, C–F triangles, metric dot-plot, distance curves	PowerPoint: Insert → Pictures → .svg, then right-click → Convert to Shape (text stays editable). Illustrator/Inkscape open SVG/PDF directly
11. Limitations
•	EDS spatial resolution is ≈1 µm (interaction volume), so a PVDF skin of tens of nanometres on a CNT is not resolved, and “mixed” pixels include CNT and PVDF overlapping in depth.
•	Corners are set per map, so class areas are relative within a map; the correlation and the F ratio are less sensitive to this.
•	Rough surfaces add a positive C–F correlation (shadows lower both signals). The F-in-CNT ratio is less affected and should be read alongside r.
•	Report images (not raw counts) are analysed; their display scaling is unknown but monotonic, which the correlation, ratio and class metrics tolerate.
•	Pixel distances depend on magnification; compare distance metrics only between maps with the same field width.
12. Methods text for the Supporting Information
SEM–EDS element maps (Oxford Instruments AZtec) were analysed with a custom Python script (NumPy, SciPy, scikit-image). C Kα and F Kα maps were resampled to the native 512 × 384 map grid and smoothed with a Gaussian filter (σ = 1 pixel). Because PVDF contains carbon, the C map alone does not locate the CNTs; in the C–F intensity plane the pixels fall inside a triangle whose vertices correspond to the background (low-signal pores and shadows), the CNTs (C without F), and PVDF (C with F). The vertices were located automatically, each pixel was expressed as a linear combination of the three end-members, and its CNT fraction, x = CNT/(CNT + PVDF), was used to classify material pixels as CNT-rich (x > 0.6), PVDF-rich (x < 0.4), or mixed. The C–F correlation within the material, the F signal at the CNT vertex relative to the PVDF vertex, the low-signal area fraction, the fraction of PVDF-rich area in the largest connected component, and path distances measured within the material were computed. All metrics were recomputed over ten analysis settings that varied the smoothing (σ = 0.7–1.5 pixels), the class thresholds (0.35/0.65 to 0.45/0.55), the material-mask threshold, and the vertex selection; reported ranges span these settings. Three independent samples of each interlayer type were mapped under identical acquisition conditions ([__] kV, [__] µm field width).
