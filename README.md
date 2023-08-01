# pnas-hurricane-regimes
Data and code repository for Cannon et al. hurricane regimes

These code, output data, and figures are associated with Cannon et al. Hurricane wind regimes for forests of North America (*PNAS*, in review).

See article for complete study description and results.

Code in these scripts requires use of the *hurrecon* library. See https://github.com/jbcannon/hurrecon for details on installation and use.

# output/

* *regimes.zip:* (contains .shp, .dbf, .prj, and .shx) shapefile indicating four-hurricane regimes boundaries delinated in Cannon et al. article including boundaries for Continental, Inland, Coastal, and Fringe regimes. See Figure 2 and Figure 3 in corresponding article 

* *hurr-risk-summary.csv:* table summarizing hurricane regime for 4 areas indentified in regimes.zip. Table includes information on 1 year occurence probability, 10-year occurence probability, 30-year occurence probability, and return period for each storm intensity level. See Table 2 in corresponding article.

* *SS_freq_tally-AL.tif:* six layer .tiff representing a tally of the number of times wind speeds were within the range for tropical storm force winds, Saffir-Simpson I, II, III, IV, and V, respectively (for each layer). Includes the 171 year period where data is available for the Atlantic Basin. 

* *SS_freq_tally-EP.tif:* six layer .tiff representing a tally of the number of times wind speeds were within the range for tropical storm force winds, Saffir-Simpson I, II, III, IV, and V, respectively (for each layer). Includes the 73 year period where data is available for the Eastern Pacific Basin.

* *SS_risk_01yr.tif:* annual occurence probability (count/years elapsed) of wind speed ranges corresponding to tropical storm force winds, Saffir-Simpson I, II, III, IV, and V, respectively (for each layer). Note that tallys for AL and EP regions counted separately (with corresponding time ranges) then combined in this layer.

* *SS_risk_10yr.tif:* 10-year occurence probability of wind speed ranges corresponding to tropical storm force winds, Saffir-Simpson I, II, III, IV, and V, respectively (for each layer). Calculated as P10 = 1 – (1 – P1)^10, where P1 is annual occurence probability (SS_risk_01yr).	

* *SS_risk_30yr.tif:* 30-year occurence probability of wind speed ranges corresponding to tropical storm force winds, Saffir-Simpson I, II, III, IV, and V, respectively (for each layer). Calculated as P30 = 1 – (1 – P1)^30, where P1 is annual occurence probability (SS_risk_01yr).

* *validation_output.csv:* Maximum sustained wind speeds as estimated from HRRR (hrrr), and hurrecon (vs), for hurricane tracks.

# figs/

* *hurrecon-wind-profile-demo_R2.png:* See Figure 1. (A) Example output of HURRECON model showing wind profile extending from hurricane eyewall using a lognormal curve to fit wind speed observations. (B) Model of sustained wind speeds (m s-1) within Hurricane Michael at landfall on the Florida, US coast at 1730 UTC 10 October 2018. 

* *cluster_map_full.png:* See Figure 2. Hurricane disturbance regimes identified for North and Central America in order of increasing hurricane activity include (A) Continental, (B) Inland, (C) Coastal, and (D) Fringe. See insets in Figure S5.

* *cluster_map_panels.png:* See Figure 3. Hurricane disturbance regimes identified for North and Central America in order of increasing hurricane activity including Continental, Inland, Coastal, and Fringe Regimes, shown for the (A) U.S. Gulf Coast, (B) U.S. Atlantic coast, and (C) Caribbean region.

* *cluster_characteristics.png:* See Figure 4. (A) 30-year occurrence probability (P30) of  hurricane winds for four identified wind disturbance regimes including probability for five levels of wind intensity corresponding to Saffir–Simpson category winds. Tropical storm (TS)-force = 17–32 m s-1. Category 1 = 33–42 m s-1, Category 2 = 42–49 m s-1, Category 3 = 49–58 m s-1, Category 4 = 58–70 m s-1, and Category 5 ≥ 70 m s-1. (B) For comparison of rare events, this panel contains the same data in A with a logarithmic y-axis.

* *cluster_analysis-scree-plot.png:* See Supplemental Figure S6. Scree plot used to determine number of clusters to use for k-means clustering algorithm showing reduction in total within-cluster sum of squares for values of k from 1 to 20.

* *data-gap-fill.png:* See Supplemental Figure S2. Results of linear models predicting rv—the radius r at which wind speed v occurs (in nautical miles, NM) based on maximum wind speed, latitude, and pressure (Eq. 2)— for post-2004 storm observations, where rv was used to predict rv for pre-2004 storm observations. Model parameters and R2 values are presented in Table S1.

* *radius_assymetry.png:* See Supplmental Figure S3. Mean relative distance from eye of hurricane to location of threshold wind speed (rv) for each of three threshold wind speeds (v = 17.5 m s-1, 25.7 m s-1, and 32.9 m s-1, corresponding to 34, 50, and 64 knots). Error bars represent ± 1 s.e. of the mean and panels are based on 12472, 4158, 4158, and 4156 observations, respectively.

* *30-yr-risk-fig_6panel.png:* See Supplemental Figure S4. 30-year probability (P30) of sustained wind speeds across the study area for tropical storm and Saffir–Simpson category windspeed ranges (A) Tropical storm (TS)-force = 17–32 m s‑1; (B) Category 1 = 33–42 m s‑1; (C) Category 2 = 43–49 m s‑1; (D) Category 3 = 50–57 m s‑1; (E) Category 4 = 58–69 m s‑1; and (F) Category 5 ≥ 70 m s‑1.For panels C-E, values of P30 < 0.0001 not shown to enhance visibility

* *validation-fig.png:* See Supplmental Figures S1. Comparison of maximum sustained wind speed estimates (vs) between the High-resolution Rapid Refresh (HRRR) model, and predictions from the HURRECON model. Blue line represents model fit based on linear regression (Table S1), compared to a 1:1 line (dashed). Zero values along x-axis are minor storms modeled by HRRR, but not captured by the HURRECON model.

# code/

* *define_roi.R:* Script to download, edit and define the region of interest used in the analysis

* *exploratory analysis.R:* Script to download and examine hurricane tracks and produce summary data, and caluclate information and models for gap filling

* *hurrecon_runs_parallel.R* Script to process individual NOAA tracks using hurrecon. Code is customizable to run in parallel for speed

* *hurrecon_cumulative_intensity.R:* Script to tally wind speed data generated from hurrecon runs and estimate P1, P10, P30

* *cluster_analysis:* Script to take P30 estimates, normalize, and use cluster analysis to generate hurricane regime clusters.

* *making_graphics.R:* Script to assemble maps, figures, etc.

* *validation_extract_tracks.R:* Data to select and extract tracks to test against HRRR runs.

* *HRR_ws_finder.R:* Script to convert data obtained from HRRR to raster and compare with outputs from *hurrecon*.

* *hrrr.py:* Python script that takes hurricane track shapefiles, downloads them from the HRRR database, and outputs wind speed and direction vectors to be processed and converted to raster using *validation_extract_tracks.R*




