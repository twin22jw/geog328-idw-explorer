# GEOG 328 Spatial Interpolation Tools

Two interactive teaching tools for GEOG 328:

- `index.html`: inverse distance weighting (IDW), including power, search neighborhood, output cell size, sampling strategy, barriers, and neighborhood shape.
- `spline.html`: compare 1D B-spline polynomial degrees 1–3, then compare a fixed global exact thin-plate spline surface with IDW using the same 12 observations. Drag the A–B line across the maps to compare profiles extracted from both 2D models, or click any map cell to inspect its interpolated values. The page explains the thin-plate spline formula and distinguishes the fixed range of all 12 sample values from values along A–B. Adjust IDW power to examine local influence.

Both pages include group discussion questions. The observations are fictional. The tools use local HTML, CSS, and JavaScript without external runtime dependencies.

Public pages: [IDW](https://twin22jw.github.io/geog328-idw-explorer/) and [Spline](https://twin22jw.github.io/geog328-idw-explorer/spline.html).

