# GEOG 328 Spatial Interpolation Tools

Three interactive teaching tools for GEOG 328:

- `index.html`: inverse distance weighting (IDW), including power, search neighborhood, output cell size, sampling strategy, barriers, and neighborhood shape.
- `spline.html`: explore exact and smoothed 1D B-spline fits with polynomial degrees 1–3, then compare a global exact thin-plate spline with IDW (fixed α = 2). Use 12 or 24 shared observations; inspect values by numbered sample ID; drag A–B; and click a map cell or profile to read predictions. Students enter their own sample minimum and maximum as graph reference lines. The page has three group discussion experiments: 1D polynomial degree and smoothing, spline versus IDW surfaces, and a 12 versus 24 observation comparison.

- `kriging.html`: build two semivariograms from independent 12-observation examples with distinct short and long fitted ranges, a fixed nugget of 10, progressive pair and lag-bin reveals, and optional range/nugget/sill annotations. A separate ordinary Kriging example compares prediction, model-based standard deviation, and IDW weights while students adjust range, nugget, sill, sampling design, and observation count.

All three pages include group discussion questions. The observations are fictional. The tools use local HTML, CSS, and JavaScript without external runtime dependencies.

Public pages: [IDW](https://twin22jw.github.io/geog328-idw-explorer/), [Spline](https://twin22jw.github.io/geog328-idw-explorer/spline.html), and [Kriging](https://twin22jw.github.io/geog328-idw-explorer/kriging.html).


