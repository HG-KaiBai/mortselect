# Mortality Forecast Explorer

The public explorer is served from this directory by GitHub Pages:
https://hg-kaibai.github.io/mortselect/

It uses the original saved forecast matrices from the five RData files in the
author's `Shiny` directory. It performs no model fitting and displays no observed
mortality data. Values are stored as unscaled rates and displayed as rate × 100.
Raw population labels are preserved in downloads. Display labels only remove
the leading `--` and shorten three total-population names.

## Update the results

From the repository root, with `jsonlite` installed:

```r
system("Rscript tools/export_forecasts.R ../Shiny docs/data")
```

The exporter checks all matrix values after JSON serialization, preserves age
and year labels, and records source-file MD5 hashes in `data/manifest.json`.
Only each result's forecast matrix, population and sex are exported. The
original workspaces and fitting objects are not published.

## Preview

```sh
node tools/preview_site.mjs docs 4173
```

Open http://127.0.0.1:4173. No build, external JavaScript, fonts, or chart service
is needed. Publishing `master:/docs` provides a static site with no R process
or session timeout. Availability remains subject to GitHub Pages service limits.

## Native Shiny version

```r
install.packages(c("shiny", "ggplot2", "jsonlite"))
shiny::runApp("shiny")
```

The Shiny app reads exactly the same exported data. The original `SingleASS`
and `MultiASS` data identifiers are retained, while the interface uses the
paper's ASB terminology. Paper references and method mappings follow the
package README.
