# Causes of range collapse and extinction in the wild of European bison

Code contributors: Julia A Pilowsky, Stuart C Brown, Damien A Fordham.

This repository contains the `R` code and some example data to run the validated ensemble of `paleopop` simulations of the European bison.

In order to run the example simulations, `paleopop v 2.0.0` or later must be installed. It is available on CRAN.

The code has a number of package dependencies listed below:

- "poems"
- "paleopop"
- "raster"
- "purrr"
- "stringr"
- "dplyr"
- "furrr"
- "data.table"
- "sf"
- "gdistance"

The scripts are designed to run in parallel across *n* sessions - please set *n* accordingly.

There are three scripts contained in the repository. You may run `bison-baseline.R`, `bison-no-hunting.R`, or both, as desired. The script `bison-summary-metrics.R` must be run last.

You may simulate the "baseline" scenario (environmental change, hunting, and land use change as drivers) using `bison-baseline.R`, and the "no hunting" scenario (only environmental change and land use change as drivers) using `bison-no-hunting.R`. The `bison_summary_metrics.R` script calculates summary metrics for the simulation outputs. These summary metrics can be compared against observed patterns for model selection by Approximate Bayesian Computation, as described in the main manuscript.

All outputs are saved in the `results/` directory.

Additional niche samples are available on reasonable request.

Directory structure:

```
european-bison/
  - Data/ # data necessary to build the models
  - k_cuts/ # niche samples required for the models
  - results/ # output directory for simulation results
```

This code is released under the following licence:

GNU GENERAL PUBLIC LICENSE
   Version 3, 29 June 2007

Copyright (C) 2007 Free Software Foundation, Inc. <https://fsf.org/>
Everyone is permitted to copy and distribute verbatim copies
of this license document, but changing it is not allowed.
