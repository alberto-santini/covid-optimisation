# Optimisation tools for the COVID-19 pandemic

## About

This repository was created by [Alberto Santini](https://santini.in/), professor of Operational Research at [Universitat Pompeu Fabra](https://upf.edu/).
It contains optimisation models and software to optimise decisions related to the COVID-19 pandemics.

### Main links

* [Interactive dashboard](https://santini.in/covid/) with results (optimisation of swab testing capacity, Italy case study).
* [Technical report](https://santini.in/files/tech-rep-swabs.pdf) on the optimisation of swab testing cpacity.

![dashboard-screenshot.png](https://github.com/alberto-santini/covid-optimisation/raw/master/images/dashboard-screenshot.png)

### Contents of the repository

* Folder `model-implementation` contains the implemented models. The software is written in C++ and uses CMake for configuration. The solver [Gurobi](https://www.gurobi.com/) is necessary to run the software. All software is released under the GPLv3 license (see file `LICENSE`).
* Folder `swab-tests-data` contains data for the problem of increasing swab tests capacity. Find a description of this problem and a mathematical model to solve it in this [technical report](https://santini.in/files/tech-rep-swabs.pdf), which is updated constantly.
    * Subfolder `italy-case-study` contains data (and the programme to generate it) relative to swabs tested in Italy during the period 01-13 April, 2020.
    * Subfolder `synthetic-data` contains simulated data used to test how the model behaves under different scenarios.

## Data Generation

This sections explains in detail how we generate data for the case studies.

### Swab tests optimisation - Italy case study

We start from the publicly available data sources:
    * A list of laboratories which are certified to test swabs for COVID-19: `italy_labs_coords.csv`.
    * A list of factories in Italy, which are certified to produce test kits: `italy_factories_coords.csv`.
    * Data about the number of swabs tested each day, region-by-region, released by the Italian Civil Defense department: `dpc-covid19-ita-regioni.csv`.

We first estimate how many swabs regions would have liked to test each day, if they had enough reagents.
This number is calculated from the actual number they tested, times some region-dependent multiplier which we estimate based on news and press releases (see `instance-generator.py` for the precise numbers).
We assume that, before the planner's intervention, each region would allocate an equal number of swabs to each of its laboratory.
Because we start our analysis mid-crisis (April 1st) we assume that there is no safety stock of reagents neither at laboratories, nor at factories.
We also assume that laboratories structural capacity (i.e., the capacity determined by staffing and test machines, but not by reagent availability) has a +25% slack compared to the number of swabs they actually tested during the period we analysed.
This hypothesis means that we think (based on news and statements released by regional health authorities) that reagent availability was the main bottleneck that prevented increasing test capacity.
The number of swabs and reagents that can be transfered to each region day by day was estimated based on geographic and economic factors peculiar to each region.
For example, Sardinian labs cannot exchange swabs with labs in other regions, because it would be hard to ship them out of the island.
Some regions, such as Lombardia, Piemonte, Emilia Romagna or Campania have good infrastructure and are well-connected, so they are able to move more material.

### Swab test optimisation - synthetic data

TODO.