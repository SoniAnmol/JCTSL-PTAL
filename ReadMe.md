# PTAL Analysis of Jaipur

This repository is dedicated to the project of mapping the Public Transport Accessibility Levels (PTAL) of Jaipur
municipality. The aim of this project is to assess and visualize the accessibility levels to public transport services
within Jaipur, leveraging data available in the public domain.

## Table of Contents

- [Introduction](#introduction)
- [Project Overview](#project-overview)
- [Data Sources](#data-sources)
- [Project Structure](#project-structure)
- [Usage](#usage)

## Introduction

This project focuses on measuring macro-level accessibility to public transport across Jaipur Municipality.
Macro-level accessibility refers to the potential opportunities for interaction between
transport service points (i.e. JCTSL bus services) at spatially segregated points (Grids of 1 square Km) of interest.
This study specifically measures access to public
transport service points from various locations within the city, without taking into account monetary or physical
constraints that may hinder access.

## Project Overview

Measurement of Public Transport Accessibility:

1. Computation of the Accessibility Index.
    1. Identification of the nearest JCTSL bus stop from each 1 sq km grid's centroid.
    1. Calculation of walking distances to the nearest bus stop.
    1. Estimation of average waiting times based on bus frequency.

1. Visualization of PTAL:
    1. Creation of PTAL maps for Jaipur, highlighting areas with different accessibility levels.
    1. Overlaying PTAL maps with population density data for a comprehensive view.

1. Policy Implications:
    1. Discussion of policy recommendations to enhance public transport accessibility in Jaipur.

## Data Sources

* **Population Data**: Sourced from the Smart Cities Mission data portal, based on the 2011 census. Note: Population counts
  are subdivided into 1 sq km grids.
* **Public Transport Data**: Includes data on the Jaipur City Transport Services Limited (JCTSL) bus routes,
  frequencies, and locations of bus stops.
* **Geospatial Data**: Road network maps used for calculating walking distances downloaded from openstreetmaps.

## Project Structure

The project directory is organized as follows:

* `` data/``: Contains population data, bus routes data, and geospatial data.
* ``figures/``: Contains all the plots and figures, visualisations and PTAL maps created during this project.
* ``notebooks/``: Contains jupyter notebooks to clean, prepare and analyse data for computation of PTAL.
* ``scripts/``: Includes scripts for data analysis, PTAL computation, and visualisation.

## Usage

* Clone the repository to your local machine.
* Install the required dependencies (listed in requirements.txt).
* Use the provided notebooks in the given order to perform PTAL analysis on your dataset.



