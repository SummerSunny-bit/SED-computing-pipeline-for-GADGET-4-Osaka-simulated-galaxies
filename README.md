# SED computing pipeline for GADGET-4 Osaka simulated galaxies
A modular postprocessing pipeline for generating synthetic spectral energy distributions (SEDs) and broadband photometry from GADGET-4 Osaka cosmological simulations using FSPS stellar population synthesis grids.

# Overview
The pipeline converts raw GADGET-4 Osaka cosmilogical simulation outputs into observable galaxy roperties by combining stellar particle inforamtion with precomputed FSPS stellar population sythesis grids.
The pipeline supports:
- extraction of stellar populations from SUBFIND subhalos,
- stellar age and formation redshift computation,
- FSPS-based SED interpolation,
- integrated galaxy SED synthesis,
- (SDSS) broadband photometry

The primary pipeline implementaion is provided in the notebook \texttt{SEDComputing_Main_Pipeline.v1.ipynb}.
The repository further includes supporting documentation and manuels, including \textttt{GADGET-4 Osaka Manuel} and the author's own documentation, \texttt{SED Computing Pipeline for Simulated Galaxies - Technical Documentation}, 
which describes the pipeline architecture, implementation details, scientific assumptions, and usage of the included code. 

# Pipeline Structure
The pipeline consists of four major stages:

1. Snapshot Processing
   Extraction of galaxy and stellar particle data from GADGET-4 outputs.
2. FSPS Grid Generation
   Construction of reusable stellar population synthesis grids.
3. SED Computation
   Interpolation and summation of stellar spectra.
4. Photometric Synthesis
   Computation of broadband photometry.

# Supported Input Data
The current interpolation assumes:
- GADGET-4 HDF5 snapshot output,
- SUBFIND group catalogs,
- stellar particles stored as PartType4,

# Main Dependencies
Required packages include:
- h5py
- astropy
- python-fsps

Recommended Python version: Python 3.10+

# Academic Context
This project was carried out within the academic exchange framework between the University of Bremen and The University of Osaka.
