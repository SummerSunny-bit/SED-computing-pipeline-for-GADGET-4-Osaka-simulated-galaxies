#===== NOTE: Imports ====#
# Standard-Bibliotheken
import os
import shutil
import time
from pathlib import Path

# Datenverarbeitung & Plotting
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from tqdm import tqdm

# Wissenschaft & Astronomie
import astropy.units as u
import fsps
import h5py
from astropy.cosmology import FlatLambdaCDM, LambdaCDM

#===== Grid Creation =====#
# This snippet will create a FSPS SSP-Grid (ages_gyr x logzsol) and stores it as a .hf/ .HDF5.
# Default parameters include MIST Isochrones, Chabrier IMF, zcontinous = 1
# NOTE: Please refer to the official documentation: https://dfm.io/python-fsps/current/stellarpop_api/#api-reference

#=== Configuration ===#
OUTFILE = "/home/rathjen/Synthesizer_Osaka/BigData/FSPS_Grids/fsps_grid_mist_chabrier.h5"
    # Please rename the Girdname for multiple grids, otherwise they will be overwritten.
Z_SUN = 0.014           # MISt ~0.014 / Padova ~0.019
ISOCHRONES = "MIST"     # "MIST" or "Padova", or set up others in your environment
IMF_TYPE = 1            # 0: Salpeter,1: Chabrier, 2: Kroupa, 3: van Dokkum, 4: Dave
ZCONTINUOUS = 1         # 1: SSPs are interpolated to the value of logzsol
SFH = 0                 # 0: SSP, 1: 6 paramter SFH, etc. 

# ! Grid Size
N_AGE = 10              # Defines the # of grid points for Age between AGE_MIN_GYR and AGE_MAX_GYR
N_LOGZ = 5              # Defines # of grid points for Metallicity, respectively
# Set Grid boarder
AGE_MIN_GYR = 1e-3
AGE_MAX_GYR = 13.8
LOGZ_MIN = -2
LOGZ_MAX = 0.5

L_SUN_CGS = 3.828e33    # For later conversion [erg/s]

#=== Store Grid ===#
DTYPE_SPEC = "f4"
COMPRESSION = "gzip"
CHUNKS = None           # Will set to (1,1,n_lambda) after the wavelenght is known
# Ensure, output dictionary exists 
os.makedirs(os.path.dirname(OUTFILE), exist_ok=True)

#=== FSPS Object ===#
sp = fsps.StellarPopulation(imf_type=IMF_TYPE,zcontinuous=ZCONTINUOUS,sfh=SFH)
#sp.params["isochrones"] = "MIST"
    
# Safety: Disable unwanted physics
sp.params["add_neb_emission"] = False
sp.params["add_neb_continuum"] = False
sp.params["dust2"] = 0.0
sp.params["dust1"] = 0.0
sp.params["igm_factor"] = 0.0
sp.params["fagn"] = 0.0

#=== Grid ===#
ages_gyr = np.logspace(np.log10(AGE_MIN_GYR), np.log10(AGE_MAX_GYR), N_AGE)
logzsol = np.linspace(LOGZ_MIN, LOGZ_MAX, N_LOGZ)

#=== Grid Creation ===#
with h5py.File(OUTFILE, "w") as f:

    # Speichere Age- und Z-Achsen
    f.create_dataset("ages_gyr", data=ages_gyr)
    f.create_dataset("logzsol", data=logzsol)

    pbar = tqdm(total=N_LOGZ, desc="Building FSPS grid", unit="Z")
    start_time = time.time()

    lam_ds = None
    spec_ds = None

    for j, lz in enumerate(logzsol):

        sp.params["logzsol"] = float(lz)

        # --- Beim ersten Durchlauf: Wellenlängenachse bestimmen ---
        if lam_ds is None:
            wave, _ = sp.get_spectrum(tage=float(ages_gyr[0]), peraa=True)
            lam_ds = f.create_dataset("wavelength", data=wave)
            n_lambda = wave.size

            # Spektren-Dataset anlegen
            spec_ds = f.create_dataset(
                "spectra",
                shape=(N_AGE, N_LOGZ, n_lambda),
                dtype=DTYPE_SPEC,
                compression=COMPRESSION,
                chunks=(1, 1, n_lambda)
            )

        # --- Block für alle Alterswerte ---
        spec_block = np.zeros((N_AGE, n_lambda), dtype=float)

        for i, age in enumerate(ages_gyr):
            wave, spec = sp.get_spectrum(tage=float(age), peraa=True)
            spec_block[i, :] = spec * L_SUN_CGS  # erg/s/Å per Msun

        # --- Block in HDF5 schreiben ---
        spec_ds[:, j, :] = spec_block.astype(DTYPE_SPEC)

        pbar.update(1)

    pbar.close()
    elapsed = time.time() - start_time

   #=== Meta Data ===#
    f.attrs["Z_sun"] = Z_SUN
    f.attrs["imf_type"] = IMF_TYPE
    f.attrs["zcontinuous"] = ZCONTINUOUS
    f.attrs["sfh"] = SFH
    f.attrs["fsps_version"] = fsps.__version__
    f.attrs["units"] = "erg/s/Angstrom per Msun"
    f.attrs["created_on"] = time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())
    f.attrs["n_age"] = N_AGE
    f.attrs["n_logz"] = N_LOGZ
    f.attrs["build_time_seconds"] = elapsed

print(f"Grid written to {OUTFILE} in {elapsed:.1f} s")