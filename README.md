# LYRA: Python Likelihood Ratio for AGN-DB

LYRA matches sources between astronomical catalogs using positional information and, optionally, magnitude priors from a likelihood-ratio analysis. It can work with overlapping or non-overlapping catalogs and with different positional uncertainties. It extends LYR (Nanni et al. 2020; Peca et al. 2021).

If you use LYRA, please cite Peca et al. (in prep.).

## Requirements

Python 3 and the following packages are required:

```bash
python -m pip install numpy scipy scikit-learn pandas matplotlib ipython
```

For FITS input, also install `astropy`. The repository does not currently provide an installable package; run `LYRA_main.py` from this directory. If your environment uses `anapython`, you can substitute it for `python` in the commands below.

## Input catalogs

The two positional catalog arguments can be CSV files or FITS tables. By default, each catalog needs columns named `ID`, `RA`, and `DEC`; use the `--ID_*`, `--RA_*`, and `--DEC_*` options when names differ. Coordinates are in degrees. The two positional-error arguments can each be a constant value or the name of a column in the corresponding catalog; errors and the angular-resolution argument are in arcseconds. Magnitude columns, when used, must be in the primary catalog.

Run `python LYRA_main.py --help` for all available options.

## Examples

Run these commands from the repository directory. Make the output directory first: the program opens the information file before it creates missing directories. Output paths should end in `/` because the program joins paths by string concatenation.

### Positional matching with constant errors

```bash
mkdir -p results
python LYRA_main.py 7.0 0.1 1.0 \
  data/primary.csv data/secondary.csv \
  --path_output results/ \
  --filename_LR positional_matches \
  --filename_outinfo positional_info.txt
```

Here `7.0` is the angular-resolution/search-radius setting, and `0.1` and `1.0` are the primary and secondary positional errors. With no `--mag_col`, LYRA performs positional matching.

### SDSS Fliri+16 with magnitude priors

This adapts a historical Stripe 82 run. Replace the sample paths with your catalog locations:

```bash
mkdir -p results
python LYRA_main.py 7.0 0.1 POS_ERR \
  data/iacstr82_Fliri+16.csv \
  data/CSC_XMM_15_comberr_CORR4_slim13.csv \
  --custom_RC_threshold 3 \
  --ID_prim sdss_iac_id --ID_sec XID_tmp \
  --RA_prim ALPHA_J2000 --DEC_prim DELTA_J2000 \
  --mag_col MAG_AUTO_u,MAG_AUTO_g,MAG_AUTO_r,MAG_AUTO_i,MAG_AUTO_z \
  --nomag 99 \
  --path_LR results/ --path_output results/ --path_images results/ \
  --filename_LR out_table_iac \
  --filename_outinfo outinfo_iac.txt \
  --save_images True
```

`POS_ERR` names a column in the secondary catalog. The magnitude columns and `99` missing-value marker refer to the primary catalog. `--mag_col` is the current option name (not `--mag`). Supply `--filename_LR` **without** `.csv`; LYRA appends the extension.

### FITS catalogs and custom coordinate columns

```bash
mkdir -p results
python -m pip install astropy
python LYRA_main.py 5.0 ERR_PRIMARY ERR_SECONDARY \
  data/primary.fits data/secondary.fits \
  --ID_prim SOURCE_ID --RA_prim RA_ICRS --DEC_prim DEC_ICRS \
  --ID_sec SOURCE_ID --RA_sec RA_ICRS --DEC_sec DEC_ICRS \
  --path_output results/ \
  --filename_LR fits_matches \
  --filename_outinfo fits_info.txt
```

The error names in this example are columns in their respective FITS tables. Change them to match your data.

## Outputs

With the example settings, LYRA writes the match table to `results/<filename_LR>.csv` and a run-information text file to `results/<filename_outinfo>`. When `--save_images True` is set, plots are written under `--path_images`. The default names are `tmp_table.csv` and `tmp_info`, but setting names and paths explicitly makes runs easier to identify.

## License and contact

This repository includes the GNU General Public License v3 in [LICENSE](LICENSE). For questions, contact alessandro.peca@miami.edu.
