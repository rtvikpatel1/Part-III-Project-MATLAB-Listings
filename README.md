# Fluid Film Height & Hydraulic Jump Analysis

This repository contains MATLAB scripts for analyzing fluid film height and hydraulic jump behavior using experimental DFI images and theoretical/simulated models.

## Contents

### Experimental & Theoretical Comparison
- **`big_plot.m`**  
  Extracts and compares experimental jump radii with theoretical predictions (Bhagat 2020) from DFI images.

- **`comparin_4_theoretical_models.m`**  
  Solves and plots four theoretical models for fluid height from an impinging jet in polar coordinates.

- **`exp_vs_theoretical_film_thickness.m`**  
  Compares film thickness profiles at various angles using experimental data and ODE-based theoretical predictions.

- **`jump_radius_as_a_function_of_angle.m`**  
  Compares experimentally measured jump radii across angles with three theoretical models in a polar plot.

- **`theroetical_integrated_vs_experimental.m`**  
  Compares experimental jump radii with those predicted by a fully integrated ODE model over 0°–180°.

### Image and Data Utilities
- **`dfi2mat.m`**  
  Loads `.dfi` files into MATLAB structures, converting Digiflow image format to standard arrays.

- **`dfi_to_png.m`**  
  Converts `.dfi` image files to `.png` for easy viewing and export.

### Geometry and Simulation Analysis
- **`heights_super_ellipse.m`**  
  Generates normalized film height profiles from a superellipse-based theoretical model for various exponents.

- **`radial_square_experimental.m`**  
  Plots angular flow quantity \( q^{3/4}(\theta) \) across different geometric exponents with experimental comparison.

- **`laminar_to_turbulent_h.m`**  
  Compares laminar, turbulent, and hybrid (transition) models for hydraulic jump radii over a full angular sweep.

- **`openfoam_jump_location.m`**  
  Analyzes jump radius from OpenFOAM and experimental data, comparing with a mixed theoretical model (Bhagat 2020).

- **`height_openfoam.m`**  
  Visual comparison of film thickness profiles from OpenFOAM simulations and experiments, validated interactively.

---

## Requirements
- MATLAB R2018a or later (R2014a minimum for some functions)
- Image Processing Toolbox
- `.dfi` files from Digiflow software
- OpenFOAM simulation CSV output (if applicable)

## License
This repository is for academic and research use only. Please contact the original author for other uses.

