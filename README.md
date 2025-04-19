# Bending Analysis of a Rectangular Plate under Uniform Load

This project performs a bending analysis of a simply supported rectangular plate subjected to a uniform load using classical plate theory. The code is written in Python and utilizes NumPy and Matplotlib for numerical computations and visualization.

## Features

- Calculates displacement, bending moments, and Von Mises stress distributions for a rectangular plate.
- Supports varying the number of Fourier terms (N) for convergence studies.
- Visualizes:
    - Displacement contours
    - Von Mises stress contours
    - Through-thickness stress profiles at the plate center
    - Strain energy as a function of N

## Problem Description

- Plate dimensions: `a = b = 0.5 m`
- Plate thickness: `h = 0.004 m`
- Young's modulus: `E = 200 GPa`
- Poisson's ratio: `nu = 0.3`
- Uniform load: `q0 = 1000 N/m²`
- Yield stress: `450 MPa`

## Code Structure

- **compute_coefficients(N):**  
    Computes the Fourier coefficients for the plate solution using N terms.

- **compute_fields(x, y, coeffs):**  
    Calculates displacement (`w`), bending moments (`Mx`, `My`), and twisting moment (`Mxy`) at a given (x, y) location.

- **Main Execution:**  
    - Sets up mesh grids for the plate surface and through-thickness direction.
    - Loops over different values of N (2, 4, 6) to study convergence.
    - For each N:
        - Computes displacement and Von Mises stress fields.
        - Plots contour maps for displacement and stress.
        - Calculates and stores strain energy.
    - Plots through-thickness stress profiles at the plate center.
    - Plots strain energy versus N.

## Usage

1. Ensure you have Python 3.x installed.
2. Install required packages:
     ```
     pip install numpy matplotlib
     ```
3. Run the script:
     ```
     python <script_name>.py
     ```
4. The script will display several plots for displacement, stress, and strain energy.

## Notes

- The code assumes simply supported boundary conditions on all edges.
- The accuracy of results improves with increasing N, at the cost of computation time.
- The script is modular and can be adapted for different plate sizes, thicknesses, or loading conditions.

## Output

- Displacement and Von Mises stress contour plots for N = 2, 4, 6.
- Through-thickness stress variation at the plate center.
- Strain energy convergence plot.

---