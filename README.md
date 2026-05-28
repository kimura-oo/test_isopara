# FEM Isoparametric Solvers

This repository contains a collection of Finite Element Method (FEM) solvers written in Fortran 90. It focuses on solving various Partial Differential Equations (PDEs) using isoparametric elements of different orders (linear P1, quadratic P2, and cubic P3) and implements the Streamline Upwind Petrov-Galerkin (SUPG) method for stabilized solutions in advection-dominated problems.

## Directory Structure

The repository is organized into several main directories based on the target equations and functionalities:

* **`adv_supg/`**
    Solvers for the advection equation using the SUPG stabilization method. It is subdivided by element order: `advec_p1/`, `advec_p2/`, and `advec_p3/`.
* **`convecdiff_supg/`**
    Solvers for the convection-diffusion equation employing the SUPG method.
* **`laplace/`**
    Solvers for the Laplace equation, including implementations for quadratic FEM and various test meshes.
* **`poisson_isopara/`**
    Solvers for the Poisson equation utilizing isoparametric elements. It includes support for linear (`poisson_p1/`), quadratic (`poisson_p2/`), and cubic (`poisson_p3/`) elements.
* **`meshconv_p1_to_p2/` & `meshconv_p1_to_p3/`**
    Mesh conversion utilities. These tools convert standard linear meshes (P1) into higher-order quadratic (P2) or cubic (P3) meshes required for higher-order solvers.
