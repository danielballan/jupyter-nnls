# A Jupyter Widget Front-end to Fortran NNLS Code

## Installation

Obtain Python 3.6+. Install Jupyter Lab and ipywidgets. Start JupyterLab:

```
jupyter lab
```

## Proof of Concept

This runs a Fortran program to add two numbers from Jupyter Widgets.

1. Compile the fortran program.

    ```
    gfortran add.f
    ```

   We accept the default output file, ``a.out``.

2. Execute the cells. Widgets will be displayed. Now we have a poor man's "web
   GUI" interface to a Fortran program.

3. Enter numbers; click Compute.

## NNLS

1. Compile the Fortran program.

    ```
    gfortran salt_study.f -o salt_study
    ```

2. Execute the cells. Widgets will be displayed. Now we have a poor man's "web
   GUI" interface to a Fortran program.

3. Enter numbers into the b and y matricies. Click compute.
