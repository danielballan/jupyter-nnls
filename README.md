## Proof of Concept

This runs a Fortran program to add two numbers from Jupyter Widgets.

1. Compile the fortran program.

    ```
    gfortran add.f
    ```

   We accept the default output file, ``a.out``.

2. Start a Jupyter notebook server.

    ```
    jupyter lab
    ```

3. Execute the cells. Widgets will be displayed. Now we have a poor man's "web
   GUI" interface to a Fortran program.

4. Enter numbers; click Compute.
