# Homework 3

The code and plots for the assignment are set up in a Jupyter Notebook file,
[_hw3.ipynb_](hw3.ipynb), along with the associated discussion.

The major components of this homework are Picard and Newton-Krylov solvers for
the _p_-Laplacian problem, along with observations on their convergence and the
number of iterations the methods use.

Issues:
* The Picard solver appears to work correctly, but gives an unexpected grid
refinement plot, easily a result of my misunderstandings or a bug.
* The Newton-Krylov solver always hits the iteration limit, although its
approximations look reasonable visually. My attempts to fix it led to it
breaking completely or outputting garbage, but lead candidates for problems
may be the setup of the convergence checks/norm use and incorrect handling of
the vectors/matrices and their shapes.
