# Homework 2

The code and plots for the assignment are set up in a Jupyter Notebook file,
[_hw2.ipynb_](hw2.ipynb), along with the associated discussion.

The main feature is an adaptive Runge-Kutta integrator--based on the
`ode_rkexplicit` function from lecture-- that uses the Bogacki-Shampine method.
Given an error tolerance, the integrator will adjust the sizes of its time steps
to keep the estimated local truncation error near the tolerance.

The integrator is tested on the Van der Pol Equation and on _u'=Au_ with the
matrix exponential solution from lecture.
