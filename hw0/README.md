# Homework 0

The code and plots for the assignment are set up in a Jupyter Notebook file, [_hw0.ipynb_](hw0.ipynb), along with more in-depth discussion.


## `diffmat(x)`

For the interior points, I used a central difference approximation for the
derivative. For the endpoints, I used forward and backward differences, which
cause the approximation as a whole to have first-order accuracy.

## `diff2mat(x)`

For the interior points, I used a central difference approximation for the
derivative. In the second derivative case, the one-sided differences produced
the same approximations as the central differences for the adjacent points,
which performed poorly, particularly for points with larger spacing.

## Testing

I tested the approximations using grid refinement on the range [-1, 1], using
u(_x_) = cos(_x_) as my test function. Error values were calculated using an
infinity norm.

When testing with unevenly-spaced points, I drew random values from a uniform
distribution over [-1, 1).
