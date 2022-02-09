.. _use_Cases:

Use Cases
********************

SAR example
============
This script (``examples/SAR_example.py``) illustrates the use of InstruPy to produce the data-metrics for SAR models with 
different possible set of configurations. The configurations shown are:

1.  L- band SAR, Dual polarization (SMAP pulse configuration), Fixed swath size (25km)
2.  L- band SAR, Dual polarization (AIRSAR pulse configuration), Full swath
3.  P-band SAR, Single  polarization, Fixed swath size (50km)

Multi Objective Optimization (MOO) of SAR
==========================================

This script is used to illustrate the use of InstruPy for the purpose of optimal instrument design.
Specifically in this example involves SAR instrument design based on a multi-objective optimization algorithm (NSGA-II), enabled by the 
``pymoo 0.4.2`` package. A baseline instrument is defined, and some of the instrument parameters are kept as optimization variables.
Multiple objectives are defined in terms of the data-metrics produced by the InstruPy package.
The result is a multi-dimensional Pareto-curve illustrating the trade-offs between different instrument designs. Any point on the
Pareto curve is optimal and a suitable point can be selected by the user for implementation.

