##Experimental results with following example: s_A(c) = 0.1 + c/2 and s_B(c) = 0.4##

1.  `winprobcalc.py`

Computes two fixed points, c⁺ and c⁻, for a given number of voters `N`, and then calculates win probabilities for two sides (A and B).
The results are plotted against `N` on a log scale.

* Binary search to find c⁺ > 0.6 such that `c = f(c)`
* Linear search to find  c⁻  minimizing `|c - f(c)|` in a given range


2.  `WinProb_Alpha.py`

Studies how alpha (a scaling parameter) and population size `N` influence equilibrium win probabilities. The parameter `c` is found using binary search for each `(alpha, N)` pair, and results are plotted.


python WinProb_Alpha.py

Output:  Plot probabilities vs. `alpha` for multiple `N`.


Dependencies:  

numpy
matplotlib
scipy

Installation: 

pip install numpy matplotlib scipy

Notes

* Both scripts are stand-alone and can be run independently.
* Parameter values (`alphaArray`, `N_values`, tolerances) can be tuned for different experiments.
* The mathematical model assumes specific definitions for `s_A(c)` and `s_B(c)`; adjust these functions for other scenarios.
