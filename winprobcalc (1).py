


import numpy as np

def find_c_plus_star(N_value, tol=1e-8, max_iter=10000):
    """
    Finds the fixed point c* > 0.6 such that c = f(c) for a given N.

    Args:
        N_value: The value of N to use in the calculations.
        tol: The tolerance for the binary search.
        max_iter: The maximum number of iterations for the binary search.

    Returns:
        The value of c* if found, otherwise None.
    """
    def s_A(c):
        return 0.1 + c / 2

    def s_B(c):
        return 0.4

    def n(c):
        return (s_A(c) + s_B(c)) * N_value

    def m(c):
        return (s_A(c) - s_B(c)) / (s_A(c) + s_B(c))

    def f(c):
        # Handle potential division by zero or invalid values
        denominator = np.sqrt(n(c)) * m(c)
        if denominator == 0 or np.isnan(denominator) or np.isinf(denominator):
            return np.inf if denominator >= 0 else -np.inf # Return infinity or negative infinity to guide search

        return 1 / denominator

    def find_fixed_point_above_0_6_internal(tol=1e-6, max_iter=1000):
        c_low = 0.6
        c_high = 1.0  # Reasonable upper bound

        for _ in range(max_iter):
            c_mid = (c_low + c_high) / 2
            f_mid = f(c_mid)

            if abs(c_mid - f_mid) < tol:
                return c_mid  # Found fixed point with desired accuracy

            # Adjust search range based on whether c_mid is greater or less than f_mid
            if c_mid > f_mid:
                c_high = c_mid
            else:
                c_low = c_mid

        return None  # No solution found within max_iter

    return find_fixed_point_above_0_6_internal(tol, max_iter)

"""##Define a Function to find C-"""

def find_c_minus_star(N_value, c_start=0.4, c_end=0.6, increment=0.0000001):
    """
    Finds the value of c > 0.6 that minimizes the absolute difference |c - f(c)|
    using a linear search.

    Args:
        N_value: The value of N to use in the calculations.
        c_start: The starting value for the linear search.
        c_end: The ending value for the linear search.
        increment: The step size for the linear search.

    Returns:
        The value of c that minimizes |c - f(c)|.
    """
    def s_A(c):
        return 0.1 + c / 2

    def s_B(c):
        return 0.4

    def n(c):
        return (s_A(c) + s_B(c)) * N_value

    def m(c):
        return (s_B(c) - s_A(c)) / (s_A(c) + s_B(c))

    def f(c):
        denominator = np.sqrt(n(c)) * m(c)
        if denominator == 0 or np.isnan(denominator) or np.isinf(denominator):
            return np.inf if denominator >= 0 else -np.inf
        return 1 / denominator

    min_diff = float('inf')
    best_c = None

    for c in np.arange(c_start, c_end, increment):
        f_c = f(c)
        diff = abs(c - f_c)
        if diff < min_diff:
            min_diff = diff
            best_c = c

    return best_c

"""## Iterate and calculate c- and c+


"""

# 1. Define a list named N_values containing the specified values:
N_values = [1000, 10000, 100000, 500000, 1000000, 2500000, 50000000, 10000000]

# 2. Create an empty dictionary to store the results
c_star_values = {}

# 3. Loop through each value in N_values.
# 4. Inside the loop, call the find_c_star() function with the current N value.
# 5. Store the returned c* value, associated with the corresponding N, in the results container.
for N_value in N_values:
    c_plus_star = find_c_plus_star(N_value)
    c_minus_star = find_c_minus_star(N_value)
    c_star_values[N_value] = [c_plus_star, c_minus_star]

# Display the results
display(c_star_values)

"""## Plot the results


"""

from scipy.stats import binom
import matplotlib.pyplot as plt
import numpy as np

# Get the N values and c* values from the dictionary
N_values = list(c_star_values.keys())
c_star_values_list = list(c_star_values.values())

probabilities_A = []
probabilities_B = []
for N, c_star in zip(N_values, c_star_values_list):
    # Calculate n(c*)
    # s_A and s_B are defined in a previous cell, assuming they are accessible
    s_A_c_star_plus = 0.1 + c_star[0] / 2
    s_A_c_star_minus = 0.1 + c_star[1] / 2
    s_B_c_star = 0.4
    n_c_star_plus = (s_A_c_star_plus + s_B_c_star) * N
    n_c_star_minus = (s_A_c_star_minus + s_B_c_star) * N
    # Calculate the probability Pr(X > n(c*)/2)
    # Ensure n_c_star is an integer for binomial distribution
    n_c_star_plus = int(round(n_c_star_plus))
    n_c_star_minus = int(round(n_c_star_minus))

    # Calculate the new probability parameter p
    denominator_plus = s_A_c_star_plus + s_B_c_star
    if denominator_plus == 0:
        p_param_plus = 0  # Handle division by zero
    else:
        p_param_plus = s_A_c_star_plus / denominator_plus
    denominator_minus = s_A_c_star_minus + s_B_c_star
    if denominator_minus == 0:
        p_param_minus = 0  # Handle division by zero
    else:
        p_param_minus = s_B_c_star / denominator_minus
    # Calculate the probability Pr(X > n(c*)/2)
    # X ~ B(n=n(c*), p=s_A(c*)/(s_A(c*) + s_B(c*)))  <-- Corrected parameter
    # We need the survival function (1 - CDF)
    k_plus = int(round(n_c_star_plus / 2))
    k_minus = int(round(n_c_star_minus / 2))
    # Pr(X > k) = 1 - Pr(X <= k) = 1 - binom.cdf(k, n=n_c_star, p=p_param)
    winProbA_plus = 1 - binom.cdf(k_plus, n=n_c_star_plus, p=p_param_plus)
    winProbB_minus = 1 - binom.cdf(k_minus, n=n_c_star_minus, p=p_param_minus)
    probabilities_A.append(winProbA_plus)
    probabilities_B.append(winProbB_minus)
# Plot the probabilities vs N (using log scale for N as before)
plt.figure(figsize=(10, 6))
plt.plot(N_values, probabilities_A, marker='o', linestyle='-', markersize=8)
plt.plot(N_values, probabilities_B, marker='o', linestyle='-', markersize=8)
plt.xscale('log')
plt.xlabel('N (log scale)')
plt.ylabel('Pr(X > n(c*)/2)')
plt.title('Probability vs. N (Log Scale)')
plt.grid(True)

