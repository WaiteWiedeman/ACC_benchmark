clear; close all; clc;

confidence = 0.99

c = 2
n = 200

upper1 = upper_bound_confidence_binomial_approximation(c, n, confidence)
lower1 = lower_bound_confidence_binomial_approximation(c, n, confidence)

confidence = 0.99

c = 12
n = 100

upper2 = upper_bound_confidence_binomial_approximation(c, n, confidence)
lower2 = lower_bound_confidence_binomial_approximation(c, n, confidence)