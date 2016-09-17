Calculates integral Isubn = integral [0,1] of (x^n e^-x dx).

Does so by using a modified algorithm and a high quality general integrator where absolute and relative error is 10^-9. (See results).

Measures time per function call of the two functions, repeated so that total running time for each is between 1 and 2 seconds.

Calculates the ratio of the general integrator to the algorithm. (High quality integrator is about 815 times faster.)
