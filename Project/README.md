# Simulation of pairwise interactive particles by Lennard Jones potential 
## Introduction
In this files there are some codes in C++ for the simulation of interactive particles using MonteCarlo techiques. The Equation Of State (EOS) is obtained using the both an NTV simulation and a NPT one
## NTV Ensemble 
We try to move displace a particle 
‘‘‘math
\vec{x_i}\Rightarrow \vec{x_i}+\vec{\xi}
‘‘‘
It is accepted either if it minimize the energy or if a random number is less than $$‘e^{-\beta \Delta}‘$$
This procedure is repeated N times (that is a Montecarlo step)
## NPT Ensemble
We either displace a particle or try to change the voume with a probability 1/N. The move is a random walk in the log(V) since it can't be negative. In doing so we need to take into account a new contribution when evaluating the enthalpy of the system
## $\mu$TV Ensemble
We fix the fugacity and the temperature. A Montecarlo Move consist either in displace the particle or change the number by increasing or decreasing it.
### Successive Umbrella Sampling 
A SUS algorithm is implemented for the evaluation of the equation of state. 