# Simulation of pairwise interactive particles by Lennard Jones potential 
## Introduction
In this files there are some codes in C++ for the simulation of interactive particles using MonteCarlo techiques. The Equation Of State (EOS) is obtained using the both an NTV simulation and a NPT one
## NTV Ensemble 
We try to move displace a particle 
‘‘‘math
\vec{x_i}\Rightarrow \vec{x_i}+\vec{\xi}
‘‘‘
It is accepted either if it minimize the energy or if a random number is less than $‘e^{-\beta \Delta}‘$
This procedure is repeated N times (that is a Montecarlo step)
## NPT Ensemble
