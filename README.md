# DPR1Eigensystem
A Mathematica implementation of an eigensolver for symmetric, positive-definite rank-one modifications of diagonal matrices.

Derived from the following paper:
https://doi.org/10.1016/j.laa.2015.09.025


# Installation and usage

Simply copy the file `DPR1Eigensystem.m` anywhere, where _Mathematica_ can find it, e.g., unto the path

     FileNameJoin[{$BaseDirectory, "Applications"}]

Run ``PR1Eigensystem`*`` to see the packages contents.
The main routine is just called `DPR1Eigensystem`. Run ``?DPR1Eigensystem`` to get its syntax information. 

Also see

    https://mathematica.stackexchange.com/a/280669/38178
    https://mathematica.stackexchange.com/a/280683/38178
    
for context and usage examples.

# Disclaimer

All this is done under the assumption that no two elements of `diag` coincide. The Cauchy interlacing theorem implies that no eigenvalue has a multiplicity. Not sure how to deal with multiplicities...
