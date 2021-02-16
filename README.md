# Gröbner Basis VooDoo

An extension of [F<sub>5</sub>](https://asdm.gmbh/2020/11/20/f5/): Given a system of multivariate polynomials, compute a [Gröbner Basis](https://asdm.gmbh/2020/08/09/introduction-to-gb-attacks-on-aoc/) and the respective vectors of origin.

## Vectors of Origin
Given input system *F* = (*f*<sub>0</sub>, …, *f*<sub>*n*-1</sub>), a Gröbner basis algorithm outputs system *G* = (*g*<sub>0</sub>, …, *g*<sub>*m*-1</sub>) spanning the same ideal &langle;*F*&rangle;. Even for mildly complex systems *F*, it is generally not at all obvious how to (weightedly) combine the *f*<sub>*i*</sub> to arrive at any given *g*<sub>*j*</sub>.

A vector *v* = (*v*<sub>0</sub>, …, *v*<sub>*n*-1</sub>) is called **Vector of Origin** (voo) for *g*<sub>*j*</sub> if the inner product of *F* and *v* equals *g*<sub>*j*</sub>. That is, Σ<sub>*i*</sub> *f*<sub>*i*</sub>·*v*<sub>*i*</sub> = *g*<sub>*j*</sub>.

Arranging all *m* vectors of origin into a matrix *V*, we have *F*·*V* = *G*. This *n*×*m* matrix *V*, where each entry is a multivariate polynomial, contains a lot of juicy information about *F*!

## Code
The [F<sub>5</sub> algorithm](https://asdm.gmbh/2020/11/20/f5/) due to Jean-Charles Faugère implicitly uses vectors of origin – this code simply makes them explicit. It is an extension of the F5 and F4 [sagemath implementation](https://bitbucket.org/malb/research-snippets/src/master/f5.py) by Martin Albrecht and John Perry.
