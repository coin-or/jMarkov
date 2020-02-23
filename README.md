# jMarkov

jMarkov is a Java framework for Markov chain (MC) modeling that provides the user with the ability to define MCs from
the basic rules that determine their dynamics. From these rules, jMarkov explores and builds the state space and the
other MC parameters, which are then used to solve the MC and determine user-defined steady-state and transient performance
measures. 

jMarkov has the following modules:

-  *jMarkov Core Module*: Permits modeling large-scale finite Markov Chains
-  *jQBD*: Modeling of Quasi Birth and Death processes (QBDs)
-  *jPhase*: Modeling of Phase type distributions
-  *jMDP*: Modeling of Markov Decision Processes

## Download

The source code of jMarkov releases can be downloaded from the following link:
 - https://github.com/coin-or/jMarkov/releases

## Required packages

The jMarkov library needs some third party libraries for
all its features to work properly. All the libraries listed here
are open source, but keep in mind that not all of them are released
under the Eclipse Public license, as it is jMarkov. We recommend to verify
that the conditions of the third party libraries licenses listed 
here fit the requirements of your project.

1. COLT: http://acs.lbl.gov/ACSSoftware/colt/
2. JAMA: http://math.nist.gov/javanumerics/jama/
3. JCOMMON: http://www.jfree.org/jcommon/
4. JFREECHART: http://www.jfree.org/jfreechart/
5. MTJ: https://github.com/fommil/matrix-toolkits-java
6. QSOPT: http://www.math.uwaterloo.ca/~bico/qsopt/
7. SSJ: http://simul.iro.umontreal.ca/ssj-2/indexe.html.

## Documentation

The following module manuals are available:
 - [Core module and jQBD](https://coin-or.github.io/jMarkov/jMarkovManual.pdf)
 - [jPhase](https://coin-or.github.io/jMarkov/jPhaseManual.pdf)
 - [jMDP](https://coin-or.github.io/jMarkov/jMDPManual.pdf)

## Maintainer

jMarkov is maintained by:

- [Julio C. Góez](https://www.nhh.no/en/research-faculty/department-of-business-and-management-science/for/cv/goez--julio-c.aspx) (jgoez at jmarkov dot org): Core Module, jQBD
- [Juan F. Pérez](https://www.researchgate.net/profile/Juan_Perez225) (jfperez at jmarkov dot org): Core Module, jPhase
- [Daniel Silva](http://www.prism.gatech.edu/~dfsi3/) (dfsilva at jmarkov dot org): jMDP


## Contributors
- Juan P. Alvarado, 
- Diego Bello, 
- Rodrigo Cáliz, 
- Marco Cote, 
- Julio C. Góez (Norwegian School of Economics, Bergen, Norway), 
- Leonardo Lozano, 
- Juan F. Pérez (Universidad del Rosario, Bogotá, Colombia), 
- Germán Riaño, 
- Andrés Sarmiento Plata,
- Andrés Sarmiento Romero, 
- Daniel Silva (Auburn University, Auburn, Alabama, USA),
- Laura Vielma
