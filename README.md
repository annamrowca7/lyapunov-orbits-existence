# Proof of the existence of vertical Lyapunov orbits in the restricted circular three-body problem.
## Project description 
The purpose of this repository is a computer-assisted proof of the existence of vertical Lyapunov orbits in the restricted circular three-body problem. Project is based on master thesis.

## Requirements
To run the program, the **CAPD (Computer Assisted Proofs in Dynamics)** library must be installed:  
[https://github.com/CAPDgroup/CAPD](https://github.com/CAPDgroup/CAPD)  

The CAPD library provides tools for:  
- interval arithmetic,
- rigorous numerical analysis of dynamical systems,  
- computer-assisted proofs in dynamical systems.

---

## Installation
Clone and build the **CAPD** library:

```bash
git clone https://github.com/CAPDgroup/CAPD.git
cd CAPD
mkdir build && cd build
cmake ..
make
```
For detailed decription on how to build the library see
http://capd.ii.uj.edu.pl/html/capd_compilation.html
