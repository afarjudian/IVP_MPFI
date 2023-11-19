# IVP_MPFI
Sound and complete interval method for solving ordinary differential equation (ODE) initial value problems (IVPs).

The code is an implementation of the method presented in the paper:

- Edalat, A., Farjudian, A., and Li, Y. “Recursive Solution of Initial Value Problems with Temporal Discretization”. In: Theoretical Computer Science (2023)
  - https://doi.org/10.1016/j.tcs.2023.114221

and it may not make much sense without reading the above article first.



### Environment
- Ubuntu 22.04 LTS

In principle, it should run on any system with a C++ compiler. No advanced feature of C++ is used in the code.

### Dependency
- libboost-all-dev
- libmpfi-dev


### Compile and Run:
On Ubuntu, the following command should compile the program:

- g++ -o ivp_mpfi *.cpp -lmpfi -lmpfr

and then run:

- ./ivp_mpfi




