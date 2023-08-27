# IVP_MPFI
Sound and complete interval method for solving initial value problems.

The code is an implementation of the method presented in the paper:

- Edalat, A., Farjudian, A., and Li, Y. "Recursive Solution of Initial Value Problems with Temporal Discretization." 2023. arXiv: 2301.03920 [math.NA]

The code is written in C++. No advanced feature of C++ is used in the code. It may not make much sense without reading the above article.

The code has only been run on Ubuntu Linux, although in principle, it should run on any system with a C++ compiler. It needs the following libraries:

libboost-all-dev
libmpfi-dev

On Ubuntu, the following command should compile the program:

g++ -o ivp_mpfi *.cpp -lmpfi -lmpfr

and then run:

./ivp_mpfi




