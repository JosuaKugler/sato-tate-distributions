# Numerical analysis of the convergence behaviour in horizontal and vertical Sato-Tate Distributions

A lot of the ideas are from https://wstein.org/talks/20071016-convergence/talk.pdf
and some of the code is taken from https://wstein.org/talks/20071016-convergence/uw/sato_tate-talk.pdf

The aim of this project is to provide numerical evidence for studying the convergence speed of the Sato-Tate-Distribution.

In horizontal_st.ipynb, we plot the convergence speed for elliptic curves with different rank and see that the convergence speed strongly depends on the rank.

In vertical_st.ipynb, we plot the convergence speed of the vertical Sato-Tate-Distribution as a function of the prime.
We conjecture that the squared L2-norm of the difference of the distribution and the semicircle is in O(p)