# Numerical analysis of the convergence behaviour in horizontal and vertical Sato-Tate Distributions

A lot of the ideas are from https://wstein.org/talks/20071016-convergence/talk.pdf
and some of the code is taken from https://wstein.org/talks/20071016-convergence/uw/sato_tate-talk.pdf

As for the horizontal case, let s(x) be the area under the arc of the semicircle and f_c(x) the area of the histogram bars for primes up to c.
Then define 
<img src="https://render.githubusercontent.com/render/math?math=\Delta(C) = \sqrt{\int_{-1}^1 (s(x) - f_c(x))^2 \mathrm{d} x}">, i.e. the L2-norm of difference between our histogram and the convergence rate.

Finally we define <img src="https://render.githubusercontent.com/render/math?math=\theta(C) = -\log_C(\Delta(C)^2)">.

In horizontal_st.ipynb, we plot theta(C) for elliptic curves with different rank and see that theta depends on the rank.