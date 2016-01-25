Finding zeros of Hankel function
================================

We had great problems to find zeros of the Hankel/Bessel function for high
orders. In the current implementation of the SFS Toolbox we solved it via
calculating them with the [Multiprecision Computing
Toolbox](http://www.advanpix.com/) and storing them under
https://github.com/sfstoolbox/data/tree/master/sphbesselh_zeros

This folder contains code that was a first try to calculate the zeros without
the Multiprecision Toolbox, see also:
http://www.mathworks.com/matlabcentral/fileexchange/22725-variable-precision-integer-arithmetic

I cannot say if the code is actual working ;)
