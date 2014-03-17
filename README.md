rfunctions
==========

some cool R functions that I either wrote, stole, or modified

## Installation

**rfunctions** is not available on CRAN, but can be installed using the R package **devtools**. **rfunctions** can be installed with the following R code:


{% highlight r %}
devtools::install_github("jaredhuling/rfunctions")
library(rfunctions)
{% endhighlight %}


## Accelerated crossprod function

A project I've been working on requires fast evaluation of $$X^TX$$ for a design matrix $$X$$. I found a great example in the paper for [RcppEigen](http://www.jstatsoft.org/v52/i05/paper) by Douglas Bates and Dirk Eddelbuettel for just such a thing. **RcppEigen** provides a simple and effective interfact between R and the blazing-fast **Eigen** C++ library for numerical linear algebra. Their example uses **inline**, a nice tool for inline C++ code in R, and I a made a proper **R** function from that. The following showcases the speed of **Eigen**. Note that since $$X^TX$$ is symmetric, we only have to compute half of the values, which further reduces computation time. 


{% highlight r %}
n.obs <- 10000
n.vars <- 100

x <- matrix(rnorm(n.obs * n.vars), n.obs, n.vars)

library(microbenchmark)

microbenchmark(crossprodcpp(x), crossprod(x), times = 25L)
{% endhighlight %}



{% highlight text %}
## Unit: milliseconds
##             expr   min    lq median    uq   max neval
##  crossprodcpp(x) 14.05 14.35  14.94 15.97 23.72    25
##     crossprod(x) 61.83 62.45  62.91 67.37 75.84    25
{% endhighlight %}



{% highlight r %}
all.equal(crossprodcpp(x), crossprod(x))
{% endhighlight %}



{% highlight text %}
## [1] TRUE
{% endhighlight %}


```crossprodcpp``` can also compute a weighted cross product $$X^T W X$$ where $$W$$ is a diagonal weight matrix


{% highlight r %}
ru <- runif(n.obs)
weights <- ru * (1 - ru)

microbenchmark(crossprodcpp(x, weights), crossprod(x, weights * x), times = 25L)
{% endhighlight %}



{% highlight text %}
## Unit: milliseconds
##                       expr    min    lq median    uq   max neval
##   crossprodcpp(x, weights)  19.01  19.5  19.76  21.0  24.0    25
##  crossprod(x, weights * x) 123.00 124.4 125.25 127.3 140.1    25
{% endhighlight %}



{% highlight r %}
all.equal(crossprodcpp(x, weights), crossprod(x, weights * x))
{% endhighlight %}



{% highlight text %}
## [1] TRUE
{% endhighlight %}



## Largest Singular Value Computation

The Lanczos algorithm is a well-known method for fast computation of extremal eigenvalues. The Golub-Kahan-Lanczos bidiagonalization algorithm is an extension of this to approximate the largest singular values of a matrix $$X$$ from below. The function ```gklBidiag``` approximates the largest singular value of a matrix. Since GKL bidiagonalization is initialized from a random vector, we can compute a probabilistic upper bound for the singular value. The following compares the speed of ```gklBidiag``` and the implementation in the popular **Fortran** library **PROPACK** found in the **svd** package 


{% highlight r %}
library(svd)

v <- runif(ncol(x))  #initialization for GKL-bidiag
opts <- list(kmax = 30L)
microbenchmark(gklBidiag(x, v, maxit = 30L), propack.svd(x, neig = 1L, opts = opts))
{% endhighlight %}



{% highlight text %}
## Error: "propack_svd" not resolved from current namespace (svd)
{% endhighlight %}



{% highlight r %}

gklBidiag(x, v, maxit = 30L)$d - propack.svd(x, neig = 1L, opts = opts)$d
{% endhighlight %}



{% highlight text %}
## Error: "propack_svd" not resolved from current namespace (svd)
{% endhighlight %}



As ```gklBidiag``` also works on sparse matrices (of the ```SparseMatrix``` class from the **Matrix** package), I can showcase another function in **rfunctions**, ```simSparseMatrix```, which unsurprisingly simulates matrices with very few nonzero values. The nonzero values can either be all 1's or generated from a normal distribution. The level of sparsity of the simulated matrix can be specified



{% highlight r %}
n.obs <- 1e+05
n.vars <- 1000

# simulate a very sparse matrix (this matrix has many zeros and few ones)
x.s.b <- simSparseMatrix(sparsity = 0.99, dim = c(n.obs, n.vars), boolean = T)
x.s.c <- simSparseMatrix(sparsity = 0.99, dim = c(n.obs, n.vars), boolean = F)
v <- runif(n.vars)

# reorthogonalization sometimes leads to higher accuracy. it helps correct
# for floating-point errors
microbenchmark(gklBidiag(x.s.b, v, maxit = 10L, 0L), gklBidiag(x.s.c, v, maxit = 10L, 
    0L))
{% endhighlight %}



{% highlight text %}
## Unit: milliseconds
##                                  expr   min    lq median    uq   max neval
##  gklBidiag(x.s.b, v, maxit = 10L, 0L) 95.61 98.09  99.94 103.6 149.8   100
##  gklBidiag(x.s.c, v, maxit = 10L, 0L) 95.44 98.96 100.40 103.6 126.0   100
{% endhighlight %}



{% highlight r %}

gklBidiag(x.s.b, v, maxit = 10L, 0L)$d
{% endhighlight %}



{% highlight text %}
## [1] 104.9
{% endhighlight %}



{% highlight r %}
gklBidiag(x.s.c, v, maxit = 10L, 0L)$d
{% endhighlight %}



{% highlight text %}
## [1] 35.37
{% endhighlight %}


##Faster Addition/Subtraction of Matrices

This may seem pointless, but I wrote functions to add and subtract matrices. It turns out my functions are faster than using the ```+``` and ```-``` operators. I'm sure someone will be quick to point out why using my ```add()``` and ```subtract()``` functions is silly and a bad idea.


{% highlight r %}
A <- simSparseMatrix(sparsity = 0.99, dim = c(n.obs, n.vars), boolean = F)
B <- simSparseMatrix(sparsity = 0.99, dim = c(n.obs, n.vars), boolean = F)

microbenchmark(add(A, B), A + B)
{% endhighlight %}



{% highlight text %}
## Unit: milliseconds
##       expr    min     lq median    uq max neval
##  add(A, B)  89.09  94.41  97.67 106.8 220   100
##      A + B 249.16 366.64 376.51 391.0 575   100
{% endhighlight %}



{% highlight r %}
microbenchmark(subtract(A, B), A - B)
{% endhighlight %}



{% highlight text %}
## Unit: milliseconds
##            expr    min     lq median    uq   max neval
##  subtract(A, B)  91.62  94.47  96.69 105.4 213.7   100
##           A - B 252.04 370.01 376.08 386.3 602.1   100
{% endhighlight %}



{% highlight r %}

all.equal(add(A, B), A + B)
{% endhighlight %}



{% highlight text %}
## [1] TRUE
{% endhighlight %}



{% highlight r %}
all.equal(subtract(A, B), A - B)
{% endhighlight %}



{% highlight text %}
## [1] TRUE
{% endhighlight %}


The ```add()``` and ```subtract()``` methods for dense matrices are slower than the corresponding operators, so they're only worth using when you have sparse matrices.


{% highlight r %}
n.obs <- 1000
n.vars <- 1000

A <- matrix(rnorm(n.obs * n.vars), n.obs, n.vars)
B <- matrix(rnorm(n.obs * n.vars), n.obs, n.vars)

microbenchmark(add(A, B), A + B)
{% endhighlight %}



{% highlight text %}
## Unit: milliseconds
##       expr   min    lq median    uq   max neval
##  add(A, B) 5.539 5.899  7.218 7.348 24.17   100
##      A + B 1.847 2.311  3.516 3.587 22.61   100
{% endhighlight %}



{% highlight r %}
microbenchmark(subtract(A, B), A - B)
{% endhighlight %}



{% highlight text %}
## Unit: milliseconds
##            expr   min    lq median    uq   max neval
##  subtract(A, B) 5.498 7.260  7.813 9.428 46.53   100
##           A - B 1.823 2.396  3.576 4.134 37.34   100
{% endhighlight %}



{% highlight r %}

all.equal(add(A, B), A + B)
{% endhighlight %}



{% highlight text %}
## [1] TRUE
{% endhighlight %}



{% highlight r %}
all.equal(subtract(A, B), A - B)
{% endhighlight %}



{% highlight text %}
## [1] TRUE
{% endhighlight %}



