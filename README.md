rfunctions
==========

some cool R functions that I either wrote, stole, or modified

## Installation

**rfunctions** is not available on CRAN, but can be installed using the R package **devtools**. **rfunctions** can be installed with the following R code:


```r
devtools::install_github("jaredhuling/rfunctions")
library(rfunctions)
```


## Accelerated crossprod function

A project I've been working on requires fast evaluation of **X'X** for a design matrix **X**. I found a great example in the [paper](http://www.jstatsoft.org/v52/i05/paper) for [RcppEigen](http://cran.r-project.org/web/packages/RcppEigen/index.html) by Douglas Bates and Dirk Eddelbuettel for just such a thing. **RcppEigen** provides a simple and effective interface between R and the blazing-fast **Eigen** C++ library for numerical linear algebra. Their example uses **inline**, a nice tool for inline C++ code in R, and I a made a proper **R** function from that. The following showcases the speed of **Eigen**. Note that since **X'X** is symmetric, we only have to compute half of the values, which further reduces computation time. 


```r
n.obs <- 10000
n.vars <- 100

x <- matrix(rnorm(n.obs * n.vars), n.obs, n.vars)

library(microbenchmark)

microbenchmark(crossprodcpp(x), crossprod(x), times = 25L)
```

```
## Unit: milliseconds
##             expr   min    lq median    uq   max neval
##  crossprodcpp(x) 11.35 11.64  11.83 12.04 16.19    25
##     crossprod(x) 50.44 50.94  52.26 52.72 57.36    25
```

```r
all.equal(crossprodcpp(x), crossprod(x))
```

```
## [1] TRUE
```


```crossprodcpp``` can also compute a weighted cross product **X'WX** where **W** is a diagonal weight matrix


```r
ru <- runif(n.obs)
weights <- ru * (1 - ru)

microbenchmark(crossprodcpp(x, weights), crossprod(x, weights * x), times = 25L)
```

```
## Unit: milliseconds
##                       expr    min     lq median     uq    max neval
##   crossprodcpp(x, weights)  15.34  15.93  16.04  16.34  18.98    25
##  crossprod(x, weights * x) 102.98 105.23 106.27 113.21 193.27    25
```

```r
all.equal(crossprodcpp(x, weights), crossprod(x, weights * x))
```

```
## [1] TRUE
```



## Largest Singular Value Computation

The Lanczos algorithm is a well-known method for fast computation of extremal eigenvalues. The Golub-Kahan-Lanczos bidiagonalization algorithm is an extension of this to approximate the largest singular values of a matrix **X** from below. The function ```gklBidiag``` approximates the largest singular value of a matrix. Since GKL bidiagonalization is initialized from a random vector, we can compute a probabilistic upper bound for the singular value. The following compares the speed of ```gklBidiag``` and the implementation in the popular **Fortran** library **PROPACK** found in the **svd** package 


```r
library(svd)

v <- runif(ncol(x))  #initialization for GKL-bidiag
opts <- list(kmax = 30L)
microbenchmark(gklBidiag(x, v, maxit = 30L), propack.svd(x, neig = 1L, opts = opts))
```

```
## Unit: milliseconds
##                                    expr    min     lq median     uq    max
##            gklBidiag(x, v, maxit = 30L)  31.61  32.08  32.43  35.12  67.78
##  propack.svd(x, neig = 1L, opts = opts) 230.24 234.09 237.05 242.29 329.95
##  neval
##    100
##    100
```

```r

gklBidiag(x, v, maxit = 30L)$d - propack.svd(x, neig = 1L, opts = opts)$d
```

```
## [1] -2.52e-11
```



As ```gklBidiag``` also works on sparse matrices (of the ```SparseMatrix``` class from the **Matrix** package), I can showcase another function in **rfunctions**, ```simSparseMatrix```, which unsurprisingly simulates matrices with very few nonzero values. The nonzero values can either be all 1's or generated from a normal distribution. The level of sparsity of the simulated matrix can be specified



```r
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
```

```
## Unit: milliseconds
##                                  expr   min    lq median    uq   max neval
##  gklBidiag(x.s.b, v, maxit = 10L, 0L) 80.49 81.61  82.05 83.10 109.9   100
##  gklBidiag(x.s.c, v, maxit = 10L, 0L) 80.42 81.45  81.91 82.97 118.6   100
```

```r

gklBidiag(x.s.b, v, maxit = 10L, 0L)$d
```

```
## [1] 104.9
```

```r
gklBidiag(x.s.c, v, maxit = 10L, 0L)$d
```

```
## [1] 35.37
```


##Faster Addition/Subtraction of Matrices

This may seem pointless, but I wrote functions to add and subtract matrices. It turns out my functions are faster than using the ```+``` and ```-``` operators. I'm sure someone will be quick to point out why using my ```add()``` and ```subtract()``` functions is silly and a bad idea.


```r
A <- simSparseMatrix(sparsity = 0.99, dim = c(n.obs, n.vars), boolean = F)
B <- simSparseMatrix(sparsity = 0.99, dim = c(n.obs, n.vars), boolean = F)

microbenchmark(add(A, B), A + B)
```

```
## Unit: milliseconds
##       expr    min     lq median     uq   max neval
##  add(A, B)  77.15  79.07  80.05  88.85 191.2   100
##      A + B 227.95 328.97 334.87 339.39 535.6   100
```

```r
microbenchmark(subtract(A, B), A - B)
```

```
## Unit: milliseconds
##            expr    min     lq median     uq   max neval
##  subtract(A, B)  77.89  79.32  80.36  88.99 191.1   100
##           A - B 217.89 250.13 338.79 342.60 360.1   100
```

```r

all.equal(add(A, B), A + B)
```

```
## [1] TRUE
```

```r
all.equal(subtract(A, B), A - B)
```

```
## [1] TRUE
```


The ```add()``` and ```subtract()``` methods for dense matrices are slower than the corresponding operators, so they're only worth using when you have sparse matrices.


```r
n.obs <- 1000
n.vars <- 1000

A <- matrix(rnorm(n.obs * n.vars), n.obs, n.vars)
B <- matrix(rnorm(n.obs * n.vars), n.obs, n.vars)

microbenchmark(add(A, B), A + B)
```

```
## Unit: milliseconds
##       expr   min    lq median    uq   max neval
##  add(A, B) 6.044 6.347  6.425 6.566 24.74   100
##      A + B 2.975 3.116  3.149 3.209 20.44   100
```

```r
microbenchmark(subtract(A, B), A - B)
```

```
## Unit: milliseconds
##            expr   min    lq median    uq   max neval
##  subtract(A, B) 6.060 6.331  6.432 6.596 22.09   100
##           A - B 2.997 3.117  3.157 3.253 21.50   100
```

```r

all.equal(add(A, B), A + B)
```

```
## [1] TRUE
```

```r
all.equal(subtract(A, B), A - B)
```

```
## [1] TRUE
```



