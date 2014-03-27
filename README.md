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
##  crossprodcpp(x) 11.53 11.84  12.23 13.02 17.37    25
##     crossprod(x) 51.39 52.70  55.06 56.78 61.04    25
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
##                       expr    min    lq median     uq    max neval
##   crossprodcpp(x, weights)  15.67  16.2  16.45  16.57  17.68    25
##  crossprod(x, weights * x) 101.16 105.7 107.18 109.19 122.84    25
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
##                                  expr   min    lq median    uq    max
##  gklBidiag(x.s.b, v, maxit = 10L, 0L) 81.31 82.32  82.81 83.71  94.72
##  gklBidiag(x.s.c, v, maxit = 10L, 0L) 80.92 82.53  83.14 84.14 106.05
##  neval
##    100
##    100
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
## [1] 35.25
```


## Fast Moore-Penrose Generalized Inverse

The speed of ```MASS::ginv()``` leaves much to be desired, as it called ```svd()``` in order to compute the Moore-Penrose generalized inverse of a matrix. I came across this [paper](), which provides a faster algorithm for the M-P inverse. I've implemented it using RcppEigen. This is useful for least squares problems when the matrix is less than full rank, like below


```r
n <- 1000
p <- 500

x <- matrix(rnorm(n * (p - 1)), n, p)
x <- cbind(x, rowMeans(x))

## compute X'X
xpx <- crossprodcpp(x)

library(MASS)

## compute generalized inverse of X'X
microbenchmark(geninv(xpx), ginv(xpx))
```

```
## Unit: milliseconds
##         expr   min    lq median    uq   max neval
##  geninv(xpx) 186.4 189.0  190.8 195.1 227.1   100
##    ginv(xpx) 501.2 510.7  517.5 523.4 550.9   100
```

```r


inv <- geninv(xpx)
inv2 <- ginv(xpx)

## check if we have computed the M-P generalized inverse
all.equal(xpx, xpx %*% inv %*% xpx)
```

```
## [1] TRUE
```

```r
all.equal(inv, inv2)
```

```
## [1] TRUE
```


## Faster Addition/Subtraction of Matrices

This may seem pointless, but I wrote functions to add and subtract matrices. It turns out my functions are faster than using the ```+``` and ```-``` operators. I'm sure someone will be quick to point out why using my ```add()``` and ```subtract()``` functions is silly and a bad idea.


```r
A <- simSparseMatrix(sparsity = 0.99, dim = c(n.obs, n.vars), boolean = F)
B <- simSparseMatrix(sparsity = 0.99, dim = c(n.obs, n.vars), boolean = F)

microbenchmark(add(A, B), A + B)
```

```
## Unit: milliseconds
##       expr    min     lq median    uq max neval
##  add(A, B)  78.31  80.31  81.74  90.2 200   100
##      A + B 210.79 238.46 335.76 340.6 362   100
```

```r
microbenchmark(subtract(A, B), A - B)
```

```
## Unit: milliseconds
##            expr   min     lq median     uq   max neval
##  subtract(A, B)  78.3  80.28  81.66  85.69 104.8   100
##           A - B 215.0 339.42 343.98 349.33 378.6   100
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
##  add(A, B) 5.008 6.355  6.453 6.628 25.57   100
##      A + B 1.734 3.093  3.172 3.271 17.04   100
```

```r
microbenchmark(subtract(A, B), A - B)
```

```
## Unit: milliseconds
##            expr   min    lq median    uq   max neval
##  subtract(A, B) 4.988 6.286  6.479 6.727 18.48   100
##           A - B 1.728 3.127  3.192 3.278 17.50   100
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



