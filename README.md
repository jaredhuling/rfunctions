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
##  crossprodcpp(x) 11.88 13.09  15.33 16.82 19.34    25
##     crossprod(x) 51.32 56.07  60.98 62.41 72.67    25
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
##                       expr    min     lq median    uq    max neval
##   crossprodcpp(x, weights)  15.93  16.34  16.68  17.2  24.61    25
##  crossprod(x, weights * x) 105.36 106.59 110.08 113.9 134.62    25
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
##                                  expr   min    lq median     uq   max
##  gklBidiag(x.s.b, v, maxit = 10L, 0L) 82.95 84.54  85.59  92.48 246.8
##  gklBidiag(x.s.c, v, maxit = 10L, 0L) 83.21 85.00  87.00 101.99 243.1
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
## [1] 35.29
```


## Fast Moore-Penrose Generalized Inverse

The speed of ```MASS::ginv()``` leaves much to be desired, as it calls ```svd()``` in order to compute the Moore-Penrose generalized inverse of a matrix. I came across this [paper](), which provides a faster algorithm for the M-P inverse. I've implemented it using RcppEigen. This is useful for least squares problems when the matrix is less than full rank, like below


```r
n <- 1000
p <- 500

x <- matrix(rnorm(n * (p - 1)), n, p - 1)
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
##  geninv(xpx) 185.9 190.4  192.1 199.6 256.2   100
##    ginv(xpx) 498.7 511.6  522.2 534.7 711.0   100
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
n.obs <- 1e+05
n.vars <- 500

A <- simSparseMatrix(sparsity = 0.99, dim = c(n.obs, n.vars), boolean = F)
B <- simSparseMatrix(sparsity = 0.99, dim = c(n.obs, n.vars), boolean = F)

microbenchmark(add(A, B), A + B, times = 25L)
```

```
## Unit: milliseconds
##       expr   min     lq median     uq   max neval
##  add(A, B) 37.59  40.08  40.93  44.32 161.2    25
##      A + B 92.69 100.69 105.72 110.32 218.2    25
```

```r
microbenchmark(subtract(A, B), A - B, times = 25L)
```

```
## Unit: milliseconds
##            expr   min    lq median    uq   max neval
##  subtract(A, B) 37.42 39.22  40.22  45.8 155.2    25
##           A - B 94.09 97.81 108.23 112.1 120.9    25
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
##  add(A, B) 4.801 6.286  6.432 7.201 24.76   100
##      A + B 1.677 3.018  3.207 3.472 19.32   100
```

```r
microbenchmark(subtract(A, B), A - B)
```

```
## Unit: milliseconds
##            expr   min    lq median    uq   max neval
##  subtract(A, B) 4.838 6.260  6.347 6.512 19.38   100
##           A - B 1.845 2.996  3.054 3.198 20.56   100
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



