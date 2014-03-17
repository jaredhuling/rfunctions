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

A project I've been working on requires fast evaluation of $$X^TX$$ for a design matrix $$X$$. I found a great example in the paper for [RcppEigen](http://www.jstatsoft.org/v52/i05/paper) by Douglas Bates and Dirk Eddelbuettel for just such a thing. **RcppEigen** provides a simple and effective interfact between R and the blazing-fast **Eigen** C++ library for numerical linear algebra. Their example uses **inline**, a nice tool for inline C++ code in R, and I a made a proper **R** function from that. The following showcases the speed of **Eigen**. Note that since $$X^TX$$ is symmetric, we only have to compute half of the values, which further reduces computation time. 


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
##  crossprodcpp(x) 14.00 14.21  14.49 14.82  23.8    25
##     crossprod(x) 61.84 63.26  64.54 71.42 193.3    25
```

```r
all.equal(crossprodcpp(x), crossprod(x))
```

```
## [1] TRUE
```


```crossprodcpp``` can also compute a weighted cross product $$X^T W X$$ where $$W$$ is a diagonal weight matrix


```r
ru <- runif(n.obs)
weights <- ru * (1 - ru)

microbenchmark(crossprodcpp(x, weights), crossprod(x, weights * x), times = 25L)
```

```
## Unit: milliseconds
##                       expr    min     lq median     uq    max neval
##   crossprodcpp(x, weights)  19.12  20.76  23.53  28.22  78.89    25
##  crossprod(x, weights * x) 125.96 129.44 142.63 175.67 319.81    25
```

```r
all.equal(crossprodcpp(x, weights), crossprod(x, weights * x))
```

```
## [1] TRUE
```



## Largest Singular Value Computation

The Lanczos algorithm is a well-known method for fast computation of extremal eigenvalues. The Golub-Kahan-Lanczos bidiagonalization algorithm is an extension of this to approximate the largest singular values of a matrix $$X$$ from below. The function ```gklBidiag``` approximates the largest singular value of a matrix. Since GKL bidiagonalization is initialized from a random vector, we can compute a probabilistic upper bound for the singular value. The following compares the speed of ```gklBidiag``` and the implementation in the popular **Fortran** library **PROPACK** found in the **svd** package 


```r
library(svd)

v <- runif(ncol(x))  #initialization for GKL-bidiag
opts <- list(kmax = 30L)
microbenchmark(gklBidiag(x, v, maxit = 30L), propack.svd(x, neig = 1L, opts = opts))
```

```
## Unit: milliseconds
##                                    expr    min     lq median     uq    max
##            gklBidiag(x, v, maxit = 30L)  36.49  36.98  38.08  42.71  61.51
##  propack.svd(x, neig = 1L, opts = opts) 350.00 357.09 364.51 400.14 567.51
##  neval
##    100
##    100
```

```r

gklBidiag(x, v, maxit = 30L)$d - propack.svd(x, neig = 1L, opts = opts)$d
```

```
## [1] -1.252e-11
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
##  gklBidiag(x.s.b, v, maxit = 10L, 0L) 94.84 96.34  97.72 105.8 183.1   100
##  gklBidiag(x.s.c, v, maxit = 10L, 0L) 94.97 96.96  98.79 114.9 230.3   100
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
## [1] 35.66
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
##       expr    min     lq median    uq   max neval
##  add(A, B)  92.39  95.04  103.7 111.9 222.3   100
##      A + B 270.58 385.46  393.9 404.0 711.0   100
```

```r
microbenchmark(subtract(A, B), A - B)
```

```
## Unit: milliseconds
##            expr    min     lq median    uq   max neval
##  subtract(A, B)  93.53  95.93  104.3 111.9 220.1   100
##           A - B 265.88 384.31  398.6 445.6 560.3   100
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
##  add(A, B) 7.031 7.307  7.779 9.599 30.01   100
##      A + B 3.368 3.533  3.676 4.498 33.05   100
```

```r
microbenchmark(subtract(A, B), A - B)
```

```
## Unit: milliseconds
##            expr   min    lq median    uq   max neval
##  subtract(A, B) 6.982 7.179  7.297 7.552 27.06   100
##           A - B 3.423 3.490  3.551 3.749 22.57   100
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



