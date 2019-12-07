3Th
================
Kalyani Cauwenberghs

This is Header 1
================

This is my work from **BIMM 143** 3Th.

``` r
# this is to demo a code chunk
plot(1:10)
```

![](3Th_files/figure-markdown_github/unnamed-chunk-1-1.png) \#\# Practice reading files (again) Here I practice reading 3 different files

``` r
table1=read.table("test1.txt", sep=",",header = T)
table2=read.table("test2.txt", sep="$",header = T)
table3=read.table("test3.txt", sep="",header = F)
```

``` r
add<-function(x,y=1){
  #sum the input of x and y
  x+y
}
```

``` r
add(x=1,y=4)
```

    ## [1] 5

``` r
add(1,4)
```

    ## [1] 5

``` r
add(1)
```

    ## [1] 2

``` r
add(c(1,2,3))
```

    ## [1] 2 3 4

``` r
add(c(1,2,3,4))
```

    ## [1] 2 3 4 5

Our function is vectorized

``` r
add(c(1,2,3))
```

    ## [1] 2 3 4

``` r
add(c(1,2,3,4))
```

    ## [1] 2 3 4 5

A new function to rescale data

``` r
rescale <- function(x) {
  rng <-range(x)
  (x - rng[1]) / (rng[2] - rng[1])
}
```

``` r
rescale(1:10)
```

    ##  [1] 0.0000000 0.1111111 0.2222222 0.3333333 0.4444444 0.5555556 0.6666667
    ##  [8] 0.7777778 0.8888889 1.0000000

test different vector:

``` r
rescale(c(1,2,NA,3,10))
```

    ## [1] NA NA NA NA NA

try to fix function:

``` r
x<-c(1,2,NA,3,10)
rng<-range(x,na.rm = T)
rng
```

    ## [1]  1 10

``` r
rescale2 <- function(x) {
  rng <-range(x,na.rm=T)
  (x - rng[1]) / (rng[2] - rng[1])
}
```

``` r
rescale2(c(1,2,NA,3,10))
```

    ## [1] 0.0000000 0.1111111        NA 0.2222222 1.0000000

``` r
rescale3 <- function(x, na.rm=TRUE, plot=FALSE) {
   rng <-range(x, na.rm=na.rm)
   print("Hello")
   answer <- (x - rng[1]) / (rng[2] - rng[1])
   print("is it me you are looking for?")
   if(plot) {
    plot(answer, typ="b", lwd=4)
   }
   print("I can see it in ...")
   return(answer)
   #by default, R will return the last thing that was calculated
   #returning will stop the function
}
```

``` r
rescale3(1:10,plot = T)
```

    ## [1] "Hello"
    ## [1] "is it me you are looking for?"

![](3Th_files/figure-markdown_github/unnamed-chunk-13-1.png)

    ## [1] "I can see it in ..."

    ##  [1] 0.0000000 0.1111111 0.2222222 0.3333333 0.4444444 0.5555556 0.6666667
    ##  [8] 0.7777778 0.8888889 1.0000000

Section 2 of Hands on Sheet
===========================

Install **bio3d** for sequence and structure analysis

``` r
#do this in console, in order to not install the package many times
#install.packages("bio3d")
```

``` r
# Can you improve this analysis code?
library(bio3d)
s1 <- read.pdb("4AKE") # kinase with drug
```

    ##   Note: Accessing on-line PDB file

``` r
s2 <- read.pdb("1AKE") # kinase no drug
```

    ##   Note: Accessing on-line PDB file
    ##    PDB has ALT records, taking A only, rm.alt=TRUE

``` r
s3 <- read.pdb("1E4Y") # kinase with drug
```

    ##   Note: Accessing on-line PDB file

``` r
s1.chainA <- trim.pdb(s1, chain="A", elety="CA")
s2.chainA <- trim.pdb(s2, chain="A", elety="CA")
s3.chainA <- trim.pdb(s1, chain="A", elety="CA")
s1.b <- s1.chainA$atom$b
s2.b <- s2.chainA$atom$b
s3.b <- s3.chainA$atom$b
plotb3(s1.b, sse=s1.chainA, typ="l", ylab="Bfactor")
```

![](3Th_files/figure-markdown_github/unnamed-chunk-15-1.png)

``` r
plotb3(s2.b, sse=s2.chainA, typ="l", ylab="Bfactor")
```

![](3Th_files/figure-markdown_github/unnamed-chunk-15-2.png)

``` r
plotb3(s3.b, sse=s3.chainA, typ="l", ylab="Bfactor")
```

![](3Th_files/figure-markdown_github/unnamed-chunk-15-3.png)
