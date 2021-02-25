Week 6 Lab
=============
  
Today we will do a short exercise to illustrate the permutation method of dealing with multiple comparisons.

First, we will simulate 10 groups of data (say, heights) from the *same* distribution using the normal distribution and put the data for each group in a list for easy access:


```r
data.all<-list()
for (i in 1:10)
  {
    data.all[[i]]<-rnorm(10)  #Note the double brackets for a list
  }
```

Now we can compare each group pairwise using a t.test.


```r
p.values<-matrix(ncol=10,nrow=10)
for (i in 1:9)
  {
    for (j in (i+1):10)
      {
        p.values[i,j]<-t.test(data.all[[i]],data.all[[j]])$p.value 
      }
  }
p.values
```

```
##       [,1]      [,2]      [,3]      [,4]      [,5]       [,6]       [,7]
##  [1,]   NA 0.2298995 0.2858893 0.2124564 0.2828324 0.04885913 0.55854000
##  [2,]   NA        NA 0.7525434 0.9881544 0.9613593 0.60438727 0.34709633
##  [3,]   NA        NA        NA 0.7517666 0.7501809 0.31474742 0.44624492
##  [4,]   NA        NA        NA        NA 0.9499370 0.57087351 0.32077318
##  [5,]   NA        NA        NA        NA        NA 0.70426570 0.41398673
##  [6,]   NA        NA        NA        NA        NA         NA 0.05008738
##  [7,]   NA        NA        NA        NA        NA         NA         NA
##  [8,]   NA        NA        NA        NA        NA         NA         NA
##  [9,]   NA        NA        NA        NA        NA         NA         NA
## [10,]   NA        NA        NA        NA        NA         NA         NA
##             [,8]      [,9]      [,10]
##  [1,] 0.60222145 0.3730972 0.11392759
##  [2,] 0.09889786 0.7275846 0.88607658
##  [3,] 0.10639723 0.9446541 0.57763690
##  [4,] 0.08484414 0.7269590 0.86559272
##  [5,] 0.14702535 0.7259625 0.94723227
##  [6,] 0.01238666 0.3420798 0.64173422
##  [7,] 0.20840332 0.5774271 0.14869129
##  [8,]         NA 0.1691917 0.03360092
##  [9,]         NA        NA 0.57643595
## [10,]         NA        NA         NA
```

Now we can see how many of these p.values are "significant". We know these are false positives, because all the data were generated from the same distribution.


```r
false.positives<-sum(p.values<0.05,na.rm=T)
false.positives
```

```
## [1] 3
```

We could correct this using the Bonferonni method:


```r
k<-45
new.threshold.B<-0.05/k
new.threshold.B
```

```
## [1] 0.001111111
```

```r
false.positives.B<-sum(p.values<new.threshold.B,na.rm=T)
false.positives.B
```

```
## [1] 0
```

We could correct this using the Dunn-Sidak method:


```r
k<-45
new.threshold.DS<-1-((1-0.05)^(1/k))
new.threshold.DS
```

```
## [1] 0.001139202
```

```r
false.positives.DS<-sum(p.values<new.threshold.DS,na.rm=T)
false.positives.DS
```

```
## [1] 0
```

We could correct this using the randomization method. This requires simulating data under the null hypothesis to generate a null distribution of p-values.



```r
p.values.all<-c()
min.p.values.all<-c()
for (k in 1:1000)
  {
    data.null<-list()
    for (i in 1:10)
      {
        data.null[[i]]<-rnorm(10)  #Note the double brackets for a list
      }
    p.values.null<-matrix(ncol=10,nrow=10)
    for (i in 1:9)
      {
        for (j in (i+1):10)
          {
            p.values.null[i,j]<-t.test(data.null[[i]],data.null[[j]])$p.value 
          }
      }
    p.values.all<-c(p.values.all,c(p.values.null)[!is.na(c(p.values.null))])
    min.p.values.all<-c(min.p.values.all,min(c(p.values.null)[!is.na(c(p.values.null))]))
  }
new.threshold.R<-quantile(min.p.values.all,probs=c(0.05))
new.threshold.R
```

```
##          5% 
## 0.001447136
```

```r
false.positives.R<-sum(p.values<new.threshold.R,na.rm=T)
false.positives.R
```

```
## [1] 0
```

If you were to do this experiment (all of the code in the preceding clock) 100 times, you should get at least 1 false positive 5 times, since we have set the threshold such that we have a 5% chance that the smallest p-value in that set of 45 comparisons will be smaller than the threshold we set.
