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
##       [,1]     [,2]      [,3]      [,4]      [,5]        [,6]       [,7]
##  [1,]   NA 0.766409 0.5035172 0.2723895 0.3279641 0.047989290 0.63474110
##  [2,]   NA       NA 0.5748155 0.2318004 0.3407202 0.006842959 0.77110354
##  [3,]   NA       NA        NA 0.6027549 0.6662664 0.085373082 0.84063500
##  [4,]   NA       NA        NA        NA 0.9916178 0.200186017 0.48812848
##  [5,]   NA       NA        NA        NA        NA 0.303770346 0.55490019
##  [6,]   NA       NA        NA        NA        NA          NA 0.07446532
##  [7,]   NA       NA        NA        NA        NA          NA         NA
##  [8,]   NA       NA        NA        NA        NA          NA         NA
##  [9,]   NA       NA        NA        NA        NA          NA         NA
## [10,]   NA       NA        NA        NA        NA          NA         NA
##            [,8]       [,9]      [,10]
##  [1,] 0.2498541 0.68403638 0.71602365
##  [2,] 0.2222855 0.84635406 0.90739645
##  [3,] 0.5401417 0.75802407 0.66549082
##  [4,] 0.8896582 0.40874033 0.30727468
##  [5,] 0.8978591 0.48596553 0.40568932
##  [6,] 0.3177327 0.04788181 0.01725296
##  [7,] 0.4408177 0.92573608 0.85220092
##  [8,]        NA 0.37107949 0.28496924
##  [9,]        NA         NA 0.92962572
## [10,]        NA         NA         NA
```

Now we can see how many of these p.values are "significant". We know these are false positives, because all the data were generated from the same distribution.


```r
false.positives<-sum(p.values<0.05,na.rm=T)
false.positives
```

```
## [1] 4
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
## 0.001766554
```

```r
false.positives.R<-sum(p.values<new.threshold.R,na.rm=T)
false.positives.R
```

```
## [1] 0
```

If you were to do this experiment (all of the code in the preceding clock) 100 times, you should get at least 1 false positive 5 times, since we have set the threshold such that we have a 5% chance that the smallest p-value in that set of 45 comparisons will be smaller than the threshold we set.
