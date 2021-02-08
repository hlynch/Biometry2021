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
##       [,1]     [,2]      [,3]      [,4]      [,5]      [,6]      [,7]      [,8]
##  [1,]   NA 0.386674 0.8731121 0.6265219 0.7037056 0.4290053 0.7824417 0.8647718
##  [2,]   NA       NA 0.2770420 0.6293179 0.1967478 0.8568907 0.2263649 0.5529290
##  [3,]   NA       NA        NA 0.4821811 0.8129037 0.2995912 0.9020609 0.7487978
##  [4,]   NA       NA        NA        NA 0.3504005 0.7206063 0.4032312 0.8097944
##  [5,]   NA       NA        NA        NA        NA 0.2081556 0.9068244 0.6059583
##  [6,]   NA       NA        NA        NA        NA        NA 0.2407097 0.6176169
##  [7,]   NA       NA        NA        NA        NA        NA        NA 0.6712492
##  [8,]   NA       NA        NA        NA        NA        NA        NA        NA
##  [9,]   NA       NA        NA        NA        NA        NA        NA        NA
## [10,]   NA       NA        NA        NA        NA        NA        NA        NA
##             [,9]      [,10]
##  [1,] 0.41352781 0.47999805
##  [2,] 0.05642677 0.09675877
##  [3,] 0.49024391 0.56269858
##  [4,] 0.11597644 0.18490667
##  [5,] 0.68616887 0.74371768
##  [6,] 0.04497730 0.09455380
##  [7,] 0.58014374 0.64813313
##  [8,] 0.36596538 0.41901739
##  [9,]         NA 0.96720154
## [10,]         NA         NA
```

Now we can see how many of these p.values are "significant". We know these are false positives, because all the data were generated from the same distribution.


```r
false.positives<-sum(p.values<0.05,na.rm=T)
false.positives
```

```
## [1] 1
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
## 0.001927143
```

```r
false.positives.R<-sum(p.values<new.threshold.R,na.rm=T)
false.positives.R
```

```
## [1] 0
```

If you were to do this experiment (all of the code in the preceding clock) 100 times, you should get at least 1 false positive 5 times, since we have set the threshold such that we have a 5% chance that the smallest p-value in that set of 45 comparisons will be smaller than the threshold we set.
