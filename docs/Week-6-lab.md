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
##       [,1]     [,2]      [,3]       [,4]      [,5]      [,6]      [,7]
##  [1,]   NA 0.696921 0.8781814 0.04189029 0.1432122 0.3887090 0.1172191
##  [2,]   NA       NA 0.8326266 0.16001412 0.5293049 0.7044614 0.4007272
##  [3,]   NA       NA        NA 0.09626254 0.3481246 0.5410667 0.2600681
##  [4,]   NA       NA        NA         NA 0.2546678 0.2847480 0.4345549
##  [5,]   NA       NA        NA         NA        NA 0.8610593 0.7125811
##  [6,]   NA       NA        NA         NA        NA        NA 0.6642354
##  [7,]   NA       NA        NA         NA        NA        NA        NA
##  [8,]   NA       NA        NA         NA        NA        NA        NA
##  [9,]   NA       NA        NA         NA        NA        NA        NA
## [10,]   NA       NA        NA         NA        NA        NA        NA
##            [,8]        [,9]      [,10]
##  [1,] 0.3147552 0.004496465 0.01260066
##  [2,] 0.7367924 0.112550863 0.08204760
##  [3,] 0.5364939 0.052139933 0.04303644
##  [4,] 0.1707156 0.866521489 0.76633395
##  [5,] 0.6859272 0.112692261 0.11171996
##  [6,] 0.9073044 0.231077491 0.15898856
##  [7,] 0.4901517 0.364098573 0.24273261
##  [8,]        NA 0.064819713 0.07021395
##  [9,]        NA          NA 0.55905635
## [10,]        NA          NA         NA
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
## 0.001696977
```

```r
false.positives.R<-sum(p.values<new.threshold.R,na.rm=T)
false.positives.R
```

```
## [1] 0
```

If you were to do this experiment (all of the code in the preceding clock) 100 times, you should get at least 1 false positive 5 times, since we have set the threshold such that we have a 5% chance that the smallest p-value in that set of 45 comparisons will be smaller than the threshold we set.
