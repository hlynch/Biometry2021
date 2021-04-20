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
##       [,1]       [,2]      [,3]      [,4]       [,5]      [,6]       [,7]
##  [1,]   NA 0.04186525 0.1869201 0.8419615 0.70022178 0.1807904 0.90863566
##  [2,]   NA         NA 0.5138496 0.1122906 0.08013699 0.3281831 0.06653403
##  [3,]   NA         NA        NA 0.3290658 0.31153291 0.8363496 0.24690475
##  [4,]   NA         NA        NA        NA 0.90310446 0.3617358 0.92671852
##  [5,]   NA         NA        NA        NA         NA 0.3286601 0.80692945
##  [6,]   NA         NA        NA        NA         NA        NA 0.25721932
##  [7,]   NA         NA        NA        NA         NA        NA         NA
##  [8,]   NA         NA        NA        NA         NA        NA         NA
##  [9,]   NA         NA        NA        NA         NA        NA         NA
## [10,]   NA         NA        NA        NA         NA        NA         NA
##            [,8]      [,9]     [,10]
##  [1,] 0.3802870 0.4366381 0.3085945
##  [2,] 0.2268334 0.3899732 0.4494167
##  [3,] 0.6137700 0.7550886 0.8689438
##  [4,] 0.5756065 0.5825422 0.4557632
##  [5,] 0.5965304 0.6126718 0.4646608
##  [6,] 0.7108376 0.8603011 0.9985171
##  [7,] 0.4716655 0.5055257 0.3748889
##  [8,]        NA 0.9224314 0.7770814
##  [9,]        NA        NA 0.8821982
## [10,]        NA        NA        NA
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
## 0.001630565
```

```r
false.positives.R<-sum(p.values<new.threshold.R,na.rm=T)
false.positives.R
```

```
## [1] 0
```

If you were to do this experiment (all of the code in the preceding clock) 100 times, you should get at least 1 false positive 5 times, since we have set the threshold such that we have a 5% chance that the smallest p-value in that set of 45 comparisons will be smaller than the threshold we set.
