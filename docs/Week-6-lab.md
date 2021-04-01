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
##       [,1]      [,2]       [,3]       [,4]      [,5]      [,6]       [,7]
##  [1,]   NA 0.0581118 0.09788518 0.90514026 0.2998377 0.5088893 0.01941372
##  [2,]   NA        NA 0.71322213 0.06173588 0.3391861 0.4284157 0.43339009
##  [3,]   NA        NA         NA 0.10650367 0.5173788 0.5636554 0.26683383
##  [4,]   NA        NA         NA         NA 0.3381724 0.5601843 0.01997926
##  [5,]   NA        NA         NA         NA        NA 0.8965854 0.11922369
##  [6,]   NA        NA         NA         NA        NA        NA 0.20964765
##  [7,]   NA        NA         NA         NA        NA        NA         NA
##  [8,]   NA        NA         NA         NA        NA        NA         NA
##  [9,]   NA        NA         NA         NA        NA        NA         NA
## [10,]   NA        NA         NA         NA        NA        NA         NA
##             [,8]       [,9]      [,10]
##  [1,] 0.92474510 0.36869921 0.84641486
##  [2,] 0.05946597 0.24140516 0.06865379
##  [3,] 0.09709067 0.38858915 0.11870985
##  [4,] 0.83269516 0.41682777 0.93913249
##  [5,] 0.27888086 0.85708351 0.37180730
##  [6,] 0.47256191 0.99724678 0.59667931
##  [7,] 0.02072630 0.07944002 0.02205608
##  [8,]         NA 0.34088106 0.77718813
##  [9,]         NA         NA 0.45722745
## [10,]         NA         NA         NA
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
##         5% 
## 0.00171418
```

```r
false.positives.R<-sum(p.values<new.threshold.R,na.rm=T)
false.positives.R
```

```
## [1] 0
```

If you were to do this experiment (all of the code in the preceding clock) 100 times, you should get at least 1 false positive 5 times, since we have set the threshold such that we have a 5% chance that the smallest p-value in that set of 45 comparisons will be smaller than the threshold we set.
