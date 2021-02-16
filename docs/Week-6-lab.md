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
##       [,1]      [,2]      [,3]      [,4]      [,5]      [,6]      [,7]
##  [1,]   NA 0.8956663 0.4245225 0.7641054 0.9082289 0.3669995 0.8352039
##  [2,]   NA        NA 0.3324649 0.8473892 0.9864740 0.2542594 0.9188663
##  [3,]   NA        NA        NA 0.2712251 0.3410228 0.9702608 0.3383040
##  [4,]   NA        NA        NA        NA 0.8355912 0.1971998 0.9487875
##  [5,]   NA        NA        NA        NA        NA 0.2641046 0.9079003
##  [6,]   NA        NA        NA        NA        NA        NA 0.2820048
##  [7,]   NA        NA        NA        NA        NA        NA        NA
##  [8,]   NA        NA        NA        NA        NA        NA        NA
##  [9,]   NA        NA        NA        NA        NA        NA        NA
## [10,]   NA        NA        NA        NA        NA        NA        NA
##            [,8]       [,9]      [,10]
##  [1,] 0.9198139 0.38704864 0.28894340
##  [2,] 0.8180883 0.41210935 0.30476229
##  [3,] 0.5045566 0.12775364 0.09910650
##  [4,] 0.7005925 0.51874338 0.38300215
##  [5,] 0.8298160 0.40667200 0.30099583
##  [6,] 0.4632491 0.08457639 0.06976178
##  [7,] 0.7698320 0.52308595 0.39246616
##  [8,]        NA 0.36386579 0.27361748
##  [9,]        NA         NA 0.77968478
## [10,]        NA         NA         NA
```

Now we can see how many of these p.values are "significant". We know these are false positives, because all the data were generated from the same distribution.


```r
false.positives<-sum(p.values<0.05,na.rm=T)
false.positives
```

```
## [1] 0
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
## 0.001587962
```

```r
false.positives.R<-sum(p.values<new.threshold.R,na.rm=T)
false.positives.R
```

```
## [1] 0
```

If you were to do this experiment (all of the code in the preceding clock) 100 times, you should get at least 1 false positive 5 times, since we have set the threshold such that we have a 5% chance that the smallest p-value in that set of 45 comparisons will be smaller than the threshold we set.
