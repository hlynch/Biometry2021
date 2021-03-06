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
##       [,1]       [,2]       [,3]        [,4]      [,5]      [,6]       [,7]
##  [1,]   NA 0.02318998 0.03921772 0.008002255 0.1641074 0.2144428 0.48174699
##  [2,]   NA         NA 0.82045512 0.650537367 0.4997475 0.5250178 0.10396891
##  [3,]   NA         NA         NA 0.908063196 0.4446041 0.4631071 0.12634121
##  [4,]   NA         NA         NA          NA 0.2943971 0.3337476 0.04179301
##  [5,]   NA         NA         NA          NA        NA 0.9787854 0.43496043
##  [6,]   NA         NA         NA          NA        NA        NA 0.49522891
##  [7,]   NA         NA         NA          NA        NA        NA         NA
##  [8,]   NA         NA         NA          NA        NA        NA         NA
##  [9,]   NA         NA         NA          NA        NA        NA         NA
## [10,]   NA         NA         NA          NA        NA        NA         NA
##             [,8]       [,9]      [,10]
##  [1,] 0.51853183 0.07978914 0.38585711
##  [2,] 0.10884960 0.93089390 0.05635670
##  [3,] 0.12748267 0.92121384 0.09752250
##  [4,] 0.04622628 0.82987934 0.01494479
##  [5,] 0.42779766 0.55552645 0.39548765
##  [6,] 0.48564147 0.56372712 0.47017346
##  [7,] 0.97098694 0.20652471 0.96639679
##  [8,]         NA 0.20483784 0.93415338
##  [9,]         NA         NA 0.18140946
## [10,]         NA         NA         NA
```

Now we can see how many of these p.values are "significant". We know these are false positives, because all the data were generated from the same distribution.


```r
false.positives<-sum(p.values<0.05,na.rm=T)
false.positives
```

```
## [1] 6
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
## 0.001776443
```

```r
false.positives.R<-sum(p.values<new.threshold.R,na.rm=T)
false.positives.R
```

```
## [1] 0
```

If you were to do this experiment (all of the code in the preceding clock) 100 times, you should get at least 1 false positive 5 times, since we have set the threshold such that we have a 5% chance that the smallest p-value in that set of 45 comparisons will be smaller than the threshold we set.
