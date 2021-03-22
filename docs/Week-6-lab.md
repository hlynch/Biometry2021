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
##       [,1]     [,2]      [,3]      [,4]      [,5]      [,6]       [,7]
##  [1,]   NA 0.799899 0.9897612 0.1008905 0.9846603 0.1513081 0.09330358
##  [2,]   NA       NA 0.8129506 0.3277230 0.8177337 0.4464267 0.27685211
##  [3,]   NA       NA        NA 0.1210558 0.9950542 0.1811853 0.10835930
##  [4,]   NA       NA        NA        NA 0.1260748 0.7222243 0.83552887
##  [5,]   NA       NA        NA        NA        NA 0.1884382 0.11220516
##  [6,]   NA       NA        NA        NA        NA        NA 0.59031925
##  [7,]   NA       NA        NA        NA        NA        NA         NA
##  [8,]   NA       NA        NA        NA        NA        NA         NA
##  [9,]   NA       NA        NA        NA        NA        NA         NA
## [10,]   NA       NA        NA        NA        NA        NA         NA
##            [,8]       [,9]      [,10]
##  [1,] 0.8612507 0.98490584 0.08726046
##  [2,] 0.7321573 0.78134999 0.27665807
##  [3,] 0.8581036 0.97485444 0.10294149
##  [4,] 0.1993277 0.07502943 0.85040088
##  [5,] 0.8553911 0.96952529 0.10691667
##  [6,] 0.2728843 0.11011848 0.59583687
##  [7,] 0.1697395 0.07438744 0.98132650
##  [8,]        NA 0.86479897 0.16878218
##  [9,]        NA         NA 0.06765135
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
## 0.001207251
```

```r
false.positives.R<-sum(p.values<new.threshold.R,na.rm=T)
false.positives.R
```

```
## [1] 0
```

If you were to do this experiment (all of the code in the preceding clock) 100 times, you should get at least 1 false positive 5 times, since we have set the threshold such that we have a 5% chance that the smallest p-value in that set of 45 comparisons will be smaller than the threshold we set.
