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
##       [,1]      [,2]      [,3]       [,4]       [,5]      [,6]       [,7]
##  [1,]   NA 0.1416562 0.7240392 0.41288398 0.44510241 0.9005428 0.18783667
##  [2,]   NA        NA 0.3687797 0.08640431 0.04723061 0.2220358 0.70460908
##  [3,]   NA        NA        NA 0.31455189 0.32618008 0.6912887 0.32497486
##  [4,]   NA        NA        NA         NA 0.80996494 0.5330303 0.09152028
##  [5,]   NA        NA        NA         NA         NA 0.6219036 0.08136982
##  [6,]   NA        NA        NA         NA         NA        NA 0.21222318
##  [7,]   NA        NA        NA         NA         NA        NA         NA
##  [8,]   NA        NA        NA         NA         NA        NA         NA
##  [9,]   NA        NA        NA         NA         NA        NA         NA
## [10,]   NA        NA        NA         NA         NA        NA         NA
##            [,8]      [,9]     [,10]
##  [1,] 0.6062913 0.5632074 0.8978224
##  [2,] 0.1540624 0.3539545 0.1417644
##  [3,] 0.4684438 0.8949575 0.6581157
##  [4,] 0.8111920 0.2382567 0.4809733
##  [5,] 0.9654899 0.2092335 0.5467816
##  [6,] 0.7197923 0.5774390 0.9848899
##  [7,] 0.1467636 0.3340517 0.1750514
##  [8,]        NA 0.3798152 0.6804934
##  [9,]        NA        NA 0.5106704
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
##         5% 
## 0.00216642
```

```r
false.positives.R<-sum(p.values<new.threshold.R,na.rm=T)
false.positives.R
```

```
## [1] 0
```

If you were to do this experiment (all of the code in the preceding clock) 100 times, you should get at least 1 false positive 5 times, since we have set the threshold such that we have a 5% chance that the smallest p-value in that set of 45 comparisons will be smaller than the threshold we set.
