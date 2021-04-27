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
##       [,1]       [,2]       [,3]       [,4]      [,5]       [,6]       [,7]
##  [1,]   NA 0.05727256 0.67062463 0.71792931 0.2757237 0.96380027 0.88902685
##  [2,]   NA         NA 0.02672971 0.03758002 0.4790781 0.09935053 0.02846682
##  [3,]   NA         NA         NA 0.96555735 0.1513888 0.67241379 0.74118823
##  [4,]   NA         NA         NA         NA 0.1797566 0.71327223 0.79176397
##  [5,]   NA         NA         NA         NA        NA 0.34533446 0.19686063
##  [6,]   NA         NA         NA         NA        NA         NA 0.86556367
##  [7,]   NA         NA         NA         NA        NA         NA         NA
##  [8,]   NA         NA         NA         NA        NA         NA         NA
##  [9,]   NA         NA         NA         NA        NA         NA         NA
## [10,]   NA         NA         NA         NA        NA         NA         NA
##            [,8]      [,9]      [,10]
##  [1,] 0.8116719 0.3316388 0.81260009
##  [2,] 0.2455220 0.4198211 0.02783865
##  [3,] 0.5863511 0.1874789 0.82733272
##  [4,] 0.6185641 0.2185703 0.87246965
##  [5,] 0.5487468 0.9183341 0.18057667
##  [6,] 0.8503919 0.4020681 0.79929955
##  [7,] 0.7304301 0.2458334 0.90851481
##  [8,]        NA 0.6074725 0.68177789
##  [9,]        NA        NA 0.22498171
## [10,]        NA        NA         NA
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
## 0.001925717
```

```r
false.positives.R<-sum(p.values<new.threshold.R,na.rm=T)
false.positives.R
```

```
## [1] 0
```

If you were to do this experiment (all of the code in the preceding clock) 100 times, you should get at least 1 false positive 5 times, since we have set the threshold such that we have a 5% chance that the smallest p-value in that set of 45 comparisons will be smaller than the threshold we set.
