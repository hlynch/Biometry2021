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
##       [,1]      [,2]      [,3]       [,4]       [,5]       [,6]      [,7]
##  [1,]   NA 0.5120896 0.4964608 0.15653492 0.49104626 0.57454348 0.8247777
##  [2,]   NA        NA 0.9979817 0.01987823 0.96677151 0.92407914 0.3553376
##  [3,]   NA        NA        NA 0.01535459 0.96239059 0.92118814 0.3367900
##  [4,]   NA        NA        NA         NA 0.01821685 0.02879134 0.2222320
##  [5,]   NA        NA        NA         NA         NA 0.89247172 0.3381493
##  [6,]   NA        NA        NA         NA         NA         NA 0.4114975
##  [7,]   NA        NA        NA         NA         NA         NA        NA
##  [8,]   NA        NA        NA         NA         NA         NA        NA
##  [9,]   NA        NA        NA         NA         NA         NA        NA
## [10,]   NA        NA        NA         NA         NA         NA        NA
##            [,8]       [,9]      [,10]
##  [1,] 0.8327511 0.28535031 0.18180084
##  [2,] 0.2742285 0.55523001 0.03709608
##  [3,] 0.2429980 0.53608650 0.03222973
##  [4,] 0.1332988 0.01301655 0.95549271
##  [5,] 0.2561190 0.57780876 0.03462578
##  [6,] 0.3418095 0.51534506 0.04805525
##  [7,] 0.9631421 0.19217286 0.24926990
##  [8,]        NA 0.14173048 0.17541173
##  [9,]        NA         NA 0.02125672
## [10,]        NA         NA         NA
```

Now we can see how many of these p.values are "significant". We know these are false positives, because all the data were generated from the same distribution.


```r
false.positives<-sum(p.values<0.05,na.rm=T)
false.positives
```

```
## [1] 10
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
## 0.001747646
```

```r
false.positives.R<-sum(p.values<new.threshold.R,na.rm=T)
false.positives.R
```

```
## [1] 0
```

If you were to do this experiment (all of the code in the preceding clock) 100 times, you should get at least 1 false positive 5 times, since we have set the threshold such that we have a 5% chance that the smallest p-value in that set of 45 comparisons will be smaller than the threshold we set.
