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
##       [,1]      [,2]       [,3]       [,4]      [,5]       [,6]       [,7]
##  [1,]   NA 0.1790671 0.49389565 0.42565891 0.4165550 0.11078940 0.51240664
##  [2,]   NA        NA 0.05806477 0.06711235 0.5397060 0.97603898 0.08609956
##  [3,]   NA        NA         NA 0.81195239 0.1468869 0.02656330 0.92770095
##  [4,]   NA        NA         NA         NA 0.1532970 0.04286506 0.89736230
##  [5,]   NA        NA         NA         NA        NA 0.45556539 0.19397010
##  [6,]   NA        NA         NA         NA        NA         NA 0.05665428
##  [7,]   NA        NA         NA         NA        NA         NA         NA
##  [8,]   NA        NA         NA         NA        NA         NA         NA
##  [9,]   NA        NA         NA         NA        NA         NA         NA
## [10,]   NA        NA         NA         NA        NA         NA         NA
##            [,8]      [,9]     [,10]
##  [1,] 0.4517271 0.9452434 0.9913079
##  [2,] 0.5938037 0.1772074 0.2184898
##  [3,] 0.1852733 0.5607196 0.5304531
##  [4,] 0.1776722 0.4755379 0.4522200
##  [5,] 0.9781845 0.4015730 0.4681665
##  [6,] 0.5285635 0.1156126 0.1563528
##  [7,] 0.2200004 0.5648290 0.5358149
##  [8,]        NA 0.4336261 0.4927799
##  [9,]        NA        NA 0.9418747
## [10,]        NA        NA        NA
```

Now we can see how many of these p.values are "significant". We know these are false positives, because all the data were generated from the same distribution.


```r
false.positives<-sum(p.values<0.05,na.rm=T)
false.positives
```

```
## [1] 2
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
## 0.001836248
```

```r
false.positives.R<-sum(p.values<new.threshold.R,na.rm=T)
false.positives.R
```

```
## [1] 0
```

If you were to do this experiment (all of the code in the preceding clock) 100 times, you should get at least 1 false positive 5 times, since we have set the threshold such that we have a 5% chance that the smallest p-value in that set of 45 comparisons will be smaller than the threshold we set.
