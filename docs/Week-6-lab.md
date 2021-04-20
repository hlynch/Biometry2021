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
##       [,1]      [,2]      [,3]      [,4]       [,5]      [,6]      [,7]
##  [1,]   NA 0.9578285 0.6032199 0.9514301 0.18850240 0.4333334 0.8146620
##  [2,]   NA        NA 0.5981183 0.9122368 0.25061193 0.4858783 0.8784925
##  [3,]   NA        NA        NA 0.6344384 0.04742805 0.2175201 0.3817419
##  [4,]   NA        NA        NA        NA 0.15349451 0.3953626 0.7529183
##  [5,]   NA        NA        NA        NA         NA 0.7863462 0.2034886
##  [6,]   NA        NA        NA        NA         NA        NA 0.5064344
##  [7,]   NA        NA        NA        NA         NA        NA        NA
##  [8,]   NA        NA        NA        NA         NA        NA        NA
##  [9,]   NA        NA        NA        NA         NA        NA        NA
## [10,]   NA        NA        NA        NA         NA        NA        NA
##             [,8]      [,9]     [,10]
##  [1,] 0.50909055 0.3552575 0.4532335
##  [2,] 0.50632588 0.4292617 0.5490245
##  [3,] 0.79944311 0.1172269 0.1166935
##  [4,] 0.53215947 0.3060788 0.3841641
##  [5,] 0.06665482 0.6643981 0.3571149
##  [6,] 0.19943544 0.9529782 0.7371220
##  [7,] 0.35049870 0.4139370 0.5382094
##  [8,]         NA 0.1335253 0.1586116
##  [9,]         NA        NA 0.6952917
## [10,]         NA        NA        NA
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
## 0.001590422
```

```r
false.positives.R<-sum(p.values<new.threshold.R,na.rm=T)
false.positives.R
```

```
## [1] 0
```

If you were to do this experiment (all of the code in the preceding clock) 100 times, you should get at least 1 false positive 5 times, since we have set the threshold such that we have a 5% chance that the smallest p-value in that set of 45 comparisons will be smaller than the threshold we set.
