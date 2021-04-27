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
##  [1,]   NA 0.370612 0.9470316 0.2588172 0.2461497 0.1371328 0.93595698
##  [2,]   NA       NA 0.3597884 0.8356682 0.9477449 0.6091322 0.21748958
##  [3,]   NA       NA        NA 0.2378272 0.2044340 0.1078748 0.84898851
##  [4,]   NA       NA        NA        NA 0.8445556 0.7620829 0.11545264
##  [5,]   NA       NA        NA        NA        NA 0.5328797 0.03673583
##  [6,]   NA       NA        NA        NA        NA        NA 0.02675870
##  [7,]   NA       NA        NA        NA        NA        NA         NA
##  [8,]   NA       NA        NA        NA        NA        NA         NA
##  [9,]   NA       NA        NA        NA        NA        NA         NA
## [10,]   NA       NA        NA        NA        NA        NA         NA
##             [,8]       [,9]      [,10]
##  [1,] 0.70232200 0.90207316 0.59221351
##  [2,] 0.16232648 0.27300300 0.12699305
##  [3,] 0.60680679 0.83223883 0.49377399
##  [4,] 0.08986505 0.17234165 0.06851527
##  [5,] 0.04829295 0.13563097 0.03521532
##  [6,] 0.02706482 0.07179311 0.02000518
##  [7,] 0.63383801 0.93650301 0.48598570
##  [8,]         NA 0.77483461 0.84318812
##  [9,]         NA         NA 0.64484157
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
## 0.001381036
```

```r
false.positives.R<-sum(p.values<new.threshold.R,na.rm=T)
false.positives.R
```

```
## [1] 0
```

If you were to do this experiment (all of the code in the preceding clock) 100 times, you should get at least 1 false positive 5 times, since we have set the threshold such that we have a 5% chance that the smallest p-value in that set of 45 comparisons will be smaller than the threshold we set.
