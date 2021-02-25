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
##       [,1]        [,2]       [,3]        [,4]         [,5]       [,6]
##  [1,]   NA 0.004037271 0.05364278 0.004178771 0.0002485922 0.04415976
##  [2,]   NA          NA 0.16137766 0.848029326 0.4122126645 0.42653670
##  [3,]   NA          NA         NA 0.195192326 0.0164320661 0.65992010
##  [4,]   NA          NA         NA          NA 0.2728652672 0.50853988
##  [5,]   NA          NA         NA          NA           NA 0.11889685
##  [6,]   NA          NA         NA          NA           NA         NA
##  [7,]   NA          NA         NA          NA           NA         NA
##  [8,]   NA          NA         NA          NA           NA         NA
##  [9,]   NA          NA         NA          NA           NA         NA
## [10,]   NA          NA         NA          NA           NA         NA
##             [,7]        [,8]       [,9]       [,10]
##  [1,] 0.01562405 0.006894504 0.14196253 0.370499195
##  [2,] 0.57947420 0.987360631 0.11457168 0.044615865
##  [3,] 0.41792633 0.199564672 0.72596984 0.361438138
##  [4,] 0.69146230 0.872013266 0.13755716 0.052484118
##  [5,] 0.16101166 0.442484122 0.01565793 0.005628277
##  [6,] 0.77618961 0.460395395 0.48241081 0.241833781
##  [7,]         NA 0.614568626 0.29329927 0.127626285
##  [8,]         NA          NA 0.14133161 0.058686179
##  [9,]         NA          NA         NA 0.591337180
## [10,]         NA          NA         NA          NA
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
## [1] 1
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
## [1] 1
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
## 0.00155177
```

```r
false.positives.R<-sum(p.values<new.threshold.R,na.rm=T)
false.positives.R
```

```
## [1] 1
```

If you were to do this experiment (all of the code in the preceding clock) 100 times, you should get at least 1 false positive 5 times, since we have set the threshold such that we have a 5% chance that the smallest p-value in that set of 45 comparisons will be smaller than the threshold we set.
