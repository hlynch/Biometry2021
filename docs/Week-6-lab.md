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
##       [,1]      [,2]      [,3]       [,4]       [,5]      [,6]      [,7]
##  [1,]   NA 0.1837722 0.6854507 0.58013437 0.14383184 0.7345513 0.5563851
##  [2,]   NA        NA 0.1110884 0.05204322 0.91594972 0.4460167 0.3844770
##  [3,]   NA        NA        NA 0.94825956 0.08714417 0.5105166 0.3408502
##  [4,]   NA        NA        NA         NA 0.03552543 0.4286779 0.2147988
##  [5,]   NA        NA        NA         NA         NA 0.3913706 0.3129061
##  [6,]   NA        NA        NA         NA         NA        NA 0.9097747
##  [7,]   NA        NA        NA         NA         NA        NA        NA
##  [8,]   NA        NA        NA         NA         NA        NA        NA
##  [9,]   NA        NA        NA         NA         NA        NA        NA
## [10,]   NA        NA        NA         NA         NA        NA        NA
##            [,8]      [,9]      [,10]
##  [1,] 0.5796148 0.6373534 0.70964933
##  [2,] 0.3962932 0.4598176 0.08754844
##  [3,] 0.3609609 0.4214938 0.93950060
##  [4,] 0.2412960 0.3263718 0.86825929
##  [5,] 0.3270711 0.3975927 0.06413793
##  [6,] 0.9197906 0.9324820 0.51934503
##  [7,] 0.9893381 0.9844508 0.31598631
##  [8,]        NA 0.9936311 0.34163171
##  [9,]        NA        NA 0.41830613
## [10,]        NA        NA         NA
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
## 0.001991964
```

```r
false.positives.R<-sum(p.values<new.threshold.R,na.rm=T)
false.positives.R
```

```
## [1] 0
```

If you were to do this experiment (all of the code in the preceding clock) 100 times, you should get at least 1 false positive 5 times, since we have set the threshold such that we have a 5% chance that the smallest p-value in that set of 45 comparisons will be smaller than the threshold we set.
