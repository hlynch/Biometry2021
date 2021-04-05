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
##       [,1]      [,2]       [,3]      [,4]      [,5]      [,6]      [,7]
##  [1,]   NA 0.3222741 0.71360490 0.5058128 0.9945450 0.8421853 0.6221801
##  [2,]   NA        NA 0.04240084 0.6332494 0.1181547 0.4060628 0.4466807
##  [3,]   NA        NA         NA 0.1094230 0.4919607 0.5040487 0.1592129
##  [4,]   NA        NA         NA        NA 0.2772699 0.6355872 0.7841793
##  [5,]   NA        NA         NA        NA        NA 0.7816359 0.3980127
##  [6,]   NA        NA         NA        NA        NA        NA 0.7786244
##  [7,]   NA        NA         NA        NA        NA        NA        NA
##  [8,]   NA        NA         NA        NA        NA        NA        NA
##  [9,]   NA        NA         NA        NA        NA        NA        NA
## [10,]   NA        NA         NA        NA        NA        NA        NA
##            [,8]       [,9]     [,10]
##  [1,] 0.5356557 0.38041217 0.8021941
##  [2,] 0.5157948 0.84684759 0.3473179
##  [3,] 0.0814106 0.05072868 0.3848797
##  [4,] 0.9017005 0.76278240 0.6021204
##  [5,] 0.2569409 0.14809297 0.6978134
##  [6,] 0.6782215 0.48074344 0.9736952
##  [7,] 0.8567384 0.55214179 0.7681860
##  [8,]        NA 0.64080500 0.6439360
##  [9,]        NA         NA 0.4240882
## [10,]        NA         NA        NA
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
## 0.001662329
```

```r
false.positives.R<-sum(p.values<new.threshold.R,na.rm=T)
false.positives.R
```

```
## [1] 0
```

If you were to do this experiment (all of the code in the preceding clock) 100 times, you should get at least 1 false positive 5 times, since we have set the threshold such that we have a 5% chance that the smallest p-value in that set of 45 comparisons will be smaller than the threshold we set.
