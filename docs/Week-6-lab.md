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
##       [,1]      [,2]      [,3]       [,4]      [,5]      [,6]      [,7]
##  [1,]   NA 0.3409058 0.6445162 0.19253712 0.8777781 0.7093718 0.9032732
##  [2,]   NA        NA 0.6662958 0.07649484 0.3681016 0.6643270 0.4648202
##  [3,]   NA        NA        NA 0.14830469 0.6201594 0.9730745 0.7675865
##  [4,]   NA        NA        NA         NA 0.3166721 0.1824301 0.2138004
##  [5,]   NA        NA        NA         NA        NA 0.6690177 0.8210520
##  [6,]   NA        NA        NA         NA        NA        NA 0.8124014
##  [7,]   NA        NA        NA         NA        NA        NA        NA
##  [8,]   NA        NA        NA         NA        NA        NA        NA
##  [9,]   NA        NA        NA         NA        NA        NA        NA
## [10,]   NA        NA        NA         NA        NA        NA        NA
##            [,8]       [,9]       [,10]
##  [1,] 0.6910762 0.03026196 0.051262891
##  [2,] 0.5947222 0.01862010 0.418668140
##  [3,] 0.9342165 0.04285557 0.206281269
##  [4,] 0.1514937 0.79552335 0.014093449
##  [5,] 0.6595325 0.13489072 0.088079928
##  [6,] 0.9676164 0.06955537 0.229221956
##  [7,] 0.8206926 0.06559597 0.111218852
##  [8,]        NA 0.03739807 0.158746215
##  [9,]        NA         NA 0.001549432
## [10,]        NA         NA          NA
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
## 0.001604466
```

```r
false.positives.R<-sum(p.values<new.threshold.R,na.rm=T)
false.positives.R
```

```
## [1] 1
```

If you were to do this experiment (all of the code in the preceding clock) 100 times, you should get at least 1 false positive 5 times, since we have set the threshold such that we have a 5% chance that the smallest p-value in that set of 45 comparisons will be smaller than the threshold we set.
