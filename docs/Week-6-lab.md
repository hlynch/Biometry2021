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
##  [1,]   NA 0.6112596 0.2207574 0.36457865 0.37206097 0.6603578 0.5210114
##  [2,]   NA        NA 0.1878838 0.92762857 0.91903446 0.4366352 0.8950418
##  [3,]   NA        NA        NA 0.05401178 0.05804099 0.4970019 0.0746872
##  [4,]   NA        NA        NA         NA 0.98562075 0.2396548 0.7177566
##  [5,]   NA        NA        NA         NA         NA 0.2448061 0.7144024
##  [6,]   NA        NA        NA         NA         NA        NA 0.3319239
##  [7,]   NA        NA        NA         NA         NA        NA        NA
##  [8,]   NA        NA        NA         NA         NA        NA        NA
##  [9,]   NA        NA        NA         NA         NA        NA        NA
## [10,]   NA        NA        NA         NA         NA        NA        NA
##             [,8]      [,9]      [,10]
##  [1,] 0.32142080 0.3761821 0.31765177
##  [2,] 0.90661844 0.9310647 0.79934642
##  [3,] 0.04367422 0.0574847 0.05713087
##  [4,] 0.96951236 0.9956749 0.81410376
##  [5,] 0.98573568 0.9816791 0.83068139
##  [6,] 0.21342953 0.2470573 0.21271933
##  [7,] 0.66926747 0.7280724 0.58375069
##  [8,]         NA 0.9655790 0.83241606
##  [9,]         NA        NA 0.81284159
## [10,]         NA        NA         NA
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
## 0.001489526
```

```r
false.positives.R<-sum(p.values<new.threshold.R,na.rm=T)
false.positives.R
```

```
## [1] 0
```

If you were to do this experiment (all of the code in the preceding clock) 100 times, you should get at least 1 false positive 5 times, since we have set the threshold such that we have a 5% chance that the smallest p-value in that set of 45 comparisons will be smaller than the threshold we set.
