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
##       [,1]      [,2]      [,3]      [,4]      [,5]       [,6]      [,7]
##  [1,]   NA 0.7025602 0.8862246 0.6704849 0.7927529 0.07540632 0.3479086
##  [2,]   NA        NA 0.6372134 0.9855583 0.8848196 0.27702157 0.6130799
##  [3,]   NA        NA        NA 0.6077012 0.7110985 0.09123174 0.3262412
##  [4,]   NA        NA        NA        NA 0.8628933 0.25637489 0.6106836
##  [5,]   NA        NA        NA        NA        NA 0.15772662 0.4893171
##  [6,]   NA        NA        NA        NA        NA         NA 0.6313965
##  [7,]   NA        NA        NA        NA        NA         NA        NA
##  [8,]   NA        NA        NA        NA        NA         NA        NA
##  [9,]   NA        NA        NA        NA        NA         NA        NA
## [10,]   NA        NA        NA        NA        NA         NA        NA
##             [,8]       [,9]      [,10]
##  [1,] 0.53225028 0.13142307 0.89273757
##  [2,] 0.38128925 0.31938623 0.63385770
##  [3,] 0.65656836 0.13365950 0.98844962
##  [4,] 0.35082268 0.30660595 0.60221300
##  [5,] 0.41261985 0.21682404 0.70886772
##  [6,] 0.03906835 0.94759894 0.07653418
##  [7,] 0.17640713 0.63834977 0.31540212
##  [8,]         NA 0.06422354 0.63343043
##  [9,]         NA         NA 0.12256473
## [10,]         NA         NA         NA
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
## 0.002171137
```

```r
false.positives.R<-sum(p.values<new.threshold.R,na.rm=T)
false.positives.R
```

```
## [1] 0
```

If you were to do this experiment (all of the code in the preceding clock) 100 times, you should get at least 1 false positive 5 times, since we have set the threshold such that we have a 5% chance that the smallest p-value in that set of 45 comparisons will be smaller than the threshold we set.
