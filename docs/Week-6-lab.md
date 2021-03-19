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
##       [,1]      [,2]      [,3]      [,4]      [,5]      [,6]       [,7]
##  [1,]   NA 0.1278338 0.5806716 0.7993356 0.6634569 0.5300475 0.03908005
##  [2,]   NA        NA 0.1803076 0.1504388 0.4017030 0.3896017 0.46153807
##  [3,]   NA        NA        NA 0.7532662 0.9705715 0.8409507 0.03622699
##  [4,]   NA        NA        NA        NA 0.8045205 0.6655547 0.03847441
##  [5,]   NA        NA        NA        NA        NA 0.9102868 0.18651684
##  [6,]   NA        NA        NA        NA        NA        NA 0.14764817
##  [7,]   NA        NA        NA        NA        NA        NA         NA
##  [8,]   NA        NA        NA        NA        NA        NA         NA
##  [9,]   NA        NA        NA        NA        NA        NA         NA
## [10,]   NA        NA        NA        NA        NA        NA         NA
##            [,8]      [,9]     [,10]
##  [1,] 0.3130959 0.3324465 0.1406420
##  [2,] 0.5181008 0.6703465 0.9983790
##  [3,] 0.5180898 0.5310915 0.2076044
##  [4,] 0.3971126 0.4173207 0.1688519
##  [5,] 0.6984804 0.6670025 0.4137721
##  [6,] 0.7513877 0.7134339 0.4072629
##  [7,] 0.1731734 0.3133822 0.4853066
##  [8,]        NA 0.9166628 0.5404521
##  [9,]        NA        NA 0.6820325
## [10,]        NA        NA        NA
```

Now we can see how many of these p.values are "significant". We know these are false positives, because all the data were generated from the same distribution.


```r
false.positives<-sum(p.values<0.05,na.rm=T)
false.positives
```

```
## [1] 3
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
## 0.001931768
```

```r
false.positives.R<-sum(p.values<new.threshold.R,na.rm=T)
false.positives.R
```

```
## [1] 0
```

If you were to do this experiment (all of the code in the preceding clock) 100 times, you should get at least 1 false positive 5 times, since we have set the threshold such that we have a 5% chance that the smallest p-value in that set of 45 comparisons will be smaller than the threshold we set.
