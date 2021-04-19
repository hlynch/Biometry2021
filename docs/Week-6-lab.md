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
##       [,1]       [,2]      [,3]      [,4]      [,5]       [,6]      [,7]
##  [1,]   NA 0.01808993 0.1872514 0.1524112 0.6449736 0.05731701 0.3311274
##  [2,]   NA         NA 0.8795648 0.9102362 0.2362569 0.89643623 0.2771366
##  [3,]   NA         NA        NA 0.8363663 0.4325390 0.81794312 0.5392502
##  [4,]   NA         NA        NA        NA 0.3437049 0.98897134 0.4240925
##  [5,]   NA         NA        NA        NA        NA 0.25625575 0.7749108
##  [6,]   NA         NA        NA        NA        NA         NA 0.3137991
##  [7,]   NA         NA        NA        NA        NA         NA        NA
##  [8,]   NA         NA        NA        NA        NA         NA        NA
##  [9,]   NA         NA        NA        NA        NA         NA        NA
## [10,]   NA         NA        NA        NA        NA         NA        NA
##             [,8]       [,9]     [,10]
##  [1,] 0.02898771 0.06191261 0.2964931
##  [2,] 0.64792324 0.88825531 0.1489846
##  [3,] 0.63947317 0.96464433 0.4504599
##  [4,] 0.83513154 0.84403679 0.3515453
##  [5,] 0.16676374 0.31961934 0.8380888
##  [6,] 0.78175915 0.81532002 0.2189324
##  [7,] 0.19539720 0.39776737 0.8915328
##  [8,]         NA 0.59626415 0.1237495
##  [9,]         NA         NA 0.2751828
## [10,]         NA         NA        NA
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
## 0.001611674
```

```r
false.positives.R<-sum(p.values<new.threshold.R,na.rm=T)
false.positives.R
```

```
## [1] 0
```

If you were to do this experiment (all of the code in the preceding clock) 100 times, you should get at least 1 false positive 5 times, since we have set the threshold such that we have a 5% chance that the smallest p-value in that set of 45 comparisons will be smaller than the threshold we set.
