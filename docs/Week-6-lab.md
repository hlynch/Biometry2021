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
##       [,1]      [,2]      [,3]      [,4]       [,5]      [,6]       [,7]
##  [1,]   NA 0.9664451 0.3781818 0.2364026 0.15083488 0.2071160 0.13520275
##  [2,]   NA        NA 0.2877851 0.1493874 0.06385648 0.1402809 0.06241044
##  [3,]   NA        NA        NA 0.7754938 0.64229684 0.6528666 0.56156784
##  [4,]   NA        NA        NA        NA 0.88046919 0.8389004 0.76558687
##  [5,]   NA        NA        NA        NA         NA 0.9182147 0.84279776
##  [6,]   NA        NA        NA        NA         NA        NA 0.96598104
##  [7,]   NA        NA        NA        NA         NA        NA         NA
##  [8,]   NA        NA        NA        NA         NA        NA         NA
##  [9,]   NA        NA        NA        NA         NA        NA         NA
## [10,]   NA        NA        NA        NA         NA        NA         NA
##             [,8]       [,9]     [,10]
##  [1,] 0.16524858 0.05308328 0.3571013
##  [2,] 0.09726411 0.02013563 0.2941001
##  [3,] 0.59137716 0.25507701 0.8574625
##  [4,] 0.78067637 0.35635805 0.9620380
##  [5,] 0.85513836 0.34079435 0.8740788
##  [6,] 0.95842577 0.53745546 0.8355642
##  [7,] 0.98705221 0.47392775 0.7876735
##  [8,]         NA 0.53968432 0.7925061
##  [9,]         NA         NA 0.4546297
## [10,]         NA         NA        NA
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
## 0.001245255
```

```r
false.positives.R<-sum(p.values<new.threshold.R,na.rm=T)
false.positives.R
```

```
## [1] 0
```

If you were to do this experiment (all of the code in the preceding clock) 100 times, you should get at least 1 false positive 5 times, since we have set the threshold such that we have a 5% chance that the smallest p-value in that set of 45 comparisons will be smaller than the threshold we set.
