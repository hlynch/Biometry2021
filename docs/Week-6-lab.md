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
##       [,1]      [,2]      [,3]      [,4]      [,5]      [,6]      [,7]
##  [1,]   NA 0.7616473 0.5189500 0.7217302 0.3740912 0.7224755 0.3039820
##  [2,]   NA        NA 0.7946591 0.5784191 0.6668577 0.9935532 0.2646924
##  [3,]   NA        NA        NA 0.3860927 0.8654872 0.7594554 0.1659147
##  [4,]   NA        NA        NA        NA 0.2794914 0.5281269 0.4894357
##  [5,]   NA        NA        NA        NA        NA 0.6108020 0.1177328
##  [6,]   NA        NA        NA        NA        NA        NA 0.2235116
##  [7,]   NA        NA        NA        NA        NA        NA        NA
##  [8,]   NA        NA        NA        NA        NA        NA        NA
##  [9,]   NA        NA        NA        NA        NA        NA        NA
## [10,]   NA        NA        NA        NA        NA        NA        NA
##            [,8]      [,9]     [,10]
##  [1,] 0.9233667 0.5829747 0.6581017
##  [2,] 0.8671960 0.4784031 0.5346781
##  [3,] 0.6675778 0.3078735 0.3607688
##  [4,] 0.7148117 0.8578583 0.9177490
##  [5,] 0.5473965 0.2171163 0.2659107
##  [6,] 0.8576345 0.4223951 0.4878627
##  [7,] 0.3427947 0.5939479 0.5661035
##  [8,]        NA 0.6033211 0.6598912
##  [9,]        NA        NA 0.9467434
## [10,]        NA        NA        NA
```

Now we can see how many of these p.values are "significant". We know these are false positives, because all the data were generated from the same distribution.


```r
false.positives<-sum(p.values<0.05,na.rm=T)
false.positives
```

```
## [1] 0
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
## 0.001937069
```

```r
false.positives.R<-sum(p.values<new.threshold.R,na.rm=T)
false.positives.R
```

```
## [1] 0
```

If you were to do this experiment (all of the code in the preceding clock) 100 times, you should get at least 1 false positive 5 times, since we have set the threshold such that we have a 5% chance that the smallest p-value in that set of 45 comparisons will be smaller than the threshold we set.
