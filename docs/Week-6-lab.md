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
##       [,1]       [,2]      [,3]      [,4]      [,5]      [,6]      [,7]
##  [1,]   NA 0.09657283 0.1871760 0.5205089 0.4636804 0.7294002 0.4350490
##  [2,]   NA         NA 0.9799941 0.4659550 0.5485559 0.1184849 0.4052536
##  [3,]   NA         NA        NA 0.5382010 0.6122307 0.1750630 0.4993235
##  [4,]   NA         NA        NA        NA 0.9219899 0.4263435 0.9958554
##  [5,]   NA         NA        NA        NA        NA 0.3845472 0.9174481
##  [6,]   NA         NA        NA        NA        NA        NA 0.3730150
##  [7,]   NA         NA        NA        NA        NA        NA        NA
##  [8,]   NA         NA        NA        NA        NA        NA        NA
##  [9,]   NA         NA        NA        NA        NA        NA        NA
## [10,]   NA         NA        NA        NA        NA        NA        NA
##            [,8]       [,9]     [,10]
##  [1,] 0.4266972 0.07870695 0.6601776
##  [2,] 0.4379227 0.98180925 0.3617543
##  [3,] 0.5265420 0.96436184 0.4376708
##  [4,] 0.9700661 0.44358916 0.8662612
##  [5,] 0.9438901 0.52711115 0.7934355
##  [6,] 0.3652769 0.10676784 0.5268119
##  [7,] 0.9702895 0.37755750 0.8465794
##  [8,]        NA 0.41107209 0.8241955
##  [9,]        NA         NA 0.3405249
## [10,]        NA         NA        NA
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
## 0.001723585
```

```r
false.positives.R<-sum(p.values<new.threshold.R,na.rm=T)
false.positives.R
```

```
## [1] 0
```

If you were to do this experiment (all of the code in the preceding clock) 100 times, you should get at least 1 false positive 5 times, since we have set the threshold such that we have a 5% chance that the smallest p-value in that set of 45 comparisons will be smaller than the threshold we set.
