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
##       [,1]      [,2]        [,3]      [,4]       [,5]       [,6]       [,7]
##  [1,]   NA 0.5499847 0.003256125 0.3182146 0.24115792 0.51365389 0.07710116
##  [2,]   NA        NA 0.043391469 0.6053387 0.75091681 0.97456811 0.27276864
##  [3,]   NA        NA          NA 0.2827111 0.03162916 0.04141813 0.41188821
##  [4,]   NA        NA          NA        NA 0.73737141 0.61818421 0.69238909
##  [5,]   NA        NA          NA        NA         NA 0.77241877 0.32567980
##  [6,]   NA        NA          NA        NA         NA         NA 0.27655807
##  [7,]   NA        NA          NA        NA         NA         NA         NA
##  [8,]   NA        NA          NA        NA         NA         NA         NA
##  [9,]   NA        NA          NA        NA         NA         NA         NA
## [10,]   NA        NA          NA        NA         NA         NA         NA
##             [,8]      [,9]      [,10]
##  [1,] 0.06137235 0.4074293 0.07604571
##  [2,] 0.32198945 0.7555310 0.33228467
##  [3,] 0.19190961 0.1592980 0.22614537
##  [4,] 0.83300480 0.8360162 0.82384713
##  [5,] 0.38193216 0.9271663 0.40015052
##  [6,] 0.32592968 0.7725987 0.33708483
##  [7,] 0.78014078 0.5091956 0.80506799
##  [8,]         NA 0.6211373 0.97979821
##  [9,]         NA        NA 0.61862372
## [10,]         NA        NA         NA
```

Now we can see how many of these p.values are "significant". We know these are false positives, because all the data were generated from the same distribution.


```r
false.positives<-sum(p.values<0.05,na.rm=T)
false.positives
```

```
## [1] 4
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
##         5% 
## 0.00167183
```

```r
false.positives.R<-sum(p.values<new.threshold.R,na.rm=T)
false.positives.R
```

```
## [1] 0
```

If you were to do this experiment (all of the code in the preceding clock) 100 times, you should get at least 1 false positive 5 times, since we have set the threshold such that we have a 5% chance that the smallest p-value in that set of 45 comparisons will be smaller than the threshold we set.
