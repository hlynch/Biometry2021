Week 2 Lab
=============

Confidence intervals
-----------------------

Before getting too far, we need to circle back and make sure we understand what is meant by a confidence interval. 

A 95th percentile confidence interval say “If I repeat this procedure 100 times using 100 different datasets, 95% of the time my confidence intervals will capture the true parameter”. It does NOT say that there is a 95% chance that the parameter is in the interval.

**Quiz time! (Don't worry, not a real quiz)**

*Important note*: This is an area where Aho is **WRONG**. Aho is correct on only one point. It is true that *once the 95th CI has been constructed*, it is no longer possible to assign a $\%$ to the probability that that CI contains the true value or not. Because that CI, once created, either DOES or DOES NOT contain the true value. However, we often talk about the interval in the abstract. **<span style="color: orangered;">When we say "There is a 95$\%$ chance that the interval contains the true value" what we mean is that there is a 95$\%$ probability that a CI created using that methodology would contain the true value.</span>**

Do not let Week 2 pass by without fundamentally understanding the interpretation of a confidence interval. 

Testing hypotheses through permutation
------------------------------------

These examples use data on the speeds of the top 20 racing pigeons from a race in Alma, GA on February 7,2021. 

**Example #1**: Use permutation methods to test whether Cock or Hen birds fly at different speeds (speeds are in meters-per-minute) (in other word: $H_{0}$: No difference in speeds between the C and H groups):

C=$\{1359.8,1355.3,1355.1,1353.0,1349.8,1348.8,1345.2\}$

H=$\{1357.5,1356.4,1355.1,1353.5,1353.2,1352.5,1350.0,1349.8,1346.2,1344.9,1344.4,1343.9,1342.6\}$

**<span style="color: green;">Checkpoint #1: Is this a one-tailed or a two-tailed test?</span>**

Make sure that you understand what is being done here, as this example is very closely related to the problem set.


**Example #2**: Using the same data, provide a 95% confidence interval for the difference in mean speed based on 1000 bootstrap samples

Note that these two approaches are very closely related. Do you see why either approach can be used to test the null hypothesis? **<span style="color: green;">Checkpoint #2: What is the null hypothesis here?</span>**

**Example #3**: Now we will do one slightly more complicated example from Phillip Good's book "Permutation tests: A practical guide to resampling methods and testing hypotheses":

Holmes and Williams (1954) studied tonsil size in children to verify a possible association with the virus \textit{S. pyrogenes}. Test for an association between \textit{S. pyrogenes} status and tonsil size. (Note that you will need to come up with a reasonable test statistic.)

<div class="figure" style="text-align: center">
<img src="Table2categories.png" alt="Data on tonsil size and S. pyrogenes status. Source: Good (1994)" width="40%" />
<p class="caption">(\#fig:unnamed-chunk-1)Data on tonsil size and S. pyrogenes status. Source: Good (1994)</p>
</div>

Now lets consider the full dataset, where tonsil size is divided into three categories. How would we do the test now? **<span style="color: green;">Checkpoint #3: What is the new test statistic? (There are many options.)</span>** What 'labels' do you permute?

<div class="figure" style="text-align: center">
<img src="Table3categories.png" alt="Fill dataset on tonsil size and S. pyrogenes status. Source: Good (1994)" width="50%" />
<p class="caption">(\#fig:unnamed-chunk-2)Fill dataset on tonsil size and S. pyrogenes status. Source: Good (1994)</p>
</div>

Basics of bootstrap and jackknife
------------------------------------

To get started with bootstrap and jackknife techniques, we start by working through a very simple example. First we simulate some data


```r
x<-seq(0,9,by=1)
```

This will constutute our "data". Let's print the result of sampling with replacement to get a sense for it...


```r
table(sample(x,size=length(x),replace=T))
```

```
## 
## 2 3 4 5 7 8 9 
## 2 1 1 2 1 1 2
```

Now we will write a little script to take bootstrap samples and calculate the means of each of these bootstrap samples


```r
xmeans<-vector(length=1000)
for (i in 1:1000)
  {
  xmeans[i]<-mean(sample(x,replace=T))
  }
```

The actual number of bootstrapped samples is arbitrary *at this point* but there are ways of characterizing the precision of the bootstrap (jackknife-after-bootstrap) which might inform the number of bootstrap samples needed. *In practice*, people tend to pick some arbitrary but large number of bootstrap samples because computers are so fast that it is often easy to draw far more samples than are actually needed. When calculation of the statistic is slow (as might be the case if you are using the samples to construct a phylogeny, for example), then you would need to be more concerned with the number of bootstrap samples. 

First, lets just look at a histogram of the bootstrapped means and plot the actual sample mean on the histogram for comparison



```r
hist(xmeans,breaks=30,col="pink")
abline(v=mean(x),lwd=2)
```

<img src="Week-2-lab_files/figure-html/unnamed-chunk-6-1.png" width="672" />

Calculating bias and standard error
-----------------------------------

From these we can calculate the bias and standard deviation for the mean (which is the "statistic"):

$$
\widehat{Bias_{boot}} = \left(\frac{1}{k}\sum^{k}_{i=1}\theta^{*}_{i}\right)-\hat{\theta}
$$


```r
bias.boot<-mean(xmeans)-mean(x)
bias.boot
```

```
## [1] -0.0298
```

```r
hist(xmeans,breaks=30,col="pink")
abline(v=mean(x),lwd=5,col="black")
abline(v=mean(xmeans),lwd=2,col="yellow")
```

<img src="Week-2-lab_files/figure-html/unnamed-chunk-7-1.png" width="672" />

$$
\widehat{s.e._{boot}} = \sqrt{\frac{1}{k-1}\sum^{k}_{i=1}(\theta^{*}_{i}-\bar{\theta^{*}})^{2}}
$$


```r
se.boot<-sd(xmeans)
```

We can find the confidence intervals in two ways:

Method #1: Assume the bootstrap statistics are normally distributed


```r
LL.boot<-mean(xmeans)-1.96*se.boot #where did 1.96 come from?
UL.boot<-mean(xmeans)+1.96*se.boot
LL.boot
```

```
## [1] 2.707397
```

```r
UL.boot
```

```
## [1] 6.233003
```

Method #2: Simply take the quantiles of the bootstrap statistics


```r
quantile(xmeans,c(0.025,0.975))
```

```
##  2.5% 97.5% 
##   2.7   6.3
```

Let's compare this to what we would have gotten if we had used normal distribution theory. First we have to calculate the standard error:


```r
se.normal<-sqrt(var(x)/length(x))
LL.normal<-mean(x)-qt(0.975,length(x)-1)*se.normal
UL.normal<-mean(x)+qt(0.975,length(x)-1)*se.normal
LL.normal
```

```
## [1] 2.334149
```

```r
UL.normal
```

```
## [1] 6.665851
```

In this case, the confidence intervals we got from the normal distribution theory are too wide.

**<span style="color: green;">Checkpoint #4: Does it make sense why the normal distribution theory intervals are too wide?</span>** Because the original were were uniformly distributed, the data has higher variance than would be expected and therefore the standard error is higher than would be expected.

There are two packages that provide functions for bootstrapping, 'boot' and 'boostrap'. We will start by using the 'bootstrap' package, which was originally designed for Efron and Tibshirani's monograph on the bootstrap. 

To test the main functionality of the 'bootstrap' package, we will use the data we already have. The 'bootstrap' function requires the input of a user-defined function to calculate the statistic of interest. Here I will write a function that calculates the mean of the input values.


```r
library(bootstrap)
theta<-function(x)
  {
    mean(x)
  }
results<-bootstrap(x=x,nboot=1000,theta=theta)
results
```

```
## $thetastar
##    [1] 4.5 4.8 5.0 4.5 2.8 4.5 4.3 3.9 5.0 4.3 4.4 5.4 6.2 5.0 4.3 4.0 5.1 4.4
##   [19] 3.1 4.7 4.3 3.4 5.0 4.3 3.6 6.3 3.0 4.4 5.4 5.2 4.2 5.4 3.7 6.4 4.7 3.7
##   [37] 6.0 4.8 4.6 5.1 4.1 3.3 3.9 5.3 2.9 3.9 4.6 4.6 4.1 5.0 5.7 5.7 3.7 3.1
##   [55] 3.1 4.8 3.7 5.1 5.2 4.4 4.2 3.5 3.3 5.3 4.2 4.8 3.4 3.3 3.9 3.2 2.6 4.5
##   [73] 3.6 5.0 5.9 5.5 3.6 5.9 4.2 4.2 5.5 4.6 3.9 5.4 4.4 5.5 4.1 3.2 5.7 4.7
##   [91] 6.2 4.3 4.9 4.6 6.4 6.2 4.6 4.7 6.4 5.3 3.9 4.9 4.2 3.7 3.6 5.5 4.4 5.1
##  [109] 4.7 5.1 4.9 4.8 5.2 4.1 4.3 3.9 3.1 4.6 6.8 4.1 3.2 4.6 4.3 4.6 4.8 3.8
##  [127] 4.5 4.7 4.3 4.2 3.0 5.5 2.8 3.9 5.1 4.6 4.8 4.1 5.0 5.2 4.4 4.2 4.8 5.5
##  [145] 4.1 5.1 3.5 4.1 3.5 4.7 5.1 4.7 3.0 4.9 4.7 3.5 3.6 4.1 4.6 3.5 3.0 3.3
##  [163] 4.7 5.3 4.3 5.7 5.2 4.3 5.2 4.3 4.2 5.8 4.9 4.2 4.1 3.7 4.7 2.7 4.3 4.7
##  [181] 4.0 3.6 4.8 5.1 4.0 5.7 5.7 4.2 5.9 3.7 5.1 5.9 4.5 3.9 5.6 3.1 4.5 4.1
##  [199] 4.1 4.3 4.8 5.0 5.2 5.0 2.9 5.0 3.0 5.8 4.0 5.3 4.2 5.1 3.6 3.6 4.5 4.1
##  [217] 6.2 4.7 5.8 2.9 5.1 3.6 3.8 4.8 3.3 4.8 3.7 4.3 4.5 4.8 3.8 4.3 4.4 5.0
##  [235] 4.2 5.2 5.1 4.8 4.4 4.9 3.6 5.1 4.0 5.2 5.2 2.8 6.2 3.8 4.2 5.1 6.0 4.7
##  [253] 6.0 5.3 6.2 4.8 4.8 5.0 4.7 4.4 4.8 6.2 3.6 5.0 3.9 4.3 5.9 4.5 2.8 3.6
##  [271] 5.5 3.5 3.7 3.9 4.3 3.5 5.5 3.6 3.6 4.5 5.5 3.7 4.5 4.6 4.8 3.9 4.4 2.9
##  [289] 4.1 5.4 4.7 5.4 4.0 3.0 5.4 3.1 2.2 4.1 5.4 5.4 5.4 3.7 3.7 5.1 2.8 4.1
##  [307] 5.3 3.5 2.8 4.7 4.7 4.3 3.9 4.2 5.0 4.0 4.7 4.1 5.6 4.9 5.3 2.8 4.4 5.9
##  [325] 4.3 4.8 7.8 3.1 3.9 6.1 4.5 4.6 3.9 4.6 6.5 5.4 3.5 5.0 4.6 4.9 5.5 5.4
##  [343] 3.6 5.8 3.1 4.3 4.7 4.8 6.5 4.5 3.2 2.6 3.6 4.4 5.9 4.4 3.9 4.5 4.5 5.2
##  [361] 3.9 3.7 6.2 3.9 3.9 4.7 3.5 3.8 3.7 6.0 5.0 5.7 4.2 3.5 3.7 3.3 3.5 4.3
##  [379] 4.2 4.8 5.5 5.0 3.3 4.8 4.5 3.7 4.9 4.9 3.3 2.8 3.6 4.4 4.4 3.9 4.6 5.8
##  [397] 3.9 6.0 3.9 2.5 3.3 6.2 5.5 3.5 4.7 4.9 4.2 3.9 3.4 5.1 3.8 3.2 3.9 2.4
##  [415] 4.1 5.0 3.6 4.2 4.4 5.6 5.6 3.7 3.9 3.3 6.4 4.8 4.1 3.8 4.1 4.2 4.5 4.2
##  [433] 5.5 4.1 3.6 3.1 3.8 3.7 5.0 4.3 5.0 6.2 4.7 4.6 3.5 5.0 4.4 4.2 5.2 4.6
##  [451] 3.7 3.9 4.1 3.6 3.2 4.6 4.1 3.5 3.8 4.0 5.2 4.9 3.7 5.0 4.9 4.4 5.6 4.6
##  [469] 3.6 3.6 5.1 5.1 3.3 2.7 5.7 3.6 5.8 4.1 4.3 5.1 2.5 4.5 4.6 4.7 5.0 4.7
##  [487] 4.3 3.5 5.3 4.2 4.2 6.8 5.6 4.6 6.1 4.2 3.2 5.9 3.9 3.1 4.5 6.4 4.4 3.4
##  [505] 4.6 3.9 4.3 5.5 5.0 4.3 3.1 4.7 5.4 5.0 4.8 5.5 4.6 4.6 4.0 2.9 3.8 4.1
##  [523] 4.3 5.1 3.3 3.7 5.0 5.8 4.1 4.1 3.8 4.6 4.0 6.2 6.1 3.5 5.9 5.0 4.7 5.6
##  [541] 6.4 3.6 3.8 4.2 3.8 3.7 3.7 4.3 4.4 4.7 3.9 3.2 5.4 4.2 4.4 4.7 4.5 4.8
##  [559] 3.9 4.1 4.2 3.9 4.3 2.5 5.9 4.8 5.3 6.1 4.1 4.5 4.3 3.0 5.0 5.4 5.6 5.0
##  [577] 4.0 4.3 4.8 4.3 3.5 4.6 5.0 4.8 4.5 4.2 3.6 5.1 3.7 4.8 5.0 5.2 3.5 4.8
##  [595] 5.6 5.3 3.6 4.6 4.6 6.0 4.9 3.7 3.9 2.2 4.8 4.5 5.0 2.9 3.5 5.7 5.4 4.2
##  [613] 5.3 6.3 4.5 3.7 4.2 3.3 4.4 4.5 3.5 3.9 3.7 4.3 4.7 4.4 5.2 4.3 5.6 5.1
##  [631] 4.5 3.3 4.9 5.9 4.5 3.7 3.8 5.6 3.8 4.2 3.7 4.5 4.7 6.0 4.6 5.6 5.6 5.7
##  [649] 3.5 4.5 5.8 3.1 4.9 4.5 6.8 4.2 4.7 5.9 5.3 4.8 5.4 5.2 4.7 4.5 5.8 3.6
##  [667] 4.2 4.1 4.9 4.4 3.7 3.5 5.1 3.4 3.3 3.4 4.4 3.7 6.0 4.8 5.1 4.8 3.6 5.4
##  [685] 3.8 3.8 3.0 3.5 2.7 4.8 4.7 3.6 4.1 4.6 5.0 5.1 4.8 5.9 4.0 4.9 5.0 4.5
##  [703] 2.9 3.5 4.0 2.9 4.7 4.9 5.6 4.6 3.9 5.6 4.2 5.3 4.9 4.3 2.8 3.4 4.2 4.2
##  [721] 5.8 4.8 4.0 4.0 5.4 4.4 4.0 4.5 5.2 5.0 6.2 4.6 4.5 2.8 4.0 4.1 5.0 5.1
##  [739] 5.8 5.0 4.6 4.4 3.8 4.8 5.8 3.5 5.2 4.8 6.1 4.1 4.1 4.9 2.0 4.9 4.5 3.8
##  [757] 5.5 4.7 4.6 4.6 5.9 3.5 3.7 5.2 6.2 4.0 3.2 4.2 5.5 3.9 4.1 4.0 4.9 5.5
##  [775] 4.9 4.8 4.5 4.3 5.1 4.3 3.6 4.2 4.0 3.8 3.7 4.5 5.0 4.7 3.0 2.0 5.6 4.2
##  [793] 4.6 5.1 4.8 4.1 5.2 6.1 4.5 4.5 5.3 4.3 5.2 4.2 5.6 5.6 4.1 4.6 4.4 3.9
##  [811] 3.2 4.1 4.0 2.5 3.7 6.2 4.6 4.1 3.5 5.1 4.0 7.9 4.3 4.4 3.2 4.1 4.8 5.5
##  [829] 4.2 5.6 5.2 4.9 4.1 4.6 5.3 5.1 3.4 5.3 4.5 5.4 5.5 2.9 4.3 3.7 3.7 3.4
##  [847] 4.3 4.4 3.8 5.3 2.8 5.3 4.1 4.3 4.2 4.6 5.7 6.5 5.0 4.8 4.2 4.6 3.1 4.8
##  [865] 4.8 3.8 6.9 3.3 4.3 4.8 2.8 4.4 3.8 4.5 4.6 3.6 4.8 4.9 5.2 6.0 4.1 4.6
##  [883] 5.1 3.7 4.2 4.7 3.0 4.7 6.1 2.5 2.0 4.2 2.0 7.8 4.7 3.7 3.5 4.9 5.3 5.8
##  [901] 4.4 4.7 4.7 4.0 5.7 3.7 4.7 4.0 3.1 5.4 4.6 4.7 4.0 5.3 4.2 6.3 4.8 4.3
##  [919] 5.5 4.8 5.8 3.6 4.2 4.9 5.8 6.9 4.0 4.9 5.8 4.3 3.6 5.7 4.3 5.5 4.8 4.1
##  [937] 5.4 4.8 4.7 5.1 3.1 4.9 3.7 3.6 5.8 6.4 4.2 3.2 3.9 4.6 3.0 3.8 4.9 3.7
##  [955] 4.7 3.9 2.9 6.0 4.9 5.9 4.7 5.5 5.4 4.6 4.4 2.9 6.0 5.7 6.0 4.5 6.0 4.9
##  [973] 4.7 5.5 5.3 4.7 2.7 4.4 4.6 4.1 3.3 4.7 2.9 5.1 3.7 3.6 5.4 5.7 3.9 3.7
##  [991] 2.5 5.0 4.6 4.2 5.4 3.9 4.3 5.1 4.5 4.4
## 
## $func.thetastar
## NULL
## 
## $jack.boot.val
## NULL
## 
## $jack.boot.se
## NULL
## 
## $call
## bootstrap(x = x, nboot = 1000, theta = theta)
```

```r
quantile(results$thetastar,c(0.025,0.975))
```

```
##  2.5% 97.5% 
##   2.8   6.2
```

Notice that we get exactly what we got last time. This illustrates an important point, which is that the bootstrap functions are often no easier to use than something you could write yourself.

You can also define a function of the bootstrapped statistics (we have been calling this theta) to pull out immediately any summary statistics you are interested in from the bootstrapped thetas.

Here I will write a function that calculates the bias of my estimate of the mean (which is 4.5 [i.e. the mean of the number 0,1,2,3,4,5,6,7,8,9])


```r
bias<-function(x)
  {
  mean(x)-4.5
  }
results<-bootstrap(x=x,nboot=1000,theta=theta,func=bias)
results
```

```
## $thetastar
##    [1] 5.8 3.1 5.4 5.7 3.0 3.4 5.4 4.3 5.0 3.1 4.8 5.1 4.0 3.9 2.2 5.1 3.9 6.2
##   [19] 4.1 3.7 4.5 7.1 4.3 5.3 4.9 2.7 3.2 4.4 4.2 4.8 3.5 4.0 4.4 3.5 5.0 3.4
##   [37] 5.5 3.8 5.3 4.0 4.5 6.0 4.2 5.3 3.8 3.8 5.8 4.9 5.5 3.8 4.3 4.4 4.2 4.0
##   [55] 5.2 4.4 5.1 2.6 5.2 6.5 3.8 6.1 2.9 5.5 4.5 4.4 5.0 5.7 5.4 3.8 5.8 3.9
##   [73] 5.3 6.0 5.6 4.6 4.5 2.6 6.3 5.7 4.1 5.9 6.1 4.6 4.7 4.3 2.8 5.2 4.0 3.8
##   [91] 6.0 6.1 5.0 5.5 5.8 5.0 4.2 4.8 4.4 2.6 5.3 3.1 4.3 6.2 5.6 5.9 5.2 4.7
##  [109] 4.4 4.0 4.7 4.5 4.6 3.6 4.8 4.1 4.4 5.1 4.1 3.6 4.6 3.8 5.6 3.8 3.4 4.7
##  [127] 2.9 4.0 4.0 5.3 3.8 4.4 4.3 3.9 3.2 4.6 4.8 5.2 4.9 3.7 6.5 5.6 3.9 4.8
##  [145] 5.2 4.3 5.0 4.2 3.8 3.9 3.1 4.5 4.0 4.3 5.3 5.0 4.7 4.1 4.8 4.5 4.7 6.9
##  [163] 4.9 5.9 3.5 4.1 5.8 4.2 4.3 4.3 4.1 4.8 4.4 3.9 4.8 4.4 4.3 4.7 5.0 5.6
##  [181] 4.4 4.4 5.2 4.9 4.2 3.3 5.2 4.5 5.4 3.8 5.1 4.5 4.5 4.5 5.3 4.9 5.2 5.3
##  [199] 3.8 4.0 3.5 4.2 4.5 4.5 3.5 3.8 4.2 6.0 4.6 4.5 4.7 5.6 3.5 3.2 4.6 4.9
##  [217] 4.2 5.5 3.8 5.5 6.3 2.7 2.7 5.5 4.3 5.9 4.0 4.5 4.8 6.6 3.3 3.5 3.7 3.4
##  [235] 5.4 3.7 4.9 3.8 4.2 5.2 6.0 6.6 4.5 3.0 4.2 4.4 4.1 3.2 3.7 4.8 4.3 5.1
##  [253] 2.6 5.6 4.2 4.4 4.2 3.8 4.0 6.8 3.6 3.4 4.7 4.4 4.6 4.4 3.9 4.6 4.9 4.3
##  [271] 4.7 4.7 4.4 4.2 5.9 5.5 5.4 6.5 5.5 4.4 3.4 4.1 4.2 4.5 5.8 4.2 4.9 5.1
##  [289] 5.4 4.7 3.9 4.5 4.6 4.9 4.9 3.6 4.4 3.8 4.9 6.3 4.8 4.9 3.8 5.4 4.6 4.4
##  [307] 4.1 5.7 5.6 4.5 5.5 3.8 4.9 3.9 4.5 5.3 5.4 4.7 4.7 3.7 4.2 4.2 3.3 4.9
##  [325] 4.6 4.6 4.9 3.7 3.8 3.3 2.9 2.4 3.4 3.4 7.0 4.3 6.2 5.9 4.4 6.1 5.6 5.4
##  [343] 4.4 4.7 4.5 5.3 3.5 4.7 6.4 5.2 2.7 5.3 4.5 4.9 5.0 3.7 4.0 4.2 5.7 5.2
##  [361] 5.7 6.0 3.2 5.7 4.3 3.0 4.5 5.6 5.5 6.1 3.6 5.2 4.9 6.0 2.8 4.6 4.1 3.9
##  [379] 6.3 5.1 4.2 4.7 5.5 3.5 4.2 4.3 5.6 4.4 4.8 4.2 5.9 4.3 3.8 4.6 5.7 3.8
##  [397] 3.8 3.5 4.3 4.7 4.0 4.8 5.6 4.8 5.1 3.8 4.2 3.2 5.1 4.4 4.5 4.9 5.7 4.5
##  [415] 4.6 5.1 5.5 5.2 4.5 4.4 4.0 4.4 4.2 3.4 5.8 4.0 4.3 5.6 3.1 3.2 3.5 3.5
##  [433] 4.6 4.5 4.1 4.6 3.9 3.9 6.1 4.1 6.2 4.7 5.0 4.1 4.5 3.4 5.3 5.3 4.6 3.4
##  [451] 5.4 4.3 4.7 4.7 3.5 4.1 5.9 4.1 4.8 5.3 4.6 4.8 5.2 6.0 4.3 5.2 3.7 6.7
##  [469] 4.3 5.5 4.1 4.9 4.1 6.6 4.5 3.6 4.3 4.5 2.0 5.9 4.1 6.0 3.9 4.2 4.3 4.0
##  [487] 2.9 4.8 4.0 6.2 3.2 4.2 4.1 5.4 5.8 3.9 4.5 4.0 4.6 3.7 3.7 3.8 4.4 7.1
##  [505] 3.0 4.9 4.7 4.9 5.3 5.7 4.3 5.8 4.8 3.6 2.5 4.3 4.7 3.9 4.0 3.7 4.3 4.5
##  [523] 3.8 4.3 5.1 5.8 5.2 3.4 5.6 4.4 4.5 5.7 5.8 3.5 5.2 3.2 4.2 4.7 4.1 3.7
##  [541] 4.3 5.0 5.0 3.6 4.0 3.9 3.4 5.0 4.6 4.8 4.7 2.4 5.6 5.0 4.9 4.5 3.1 3.2
##  [559] 4.9 5.2 4.3 4.9 5.7 3.9 4.0 5.3 4.6 4.6 4.2 4.8 4.8 3.7 5.4 6.2 4.6 3.7
##  [577] 4.1 5.4 5.3 2.7 5.7 4.8 4.6 2.4 4.3 4.1 5.5 5.9 5.0 4.3 6.6 4.2 3.0 4.6
##  [595] 3.1 5.1 3.5 3.6 4.7 3.2 4.5 3.6 6.0 3.9 5.2 4.3 6.0 3.8 4.3 4.8 4.4 2.7
##  [613] 4.4 3.6 5.5 4.9 4.2 5.9 5.3 6.0 3.9 2.5 4.4 3.2 5.1 3.3 4.3 4.2 6.5 3.8
##  [631] 6.4 4.9 4.0 3.6 5.3 7.2 5.1 3.2 3.2 4.2 3.4 4.0 4.8 4.7 4.3 5.5 4.3 3.7
##  [649] 4.8 4.3 3.6 2.8 4.6 4.9 3.4 4.9 5.4 3.6 5.1 4.3 5.0 6.6 4.4 4.0 6.0 4.1
##  [667] 3.5 3.3 3.5 4.5 4.9 3.8 4.5 4.5 4.7 4.2 4.8 3.1 3.9 5.7 4.8 4.4 4.8 4.7
##  [685] 5.1 5.1 4.7 4.3 4.5 3.4 5.3 4.3 4.8 3.5 4.6 4.2 4.7 4.1 4.5 4.5 4.7 4.0
##  [703] 4.7 5.2 4.8 4.0 3.2 4.3 5.0 5.3 5.6 4.2 4.4 5.0 4.6 4.1 3.9 5.3 3.6 3.5
##  [721] 3.9 4.9 5.7 5.8 4.5 2.3 4.4 4.2 4.1 5.7 5.6 4.7 4.8 4.1 4.6 3.1 4.4 3.9
##  [739] 4.7 6.0 4.4 5.0 4.7 4.4 4.7 5.0 4.6 3.9 4.3 6.3 4.7 5.4 4.2 3.7 4.1 4.1
##  [757] 2.9 4.9 4.5 3.1 5.4 4.6 4.6 6.1 5.9 4.6 4.2 3.8 5.1 5.5 4.5 4.5 4.8 4.7
##  [775] 3.9 4.4 5.3 5.4 5.0 4.1 4.4 3.5 5.0 4.3 2.6 4.1 3.4 6.0 6.2 4.4 5.2 3.7
##  [793] 5.1 2.6 5.4 3.1 4.2 4.0 3.6 4.0 4.6 4.0 2.6 4.8 4.0 4.4 4.3 4.8 3.8 5.3
##  [811] 3.2 5.0 3.3 4.3 3.3 4.7 4.4 5.1 6.0 4.1 4.3 4.4 4.4 4.6 4.4 3.2 5.0 3.3
##  [829] 4.2 5.6 5.1 4.9 3.6 3.7 4.8 4.9 3.8 3.7 4.2 5.2 4.2 3.4 5.1 3.2 4.6 4.6
##  [847] 4.4 3.3 5.2 5.5 5.3 4.1 5.8 3.7 3.1 4.0 3.3 3.7 5.7 3.4 4.1 4.5 4.8 3.9
##  [865] 3.3 4.3 3.2 5.3 3.7 5.0 5.2 2.8 5.8 4.2 3.2 3.6 4.4 3.6 4.1 2.4 3.0 5.3
##  [883] 3.4 5.7 4.5 4.3 4.0 3.9 4.9 5.9 4.5 3.9 3.7 4.6 6.1 5.0 4.5 5.2 5.3 4.5
##  [901] 4.9 3.8 3.8 3.3 4.9 6.1 5.2 4.7 6.0 6.1 3.9 3.7 4.8 5.6 3.9 3.9 4.6 5.5
##  [919] 4.0 4.7 4.6 4.9 4.5 4.7 4.5 6.3 4.6 4.9 4.5 4.1 4.5 4.9 5.9 5.8 4.5 3.1
##  [937] 4.2 5.5 5.8 3.9 4.4 5.0 3.5 5.7 2.6 3.9 3.6 3.9 4.0 4.7 5.2 4.9 5.4 4.9
##  [955] 3.9 4.2 4.0 4.2 2.9 4.4 3.4 4.7 3.7 4.8 6.3 5.4 5.2 3.6 3.8 5.1 4.0 5.7
##  [973] 3.2 4.2 5.0 4.7 3.6 5.9 3.5 4.9 4.8 4.9 4.6 4.7 4.4 4.2 4.7 5.0 4.3 4.5
##  [991] 3.6 2.7 4.5 4.4 5.0 4.6 6.0 3.1 4.9 3.6
## 
## $func.thetastar
## [1] 0.0182
## 
## $jack.boot.val
##  [1]  0.52123894  0.41144578  0.26724138  0.20504202  0.07316384 -0.03473389
##  [7] -0.14732620 -0.21744548 -0.31805556 -0.46212121
## 
## $jack.boot.se
## [1] 0.917722
## 
## $call
## bootstrap(x = x, nboot = 1000, theta = theta, func = bias)
```

Compare this to 'bias.boot' (our result from above). Why might it not be the same? Try running the same section of code several times. See how the value of the bias ($func.thetastar) jumps around? We should not be surprised by this because we can look at the jackknife-after-bootstrap estimate of the standard error of the function (in this case, that function is the bias) and we can see that it is not so small that we wouldn't expect some variation in these values.

Remember, everything we have discussed today are estimates. The statistic as applied to your data will change with new data, as will the standard error, the confidence intervals - everything! All of these values have sampling distributions and are subject to change if you repeated the procedure with new data.

Note that we can calculate any function of $\theta^{*}$. A simple example would be the 72nd percentile:


```r
perc72<-function(x)
  {
  quantile(x,probs=c(0.72))
  }
results<-bootstrap(x=x,nboot=1000,theta=theta,func=perc72)
results
```

```
## $thetastar
##    [1] 4.0 4.6 5.0 3.4 5.0 3.1 3.5 3.4 4.2 3.2 4.2 4.7 3.3 4.1 4.5 4.5 5.3 5.9
##   [19] 3.9 4.0 4.7 5.3 4.8 4.3 3.2 2.5 3.6 5.9 5.6 4.0 4.9 4.1 4.3 4.0 5.6 4.8
##   [37] 4.5 4.0 4.5 4.8 4.1 4.6 5.4 3.6 5.5 3.3 4.1 4.5 3.0 6.1 4.1 4.3 5.0 3.3
##   [55] 6.0 4.6 4.1 5.0 3.3 3.9 3.5 3.9 4.3 5.6 4.3 3.6 3.8 3.8 3.4 5.2 4.3 4.7
##   [73] 5.4 5.2 4.0 4.5 3.3 3.8 6.0 4.6 3.9 2.7 3.6 4.2 5.0 3.7 4.7 4.2 4.3 4.8
##   [91] 5.3 2.8 5.4 3.7 4.8 5.5 3.3 5.0 6.1 3.7 5.1 4.7 4.0 3.9 4.3 5.2 4.8 5.9
##  [109] 5.4 3.9 3.8 2.7 4.2 4.0 4.1 3.3 3.3 5.9 5.5 3.3 3.8 5.3 4.4 5.3 3.9 4.9
##  [127] 4.1 3.6 5.4 4.1 5.8 5.0 2.8 4.7 3.8 5.3 4.1 4.9 3.2 4.9 4.0 5.0 4.0 5.8
##  [145] 4.8 5.6 4.9 4.5 5.0 3.7 5.1 3.6 5.6 5.7 5.6 4.2 5.7 4.9 5.2 2.4 3.2 3.7
##  [163] 5.8 3.6 4.3 5.7 2.9 4.2 5.1 4.0 4.5 4.3 4.5 4.9 4.1 5.7 3.4 5.1 4.8 5.1
##  [181] 4.4 4.5 3.4 5.3 3.9 2.7 4.6 4.5 3.8 4.3 5.4 5.0 5.1 5.6 2.9 4.2 4.9 5.2
##  [199] 6.3 3.5 4.9 4.4 4.7 4.9 4.5 2.9 4.3 5.4 4.6 4.4 5.9 3.9 3.7 4.8 3.5 3.7
##  [217] 5.4 4.1 3.2 3.5 4.0 3.8 4.8 3.5 5.8 2.6 4.4 2.6 4.2 4.8 3.5 3.2 3.7 4.5
##  [235] 5.3 5.6 3.8 5.6 4.2 4.3 5.5 4.4 4.9 5.2 3.5 3.3 5.3 4.6 5.3 3.9 4.4 3.8
##  [253] 4.6 5.6 4.4 4.3 3.3 4.8 5.0 5.3 3.0 4.6 3.9 3.3 2.7 4.6 6.2 3.6 6.3 5.5
##  [271] 5.4 5.9 4.9 5.1 4.5 2.9 4.8 5.0 4.0 4.4 6.2 4.1 3.2 1.7 6.7 5.3 5.8 4.4
##  [289] 4.0 4.8 5.1 2.6 3.3 3.1 3.6 4.4 5.7 3.1 4.0 3.4 3.6 5.4 3.5 4.7 3.6 2.8
##  [307] 4.9 5.2 5.6 4.9 3.5 5.3 4.0 2.9 3.1 4.1 4.9 3.5 6.9 4.7 6.5 5.4 5.0 3.5
##  [325] 5.4 5.0 4.7 4.5 6.2 4.0 4.8 3.9 2.4 5.3 4.2 5.1 1.8 3.1 4.9 5.1 5.5 4.9
##  [343] 6.2 4.8 5.0 4.6 4.6 3.6 4.7 3.0 5.7 3.5 4.4 5.2 4.0 4.8 5.5 5.3 4.5 4.5
##  [361] 4.3 4.7 4.1 4.2 4.0 5.4 2.8 3.8 3.5 3.4 4.7 3.7 4.6 4.0 4.3 5.3 6.4 3.2
##  [379] 3.8 5.5 4.4 3.5 7.4 3.9 3.4 4.5 4.9 4.2 5.1 4.8 5.6 3.9 4.7 4.0 4.9 5.3
##  [397] 4.6 2.7 5.0 3.4 2.3 4.4 4.6 3.7 3.8 2.9 4.1 4.0 5.3 3.4 4.5 4.7 4.7 4.8
##  [415] 3.5 5.5 4.5 4.9 3.4 7.0 3.5 4.2 5.1 4.3 4.8 3.5 5.5 5.4 4.8 4.6 3.5 6.5
##  [433] 4.7 3.6 4.4 4.1 5.0 5.4 3.3 4.7 3.9 3.9 4.3 5.9 3.3 3.6 4.6 4.2 3.8 4.3
##  [451] 4.4 5.1 4.5 5.2 4.0 4.8 4.5 4.5 4.8 3.8 3.7 4.5 4.5 4.9 5.4 4.9 3.9 4.7
##  [469] 4.2 3.0 5.7 5.9 5.5 6.5 2.7 5.2 3.3 5.6 6.7 4.9 4.5 4.3 4.3 5.1 4.2 5.6
##  [487] 4.9 4.7 5.6 2.9 3.8 5.8 3.8 4.4 4.2 4.5 5.0 4.4 3.9 5.2 5.3 4.3 4.9 4.2
##  [505] 6.0 3.8 5.4 4.9 3.5 3.9 5.6 4.1 4.1 4.2 3.3 4.0 3.0 5.3 5.6 4.7 4.7 5.6
##  [523] 3.7 5.7 4.7 3.9 5.2 4.6 3.6 5.4 3.0 5.1 6.0 4.0 3.6 6.4 4.3 5.5 2.5 5.2
##  [541] 5.7 5.9 3.8 4.4 4.5 4.3 5.1 4.6 6.3 4.5 5.6 3.4 5.3 3.5 4.4 5.6 4.1 3.8
##  [559] 4.8 3.9 4.2 4.4 3.8 5.5 4.5 5.1 6.0 4.7 5.4 4.7 2.4 3.1 3.4 6.1 4.2 5.1
##  [577] 4.7 6.1 2.9 4.2 5.3 3.4 4.6 5.7 4.9 4.1 4.9 5.7 3.7 2.9 4.0 4.4 5.8 3.8
##  [595] 4.9 3.4 5.2 2.7 3.2 3.9 5.5 4.7 3.4 3.3 2.5 2.9 4.9 2.9 2.4 5.1 4.8 3.8
##  [613] 4.5 3.3 2.7 3.6 4.2 4.0 4.5 3.5 5.5 4.7 3.8 4.5 4.3 4.5 3.7 4.6 5.1 3.7
##  [631] 3.6 5.7 5.0 5.1 5.9 4.6 5.5 4.4 3.3 3.5 4.5 5.9 6.1 4.7 4.8 4.9 4.7 3.4
##  [649] 5.2 3.0 4.3 6.5 4.1 4.4 5.6 5.0 4.5 5.2 6.0 5.2 5.3 3.5 5.7 3.5 5.0 3.9
##  [667] 5.5 4.8 4.4 6.0 5.0 4.0 3.0 5.4 5.2 3.0 5.4 5.2 6.3 2.9 3.7 3.5 3.2 4.0
##  [685] 4.3 5.0 5.0 3.8 4.9 5.2 3.3 4.4 5.2 5.7 3.2 4.4 4.2 3.0 4.5 5.0 4.6 4.2
##  [703] 4.3 4.5 5.7 5.0 5.5 4.7 4.1 4.8 5.0 3.2 2.7 3.3 3.8 6.0 3.4 5.4 3.0 5.9
##  [721] 3.7 2.5 4.8 5.2 4.5 7.0 2.9 3.7 3.5 2.3 3.8 4.4 3.9 3.1 5.1 5.6 3.5 3.9
##  [739] 6.0 4.1 3.8 4.2 5.1 3.3 3.6 3.8 4.7 5.2 3.7 4.1 2.6 3.8 5.0 4.8 5.8 5.1
##  [757] 3.0 4.2 5.2 4.1 5.5 4.8 4.0 4.6 4.1 4.5 4.4 3.9 6.1 5.0 6.0 3.7 4.7 5.3
##  [775] 5.9 4.9 2.2 3.8 5.6 4.8 3.7 5.0 5.4 5.8 5.7 3.4 5.7 5.1 5.4 4.3 4.3 4.5
##  [793] 4.1 5.6 3.3 5.1 3.8 4.6 3.7 5.2 3.8 4.9 5.4 4.5 4.8 5.7 5.8 4.4 3.7 4.0
##  [811] 3.7 4.1 4.7 4.8 3.8 3.4 6.0 4.0 5.5 4.2 3.3 5.4 4.4 5.4 4.3 4.3 4.7 5.5
##  [829] 4.1 4.8 5.2 5.0 5.0 5.1 3.8 4.2 5.6 3.6 5.1 5.4 4.0 4.5 4.0 4.5 5.4 4.1
##  [847] 5.1 6.6 4.6 3.7 4.6 5.4 4.1 3.8 5.1 5.9 5.6 5.9 5.3 3.3 6.0 4.6 5.4 4.5
##  [865] 4.8 4.0 3.8 3.3 3.6 5.5 5.0 3.4 2.7 5.2 4.6 5.1 5.0 4.0 5.0 5.7 5.3 3.2
##  [883] 4.0 2.4 5.4 4.9 4.0 3.6 4.6 5.8 4.7 3.9 4.9 5.8 4.4 5.6 4.5 5.6 4.9 4.1
##  [901] 3.8 4.4 3.8 5.0 4.0 4.5 5.0 3.9 5.5 3.2 3.3 4.3 3.6 4.3 4.2 4.1 5.1 5.8
##  [919] 5.9 4.0 3.6 4.0 4.6 5.1 5.7 2.9 5.8 4.6 5.4 5.4 4.8 4.1 4.4 4.8 4.9 5.9
##  [937] 4.2 5.3 4.0 5.1 4.5 4.4 3.5 3.9 4.8 5.2 4.6 4.4 5.0 5.5 3.6 5.9 4.2 3.1
##  [955] 4.8 5.1 3.9 4.1 3.6 4.4 4.0 3.9 4.9 3.9 5.7 5.2 2.5 3.7 4.5 4.7 5.1 5.1
##  [973] 5.2 2.9 4.8 4.0 4.1 3.9 4.6 3.3 4.2 1.9 4.4 4.0 4.0 4.6 3.1 4.6 4.2 5.1
##  [991] 3.7 4.9 3.9 5.0 4.5 5.6 5.8 5.6 5.7 3.8
## 
## $func.thetastar
## 72% 
##   5 
## 
## $jack.boot.val
##  [1] 5.400 5.444 5.300 5.200 5.064 5.000 4.900 4.700 4.500 4.400
## 
## $jack.boot.se
## [1] 1.039397
## 
## $call
## bootstrap(x = x, nboot = 1000, theta = theta, func = perc72)
```

On Tuesday we went over an example in which we bootstrapped the correlation coefficient between LSAT scores and GPA. To do that, we sampled pairs of (LSAT,GPA) data with replacement. Here is a little script that would do something like that using (X,Y) data that are independently drawn from the normal distribution


```r
xdata<-matrix(rnorm(30),ncol=2)
```

Everyone's data is going to be different. With such a small sample size, it would be easy to get a positive or negative correlation by random change, but on average across everyone's datasets, there should be zero correlation because the two columns are drawn independently.


```r
n<-15
theta<-function(x,xdata)
  {
  cor(xdata[x,1],xdata[x,2])
  }
results<-bootstrap(x=1:n,nboot=50,theta=theta,xdata=xdata) 
#NB: xdata is passed to the theta function, not needed for bootstrap function itself
```

Notice the parameters that get passed to the 'bootstrap' function are: (1) the indexes which will be sampled with replacement. This is different that the raw data but the end result is the same because both the indices and the raw data get passed to the function 'theta' (2) the number of bootrapped samples (in this case 50) (3) the function to calculate the statistic (4) the raw data.

Lets look at a histogram of the bootstrapped statistics $\theta^{*}$ and draw a vertical line for the statistic as applied to the original data.


```r
hist(results$thetastar,breaks=30,col="pink")
abline(v=cor(xdata[,1],xdata[,2]),lwd=2)
```

<img src="Week-2-lab_files/figure-html/unnamed-chunk-17-1.png" width="672" />

Parametric bootstrap
---------------------

Let's do one quick example of a parametric bootstrap. We haven't introduced distributions yet (except for the Gaussian, or Normal, distribution, which is the most familiar), so lets spend a few minutes exploring the Gamma distribution, just so we have it to work with for testing out parametric bootstrap. All we need to know is that the Gamma distribution is a continuous, non-negative distribution that takes two parameters, which we call "shape" and "rate". Lets plot a few examples just to see what a Gamma distribution looks like. (Note that the Gamma distribution can be parameterized by "shape" and "rate" OR by "shape" and "scale", where "scale" is just 1/"rate". R will allow you to use either (shape,rate) or (shape,scale) as long as you specify which you are providing.

<img src="Week-2-lab_files/figure-html/unnamed-chunk-18-1.png" width="672" />


Let's generate some fairly sparse data from a Gamma distribution


```r
original.data<-rgamma(10,3,5)
```

and calculate the skew of the data using the R function 'skewness' from the 'moments' package. 


```r
library(moments)
theta<-skewness(original.data)
head(theta)
```

```
## [1] 0.9816691
```

What is skew? Skew describes how assymetric a distribution is. A distribution with a positive skew is a distribution that is "slumped over" to the right, with a right tail that is longer than the left tail. Alternatively, a distribution with negative skew has a longer left tail. Here we are just using it for illustration, as a property of a distribution that you may want to estimate using your data.

Lets use 'fitdistr' to fit a gamma distribution to these data. This function is an extremely handy function that takes in your data, the name of the distribution you are fitting, and some starting values (for the estimation optimizer under the hood), and it will return the parameter values (and their standard errors). We will learn in a couple weeks how R is doing this, but for now we will just use it out of the box. (Because we generated the data, we happen to know that the data are gamma distributed. In general we wouldn't know that, and we will see in a second that our assumption about the shape of the data really does make a difference.)


```r
library(MASS)
fit<-fitdistr(original.data,dgamma,list(shape=1,rate=1))
```

```
## Warning in densfun(x, parm[1], parm[2], ...): NaNs produced
```

```r
# fit<-fitdistr(original.data,"gamma")
# The second version would also work.
fit
```

```
##     shape       rate  
##   4.105315   7.505225 
##  (1.766425) (3.435208)
```

Now lets sample with replacement from this new distribution and calculate the skewness at each step:


```r
results<-c()
for (i in 1:1000)
  {
  x.star<-rgamma(length(original.data),shape=fit$estimate[1],rate=fit$estimate[2])
  results<-c(results,skewness(x.star))
  }
head(results)
```

```
## [1]  1.4211517  1.2891195 -0.2071598  0.9740248 -0.1989893  0.2485678
```

```r
hist(results,breaks=30,col="pink",ylim=c(0,1),freq=F)
```

<img src="Week-2-lab_files/figure-html/unnamed-chunk-22-1.png" width="672" />

Now we have the bootstrap distribution for skewness (the $\theta^{*}$ s), we can compare that to the equivalent non-parametric bootstrap:


```r
results2<-bootstrap(x=original.data,nboot=1000,theta=skewness)
results2
```

```
## $thetastar
##    [1]  1.6057449630  0.6788923883 -0.3034483165  0.4885899677  1.0494095697
##    [6]  0.9226380447  1.9690493012  0.5416395242 -0.0136641221  0.5899723017
##   [11]  1.7214367328  0.5715112084 -0.3004428052  0.7751950311  1.2708355127
##   [16]  1.3152853038  0.6406762079  0.8744481401  0.9735550939  1.1970255892
##   [21]  0.7616201881  1.2850922379  0.0255238006  0.3718577925  1.1765539811
##   [26]  1.1808862158  0.6857305293  0.4838681525  1.2415566885  0.6958492112
##   [31]  1.0672763539  0.0005921511  0.3354901757  0.5274104708  1.1654197656
##   [36]  1.0356947120  0.8698465472  0.6602917769  1.4472088650  0.6156482056
##   [41] -0.5476869736 -0.0147997228  0.9697107810  1.2893461241  0.6969687914
##   [46]  0.0211744119  1.0367691853  1.4183792449 -0.5797850683 -0.0760219915
##   [51]  0.5740688555  1.1020423815  1.4644064310  1.9405901589  2.0970757261
##   [56]  0.5238878679  1.3716402442  1.4557988963  0.7334089497  1.1558826982
##   [61]  1.0767264635  0.4326118186 -1.7102613402 -0.4104082714  0.7520364318
##   [66]  0.7349106101  0.6312108387  0.5599206062  1.0998943413  1.3144289763
##   [71]  1.1627482600  1.7954602441  0.9092910920 -1.4382207109  0.2721368051
##   [76]  1.7672347853  2.0437133191  0.7295270555  1.2263798091  0.2237626838
##   [81]  1.7972288461  0.5468607367  0.7604182112  1.0983759603  0.8363733711
##   [86]  1.1315504042  0.5247326790  1.0622333732  1.4437932300  1.2650783582
##   [91]  0.5295092594  1.0220311371  1.1322133340  1.0052827606  1.9193717065
##   [96]  2.1182372429 -0.3345385785  0.6912420835  0.1531540648  0.8784297162
##  [101]  0.4881581149  0.8880875061 -0.3490144474  1.4474474742  2.1954182634
##  [106]  1.0394168693 -1.0368115356 -0.2678196961  0.2177533886  0.5746490101
##  [111]  0.6624857657  0.6620723414  1.2609002600  0.1092702577  0.6279645894
##  [116]  0.5449297813  0.1500319630  0.9435973838 -0.4969925953  0.1580885612
##  [121]  0.6767600965  0.5672759009  0.9502400863  1.3492657952  1.0545370197
##  [126]  1.3473727910  1.3669016079  1.0919573130  0.8360924163  0.5576147216
##  [131]  0.3907948498  1.2893064773  1.8389331942  0.9735550939  0.6000017483
##  [136]  2.0513187314  1.0322791029  1.9672595525  1.4704403380  0.9439628885
##  [141]  0.9899534869  1.6841099695  1.0329890325 -0.1308978760  0.5853859513
##  [146]  2.2848658879  0.6373394087  1.4029986948  0.2345864150  0.0864283239
##  [151]  1.4903172895  1.0924539096  0.9574927646 -0.2621017254  1.8470315721
##  [156]  0.7819918525  0.5921621377  1.0409377679  1.8705900579 -0.0973643541
##  [161]  1.0756856034  1.1953480500  1.3143431370  1.3312351209  0.9661297832
##  [166]  0.4640403323  1.1371573792  0.6246406763  0.8564510159  2.1002933480
##  [171] -0.6754386813  0.0858748187  1.4389752924  1.0954239939 -0.6057728497
##  [176] -0.1128487700  0.6930675785  1.0038030566  0.9120639768 -0.6959114811
##  [181]  0.4172453634  1.8820637113  0.9748246459  0.9899534869 -0.4331067455
##  [186]  1.7170481712  1.1844863861  0.8192906643  0.7983461878  0.9299907213
##  [191]  0.3741759117  1.4755041192  1.6892678993  0.1569712094  1.0613493972
##  [196]  1.6206281378  1.0284834533  1.0412663477  0.1681782095  1.5805624421
##  [201]  1.1977689567  0.3612404218  1.1653724369  0.6453547876  0.5801476348
##  [206]  0.9995573551  1.1755361624  1.3464158889 -0.3706207461 -0.6293415439
##  [211]  0.9638659804  0.5763135742  1.9408675720  0.6420730364  0.7846957476
##  [216]  0.6425868736  0.7694148514  1.7295158381  0.9389613259  1.9110780547
##  [221]  0.0195275321  0.3948675410 -0.1916402459  1.1751766327 -0.0813611293
##  [226]  0.6592723242  0.9183676442  1.5345313508  1.3727771237  1.1705373606
##  [231] -0.5450671894  0.3264213435 -1.8563770625  1.8608575167  1.1431177671
##  [236]  1.1560156145  1.6576074280 -0.7445881354 -0.1941682280  1.2016914819
##  [241]  2.1163364219  0.5007316359  1.1082172049  0.1256421462 -0.4493942591
##  [246]  0.8892182404  0.5069140700 -0.5571761570  1.2540154945 -0.4137558223
##  [251]  1.0589539800  1.7569310554  0.9778526302  1.1256810954  1.0332585242
##  [256]  1.1490749273  1.0206317342  1.2567328904  0.9344975291  1.0361593279
##  [261]  0.6801911422  1.2202316048  1.6528210327  1.4510399573 -0.4549585562
##  [266]  0.2432236270  1.0760354884  1.0768582482 -0.0455861548 -0.1411646936
##  [271]  0.2179465649 -1.1710379206  1.0248839768  0.2731282563  1.0328736212
##  [276]  0.6552031437  1.0406688813  0.3937124920  0.6292323919  1.3358696056
##  [281]  0.9254340799  1.3180306239  0.1110533829  0.0357809519  0.6492086498
##  [286]  0.4726977627  1.0974384504  1.3758139327  1.0723286709  2.3039116604
##  [291]  0.6179659982  1.5734053422  2.0884720783  0.4252815441  0.5256608322
##  [296]  0.6514111512  0.2004690176 -0.5507345969  1.5915531797 -0.4686935633
##  [301]  0.4908446956  0.6082146358 -0.5868151514  1.8086943897  1.1623156759
##  [306]  0.5590277935  0.6784598389  0.7325190177  1.3235711977  1.0290871356
##  [311]  1.0529541061  0.2083625530  1.0742196984  0.2068877356  0.4264827938
##  [316]  1.6737468234  0.9618344516  0.2796290333 -1.9168367563 -0.0321296940
##  [321] -0.6412316227  1.1484965115 -0.3254985145  0.8295548376  0.5842271504
##  [326]  1.8743328901  0.8841254316  1.1694007345  0.5771905520 -0.2948516873
##  [331]  0.7128370213  1.5435285614  0.5117346328  0.5298802671  0.9692043661
##  [336]  0.9812876941  1.2128389816  2.0567667902  0.6185290092 -0.0865952698
##  [341]  1.8642615338  0.8572521819  1.8036222289  0.7946256674  0.7984507137
##  [346]  0.7038294605  0.8426377683  0.1571033067  1.3555854310  0.9272604815
##  [351]  0.6616263497  1.0695251436  0.4235001793  1.0239678973  0.5655562595
##  [356] -0.3436082539  1.4503013712  1.2358573874  1.3512281057  1.5551449585
##  [361]  1.1575134669  0.8111000241  0.6006418959  0.8549037796  0.7203354091
##  [366]  0.8672654488  1.1967319088  1.2427085074  0.5516829613  0.6105696076
##  [371]  0.9591241499  0.8910444327  0.6929320998  0.2368638800 -0.1989101174
##  [376] -0.1698169731  0.8487527627  1.2619177127  1.0943122659  0.2248525153
##  [381]  1.0990979056  0.2422901713  0.6025491298  1.0443395301  0.4771682164
##  [386]  0.7892430651 -0.5571494222  0.3747266897  1.3805646824  1.0426628875
##  [391]  1.8926765324 -0.1453377618  0.6494015935  1.0293468878  1.1075180112
##  [396]  0.7627078608  0.1022496834  0.6893995610  0.5940796471 -0.5074924296
##  [401]  1.2711436881  1.1064236995  0.6422253470  1.2112998474  1.6373020191
##  [406]  0.3747266897  1.5452312580  0.7851293484  1.0936485354  0.6359928140
##  [411]  0.6264051307  0.1306392434  1.0950681323  1.0989068124  1.0097259080
##  [416]  2.1621839509 -0.0472418730  0.2865268883 -0.2455900847  0.2944244785
##  [421]  0.5279914255  1.1653724369  1.6192139542  0.7890172638  1.7279996365
##  [426]  1.1774666359  0.1284124869  2.1167465116  1.2509463363 -0.2157768580
##  [431]  0.9816691257  0.4952399419  0.6359928140  1.3290841674  0.8196450469
##  [436] -0.2266173220  0.9932685874  1.0613293211  0.6761283279  1.1804111264
##  [441]  1.6598818358  1.3471945160  0.1996818726 -0.9096340649  0.8920595643
##  [446]  0.5899938791  1.1562355446  1.8045997469  0.3920477074  1.1752183407
##  [451]  1.2388874693  2.5219990684  1.0174486195  1.4705732258  0.6658114343
##  [456] -0.9892756007 -0.4708897074  0.9446191245  0.8438408036 -0.4578339012
##  [461]  1.0971901065  1.5945919147  0.5944609843  2.5489790044  0.5836657023
##  [466]  1.2304906244  0.5596958901  1.1036096722  1.2176164262  1.5302172321
##  [471]  0.5620884572 -0.1960213285  1.0871253013  1.9049245930  0.7035789060
##  [476]  0.4653967804  0.5399088625  1.0372061141  0.5420371794  0.6620645331
##  [481]  0.4680435001  1.5760236088  1.0797138976  0.6889486441  0.6065084219
##  [486]  1.1110895952  0.6350145015  0.9220230741  1.6668354350 -0.0508582050
##  [491]  0.6198437945  2.5146466718  1.1474467554  0.6037732532  1.1822970646
##  [496]  1.5996500967  1.0642519854 -0.7205074325  1.0581428386  0.9803557795
##  [501]  0.9235836137  1.5729443834  0.5722874558  0.6719660809  1.4941928741
##  [506]  1.9121249230  1.5770057168  1.0294743176  0.6767600965  1.7232251668
##  [511]  0.1198988816  0.9700246989  1.2475644655  1.2370743701  1.0035030216
##  [516]  1.1715412888  0.1905936957  1.3914033040  0.5624987547  0.2940434062
##  [521]  1.1537601670  1.2635945107  1.9667380951  1.2549373926  1.1525773121
##  [526]  1.2194882080  0.0303457829 -0.3478930926  0.7590618156 -0.1872260125
##  [531]  0.4383173911 -0.1767792850  0.4674490882  0.5845198035 -0.0064556874
##  [536]  0.3839882031  1.2954870462  1.3503516859  0.3690225789  1.5287767551
##  [541]  0.5402947473  0.5662729808  0.7984507137  1.0203259927  1.0993679413
##  [546]  0.0754291102  1.2369145321  1.0808813215  0.3172978854  0.9860267297
##  [551]  1.0840341262  2.1091672273  1.2561426739  0.9293995688  0.9280515068
##  [556]  1.5449152177  0.3482861541 -0.2305127202  0.4694697048  0.2566713015
##  [561]  0.7130383477  2.5007545580  1.2354875501  1.0851986313  0.6244957500
##  [566]  0.6845685314  1.1744749926  1.1046874033  0.2208727471  1.8625618607
##  [571]  0.9333701960  1.0140976989  1.2515883205  1.0595120888  0.5272027642
##  [576]  0.9168800795  1.6049806246  1.7668809219  0.9231246994  0.4889291976
##  [581]  1.4457059993  0.8296483026  1.5928536857  1.7377036122  0.4903208697
##  [586]  1.4954268698  1.3341474035  0.6659506694  0.6927137076  0.7887070102
##  [591]  0.2749512569  0.4225180469  0.5197248468  1.0144686735  1.2244681944
##  [596]  1.5771121861  1.1169485085 -1.0445006778  1.0735508153  0.1613822449
##  [601]  1.2484416147  1.0851651175  0.6382064312  0.5303047685  0.6057974628
##  [606] -0.5447597977  0.8138971006 -0.1947955606  1.2790917589  0.9448867949
##  [611]  0.9971117310  1.0688807562  1.5750259918  1.4559072303  0.0012532715
##  [616]  0.7080348117  0.5254089357  0.3847112566  1.1460703688  2.1273564697
##  [621]  0.5187300799  0.9221025803  1.1332186739  1.8762337080 -0.4704967858
##  [626]  1.4131686444  0.7798187330  0.8403220508 -0.6060743265 -0.8338513014
##  [631]  1.5031536918  0.6302704903  1.1599625374  1.0993436679  1.6236819851
##  [636]  1.2881328625  1.5355901651  1.5537233051  1.0301736954  0.8708426655
##  [641]  1.3712578003  0.8586426068 -0.9392552214  1.8998111396  0.6430147870
##  [646]  0.8228483461 -0.4247265793  0.8296483026  0.8926049027  0.5364528392
##  [651]  0.9565064813  1.9244650431  0.7976345452  1.2180602913 -0.4524357702
##  [656]  1.4469727929  1.3598188846  1.3555038232  0.6769906233  1.0450254610
##  [661]  1.1959042583  1.1605082817  1.6512223779 -0.4880807606  0.9948794962
##  [666]  0.9641779773  0.2750517211  1.8107975934  0.1329350572 -0.2850830981
##  [671]  1.1843097083  0.7616522063  1.5086154331  1.5486658383  1.5783293515
##  [676]  1.1243895518  0.7870885759  0.8453556913  1.1415905329 -0.0274204953
##  [681]  0.8008682543  1.0758138096  0.6931582085  0.6905661881  2.0394048369
##  [686]  1.0443395301  0.5215162281  1.1836593207  1.6253231670 -0.5451770738
##  [691]  1.0535972848  0.4194682711  0.1863113352  1.3884493663  0.1479788123
##  [696]  1.0347925337  0.8698807789  1.5076562016  0.2046046881  0.5521749654
##  [701]  0.1482892628  0.9070974028  0.0222366432 -0.2758546317  0.9487016143
##  [706]  1.3719532325  0.5763011817  1.0206317342  1.3690135316  1.0807890006
##  [711]  0.2159585215 -0.8198465629  0.6341964373  1.0116174915 -0.7845036945
##  [716]  0.9953257238  1.7023190228  1.1139002318  0.2890376928  1.3009897983
##  [721]  0.6759429576  0.8545004615  1.7948440879 -0.0032505713 -0.1207003423
##  [726]  0.9235240979  1.6843989553  1.0743943798  1.0870974420  0.9674904627
##  [731]  0.8118770707  0.4978143460  2.1676738794  0.6591635409  1.7785927897
##  [736]  0.9651429760 -0.5509159051  0.3739637760  0.2395164478  1.6405404651
##  [741]  0.4305651814  1.2358565318  1.6903496091  0.5314139949  1.2992311150
##  [746]  1.2851553107  1.4654567835  1.7891119917  1.2652292450  0.6900815762
##  [751]  1.1534530593  0.7731230601  1.8149009501  0.6434928005  1.6904191919
##  [756]  0.4769015862  1.4211666421  0.1877570609 -0.2111006513  1.9902607139
##  [761]  1.1724689157  0.6537036497  2.0694324039  0.5229380877  0.3650194451
##  [766]  1.9049245930  0.8253903961  1.0236111462  0.8560079781  0.8384350655
##  [771]  0.5854288619  0.4713605305  0.2047134961 -0.9175050363 -0.0155078481
##  [776]  1.2257451497  1.0983759603  0.2647061627  1.2021338084  0.2860457531
##  [781]  1.9008120201  0.0406103453  1.3637585671  0.9793142099  1.0023488352
##  [786]  0.7774915403  1.0358230561  1.1669900992  0.6620723414 -0.5915559418
##  [791]  0.9853517283  0.3300973543  0.5127859513  0.8130193721  0.5977179862
##  [796]  1.2112591752  0.3730384960  0.2644572305  0.0435212657  0.6691339542
##  [801]  1.9608186658  0.4717218383  0.4963511277 -0.1909196372  1.0794370119
##  [806]  1.5714369626  0.5135801602  1.2015746290 -0.0727665922  1.2729084806
##  [811]  0.2153081659 -1.5969532682  1.1601101062  1.7900691988  0.9297535937
##  [816]  1.3140704081  0.9258427323  2.5129603833  0.5401242508  1.1921004667
##  [821]  1.1410229763  0.3229695335  1.0114253406  1.3192128136 -0.5138939480
##  [826]  1.6148872990  0.1201743892  1.2663931597  0.4610218918  0.6887271454
##  [831]  0.9467742296  1.6013061325  1.7288264797  1.3559672545  1.1395241680
##  [836]  0.5195310081  1.9514052160 -0.1065328481  1.2285730140  1.5166613403
##  [841] -0.0752791219  0.7080348117 -0.0981984085  1.6330586615  0.4976851769
##  [846]  0.4189708393  1.3961473416  0.9776131262  0.4820145746  1.6594980524
##  [851]  0.0616181557  1.0990461346  1.0417950205  1.5793831634  0.2677226265
##  [856]  1.6110037450  1.1333851500  1.5739673531  0.8571397953 -0.0606905393
##  [861]  0.7596117930  0.6124809356  0.2129968175  0.6143488372  0.6780760715
##  [866]  0.9602933967  1.0441026729  0.3318209582  0.9884174629  1.5829985544
##  [871]  0.3280330864  0.9964405138  1.6904191919  1.5209576771 -0.0552491412
##  [876]  1.5292941101  1.1527023854  1.4444398877  0.9481858402  1.9568732017
##  [881]  0.5830315124  0.2125454767  0.4467454735  1.6238618117  0.6818468119
##  [886]  0.9318242935  1.2647056423  0.6705297183 -0.4807990906  1.0251946541
##  [891]  0.2920153398  0.7401326572  0.7851293484  1.1398045319  0.7411918502
##  [896]  1.5822350306  0.2492181355 -0.5002075331  0.4490405621  1.3392748258
##  [901]  0.4983270796  0.6372560292 -0.4617015214 -0.2022429069  1.0541167691
##  [906]  1.0810696236  0.4002190001  0.2347850011  0.9837928309 -0.4463544158
##  [911]  0.1643119819  0.8389607437  1.2973247978 -0.2752532196  1.0866569982
##  [916]  0.3911804724  0.5354194454  1.0522656420  0.2459669512  1.0945221337
##  [921]  1.1177949747 -0.1264730996 -1.0748319040  1.6121049332  0.9328894796
##  [926]  0.9145100064  0.7160494752  1.2474879424  1.0214920681  0.5585673862
##  [931]  2.0101590079  0.6822143855  1.1874557726  1.5857569538  0.7272732038
##  [936]  1.4523266261  1.0689856491  1.0322791029  0.2422901713  1.7490786993
##  [941]  0.7090187992  2.0022909813  1.1421071740  1.5302172321 -0.3326107688
##  [946]  0.9587126489  0.4130363605 -0.8322038574  2.0516274784  0.7233727983
##  [951]  0.5018127940  0.5033686235  1.3079988491  0.9448665176  1.0054526316
##  [956]  1.4114146945  1.0545567056  0.8700442932 -0.0827408294  0.7066232810
##  [961]  0.9012301385  0.5917288520  0.5521085697  1.5806025286  0.4851663322
##  [966]  1.1047156984  1.0523811761  0.5356615685  0.5538415883 -0.6507182195
##  [971]  0.6300020234  1.4045616824  1.4833242145  0.6259663823  0.6244957500
##  [976]  1.1865764593  1.6226060428  0.1503675509 -0.3688632920 -0.3620763576
##  [981]  1.1342411151  0.9096740752  0.6818913205  1.0094881866  1.4142546829
##  [986]  1.3057251575  0.5819476873 -0.1693064911  0.9007718722  0.4092438521
##  [991]  1.1053344581  0.5560203135  1.0442611065  1.8421673830  1.7398849774
##  [996]  1.0261621624  1.4753430683  0.5235168230  1.4488249064  1.4426890759
## 
## $func.thetastar
## NULL
## 
## $jack.boot.val
## NULL
## 
## $jack.boot.se
## NULL
## 
## $call
## bootstrap(x = original.data, nboot = 1000, theta = skewness)
```

```r
hist(results,breaks=30,col="pink",ylim=c(0,1),freq=F)
hist(results2$thetastar,breaks=30,border="purple",add=T,density=20,col="purple",freq=F)
```

<img src="Week-2-lab_files/figure-html/unnamed-chunk-23-1.png" width="672" />

What would have happened if we would have fit a normal distribution instead of a gamma distribution?


```r
fit2<-fitdistr(original.data,dnorm,start=list(mean=1,sd=1))
```

```
## Warning in densfun(x, parm[1], parm[2], ...): NaNs produced

## Warning in densfun(x, parm[1], parm[2], ...): NaNs produced

## Warning in densfun(x, parm[1], parm[2], ...): NaNs produced

## Warning in densfun(x, parm[1], parm[2], ...): NaNs produced

## Warning in densfun(x, parm[1], parm[2], ...): NaNs produced

## Warning in densfun(x, parm[1], parm[2], ...): NaNs produced
```

```r
fit2
```

```
##       mean          sd    
##   0.54699244   0.28653547 
##  (0.09061047) (0.06406753)
```

```r
results.norm<-c()
for (i in 1:1000)
  {
  x.star<-rnorm(length(original.data),mean=fit2$estimate[1],sd=fit2$estimate[2])
  results.norm<-c(results.norm,skewness(x.star))
  }
head(results.norm)
```

```
## [1]  0.383678911 -0.469041908 -0.175681961 -0.133035677 -0.296936961
## [6] -0.006750163
```

```r
hist(results,breaks=30,col="pink",ylim=c(0,1),freq=F)
hist(results.norm,breaks=30,col="lightgreen",freq=F,add=T)
hist(results2$thetastar,breaks=30,border="purple",add=T,density=20,col="purple",freq=F)
```

<img src="Week-2-lab_files/figure-html/unnamed-chunk-24-1.png" width="672" />

All three methods (two parametric and one non-parametric) really do give different distributions for the bootstrapped statistic, so the choice of which method is best depends a lot on the situation, how much data you have, and what you might already know about the underlying distribution.

Jackknifing is just as easy at bootstrapping. Here we will do a trivial example for illustration. We will write a little function for the mean even though you could put the function in directly with 'jackknife(x,mean)'


```r
theta<-function(x)
  {
  mean(x)
  }
x<-seq(0,9,by=1)
results<-jackknife(x=x,theta=theta)
results
```

```
## $jack.se
## [1] 0.9574271
## 
## $jack.bias
## [1] 0
## 
## $jack.values
##  [1] 5.000000 4.888889 4.777778 4.666667 4.555556 4.444444 4.333333 4.222222
##  [9] 4.111111 4.000000
## 
## $call
## jackknife(x = x, theta = theta)
```

**<span style="color: green;">Checkpoint #6: Why do we not have to tell the 'jackknife' function how many replicates to do?</span>**

Let's compare this with what we would have obtained from bootstrapping


```r
results2<-bootstrap(x,1000,theta)
mean(results2$thetastar)-mean(x)  #this is the bias
```

```
## [1] 0.0223
```

```r
sd(results2$thetastar)  #the standard deviation of the theta stars is the SE of the statistic (in this case, the mean)
```

```
## [1] 0.8867915
```


Everything we have done to this point used the R package 'bootstrap' - now lets compare that with the R package 'boot'. To avoid any confusion (a.k.a. masking) between the two packages, I recommend detaching the bootstrap package from the workspace with


```r
detach("package:bootstrap")
```


The 'boot' package is now recommended over the 'bootstrap' package, but they give the same answers and to some extent it is personal preference which one prefers to use.

We will still use the mean as the statistic of interest, but we will have to write a new function for it because the syntax of the 'boot' package is slightly different:


```r
library(boot)
theta<-function(x,index)
  {
  mean(x[index])
  }
boot(x,theta,R=999)
```

```
## 
## ORDINARY NONPARAMETRIC BOOTSTRAP
## 
## 
## Call:
## boot(data = x, statistic = theta, R = 999)
## 
## 
## Bootstrap Statistics :
##     original      bias    std. error
## t1*      4.5 -0.02882883   0.9375007
```

One of the main advantages to the 'boot' package over the 'bootstrap' package is the nicer formatting of the output.

Going back to our original code, lets see how we could reproduce all of these numbers:


```r
table(sample(x,size=length(x),replace=T))
```

```
## 
## 1 2 5 7 8 
## 2 3 1 1 3
```

```r
xmeans<-vector(length=1000)
for (i in 1:1000)
  {
  xmeans[i]<-mean(sample(x,replace=T))
  }
mean(x)
```

```
## [1] 4.5
```

```r
bias<-mean(xmeans)-mean(x)
se.boot<-sd(xmeans)
bias
```

```
## [1] 0.0219
```

```r
se.boot
```

```
## [1] 0.8872416
```

Why do our numbers not agree exactly with those of the boot package? This is because our estimates of bias and standard error are just estimates, and they carry with them their own uncertainties. That is one of the reasons we might bother doing jackknife-after-bootstrap.

The 'boot' package has a LOT of functionality. If we have time, we will come back to some of these more complex functions later in the semester as we cover topics like regression and glm.

