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
## 0 1 2 4 5 6 
## 2 2 2 1 2 1
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
## [1] -0.0243
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
## [1] 2.691122
```

```r
UL.boot
```

```
## [1] 6.260278
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
##    [1] 4.2 3.9 4.4 2.9 2.5 2.7 4.5 4.3 4.3 5.5 3.1 3.8 4.0 6.4 4.7 5.3 5.2 3.2
##   [19] 4.8 4.7 3.9 4.5 6.4 4.1 5.3 4.1 5.6 2.8 4.2 4.5 5.7 4.5 5.7 4.8 4.5 5.2
##   [37] 4.6 5.4 2.4 5.4 4.9 5.8 5.3 6.3 4.4 4.1 4.8 3.7 5.3 4.3 4.3 2.7 3.5 5.3
##   [55] 2.8 3.0 5.7 3.1 3.5 4.5 4.1 3.6 3.5 3.6 3.3 3.5 5.1 5.5 4.5 4.0 4.0 3.7
##   [73] 6.1 5.4 4.5 4.3 4.8 3.4 6.2 3.8 4.8 3.8 5.6 4.7 3.0 3.7 5.2 4.2 5.0 4.6
##   [91] 5.0 4.5 3.8 5.4 4.1 3.1 4.4 5.1 4.4 7.5 6.1 2.8 4.5 3.7 3.0 4.4 4.4 4.5
##  [109] 3.7 6.5 6.6 4.1 4.2 3.5 5.2 4.0 5.5 4.9 4.0 4.1 4.3 4.1 3.8 4.0 3.4 6.6
##  [127] 5.2 6.0 6.0 3.5 4.6 3.8 3.4 5.9 4.8 5.2 6.2 3.7 3.1 4.0 4.9 4.4 3.6 3.4
##  [145] 5.0 4.2 6.4 5.3 3.6 3.4 4.4 4.5 3.6 5.2 7.0 4.8 5.5 4.8 4.3 6.5 5.3 4.5
##  [163] 4.4 3.7 4.7 4.8 4.2 4.7 3.1 3.1 5.0 4.7 5.3 4.8 4.3 3.8 4.1 5.0 3.9 3.1
##  [181] 3.8 4.5 4.3 5.5 4.6 3.5 5.2 3.3 4.3 4.3 3.7 4.2 5.0 4.8 4.7 6.1 4.3 2.7
##  [199] 3.0 3.6 5.1 4.4 4.4 3.9 2.9 4.1 4.2 5.0 3.3 7.0 4.0 5.7 3.5 5.0 4.4 5.1
##  [217] 4.9 4.1 5.6 5.0 4.2 3.4 5.9 6.1 5.1 4.7 3.5 4.6 4.6 3.4 5.3 5.1 3.3 5.2
##  [235] 4.9 5.0 3.1 4.9 3.9 5.1 3.5 4.8 4.5 4.5 4.6 5.3 2.8 5.0 3.6 2.5 6.8 5.0
##  [253] 3.3 6.1 3.0 3.7 4.1 4.3 5.9 2.3 5.6 5.5 5.2 5.8 3.3 4.3 4.9 4.4 5.4 3.4
##  [271] 5.2 4.6 5.2 4.5 3.0 4.0 5.7 3.7 5.7 4.9 5.2 4.4 4.1 4.2 3.8 3.8 6.4 4.1
##  [289] 5.6 3.9 5.2 3.7 3.2 3.9 3.6 6.2 6.1 5.2 4.3 4.2 5.7 5.6 3.8 5.0 4.8 6.4
##  [307] 4.6 4.7 5.1 4.1 3.9 4.2 3.3 4.9 5.3 3.2 4.1 4.8 4.5 5.7 4.4 4.1 3.6 4.4
##  [325] 4.7 3.1 4.1 4.4 5.2 4.3 4.4 3.5 4.3 4.7 4.1 3.7 4.4 6.2 6.0 4.1 4.9 4.9
##  [343] 5.7 5.7 5.6 3.2 4.2 4.2 5.4 3.6 5.4 1.9 6.2 4.7 4.5 5.6 5.7 3.8 3.9 5.2
##  [361] 5.0 5.8 6.0 3.3 5.6 5.0 4.6 4.9 5.8 3.6 5.3 4.4 4.4 3.9 6.0 2.9 5.3 4.1
##  [379] 6.0 4.0 3.6 3.5 3.9 2.4 4.7 4.6 2.8 6.8 4.8 3.6 5.2 2.8 4.6 3.8 6.5 3.6
##  [397] 3.2 3.1 4.3 5.2 5.0 3.7 4.1 4.6 3.6 6.0 4.5 5.2 2.5 3.3 4.3 6.1 5.1 5.0
##  [415] 4.2 3.7 5.8 5.3 3.6 4.0 4.9 5.9 3.0 4.9 5.4 6.1 3.6 4.6 4.5 6.2 3.8 5.2
##  [433] 3.7 4.0 2.9 4.5 6.8 2.9 4.6 5.5 3.6 3.5 4.0 6.0 4.2 2.9 5.5 5.2 2.2 4.6
##  [451] 3.7 3.8 4.8 3.8 5.0 3.9 4.3 5.4 4.9 3.4 5.3 4.2 5.1 4.5 4.5 4.4 3.7 5.3
##  [469] 3.0 4.0 5.0 5.4 4.3 4.8 4.2 5.1 3.9 5.5 4.1 2.6 4.7 4.8 2.9 3.8 4.5 5.0
##  [487] 5.5 2.5 2.8 5.4 4.3 5.1 5.1 4.2 3.4 5.8 6.1 5.5 4.5 2.6 3.7 5.0 4.1 5.1
##  [505] 5.6 4.9 5.1 5.3 3.6 3.6 5.2 3.3 5.1 3.4 5.1 5.1 4.4 5.3 3.5 5.7 6.0 5.1
##  [523] 4.8 3.1 4.6 4.1 4.3 4.2 3.8 3.9 5.3 7.0 3.9 4.3 4.2 6.1 5.6 4.1 4.1 5.0
##  [541] 3.2 4.2 6.0 4.7 5.8 5.2 3.8 3.2 3.2 5.5 4.4 5.9 4.5 3.8 2.7 6.6 5.2 5.5
##  [559] 3.0 6.1 6.0 2.6 4.8 5.5 6.4 3.3 4.1 5.1 3.9 6.0 4.8 3.8 6.4 4.4 5.2 6.2
##  [577] 4.7 5.7 5.3 5.1 5.3 3.7 4.0 3.6 3.6 2.9 4.8 5.6 4.6 4.2 4.4 4.2 4.7 5.8
##  [595] 4.3 4.4 4.8 4.2 3.8 5.5 4.7 5.1 4.9 5.0 3.5 3.1 4.4 2.7 3.5 4.5 4.1 3.6
##  [613] 4.0 4.3 5.3 4.8 3.8 4.3 3.2 3.8 3.5 5.2 5.6 5.3 4.2 4.7 5.1 4.3 5.1 6.1
##  [631] 5.7 3.8 5.1 3.8 4.7 5.0 4.4 5.1 5.8 3.9 3.2 5.5 4.6 4.8 4.0 3.6 6.3 4.6
##  [649] 5.3 4.4 4.7 3.3 4.5 5.4 3.5 3.9 3.5 4.8 5.2 5.0 4.4 4.3 4.3 5.0 2.7 5.5
##  [667] 3.5 2.7 4.8 6.6 3.7 5.3 6.1 2.8 3.1 4.2 5.7 2.6 5.0 4.4 4.9 5.1 4.2 4.0
##  [685] 4.8 6.1 4.0 5.6 4.5 5.0 5.2 3.4 4.0 4.5 3.7 4.5 3.7 5.0 5.8 5.0 4.6 4.2
##  [703] 4.5 4.0 5.4 4.6 3.9 4.6 5.7 4.3 6.9 4.6 3.1 4.6 4.6 3.7 4.8 6.6 4.7 4.9
##  [721] 4.0 5.2 3.6 2.9 5.3 4.6 5.9 3.1 3.9 4.1 5.5 5.0 3.6 4.2 4.7 4.5 4.8 4.2
##  [739] 4.7 5.3 5.9 4.4 2.9 4.9 4.0 4.3 3.6 5.5 4.6 6.1 4.0 5.5 3.2 4.7 5.4 4.2
##  [757] 5.3 3.3 6.1 4.4 5.6 6.2 5.4 5.5 5.3 4.2 2.9 3.2 4.2 3.0 3.5 6.4 4.3 3.7
##  [775] 4.1 5.5 5.6 2.2 5.1 3.3 5.4 4.2 4.3 5.1 5.6 4.8 4.1 4.9 4.8 4.3 2.7 5.1
##  [793] 5.3 5.5 6.0 4.7 5.7 3.9 3.7 4.7 5.2 5.2 3.1 4.7 5.4 5.5 5.4 4.7 5.0 3.6
##  [811] 5.2 4.1 4.0 4.5 4.6 3.1 5.3 3.4 4.1 4.1 5.1 5.9 4.9 4.9 4.2 5.1 4.2 3.7
##  [829] 4.9 5.2 6.0 5.0 5.4 6.6 4.3 4.8 4.2 3.5 3.5 5.4 5.3 4.0 4.9 4.8 4.5 3.7
##  [847] 5.4 6.1 4.0 5.3 3.2 4.4 4.9 4.3 4.5 4.5 4.6 5.4 4.6 5.4 4.3 4.9 5.8 4.4
##  [865] 4.1 2.2 4.5 3.9 3.2 4.0 4.4 4.2 4.4 5.1 4.5 3.5 5.4 5.3 4.0 5.5 5.5 5.9
##  [883] 5.6 3.7 5.1 2.3 4.5 4.8 3.2 4.0 4.6 6.5 4.1 4.1 3.2 4.7 4.3 6.2 3.4 3.8
##  [901] 4.2 4.3 5.7 4.5 6.2 3.9 4.7 4.8 4.9 4.1 3.4 5.3 5.0 4.5 4.1 4.1 3.8 4.7
##  [919] 5.1 4.3 6.7 5.1 4.9 4.9 5.0 3.9 3.3 5.0 5.6 3.6 4.4 4.1 3.6 5.8 4.1 3.5
##  [937] 4.7 3.4 4.4 4.6 5.2 3.9 3.4 5.7 5.2 4.7 3.1 4.5 5.0 5.8 3.3 5.5 2.8 6.5
##  [955] 4.4 2.8 4.4 4.7 4.7 4.0 5.3 3.5 6.2 4.4 4.0 3.6 4.4 4.1 5.2 2.5 4.6 4.5
##  [973] 4.1 4.6 5.5 3.5 3.8 5.7 4.8 4.5 6.8 4.3 4.7 4.7 5.4 4.0 4.1 3.0 3.6 2.3
##  [991] 4.4 6.4 3.2 4.3 4.4 5.6 4.0 5.0 4.4 5.6
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
##   2.7   6.4
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
##    [1] 5.2 5.5 5.5 4.4 4.3 4.8 3.8 5.1 3.9 4.0 3.4 5.4 5.6 5.2 2.6 6.0 3.6 4.3
##   [19] 3.7 4.7 4.3 4.0 5.1 4.2 3.7 2.8 5.7 3.9 5.4 4.8 4.3 4.0 4.0 3.9 5.6 5.4
##   [37] 3.3 4.8 5.2 2.6 3.8 4.7 3.8 2.9 2.4 4.5 5.1 5.7 5.0 6.1 5.8 3.9 3.6 4.3
##   [55] 4.3 3.7 5.3 5.3 4.6 2.6 6.0 5.3 6.0 4.6 5.6 5.5 5.9 4.9 6.2 4.7 5.6 4.9
##   [73] 3.8 4.7 4.3 4.2 4.0 5.0 4.9 2.7 3.8 4.6 4.7 5.8 5.8 5.2 2.8 4.6 4.0 5.6
##   [91] 4.3 5.6 4.4 4.4 3.8 4.5 5.4 5.4 4.7 4.0 3.6 5.8 5.4 3.5 4.1 5.0 4.7 5.6
##  [109] 4.6 4.2 4.7 5.2 5.5 2.8 3.9 4.6 5.9 5.0 3.4 3.7 4.1 5.6 4.7 4.5 4.4 4.1
##  [127] 5.0 4.4 4.9 4.5 3.6 3.4 4.0 3.9 4.0 3.1 5.0 5.5 4.3 3.5 3.7 4.8 5.7 4.8
##  [145] 5.8 5.9 5.4 3.7 4.9 4.7 3.2 6.3 3.4 4.9 5.1 3.7 4.8 5.3 4.3 5.3 4.3 7.3
##  [163] 5.1 4.2 3.8 5.6 4.7 3.7 3.7 3.3 3.6 2.7 6.1 4.5 4.2 3.6 4.7 5.5 5.2 3.6
##  [181] 4.4 3.6 4.0 3.7 3.1 3.7 4.4 5.8 5.1 5.1 3.8 5.5 4.5 2.6 4.5 4.7 6.3 4.8
##  [199] 4.3 3.7 5.1 5.4 4.5 3.5 5.2 3.8 3.8 5.0 5.2 6.0 3.4 3.9 4.7 4.2 4.0 4.5
##  [217] 5.7 3.8 5.7 4.9 3.4 4.6 6.6 4.0 5.2 4.8 6.3 3.4 3.4 4.9 3.8 5.2 3.7 4.2
##  [235] 6.3 4.2 5.8 5.3 4.7 5.7 5.1 5.3 3.2 5.1 6.2 3.7 5.6 2.8 2.8 6.5 4.5 4.2
##  [253] 5.1 4.8 5.5 4.2 4.4 4.5 5.3 4.3 3.3 4.8 5.8 4.7 5.3 5.4 3.9 5.6 4.4 5.1
##  [271] 4.2 4.9 4.3 5.8 4.8 5.9 3.3 3.6 4.7 5.0 4.7 5.3 5.1 3.8 4.0 4.5 2.7 3.8
##  [289] 2.9 3.8 4.8 5.1 5.1 3.9 3.6 3.5 3.7 5.5 5.1 4.9 4.8 6.2 4.8 5.3 4.6 5.2
##  [307] 5.9 5.7 3.0 5.2 4.9 4.1 4.5 4.0 5.3 4.5 3.8 3.6 4.4 4.6 4.3 4.3 4.0 5.2
##  [325] 5.1 6.7 3.0 3.5 4.5 5.5 3.8 4.1 4.7 4.4 4.7 3.9 4.3 4.1 3.4 5.5 3.4 5.6
##  [343] 4.0 4.6 5.0 2.9 4.1 4.2 5.5 4.0 5.2 3.2 5.1 5.0 4.0 4.6 4.9 3.0 4.1 5.7
##  [361] 5.3 5.9 4.7 3.4 4.3 4.3 4.4 4.5 3.9 4.2 5.3 4.1 3.2 3.9 3.8 4.1 5.0 3.6
##  [379] 4.0 3.6 3.9 5.0 4.1 5.6 4.0 5.2 4.6 3.1 4.7 4.0 4.5 3.4 2.9 3.6 4.8 3.7
##  [397] 3.8 5.7 3.7 4.3 5.4 4.5 3.1 5.8 4.9 3.6 6.1 5.2 6.9 5.2 4.7 6.0 4.6 4.4
##  [415] 5.6 4.7 3.2 3.3 5.6 3.4 5.2 4.2 5.0 4.4 4.2 4.4 5.4 4.1 4.8 4.6 5.0 4.4
##  [433] 5.9 4.8 3.4 4.7 4.8 4.2 5.0 5.7 5.3 3.2 3.4 3.5 5.2 6.0 5.5 5.6 5.2 6.2
##  [451] 3.8 5.4 5.8 4.7 3.0 4.3 4.1 4.6 4.6 4.8 2.6 4.6 4.5 4.2 4.0 4.7 4.7 5.1
##  [469] 3.0 3.7 5.6 5.3 6.7 5.7 5.0 5.5 4.1 4.1 4.6 5.1 5.5 3.9 3.9 3.7 2.6 2.9
##  [487] 4.4 4.5 5.7 4.3 2.5 4.4 6.2 4.8 3.0 6.0 4.4 4.4 2.6 4.1 3.2 5.1 4.5 5.7
##  [505] 3.5 3.7 4.6 4.4 5.0 5.5 5.4 4.7 5.1 5.0 5.3 5.2 4.4 6.2 3.2 5.0 4.5 3.6
##  [523] 2.6 5.4 4.7 4.2 4.1 3.7 6.0 2.2 4.4 2.8 4.4 5.1 4.4 3.9 2.9 4.8 6.3 4.2
##  [541] 4.7 3.0 3.9 3.5 3.2 4.7 3.9 3.5 2.7 5.0 4.0 4.1 3.5 4.7 4.6 2.4 5.4 4.8
##  [559] 4.5 4.2 4.0 4.2 3.8 5.5 7.0 3.3 4.2 4.4 1.7 3.9 6.0 5.2 6.2 4.4 4.5 3.0
##  [577] 4.7 5.3 4.5 2.9 5.2 4.2 4.6 5.5 2.4 4.8 4.5 5.5 6.4 4.0 5.3 4.9 4.2 7.2
##  [595] 3.5 4.2 4.6 3.4 4.8 5.6 4.4 3.7 4.1 4.6 4.6 4.2 3.8 4.1 3.4 3.6 3.2 4.7
##  [613] 3.7 3.8 5.7 4.8 6.7 4.6 4.2 3.8 2.9 4.2 6.0 2.6 4.7 5.2 6.0 4.1 2.2 4.8
##  [631] 4.5 3.4 3.3 4.4 5.6 3.7 4.0 4.4 4.8 3.9 4.1 4.3 3.6 4.4 3.9 5.4 5.3 3.9
##  [649] 5.2 4.7 4.3 2.7 3.6 5.2 5.6 4.7 4.3 2.9 5.2 3.7 6.0 5.5 4.2 3.4 5.2 4.7
##  [667] 2.5 3.0 3.5 4.9 4.7 5.2 4.8 5.6 3.3 5.9 3.4 4.1 3.4 3.6 4.7 4.9 4.0 5.2
##  [685] 4.1 5.4 5.9 3.8 4.7 5.5 4.6 5.8 5.9 2.4 4.9 4.3 4.1 4.6 3.7 5.5 5.2 2.9
##  [703] 5.4 3.2 2.1 3.4 4.8 3.7 5.3 4.7 5.0 3.6 3.7 5.5 4.9 5.6 3.6 5.8 6.0 2.9
##  [721] 2.4 5.4 6.5 3.6 4.8 2.7 4.8 4.8 4.2 4.7 4.8 4.5 2.9 5.3 3.5 4.6 4.6 2.8
##  [739] 4.5 5.2 4.5 5.8 5.0 4.2 4.2 7.1 4.5 4.8 3.3 5.4 3.3 4.9 5.9 4.4 3.6 4.7
##  [757] 3.5 3.2 5.1 3.9 3.3 3.9 5.6 4.3 3.9 5.3 3.3 5.1 3.5 4.8 3.1 5.3 5.4 2.8
##  [775] 3.6 4.4 4.1 5.3 5.8 5.2 5.5 3.5 3.8 4.9 4.2 6.1 4.3 3.5 2.9 5.1 4.2 4.3
##  [793] 5.5 3.7 5.2 4.2 5.3 5.5 3.9 5.4 5.6 4.6 5.5 5.5 5.9 5.0 2.3 5.3 5.8 4.1
##  [811] 3.9 4.0 5.6 2.7 6.3 5.1 3.8 4.6 3.1 4.9 4.9 4.0 3.8 5.1 5.4 5.6 5.6 2.9
##  [829] 3.6 5.2 2.0 4.2 4.2 5.0 4.2 4.1 5.7 4.4 3.9 4.5 4.6 4.4 4.8 4.3 4.0 4.2
##  [847] 4.8 3.2 3.6 3.9 4.1 5.0 4.2 3.4 5.2 4.8 5.8 3.6 5.1 4.6 6.2 3.2 5.6 3.4
##  [865] 3.8 4.1 5.0 2.4 3.6 5.3 4.0 6.3 4.2 4.1 4.7 5.4 6.3 4.5 5.2 3.8 3.8 4.9
##  [883] 4.9 4.6 5.7 3.9 5.1 3.0 4.8 4.1 4.5 4.6 6.1 5.1 4.9 3.5 4.9 4.6 5.8 5.2
##  [901] 5.7 3.6 3.8 4.2 5.0 2.5 3.9 5.3 4.4 5.6 4.6 3.8 5.1 2.6 4.5 4.4 5.5 5.3
##  [919] 4.3 5.5 4.8 3.8 3.8 6.2 2.9 2.9 4.7 4.3 5.4 4.2 3.8 3.3 4.3 5.0 5.6 4.1
##  [937] 3.6 4.0 4.5 5.0 3.7 4.3 2.7 3.3 3.9 4.2 3.9 5.0 2.9 4.4 5.1 3.4 4.3 3.7
##  [955] 5.1 5.1 2.0 5.0 5.8 4.9 5.2 4.5 4.5 4.4 4.1 4.4 6.8 4.2 4.3 5.5 3.9 4.3
##  [973] 4.1 4.4 3.2 4.3 4.2 3.8 4.7 4.8 3.7 5.5 4.8 3.2 4.6 3.7 6.6 3.9 3.7 5.5
##  [991] 4.9 4.9 4.3 4.8 6.1 5.5 4.2 4.1 4.0 4.0
## 
## $func.thetastar
## [1] -0.0078
## 
## $jack.boot.val
##  [1]  0.49241983  0.36581197  0.35451807  0.13653333  0.03011696 -0.02660819
##  [7] -0.10395480 -0.41142857 -0.37454068 -0.50596591
## 
## $jack.boot.se
## [1] 0.9942688
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
##    [1] 5.4 4.7 3.9 3.2 4.8 4.1 3.3 5.1 4.7 1.7 4.5 4.8 5.8 5.5 2.5 4.8 4.1 4.6
##   [19] 3.9 4.6 5.0 4.0 3.2 2.9 4.4 4.2 4.5 6.1 4.5 3.4 4.8 3.6 4.3 4.0 5.8 6.0
##   [37] 3.0 3.6 5.6 5.7 3.6 5.0 4.5 5.1 4.9 5.6 5.6 3.2 5.3 5.0 3.4 6.0 3.3 5.9
##   [55] 4.7 2.9 4.9 3.2 5.0 5.2 4.6 4.3 5.4 3.8 3.3 4.7 4.1 5.8 3.6 4.0 3.4 3.9
##   [73] 4.9 5.3 5.3 3.9 3.8 4.2 4.9 4.4 3.6 6.1 5.1 3.8 5.2 4.0 3.4 4.4 5.7 5.7
##   [91] 5.7 4.1 4.4 4.8 4.0 6.0 3.6 4.4 3.3 5.4 4.5 3.9 4.0 5.0 3.7 5.6 5.2 4.4
##  [109] 3.5 3.7 4.2 4.8 6.4 4.5 4.5 4.6 4.8 2.9 5.5 5.8 4.2 4.0 3.9 3.5 4.0 4.5
##  [127] 3.9 4.1 3.1 5.6 5.5 6.4 5.3 6.0 4.6 4.9 3.9 3.1 5.3 5.3 3.9 4.2 2.8 5.0
##  [145] 4.7 4.2 5.2 5.3 5.3 4.5 4.4 4.4 4.2 4.1 4.9 4.2 3.1 6.3 4.2 5.5 4.1 5.2
##  [163] 3.1 5.3 4.4 4.9 2.4 4.1 4.1 4.6 5.5 2.7 5.9 6.0 4.6 5.5 4.2 3.5 4.1 5.8
##  [181] 4.7 4.5 4.6 4.9 3.0 5.0 4.8 3.8 5.7 4.0 2.4 3.2 4.6 3.9 4.4 2.8 4.7 3.4
##  [199] 3.3 4.8 4.3 3.8 5.1 5.6 4.9 5.3 4.1 5.7 4.5 6.0 4.2 3.8 4.3 4.2 3.4 4.5
##  [217] 3.3 5.3 5.5 4.8 2.4 2.3 5.5 5.4 4.8 4.0 4.5 5.3 5.0 5.2 3.0 4.5 4.2 5.5
##  [235] 5.4 3.8 5.7 5.1 5.3 4.3 5.5 4.6 5.2 4.3 5.4 4.6 3.6 4.9 4.3 5.9 3.8 4.8
##  [253] 5.0 4.8 4.7 3.1 5.2 4.6 4.0 3.7 3.7 5.7 3.6 3.7 5.5 4.8 4.0 4.0 5.0 4.9
##  [271] 4.1 4.1 4.0 4.9 3.6 4.7 3.8 3.5 3.8 5.4 5.3 4.1 5.1 5.1 4.9 3.2 3.9 4.4
##  [289] 3.9 4.1 4.2 3.1 5.5 3.0 5.2 4.3 4.6 4.6 4.5 5.2 5.2 4.3 3.5 4.4 4.3 3.2
##  [307] 3.6 4.8 3.8 6.2 3.1 5.0 3.9 5.1 4.6 4.4 3.7 4.1 4.4 5.6 3.3 5.5 5.6 4.6
##  [325] 3.5 5.7 5.2 4.5 5.1 4.3 5.2 5.1 6.0 3.7 4.6 5.3 3.2 4.3 5.6 4.6 4.2 3.8
##  [343] 4.8 5.2 4.2 4.0 3.9 3.7 4.4 4.1 3.9 4.3 5.2 3.8 5.5 5.9 4.1 4.9 4.6 2.7
##  [361] 4.2 5.1 3.1 4.9 4.3 4.4 4.6 2.2 6.0 4.5 4.9 5.0 3.5 4.0 6.2 4.0 6.0 5.8
##  [379] 3.0 5.5 4.0 5.3 2.2 4.7 5.4 5.5 6.1 3.9 4.8 5.7 4.7 4.5 6.2 4.6 3.6 5.2
##  [397] 4.6 4.6 3.5 4.9 3.6 3.5 4.5 5.4 3.9 4.0 5.3 2.8 5.2 4.0 5.4 2.6 4.1 4.6
##  [415] 4.8 4.7 4.1 4.4 4.5 4.0 2.6 3.6 3.9 4.4 4.8 4.5 4.6 4.8 4.2 5.9 3.7 5.3
##  [433] 6.0 4.7 3.5 4.3 2.9 4.5 4.6 5.5 4.5 4.4 3.9 3.9 3.1 4.1 3.9 5.1 2.3 5.7
##  [451] 3.7 3.9 6.0 3.3 4.0 5.5 5.2 4.7 3.5 3.6 4.7 5.2 5.0 3.7 5.5 5.6 3.8 3.0
##  [469] 3.3 4.4 2.5 3.6 6.1 4.2 5.1 4.4 3.6 6.0 4.9 6.5 6.4 4.8 5.9 5.0 4.9 4.4
##  [487] 3.5 6.6 4.3 3.9 4.9 3.3 2.7 4.8 5.3 4.8 4.8 3.9 6.3 3.8 3.9 5.7 3.5 3.6
##  [505] 3.9 4.0 4.5 4.8 4.6 5.3 5.0 4.7 5.3 4.1 5.0 4.8 5.9 4.9 3.2 5.4 3.3 4.1
##  [523] 4.8 4.4 4.9 4.4 4.1 3.9 4.5 5.2 5.2 4.5 6.1 4.4 3.9 4.3 3.1 5.6 5.2 4.7
##  [541] 4.6 3.7 5.4 4.0 4.0 3.5 5.5 2.9 4.7 4.9 5.3 5.3 4.4 4.6 4.8 4.2 4.1 4.4
##  [559] 3.9 4.3 3.6 4.7 5.5 4.6 3.8 4.6 4.7 4.9 4.8 4.1 2.7 5.6 4.3 5.3 5.1 5.2
##  [577] 5.3 4.6 5.3 4.6 2.9 2.9 5.6 3.1 4.7 5.0 3.2 4.0 5.9 4.8 4.2 3.7 3.0 5.1
##  [595] 4.2 5.5 3.9 5.9 4.8 2.7 5.8 3.0 6.0 4.7 4.4 4.8 3.9 4.2 4.3 2.2 3.6 4.8
##  [613] 3.7 4.4 4.6 5.1 3.7 4.2 6.2 5.6 4.9 4.5 3.6 5.0 4.1 6.2 5.3 3.4 2.7 4.3
##  [631] 2.7 5.2 4.6 6.2 3.3 4.4 5.9 4.6 4.1 4.6 3.8 5.8 4.8 4.5 2.9 5.4 4.3 4.2
##  [649] 5.0 5.8 4.1 4.0 4.4 3.1 5.6 3.9 6.0 3.6 5.6 5.2 2.2 3.3 5.3 4.9 4.7 5.2
##  [667] 5.8 2.3 6.0 4.5 5.6 3.2 4.8 5.6 3.9 4.1 4.5 4.5 6.3 5.5 4.6 3.9 5.3 3.8
##  [685] 6.1 4.8 7.0 5.8 5.0 5.2 4.3 4.3 4.7 3.4 5.4 4.6 4.4 4.1 4.8 4.8 6.0 5.6
##  [703] 2.2 3.4 4.5 4.9 5.2 3.9 4.5 7.1 4.1 4.9 2.8 4.4 4.1 5.4 3.8 4.6 5.1 5.4
##  [721] 6.1 4.7 3.6 3.9 4.2 4.0 2.5 4.3 4.1 3.2 4.9 3.2 3.1 4.4 3.9 3.5 3.7 5.3
##  [739] 4.0 4.4 5.3 4.8 3.8 4.2 5.1 3.6 4.9 5.7 6.4 4.1 3.8 4.6 5.0 4.8 4.6 3.3
##  [757] 3.4 4.5 4.0 5.1 4.9 6.2 5.8 4.0 3.4 4.3 2.9 4.8 5.3 5.2 6.0 4.4 3.2 4.3
##  [775] 3.0 4.1 4.9 2.7 4.5 3.0 6.2 4.8 4.1 5.1 4.8 4.8 5.1 4.5 4.2 3.0 3.7 4.2
##  [793] 5.3 3.2 5.7 4.1 4.4 5.3 3.6 3.3 3.6 4.6 5.6 6.0 4.3 5.5 3.6 4.9 4.1 2.8
##  [811] 3.8 3.6 2.9 3.6 4.4 4.1 5.7 4.4 3.4 5.0 4.9 3.8 4.8 5.1 5.0 3.8 5.5 5.9
##  [829] 4.1 5.0 4.5 5.1 4.8 3.6 4.4 3.6 4.5 5.1 5.6 5.4 3.2 6.4 3.6 4.5 5.3 4.0
##  [847] 5.9 4.5 3.6 5.2 3.0 3.5 4.1 6.1 3.9 4.6 5.2 4.4 4.1 4.1 3.9 3.9 4.1 6.3
##  [865] 4.7 3.7 4.0 3.1 4.0 3.7 4.2 4.6 5.4 2.5 5.2 5.2 3.1 3.0 4.7 3.9 6.5 5.1
##  [883] 4.2 5.2 5.4 4.4 4.2 5.1 3.3 5.1 4.8 3.2 3.7 5.1 4.9 5.1 5.0 3.3 5.0 2.6
##  [901] 5.5 6.5 4.6 4.8 3.4 4.2 5.8 3.7 4.8 5.4 2.9 4.2 4.3 4.6 4.8 4.3 5.5 5.3
##  [919] 5.8 3.1 4.3 4.4 6.5 3.7 5.1 4.9 4.1 6.4 5.1 3.7 5.1 4.5 4.4 3.3 3.7 3.5
##  [937] 4.5 5.1 2.8 2.8 5.3 4.6 5.2 5.2 4.5 4.6 3.1 3.3 4.4 5.3 6.4 4.7 3.4 4.1
##  [955] 4.8 3.1 4.0 3.2 5.1 4.7 4.9 2.7 5.6 3.6 3.7 3.4 4.8 5.9 4.7 5.2 4.2 4.4
##  [973] 5.1 4.9 4.6 3.0 3.2 3.8 4.4 3.7 3.8 4.2 4.3 3.9 5.5 5.3 4.9 4.1 5.2 5.3
##  [991] 4.4 3.6 4.6 3.7 3.3 5.8 3.8 5.6 5.5 4.1
## 
## $func.thetastar
##   72% 
## 5.028 
## 
## $jack.boot.val
##  [1] 5.4 5.4 5.3 5.2 5.2 5.1 4.8 4.8 4.6 4.5
## 
## $jack.boot.se
## [1] 0.9396276
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
## [1] 0.7974202
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
##   3.630107   6.991499 
##  (1.554610) (3.211142)
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
## [1]  0.955872849 -0.256703950  1.659703265 -0.694226178  0.009640977
## [6]  0.433606382
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
##    [1]  1.814576322  0.477980710  0.572126497  1.444639968 -0.486507536
##    [6]  0.182709374  0.858830136  1.473220409  1.214658053  1.360133983
##   [11]  0.463816941  0.950411688  0.572093179  0.091664199  0.398676856
##   [16]  1.400130191  1.392268376  0.302826544  0.933539603  0.774867305
##   [21]  1.473275478  0.521880116  0.825842855  2.302384393  0.635197526
##   [26]  0.571728455  0.915912780 -0.037844170  0.510048088 -0.934314954
##   [31]  0.486544487  1.508009343  0.780202210 -1.019850273  0.793148617
##   [36]  0.834754394  1.300003678  2.425257158  2.313642935  1.052748746
##   [41]  0.891032611 -0.046739374  0.982785698  1.059723207  2.020625097
##   [46]  0.513969075  0.776154704  1.271359119  0.936474366 -0.025834627
##   [51] -0.241450460 -0.006397086  0.687547296  0.181092021  1.327543399
##   [56]  0.246880156  0.850199095  0.935034953  0.444314793  0.774033644
##   [61]  0.487944945  1.751926498  0.572464194  1.177391142 -1.493227869
##   [66]  2.084943962 -0.766546590  1.359200681  1.265485968  0.007763855
##   [71]  1.025577942  1.462806178 -0.175985056  0.863160113  0.466534974
##   [76]  0.230945072  0.144692511  0.160829064 -0.808138763  0.725969409
##   [81]  0.937641039  0.948717623 -0.023764178  0.810799254 -1.009441987
##   [86]  2.138338522  0.584353735  1.187248334  1.440588117  1.680895545
##   [91] -0.015909756  0.280272831  0.227133334  0.759493249  1.108731150
##   [96]  0.918634726  0.916504409  0.907207094  0.695138722  0.054497459
##  [101]  0.222238621  0.268072629  1.367785640  0.295450146  1.049020349
##  [106] -0.524408162  0.954187854  0.863259602  0.743796328  0.258251190
##  [111]  0.446879197  2.082669681  0.763539372  0.147537167  0.856273036
##  [116]  0.931335264  0.030958051  0.675008498  0.080286809  0.021042365
##  [121] -1.596250362  0.125005361  0.093624711  1.392581779  0.291731589
##  [126]  0.479733386  0.681002726  0.763167955  0.429780056  2.174230730
##  [131]  1.547205082  0.623758361  1.319073564  1.030926431  0.809181696
##  [136]  1.374472988  0.779073309  0.009653690  1.383006354  0.447284338
##  [141]  0.470668696 -1.076527748 -1.389381954  0.437027021  0.044479793
##  [146]  0.797647206  0.451823370  0.744858236  1.352746170  0.286471545
##  [151]  0.868125423  1.364534697  1.245024343  0.395297576  0.958761371
##  [156] -0.834757244  1.009651991  1.279229423  1.630768785  0.133225245
##  [161] -0.052108524  0.388211402  0.961362895 -0.101662233 -0.422153815
##  [166]  0.582107876  0.720663317  1.386672915  0.484022970  1.388994956
##  [171]  1.384080422  1.940094224  0.291316943  1.080092097  0.559841173
##  [176]  0.332411180  2.070689934  0.278079040  0.718967345  1.674207329
##  [181]  0.583108295  0.840416493 -1.283553604  0.486514207 -0.322011497
##  [186]  0.903430115  0.645530981  0.715756858  1.446896256  0.618898413
##  [191]  0.552884333  1.605007548 -1.591760996  1.235855607  0.411793268
##  [196]  1.474015651  0.647366490  2.106665428  0.598104872  0.975205276
##  [201]  0.232898295  0.924392379  0.021481493  2.230905146  1.915484034
##  [206]  0.840288884  0.856788488  1.371012402  1.430505042  0.803443920
##  [211]  0.685566428  0.512896331  0.802712329  0.204153495 -1.617766871
##  [216] -0.789747695  0.444053370  1.140178106  0.067805477 -0.219142123
##  [221]  1.403427198  1.153188125  0.844356432  0.940303363  2.246018266
##  [226]  0.905744069  1.179060861 -0.900432947  0.421282243 -0.356632602
##  [231]  1.374973618  0.585846735  0.950195890  0.686600744  0.617644874
##  [236]  0.353356140  1.355266225  2.188905530  1.021853588  0.784297032
##  [241] -0.118513802  0.335431461  1.494038915  1.998070490 -1.381454741
##  [246]  1.321497862  0.877334779  0.175633812  1.300209325  0.689289377
##  [251]  0.873631972 -0.303197737 -0.330602743  0.683524985  1.981874988
##  [256]  0.033840976  0.524175299  2.097507589  0.229834083  1.013129785
##  [261]  0.938047393 -1.173039936  0.731985343 -0.069592173  0.553889851
##  [266]  0.127662686  0.430861146  1.066509717  1.677192182  0.233843971
##  [271]  1.070115331  0.855008307 -1.137333525 -0.529520107 -0.177472662
##  [276] -0.043300079  1.824989185  1.505895043  0.672903759  0.674591798
##  [281] -0.225847688  0.983757833  0.177152323  0.704955775  0.613861246
##  [286]  1.449280016  1.088118639  2.448388122  0.483327658  0.282841852
##  [291]  0.760682784  0.801166439 -0.397128598  0.780283728  0.538881227
##  [296]  0.460847221  0.255361547  0.795171850  0.827377276  0.946953883
##  [301]  1.036314876  0.197649933  2.265339176  2.176928306 -1.326346099
##  [306] -1.187397465  0.943781460  0.817429534  0.341407227 -0.877163713
##  [311]  0.128952068  0.552180947 -1.295305273  1.499297068  0.619547873
##  [316]  1.350886803  2.106749052  0.007257150  0.216832973  0.645530981
##  [321]  1.011625112  0.013423549  0.879617361 -0.177365186  1.436280620
##  [326]  0.093791571  0.939420023  0.765543118  1.023981278  0.100281104
##  [331] -1.105054649  0.903278254  1.016117832 -1.698879059  0.363117896
##  [336]  0.318771378  0.583131522 -0.798593433  2.049680879  0.881688944
##  [341]  0.628443746  0.767018048  0.956578148  2.333817853  1.477608259
##  [346]  0.634206457  1.067665066  0.818914524  0.988398637  0.954077970
##  [351]  0.573022869  0.010393629  0.946341043  1.378297340  1.241194683
##  [356]  1.756600287  0.019912032  0.959489403  0.902802780  0.769613067
##  [361]  0.004642948  0.911270502  0.408842606  1.222705520 -0.075147864
##  [366]  0.770632081  0.108564020  0.333994556  0.269133508  0.850225430
##  [371] -0.161575971  2.259955359  1.150188647  0.206091494 -0.832268111
##  [376]  0.179232212  1.145549749  1.216044207  0.128512913  0.093609787
##  [381]  0.468916472  0.477824259  0.543439875  0.943907159  0.497061961
##  [386]  1.468413914  0.152684590  0.391483941  0.947010316  0.458289850
##  [391]  0.006767717  1.302630496  0.480327401 -0.122589964  1.450228325
##  [396]  0.212273941  0.895256112  0.518423030  0.754867134 -0.150233216
##  [401] -1.589637487  1.777381219  0.469995208  1.674927599  2.105454797
##  [406] -0.488542750  0.718146605  0.538156911  0.568023322  1.237399438
##  [411] -0.618906332  0.149371069  0.700743202  0.213174092  1.406370855
##  [416] -0.028963354  0.734015407  1.744689015  0.688815308  0.896486344
##  [421]  0.500561118  0.050560501  0.423092400  0.555418832  0.367002914
##  [426]  0.179085360  0.730274714  0.652945345  0.857122524  0.147532372
##  [431]  0.922376488  0.697998193  0.308871406  0.462264587  0.880500577
##  [436]  0.149275638  0.357810076  1.499707892 -0.707884979  0.351128043
##  [441]  0.847336205  1.569743539  0.373700198  0.367002914  0.673629337
##  [446]  0.522808208  0.687342079  0.127981593  0.492005884  0.858912168
##  [451]  0.938094947  0.434720620  0.385621723  2.009046359 -1.547125539
##  [456]  1.565369795  0.296925263  1.766803453  0.384578850  1.449132230
##  [461]  0.416132010  0.660553983  0.640199457  1.200889915 -1.643984169
##  [466]  0.889692880 -1.435186687  1.486800010  0.951242325  0.275928145
##  [471]  0.687473055  0.038639315 -1.315957641  1.328376803  0.112805315
##  [476]  0.821352152  0.396850979 -1.132851055 -0.844756759  0.799679921
##  [481]  0.642064152  0.267612525  1.262341967 -0.044436465  0.137004453
##  [486]  1.491797861  0.072109416  1.243062823 -1.956993918  0.532982000
##  [491]  1.519626897 -1.276325207  0.998021476  0.244900858  0.070391621
##  [496] -0.207325162  0.364727628  1.273852730  1.286282853  2.180189511
##  [501]  2.314838231  0.163262328 -1.044049419  2.320659011 -0.014712143
##  [506]  0.379061120  1.469360403 -1.130750924 -0.809557584  0.207731711
##  [511]  0.757495438  1.369827467  0.393235513  0.972312610  1.170517518
##  [516]  0.104374893  0.169783873  1.773143605  1.500524929  0.670505876
##  [521]  1.449318989  0.404441107 -0.061377046 -1.247146123  0.914453885
##  [526]  0.139885955  1.546373358  0.156458434  0.075362690  1.092984268
##  [531]  0.856476657  1.302495847  2.436196091  0.074409651  0.530777744
##  [536]  0.735093639  1.386215285 -0.714502747  2.132434689  2.195638212
##  [541]  0.181719159  0.387856465  0.553306758  0.899126280  0.813443425
##  [546]  1.076771274  0.708829827  0.813673742  0.670288088  1.490987107
##  [551] -0.103398358  0.322483227  0.236009787  1.057336761  0.894132974
##  [556] -0.188636091 -0.869310213  0.796987684  0.954792822 -0.674940404
##  [561]  1.799933272  1.335203422  1.064068755  0.307079505  1.413681458
##  [566]  0.271692519  1.074438266  1.684803723 -1.293274347  0.049904767
##  [571]  0.339049411  1.259270691  0.268058201  1.109927550  0.716646301
##  [576]  0.590316238  0.782717302  0.536799044  0.360614501 -0.023505691
##  [581]  0.718558384 -0.301630310  0.389682136 -1.565070367  1.416887757
##  [586]  0.907947253  1.683646864  2.062317512 -1.619745299  2.292827335
##  [591]  1.361005764 -0.005876721  0.247336582  1.437066664  2.202722262
##  [596]  1.289177920 -1.174803385  0.006173141  0.705834632 -0.019618056
##  [601]  0.350807488  0.688177057  1.447368303 -0.709362938  1.323040930
##  [606]  0.610842739 -1.709438925  0.194753341  1.980459445  1.176575080
##  [611]  1.519828406 -0.225006974  0.542810283  0.962522232  0.962289153
##  [616]  0.328129604  0.050826719 -0.834918967  0.654166004 -0.946261440
##  [621] -0.946691085  1.312159756  0.377881534  0.395474929  0.874647554
##  [626]  0.932184587  1.444639968  0.947883944  0.354999258  0.431296875
##  [631]  1.173566246  0.846101150  2.106028447  0.733643189  1.042840297
##  [636]  0.808079777  0.975805096  0.458257651  1.448290027  0.889551542
##  [641]  0.185256466  1.094422131 -0.119918107  1.555641051  0.277409390
##  [646]  0.654972395  1.269577964  0.805022203  1.135078417  1.357685964
##  [651]  0.868475958  0.417784807 -1.618726748 -0.046026211  0.418676019
##  [656]  1.278740073  0.652831659  0.502063510  0.403540302  0.173805351
##  [661]  0.441226628  0.195889883  1.222227734  1.560139033  0.383291356
##  [666]  1.243615853  1.491768346  0.891018893  1.387643428  0.395254660
##  [671] -1.998145386  0.260851595  0.857766552  0.360879196  0.767839191
##  [676]  0.967575567  0.623695063  2.284907192  1.231238222  0.064083957
##  [681]  0.890127266 -2.050819386  0.317838352  0.460739156  0.530456645
##  [686]  0.192553744  0.309154204  2.060113780  0.784297032  1.385051157
##  [691]  0.306332099  1.421617362  0.376708716  0.386175560  0.465287227
##  [696]  0.454266918  0.919401230  1.258710355 -0.860931016  2.194724970
##  [701]  0.254937567  1.426829356  1.450302835 -0.243857190  0.449139886
##  [706]  0.734723096  1.419183578  0.594343505  0.523806898  0.385195226
##  [711]  1.415254174  0.574779256  0.848113798  0.954316756  0.691028650
##  [716] -1.349218549  0.601918435  2.376595897  0.765963587 -0.883580506
##  [721]  1.765572130  1.222193643  1.067665066  0.146551056  0.898492369
##  [726]  0.316089668  0.528083481  0.118735740  0.458754536  0.369146705
##  [731]  0.318330797  0.761710676  0.391559252  0.994514047  2.419645526
##  [736]  0.551952446  0.622900722  0.784857631  0.674828930  1.355736110
##  [741]  1.059446782  0.562384317  1.044490395  0.152342230  0.939072665
##  [746]  1.240781682  0.616201825  0.309025284 -0.170792810  0.127095440
##  [751]  1.458153088  0.730528685  0.577849222  0.396531942  0.451886700
##  [756]  0.464422144  0.493309053  0.775146977  0.574094814  2.188702454
##  [761]  1.203520754  1.437242678  1.186752960  1.502255066  0.404598339
##  [766]  0.932531939  1.923514306  0.811180870  1.552085399  0.216886969
##  [771]  0.912284300  0.084208964  0.489583262  1.318105665  1.699100409
##  [776] -0.006343545 -0.084948222  1.298301919  1.324717835 -0.051299316
##  [781]  1.381787085  1.453451391  1.426691779 -0.263419958  0.149833664
##  [786]  0.887129121  0.539499102  0.441462721  1.582976921  1.103059415
##  [791]  0.434436030  0.753574632  0.516346309  0.202263993  0.499144284
##  [796]  1.246055338 -0.018574681  0.080089378  2.137338309  0.619871743
##  [801]  0.072573665  0.322579416  0.589111912  1.794759222  0.450822719
##  [806]  0.582022214  0.169737712  0.746453544 -0.097726187  1.201283860
##  [811]  0.988272039  0.559445682  0.497456796  0.303674241  0.345618569
##  [816]  0.668205316  0.451660071  0.431320756  0.125604131  0.418845292
##  [821]  1.400945240  1.357343676  0.481777662  1.223316406  0.458489791
##  [826] -0.719236818  0.523223321  1.671278377  0.567524972 -0.647102164
##  [831]  1.705247098  0.127173820  0.405689713  1.302495847  0.460905870
##  [836]  1.636265243  0.555537082  2.053257828  0.357220855  2.073680004
##  [841]  0.830259039  0.810403439  0.382446069  0.939420023  0.109907738
##  [846]  1.464263294  0.765139851  0.887691372  0.512073821  0.432219745
##  [851]  0.203258115  1.029877961  0.063742027  0.782137363  1.353069000
##  [856]  1.044695285 -1.482280552  0.927731436  0.435481240  0.376683328
##  [861]  1.063007978  0.724692040  0.428909431  0.724759823  1.644403340
##  [866]  0.739553558 -1.879315434  0.785423552  1.046380213  0.189382740
##  [871]  0.391878779  0.806076194  0.884292286  0.548690223  0.671082519
##  [876] -1.048924908 -0.273666874  0.016281215  1.807956959  0.483327658
##  [881]  0.537257515  0.120455815  0.442262966  1.225336684  0.092080534
##  [886]  0.732549262  1.056334842 -1.289570803  0.644608768  0.266796727
##  [891]  0.796263245 -1.074140636  0.737334764  0.883463768  0.722941742
##  [896] -0.522969367  0.845356302  0.911621358 -1.329429349  0.613290585
##  [901]  0.961694092  0.012084593  1.455570150  0.837718974  0.859883002
##  [906]  0.750009957  1.195348050  0.163051475  1.179531928  0.108002303
##  [911] -0.111967332  0.478949673  0.412126157  0.407480816  1.704189898
##  [916] -1.208850881  0.236202802  2.424155949  0.315213023  1.262485662
##  [921]  0.021577299  0.326253296  0.066028759  0.683524985  0.176559579
##  [926]  0.866748260  0.256801788 -0.948074520  0.428670706  1.351375449
##  [931]  1.340083843  0.423746153  0.601689890  0.283025285  2.278155683
##  [936]  0.144090099  1.069922196  0.381899915  0.664041667  0.436550260
##  [941]  0.463721389  0.813443425  0.312102698  0.523062637  1.173299325
##  [946] -2.056940624  0.194066359  1.183955297  0.247090308 -0.074967504
##  [951] -0.179147933  1.416550996 -0.068518098  1.901083600  2.188847928
##  [956]  0.702561175  1.525881532  0.508501070  0.508845609  1.022653828
##  [961]  0.345618569  1.462743651  0.028101418  0.578665126  0.535662768
##  [966]  0.439528280  0.294648265  0.950541833  0.633538076  0.199031428
##  [971] -0.675795148  0.950546128 -0.855390636 -1.575184383  0.662050946
##  [976]  2.400648704  0.919322134  1.313749894  0.302049920  0.268072629
##  [981]  0.873392642  0.497698609  0.430454220  1.233823815  1.529245308
##  [986]  1.943396053  1.317522514  1.040863844 -1.288494224  1.059082969
##  [991]  0.718188573 -0.196083031  1.290754365  2.170808918 -0.539178175
##  [996] -0.240311115 -1.648028803  1.393765819  0.315213023  0.357601095
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

## Warning in densfun(x, parm[1], parm[2], ...): NaNs produced
```

```r
fit2
```

```
##       mean          sd    
##   0.51921757   0.25915418 
##  (0.08195175) (0.05794532)
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
## [1] -0.2630822 -0.4537855  0.3908646  0.4776744  0.5502591  0.2420518
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
## [1] 0.0642
```

```r
sd(results2$thetastar)  #the standard deviation of the theta stars is the SE of the statistic (in this case, the mean)
```

```
## [1] 0.9031578
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
## t1*      4.5 0.002002002   0.9163927
```

One of the main advantages to the 'boot' package over the 'bootstrap' package is the nicer formatting of the output.

Going back to our original code, lets see how we could reproduce all of these numbers:


```r
table(sample(x,size=length(x),replace=T))
```

```
## 
## 0 3 6 7 8 9 
## 1 1 3 2 2 1
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
## [1] 0.0097
```

```r
se.boot
```

```
## [1] 0.9120349
```

Why do our numbers not agree exactly with those of the boot package? This is because our estimates of bias and standard error are just estimates, and they carry with them their own uncertainties. That is one of the reasons we might bother doing jackknife-after-bootstrap.

The 'boot' package has a LOT of functionality. If we have time, we will come back to some of these more complex functions later in the semester as we cover topics like regression and glm.

