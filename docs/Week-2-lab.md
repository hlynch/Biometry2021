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
## 2 3 4 5 6 7 8 9 
## 1 1 2 1 2 1 1 1
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
## [1] -0.0499
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
## [1] 2.651292
```

```r
UL.boot
```

```
## [1] 6.248908
```

Method #2: Simply take the quantiles of the bootstrap statistics


```r
quantile(xmeans,c(0.025,0.975))
```

```
##   2.5%  97.5% 
## 2.7000 6.1025
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
##    [1] 4.0 5.3 4.3 4.7 4.3 5.0 4.1 4.2 6.2 5.2 5.4 4.5 5.1 5.1 4.5 5.0 4.2 4.9
##   [19] 4.0 6.7 3.2 3.9 5.5 4.5 4.6 3.2 6.2 3.5 3.9 5.2 4.5 4.4 5.2 3.9 4.0 5.8
##   [37] 4.8 3.8 4.8 3.3 5.7 5.1 7.2 3.4 5.3 3.8 3.0 4.8 5.0 3.1 4.0 3.3 5.8 4.0
##   [55] 3.8 5.3 5.2 5.9 3.8 3.1 4.5 5.1 3.4 4.5 2.4 6.2 4.7 3.6 5.8 4.8 3.6 4.6
##   [73] 3.9 4.2 2.8 3.9 3.8 5.9 6.1 2.1 4.3 5.0 5.1 3.4 2.9 3.4 5.4 1.8 4.3 4.9
##   [91] 4.9 5.9 3.6 5.2 4.1 6.2 5.5 3.8 5.4 3.8 3.7 6.4 4.5 4.4 3.4 4.6 4.7 3.3
##  [109] 4.4 5.5 4.1 6.0 4.5 5.1 4.0 3.7 4.7 5.8 3.6 4.2 5.1 4.0 3.7 4.7 3.7 5.6
##  [127] 3.7 4.9 4.2 3.9 4.6 4.8 5.5 3.8 4.2 4.6 4.4 4.4 2.8 4.7 5.1 5.3 4.9 5.8
##  [145] 5.1 4.1 4.6 5.0 5.0 5.8 5.1 3.9 5.4 5.5 4.3 4.7 5.4 5.9 4.4 4.9 3.9 3.8
##  [163] 5.1 3.5 4.2 3.8 4.6 4.3 5.5 5.5 4.3 4.5 4.8 3.4 6.3 3.3 4.9 5.2 4.5 4.2
##  [181] 5.1 3.8 5.0 4.5 3.8 5.8 4.8 5.4 5.4 4.6 3.1 2.7 4.3 4.2 5.7 3.9 5.3 5.8
##  [199] 3.9 4.5 3.6 3.2 6.3 5.8 2.8 5.2 4.6 5.1 3.6 4.0 4.3 3.9 4.7 4.1 4.4 2.9
##  [217] 3.8 3.4 6.3 3.5 4.4 6.8 4.6 3.8 4.8 3.6 4.1 4.6 5.3 3.9 3.2 5.7 4.3 5.5
##  [235] 3.7 4.4 6.1 3.9 5.2 4.2 3.7 4.5 4.9 5.3 4.8 5.8 3.7 4.2 3.6 4.3 3.1 5.2
##  [253] 4.4 5.7 3.7 1.9 5.4 4.7 4.9 5.2 5.0 3.5 4.5 4.6 4.9 3.8 4.9 3.7 5.2 3.5
##  [271] 4.0 4.6 4.0 4.5 5.5 4.5 3.0 4.5 4.2 3.7 5.1 5.2 2.9 4.7 5.3 3.8 3.6 5.0
##  [289] 4.5 4.0 4.9 3.8 4.8 3.9 3.2 4.7 4.6 5.7 5.9 4.9 2.9 3.7 6.0 4.9 4.6 4.3
##  [307] 4.3 4.1 6.2 4.4 5.1 6.8 5.0 2.8 2.7 4.7 4.9 5.2 3.2 4.7 4.5 5.0 4.4 4.4
##  [325] 3.6 4.3 4.2 3.5 4.0 3.7 3.6 5.6 6.2 4.9 3.6 3.7 2.8 3.1 3.2 4.7 3.8 3.6
##  [343] 6.1 4.4 5.1 4.7 4.4 5.0 4.6 5.3 4.0 4.0 4.3 4.3 6.3 3.2 3.8 4.5 5.1 3.7
##  [361] 3.5 5.4 5.3 4.8 3.4 5.2 5.5 3.8 5.0 3.8 4.4 5.9 4.6 3.4 5.4 5.0 3.5 2.2
##  [379] 5.2 3.1 4.7 4.7 4.6 4.5 4.0 3.9 4.1 3.7 4.9 3.2 2.4 4.4 5.4 5.0 5.8 4.8
##  [397] 6.1 3.7 5.1 3.5 3.8 4.8 4.9 4.4 5.1 3.8 5.4 5.5 5.1 4.0 5.1 3.6 4.6 5.6
##  [415] 4.9 5.5 4.1 4.8 6.1 6.0 3.8 5.6 5.2 5.9 3.9 3.1 4.6 3.1 4.4 4.3 3.9 4.6
##  [433] 3.5 4.1 6.6 4.9 4.3 3.5 5.7 4.4 5.6 6.1 4.3 3.8 4.9 3.8 5.3 4.0 4.2 6.9
##  [451] 5.6 4.9 4.8 4.8 4.8 2.8 4.0 4.8 4.3 3.5 5.8 5.7 4.0 5.5 4.3 4.4 5.4 4.9
##  [469] 5.0 2.7 4.9 4.0 3.7 5.5 3.9 5.1 2.9 4.6 5.1 4.4 4.2 3.7 5.0 4.4 4.7 4.4
##  [487] 4.8 3.3 5.5 4.8 4.4 4.8 5.3 4.5 5.4 4.3 6.1 4.2 3.8 6.4 3.5 5.1 6.5 5.8
##  [505] 3.6 3.5 4.5 3.8 4.0 3.4 2.9 4.3 5.3 5.6 3.4 5.8 4.3 3.9 4.1 5.3 4.9 4.3
##  [523] 5.1 4.6 4.0 4.5 5.4 4.1 5.0 3.5 6.2 6.5 3.9 3.5 5.7 4.6 3.0 4.9 4.7 5.9
##  [541] 5.0 5.4 4.4 5.8 3.6 4.3 4.3 5.3 5.7 5.2 3.4 4.7 4.9 4.6 2.5 2.6 5.7 3.0
##  [559] 5.2 3.8 4.3 4.2 5.0 5.7 3.9 4.3 4.4 4.7 4.3 4.5 4.4 5.4 4.7 5.5 2.9 3.0
##  [577] 5.3 4.3 4.1 2.9 3.6 6.1 3.9 4.3 3.5 4.7 5.4 3.2 6.1 5.1 4.9 4.2 5.5 5.4
##  [595] 5.3 5.1 4.8 3.6 4.5 4.6 4.4 4.3 4.7 3.7 3.8 3.7 5.0 4.8 5.8 3.8 2.7 5.6
##  [613] 3.4 4.6 4.4 4.3 3.8 3.8 3.0 5.4 2.8 3.9 4.9 3.2 3.9 4.6 3.9 4.1 4.1 2.1
##  [631] 4.7 3.7 3.4 3.3 5.4 6.7 5.6 6.1 4.4 5.0 4.2 3.5 3.1 4.3 5.6 2.7 4.6 3.9
##  [649] 4.6 3.8 6.3 5.2 5.3 2.9 4.4 5.0 3.9 4.5 4.2 4.0 3.5 2.5 5.9 4.8 4.2 4.0
##  [667] 5.9 4.1 4.8 4.7 4.6 4.3 5.7 4.5 4.6 3.0 3.7 3.4 3.9 6.2 4.8 4.9 4.8 7.0
##  [685] 4.4 3.9 4.7 4.9 3.3 3.3 5.1 3.4 3.7 4.1 4.3 4.4 5.3 3.0 3.0 5.1 5.5 5.2
##  [703] 3.0 4.0 4.5 4.5 4.5 4.2 4.5 3.4 3.7 4.0 2.9 5.3 4.8 4.6 3.3 4.7 3.3 7.1
##  [721] 3.3 3.9 2.8 5.3 3.4 4.3 6.1 4.5 5.1 6.4 3.3 3.9 3.6 5.3 6.0 5.3 2.7 5.0
##  [739] 4.4 5.8 5.4 3.9 5.1 5.8 4.4 5.3 3.8 3.5 4.0 4.8 4.7 6.0 4.5 4.3 6.3 3.4
##  [757] 4.4 4.0 3.6 4.8 5.2 5.0 5.0 5.0 4.7 5.9 3.8 3.8 5.6 4.5 5.9 4.0 3.8 4.0
##  [775] 5.3 3.8 4.1 5.2 5.7 5.0 5.5 5.6 3.5 4.3 2.3 5.5 3.1 4.0 4.7 3.6 4.7 4.7
##  [793] 5.8 3.6 4.8 4.1 6.1 4.0 3.7 4.2 5.3 3.2 2.4 6.6 4.8 5.8 5.2 6.4 5.0 4.2
##  [811] 4.1 4.3 3.5 4.3 5.6 4.4 4.3 2.9 4.7 4.7 3.5 4.3 4.3 4.9 2.2 5.6 3.2 2.6
##  [829] 4.3 4.5 5.4 5.6 4.7 3.8 3.0 3.9 4.9 4.8 3.8 4.7 3.6 4.2 4.8 4.2 6.1 3.1
##  [847] 3.7 6.5 4.9 4.4 4.9 2.5 3.1 5.0 5.2 3.6 4.3 4.0 4.7 6.3 3.8 4.3 4.6 3.9
##  [865] 3.6 3.7 4.9 4.0 4.7 5.3 5.0 3.8 4.0 5.4 4.9 5.4 5.8 5.1 3.7 4.0 5.2 2.2
##  [883] 4.9 5.4 4.1 4.0 4.3 5.8 2.8 4.3 4.2 4.9 5.6 5.1 4.9 3.7 4.7 4.5 3.9 4.2
##  [901] 5.8 4.6 5.9 6.6 3.6 4.5 3.8 3.8 5.6 4.3 5.2 5.3 4.4 6.3 4.7 4.6 4.8 4.5
##  [919] 5.7 3.7 5.0 4.5 4.6 3.3 3.9 4.4 3.7 3.3 3.8 4.2 4.4 4.4 4.0 3.8 4.4 3.3
##  [937] 4.1 4.4 5.3 3.6 4.3 4.8 5.3 4.4 3.7 4.5 3.4 5.6 4.5 4.7 5.2 4.8 5.3 3.9
##  [955] 4.1 3.7 4.4 3.8 3.7 4.3 4.4 3.7 4.3 3.7 3.5 5.6 3.2 4.6 4.6 3.8 3.1 4.3
##  [973] 4.7 4.7 5.1 2.6 5.0 4.7 6.4 5.3 3.0 4.7 2.9 4.8 6.1 4.8 4.2 4.6 5.5 3.6
##  [991] 4.6 5.1 2.4 3.5 3.1 6.9 3.3 5.1 6.5 4.2
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
##   2.8   6.3
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
##    [1] 4.9 5.4 3.2 4.8 3.8 4.9 2.8 4.3 3.5 5.0 4.7 6.3 6.4 3.7 4.8 4.6 5.3 3.4
##   [19] 4.5 4.5 5.8 4.6 5.3 5.0 4.3 5.3 3.2 4.5 5.5 4.1 4.7 4.8 5.0 4.9 4.8 3.7
##   [37] 2.4 5.9 5.6 4.1 3.2 4.2 4.5 3.6 2.9 5.4 3.4 4.7 5.5 3.0 5.0 4.2 5.3 4.3
##   [55] 2.7 4.8 4.7 3.9 3.7 5.5 4.5 4.4 5.3 4.6 4.7 5.5 3.5 5.1 4.2 4.5 5.1 4.0
##   [73] 3.8 4.0 4.1 4.8 2.5 4.4 4.9 4.8 3.6 4.0 4.7 4.2 5.1 3.8 5.1 4.3 4.5 5.5
##   [91] 5.8 3.7 4.0 3.9 5.5 5.1 3.9 2.9 4.6 4.7 4.0 6.1 2.3 4.5 5.6 3.8 3.6 4.7
##  [109] 4.6 5.2 4.4 6.2 5.9 5.1 6.6 4.5 5.1 3.4 4.1 4.5 2.7 3.2 5.6 3.5 6.3 3.5
##  [127] 3.9 2.2 5.5 4.0 5.6 6.3 4.0 3.9 4.4 4.3 3.0 2.8 6.1 4.0 4.4 5.2 3.6 4.2
##  [145] 4.4 4.7 2.9 5.1 4.2 4.7 4.3 4.5 5.0 6.3 4.7 4.6 3.1 4.4 4.1 5.1 2.7 3.7
##  [163] 4.4 3.5 3.0 3.8 4.5 3.1 5.8 4.2 4.4 4.1 4.9 4.3 4.2 5.4 4.9 5.7 4.7 3.4
##  [181] 4.0 3.3 4.5 4.9 3.5 5.3 3.7 6.3 3.8 4.3 5.2 4.0 6.2 3.1 3.3 3.0 4.1 4.7
##  [199] 5.4 5.1 4.2 4.2 4.3 3.5 5.3 5.3 3.3 3.3 5.1 3.3 4.9 3.6 5.0 5.2 4.5 4.7
##  [217] 3.0 3.5 5.2 4.5 2.6 4.6 3.4 1.7 3.8 4.7 5.4 4.2 4.3 5.3 5.1 4.4 5.0 3.4
##  [235] 3.9 3.9 3.0 6.3 3.5 3.1 5.0 5.0 4.2 5.1 3.3 4.8 4.0 5.6 5.6 3.9 4.5 4.5
##  [253] 6.0 3.2 4.7 3.4 3.0 4.4 6.3 4.3 3.9 4.6 3.2 4.5 4.4 5.1 6.1 4.3 4.0 4.3
##  [271] 4.4 6.5 3.0 4.0 3.7 5.1 4.9 6.6 5.2 3.9 4.4 4.4 3.6 4.5 5.9 4.4 4.9 5.2
##  [289] 3.6 4.5 5.0 4.8 6.2 3.2 3.3 4.2 4.8 4.3 5.1 4.0 3.9 4.5 4.9 4.6 5.1 5.7
##  [307] 4.5 4.1 4.6 5.8 4.5 5.4 4.7 5.3 3.8 5.0 4.4 5.8 3.8 6.2 3.6 5.3 4.1 4.3
##  [325] 3.8 4.2 5.0 4.3 3.7 4.9 4.2 3.7 5.2 5.6 5.1 3.4 3.6 6.1 3.9 5.2 5.0 4.4
##  [343] 4.0 3.6 4.6 3.7 3.6 4.9 3.8 4.7 4.1 4.5 4.9 3.9 3.4 3.2 4.6 4.7 3.7 6.0
##  [361] 4.1 3.9 3.6 5.3 4.5 4.3 3.4 3.2 4.5 4.4 5.8 4.4 5.3 3.0 5.3 4.2 4.7 5.1
##  [379] 4.8 4.3 4.7 4.6 5.1 3.0 4.3 5.2 5.9 5.1 4.8 4.4 2.3 4.7 5.8 3.4 4.1 4.4
##  [397] 3.6 4.8 4.6 5.3 6.1 5.2 4.2 2.5 5.1 4.3 5.2 4.0 5.4 5.6 4.1 5.3 3.7 3.6
##  [415] 6.1 5.3 3.9 3.5 4.3 4.0 3.9 4.6 3.7 3.9 4.3 3.9 4.6 3.5 4.9 6.2 4.7 5.3
##  [433] 4.6 4.6 3.8 4.5 3.4 5.1 6.1 4.1 5.7 5.1 4.2 4.3 4.2 4.9 4.9 5.0 5.7 4.3
##  [451] 2.9 4.1 2.9 3.9 4.3 5.7 4.0 4.5 4.7 5.0 5.3 5.0 5.2 4.2 4.7 4.7 4.6 3.6
##  [469] 5.1 4.4 7.2 3.8 4.9 5.9 4.9 4.4 3.8 4.4 5.2 3.2 4.4 4.5 5.0 4.0 5.4 5.1
##  [487] 5.3 4.4 5.9 3.9 5.7 3.9 4.7 3.6 4.9 4.8 4.7 3.4 4.8 5.9 4.4 4.9 3.3 3.7
##  [505] 4.9 3.1 4.0 4.1 4.3 5.1 3.7 4.0 3.4 4.8 6.4 3.4 4.8 4.2 4.9 5.0 4.8 4.8
##  [523] 6.3 5.0 4.3 2.2 5.0 4.5 4.4 5.1 3.7 5.3 5.9 4.4 5.2 5.1 5.2 5.8 3.9 5.5
##  [541] 5.3 4.9 3.5 2.8 3.3 3.0 4.4 5.0 5.4 5.2 3.4 3.8 4.7 5.0 2.7 5.3 5.4 5.4
##  [559] 4.4 3.8 6.3 3.2 6.9 4.3 3.7 3.5 4.1 4.1 4.0 4.2 3.2 4.5 4.2 4.8 3.7 4.4
##  [577] 4.1 3.1 4.9 4.3 5.2 3.9 3.4 5.3 5.4 3.4 6.2 6.2 6.3 4.6 4.7 5.1 5.6 4.7
##  [595] 4.5 5.1 4.8 4.0 4.4 3.8 4.1 5.1 4.6 4.2 5.7 4.9 4.7 3.1 3.9 4.7 5.3 6.5
##  [613] 3.7 3.7 2.3 4.9 5.7 4.9 5.6 3.8 5.4 2.8 3.7 4.8 3.8 5.7 5.9 4.3 4.6 4.0
##  [631] 3.8 4.1 3.9 4.4 4.9 5.2 6.1 2.9 6.7 5.0 4.3 4.2 4.3 5.6 4.0 4.2 3.8 4.2
##  [649] 3.8 3.1 3.7 5.6 5.7 4.2 4.3 3.7 3.7 3.4 4.0 4.9 4.3 4.0 4.4 4.9 3.6 4.4
##  [667] 5.1 5.3 3.8 3.9 4.6 2.7 3.1 3.9 4.2 4.1 3.8 3.7 3.4 4.5 4.0 6.2 4.0 5.1
##  [685] 4.2 4.7 4.6 3.5 4.5 4.5 5.4 5.4 5.9 4.1 6.3 4.5 4.6 6.0 5.9 4.0 4.3 4.5
##  [703] 5.3 4.7 3.9 3.6 4.7 4.4 6.6 1.6 4.1 5.1 5.5 4.5 3.3 3.1 4.9 4.0 2.6 4.1
##  [721] 4.2 5.1 4.4 5.4 4.5 3.8 3.9 2.6 4.4 3.9 3.4 6.0 4.1 5.6 3.2 5.0 6.1 4.2
##  [739] 4.6 3.7 5.3 4.1 4.7 5.0 3.1 4.1 3.4 4.1 5.5 4.7 4.3 5.0 4.3 4.2 5.8 5.0
##  [757] 4.3 5.3 4.5 4.1 5.8 5.5 5.5 5.3 4.3 5.0 5.2 3.7 5.3 4.5 5.2 5.7 4.0 4.0
##  [775] 3.2 4.3 5.0 3.1 4.1 4.0 2.5 4.6 4.2 6.1 4.9 4.6 4.0 3.5 3.6 5.6 4.9 4.7
##  [793] 4.8 5.3 2.7 5.6 4.0 3.7 2.8 4.9 3.5 3.5 2.9 2.9 4.1 4.6 2.2 3.5 4.3 5.3
##  [811] 3.0 4.3 4.2 5.1 5.4 5.9 3.7 4.1 3.8 3.8 4.2 4.2 3.0 5.1 3.8 3.4 4.3 2.9
##  [829] 4.3 3.8 5.1 5.6 5.2 5.0 4.8 2.9 4.3 4.8 4.8 5.4 4.1 3.7 5.2 4.3 3.4 5.1
##  [847] 3.0 4.6 4.5 6.3 3.5 4.0 5.6 3.5 5.2 6.6 5.3 5.2 2.3 4.3 3.6 5.2 4.9 6.9
##  [865] 4.9 2.1 4.4 4.9 3.2 4.2 6.7 2.8 5.5 4.6 2.8 3.4 3.9 5.2 2.3 3.5 5.1 3.3
##  [883] 4.4 5.4 4.6 5.4 3.6 4.6 5.9 4.1 4.5 4.2 4.7 5.1 4.9 3.8 5.1 5.1 5.0 3.7
##  [901] 3.3 4.6 4.0 4.9 3.8 5.5 4.1 2.7 4.1 3.6 5.8 5.4 4.9 4.3 4.2 4.1 4.6 4.9
##  [919] 4.1 5.2 5.3 4.5 4.6 3.3 3.4 5.8 5.2 4.8 4.5 4.1 4.9 5.1 4.5 5.6 3.0 4.5
##  [937] 5.1 4.1 3.8 5.2 5.7 5.3 3.7 5.1 3.7 3.6 3.8 3.8 3.4 3.8 5.3 4.6 4.1 4.6
##  [955] 5.0 5.1 3.9 3.2 4.8 4.7 5.1 4.7 4.5 3.2 4.8 4.9 4.1 4.0 4.8 4.5 5.6 5.6
##  [973] 4.4 4.3 5.1 4.7 4.8 5.3 4.6 4.6 3.3 5.6 4.5 6.0 5.4 4.3 4.2 4.7 4.0 4.5
##  [991] 4.5 4.3 5.8 4.1 3.9 3.8 4.4 4.7 4.6 3.8
## 
## $func.thetastar
## [1] -0.0359
## 
## $jack.boot.val
##  [1]  0.43052326  0.31930836  0.28018293  0.19882698  0.01607143 -0.08179487
##  [7] -0.20027855 -0.32507042 -0.37767584 -0.51934605
## 
## $jack.boot.se
## [1] 0.9313113
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
##    [1] 4.8 6.9 4.3 5.3 3.8 3.8 4.4 4.0 4.3 5.7 3.4 4.2 2.2 2.7 4.7 4.4 5.7 5.0
##   [19] 3.6 5.6 2.9 3.6 2.3 3.0 3.9 3.5 4.0 3.8 3.7 3.7 4.8 4.8 5.5 4.6 4.4 5.6
##   [37] 4.3 3.9 4.7 5.2 3.3 5.7 4.9 3.5 4.8 3.8 4.3 3.9 5.6 6.1 4.9 4.1 4.0 3.3
##   [55] 5.1 5.3 3.7 3.5 6.3 4.5 4.0 3.9 4.6 4.9 5.7 4.0 3.6 3.8 4.0 4.2 5.7 2.0
##   [73] 6.5 4.7 6.1 3.6 4.7 3.0 3.6 5.2 4.6 5.9 4.7 4.3 6.1 4.3 3.9 4.5 4.9 4.4
##   [91] 4.4 4.9 5.5 4.2 5.3 3.3 5.1 3.8 4.1 3.8 5.1 4.1 5.5 5.5 2.3 3.6 5.6 5.4
##  [109] 5.2 4.2 3.6 4.5 3.4 3.2 4.6 6.3 6.1 5.2 3.7 5.0 5.7 3.6 4.8 4.2 4.5 5.5
##  [127] 2.9 5.5 4.4 3.3 6.0 4.3 5.3 3.2 4.9 6.2 6.1 3.6 6.5 3.8 4.5 5.2 4.3 5.6
##  [145] 5.2 3.5 5.0 3.2 2.8 3.3 4.5 5.4 4.1 4.5 6.8 4.7 2.9 3.9 5.3 5.7 4.7 5.6
##  [163] 3.2 4.4 4.5 4.6 3.2 5.2 4.0 3.4 3.6 4.8 6.1 4.2 3.9 5.7 5.0 6.1 4.0 5.0
##  [181] 4.1 5.1 4.7 4.7 6.1 3.9 4.8 5.4 4.6 5.7 3.1 4.3 4.5 4.0 4.6 4.3 4.8 3.3
##  [199] 6.2 4.6 3.5 5.6 3.2 4.6 5.9 4.2 3.4 6.0 6.3 3.1 6.0 3.8 3.3 4.2 4.5 4.5
##  [217] 2.0 4.1 4.6 4.3 4.0 4.4 3.5 5.8 4.4 3.8 6.3 4.1 5.0 3.8 5.3 5.5 4.5 5.5
##  [235] 5.9 2.1 2.1 5.8 4.0 3.4 4.0 5.3 5.2 4.9 4.3 3.5 3.4 3.9 4.9 2.6 5.3 5.0
##  [253] 5.4 3.7 6.6 4.7 5.8 3.5 4.4 5.4 4.1 4.2 5.4 4.5 4.5 4.3 5.5 4.3 4.0 4.9
##  [271] 5.3 5.8 5.3 4.5 4.5 5.3 6.0 5.3 3.7 4.1 4.8 4.9 4.2 4.1 4.2 3.9 4.6 4.8
##  [289] 5.4 4.7 3.9 4.7 4.8 4.5 5.9 4.4 6.7 5.3 3.5 5.5 2.9 3.9 4.5 5.1 3.2 5.1
##  [307] 4.2 5.9 4.3 3.6 4.1 4.8 5.7 4.5 4.0 3.8 4.9 5.0 4.0 3.7 5.1 5.2 3.8 4.6
##  [325] 5.0 4.7 3.0 3.3 4.7 5.1 5.5 3.4 3.9 4.9 4.9 4.3 5.5 4.7 3.4 3.3 4.8 3.7
##  [343] 3.8 3.2 4.5 4.0 4.8 6.6 5.1 5.6 5.2 4.3 3.9 4.5 3.0 6.1 5.3 3.4 4.9 4.2
##  [361] 5.3 2.8 4.0 4.5 5.5 3.7 5.1 4.8 4.7 4.2 2.6 5.1 3.9 4.9 3.4 3.6 5.5 4.1
##  [379] 3.9 5.5 5.8 5.1 5.1 3.7 6.0 3.5 5.8 4.6 4.8 4.2 3.6 3.9 7.0 3.9 5.2 3.7
##  [397] 5.4 3.4 5.1 2.6 2.9 3.5 4.6 4.3 4.3 5.6 5.7 5.2 4.4 4.5 4.4 4.2 5.0 3.9
##  [415] 3.5 6.0 3.0 5.4 5.2 4.9 4.8 5.7 3.8 3.9 4.1 4.4 3.1 5.5 3.3 3.7 3.7 4.9
##  [433] 5.5 5.2 4.0 5.5 5.6 4.8 6.1 3.8 6.2 3.0 3.9 5.8 5.4 2.6 4.2 4.8 4.5 5.1
##  [451] 4.2 4.5 5.7 5.9 3.2 5.6 3.2 5.4 4.8 5.6 4.1 3.6 3.7 4.8 3.6 4.9 4.7 4.1
##  [469] 4.3 4.5 5.0 4.7 5.1 4.8 5.0 5.9 4.1 5.8 4.1 6.0 4.1 4.6 2.8 4.4 4.8 3.9
##  [487] 6.0 3.0 2.9 4.2 4.2 3.6 4.3 4.8 4.9 6.0 4.0 5.6 3.0 4.2 5.0 4.1 4.1 5.8
##  [505] 5.2 4.7 5.2 4.0 3.7 5.1 4.3 3.4 5.1 4.1 4.7 4.8 3.7 4.4 5.5 4.8 4.9 2.4
##  [523] 5.2 5.6 5.7 2.3 3.7 5.0 5.1 3.7 4.4 4.4 5.4 4.4 4.7 4.6 3.5 5.2 3.9 4.1
##  [541] 4.6 2.8 4.2 2.7 6.2 3.7 4.6 5.4 3.8 5.3 3.9 3.4 5.5 4.9 3.3 4.8 3.6 4.9
##  [559] 4.9 5.8 3.5 5.4 3.9 3.4 4.6 3.4 4.8 3.9 4.0 5.8 3.4 6.5 4.1 3.4 3.7 3.3
##  [577] 4.8 4.9 3.0 5.1 2.8 6.5 3.8 4.9 4.9 6.4 5.3 4.1 4.6 4.5 5.2 2.5 3.9 4.8
##  [595] 3.3 5.8 3.6 4.9 3.9 4.3 3.3 4.6 3.6 5.1 4.7 5.7 3.4 5.6 4.8 5.0 3.6 4.2
##  [613] 4.2 4.2 4.5 4.7 5.2 4.5 2.9 5.2 5.2 5.7 6.1 3.1 4.4 5.4 4.5 4.0 4.7 6.0
##  [631] 5.0 5.8 4.8 5.1 5.4 4.8 5.0 5.8 3.8 4.8 4.0 4.5 5.1 5.0 5.0 4.1 3.9 5.1
##  [649] 4.5 4.5 5.2 5.5 5.5 4.3 3.8 3.7 4.4 4.6 4.6 3.8 4.2 4.3 5.9 5.4 6.5 3.9
##  [667] 6.0 5.6 3.3 5.4 4.3 4.6 5.5 5.0 3.5 3.9 3.8 4.4 4.6 3.5 3.5 5.8 5.2 5.7
##  [685] 4.1 5.7 4.6 5.8 3.6 3.8 4.6 4.7 5.8 4.8 5.3 4.2 3.3 3.7 4.3 5.1 3.2 3.3
##  [703] 2.9 3.3 6.6 4.9 5.5 5.2 6.1 3.3 4.5 3.9 4.1 3.4 3.6 3.8 3.8 5.6 4.3 3.2
##  [721] 3.6 5.4 5.0 4.0 3.4 4.8 3.4 4.8 4.8 5.5 4.4 5.2 5.0 4.9 5.2 5.2 4.9 4.5
##  [739] 3.5 5.3 6.8 4.6 4.4 3.8 7.0 5.4 5.6 4.0 4.3 5.5 3.7 4.9 3.7 5.7 3.9 5.2
##  [757] 5.2 5.8 5.3 4.6 5.9 5.1 4.9 5.3 4.5 3.9 2.9 4.5 3.2 5.3 2.8 6.3 4.8 4.6
##  [775] 3.1 2.1 5.1 5.7 4.1 4.3 4.2 5.5 2.8 4.0 3.3 4.8 4.3 4.2 4.3 5.2 6.0 3.3
##  [793] 2.9 3.9 3.9 4.6 5.2 5.1 5.6 3.9 6.5 5.3 4.0 4.1 4.8 5.1 4.5 5.2 6.3 3.7
##  [811] 4.3 6.2 5.8 5.1 6.2 5.3 3.8 4.3 2.9 4.0 3.1 4.6 4.5 4.5 4.9 6.2 6.7 3.8
##  [829] 6.4 4.5 4.6 4.9 4.1 6.5 5.2 3.9 4.0 3.9 5.1 4.7 4.9 4.0 4.5 4.7 3.2 5.2
##  [847] 5.0 2.8 4.0 4.6 5.9 4.9 2.6 4.2 4.7 4.2 3.9 6.2 5.3 3.5 4.3 4.2 5.0 6.0
##  [865] 5.0 2.4 5.4 2.7 5.3 5.6 4.4 4.3 5.3 2.8 4.6 4.0 3.8 3.0 3.6 5.9 4.8 4.6
##  [883] 5.9 4.1 4.8 3.3 4.9 5.0 4.5 3.2 4.4 5.7 4.7 5.6 4.0 4.6 4.4 3.8 5.0 5.4
##  [901] 4.0 4.7 3.6 3.8 5.5 5.5 4.8 4.8 4.0 4.4 4.4 5.6 4.3 4.8 5.6 3.2 4.3 4.6
##  [919] 3.2 3.8 4.1 4.2 4.9 4.6 5.2 3.8 3.8 5.4 3.4 4.3 4.8 5.9 5.0 3.1 3.8 3.9
##  [937] 3.6 4.2 4.0 2.1 3.3 4.7 4.1 5.4 3.8 4.1 4.2 4.0 5.1 3.9 4.3 4.4 5.1 5.1
##  [955] 4.3 3.8 2.5 4.5 5.0 5.8 4.9 4.6 4.4 5.8 4.1 3.7 5.0 4.6 3.2 5.1 4.8 4.8
##  [973] 4.3 5.7 5.2 4.3 3.4 4.6 4.7 4.4 4.3 3.6 5.3 6.2 5.1 3.6 4.7 4.5 4.4 4.8
##  [991] 5.6 3.5 3.4 3.4 3.8 3.6 5.4 4.5 4.1 4.8
## 
## $func.thetastar
## 72% 
## 5.1 
## 
## $jack.boot.val
##  [1] 5.5 5.4 5.3 5.4 5.2 5.0 5.0 4.8 4.6 4.5
## 
## $jack.boot.se
## [1] 0.9954396
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
## [1] 0.4375835
```

What is skew? Skew describes how assymetric a distribution is. A distribution with a positive skew is a distribution that is "slumped over" to the right, with a right tail that is longer than the left tail. Alternatively, a distribution with negative skew has a longer left tail. Here we are just using it for illustration, as a property of a distribution that you may want to estimate using your data.

Lets use 'fitdistr' to fit a gamma distribution to these data. This function is an extremely handy function that takes in your data, the name of the distribution you are fitting, and some starting values (for the estimation optimizer under the hood), and it will return the parameter values (and their standard errors). We will learn in a couple weeks how R is doing this, but for now we will just use it out of the box. (Because we generated the data, we happen to know that the data are gamma distributed. In general we wouldn't know that, and we will see in a second that our assumption about the shape of the data really does make a difference.)


```r
library(MASS)
fit<-fitdistr(original.data,dgamma,list(shape=1,rate=1))
```

```
## Warning in densfun(x, parm[1], parm[2], ...): NaNs produced

## Warning in densfun(x, parm[1], parm[2], ...): NaNs produced
```

```r
# fit<-fitdistr(original.data,"gamma")
# The second version would also work.
fit
```

```
##     shape       rate  
##   1.734785   4.164349 
##  (0.713731) (1.983705)
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
## [1]  1.1201255  1.0748725 -0.2371960  0.8685136  2.3226662  0.9587821
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
##    [1]  0.2898207519  0.4675114931 -0.0064623872  0.2351764045  0.6556261213
##    [6]  0.5090491514  0.9432274766  0.6509536470  0.2386468056  1.6188491829
##   [11]  1.4624191290  0.5041831257  0.0857884113  0.2726716582  0.7539328235
##   [16] -1.0134420768  0.1175200926  0.9569631065  0.2040324337  0.5701825033
##   [21]  0.7215626290  0.5041220624 -0.2617791842  0.4177754708  0.0824297772
##   [26]  0.3999194283  0.6167650489  0.0373203656 -0.0540478503 -0.2843223963
##   [31]  0.0661514864  0.7231751746  1.1016924218  0.4295332900  0.1573273378
##   [36] -0.6137699938  0.4430140973  1.1838519024  1.0178683042  1.0212936694
##   [41] -0.7940797335  0.3592786613  1.1570404724  1.5147526291 -0.0626573107
##   [46]  0.8278364536  0.0926914043  0.2750514897  0.7744017507  0.5790789133
##   [51]  0.4810460286  0.4391883392  0.2562864802 -0.2309776749  1.2201273245
##   [56]  0.5269914298  1.0412251182  0.6863465014  0.9844520959  0.6995299584
##   [61]  0.6420348784  0.6525527016 -0.6982662793  0.6293394674 -0.0481076615
##   [66]  0.5013267856  0.2777814685  0.3791750728 -0.3699381238  0.0676772327
##   [71]  0.8549313411  1.3576499240  0.7898750744  0.0095034679  0.4646620754
##   [76]  1.5986997783  0.1670818460  1.3438413174  0.0700730058  0.1123875696
##   [81]  0.2656780148  0.9367392925  0.4641503963  0.9447470266  0.3423607184
##   [86] -0.0335878821  0.8598236445 -0.1494598981 -0.2119875337  1.2867957225
##   [91]  0.6092653534  0.5594022949 -0.1005566694  1.2677284115  1.3100869781
##   [96]  0.7353488665  1.3788381630 -0.4149332820  0.5752915949 -0.2350610350
##  [101]  0.4332237610  1.0443223124  0.1387780226 -0.3943885631  0.7604390203
##  [106]  0.7292792064  0.6579792541 -1.1577997707  0.7975070170  0.0531097638
##  [111] -0.0926603907  0.3997433278  1.9244003654  0.6322342438  0.2294276471
##  [116]  0.4257075094  0.8159331605 -0.2944511353  0.4831355697  0.6627447202
##  [121] -0.0265865387  1.2980202734 -0.1007318329  0.6332690734  0.8492311026
##  [126] -0.1488892441  1.8016688688  1.3619967211 -0.2133084487  1.1248705770
##  [131]  0.6736371824 -0.2502678645  1.3917222200  0.2045915295  1.9130068706
##  [136]  0.8887502652  1.4254269002  1.1961224427  0.3868935286 -0.1204736022
##  [141]  0.0103406398 -0.2804506966  0.7635645451  0.6470529121  0.7425818793
##  [146] -0.8743106762 -0.0823509070  1.2593180156  0.2173299668  0.4683948796
##  [151]  1.1519303735  2.2476906738  0.4412581175 -0.1455146733 -0.2266239750
##  [156]  0.0456870091  0.1689259159  0.4256066086  0.6885134970  0.1213225566
##  [161]  0.2994977776  0.4799316749  0.0188954387 -0.5066921610  0.8191385239
##  [166]  0.0702192428 -0.6414765346  1.2135674801  0.2378389910  0.2613658274
##  [171]  1.5143033882  1.1524146658  0.6596057104 -0.1466254579 -0.3428357347
##  [176] -0.4944956351  1.1862732420  1.0594787052 -0.0797800386  0.0482954499
##  [181] -0.0836066927  0.5074026557  1.3524337630  0.1050854612  0.7517042882
##  [186] -0.7844545863  0.0864985595  2.1150129359 -0.1552428027  0.8921353682
##  [191]  0.8526798458  0.1896911633  0.3442373427  0.4857586679  0.3296811873
##  [196]  0.0162700336  0.1945725495  0.4344197998  0.6146718838  0.2012476947
##  [201]  1.0245554474  0.0088743622  0.4125740394  0.6128175156  0.5654592079
##  [206]  0.8100934088  1.5423904947  0.4542754858  1.5205073805  1.1376227243
##  [211]  0.1626372462  0.5091330955 -0.1142387630 -0.0375154651  0.6624065152
##  [216]  0.4195920197  0.8838244408  0.8282687338 -0.0349636696  1.1567727999
##  [221] -0.3032021771  1.0274707243 -1.0734018270  0.4010222715 -0.4010485525
##  [226]  0.8561484032 -0.3790693001  1.0958290350  0.9298384182  0.8550770601
##  [231]  0.3960751551  0.1854232924  0.7491849729  1.4438935622  0.8945833861
##  [236]  0.8300383381 -0.0761817576 -0.3428357347  0.5873376817  0.0194365016
##  [241]  0.5387655121  1.0504538503  0.1038738148  0.7367151939  0.7740297415
##  [246]  0.6983220272 -0.3207402871  0.1453280919  1.1772593277  0.1630933205
##  [251]  1.1788142541  0.5140067794 -0.1626150275  1.5039074667  1.1304908809
##  [256]  1.6934904947 -0.5558749015 -0.0042424638  0.2439498353  0.3408939162
##  [261]  0.9053489637  0.5573911456  0.8564007051  0.6142925006 -0.1374803957
##  [266]  0.4805539007  0.9290825730  1.2995097598 -0.1491816991 -0.5896725315
##  [271]  0.3345829500 -0.0711864387  0.9059893882  1.0190257590  0.1330881032
##  [276]  1.0836484843  0.0868160970  1.3619619308  0.0973767862  0.8853286636
##  [281] -0.2881830560  0.1828819549  0.4941516204 -0.3665027877  0.0747763355
##  [286] -0.6812001204  0.4741577612 -0.1162385277  0.2668723425  0.7870853983
##  [291]  0.9602688871  0.4097197656  0.4569489691  0.7644599945  0.4106876682
##  [296]  1.1036516233  0.8379942393  0.1216344638  0.1013965365 -0.2760631531
##  [301]  0.7565058127  0.0539526956  1.0886324947  1.6395165248  0.6850189310
##  [306] -0.1126938562 -0.1608862298  1.3724461818 -0.5828578381 -0.3311050203
##  [311] -0.0009801093  0.3370014443  0.8536977310 -0.1059191137 -0.2477897131
##  [316] -0.0088191721  0.9580400779  0.2429664262 -0.1524790860  1.5125221541
##  [321] -0.2050718016  0.5806701946  0.6955259246 -0.6544301465  0.0248490846
##  [326]  0.9553069004  0.0796694447  1.3696910573  1.2574977836  0.1923233076
##  [331]  0.1319730630  0.7081598384 -0.5541301243  0.8216082554  0.6717966054
##  [336]  0.9528306868  0.9189385880  0.6717463226  0.5958999328  0.2077823021
##  [341] -0.6192398866  0.6924647322  0.4323312299  0.8579700353  0.7447049068
##  [346] -0.1210211863 -0.2288972471  0.1830442118  0.1733769422  0.5857283269
##  [351]  0.4217926249  0.5210227154  0.0314250649 -0.5756333504  1.1440239158
##  [356]  0.7802345985  0.4999260626  0.8580495834  0.8811148329  0.3842433999
##  [361]  0.4162497467  0.6172771997  0.2155016409  1.2265717255 -0.3981676972
##  [366]  0.2866018941  0.6637393282  0.0768770837  0.1377367068  0.4011603292
##  [371]  0.2958055680  1.5783762257 -0.3881629715  1.1372568543  0.5685316443
##  [376]  0.6043296428  1.3917122058 -0.5602876981  0.7410094823 -0.0560474749
##  [381]  0.4337054377  0.8461810423  0.2439458906  0.0760220273 -0.6953534305
##  [386] -0.3806006509  0.5976010285 -0.0066043135 -0.1811562542  0.7175270929
##  [391] -0.0045226339  1.1292298887  1.3151865796  1.4230397112  0.3209326074
##  [396]  0.3894877829  0.6098639227 -0.3440495417  0.0597076986  0.2109079193
##  [401]  0.8858007898  1.8827526361  0.7760282252  0.0597588070  1.7459080005
##  [406]  0.0332971487  0.9560326731 -0.1514452741  1.1506641170  0.2951427869
##  [411] -0.1733464847  0.7337134639  0.6349542263  0.6466043584  1.0571571177
##  [416]  0.2922458605  1.2376697012  0.4885918039  0.8613638041 -0.2581628870
##  [421]  1.0048688471  1.0686705905  0.3690405127  0.0242744499  1.2361443585
##  [426]  0.3345758194  1.2382242601  0.5654847288  0.2578940765 -0.7716655768
##  [431]  0.5289494602  0.3091622084 -0.0567476812  1.1052407032  0.3517885660
##  [436] -0.3357549023 -0.2195433962  1.0141050889  0.2094484178 -0.4208687410
##  [441]  0.5991464774  0.1509582836  0.4806686157  0.5904708539  0.0219489170
##  [446]  0.3603806628  0.0682972083  0.4268480484  0.6466445415  0.2306228746
##  [451]  0.4734048318  1.0429696276 -0.0216828618  1.2821578944  0.9977234759
##  [456] -0.3291211472  0.2449033999 -0.2804999722 -0.3030716965  0.2053600061
##  [461]  1.0299308106  0.5652859172  0.8091135112  0.4544364426 -0.1785038269
##  [466]  0.1552429783  0.0750770209  0.4132026268  0.1880831683  1.0451776259
##  [471] -0.0992543938  0.0304329409  0.1349143822  0.6605563500  0.4103183383
##  [476]  0.5463281249  0.6993249084  0.0976001077 -0.3590686616  1.6762863419
##  [481]  0.6329568227  0.3822448892  1.6664757784 -0.3110197946  0.3869336833
##  [486]  0.5779181678  0.4541419233 -0.1998073120  0.2762884022  1.3885840474
##  [491]  0.6050073265  0.5638591652  0.4543040996  0.5210654426  0.3950159213
##  [496] -0.1708686000  0.9623517224  0.3717051449  0.0175710653  0.0002075923
##  [501]  0.0935082167 -0.4024274141  0.5776768502  0.3236521225  0.2999928875
##  [506]  0.8315042650  0.4043991377  0.3167140639  0.8128071120 -0.9558606296
##  [511]  1.1219842077  0.6336357662  1.1874337750 -0.5798444207  0.1915194299
##  [516]  0.5341446207 -0.2660583209  0.0610415591  0.3513853145 -0.7263154441
##  [521] -0.0990773820  0.9514708417  0.3846729463  1.3187443088  0.1749687763
##  [526]  1.1464403827 -0.6813046280  1.0854094949  0.0593477908  1.0304959974
##  [531] -0.6507340634  1.0153054579  0.3637111964 -0.0576855190  0.7736071346
##  [536]  1.2438563567  0.5876997572  0.6636659948  0.1287596052  0.3496876582
##  [541]  0.3322526707  1.1681026488  0.2268108871 -0.5737518885 -0.0530856225
##  [546]  0.7983068507 -0.1313767905  0.7970491014  0.5222655910  0.4807313551
##  [551]  0.6749344928  0.6733336042  0.0786002540  0.3626368759  0.7357289849
##  [556]  1.0478839796  0.2046968553 -0.1285599364  0.6877205578  0.7293830874
##  [561]  0.2192679102  0.0432574750  1.0529545016 -0.2486121749  0.9206766289
##  [566]  0.1181055202  0.2220751189  0.3447782882  0.5040534441  0.0554772257
##  [571]  0.5527047359  1.1373813010  0.6436273683  0.6018539495  0.6973034419
##  [576]  0.8821140711  2.2086665532  0.2174976462 -0.0345908635  0.8128071120
##  [581] -0.2954316631  1.3968796931  0.4690193142  0.7897201017  0.8621051590
##  [586]  0.4854756804  0.3261378504  1.0477801957  0.1201006608  0.7150725193
##  [591]  0.2192679102 -0.3545530488 -0.5029368753  0.1177283807 -0.1596347669
##  [596] -0.6956040407  0.2624188012  0.9777278270  0.7117081562  1.5692743255
##  [601]  0.2831804894  0.2165565255  1.0419578237  0.5886271850  0.9442722110
##  [606]  0.0957066468  0.1283229881  0.4164775640 -0.0035358955  0.4534361991
##  [611]  0.5950542642  0.8161349765  0.8326575674 -0.3513642762  0.1564680698
##  [616] -0.0381455502 -0.4577224828  1.3553720005 -0.1549088566  0.7614730475
##  [621] -0.2387423636 -0.5355680104  0.3207138800  0.4210931332  0.7132226966
##  [626]  0.7984881109  0.1168153640  2.2086665532  0.4135578646  0.4959598794
##  [631]  0.4028050632  1.0053976022  0.6325582879  1.5509376263  0.4272499680
##  [636]  2.0399198448 -0.1368804632  0.5856937744  0.0283461061  0.5471934760
##  [641]  0.0213723180  0.6266747903  0.6926127525  0.6167650489  0.4993716885
##  [646]  0.5737033163  0.4091616366 -0.1059256336  0.5798254645  0.4609386417
##  [651]  0.4415253539  0.7384102781  0.3967846409 -0.4577874916 -0.1436122874
##  [656] -0.4632218347  0.7226457767  1.0288202979  0.3465200820  0.4158879931
##  [661] -0.1730473710  0.5713250681 -0.3587476196  0.9075904507  0.3057571005
##  [666] -0.0913409859 -0.1514657539  0.8609006023 -0.1250108941  0.6397074429
##  [671]  1.6515224892  0.7661159045  1.2878499706  0.7580762497  0.4844657456
##  [676]  0.0217521676  0.5849978563 -0.0573965079  0.5655436654 -0.4032289594
##  [681] -0.2866566241  0.2340099530  0.8634493158  0.5102282923 -0.3410792176
##  [686]  0.2956173363  1.0110007829  0.2826346546  0.3458151790  0.3638153721
##  [691] -0.0987123586  0.8635156859  1.0536344424  0.9708832079  0.3890876259
##  [696] -0.0027469888  0.2208832362  0.7821852253  0.2438861978  0.9570465839
##  [701]  1.0584122285  0.3158183466  0.5856706994  0.5790429792  0.2586429522
##  [706]  0.5471768418  0.6228552851  0.2232081127  0.9794354153  0.1284714400
##  [711]  0.5832372489 -0.0881433718  1.1044928771  0.3308702187 -0.1882065319
##  [716]  0.4274885956  0.2437277049 -0.2752883923  0.1946305892  0.9022246360
##  [721]  0.9388890490  0.0157626915  0.9958933946  0.2825478941  0.7764898752
##  [726]  0.9716237684  0.1145541042  0.9343283160  0.3103749683  0.2763647935
##  [731]  0.2919094695 -0.1716965215  0.8620737266  0.5041220624  0.6956881361
##  [736]  0.0871194350 -0.7562008578  0.6182767669 -0.5140756940  0.5055829538
##  [741]  0.4894150362 -0.2962131323  0.4011603292  0.9249804922  1.0913278773
##  [746]  0.4671724143 -0.5990452350  1.0682383436  0.4478586148  0.1224982776
##  [751]  1.1004256297  1.1099084898 -0.5251873621  0.2414847607 -0.0590069500
##  [756] -0.4114112046  0.1354464351  0.7746237650 -0.3788097090  0.8881872997
##  [761]  0.9001317724  1.2514060728  0.5681333595 -0.4083596387  0.8831855230
##  [766]  0.3810448178 -0.3665380202 -0.0020100588  0.6595888443  0.8898409613
##  [771]  0.6710529547  1.4665431701  0.5405675527  0.3516004302  0.7309414947
##  [776]  1.2538024927  0.3064087599  0.6602186812  1.1905713772  0.2727936171
##  [781]  1.4873372063  0.5213052980  0.4034780249 -0.0249346630  0.0192679612
##  [786]  0.1873601914  1.7702495597  0.4252542795 -0.5944074501 -0.0389929983
##  [791]  0.0323891478  0.1786528296  0.6069681488  0.3964806508  0.5142905370
##  [796]  0.3607637576  0.2137052382  0.4393888170  0.4339910769 -0.2226686877
##  [801] -0.0024998232  0.4063213748  0.6531310422  0.3705704100  0.7559551692
##  [806]  0.0408743096  0.3011137173 -0.4030061214 -0.3703854117 -0.3633247817
##  [811] -0.4337908140 -0.0558972415  0.2989570101  1.1513392138 -0.3502302670
##  [816]  0.7857627440  0.4609386417  0.6014560920  0.2648421005  0.4554352074
##  [821]  0.4299923179  0.8902238723  0.7823903650 -0.5134919814  0.3868577597
##  [826]  0.0729564399  0.1167646578  0.4600277206  0.3303849354 -0.9616584967
##  [831]  1.0923204030  0.9650446730  0.1385826437  0.7597783734  0.4764560240
##  [836]  0.2218102192 -0.5484772013 -0.2924198075 -0.0520088335  0.0613158337
##  [841]  0.8640125695  0.3653600983  0.6650329996  1.2664877584  0.7629686895
##  [846]  0.4210070880 -0.1698841112  0.4057718927  0.1223803072  0.1699040544
##  [851]  0.2765393075  0.4185987951  1.4589681613  0.4058997449 -0.4782543503
##  [856]  0.4527514749  0.2148956622 -0.5159317363  0.2385827759  0.4260146083
##  [861]  0.2991304963  0.6493825974  0.2623906686 -0.3251116728  1.0652967912
##  [866]  0.9538029571  1.0394110476  0.2345615290 -0.1248359819  0.9244196293
##  [871]  0.4369208800  0.0860712596  0.8527689096  0.5233962211  0.6215431554
##  [876]  0.4789185070  0.0194566828 -0.3293665039  1.8997973680  0.5709609462
##  [881]  1.2189981771  0.0299823361  0.8265257364 -0.0149905403  0.4426514720
##  [886]  0.3654278716  0.2518421887  1.2575190152  0.3198174502  0.8754650979
##  [891]  0.4456649032  0.1796965631 -0.0319165461  0.6153023081  0.1674926100
##  [896]  1.1872471415  0.5929972615  1.1733147110  0.4131683270  0.1443065410
##  [901]  0.7287031763  0.3098101020 -0.1735261134  0.1960662325  0.5755441493
##  [906]  0.5500910467 -0.1684427193 -0.1714092207 -0.2593246099  0.4487398552
##  [911]  0.2097672181  0.6459102906  0.0299331935  0.2593513979  0.8718487845
##  [916]  0.3302851129  0.6833519841 -0.0233934576  0.1275959760 -0.2806333689
##  [921]  0.5203171505  0.0453517092  0.7948234486  0.6885586663  0.3921791575
##  [926]  0.4664383650  0.5341446207  0.2989572735 -0.0238859682  0.6573411826
##  [931]  0.4340189448 -0.2433811080  1.2100736634 -0.4925361293  1.0181114401
##  [936]  0.1904221462  0.5529962388  0.3311055386  0.6089181181 -0.1395349493
##  [941]  1.0368894498 -0.0786324511  0.0532864989  0.3945648733  0.5722144138
##  [946]  0.8126968385  0.4565346389 -0.7135422004 -0.0093077476  0.6909472472
##  [951]  0.6957990507 -0.3481985214  0.5575311945  0.3096917294  0.2726927567
##  [956]  0.2550176082 -0.3411223746  0.1973741703  0.1420198431 -0.5979857859
##  [961]  0.1492755255 -0.1600199888  0.3949941027  1.4658663317 -0.1902974182
##  [966]  1.7807597986  0.3838156559  0.0854393514  0.2344298360  1.2376930912
##  [971]  0.2187409570  1.5032766488  0.4123676319  0.2252911681  0.6517362327
##  [976]  0.0383701562  0.2594670553  0.0282342730  0.1832386670  0.3199151757
##  [981]  1.4567386839 -0.0567374431  0.5065018036  0.2290611909  1.1511568693
##  [986]  0.2540962848  0.2170560945  0.1658713244 -0.0379798497 -0.4206149545
##  [991]  0.9207101740  0.2809299709 -0.0761487983  0.7137382631 -0.0836066927
##  [996] -0.4262607843 -0.1670427618  0.0300300606  0.2263574332 -1.3232534284
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
##   0.41657945   0.29322927 
##  (0.09272724) (0.06556487)
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
## [1]  0.31434260 -0.74705582 -0.05282655 -0.13679328  0.73314574 -0.63356713
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
## [1] 0.0103
```

```r
sd(results2$thetastar)  #the standard deviation of the theta stars is the SE of the statistic (in this case, the mean)
```

```
## [1] 0.8882076
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
##     original       bias    std. error
## t1*      4.5 -0.003203203   0.9145724
```

One of the main advantages to the 'boot' package over the 'bootstrap' package is the nicer formatting of the output.

Going back to our original code, lets see how we could reproduce all of these numbers:


```r
table(sample(x,size=length(x),replace=T))
```

```
## 
## 0 1 3 4 5 7 8 
## 1 3 1 1 1 1 2
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
## [1] -0.0316
```

```r
se.boot
```

```
## [1] 0.9097852
```

Why do our numbers not agree exactly with those of the boot package? This is because our estimates of bias and standard error are just estimates, and they carry with them their own uncertainties. That is one of the reasons we might bother doing jackknife-after-bootstrap.

The 'boot' package has a LOT of functionality. If we have time, we will come back to some of these more complex functions later in the semester as we cover topics like regression and glm.

