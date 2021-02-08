Week 2 Lab
=============

Confidence intervals
-----------------------

Before getting too far, we need to circle back and make sure we understand what is meant by a confidence interval. 

A 95th percentile confidence interval say “If I repeat this procedure 100 times using 100 different datasets, 95% of the time my confidence intervals will capture the true parameter”. It does NOT say that there is a 95% chance that the parameter is in the interval.

**Quiz time! (Don't worry, not a real quiz)**

*Important note*: This is an area where Aho is WRONG. I will not repeat Aho's interpretation here because I think he's just wrong. Aho is correct on only one point. It is true that ONCE THE 95th CI HAS BEEN CONSTRUCTED, it is no longer possible to assign a $%$ to the probability that that CI contains the true value or not. Because that CI, once created, either DOES or DOES NOT contain the true value. However, we often talk about the interval in the abstract. When we say "There is a 95$%$ chance that the interval contains the true value" what we mean is that there is a 95$%$ probability that a CI created using that methodology would contain the true value.

Do not let Week 2 pass by without fundamentally understanding the interpretation of a confidence interval. 

Testing hypotheses through permutation
------------------------------------

We'll start off by working through two examples from Phillip Good's book "Introduction to Statistics Through Resampling Methods and R/S-PLUS":

**Example #1**: Use permutation methods to test the null hypothesis that the treatment does not increase survival time (in other word: $H_{0}$: No difference in survival between the treated and control groups):

Survival.treated=$\{94,197,16,38,99,141,23 \}$

Survival.control=$\{52,104,146,10,51,30,40,27,46 \}$

(Is this a one-tailed or a two-tailed test?)

Make sure that you understand what is being done here, as this example is very closely related to the problem set.


**Example #2**: Using the same data, provide a 95% confidence interval for the difference in mean survival days based on 1000 bootstrap samples

Note that these two approaches are very closely related. Do you see why either approach can be used to test the null hypothesis? (What is the null hypothesis here?)

Now we will do one slightly more complicated example from Phillip Good's book "Permutation tests: A practical guide to resampling methods and testing hypotheses":

Holmes and Williams (1954) studied tonsil size in children to verify a possible association with the virus \textit{S. pyrogenes}. Test for an association between \textit{S. pyrogenes} status and tonsil size. (Note that you will need to come up with a reasonable test statistic.)

<div class="figure" style="text-align: center">
<img src="Table2categories.png" alt="Data on tonsil size and S. pyrogenes status. Source: Good (1994)" width="40%" />
<p class="caption">(\#fig:unnamed-chunk-1)Data on tonsil size and S. pyrogenes status. Source: Good (1994)</p>
</div>

Now lets consider the full dataset, where tonsil size is divided into three categories. How would we do the test now? What is the new test statistic? (There are many options.) What 'labels' do you permute?

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
## 0 1 3 6 7 8 9 
## 1 1 2 1 2 1 2
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
## [1] -0.0157
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
## [1] 2.748413
```

```r
UL.boot
```

```
## [1] 6.220187
```

Method #2: Simply take the quantiles of the bootstrap statistics


```r
quantile(xmeans,c(0.025,0.975))
```

```
##  2.5% 97.5% 
##   2.7   6.2
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

Does it make sense why the normal distribution theory intervals are too wide? Because the original were were uniformly distributed, the data has higher variance than would be expected and therefore the standard error is higher than would be expected.

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
##    [1] 4.7 5.4 3.9 3.0 4.8 4.3 3.5 5.8 5.2 3.7 2.6 4.8 4.7 4.9 3.5 4.9 2.9 4.6
##   [19] 4.0 1.5 6.5 4.2 4.5 4.4 5.3 5.2 3.3 4.9 4.8 4.2 3.4 5.7 6.3 5.5 3.6 4.6
##   [37] 3.9 4.9 4.2 5.0 4.5 6.1 5.2 4.4 4.9 6.1 4.7 3.4 5.4 4.8 5.3 5.0 3.9 4.3
##   [55] 3.7 3.3 4.7 4.0 4.9 4.6 5.3 3.0 5.7 4.5 3.9 6.6 3.7 4.5 5.1 5.2 5.6 3.6
##   [73] 4.0 4.7 5.8 4.2 4.1 5.4 3.5 4.4 2.8 4.2 4.9 4.4 5.1 4.3 5.5 3.8 5.6 4.4
##   [91] 3.0 4.1 5.3 4.0 5.6 5.8 5.3 5.3 4.5 4.3 5.4 3.6 4.3 4.0 5.4 3.2 5.8 4.7
##  [109] 4.8 5.5 4.7 3.9 6.3 4.4 2.7 4.7 4.6 1.5 4.1 3.2 5.1 3.5 5.2 3.1 5.9 3.9
##  [127] 5.2 3.5 3.2 5.4 5.8 4.8 5.0 5.1 6.7 4.4 4.4 3.8 4.5 3.4 4.9 6.0 5.1 5.0
##  [145] 3.5 6.1 3.7 5.0 4.4 4.2 4.4 4.7 4.9 5.5 4.4 4.9 5.4 4.1 4.6 5.0 5.4 3.6
##  [163] 4.1 4.0 4.5 3.6 5.0 4.5 4.0 4.3 3.2 4.5 4.7 4.0 3.2 4.8 2.2 4.2 4.1 4.2
##  [181] 4.4 4.6 5.1 4.4 4.8 5.6 2.5 4.9 4.7 5.1 4.9 5.0 6.7 4.7 5.2 4.3 3.6 4.8
##  [199] 3.5 3.7 2.3 4.4 5.4 4.8 4.2 4.5 3.8 5.9 6.0 5.6 5.3 4.1 3.7 4.7 3.9 3.5
##  [217] 4.5 4.0 5.0 4.6 4.9 3.9 4.7 3.4 4.6 3.6 3.6 4.4 4.0 5.0 4.2 5.8 3.7 3.7
##  [235] 4.4 5.8 5.8 4.3 5.0 5.4 4.8 4.0 4.8 5.5 4.5 7.0 4.5 4.2 4.3 4.5 4.8 6.1
##  [253] 5.1 3.5 3.9 4.1 3.5 4.7 7.4 4.8 5.6 4.8 3.0 4.2 5.6 4.8 5.3 2.3 4.3 3.0
##  [271] 4.6 3.7 4.8 4.5 4.4 5.8 3.2 4.2 5.1 5.1 3.7 5.3 5.1 5.1 3.6 3.9 4.4 5.1
##  [289] 5.6 4.6 4.9 3.7 3.9 4.3 4.0 5.4 5.4 4.4 5.6 5.3 4.8 5.6 4.3 4.1 5.5 5.0
##  [307] 5.8 3.1 4.7 3.2 4.1 4.4 4.8 3.9 4.2 4.9 2.3 5.9 5.1 4.3 3.1 4.7 3.9 4.7
##  [325] 3.7 4.9 3.6 5.0 3.7 3.6 3.2 4.5 2.7 4.8 4.9 4.4 5.4 5.3 3.5 4.7 4.8 5.4
##  [343] 3.9 3.8 4.4 6.1 3.1 4.4 3.7 5.2 5.1 5.1 4.7 4.4 5.7 4.9 5.1 4.2 5.1 4.1
##  [361] 4.2 5.1 3.7 4.2 4.5 5.5 5.2 5.6 5.4 5.5 3.7 5.4 5.6 3.7 2.3 6.1 3.7 4.6
##  [379] 3.7 4.3 4.8 4.9 4.2 5.1 4.6 2.3 3.4 5.6 4.5 4.4 3.9 4.8 4.1 4.1 5.0 4.7
##  [397] 4.0 5.5 4.0 4.3 4.3 5.6 4.6 4.3 4.1 3.3 4.7 4.2 3.5 3.5 5.2 5.1 5.0 4.1
##  [415] 4.1 5.0 3.2 3.0 5.3 5.1 5.1 4.5 5.2 3.0 5.6 5.3 3.5 4.1 4.8 4.1 2.7 3.1
##  [433] 5.2 4.3 2.0 3.6 5.1 4.7 3.1 4.0 5.3 4.0 3.7 5.6 4.7 6.0 5.5 5.4 3.8 2.6
##  [451] 6.5 3.8 5.9 4.6 4.6 4.8 5.6 5.9 4.2 4.6 3.0 4.8 5.5 4.2 4.2 4.1 4.5 5.7
##  [469] 5.0 4.1 4.5 6.1 5.6 5.4 6.2 4.7 4.3 5.5 5.9 4.4 4.9 4.1 3.6 5.3 4.9 5.6
##  [487] 4.3 5.1 3.9 5.1 4.7 4.4 4.3 5.4 3.6 5.5 5.2 4.5 4.2 3.6 4.1 5.8 5.3 4.5
##  [505] 4.9 4.5 4.6 3.6 7.2 3.7 4.1 4.6 3.9 5.0 5.5 5.8 6.4 5.3 4.1 4.4 5.5 5.4
##  [523] 5.3 6.3 2.7 4.4 5.3 4.0 4.9 4.3 4.6 6.0 3.0 3.0 5.7 3.8 5.3 5.2 4.4 3.4
##  [541] 4.3 4.1 4.7 4.6 4.1 4.3 3.4 4.1 5.1 3.3 5.9 3.7 4.2 2.5 5.9 4.6 4.0 5.1
##  [559] 5.9 4.6 5.7 2.9 5.1 4.2 4.9 5.9 4.6 5.4 5.6 3.0 4.1 4.9 4.7 4.1 5.3 3.8
##  [577] 4.8 2.6 4.9 5.1 7.2 4.4 3.0 4.2 3.1 5.9 5.5 3.9 4.2 4.8 5.4 3.9 5.2 3.9
##  [595] 3.8 4.7 3.1 5.1 4.2 5.5 3.7 3.6 3.8 4.2 4.2 4.0 4.0 5.9 6.7 3.1 4.8 3.7
##  [613] 5.5 5.0 4.5 4.5 4.4 4.4 3.6 5.0 4.3 5.3 4.9 3.2 6.2 4.6 3.5 4.9 5.0 3.0
##  [631] 5.2 4.8 5.2 3.5 4.1 4.4 4.7 4.7 4.8 4.4 4.0 4.1 4.5 6.5 4.6 4.0 2.1 4.6
##  [649] 4.5 5.1 5.1 4.2 4.7 4.1 4.7 3.2 2.9 6.5 4.8 4.8 3.2 5.5 2.9 4.6 4.4 4.8
##  [667] 3.4 3.5 3.8 5.1 3.2 4.5 4.2 6.0 5.5 4.5 3.5 5.1 5.7 4.3 4.8 6.1 3.3 2.9
##  [685] 4.0 3.6 5.4 5.3 5.5 5.0 3.4 3.2 5.0 4.1 5.9 3.7 3.9 3.3 3.9 4.2 3.0 3.5
##  [703] 3.9 5.7 5.7 4.7 5.1 2.7 5.4 4.2 4.3 4.3 5.1 5.4 6.6 4.0 3.6 2.8 4.9 3.9
##  [721] 4.6 4.4 4.9 4.7 3.6 2.9 4.1 3.2 5.7 2.4 4.5 3.3 5.3 3.9 3.9 4.5 4.1 2.7
##  [739] 4.9 5.1 4.7 4.0 4.1 4.8 5.2 3.0 4.6 4.8 4.0 3.6 5.2 3.6 6.0 3.0 4.5 4.5
##  [757] 5.7 5.2 4.0 3.7 6.1 4.8 4.2 5.8 5.3 4.5 4.3 5.1 3.1 2.5 4.9 4.8 4.6 6.3
##  [775] 6.3 5.6 4.0 4.7 4.9 5.3 4.1 5.0 5.2 3.5 4.1 5.8 4.3 3.2 5.3 3.2 4.8 4.7
##  [793] 5.0 4.8 3.6 5.1 4.8 4.8 5.1 2.5 3.3 5.0 2.9 5.6 4.4 5.2 4.9 4.5 4.4 4.0
##  [811] 6.0 4.2 5.5 4.0 5.3 3.6 4.0 5.3 3.6 2.5 5.2 4.8 3.8 2.9 5.3 5.4 3.5 3.3
##  [829] 4.7 3.3 3.5 4.3 3.0 4.4 3.4 5.9 4.4 6.0 6.2 2.2 4.6 4.0 3.9 5.6 3.5 3.8
##  [847] 5.9 5.1 3.9 5.9 4.9 3.5 3.3 5.1 3.7 5.6 3.9 4.8 5.8 2.9 6.1 5.2 5.3 5.3
##  [865] 3.9 5.7 4.0 3.9 3.7 5.5 3.7 3.2 4.0 5.8 5.7 5.6 4.3 4.6 6.1 5.3 3.6 3.9
##  [883] 6.1 4.7 3.4 4.8 3.7 2.8 4.6 4.3 4.1 4.6 5.0 4.7 4.6 3.6 4.9 5.5 4.3 3.1
##  [901] 4.3 5.2 3.1 3.5 5.6 4.5 3.5 4.6 5.9 2.0 6.0 4.9 4.4 4.3 3.9 5.4 4.8 4.8
##  [919] 5.4 4.3 3.5 4.5 2.7 7.0 5.4 3.6 3.7 3.5 3.1 3.8 4.2 6.5 4.3 4.7 5.4 4.9
##  [937] 5.5 3.6 2.8 4.4 3.7 3.7 4.1 5.7 5.5 5.4 3.5 4.8 3.9 4.0 4.4 4.1 4.3 4.9
##  [955] 3.6 3.6 4.0 4.0 4.0 6.3 5.6 4.0 4.2 3.3 3.4 5.4 3.8 3.9 4.8 5.5 4.3 6.8
##  [973] 4.2 5.5 4.3 4.7 5.6 5.0 5.3 5.0 3.7 4.2 4.0 5.8 3.8 6.0 5.0 5.4 4.9 5.0
##  [991] 3.4 5.3 5.2 3.8 2.6 3.6 5.0 4.1 5.3 4.9
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
##   2.7   6.2
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
##    [1] 2.6 4.1 3.6 4.9 4.7 5.7 5.6 3.9 5.2 3.1 4.4 4.8 6.2 4.1 4.0 5.7 3.9 5.6
##   [19] 5.6 5.2 3.5 4.3 5.3 4.1 4.0 4.4 4.8 6.0 4.4 5.3 4.6 3.8 5.5 3.7 2.8 6.8
##   [37] 3.8 6.1 3.7 4.2 5.7 2.7 2.8 4.2 3.7 5.3 2.4 6.7 4.3 4.8 4.4 7.2 3.3 4.1
##   [55] 4.5 4.0 4.2 4.8 4.9 4.4 3.8 3.1 2.7 3.6 4.2 3.5 3.9 4.8 4.6 3.6 2.9 4.8
##   [73] 4.7 4.8 3.9 4.1 6.3 4.5 3.4 6.7 4.6 5.4 4.4 5.6 3.9 4.4 5.0 3.8 3.2 3.6
##   [91] 5.7 3.3 5.3 5.5 5.6 4.9 4.9 5.4 4.2 5.1 5.9 3.9 3.7 3.9 5.2 5.3 5.3 3.9
##  [109] 4.2 4.8 3.0 4.7 4.4 3.0 5.3 3.9 4.0 5.5 4.3 4.7 4.8 2.7 3.7 7.2 3.5 3.3
##  [127] 5.8 5.7 4.0 5.2 4.4 3.5 3.1 4.5 3.7 5.5 6.3 5.0 4.8 2.8 3.7 3.9 5.2 5.4
##  [145] 4.2 4.9 4.2 6.0 4.2 3.9 4.7 4.2 4.8 4.4 4.6 5.3 4.6 5.0 5.2 5.9 3.5 5.2
##  [163] 2.9 4.0 5.1 4.5 5.2 5.9 3.9 2.3 4.1 3.4 3.4 4.4 4.5 4.8 4.1 5.1 5.9 5.4
##  [181] 4.7 5.1 5.2 4.8 4.1 3.7 4.4 5.4 4.5 3.5 4.2 5.0 6.0 3.9 3.2 5.4 4.7 4.1
##  [199] 6.2 3.7 4.8 4.2 3.3 4.7 3.2 3.3 4.9 4.0 5.3 4.9 4.0 3.3 5.7 5.0 4.7 3.5
##  [217] 3.8 3.9 4.5 4.2 4.7 5.0 6.1 4.0 6.0 5.7 4.5 4.4 3.0 4.8 5.8 4.4 5.4 4.7
##  [235] 3.3 4.5 4.2 4.7 5.0 4.0 4.5 3.6 2.0 2.7 4.8 3.5 4.6 5.0 6.2 3.1 4.1 3.9
##  [253] 5.1 4.7 5.3 6.5 5.6 5.6 4.2 5.0 5.1 5.9 4.4 6.1 3.7 3.5 4.9 4.9 5.1 3.6
##  [271] 5.6 5.1 4.8 3.6 3.9 4.8 3.9 2.6 4.2 3.9 3.5 3.5 4.0 5.1 3.9 3.3 4.7 4.4
##  [289] 3.6 3.0 5.6 4.4 3.4 4.2 4.9 4.1 5.3 5.3 3.1 2.8 6.5 5.2 3.1 3.8 4.0 4.5
##  [307] 4.6 6.2 4.8 4.1 3.5 4.1 4.7 5.2 4.7 5.6 4.8 5.0 5.5 4.4 3.7 4.2 5.3 4.2
##  [325] 5.4 6.2 4.2 1.8 4.7 4.6 4.5 5.2 3.7 2.2 4.4 4.0 4.8 4.7 3.7 6.9 3.6 3.6
##  [343] 5.2 4.7 4.9 4.2 4.7 4.3 5.3 4.9 5.7 3.2 4.5 6.1 5.3 4.0 2.0 3.1 4.7 5.6
##  [361] 4.0 3.8 3.4 4.5 4.7 4.0 6.0 4.1 4.5 5.4 3.3 6.5 5.6 5.3 5.0 5.7 4.4 4.2
##  [379] 4.6 5.3 5.3 4.7 5.3 5.4 4.1 5.8 4.1 4.3 5.2 4.7 4.3 5.2 5.0 5.0 5.9 4.2
##  [397] 4.1 3.4 4.2 5.8 5.0 4.2 5.2 4.9 4.4 4.5 3.5 3.5 4.4 5.3 5.9 3.5 4.4 2.9
##  [415] 5.1 7.1 2.4 4.0 5.6 3.1 4.2 4.9 4.7 5.1 4.4 5.5 5.0 4.2 5.0 4.9 3.9 3.9
##  [433] 3.7 4.6 3.2 4.6 2.7 5.5 5.7 5.0 5.9 5.0 5.8 4.8 3.1 3.9 3.5 3.2 2.9 4.0
##  [451] 4.4 5.5 5.3 5.6 3.1 5.8 5.5 4.3 5.1 5.1 3.6 3.5 4.6 3.4 4.9 4.4 4.5 5.0
##  [469] 3.1 3.8 5.9 5.3 4.0 5.1 3.1 5.8 4.0 5.7 5.0 5.6 5.5 5.2 3.9 4.5 5.8 4.2
##  [487] 3.9 4.1 4.7 5.6 5.1 3.2 4.5 3.4 2.8 2.1 4.9 3.6 5.0 4.8 4.0 3.7 3.5 5.3
##  [505] 6.4 5.3 4.3 4.5 4.8 5.2 3.8 3.8 4.4 5.2 4.0 4.9 4.8 5.2 4.4 5.3 4.3 5.1
##  [523] 5.7 5.4 4.6 4.4 5.2 4.1 2.2 4.3 3.3 6.3 5.3 5.0 3.3 6.0 5.3 5.9 4.4 2.8
##  [541] 5.0 3.5 4.2 6.4 3.9 4.7 5.2 3.0 4.5 3.5 5.5 2.5 5.1 5.0 3.7 5.9 4.2 5.5
##  [559] 5.9 4.3 4.1 4.4 6.1 4.2 5.1 6.2 4.4 3.3 4.2 4.2 5.5 4.2 3.5 4.0 3.9 2.1
##  [577] 6.1 4.1 3.0 2.4 4.5 4.5 6.2 6.3 4.4 5.4 2.9 4.3 4.6 3.0 4.0 4.3 2.8 3.6
##  [595] 3.8 4.4 2.7 5.2 4.1 1.9 4.7 2.6 4.0 4.1 4.8 4.2 5.0 4.3 5.0 2.8 4.4 3.2
##  [613] 5.4 4.3 5.9 4.2 5.0 3.3 4.2 3.3 4.8 3.6 6.1 4.3 2.3 4.8 5.8 4.8 4.2 3.5
##  [631] 3.7 4.1 4.0 4.9 4.2 4.1 4.1 3.2 4.4 4.6 5.4 5.0 4.8 6.0 5.5 3.2 3.9 3.9
##  [649] 4.0 4.3 5.2 6.6 4.1 5.1 4.9 5.4 3.0 3.9 5.2 4.3 4.3 5.0 3.2 3.8 5.4 5.1
##  [667] 2.8 2.6 5.0 4.8 4.4 3.6 5.0 5.7 5.4 3.4 5.8 4.7 4.4 3.3 3.9 3.8 4.9 5.8
##  [685] 5.4 3.9 4.3 6.0 6.4 4.4 2.9 5.7 3.9 4.5 4.4 4.9 4.6 5.3 5.5 4.1 5.3 5.3
##  [703] 4.1 4.5 4.0 3.4 3.2 4.1 5.4 3.1 6.1 2.5 4.0 5.9 3.4 3.7 4.1 3.7 5.9 4.4
##  [721] 4.7 5.0 5.2 3.7 5.3 3.9 4.9 3.0 2.7 4.0 4.7 5.8 2.7 4.6 5.4 3.4 3.4 5.8
##  [739] 4.3 5.3 6.1 5.5 2.8 4.7 5.6 4.8 5.2 6.1 3.8 4.8 4.9 5.5 3.8 5.3 4.2 5.3
##  [757] 4.5 4.5 4.1 4.6 3.3 4.1 3.1 3.7 6.0 4.6 4.7 2.7 3.3 5.3 4.6 5.8 5.1 3.0
##  [775] 5.3 4.1 2.6 5.7 4.0 4.2 5.2 3.4 6.5 4.9 4.5 3.3 4.5 4.9 3.3 4.9 3.9 4.7
##  [793] 4.6 5.2 2.3 5.3 3.6 5.0 5.6 3.7 6.4 4.1 2.8 4.8 3.1 4.5 4.8 4.3 4.8 4.9
##  [811] 6.4 4.1 3.6 5.3 3.1 6.3 4.3 4.6 5.8 4.1 3.3 2.4 3.4 4.7 4.7 5.8 4.7 4.0
##  [829] 4.6 5.4 3.3 3.1 4.3 5.2 5.1 5.7 5.7 3.1 4.5 5.3 2.0 4.9 3.5 3.5 4.5 4.3
##  [847] 4.1 3.2 2.7 4.0 3.8 5.6 4.0 4.8 3.1 4.7 5.2 3.9 4.0 5.5 3.4 4.3 4.1 4.4
##  [865] 4.2 3.6 4.7 4.1 4.3 4.5 3.4 3.8 3.7 3.2 4.7 5.9 4.1 2.0 4.9 4.1 3.8 4.4
##  [883] 4.0 2.3 5.4 5.6 4.3 5.4 3.8 4.4 7.0 3.9 6.4 5.2 4.2 3.8 2.0 2.6 3.2 4.5
##  [901] 4.1 4.0 5.5 4.3 3.1 6.5 5.8 3.4 4.1 2.9 5.1 4.9 5.5 3.8 3.9 5.3 4.7 3.3
##  [919] 4.7 4.3 6.4 4.0 2.3 6.1 4.6 3.7 4.1 2.3 5.1 5.1 2.7 4.7 3.5 1.8 4.0 4.8
##  [937] 4.5 5.1 4.0 6.0 2.8 3.4 4.1 4.6 4.0 3.8 4.0 3.9 4.8 6.0 3.6 5.4 4.0 4.1
##  [955] 3.4 5.8 2.7 4.2 4.9 3.8 5.1 4.8 5.0 4.8 3.7 3.9 5.9 4.0 4.7 4.6 3.9 3.3
##  [973] 4.9 4.8 4.7 3.8 3.3 4.5 5.7 4.9 3.0 4.7 3.8 3.4 4.7 3.4 5.8 4.8 5.7 4.4
##  [991] 6.4 4.8 4.2 4.2 4.9 3.6 3.9 4.7 3.8 4.7
## 
## $func.thetastar
## [1] -0.044
## 
## $jack.boot.val
##  [1]  0.48798799  0.36768802  0.27189349  0.19482289  0.03665689 -0.09745763
##  [7] -0.21820652 -0.34687500 -0.48513514 -0.59387755
## 
## $jack.boot.se
## [1] 1.055302
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
##    [1] 3.3 5.5 4.5 4.6 4.6 5.0 3.7 3.7 4.7 4.2 6.2 2.7 4.1 4.4 5.5 5.0 4.6 4.9
##   [19] 4.2 3.5 4.4 5.4 5.8 2.7 4.7 3.0 3.5 4.7 4.2 5.6 6.0 5.1 3.6 5.3 4.3 4.0
##   [37] 3.5 4.4 5.5 5.3 3.9 2.9 4.1 4.4 4.0 3.7 5.1 3.9 4.8 4.1 4.2 4.4 4.7 3.9
##   [55] 4.7 3.6 5.5 4.9 5.3 4.4 4.5 5.0 5.3 3.4 5.7 4.0 4.2 2.5 4.0 5.4 3.9 6.4
##   [73] 4.2 3.6 3.8 5.5 5.8 3.9 3.8 4.8 3.5 2.8 4.0 5.6 3.9 4.6 4.4 4.2 4.2 4.0
##   [91] 5.4 5.9 4.5 4.1 6.3 5.7 4.7 3.0 5.1 3.0 3.5 5.6 4.3 3.4 5.4 4.2 5.1 6.0
##  [109] 3.7 3.3 3.9 4.5 4.8 4.6 3.7 5.7 6.0 4.0 2.9 3.5 2.8 6.6 4.0 4.4 3.4 5.4
##  [127] 3.8 5.2 5.7 4.5 5.8 5.8 4.6 4.0 3.6 4.6 4.6 3.9 2.7 3.6 5.9 6.2 5.2 4.5
##  [145] 4.6 4.9 4.9 4.6 2.9 5.9 3.6 4.7 4.1 4.1 4.1 2.7 5.3 3.3 5.2 4.4 4.7 3.9
##  [163] 4.5 4.1 5.4 5.2 4.4 4.9 3.3 4.7 4.4 4.1 4.3 4.1 5.6 4.5 3.9 5.5 4.6 6.6
##  [181] 3.1 4.6 4.0 6.5 5.0 4.6 3.5 4.6 4.0 5.7 4.2 5.5 5.3 5.0 4.3 3.7 5.3 1.6
##  [199] 4.4 5.8 4.9 6.3 4.3 4.3 3.5 4.9 4.2 3.9 3.7 6.1 5.3 4.5 4.8 3.1 3.8 5.2
##  [217] 4.4 4.6 5.7 4.2 5.3 4.5 4.9 3.9 4.8 4.3 4.8 5.3 3.7 3.4 3.9 5.4 5.2 3.2
##  [235] 5.9 4.8 5.0 4.1 4.8 4.1 5.2 4.4 6.0 5.8 5.4 4.8 5.5 5.3 2.2 5.2 5.2 3.2
##  [253] 3.5 5.0 4.7 4.0 5.4 3.9 4.6 4.1 3.5 3.6 3.1 5.2 2.9 4.9 5.3 4.7 4.8 5.7
##  [271] 3.6 5.5 4.9 3.9 4.3 4.5 5.5 3.4 4.9 6.0 3.0 6.0 4.6 5.1 3.7 3.3 3.4 4.9
##  [289] 5.5 3.5 5.3 4.7 3.4 4.4 4.7 6.0 5.1 4.6 5.3 5.6 3.0 5.0 4.4 5.9 5.5 3.7
##  [307] 5.3 5.1 4.7 6.3 4.2 5.5 4.3 4.7 4.6 3.4 4.1 5.2 3.9 4.5 5.6 5.5 4.9 2.6
##  [325] 5.1 3.9 4.6 4.1 4.8 3.8 3.9 4.9 4.2 3.7 5.3 6.6 4.3 5.9 3.6 5.0 4.1 4.2
##  [343] 4.9 5.3 3.8 3.1 5.2 3.8 3.9 5.6 5.4 4.0 4.6 3.4 3.6 5.2 3.5 5.8 5.9 5.2
##  [361] 2.8 4.0 5.3 3.5 5.5 3.8 5.7 4.9 4.0 3.4 5.2 4.5 5.8 4.5 4.9 3.7 3.7 4.5
##  [379] 4.8 4.2 5.6 4.0 2.7 5.1 4.5 4.5 4.1 5.3 3.4 4.7 2.8 5.0 4.5 5.1 4.9 5.2
##  [397] 5.4 3.3 4.6 4.7 3.6 4.9 5.8 3.9 4.4 4.2 4.3 5.8 4.0 4.5 4.7 2.8 4.6 4.4
##  [415] 6.4 4.3 4.7 4.3 5.2 4.3 5.1 3.3 4.2 4.5 4.8 5.2 4.2 4.9 5.2 4.3 3.4 5.2
##  [433] 3.5 2.9 6.0 4.8 4.7 3.0 3.7 4.1 4.1 5.6 5.6 3.0 4.1 3.1 6.0 3.8 3.6 4.2
##  [451] 4.3 1.8 5.3 6.9 6.7 3.8 5.1 4.1 5.6 2.9 5.3 6.8 4.2 3.3 4.1 5.5 5.7 4.6
##  [469] 4.1 4.5 4.5 4.5 5.1 4.2 5.0 4.9 3.0 5.5 5.1 5.7 6.1 4.7 5.0 3.5 3.4 4.1
##  [487] 6.5 4.2 5.8 4.6 4.6 5.6 4.4 5.1 3.5 3.6 2.6 5.7 5.3 4.2 4.2 5.6 2.8 3.3
##  [505] 4.3 2.4 4.4 3.4 3.6 4.9 4.8 3.8 2.8 4.4 5.6 4.1 4.1 3.7 3.6 6.0 4.2 4.6
##  [523] 3.5 5.4 4.6 5.4 4.2 4.0 4.3 5.3 4.6 5.4 2.7 4.5 4.9 3.8 3.2 4.8 3.9 4.3
##  [541] 5.0 5.1 5.6 4.1 5.3 3.6 4.9 4.0 4.7 4.6 3.2 3.9 4.8 4.6 4.2 4.1 6.3 3.7
##  [559] 3.7 3.7 4.2 5.9 5.4 3.9 4.4 5.9 4.0 4.9 4.6 5.5 5.4 5.6 4.7 4.6 6.4 4.6
##  [577] 4.9 4.1 5.9 4.6 4.5 5.4 4.7 4.9 3.9 4.9 5.1 4.4 5.3 4.7 3.4 4.8 3.1 4.9
##  [595] 4.9 3.6 3.7 4.9 5.4 5.3 6.0 4.1 5.5 4.3 5.7 5.9 3.8 4.4 6.5 3.7 5.2 4.5
##  [613] 4.2 2.9 3.9 5.2 2.9 4.8 4.4 2.6 4.2 4.8 5.4 4.5 3.9 4.2 4.1 4.7 3.8 4.5
##  [631] 3.1 4.6 4.8 4.2 4.0 3.3 3.8 5.4 4.9 4.9 5.8 4.9 2.8 4.0 4.7 4.4 4.1 3.3
##  [649] 4.1 3.6 4.9 3.7 4.1 4.4 4.1 6.6 3.7 3.6 6.1 4.0 4.8 4.5 3.6 3.6 4.6 4.0
##  [667] 4.8 3.9 4.5 5.2 4.8 3.5 4.7 5.7 4.4 5.4 3.8 4.4 4.7 4.6 6.7 3.8 5.7 3.0
##  [685] 4.0 4.9 3.8 4.8 1.9 3.7 3.7 4.1 4.1 2.6 3.8 5.8 4.3 5.1 4.9 4.9 5.3 4.0
##  [703] 5.6 3.6 5.6 5.2 4.7 5.1 3.8 3.9 4.6 3.9 4.5 5.8 4.9 5.2 5.1 4.9 5.9 4.9
##  [721] 3.6 4.1 6.3 5.2 4.4 5.0 4.7 4.8 4.2 4.0 5.1 4.3 6.3 4.3 3.8 4.1 4.0 3.9
##  [739] 4.7 4.2 4.3 5.0 5.6 4.2 6.0 4.6 3.9 5.3 3.4 3.7 5.4 3.8 6.0 4.0 4.5 6.2
##  [757] 2.5 4.6 3.4 7.1 3.7 2.9 4.2 4.5 5.8 4.4 4.8 4.2 4.4 4.5 3.1 5.4 4.8 5.0
##  [775] 5.1 3.2 4.3 4.0 4.9 2.9 3.5 4.2 5.8 5.7 5.3 4.5 5.8 4.7 5.1 5.5 4.7 5.2
##  [793] 5.3 4.7 4.6 3.2 4.3 5.1 4.1 4.6 3.3 3.7 4.5 5.5 6.0 5.5 5.7 4.0 4.4 5.9
##  [811] 2.5 5.1 4.6 6.1 5.9 4.5 3.8 5.4 4.3 4.6 5.1 4.2 2.7 3.7 4.2 4.6 5.7 5.4
##  [829] 4.0 4.0 4.0 3.6 3.9 4.5 5.5 6.1 5.1 5.5 2.9 3.3 3.9 3.0 4.4 4.5 5.6 4.0
##  [847] 5.7 4.6 3.8 4.5 6.6 4.9 3.1 4.2 4.4 3.0 3.8 5.4 4.4 3.9 5.2 4.0 3.8 4.6
##  [865] 3.3 5.0 5.7 3.5 3.2 5.6 3.5 4.6 5.0 3.4 2.3 4.0 3.4 4.8 4.2 4.4 3.5 5.0
##  [883] 4.0 5.0 3.5 3.7 3.9 5.4 5.2 4.2 4.4 4.5 4.1 4.8 3.6 4.2 3.4 3.9 4.8 4.6
##  [901] 5.2 3.4 3.7 3.9 4.1 4.7 4.6 4.9 4.4 4.4 3.5 4.0 4.0 3.9 4.2 5.9 4.4 2.7
##  [919] 3.4 5.2 4.6 4.2 4.0 3.8 5.7 5.1 4.3 2.5 5.3 3.1 5.7 5.3 3.7 4.3 4.3 5.0
##  [937] 4.3 4.2 3.5 4.9 4.8 3.8 3.2 6.4 4.9 6.0 6.6 4.0 3.5 4.7 4.9 5.0 3.3 4.7
##  [955] 4.0 5.4 4.5 4.9 4.3 4.2 3.2 4.9 6.0 4.9 4.5 5.9 6.3 3.7 4.8 4.3 3.3 5.6
##  [973] 3.2 4.8 3.5 6.7 4.1 5.7 4.7 7.3 4.2 5.2 4.9 6.7 4.3 4.6 4.3 5.6 4.3 4.7
##  [991] 5.0 5.0 5.3 3.9 6.0 4.9 4.4 5.6 4.8 5.3
## 
## $func.thetastar
## 72% 
##   5 
## 
## $jack.boot.val
##  [1] 5.5 5.4 5.4 5.2 5.2 5.1 4.9 4.8 4.6 4.4
## 
## $jack.boot.se
## [1] 1.041393
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
## [1] 0.3632603
```

What is skew? Skew describes how assymetric a distribution is. A distribution with a positive skew is a distribution that is "slumped over" to the right, with a right tail that is longer than the left tail. Alternatively, a distribution with negative skew has a longer left tail. Here we are just using it for illustration, as a property of a distribution that you may want to estimate using your data.

Lets use 'fitdistr' to fit a gamma distribution to these data. This function is an extremely handy function that takes in your data, the name of the distribution you are fitting, and some starting values (for the estimation optimizer under the hood), and it will return the parameter values (and their standard errors). We will learn in a couple weeks how R is doing this, but for now we will just use it out of the box. (Because we generated the data, we happen to know that the data are gamma distributed. In general we wouldn't know that, and we will see in a second that our assumption about the shape of the data really does make a difference.)


```r
library(MASS)
fit<-fitdistr(original.data,dgamma,list(shape=1,rate=1))
# fit<-fitdistr(original.data,"gamma")
# The second version would also work.
fit
```

```
##     shape       rate  
##   5.395744   7.719333 
##  (2.342256) (3.511830)
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
## [1]  0.03686457  0.20648673  0.80775623 -0.09380957  0.07250676 -0.45153720
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
##    [1]  0.305076735 -0.672316609  0.568188852  0.645422419  0.791274604
##    [6]  0.862685357  0.717327318  1.079080235  0.294211267  0.856934937
##   [11] -0.164560987  0.804536436  0.847548951 -0.695269019  0.720089607
##   [16]  1.077713130  0.817661558  0.266918652  0.176493428  1.273443001
##   [21] -0.598372877 -0.114019182  0.108460024  0.678721490  0.422253889
##   [26]  0.623303521  0.223957084  0.620024354  0.010289411  0.041447517
##   [31]  0.655381415  0.677274592  0.102307261  0.827076901  0.556350749
##   [36]  0.437663650  0.474530017  0.176505161  0.324222562  0.325207526
##   [41]  1.352116631  0.165257364 -2.068097607  0.590987761  0.431098254
##   [46] -0.164061420  0.061340235  0.465163946  0.082899728  0.619539200
##   [51]  1.027093493  0.606368026 -0.030398500  0.230871621  0.657714287
##   [56] -0.005368417  0.386722335  0.682927782  0.340186069 -0.221265349
##   [61]  1.360827354 -0.017056538  0.584241226  1.104678754  0.913209752
##   [66]  0.412102292  0.708949082  0.851619611  0.490096012  0.276686830
##   [71]  0.813170622  0.009233673 -0.384860292  1.225055280  0.500916058
##   [76] -0.183952347  0.733629100  0.234118711  0.024496136  0.520107755
##   [81] -0.159173898  0.045241616  0.970863644 -0.373352056  0.538382254
##   [86]  0.649397489  0.109195590  0.876293847  0.545712841  0.418252922
##   [91]  1.197308586 -1.035983304  0.345228173  0.545391114  0.791756432
##   [96]  0.640586483 -0.732252238  0.583761994  0.132805921 -0.325958049
##  [101]  0.298065417  0.885280876 -0.069107545  0.510224072  0.045806285
##  [106]  0.335936293  0.318103793  0.405622123  1.000918873  0.451172986
##  [111]  0.261494381  0.784129753  0.124303036 -0.119105117  0.351611784
##  [116] -0.496724970 -0.110426564  0.366791790  0.546419222  0.043966503
##  [121]  0.193271396  0.044776956  0.258032709 -0.038604715  1.020612276
##  [126]  0.769695181 -0.258441615  0.294176441  1.335295912  0.324542180
##  [131]  0.948108506  0.802856643  1.027146052 -0.280663555  0.440747970
##  [136]  0.363747324  0.714744352  0.353492773  0.299639859  0.250444024
##  [141]  1.879499122  1.120421180  0.466949503  0.335560454  0.460830765
##  [146]  0.423796324  0.349076458  0.451889317  0.158136491 -0.156513704
##  [151]  1.128430347  0.137274011 -0.242671912  0.856370743  0.587288682
##  [156] -0.954574894 -0.989704211  0.435633586  0.480713031  0.715803337
##  [161] -0.091128173  1.074563249 -0.135610039  0.544623749  0.764441727
##  [166]  0.171166187  0.114712088  0.012240312 -0.133529086 -0.653751699
##  [171]  0.672405261  0.431472205 -0.196641202  0.474575683  0.236105636
##  [176]  0.185217033  0.305267927  0.126071234  0.769695181  1.053687088
##  [181]  0.501870943 -0.256582209  0.124655984  0.283073904  0.334651778
##  [186]  0.814168131  0.615272319  0.575110534  0.485103666  0.736171996
##  [191]  0.789205258  1.519173390  0.037840601  0.512453760  0.453984756
##  [196] -0.106830029  0.363079456  0.056216306  0.740010608  0.984810287
##  [201]  0.927177064 -0.243021239  0.209863146  0.546686113  0.731875003
##  [206]  0.541314547  1.284439012  0.087247081  0.933171574 -0.297919565
##  [211]  0.663554786  0.743636296  0.123622417  0.580092119 -0.213718964
##  [216]  2.619068898  0.373865505 -0.208434452  0.691493722  0.817036940
##  [221]  0.135549384 -0.728204205  0.028979066 -0.179733983  0.522291674
##  [226]  0.048849380  0.536633420  0.388425318  0.350278785  0.433694082
##  [231]  0.427057509  0.628807470  0.040463472  0.112417040 -0.883200458
##  [236]  0.035422737  1.101576346  0.191108530  1.327136256 -0.630755817
##  [241]  0.707599775  1.256060532  0.009370807  0.123888280  0.283608226
##  [246]  0.684630788  0.650678945 -0.498263884  0.149224952  0.385031612
##  [251] -0.303628612  0.347208802  0.449710859  0.598385364  0.555356084
##  [256]  0.445135747 -0.051689593 -0.336995283  0.474438838 -0.545723405
##  [261]  0.091621946  0.617674774  0.465035459  0.333975476  0.251372583
##  [266]  0.086829987  0.467020991  0.018372531 -0.267309352  0.049576550
##  [271]  0.106842378  0.467191332  0.817620260  1.743917731  1.293900781
##  [276]  0.970857014  0.440412577  0.451480268  0.468051640 -0.165418042
##  [281]  0.587800756  0.609917102  0.141805856  0.340625184  0.789268296
##  [286]  0.541109983  0.497670067  0.051684897  0.645117373 -0.179521599
##  [291]  1.433765500  0.959848760 -0.604179475  0.494090939  0.288584334
##  [296]  0.465961008  0.889566762 -0.316522239  1.471940645 -0.545501963
##  [301]  0.168739996  0.584148625  0.054181295  0.391067574  0.835101047
##  [306] -0.346159731 -0.272160791  1.295767160  0.443075474  0.655034517
##  [311]  0.689874676  0.952006870  0.076859338 -0.036316507 -0.050148908
##  [316] -0.013622065  0.111469248 -0.977953509  0.816103747  0.128860256
##  [321]  0.484557608  0.867368048 -0.264869832  0.798261726 -0.155584063
##  [326]  0.904309918  0.499297377  0.438719183  0.765094549 -0.542905649
##  [331]  0.464576113  0.148807939  0.149140691  0.836097120 -0.026415845
##  [336]  0.473279938  1.668331179 -0.726842583 -0.402747760  0.651839320
##  [341] -0.044751169  0.500669973  0.110706554 -0.180229127  0.776235334
##  [346]  0.717456531  0.165666195  0.466497911  0.221340001 -0.226089038
##  [351]  0.187505918  0.136075006  0.078648367  0.957949657  0.040546727
##  [356]  0.257935392  0.215991424  1.337407255  0.726385644  0.438663726
##  [361]  0.670281036  0.775025276 -0.007859131  0.009336604  1.088277617
##  [366]  1.197511693  0.132155247  0.959128279 -0.596875632  0.637507202
##  [371] -0.070053338  0.764503211  1.077868076  0.130564620  0.822578417
##  [376]  0.301531169 -0.045290906  0.793740693  0.431833019  0.829344039
##  [381]  1.509540364  1.583816541  0.586957693  0.222240316 -0.257925392
##  [386]  0.782138240  0.505333820  0.884376807  0.002765248  1.036401696
##  [391]  1.033471889  0.912570148 -0.279275302 -0.320686266  0.782760733
##  [396]  0.659236247  0.889951912  0.842876834  1.045561537  0.672529109
##  [401]  0.477326706  0.716689460  0.542843077 -0.196249414  1.284810478
##  [406] -1.009279994 -0.323745720 -0.555722232  0.767346732  0.195066570
##  [411]  0.377863973 -0.237612419  0.052884067  0.334585694 -1.140281219
##  [416]  0.546839393  0.101369704 -0.383376337  1.472645193  0.266991178
##  [421]  0.313998482  0.058024103  0.176072910  1.235582456  0.503163643
##  [426] -0.525724248  1.079142792  0.708313410 -0.065556813  0.421870578
##  [431]  0.163404590 -0.276693085  0.420832971  0.161475267  0.393517544
##  [436]  0.762807670  0.435205629  0.067468006 -0.063189632  0.812117714
##  [441]  0.084488069  0.493582752  0.091167159  0.617326367  0.094765466
##  [446]  1.332929534  0.700201529  0.111528320  1.312068508 -0.560811915
##  [451]  0.539187117  0.166100683  0.150205054  0.340074580  0.825925660
##  [456]  0.007161255  0.220904929  0.214491362  0.317201394  0.759465448
##  [461]  0.093586852  0.722093292  0.412843224 -0.161327501 -0.308786620
##  [466]  0.067342842  0.420818193 -0.074231239  0.821574265  0.177937113
##  [471]  0.134516340 -0.014318669 -0.241951577  1.470840094  0.619325686
##  [476] -0.073201654  0.827436318  0.491618104 -0.576022093  0.328818312
##  [481]  0.849638152  0.356126504 -0.100382047  0.115217496  0.315876129
##  [486]  0.819676064  0.661821387  0.662093026  0.298089872  1.060671266
##  [491]  0.113459089 -0.174268114  0.186714254  0.904066203 -0.095896678
##  [496]  0.695106099  1.273882947  0.110959521  1.290893678  0.112636201
##  [501]  0.101039158 -0.330204962  0.478734221  0.065754754 -0.225269449
##  [506]  0.324595848  0.507724585  0.402891890  0.110864260 -0.207413827
##  [511] -0.258987152  0.611540625  0.701242457  0.584303639  0.141503926
##  [516]  0.072886685  0.347050879 -0.687460633  1.067473641  0.863546396
##  [521]  0.952988135 -0.125461256  0.026569942  0.261454330  0.479942446
##  [526]  0.460184238 -0.005238640  0.364389143  0.970202888  0.664255248
##  [531]  0.421699062  0.648406058 -0.027417100  1.162062002 -0.074730392
##  [536]  0.038246020  0.724547160  0.169589775  0.076131570  0.321264180
##  [541] -0.104397177  0.861023899  0.672908214  0.069859519  0.082882349
##  [546]  0.711108775  0.890635714  0.027805613 -0.055986720  0.369298821
##  [551]  0.293404314  0.302130932  0.394030455  0.187500303  0.054101756
##  [556] -0.004258250  0.194651428  0.198148229  0.148400318 -0.161809488
##  [561]  0.398695529  0.449134117  0.611772806 -0.257065843  0.070860097
##  [566]  0.446396209 -0.183670512  0.303269859  0.095083553 -1.006140117
##  [571] -0.608187525  0.840621518  0.590110826  0.259919639 -0.544926655
##  [576]  0.421821149  0.783434243 -0.614307282  1.298482246  0.606399414
##  [581] -0.172447706  0.478426189  0.624041509 -0.710047691 -0.215498025
##  [586]  0.209293091  0.409130650  1.549011483  0.634732462  0.804032096
##  [591]  0.410465193  0.590870123 -1.073778051  0.859012648 -0.224960549
##  [596]  0.158084705  1.249670508  0.209944373  0.395795798  0.098629500
##  [601]  0.736470988  0.043920484 -0.198090064 -0.152849847  0.165373119
##  [606]  0.477326706  1.007095485  0.414373817 -0.408376089  0.464862693
##  [611]  0.355042265 -0.205608628  0.074636716 -0.838528108  0.533740339
##  [616]  0.286053840  0.408190015  0.424782859  0.850857019  0.035095718
##  [621] -0.022926388  0.973648872  0.438088140  0.231522362 -2.494694017
##  [626]  0.906262054  0.477342320  0.405292354  0.396295272 -0.193643006
##  [631]  0.409907912  0.679089545  0.457025677  0.482302224  0.956057753
##  [636]  0.925546549  0.365989995  0.450635213 -0.030398500 -0.272911284
##  [641] -0.451913228  0.142406109  1.253373833 -0.173585051  0.623734323
##  [646]  0.452473887  1.100446294  0.072790598 -1.422309031  0.469993273
##  [651]  0.058032129  0.418252493  0.434532013  0.156846697  1.599664195
##  [656]  0.426908813 -0.391405748  0.630735598 -0.025218090  0.255019528
##  [661]  0.763259259  0.737952416  1.053209092  0.628807470 -0.525024493
##  [666]  0.597923164  0.036754086  0.777250261  0.463883489  0.821654463
##  [671]  0.281108237  0.263177747  0.504611641  0.347815945  0.735156857
##  [676]  0.399415259 -0.009090478 -0.368996329  0.312621413 -0.291350597
##  [681]  0.717515880  0.181045609  1.128626313  0.166056554  1.150586183
##  [686]  0.410280999  0.716989540  0.216932265  0.737149802  0.696725635
##  [691] -0.362283848 -1.391726618 -0.109987357  0.883242941  0.470098914
##  [696]  1.271138778  0.092011813  0.884929097  0.303920519  0.123048385
##  [701]  1.000606959  0.439975636  0.751446888  0.374515406  0.524877046
##  [706]  0.153714378  0.254774852 -0.221265349  0.466238491  0.291254526
##  [711]  1.101475114  0.471923786 -0.128578665  0.091735796  1.054234850
##  [716]  0.540004599  0.449408267  0.937832953  0.465901227  0.904747697
##  [721]  0.415026705  0.282215656  0.717512780  0.568612440  0.095073094
##  [726]  0.822284292  0.836097120 -1.007231667  0.356024704  0.358653429
##  [731]  0.038706462  0.139809236  0.465385716  0.152789736  0.522279104
##  [736]  0.205138083  0.309571322  0.358539426  0.223984270  0.631189268
##  [741]  0.186396533 -0.445845701 -0.137965147  0.260488836  0.414515186
##  [746] -0.039846352  0.179079766  0.492859109  0.468874483  0.241577302
##  [751] -0.264160681  0.522279104  0.099428352  0.477324283 -1.557555678
##  [756]  1.308726650  0.433079177 -0.235763584 -0.373980043  0.591359118
##  [761]  0.536652725  1.161219950 -0.320351192 -0.054229571 -0.539549771
##  [766]  0.562400612  0.882699231  0.988857380  0.250616505  1.165114688
##  [771]  0.446672757 -0.129748196  1.438105418  0.735957361  0.747443354
##  [776]  0.377226779  0.313017875  0.726563694  0.791517131  0.454420273
##  [781]  0.778484294  0.300611674  0.361877909  0.410550933  0.318504842
##  [786]  0.298440547  0.648051641  0.389127809  0.094820301  0.222785077
##  [791]  1.102944111  0.225705571  0.017617772  0.681188610  0.559202222
##  [796]  0.359672019  0.317328573  0.117407513  0.449371047 -0.075482384
##  [801]  0.842132994  0.719648240  0.222785077  0.114706755 -0.050700769
##  [806]  0.603619215  0.179897144  0.595180494  0.396356513 -0.014296960
##  [811]  0.248487697  0.338041919  0.305578530  0.158032651  0.446983612
##  [816]  0.326600597  0.404585027  0.763532928 -0.209538783  0.191573749
##  [821]  0.127939072  0.051963844  0.827023259  0.107411666  0.212811714
##  [826] -0.403024746  0.518753734  0.943704070  0.490096012  0.235243455
##  [831]  0.137761952  0.034974461  0.315668480  0.048656982  0.521654058
##  [836]  0.171499507  0.830973494  0.409882087  0.022885114  0.489295093
##  [841]  0.744408172  0.718035350  0.390178665  1.009042346  0.529022316
##  [846]  0.168176526  0.419983276  1.276396507  0.709392957  0.132538205
##  [851]  0.329222166  0.097989191  0.958838136 -0.542794489  0.893745292
##  [856]  0.514048266  0.353766531  0.240166747  0.821465524  0.269318218
##  [861]  0.630297421  0.484526595  0.223490672  0.339077940  0.642358005
##  [866]  0.106674573  0.455262298  0.234158807  1.901662995  0.783756076
##  [871]  0.457468939  0.207837829 -0.195323874  0.025632770  0.448853123
##  [876]  1.158386629 -0.145411760  0.386519433  0.991091693  0.775713052
##  [881]  0.041079177  0.738233692  0.006269967 -0.282679158 -0.576125176
##  [886]  0.273430244  0.081989766  0.676704832  0.992540587 -0.070673923
##  [891]  0.758179797  0.351438843  0.289500564  0.277791351  0.593149087
##  [896]  0.756632927  0.547482176  0.427587855  0.106247461  0.045536099
##  [901]  0.294110362  0.081216303  0.611319762  0.741378815  1.079100605
##  [906]  0.455804016  0.275439890  0.654633664 -0.732323474 -0.001989239
##  [911] -0.027999664 -0.506983928  1.322911844  0.474845163 -0.019830419
##  [916]  1.114920407  0.304570307  0.471309309  0.233896041  1.786222433
##  [921]  0.135174474 -0.145880068 -0.187392392  0.820565561  0.259420709
##  [926]  1.449976614  0.119406862  0.707785657  0.054919556  0.067072476
##  [931]  0.156846697 -1.559896634 -0.168583277  0.162850747  0.073955253
##  [936]  0.678787512  0.450582646  0.472840262  0.496733125  0.743347934
##  [941]  0.190921285  0.297904248  0.718870995  0.233253787 -0.250027495
##  [946] -0.542905649  0.716920149  0.626694606  0.178276190  0.182912203
##  [951]  0.112550205  0.045196319  0.126569024  0.129463196 -0.322787666
##  [956]  0.130540035  1.388155028  0.458108951  0.051392435  0.194071759
##  [961]  0.196373735 -0.188147380  0.734269193  0.577301082  1.007652053
##  [966] -0.182590946  0.522293033 -0.664368134  0.579626699  0.822877529
##  [971]  0.408912260  0.458703944 -0.620828966  0.441162324 -0.251299033
##  [976]  1.568299555  0.687058833  0.478405543  0.630739218  1.428791654
##  [981]  0.409063316  0.022405759  0.238588142  0.367385444  0.160117555
##  [986]  0.323471750 -0.425178291  0.297963350  0.616455638 -0.582614267
##  [991]  0.609051316  0.522800931  0.802466368  0.106310582  0.220923065
##  [996]  0.424523566 -0.114709097  0.905672092  0.384093345 -0.238180645
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
```

```r
fit2
```

```
##       mean          sd    
##   0.69898817   0.29454813 
##  (0.09314430) (0.06585886)
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
## [1]  0.1695978  0.4392654 -0.6155409 -1.0684915 -0.1099665  0.1772334
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

Why do we not have to tell the 'jackknife' function how many replicates to do?

Let's compare this with what we would have obtained from bootstrapping


```r
results2<-bootstrap(x,1000,theta)
mean(results2$thetastar)-mean(x)  #this is the bias
```

```
## [1] 0.0228
```

```r
sd(results2$thetastar)  #the standard deviation of the theta stars is the SE of the statistic (in this case, the mean)
```

```
## [1] 0.8885661
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
##     original     bias    std. error
## t1*      4.5 0.03093093   0.8995848
```

One of the main advantages to the 'boot' package over the 'bootstrap' package is the nicer formatting of the output.

Going back to our original code, lets see how we could reproduce all of these numbers:


```r
table(sample(x,size=length(x),replace=T))
```

```
## 
## 4 5 6 8 9 
## 2 4 1 2 1
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
## [1] -0.0629
```

```r
se.boot
```

```
## [1] 0.9124835
```

Why do our numbers not agree exactly with those of the boot package? This is because our estimates of bias and standard error are just estimates, and they carry with them their own uncertainties. That is one of the reasons we might bother doing jackknife-after-bootstrap.

The 'boot' package has a LOT of functionality. If we have time, we will come back to some of these more complex functions later in the semester as we cover topics like regression and glm.

