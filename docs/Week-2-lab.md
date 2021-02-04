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
## 0 2 4 5 6 7 8 9 
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
## [1] 0.0176
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
## [1] 2.741997
```

```r
UL.boot
```

```
## [1] 6.293203
```

Method #2: Simply take the quantiles of the bootstrap statistics


```r
quantile(xmeans,c(0.025,0.975))
```

```
##   2.5%  97.5% 
## 2.7975 6.3025
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
##    [1] 4.2 4.6 2.8 4.8 5.2 3.7 6.1 5.5 3.8 4.8 4.6 5.5 4.1 4.5 6.2 4.6 5.1 4.8
##   [19] 3.4 4.3 4.4 5.4 3.9 4.7 4.8 4.9 3.5 5.0 3.7 3.1 4.9 6.8 5.8 5.3 2.6 5.5
##   [37] 2.7 6.0 4.5 4.2 6.2 4.6 3.0 4.1 3.9 4.0 5.0 4.7 3.3 4.4 3.4 3.7 4.0 5.4
##   [55] 6.2 3.7 5.3 4.7 5.4 4.4 4.8 4.3 4.2 4.1 4.0 3.7 4.7 6.3 4.5 4.7 4.3 3.9
##   [73] 4.8 3.3 4.3 4.1 6.5 4.6 4.1 4.6 6.6 3.3 3.7 5.4 3.5 4.8 4.9 6.0 7.1 5.5
##   [91] 2.8 3.6 3.8 2.8 4.7 3.8 6.2 5.1 3.1 6.2 5.1 5.1 4.5 4.7 5.4 4.9 3.5 7.2
##  [109] 3.2 3.2 5.0 4.3 4.6 3.5 5.3 4.6 4.9 3.1 5.2 4.8 4.6 4.1 6.1 5.4 4.5 5.2
##  [127] 4.7 3.3 4.1 5.1 3.7 3.3 4.9 4.9 4.7 4.9 4.5 4.3 4.3 3.1 3.8 2.6 5.2 5.4
##  [145] 4.5 4.2 3.9 4.9 5.2 4.6 3.9 5.5 2.3 4.1 4.4 4.8 5.3 4.7 4.4 5.3 5.9 3.2
##  [163] 5.1 3.9 4.0 4.6 4.4 3.7 3.7 4.8 5.1 3.8 4.1 5.1 3.9 4.3 3.2 3.5 6.3 5.2
##  [181] 4.9 3.0 4.4 5.3 5.8 4.8 5.6 3.6 3.9 3.3 4.5 4.1 2.6 3.5 5.2 3.6 4.8 4.8
##  [199] 4.4 3.2 4.8 5.7 4.4 4.4 4.5 4.0 5.1 4.8 4.1 5.0 3.3 3.3 5.9 4.0 4.4 4.7
##  [217] 4.6 3.4 4.8 5.2 4.7 4.3 5.0 5.0 3.0 2.9 4.2 3.7 3.4 3.6 4.8 6.1 3.3 4.4
##  [235] 4.5 6.2 5.1 4.0 5.6 3.1 4.6 3.9 4.9 3.4 5.6 4.4 3.5 3.9 5.1 4.3 4.1 5.6
##  [253] 5.1 4.6 4.8 5.2 3.7 4.7 2.8 3.2 5.2 6.5 6.3 2.5 5.3 3.5 4.8 5.0 6.6 5.3
##  [271] 5.3 5.5 4.0 4.9 4.3 3.8 4.5 5.1 4.4 3.3 3.8 4.3 5.1 3.6 4.6 3.7 4.8 6.7
##  [289] 4.9 5.2 4.1 5.4 6.9 3.9 5.8 3.1 3.4 2.3 4.1 3.7 3.7 5.2 4.5 3.9 6.2 5.3
##  [307] 5.0 5.7 3.9 6.1 5.8 3.1 3.4 5.2 5.5 5.6 5.3 3.2 4.0 4.9 4.3 2.9 3.5 5.3
##  [325] 3.7 4.7 3.8 3.6 4.3 4.5 5.2 4.5 4.4 6.4 4.9 5.2 6.2 5.1 4.1 3.3 5.1 5.6
##  [343] 4.7 4.6 3.1 3.6 3.1 3.5 4.8 4.3 5.7 3.0 5.0 3.6 5.3 3.7 5.3 5.8 5.5 3.3
##  [361] 3.8 5.7 3.8 5.8 5.8 3.5 4.2 4.2 3.0 3.3 5.5 4.3 3.3 2.6 2.0 3.9 5.0 3.4
##  [379] 3.3 4.5 4.8 2.6 4.4 4.8 5.3 4.1 5.2 4.9 4.8 4.8 3.5 3.2 5.3 5.6 4.3 4.8
##  [397] 4.2 5.3 3.9 4.5 4.7 6.1 4.0 3.9 4.0 5.5 2.4 4.8 6.2 4.7 4.4 5.1 4.1 5.1
##  [415] 4.4 5.3 5.4 5.6 4.2 3.0 2.8 4.2 4.4 3.6 5.0 6.1 5.4 5.2 5.6 5.2 3.8 3.5
##  [433] 4.3 4.3 4.6 4.9 4.5 3.0 5.5 3.7 4.8 3.2 4.8 5.0 5.7 3.6 5.6 3.4 4.0 5.9
##  [451] 4.0 5.3 3.7 3.2 4.3 4.5 3.6 5.1 4.5 3.3 4.7 4.3 5.8 5.4 3.0 3.8 4.2 2.9
##  [469] 3.9 4.8 3.3 5.1 4.4 5.0 4.0 4.1 3.8 4.8 4.0 7.1 3.6 4.8 2.4 5.4 4.7 3.3
##  [487] 4.6 4.0 4.2 3.1 3.8 3.7 5.9 3.8 5.4 5.5 3.6 4.6 4.8 5.6 5.3 4.8 5.8 4.1
##  [505] 6.2 5.3 5.0 4.8 6.2 4.9 3.8 4.2 5.4 5.0 5.8 3.9 2.7 5.0 4.8 3.2 3.3 4.3
##  [523] 4.5 4.6 5.0 4.5 4.3 4.6 1.9 4.7 4.7 5.6 4.7 4.3 3.8 4.1 5.5 4.6 5.6 3.3
##  [541] 4.2 6.2 4.4 4.0 4.6 4.4 4.0 4.6 5.6 2.4 4.2 4.5 4.0 3.7 4.6 4.7 3.7 4.6
##  [559] 4.4 3.9 4.0 5.4 3.6 3.0 5.7 3.6 3.8 3.1 5.1 2.8 5.2 4.6 6.3 4.0 5.7 7.1
##  [577] 6.5 4.5 4.9 4.8 4.5 4.4 4.5 3.8 3.6 4.5 5.4 4.2 5.5 5.1 4.1 5.6 5.3 5.1
##  [595] 4.2 5.4 4.3 5.5 3.9 3.8 5.7 5.0 3.5 3.9 4.0 3.5 3.9 5.6 4.5 4.7 3.8 5.3
##  [613] 4.5 5.5 3.3 4.0 3.1 4.5 5.0 3.8 5.0 5.0 4.7 3.9 4.4 4.0 4.3 4.9 4.0 4.3
##  [631] 4.1 4.0 5.1 3.1 4.5 3.7 3.2 3.7 4.1 3.9 4.3 3.2 3.5 3.7 4.6 3.0 4.2 5.5
##  [649] 4.2 4.6 4.8 4.5 2.6 4.6 6.0 5.4 4.8 5.6 4.1 6.3 5.5 6.1 5.3 3.5 5.2 4.2
##  [667] 5.6 4.5 3.0 2.5 4.8 4.1 5.6 3.8 3.8 5.6 4.1 4.4 4.2 5.1 3.4 5.0 6.1 3.3
##  [685] 3.1 4.5 4.2 3.8 3.1 5.5 2.7 7.1 3.8 5.8 4.0 5.6 4.6 3.8 5.0 5.6 5.3 2.4
##  [703] 2.7 5.1 3.5 4.3 5.4 4.9 4.6 6.1 4.9 5.6 4.4 4.2 3.9 3.4 4.3 5.3 6.1 4.5
##  [721] 5.4 5.4 2.8 4.2 1.7 4.1 3.1 4.8 4.7 4.9 4.8 4.1 1.9 3.8 4.6 5.5 3.3 3.5
##  [739] 5.1 6.7 5.4 3.6 5.1 4.7 3.8 5.2 4.2 4.9 4.5 4.2 5.6 5.6 3.4 4.2 5.0 3.4
##  [757] 4.8 4.6 4.8 5.9 5.5 6.1 3.8 3.4 5.6 5.7 3.6 5.6 3.7 4.9 3.4 4.7 4.7 4.1
##  [775] 3.9 5.2 4.1 4.6 3.9 5.9 4.9 5.0 3.4 4.5 3.8 5.2 2.9 4.4 6.2 2.6 3.6 4.9
##  [793] 3.1 6.6 6.5 4.2 5.5 4.7 3.0 4.8 3.5 4.3 6.1 5.5 3.6 3.7 3.4 2.7 4.6 6.0
##  [811] 4.1 3.4 3.7 4.8 3.3 3.8 6.0 6.4 3.0 5.6 5.2 5.7 5.0 5.0 5.7 3.7 5.6 4.7
##  [829] 3.8 5.8 4.7 5.4 4.7 5.8 5.0 5.6 4.4 3.9 3.4 3.1 4.7 4.9 3.6 6.9 2.8 5.4
##  [847] 2.4 4.5 4.2 5.5 5.0 4.9 5.6 4.6 6.3 4.0 4.0 4.6 4.9 4.3 6.5 3.6 5.7 4.6
##  [865] 5.0 5.8 3.3 3.7 5.7 5.0 3.5 3.5 3.7 5.3 3.8 5.4 5.8 4.0 4.4 3.8 4.4 4.9
##  [883] 2.6 6.1 5.8 4.4 4.1 3.9 4.5 3.5 3.2 4.8 4.2 4.2 4.2 4.2 4.8 4.8 4.3 4.1
##  [901] 4.5 3.9 3.3 4.1 6.5 5.0 4.2 4.5 3.4 3.7 5.5 4.2 5.4 5.0 3.0 3.6 6.1 2.5
##  [919] 3.4 3.9 5.9 3.5 3.3 5.5 3.3 3.1 4.5 5.0 3.1 5.6 3.6 4.5 4.0 2.8 7.2 4.7
##  [937] 5.4 2.0 4.9 4.5 5.3 5.1 5.3 4.2 5.6 3.4 3.6 4.1 4.5 1.4 7.0 5.1 4.4 4.4
##  [955] 3.9 3.6 4.9 4.5 6.2 5.3 3.9 5.4 3.3 4.7 3.5 3.8 5.5 6.1 4.1 3.7 4.5 4.9
##  [973] 3.4 5.2 5.8 4.2 4.4 4.7 3.4 5.5 6.1 3.8 5.4 3.1 3.9 4.8 4.3 4.0 4.6 6.1
##  [991] 5.3 4.8 3.9 3.4 5.2 3.2 4.8 2.9 5.7 4.7
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
##   2.7   6.3
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
##    [1] 3.9 4.0 3.0 4.5 3.9 2.6 3.2 5.3 5.0 5.3 4.0 4.4 3.5 3.9 3.7 3.1 4.2 4.4
##   [19] 4.5 3.2 4.3 4.2 3.2 5.8 5.6 4.0 4.9 3.8 4.2 5.2 3.9 3.6 6.8 4.2 5.0 4.0
##   [37] 3.9 4.1 3.4 5.2 4.0 3.3 5.1 3.8 5.8 4.6 4.6 2.2 5.3 5.4 3.2 6.7 6.5 5.7
##   [55] 3.8 5.0 6.1 3.1 4.8 3.0 4.2 4.9 4.2 4.8 4.2 4.8 5.4 4.8 4.4 5.5 5.6 4.3
##   [73] 5.5 3.5 5.6 5.3 4.4 5.7 5.1 5.8 4.5 2.4 3.9 4.4 5.8 3.1 4.9 3.9 4.8 4.7
##   [91] 2.5 2.9 3.9 4.0 5.1 3.8 4.0 4.5 4.8 5.0 4.2 4.4 4.7 4.4 3.2 5.4 4.4 4.0
##  [109] 4.7 5.6 4.4 3.8 3.1 5.4 3.6 6.5 4.6 3.6 5.5 5.8 3.0 3.4 4.5 6.1 4.4 5.5
##  [127] 3.1 4.8 5.3 4.6 3.4 3.5 5.8 4.4 4.8 6.8 5.4 3.6 5.0 5.6 4.1 4.6 4.5 4.0
##  [145] 4.6 3.8 3.9 4.1 5.0 4.5 4.9 4.4 2.9 3.3 4.2 4.4 3.9 3.3 4.7 4.9 3.3 5.3
##  [163] 5.0 6.2 5.3 4.4 4.8 4.3 4.0 5.8 5.7 3.7 6.1 3.6 4.9 4.1 4.7 4.1 5.0 4.0
##  [181] 5.6 3.9 4.6 4.5 4.4 4.3 4.0 4.8 4.0 4.6 4.2 4.5 4.9 3.8 5.8 5.7 4.9 4.2
##  [199] 6.2 3.5 4.7 4.1 4.4 3.2 4.2 4.4 2.5 5.6 3.7 4.9 4.3 4.5 3.4 3.6 6.4 5.6
##  [217] 4.0 4.5 4.0 4.0 3.7 4.3 4.2 4.7 5.8 4.2 4.2 4.3 5.4 4.6 3.5 4.2 4.0 5.2
##  [235] 4.9 3.8 4.2 4.2 5.2 4.2 4.7 5.3 4.9 3.5 5.1 4.9 4.8 3.8 4.6 3.4 4.8 5.1
##  [253] 4.8 3.7 4.4 4.0 5.3 4.7 6.2 4.0 4.1 4.4 3.8 4.6 5.7 3.1 3.8 2.9 4.1 2.8
##  [271] 4.8 3.9 4.1 4.0 3.9 3.0 4.7 5.8 3.6 3.9 4.0 4.4 3.8 6.2 4.2 5.3 3.1 3.7
##  [289] 5.6 4.8 5.3 6.5 4.1 6.0 3.6 4.0 4.3 3.8 4.1 4.4 2.8 3.4 2.9 4.4 4.6 3.9
##  [307] 4.0 4.6 4.9 5.4 4.5 3.0 3.8 3.5 3.0 3.4 6.4 4.2 5.5 5.5 5.4 6.6 4.9 5.6
##  [325] 4.2 5.6 5.5 5.4 3.4 4.5 6.6 5.3 4.8 3.6 7.3 4.3 4.4 3.9 4.2 4.6 6.1 4.4
##  [343] 3.3 4.3 4.3 5.2 3.7 4.6 3.7 4.8 3.7 4.3 5.7 5.4 3.5 4.8 3.7 4.9 4.0 4.8
##  [361] 3.9 4.6 2.8 3.0 5.0 4.7 5.0 4.6 4.5 4.3 4.8 4.5 4.7 4.2 4.9 3.9 5.2 5.1
##  [379] 2.8 4.8 5.4 5.0 4.9 5.7 3.0 4.9 4.6 4.7 3.8 5.5 3.5 4.8 4.5 4.1 2.9 4.4
##  [397] 3.9 4.0 4.8 4.3 5.0 4.2 4.5 4.3 4.8 4.8 3.0 4.6 4.8 4.6 3.4 3.9 4.3 2.2
##  [415] 4.7 2.9 3.4 4.5 3.2 5.6 5.6 5.1 3.1 5.5 4.9 3.4 5.1 5.4 5.0 4.6 5.6 4.1
##  [433] 4.8 3.5 3.2 5.8 4.1 4.5 5.5 4.1 6.0 5.3 6.2 3.6 3.5 4.3 4.2 4.6 4.3 3.7
##  [451] 4.9 4.3 6.2 5.3 3.7 3.5 6.7 6.1 4.6 3.6 3.0 4.1 6.4 4.2 3.3 4.4 3.0 4.9
##  [469] 4.5 3.6 4.0 5.8 4.3 5.1 3.5 5.1 3.0 3.8 3.6 3.9 5.2 3.4 5.3 4.7 5.3 6.1
##  [487] 3.5 5.0 4.9 6.3 3.8 5.6 4.2 4.7 4.6 3.8 4.4 4.7 3.7 3.4 6.2 3.6 3.9 4.3
##  [505] 3.8 6.0 6.3 3.5 3.5 4.8 4.7 3.8 3.9 5.8 5.6 3.5 6.6 2.6 4.5 1.7 4.3 3.9
##  [523] 2.7 3.8 5.3 4.1 3.9 3.2 4.6 5.1 3.6 3.8 3.8 5.6 6.1 3.6 3.4 5.6 4.8 3.1
##  [541] 3.9 4.7 3.7 4.6 4.5 4.3 4.2 4.5 4.7 4.8 4.1 2.3 5.5 3.9 4.5 5.8 5.7 4.5
##  [559] 5.0 6.1 5.2 4.3 6.2 6.1 2.9 5.1 4.2 5.2 4.2 5.1 5.7 3.9 4.4 3.3 4.6 6.4
##  [577] 3.3 3.5 5.0 3.5 3.0 4.0 3.0 5.8 4.6 4.1 3.5 5.6 4.7 3.2 6.4 4.7 3.6 4.4
##  [595] 2.6 6.3 4.0 5.7 4.3 3.8 5.0 4.1 3.0 4.2 5.4 2.6 5.9 5.7 4.8 4.6 3.5 5.4
##  [613] 4.7 5.7 4.6 4.9 2.8 5.3 4.5 4.1 4.8 3.6 4.3 5.6 4.6 4.9 5.4 4.4 3.0 5.2
##  [631] 6.0 5.2 4.0 2.9 4.7 4.7 1.9 4.6 5.8 1.8 5.3 7.0 4.6 2.5 3.9 5.8 3.7 4.0
##  [649] 4.4 5.4 5.2 2.4 3.7 3.5 3.7 5.0 3.9 4.8 3.8 5.1 4.7 6.2 4.8 1.7 5.3 4.7
##  [667] 5.2 4.5 6.0 3.6 4.5 4.3 4.1 4.1 4.8 5.3 4.8 3.8 3.5 4.8 5.0 4.1 4.3 2.7
##  [685] 3.8 4.6 5.2 3.7 3.4 6.0 4.2 5.1 3.9 6.0 4.5 4.5 4.7 4.2 4.1 3.4 5.8 5.1
##  [703] 4.1 4.7 3.3 4.3 6.1 3.6 5.7 4.4 5.2 4.9 4.6 5.5 4.2 3.4 4.7 5.0 3.2 4.3
##  [721] 2.7 2.0 5.3 5.5 4.0 5.9 6.8 4.4 4.6 3.1 5.4 3.9 6.2 4.8 3.9 4.1 3.8 3.4
##  [739] 4.3 6.0 2.7 3.9 5.2 5.6 3.3 3.9 4.5 5.7 4.0 4.6 3.3 4.4 4.2 4.5 4.9 5.0
##  [757] 1.9 3.6 4.2 4.7 4.3 3.5 4.9 4.7 4.1 3.2 4.2 3.3 5.5 6.6 2.6 2.4 5.5 4.8
##  [775] 5.3 3.6 4.3 4.7 5.0 3.5 3.2 4.3 4.5 2.7 4.1 5.2 4.1 5.3 4.7 5.2 6.2 4.7
##  [793] 3.4 3.7 4.1 5.1 4.8 3.0 6.0 4.6 4.9 5.3 4.8 4.5 4.8 4.1 4.9 3.5 4.3 5.7
##  [811] 4.4 3.7 3.8 5.1 4.3 4.8 5.3 4.6 5.7 5.8 5.5 5.4 5.2 4.3 4.7 5.1 3.5 5.5
##  [829] 4.4 5.8 4.1 4.3 4.8 2.0 4.9 3.9 4.2 4.2 5.5 6.4 5.6 4.5 3.8 4.7 5.2 5.8
##  [847] 5.2 4.2 4.9 4.2 5.4 3.4 3.2 3.2 4.8 4.9 5.2 4.5 5.2 5.3 6.0 2.7 4.6 4.3
##  [865] 3.5 4.0 4.2 3.7 4.3 4.0 5.3 5.9 6.7 5.6 3.9 5.4 5.8 4.4 3.7 4.7 3.8 4.6
##  [883] 4.6 3.6 4.1 5.4 3.1 5.4 6.0 3.7 3.3 3.5 4.4 5.1 3.6 3.6 5.6 3.6 4.4 4.5
##  [901] 4.3 4.1 3.3 3.5 4.1 6.3 3.6 5.1 4.7 3.6 4.2 3.3 4.5 3.7 5.3 5.7 4.1 3.8
##  [919] 5.3 4.1 4.4 4.4 3.6 5.4 3.8 4.0 6.0 5.8 5.8 3.5 4.2 4.9 3.5 4.9 2.1 4.8
##  [937] 3.3 5.6 4.2 3.5 5.0 4.9 6.6 4.3 5.4 3.8 5.8 4.3 3.8 6.5 2.1 6.2 4.1 3.6
##  [955] 4.1 6.2 3.6 2.7 4.5 3.6 4.1 3.1 3.9 4.9 5.8 3.6 5.8 3.8 3.3 3.8 4.6 3.3
##  [973] 3.8 3.5 6.2 3.9 4.6 5.0 3.1 4.0 4.4 5.3 4.5 4.0 5.7 4.2 4.2 3.3 4.8 4.4
##  [991] 4.2 4.6 4.1 6.0 5.1 4.8 3.8 5.4 4.6 3.7
## 
## $func.thetastar
## [1] -0.0397
## 
## $jack.boot.val
##  [1]  0.4585507  0.3815341  0.2731928  0.1436747  0.0750000 -0.0960961
##  [7] -0.2038251 -0.3513120 -0.4571031 -0.5626781
## 
## $jack.boot.se
## [1] 1.015132
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
##    [1] 5.8 4.6 4.2 3.7 3.8 5.4 4.2 6.3 4.8 3.0 3.7 5.7 4.1 2.8 2.3 4.6 4.4 5.4
##   [19] 3.4 4.7 4.5 4.5 5.8 3.5 3.7 2.7 3.4 4.4 5.1 5.1 6.3 6.4 3.8 3.6 5.1 3.9
##   [37] 5.2 3.8 3.6 3.6 2.1 3.2 4.2 5.0 4.9 7.3 4.0 3.8 5.7 3.6 5.5 3.8 5.2 4.1
##   [55] 3.8 5.6 5.6 3.6 5.2 4.7 3.9 4.5 3.3 4.8 4.8 4.5 3.6 3.7 3.3 4.7 5.4 5.3
##   [73] 3.6 5.5 5.6 4.1 4.4 3.7 5.8 3.1 3.8 4.3 5.1 5.6 4.7 5.2 4.4 3.8 5.5 5.0
##   [91] 4.6 4.2 5.2 3.6 4.0 3.1 4.2 4.2 4.7 6.5 3.3 3.5 3.3 3.2 5.1 4.2 4.8 4.8
##  [109] 4.3 5.2 4.7 5.2 4.4 5.3 4.6 5.7 4.3 5.3 3.2 4.8 6.5 5.5 3.1 4.3 4.4 4.8
##  [127] 5.6 3.8 4.0 4.4 4.0 4.3 3.6 4.8 3.8 5.3 4.0 4.2 4.1 4.9 4.1 6.0 3.2 3.7
##  [145] 5.5 6.0 5.6 4.9 4.3 4.1 2.4 6.2 3.8 2.9 5.4 3.8 4.1 5.8 3.6 4.7 5.4 5.2
##  [163] 5.8 4.1 3.6 5.5 6.5 3.9 6.2 3.6 3.6 4.4 3.4 4.2 4.7 5.5 6.1 3.5 4.1 4.9
##  [181] 5.2 3.9 4.3 5.2 6.6 4.1 4.6 3.4 4.8 4.3 4.6 4.2 5.1 4.3 5.2 3.5 5.0 5.3
##  [199] 3.8 4.2 4.5 3.2 4.2 6.5 5.9 3.7 3.1 3.1 4.3 6.5 3.9 5.7 4.9 3.6 4.1 3.8
##  [217] 6.2 4.0 6.1 3.5 4.2 4.3 6.9 5.5 5.3 4.4 4.2 4.7 5.6 6.0 3.9 5.3 5.0 3.9
##  [235] 6.5 5.6 5.9 4.7 4.6 3.1 3.2 2.7 5.1 4.9 4.1 5.4 4.5 4.1 4.1 4.9 4.4 4.3
##  [253] 5.2 4.6 3.4 4.2 4.8 4.8 4.3 4.5 2.3 4.2 3.4 5.4 4.9 4.1 3.9 4.6 6.4 5.1
##  [271] 5.0 4.5 7.0 4.4 6.5 4.2 5.2 2.9 5.1 4.4 6.1 3.7 4.8 4.7 5.0 3.2 4.0 3.7
##  [289] 6.7 5.8 3.8 3.8 4.2 2.9 4.1 5.8 6.0 4.2 4.5 5.9 3.6 3.6 5.2 4.5 4.1 5.9
##  [307] 5.2 5.1 4.0 4.8 4.8 4.4 3.2 5.6 3.7 5.2 4.6 4.3 3.6 4.1 2.0 4.2 4.1 4.8
##  [325] 4.6 4.1 4.4 4.4 4.4 4.9 2.8 5.1 5.6 3.4 4.0 5.8 3.6 2.9 3.7 6.4 3.5 4.4
##  [343] 5.5 6.3 4.8 3.6 3.3 2.1 4.0 4.1 5.6 4.8 2.6 4.0 5.1 3.0 6.0 3.5 3.3 5.0
##  [361] 5.3 4.6 4.3 3.4 4.8 3.6 3.4 4.6 5.2 6.0 3.7 5.2 5.1 3.3 4.6 4.0 4.5 3.5
##  [379] 6.5 4.1 4.3 5.6 2.7 4.3 4.7 4.3 5.1 3.9 5.0 4.5 4.5 3.4 3.8 3.3 5.3 4.5
##  [397] 5.3 4.7 4.4 3.8 5.0 4.6 2.7 3.3 3.7 4.2 4.2 3.6 4.2 6.2 4.1 3.2 5.7 3.1
##  [415] 4.5 4.4 7.0 4.7 3.8 4.2 3.1 4.6 5.0 5.2 4.3 6.1 3.5 4.3 5.5 3.5 5.1 2.8
##  [433] 4.6 4.6 5.6 4.9 3.3 4.6 4.9 3.7 3.6 2.7 5.3 5.9 4.1 4.9 4.1 4.5 4.6 4.4
##  [451] 3.6 6.3 4.8 4.4 4.0 4.8 7.2 3.9 3.8 3.9 3.3 5.4 4.9 6.0 2.3 5.6 3.8 5.3
##  [469] 4.5 5.9 3.9 4.9 6.5 4.6 5.0 5.9 3.0 3.8 4.0 4.0 4.2 3.8 4.5 4.4 6.5 3.0
##  [487] 4.4 3.8 4.4 3.5 4.7 4.3 4.3 5.4 4.7 3.6 4.9 3.7 5.2 5.0 3.0 4.9 5.9 4.4
##  [505] 3.0 4.8 3.7 5.2 3.2 4.0 4.2 4.4 3.6 4.2 4.1 4.3 3.9 5.4 4.2 4.6 4.9 4.5
##  [523] 3.9 3.8 4.5 4.5 4.2 5.1 5.4 5.3 5.3 4.7 3.5 3.9 3.5 4.2 5.2 5.0 4.7 3.8
##  [541] 4.7 5.8 4.9 4.2 5.2 5.7 3.9 5.6 4.1 2.5 5.6 5.3 4.4 5.7 3.4 4.0 4.8 5.4
##  [559] 5.6 6.7 2.8 5.1 4.2 3.3 4.1 5.8 1.5 6.2 5.8 6.1 5.0 4.4 4.4 5.8 4.2 4.4
##  [577] 6.1 3.5 6.2 4.9 4.1 4.2 5.6 5.0 3.5 4.4 3.6 4.1 3.5 4.1 6.3 4.8 3.2 5.0
##  [595] 5.1 4.1 4.9 4.8 2.8 4.1 3.6 4.8 4.1 4.5 4.2 3.9 4.7 5.5 5.3 5.2 3.9 4.7
##  [613] 4.0 5.3 4.5 5.3 5.1 3.7 5.6 3.7 5.8 2.3 3.7 4.0 3.6 5.7 3.8 5.2 4.4 5.8
##  [631] 3.4 6.6 3.3 3.9 4.4 4.6 3.8 5.7 3.0 7.0 4.1 3.9 5.3 4.9 4.6 3.2 3.9 3.5
##  [649] 5.0 4.3 3.4 5.1 4.9 3.0 6.1 5.1 3.6 5.0 4.8 4.5 3.7 3.7 4.0 4.5 4.7 5.2
##  [667] 3.9 5.5 3.6 4.0 3.7 3.0 5.1 4.5 4.0 4.6 3.5 3.6 4.3 5.1 6.3 3.6 5.9 3.6
##  [685] 5.2 5.5 4.8 3.2 3.6 2.7 5.7 3.6 5.2 3.3 5.2 5.4 4.6 3.3 3.4 4.7 6.3 3.6
##  [703] 4.5 3.4 4.6 4.2 3.8 2.3 5.0 4.0 4.5 4.7 3.1 3.4 3.9 5.6 3.8 4.3 5.2 4.9
##  [721] 4.7 4.1 3.1 3.1 3.5 4.1 5.4 3.3 4.6 4.2 6.0 4.1 4.8 2.3 3.2 5.8 3.4 6.4
##  [739] 5.3 3.5 5.1 4.4 3.9 6.8 6.2 5.7 5.4 4.5 4.4 3.6 4.9 6.1 5.0 5.8 5.0 3.9
##  [757] 4.3 3.5 5.6 2.5 3.9 3.2 3.7 4.4 3.5 5.3 4.7 5.5 3.6 4.4 3.8 4.2 3.9 2.0
##  [775] 4.7 3.8 3.9 4.8 6.8 4.6 4.0 3.8 4.2 4.6 4.3 2.5 5.5 6.8 4.6 5.2 4.1 4.7
##  [793] 5.4 4.9 6.0 4.5 4.6 5.0 3.2 4.5 4.8 3.2 4.1 4.7 4.5 5.3 5.2 6.2 3.1 5.2
##  [811] 2.6 6.0 3.7 3.5 4.7 4.3 4.7 4.4 3.3 4.1 3.4 5.3 4.8 4.7 4.8 5.7 4.5 4.3
##  [829] 4.9 4.8 5.1 3.6 5.6 5.1 4.5 4.4 5.1 4.4 4.2 4.5 4.2 5.3 4.5 4.2 3.9 3.6
##  [847] 3.2 4.6 4.2 2.8 5.9 4.3 3.9 2.9 5.1 4.8 4.9 4.0 6.1 3.5 5.1 3.4 3.1 5.6
##  [865] 5.8 6.2 4.5 3.0 4.8 5.3 4.4 3.7 5.2 4.3 3.1 5.0 3.6 3.3 3.2 3.7 4.3 4.0
##  [883] 6.2 3.9 4.7 5.2 5.6 5.2 3.9 4.0 5.6 5.3 4.3 3.5 5.4 4.2 4.6 6.8 3.4 4.6
##  [901] 4.8 5.5 4.3 5.9 2.7 4.8 4.4 4.3 5.8 3.5 4.0 4.0 5.3 5.8 4.9 4.3 4.4 4.8
##  [919] 4.1 2.9 3.7 5.7 6.0 5.1 4.3 4.1 4.6 5.3 3.9 6.0 6.1 4.6 3.7 4.1 3.8 3.7
##  [937] 6.5 5.1 5.2 5.8 5.4 5.8 4.4 5.0 5.0 5.2 5.7 4.6 5.7 6.0 4.6 3.0 4.4 4.6
##  [955] 4.0 5.8 4.5 5.5 4.8 3.3 4.6 4.7 4.5 3.4 5.4 7.3 4.5 3.9 4.2 3.5 4.0 2.9
##  [973] 3.3 4.9 2.8 5.2 5.7 4.1 3.8 3.2 6.8 3.1 5.5 5.1 5.6 3.1 3.6 6.5 3.4 4.6
##  [991] 4.8 3.6 2.9 4.5 3.2 4.4 5.6 4.7 3.5 4.3
## 
## $func.thetastar
## 72% 
## 5.1 
## 
## $jack.boot.val
##  [1] 5.6 5.5 5.4 5.2 5.2 5.1 4.9 4.7 4.6 4.4
## 
## $jack.boot.se
## [1] 1.14
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
## [1] 0.4924596
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
##   2.960997   5.206098 
##  (1.256766) (2.407907)
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
## [1] 0.6395853 0.9300531 0.9918624 0.9433125 0.3380203 0.9737417
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
##    [1]  0.1311027604  2.4019657990  0.8924412701  0.8065237708  0.4506056979
##    [6]  0.6736379973  0.8386546561  0.5503895795  2.2639747003 -0.4575596753
##   [11]  0.9270510992  1.3323152189 -0.6345368949  2.0228430166  0.8549495907
##   [16]  2.1211196219  1.2721963276 -0.9514243186  0.9722947340 -0.1641160129
##   [21]  1.5404819193 -0.2363861246  0.8606763542 -1.4263513689  0.1188902189
##   [26] -0.1982403875  0.5637625371 -0.0235611789  0.4587410930  1.5146351157
##   [31]  0.1474460128  0.1620602316  1.1971454471  1.2209167889  0.4962548982
##   [36]  0.5321188893  0.3161391066  0.6736379973  1.5158316693  0.4931340355
##   [41]  0.9218508167  0.7076678370 -0.0141371779  2.4393664643 -0.1787045785
##   [46]  1.1492869720  0.5622496142  1.2907864106  0.4926706178  0.3132066038
##   [51]  1.2912605834 -0.2292695521  0.5449944654 -0.2394571408  1.2346427017
##   [56] -0.3013361185 -0.0231921629  2.5955061708  0.4754205235  0.5426922594
##   [61]  1.4112363262  0.1114715324  1.2664945639  0.8420417258  0.5119183873
##   [66]  0.1685054347  0.1191942806  2.1238408917  0.0983181219  0.4806579879
##   [71]  1.5535123304  0.4991073812 -0.4009281148  0.9539903223  0.4722412765
##   [76]  1.2817229853 -0.2912929872  0.5684835988  0.4883573210  0.4085537858
##   [81]  0.5622553908  0.1380390275  0.1615757192  0.4943680307  0.6946967391
##   [86]  0.4882589031  0.4008975459  0.0428271667 -0.1904801827 -0.5359464184
##   [91] -0.2145264507  0.1955433741  0.8225948060  1.8350751728  0.7387827958
##   [96]  2.0749615192  1.2715617777 -0.1291025924  0.4345745582  0.7602368271
##  [101]  0.7737257725  0.9098974815 -0.1921260659 -0.8936001713  1.7226023878
##  [106] -0.5496354198  0.5991676577  0.5451555269  1.0511824101 -0.0795532698
##  [111]  0.1316272246  0.7186897559  1.1641786422  0.8456273662  0.5363855581
##  [116]  0.2076510737  1.4280520598  1.9315586435  0.8713280774  0.3292870232
##  [121]  0.9540354194  0.2450750016 -0.0871363497  1.5856855302  0.4684793979
##  [126] -0.6030351201  0.1199297226  1.5137762807  1.1562849869  0.5033200433
##  [131] -0.0141347778  0.5370639148  0.5362593013  0.1929235688 -0.1669081578
##  [136]  1.4146592101 -0.4370715361 -0.0554813414  0.3348648191  1.4586975043
##  [141]  0.5524115765  0.1268744843  1.0010351947  1.5049007368 -0.1926297335
##  [146] -0.5841580631  0.4898617576  0.8875991641 -0.2946197438 -0.1890213556
##  [151]  0.8965191975  0.1629024891  0.4605413344  0.1036201332  1.9212859826
##  [156]  1.6089732178 -0.1575545200 -0.6389697033  2.1316753067 -0.4180412126
##  [161]  0.8887832909 -0.2546868072 -0.2745147195  0.9416113059  0.4819250321
##  [166]  0.3868979717 -1.0449585096  0.8454042519  0.4453885480  0.8929168044
##  [171]  2.1015222260  0.9102116042  1.7881318047 -0.3074695799  1.1852076769
##  [176]  0.9292048547  0.9913434255  0.8409446060  0.8066400069  1.6074284784
##  [181]  0.5468178225 -0.2315916053  0.4718995578  0.1220667961  1.1821970810
##  [186]  0.9068573259  0.5413958810  0.4473834003  1.5149455743  0.5126363085
##  [191]  0.3887024272  1.4017171714 -0.1956565976  0.5215522097  0.5634401262
##  [196] -0.1672915159  0.8111304071  0.5058310905 -0.6131095581  0.1775192935
##  [201]  0.1645575888 -0.6175838910  0.5086544797  0.4579412246  1.1147724862
##  [206]  0.5745832415  0.1999931760  1.5466272607  2.2135966365  0.0941303246
##  [211] -0.9359936511  0.5594292792  0.0807597337  0.1381777322  1.9224240858
##  [216] -0.1888006037  0.1803558748  0.8878998064  0.3903694442  0.3031640379
##  [221]  0.5373405625 -0.0057830748 -0.0304988393  0.1384693974  0.1447109900
##  [226]  0.4711363247  0.7737257725 -0.2713535081  1.2949773657 -0.0922270877
##  [231] -0.0672477442 -0.1657013977  0.7625764519  0.1425765855  0.4450851214
##  [236]  0.8480914193  0.1226884160 -0.2090934424  0.4748900572 -0.0653907811
##  [241]  0.8695528664  0.4530108239  1.4174500837  0.8437905504 -0.2311951523
##  [246]  0.4396939789  0.5257947232  0.4296607407 -0.6895024607 -0.2797791574
##  [251]  0.5331415997  0.8887832909  0.9327074622  0.9080834083  0.3835665081
##  [256]  0.2345808251 -0.4008087840 -0.5349110792  1.0877766693  0.4978084045
##  [261] -0.0201027763  0.2544523558  0.1671104955  1.5863846789  0.7348233128
##  [266]  0.4727003645  0.1614989837 -0.1212029325 -0.0972807012  0.9268771672
##  [271]  0.7638806592  2.3640138905  0.1103573360 -0.5263143517  0.7677970557
##  [276]  0.1040178634  0.9966072877  0.1582694026 -0.1830767434  0.4924596270
##  [281]  2.3608312817 -0.5336407890  0.7245105583  2.2563208327  0.4609822802
##  [286]  0.5054852831  0.5775833656  0.1490955448  0.2416153470 -0.0428815631
##  [291]  1.9297778200  0.1185993257  1.4646465276  0.2918803589  0.5067276523
##  [296]  1.1774152909  0.0897712427 -0.5672319833  0.1136327714  0.8565247369
##  [301] -0.2063840409 -0.2549887492  1.1668196270  1.0016028627  0.1988546564
##  [306]  0.5332846750  0.2273843748  0.3222750778  0.4641253623  0.3154642557
##  [311] -0.0943127313  1.0143680824  0.5237097278 -0.5112229412  0.9513030293
##  [316]  0.1119265461 -0.0394565789  1.0143680824  0.9468240617  1.5734829738
##  [321]  0.8680682382  1.4141209096  0.0521115749  0.1238562405  0.8234534080
##  [326]  0.4580717603  1.2462372499  0.1748867654 -0.2841369682  0.9335656684
##  [331]  0.5190887602  1.3037420212  0.0784732417  2.1020850047  0.9615451116
##  [336]  0.4921194210  0.4429870982  0.1880198694  0.5443043841  0.8551351680
##  [341]  1.4923898767  0.8847430475  0.6860705247 -0.4117826452  0.5416423558
##  [346] -0.5687812894  0.8465651451  0.3130987897  0.3156797344  0.8357561094
##  [351] -0.0527613516  1.0415202967  0.7459911603  0.1308092883 -0.1846805861
##  [356]  0.3179505784 -0.4734195633  0.3256654946  0.7433615756  0.5325248797
##  [361]  0.9602901144  0.4093212496  0.3203334551  1.4106957027 -0.4099270165
##  [366]  0.2942378739  0.5380210951  1.3992616187  0.0837824467  0.8234762724
##  [371]  0.9668985048  0.5294574298  0.9421217836  1.4195837378  0.4676468608
##  [376]  0.8305842107  0.1785064962  0.8548337312 -0.8909178500  0.4370436881
##  [381]  0.9208874358  0.1638211415  0.8277876350  0.4258574978  0.3845309820
##  [386]  0.5065352155  1.4011715635  1.4466866321 -0.4354715005  0.1474344319
##  [391]  0.7475319382  1.5770480415  0.9337797236 -0.6136139833  0.0884173599
##  [396]  0.4571156972  0.0879298341  0.9367046287  0.4672694316 -0.0206877590
##  [401]  0.3112403102  2.4846319215  0.3622142456 -0.1398474685  1.2992791206
##  [406]  0.5510748816  0.0250420463 -0.0161202057  2.3450157767  0.1581347410
##  [411]  2.5526633610  1.4071432412 -0.1832087834  2.4213668897  1.3767762097
##  [416]  1.6242529363  0.4717792262  0.1729269937  1.2726211672  0.1999853212
##  [421]  1.4027558213  0.1421009726 -0.4731837592  0.6416152229  0.4292282664
##  [426]  1.2176302938  1.4222781318  0.5054891680  0.3329089692  1.1912029539
##  [431] -0.0938368424  0.1216704455  0.5176250331  0.0891934926  0.5818262270
##  [436]  0.2972489866 -0.6162206427  0.2458108118  0.4950066269  0.7311513785
##  [441]  1.1501697239  0.5520456580  0.7630446599  0.9226110035  1.5975117541
##  [446]  0.5116564985 -0.2203226229  0.9069774152  0.8207448788  0.1365565569
##  [451]  0.0925725541  0.6055388678  1.4796622477  0.9104818790  0.8843568493
##  [456]  0.3264400641 -0.2206275165  0.5560841216 -0.6684557913  0.0870993437
##  [461]  1.0304627514 -0.4994483490  0.8482964499  0.6837052826 -0.2910711494
##  [466]  0.4801652957  1.1393616237  0.1367718605  0.2385503212  0.7045830771
##  [471]  0.1833943677  1.0436542554  1.5729180145  1.1839450670  0.3209667841
##  [476]  0.1740095741  1.0110734582  1.4169104578  0.4784464065  0.5083584046
##  [481]  0.7403034331 -0.0689468618  0.4688059677  0.1000306016  0.9532678239
##  [486]  1.0016611019  0.9604111249  0.1328878129  0.7459786534  2.1435820179
##  [491]  0.9421987156  0.5042616781 -0.0624605262  0.4372613003  0.1001209789
##  [496]  0.3947085853  0.9735501934  0.3149596564  0.9296803905  0.9212346095
##  [501] -0.6620747687  1.1720063299  0.3308058459  0.7974865854 -1.4847539593
##  [506]  0.7464342902  0.4870483306  0.6475654316  2.2662391261  0.4363618364
##  [511]  0.4806998312  0.9752906354  1.3156442875  0.1327534029  0.1947718019
##  [516]  0.7433853822  0.1077427243 -1.4931680760 -0.5962347991  0.3204639293
##  [521]  1.4235493393  1.4139370942  0.6952394404 -0.4828884190  0.8880548454
##  [526]  0.1902217298 -0.4328735226  0.4370066752 -0.6477506747  0.8687037021
##  [531] -0.1582429306  0.9830088762 -0.2559501502  0.4732965010  0.4462597969
##  [536]  0.5312049022  1.3992615892 -0.0179690055  0.8338221381  0.2505467927
##  [541]  0.2968923925  0.5151990844  0.9638405399 -1.0269963235  0.5523703738
##  [546]  0.1649205558  0.1847968609  0.9065865440  0.8617182365  0.9898870059
##  [551]  0.0071111079  1.2987454941  0.3022357575  0.3880931406  0.1455603524
##  [556]  0.3434365864  0.5157245634  0.8168891853  0.5033524567 -0.2210221794
##  [561]  0.3359819551  0.7661848462  0.4599251290 -0.2314794418  0.8189629527
##  [566]  0.8925653662  0.1328878129  0.9222632318  0.1489576822 -0.2674936678
##  [571]  0.5034326456  1.5645213554  0.6848370353  0.4085330452  0.3892420480
##  [576]  0.5301271919  0.4610835597  1.0973202964  0.2031856424  0.1275868448
##  [581]  0.7346325271  0.8322359560  1.4541267660  0.1020729941 -0.7888876004
##  [586]  0.3105476217  1.4477941211  0.1975895241 -0.5142949016  1.4625440134
##  [591]  0.9344594350  0.5177715515  0.9295825878  2.4331988676  0.5211406568
##  [596]  0.5442854297  0.2304111893  1.5585702396  0.5081203151  0.5552625330
##  [601] -0.1095544951  0.0486572890  0.1129818661  1.6742069278 -0.0021569989
##  [606]  0.1170885952  0.9376075005  0.8222327550 -0.0625106414 -0.5118651224
##  [611]  0.1139196618  0.8942682830  0.1847968609  0.4858822588  0.4167438288
##  [616]  0.8444023640  0.9917465712  0.1361434266  0.4925770313  0.5201454769
##  [621]  0.5102739240  0.8314370867 -0.6603116773  0.0705130900  1.0046638639
##  [626]  0.5188908912  0.9643045850  0.3073708785  0.8303426661  1.1421716766
##  [631]  0.8906137725 -0.2205115435  0.4774052217  0.8313301444  0.5022682723
##  [636]  0.1752993780  0.4502925596  0.5758423481 -0.8432274014  0.3051090441
##  [641]  0.0841118733  1.0465537871  0.3757985623  0.5296501650  1.2851969361
##  [646]  0.4121731284  2.4827155769 -0.3583878972  0.9820768875 -0.2515707746
##  [651]  0.2255262579 -0.0070808470 -0.1794230192  0.3225673788  2.1747556211
##  [656]  1.4594866829 -0.2535118523 -0.0474965463 -0.3275846319  0.4187630567
##  [661]  0.2219598660  0.0792193458 -0.0891050164  0.4611112104 -0.0162643254
##  [666]  0.5520577601 -0.0767773988 -0.0781849353  0.5161800598  0.4740318167
##  [671]  0.1446581836  0.2729873865  0.9879323869  0.4839994920  0.9942328044
##  [676]  0.9585204724  0.1181185261  0.3360844158  0.4009372008  1.1669761259
##  [681]  0.8313301444  0.8301219072  0.8762017172 -0.3230274822  1.3928688274
##  [686]  0.1785957095 -0.0721605704  0.3081285578  1.5065297606  0.1545683192
##  [691] -0.6438690320  0.9109593404  0.4756511697  0.4737355848 -0.8890164768
##  [696]  0.7704892907  0.5221539622  0.5078663418 -0.5677575332  1.0174116418
##  [701]  0.3268363811  2.3017770850 -0.2472736756  1.4687354139  0.1339049574
##  [706] -0.0879994752  0.1115972096  0.9339848428 -0.1660378256  0.1070521452
##  [711]  0.9679803838 -0.1375984602  0.6714122067  0.3275684339  0.5161141240
##  [716]  0.8091213234  1.4513003484  0.1471400920 -0.0390064758  0.7096670546
##  [721]  1.3931222811  0.4577259859  0.9345977314  2.3425173366  0.3102024183
##  [726]  1.6120089846  0.8889975287  1.3239643153  0.7284038116  0.2946073815
##  [731]  1.8771342158  0.1401841687  0.1643476581  1.1572864698  0.1227307454
##  [736]  0.4535181643  0.3069776398  0.4346548213 -1.0401228816  2.5328558264
##  [741]  0.2968900230  0.9691979855  0.3849034335  0.2120873135  1.0292739548
##  [746]  0.1674007390 -0.1580456927  0.8773300371  1.1043449776 -0.5005995535
##  [751]  0.3726219603 -0.0824962070  1.2322844355  0.1464669636  2.3888536672
##  [756]  1.3310307752  0.7345662396  0.7989006134  0.3002117115  0.1321690549
##  [761] -0.0006265342  0.4975440620 -0.5956423190  0.3028084805  0.1594144593
##  [766]  0.4274278057  0.5053772669  0.4827562496  0.4586772850  0.1639173980
##  [771]  2.4846319215  0.4922729830  0.5393517138  1.2424194915 -0.3908547692
##  [776]  0.4444734178  0.8971291181  0.8753122315  1.0408379651 -0.8761037320
##  [781]  1.1755710475  0.9267620126  0.5676480262  0.9197577376  0.9620867164
##  [786]  0.8247007261 -0.0455089032 -0.5590363268  0.4845545283  0.8805088025
##  [791]  1.4715641015  0.9449992560  0.1498746164  0.4955786356  1.2916486307
##  [796]  0.4730962313  1.2973833709  0.8364744953  0.4243705771  0.9389503424
##  [801]  0.1261135276  0.6736382094  0.3322112808  0.1490804255 -0.9183108200
##  [806] -0.2349945114  0.7518176097 -0.6486933443 -0.5372590175  0.9054202061
##  [811] -0.0873742075  1.1439355565  0.9421217836  0.5981237695 -0.5285619581
##  [816]  0.1069435038  0.8232995696  0.5266758384  0.9520223978  0.3852426631
##  [821]  1.4330334416  0.1178952555 -0.4685564242  0.6566257458  1.4401853864
##  [826] -0.5580845433  2.0854749591  0.8539963451  1.0060986809  0.1242100905
##  [831]  1.4837917263  0.7000412963  0.2420991674  0.8616976704  0.3246328191
##  [836]  0.1971460229  0.0055782370  0.9243289343  0.1122703546  0.5642533904
##  [841] -0.7023497206  0.4915644233  0.1706140235  0.3216601204  2.4358042447
##  [846]  0.1578869404  0.5296501650  0.5086544797  1.2888923404  0.8774734567
##  [851] -0.2518423021 -0.2036795404 -0.4917348834  0.1616766688  0.1831605889
##  [856]  1.3885768112  1.0435425703 -0.9964313135  0.1494511129  2.2307820157
##  [861] -0.2118015504  1.2215068893 -1.4360103450 -0.5467979728  0.5493260862
##  [866]  0.3796014670  0.4396939789  0.8465114884  0.5105240092  0.4683879054
##  [871]  0.2286619285  1.2749667352  0.6887322134 -0.4237028978  0.0496465122
##  [876]  0.7486637545  1.9595589111  1.2657342316  0.9191110146 -0.3190812178
##  [881]  0.4292754636  1.2968193571  0.8900500617  0.9448979411  0.9730302295
##  [886]  1.4712322468  1.4657223738  0.5348448723  1.2390071504  1.4398812165
##  [891]  0.5434148317  0.4241574094  1.3858666335  0.5014186785  0.0097433312
##  [896]  0.2877699606 -0.6575396420  1.5486134560  0.8135393787  0.1706992976
##  [901]  2.4624606570  0.2560464413  2.2506206608  0.5349125921  0.1054462579
##  [906] -0.1944438939 -0.6555986103  1.2000950706  1.4059209996 -0.2387431065
##  [911]  0.7111007410  0.6659647076  0.9636306741  1.3953392407  0.1559535564
##  [916] -0.0860420918  0.0120605027  0.1397902509  0.0621338164 -0.5156322704
##  [921] -0.5472716496  0.6680904368  0.6515386408 -0.2222259872 -0.2757696123
##  [926] -0.4178730414  0.3835019612  2.4753237481  0.3880372885  0.2905944349
##  [931]  0.5604434181  0.5033501042  0.8329436498  0.8467720643  1.6245856632
##  [936]  0.1971390683  0.1989662241  0.1596176884  0.6652551627  0.8225054363
##  [941]  0.4555949631  0.4802333292  0.5081203151  2.0719848432  0.8339579785
##  [946]  0.2508442543  0.1531455901 -0.2512661785  0.9495169033  0.5515057635
##  [951]  0.6867580636  1.3185147148  0.3031640379  0.4668938793 -0.1744277199
##  [956]  0.1508058278  0.2993323914  0.7939105346  1.1987377388  0.4407012163
##  [961]  0.4055121595  0.4617547336  1.8650855190  0.1252918300 -0.2159410115
##  [966] -0.2814590714  0.1475370184  0.1449323635 -0.2314710761  0.0585648413
##  [971]  0.4326561211  1.1765912685  0.3816218868  0.8733030492  0.5178265025
##  [976] -0.0249475839  1.1076678841  0.1928292912  0.1963603432  1.3891156563
##  [981]  0.0577214073  0.7646150795  0.4405300861  0.4683248040  0.9420628725
##  [986]  0.1205554516  0.1713898801  0.0604711512 -0.1722087458  0.1789741713
##  [991]  0.2811294154 -0.2040973367  0.7678047439  0.5165497451  0.3739245630
##  [996]  0.8642183680  0.9068573259  0.7263593908 -0.2932005045  0.1915344435
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
##   0.56875480   0.33322456 
##  (0.10537486) (0.07450866)
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
## [1] -0.07255453  0.61874953  0.34910445  0.87354255 -0.08915397  0.10398248
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
## [1] 0.0222
```

```r
sd(results2$thetastar)  #the standard deviation of the theta stars is the SE of the statistic (in this case, the mean)
```

```
## [1] 0.9173706
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
## t1*      4.5 0.01351351   0.8938711
```

One of the main advantages to the 'boot' package over the 'bootstrap' package is the nicer formatting of the output.

Going back to our original code, lets see how we could reproduce all of these numbers:


```r
table(sample(x,size=length(x),replace=T))
```

```
## 
## 2 3 4 5 8 
## 1 3 3 2 1
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
## [1] -0.0215
```

```r
se.boot
```

```
## [1] 0.891225
```

Why do our numbers not agree exactly with those of the boot package? This is because our estimates of bias and standard error are just estimates, and they carry with them their own uncertainties. That is one of the reasons we might bother doing jackknife-after-bootstrap.

The 'boot' package has a LOT of functionality. If we have time, we will come back to some of these more complex functions later in the semester as we cover topics like regression and glm.

