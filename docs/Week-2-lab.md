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
## 1 3 4 5 6 8 
## 2 3 1 1 1 2
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
## [1] 0.0222
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
## [1] 2.666997
```

```r
UL.boot
```

```
## [1] 6.377403
```

Method #2: Simply take the quantiles of the bootstrap statistics


```r
quantile(xmeans,c(0.025,0.975))
```

```
##   2.5%  97.5% 
## 2.6975 6.5000
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
##    [1] 4.6 6.1 4.6 2.2 5.9 3.1 6.2 4.8 4.8 4.5 3.9 4.1 7.0 4.4 4.7 4.1 5.3 4.3
##   [19] 5.4 3.3 3.2 4.3 2.2 5.5 3.3 3.8 5.0 4.6 3.2 4.7 5.0 4.8 4.9 5.2 4.2 4.9
##   [37] 5.5 3.8 4.5 3.9 4.0 4.7 3.1 5.7 5.0 2.4 3.6 4.7 5.1 4.0 5.3 3.9 4.9 2.8
##   [55] 5.1 5.0 3.7 4.1 3.8 3.3 4.9 4.9 4.2 6.2 4.6 4.1 5.4 4.9 5.3 4.6 4.3 2.9
##   [73] 2.8 5.5 3.0 4.0 5.0 6.4 4.9 5.4 3.4 4.2 4.4 4.5 6.1 4.5 3.9 5.6 4.3 3.9
##   [91] 3.7 4.7 4.6 6.1 3.4 5.3 3.4 3.7 4.4 4.9 5.1 5.5 3.9 4.3 5.3 5.5 4.0 4.0
##  [109] 4.0 3.7 2.6 3.9 4.9 5.5 5.3 5.3 5.4 4.2 2.4 4.3 4.2 3.4 3.7 3.2 4.1 4.7
##  [127] 5.8 3.7 3.3 4.1 2.3 4.3 6.2 5.3 5.9 4.9 4.7 5.1 5.6 4.7 3.9 4.9 3.7 3.6
##  [145] 4.0 4.0 4.3 4.9 3.0 4.4 4.4 4.7 2.5 6.5 3.2 3.0 4.3 4.4 5.2 3.9 2.7 3.6
##  [163] 6.3 5.9 1.3 4.4 1.6 6.6 3.7 4.0 3.3 5.3 3.5 4.5 4.1 4.4 5.0 5.6 4.5 5.0
##  [181] 4.2 4.7 5.7 4.1 5.0 4.6 4.2 2.8 5.8 3.3 5.7 5.5 3.8 4.1 5.1 3.8 4.9 3.5
##  [199] 4.7 4.8 4.1 5.8 3.9 4.2 5.3 5.6 3.4 4.6 5.3 3.7 5.3 3.4 3.2 5.7 5.7 4.9
##  [217] 4.5 3.2 4.5 4.3 5.1 4.4 4.1 4.8 3.5 5.0 4.4 5.0 4.0 3.8 4.9 4.6 5.7 4.3
##  [235] 5.2 3.8 4.0 3.8 4.7 4.4 4.5 5.1 4.4 4.2 4.8 5.9 4.3 5.4 5.8 3.2 7.0 6.7
##  [253] 5.0 6.3 3.0 1.6 5.2 4.2 5.4 4.0 5.2 6.4 5.7 5.8 4.2 3.8 4.7 3.5 5.1 3.4
##  [271] 3.5 6.3 4.3 3.2 5.4 5.1 4.9 5.6 5.0 5.2 4.9 4.1 4.2 4.6 3.2 3.6 3.3 4.7
##  [289] 6.0 3.3 4.6 2.8 3.3 4.3 4.1 4.1 5.1 4.5 4.9 4.1 2.9 4.4 4.0 3.7 5.3 3.2
##  [307] 2.7 2.7 4.2 5.3 3.0 3.7 5.0 4.6 5.1 4.6 4.9 2.6 3.9 3.8 5.0 3.2 4.4 5.5
##  [325] 3.7 4.7 5.7 5.1 5.1 5.9 5.5 4.6 5.2 4.8 4.5 4.4 5.9 4.8 4.5 3.5 5.5 6.4
##  [343] 4.6 3.5 4.2 5.4 6.3 5.1 3.2 3.5 2.9 5.9 4.6 5.1 4.6 5.4 4.5 3.3 3.9 4.7
##  [361] 6.2 4.0 4.6 4.3 5.3 3.0 4.8 3.8 6.2 4.1 4.4 4.3 5.3 4.0 4.4 2.8 5.9 4.4
##  [379] 6.5 4.8 4.8 4.5 5.1 4.8 4.3 5.4 5.3 3.8 5.9 5.9 4.8 4.5 4.4 5.2 3.7 4.2
##  [397] 4.8 3.7 5.1 5.0 4.4 5.2 4.5 2.9 5.1 4.4 6.2 3.9 3.9 4.7 4.1 4.4 4.5 3.4
##  [415] 4.8 4.8 5.2 4.3 5.8 4.0 5.1 3.1 6.1 4.7 4.7 3.4 4.9 2.5 4.5 3.2 4.9 4.0
##  [433] 5.2 3.9 3.7 2.3 5.6 6.6 3.9 3.5 3.2 5.0 5.6 4.4 4.7 5.9 5.0 3.7 5.4 4.6
##  [451] 5.6 3.3 4.6 4.3 3.8 4.2 5.5 4.6 6.1 3.2 5.3 4.5 6.2 4.3 5.8 5.1 5.5 3.9
##  [469] 4.0 4.9 2.9 4.5 3.5 4.5 5.2 5.6 4.3 3.1 5.9 4.0 4.5 6.5 4.3 4.1 3.7 4.9
##  [487] 4.0 5.2 4.8 4.6 4.0 3.4 4.8 6.8 5.8 4.9 5.7 4.0 5.8 6.3 5.4 4.7 4.1 2.8
##  [505] 4.5 4.7 3.5 5.5 4.8 5.9 4.7 4.5 5.0 6.4 4.2 4.3 3.3 5.0 2.4 5.7 4.2 4.0
##  [523] 3.4 4.0 5.0 5.3 4.1 2.6 4.4 4.4 5.7 5.5 3.5 2.1 4.1 4.5 3.8 4.5 4.2 4.4
##  [541] 6.1 4.3 4.0 2.8 5.7 4.7 4.5 6.5 3.8 4.4 4.4 4.5 4.0 4.8 4.1 5.3 3.3 5.1
##  [559] 5.5 5.5 4.6 4.9 5.4 4.7 5.7 3.8 5.9 6.0 4.7 5.0 5.0 4.7 4.5 4.5 4.0 4.3
##  [577] 4.6 4.6 4.1 4.5 6.0 3.3 5.4 3.2 3.1 4.6 5.1 5.0 5.5 4.1 4.7 4.4 6.5 4.8
##  [595] 4.2 4.0 4.6 4.5 3.2 3.4 6.8 2.4 4.3 4.1 4.3 4.8 5.5 4.7 5.4 4.6 5.1 3.0
##  [613] 4.9 5.6 6.6 6.4 3.5 3.6 3.9 5.3 3.8 3.8 6.6 4.6 3.7 3.6 4.4 3.1 3.4 5.9
##  [631] 4.6 4.8 3.1 5.6 4.7 3.8 5.3 4.2 3.7 3.4 4.9 3.7 3.1 3.8 3.9 3.5 4.0 6.1
##  [649] 3.9 5.7 5.3 4.5 2.7 3.9 4.1 6.7 5.3 3.9 5.3 4.7 5.9 4.6 5.4 2.9 5.8 5.0
##  [667] 5.7 3.9 5.6 3.6 5.1 4.6 3.1 4.3 4.2 5.0 4.9 4.4 5.0 4.7 5.0 4.9 3.1 3.9
##  [685] 3.4 3.5 4.8 5.3 4.3 6.3 5.5 5.1 4.7 4.2 5.6 3.9 3.4 4.2 5.5 5.1 5.1 4.8
##  [703] 4.1 5.2 4.5 6.0 4.9 5.2 5.0 2.6 4.3 4.1 3.9 3.1 5.1 4.7 5.9 4.1 5.9 4.1
##  [721] 3.6 4.4 2.4 4.9 4.4 5.3 3.8 3.4 4.4 5.0 4.8 4.3 4.8 4.9 3.1 3.5 5.1 5.1
##  [739] 5.9 5.6 4.2 4.7 4.2 5.0 3.9 4.3 4.1 5.8 5.1 4.8 2.9 3.9 4.6 3.1 2.8 4.0
##  [757] 5.1 3.6 3.3 5.4 4.5 4.6 4.0 4.3 5.3 5.1 4.1 5.3 5.3 4.0 6.1 6.1 3.9 2.0
##  [775] 5.4 5.9 4.3 5.5 3.5 3.6 4.6 3.4 4.4 3.7 4.9 3.8 5.4 4.6 3.8 4.2 5.5 4.0
##  [793] 4.6 3.8 5.9 3.4 4.9 5.0 4.0 3.9 5.4 4.9 4.3 5.2 5.9 4.7 4.7 3.9 3.1 4.0
##  [811] 4.2 5.1 5.2 4.2 4.4 4.9 3.4 6.6 4.6 4.9 5.0 4.9 6.6 5.9 5.9 2.8 3.1 5.7
##  [829] 4.2 5.5 3.8 6.8 4.8 4.8 3.3 3.4 4.9 3.2 4.4 6.3 3.2 4.3 4.2 3.5 3.3 4.5
##  [847] 3.6 3.9 4.7 4.4 3.5 4.6 5.0 5.0 4.3 4.5 5.0 3.2 4.0 3.9 5.3 3.7 4.7 5.2
##  [865] 4.8 3.9 4.9 2.3 5.7 4.4 4.6 5.2 3.2 4.7 5.3 5.5 4.7 4.8 5.4 6.1 3.7 3.4
##  [883] 4.5 3.3 3.2 4.8 2.8 2.2 5.3 5.0 5.3 5.4 3.9 4.3 4.2 4.3 4.4 4.8 3.9 6.9
##  [901] 4.0 5.7 4.3 3.5 3.4 4.8 1.5 5.2 4.5 3.4 4.3 4.6 5.1 5.3 4.7 3.5 4.6 2.4
##  [919] 4.5 5.5 5.0 2.6 4.3 5.3 5.3 4.5 3.9 3.7 4.9 4.3 4.0 4.6 3.8 4.2 4.0 1.8
##  [937] 4.9 5.5 5.9 3.9 5.8 4.4 5.5 3.1 5.7 4.3 4.7 6.3 5.2 5.4 2.9 3.9 5.1 5.5
##  [955] 4.0 4.3 4.3 5.0 5.1 4.8 5.3 3.5 5.0 4.4 4.6 5.3 6.3 3.7 2.9 5.3 4.3 3.8
##  [973] 6.3 3.0 5.1 4.3 3.8 5.0 3.7 4.9 4.5 5.7 5.1 3.8 4.9 4.6 3.6 4.9 5.3 5.4
##  [991] 4.6 4.3 3.3 3.7 4.5 4.4 4.9 4.0 5.8 5.5
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
##   2.6   6.3
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
##    [1] 4.2 4.3 5.1 4.4 4.6 5.0 2.8 3.9 4.3 4.9 5.3 4.7 5.0 3.2 3.4 6.3 2.8 3.3
##   [19] 4.6 3.2 3.9 3.4 5.4 5.8 4.3 4.2 3.8 4.7 5.9 4.1 3.8 3.7 3.1 5.2 4.9 4.6
##   [37] 5.2 3.0 5.9 5.5 3.7 4.5 4.1 3.6 4.1 4.0 5.0 4.3 5.2 3.9 6.6 2.8 4.3 4.7
##   [55] 5.2 3.7 4.7 2.4 4.5 4.8 4.5 5.0 3.5 6.0 3.6 5.8 4.4 3.5 6.7 2.9 4.3 4.0
##   [73] 5.6 3.8 4.8 5.3 3.4 5.0 4.4 4.7 5.4 5.4 4.5 5.2 6.0 4.1 4.9 5.4 6.0 4.3
##   [91] 4.7 3.3 6.2 4.2 3.3 3.5 4.4 3.9 4.7 5.0 4.8 4.2 4.7 5.3 3.3 4.4 4.4 3.9
##  [109] 4.0 5.7 4.8 4.0 4.4 4.5 3.7 4.0 3.4 4.6 5.5 5.5 4.8 3.9 5.5 4.6 4.6 4.0
##  [127] 3.8 4.7 5.0 4.8 4.3 4.0 5.9 4.4 4.6 4.6 3.5 4.5 4.9 3.2 4.3 5.6 4.7 5.2
##  [145] 5.0 4.2 5.4 6.2 4.0 4.6 4.0 4.9 4.2 3.5 4.5 5.0 4.7 5.2 4.2 3.8 4.4 5.1
##  [163] 5.3 6.1 3.3 4.9 5.3 3.8 3.2 3.5 2.1 5.8 4.4 5.1 2.2 4.0 5.6 3.7 4.1 4.2
##  [181] 4.9 5.0 5.0 5.2 4.1 6.0 5.6 3.6 4.5 4.2 5.1 3.1 3.7 3.8 4.0 6.1 5.6 6.4
##  [199] 3.7 3.7 4.5 4.8 4.3 4.5 5.0 4.6 3.5 5.3 3.9 6.1 3.7 4.0 3.9 4.6 6.2 4.0
##  [217] 5.2 3.3 3.9 4.7 5.1 5.9 4.7 5.3 4.6 2.4 5.2 5.3 6.0 3.9 3.9 4.5 3.7 4.3
##  [235] 4.2 4.3 2.6 6.3 2.4 3.6 5.5 4.9 5.0 4.5 4.1 5.6 4.6 3.2 3.0 4.5 4.1 4.8
##  [253] 4.0 3.5 4.3 6.0 5.4 5.0 4.1 4.3 4.7 6.0 3.6 5.1 4.7 4.8 5.7 5.9 3.9 4.9
##  [271] 5.4 4.1 2.9 5.0 5.4 4.5 5.4 4.2 4.6 2.6 4.1 4.5 4.9 4.3 3.3 4.1 4.8 4.9
##  [289] 3.0 4.2 4.6 5.5 2.9 4.3 6.3 4.2 4.5 5.4 3.8 4.3 4.2 4.4 5.8 4.9 5.2 4.5
##  [307] 3.7 3.6 4.6 3.7 3.6 4.6 6.2 5.1 4.5 4.0 4.9 5.5 4.9 5.7 4.6 4.3 4.8 4.4
##  [325] 5.7 3.5 4.0 3.4 4.8 5.1 4.9 4.9 4.7 2.8 3.9 3.5 4.5 5.0 4.1 3.8 5.4 2.7
##  [343] 4.4 4.4 3.7 4.9 3.7 4.4 4.1 4.7 4.3 3.3 4.3 4.2 4.5 4.5 3.4 4.9 5.0 4.0
##  [361] 5.0 4.6 5.0 3.5 5.3 6.3 5.0 3.3 4.5 4.0 4.0 4.6 3.8 4.3 5.1 5.8 5.8 6.2
##  [379] 3.7 5.0 5.6 4.1 4.0 5.0 4.0 2.7 4.6 5.8 4.7 3.2 5.9 4.3 3.5 3.8 4.4 4.5
##  [397] 3.6 3.8 3.9 3.9 5.9 4.2 4.7 4.8 3.6 5.3 6.8 5.0 4.6 5.0 4.3 6.7 4.1 6.0
##  [415] 4.3 4.5 4.5 4.2 3.3 4.8 5.1 4.3 4.4 4.5 5.5 4.6 5.3 4.1 5.7 4.9 4.1 2.4
##  [433] 4.3 5.0 3.4 4.5 3.6 4.8 3.3 5.1 3.5 5.4 4.5 3.0 5.1 3.8 4.2 4.0 4.6 4.0
##  [451] 2.5 5.6 5.7 2.4 5.9 3.1 3.3 5.7 3.7 5.8 4.3 2.7 3.5 4.5 6.0 3.7 4.6 5.5
##  [469] 3.8 3.4 3.7 3.4 4.4 4.9 3.9 4.4 5.0 3.7 3.1 2.9 4.8 5.5 5.9 4.9 5.6 4.3
##  [487] 4.6 4.3 5.5 4.9 3.8 5.9 4.9 5.0 2.9 4.2 3.6 3.2 5.0 4.1 5.5 4.2 2.9 3.3
##  [505] 4.2 3.6 6.0 5.4 5.7 5.2 5.3 5.3 4.4 4.2 3.9 4.5 3.6 4.4 3.7 4.4 4.1 4.1
##  [523] 6.4 5.0 3.7 4.5 4.5 4.0 3.9 3.5 5.7 7.0 5.4 6.0 3.0 4.4 5.2 5.6 3.5 4.4
##  [541] 5.2 4.3 5.2 5.6 3.8 6.0 4.3 4.0 3.0 4.0 5.3 3.8 4.3 4.2 4.9 3.7 4.8 5.3
##  [559] 3.4 4.2 5.6 4.2 3.1 4.6 5.3 3.3 3.7 4.8 4.3 5.2 5.2 3.6 3.6 3.9 3.0 3.3
##  [577] 5.8 5.1 4.0 5.1 4.5 3.5 4.5 6.8 4.9 3.8 4.8 3.8 4.7 5.3 5.2 2.8 5.4 5.2
##  [595] 3.4 3.8 5.2 4.3 3.3 4.0 3.9 6.2 4.3 4.8 5.0 4.7 6.5 4.6 5.0 4.3 5.2 4.3
##  [613] 4.5 5.5 5.1 4.4 2.7 3.4 6.4 4.4 3.9 4.1 4.2 4.1 5.7 5.0 3.1 5.4 3.6 4.3
##  [631] 4.2 3.9 3.7 4.9 4.3 5.7 5.2 5.8 3.1 5.7 4.2 4.1 3.4 4.4 4.3 3.6 4.5 4.1
##  [649] 3.2 5.0 3.8 4.5 5.0 2.7 5.5 4.0 2.6 3.6 4.7 3.6 5.7 4.9 4.4 5.6 3.9 5.5
##  [667] 4.3 4.3 4.8 4.1 3.5 5.2 6.4 3.3 2.8 5.7 4.4 6.1 4.2 4.6 4.9 5.4 4.4 4.0
##  [685] 2.6 4.8 3.6 6.6 5.0 3.3 4.6 3.2 5.8 5.2 5.9 4.2 4.6 5.6 3.8 4.3 4.7 6.5
##  [703] 5.1 4.5 4.4 3.2 3.8 4.3 6.3 3.3 5.5 4.9 3.5 4.7 3.8 4.6 3.8 4.2 3.2 5.2
##  [721] 3.2 6.2 4.8 5.7 4.7 4.6 5.0 3.8 3.3 5.5 6.5 4.2 5.2 5.3 3.8 4.5 4.7 3.9
##  [739] 4.4 4.1 5.0 3.9 5.2 3.6 3.2 4.2 4.7 5.2 6.3 5.5 4.4 4.7 5.2 3.7 3.9 4.5
##  [757] 3.5 3.9 4.5 3.4 4.7 5.3 4.6 5.1 4.1 3.8 4.5 4.7 4.1 4.8 3.1 4.3 4.6 3.5
##  [775] 4.8 5.4 6.1 4.4 5.3 6.6 3.8 4.3 5.2 5.9 3.8 3.5 5.0 4.5 5.1 3.3 4.8 2.7
##  [793] 4.9 4.7 3.6 4.1 5.2 3.8 5.8 5.2 3.4 4.0 4.6 6.5 3.5 4.7 6.2 5.2 6.9 5.4
##  [811] 2.8 3.2 3.8 4.0 4.1 4.7 4.6 4.0 4.2 5.9 4.5 3.9 5.1 4.2 5.4 4.8 4.9 5.5
##  [829] 5.1 4.0 4.4 4.8 3.0 5.5 4.5 3.8 3.2 5.3 4.3 5.5 5.8 5.6 4.8 4.2 3.2 3.2
##  [847] 4.7 4.5 5.0 6.3 4.3 2.7 3.1 4.0 4.2 4.7 3.7 4.9 2.9 4.6 4.9 6.2 3.9 3.9
##  [865] 4.3 4.2 4.8 4.1 4.8 3.5 4.0 3.5 4.2 2.8 3.8 3.1 4.8 3.8 4.8 5.8 5.6 5.1
##  [883] 6.0 5.1 5.7 5.3 5.3 3.6 4.5 5.0 5.2 6.0 5.2 3.7 4.1 4.0 5.1 5.8 4.9 3.2
##  [901] 3.7 4.4 3.7 5.5 4.2 3.3 4.5 3.9 5.6 5.7 5.8 4.9 3.8 4.2 4.3 4.5 4.8 5.0
##  [919] 3.5 3.7 3.6 3.4 4.6 3.8 4.7 4.2 4.2 5.2 3.8 4.1 4.9 4.1 4.4 5.2 4.8 5.9
##  [937] 4.8 5.3 7.1 5.5 3.7 4.8 4.6 6.1 2.8 3.3 5.4 3.5 4.9 3.7 5.3 4.5 4.0 2.6
##  [955] 4.5 6.1 4.5 5.3 5.0 4.6 3.8 4.8 5.8 5.9 4.3 4.9 5.1 5.7 5.3 5.9 4.7 3.3
##  [973] 4.3 3.3 4.9 4.8 4.5 5.1 4.2 4.2 6.6 4.3 5.2 6.4 4.9 3.2 4.0 4.7 3.7 4.9
##  [991] 4.7 5.6 5.9 4.5 4.7 3.6 5.0 4.6 5.8 2.5
## 
## $func.thetastar
## [1] 0.0088
## 
## $jack.boot.val
##  [1]  0.466268657  0.373387097  0.277714286  0.181216931  0.005780347
##  [6] -0.063128492 -0.110588235 -0.244347826 -0.366944444 -0.489940828
## 
## $jack.boot.se
## [1] 0.9086322
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
##    [1] 5.8 3.0 4.5 4.8 4.3 3.9 3.5 4.7 3.2 3.8 6.5 4.7 4.6 4.2 4.9 4.7 5.5 3.8
##   [19] 3.7 5.0 6.0 4.9 4.2 4.2 4.7 4.5 3.9 4.8 3.3 4.4 3.7 5.1 5.2 4.2 3.6 3.5
##   [37] 3.1 4.1 4.7 5.1 3.6 3.8 3.7 6.0 2.6 5.1 5.7 5.2 5.3 5.4 3.1 5.6 2.9 6.1
##   [55] 6.4 3.1 5.1 5.0 4.5 3.5 4.0 4.0 4.2 4.6 6.6 3.5 3.6 3.2 5.1 3.9 3.7 5.8
##   [73] 3.8 3.8 4.8 4.1 4.4 4.5 5.1 4.2 4.4 4.9 3.0 3.3 3.1 5.7 2.8 5.4 6.0 4.1
##   [91] 5.5 4.6 6.5 5.3 4.0 4.6 6.0 3.8 4.6 3.6 5.4 4.4 5.0 3.4 5.0 5.6 4.6 5.1
##  [109] 5.1 5.8 5.0 3.3 5.7 4.2 6.2 3.5 5.0 5.4 5.1 4.6 5.4 7.1 4.3 6.0 4.4 4.8
##  [127] 4.7 3.6 2.7 4.2 3.5 4.8 4.8 4.3 6.2 4.5 4.0 4.7 5.7 5.6 4.2 4.3 4.0 4.6
##  [145] 5.2 4.3 2.8 4.3 2.8 5.2 3.4 4.3 3.7 5.5 2.6 4.4 3.2 5.1 6.2 4.7 3.9 4.7
##  [163] 5.0 4.3 3.3 4.0 5.7 3.3 6.7 4.4 6.2 5.2 3.1 5.3 2.9 6.4 4.0 5.2 4.3 3.8
##  [181] 4.5 2.7 2.0 5.5 3.9 3.6 5.3 7.1 3.8 5.7 4.0 4.7 4.4 3.7 4.5 4.8 5.0 4.2
##  [199] 5.5 4.3 2.4 4.3 4.0 3.9 4.9 5.3 3.9 2.1 5.1 4.4 5.3 5.6 4.3 3.1 6.3 3.0
##  [217] 3.6 4.9 3.1 4.6 4.4 4.3 3.3 4.5 4.0 4.6 4.9 5.0 4.8 3.5 4.9 3.0 4.3 5.9
##  [235] 3.6 4.9 4.4 4.0 4.9 4.7 5.1 3.7 4.2 3.9 4.2 3.7 5.1 4.5 3.4 3.4 4.2 3.5
##  [253] 3.7 5.2 4.2 4.3 4.0 3.8 3.7 3.9 4.2 3.3 4.3 2.9 4.3 5.8 3.4 5.4 5.1 4.2
##  [271] 4.4 3.5 5.1 4.5 2.8 4.3 4.8 4.3 4.7 4.1 4.1 6.0 4.3 5.3 4.5 5.6 5.0 6.0
##  [289] 4.8 4.2 2.5 4.3 3.8 3.4 3.6 5.7 3.8 4.7 4.4 3.7 4.5 4.7 5.2 4.4 5.9 6.2
##  [307] 5.7 4.2 4.0 5.5 4.1 4.2 3.7 3.1 4.7 5.0 4.3 4.5 6.1 4.5 5.5 4.6 4.7 5.7
##  [325] 6.0 4.6 4.6 3.8 5.4 4.5 3.0 4.9 4.7 5.4 3.9 4.6 5.0 3.3 3.9 4.1 4.3 5.6
##  [343] 5.1 3.9 4.1 3.9 7.0 4.2 4.4 5.7 3.7 5.1 4.7 4.5 5.3 5.5 5.2 5.8 3.6 6.0
##  [361] 2.7 5.4 6.5 4.2 4.2 4.5 3.1 4.9 5.2 5.0 6.5 3.5 3.5 5.5 3.9 5.1 4.0 3.9
##  [379] 3.2 3.7 3.9 4.7 5.2 5.0 4.4 3.9 3.9 4.6 3.6 3.2 3.9 3.4 2.9 3.7 5.2 5.2
##  [397] 2.3 4.5 5.6 5.8 4.1 4.3 4.9 2.8 4.8 3.9 4.4 3.3 3.8 4.3 5.4 3.9 4.4 4.3
##  [415] 6.6 4.6 5.3 3.4 4.7 2.8 4.6 5.0 5.2 4.0 5.1 6.3 5.8 2.2 5.1 4.3 4.3 5.2
##  [433] 4.2 5.4 4.5 4.1 4.3 5.1 4.2 4.0 4.1 4.3 3.7 5.1 3.3 3.0 3.2 4.8 4.4 5.0
##  [451] 5.0 4.5 5.6 5.5 3.4 3.9 4.8 5.1 4.0 5.5 6.4 5.9 5.2 4.7 4.2 4.7 3.1 4.9
##  [469] 4.4 5.1 4.3 4.7 5.6 3.8 6.2 2.8 4.7 3.7 3.3 5.6 2.4 4.8 4.4 4.2 3.9 4.3
##  [487] 4.4 5.6 4.6 3.2 5.7 4.4 4.9 6.2 4.3 4.6 5.3 4.4 5.3 4.7 3.0 3.7 3.4 6.1
##  [505] 4.4 5.3 4.0 6.6 5.9 5.3 2.9 4.0 4.9 3.3 3.7 4.8 5.0 3.6 5.1 4.9 4.2 4.2
##  [523] 3.7 3.8 5.0 4.4 5.7 4.1 3.7 5.4 5.0 5.3 5.4 4.5 4.8 3.9 3.9 3.2 6.7 3.8
##  [541] 5.3 4.1 3.8 5.3 5.1 5.1 4.6 4.8 4.1 3.4 6.9 6.1 5.0 2.9 4.8 6.3 6.0 4.8
##  [559] 5.5 4.7 5.7 3.9 5.2 4.8 5.4 2.7 5.2 5.3 5.2 5.8 5.1 6.8 4.8 6.0 4.6 4.5
##  [577] 4.3 4.8 4.1 3.6 3.9 3.6 4.5 4.5 3.6 4.6 4.5 2.7 4.3 5.0 4.5 3.1 5.5 3.5
##  [595] 6.4 4.7 4.2 3.1 3.2 5.5 6.4 6.5 4.0 4.0 3.8 3.4 2.9 3.7 4.8 3.4 2.9 4.6
##  [613] 3.7 5.6 4.6 2.1 6.5 5.3 3.5 4.9 4.5 4.5 5.6 5.0 4.6 5.5 4.3 4.8 5.4 4.7
##  [631] 3.9 3.4 5.6 6.0 5.1 4.6 4.4 3.3 4.3 4.8 2.8 5.5 4.5 4.1 4.7 3.9 5.0 4.3
##  [649] 5.2 5.3 6.3 4.9 5.3 5.4 3.4 5.9 4.8 5.1 4.8 2.3 3.2 5.7 4.6 4.5 4.4 4.6
##  [667] 5.5 6.0 5.1 5.7 3.8 5.8 2.6 4.8 5.4 4.6 3.2 5.4 5.1 6.0 5.1 5.7 4.2 4.4
##  [685] 4.0 4.9 3.9 5.2 5.9 6.1 4.9 4.2 4.9 5.1 3.3 5.1 3.5 6.1 4.0 4.2 6.2 2.4
##  [703] 4.8 5.1 5.1 6.0 4.4 5.2 4.9 5.1 3.3 2.6 3.8 3.6 4.5 5.4 5.2 3.7 3.2 4.6
##  [721] 4.0 2.7 3.5 3.6 3.9 5.2 3.6 5.0 4.7 4.8 5.2 4.3 4.8 4.5 4.2 5.9 5.7 3.5
##  [739] 3.7 6.1 2.7 5.4 3.9 4.6 4.9 5.4 5.7 4.4 6.4 3.8 1.1 3.7 5.9 4.1 5.8 3.2
##  [757] 3.6 3.6 4.7 4.6 5.1 3.7 3.4 4.8 6.4 4.5 4.3 5.1 4.2 5.1 5.3 3.8 5.1 3.4
##  [775] 5.8 4.6 2.1 5.9 3.8 4.6 4.8 5.3 5.3 3.4 4.8 5.8 6.0 3.2 5.3 4.1 5.8 3.6
##  [793] 4.7 5.0 5.5 4.0 4.7 5.0 5.3 3.1 6.1 4.0 4.3 6.2 5.5 4.2 4.0 4.6 4.1 4.6
##  [811] 5.5 5.7 5.1 3.0 5.2 5.9 6.0 6.3 3.4 5.2 3.3 2.9 3.6 3.5 4.6 4.6 6.2 5.0
##  [829] 3.7 4.5 2.9 3.7 4.9 3.7 6.9 4.5 4.7 4.3 4.8 5.8 2.9 5.2 4.6 4.1 4.6 5.4
##  [847] 4.0 4.5 5.5 2.7 5.3 5.4 3.6 3.2 3.8 4.2 4.1 5.5 5.6 4.2 5.9 3.5 4.8 6.3
##  [865] 3.6 3.6 4.4 5.0 4.3 3.8 3.4 6.1 4.6 5.3 4.6 4.5 5.1 5.4 4.5 4.0 4.8 4.0
##  [883] 4.2 6.1 5.4 4.2 5.3 4.7 4.3 5.3 5.3 5.5 5.5 4.9 4.5 4.4 5.6 3.6 5.7 4.5
##  [901] 6.1 4.5 3.9 5.4 5.9 5.3 3.7 4.4 5.8 6.2 5.4 3.8 4.6 3.0 5.3 5.0 4.9 4.3
##  [919] 4.4 3.3 3.4 5.3 4.3 3.7 4.0 3.2 4.8 5.6 3.9 4.7 4.8 2.9 5.7 4.7 3.6 3.9
##  [937] 5.1 5.2 6.0 5.0 5.7 6.9 5.1 4.3 3.7 4.1 4.3 4.3 4.3 4.3 4.9 5.0 5.7 5.3
##  [955] 4.3 4.5 6.1 5.5 4.4 6.7 4.2 4.2 5.3 4.7 5.7 5.2 3.3 5.6 2.5 5.4 4.9 5.5
##  [973] 4.4 5.0 4.0 3.6 3.6 4.1 3.4 3.0 3.5 5.5 3.3 4.6 4.1 5.3 4.5 5.2 4.0 3.7
##  [991] 2.6 3.3 4.4 4.1 5.9 5.3 3.8 3.4 3.8 3.4
## 
## $func.thetastar
## 72% 
## 5.1 
## 
## $jack.boot.val
##  [1] 5.6 5.5 5.3 5.3 5.2 5.1 5.0 4.8 4.6 4.6
## 
## $jack.boot.se
## [1] 0.9949874
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
## [1] 1.362973
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
##   2.425850   4.467876 
##  (1.019118) (2.084698)
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
## [1]  0.6677042  0.6068928  0.6834122  0.1075979  0.1765846 -0.1267794
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
##    [1]  0.400776962  1.296728050 -0.208674880  2.131568855  0.814066990
##    [6]  1.537956715  0.632422072  0.912311693  2.420744905  1.352825044
##   [11]  1.429926530  1.245782841  2.152177656  1.323500595  1.173032920
##   [16]  1.413202837  0.907027761  0.735886784 -0.170853916  2.116421831
##   [21]  0.089105704  0.654184419 -0.428925376  1.351731918 -0.795938344
##   [26]  2.421185216  0.883546510 -0.162394682  1.382542327 -0.152881653
##   [31]  0.927629495 -0.664192934  1.299312948  0.907027761  1.826487328
##   [36]  1.897186429  1.421809561  0.833549888  0.853637680  0.907027761
##   [41]  0.782629038  0.699098963  1.366327922  0.850124464  0.817340178
##   [46]  0.924860914  2.097701388  2.516234806  0.747859886  0.523512337
##   [51]  2.167534618  1.080957647  1.712356706  1.419767744  1.282371663
##   [56] -0.064620602  1.356589937  1.370058699  0.367541711  0.448047287
##   [61]  1.200410246  0.287277487  0.885210140 -0.808925585  2.011189327
##   [66]  2.198491028  0.784144889  1.317791810  2.432174734  1.327227623
##   [71]  0.105826592  1.342571624  1.722253238 -0.272253752  1.655502221
##   [76]  2.241706414  1.650014173  0.855791156  0.888021437  1.600035014
##   [81]  1.110479905  2.110193543  1.260066337  2.310367379  1.403881101
##   [86]  1.456829411  1.329405256  0.795387384 -0.372590869  1.279123745
##   [91]  1.232168489  1.465426961  1.443083952  0.944087750  1.660849856
##   [96]  1.243499427  2.360369789 -0.039597846  1.941326900  0.910403754
##  [101]  1.962376126  1.313774995  0.908894148  2.047296449  0.367633487
##  [106]  0.838642257  2.147518781 -0.277335131  0.434102850  0.394522721
##  [111]  1.240815970  1.377923356  1.079056270  0.860076365  1.339812954
##  [116]  0.796708517  2.391111153  1.436309801  0.212399938  1.257495776
##  [121]  2.022711297  2.197016843  1.365915641  1.379531855  0.125993038
##  [126]  0.888373471  1.348477768  2.061102884  1.475166870  0.784380052
##  [131]  2.275569329  1.540700173 -0.556456261  1.780785685  1.466639286
##  [136]  1.301807052  1.309956018  0.793850363  0.402084684  1.329576584
##  [141]  1.162702314  2.090768200  0.386020611  1.833892285  2.209780016
##  [146]  2.247841591  0.143941614  1.807584495  2.093338433  0.417443337
##  [151]  1.378343877  1.388284242  1.163946181  0.381543369  1.866308670
##  [156]  0.408722358  0.808587877  1.683480973  0.801116567 -0.030355185
##  [161]  1.437783880  1.239276645  2.238022763  2.212586766  1.307066897
##  [166]  0.786175260  1.916980825  2.041755766  1.397055964  1.346095630
##  [171]  1.336720237  1.309023290  1.808421545  2.085589269  0.362631180
##  [176]  2.193299007  0.849836063  0.830370883  2.044815151  1.778985061
##  [181]  1.798349023 -0.310208758  1.300272777  0.496268871 -0.132085500
##  [186]  2.432174734  1.358709179  2.247324133  0.830632669  1.156422607
##  [191]  1.847296824  0.791564018  0.811495217  1.406035995  2.119409978
##  [196]  1.087884004  2.553908322  1.340019087  1.322731846  1.086844862
##  [201]  2.303693522  2.186107935 -0.620665452  0.101925307  0.795727074
##  [206]  0.520012558  1.137509893 -0.548135308  1.379393576  1.333495827
##  [211]  2.085250824  1.358969784  1.370815425  0.358439132  1.489387720
##  [216]  0.365197891  0.803912937  0.874284412  0.447986847  1.669795962
##  [221]  1.847201799  2.222141842  2.028012863  1.368258747  1.426002817
##  [226] -0.186937654  1.387454245 -0.487435276  1.529576642  1.382134584
##  [231]  0.401705776  0.404485465  1.304232168  0.956159048  1.327195257
##  [236]  1.199421432 -0.141154815  1.349984678  1.373052069 -0.371905768
##  [241]  1.460276220  0.079045117  1.193521933  0.430527148  0.782841734
##  [246]  1.829300315  1.384949304  2.008258840  2.344362505  0.861758191
##  [251]  2.165938141  1.287545050  1.971011995  0.846064307  1.453190480
##  [256]  1.404240528  1.386472879  1.338504911 -0.361825879  0.134156697
##  [261]  1.328448881  1.908327173 -0.057153117  2.166371315  1.398573123
##  [266] -0.231948316 -0.422933129  0.838669234  1.937270017  0.888787949
##  [271]  0.827318940  2.371169164  2.306896327  1.329918426  0.866646968
##  [276]  0.235361375  2.351644966 -0.027539427  1.353417160  1.977248244
##  [281]  1.214219811  1.701027652  0.388499262  2.103042736  1.195222624
##  [286]  0.798800671  1.306386755  2.149432257  1.423538795  1.910749720
##  [291]  2.143263121  0.791251301  0.764402200 -0.447419036  1.991120556
##  [296] -0.072147164 -0.374295253  0.146351182 -0.367090604  2.343304949
##  [301] -0.403185895 -0.625457148  2.292710258  0.783443379  1.480061842
##  [306]  1.526943124  1.199970724 -0.118082248  2.118801922  1.421137983
##  [311]  1.447780702  1.417456650 -0.587631185  2.310360067  1.328685815
##  [316]  1.417525442  2.047261130  1.300613799  0.451927709 -0.316303151
##  [321]  0.768265217  0.889026610  1.885252029  0.064125136  0.318594517
##  [326]  1.628068331  0.244672393  1.343350000 -0.351710362  0.414742190
##  [331]  0.872247457  1.331473898  1.442484855 -0.986850129  0.414578833
##  [336]  2.249177367 -1.050554496  2.082781047  1.100516085  1.725260305
##  [341]  0.438812527  1.358773495  1.912348801  1.368428957  0.411633168
##  [346]  2.292826536  1.199992858  2.329846917  1.346095630  0.475256931
##  [351]  1.356727357 -0.320356021  1.404240528  0.884301042  0.460287573
##  [356]  2.134638531  0.525400252 -0.394171325  0.479621476  0.834588808
##  [361]  2.220071604  0.832862002  2.084512836  0.377934878 -0.711242081
##  [366]  0.622428870  1.296698269  0.799203282  1.393717550  1.362970581
##  [371]  2.184109985  1.352651070  0.927093657  0.848076375  2.258852485
##  [376]  2.178088706  2.118951119  2.466830150  0.052025531  0.781077844
##  [381]  1.315149744  1.197371319  1.390508794  2.310671079  1.376239677
##  [386]  1.367468634  1.714815448  1.857501839  2.327921538  2.147588313
##  [391]  0.834588808  2.108946046  0.629968287  1.134962533  2.262832581
##  [396]  0.760552008  2.312351599  0.828533506  0.420915431 -0.180973450
##  [401]  1.362070694  1.901912209  0.918200789 -0.213153857  0.833623525
##  [406]  1.343465372  0.469354969  2.378239462  0.417773544  2.055157847
##  [411]  1.291474709  0.849526852  1.384547221  1.323101829  1.208350147
##  [416]  2.368743410  1.128588918  1.375711323  0.008538915  2.125767695
##  [421]  1.360790840  0.869932900  1.516030065  1.296005640  1.433434570
##  [426]  0.820410170  1.165272803  0.328687499  0.848742263  1.446303987
##  [431]  0.817945547  1.958398048 -0.299658742  0.362104811  1.350345261
##  [436]  1.745034947  0.640578025 -0.322018299  0.166631181  2.037379498
##  [441]  1.758579055  2.166555026  2.020618488  0.926258216  1.170780206
##  [446]  1.441094023  1.244375682  1.467285773  2.258337675  2.325818716
##  [451]  1.401101014  0.833565128  1.340144285  1.230761217  1.061206303
##  [456] -0.824974587  1.844979943 -0.207301801 -0.168505663  2.349490476
##  [461]  1.419133744  2.110193543  1.291649292  1.371052630 -0.295655037
##  [466]  0.806763793  1.761233516  0.384770929  2.109664223  1.364607746
##  [471]  1.433972464 -0.145045674 -0.313740042  1.797893935  0.880374911
##  [476]  1.384588561  0.022152276  0.766365346  2.159385424  2.146288722
##  [481]  1.458311681  0.832400679  1.464482600  1.327208024  1.554385843
##  [486]  2.381644594  1.471520949  2.461622533  1.455050369  2.254574308
##  [491]  1.926590133  0.380890319  1.259437798  0.984137413  0.809610842
##  [496]  0.060632864  0.792389123  1.874826372  0.673310545  1.503028886
##  [501]  0.912947478  1.391216572  2.220011153  2.254524059  2.263257990
##  [506]  0.851683059  1.376262905  2.039494782  1.394729215  1.948391948
##  [511]  0.689390579  2.210272677  0.610821618  1.115243045  1.362975371
##  [516]  0.092994886  2.453294304 -0.427421293  0.636021839  0.469552861
##  [521]  1.256632474 -0.183677071  1.270583884  1.173032920  1.279543009
##  [526]  1.299150168  1.323099808  0.440775796  1.691681978  0.803971471
##  [531] -0.448196201 -0.713501534  1.235219037  0.781825943  0.813909132
##  [536]  0.823192631  1.912288813  2.090009785  2.001424748  0.844790412
##  [541] -0.725369992  0.257875982 -0.319726560  0.747749868  1.166356033
##  [546]  1.591324145  2.385580868 -0.384668665  1.096195465  1.930196014
##  [551]  0.830632669  0.791116819 -0.420423202  1.763734375  0.820095203
##  [556]  1.319326716  0.318905949  2.099449226  1.462466641 -1.183961756
##  [561]  2.009819208  2.034859607 -0.787554535  2.076268151  1.320576443
##  [566]  1.456145589  0.409826641 -0.088198439  1.160249650  0.841069399
##  [571]  1.341550042  0.841962084  0.372083243 -0.581021529  1.178980322
##  [576]  2.301655542  0.894243898  1.312914301  1.782700901  0.438832618
##  [581]  2.093338433  1.306578509 -1.298737249  0.428513991  1.834903184
##  [586]  1.495973085  1.538261162  2.216591680  2.011281083  2.179612479
##  [591]  0.080936291  0.880106143  2.158198591  0.963412624  2.261537348
##  [596]  0.913705462  0.511915221  0.799999046  0.426382291  2.225823017
##  [601]  1.237661272  2.271878774 -0.193959979  2.272481627  2.534636395
##  [606]  1.248004629  1.436181809  2.318749522  2.304771623  1.334712201
##  [611]  1.329624817  0.486767018  0.509445867  1.554733347  1.447445384
##  [616]  0.661600223  1.326676927  0.394824485  1.340173479  1.386560371
##  [621]  2.153871178  2.442862527  1.266029783  0.687079172  0.594781732
##  [626]  0.856309716  1.047627676  1.457239275  0.361612470  1.204439226
##  [631]  1.977120542  0.853049597 -0.199222856  2.456014671 -0.250910909
##  [636]  1.239020291  0.453835405 -0.431458928  2.497537047  0.892182723
##  [641]  2.493613081  1.374139696  1.242993663  0.852815902  1.144995292
##  [646]  0.784849375  2.034020551  0.096418227  1.304948493  1.909534014
##  [651]  0.813368112  1.531828072 -0.622961170  1.789220566  1.430913824
##  [656]  0.835266155  0.396275764  1.847691859  1.659119859  0.748887004
##  [661]  1.430546608 -1.112060545  0.803355945  2.095150069  0.888624900
##  [666]  0.408232241 -0.514661687  0.796070862  1.897167794  0.832400679
##  [671]  0.085545942  0.835445199 -0.011632348  1.413689895  0.857012897
##  [676] -0.719237650  1.290530948  1.227205500  2.266763844  1.316424733
##  [681]  1.682319970  1.519461937  1.428121738  1.999102945  2.502376417
##  [686]  1.354898043  0.864640001  2.292473888  0.692153649  2.221982297
##  [691]  0.855828824  0.785486114  2.288544516  2.168751493  1.287450477
##  [696]  0.333331554  2.235826128  0.778296882  1.966412401  0.884107527
##  [701]  1.359147104  2.199118953  2.424862102  1.340510365  1.305178715
##  [706]  2.059668609  0.783605915  0.673844295  0.649944904  0.425815297
##  [711]  1.354500543  2.222141842  0.789525550  2.331151587  1.262283346
##  [716]  1.845423522  0.112943667 -0.200339322  2.113163436  0.838368271
##  [721]  0.837916703  0.895511259  0.767528210  1.317337304 -0.147303359
##  [726]  0.157934145  0.845829495  0.846704185  1.368258747  0.787422211
##  [731]  1.943789794  1.387518430  1.341832294  2.175013360  1.445673777
##  [736]  2.330506103  1.333963365  1.383741178  2.071543638 -0.474534566
##  [741]  1.486201898  2.173476099  2.187063904  2.217332142  1.344590411
##  [746] -0.814533836  0.831520276  1.308889241  1.446667717  1.439493902
##  [751]  1.260676070  2.560820676  0.698061981  1.168209992  1.821186358
##  [756]  2.227285966  0.825647726  1.189020685  1.124987416  0.511866389
##  [761]  0.824887785  1.303646102  1.318164182  1.403839803  0.773919510
##  [766]  0.883288384  0.755620128  1.379426115  0.396059472  0.396297349
##  [771]  1.949151513  0.816064091  0.470899470  2.202588537  0.945914234
##  [776]  2.244456641  2.462851383  1.230468112  0.804292457  1.937269364
##  [781] -0.056644387 -0.764551969  0.524895206  2.276558603  1.634062453
##  [786]  1.753172841  0.776746031 -0.353347990 -0.113419937  1.870512763
##  [791]  1.099550369  1.912860787  2.180226370  0.823205509  1.142013485
##  [796]  0.928700453 -0.894473417  1.327836939  1.348500063  1.262676184
##  [801]  0.919289110  1.991966924  1.541508862  1.781248880  2.128582705
##  [806]  0.508153831  2.221522749  1.518942785  0.831079954  0.390338942
##  [811]  2.089946588  0.512479471  2.309649710  1.361083614 -0.698179244
##  [816] -1.518429332  2.127977410  0.439870900  0.815020926  1.415084006
##  [821]  1.338822634  0.055317772  0.868913447  1.265531498  1.520854729
##  [826]  0.273453711  0.461555500  0.170920164  2.351516407  1.394673681
##  [831]  0.416999925  0.660668934  2.298282571  0.429948779  1.396526954
##  [836]  1.422063643  0.868691870  2.016530427  2.009252347  1.265239938
##  [841]  1.442520858  0.869335263 -0.278659160  1.376077845 -0.626901156
##  [846]  1.815741882  0.851743780  1.158619507 -0.239304962  2.372210566
##  [851]  1.387827592  0.892182723  0.793695212  0.886587186  1.808214537
##  [856]  1.481682473  1.142524985  0.133828033  2.040363031  0.341082415
##  [861]  1.347917932  1.436647199  1.306263278  2.255085777  2.426685520
##  [866]  2.187132634  1.440970066  2.337205057  0.902340952  0.242090490
##  [871]  1.609363843  1.329163604  1.260676070  1.478248761  2.234189548
##  [876]  2.238022763  0.506929837  1.893451537  1.334201925  0.722137333
##  [881]  2.149411897  2.342805843  2.183778585  1.314723255  1.387107493
##  [886]  0.738368815  0.967576817  1.871481556  1.256683306  1.784614897
##  [891]  2.419149052  2.478432878  2.148331027 -0.389873822 -0.284129588
##  [896]  0.787128687  1.389736626 -0.199749437  1.398576954  0.796347761
##  [901]  1.265446914  1.299469764  0.117584783  1.820481621 -0.319648707
##  [906]  0.401692817  1.559023447  1.562621631  0.122354361  0.709070729
##  [911]  0.779433280  1.368135096  1.980458900  0.465429971  0.361304790
##  [916]  0.655254867  0.090853232  0.417129571  0.162933660  2.173722694
##  [921] -0.144088045 -0.929387956  1.537801345  2.258592125  0.387919167
##  [926]  2.188652726  1.403674044  0.853706396  0.861488236  0.808922468
##  [931] -0.933298056  0.744367269  1.870512763  1.436302372  0.401936989
##  [936]  0.420099608  2.169007562  1.711390609  1.278930660  2.279571997
##  [941]  2.335318817  0.441628885  1.286614146  1.521141579  0.848840578
##  [946]  1.386984286  0.921965825  0.736965748  1.501964602  0.175547861
##  [951]  1.938603704  1.419702824  2.207055302  1.655231817  2.223858210
##  [956]  0.455382112  1.393717550  1.500797569 -1.637170070  1.367187590
##  [961]  1.303115179  2.341433556  0.984942756  1.178708039  2.230319529
##  [966]  1.200242803  0.399824716  1.307066897 -0.272521138  0.876930638
##  [971]  0.866555480  1.530824938  0.129613556 -0.182669365  0.931542669
##  [976]  0.396289722  1.398414367  2.417739793 -0.027794989  2.424000495
##  [981]  1.315591647  2.202895242  2.401898144 -0.528800287  2.058569690
##  [986]  0.807501497  0.644345970  1.335161804  1.439493902  1.556290575
##  [991]  0.177191833  1.564139172  0.458236968 -0.208909581  2.133875092
##  [996]  0.059286984  0.458258505  0.516020112  1.186367615  2.092648944
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
##   0.54295380   0.39880616 
##  (0.12611358) (0.08917388)
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
## [1] -0.7063895 -0.2634284 -1.0683424  0.3274153 -0.2184052  0.3752975
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
## [1] -0.015
```

```r
sd(results2$thetastar)  #the standard deviation of the theta stars is the SE of the statistic (in this case, the mean)
```

```
## [1] 0.920849
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
## t1*      4.5 0.02112112   0.9315686
```

One of the main advantages to the 'boot' package over the 'bootstrap' package is the nicer formatting of the output.

Going back to our original code, lets see how we could reproduce all of these numbers:


```r
table(sample(x,size=length(x),replace=T))
```

```
## 
## 0 1 4 6 8 9 
## 2 2 1 2 1 2
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
## [1] 0.0234
```

```r
se.boot
```

```
## [1] 0.9027998
```

Why do our numbers not agree exactly with those of the boot package? This is because our estimates of bias and standard error are just estimates, and they carry with them their own uncertainties. That is one of the reasons we might bother doing jackknife-after-bootstrap.

The 'boot' package has a LOT of functionality. If we have time, we will come back to some of these more complex functions later in the semester as we cover topics like regression and glm.

