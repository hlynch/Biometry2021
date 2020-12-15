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
## 1 2 3 5 7 9 
## 3 1 2 1 2 1
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
## [1] 0.0526
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
## [1] 2.755296
```

```r
UL.boot
```

```
## [1] 6.349904
```

Method #2: Simply take the quantiles of the bootstrap statistics


```r
quantile(xmeans,c(0.025,0.975))
```

```
##   2.5%  97.5% 
## 2.8975 6.3025
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
##    [1] 6.6 5.2 3.3 4.6 3.4 2.3 5.7 6.0 3.9 5.8 3.4 5.1 4.7 4.5 3.3 6.3 4.5 4.5
##   [19] 5.0 4.2 4.1 3.6 4.5 4.7 5.1 3.8 4.7 4.2 4.7 3.3 4.1 4.9 4.3 5.5 5.8 4.6
##   [37] 4.1 3.9 5.5 5.3 4.7 5.3 6.5 4.2 6.0 4.7 4.7 3.1 3.4 3.7 4.7 4.8 4.4 5.1
##   [55] 4.7 5.9 5.0 4.5 4.3 4.9 2.4 3.4 4.2 3.3 4.3 6.7 2.4 4.0 4.6 3.2 5.3 5.8
##   [73] 4.3 5.7 4.5 4.9 4.1 4.4 4.7 4.9 3.7 3.7 4.4 4.4 5.5 3.5 5.0 4.0 3.3 5.4
##   [91] 3.9 5.3 6.4 4.0 4.8 4.9 4.2 5.4 4.2 5.5 3.5 5.4 5.5 5.0 4.2 5.2 5.0 3.6
##  [109] 5.4 4.0 4.9 3.7 4.3 3.9 4.2 5.2 3.1 5.6 4.3 4.7 4.1 4.7 4.1 6.9 4.3 3.6
##  [127] 4.8 2.7 6.5 5.9 5.1 4.9 2.7 4.0 4.1 4.7 4.8 4.2 5.6 5.5 5.5 3.0 4.3 4.7
##  [145] 4.2 4.4 4.6 3.4 4.4 2.4 5.6 5.3 5.0 4.6 5.7 3.4 5.6 4.7 4.2 4.7 4.3 4.7
##  [163] 4.5 5.7 6.0 4.2 5.2 5.2 5.6 3.8 6.0 4.6 4.9 5.2 4.7 5.5 5.3 4.7 4.0 4.9
##  [181] 5.1 5.6 5.7 6.9 5.6 2.8 4.2 5.0 2.9 5.0 4.6 6.8 3.8 4.5 3.7 5.2 5.4 4.7
##  [199] 5.8 4.9 2.9 5.5 4.8 3.7 4.8 5.9 3.9 6.4 3.1 5.2 4.4 5.9 3.5 4.9 4.3 4.3
##  [217] 5.7 4.4 5.8 4.6 5.1 4.1 5.0 4.4 2.7 4.5 3.7 4.0 5.1 5.1 4.2 5.1 3.2 6.0
##  [235] 4.9 6.3 5.2 3.3 5.2 4.7 4.3 4.5 4.6 4.0 4.0 3.4 6.3 4.6 3.3 4.6 3.7 4.3
##  [253] 3.7 4.9 4.4 3.2 4.0 3.2 3.1 3.7 6.6 5.1 3.9 3.2 5.8 4.8 3.5 4.3 3.1 6.6
##  [271] 4.5 5.4 2.8 4.7 3.7 4.2 4.6 4.0 5.7 5.9 4.5 2.8 5.8 4.2 5.1 4.9 4.5 3.5
##  [289] 4.0 4.3 4.7 3.9 5.2 4.7 4.3 5.3 5.8 6.2 6.3 3.2 4.4 5.5 6.8 3.5 5.1 3.5
##  [307] 4.7 3.7 5.1 5.6 4.5 5.8 5.1 2.9 4.7 3.8 6.6 3.2 3.9 2.3 3.0 5.8 3.2 5.3
##  [325] 4.5 3.1 4.2 4.3 3.6 5.9 3.9 5.5 4.7 5.3 4.0 6.0 3.1 4.0 4.4 4.3 6.0 4.2
##  [343] 5.2 4.8 4.2 4.5 5.7 3.1 4.8 4.7 4.8 5.4 5.4 4.6 4.5 4.8 4.1 6.1 5.3 5.1
##  [361] 3.5 3.5 4.7 3.3 5.3 4.9 3.4 4.5 3.4 3.3 5.5 4.8 3.6 4.1 5.2 5.7 4.7 3.9
##  [379] 4.1 4.3 3.9 4.1 4.3 3.4 4.5 4.8 4.4 4.2 3.1 3.7 3.4 4.4 6.0 3.3 6.4 3.1
##  [397] 4.5 5.1 3.8 5.4 3.0 3.9 4.6 3.8 5.5 4.0 4.9 5.9 3.9 4.0 3.7 3.1 3.8 4.7
##  [415] 4.3 2.9 4.4 3.9 2.4 2.7 4.8 3.3 6.1 4.4 4.8 3.4 4.3 1.7 4.8 4.6 4.4 5.2
##  [433] 5.6 4.3 6.0 3.2 5.5 4.7 5.3 5.1 4.0 4.2 6.9 3.7 4.7 4.3 4.9 2.9 4.2 4.0
##  [451] 4.8 3.0 5.8 3.4 5.0 4.9 4.6 5.0 3.2 4.5 3.7 4.2 3.7 4.2 5.4 3.2 6.0 3.1
##  [469] 5.4 5.1 5.3 4.5 3.3 3.6 4.9 5.8 4.6 5.1 5.2 4.9 4.2 3.4 4.2 5.8 3.7 4.5
##  [487] 4.0 3.5 4.8 4.7 5.2 4.5 4.7 4.8 4.2 4.4 6.1 3.4 4.5 4.1 4.3 4.6 5.4 4.6
##  [505] 5.2 5.8 5.7 4.4 4.0 4.6 5.5 4.1 4.9 3.9 4.6 6.4 4.1 2.5 3.9 5.8 5.8 4.0
##  [523] 5.0 5.8 4.0 4.1 5.6 4.8 3.4 3.5 5.7 5.7 5.2 4.4 5.0 5.1 4.0 4.9 5.2 5.3
##  [541] 5.2 3.9 3.9 5.3 2.6 3.6 3.7 2.8 4.0 3.7 2.7 4.0 4.4 5.3 4.0 5.1 5.6 6.1
##  [559] 5.3 4.9 4.6 4.9 4.5 3.9 5.4 5.7 4.1 3.8 5.5 6.0 4.9 6.1 3.6 5.3 6.5 5.0
##  [577] 5.2 5.3 1.9 2.1 5.0 4.9 3.7 5.3 5.0 4.0 5.0 5.5 5.2 4.5 4.8 4.5 4.9 6.5
##  [595] 6.7 5.2 4.6 3.9 4.3 4.3 6.3 5.2 5.6 4.1 5.7 4.7 4.4 4.5 4.0 3.8 4.9 4.9
##  [613] 5.6 4.7 5.4 6.9 2.6 4.9 4.4 4.1 6.1 4.7 4.3 3.5 5.4 3.9 5.9 3.6 4.1 4.8
##  [631] 3.8 3.2 2.7 4.5 4.5 3.6 5.0 7.3 5.2 4.7 4.8 4.1 5.5 4.8 4.5 4.1 4.1 3.3
##  [649] 3.5 4.9 5.0 4.2 5.5 4.6 4.9 3.7 5.9 4.7 3.7 4.4 5.6 3.4 2.6 4.2 3.9 4.2
##  [667] 6.0 3.0 3.1 4.4 4.1 5.3 4.2 4.8 4.4 4.3 4.7 4.6 4.2 4.6 5.5 4.8 3.9 4.5
##  [685] 4.4 5.9 6.3 3.7 2.8 3.7 4.9 5.5 5.3 4.6 4.4 3.6 6.1 5.8 3.7 4.2 4.1 5.4
##  [703] 4.9 3.2 4.3 4.6 4.7 4.0 5.4 5.2 4.1 4.4 4.6 4.0 4.3 5.2 4.6 5.1 3.1 4.5
##  [721] 3.8 3.1 5.5 4.6 5.0 5.5 3.9 5.1 4.6 5.2 4.8 4.9 4.1 4.1 5.9 4.8 4.4 3.7
##  [739] 5.2 3.5 4.3 2.9 3.0 5.0 5.1 5.3 4.5 3.7 3.9 3.9 3.8 5.1 6.2 4.3 3.3 5.9
##  [757] 4.8 3.8 5.1 5.6 6.3 4.9 5.5 4.5 4.9 3.1 3.6 3.0 5.2 5.7 4.8 3.6 4.2 5.5
##  [775] 4.7 3.5 3.7 3.7 4.0 3.3 5.3 4.0 4.5 1.8 4.3 5.3 6.7 3.9 5.1 4.9 4.5 2.4
##  [793] 3.8 5.0 4.7 4.9 4.7 5.3 5.7 5.3 4.9 5.2 4.8 5.3 4.0 3.9 1.6 5.8 3.5 4.5
##  [811] 5.2 5.7 4.6 4.6 3.3 4.6 5.1 3.3 5.3 4.6 4.4 5.7 3.7 5.6 4.3 4.9 3.9 2.6
##  [829] 2.7 3.7 5.3 4.6 3.5 5.4 5.1 5.0 5.6 4.4 4.4 5.9 4.4 3.7 3.0 4.3 5.1 4.2
##  [847] 6.3 4.9 4.9 5.8 4.9 7.5 5.2 4.3 4.8 5.7 5.8 5.4 3.0 3.3 4.0 5.1 5.2 5.8
##  [865] 4.3 5.3 4.0 5.1 4.9 4.3 3.3 3.9 6.0 4.7 6.2 4.2 6.0 5.6 4.9 4.2 4.5 4.3
##  [883] 4.5 4.9 5.6 5.5 5.3 3.4 3.4 3.6 3.7 5.5 6.0 4.4 3.9 5.3 4.7 5.8 3.0 3.9
##  [901] 4.5 4.5 5.4 4.0 3.5 4.5 5.4 4.2 5.3 3.2 1.8 3.5 5.8 3.1 4.2 5.0 4.5 4.9
##  [919] 4.1 5.1 4.6 5.7 5.7 3.0 5.3 6.2 2.7 5.5 5.2 3.6 5.9 1.6 3.0 3.4 4.1 4.5
##  [937] 5.2 4.9 4.3 3.0 3.8 3.6 5.0 3.5 5.2 2.9 4.7 4.6 5.0 5.1 3.6 6.1 4.9 4.5
##  [955] 4.8 5.4 4.4 5.8 5.7 5.7 4.2 5.0 3.8 3.6 2.9 3.1 3.7 5.0 2.8 6.1 5.3 5.6
##  [973] 4.7 4.4 4.3 3.0 4.5 5.4 5.3 5.2 5.5 4.8 4.8 4.5 4.0 4.3 5.3 4.9 3.4 4.2
##  [991] 4.8 5.2 3.5 4.8 2.9 6.0 4.3 4.5 3.9 4.7
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
##    [1] 5.0 3.6 5.4 4.8 5.1 5.3 5.6 3.7 6.6 4.9 2.9 4.8 4.0 3.9 4.8 3.9 2.8 3.6
##   [19] 3.6 3.9 3.9 4.7 5.0 5.7 4.8 5.4 3.6 3.9 4.6 4.5 4.6 3.5 4.4 5.2 4.6 4.9
##   [37] 6.1 3.1 4.6 4.5 5.7 4.8 5.8 3.8 3.2 4.3 3.4 4.4 4.7 5.5 2.9 5.4 5.3 2.9
##   [55] 4.1 5.7 6.3 3.7 5.5 4.0 3.7 3.4 4.5 3.3 5.3 4.7 3.8 3.6 4.9 5.3 4.8 4.0
##   [73] 5.1 5.0 3.2 5.3 3.0 4.1 5.1 4.0 5.0 3.6 5.8 4.1 3.6 4.6 5.2 3.9 2.4 4.7
##   [91] 3.9 5.2 3.5 4.8 5.9 3.4 5.6 4.6 2.6 6.2 2.9 3.2 5.4 6.2 4.5 5.2 5.4 4.5
##  [109] 4.9 4.8 5.2 4.1 5.1 4.8 4.2 5.3 4.6 3.6 5.5 4.3 5.3 3.8 4.7 5.3 4.0 2.5
##  [127] 4.8 3.3 4.3 3.4 6.4 3.3 4.2 5.2 4.0 3.4 5.8 7.0 5.0 5.6 3.9 4.5 4.3 4.7
##  [145] 3.7 3.8 4.5 4.3 6.2 4.3 4.5 6.7 3.9 3.2 4.8 5.4 3.5 4.1 5.3 2.8 3.8 5.6
##  [163] 4.1 5.1 3.5 6.3 4.4 3.8 5.6 2.9 4.5 4.2 3.7 5.1 2.7 2.3 4.3 5.2 3.1 4.8
##  [181] 4.7 3.7 4.4 5.4 4.9 5.7 4.9 3.7 4.2 3.2 3.7 4.0 4.6 5.6 5.3 4.4 5.5 5.1
##  [199] 6.2 4.7 4.9 5.6 3.4 6.3 4.8 3.2 5.1 5.1 4.7 4.8 5.0 3.8 4.3 4.9 5.3 5.1
##  [217] 5.2 3.1 3.2 3.7 5.8 4.9 4.7 5.7 3.7 4.1 6.0 5.8 3.7 5.7 4.5 4.1 5.2 5.0
##  [235] 3.7 5.3 6.4 5.5 6.6 4.8 2.8 3.9 5.3 2.6 4.2 4.5 4.0 3.4 4.3 4.3 5.4 5.4
##  [253] 3.7 5.8 3.8 4.3 6.7 5.6 4.7 5.5 4.4 5.4 3.2 5.7 4.9 4.7 5.5 4.6 6.1 3.4
##  [271] 4.6 6.1 4.4 4.5 4.8 2.6 4.6 4.9 3.0 3.1 5.8 4.0 5.7 3.8 3.9 7.0 2.6 4.0
##  [289] 4.1 4.5 2.5 3.6 4.0 6.4 5.9 5.0 5.3 2.8 3.0 4.0 4.2 4.0 5.0 2.8 5.0 5.3
##  [307] 4.0 4.7 4.8 6.1 4.4 4.7 5.3 4.6 4.5 3.1 3.7 5.0 5.5 4.7 5.1 5.4 5.2 5.9
##  [325] 5.6 2.8 5.4 3.9 4.2 4.8 4.3 5.0 6.3 5.4 3.7 5.4 6.0 5.3 4.0 3.6 4.9 5.0
##  [343] 3.0 3.4 3.3 5.4 4.6 5.5 3.8 4.7 4.9 4.4 4.1 4.3 3.3 5.0 4.6 6.9 4.5 4.9
##  [361] 4.6 5.2 3.4 4.6 4.1 4.8 4.0 4.5 6.4 4.8 4.1 2.1 5.6 4.5 5.0 4.4 5.9 4.7
##  [379] 3.7 4.2 3.2 4.5 4.8 5.6 3.5 4.8 3.8 4.2 5.1 6.5 3.0 4.0 3.7 5.6 5.3 3.7
##  [397] 3.5 4.0 4.8 4.2 5.3 4.0 3.9 6.3 4.3 4.7 3.5 4.2 6.2 4.1 4.1 4.8 4.5 5.1
##  [415] 3.6 3.9 6.3 4.0 3.5 5.2 4.4 3.0 4.3 4.1 4.5 5.0 5.5 5.7 4.2 4.9 4.6 4.8
##  [433] 3.2 5.9 3.9 4.2 4.1 4.2 5.0 4.4 1.7 3.8 5.9 4.9 4.7 5.2 5.2 5.4 2.7 3.6
##  [451] 4.8 3.6 5.0 4.9 3.2 4.7 5.2 4.9 6.9 5.2 3.6 2.8 4.9 3.9 4.2 4.5 5.6 5.6
##  [469] 4.4 5.1 4.1 4.2 3.4 4.0 4.6 3.1 3.8 3.5 4.6 3.3 4.9 4.4 4.3 5.6 4.3 3.1
##  [487] 4.3 6.2 6.7 4.2 3.6 2.3 5.6 4.1 3.8 4.6 5.2 3.9 2.2 4.0 4.1 5.0 4.3 4.7
##  [505] 4.6 3.1 3.7 4.9 3.8 4.0 3.5 5.1 4.7 4.9 3.6 5.7 4.7 4.7 4.2 4.5 5.3 4.4
##  [523] 4.3 5.2 4.2 4.7 5.8 4.6 4.7 6.2 5.4 4.6 4.9 4.2 4.8 4.0 5.8 5.1 4.6 5.2
##  [541] 4.5 3.8 4.9 4.0 4.1 2.0 5.1 4.1 4.1 3.9 3.9 4.5 5.5 4.4 3.7 4.2 3.3 5.3
##  [559] 4.9 3.9 4.5 3.2 6.1 2.3 5.5 3.9 3.5 5.2 4.8 4.6 5.9 5.2 5.2 3.8 5.6 4.4
##  [577] 3.1 5.3 3.4 3.3 5.4 4.6 5.4 6.3 5.0 4.5 4.6 5.3 3.4 5.0 3.8 6.0 3.8 5.1
##  [595] 6.8 5.3 3.4 2.9 5.0 5.4 2.4 4.1 4.8 2.4 3.9 2.9 4.3 5.2 4.7 5.7 3.1 5.3
##  [613] 5.6 4.6 5.1 5.4 4.8 3.4 4.4 5.5 3.4 3.6 5.4 4.2 5.3 3.8 5.5 5.1 5.0 4.3
##  [631] 6.8 4.4 5.4 4.6 4.1 3.0 4.3 4.2 4.8 4.4 4.6 5.1 3.1 5.5 4.5 4.1 3.8 4.5
##  [649] 5.0 3.2 4.9 3.3 3.9 4.9 3.9 4.2 3.8 4.7 3.2 4.3 4.3 4.0 4.8 4.0 5.0 4.6
##  [667] 4.8 3.3 4.4 6.3 3.2 6.7 4.1 3.5 4.0 4.8 3.9 5.2 4.4 4.1 5.2 2.5 5.6 3.7
##  [685] 5.1 3.5 4.2 4.4 4.2 2.5 4.5 6.1 3.8 4.0 4.0 5.6 5.1 4.2 4.1 3.6 5.4 3.9
##  [703] 3.7 3.2 2.3 3.9 3.0 4.2 5.2 6.2 3.1 4.9 3.6 3.3 3.2 3.9 4.6 4.3 4.9 4.8
##  [721] 3.0 3.8 4.2 4.1 4.3 5.6 4.1 5.6 4.4 5.0 4.6 6.0 4.0 5.5 5.2 4.3 6.7 5.0
##  [739] 3.2 2.5 3.8 3.5 3.2 4.0 3.0 5.2 4.0 4.4 4.4 4.5 3.8 5.0 4.2 5.1 3.3 5.0
##  [757] 4.6 5.4 4.2 3.1 4.9 6.3 4.9 4.5 2.7 3.6 3.4 3.7 3.9 5.0 4.3 5.8 4.6 3.0
##  [775] 3.6 4.7 3.9 5.0 5.1 5.4 5.9 3.4 6.9 4.6 5.5 4.7 4.4 5.1 4.4 4.4 2.4 4.9
##  [793] 5.4 2.4 4.8 2.9 3.6 4.5 6.8 5.3 4.9 3.5 3.2 4.4 4.5 5.8 5.2 6.0 5.8 4.9
##  [811] 2.9 5.3 5.6 4.6 4.4 4.4 4.4 4.9 4.2 5.4 4.1 5.1 4.0 3.9 5.1 3.5 4.5 4.1
##  [829] 5.1 4.4 3.9 3.2 4.9 5.8 5.2 4.7 3.5 4.9 4.8 3.6 5.3 3.4 2.9 4.4 5.5 4.5
##  [847] 4.5 3.5 4.8 5.3 4.3 4.3 4.0 5.9 5.1 4.6 5.3 4.8 3.9 4.9 3.5 5.2 5.6 3.8
##  [865] 4.3 2.0 6.2 6.1 4.4 3.5 3.7 3.7 4.1 4.0 4.4 6.8 4.0 4.8 4.7 5.8 2.9 5.7
##  [883] 4.2 3.7 5.5 4.1 5.8 2.6 5.9 4.2 4.6 5.0 5.1 3.4 4.6 4.0 4.1 4.7 5.4 3.6
##  [901] 3.4 6.6 5.3 2.1 4.1 5.6 4.5 5.4 3.9 3.7 5.8 6.2 5.0 6.2 3.0 3.7 3.4 5.2
##  [919] 4.8 4.2 4.0 5.4 3.4 3.8 7.1 6.3 3.6 3.9 4.2 4.8 3.9 5.3 3.8 3.8 5.8 5.0
##  [937] 5.9 3.5 4.6 2.9 5.9 5.3 5.5 3.8 4.0 5.2 4.6 4.2 3.6 2.8 4.2 4.5 4.8 3.1
##  [955] 4.0 2.8 3.9 5.6 3.4 3.4 5.4 3.6 4.3 4.4 4.7 5.2 3.3 5.2 4.9 4.2 5.8 6.7
##  [973] 4.5 4.6 4.8 5.7 3.5 4.9 2.7 4.0 2.9 5.6 2.5 3.1 3.9 5.8 4.5 5.1 4.8 4.1
##  [991] 5.8 3.7 4.7 4.6 5.1 4.1 3.3 4.5 5.5 4.2
## 
## $func.thetastar
## [1] -0.0064
## 
## $jack.boot.val
##  [1]  0.5491620112  0.3938775510  0.3584070796  0.1843373494  0.0005952381
##  [6] -0.0943620178 -0.1529891304 -0.3038356164 -0.4005882353 -0.5780281690
## 
## $jack.boot.se
## [1] 1.055633
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
##    [1] 4.7 4.9 2.9 4.1 3.4 2.7 5.0 3.7 3.9 3.9 4.8 4.4 4.3 3.5 3.4 4.9 4.3 5.4
##   [19] 4.8 5.7 4.2 2.9 6.7 5.0 6.1 4.6 5.3 4.3 4.8 2.9 5.5 5.1 3.5 3.2 3.7 5.7
##   [37] 5.1 3.8 3.6 5.4 5.6 2.6 4.8 5.8 6.2 4.7 3.9 4.8 4.1 5.8 3.2 3.6 4.3 3.6
##   [55] 4.9 4.7 4.5 4.1 3.0 2.2 5.3 3.6 4.0 4.8 4.9 5.1 4.4 3.4 3.3 3.5 5.8 3.3
##   [73] 3.4 6.2 4.5 4.8 6.0 5.2 4.0 4.0 4.7 4.9 4.9 3.7 4.6 4.3 4.6 4.3 4.6 5.0
##   [91] 4.7 4.0 3.2 6.0 3.5 4.9 1.8 4.7 4.4 5.1 3.5 4.7 3.7 5.6 5.3 5.9 3.7 4.5
##  [109] 2.8 5.3 4.4 5.0 7.3 3.0 4.0 5.8 4.0 4.7 5.4 4.6 4.3 4.7 4.4 5.4 3.9 5.9
##  [127] 4.4 5.3 3.5 4.2 3.6 4.1 4.4 4.5 3.9 4.5 4.8 4.7 5.7 4.3 5.3 3.7 6.1 3.4
##  [145] 5.0 6.6 3.0 4.5 4.0 5.3 4.3 4.0 3.6 2.9 5.5 5.2 2.6 4.5 4.5 3.4 5.2 4.2
##  [163] 3.9 4.0 4.4 5.1 3.7 3.1 5.4 5.3 5.0 3.8 6.2 3.5 4.4 5.1 4.4 4.4 4.6 5.2
##  [181] 5.1 5.5 4.2 3.4 5.2 4.8 4.3 4.3 6.2 3.1 2.5 6.8 5.6 4.0 3.5 4.2 3.1 4.0
##  [199] 4.3 5.4 4.6 4.6 4.3 3.2 3.6 4.7 5.0 5.0 5.2 2.6 4.2 4.8 4.1 3.8 4.4 3.5
##  [217] 3.4 4.5 4.9 2.9 5.3 4.5 4.1 3.8 5.3 5.4 5.6 4.4 5.5 3.3 4.6 5.4 3.1 5.5
##  [235] 4.7 2.6 3.7 5.2 4.2 3.3 4.6 3.6 3.7 4.5 4.9 6.1 5.0 3.4 6.4 4.9 3.7 5.2
##  [253] 5.2 5.3 3.7 4.2 4.1 3.3 5.7 2.9 5.5 4.5 3.1 3.1 3.8 4.2 3.9 4.1 4.0 5.0
##  [271] 4.5 5.8 3.1 4.5 5.9 4.4 5.4 5.2 4.7 5.5 5.0 6.2 4.6 4.7 4.4 3.3 3.6 5.3
##  [289] 4.9 6.0 3.9 4.3 4.6 4.9 4.7 5.3 3.5 5.1 5.3 6.1 4.0 3.6 4.1 4.8 6.6 4.1
##  [307] 5.3 5.4 2.9 4.4 3.9 4.8 4.3 3.6 4.8 4.0 4.2 4.2 3.1 4.8 5.5 5.7 3.5 6.0
##  [325] 4.8 5.2 4.0 7.4 3.7 3.4 4.6 4.2 3.4 4.6 5.5 3.6 3.5 5.2 5.8 6.0 2.7 2.5
##  [343] 5.8 3.2 3.4 4.8 4.3 6.0 5.1 5.4 4.0 3.8 4.2 4.9 4.7 4.1 4.2 6.2 4.4 4.5
##  [361] 4.5 5.0 4.8 4.5 5.4 6.0 5.6 5.6 2.8 3.3 5.0 4.6 4.2 4.0 2.6 6.0 4.4 3.5
##  [379] 3.6 5.6 3.3 4.0 3.8 4.0 4.4 4.7 3.4 5.5 4.4 4.3 4.1 3.6 4.4 4.3 2.9 6.0
##  [397] 3.4 4.6 3.1 3.8 4.6 3.7 5.6 4.0 3.8 4.4 6.1 4.1 3.9 3.4 5.9 3.4 3.2 5.3
##  [415] 4.3 5.0 5.5 3.0 5.9 2.2 5.0 4.1 4.3 4.0 4.0 4.8 3.4 5.0 3.6 4.3 4.3 3.8
##  [433] 5.4 2.2 5.3 4.4 4.7 3.6 4.7 4.8 2.5 4.3 3.4 5.3 5.0 4.0 3.8 4.1 1.9 3.7
##  [451] 5.1 3.2 4.3 3.7 6.1 4.2 4.8 6.0 3.5 4.1 1.7 3.1 4.2 3.9 5.0 4.7 2.9 3.9
##  [469] 3.1 5.2 4.0 4.4 3.7 5.1 5.9 3.2 4.1 3.4 4.0 3.2 5.0 6.6 5.4 4.1 3.5 4.3
##  [487] 4.8 3.0 5.4 2.9 4.1 4.6 5.9 5.2 3.7 3.9 3.6 4.7 5.5 3.6 5.4 5.2 4.4 3.5
##  [505] 3.6 4.0 5.1 4.2 4.3 4.8 4.1 5.0 3.6 4.8 3.7 5.8 5.0 5.0 5.5 3.2 4.7 3.8
##  [523] 3.2 3.8 4.8 3.2 4.3 5.0 5.1 5.8 3.3 4.2 5.8 3.0 5.6 4.5 2.8 3.2 5.6 4.0
##  [541] 6.1 5.1 3.7 2.9 4.3 4.4 4.0 5.6 4.4 4.7 4.6 2.6 5.2 5.6 4.2 5.5 4.0 5.9
##  [559] 3.3 5.1 3.5 1.6 5.8 4.6 4.5 4.8 2.7 4.1 6.3 5.8 2.9 3.2 5.2 5.6 5.1 4.0
##  [577] 4.7 4.4 4.3 4.0 6.2 4.0 3.4 3.0 5.0 4.7 5.4 5.9 4.6 4.8 5.2 4.5 4.5 4.0
##  [595] 3.8 3.3 3.4 4.3 5.0 4.7 5.2 4.4 5.7 5.7 5.5 4.4 3.1 5.3 4.9 4.9 4.9 4.3
##  [613] 4.8 5.1 3.2 4.6 3.8 4.2 5.4 5.7 4.6 3.2 3.4 5.4 3.0 3.5 5.2 4.6 5.6 4.1
##  [631] 4.9 6.3 4.2 4.0 4.2 4.6 4.9 3.4 5.1 4.4 3.2 4.5 5.8 3.8 4.9 3.9 4.2 3.8
##  [649] 2.9 5.1 5.4 4.8 2.8 5.1 3.9 5.6 4.3 4.4 4.0 6.8 3.9 6.0 5.6 3.3 4.1 5.5
##  [667] 4.7 5.1 4.2 5.2 5.0 3.5 4.4 4.3 6.0 7.0 4.8 4.5 4.0 3.2 4.7 3.4 4.1 4.1
##  [685] 4.6 4.0 3.9 5.0 4.4 4.0 5.3 5.0 4.0 5.3 3.8 6.6 3.9 4.1 3.7 2.9 4.0 4.0
##  [703] 4.9 2.9 3.1 4.2 3.9 4.1 3.9 4.5 4.7 4.7 5.1 3.3 4.4 5.3 4.2 4.5 4.7 3.1
##  [721] 4.9 5.3 6.1 4.9 5.5 4.0 4.5 3.8 5.9 4.0 4.6 3.2 5.4 5.5 4.0 4.3 5.7 4.5
##  [739] 4.3 3.8 5.3 5.0 5.4 3.1 4.3 5.0 4.6 5.1 4.3 4.9 5.0 4.6 4.5 4.2 4.2 4.2
##  [757] 6.3 2.9 4.6 4.1 4.1 4.9 4.9 3.5 5.3 4.3 5.7 5.4 4.2 3.2 3.9 6.6 5.8 3.4
##  [775] 3.8 3.9 4.8 3.3 4.4 5.1 6.5 4.6 5.7 4.5 4.9 5.7 4.2 4.4 3.7 3.3 3.3 4.5
##  [793] 4.2 4.9 3.0 3.6 4.5 3.6 4.7 2.9 4.2 3.3 3.7 4.2 5.0 5.3 5.7 5.0 2.8 6.0
##  [811] 4.2 3.6 4.3 3.7 4.0 4.8 4.4 4.3 4.4 5.9 4.7 3.6 4.1 4.7 4.8 5.2 5.7 3.7
##  [829] 5.3 4.4 5.1 5.8 4.0 5.1 4.4 3.3 4.0 4.3 3.3 2.6 6.6 5.2 3.9 4.8 5.1 4.7
##  [847] 3.3 4.4 4.7 4.5 5.1 5.1 4.2 5.8 5.1 3.7 4.5 3.7 4.6 5.5 5.1 4.6 3.8 3.8
##  [865] 4.4 3.9 5.1 3.9 3.6 3.9 5.4 3.5 5.2 4.1 3.2 4.1 4.3 5.9 4.9 3.4 4.4 5.1
##  [883] 5.2 4.9 4.6 4.6 4.9 5.3 4.7 5.1 5.8 6.2 4.3 3.4 5.0 4.7 5.2 6.5 5.3 4.1
##  [901] 3.1 4.9 3.8 4.7 4.8 4.2 4.3 5.7 2.8 3.8 4.3 4.3 3.9 4.0 3.1 5.3 3.8 4.6
##  [919] 4.7 5.1 4.7 5.1 5.7 3.3 2.5 3.6 4.1 2.5 5.7 5.3 5.7 5.3 5.6 4.6 5.6 4.4
##  [937] 5.3 4.7 4.9 4.4 5.2 3.7 3.5 3.7 4.1 4.8 5.0 5.0 3.9 3.9 4.1 6.0 4.6 5.6
##  [955] 3.5 4.5 4.5 3.9 4.1 5.8 4.3 4.5 6.2 4.1 4.1 3.6 4.1 4.9 6.0 2.5 4.4 6.1
##  [973] 3.0 4.6 3.7 4.5 4.2 4.5 3.5 5.3 4.1 4.9 4.3 4.2 4.1 3.9 4.7 3.3 4.9 3.7
##  [991] 4.1 5.2 5.3 4.6 3.3 3.5 4.2 4.5 4.4 3.4
## 
## $func.thetastar
## 72% 
##   5 
## 
## $jack.boot.val
##  [1] 5.500 5.308 5.200 5.100 5.000 5.100 4.800 4.700 4.600 4.400
## 
## $jack.boot.se
## [1] 0.9704143
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
## [1] -0.8332946
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
##      shape       rate   
##    6.653137   12.415205 
##  ( 2.903834) ( 5.628466)
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
## [1] -0.06213694 -0.15933199  1.22009125  1.22336798 -0.26079384  0.53783690
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
##    [1] -0.4243650014 -1.4518611925 -1.0031734548 -1.2926478662 -0.4900763739
##    [6] -0.4542992189  0.9553811946 -1.0548726732 -0.4638488215 -0.7918012508
##   [11] -1.6469702884 -1.0584782036 -0.4507486093  0.1538163082 -0.1180211898
##   [16] -1.3423282000 -0.0369443500 -0.3138139585 -0.5137697477 -0.5918310634
##   [21] -0.0200908023 -0.8881330343 -0.6303249788 -0.5403067868 -0.8053663586
##   [26] -0.5375785477  0.0703243779 -1.3455827965 -1.0924513725 -0.6285531558
##   [31] -0.5502538565  0.4579503463 -0.0278232908 -0.6866490839 -0.5049967369
##   [36] -0.8671376674 -0.9125995554 -0.2297751829 -1.2182205031 -0.3556861956
##   [41] -0.0665855204 -0.6647981630  0.0928686931 -0.9937063227 -0.3904852211
##   [46]  0.9671478522  0.3456292017  0.1693642689 -0.0785719066 -0.6618329387
##   [51] -1.5755387986 -0.8623305895 -0.6426735919 -0.6169861591 -0.5479126818
##   [56] -1.3606207889 -0.6227914751 -0.5681838922 -0.1509703086 -0.1674699065
##   [61] -0.1791899823 -0.4803574255 -1.3675819689 -0.9125995554 -1.4492430142
##   [66] -0.0562382790 -0.2595550455 -0.3278616060 -0.7309440945 -0.6732571479
##   [71] -1.2286700370 -0.8277061388 -0.1859683199 -0.6349611433  0.1200037539
##   [76] -0.8815451709 -0.9712910859  0.1275268605 -0.6474087879 -0.9017287709
##   [81] -0.2029995939  0.0593998424  0.3644347456 -1.0044153638 -0.6939607665
##   [86] -1.2858079479  0.6234878826 -0.7979590047 -0.4601553212 -0.9061263917
##   [91]  0.1597763613  0.3704008382 -0.5558617898 -1.2349039165 -1.1151663544
##   [96] -1.1873457518  0.8869858005  0.1422793901 -1.6683611987 -0.6979142953
##  [101] -0.5039542661 -1.8196443409 -0.3173871623 -0.6710527755 -0.5258254055
##  [106] -0.3181309028 -1.3709222356 -1.1765172877 -1.2093698788 -1.3935067251
##  [111] -0.6720304970 -0.7175040013 -1.0684575334 -0.9005310878 -0.6164047967
##  [116]  0.1107524002 -1.1190281693 -1.2031539337 -0.5489991771 -0.7175278590
##  [121]  0.1905420472 -0.9995174539 -0.7903668029 -0.7694654710  0.5633914789
##  [126]  0.3817484672 -1.0425870818 -0.3370388330 -0.1339136364 -1.3665005148
##  [131] -1.0553226038 -0.3931629663 -0.7805022106 -0.1786370818 -0.6363122464
##  [136]  0.1971472336 -0.6727204783 -1.1306799792 -0.4419855048  0.1338479511
##  [141] -0.5113143463 -0.6017057631 -0.3039392798 -1.2083937953 -1.4835762841
##  [146]  0.3851209561 -0.5414141280 -1.1873457518 -1.7757385459  0.1023987849
##  [151]  0.2109789611 -1.2433676794  0.5186100440 -0.7629684617 -0.9428697518
##  [156] -1.0030712875 -0.0234494128 -1.0813201846 -0.2355564668 -0.6789977639
##  [161] -0.4298776419  0.2758096410  0.3907136492 -0.9931692548 -0.7744891352
##  [166] -1.7309069996 -0.8264460718 -0.9421328419 -0.8500403375 -0.4543344101
##  [171] -0.4607494536 -0.4733195991 -0.6322090907  0.0291233595 -0.4270854419
##  [176] -1.0689768898  0.2131948569 -1.2979508914  0.2551888225 -0.0795035951
##  [181] -0.2532494233 -0.3087037192 -0.5288104173 -0.9001142986 -0.8851239333
##  [186]  0.0294910026 -0.7673687264 -0.2211363815 -1.1957038384 -0.8379118413
##  [191] -1.0485938444  0.4836956408 -0.5278654946 -0.4756795714 -1.1552223414
##  [196] -1.0820876444 -0.5357536416 -0.4676138617 -0.9608667096  0.1239432431
##  [201]  0.2444866231 -0.9942897274 -1.1401420137 -0.3622799087  0.2273471669
##  [206] -0.6163194999 -0.6225621064 -0.7573073177 -1.4083573486 -1.3292925098
##  [211] -0.4947196809 -0.1639615377 -0.3464225532  0.1938279423 -0.4565467802
##  [216]  0.3223568006 -0.3168358866 -1.6342947172 -0.2714817495 -0.5387315536
##  [221] -1.0749961765 -0.9190910868  0.4663123815 -1.3325239959 -1.4473306125
##  [226] -0.5715647259 -0.3773541051 -0.6332663779 -0.2724381748 -0.9544901327
##  [231] -1.2592532386  0.0931579470 -0.0416131351  0.5042734079 -0.1157650547
##  [236]  0.5123752959 -0.3976977044 -1.1370039508 -0.5812108442 -0.0365836622
##  [241] -1.6154238524 -0.4374481006 -1.1079907421 -0.3673975616 -0.4773047362
##  [246] -0.8550730038 -0.5270061898  0.0572047115 -1.2396401310 -0.5145872713
##  [251] -0.5758329207 -1.1153074444 -0.5559781885 -0.4717983819 -0.7138757568
##  [256]  0.3402923236  0.2618503875 -1.0888535365 -0.2922795852 -0.7507245408
##  [261]  0.1061197529 -0.9296145904 -0.4677596044 -0.8249144923  0.0427928550
##  [266] -0.5293886194 -0.0950360875 -0.4995402122 -1.4451648820  0.2387940055
##  [271]  0.8823960752  0.0027492937 -1.1678075539 -0.4858655377 -0.1081796746
##  [276] -0.9454432219 -0.3202529579 -0.4237686849  0.7817774125 -0.9535296430
##  [281] -0.0818546400 -0.8363948946 -1.1148134720 -0.7950936311 -1.4126318564
##  [286] -0.6382334321 -0.7486391813 -0.5460146926 -0.6096984853 -0.0674581589
##  [291] -1.2061792712 -0.4011359963 -0.2364966244 -0.7565279940 -0.8223082470
##  [296]  0.3015428833 -0.9710783307 -0.3681489427  1.3986862075 -1.0752878992
##  [301] -0.6052863083 -0.8088322517  0.0590418887 -0.3749876109 -1.3188784010
##  [306] -0.6445275912 -0.3649394247 -0.4078675458 -0.6925339534 -0.1397911889
##  [311] -0.7540418203 -0.0183212449 -1.2414058133 -0.8888365772 -0.9130297403
##  [316] -0.3043280714 -0.1510148970 -0.4076721508 -1.0250680654 -0.4684976763
##  [321] -1.2599939204 -0.3295456415 -1.0457759107 -1.1412024769 -1.1558193368
##  [326] -0.0863486002 -0.6194340529 -1.5507383823 -0.1173352785  0.1234617211
##  [331] -0.0329908573 -0.5536611540 -1.2020638091  0.3082843317 -0.5977008721
##  [336] -1.2720437368 -0.6949680343  0.1354975530 -1.1496510901 -0.5546329094
##  [341] -0.4720509110 -1.0112082389 -1.0472648650 -0.4204676568 -0.4182919487
##  [346] -1.2345089316  0.5211928573 -0.9511389919 -0.0145497306 -0.2491253423
##  [351] -0.5450084163 -0.1394303336 -0.1724967357 -0.3409925731 -0.1830631392
##  [356]  0.2677021432  0.0641104968 -0.6616881651 -0.9441337805  0.5327236055
##  [361] -0.2891890258  0.3143531701 -1.4933035850 -0.5461180480 -0.9649578983
##  [366] -0.5619955414 -0.4155620027  0.1390085213 -0.8053692984  0.6500925127
##  [371] -0.7926868750 -1.5139227269 -0.7618968779 -0.2072375237  0.3590564803
##  [376] -0.6960581468 -0.4986625885 -0.4804385178  0.6682184844 -0.3706252165
##  [381] -0.4256778619 -0.7348896187 -1.0428432317 -0.8542714240 -1.2073716455
##  [386] -0.6114763261 -0.4947735313 -1.0924608128 -0.5659528626 -0.4824443829
##  [391] -0.9736051848 -0.9411945426 -0.2284661808 -0.4934874496 -0.1842139992
##  [396] -0.5841726709  0.8115280869  0.2476061075 -1.1755063582 -0.9596996695
##  [401]  0.3611866111 -0.4260118202 -1.2620389739 -1.1295649752 -0.4989098760
##  [406]  0.5882624965 -0.4933099391 -0.9541722033  0.1355628543 -0.7948254463
##  [411] -1.0417462792 -0.4607494536  0.3728478301 -0.9602205567  1.0425819637
##  [416] -0.8961027730 -0.1135490333 -0.8122994205 -1.4504379999 -1.0908449502
##  [421] -0.8556181083  0.0576654062 -1.4551361299 -0.5842966036 -0.9753948320
##  [426] -1.1623925366 -0.4981958941  0.2182152178 -1.8656935917 -0.6285531558
##  [431] -1.4777869163 -0.6859117209 -0.6058200544 -0.1779301201 -0.8725284462
##  [436] -1.5735464066 -0.6212512541 -1.3201322731 -1.5039097172 -0.9892647815
##  [441] -0.6542882276 -0.8207367522 -0.2660716453 -1.0846174896 -0.3149262261
##  [446] -0.9436700374  0.2348337045  0.3849293860  0.4473502702 -0.8462416588
##  [451] -0.8283237494 -1.0607361191 -1.0572382688 -1.5030424667 -0.8332946346
##  [456] -0.4961998579 -0.8107072171 -0.8338222391 -0.6702904969 -1.0036640522
##  [461] -0.5601425141 -0.7415304773  0.1378795493 -0.7549690977 -0.9531035393
##  [466] -1.0828804848 -0.9918868515 -0.8090473396 -0.3733173677 -1.4217931897
##  [471] -1.0218742921 -0.2148841093 -0.8841734774 -1.5329974391 -0.9636564860
##  [476] -0.7755702716  0.3731938570 -0.3468797979 -0.8952078997 -0.3431079581
##  [481] -1.0919489329 -1.0818400018 -0.2425030767 -0.5579013463 -1.2863750447
##  [486]  0.3130000178 -0.8733346564 -0.4084181796  0.0121347500 -2.2410870156
##  [491] -0.2133852638 -0.6018246341 -0.4209562520  1.2050795419 -0.3844937605
##  [496] -0.7136374627  0.1435474343 -1.0670472969 -0.6186048153 -0.8956071500
##  [501] -1.1916467136 -1.6849944165 -0.6549695141 -1.4632533051  0.1477733944
##  [506] -0.2603361034 -0.5137545974  0.2448341858 -1.5356896362  0.2352833876
##  [511] -0.5223867006 -0.2390494482 -0.3898137654  0.1963813029 -0.6840504890
##  [516] -0.7410181943 -0.1574230930 -0.0913165784 -1.1827188965 -0.5167295185
##  [521] -0.9431996255  0.2024223214 -0.7731579032  0.5032356175 -0.5968765045
##  [526] -0.7189973745  0.1017290845 -0.1909126823  0.4879416891  1.1510973744
##  [531] -0.6514645894 -1.2504425471 -1.0239645552 -1.0498372635  0.1073631460
##  [536] -0.6470129704 -0.5322920109 -0.5807275646 -0.9828727050 -0.8726132055
##  [541] -0.3623539901 -0.6857571747 -0.8708941777 -1.4042247193 -1.4310573956
##  [546]  0.2202006580 -1.0888535365 -0.5167023529  1.1951169381 -1.1707420667
##  [551] -0.2096060165 -0.9084845738 -0.5627793604  0.2582607958 -0.9433909392
##  [556]  0.2077290495 -0.0738442270 -0.0653552043 -0.9071342067 -0.5979613128
##  [561] -0.7480317397 -0.7219619753 -0.5384094965 -1.5665986260  0.3017733603
##  [566] -0.9485975878  0.2143891717  0.0957383143 -0.7069294257 -1.0706982535
##  [571] -0.0923573247  0.8275090173 -0.5367336538 -0.8233855772 -0.3564084418
##  [576] -1.1318807120 -0.9887225253 -0.2458311206  0.1363493262 -0.7405909843
##  [581]  0.3171956449 -1.1732161685 -1.7883165673 -1.4970989730 -0.1777445963
##  [586] -0.1209687671 -1.1786526787 -1.0492088498 -0.9687955559 -0.8537046682
##  [591] -0.0830295616 -1.5670040061 -0.8231878078 -1.3018679419 -0.0627286208
##  [596]  0.0266310419 -0.8305690671 -1.1106620681 -1.4415025982 -1.4446800736
##  [601]  1.4218750790 -0.8682045946 -0.6573368638 -0.9080209918  0.2023653729
##  [606] -0.2021482597  0.1704674511  0.0883632883 -1.7254144405  0.0290003725
##  [611]  0.2299617125  0.1518950457 -0.6279973755  0.3274783710  0.7027618172
##  [616] -1.8896899110 -0.4228248026 -1.2501468926 -0.8158236515 -0.7860658470
##  [621] -0.0640097267  0.1275868418 -1.1230210085  0.0108122036  0.0320248745
##  [626] -0.4353264677 -0.8966976141 -0.5288104173 -0.3012459590 -1.2031830074
##  [631] -0.1242721765 -0.9992240145 -1.0931148435  0.6098917974  0.5986387732
##  [636] -0.6647513811 -0.7505279098  0.1211790558  0.0438927450 -0.1063169270
##  [641] -1.1287147790 -0.2311732060 -0.2956394538 -1.2710327194  0.4181392559
##  [646]  0.2585483323 -0.3822123953 -0.5725455987 -0.3928213536 -0.3813470099
##  [651] -1.0151656723 -0.4474336511 -1.1452456655 -0.3857864091 -0.6665358144
##  [656] -1.2844821623 -0.3790068258 -0.8264460718 -0.0398724627 -0.7809368879
##  [661] -0.7863772961 -0.9639626909 -1.2899328810 -1.0177382777 -0.8850281892
##  [666] -0.9706125904 -0.1892508918 -0.3338413329 -0.9658613330  0.1625029225
##  [671] -1.2310221783 -0.2868343464 -0.2591622122 -0.6483713923 -1.0684575334
##  [676] -0.9480908039 -0.6831795654 -1.0580437919  0.7381361575 -0.6559777190
##  [681] -0.6781920480 -0.1289202004  0.3046977582 -0.7941589835  0.1277509473
##  [686]  0.4432430548  0.5472889656 -0.0175686434 -0.0612721025 -0.2894911252
##  [691] -0.3717758368 -1.7242421017 -0.0773943765 -1.4407829210 -0.4989098760
##  [696]  0.1905420472 -0.0255956977 -0.2007332731 -0.8423889626 -0.8638842382
##  [701] -0.1260035424 -0.1695778666 -0.8067699143 -1.2131667231  0.1795625688
##  [706] -0.6974321268 -1.4586114883  0.2338096146  0.1595914409 -0.8719258318
##  [711] -1.3072570184 -0.3564854054 -0.0938979010 -0.0868363795 -1.0539283567
##  [716] -0.8379715497  0.1307777762  0.4774750845 -0.1265287098 -0.7227286971
##  [721] -0.7900102512 -1.2800919451  0.1118698902 -0.4298179588 -0.2928958169
##  [726] -0.7674928615 -0.2955970692  0.5083794690 -1.2930680043 -1.3892502006
##  [731] -1.3016123896 -0.3512004711  0.5001102565 -0.2809104781 -0.6959377835
##  [736] -0.5755331847 -1.0014168011 -1.2303913787 -1.4989456378 -1.5511215114
##  [741] -0.2418159555 -1.3707563663  0.5156829069 -0.7881927541  0.5826591084
##  [746] -1.0594626839 -0.9897213836 -1.1948356489  0.2752005593 -0.8783867582
##  [751] -0.3076788370  0.4360849280  0.0529646740  0.2662984633 -0.4981275988
##  [756] -0.6142004457 -1.0302073796  0.3057731290  0.6073881866 -0.8022315031
##  [761] -0.2826700661 -0.1448737833  0.1072555268 -1.1016284831 -0.5910915121
##  [766] -1.2461054345 -0.3956363339  0.3166988144 -0.3355898783  0.2601646795
##  [771] -0.2533589503 -0.8973228835 -1.1011651067  0.1504632208 -1.4406016659
##  [776] -0.6070066149 -0.5513110607 -0.5237569968  0.0794725470  0.1854394021
##  [781] -0.3291938554 -0.3630740223 -0.2273180260 -0.1826445172 -0.7058197950
##  [786] -0.6346414572 -1.6617309859 -0.7609774775 -1.3022883976 -0.1456061815
##  [791] -0.1684646846 -0.8312539778 -0.5725836599 -0.6722075082 -1.0377356333
##  [796] -0.1953936740 -0.3693898167  0.1061501459 -0.8754168257  0.5469701180
##  [801]  0.3360147493 -1.2476149969 -0.3000522934 -0.7226536884 -0.4357653899
##  [806] -0.1153576798 -0.8569075446  0.6564724172  0.3682898772 -1.1103119019
##  [811] -1.2618531949 -0.7865278695  0.1820185982 -1.6060717132 -0.4156368271
##  [816] -0.9271523758 -0.5053470524 -1.2052611036 -0.3419320833  0.2666195643
##  [821]  0.0003379775  0.2442159428  0.2228302587 -1.1656249014 -1.0528209298
##  [826] -0.5913987468 -0.1859435975 -0.9968299389 -1.3821479460 -0.1953936740
##  [831]  0.0263222246 -0.8592958219 -0.2224981266  0.0201127373 -0.6380105541
##  [836] -0.6640207592 -0.4995552272 -0.9913871912 -0.7256839387 -0.8469402695
##  [841] -1.3974643904 -1.1496970455  0.1625029225 -1.0050932295 -0.7130350189
##  [846]  0.1986885357 -1.4886067354 -0.6636982438 -1.2473263902 -0.4281591344
##  [851] -0.9587235708 -0.6951077766 -0.7315566581 -0.6127271739 -1.2550907974
##  [856] -0.3174501701 -1.6509022564 -0.9951264626  0.2866608773 -0.7610230689
##  [861] -0.5546498858 -0.3528419930 -0.7959510588 -0.2711030332  0.2273471669
##  [866] -0.6144717940 -1.2672471266 -0.3855286382 -0.9587485025  0.9276927508
##  [871]  0.0235917776 -0.4110112663 -0.9418823942 -0.3028139436 -0.6237008680
##  [876]  0.5133378142 -0.2386624412 -0.9906781992  0.3978331326 -0.9642020035
##  [881] -0.0683220252 -0.7365351049 -0.0541776563 -0.1824109682 -1.8079236450
##  [886] -0.2827217117 -0.3651926149 -0.5160357753 -0.9519922085 -0.9505001578
##  [891]  1.0708290614 -0.7835972484 -0.8426295121 -1.5756560966 -0.8192527837
##  [896]  0.9449622648 -0.7800876362  0.3129167433 -0.8183159725 -0.5459573012
##  [901] -0.8282520426  0.2416566305 -0.9382666490 -1.2196627440 -1.2289144247
##  [906]  0.1888582619 -0.8332946346 -0.3854441291 -0.5016617260 -1.2396631033
##  [911] -0.4780794158 -0.4721982406 -1.7583510349 -0.6126894748 -0.8094865143
##  [916] -1.1907418621 -0.1659816354 -0.6774825817 -1.6133716995 -0.8978237023
##  [921] -0.4633648317  0.3969553675 -1.5858491076 -1.1351324955 -1.2060997177
##  [926] -1.2086183660 -0.5292737308 -0.6859723486 -0.6032943801  0.3726292833
##  [931]  0.2643338999 -0.9182218728 -1.5926796072 -1.0986298682 -0.5710344431
##  [936] -0.7835177999 -0.2207735373 -0.3417578986 -1.0772680722 -0.6402490020
##  [941]  0.4129713293 -0.3223845101  0.4238631423 -0.9854082256 -0.6484276481
##  [946]  0.1566629575 -0.2012896949 -0.9623738326 -1.1899266742 -0.6364072728
##  [951] -1.0819893365  0.1999694729  0.6818045820  0.0957051247 -1.5171288416
##  [956] -0.5341945017 -0.1021608183  0.5419275704 -0.7464210040 -0.9023248777
##  [961] -1.7823624440 -0.1509947272 -0.5341790361 -0.1103996373 -1.4220847346
##  [966] -1.4006017187 -0.5750219204 -0.8347444519 -0.1699912788 -0.6478780782
##  [971] -0.6478983264 -0.8614122973 -1.2728278373 -0.6514645894 -0.3933901033
##  [976]  0.2016490623 -0.2600371233 -0.7871772447 -1.1836506228 -0.9856744928
##  [981] -0.7237514878 -0.6377058675 -1.0928671914 -0.6189159618 -0.6656983939
##  [986] -1.9816865228 -0.1260903346 -1.5838372886  0.5346782475 -0.5461401309
##  [991] -0.4907329739 -0.7357815534 -1.5433249748 -0.3257953951 -0.7837598912
##  [996] -0.3067148905 -2.2286549087 -0.6185461322  1.4053414887  0.0076972096
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

## Warning in densfun(x, parm[1], parm[2], ...): NaNs produced
```

```r
fit2
```

```
##       mean          sd    
##   0.53588514   0.16954369 
##  (0.05361442) (0.03790678)
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
## [1]  1.62767929 -1.15299840 -0.02788475  0.39483466  0.95825966 -1.13588692
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
## [1] 0.0117
```

```r
sd(results2$thetastar)  #the standard deviation of the theta stars is the SE of the statistic (in this case, the mean)
```

```
## [1] 0.9107374
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
## t1*      4.5 -0.01271271   0.8559052
```

One of the main advantages to the 'boot' package over the 'bootstrap' package is the nicer formatting of the output.

Going back to our original code, lets see how we could reproduce all of these numbers:


```r
table(sample(x,size=length(x),replace=T))
```

```
## 
## 1 3 4 5 6 7 8 
## 1 1 1 1 1 2 3
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
## [1] -0.015
```

```r
se.boot
```

```
## [1] 0.918062
```

Why do our numbers not agree exactly with those of the boot package? This is because our estimates of bias and standard error are just estimates, and they carry with them their own uncertainties. That is one of the reasons we might bother doing jackknife-after-bootstrap.

The 'boot' package has a LOT of functionality. If we have time, we will come back to some of these more complex functions later in the semester as we cover topics like regression and glm.

