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
## 1 2 4 6 8 
## 2 4 2 1 1
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
## [1] 0.0234
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
## [1] 2.733322
```

```r
UL.boot
```

```
## [1] 6.313478
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
##    [1] 5.1 3.4 4.6 3.0 4.6 5.3 4.5 5.0 4.0 4.6 3.5 5.0 5.9 4.3 5.1 4.3 5.5 3.5
##   [19] 5.2 4.9 4.7 6.0 3.6 4.9 5.5 3.9 4.5 4.0 4.0 6.3 4.2 4.8 4.4 4.6 6.4 3.9
##   [37] 3.8 4.2 4.1 3.6 4.4 3.6 3.9 4.2 5.2 3.3 3.4 3.8 4.5 4.7 6.7 3.8 4.8 4.1
##   [55] 4.1 4.3 3.2 4.9 2.8 5.8 5.2 5.9 3.6 2.7 4.2 4.8 3.0 4.7 4.5 3.4 2.7 3.9
##   [73] 6.8 4.5 3.6 5.5 3.2 6.4 3.7 3.9 3.3 4.7 3.8 5.1 4.4 4.7 4.9 4.4 4.6 5.6
##   [91] 4.9 3.7 4.7 5.1 4.8 5.1 5.6 3.6 4.0 4.6 2.8 5.0 3.2 5.7 4.4 4.7 2.6 3.8
##  [109] 4.7 5.4 4.2 4.1 5.0 6.0 3.7 3.4 4.6 5.1 4.9 4.7 3.8 4.5 6.3 3.1 5.4 5.5
##  [127] 4.4 4.7 6.4 4.9 4.9 5.0 3.5 4.3 5.1 3.7 3.7 4.8 5.7 5.8 3.9 4.6 4.6 5.0
##  [145] 4.2 3.6 3.6 3.3 4.3 6.3 3.2 4.9 4.1 5.7 4.4 3.3 4.7 5.4 5.3 4.0 2.8 6.6
##  [163] 5.0 5.8 3.2 5.3 6.2 6.4 4.3 4.5 4.6 4.8 4.5 5.7 5.1 4.8 3.7 4.9 4.9 5.5
##  [181] 5.0 4.2 3.7 4.0 3.2 5.5 6.5 3.0 4.9 4.2 4.0 4.6 5.1 4.2 3.9 3.4 3.9 3.3
##  [199] 5.5 4.8 3.9 2.6 4.0 3.2 3.0 5.4 4.1 4.5 4.1 5.0 3.9 3.7 4.4 5.6 5.8 3.9
##  [217] 3.5 4.1 4.1 5.6 3.4 4.4 4.9 3.7 4.8 4.1 4.7 3.7 4.9 4.9 4.0 4.1 4.3 2.8
##  [235] 3.6 5.2 4.3 3.1 4.5 4.1 4.1 3.8 2.8 4.2 5.4 6.1 5.9 5.1 5.1 5.2 3.9 4.6
##  [253] 4.9 4.5 4.4 5.7 2.6 3.9 4.1 3.9 2.9 4.5 5.6 4.9 5.7 3.1 5.4 5.1 5.8 5.3
##  [271] 4.1 6.9 3.4 4.9 3.2 4.7 4.8 3.4 3.7 4.1 5.3 4.3 4.1 2.5 4.5 4.0 5.6 5.5
##  [289] 3.1 4.4 3.9 4.2 4.3 1.9 5.0 3.8 4.9 3.9 4.4 4.7 3.4 1.2 4.7 4.9 3.3 3.6
##  [307] 4.5 4.5 5.4 4.1 3.9 3.9 6.7 5.0 4.2 4.4 4.5 4.8 5.8 3.0 5.2 4.4 5.6 5.5
##  [325] 4.9 3.6 4.6 3.1 5.1 4.3 5.8 3.0 5.6 3.1 4.6 4.8 4.5 5.2 3.7 2.8 2.9 3.5
##  [343] 3.4 4.1 4.1 4.3 4.1 4.6 2.6 3.6 5.2 3.2 5.8 5.8 4.8 3.9 4.5 3.6 4.3 6.1
##  [361] 4.3 3.7 3.0 5.1 2.9 4.5 4.4 2.5 5.4 5.6 5.7 4.5 4.8 5.7 4.3 5.1 4.7 5.0
##  [379] 5.2 3.9 4.3 5.0 3.3 6.2 4.6 2.8 3.8 4.1 4.7 4.8 5.2 4.2 3.0 4.4 6.5 3.8
##  [397] 5.0 5.3 5.2 3.4 4.6 2.7 3.5 4.5 4.2 3.9 3.4 3.9 4.5 4.9 5.4 5.1 4.6 3.5
##  [415] 5.3 5.6 5.9 5.3 5.4 4.8 4.0 3.8 5.1 3.5 3.8 6.2 5.4 3.0 5.3 4.7 4.0 5.7
##  [433] 4.4 5.1 5.1 4.5 6.0 5.2 4.9 5.9 5.7 3.3 5.0 2.8 4.0 5.3 3.9 4.7 4.4 3.5
##  [451] 3.8 3.7 3.4 3.6 4.4 3.6 5.9 3.5 3.8 5.2 3.7 4.1 4.7 4.0 4.1 4.9 5.3 3.3
##  [469] 6.1 3.1 3.1 3.5 5.1 5.0 6.0 3.9 4.0 4.0 4.2 4.2 4.9 4.7 5.2 3.9 4.3 5.4
##  [487] 2.9 3.0 3.6 4.4 6.1 4.2 5.5 4.6 6.2 4.6 4.2 4.4 4.6 4.6 4.3 4.3 4.6 5.0
##  [505] 5.3 5.5 4.9 4.2 4.5 5.0 4.3 4.0 4.6 4.4 4.6 5.6 3.6 4.7 3.7 4.0 4.8 5.1
##  [523] 4.4 4.3 4.6 4.1 2.9 4.9 2.3 5.5 4.7 2.0 4.6 2.3 4.7 4.9 4.3 4.6 4.8 3.8
##  [541] 4.6 4.0 5.5 5.6 3.6 3.3 4.9 3.4 5.4 2.3 2.9 3.5 5.5 5.3 4.9 4.3 3.9 5.0
##  [559] 5.5 3.5 3.9 4.3 6.0 5.5 4.4 6.0 6.7 3.2 4.2 4.8 4.9 4.2 4.9 4.2 4.6 5.2
##  [577] 4.9 4.9 4.7 2.5 4.1 5.5 6.2 4.8 5.4 5.1 5.8 2.6 4.4 4.4 4.0 3.2 3.2 3.6
##  [595] 4.5 5.5 4.7 3.1 5.1 4.5 4.3 4.1 6.1 5.2 3.0 5.8 4.0 4.1 5.1 4.7 2.5 3.7
##  [613] 2.4 4.8 5.3 3.2 5.0 4.0 4.6 4.1 5.0 4.8 6.0 4.4 4.2 5.0 5.5 4.3 3.9 4.2
##  [631] 5.0 3.1 4.6 5.3 4.6 3.8 5.6 4.4 5.6 3.4 5.0 6.4 3.5 4.4 4.8 4.2 5.9 4.2
##  [649] 2.4 4.5 4.9 4.1 5.8 3.7 4.5 3.4 6.3 2.7 6.8 3.3 4.1 6.6 4.3 3.0 4.3 3.2
##  [667] 4.8 5.0 4.1 6.7 5.5 4.0 6.9 6.5 4.7 3.5 4.6 5.0 3.6 4.7 4.7 5.5 4.0 5.7
##  [685] 6.4 2.3 5.1 4.6 3.6 4.7 5.5 4.9 3.4 2.9 4.9 3.9 3.6 4.8 5.1 4.2 5.0 4.4
##  [703] 4.1 3.0 3.6 4.1 5.2 3.5 4.2 4.0 3.1 4.5 3.3 4.6 4.5 4.3 3.7 5.6 4.6 3.9
##  [721] 4.7 5.0 3.8 5.8 5.5 6.3 3.3 4.0 4.2 5.5 5.6 5.1 2.6 4.6 4.8 5.4 3.0 5.8
##  [739] 4.3 4.0 6.3 3.9 4.0 3.7 4.0 2.7 4.6 4.2 3.6 5.3 4.6 3.7 5.7 4.2 3.8 3.9
##  [757] 4.2 3.9 4.2 4.0 3.2 3.7 4.5 3.5 5.1 4.5 5.0 4.7 5.6 5.1 4.7 4.7 3.9 3.8
##  [775] 5.0 6.1 3.6 4.6 5.9 3.3 2.7 3.8 4.4 4.5 5.4 4.6 4.4 3.6 5.2 4.7 4.9 4.8
##  [793] 3.7 4.3 5.6 5.2 5.7 4.3 4.9 5.0 5.1 3.4 4.5 4.8 5.4 4.4 4.3 5.6 5.0 5.5
##  [811] 5.0 4.7 5.2 3.7 6.3 3.0 6.1 5.0 5.3 5.1 3.8 5.1 3.5 3.3 6.1 3.8 5.5 4.6
##  [829] 6.1 4.8 5.4 4.1 3.8 4.7 4.8 3.8 4.3 6.7 4.3 3.8 4.2 3.5 4.2 3.2 3.5 3.4
##  [847] 4.9 3.5 4.6 5.4 4.2 4.2 4.3 3.8 4.0 5.6 5.0 5.0 4.5 4.5 4.4 5.1 4.4 6.1
##  [865] 4.9 2.8 4.8 3.2 4.7 3.1 4.5 4.0 3.5 4.6 4.0 2.6 5.2 4.8 4.2 5.7 2.9 3.9
##  [883] 4.5 5.7 3.7 4.0 3.6 5.1 5.3 5.3 3.8 4.7 5.2 5.5 4.8 3.5 2.7 5.4 5.2 3.6
##  [901] 3.4 3.0 4.1 4.9 5.0 3.1 4.0 3.2 5.0 4.4 4.9 3.7 4.6 4.1 3.1 6.1 6.4 3.6
##  [919] 4.4 4.4 4.5 4.9 5.4 5.4 2.9 4.2 5.3 5.5 4.5 4.8 4.8 5.2 5.2 4.2 4.4 4.6
##  [937] 3.0 4.4 3.7 5.3 4.3 4.4 5.3 5.1 5.8 5.3 3.7 4.8 3.2 4.0 4.2 5.9 2.4 3.7
##  [955] 4.5 4.4 3.8 3.8 6.3 3.9 3.2 5.2 6.0 3.7 3.0 4.5 3.9 4.3 5.3 4.6 5.1 4.9
##  [973] 4.1 3.8 2.5 4.2 2.8 4.1 4.8 3.9 4.6 4.6 5.3 4.9 4.0 4.0 6.2 5.0 3.0 3.7
##  [991] 4.7 2.9 4.3 4.5 4.2 4.8 6.9 5.6 3.0 4.5
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
##    [1] 4.8 4.3 2.9 4.8 5.5 4.7 6.2 6.0 5.7 3.4 3.3 4.4 4.5 4.0 4.8 5.8 4.0 7.2
##   [19] 4.6 5.1 3.5 4.0 5.3 3.5 4.3 6.1 6.0 4.3 3.8 4.0 4.5 4.0 3.4 4.2 5.2 4.8
##   [37] 5.0 2.3 5.0 4.0 3.3 4.4 4.6 5.9 5.7 3.0 5.4 4.9 3.5 5.8 4.3 3.7 2.7 5.2
##   [55] 6.1 4.9 4.5 4.1 4.5 4.4 3.6 3.3 4.3 4.9 3.4 5.1 5.8 4.3 3.6 5.1 5.3 5.5
##   [73] 5.1 3.8 4.8 4.4 4.6 4.8 3.5 6.1 4.6 5.6 4.9 5.6 4.8 4.6 5.3 4.7 4.7 5.5
##   [91] 3.3 3.8 2.6 4.8 4.4 3.7 3.7 4.8 4.1 3.8 2.7 2.9 5.6 5.0 5.8 4.1 4.9 5.6
##  [109] 5.2 4.7 5.2 3.9 4.0 4.2 4.6 4.3 4.0 2.8 3.9 4.5 4.8 5.7 4.4 2.8 4.6 5.9
##  [127] 4.4 3.7 3.8 5.0 2.1 6.2 4.6 4.4 4.3 2.9 5.3 4.5 5.5 6.1 2.7 3.8 4.3 4.5
##  [145] 3.2 5.2 4.8 3.2 4.3 4.5 3.0 6.1 5.6 5.7 5.2 3.4 2.3 2.7 5.3 4.4 4.7 4.0
##  [163] 4.8 5.9 4.2 4.4 5.8 3.0 2.6 4.7 2.7 4.3 4.5 3.8 4.4 3.1 6.1 4.6 5.0 5.8
##  [181] 3.7 3.4 5.3 5.3 4.1 3.1 3.8 5.0 3.5 3.3 6.6 4.0 5.5 6.0 4.1 5.2 5.3 4.7
##  [199] 5.2 5.6 4.2 5.0 4.7 2.8 5.5 4.1 4.4 4.4 3.9 4.9 3.1 4.7 5.5 3.9 4.1 5.4
##  [217] 4.3 3.8 4.4 3.8 5.3 4.5 4.3 4.9 5.2 3.8 4.2 6.3 4.5 4.1 4.7 5.0 3.9 3.3
##  [235] 5.8 5.1 3.0 5.6 5.5 4.6 4.7 3.8 4.2 5.8 3.0 4.5 3.9 3.7 4.7 4.4 4.0 4.9
##  [253] 5.3 4.7 4.6 3.7 5.3 3.3 5.2 5.1 5.7 2.9 3.4 4.8 3.6 3.3 5.9 5.4 5.0 5.8
##  [271] 4.9 3.4 3.4 5.7 3.2 4.0 5.2 3.8 4.8 5.1 5.6 4.6 3.8 5.0 3.8 3.6 4.3 4.6
##  [289] 5.1 3.0 2.6 3.2 3.5 5.2 5.6 3.0 5.3 4.4 4.9 3.3 4.1 4.4 5.1 6.1 4.9 3.7
##  [307] 3.5 4.2 4.8 4.2 5.5 4.2 4.2 4.8 5.3 5.7 5.2 5.4 3.2 5.2 5.2 4.2 3.2 5.4
##  [325] 4.8 4.2 4.4 5.5 6.3 4.6 3.7 3.0 5.0 3.6 3.0 4.5 6.7 4.8 4.2 5.1 3.3 3.1
##  [343] 4.0 4.5 5.2 4.2 5.5 7.0 3.7 3.8 6.6 5.6 5.4 4.2 4.3 4.1 4.0 4.4 5.0 3.5
##  [361] 5.0 4.6 5.2 5.3 5.0 3.1 4.9 4.1 4.4 5.7 5.1 2.5 3.5 4.9 5.0 3.6 4.2 3.7
##  [379] 5.3 4.4 6.1 4.7 4.1 5.2 6.6 5.0 4.9 4.8 4.8 4.3 4.5 4.1 5.2 5.0 5.3 5.1
##  [397] 5.5 5.9 3.5 3.1 4.7 4.0 5.6 5.4 6.0 5.0 3.7 4.8 4.2 3.4 4.3 4.1 6.4 5.0
##  [415] 5.2 4.6 4.5 3.7 4.3 4.9 4.6 4.1 4.0 4.6 5.7 4.5 4.1 3.4 5.1 6.3 5.7 4.1
##  [433] 4.5 5.5 3.8 3.0 4.8 3.7 5.5 6.1 5.0 4.2 3.9 4.3 4.5 5.5 4.3 5.4 5.3 4.6
##  [451] 6.0 4.0 3.5 3.8 4.2 4.1 5.3 5.5 5.4 4.3 3.5 3.8 4.3 4.2 3.0 6.4 4.2 3.5
##  [469] 3.6 4.8 4.6 5.1 5.5 3.4 3.7 2.5 3.0 4.6 3.5 4.5 3.9 5.9 4.3 6.2 3.7 5.3
##  [487] 3.8 4.6 3.0 5.3 4.3 5.1 6.1 6.5 5.4 2.0 4.9 4.0 4.4 3.4 6.4 3.6 5.7 6.0
##  [505] 3.6 2.6 5.3 4.8 4.4 4.9 4.8 4.4 5.7 3.1 4.0 4.0 5.5 4.2 3.8 5.1 5.8 6.2
##  [523] 3.4 4.1 4.1 3.6 4.3 4.9 5.2 5.8 3.7 5.0 5.4 5.0 4.3 5.4 4.9 4.5 5.7 4.1
##  [541] 4.5 5.2 5.4 5.3 5.3 3.9 6.1 4.3 4.6 6.4 4.7 3.2 4.9 6.6 3.0 3.7 5.4 3.7
##  [559] 3.8 4.2 4.3 4.7 4.9 5.1 4.4 3.4 3.8 2.9 3.6 3.9 3.2 4.5 4.8 5.3 4.1 4.0
##  [577] 4.3 3.6 4.8 4.2 3.1 3.5 4.2 3.7 3.5 4.0 4.2 3.1 2.7 5.5 5.0 3.0 4.9 5.6
##  [595] 3.9 6.0 4.8 4.3 2.3 3.6 4.7 2.8 4.4 4.1 4.6 4.0 4.3 4.5 3.8 5.8 4.6 4.0
##  [613] 5.2 4.0 5.0 3.7 4.9 4.3 5.1 4.6 3.4 3.8 5.0 4.9 5.8 4.5 4.9 4.6 3.4 4.7
##  [631] 4.9 4.4 5.2 4.2 3.6 4.5 2.6 4.5 4.5 4.7 3.0 3.9 3.0 3.3 3.9 3.8 4.6 3.5
##  [649] 4.9 3.5 5.5 4.0 3.4 4.2 4.9 5.4 4.5 4.0 3.5 3.4 4.9 5.2 5.3 4.1 5.8 3.8
##  [667] 4.2 6.2 6.0 4.4 5.0 4.8 3.4 4.1 3.7 4.1 5.1 4.8 4.2 3.2 4.5 3.6 6.0 5.2
##  [685] 6.0 4.6 4.7 5.0 4.5 4.7 3.5 4.5 6.0 2.7 4.8 3.5 3.5 4.7 4.7 4.6 4.6 3.8
##  [703] 4.6 5.7 4.8 4.2 3.5 4.0 3.3 3.8 4.9 4.8 3.7 4.8 2.9 5.8 5.3 4.1 4.0 4.6
##  [721] 4.9 3.8 3.0 3.6 5.1 5.4 3.0 6.4 4.4 5.7 1.7 4.2 2.8 5.1 4.7 4.2 3.1 2.8
##  [739] 3.8 4.8 4.5 3.1 3.4 4.2 4.8 4.7 3.2 5.2 3.8 4.7 4.5 4.6 5.9 4.4 4.4 4.5
##  [757] 3.8 4.6 4.9 3.7 2.7 3.7 2.1 5.3 4.2 3.9 4.0 5.9 4.8 3.4 5.8 5.3 3.5 5.8
##  [775] 2.5 3.3 3.9 3.9 4.0 5.0 4.0 5.0 3.8 4.0 3.6 4.5 4.3 4.1 4.9 5.9 6.3 6.7
##  [793] 3.9 5.6 3.6 4.7 4.3 5.1 5.5 4.0 5.8 4.0 4.9 5.2 3.2 4.6 5.5 5.4 6.0 5.1
##  [811] 4.3 4.2 5.9 5.1 5.5 6.4 6.4 5.3 3.1 3.7 2.9 3.0 4.0 3.3 4.5 4.3 4.3 4.3
##  [829] 5.6 4.7 4.0 3.8 4.7 4.8 3.7 4.7 3.9 5.0 5.4 4.2 4.0 4.9 3.6 4.6 4.1 5.0
##  [847] 3.8 4.6 5.1 6.0 5.0 2.9 5.1 5.2 4.0 4.3 3.8 4.8 3.7 5.7 5.8 4.2 4.7 3.2
##  [865] 4.0 5.2 4.7 4.4 4.5 4.1 5.3 4.7 4.3 4.2 4.8 3.6 3.5 3.9 5.4 4.7 4.5 5.7
##  [883] 3.6 3.0 4.5 3.5 5.8 4.1 3.9 5.2 4.0 5.5 6.4 4.2 4.8 5.6 4.6 5.0 3.1 4.3
##  [901] 5.3 5.6 5.1 4.9 4.8 3.6 5.8 4.3 7.3 5.0 4.6 3.5 4.2 4.7 5.3 4.2 4.5 4.5
##  [919] 5.8 4.7 4.3 4.3 4.0 3.1 5.1 5.5 4.1 5.3 5.5 5.0 5.9 5.7 4.5 5.5 3.3 4.6
##  [937] 4.2 4.3 4.9 3.9 3.1 3.5 4.3 5.9 5.5 5.1 4.2 5.0 5.8 4.5 4.6 4.5 5.4 3.0
##  [955] 3.6 4.0 4.8 5.4 6.9 3.5 5.9 4.3 5.7 4.0 4.3 5.5 4.6 4.2 4.3 4.6 4.3 4.7
##  [973] 4.5 4.5 4.7 4.0 3.5 5.2 6.3 4.7 3.5 4.5 4.2 3.2 4.4 5.8 4.3 5.3 2.6 4.5
##  [991] 5.4 3.0 4.0 5.2 4.8 4.3 5.8 4.5 4.2 4.0
## 
## $func.thetastar
## [1] -7e-04
## 
## $jack.boot.val
##  [1]  0.52594595  0.38252149  0.31682243  0.15301205  0.03142077 -0.07178771
##  [7] -0.16611570 -0.31120219 -0.34545455 -0.52619048
## 
## $jack.boot.se
## [1] 0.9829366
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
##    [1] 3.8 3.3 4.3 3.3 5.0 4.0 1.9 5.7 5.4 3.7 3.4 4.9 5.5 4.0 5.7 5.0 5.7 4.8
##   [19] 5.8 4.3 5.0 4.9 3.5 4.4 3.8 4.2 3.5 4.8 5.5 3.5 5.9 3.5 5.4 3.9 3.8 3.6
##   [37] 3.9 5.0 4.3 5.0 4.2 4.5 5.4 5.3 4.1 7.2 3.8 4.4 3.8 5.7 4.4 2.4 3.8 5.4
##   [55] 5.0 4.8 4.3 4.9 6.1 3.7 4.1 4.8 5.0 4.2 5.1 4.1 5.6 4.4 3.8 5.2 3.1 2.6
##   [73] 3.6 4.3 4.8 6.0 4.2 4.2 4.5 4.2 3.2 5.3 3.8 5.4 4.4 4.3 6.0 3.0 4.5 3.8
##   [91] 3.7 4.7 4.8 3.7 3.8 3.9 3.7 4.5 5.3 4.8 2.7 4.9 4.0 5.5 5.9 4.0 5.1 4.7
##  [109] 4.1 6.6 5.0 3.9 4.9 4.8 4.4 2.7 5.6 4.4 4.8 4.5 4.3 2.6 5.9 4.6 5.6 4.1
##  [127] 4.2 6.2 3.9 4.6 4.4 3.4 4.5 2.4 3.1 3.6 4.0 5.7 5.7 3.4 5.1 5.5 5.1 4.2
##  [145] 4.2 3.8 4.4 5.3 4.3 5.9 6.5 3.8 4.3 4.4 4.4 5.4 5.0 4.0 4.7 6.7 4.4 5.1
##  [163] 4.9 3.1 5.4 4.7 3.6 2.8 3.1 5.6 3.9 4.7 2.8 5.1 4.7 3.6 4.7 3.9 5.1 5.0
##  [181] 4.9 4.3 3.7 3.8 5.1 4.3 6.0 4.3 4.8 3.4 5.0 3.2 5.1 3.6 5.6 4.8 2.5 4.0
##  [199] 4.2 4.8 3.6 5.9 5.4 5.4 5.1 3.9 4.7 3.8 3.7 4.5 5.9 5.5 4.1 6.1 6.0 5.1
##  [217] 5.4 4.8 4.4 5.3 1.8 4.3 6.0 5.8 4.2 3.7 5.0 4.6 4.6 5.1 6.5 3.2 5.0 2.5
##  [235] 5.1 3.8 5.5 4.5 4.1 3.9 4.0 4.3 6.1 4.4 4.1 5.8 5.8 4.3 2.5 3.3 4.2 5.9
##  [253] 4.1 3.6 3.7 5.8 5.1 5.8 3.8 6.0 4.7 3.9 4.1 5.9 5.9 2.9 5.0 4.2 4.1 4.3
##  [271] 4.0 5.1 4.5 4.6 4.1 5.6 2.5 5.4 4.8 4.1 4.8 3.3 4.6 3.6 4.8 4.8 4.3 4.5
##  [289] 3.4 4.2 4.6 3.8 4.4 5.6 6.0 3.7 5.3 4.0 5.1 6.4 3.3 3.8 5.8 4.3 4.2 5.5
##  [307] 3.4 3.6 3.0 4.3 4.6 5.6 5.6 3.8 3.3 5.2 6.2 4.3 5.3 4.3 4.5 5.9 4.4 5.4
##  [325] 5.9 4.9 4.9 4.8 4.7 5.7 4.2 5.6 4.3 6.4 1.9 7.1 4.1 2.4 4.3 5.9 5.2 3.1
##  [343] 5.3 3.9 4.6 5.2 4.7 5.2 3.7 5.7 5.3 4.3 4.6 5.8 5.0 3.8 5.3 4.7 5.8 6.6
##  [361] 2.9 4.6 3.6 3.2 4.7 3.2 4.5 5.0 5.0 4.6 4.8 4.3 3.6 3.0 4.5 3.7 5.9 2.6
##  [379] 4.4 4.3 4.5 4.2 3.8 4.3 3.7 4.1 2.7 5.8 5.5 5.2 5.2 5.4 4.9 3.7 3.6 4.5
##  [397] 4.4 5.2 5.0 6.0 2.9 3.0 4.3 5.5 6.1 5.1 3.4 4.8 2.7 5.9 4.0 4.9 4.9 4.9
##  [415] 4.5 4.7 5.2 4.0 4.5 4.3 4.5 3.8 4.2 4.7 5.2 2.9 4.1 4.0 5.4 4.1 4.7 5.0
##  [433] 4.1 4.6 4.3 5.4 3.5 2.7 4.0 4.8 3.4 3.7 5.1 3.3 4.4 4.9 3.9 5.0 6.4 5.3
##  [451] 6.0 4.9 4.3 4.7 5.2 3.4 5.4 5.9 3.8 4.3 4.2 5.5 5.7 3.2 4.3 2.8 5.4 6.2
##  [469] 5.2 4.4 4.9 5.6 6.0 4.8 5.4 5.3 4.5 4.2 5.4 5.5 4.5 4.0 3.9 5.7 4.2 3.7
##  [487] 4.0 4.9 4.0 4.9 5.1 5.6 3.8 3.1 4.9 4.8 4.2 5.2 3.8 4.4 5.3 5.9 3.5 4.7
##  [505] 4.6 3.7 5.5 4.7 5.2 4.4 5.2 4.4 5.5 4.2 5.4 5.0 5.0 4.2 5.2 5.5 4.8 4.9
##  [523] 4.9 4.3 6.0 3.4 5.4 4.3 5.6 4.8 4.2 3.9 5.1 4.4 3.2 2.6 4.6 4.8 5.3 3.7
##  [541] 4.4 4.7 3.3 5.7 5.6 4.4 5.4 3.2 4.2 4.2 4.7 6.1 4.8 4.2 4.7 4.1 4.0 3.1
##  [559] 5.1 4.5 3.9 5.1 4.0 3.8 4.6 4.8 4.2 4.0 4.5 5.6 5.5 5.6 4.2 5.2 4.3 2.8
##  [577] 5.7 4.2 2.2 4.2 4.0 5.0 3.5 3.5 5.6 4.6 3.4 3.7 3.4 4.5 4.3 4.9 3.6 4.7
##  [595] 6.0 4.2 4.2 4.8 5.2 5.3 3.3 5.0 5.8 4.4 6.1 5.6 5.4 6.1 5.3 5.0 3.9 4.3
##  [613] 4.1 5.4 5.0 5.2 4.3 3.6 3.9 5.7 4.9 4.7 4.6 4.3 4.6 4.5 2.4 3.1 5.0 4.1
##  [631] 5.2 5.7 5.2 5.3 5.8 3.5 5.2 6.0 4.1 3.0 5.1 4.4 3.8 3.5 5.6 4.6 5.0 5.4
##  [649] 4.3 3.4 5.5 4.7 4.6 4.2 5.0 3.4 4.5 4.3 5.2 4.7 5.2 3.6 3.8 3.5 3.9 5.0
##  [667] 4.5 3.7 4.5 5.2 3.2 5.6 5.3 6.9 4.0 2.9 3.6 6.5 5.5 4.2 5.0 4.7 3.2 4.8
##  [685] 4.1 5.7 4.4 4.8 3.9 5.6 3.8 4.9 2.7 4.3 4.7 4.9 2.7 3.7 4.2 4.6 4.8 5.5
##  [703] 4.1 3.8 3.5 4.2 3.1 4.3 4.1 5.2 3.9 5.6 4.7 4.5 4.7 3.7 5.4 5.4 4.9 3.9
##  [721] 3.9 5.5 2.8 6.1 4.4 4.6 5.4 5.6 4.7 4.7 5.3 4.4 4.3 4.6 5.4 4.2 3.5 3.3
##  [739] 4.3 3.7 4.7 4.4 3.2 5.6 3.7 2.9 3.6 4.5 4.7 5.3 3.9 4.7 4.4 3.8 5.6 5.4
##  [757] 4.2 5.3 5.1 6.9 4.9 4.7 5.3 4.4 4.4 4.4 5.2 4.6 5.3 5.3 4.3 4.1 4.6 4.4
##  [775] 5.1 4.4 5.1 4.8 4.6 4.7 3.3 7.5 2.2 3.8 3.5 4.4 5.7 3.9 4.5 5.4 3.8 3.3
##  [793] 3.0 4.0 4.6 4.2 4.1 5.9 3.9 3.4 3.2 1.6 3.0 4.5 6.1 4.6 4.0 3.9 4.7 3.2
##  [811] 4.3 4.8 4.9 5.0 5.0 4.2 6.1 4.3 2.5 3.2 3.7 4.7 2.6 5.2 3.1 4.2 5.0 5.1
##  [829] 4.3 4.3 4.9 4.6 5.2 3.3 4.7 4.2 3.6 4.4 4.2 3.6 4.2 2.9 4.5 5.5 3.9 2.0
##  [847] 4.6 4.4 4.6 3.6 5.0 5.6 4.4 5.9 4.2 5.1 4.4 3.0 4.5 3.4 2.9 5.1 4.4 5.0
##  [865] 3.5 4.3 3.3 5.2 4.3 4.5 3.9 5.3 5.3 4.0 3.5 3.2 5.3 5.1 5.8 4.7 4.2 4.0
##  [883] 3.1 5.5 4.9 2.8 3.4 5.2 4.1 4.1 4.2 4.9 4.6 5.0 3.8 4.8 4.8 4.8 5.1 4.8
##  [901] 4.5 3.9 4.1 5.4 4.3 4.3 5.1 5.4 4.5 3.3 4.7 3.2 5.0 2.4 3.0 4.7 3.1 5.3
##  [919] 4.4 4.6 4.5 4.8 3.7 3.8 5.9 6.3 4.3 4.7 4.6 5.2 3.2 5.0 5.3 3.5 6.0 4.2
##  [937] 5.2 4.8 5.3 5.5 2.4 4.6 5.1 4.6 5.0 3.5 5.8 4.2 4.1 3.5 3.9 4.6 3.8 2.7
##  [955] 4.1 3.7 5.8 4.4 4.3 6.7 4.8 5.1 4.5 4.7 4.5 5.0 4.0 4.9 3.9 4.2 6.0 5.3
##  [973] 5.0 2.9 4.0 5.4 3.7 5.2 4.8 5.4 3.8 4.5 3.9 3.3 4.7 4.6 3.5 3.3 4.6 3.7
##  [991] 4.8 4.0 3.8 5.7 4.2 4.0 4.9 4.9 3.7 4.6
## 
## $func.thetastar
## 72% 
##   5 
## 
## $jack.boot.val
##  [1] 5.4 5.3 5.4 5.4 5.2 5.0 4.8 4.7 4.7 4.5
## 
## $jack.boot.se
## [1] 0.9785704
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
## [1] -0.01784702
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
##      shape       rate   
##    7.580041   13.043926 
##  ( 3.317973) ( 5.902945)
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
## [1]  0.4319327  0.8480178 -0.6252574  0.6734304  1.3599301  0.2771909
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
##    [1]  0.2032292567 -0.6054398347 -0.9670463531 -0.5999067330 -0.1419342109
##    [6]  0.6802170999 -0.0268460396 -0.0918066170  0.3936369769 -0.2036518634
##   [11]  0.6412542752 -0.4294029667  0.6180758075  0.2536043849  0.1113136204
##   [16]  0.2978816546  0.2699325133  0.4244406629 -0.0293828124 -0.4354835647
##   [21]  0.4062037621  0.6788182293 -0.1819524150  0.3957932308 -0.2116499555
##   [26] -0.6032687375 -0.6978305729  1.0665185269 -0.0226794680  0.2016074463
##   [31]  0.6663690617  0.4082721528  0.6306997824 -0.2232902304 -0.2467828915
##   [36] -0.4534889765  0.4414502565 -2.1324083363 -0.4399675138 -0.2390890957
##   [41] -0.1270262155 -0.4208727009  0.4944537457 -1.4658719752  0.4090137866
##   [46] -0.2538279523  0.5240732403 -0.1130164416  0.3164266358 -1.0797105065
##   [51] -0.4699280538  0.8740594506 -0.1048950639 -0.6673104285  0.1073413761
##   [56]  0.0177653430 -0.2149021167 -1.3052962020 -0.4991177220 -0.6817361742
##   [61]  0.3046362163  0.3746576550 -0.7960192947  0.5173777857 -0.9719033082
##   [66] -0.8390479476 -0.5074187435 -0.7581982995  0.0392562027 -0.9043970617
##   [71]  0.3735194881  0.0145822954 -0.2461309201  0.6013646246 -0.5349705179
##   [76] -0.7539352155 -0.9165495344  1.0566407937 -0.0984852525  0.0907584023
##   [81]  0.0749121061 -0.3943647556 -0.5943146199 -0.8420963076 -0.8347652737
##   [86] -0.7838031888  1.0546195750  0.0268831427 -0.4020530507 -0.2233735367
##   [91] -0.0942433615 -0.0578820085  0.8578773492 -0.7744985325 -0.1788931501
##   [96] -0.1128808148  0.2639150842 -0.2000599923  0.1532009304  0.5058551770
##  [101]  0.0064662837 -0.0456902276  0.3598576595 -0.3659238543  0.7047733216
##  [106] -0.2507971228  0.3792966066 -0.4178414424 -0.4902014097 -0.5946849733
##  [111] -0.2305511671 -0.4575076696 -0.0072972873 -0.0856359461  1.0546572337
##  [116]  0.5535846483 -0.1597322597 -0.4039092226 -0.5608807037  0.4575988264
##  [121] -0.5751721961 -0.7761612729  0.2675745532 -0.3688103210 -0.0097454729
##  [126]  1.1858391711 -1.2652859423  1.2844329425  0.5757546717 -0.0997857366
##  [131] -0.0670869629 -0.1006245663 -0.3582957817  0.3232261069  0.2438214188
##  [136]  0.5980204911  0.1060195898 -0.0896564216  0.0705718931 -0.2333764735
##  [141]  1.2402428321  0.2761343507 -0.2683973582 -0.1727141720  0.2935765359
##  [146] -0.5178510479 -0.5700103321  0.0542173630 -0.3388224592 -1.7459451819
##  [151] -0.1066662832  0.1266762969  0.2856114580 -0.5010758825 -0.3092253646
##  [156] -0.2289836500 -0.0624700536 -0.3995806198  0.6073609556 -0.3388245952
##  [161] -0.3624424117  0.2126463519 -0.0448183831 -0.0456298027 -0.1523688663
##  [166] -0.0818799255 -0.3863811041  0.3843182459  0.0613873829  0.5655399181
##  [171] -0.2285810342 -0.0007977313  0.0516043931 -0.4901581151 -0.0171411203
##  [176] -0.0260840491 -0.2513321866  0.0418270199 -0.6302158974  0.3408679853
##  [181] -0.0445097558 -0.7593367949  0.6738523055 -0.1785238413 -0.4459371558
##  [186] -0.7276467320 -0.5449609767 -0.1145586386 -0.1743453061 -0.4192852206
##  [191]  0.1218404008 -0.6009967193 -0.0375202056  0.4037035849 -0.0711912793
##  [196]  0.2176407652 -1.5146583291 -0.3997453042  0.7862566217 -0.2509164323
##  [201] -0.1670366397  0.2292161280 -0.4687487200 -1.2021537081 -0.1294123065
##  [206] -0.5117639794 -0.3578690860 -0.6440103615 -0.3851159449 -0.4901759088
##  [211]  0.0658630234 -1.4404692988 -0.3773647805 -1.0610543930 -0.3968724125
##  [216]  0.2130965722  0.0467227705 -0.0993300313 -0.8535734211 -0.6140977351
##  [221] -0.5541626382 -0.6009967193  0.4435719033 -0.0587108781 -1.0911905165
##  [226]  0.5625542658  1.1295793695  0.6235823175  0.1644312931 -0.7160788967
##  [231] -0.0956408867 -0.1470194410  0.0338170698  0.3862709041  0.0554639713
##  [236]  0.5101856734 -0.6707748016  0.6380192892  0.2768219064 -0.1495686342
##  [241]  0.0591258923  0.5420034194  0.0828259352 -0.6293291592  0.0602692237
##  [246]  0.1417506301 -0.4201337867  0.2474178408  0.3021499966 -0.8593165415
##  [251] -0.9083825529 -0.1882115979  0.5431570203  0.4271536994  0.2662258324
##  [256] -0.5073443851  0.7365587549 -0.1167645664 -0.1430312582 -0.0609821444
##  [261]  0.3152528030  1.1687075317  0.1161970972 -0.2491227883 -0.1484983504
##  [266] -0.2370656466 -0.7674581941  0.8537301735 -0.1507576030  0.5661711837
##  [271] -0.1208416996 -0.1943658292  0.1071179529  0.1251146806  0.6561190930
##  [276] -0.0610584234 -0.1738913110 -0.3092459942 -0.5862357398 -0.7887644515
##  [281] -0.4925175453 -0.9159969692  0.6224645821  1.2163857983 -0.3978328894
##  [286] -0.6277987097  0.1269960469 -0.3657184636  0.5535618359 -0.1451127457
##  [291]  0.2136442143  0.4164921840 -0.1137879432 -0.5712737586 -0.0750827802
##  [296]  0.5897083675 -0.7312432207 -0.3931715655 -0.2629213820  0.8807291710
##  [301]  0.2444866982 -0.3414247565 -0.4510638908 -0.5002790897 -0.8462016140
##  [306] -0.6501314685  0.3907438078  0.3938622574  0.2432762413  0.0460175449
##  [311] -2.0245689868 -0.8564850891 -0.3994095749 -0.2989664495 -0.1480993595
##  [316]  0.9850907431 -0.3218818436 -0.0339673315  0.1781130061  0.2629530661
##  [321] -0.9129696006 -0.4957888551 -0.8977352300  0.4265822337 -0.2664036949
##  [326] -0.8261798735 -0.2138564918  0.7645488548  0.0455914120 -0.6587221633
##  [331]  0.7304878781 -0.1961297593 -0.3321428580 -0.0984852525 -0.1165357222
##  [336]  0.0081820559 -0.3026876811 -0.0097857126 -0.0545212514  0.0833339709
##  [341]  0.0925705778 -0.9473480288  0.3383068428  0.2347725565 -0.1239430427
##  [346] -0.2531951832 -0.0920616669 -0.1408143077 -0.1785238413 -0.0394768580
##  [351] -0.1123162683 -0.1610801571 -0.8155319458 -0.4508246801 -0.4876153725
##  [356]  0.4837457434 -0.1203814971  0.5826377357 -0.5875221292 -0.2629213820
##  [361] -0.3097562627 -0.4415077998 -0.7497898235 -0.4887873713  0.3966665689
##  [366] -0.0752361889  0.2633936972  0.1340905063  0.2213532601  0.0025883454
##  [371]  0.0341575031 -0.2278910446  0.2714689337 -0.3081560782  0.6229696433
##  [376] -0.3598748174  0.2622615364  0.1181266070 -0.4084777632 -1.0599070346
##  [381] -0.0281682520 -0.3119043615 -0.5371825724 -0.2477612976  0.2310048487
##  [386] -0.8383773650 -0.2346358927  0.6539695994 -0.0157213177  0.3050623702
##  [391] -0.3504338818  0.0005826699  0.3518871223  0.0121850106 -0.0073597108
##  [396]  0.1662644267  0.0233931110 -0.2285561785  0.0635339469 -0.0996311614
##  [401] -0.6839566239 -0.9582224046 -0.4753688628  0.7432820121 -0.0668378520
##  [406] -0.4795399484 -0.4287169150  0.0138062124 -0.2688453680  0.7725890217
##  [411]  0.3890483610 -0.1869954656 -0.6404872463  0.7583577599  0.2894742550
##  [416]  0.0551219463  0.0408112892  0.1335799860  0.8735991644 -0.1426524096
##  [421] -0.9211182529 -0.1775026528 -0.1927504256 -0.8282872805  0.6252034614
##  [426] -0.0121347436 -0.8850143161  0.3191236142  0.2395930851  0.4801226917
##  [431] -0.5213112424  0.0704517616  0.4980528802 -0.4503071604  0.9479417902
##  [436]  0.0355426289 -0.1201151912  0.1665981301  0.2664597913 -0.0854124668
##  [441] -0.0062051681 -0.2997719628 -0.2279059952  0.6147431285 -0.8953235549
##  [446] -0.0569325990 -1.4070319785  1.5444568964 -0.0317276554 -0.3449734093
##  [451]  0.0777424584  0.2057337506 -0.0436036055  0.2111776573  0.2477062634
##  [456] -1.3363049693 -0.0155926676 -0.1852144582  0.6737641719  0.7739305015
##  [461] -0.4666663026 -0.0083063876 -0.3759261287 -0.1764671562  0.1621861693
##  [466]  0.6436768505  0.0996993252 -0.8389928164  0.1619380641 -0.4901872081
##  [471] -0.7957260841  0.1192741928  0.4530013632  0.3280384438  0.3188903070
##  [476]  0.1866541621 -0.0766030696  0.1568748260 -0.0692547241  0.2205805987
##  [481]  0.6803373794 -0.2690137334 -0.2468703064 -0.2110054783 -0.1773041813
##  [486] -0.6294955952 -0.1895265422 -0.1179750076 -1.0574336534  0.1468505575
##  [491]  0.1652560774  1.3453603090 -0.0072292678 -2.1454423115 -0.0533027165
##  [496]  0.0551865886  0.0893512535  0.0857570848  0.1437628542 -0.5223056928
##  [501]  0.5314026705 -0.2271913028  0.3227182104  0.2483292642 -0.2424198678
##  [506]  0.1180784005 -1.2480289742 -1.1914208058  0.4239355929 -0.1296443681
##  [511]  0.5076247881 -0.0600691347 -0.4169173938  0.7438425445  1.2389851433
##  [516] -0.3269132173  0.1248185295  0.5480643277 -0.6544788911 -0.3457837867
##  [521] -0.4708525601 -0.7215139383 -1.2410373427 -0.0848629273 -0.0049235372
##  [526]  0.2158884470 -1.3284908626  0.0210119898  0.3982100939  0.0031691532
##  [531]  0.7549820048  0.2008129279 -0.4714800189  0.2127850521  0.3097868208
##  [536] -0.7404819105  0.0650721649  0.1567158376 -0.0302456948 -0.2810371804
##  [541] -0.0673004265 -1.2995718727  0.0730658760 -0.0429628267  0.4242461295
##  [546]  0.0361715097 -1.1698993207 -0.7402803543  0.6137429545 -0.4425349779
##  [551] -0.9506988101  0.2443986601 -0.4392038927  1.1007446932  0.0655775874
##  [556] -0.4253358989  0.2906370137 -0.1357349742  0.2579954968 -0.1861120841
##  [561]  0.2562613364 -0.1160731309  0.1349832048 -0.1404191903  0.5595758564
##  [566] -0.6556975011  0.1532222802 -0.3546162834 -0.1044639277 -0.1046003949
##  [571] -0.8561089760 -0.2242000166  0.6319613664  0.2548294022  0.1670182257
##  [576] -0.1209508131  0.5065833929  0.4882659701 -0.3968050414 -1.2747340468
##  [581]  1.1677252820 -0.4719824971 -0.4739665172 -0.3394042464 -0.5468662106
##  [586]  1.6596541557 -0.6055774983 -0.3310087468 -0.1865692622 -0.1467735752
##  [591] -0.0654449071 -0.0954423778  0.2032658555 -1.0303779200  0.1054436617
##  [596] -0.2804060663 -0.0329165092  0.3501033621  0.3680289470 -0.0552625182
##  [601] -0.3678571010  0.1072312405 -0.0744423641  0.6848554001 -0.2672387762
##  [606]  0.3944916602 -0.5611848046  0.1487312929 -0.1130232554 -1.2719753220
##  [611] -0.2492582417  0.2088922842  0.3359547509  0.8313790067  1.0553271317
##  [616] -0.0284477411 -0.2859013146  0.3687183139 -0.3170529613 -0.3219366252
##  [621]  0.0371118032  0.2435437543 -0.3517396489 -0.1878824619  0.2881765104
##  [626] -0.1766801219  0.1220532158 -0.1474138844  0.1486990192  0.6697130228
##  [631]  1.0374917989 -0.1396591010  0.1055089556 -0.1981729854  0.2601991865
##  [636] -1.1725140065  0.2489939521  0.4835839306  1.9993272364  0.2170339068
##  [641]  0.7553917782  0.0554639713  0.8200278376  0.4935727070 -0.0642666707
##  [646]  0.0784101794  0.5779102745 -0.1140972848 -0.1046138814  0.2333731520
##  [651] -0.0835217752 -0.0703578100  0.1114231893 -0.8114977055 -1.4027644606
##  [656]  0.6366165410  0.2186258553  0.0409741072  0.9041921039  1.3424689233
##  [661] -0.0893121778 -0.6632763757 -0.0938052656  0.8333653179  0.0179312496
##  [666] -0.0764690553  0.6041900278  0.3816347831  0.1720721147  0.1080883590
##  [671] -0.2492582417 -0.3098749163 -0.2672703907 -0.5202062402  0.0390000149
##  [676]  0.2694156266  0.2783843350 -0.6106691508  1.5386342043 -0.1214429399
##  [681] -0.6878273920  0.1702054958  0.7589106454 -0.7971526808  0.3124631559
##  [686] -0.4933296397  0.3690671383  0.1192741928  0.4098933071 -0.3182628920
##  [691] -0.9236224406 -0.0547298461  0.2286171620  0.4438346804  0.0429088209
##  [696] -0.9449413611 -0.0385119925  0.1359422329 -0.3485206116  0.2700946856
##  [701]  0.1581288570  0.5936645741 -0.2932241865 -0.0947993560 -0.9676543205
##  [706] -0.8736697699  0.5788819466  1.0427368754  0.6976197916  0.2351462192
##  [711] -0.6453209622 -0.5513645897  0.0283908261  0.1925638024  0.3167874463
##  [716]  0.1585038943 -0.3227769641 -0.0296265666 -0.0311607870 -0.3234038338
##  [721] -0.2278893098 -0.6067546997 -0.1075783749  0.1970800553 -0.3585671148
##  [726] -0.1265662550 -0.7186335628 -1.0568145808 -1.2241462207 -0.4965518930
##  [731] -0.0945581340 -0.2149021167 -0.6905196696 -0.7080543208  0.6687331819
##  [736]  0.1072871278  0.3598350689  0.0704517616 -0.2216306334  0.4302120242
##  [741]  0.0127557866 -0.8969798084  0.0293060293 -0.4253358989  0.4834984612
##  [746] -0.4264255030 -0.2231265906 -0.1541145160  0.1924565125 -0.9116887937
##  [751] -0.3738788141  0.4587759999 -0.9618807728  0.6106594832  0.1454545078
##  [756] -1.4124466588 -0.5286716519 -0.3961076822  0.2126813728  0.2433770045
##  [761] -2.0071988332 -1.2182148447  0.6202945512 -0.3488237930 -0.1928499305
##  [766]  0.0060547540 -0.1403947903  0.2008129279 -0.0684454473  0.2029441949
##  [771] -0.2127892552  0.4625379936  0.2548163999  0.4878603301  0.1459768384
##  [776]  0.1817100953  0.3533688589 -0.1869902888 -1.3928907411  0.5264362022
##  [781]  0.0126170851  0.3363627156 -0.8835872572 -0.1320842942 -0.7433339856
##  [786] -0.1951900526 -0.4016496610 -0.0069449459  1.0050482316  0.4026611782
##  [791] -0.4740551157  0.2658304560 -0.3022299312 -0.3280486704  0.1180409013
##  [796] -0.4583065389 -0.6909147934 -0.5547292882 -0.1879890458 -0.4423356842
##  [801]  0.1764645522 -0.6706866291  0.4016701113  0.3511867026 -0.4998560246
##  [806]  0.1853018916 -0.1014275630 -0.0632466753  0.4065593205  0.5573107746
##  [811] -0.3685680668 -0.2579633185 -0.6481057707  0.1124145251 -0.4901945769
##  [816]  0.2876463272 -0.0808192145  1.2871690063 -0.4953602707  0.4904216321
##  [821] -0.3993696916  0.7049962105 -0.5455039773 -0.4127682926 -0.3414386239
##  [826]  0.2851936386 -0.1314045961 -0.5708752783  0.2835115003 -0.1681888276
##  [831] -0.3409118122 -1.3914569653  0.0639654167 -0.8578144397  0.1976230644
##  [836] -0.4330138788 -0.7212929744 -0.4453625033  0.2747123717 -0.7990418330
##  [841] -0.0034634252  0.0990250532 -0.9158313356 -0.2703775094 -0.0823321982
##  [846]  0.2026944005  0.5172577784 -0.8283323636  0.3739597082 -1.5573854264
##  [851] -1.5259478853  0.1977151617 -0.2237963902  0.0285775993  0.5290357588
##  [856] -0.1064522318  0.1773378583 -0.0111940638 -0.8006154983 -0.0144113393
##  [861]  0.0278177902  0.0966283963  0.0248631497 -0.7379614119  0.3877921366
##  [866]  0.5229112756  0.3224403398 -0.0164299531  0.5771229364 -0.0271428039
##  [871] -0.4352828636 -0.5407757225 -0.1164554091 -0.2349821320  0.6966332931
##  [876]  0.8886371474 -0.1310245912 -0.8880192916 -0.3434771908 -0.4856936908
##  [881] -0.3612472646 -0.0119649462 -0.4986818922 -0.0793279986 -0.7580961723
##  [886]  0.2966235070  0.6229696433  0.2126463519 -0.7910644542 -0.3711912060
##  [891] -0.9222009406  0.5634755949 -0.8068762718  0.3549688735  0.2505887833
##  [896] -0.0383254936 -0.1270601849 -0.2719193466 -0.3534899305 -0.8562323266
##  [901] -0.1998791329 -0.0878616302  0.1227401420 -0.1232570421  0.3784199901
##  [906]  0.3744873852  0.3549649411 -0.0182392653 -0.0390673654  0.3961782034
##  [911]  0.2224514430  0.6945619499  0.0204782749 -0.8746694922  0.7270880096
##  [916] -0.8069370830 -0.3787171020  0.2328037228 -0.4644006445  1.0430829507
##  [921]  0.1764645522  0.1862145391  0.1102294178  0.6831704259 -0.1492326026
##  [926] -0.2564151784  0.5865327545  0.6135243637  0.2565082419 -0.0372941164
##  [931] -0.6951896462  0.2058427185 -1.1619920767  0.0864212142  0.1868472587
##  [936] -0.4710173195 -0.0179835765 -0.0962150936 -0.3498610559 -0.0387160273
##  [941] -0.4943096571 -0.2544422101 -0.5402138951  0.4835677717  0.1644355054
##  [946]  0.2967057538  0.2383947409  0.3901013879  0.0388333736 -0.8095888950
##  [951]  0.1144875751 -0.0343634395  0.3957932308 -0.1137049769 -1.2633662700
##  [956] -0.3900358901 -0.6493042332 -0.9128154312 -0.3375157619 -0.1887846107
##  [961] -1.4296584923  0.4257012777  0.3584082176 -0.3054445013 -0.1881862433
##  [966]  0.4078758644 -0.1829245090  0.3124631559  0.0823350116 -0.1561268571
##  [971] -0.3380630838 -0.8358951292 -0.4104453350  0.0893204947 -0.7039083385
##  [976] -0.2461755605 -0.1333410439  0.9160907587  0.1777733767 -0.3064904729
##  [981] -0.1118540576 -0.1711498784  0.3834961536 -0.2407531021 -0.1189468115
##  [986] -0.6577298212 -0.0463320449  0.4739142604 -0.1064522318  0.5660041771
##  [991] -0.0777302639 -0.5186554736  0.2530542914 -0.2237668834 -0.3519348883
##  [996]  0.3604960075  0.2325274616  0.5451069959 -0.5540072158 -1.1211022895
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
##   0.58111621   0.20033308 
##  (0.06335088) (0.04479198)
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
## [1] -0.3055007  0.4195328 -0.6038593 -0.1746850  1.0321422  0.7404767
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
## [1] -0.0089
```

```r
sd(results2$thetastar)  #the standard deviation of the theta stars is the SE of the statistic (in this case, the mean)
```

```
## [1] 0.8709704
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
## t1*      4.5 -0.05825826   0.9215791
```

One of the main advantages to the 'boot' package over the 'bootstrap' package is the nicer formatting of the output.

Going back to our original code, lets see how we could reproduce all of these numbers:


```r
table(sample(x,size=length(x),replace=T))
```

```
## 
## 1 2 3 6 7 9 
## 1 2 2 1 1 3
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
## [1] 0.0362
```

```r
se.boot
```

```
## [1] 0.925541
```

Why do our numbers not agree exactly with those of the boot package? This is because our estimates of bias and standard error are just estimates, and they carry with them their own uncertainties. That is one of the reasons we might bother doing jackknife-after-bootstrap.

The 'boot' package has a LOT of functionality. If we have time, we will come back to some of these more complex functions later in the semester as we cover topics like regression and glm.

