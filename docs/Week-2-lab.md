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
## 0 2 3 4 6 7 8 9 
## 1 1 1 1 1 1 1 3
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
## [1] -0.0094
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
## [1] 2.752401
```

```r
UL.boot
```

```
## [1] 6.228799
```

Method #2: Simply take the quantiles of the bootstrap statistics


```r
quantile(xmeans,c(0.025,0.975))
```

```
##  2.5% 97.5% 
##   2.8   6.2
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
##    [1] 4.0 4.0 2.5 3.5 3.4 3.3 2.7 5.3 3.8 3.7 5.7 3.5 4.7 4.8 4.0 5.1 4.4 3.9
##   [19] 5.2 3.3 4.2 4.1 2.9 3.2 4.1 4.5 4.5 5.2 5.0 4.8 6.3 4.2 3.9 4.7 4.8 3.6
##   [37] 4.7 3.4 3.7 5.2 5.0 4.8 6.1 2.8 3.9 4.5 6.2 5.3 5.4 4.7 5.3 4.2 5.0 4.0
##   [55] 6.1 5.9 6.6 4.6 4.0 5.4 4.5 6.0 5.2 4.5 4.2 4.0 4.6 4.2 6.6 2.9 5.4 3.6
##   [73] 3.9 3.4 3.7 5.9 3.4 4.3 4.6 4.9 5.8 3.7 4.1 4.2 5.3 3.8 3.2 4.5 3.3 4.7
##   [91] 4.0 4.8 5.9 3.8 4.2 4.4 5.7 3.9 4.5 5.7 5.4 5.1 3.3 5.3 3.4 4.8 4.7 2.6
##  [109] 4.5 5.7 5.3 5.2 5.4 3.1 5.2 3.1 4.6 5.5 5.1 4.2 4.4 3.8 5.3 4.4 4.2 6.7
##  [127] 4.2 3.9 3.7 4.8 5.6 4.8 4.6 3.6 1.6 5.2 4.9 5.9 5.6 2.9 6.2 5.3 6.7 4.9
##  [145] 3.4 3.6 4.8 4.7 4.5 4.3 4.5 3.5 4.2 5.8 4.5 3.9 4.9 5.0 3.6 3.7 5.6 5.2
##  [163] 4.4 5.8 4.1 4.5 3.5 3.0 4.4 2.4 4.8 4.9 4.5 4.0 4.1 4.4 4.4 3.7 5.4 6.2
##  [181] 4.8 6.3 3.8 2.7 4.3 4.1 5.5 5.0 4.0 4.9 3.6 4.0 5.0 5.2 6.4 3.0 4.0 5.3
##  [199] 4.6 6.0 4.1 4.9 5.2 4.0 3.7 4.2 5.0 4.9 3.9 4.9 3.0 6.1 5.1 3.1 4.3 4.1
##  [217] 3.9 5.7 4.3 4.1 4.8 4.8 4.4 3.8 3.8 2.5 3.4 5.1 4.4 5.6 3.8 5.0 5.6 3.7
##  [235] 4.0 4.6 5.8 3.4 3.9 3.6 6.0 4.6 4.3 3.9 5.6 4.3 3.2 5.6 5.3 5.0 2.3 4.1
##  [253] 4.7 3.4 5.1 4.1 5.8 4.5 6.1 4.4 3.3 5.0 4.8 5.9 4.8 3.6 4.6 3.6 5.1 3.1
##  [271] 2.4 3.1 5.6 3.5 4.9 3.3 4.9 5.1 3.0 5.4 6.5 5.3 3.7 3.5 5.2 4.0 3.6 5.8
##  [289] 2.7 6.1 3.5 5.0 5.4 6.0 6.0 5.0 4.2 5.1 3.8 3.0 4.8 4.7 2.6 4.1 4.8 5.6
##  [307] 5.1 6.5 4.8 4.5 4.3 3.7 4.3 5.6 4.8 4.6 3.7 4.4 4.0 3.7 4.2 5.0 4.2 4.9
##  [325] 4.9 6.2 2.7 4.1 3.9 4.1 4.0 4.6 4.1 5.2 5.2 4.2 5.1 4.8 4.4 5.0 5.8 4.5
##  [343] 3.3 5.0 5.9 3.3 6.7 3.5 6.5 4.8 5.3 4.1 4.4 3.5 6.3 3.8 3.0 3.9 5.1 4.3
##  [361] 5.5 3.6 6.2 2.8 5.1 2.8 6.6 5.4 4.1 2.4 3.6 4.4 4.1 4.6 4.3 4.0 3.4 4.8
##  [379] 4.5 5.1 4.2 5.5 3.8 4.8 4.6 4.3 3.8 4.6 6.0 5.1 5.2 1.5 4.6 3.1 4.0 4.3
##  [397] 6.2 4.0 4.6 5.7 5.1 4.0 4.3 5.1 5.8 4.0 4.4 3.9 6.0 3.2 3.8 4.2 4.1 3.3
##  [415] 3.5 5.2 4.0 4.2 5.3 4.5 4.6 3.2 5.8 5.3 4.1 4.4 4.2 4.1 3.3 3.1 3.7 4.1
##  [433] 3.9 3.5 5.1 3.0 4.1 5.2 4.4 5.5 2.6 5.0 4.2 5.7 2.4 4.4 4.3 3.6 4.5 4.8
##  [451] 4.3 5.4 2.8 3.9 3.9 5.1 3.2 4.4 4.1 4.4 5.0 3.1 3.3 5.3 3.9 3.6 4.7 3.7
##  [469] 4.7 4.5 3.8 3.9 4.1 4.7 4.0 1.7 6.0 3.3 4.5 4.0 3.7 3.9 6.3 4.6 5.1 4.3
##  [487] 4.4 4.7 4.0 5.1 3.7 5.7 4.5 2.9 3.4 5.4 4.9 3.6 3.7 4.7 4.1 4.4 3.9 6.0
##  [505] 3.4 4.2 5.6 5.0 6.5 4.5 5.9 5.6 3.7 4.1 3.8 5.3 3.9 4.4 4.8 4.6 4.1 5.0
##  [523] 4.7 4.1 4.5 4.7 5.1 4.3 4.5 4.9 5.4 4.6 4.4 6.3 3.1 4.1 5.0 3.3 5.1 5.6
##  [541] 6.7 5.1 5.3 5.0 4.5 4.7 3.8 4.1 3.3 4.7 4.9 4.8 4.0 3.3 4.5 3.6 3.7 4.5
##  [559] 4.4 4.7 5.0 3.0 4.5 3.2 4.1 5.4 4.1 4.0 3.6 6.4 5.4 3.9 4.4 2.9 4.2 3.9
##  [577] 6.1 4.8 4.3 4.5 3.7 4.5 3.4 4.5 4.9 4.2 5.6 4.8 4.3 6.2 5.6 4.2 3.5 6.6
##  [595] 5.8 3.4 6.2 4.3 2.7 4.0 5.1 2.6 5.8 3.7 5.0 5.9 5.4 3.1 4.8 5.2 6.1 3.4
##  [613] 3.9 5.7 3.4 3.8 5.2 3.7 4.3 3.5 4.8 3.9 5.1 3.1 3.3 4.3 2.7 3.8 4.1 4.6
##  [631] 4.5 3.4 5.6 3.9 5.3 3.5 2.7 4.0 4.9 3.8 4.0 3.0 3.8 3.9 4.3 5.0 4.1 3.8
##  [649] 4.9 4.4 4.4 3.5 4.6 5.0 5.3 3.5 5.0 3.4 5.7 4.3 5.6 3.7 4.6 4.8 4.8 4.7
##  [667] 3.3 5.5 5.4 5.8 6.1 5.4 4.1 3.8 4.8 4.6 5.4 5.1 4.5 5.7 4.1 5.4 5.2 5.9
##  [685] 3.5 4.3 5.6 5.0 3.5 3.8 4.6 4.0 3.5 4.1 3.5 4.2 6.2 4.1 6.2 5.1 5.0 4.6
##  [703] 5.1 5.2 5.2 5.1 3.9 4.0 4.3 3.9 6.5 4.2 3.4 4.6 3.4 5.3 5.4 5.2 4.7 4.8
##  [721] 4.2 5.2 4.0 6.1 4.3 4.5 3.5 4.7 4.0 5.4 5.5 4.6 3.1 5.6 5.7 5.5 6.5 4.8
##  [739] 4.4 5.2 4.1 2.9 4.9 4.4 4.9 5.9 3.5 4.3 5.4 4.7 3.9 5.6 3.4 3.7 5.1 4.8
##  [757] 3.3 3.5 3.8 5.7 5.3 3.8 4.4 4.1 4.0 4.8 4.6 3.8 3.7 5.9 3.8 5.1 4.3 6.0
##  [775] 3.6 2.9 6.8 4.6 4.8 5.1 3.2 3.0 4.1 4.9 5.0 6.3 3.1 3.7 3.9 5.3 3.4 4.7
##  [793] 4.1 5.7 4.4 4.4 6.2 5.2 5.1 5.7 5.0 4.5 4.7 4.5 4.4 5.8 5.5 5.6 3.6 4.2
##  [811] 4.8 6.2 5.6 5.9 4.0 5.0 3.5 6.7 5.9 4.8 4.7 5.5 3.9 5.2 4.6 3.9 6.6 6.5
##  [829] 4.4 4.8 4.9 5.0 4.7 4.9 5.1 2.7 6.4 4.1 5.2 5.0 5.7 4.5 4.7 4.6 5.2 4.7
##  [847] 5.4 2.9 5.4 4.3 3.3 6.0 3.1 4.0 3.7 3.4 5.9 3.3 4.9 5.0 5.0 6.2 4.3 4.9
##  [865] 5.0 3.9 4.0 4.5 5.7 4.1 3.8 6.0 4.4 4.1 3.6 4.8 4.8 3.3 6.0 4.4 5.1 4.2
##  [883] 2.7 4.5 5.6 4.4 4.3 3.1 4.0 4.6 3.8 4.6 4.0 5.1 3.5 4.4 6.6 5.3 3.5 6.3
##  [901] 5.4 4.1 4.4 4.1 3.6 5.1 3.9 4.0 5.4 5.2 4.2 2.8 4.1 3.2 4.1 4.0 4.8 4.9
##  [919] 4.3 4.5 4.6 5.4 3.4 4.8 3.8 4.8 5.7 4.0 6.1 6.3 4.0 3.6 4.5 4.8 4.9 4.2
##  [937] 4.6 4.3 5.6 2.9 4.3 5.4 2.6 5.0 5.1 5.1 5.0 3.5 4.8 5.2 4.9 4.7 4.3 4.4
##  [955] 3.0 4.4 5.1 6.1 4.3 3.9 5.6 4.9 4.3 3.8 4.6 4.9 4.1 5.1 4.2 3.2 4.6 2.9
##  [973] 5.5 5.5 5.4 3.2 3.7 5.4 5.1 3.5 3.3 3.5 3.1 5.0 4.6 4.7 4.8 6.5 5.1 3.5
##  [991] 4.4 5.1 4.5 5.4 4.1 3.8 4.0 5.1 4.5 5.1
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
##    [1] 4.1 4.0 4.7 4.1 6.3 4.2 4.8 3.5 5.7 3.2 5.3 2.6 5.4 3.9 4.3 2.4 5.4 4.6
##   [19] 4.2 4.2 5.6 4.8 4.8 4.4 5.1 4.6 2.1 4.0 3.4 5.5 5.1 6.0 4.7 3.4 2.8 3.7
##   [37] 5.1 5.0 3.6 4.1 4.6 5.7 4.4 4.9 4.2 3.2 3.9 4.0 5.3 6.4 3.1 5.4 5.5 4.7
##   [55] 4.7 3.1 4.2 4.3 3.5 3.4 5.2 4.3 5.0 4.2 5.0 5.5 4.4 3.3 3.3 4.7 3.0 3.1
##   [73] 5.5 4.8 3.9 3.6 4.7 4.5 3.8 4.4 4.8 4.1 4.3 6.7 2.4 5.4 5.9 5.5 5.5 5.5
##   [91] 3.3 3.6 6.2 6.0 4.3 4.9 4.2 3.5 5.9 2.6 3.5 5.0 3.9 4.3 5.6 3.9 4.7 5.6
##  [109] 4.0 4.7 4.9 5.7 2.7 3.4 4.1 3.5 4.6 5.5 3.9 4.8 4.4 4.6 5.4 3.5 4.6 2.8
##  [127] 5.3 5.4 4.7 4.7 4.7 4.1 4.0 3.7 4.1 4.9 5.7 3.2 5.9 4.8 5.0 5.3 5.0 4.5
##  [145] 5.5 4.1 3.2 4.8 3.9 5.9 3.6 2.6 5.6 5.1 6.4 3.8 3.6 4.3 4.5 4.6 4.8 5.7
##  [163] 5.4 6.2 3.1 2.9 4.1 3.8 3.3 4.8 4.1 5.6 6.4 5.5 4.2 3.7 3.6 5.4 5.8 4.4
##  [181] 4.7 4.2 2.8 6.8 3.1 4.8 4.2 6.0 5.8 5.3 5.0 4.3 3.6 2.3 3.9 4.0 5.4 3.9
##  [199] 5.2 3.9 5.3 3.5 4.1 5.1 4.7 4.9 2.8 5.1 3.5 5.2 3.2 4.1 3.4 4.9 5.0 4.3
##  [217] 5.3 5.0 4.2 4.8 4.3 4.7 4.8 4.0 5.9 5.5 4.4 4.3 4.4 3.6 5.5 3.9 6.1 4.9
##  [235] 3.0 5.0 3.2 5.7 3.8 4.3 4.4 4.4 2.9 4.4 5.7 5.2 4.6 2.5 4.7 3.4 4.6 4.9
##  [253] 4.9 3.2 5.4 3.1 4.3 4.0 4.3 2.6 4.8 4.6 4.0 4.6 4.1 5.2 3.7 3.4 4.0 3.6
##  [271] 3.7 3.0 4.1 4.0 3.4 5.5 5.2 5.7 5.9 3.3 3.1 4.6 5.8 4.9 3.8 3.5 4.8 4.2
##  [289] 3.6 4.3 5.1 3.5 3.8 4.2 4.7 4.5 3.6 4.8 2.7 4.6 5.3 6.8 6.5 3.8 4.6 5.4
##  [307] 5.6 5.3 5.2 3.0 4.4 3.9 5.7 3.6 3.2 4.8 5.2 6.3 4.1 3.5 3.8 4.9 4.4 2.8
##  [325] 4.0 4.5 5.1 3.9 4.4 5.7 4.7 5.1 4.8 3.5 5.7 4.5 5.3 5.3 3.2 3.4 4.2 4.0
##  [343] 5.4 4.6 2.1 4.5 3.4 3.8 5.2 3.8 4.1 4.3 4.6 6.0 5.1 3.6 4.7 4.2 3.8 3.2
##  [361] 3.8 5.0 5.0 3.0 3.4 4.3 4.0 4.7 4.7 4.0 2.9 5.3 2.9 3.8 4.1 2.8 5.6 4.4
##  [379] 3.8 6.3 3.4 4.7 3.8 4.3 3.8 5.1 4.0 5.0 4.4 3.4 3.4 3.8 3.9 6.0 4.9 3.9
##  [397] 5.1 4.9 3.7 4.7 3.9 3.5 3.5 6.0 4.7 3.7 2.1 3.5 3.9 5.3 4.6 4.2 5.4 5.1
##  [415] 4.7 4.7 4.8 3.7 2.9 4.1 4.6 3.9 5.0 4.7 4.4 5.2 4.0 4.7 4.1 3.2 5.4 6.9
##  [433] 5.5 3.5 4.8 3.5 4.4 3.5 6.2 4.8 5.2 5.4 2.8 4.6 4.4 5.9 4.2 4.5 5.2 4.2
##  [451] 2.3 4.5 4.5 4.8 5.4 6.0 5.4 4.7 3.7 4.3 6.5 6.2 4.9 5.6 3.9 4.4 3.9 4.8
##  [469] 5.3 3.6 5.5 4.6 5.0 4.0 3.7 3.1 3.9 3.8 5.6 5.4 5.8 5.2 5.9 2.6 4.3 3.9
##  [487] 4.4 3.7 5.3 3.6 5.9 5.5 3.1 4.3 3.9 3.1 4.8 4.2 5.9 3.5 4.4 4.9 3.4 3.1
##  [505] 3.6 4.6 5.3 4.7 5.1 3.6 4.0 3.9 4.8 3.6 5.2 4.5 4.0 2.4 5.0 4.4 2.3 4.3
##  [523] 5.5 4.7 4.8 5.9 3.2 5.0 3.2 4.1 3.9 3.9 4.5 4.1 4.4 3.4 4.6 3.1 3.5 4.4
##  [541] 6.0 6.0 4.1 5.5 4.6 4.6 5.9 4.8 3.4 4.2 4.6 3.8 5.0 6.0 4.9 3.4 3.4 4.1
##  [559] 5.0 7.2 5.0 4.8 5.4 4.0 4.3 3.1 4.7 5.0 4.8 4.9 2.8 3.5 3.8 5.0 5.5 4.9
##  [577] 4.8 4.0 5.0 5.0 6.3 3.6 4.9 6.3 3.0 4.4 5.9 4.6 6.8 5.1 4.8 3.8 3.5 5.3
##  [595] 5.2 4.3 4.0 5.4 4.7 5.5 4.6 4.9 5.5 4.8 2.9 4.6 4.8 3.4 5.6 2.9 5.0 4.8
##  [613] 3.8 2.9 4.4 6.2 4.6 4.8 4.9 4.8 4.1 4.5 4.9 4.6 2.9 4.9 5.3 5.7 5.4 5.2
##  [631] 4.9 3.0 4.6 5.3 3.4 4.0 3.6 5.7 3.5 4.3 5.1 2.8 3.2 4.1 5.6 6.4 4.1 4.3
##  [649] 3.9 6.0 5.0 5.7 4.5 2.8 3.9 3.4 5.2 4.5 4.6 4.0 4.1 2.7 2.1 4.9 3.7 6.1
##  [667] 5.6 3.8 4.2 5.6 3.2 3.7 3.5 5.6 3.5 2.7 4.5 4.1 5.6 4.0 5.2 4.7 3.8 5.2
##  [685] 5.0 3.4 3.6 4.1 3.8 6.3 4.5 4.5 3.8 2.9 4.8 5.0 4.0 3.8 4.2 4.7 5.3 4.4
##  [703] 5.1 5.2 5.5 4.2 3.7 3.8 4.3 5.2 4.4 3.9 4.1 4.8 3.6 2.9 5.0 5.1 5.0 5.6
##  [721] 5.1 3.4 3.7 5.0 5.2 6.3 4.3 7.2 4.6 3.4 2.9 3.1 3.6 5.4 5.0 5.3 3.9 3.4
##  [739] 4.2 4.4 5.8 5.0 4.7 4.0 4.7 4.2 4.6 4.5 5.1 6.0 3.6 4.3 4.9 3.8 6.9 4.4
##  [757] 3.6 4.3 5.4 4.1 4.4 5.0 5.7 5.6 5.2 4.4 4.3 5.5 4.3 5.6 4.0 2.9 5.4 5.7
##  [775] 3.5 4.6 5.3 3.8 5.8 3.2 4.2 5.2 3.3 4.8 4.5 6.1 4.6 3.8 4.9 3.8 4.9 4.9
##  [793] 5.7 2.2 5.4 4.7 5.2 4.1 5.0 5.5 6.3 4.5 3.7 4.3 3.5 5.7 3.8 5.7 6.1 4.7
##  [811] 5.8 3.9 4.3 3.6 4.5 4.1 4.8 4.7 3.7 4.1 5.3 4.9 4.4 5.0 4.6 5.3 5.7 3.3
##  [829] 4.7 5.1 5.0 2.7 6.4 3.7 3.8 4.1 5.0 4.9 4.0 2.5 4.3 4.3 4.5 3.5 4.5 4.6
##  [847] 4.5 4.2 4.9 4.1 3.4 4.2 5.4 4.0 4.7 4.4 4.1 4.5 4.7 5.5 5.8 5.5 5.1 5.2
##  [865] 3.4 4.6 5.1 4.6 2.5 3.2 4.5 4.0 3.6 3.6 4.7 4.6 3.2 3.9 3.6 3.8 3.5 6.4
##  [883] 3.6 4.2 4.7 4.2 6.0 4.2 3.6 4.5 4.5 4.4 5.3 5.0 4.5 6.3 4.8 4.2 4.4 4.1
##  [901] 4.2 5.2 5.5 4.5 4.9 5.4 4.7 4.4 3.1 4.7 5.4 4.5 4.3 4.4 3.7 6.2 6.1 6.2
##  [919] 2.5 3.9 5.3 5.7 2.7 4.7 4.8 3.1 2.6 5.1 4.4 7.1 3.5 5.5 2.8 4.8 3.2 4.1
##  [937] 6.1 3.6 4.4 5.9 6.2 5.0 4.5 3.7 3.7 5.4 4.7 4.2 4.2 5.5 3.4 4.6 4.3 4.1
##  [955] 4.9 5.6 4.6 4.8 3.9 3.1 3.4 3.9 4.7 4.4 4.5 6.1 4.6 4.2 5.0 4.2 4.1 5.8
##  [973] 4.1 5.7 5.4 3.7 4.1 3.9 3.4 5.4 4.2 5.5 3.2 4.7 5.3 5.1 3.4 2.5 3.5 4.8
##  [991] 4.6 2.6 3.5 3.3 3.6 5.2 3.4 4.7 5.3 4.1
## 
## $func.thetastar
## [1] -0.0332
## 
## $jack.boot.val
##  [1]  0.48739496  0.34334365  0.31295181  0.13050847  0.02208589 -0.01661442
##  [7] -0.17595308 -0.31505376 -0.40817439 -0.53019391
## 
## $jack.boot.se
## [1] 0.9706578
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
##    [1] 4.0 4.4 4.0 4.8 2.9 5.3 5.7 3.0 6.0 5.1 3.6 4.0 4.3 5.5 5.4 5.2 4.3 4.8
##   [19] 4.9 3.8 3.1 4.2 2.7 3.5 4.5 4.5 4.2 3.8 5.0 4.0 5.5 5.2 4.0 3.3 3.9 3.9
##   [37] 3.3 4.7 3.5 4.9 5.7 4.8 4.0 4.8 5.0 6.3 4.0 3.8 4.3 4.1 4.0 4.8 4.7 3.5
##   [55] 5.6 4.3 2.8 4.1 3.7 6.0 5.7 4.5 4.6 5.7 6.6 4.3 5.1 5.1 3.9 3.3 5.6 4.8
##   [73] 5.5 4.4 4.3 3.0 5.5 4.0 4.0 4.0 3.3 4.3 5.1 5.3 3.6 6.2 4.1 4.9 4.6 4.4
##   [91] 6.2 2.5 4.3 3.5 5.2 4.6 5.5 4.0 4.2 5.5 5.1 3.2 6.0 3.9 2.6 3.7 5.4 4.8
##  [109] 4.3 4.0 4.6 5.6 4.8 6.3 3.5 3.6 4.9 5.6 6.4 2.8 4.2 3.5 5.1 4.6 4.7 3.3
##  [127] 4.1 3.4 5.2 2.9 3.5 4.1 4.8 3.9 5.0 3.4 3.6 5.1 3.4 5.7 3.8 4.8 4.5 4.8
##  [145] 3.9 5.1 4.8 5.1 4.4 5.2 3.3 2.9 4.1 5.7 3.9 4.8 5.1 4.9 4.1 4.8 5.3 4.0
##  [163] 5.1 4.2 5.2 3.7 4.4 4.6 4.1 4.1 5.8 4.1 4.9 3.1 6.5 5.0 6.4 5.1 5.4 5.5
##  [181] 5.3 4.1 5.7 4.2 5.5 4.4 5.5 2.8 6.2 4.8 5.9 4.9 5.1 6.0 4.7 4.9 3.7 6.8
##  [199] 3.3 5.0 3.6 3.6 5.9 1.8 3.0 4.2 4.8 4.1 4.5 7.1 3.1 4.2 4.8 4.3 6.6 3.9
##  [217] 3.3 3.7 5.1 4.9 4.8 4.9 4.5 4.4 3.4 5.4 4.8 6.0 4.3 5.0 3.5 5.1 4.3 4.0
##  [235] 3.1 3.3 5.5 4.2 3.9 3.8 4.5 3.5 5.6 4.4 4.0 4.0 3.8 3.7 2.6 3.9 5.7 4.0
##  [253] 4.8 5.7 4.4 2.9 4.0 3.5 4.2 3.4 5.7 5.8 2.1 5.1 2.9 4.0 3.5 5.2 3.8 3.6
##  [271] 4.4 5.3 2.4 5.0 3.8 4.7 4.6 4.8 5.9 4.2 3.0 3.6 4.5 4.3 5.5 5.0 5.3 3.7
##  [289] 4.5 4.6 4.1 2.9 3.2 4.6 3.9 6.4 4.7 4.7 4.9 5.2 3.7 4.6 4.5 4.0 3.5 6.9
##  [307] 3.7 3.3 3.0 5.2 4.5 4.4 4.6 4.5 5.1 3.8 2.8 5.8 3.6 5.5 5.8 5.5 3.0 3.1
##  [325] 4.7 3.1 4.7 3.9 3.7 5.2 5.2 4.5 4.1 4.4 4.2 5.3 5.1 3.8 5.0 4.5 5.8 3.7
##  [343] 4.8 4.6 4.9 3.2 4.1 2.6 3.9 4.4 4.7 6.6 7.2 4.4 5.9 3.3 5.1 4.1 4.3 3.7
##  [361] 3.7 3.9 4.5 3.7 3.0 3.9 5.9 5.0 3.5 5.7 4.3 4.2 3.2 5.0 6.1 5.4 4.1 4.0
##  [379] 4.5 4.2 5.8 5.0 4.8 4.3 3.8 4.6 4.7 5.0 5.6 5.1 4.4 3.0 4.2 3.6 4.3 4.0
##  [397] 5.3 4.4 4.3 4.6 2.7 4.7 2.4 5.3 5.4 3.2 2.6 4.1 3.6 5.1 5.8 4.7 2.8 5.0
##  [415] 4.3 4.9 5.7 4.9 4.7 4.2 4.8 5.9 3.3 4.4 5.7 5.8 4.2 4.9 3.5 4.6 5.8 4.7
##  [433] 4.0 5.2 4.3 3.5 4.0 4.5 3.0 3.9 4.7 3.2 6.6 4.8 5.9 3.2 3.3 3.6 3.6 5.2
##  [451] 4.1 4.0 5.4 4.8 4.8 4.0 3.8 3.5 3.6 6.0 3.8 4.1 3.9 4.3 4.3 4.8 5.1 5.4
##  [469] 3.9 4.6 5.0 5.1 4.9 3.9 5.7 3.9 3.7 5.4 4.7 7.0 3.7 4.2 6.7 3.8 4.0 5.3
##  [487] 3.8 4.8 4.8 5.1 3.2 3.8 3.5 4.0 4.1 4.0 5.3 4.1 5.1 4.3 3.8 5.9 4.6 4.0
##  [505] 5.0 4.2 3.8 3.7 5.1 3.2 5.2 5.0 5.2 4.6 3.8 4.7 4.5 5.9 4.1 4.9 3.0 4.8
##  [523] 4.4 4.7 3.6 3.4 5.4 4.5 4.1 4.9 2.7 3.7 3.8 2.6 6.0 4.3 3.5 6.1 5.2 4.0
##  [541] 5.8 4.8 4.2 4.7 5.3 4.1 2.7 3.6 5.0 4.2 5.0 5.7 3.4 4.3 2.7 2.6 3.6 5.7
##  [559] 6.7 5.3 5.2 4.0 3.9 3.9 4.3 2.8 7.6 4.7 5.8 4.8 5.7 4.1 4.5 4.2 3.9 5.3
##  [577] 3.8 2.8 4.0 4.7 3.1 4.9 4.4 3.8 4.9 2.1 4.4 4.1 4.1 5.3 4.4 3.4 5.8 6.0
##  [595] 4.7 3.2 5.4 3.9 3.9 2.4 4.2 4.4 5.1 4.0 4.6 3.9 4.9 5.6 3.7 3.1 4.8 3.8
##  [613] 4.2 3.8 3.5 2.1 2.0 5.2 3.9 6.2 5.5 4.1 2.7 5.1 4.8 4.1 5.6 5.8 2.5 5.3
##  [631] 3.5 5.4 3.4 4.5 4.0 3.3 5.5 6.3 4.8 4.3 5.0 4.3 3.2 4.0 5.0 4.9 4.5 5.0
##  [649] 5.3 5.0 5.1 4.6 6.0 4.8 4.8 3.4 3.3 2.6 4.9 5.5 4.7 3.0 6.6 2.5 5.0 4.4
##  [667] 5.8 4.3 4.2 4.4 3.8 4.7 3.3 4.4 5.2 4.4 4.1 5.5 4.6 3.8 3.7 4.5 4.0 4.1
##  [685] 5.7 3.7 5.4 2.9 4.7 5.3 3.8 4.1 5.8 4.0 5.0 2.9 3.1 3.8 5.6 3.8 5.4 4.7
##  [703] 4.4 4.5 4.4 4.3 5.8 4.6 3.4 4.8 5.1 4.6 4.4 4.0 6.3 4.5 5.0 3.6 3.4 3.5
##  [721] 4.4 4.4 3.6 5.5 5.5 5.0 5.2 5.2 4.1 4.9 5.2 5.9 4.2 3.9 5.1 4.0 3.4 3.8
##  [739] 4.4 1.9 4.6 4.5 4.9 4.9 4.7 2.9 5.4 4.6 3.5 6.0 3.7 3.9 4.8 2.9 4.1 3.8
##  [757] 5.2 4.5 4.3 5.3 3.7 3.8 3.2 6.0 5.7 4.3 3.8 3.7 4.0 4.6 5.1 4.0 4.2 4.4
##  [775] 5.2 3.6 4.9 3.7 4.6 3.9 4.8 3.5 4.9 4.7 4.8 4.3 5.0 4.3 4.6 4.5 3.1 2.7
##  [793] 4.2 3.6 5.2 4.6 5.2 5.6 4.1 4.0 4.3 4.7 4.4 4.9 5.1 5.1 5.3 3.2 4.5 5.5
##  [811] 4.7 4.7 2.8 5.1 4.7 6.3 5.0 5.5 3.6 4.3 5.8 4.3 4.9 3.4 4.6 3.0 3.8 4.0
##  [829] 4.8 3.5 4.8 3.3 5.6 4.2 5.2 3.7 4.3 5.0 2.7 6.2 3.9 5.0 6.2 3.3 3.7 5.3
##  [847] 4.7 3.7 4.5 6.7 6.4 5.4 4.1 5.3 4.6 4.4 4.1 4.8 6.7 4.8 4.2 2.6 5.2 4.0
##  [865] 3.5 4.5 3.8 4.3 3.9 3.7 3.2 4.1 4.7 3.0 4.2 4.5 4.9 3.4 3.2 5.4 5.9 4.5
##  [883] 3.6 6.0 4.1 4.0 5.5 4.6 4.8 4.6 5.3 4.4 4.7 5.0 4.9 4.7 4.5 5.4 6.5 4.7
##  [901] 4.4 4.5 4.3 3.7 4.1 5.2 3.9 5.3 4.0 3.6 2.8 4.8 4.4 2.8 5.7 5.1 4.5 4.3
##  [919] 4.5 4.9 4.1 2.8 4.1 5.4 4.8 5.3 3.3 5.8 5.1 4.6 4.0 4.9 3.9 4.6 4.4 4.0
##  [937] 2.7 4.1 2.9 3.2 5.8 5.5 3.7 4.6 4.5 4.1 3.9 3.2 5.9 4.8 5.6 2.8 4.6 3.4
##  [955] 4.4 5.1 4.6 5.1 4.4 3.9 3.8 3.4 3.9 4.9 3.0 6.3 4.2 5.0 2.8 5.2 3.6 4.6
##  [973] 5.2 4.4 3.4 4.4 4.1 5.2 5.3 4.2 4.9 3.5 5.0 3.4 5.4 3.7 3.1 4.2 4.0 4.3
##  [991] 5.4 5.0 6.0 6.0 4.4 4.4 5.8 5.4 5.0 3.6
## 
## $func.thetastar
## 72% 
##   5 
## 
## $jack.boot.val
##  [1] 5.5 5.3 5.3 5.2 5.0 4.9 4.8 4.7 4.5 4.4
## 
## $jack.boot.se
## [1] 1.040961
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
## [1] 0.480138
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
##   3.463087   5.333432 
##  (1.480212) (2.453196)
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
## [1]  1.13257088  0.74317030 -0.03207377  1.30014174  0.19473496  1.82777550
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
##    [1] -0.0974663068  0.7912440668  0.2476980168  0.7728846956 -0.0264795898
##    [6]  0.1131262327 -0.0287531353  0.0198908254  0.8015438015  1.0662884381
##   [11]  0.2919744351  0.4525437846  0.3534457257  0.4853202213  0.7158455128
##   [16]  0.3658881226  1.3424960696  0.4524821552  0.5315581972 -0.1611922154
##   [21]  0.7191036823  0.2244967255  0.3344593775  1.3578605737  0.3384480433
##   [26]  0.1471438641  0.3264148888  0.0868480277  0.1972079488 -0.3255091975
##   [31]  0.4468687505  0.8077300916 -0.0252974057 -0.0122141067  0.0399317508
##   [36]  1.3025710111  0.2137995149  1.4065596744 -0.1330874600  0.0321754817
##   [41]  0.3950199248  0.0824177267 -0.2083547751  0.4707863365 -0.1018473644
##   [46]  1.1259018658  0.0425281595  0.0900681726  0.6526896724  0.8533227395
##   [51]  0.4658061296  0.2911820178  0.8146322272  1.1079612633  0.2999552384
##   [56]  0.5377453046  0.7364000469  0.1838177188  0.1231427267  0.8181966816
##   [61]  0.8564724825  0.6770084599  1.0221481087 -0.0391201796 -0.3083445221
##   [66]  0.3784378043  0.8290318669 -0.1821928175  0.4866236643 -0.1014367037
##   [71]  0.1622600854  0.5144748273 -0.4612539710  0.2919744351  0.3181884509
##   [76]  0.3609561674  0.8419285178  1.9074030974  0.7079362827  0.8172585739
##   [81]  0.4985257047  1.5097490888  0.5330482551  0.0590671780  0.2445252236
##   [86]  0.1781998045  0.1683682681  0.0474612132  0.7530154559  0.2589251525
##   [91]  0.7078659195  0.1290389041  0.8741684099 -0.3349799142  0.8029059966
##   [96]  0.0693171015  0.3730270405  0.5720138597 -0.0402452963  0.4360581598
##  [101]  0.6775603805  0.6535821737  0.8576059872  1.1715138338 -0.2302374385
##  [106] -0.3969969063 -0.1823529147  0.5682810568  0.9010484144 -0.0776770959
##  [111]  0.7988031240  0.6180154921  0.4163065657  0.8801987669  0.5316066559
##  [116] -0.0574452779  0.3572151561  1.1246282663  1.4295849908  0.3849852240
##  [121]  1.9068334786  0.6677149625  1.8993840511  0.1780857365  0.2631957079
##  [126]  0.1542773327  0.7543141488 -0.0142730044  0.7870151999 -0.0959746405
##  [131]  1.2541687093 -0.0704382160  0.7526961873  0.2959714444 -0.0210122983
##  [136]  0.7823115967  0.4964970197  0.1519644308  0.0354238292  0.4680422370
##  [141]  0.9580947471  0.3319168933  0.4981563236 -0.0670342854  0.1823954118
##  [146] -0.1012071872  0.9870890119  0.2883546625  0.3441914948  0.2021189369
##  [151]  0.4211149403  1.0524715001  0.1031695255 -0.4429179797  0.7225673184
##  [156]  0.2657649235  0.2112307932  0.8357963285  0.4131830498  0.1330150327
##  [161]  0.1518078475  0.4934244790  0.1056254401 -0.0653439021  0.0917512746
##  [166]  0.9726883986  0.3101404768  0.2814979114  0.3205743196  0.1298251389
##  [171]  0.0405753356  0.8153022012  0.2195868241 -0.1682125810  0.4777214623
##  [176]  1.4726683215 -0.7018032477  0.1194010686 -0.0406413264  0.0066582513
##  [181]  0.7064294376 -0.2953858652  0.1742097655  0.9978313207 -0.0649218951
##  [186]  0.5947021184  0.0303781894  0.5023341762  0.4887958008  0.0447930380
##  [191]  1.7690603159  0.3085295041  0.7419929857  0.3174871999  1.0702577872
##  [196]  1.0901545993  0.6817049722 -0.2845060801  0.1643352902  0.6137723326
##  [201]  1.1958285690  0.3059602925  1.5849102313  0.3697964235 -0.0755445378
##  [206]  0.2702845749 -0.1724393949  1.1216509617  0.8787328978  0.6625132965
##  [211]  1.2092920998  1.2203747369  0.4057742481  0.8182577588  1.0154454178
##  [216]  0.8534129358  0.3450469958  0.5294251304  0.3848529140  0.7868143759
##  [221] -0.8811618503 -0.2256969946  0.8784981854  0.3609362779  0.3891281964
##  [226]  1.3074460465  0.4291406988  0.6788547344  0.0712269688  0.4069435249
##  [231]  0.3010173594  0.4997987074  1.4614827063 -0.6540692977  1.8329269576
##  [236]  0.2480965608  0.5289481924  0.3414423777  0.5962226795  0.7661870263
##  [241]  0.9134194194  0.3308157014  0.8051945319  0.7348305298  0.0198908254
##  [246]  0.0617538877  0.2637966015  0.6554255842  0.5951928100  0.1583156043
##  [251]  0.7428854714 -0.0885641326 -0.0811317223 -0.0445283464  0.6754845411
##  [256] -0.5651984949  0.2128044438  0.5768507531  0.7032959805  1.4039747639
##  [261] -0.1055991076  0.5384506689  0.0837933688  0.6385121269  0.8375484850
##  [266] -0.4312666128  0.7351983766  0.4000973145  0.6228668243 -0.4755644771
##  [271] -0.1359848282  0.7868079866  0.5098647102 -0.2598002447  0.1581284957
##  [276] -0.1638370072  0.4177738198  0.6739684970  0.8571492658  0.7614047076
##  [281]  0.4726897756  0.5885260357  0.3142845689  0.0980162317  0.8857226682
##  [286]  0.5044458524  0.5468832199  0.2226090630  0.2146922511  0.2160920048
##  [291]  0.0591795289  0.0714634103  0.2978339070  0.6800398191  0.3407164604
##  [296]  0.3886003065  0.7240906681  0.2779519076  1.3249898829  0.6810687924
##  [301]  0.1891774299  0.4428332532  0.1809851977  0.5806729791  0.1770306153
##  [306]  0.3387193214  0.5872650950  0.8098860257  0.3073784499  0.4234287673
##  [311]  0.9563803398  0.7264646172  0.4982839786  0.2675748239  0.4841705997
##  [316]  0.6703669738  1.2086449832  0.6785604463 -0.5440773048  0.4753957130
##  [321] -0.0832634895  0.5229737129  0.4069401694  0.4076363170 -0.5560484015
##  [326]  0.6591464236  0.4626057280  0.3215284257  0.1496948213  0.7888737656
##  [331]  0.6239740178 -0.3035207952  0.6543939414 -0.0812969942  0.8279698263
##  [336]  0.5621461308 -0.1998025688  0.2080631438 -0.7417411674  1.1362863137
##  [341] -0.1779869308  0.7526803791  1.3633821529  0.0714634103  1.0592216395
##  [346]  0.3738926554  0.8144693256  0.1731513310  0.7178693494  0.2114930479
##  [351]  0.3895289168  0.6119805486 -0.1201309981  1.3147841842  0.0985433811
##  [356]  0.1024474873  0.9432562291  0.3784056443  0.0951760818  0.6722119382
##  [361]  0.5838591672 -0.1769179624 -0.4013556732 -0.0511613834  0.6224777769
##  [366]  0.8063508807  0.1028206688  1.0380365572  0.1636611363  0.3919359102
##  [371]  0.6139220524  0.6446622627  0.3176955913  0.8674723886 -0.2659547779
##  [376]  0.4141856244  1.8321653153  1.0846452978  0.4935227855 -0.3866550865
##  [381]  0.4943232822  0.9877507694  0.7023419544 -0.0827347991  0.1418781450
##  [386]  0.5481611399  0.2585567856  0.3864011522  0.6656219511  0.9050034603
##  [391] -0.1843289096  0.4952023537  0.3685887543 -0.1483719004  0.3904143875
##  [396]  1.2958168967  0.3549214975 -0.3066235697  1.0851120989  0.4260618904
##  [401] -0.0663828303 -0.0601223645  0.5020728893 -0.0752383691  0.8458439389
##  [406]  0.8912271178  0.3595900842  1.3242302287  0.4465862074  0.1044291664
##  [411]  0.0122267945 -0.1590012350  0.4935869924  0.7223675356 -0.2029727654
##  [416]  0.4184162659  0.4669909749  1.4844791481 -0.0239950337  0.5780785174
##  [421]  0.3581064383  1.2127959723 -0.4745759320 -0.1660669877 -0.2539319423
##  [426]  0.4757223812  0.1036368283  1.1822568036 -0.3315650448  1.0756930545
##  [431]  0.7373359399 -0.3999211768  0.2501317519  0.9868028309  0.3977127756
##  [436]  0.7437288782 -0.7519431744  0.4673977780  0.5744542347  0.0455115055
##  [441]  0.3844432408  0.6230697705  0.6608226013  0.3956224060  0.1778485783
##  [446]  0.8147695724 -0.0624020678  0.5315020354  0.8846398462  0.6174307506
##  [451]  0.6191461438  0.2584557959  0.3967190898  0.9148598620  1.0521947276
##  [456] -0.1900073041  0.4518809031  1.0111984150  0.0646930339  1.3560662420
##  [461]  0.4511519943  2.0249765558  1.0653695428  0.3663008900  0.7091195725
##  [466] -0.0491704240  1.3186856156  1.0766970376  0.4827642624  1.1508220576
##  [471]  0.3745213850  0.1317138412  0.4584030373  1.1637682324  0.1297397698
##  [476]  0.1783936558 -0.0475615222  1.1365373328  0.3573415492  0.5256816887
##  [481]  1.4924882980  1.8133020458  0.9461775079  0.7323439611  0.4672781515
##  [486] -0.0244053906  0.2072384779 -0.1206881481 -0.0445406573  0.5164162497
##  [491]  0.3532278713 -0.2626901778  0.0352266117 -0.1436300630  0.5514816210
##  [496] -0.0391242159  0.4753774287  0.4562999960  0.8274315305  0.3208633620
##  [501]  0.6301794408  0.4100217741  0.7319772367 -0.0424332967  0.5059545996
##  [506]  0.7509110307  0.6106281328  0.2991140348  1.7347315651  0.6495347126
##  [511]  0.7241343704 -0.5690787879  0.5052995692  1.0915098529  0.7362466735
##  [516]  0.0259604533 -0.3958325265  0.6757706645  0.5807118551  0.4140351961
##  [521]  0.7542768246  0.6380425387 -0.0748841887  1.1305941003  0.1786482878
##  [526]  0.1239394366  0.3767749895 -0.0030381794  0.5077732483  0.0357728666
##  [531] -0.4156548223  0.6736794711 -0.2605762321  0.5806729791 -0.1165748299
##  [536]  0.6959586518  0.3653639596  0.3041867684  1.2125373671  1.0242174466
##  [541]  0.9805232324  0.9453813741  0.3512685210  1.1602857189  0.6490967222
##  [546] -0.2160046722  1.1684110624  0.3296768154  0.4307850989 -0.3460955252
##  [551]  0.3781182064  0.1968398105  0.8298642081  0.5494128382  0.5256718277
##  [556] -0.1587907685  0.3039317890 -0.2498233668  0.3589659401  0.4320012306
##  [561] -0.3429873339  0.2996521465 -0.3418562322  0.5192434822  1.1536186592
##  [566]  0.5746678791 -0.0550140415  1.1186155292 -0.3492985701  1.1564505144
##  [571]  0.7980247116  1.1928942348  0.1020774450  1.0448738793  0.9121513038
##  [576]  0.4750172130 -0.0007216855  0.6050501947  0.3795467960  1.3439145921
##  [581]  0.6394620697  0.4115573976  0.5619013743  1.2859137556  0.7833696704
##  [586] -0.0674254728  0.4705569964  0.6587151172  0.4913004018  0.4036631108
##  [591] -0.3799207590 -0.3977032287  0.9643899514  0.4368209325 -0.3068867098
##  [596]  0.1824173303  0.4098557865  0.5873072523  0.5150612570  0.9675692989
##  [601] -0.0599989498  0.2290638653 -0.1031130625  0.2409093826 -0.0208066908
##  [606]  0.2438647252  0.2647364924  1.2020438732  0.2495409534  0.7693520989
##  [611]  1.3080414126  0.1075942003  0.0663149048  0.7526366229  1.3218009643
##  [616]  0.7749908786  0.5481610802  1.5232043045  2.0217601267  0.3447290575
##  [621]  0.7730064381  0.5071425815  0.8860848662  0.7373618378  0.1067743073
##  [626]  0.8266078516  1.4343154143 -0.0484932682 -0.0108797813  0.4271184973
##  [631] -0.3454373675 -0.0001120349  0.5953520067  0.4271928221  0.4896585351
##  [636]  1.1112838624  0.6776325691  0.6487049691 -0.4761412683  0.1644317707
##  [641]  0.6144521090  0.2382894268  0.2717490698  0.5324984244  0.7139100229
##  [646]  0.4253558731  1.1841182321  0.4872864133  0.5965185248  0.9506857869
##  [651]  0.3612464131  0.9771148475  0.5036397565 -0.1472714129 -0.1504190892
##  [656]  0.6871430693  0.0316475351  0.6656219511 -0.2946755152  0.5274183060
##  [661]  0.6659584939 -0.5688204367  0.7630928952  0.9577414760 -0.1574733775
##  [666]  0.8303291675  0.0445547576  1.0440576131 -0.1946145094 -0.3171654653
##  [671]  0.4554892940  0.1626500243  0.2885839690  0.1466292112  0.8056812372
##  [676]  0.7965640537  0.3374748419  0.3626505362  0.0813168724  0.1856963026
##  [681]  0.0105639478  0.8120773319  0.2603171841 -0.0742639811  0.9266384078
##  [686]  0.4316361804  0.3502809769 -0.2352020159  1.3894825310 -0.4487799757
##  [691]  0.3543777017  1.1363496210  0.8244604722  0.5839654034 -0.0503746796
##  [696] -0.1822047190  0.6414527507  0.7160756837 -0.1546417714  0.3642100377
##  [701]  0.6501214031  0.3833690936  0.1982006305  0.5990643716  0.6724623294
##  [706]  0.9979088094  0.3461281508  1.4819394021  0.9445956863  0.6239671478
##  [711]  0.4972504959  0.1827047867 -0.1412170350  0.7624057836  0.6858381626
##  [716]  0.5944610240  0.9687323397 -0.1728901937  0.2517869683 -0.2687183584
##  [721]  0.5330353682  0.2181945508  0.3683973702  0.5479978712 -0.2032069378
##  [726]  0.3275005126 -0.1411636809  0.4801984808  0.1564094602  0.0617538877
##  [731]  0.3919310641  0.4342120350  0.5276368962  0.0127515071 -0.3133442833
##  [736] -0.4599980289  0.1837061961  0.5201204456 -0.0229992979  0.3038076023
##  [741]  0.2133651965  0.3276103621 -0.1691741649  0.8378755867  0.1262694449
##  [746]  0.6016327926  0.4388328712 -0.2092814761  0.5326872559  0.4467081692
##  [751]  0.5385143842  0.8248222111  0.3627222144  0.5587028703  1.3148577269
##  [756]  0.4013272747  0.6134658446  0.2198004263  0.6329903343  0.2927204480
##  [761]  1.1778060392  0.3160428900  1.2807583674 -0.0437611395  0.4658264970
##  [766]  0.2949556458 -0.7772678484  0.8677816556  0.8469115682  0.7831012771
##  [771]  0.3484917882  1.1176717265  0.3852604241  0.5012708459  0.5187989556
##  [776] -0.1882134127  0.3188628824  0.2946101710  0.7759057573  0.7090483108
##  [781]  0.9474217368  0.3656556558  0.3654159568  0.0731221406  0.4766553167
##  [786]  0.6955189249  0.0368533879  0.9019074821  1.0036594844 -0.0613980692
##  [791]  0.8351569243  0.3061133108  0.1656339638  1.1087513407  0.6209606977
##  [796]  0.0155062514  0.7960012164  1.2216415718  0.5642279022  0.6871218557
##  [801]  0.4487982827  0.3620859267  0.5091609394  0.2006671339  0.4936960195
##  [806] -0.0316729197  1.0161020695  0.9841291148  0.0419300736  0.0835800290
##  [811]  0.2576733623  0.3151002559  0.6115161119  0.2760949391 -0.2558056194
##  [816] -0.0574688783  0.5520609253  1.5724783538  0.2109315034  0.3264148888
##  [821]  0.4146708000  1.1828990081  0.6938330921  0.5913800098  0.7777330658
##  [826]  0.8053682119  1.7132857784 -0.1376230035  0.2999359407  0.0811935381
##  [831]  0.5596535502 -0.0118568475  0.9776255459  1.4305452393  0.4962915671
##  [836]  0.4504529126 -0.2676748987  0.3924395967  0.3774882395  0.5479312070
##  [841]  0.1337723072  0.2469175241  1.3590923163  0.0961845990  0.5717905100
##  [846]  0.9418415001  0.1044291664  1.2692752076 -0.3449314385  0.5061600041
##  [851]  0.9248097607  0.3694604711  0.0002208574 -0.1383117147  1.2913198953
##  [856]  0.7136728055  0.5349583427  0.3921278901  1.7652102280  1.1747454097
##  [861]  0.8854179170 -0.2906379587  0.6145573089  0.4731419848  0.4068213813
##  [866]  0.0333883695 -0.0344337153  0.3044787406 -0.6047787888  1.0865119579
##  [871]  0.6723387656  0.2355223808  1.0682175320 -0.1378309249  0.1351777169
##  [876]  0.3739027428  0.1131393043  0.4800956170  1.0102456170  0.2790706002
##  [881]  0.6363831779  0.3782165589  0.7353429472  0.4130000130  0.9622179746
##  [886]  0.3297781666  0.4049740557  1.0341906256  0.2311274952  1.0092950370
##  [891]  0.7021009486  0.3206488915  0.9657396892  0.2804982945  0.3785642366
##  [896]  0.0466933930  0.9358207575  0.5720883015  0.3885298148 -0.0973874473
##  [901] -0.3699469626  0.2413404766  0.3413728637  0.0610810805  0.0617724030
##  [906]  1.8436544252  0.7719173044  0.3077776091  0.8039278070  0.9446926440
##  [911]  0.3755455246  0.9790769405  0.3359813191  0.7749908786  0.4994758400
##  [916]  0.1844344866  0.2586712676  0.8676156834  0.9282430804  0.2390358920
##  [921]  0.4022540663  1.0839592389  0.3528358263  0.6518468579  1.0640507982
##  [926]  0.7008112125  0.7769347048  0.1195235692  0.3718199822  0.7449240059
##  [931]  0.4257993421  0.8268061961  0.7112131937  0.2103916634  0.4297394398
##  [936]  0.4748369464  0.3646385822  0.1847742751  0.6054589306  0.1872759238
##  [941]  0.7646038239  0.5069164676  0.5412805994  0.0264429540  1.1647589307
##  [946]  1.1341996423 -0.3595247644 -0.2482666895  0.0857349957  0.6497824694
##  [951] -0.0797954264  0.3129320003 -0.2539319423  1.2042716753 -0.2498233668
##  [956]  0.7679224805  0.2889786901  0.5687477615  0.7123285638  1.1509366784
##  [961]  0.2135237515  0.7028115810  0.1015403851  1.3071687728  0.2124569546
##  [966] -0.0383183986  0.4678279401  0.6474584339  0.4122984825  0.0799386393
##  [971]  0.5257413739  1.3375994198  0.8090125234  1.5236271002  0.6023422749
##  [976]  0.3327297893  0.1239114328  0.3814858719  0.7001924503  0.0223561503
##  [981]  0.9119495223 -0.1371574274  0.1029486606 -0.1558870872  1.4928656178
##  [986]  0.4472666198  0.8963045198 -0.0872526544  0.8532643961  0.4593096680
##  [991]  1.0600596208 -0.3062574107  1.3233159058  0.2039436072  1.1237068366
##  [996]  0.1823954118  0.3350669907  0.5455488042  0.3155002688 -0.6124976073
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
##   0.64930197   0.33974168 
##  (0.10743575) (0.07596863)
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
## [1] -0.14142598  0.58653826  0.24147043  0.07511422 -0.61669404 -1.64230114
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
## [1] 0.0092
```

```r
sd(results2$thetastar)  #the standard deviation of the theta stars is the SE of the statistic (in this case, the mean)
```

```
## [1] 0.8995914
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
## t1*      4.5 0.004104104   0.8865919
```

One of the main advantages to the 'boot' package over the 'bootstrap' package is the nicer formatting of the output.

Going back to our original code, lets see how we could reproduce all of these numbers:


```r
table(sample(x,size=length(x),replace=T))
```

```
## 
## 0 2 5 6 7 9 
## 2 1 1 1 2 3
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
## [1] -0.0201
```

```r
se.boot
```

```
## [1] 0.8844253
```

Why do our numbers not agree exactly with those of the boot package? This is because our estimates of bias and standard error are just estimates, and they carry with them their own uncertainties. That is one of the reasons we might bother doing jackknife-after-bootstrap.

The 'boot' package has a LOT of functionality. If we have time, we will come back to some of these more complex functions later in the semester as we cover topics like regression and glm.

