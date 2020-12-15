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
## 0 2 4 6 8 9 
## 2 3 1 2 1 1
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
## [1] 0.0361
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
## [1] 2.83849
```

```r
UL.boot
```

```
## [1] 6.23371
```

Method #2: Simply take the quantiles of the bootstrap statistics


```r
quantile(xmeans,c(0.025,0.975))
```

```
##  2.5% 97.5% 
##   2.9   6.2
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
##    [1] 5.1 6.6 4.4 4.4 4.6 4.3 5.8 6.5 5.2 3.3 5.8 3.8 4.4 4.0 5.3 4.4 4.2 4.1
##   [19] 3.2 3.5 3.7 4.3 4.3 4.4 5.4 5.6 5.3 5.4 4.2 4.9 3.2 4.8 4.3 5.3 4.9 4.2
##   [37] 7.1 4.3 4.5 4.1 4.4 3.6 4.1 4.5 4.0 5.1 5.5 4.4 4.8 5.8 4.3 5.3 3.9 3.7
##   [55] 3.0 3.1 2.8 5.1 2.9 5.5 4.7 4.6 5.6 4.3 5.0 5.3 4.4 5.2 5.1 5.7 4.0 3.6
##   [73] 3.8 3.9 5.9 4.6 4.1 5.6 5.1 4.4 5.5 4.9 5.3 2.9 4.6 6.0 3.8 2.7 5.0 4.5
##   [91] 4.1 5.6 4.8 4.7 5.4 3.9 4.4 4.6 4.5 5.0 4.4 4.4 5.7 4.7 5.3 4.6 3.8 5.2
##  [109] 4.5 4.1 3.2 3.6 4.2 5.3 4.6 6.7 5.2 3.0 3.7 4.1 3.9 4.1 2.6 3.6 5.3 4.2
##  [127] 4.6 4.0 6.7 5.4 6.3 3.4 5.6 5.0 4.6 4.6 3.9 6.0 3.1 5.0 3.5 5.7 3.9 5.2
##  [145] 4.2 3.5 3.4 5.0 3.6 4.6 5.0 4.6 4.6 4.1 3.6 4.2 3.3 5.1 4.6 4.1 4.2 3.8
##  [163] 4.6 4.5 5.4 5.2 3.6 2.7 3.9 2.8 4.1 4.9 4.3 4.4 5.4 4.0 5.1 4.6 3.9 4.7
##  [181] 4.1 5.0 4.3 3.4 6.5 6.1 5.8 4.5 4.7 6.1 6.6 4.4 5.9 3.2 4.2 4.5 4.4 4.7
##  [199] 5.5 4.2 5.3 3.9 3.5 3.1 4.3 5.2 3.7 4.3 5.0 2.9 5.6 2.8 4.8 5.2 5.7 3.1
##  [217] 5.3 5.4 5.5 3.2 3.7 4.2 4.0 4.8 3.3 5.4 5.4 4.6 5.2 4.9 4.1 5.5 5.4 4.4
##  [235] 3.4 5.9 5.4 4.8 4.4 4.0 4.6 4.7 5.6 2.8 6.7 4.0 3.4 4.7 4.4 2.7 4.7 4.3
##  [253] 5.2 4.2 4.6 4.5 2.8 5.1 2.4 5.3 5.6 5.5 6.0 3.7 3.0 3.0 4.9 4.8 5.9 4.7
##  [271] 4.1 4.5 4.7 3.5 3.9 4.9 4.5 3.2 3.8 4.3 3.3 4.2 5.6 3.1 5.1 4.1 5.4 5.3
##  [289] 4.6 5.8 5.1 4.4 4.9 3.3 3.8 4.0 3.6 3.6 3.8 4.9 5.5 3.2 4.2 5.2 5.6 3.9
##  [307] 4.2 4.3 5.3 2.4 4.4 5.6 4.6 5.4 4.9 4.2 3.6 3.5 4.9 4.6 3.7 4.1 4.3 3.7
##  [325] 6.0 4.8 5.2 4.5 3.8 5.5 5.0 4.6 4.8 3.6 4.6 5.3 3.8 5.3 4.8 4.4 5.4 4.5
##  [343] 4.5 4.8 3.7 3.3 3.0 2.9 5.2 5.6 4.1 4.4 4.6 4.3 3.7 2.5 3.4 4.3 4.9 3.7
##  [361] 4.1 5.2 4.4 3.8 4.3 4.8 6.4 4.7 4.1 2.6 5.3 3.6 6.3 4.0 5.5 5.8 4.5 5.8
##  [379] 6.2 4.4 3.9 5.3 2.9 4.0 4.4 4.9 4.1 3.8 5.0 5.4 4.0 5.3 4.6 4.1 4.0 6.0
##  [397] 4.2 3.7 5.3 4.4 5.0 5.4 2.6 6.5 3.2 4.0 5.4 4.2 3.0 3.3 5.2 3.1 2.6 3.1
##  [415] 3.6 4.7 4.7 3.6 4.0 4.0 4.5 4.1 5.7 3.4 4.6 4.7 3.3 4.9 5.2 4.0 5.5 3.6
##  [433] 6.1 4.9 3.7 4.5 5.9 4.8 5.7 4.4 5.4 3.3 4.5 5.3 4.4 4.8 4.1 4.8 4.1 5.6
##  [451] 5.4 4.3 5.3 4.1 3.7 3.2 3.9 4.7 5.3 5.3 4.6 3.5 4.1 4.2 4.1 3.6 3.9 4.6
##  [469] 3.8 4.9 5.7 4.7 3.7 4.4 5.0 3.8 2.4 3.5 3.5 3.6 4.9 4.3 5.6 5.0 3.6 4.6
##  [487] 5.5 4.3 5.0 3.5 5.1 4.1 5.8 3.1 5.2 4.9 4.6 2.8 4.0 4.0 4.8 5.0 4.7 4.3
##  [505] 5.5 4.0 5.0 4.9 4.1 5.6 2.9 4.5 3.5 4.6 6.0 3.8 5.9 3.6 3.9 4.2 4.2 3.9
##  [523] 3.6 5.2 3.8 4.8 3.9 4.4 5.7 5.0 3.1 4.5 3.0 6.0 4.6 4.0 5.5 3.5 5.4 4.9
##  [541] 3.5 4.5 3.4 3.7 5.5 4.1 5.2 4.9 4.2 4.5 4.1 2.5 4.5 3.6 4.5 2.1 4.2 4.0
##  [559] 4.3 3.7 5.0 5.0 4.5 5.3 4.4 4.6 4.9 4.0 3.1 4.5 6.0 4.8 4.8 4.0 3.3 5.0
##  [577] 3.2 3.6 4.7 3.4 4.5 5.6 6.0 5.7 3.3 4.1 3.9 3.8 6.6 3.4 5.3 4.4 4.0 3.5
##  [595] 4.6 3.2 3.8 4.7 4.8 3.9 4.6 3.2 4.9 6.8 6.1 4.9 4.2 6.0 3.4 4.1 4.6 4.1
##  [613] 4.3 4.1 5.4 6.0 3.6 3.5 6.0 4.3 4.2 4.5 4.9 4.1 3.8 5.1 3.7 3.0 4.1 4.0
##  [631] 4.5 3.5 3.5 4.7 5.3 3.9 3.8 5.9 6.1 4.3 5.1 5.8 3.1 5.5 2.2 5.2 4.0 4.4
##  [649] 4.9 4.4 4.3 3.6 5.6 5.7 4.2 3.0 3.8 5.1 4.9 4.9 3.8 5.1 3.2 4.7 5.2 5.2
##  [667] 5.2 5.3 4.5 5.4 2.4 4.1 4.9 4.0 5.4 4.0 5.9 3.5 4.6 5.1 5.3 4.6 4.2 3.3
##  [685] 3.1 4.1 4.7 4.7 5.2 4.7 5.0 4.4 5.3 7.2 4.4 4.2 4.3 4.8 4.5 4.6 2.5 5.5
##  [703] 5.0 5.2 4.9 3.1 6.3 5.4 6.1 3.9 4.7 3.9 4.6 5.3 6.2 3.6 3.9 4.5 5.4 4.8
##  [721] 4.9 4.8 5.5 5.1 5.6 5.7 5.4 4.0 5.0 4.2 3.7 4.7 4.6 5.2 5.1 5.0 5.4 4.7
##  [739] 3.0 5.4 2.8 4.2 3.1 4.4 3.2 4.4 4.7 4.6 3.2 3.3 3.7 6.1 3.7 3.4 4.3 3.9
##  [757] 4.6 3.4 4.7 3.7 5.4 5.5 5.0 3.7 3.9 4.3 4.5 4.8 3.4 5.9 4.0 3.1 6.0 3.5
##  [775] 3.8 4.0 2.7 4.6 5.0 5.4 4.9 3.3 3.5 5.7 5.1 5.1 4.6 3.9 3.3 3.9 3.0 4.1
##  [793] 3.8 2.3 3.7 4.4 4.4 5.2 3.3 3.2 2.9 5.8 4.1 3.7 5.9 4.5 3.4 5.0 3.5 4.7
##  [811] 4.4 5.0 4.1 5.3 5.7 4.0 5.1 6.7 4.2 6.2 4.0 3.9 3.6 4.1 4.8 4.1 5.6 4.2
##  [829] 4.3 3.0 5.0 4.5 4.6 3.9 4.9 5.3 5.7 4.3 4.6 4.3 5.1 3.7 4.3 5.3 4.1 4.7
##  [847] 4.5 5.5 5.0 4.7 4.5 2.8 4.2 4.1 5.6 2.4 3.5 4.9 3.4 4.5 4.4 4.8 5.7 3.1
##  [865] 5.6 4.8 4.7 5.2 5.1 1.8 3.0 5.3 2.9 6.4 3.5 4.8 4.9 4.2 4.1 3.0 6.2 5.8
##  [883] 4.7 4.1 4.8 4.8 4.2 4.8 5.5 4.2 5.8 3.4 4.5 4.8 3.2 4.5 4.5 5.7 3.8 3.8
##  [901] 3.5 5.3 4.3 6.3 4.6 5.2 3.2 6.2 4.3 4.4 6.0 5.1 6.6 4.7 4.0 4.5 4.5 4.9
##  [919] 4.4 4.1 4.3 5.8 4.3 4.8 3.5 4.8 4.7 5.4 4.9 3.7 5.1 5.1 3.8 3.8 6.1 2.7
##  [937] 3.9 4.8 2.5 4.1 4.5 4.0 3.7 4.5 2.9 5.3 4.1 5.7 5.3 3.8 5.1 5.6 4.3 4.8
##  [955] 4.2 2.2 3.4 4.5 6.0 4.7 4.1 3.7 4.7 5.6 3.9 4.7 4.8 2.9 4.3 3.7 5.4 4.8
##  [973] 5.7 6.3 5.5 4.1 4.9 4.5 4.9 2.3 5.7 5.8 5.8 3.7 5.2 5.0 3.9 5.5 5.8 5.4
##  [991] 3.6 4.5 5.5 4.3 3.4 4.8 4.9 3.7 5.3 4.9
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
##    [1] 3.9 3.5 4.4 4.5 5.2 4.2 5.2 3.3 3.7 4.1 4.7 5.5 5.1 4.4 5.4 5.6 3.5 3.4
##   [19] 2.8 5.4 4.4 4.7 5.3 3.2 5.4 5.4 3.0 4.2 3.6 3.5 4.8 3.7 4.9 5.6 5.6 3.6
##   [37] 4.6 3.2 6.5 6.3 2.1 3.5 4.8 2.7 4.0 5.2 4.4 5.8 4.8 5.4 4.4 4.1 5.8 3.9
##   [55] 3.8 3.8 5.1 3.3 3.8 4.8 3.0 4.9 4.5 4.4 5.6 2.5 3.3 5.3 4.3 3.1 5.8 4.3
##   [73] 4.3 3.9 5.9 4.8 5.7 5.4 4.2 2.8 4.9 4.8 4.6 3.8 4.7 4.0 5.5 4.6 4.6 6.0
##   [91] 3.3 5.0 4.5 2.6 3.4 4.7 5.8 5.0 4.3 4.0 4.7 4.6 4.8 5.5 4.8 5.1 5.4 5.5
##  [109] 3.3 4.0 5.7 6.0 5.3 3.5 4.2 3.9 4.4 4.4 3.7 4.1 3.8 5.5 4.1 6.0 4.6 4.2
##  [127] 4.8 4.4 5.3 5.1 5.0 4.7 5.1 4.3 4.5 5.1 4.7 4.1 5.3 4.4 5.2 3.6 3.1 4.7
##  [145] 5.1 6.5 2.5 4.6 5.5 4.2 4.2 4.6 3.5 2.0 5.6 4.5 4.9 4.0 3.8 4.4 4.8 5.4
##  [163] 5.8 3.6 4.9 2.6 4.5 3.2 3.6 3.7 4.1 4.3 4.7 4.6 5.0 3.3 5.4 3.9 4.3 4.6
##  [181] 6.7 3.8 4.6 5.1 4.5 4.9 3.4 4.3 5.9 5.7 4.1 4.2 4.0 5.2 5.5 5.7 4.9 3.5
##  [199] 5.1 4.6 4.1 3.6 5.3 3.8 5.3 4.1 5.7 3.6 5.0 3.8 3.7 3.6 3.3 5.1 4.3 5.8
##  [217] 4.0 4.9 3.5 4.6 5.8 5.0 3.8 4.2 2.9 4.7 3.4 4.2 4.8 4.7 5.3 4.9 6.1 2.4
##  [235] 5.4 5.9 4.3 5.3 3.6 3.3 3.9 2.8 5.5 6.4 5.0 4.5 5.0 3.9 5.6 5.2 4.1 4.1
##  [253] 4.7 3.2 6.1 4.1 5.0 5.0 3.9 4.9 4.3 5.1 5.0 4.2 5.3 4.8 5.0 4.3 4.4 6.3
##  [271] 4.9 4.4 4.7 4.1 5.7 5.1 4.2 4.5 5.3 5.9 5.5 5.3 4.8 6.6 5.0 4.3 4.5 5.7
##  [289] 4.9 2.4 5.7 5.4 4.2 3.6 5.2 5.1 5.2 4.8 2.5 4.5 4.5 4.2 4.4 4.9 3.3 3.1
##  [307] 4.2 4.6 4.8 4.2 4.8 5.5 3.7 3.6 5.0 3.6 4.5 4.9 3.9 5.9 5.7 3.6 3.9 4.6
##  [325] 3.8 4.9 5.8 6.0 3.7 4.2 4.2 4.3 4.1 4.4 4.1 4.4 5.0 6.3 5.2 3.6 4.6 5.0
##  [343] 4.9 4.3 5.2 3.4 4.6 5.8 4.3 5.9 3.2 5.1 3.9 5.4 4.0 3.9 4.3 3.8 5.5 4.6
##  [361] 4.1 5.0 2.8 5.0 4.0 4.6 5.1 3.5 4.4 5.2 4.9 4.6 4.6 5.6 3.9 4.8 4.6 4.3
##  [379] 4.2 2.5 5.0 4.7 5.6 4.8 4.5 5.3 5.0 5.0 4.8 4.9 4.3 4.0 3.7 5.4 5.7 5.3
##  [397] 2.5 4.8 6.7 4.4 4.6 4.9 5.3 4.4 5.5 4.7 5.5 3.7 5.5 4.6 3.3 5.6 4.5 4.8
##  [415] 4.6 3.9 4.7 5.3 2.6 4.8 5.7 4.2 5.6 4.6 5.4 5.7 4.7 4.8 4.0 4.4 3.7 3.9
##  [433] 3.8 6.3 5.8 5.7 4.0 4.5 4.5 4.8 2.7 3.4 3.8 5.9 4.3 5.8 4.4 4.5 3.7 5.7
##  [451] 3.5 4.4 5.1 4.3 4.6 3.5 6.2 4.7 4.7 5.4 3.5 4.4 4.3 4.8 4.0 5.1 4.8 4.7
##  [469] 5.6 4.0 5.4 4.2 5.7 4.6 3.2 3.5 2.9 3.8 6.0 5.5 4.2 3.7 3.9 4.3 4.2 4.7
##  [487] 3.5 3.5 3.7 4.5 4.6 5.7 4.3 5.4 4.0 3.0 4.7 4.4 4.6 5.4 3.5 3.3 7.2 4.7
##  [505] 4.8 4.8 5.5 4.3 4.6 5.4 5.0 5.6 4.3 5.9 4.4 3.9 4.5 4.6 5.1 4.3 5.8 4.0
##  [523] 5.6 4.8 6.4 5.9 5.3 5.2 4.1 3.9 2.7 5.2 5.4 4.0 2.1 5.3 5.4 3.5 3.9 4.2
##  [541] 5.0 4.1 4.4 4.5 3.6 3.8 5.8 4.5 5.8 5.1 4.4 2.5 5.4 5.4 5.3 3.9 4.6 4.7
##  [559] 5.0 4.0 4.8 3.3 5.3 4.5 5.2 2.7 5.3 3.9 4.8 4.5 4.9 3.9 3.5 3.8 4.2 4.1
##  [577] 2.8 5.7 5.2 3.8 3.5 5.2 5.1 5.5 4.4 3.3 4.3 7.0 4.4 6.3 4.3 2.8 3.7 3.8
##  [595] 4.1 5.3 3.2 3.8 4.2 5.6 5.5 4.2 4.4 5.5 4.3 6.1 3.3 4.6 2.9 4.4 4.2 4.7
##  [613] 4.2 2.8 2.6 5.0 2.1 5.0 3.1 3.2 5.8 5.3 4.4 4.1 5.9 4.8 4.2 4.7 3.4 6.7
##  [631] 4.5 3.7 3.3 5.2 4.3 4.2 4.8 4.9 4.1 5.0 4.6 3.9 3.2 3.8 5.4 3.6 4.9 6.5
##  [649] 4.0 5.1 4.6 5.0 5.0 4.6 4.9 4.9 3.1 3.6 5.4 3.6 5.0 4.9 2.8 4.2 5.4 5.0
##  [667] 4.9 5.0 3.4 4.8 4.2 3.5 2.8 5.3 3.1 3.9 3.4 4.7 4.7 4.7 5.4 5.6 5.5 3.4
##  [685] 4.6 5.1 4.8 6.0 5.8 5.2 4.3 3.6 4.6 5.6 4.7 3.6 4.9 5.4 5.0 5.8 3.2 4.8
##  [703] 3.7 6.4 4.3 5.8 4.6 5.0 3.6 4.4 4.5 6.7 5.4 5.0 5.6 4.5 4.8 5.0 3.3 5.4
##  [721] 5.0 5.3 4.1 3.8 5.0 5.2 3.9 4.0 4.5 4.6 4.7 4.6 5.5 4.9 3.9 4.2 2.6 4.6
##  [739] 5.5 5.6 4.9 3.7 5.3 3.7 4.7 6.8 5.8 5.3 5.1 2.5 5.7 5.0 5.4 4.3 4.2 3.9
##  [757] 3.6 3.6 5.5 5.8 4.3 4.1 5.5 6.2 4.3 4.4 5.7 4.0 3.9 4.5 6.7 5.4 3.9 4.4
##  [775] 6.7 3.6 4.5 4.8 5.0 5.0 4.0 3.7 3.6 5.3 4.5 5.2 3.7 5.2 3.9 2.6 4.0 6.0
##  [793] 3.1 2.9 1.9 6.9 3.5 5.2 3.6 4.9 4.4 4.6 2.6 2.9 3.1 4.4 3.5 5.2 5.9 5.8
##  [811] 3.6 5.3 6.5 5.3 5.8 4.1 4.7 4.1 4.6 5.3 3.9 2.7 3.6 5.1 4.5 4.4 4.1 4.5
##  [829] 4.5 3.3 4.8 4.9 4.0 4.4 5.4 4.3 4.6 4.7 3.5 3.2 4.7 4.6 4.1 2.9 2.7 5.4
##  [847] 3.1 3.0 3.2 4.3 5.6 4.1 3.5 4.4 5.2 5.7 5.0 4.2 6.6 4.6 3.5 5.0 3.8 6.0
##  [865] 4.1 5.3 3.7 3.8 4.5 5.4 4.6 3.3 3.1 3.5 3.9 3.5 4.2 3.1 4.1 5.1 5.5 4.0
##  [883] 4.4 5.9 5.3 4.6 5.1 5.2 5.0 5.2 4.9 5.6 5.4 5.8 2.2 3.2 4.8 3.8 4.2 5.5
##  [901] 4.8 4.1 4.1 6.6 2.7 3.6 4.0 4.9 4.5 4.2 5.6 4.7 3.8 4.7 4.3 4.2 4.1 4.8
##  [919] 5.3 4.8 5.1 4.1 4.0 4.5 4.0 5.6 3.5 4.9 6.1 4.6 5.0 3.3 4.1 4.3 3.6 5.0
##  [937] 5.4 6.7 6.5 6.0 5.3 3.9 4.6 5.1 5.7 3.6 3.9 4.9 3.4 4.7 3.3 6.6 5.7 4.2
##  [955] 5.2 3.7 4.3 5.2 5.1 6.2 4.4 5.0 4.0 5.1 2.1 4.9 4.8 6.3 4.5 4.1 4.6 6.1
##  [973] 4.2 4.4 5.9 4.5 4.2 6.8 4.0 5.0 5.1 4.9 5.1 5.8 4.5 4.2 4.5 3.6 2.9 5.1
##  [991] 5.7 4.6 5.8 3.8 5.3 3.4 5.2 5.8 5.1 4.5
## 
## $func.thetastar
## [1] 0.0544
## 
## $jack.boot.val
##  [1]  0.47206704  0.36345609  0.33156425  0.20420420  0.14000000 -0.01356383
##  [7] -0.10802139 -0.15394322 -0.37701149 -0.45676471
## 
## $jack.boot.se
## [1] 0.8981708
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
##    [1] 4.9 5.6 5.2 2.6 3.6 4.5 5.0 5.2 6.0 6.3 4.2 5.0 4.7 3.7 4.3 3.5 4.3 4.1
##   [19] 4.5 4.0 2.8 4.4 3.6 4.0 5.1 4.0 4.5 2.2 4.7 3.7 2.6 4.3 3.7 5.7 2.7 4.2
##   [37] 5.1 4.9 5.3 4.2 7.0 4.7 5.3 4.8 5.9 5.3 4.8 3.3 4.2 3.6 3.7 4.1 5.5 4.0
##   [55] 4.6 4.2 5.1 4.9 4.4 4.4 4.8 5.7 4.2 5.1 4.0 5.2 5.6 4.5 3.4 3.6 4.0 4.7
##   [73] 4.0 3.2 4.9 4.7 4.2 4.1 4.7 5.4 3.5 3.8 3.0 4.8 3.5 4.2 5.0 5.3 3.9 4.3
##   [91] 5.1 3.3 3.9 3.7 5.1 6.7 3.9 3.1 3.9 5.2 3.9 4.9 3.4 4.9 5.1 5.6 3.9 6.0
##  [109] 5.5 5.2 3.3 5.0 1.8 3.8 3.7 4.7 4.6 3.4 4.5 5.2 4.4 5.5 3.7 4.2 4.8 5.1
##  [127] 5.0 3.7 3.1 4.7 4.9 4.0 4.5 4.7 5.7 4.7 5.1 2.9 4.5 3.9 4.6 4.3 3.8 4.6
##  [145] 4.6 4.4 3.9 2.8 3.4 4.1 5.0 5.0 5.1 4.3 4.4 6.3 4.4 4.4 4.7 5.3 3.7 5.0
##  [163] 4.7 4.5 3.9 4.8 5.9 4.8 4.8 3.8 3.8 4.8 3.6 4.1 4.9 5.1 5.1 4.9 5.0 3.8
##  [181] 4.1 6.2 3.6 4.1 4.0 6.1 2.9 4.1 3.0 3.6 4.7 6.5 5.1 5.0 4.6 4.1 4.9 6.2
##  [199] 4.2 5.8 3.0 3.1 4.5 2.8 2.8 3.7 3.5 4.1 6.5 5.7 4.2 4.5 4.4 4.4 4.3 4.1
##  [217] 3.5 5.2 3.9 3.7 4.7 3.7 4.4 4.3 4.9 2.5 5.2 3.7 5.4 4.0 4.1 5.3 4.8 4.9
##  [235] 3.9 3.5 2.9 4.6 3.3 4.3 5.0 3.4 3.3 4.3 6.1 5.4 4.7 4.7 3.5 3.8 3.2 5.1
##  [253] 4.6 4.2 4.2 5.6 2.3 3.7 6.6 5.7 5.3 3.8 2.6 4.0 3.9 4.8 5.0 3.9 3.4 5.3
##  [271] 3.2 5.2 5.3 6.1 3.8 4.0 5.4 4.1 4.2 4.6 6.6 4.8 5.1 4.1 3.9 2.4 4.6 4.7
##  [289] 5.1 5.0 6.4 4.1 6.3 5.5 3.6 3.7 5.6 4.2 6.3 5.4 5.7 5.2 3.8 5.0 5.4 4.2
##  [307] 3.3 4.3 3.8 4.4 4.5 4.4 5.4 4.9 3.8 4.0 4.5 5.2 4.6 5.3 4.6 4.4 5.8 5.1
##  [325] 3.9 5.5 5.6 4.7 3.8 4.8 3.5 4.1 7.0 5.3 3.9 5.7 3.6 2.5 4.4 5.8 4.2 4.6
##  [343] 5.5 6.2 4.5 5.3 6.4 4.0 5.3 3.6 2.3 3.9 3.6 2.9 5.5 5.0 4.2 5.7 4.4 5.2
##  [361] 4.1 4.1 3.8 4.5 4.3 3.7 5.2 5.4 5.8 3.9 3.9 4.9 3.8 4.7 5.8 4.2 4.8 4.0
##  [379] 4.7 4.7 3.3 4.7 3.3 4.6 6.0 3.2 1.9 5.4 4.3 5.3 3.7 4.2 4.7 4.9 4.8 5.8
##  [397] 3.4 5.5 4.9 5.3 4.6 4.9 5.4 6.5 4.4 4.7 4.8 3.5 4.8 5.4 2.9 4.2 4.7 3.6
##  [415] 4.3 4.4 4.4 4.7 4.0 4.5 4.9 4.1 4.0 6.1 4.4 4.3 2.4 4.0 2.1 5.1 4.9 3.6
##  [433] 4.9 3.3 3.8 4.5 4.3 4.3 4.5 4.2 2.7 4.3 5.6 3.6 2.7 3.9 3.3 5.2 3.3 4.2
##  [451] 4.4 3.9 4.1 3.8 3.2 4.8 4.9 4.9 3.6 4.4 4.8 3.4 3.3 4.4 3.9 4.1 5.8 5.9
##  [469] 4.0 3.7 4.2 3.7 4.0 3.7 4.0 3.6 6.1 3.6 4.1 3.6 6.0 4.4 5.0 4.8 4.2 3.7
##  [487] 4.8 4.9 5.3 5.8 4.0 4.1 3.8 3.9 4.1 3.7 5.4 3.8 4.8 4.8 4.3 5.3 3.5 4.3
##  [505] 2.8 4.8 4.0 4.5 2.1 4.6 4.6 4.5 3.3 5.3 4.2 4.0 5.5 3.5 4.5 4.2 4.4 5.3
##  [523] 4.2 4.7 2.8 3.4 6.1 2.6 4.1 4.0 4.2 3.9 4.0 5.6 5.7 3.8 4.7 2.6 6.8 5.3
##  [541] 3.9 4.8 4.5 4.6 4.7 4.2 3.4 4.7 4.2 5.0 4.3 4.0 4.4 3.0 4.7 3.5 4.4 5.2
##  [559] 3.8 6.0 3.9 2.8 4.8 4.5 5.7 3.6 3.2 3.2 1.8 5.6 5.0 7.4 4.8 4.5 5.1 3.9
##  [577] 4.4 5.9 5.7 2.9 4.4 4.5 4.6 4.0 4.4 4.8 3.5 5.0 3.2 3.6 5.1 6.0 2.6 4.2
##  [595] 4.8 5.0 2.7 5.8 2.6 3.3 2.6 4.1 6.1 3.5 6.0 3.0 4.2 3.2 5.9 3.6 3.7 6.4
##  [613] 3.3 4.3 4.7 4.8 4.1 3.5 6.5 4.0 2.1 3.7 5.0 4.6 3.9 4.7 5.2 3.4 3.4 4.4
##  [631] 5.5 4.5 4.0 5.3 4.7 4.0 4.4 5.1 3.3 4.3 4.2 3.7 3.9 4.1 4.5 5.0 5.1 4.3
##  [649] 4.9 6.0 5.8 3.9 5.1 2.9 5.0 4.7 5.3 4.3 4.8 5.1 4.0 3.9 4.0 5.7 5.3 5.9
##  [667] 2.4 4.7 3.9 4.5 3.6 4.2 4.8 5.2 3.3 4.1 5.0 4.0 4.6 5.7 3.5 3.3 4.3 5.2
##  [685] 4.0 4.2 4.2 6.0 4.4 4.7 5.0 3.5 4.6 4.1 3.2 4.6 2.6 4.8 5.2 5.6 6.5 2.9
##  [703] 3.2 4.3 7.3 3.2 6.1 5.0 5.0 2.5 4.9 5.4 4.1 4.7 5.4 4.7 5.4 4.7 4.1 4.7
##  [721] 6.0 4.4 4.1 3.4 3.7 4.3 3.5 3.3 3.2 6.6 4.7 5.0 4.6 4.3 4.9 5.3 4.2 4.1
##  [739] 4.3 5.8 5.4 6.3 5.5 5.2 5.0 4.5 5.3 5.3 3.6 4.1 4.4 2.9 3.9 4.8 5.0 4.0
##  [757] 4.0 3.8 5.2 4.4 4.4 4.9 3.3 4.4 4.7 4.9 5.1 6.2 6.3 4.3 4.2 2.8 4.8 3.0
##  [775] 3.9 3.8 4.9 5.2 4.3 4.5 3.4 4.8 5.3 3.6 4.4 4.4 3.5 3.5 3.5 6.0 4.1 3.4
##  [793] 5.5 3.9 5.6 4.1 4.7 3.9 4.0 4.5 5.5 6.1 4.9 4.6 4.8 3.8 4.8 4.4 5.0 4.1
##  [811] 5.0 5.4 6.0 4.9 4.9 3.9 5.4 4.5 3.7 5.0 3.9 3.9 4.2 3.3 4.7 4.8 3.3 4.2
##  [829] 3.8 4.3 5.5 4.7 3.7 4.7 4.3 4.6 4.0 4.2 4.9 4.5 4.1 4.0 5.0 5.5 4.5 5.0
##  [847] 5.5 5.0 3.9 5.3 3.8 5.5 4.4 4.1 4.5 4.4 3.7 3.3 4.3 4.2 5.0 3.3 4.2 3.6
##  [865] 4.1 4.4 5.1 5.0 5.1 4.4 4.7 4.0 6.1 4.2 3.5 4.9 3.3 5.5 5.3 4.0 4.4 5.2
##  [883] 4.1 6.3 4.7 4.4 5.6 4.2 5.1 5.0 3.4 4.3 5.6 4.4 4.5 3.6 3.6 4.1 3.9 5.5
##  [901] 4.6 4.3 5.6 5.5 3.7 2.8 5.1 3.9 4.8 3.7 5.9 4.6 4.9 5.2 3.7 5.8 2.8 5.5
##  [919] 4.8 3.9 5.9 3.2 4.7 4.6 5.7 4.9 4.3 5.3 3.7 5.1 5.6 4.0 5.3 2.8 4.0 5.2
##  [937] 3.9 5.5 5.9 5.5 4.8 5.2 5.3 4.4 5.3 3.4 2.6 4.4 5.7 5.0 3.5 5.3 3.6 4.6
##  [955] 5.2 5.6 4.5 4.2 3.6 4.1 4.0 3.7 5.1 4.8 5.1 4.7 5.1 4.5 4.5 4.5 4.2 4.5
##  [973] 5.2 3.7 4.3 2.6 4.6 4.5 5.8 6.3 5.1 6.6 3.0 4.1 6.5 4.7 4.6 3.9 6.0 4.3
##  [991] 4.2 3.8 4.5 5.3 4.5 5.4 3.7 5.3 4.6 6.1
## 
## $func.thetastar
## 72% 
##   5 
## 
## $jack.boot.val
##  [1] 5.4 5.3 5.2 5.1 4.9 4.9 4.8 4.7 4.6 4.4
## 
## $jack.boot.se
## [1] 0.9104395
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
## [1] 1.216774
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
##   3.128156   4.627659 
##  (1.331117) (2.135961)
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
## [1]  0.2908519  0.8158708  1.0479955  0.9007109  0.6333657 -0.4447578
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
##    [1] -0.312727468  1.067808217  0.361835581  0.248570919  1.150259185
##    [6]  1.210223988  1.177628432  1.324140690  0.314311484 -1.027002925
##   [11]  1.244096850  1.543889064  1.667209947  1.053308355  1.143360856
##   [16] -0.169184427 -0.337683101  0.221585154  1.048377227  0.701027654
##   [21]  0.641016790  0.249140933  0.632591093 -0.033395891  0.591480329
##   [26]  1.550927803  0.949553367 -0.839821931  0.077031505  0.394923689
##   [31]  1.055476454  1.328997237  1.568190147  0.145383858  1.042935438
##   [36]  0.427715139  1.447619512  0.970492654  0.404435791  0.297834840
##   [41]  0.642418299 -0.448748934  0.985879519  1.124064965  1.080699737
##   [46]  1.626005972  0.450070546  1.267298184  0.758587014  0.888759105
##   [51] -0.036736408  0.231972968  0.189693801  1.303309198  0.280357548
##   [56]  0.382992089  1.006954805  0.980154063  0.226902103  1.144323015
##   [61]  1.156978770  0.434956637  1.921235920  1.059983559  0.558355099
##   [66]  0.689785707  0.961281705  1.585171777  1.270353850  0.005725954
##   [71]  1.055567965  0.273820445  0.831370679  1.240832682  1.219957646
##   [76]  0.028858557  1.117795140 -0.039816279  1.020334094  1.903607961
##   [81]  1.086059650  0.685615416  1.131456114  0.925208299  0.278641221
##   [86]  0.595495549  1.759756175  0.058044257 -0.154717935  0.161889380
##   [91]  1.416895769  0.195007497  0.860633881  0.310565782  0.927887444
##   [96]  0.513864396  0.927786334  1.387152113 -0.903442265  1.155696630
##  [101]  0.297031246  1.426824970  0.948791539  0.654318326  0.410950064
##  [106]  1.372851485  0.210424145  1.259168849  1.053543765  0.865970925
##  [111]  2.131895096  0.682735518  0.773957466  0.947383722  1.407998151
##  [116]  1.458457195  1.540909011  0.542830902  0.427054943  1.612178078
##  [121]  0.620274285  0.014946606  1.710094845  1.676257903  1.636158598
##  [126]  1.269990453  0.757314994  0.300462283  0.155044986  1.130254609
##  [131]  1.499933681  1.441329014  0.397540133  1.778865266  1.368344177
##  [136]  0.793283703  0.611755874  1.788577831  0.687563556  0.960565338
##  [141]  1.183284128  0.726419836 -0.145957472 -0.484741468  1.408723666
##  [146]  1.078324190  1.457239827  0.974890657  0.816033639  0.802333068
##  [151]  1.133932431  0.660600855  0.990023776  0.919972769  1.382543827
##  [156]  0.976553655  1.881009440  1.275026192 -0.144543025  1.378000049
##  [161]  0.829386127  1.909675463  0.870571928  0.657979571  1.039976880
##  [166]  0.823015525  0.173411962  0.875220153  0.451263862  1.026079634
##  [171]  1.515801909  0.040402235  0.569215113  1.386184735  0.853529854
##  [176]  1.069287556  0.826836269 -0.240575349  0.491261119  0.005744811
##  [181]  1.775212755 -0.716544303 -0.207541211  1.353746739  1.430282960
##  [186] -0.148842198  0.763503782  0.905844447  0.630764862 -0.622215721
##  [191]  1.245535275  1.210058691  0.959628698  0.599233562  0.632416579
##  [196]  0.644101964  1.090379550 -0.315786276  1.592760362  0.783563987
##  [201]  1.057579153  0.828197402  1.820063370  0.822912123  0.919291774
##  [206] -0.833289603  0.599208312  0.726204352  1.508086332  0.098816641
##  [211]  0.425282996  1.010981007  2.073888539 -0.109919901  1.885939532
##  [216]  1.367219980  0.425897334 -0.227088721 -0.157980956  1.004145129
##  [221]  0.442670371 -0.105055337  0.910439445  1.751367458  1.011893050
##  [226]  0.392648710  0.791111088  1.631690473 -0.102361079 -0.451608297
##  [231]  0.753225754  0.737570395  2.068032625  0.795884714  1.098797972
##  [236]  1.095453621  1.319575224  1.298387844  0.262636812  0.745697017
##  [241]  0.760968374 -0.624756718  0.420140842  1.416022447  0.843728907
##  [246]  1.051373558  0.557494995  1.282821864  1.614833611  1.518963499
##  [251]  0.767364334  0.745821331  1.161397262  1.374712767 -1.379945718
##  [256]  0.205315770  1.693470044  0.335683576  1.106891793  0.838381015
##  [261]  2.106841788  0.451654168  1.157964009  1.487464250  1.122680380
##  [266]  0.637250211  0.182918365  0.377083699  1.036852800  1.130492659
##  [271] -0.453557069  1.573704572  0.905621822  1.815227233  1.501747102
##  [276]  0.710332773  1.400515876  1.512667235  1.017475766 -0.162856402
##  [281]  0.319008786 -0.132761200  0.139038265  0.319246134  0.119101948
##  [286]  0.901288545 -0.257590665  0.495618685  0.459607341  1.413387902
##  [291]  0.920927317  0.773270901  0.142098589 -0.228012129 -0.118214319
##  [296]  1.689748230 -2.532878960 -0.601729678  1.341057472 -0.281844903
##  [301]  1.068190601  1.386393988  0.503747874  1.212579603  1.378767773
##  [306]  1.326714879 -0.440148652 -0.410139686  0.930973869  1.043753184
##  [311]  1.292963124  0.953187577 -0.179739560  0.403504374 -0.480579062
##  [316]  1.031062898 -0.006969171  0.665640420  1.389146987  0.442480899
##  [321]  1.050074676  0.002653965  1.625321804  0.317534811  0.752028207
##  [326]  1.236734467 -0.021924313  0.983138679  0.654053846  0.495382392
##  [331]  1.420432130 -0.514983725  0.422508632  1.071463607  0.739961459
##  [336]  0.200039137  0.876522518  1.986396411  1.283705321  1.519916309
##  [341]  0.830542783  0.549210766  0.523893460  0.335027924  1.051662718
##  [346]  0.353570492  1.027808067  1.717336393 -0.286708090  0.734290140
##  [351]  0.980154063  1.609877242  0.097706866  0.754445174  0.936080497
##  [356]  0.371450879  0.641429357  1.265541936  1.616672020  1.055360332
##  [361]  1.362326921  1.224691619 -0.195974144  0.993475684  0.815442047
##  [366]  0.584809167  0.362576667  0.381039643  0.352225577  0.404050524
##  [371] -0.147254991  1.111142171  0.823574027  1.084192722  0.957182283
##  [376]  0.998612485  1.409491969  1.248379935  1.644226846  1.249490242
##  [381]  1.072964480  1.526110298  0.862645125  1.177625076  0.053952565
##  [386]  1.647877576  1.609234869 -0.076169633  1.781370226  1.299373319
##  [391]  1.729645763  1.363291795  0.263027733  0.597852031  0.003553836
##  [396]  1.228432914  0.884061097  0.862423972 -0.006035465  0.255733883
##  [401] -0.072569912 -0.514690498  0.818276733 -0.074413644 -0.092784396
##  [406]  0.703206862  1.127472037 -0.135287985  1.216574919  1.033019253
##  [411]  0.016564640  1.254706478  0.433632422  1.166255163  0.103713350
##  [416]  1.631341778  0.869100426  0.654127728 -0.481859229  1.406375894
##  [421]  0.689481118  0.576891015 -0.224648385 -0.184526229  1.730209748
##  [426]  1.247830960  1.923257404  0.284306061  1.055526300  0.148819986
##  [431]  0.462767580  0.446784463  0.241038520  0.646756960  0.783937661
##  [436]  1.123064523  0.171255304  1.270714170  0.461708931  0.766249716
##  [441] -1.337074080  0.853121942 -0.841208874  1.499183076  0.177168261
##  [446]  1.359114197  0.213450121 -0.073497812  1.684997598  0.090403255
##  [451]  0.432645023  0.801030406  1.303437493  0.970387790  0.436835722
##  [456]  1.463568686  0.941560182  1.900659324  0.982190985  0.887158716
##  [461]  1.726525236  0.150037897 -0.047820852  0.288363180  1.319772534
##  [466]  0.531285226  0.625897714 -0.440631407  0.967007503  0.664908895
##  [471]  0.733889368  0.717113739 -0.641352714  0.125274033  1.202994129
##  [476]  1.454779628  0.624166792  1.682091934  0.716776574  0.876658834
##  [481]  1.082625764  0.264414638 -0.210417352  0.486950386 -0.178444830
##  [486]  1.046955134  1.077846587 -0.893968589 -0.140909549  0.599086098
##  [491]  0.105797905  0.829798872  0.446863620  1.077983661  1.785493429
##  [496] -0.482996002  1.591758881  0.180493837  0.405528132  0.561794992
##  [501] -0.206699621  0.613036353  0.734322603  0.175768635  0.939267970
##  [506]  1.626173570  0.877864866  1.125840160  1.732849008  0.863858408
##  [511]  1.293003999  0.049237571  0.868753865  2.206596452 -0.423316860
##  [516]  1.821471889  1.881279148  0.388358744  0.848314290  0.867514870
##  [521]  0.187402325  0.179084288 -0.081810163 -0.654588537 -0.965759047
##  [526]  0.557584862  0.983128504  0.883940315  1.677423344  0.385232023
##  [531]  1.411954447 -0.387912131  0.906128249  1.191060379  0.379512316
##  [536]  1.590970766  0.494688439  0.428557057  0.745049816  1.692025439
##  [541]  0.976451681  1.006321252  0.956200361  0.373287015 -0.306737775
##  [546]  1.672762184  0.853306340  1.061862511  1.595502383  1.862881840
##  [551]  1.748482403  0.495802480  0.788151895  0.678707949  1.168920319
##  [556]  0.759203284  1.134960996  1.260837796  1.040832322  0.901219094
##  [561]  1.609403828 -0.094692311  1.431659803  0.990004978  1.029883793
##  [566]  1.987874709  1.738426430  1.065127204 -0.086287835  2.143179411
##  [571]  1.136515481  0.944814332 -0.337729978  1.083477360  1.764316820
##  [576]  1.470021282  0.187402325  1.006440348  1.124743004  0.771624588
##  [581]  0.221665443  1.644504206  0.006636287  1.254464516  0.829572233
##  [586]  1.283243020  0.145226526 -0.318416093  0.985896484  0.411435517
##  [591]  0.646266107  0.861532877 -0.632492517  0.610861059  0.538509160
##  [596]  1.831441115  1.404018475  1.059215636  1.069548110 -0.343924632
##  [601]  0.264221491 -0.090066491  0.037011611  0.858649179  0.359568329
##  [606]  1.105727055  1.243382835  0.978171182  0.899635072  0.536591752
##  [611]  1.261908930  0.806785395  0.475827638 -0.279245446  1.081182684
##  [616]  0.020183445  1.388304656  0.591863179  1.130400664  0.349569496
##  [621]  0.137514271  0.380639377 -0.494311182  0.837208810  1.940096052
##  [626]  0.381810416  1.071712808  0.914689493  0.473964568  0.010726647
##  [631]  0.991288194  0.779068184  0.515089244  0.273316656  1.737142322
##  [636]  0.514209061  0.342572678  1.434644220 -0.157435642  1.264117612
##  [641]  0.816708405  0.191222518  0.744279011  0.882283099  0.814142932
##  [646]  0.174295775  0.817404938  0.100740855  0.631976803  0.365684881
##  [651]  0.916348925  0.803230008  1.283552008  1.021258008  0.988298973
##  [656]  1.082438179 -0.952116711  0.576766602  0.017630321  1.046466679
##  [661]  1.636526170  0.241752392  0.150651862  0.398735158  1.112679888
##  [666]  1.771386961  1.630025304 -0.479475716  0.782057894  1.187980042
##  [671]  0.471247384  0.176198218  2.084910021  0.961281705  1.182807664
##  [676]  1.778219044  0.404569567  0.630687312  0.589898833  1.518827931
##  [681]  1.096354981  2.060223657  0.836212806  0.568769817  0.369548834
##  [686]  0.695408394  0.434736157  0.817150457  1.013461182  0.830114109
##  [691] -0.553376953  0.537932587 -0.304869735  1.342605795  0.647956113
##  [696]  0.888759105  0.286506003  0.899542578  0.103713350  0.440619433
##  [701]  0.669006939  1.742788618  1.283125356  1.061601829  0.484911522
##  [706]  0.182106152  1.327504959 -0.532813555  1.698321092  0.782111462
##  [711]  0.539051911 -0.071131168  0.739201459  0.746961816  0.154105734
##  [716]  1.689748230  0.417297829  0.139963483  0.522667702  0.383537286
##  [721] -0.092609577  1.288497907  1.997511524  0.673694766  0.154814754
##  [726] -0.393337767  0.507481112  1.004745888  1.292677589  0.187326409
##  [731]  0.619424433  1.607200420 -1.359461003  0.153875667  0.678502211
##  [736]  1.274491034  0.558622840  1.591633586  0.426130732  1.238947676
##  [741]  0.486560900  0.995476281  0.467641308  0.257204396  0.357809639
##  [746]  0.381662695  1.671096624  0.873724497  0.686115833  0.161629152
##  [751]  0.742664351  1.491187936  1.240871889  1.626449526  1.046709075
##  [756]  1.333709022  1.614416233  0.417058052  0.241895462  1.118206625
##  [761]  0.258432571  0.503364852  1.879094797  1.316112186  0.506569104
##  [766]  1.212756247  1.225572613  1.875667022  1.234389371  0.738934709
##  [771]  0.775091050  0.838807870  1.174198440  0.850402034  0.958979461
##  [776]  0.455437241  0.689772806  0.410796550 -0.248616010  0.965873699
##  [781]  1.766761401  0.961885502  0.965056913  1.907594192  1.568108616
##  [786]  0.237729390  0.026650972  0.883596646  1.363460701  1.092260208
##  [791]  1.430621810  1.551710891  1.140749807  0.849853092  1.103984301
##  [796]  0.766766876  1.033019253  1.387736348  1.307234998 -0.116831556
##  [801] -0.848906923  1.107730452 -0.049732965  0.525558342  1.301138508
##  [806] -0.048062473 -0.953236642 -0.526900450  1.201967016  1.238059994
##  [811]  1.729447850 -0.083712890  1.679446120  0.750938166  1.324140690
##  [816]  1.098094226  1.340079214 -0.776433223  0.267286174  0.552585551
##  [821]  1.105983789  1.213842900  1.795829350  1.063092729  1.083506154
##  [826] -0.226883386  0.703736194  0.937919616  0.589237286  1.007225811
##  [831]  0.363430834  0.720612987  1.038504868  0.182032695  0.774994023
##  [836]  0.983566785  0.836432607  1.313115721  1.342700575  0.319246134
##  [841]  0.176421104  0.423475843  0.948941473  0.056028533  0.642418299
##  [846]  1.328510454  1.348790325  0.808308729  1.373009984  0.229683968
##  [851]  0.957894674  0.550102543  1.130400664  0.829843540  1.347628116
##  [856] -0.044105419  0.226565109  0.463896197 -0.271474553  1.601023259
##  [861]  1.382546714  0.925544131  0.032398797  0.886602194  1.710476931
##  [866]  1.663456084 -1.005081446  0.859607252  0.929149944  0.696543323
##  [871]  0.605646300 -0.137882890  0.597731149  1.086769906  1.850333440
##  [876]  1.056654053  1.356002013  0.665762880  0.735054404 -0.040702922
##  [881]  1.230466577 -0.251867746  0.472397744  1.035829554  1.079343157
##  [886]  1.031400424  0.645900474  0.781798103  0.868060690  1.148320806
##  [891] -0.266396172  0.776281346  0.273699961  1.190397693  1.542767886
##  [896]  0.944817948  1.303450552  0.255960046  1.152421833  1.536884807
##  [901]  0.878255250  0.316635410  0.570203032  1.258378762  0.425727848
##  [906]  0.931653546  0.073612565  1.558665195  1.405706006  1.033925736
##  [911]  0.865229625  0.266947592  1.271498995  1.183067096  0.292106664
##  [916] -1.027270431  0.950237196  1.820938104  0.792136097  1.726669221
##  [921]  0.125216538  1.303564197  0.803771923  1.254181652  0.043302394
##  [926]  0.230419524  1.389911453  1.437204326  1.727521848  1.500291173
##  [931]  0.959363657  0.818702603  0.076143625  1.106237062  1.336026178
##  [936]  1.248628364  1.106891793  0.958203475  1.708568654  1.785719293
##  [941]  1.102232663  2.135880291 -0.096641969  0.809941860  0.026120531
##  [946]  1.152830102  0.546084956  1.598797414  1.070710749  0.825917350
##  [951]  0.478533128  1.318321747  0.578258916  0.224801411  0.935996259
##  [956]  0.989646925  0.788359180  0.109687361  1.316441853  1.243282840
##  [961] -0.088721136  0.452842843  1.297742384  1.932584063  0.152527217
##  [966]  1.355695961  0.021157793  1.303437493  0.179965227  0.300862761
##  [971]  1.343386895  0.588615860  0.156059733  1.668028361  0.941649107
##  [976]  0.097234132 -0.014937090 -0.015007508  0.302001758  1.291800740
##  [981]  1.163087501  0.679411310 -0.456936439  1.151974689  1.981028015
##  [986] -1.292395778  1.647744361  0.983605135  0.976606808  1.249602088
##  [991]  0.813804825 -0.031281875 -0.193281857  0.370068295  1.091670143
##  [996]  1.509997209  0.659337074  1.640784795  0.783466501 -0.576654177
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
##      mean         sd    
##   0.6759609   0.4059710 
##  (0.1283793) (0.0907758)
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
## [1]  0.15212547  0.17529847 -0.39507904 -0.14098788 -0.07128585  0.25651132
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
## [1] 0.0506
```

```r
sd(results2$thetastar)  #the standard deviation of the theta stars is the SE of the statistic (in this case, the mean)
```

```
## [1] 0.94387
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
## t1*      4.5 -0.01611612   0.9173623
```

One of the main advantages to the 'boot' package over the 'bootstrap' package is the nicer formatting of the output.

Going back to our original code, lets see how we could reproduce all of these numbers:


```r
table(sample(x,size=length(x),replace=T))
```

```
## 
## 0 2 4 5 7 8 9 
## 1 1 2 2 2 1 1
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
## [1] -0.029
```

```r
se.boot
```

```
## [1] 0.930766
```

Why do our numbers not agree exactly with those of the boot package? This is because our estimates of bias and standard error are just estimates, and they carry with them their own uncertainties. That is one of the reasons we might bother doing jackknife-after-bootstrap.

The 'boot' package has a LOT of functionality. If we have time, we will come back to some of these more complex functions later in the semester as we cover topics like regression and glm.

