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
## 0 1 3 4 5 6 7 
## 1 2 2 1 1 1 2
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
## [1] -0.0338
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
## [1] 2.65897
```

```r
UL.boot
```

```
## [1] 6.27343
```

Method #2: Simply take the quantiles of the bootstrap statistics


```r
quantile(xmeans,c(0.025,0.975))
```

```
##  2.5% 97.5% 
##   2.6   6.3
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
##    [1] 4.4 3.4 5.1 5.0 5.7 4.5 5.1 6.3 4.1 5.4 3.4 6.0 5.7 6.2 4.4 5.1 3.7 4.8
##   [19] 4.8 5.2 6.6 5.7 5.0 5.8 6.2 3.0 6.7 5.8 4.5 4.1 4.4 4.3 4.2 4.5 5.0 3.1
##   [37] 4.2 3.7 2.5 5.9 3.5 5.1 4.6 3.7 4.0 4.4 5.0 5.7 4.7 4.5 4.4 4.9 4.2 3.6
##   [55] 5.9 5.4 3.6 4.2 4.4 4.3 4.6 4.2 4.0 5.0 2.8 3.5 4.5 5.8 4.8 5.0 4.6 4.6
##   [73] 3.8 5.5 3.9 6.8 3.5 5.3 4.2 4.3 5.4 5.3 3.0 4.1 4.9 4.3 6.7 4.7 3.8 4.9
##   [91] 5.9 4.1 4.3 5.6 4.5 4.7 3.9 3.0 3.5 4.3 5.6 3.9 3.2 3.7 4.6 3.7 3.7 4.4
##  [109] 6.4 4.4 3.4 4.7 4.0 4.5 5.9 4.8 3.3 6.3 4.8 4.9 4.7 5.4 4.0 3.7 4.9 4.3
##  [127] 3.5 3.9 4.0 4.9 4.8 5.8 6.5 2.9 3.9 4.5 3.6 4.0 3.2 4.6 4.3 4.0 3.4 4.5
##  [145] 4.4 3.5 5.0 4.0 5.9 5.0 3.2 2.9 4.2 4.1 4.6 6.1 3.6 3.3 4.0 4.4 5.0 3.5
##  [163] 4.7 4.4 6.1 4.8 5.0 4.7 4.1 3.1 3.9 5.5 5.2 4.9 2.9 5.1 4.2 3.2 4.0 3.4
##  [181] 4.1 6.2 4.8 5.6 5.0 6.1 4.7 5.7 4.4 4.2 4.5 4.6 3.9 4.3 3.5 5.4 3.1 5.0
##  [199] 5.3 5.0 5.8 3.7 5.0 3.5 5.3 4.0 5.2 2.8 5.3 5.7 3.7 6.0 4.8 5.0 4.1 4.7
##  [217] 2.7 5.5 4.5 5.2 4.6 4.7 5.9 3.4 5.9 4.5 4.6 6.1 2.9 5.2 5.9 3.9 5.2 4.8
##  [235] 4.4 3.9 3.5 4.4 4.4 5.6 5.0 5.0 4.6 4.0 3.6 3.2 3.1 3.7 4.2 5.8 4.6 4.7
##  [253] 5.2 5.0 5.6 4.8 4.8 5.6 5.3 4.6 4.6 3.8 3.0 3.6 5.6 5.3 3.1 4.0 4.1 4.6
##  [271] 3.6 3.9 4.2 3.4 5.5 4.9 4.6 5.8 4.2 6.0 5.9 4.9 6.3 5.3 4.5 3.5 2.6 6.1
##  [289] 4.9 4.7 4.8 3.8 5.3 3.6 4.3 6.4 4.6 3.0 3.9 3.3 6.0 3.7 4.6 4.5 4.1 4.1
##  [307] 5.2 6.2 6.1 4.4 4.6 6.0 4.1 4.2 4.6 4.5 4.7 3.9 4.6 2.4 3.7 4.5 4.6 4.0
##  [325] 4.6 4.5 3.8 4.1 5.7 4.4 5.1 5.9 4.9 3.6 4.2 5.1 3.9 5.2 2.9 4.2 4.2 4.4
##  [343] 4.8 3.4 3.9 4.6 4.7 4.0 3.9 4.2 3.5 4.7 3.9 3.7 5.8 4.8 4.3 4.2 5.2 3.8
##  [361] 4.9 4.2 3.6 3.3 6.3 5.7 4.6 5.2 4.9 4.6 5.2 3.4 4.0 6.1 3.4 4.0 2.9 3.6
##  [379] 6.1 5.2 4.2 4.8 3.8 3.3 5.8 3.6 5.1 6.2 5.0 3.0 4.5 4.2 5.9 3.3 4.9 3.4
##  [397] 2.8 5.3 4.7 4.9 4.0 5.0 4.3 4.2 6.6 3.8 4.9 4.5 5.5 3.7 3.3 4.0 4.2 4.4
##  [415] 4.3 5.8 3.5 4.2 4.8 4.8 3.9 4.4 5.3 4.6 3.7 3.8 4.2 3.7 3.8 4.3 4.5 3.8
##  [433] 4.7 5.8 4.8 5.4 3.2 5.1 4.3 4.1 3.7 6.0 3.6 5.4 5.2 5.6 4.9 4.7 4.4 3.1
##  [451] 4.7 3.6 2.7 4.6 4.1 4.2 4.8 4.8 2.5 3.3 3.9 5.7 4.8 4.5 4.6 5.6 4.6 5.0
##  [469] 3.2 7.2 6.4 3.5 4.7 4.7 5.3 6.5 2.7 2.4 4.6 6.0 4.2 5.3 4.0 4.6 3.0 4.0
##  [487] 4.5 3.2 5.5 5.2 1.9 4.2 4.3 3.9 3.8 5.2 2.9 4.2 4.8 4.2 3.6 3.1 4.3 3.0
##  [505] 3.0 2.8 4.9 3.5 5.0 2.9 4.9 5.5 4.3 4.7 3.6 4.1 3.1 4.4 4.0 5.0 3.0 3.6
##  [523] 5.0 4.7 6.4 4.2 4.3 4.6 5.8 6.6 4.9 3.5 5.4 5.5 3.7 5.9 4.2 4.8 6.2 4.9
##  [541] 3.7 5.6 4.3 4.1 6.0 5.7 5.6 5.2 4.4 5.6 3.7 4.0 4.0 4.8 3.5 3.0 3.8 4.4
##  [559] 3.4 4.6 5.0 5.1 6.4 3.0 5.9 2.0 3.7 4.8 3.3 3.2 5.1 5.8 4.5 4.6 3.9 5.5
##  [577] 3.1 4.8 4.0 2.5 4.3 4.8 4.8 5.4 4.8 5.0 5.5 5.5 5.3 6.5 4.6 5.7 4.9 5.7
##  [595] 5.0 4.6 3.6 4.6 5.0 5.6 5.9 4.9 4.0 6.1 5.8 4.3 3.6 3.3 5.3 2.8 4.0 5.6
##  [613] 2.8 4.4 5.2 5.4 3.9 3.2 3.7 4.8 3.9 4.5 4.0 4.5 4.6 5.0 5.7 4.8 4.3 4.8
##  [631] 5.7 6.3 5.6 5.6 5.1 2.9 6.4 4.5 3.7 4.8 5.4 4.7 4.2 6.1 4.9 6.3 5.5 4.0
##  [649] 4.0 4.9 5.3 4.3 4.6 4.4 5.2 4.4 4.0 3.4 4.0 4.5 3.8 3.7 3.6 5.8 3.9 5.3
##  [667] 2.8 4.1 3.5 5.6 3.6 3.4 3.2 4.0 4.5 3.9 5.1 6.0 4.9 5.0 4.5 4.6 4.4 5.2
##  [685] 5.3 3.7 5.0 2.6 4.1 6.1 4.2 5.0 5.5 3.9 5.3 4.5 3.4 5.6 5.1 4.9 4.6 3.3
##  [703] 4.9 5.8 4.6 5.5 4.3 4.4 6.5 3.2 2.7 4.9 4.4 3.9 6.3 5.4 5.2 6.1 5.2 4.1
##  [721] 4.4 4.5 3.9 3.8 4.6 3.7 6.7 4.8 5.2 5.2 4.3 5.2 5.2 4.6 4.5 2.9 3.6 4.2
##  [739] 3.7 3.4 5.7 4.6 6.4 4.2 5.1 4.5 3.9 6.2 4.4 4.1 4.9 3.0 5.6 4.6 3.5 4.9
##  [757] 3.8 3.0 2.9 4.9 5.8 5.3 4.6 3.2 5.1 5.1 4.6 3.1 5.8 4.4 5.8 3.2 5.5 6.0
##  [775] 3.3 5.0 4.1 4.6 4.5 5.1 5.8 3.8 6.0 5.7 5.9 5.2 2.9 4.3 4.6 4.7 5.1 4.1
##  [793] 3.9 5.5 5.1 5.5 3.9 3.4 4.1 4.6 2.5 4.0 4.9 5.1 5.2 4.5 5.4 4.3 4.4 3.7
##  [811] 3.7 4.4 4.6 6.1 3.5 4.2 4.7 5.2 3.1 4.9 5.1 3.9 3.9 4.5 3.2 4.8 6.0 4.9
##  [829] 4.0 4.3 4.1 4.8 4.6 5.2 4.5 5.5 4.8 5.4 4.6 4.5 4.7 4.6 4.5 4.3 3.8 4.6
##  [847] 3.5 3.5 4.4 5.0 3.8 4.8 4.2 4.0 4.4 4.3 4.6 4.0 5.4 4.7 4.7 3.1 5.2 6.4
##  [865] 3.5 4.2 3.4 3.5 2.8 3.4 5.8 3.5 4.1 5.0 3.0 4.6 4.7 4.2 3.5 5.9 5.4 4.7
##  [883] 4.9 5.5 3.9 4.7 4.4 3.5 5.3 3.9 3.8 3.9 6.4 5.4 6.4 6.2 4.4 3.4 5.8 5.1
##  [901] 4.6 3.9 5.7 4.7 3.6 5.6 4.2 6.0 3.6 3.9 3.0 3.0 6.1 4.5 5.2 4.8 6.6 3.6
##  [919] 3.5 4.5 4.4 5.4 4.2 3.4 4.9 4.5 5.8 4.3 4.8 5.5 5.3 2.6 5.8 4.9 3.3 6.7
##  [937] 5.4 3.5 5.4 4.8 4.9 4.0 4.5 5.1 4.0 3.6 5.1 3.9 5.6 4.7 3.3 3.8 4.1 5.6
##  [955] 3.4 3.4 3.7 4.4 3.4 3.8 5.4 3.3 4.2 5.9 2.7 4.0 4.4 3.6 3.7 4.7 5.3 5.1
##  [973] 3.9 3.7 4.8 3.6 5.7 4.1 4.7 5.0 2.5 4.0 4.6 3.1 4.8 2.5 3.4 3.4 4.4 3.4
##  [991] 3.4 6.1 4.0 4.1 4.0 4.6 3.1 6.6 5.6 4.9
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
##   2.5%  97.5% 
## 2.8000 6.3025
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
##    [1] 5.7 3.8 2.9 4.4 4.1 4.1 4.0 5.1 4.0 2.6 3.0 4.3 5.4 3.2 6.2 4.0 5.4 5.3
##   [19] 5.4 3.4 5.5 3.9 5.2 2.6 4.6 5.2 4.0 3.7 2.9 4.1 4.6 5.6 3.6 3.7 3.8 3.4
##   [37] 4.7 4.5 4.6 4.4 5.2 4.2 3.3 5.5 4.9 5.8 4.9 5.2 5.8 5.9 2.7 4.9 5.3 3.9
##   [55] 5.4 5.3 5.3 5.0 5.0 5.1 3.6 6.2 4.0 4.9 2.8 3.7 6.1 4.4 3.5 3.4 4.8 4.5
##   [73] 4.4 3.4 4.6 4.6 5.7 4.1 5.8 3.7 2.6 3.1 4.8 5.1 5.1 3.7 5.6 4.8 7.1 4.9
##   [91] 3.9 5.0 2.8 5.8 3.1 3.7 5.4 4.5 4.5 5.3 4.0 2.9 4.8 3.2 3.4 4.2 5.1 4.1
##  [109] 4.7 3.6 4.0 4.6 5.1 5.6 4.1 4.5 3.5 2.9 3.4 4.5 3.8 4.6 4.0 2.3 4.4 4.4
##  [127] 5.5 2.7 4.8 4.9 4.4 4.0 3.8 3.5 5.8 3.1 5.1 3.6 5.9 4.7 5.6 4.7 4.9 3.9
##  [145] 3.9 4.5 5.1 3.8 5.0 4.2 3.8 3.8 5.5 3.4 5.5 3.8 6.2 4.6 2.5 2.9 5.0 5.3
##  [163] 3.0 3.6 4.1 4.8 3.5 3.4 4.0 4.1 4.7 4.8 3.4 2.9 4.0 4.4 3.2 4.7 5.8 4.5
##  [181] 5.0 5.3 5.5 4.8 6.0 2.3 4.3 4.1 4.3 5.1 5.5 5.7 4.9 4.9 3.8 5.0 3.9 3.9
##  [199] 2.8 4.6 5.3 3.8 3.5 3.8 4.7 4.0 4.7 4.0 5.7 4.3 6.0 2.9 4.8 4.9 4.7 4.4
##  [217] 3.8 4.3 5.5 5.5 4.3 5.1 2.4 6.1 4.0 4.5 5.0 4.5 5.5 5.2 3.6 4.1 4.6 4.2
##  [235] 4.0 3.6 4.4 5.7 5.7 4.9 4.0 5.7 3.7 3.5 5.2 4.6 4.4 3.6 4.5 2.2 4.9 4.9
##  [253] 4.6 4.1 5.1 4.5 4.5 3.5 3.3 3.9 5.5 3.3 3.6 4.8 4.1 5.0 3.7 3.2 3.2 5.4
##  [271] 4.3 3.3 4.5 4.1 5.2 3.7 5.4 3.9 2.9 4.1 6.2 6.1 4.2 4.2 4.6 3.4 5.3 5.4
##  [289] 3.7 3.0 3.7 4.6 4.9 5.0 6.8 6.3 5.2 3.1 4.6 3.8 3.6 2.5 3.5 5.4 5.3 3.3
##  [307] 2.8 5.7 3.4 5.8 4.5 5.1 5.2 5.1 6.0 5.3 3.9 3.2 4.1 2.7 5.0 4.8 6.1 4.4
##  [325] 5.1 5.5 7.4 3.8 4.5 3.6 4.8 4.9 4.7 4.7 5.6 4.5 3.3 2.2 3.4 3.6 3.8 4.9
##  [343] 4.9 2.6 6.9 4.1 5.7 5.3 4.2 4.4 3.8 4.4 3.6 5.7 4.0 3.4 3.5 4.1 4.2 4.1
##  [361] 4.8 5.4 4.0 4.7 6.2 5.1 5.5 3.3 6.4 4.3 6.0 6.5 4.2 5.4 4.7 5.8 5.7 3.9
##  [379] 4.8 5.0 5.0 4.4 4.5 5.3 5.3 4.9 3.7 4.7 2.8 6.7 4.2 4.3 5.6 4.5 4.1 4.6
##  [397] 5.7 3.7 4.1 5.8 5.0 4.6 4.0 4.7 5.3 4.3 3.4 5.0 5.4 6.0 5.0 3.5 3.8 4.7
##  [415] 3.4 5.6 3.6 4.0 3.8 5.6 5.7 6.1 4.1 2.8 4.2 4.1 4.9 3.0 3.6 5.8 4.8 3.8
##  [433] 4.2 5.9 4.8 4.4 4.0 4.9 5.4 4.4 4.9 4.9 4.4 4.5 3.8 4.1 5.2 4.5 5.2 4.7
##  [451] 3.6 3.8 5.0 5.0 3.8 3.7 4.5 3.5 3.8 6.4 3.5 4.5 4.6 2.7 4.3 5.2 5.4 3.5
##  [469] 3.0 3.3 4.8 4.8 4.0 5.1 5.0 4.4 5.2 2.5 3.1 5.7 4.0 4.8 4.3 5.1 3.9 4.4
##  [487] 5.6 5.6 4.1 3.3 4.9 3.9 3.0 5.9 3.8 4.3 4.3 4.5 2.8 3.3 5.1 3.6 5.8 5.4
##  [505] 4.0 5.1 5.6 5.5 5.1 4.5 6.3 4.9 4.0 5.4 5.7 4.1 5.4 5.4 4.4 5.0 5.0 4.1
##  [523] 4.9 5.6 4.9 4.5 4.0 5.9 4.9 4.5 3.8 4.9 4.2 3.3 4.8 5.4 5.4 6.1 3.9 5.5
##  [541] 4.1 5.7 4.7 3.5 5.0 4.5 5.0 2.8 3.7 5.5 5.7 4.3 5.1 3.9 4.3 5.1 5.0 4.5
##  [559] 5.0 4.7 3.4 3.8 5.2 5.0 3.4 4.7 3.8 3.7 4.8 3.6 6.1 4.5 4.3 4.8 5.2 3.2
##  [577] 4.7 5.4 5.7 3.7 5.7 5.9 4.8 3.1 4.9 5.9 4.8 4.2 4.4 4.9 5.0 4.6 4.8 4.2
##  [595] 3.7 4.8 4.5 4.0 4.9 5.1 4.6 3.4 5.3 5.4 4.2 4.5 3.7 5.3 3.9 4.4 6.3 4.6
##  [613] 3.6 3.4 5.0 4.1 5.0 4.9 5.8 4.0 4.0 4.9 4.1 4.5 4.5 3.9 4.7 3.5 3.1 3.8
##  [631] 3.9 4.8 2.8 5.9 4.9 3.9 6.1 6.2 5.2 5.2 4.5 3.8 4.3 4.7 4.4 4.3 4.5 3.9
##  [649] 5.2 4.6 4.8 3.1 4.8 2.8 3.8 5.0 3.5 5.3 3.7 4.4 3.5 3.2 4.8 4.1 6.3 4.2
##  [667] 4.1 4.7 4.9 3.8 4.7 3.7 4.6 2.8 4.9 4.6 3.7 5.0 3.5 3.9 4.6 2.9 4.8 6.0
##  [685] 5.6 5.1 5.2 4.4 6.3 5.1 4.6 2.4 5.3 6.8 4.2 3.0 5.1 3.8 5.3 5.1 4.0 5.5
##  [703] 4.6 4.4 4.2 5.5 4.1 5.0 4.2 5.7 4.9 4.8 6.1 2.5 3.4 5.1 6.4 4.5 5.4 3.2
##  [721] 3.8 4.8 3.9 6.0 3.3 3.9 3.1 5.8 4.8 5.0 3.3 3.0 3.3 4.4 4.6 2.9 3.8 6.0
##  [739] 4.9 3.9 4.1 7.1 5.8 3.3 4.6 4.9 5.1 5.0 3.6 4.9 4.7 4.8 5.4 5.4 4.2 4.7
##  [757] 5.5 4.5 5.6 4.9 4.1 4.8 3.2 4.4 4.9 5.7 3.8 3.5 2.9 4.6 4.1 3.3 3.1 5.2
##  [775] 3.9 2.6 3.1 5.1 3.6 4.0 5.2 5.2 4.2 5.2 3.9 3.3 3.3 3.8 4.5 3.5 3.8 4.2
##  [793] 3.9 3.9 5.2 4.9 3.6 4.8 4.8 5.1 2.5 5.4 4.7 5.0 5.2 4.1 4.2 4.7 5.2 5.6
##  [811] 3.6 5.2 5.3 3.8 4.2 4.4 4.5 5.3 5.5 5.7 4.2 4.5 4.2 4.2 4.1 5.1 4.0 3.4
##  [829] 5.2 5.2 5.7 4.6 5.5 3.6 4.1 6.2 3.8 5.3 6.2 3.0 5.0 6.0 4.5 4.4 4.2 4.0
##  [847] 3.8 4.0 4.3 3.4 5.3 4.4 3.9 6.0 4.4 5.5 3.6 5.7 3.7 4.7 3.7 3.8 4.1 5.5
##  [865] 4.4 5.4 3.7 5.0 5.0 5.5 2.8 4.4 4.9 3.6 4.0 3.8 3.9 2.7 4.3 4.8 4.8 5.6
##  [883] 4.7 5.7 4.7 5.0 3.1 3.6 5.8 5.3 5.3 3.3 3.1 5.0 5.1 3.3 5.5 5.9 5.0 3.7
##  [901] 3.6 4.2 5.3 4.3 4.0 4.8 4.4 3.2 3.7 4.8 5.5 3.6 4.9 3.2 5.3 4.1 4.2 4.0
##  [919] 4.2 5.3 3.8 4.6 3.0 3.8 4.8 5.5 5.9 4.0 5.1 4.3 5.5 5.4 5.4 2.7 3.8 4.2
##  [937] 4.0 3.2 3.8 3.8 5.6 3.0 4.9 5.2 4.4 4.3 4.6 4.0 4.4 4.2 5.3 4.0 3.0 4.5
##  [955] 5.0 4.8 3.8 3.3 3.4 4.6 4.6 4.6 4.5 4.3 5.4 5.4 5.8 4.3 5.1 4.0 4.2 3.9
##  [973] 5.8 4.3 4.5 5.0 4.3 5.1 5.1 2.7 3.4 3.2 3.9 4.5 4.7 3.2 6.6 4.6 5.3 3.7
##  [991] 4.2 3.7 5.2 3.7 4.3 4.6 6.4 5.2 4.5 4.6
## 
## $func.thetastar
## [1] -0.0162
## 
## $jack.boot.val
##  [1]  0.531987578  0.389373297  0.317611940  0.075428571 -0.009171598
##  [6] -0.050746269 -0.138108883 -0.249858357 -0.409859155 -0.519335347
## 
## $jack.boot.se
## [1] 0.9779506
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
##    [1] 4.2 3.8 2.6 4.1 4.1 4.5 3.5 5.3 4.3 4.1 5.1 5.6 3.9 5.4 4.2 4.5 5.0 5.7
##   [19] 4.1 5.4 5.4 3.0 3.9 1.8 5.7 3.5 4.7 3.7 4.4 2.4 3.1 5.2 5.6 3.7 4.4 4.1
##   [37] 5.6 4.0 3.9 4.5 4.5 4.5 4.3 5.9 4.9 4.2 3.7 4.5 4.7 4.6 5.5 3.9 4.2 5.3
##   [55] 3.4 4.3 4.5 2.8 4.8 3.0 5.8 4.5 4.9 5.1 4.1 2.6 4.4 4.9 4.8 5.5 4.5 3.0
##   [73] 4.9 4.6 5.5 4.4 4.9 5.1 3.7 4.2 3.9 4.7 5.0 4.4 2.8 4.1 4.1 4.9 5.6 3.9
##   [91] 5.5 4.0 2.6 4.2 4.5 5.9 3.4 4.9 5.3 5.9 4.2 4.2 3.9 5.7 5.5 3.4 4.5 3.1
##  [109] 4.8 4.8 2.6 3.9 5.2 4.2 4.2 4.3 4.9 6.1 5.3 4.8 5.0 4.9 5.0 2.9 4.5 5.0
##  [127] 4.2 5.2 3.4 3.1 4.3 4.0 3.8 5.7 4.4 4.1 5.1 4.7 3.8 3.8 6.0 4.0 4.2 4.5
##  [145] 3.5 3.6 4.6 3.6 3.0 3.7 5.0 3.1 4.7 3.2 4.8 5.4 3.7 4.4 4.0 6.0 5.1 4.1
##  [163] 4.2 3.4 3.8 4.1 4.1 5.1 5.3 2.6 4.8 5.3 6.7 6.0 4.3 4.5 3.8 4.0 3.8 5.8
##  [181] 3.3 6.6 4.6 3.8 4.7 3.7 3.7 3.5 3.1 4.8 5.4 5.1 3.1 2.4 4.9 5.0 3.8 3.3
##  [199] 4.5 4.6 3.8 4.8 3.6 5.2 3.9 5.9 5.3 3.2 3.5 3.9 4.3 5.9 4.5 4.9 5.2 3.7
##  [217] 5.5 4.4 7.3 4.7 3.8 3.0 5.1 3.1 5.3 6.0 3.8 6.0 5.1 4.5 5.0 5.3 5.0 6.7
##  [235] 4.1 5.3 6.2 3.3 3.8 2.8 3.9 3.9 5.2 5.6 5.0 4.1 6.4 4.9 4.1 4.7 5.0 2.4
##  [253] 5.3 5.5 2.9 4.7 3.6 4.2 2.9 4.8 4.9 5.7 3.2 4.6 3.2 4.3 6.9 3.7 3.8 3.9
##  [271] 4.5 3.5 4.8 4.0 3.5 2.9 4.8 4.8 4.1 4.4 4.2 5.9 4.4 3.8 4.0 5.8 3.6 5.1
##  [289] 4.0 5.4 2.8 4.5 4.6 4.3 5.3 4.1 3.9 4.0 3.1 5.2 5.5 4.1 5.9 5.7 4.1 4.3
##  [307] 5.1 4.9 3.6 3.3 3.9 3.0 4.7 3.9 4.2 5.4 4.4 5.7 5.3 4.9 4.7 4.4 3.7 4.9
##  [325] 4.3 3.1 5.7 4.6 3.5 4.1 3.2 5.4 5.2 3.8 4.6 4.5 4.9 6.0 5.6 2.8 4.6 4.7
##  [343] 4.1 4.4 2.2 5.6 4.2 4.4 5.1 3.2 4.8 5.7 4.9 4.2 3.4 3.4 5.2 5.1 5.6 4.1
##  [361] 3.6 4.3 4.6 4.6 5.6 3.8 5.3 3.9 3.5 4.7 3.7 3.9 4.1 6.0 2.7 5.9 4.1 5.4
##  [379] 4.1 3.4 4.7 4.6 4.5 4.7 3.7 4.9 5.2 4.3 5.3 5.0 4.3 3.2 3.1 3.4 5.5 2.7
##  [397] 4.0 5.1 4.8 4.6 4.6 3.3 3.6 4.6 3.5 4.2 5.7 5.3 6.5 4.9 3.2 4.3 2.5 4.6
##  [415] 4.7 5.6 3.1 4.9 4.0 3.6 3.9 4.4 5.9 5.1 4.3 4.8 4.2 5.3 4.5 5.0 5.6 3.9
##  [433] 6.1 4.7 4.4 5.7 3.4 5.4 4.4 3.6 4.5 3.4 4.1 6.1 5.1 4.4 4.7 4.2 4.7 3.3
##  [451] 6.2 2.8 4.6 2.9 5.1 4.5 3.8 4.1 4.8 3.5 5.3 5.1 4.9 5.9 2.8 4.5 3.9 5.5
##  [469] 4.3 5.0 5.9 4.7 6.2 4.7 5.1 4.2 3.7 6.0 3.5 4.3 5.3 4.6 5.8 5.8 4.2 2.5
##  [487] 4.9 3.5 4.8 2.8 5.0 3.7 3.2 2.5 4.6 5.2 5.5 5.7 3.7 6.3 4.8 4.3 3.0 5.6
##  [505] 4.2 4.8 6.3 4.6 4.1 4.6 3.4 3.6 3.7 3.7 4.2 4.6 5.3 4.9 4.6 2.9 4.8 4.5
##  [523] 3.2 5.1 4.9 4.8 3.7 3.8 5.6 5.0 5.5 3.6 4.1 4.3 5.7 5.3 5.5 3.7 4.5 4.7
##  [541] 4.2 4.0 5.3 5.5 5.4 3.3 3.8 4.9 5.0 5.7 6.2 5.5 4.5 4.4 4.4 4.6 6.4 5.5
##  [559] 3.5 5.4 4.2 5.6 4.9 3.3 4.6 4.1 4.1 5.0 3.9 4.6 5.4 5.8 4.9 3.9 4.7 3.3
##  [577] 5.1 4.9 6.0 4.2 5.0 4.8 4.1 3.4 5.1 3.9 3.8 3.7 4.1 2.9 5.4 3.7 5.7 5.0
##  [595] 2.9 5.8 5.3 3.3 3.6 4.1 4.3 5.4 3.8 5.2 4.6 3.4 3.9 4.5 5.9 5.4 5.7 3.0
##  [613] 3.7 3.9 2.9 3.1 4.3 4.0 4.7 4.9 4.3 6.1 3.2 4.2 5.8 6.4 3.7 2.8 4.1 4.5
##  [631] 4.4 4.7 4.9 5.2 4.0 4.2 5.3 4.5 4.9 4.3 5.6 3.0 4.3 3.6 4.3 6.7 2.5 4.6
##  [649] 4.2 4.5 4.6 5.0 5.2 4.4 3.7 5.0 5.0 4.5 2.1 2.6 5.6 4.7 3.8 5.2 5.1 5.8
##  [667] 4.3 5.6 3.5 2.7 3.5 5.4 4.3 6.6 3.4 5.5 3.0 6.3 4.1 3.5 5.0 2.4 3.7 5.5
##  [685] 6.1 5.5 5.3 4.5 5.2 2.7 4.2 3.7 3.8 4.0 3.1 4.0 5.8 5.9 3.5 4.1 5.0 5.5
##  [703] 6.2 5.7 4.8 4.2 5.2 4.0 6.1 3.2 4.4 5.5 3.6 3.6 5.7 3.8 5.8 3.9 4.6 5.5
##  [721] 4.7 3.7 2.9 4.0 6.2 6.1 3.3 4.2 3.5 5.7 5.5 4.1 6.1 4.7 4.3 4.2 3.3 5.5
##  [739] 5.1 5.0 5.1 4.7 4.5 4.2 4.4 3.8 4.4 3.3 5.6 5.4 5.4 6.6 2.1 3.9 5.5 5.4
##  [757] 4.7 3.5 5.7 3.1 4.7 3.9 5.9 3.9 4.5 5.7 5.5 5.1 4.2 5.1 5.0 4.9 4.2 5.3
##  [775] 5.6 4.7 4.7 4.6 4.0 4.5 5.3 6.1 5.5 5.7 4.2 3.7 3.9 3.6 5.1 3.3 3.8 5.0
##  [793] 4.6 3.8 4.3 2.0 5.0 5.5 4.3 3.0 2.8 4.8 3.5 4.2 4.7 5.6 6.6 5.1 3.6 4.5
##  [811] 4.7 5.5 3.4 5.2 3.4 4.5 2.9 4.5 3.6 4.7 2.8 4.9 4.5 4.2 3.7 3.5 5.2 3.2
##  [829] 3.5 3.9 3.8 3.7 4.9 4.5 4.2 4.9 4.6 3.3 4.6 4.5 5.6 3.1 4.7 4.3 3.7 6.1
##  [847] 4.5 3.9 2.2 6.0 3.7 3.6 2.8 5.0 5.1 6.3 4.4 3.0 3.8 5.7 5.6 6.2 6.2 4.5
##  [865] 4.1 5.0 5.6 4.4 3.0 2.4 3.0 4.5 4.3 2.9 5.2 4.3 5.5 4.6 5.5 3.9 5.0 6.3
##  [883] 4.8 3.9 4.6 4.6 3.6 4.0 3.5 5.0 6.0 3.7 5.9 3.1 3.8 4.7 4.6 5.7 5.3 3.3
##  [901] 4.7 5.1 5.6 4.0 4.9 4.2 4.3 5.0 4.2 4.1 3.9 3.3 5.6 5.0 4.8 4.9 5.8 4.8
##  [919] 6.3 5.0 3.9 4.9 4.2 4.6 3.8 4.3 5.0 5.1 3.2 4.4 4.0 6.3 4.2 4.2 3.1 4.3
##  [937] 3.8 4.2 4.0 2.6 6.6 3.5 5.3 4.2 4.7 4.1 3.3 4.3 4.2 5.1 3.5 4.6 5.4 4.3
##  [955] 3.7 3.3 3.8 4.9 3.3 4.3 5.9 3.1 3.5 5.0 5.8 5.6 4.2 3.5 4.0 5.6 5.6 4.7
##  [973] 5.7 3.5 4.9 3.9 4.2 4.9 4.8 5.0 5.1 6.0 5.4 4.2 6.9 6.1 3.0 4.4 3.1 5.5
##  [991] 6.3 4.4 5.0 3.7 3.1 3.3 6.2 5.4 4.3 5.5
## 
## $func.thetastar
## 72% 
##   5 
## 
## $jack.boot.val
##  [1] 5.600 5.400 5.400 5.300 5.264 4.944 4.900 4.700 4.700 4.500
## 
## $jack.boot.se
## [1] 1.053297
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
## [1] 1.057922
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
##   5.876498   9.153334 
##  (2.556931) (4.157860)
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
## [1]  0.5879761  0.3626036 -0.6378583  1.5214438  0.2570424  1.9030374
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
##    [1]  1.4285888599  0.7915836930  0.3468592545  0.6565424745  0.7159270030
##    [6]  0.7735869345  0.3045928836  1.6289806284  0.4793420451  1.3353607436
##   [11]  1.4871081016  1.9366375106  0.2121152189  0.3688622052  1.6556055489
##   [16]  2.0299771473  1.1612723557  0.3653510419 -0.3768904427  1.1225753866
##   [21]  1.7437914539  0.4813048354  0.3348607993  0.2283295403  0.7893421307
##   [26]  1.2153237221  1.0225748693  1.2870710950  0.2516031311  1.2223630076
##   [31]  0.6194238201  0.7066567961  1.3313342298  0.2318565453  0.6373501799
##   [36]  0.4372882529  1.4474126317  1.7900403807  1.0656808081  0.3035527281
##   [41]  1.9675116860  0.6521788228  1.0683563519  1.2025991133 -0.0663281383
##   [46]  1.0805220404  1.1648430976  1.4983491953  0.5684397270  1.0949053282
##   [51]  0.5411903599  0.0090377548  0.3517826314  1.6263447807  1.0841898373
##   [56]  1.2782448086  1.8677167420  1.1939558299  1.0119618284  0.9045739811
##   [61]  0.2172987370  2.1892791740  0.7033335767  1.1361712731  0.5902125103
##   [66]  0.2866631652  1.6802167331  2.0011808687  1.5620058008  0.3202926334
##   [71] -0.1978267795  1.0243750426  0.7730584125  0.6220910456  1.0935962184
##   [76]  1.6944410235  2.2184509087  1.9531824882  0.1833811817 -0.4962935162
##   [81]  0.3115293260  0.7169095555  1.3253340592  1.5561755407  1.5817789242
##   [86]  0.8917689603  1.6661807125  1.3533651998 -0.1668688029  0.6918591245
##   [91]  2.0623054062  1.3967797332  1.2455612737  1.5657603134  0.3092182522
##   [96]  1.3757552023  0.5803113622  0.6857190868  0.6135699467  1.0675822365
##  [101]  0.1520718500  0.3766847393  1.5640071942  0.3086774151 -0.1588213953
##  [106] -0.3014373802  0.1204266730  1.7421165390  2.2117906109 -0.3996975557
##  [111]  0.7080758026  0.1599453585  0.4751571537  2.0582684738  1.5326879773
##  [116]  0.7158294793  0.0655875838  0.3040415125  1.6846262899  0.6553259225
##  [121]  1.1536574817 -0.0435909925  1.2887666861  0.9123148091 -0.0222620668
##  [126]  0.5322705766  0.7705987464  1.7385114407  1.1384736155  0.6023987751
##  [131]  0.9167130447  0.1800528162  1.7654634764  1.9466728815 -0.4084829519
##  [136]  1.2319022342  0.4717606843  0.4529766850  0.6930591915  1.2310175204
##  [141]  0.6926179578  1.2772976109  0.8668654898  0.6555691080  0.6185547627
##  [146]  1.2261570302  1.3041407213  0.3337751648  0.8792481180  0.4008670752
##  [151] -1.0621513038  1.5929714978  1.3879209999  1.6784775681  1.2299802710
##  [156]  0.1683666405  0.6961652272  0.7416356174  1.0403971142 -0.7191785457
##  [161]  1.3181814784  1.0766801233  0.5106841838  1.1756226290  1.9289866676
##  [166] -0.4795055935  0.8836933874  1.5167997617  0.7472103385  0.2459318221
##  [171]  1.7490927940 -0.6064761012  0.7684999439  0.9323599626  1.0800896309
##  [176]  2.0474421629  0.9063017454  1.3338514711  0.7457919744  1.2370201023
##  [181]  1.7700749213  0.7123114931  0.6301874234  1.0306360244  0.6636777449
##  [186]  1.6754406992  1.3019658160  1.7440772178  1.5522706280  0.3506999721
##  [191]  1.9670536325  0.3181700683  0.9101177658  0.5953630746  1.8162504007
##  [196]  0.7626651786  1.8478740217  0.9932464253  0.4028633026  1.1531114189
##  [201]  1.1965635913  1.1387682356  2.0455880440  0.6261428293  0.1419473322
##  [206]  1.7246584991  0.6644597039 -0.5433438161  0.7431770194  1.2534857813
##  [211]  1.7795762385  1.1585666122  1.0967493994  0.7913940006  1.2800879156
##  [216]  1.5395719765  0.5075616325 -0.0777141434  1.8950301980 -0.0040656671
##  [221] -0.0273439681  0.7407447802  1.6057436835  1.1921991867  0.3344765008
##  [226]  1.1387682356  1.2858499468  0.8838313838  0.9345416169  0.0396284284
##  [231]  0.7554770751  0.5411903599  0.5272379038  0.9640916155  1.3482149847
##  [236]  0.3536581267  1.0015755185  0.5460626030  1.0682843782  0.6727990410
##  [241]  1.0860959078  0.8167158324  1.0585589600 -0.1972557976  2.4707393795
##  [246]  0.7675759181  0.9031233773  1.3781615659  0.2428866594  1.1624233700
##  [251]  0.9403887484  0.4300745310  0.5583743740  0.2243206011  1.2398899681
##  [256]  0.9805005978  0.7314664587  1.2848587393  0.8080306580  1.0982592961
##  [261]  1.2682483893 -0.5374663468  0.2677941320  1.0998032914  0.8575835084
##  [266]  0.5771148398  1.1821109391  1.2863074521  0.9981008932  1.0781540300
##  [271]  2.0546461682  1.1585999571  0.6497941401 -0.0865446163  2.2107926648
##  [276]  0.4694545138 -0.2664452148  1.0854976011  0.7630654944 -0.4169069132
##  [281]  0.0550416797  0.8834380518  0.4054827808  0.9305060119  1.1210791230
##  [286]  0.9328145449  1.3887872597  0.6280051581 -0.8800236901  1.8180767771
##  [291]  0.3484511693  0.8020965020  2.1599040324  0.5790197683  0.0440413924
##  [296]  0.6262959282  0.4212522810  1.0877859842  0.4424294259  1.3115854753
##  [301]  0.4441127709  1.6456126626  1.3108419061  1.7636886546 -0.6659442570
##  [306]  1.1263183689  0.6570503902  1.2461699835  0.9119744886  0.5332168114
##  [311]  0.1703312960  1.0501425893  1.4288671791  1.3547187769  0.6518877395
##  [316]  2.3099768076  1.7330850990  0.6038226367  0.8932733484 -0.0004623768
##  [321]  1.6868632932  1.0967029415  1.1106392426  0.5854181120  0.0695137783
##  [326]  0.6118499391  1.1556681964  0.6534910958  0.8737987342  0.9361067669
##  [331]  1.5069511366  0.6184006051  2.0042704541  1.5657603134  0.1128063217
##  [336] -0.2772639063  0.7160428143  2.0958085649  0.5349971957  0.6126508586
##  [341]  1.6580822278  0.9925182855  0.8831230413  1.1528199896  0.6010528941
##  [346]  2.0004062943  1.1286475356  0.9327041356  0.8180306460  0.4167864518
##  [351]  0.2742321343  2.0515689023  1.1673614335  1.9906827531  1.1907894652
##  [356]  0.1559326786  1.1594628805  0.5055114193  1.3369154552  1.1742875265
##  [361]  1.5217336362  1.1013367598  1.9096168491  0.3599598980  1.1964940163
##  [366]  1.2184961516  1.0646469912  1.1675946818  1.4237884292  1.2592671804
##  [371]  1.9149402123  0.7963101513  2.0659381043  1.2161366510  0.4669400494
##  [376]  1.3393464808 -0.4541550059 -0.1207262821  1.3624749458  1.3207478928
##  [381]  2.1929508104  0.3924778456 -0.2815447040  1.0920264399  1.2191459647
##  [386]  0.7919984865  1.2744518121  1.1900742252  1.7497297242  0.2488926879
##  [391]  1.4508781328  0.7242377372  0.7019191807 -0.2423693017  1.2315041054
##  [396]  1.0003458587  0.2862147266  0.7489518915  1.3061285805  0.5388024272
##  [401]  0.8657292226  0.6071952632  1.4225868539  0.3591963498  1.4538635390
##  [406]  1.0667023648  0.6910116673  0.1682518404  1.7638765547  0.9904337151
##  [411]  0.3823446108  1.5555082236  1.0280949189  1.1756065739  1.1916450663
##  [416]  1.8238591812  0.7684575194  0.9710270645  0.3597698995  0.7566003783
##  [421]  1.0204836175  0.6555691080  1.3544432990  0.7359331983  1.9148324787
##  [426]  0.7041542702  1.5347508695  0.9787338273 -0.6556648421  0.7811037760
##  [431]  0.6411914904  0.2248558646  0.8348162875  1.4082194402  0.1712438559
##  [436]  1.0921337434  1.5931484307  0.1237740513  1.3035993887  1.7533791340
##  [441]  0.2671146766  1.0469258825  1.1861746781  1.0655453336 -0.5382503556
##  [446]  0.9286874722 -0.1739864257  1.0271366966  0.8761898453  1.1364832321
##  [451]  1.0383553825  0.6747872964  0.6691905584  1.0418963261  1.6839647717
##  [456] -0.0650080451  0.1134647836  1.9047458014  1.0018497201  1.2727566702
##  [461]  1.0479028258 -0.0231561707  0.6440413989  1.5774573961  1.5756017852
##  [466]  0.8792516992  1.3591003374  0.8254004192  1.4589295643  0.8434764935
##  [471]  1.6427700409  1.3648703901  0.7651697396 -0.3502291764  1.1808948992
##  [476]  0.3092071977 -0.4079106163  1.5584427749  0.4859898686  1.7285285133
##  [481]  0.3053457104  1.1406424684  0.5690898842  1.1538817281  1.8795748533
##  [486]  0.9960103820  0.9961896222  1.6370202576  1.5666375465  1.6913811048
##  [491]  1.1080117700  0.7160363804  0.6747602668  0.1500166262  1.2577172517
##  [496]  0.8005518479  2.0725308100  0.9693629448  0.4306838029  1.9215638146
##  [501]  0.7605330039  1.0985994574  1.6203376705  1.3743661678  0.3540717910
##  [506] -0.8495762682  0.6804143592  0.5559866502  1.7828072425  0.3814563804
##  [511]  0.9380814024  1.9142936567  0.6868958683  1.2172254825  0.0723137987
##  [516]  0.6281671380  0.8837786130  1.8111112492  0.9886369337 -0.0584637652
##  [521]  0.1966759413  1.2920790965  0.5374155314  1.7187913169  0.6636777449
##  [526]  0.9745569717  1.8128821961  0.4753639727  0.6898044787  1.1717276914
##  [531]  0.5004167816  1.2281059712 -0.1509234547  0.6466497785  1.0971339226
##  [536]  1.8012017255  0.9905450576  1.1323764560  1.2093254577 -0.0476620804
##  [541]  2.2348961648  1.9450213092  0.8216959973  1.7836187360  0.6522515778
##  [546]  1.2281059712  1.5454512998  1.3174449484  0.7076164577  0.5933569041
##  [551]  0.8489216256  1.1818466893  1.3722142257 -0.5416202262  1.0901434535
##  [556]  0.3972613581  0.7110602903  0.1962958234  1.4897529862  0.4494437055
##  [561]  0.8176508769  0.3711406555  0.2761780272  1.9165669050  0.8009612514
##  [566]  0.3237941304 -0.0701258393 -0.2387890046  1.6406383685  0.5334771876
##  [571]  0.5070779263  0.9387368837  0.5920470104  1.3393703585  0.6563201821
##  [576]  0.9537089832  0.3104181981 -0.3320266351  0.1845507779  0.4335663193
##  [581]  1.0113698064  0.5834556952  0.2867860628  1.2971728202  1.6201168655
##  [586]  0.9378874899  0.7649534813  0.2987994262  1.2044568014  0.8756900283
##  [591]  1.9682543936  1.0933438987  1.4115489457  1.1159753137  0.9397995122
##  [596]  1.9238713308  1.6177461118  1.3238931123  1.3226416616 -0.1780737164
##  [601]  0.2874516285  1.1304292800  1.5754051890  0.3917144958 -0.0413032213
##  [606]  1.2524083888 -0.0346838887  0.2817732362  0.8724420739  0.7168063143
##  [611]  0.8874706264  0.6526823951  0.7037350304  0.0150717647  1.2876753329
##  [616] -0.6026462786  1.2992696883  0.9248718575  2.1741465219  1.2460871158
##  [621]  2.3038436833  1.2846504809  0.3307211665  0.8727480717  0.2483877982
##  [626]  1.6409947042  0.3409514656  0.8385936563  0.6999672984  1.1860871515
##  [631]  1.4384641811  0.8666883341  1.2299451346  1.1218827005  1.1190840404
##  [636]  0.7995495457  1.7670085983  2.0818939632  1.3209350702  2.0275972634
##  [641]  0.6743401267  1.5144457785  0.5529641421  0.5571807165 -0.5396887565
##  [646]  0.8158287288  0.7789137246  0.1212109633  0.5198390375 -0.1924531389
##  [651]  1.0194495096  0.1885423590  0.5311161498  1.3040839771  0.7746327085
##  [656]  0.3872192323 -0.2205983515  0.9210322082  0.2487828082  0.6850882709
##  [661]  1.9434724281  0.4929567069  1.2159320177 -0.5681435552  0.5354883248
##  [666] -0.0921814600 -0.3194461528  0.9348301482  1.4604418883  1.4878190497
##  [671] -0.0604021594 -0.6680835626  0.6460341552  1.0767596194  1.5214411569
##  [676]  2.2097658246  0.2702134092  1.0520280033  1.2208713090  1.3051596526
##  [681]  1.2948734950  1.6888815611  1.8099085297  1.0919955756  0.6083247318
##  [686]  1.0658694564  1.1755402242  0.6067009024  1.2592615219  0.7652573934
##  [691]  1.8672434293  1.1842145006  1.0143322068  1.0212185708  0.6352578078
##  [696]  0.3102939324  1.1585999571  0.3440454362  0.9833721281  0.2860258541
##  [701]  1.6702104629  1.0608419355  2.0561561932  0.5253681604  0.6103436344
##  [706]  1.2290783776  0.3859531065  0.1730466158  0.1639754675  0.6891850932
##  [711] -0.5246147391  1.7136338769  2.2498456417  0.4744525946  1.5313084862
##  [716]  0.7519151570  1.2193902786  1.3352992461  1.6895414158  1.6145067978
##  [721]  1.3078711465  0.6476120341  1.8965967886  0.5901809076  1.4274094895
##  [726]  1.3072435824  0.5204094767  1.0077448704  1.0125859397  1.2503042652
##  [731]  1.2306719464 -0.1729884283  1.9542829016  0.2872530800  0.1703312960
##  [736]  0.7958126860  0.9861489356  0.2671725828  0.1179001304  0.7211661120
##  [741]  0.7654339522  1.2809581616  0.3581906899  1.1359686119  1.2022195198
##  [746]  1.6043091762  1.0698758748 -0.8759602598  1.6480843512  1.1156086130
##  [751]  1.2440897593  0.7760834474  1.3154358021  0.7240340784  0.2751072406
##  [756]  1.0095502510  1.1424465236  1.0225748693  0.5709875022  0.9178074842
##  [761]  0.8587957881  0.2480611225  0.2434511597  2.0649738110  0.6304776820
##  [766] -0.1763590304  0.0797794782  0.1512232470  1.8232280123  0.0990975305
##  [771]  0.9684778201  0.7804412100  0.5318952637  0.6438667319  1.1069992783
##  [776]  0.2056980907  0.6514618452  0.1437319326  1.1436080332  1.1900806780
##  [781]  0.6960125939  1.5522706280  0.6825065817  0.8466571617  0.4244848033
##  [786]  1.1888784966  0.9697959406  1.0307384460  0.9210898516 -0.5581703354
##  [791]  0.7435629802  0.8591727849  1.2296769885 -0.0327096472  1.0792107565
##  [796]  0.8849957631  1.8537972096  1.3537089207  0.4719392928 -1.5605266859
##  [801]  1.5545984731  1.3156309401  0.2999018397 -0.3829428979  0.6364108266
##  [806]  0.8725496848 -0.1294886340  2.1648074462  1.2352375326  0.5743568719
##  [811]  1.1984881612  1.6658866765  0.3712379952  0.3247546481  1.3302833609
##  [816]  0.6039456379  1.3246945044  0.6536685129  0.3219580677  0.7407813628
##  [821]  0.3983555423  1.4991029940 -0.5524687569  0.7769681690  2.2791969334
##  [826] -0.0302994366  1.1234224138  1.4819342904  1.3616299192 -0.2102490984
##  [831]  0.3942661715  1.4095591409  0.7860346367  1.3524117400  0.7430813863
##  [836]  1.4129968186  1.4176846089  0.3823790758  1.7008432151  0.6644875308
##  [841]  1.0791710496 -0.1818892480  0.2165450950  1.3917410152  0.3712890300
##  [846]  0.4753019881  0.8894901346  1.7917121284  0.7300044390 -0.2718627816
##  [851]  1.9757733837  0.8602904240  0.7099005938  1.3420872697  1.8418930640
##  [856] -0.4430438724  0.6415873263 -0.3521004037 -0.0032361360  0.8093844220
##  [861] -0.1485592029  0.6397174488  0.6857489280  0.6578660623 -0.7661054588
##  [866] -0.1596375565  0.9132884879  1.1573705989 -0.0061076176  1.1890908268
##  [871] -0.0114525477  1.8335206104  1.8585147176  1.1261639878  0.5984599496
##  [876]  0.8518446171  1.1694657082  1.8128010889  1.6924098001  1.0755975264
##  [881]  1.9896864084  1.4814002218  0.1122631813  0.8870375260  1.8956312936
##  [886]  1.8022941662  1.5640071942  1.3048669390  1.3216971518  0.7122084735
##  [891]  0.8091295306  0.5428546426  0.0261238213  0.1395824697  0.7658505266
##  [896]  1.3949536533  2.0141565507  1.5507488628  2.2171017525 -0.0154245153
##  [901]  2.1270337403  0.2423280566  0.8309223246 -0.3102321237  1.6426457422
##  [906]  0.7220612673  0.5399739187  0.6866280749  0.2484244746  1.1575498508
##  [911]  1.3258896829  0.0401020914  1.2239182126  0.7915366938  1.1276756614
##  [916] -0.1429741327  1.3409780971 -0.1683758543  1.2281956320  0.7561645551
##  [921]  2.3840160529  0.5204094767  1.7572544140  1.7764259933 -0.2717285796
##  [926]  0.7197681352 -0.3155281354  0.7467891890  1.8434400521  1.2203670605
##  [931]  0.5851489400 -0.3136528651  1.5209124614  1.2616270065 -0.2476547308
##  [936]  1.1542761642  0.5870610553  0.6083768374 -0.0427434990  2.1900505572
##  [941]  1.0350560537  0.6180226468  1.2980966393  0.2836566862  1.8009346140
##  [946]  1.5685850106  0.5799832996  1.3547187769  1.1423126167  1.5548178207
##  [951]  1.2513428658  1.5380046947  0.5029334272  1.2766811427  1.6596104838
##  [956]  0.6467751579  1.7490927940  0.7329574834  0.8423768496  0.9293738704
##  [961]  2.2016446180  1.2681539844  0.6466103077 -0.1535701036  0.7394532743
##  [966]  0.4276318890  1.0629145513  0.6944984702  0.7843119101  1.1606210214
##  [971]  0.5597783162  2.2123971780  2.0853622841  1.1472563179  0.2422641420
##  [976]  0.2434329602  0.5715886911  1.9540490691  0.8601128557  0.3892647089
##  [981]  1.3889164797  0.8153688615  0.9904337151  0.8627948496  0.1093407245
##  [986]  1.6486488924  1.0809963857  1.1362439759  0.8414291588  1.2816617652
##  [991] -0.6647766120 -0.2815504896  0.9030080773  1.4684686601  1.2787426473
##  [996]  0.6486892614  0.5870640606  1.1340988029  0.4158618322  1.3933337260
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
##   0.64200863   0.28647901 
##  (0.09059262) (0.06405534)
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
## [1] -0.7862117794  0.0734955117  0.3335817669  0.0001333472 -0.6958068487
## [6]  0.3454233552
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
## [1] -0.0131
```

```r
sd(results2$thetastar)  #the standard deviation of the theta stars is the SE of the statistic (in this case, the mean)
```

```
## [1] 0.888306
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
## t1*      4.5 -0.03003003   0.9203767
```

One of the main advantages to the 'boot' package over the 'bootstrap' package is the nicer formatting of the output.

Going back to our original code, lets see how we could reproduce all of these numbers:


```r
table(sample(x,size=length(x),replace=T))
```

```
## 
## 0 2 4 5 7 9 
## 1 3 1 3 1 1
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
## [1] 0.008
```

```r
se.boot
```

```
## [1] 0.9142384
```

Why do our numbers not agree exactly with those of the boot package? This is because our estimates of bias and standard error are just estimates, and they carry with them their own uncertainties. That is one of the reasons we might bother doing jackknife-after-bootstrap.

The 'boot' package has a LOT of functionality. If we have time, we will come back to some of these more complex functions later in the semester as we cover topics like regression and glm.

