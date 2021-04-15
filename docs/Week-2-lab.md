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
## 0 1 3 4 5 7 8 9 
## 2 1 1 2 1 1 1 1
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
## [1] -0.0169
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
## [1] 2.749295
```

```r
UL.boot
```

```
## [1] 6.216905
```

Method #2: Simply take the quantiles of the bootstrap statistics


```r
quantile(xmeans,c(0.025,0.975))
```

```
##   2.5%  97.5% 
## 2.8000 6.1025
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
##    [1] 5.5 4.7 4.0 3.8 4.9 6.0 3.9 4.5 5.2 4.3 4.4 4.6 3.0 5.4 4.0 6.4 5.1 4.3
##   [19] 3.1 4.8 3.5 4.9 4.7 5.4 4.5 5.1 6.0 4.7 3.4 5.4 5.2 4.7 3.5 5.1 4.8 5.5
##   [37] 4.3 2.8 3.9 4.4 4.9 3.8 4.5 3.9 5.5 4.9 4.4 5.7 3.1 3.3 4.2 5.4 4.3 3.5
##   [55] 5.5 3.6 4.1 5.6 4.4 4.8 4.2 4.0 4.9 5.6 4.8 4.4 2.0 2.9 3.1 4.6 4.9 4.8
##   [73] 6.5 4.7 5.2 3.1 4.6 5.0 4.1 5.0 6.0 6.5 4.6 1.9 4.4 4.6 3.8 4.4 4.4 4.3
##   [91] 4.2 4.1 4.4 6.2 4.4 5.9 4.7 3.9 4.4 6.2 4.7 5.3 4.4 3.3 4.7 2.6 4.6 5.0
##  [109] 4.4 5.3 5.1 3.9 4.2 4.9 5.4 6.2 3.5 4.8 4.1 4.7 4.7 5.3 4.2 4.5 5.2 5.4
##  [127] 5.1 4.8 5.3 4.0 4.1 4.6 5.4 4.7 4.8 4.2 4.7 4.9 5.7 3.3 3.9 4.3 6.2 5.3
##  [145] 5.2 4.0 5.5 4.3 5.7 5.0 4.4 4.9 5.1 4.1 4.5 4.9 5.1 4.3 4.8 4.7 4.8 5.1
##  [163] 4.5 5.6 4.9 3.6 4.3 5.0 5.6 4.4 4.5 3.9 4.6 4.1 4.5 5.7 5.3 5.1 4.1 5.4
##  [181] 4.2 3.3 4.4 4.9 5.5 4.5 6.4 4.5 4.1 4.4 4.1 3.6 5.0 6.4 4.2 5.2 3.6 3.2
##  [199] 3.6 3.9 6.7 4.4 4.3 5.0 5.4 2.9 5.3 3.4 4.8 4.4 2.3 3.4 5.8 4.3 3.4 5.5
##  [217] 4.6 4.0 3.5 3.6 4.2 3.9 6.1 5.5 5.2 3.4 6.1 4.7 5.2 5.4 3.6 6.2 3.9 3.8
##  [235] 5.5 3.8 5.0 4.4 4.2 4.0 2.7 3.3 4.4 4.4 4.9 5.6 5.0 3.5 5.6 4.2 4.2 5.6
##  [253] 4.8 4.1 4.6 4.3 4.8 5.8 3.9 5.8 3.4 6.0 2.8 4.0 3.4 4.5 3.6 4.8 5.7 3.4
##  [271] 3.6 4.3 5.8 3.7 5.7 4.5 4.7 5.1 4.4 5.6 6.6 4.9 5.4 4.7 3.7 3.3 5.5 4.2
##  [289] 3.3 5.3 3.3 5.9 3.8 3.4 2.6 4.5 4.9 3.9 4.5 4.7 4.8 4.8 2.7 5.5 3.9 3.8
##  [307] 5.6 4.2 3.7 4.3 5.2 5.7 4.6 4.6 4.3 3.4 5.6 4.4 5.0 4.4 5.4 4.4 3.8 5.3
##  [325] 5.4 4.1 5.1 3.4 5.4 3.4 4.4 5.3 4.7 4.8 5.4 5.2 4.2 3.2 5.2 4.7 3.2 4.1
##  [343] 4.0 4.8 4.3 4.5 5.8 5.0 4.6 4.4 2.3 4.7 5.8 2.6 5.3 3.5 3.9 4.8 3.9 4.9
##  [361] 4.2 5.1 4.9 5.6 5.1 4.1 6.3 4.9 3.1 4.2 4.3 3.1 5.3 4.7 5.5 4.5 4.6 4.0
##  [379] 4.3 4.8 4.7 4.6 4.7 3.1 3.3 5.2 3.3 3.1 4.6 2.6 3.9 5.0 5.2 4.9 3.3 5.8
##  [397] 5.2 5.2 5.4 2.9 4.2 3.7 3.9 5.2 4.4 3.6 5.3 4.0 4.2 5.0 4.0 5.3 5.0 4.7
##  [415] 4.0 3.2 6.1 2.7 4.7 4.6 5.4 4.2 6.0 4.9 4.1 4.5 2.6 4.4 4.0 3.9 5.0 2.9
##  [433] 3.5 4.9 4.1 4.1 6.2 5.5 4.5 5.2 4.2 3.2 3.8 4.3 4.9 4.3 5.2 2.4 3.5 2.1
##  [451] 3.8 4.5 5.6 5.1 4.7 3.8 4.6 4.1 4.2 2.7 3.7 3.7 4.7 5.5 4.9 3.9 2.5 4.5
##  [469] 3.4 4.5 4.1 4.3 4.0 3.5 4.9 5.0 4.5 4.5 4.5 5.1 3.1 3.9 4.6 4.4 5.4 5.3
##  [487] 5.0 5.0 4.2 5.5 5.9 4.8 4.4 3.4 6.3 2.8 4.9 2.9 4.2 5.4 4.5 4.9 6.4 5.0
##  [505] 4.1 3.7 4.9 3.6 4.6 3.6 3.7 4.0 4.0 3.9 5.7 4.7 5.4 4.3 5.3 4.4 4.9 5.4
##  [523] 5.3 4.1 3.3 3.3 2.7 6.1 3.2 3.9 3.9 5.9 3.4 4.6 3.3 5.0 4.9 4.2 4.4 4.5
##  [541] 4.6 5.0 5.7 4.8 4.8 5.0 5.6 4.6 4.0 4.6 4.6 4.5 5.5 5.0 4.0 4.5 3.0 4.5
##  [559] 5.2 4.4 3.3 3.7 4.4 3.9 3.5 4.7 5.3 5.2 5.6 5.1 4.2 4.6 4.6 3.4 4.5 4.2
##  [577] 3.7 4.3 4.8 6.7 5.3 3.7 3.4 4.0 3.9 5.0 3.9 4.7 5.2 5.0 4.0 3.8 4.0 5.4
##  [595] 4.3 4.3 6.5 5.4 4.5 4.2 3.2 5.6 6.7 4.9 4.9 4.2 5.4 4.5 5.2 5.6 4.0 6.1
##  [613] 4.8 5.5 6.0 3.8 3.2 4.2 4.2 5.3 3.9 3.9 5.9 5.4 3.8 5.2 3.1 4.6 4.0 4.5
##  [631] 5.8 4.1 5.7 4.3 4.2 4.5 5.0 4.1 6.7 5.1 4.6 3.6 5.2 4.4 3.8 3.5 6.0 5.5
##  [649] 4.8 4.5 3.8 3.2 4.5 3.3 4.9 4.3 3.1 3.9 4.4 4.4 4.9 4.2 4.0 5.8 5.9 4.0
##  [667] 4.4 3.8 3.8 4.4 5.5 4.6 4.7 4.8 3.8 2.8 5.3 4.8 4.9 4.0 5.6 4.7 4.3 5.6
##  [685] 4.2 2.9 4.5 4.6 5.7 3.0 6.1 4.4 3.8 5.1 4.4 4.5 5.9 3.1 3.6 2.1 2.9 2.5
##  [703] 5.2 3.3 5.0 5.6 3.9 2.0 4.8 4.0 4.5 4.6 4.3 5.3 2.5 5.2 2.7 4.0 5.7 4.1
##  [721] 6.4 5.3 4.0 4.1 4.1 4.7 4.2 3.7 2.8 4.5 4.3 4.8 3.6 5.6 3.5 5.4 5.1 3.9
##  [739] 4.8 6.0 4.5 5.7 3.8 4.9 6.0 6.3 4.2 4.2 2.6 3.0 3.9 5.5 4.2 4.5 5.0 3.8
##  [757] 5.2 5.5 4.7 4.4 5.8 3.4 3.0 4.2 4.9 2.6 4.6 5.0 2.6 5.3 4.7 4.2 5.8 3.1
##  [775] 4.7 4.6 3.5 3.2 4.5 5.8 5.2 4.3 6.0 3.7 5.5 5.7 5.1 4.8 4.7 4.6 3.4 4.4
##  [793] 3.8 4.5 5.4 4.9 4.3 3.4 4.4 6.2 3.0 4.5 4.7 3.2 4.7 2.5 4.0 5.0 3.6 4.7
##  [811] 4.1 5.0 5.7 5.7 5.7 4.0 5.4 4.7 4.7 4.7 4.9 5.9 4.6 3.8 5.8 6.3 3.9 5.0
##  [829] 3.3 4.4 4.6 5.8 2.8 4.2 4.8 2.8 4.0 5.1 5.5 4.6 3.9 5.1 3.5 4.3 5.2 5.7
##  [847] 4.6 3.2 4.9 4.8 4.6 4.0 4.2 4.8 4.7 4.9 5.1 3.8 4.6 3.7 4.5 3.3 4.8 5.4
##  [865] 4.5 3.9 5.4 5.5 5.2 4.2 5.9 4.5 4.6 5.1 3.8 6.5 4.5 4.0 4.7 5.0 4.6 3.0
##  [883] 3.8 4.7 5.7 3.7 5.2 5.2 3.5 4.0 4.9 4.4 4.4 3.4 5.1 5.7 5.3 4.8 4.6 2.8
##  [901] 2.8 4.9 5.2 6.2 4.3 5.6 3.4 5.8 6.0 4.9 3.5 6.1 5.3 5.1 3.9 4.8 4.6 4.5
##  [919] 5.2 3.9 4.1 3.4 4.0 5.5 5.4 6.2 5.9 3.6 4.1 3.1 5.4 4.4 3.9 5.5 4.2 4.6
##  [937] 5.6 4.8 4.8 4.4 4.9 4.6 4.6 5.1 3.8 6.2 3.9 4.3 5.5 3.1 4.7 6.0 3.6 3.3
##  [955] 3.6 4.3 5.3 4.0 5.5 3.8 2.1 4.8 3.5 4.5 5.3 2.8 3.5 4.3 5.6 4.9 3.7 3.5
##  [973] 4.9 5.0 3.8 5.6 3.2 3.5 2.3 3.7 4.2 3.7 5.7 3.4 4.9 5.5 3.8 3.6 6.4 4.9
##  [991] 3.4 4.8 5.3 5.0 4.4 6.0 4.7 5.6 4.4 4.6
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
##    [1] 4.8 4.9 2.7 3.0 5.1 3.8 5.0 4.4 3.1 3.1 4.9 3.5 5.3 4.3 3.4 2.8 3.8 5.4
##   [19] 5.2 4.9 4.1 3.8 5.3 4.6 5.6 4.9 2.7 5.8 4.6 4.7 4.0 4.8 3.6 6.9 4.7 3.9
##   [37] 5.4 4.1 5.3 5.2 6.0 3.0 4.8 5.5 5.0 4.7 4.7 5.7 4.7 5.0 3.0 5.1 2.9 5.3
##   [55] 3.3 5.2 4.2 4.3 5.6 5.2 5.4 3.9 4.8 3.7 2.6 4.6 4.9 4.8 5.3 4.7 4.3 2.8
##   [73] 4.9 5.7 5.3 4.8 4.5 4.3 3.6 3.3 4.7 4.8 4.0 4.7 4.4 4.7 4.9 4.2 4.3 3.5
##   [91] 5.6 4.1 4.5 2.7 4.2 4.1 2.8 6.2 5.5 3.5 3.6 5.3 6.4 5.2 6.7 3.4 5.6 5.7
##  [109] 5.5 2.9 3.9 5.5 4.8 4.2 5.2 4.5 3.7 4.4 3.9 3.9 4.1 3.5 4.2 3.0 4.3 5.8
##  [127] 5.7 4.9 5.3 5.5 3.1 4.4 5.6 5.1 4.7 4.3 4.0 4.0 4.3 4.1 3.8 3.7 5.3 3.3
##  [145] 4.9 5.0 4.1 4.0 6.3 5.1 4.1 3.5 7.0 4.5 3.7 3.1 3.6 4.8 5.8 5.4 5.6 4.0
##  [163] 5.3 5.1 4.9 6.7 6.5 5.7 4.0 5.2 5.7 4.2 4.1 4.5 3.2 3.9 4.0 4.9 4.8 4.3
##  [181] 5.2 4.9 5.0 5.5 4.1 4.4 4.0 4.6 4.9 4.2 3.7 6.8 4.5 2.8 3.3 5.2 4.5 3.3
##  [199] 5.0 4.8 5.8 4.4 4.5 3.6 3.5 5.2 3.3 5.0 6.0 4.2 3.7 4.4 6.4 4.9 3.7 4.0
##  [217] 4.9 3.4 5.8 4.2 3.8 4.2 4.2 5.3 5.1 3.7 4.7 4.4 5.6 5.3 3.9 2.1 4.6 5.7
##  [235] 5.0 3.9 5.3 3.0 4.7 5.4 4.6 5.0 4.0 4.5 5.2 5.4 5.7 5.4 4.4 4.7 6.4 5.0
##  [253] 3.2 3.0 3.2 4.3 2.7 5.1 2.6 3.3 4.3 2.9 3.5 4.2 4.3 4.7 6.1 4.9 5.2 3.6
##  [271] 5.5 3.2 4.1 6.2 5.0 4.0 5.1 4.3 4.6 3.5 5.2 3.7 6.1 6.2 4.6 5.4 4.9 4.1
##  [289] 6.2 6.0 4.0 4.2 5.3 4.0 6.5 4.0 4.2 2.7 2.6 5.2 5.4 4.7 4.6 4.9 5.1 3.9
##  [307] 5.2 5.8 4.7 4.9 4.3 5.5 3.1 3.9 6.4 5.3 4.0 5.6 5.9 4.8 5.0 5.2 5.5 3.8
##  [325] 2.4 5.9 4.3 4.3 3.3 5.2 3.2 5.1 4.6 4.4 4.0 4.0 5.4 3.3 4.3 3.3 3.8 5.2
##  [343] 5.9 2.8 5.1 4.5 4.0 4.0 6.1 3.4 4.8 2.8 4.9 4.1 5.9 2.5 2.9 4.5 4.5 4.3
##  [361] 3.8 5.2 2.7 5.9 4.4 2.7 3.9 4.3 4.9 4.5 4.9 5.0 4.6 5.4 6.1 5.1 5.7 4.7
##  [379] 4.4 5.4 4.6 4.7 5.1 3.8 4.2 4.3 4.1 4.5 5.6 3.5 4.6 4.7 3.8 3.8 3.3 4.5
##  [397] 5.2 4.4 4.8 2.3 4.3 5.4 5.0 3.1 5.0 3.8 3.7 4.2 4.8 4.8 3.2 5.6 6.4 5.4
##  [415] 3.6 4.6 5.3 3.1 5.4 3.5 4.8 5.4 6.0 4.5 4.3 5.1 5.2 5.4 4.2 3.2 3.8 3.7
##  [433] 4.8 5.9 4.3 4.3 4.3 3.6 4.8 6.1 4.1 4.0 4.6 4.3 4.1 4.0 4.6 5.0 2.9 4.1
##  [451] 5.0 6.5 5.8 5.0 4.2 4.6 4.3 6.3 4.9 4.6 4.7 4.8 5.9 6.4 4.6 5.6 5.2 5.4
##  [469] 3.6 4.7 6.9 5.8 6.1 4.1 4.2 4.5 4.7 4.2 4.3 3.1 3.7 4.9 2.3 4.8 3.9 4.3
##  [487] 5.4 4.2 5.5 5.1 3.3 3.8 5.0 4.5 5.3 4.5 3.5 4.2 4.7 4.7 3.9 4.9 5.5 4.4
##  [505] 5.2 5.3 4.3 4.2 5.2 4.8 5.8 5.4 4.3 4.1 4.4 5.0 4.9 5.2 2.6 5.4 5.1 5.7
##  [523] 5.0 4.7 5.3 2.7 3.4 4.8 5.5 3.9 4.3 4.5 7.6 4.5 5.6 3.6 4.9 5.2 3.5 5.2
##  [541] 4.5 3.1 5.1 4.5 5.0 3.8 5.0 5.1 3.0 3.9 3.6 4.9 5.8 4.1 4.7 3.9 5.7 3.1
##  [559] 2.8 3.6 3.1 4.1 4.0 5.1 5.6 4.1 4.1 4.6 4.2 3.4 3.1 4.9 4.6 6.4 4.3 3.6
##  [577] 5.2 4.9 2.1 4.5 3.5 3.4 4.3 4.0 4.6 6.2 4.3 5.1 6.2 4.9 5.0 3.6 4.4 2.9
##  [595] 5.2 5.5 4.2 5.4 4.9 4.4 4.6 6.2 4.4 4.2 5.0 2.5 4.3 3.5 4.3 6.4 6.1 6.9
##  [613] 5.0 4.0 5.9 3.6 4.0 5.5 2.8 4.6 3.8 6.5 4.4 5.4 3.9 3.8 4.9 3.9 4.8 4.6
##  [631] 4.7 6.0 6.4 3.7 3.1 4.8 3.9 4.0 5.8 3.0 4.4 5.3 4.8 4.4 3.9 4.2 5.1 5.4
##  [649] 5.1 4.9 3.8 4.8 4.1 4.4 3.4 5.2 4.6 5.3 3.3 4.2 5.5 3.6 3.9 4.5 4.3 5.8
##  [667] 5.0 4.0 4.7 4.5 4.2 4.2 3.0 5.3 5.8 4.2 5.3 4.7 6.8 5.4 3.5 5.8 4.9 5.4
##  [685] 5.2 5.0 4.8 3.9 4.4 4.7 2.6 3.8 4.4 3.9 4.7 3.7 4.5 2.6 3.3 5.3 3.8 4.4
##  [703] 4.9 4.7 4.7 5.4 5.8 3.6 5.5 3.4 6.0 4.0 5.2 3.7 5.2 3.2 3.3 4.1 5.2 5.3
##  [721] 5.2 4.3 4.4 5.5 2.8 4.2 3.7 3.9 6.2 4.7 5.2 4.9 3.9 4.4 5.1 5.1 3.3 5.3
##  [739] 4.6 2.8 4.2 4.8 4.3 5.9 4.8 5.7 5.8 4.4 4.8 2.8 4.6 4.8 3.8 4.4 5.2 3.7
##  [757] 4.7 3.9 3.3 3.9 5.4 3.8 4.7 4.0 4.4 5.8 4.6 6.4 3.0 3.0 3.6 4.5 3.6 5.1
##  [775] 4.2 6.4 4.2 5.1 4.6 4.3 3.9 5.5 5.1 4.9 4.8 5.2 5.8 3.7 4.7 3.9 4.4 3.7
##  [793] 5.7 3.8 5.1 3.9 5.1 5.1 5.8 4.5 3.0 3.7 6.1 4.9 4.1 4.7 5.0 3.4 5.2 3.8
##  [811] 3.8 5.3 2.8 5.9 4.9 5.7 3.4 4.7 3.9 4.0 5.3 4.7 4.7 3.5 4.9 5.1 5.1 5.4
##  [829] 4.9 3.2 4.2 5.2 5.2 3.8 4.3 4.5 5.7 4.9 5.3 5.4 4.7 3.5 4.6 3.9 5.1 5.4
##  [847] 4.1 5.0 3.0 4.2 5.2 6.1 2.6 4.3 4.7 5.2 5.6 4.0 3.5 2.4 4.0 5.4 5.4 5.8
##  [865] 4.6 6.1 4.5 6.4 4.9 5.5 4.2 4.3 4.7 5.4 3.7 3.4 5.0 4.9 3.1 5.4 4.5 4.1
##  [883] 4.9 5.0 5.1 5.0 3.7 3.9 5.0 3.0 3.6 3.3 5.3 5.7 4.2 4.9 4.2 3.6 3.7 4.0
##  [901] 4.2 4.6 4.5 3.9 5.4 4.3 5.5 2.8 5.9 3.6 3.9 4.2 3.0 3.7 2.9 2.9 5.9 3.3
##  [919] 4.6 6.0 4.1 5.8 4.1 5.0 4.8 3.9 5.2 4.7 4.9 4.9 3.0 4.7 6.9 4.1 4.2 5.3
##  [937] 4.5 4.9 5.1 5.9 4.3 4.2 4.0 3.8 4.3 5.0 5.2 3.7 5.0 4.6 4.2 4.3 4.2 4.7
##  [955] 2.9 4.8 6.1 5.2 3.4 3.2 5.0 5.4 5.4 5.0 4.3 6.1 2.6 4.1 4.9 3.9 4.2 5.1
##  [973] 4.8 4.0 4.6 3.5 3.2 3.8 4.8 5.6 3.1 3.4 3.6 3.6 4.2 3.0 3.4 3.5 5.0 5.8
##  [991] 4.3 4.7 4.9 4.6 4.5 3.8 3.7 4.4 4.3 4.1
## 
## $func.thetastar
## [1] 0.0349
## 
## $jack.boot.val
##  [1]  0.5248587571  0.3901639344  0.3172910663  0.2605882353  0.1051575931
##  [6]  0.0005434783 -0.1884507042 -0.2157480315 -0.3727554180 -0.4780780781
## 
## $jack.boot.se
## [1] 0.9698425
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
##    [1] 4.9 4.5 5.2 4.1 6.1 5.0 4.4 4.2 4.7 4.1 5.3 4.0 3.9 4.8 4.1 4.0 4.9 5.5
##   [19] 4.9 3.8 4.5 4.7 3.5 5.2 3.2 4.5 5.3 4.0 4.6 4.3 4.6 4.8 4.0 2.9 6.5 6.4
##   [37] 4.6 4.2 4.5 4.5 6.0 4.1 3.1 5.9 5.0 4.0 5.0 3.8 4.9 4.8 5.7 5.6 3.1 2.4
##   [55] 5.7 3.3 3.6 3.8 4.0 3.1 4.8 5.4 4.6 4.7 4.9 4.1 4.7 3.5 3.5 5.0 5.6 5.6
##   [73] 5.5 3.1 5.4 6.0 4.8 3.6 3.9 4.7 5.6 6.0 3.0 4.3 4.1 3.4 4.3 5.1 3.8 5.1
##   [91] 4.6 3.1 5.0 3.4 4.5 4.8 4.1 5.6 4.6 3.7 5.1 4.7 5.9 2.9 2.9 5.7 3.3 5.7
##  [109] 3.6 4.0 4.9 5.6 4.0 3.5 3.9 6.6 4.4 3.7 5.2 2.4 3.9 4.1 5.6 5.2 4.3 5.4
##  [127] 3.9 4.7 5.1 3.8 4.3 3.9 3.5 2.8 3.6 4.3 3.0 5.2 3.2 2.6 4.3 4.9 4.3 4.1
##  [145] 5.0 2.8 5.2 3.6 4.3 4.2 3.6 3.8 4.3 5.5 5.4 4.6 5.3 5.2 4.4 4.4 4.9 2.7
##  [163] 3.1 4.0 5.4 5.1 4.1 4.5 4.5 5.1 4.4 4.1 4.4 6.3 5.3 3.5 5.0 4.2 3.0 4.5
##  [181] 5.6 5.7 4.1 4.7 5.4 5.6 4.4 5.0 4.3 4.7 4.0 4.6 4.3 4.2 5.3 4.8 4.4 3.0
##  [199] 3.3 5.5 3.2 4.7 4.5 3.2 3.9 3.3 4.1 5.2 3.9 3.2 5.4 5.2 5.5 5.9 4.2 3.5
##  [217] 3.5 4.9 3.7 4.8 3.8 5.5 4.9 3.8 4.2 3.6 3.9 4.9 5.4 5.1 5.4 4.5 3.6 4.1
##  [235] 5.3 3.6 3.6 3.5 3.6 5.5 5.6 2.7 3.8 4.4 4.7 5.0 2.2 3.7 4.7 5.1 4.5 4.4
##  [253] 4.0 3.1 4.9 4.9 3.3 5.3 5.2 5.4 3.6 6.2 5.6 5.2 4.6 5.3 4.5 5.4 5.7 3.1
##  [271] 3.7 6.4 2.6 4.6 3.9 4.2 4.4 5.2 3.1 5.1 5.0 4.1 4.2 3.5 5.3 4.5 5.4 3.9
##  [289] 4.0 3.9 4.6 5.0 4.0 4.1 6.9 5.0 4.3 5.4 5.0 4.8 3.3 4.8 4.9 4.6 4.8 3.9
##  [307] 4.2 3.2 4.5 5.6 4.4 3.0 4.7 6.0 4.3 6.1 2.8 5.0 3.8 5.5 3.2 4.4 3.2 3.6
##  [325] 2.5 4.3 4.4 4.7 4.8 5.0 2.7 3.8 5.4 5.0 4.1 4.6 4.1 4.8 6.0 4.1 6.1 4.5
##  [343] 3.1 3.5 3.3 4.7 5.0 5.5 4.1 4.7 3.7 3.6 4.3 5.7 4.4 2.9 4.0 4.0 4.6 5.3
##  [361] 3.3 2.7 4.2 3.5 6.0 4.7 4.6 3.4 5.0 5.4 5.8 5.3 4.4 3.2 2.9 3.3 5.0 5.0
##  [379] 5.2 3.7 5.0 5.1 4.8 5.4 4.9 3.2 3.2 4.3 5.1 4.8 6.2 5.1 4.0 3.4 4.7 4.0
##  [397] 3.3 4.5 5.4 5.1 5.5 4.5 4.0 4.4 5.1 5.7 3.6 5.3 3.1 5.1 5.3 4.2 5.0 3.9
##  [415] 6.3 4.5 3.4 3.7 3.5 5.6 4.1 5.0 3.3 4.3 6.2 5.0 4.6 4.2 5.4 6.3 5.5 3.9
##  [433] 5.9 4.4 4.8 3.5 5.1 4.1 4.5 4.5 4.7 4.1 7.4 4.2 4.8 7.9 5.8 3.0 5.7 5.0
##  [451] 5.2 3.9 4.4 5.3 4.3 4.9 5.0 4.0 5.0 4.9 2.9 4.4 4.1 5.5 5.4 4.6 4.3 4.4
##  [469] 3.8 4.2 5.2 4.3 5.1 4.8 6.6 3.8 2.8 5.1 5.4 4.3 3.6 5.3 4.9 4.7 4.4 4.1
##  [487] 6.6 5.3 4.0 4.4 4.5 5.3 3.7 6.1 4.8 4.1 4.7 5.6 3.9 4.2 4.9 3.9 3.8 4.7
##  [505] 5.8 7.0 4.6 5.1 4.6 4.1 4.4 3.1 5.5 4.4 5.7 1.6 3.0 6.1 3.5 4.6 3.6 3.7
##  [523] 6.4 6.0 3.1 7.1 4.6 6.0 4.8 4.0 5.1 4.0 4.1 4.7 3.2 4.4 3.5 3.8 5.0 4.7
##  [541] 5.6 4.2 4.9 5.5 4.9 4.5 2.9 5.6 5.3 4.1 4.1 4.1 3.0 4.8 5.0 3.4 4.8 4.4
##  [559] 6.6 3.9 3.3 3.4 4.0 3.1 4.6 5.2 4.5 4.7 3.2 2.7 4.9 3.5 5.1 4.5 4.6 3.9
##  [577] 4.5 4.9 2.9 4.6 4.1 5.1 5.7 3.4 3.9 3.9 4.7 3.5 3.4 5.0 2.7 3.1 4.4 4.1
##  [595] 3.0 3.1 4.1 5.6 5.6 4.1 3.2 4.9 3.0 4.0 4.1 4.0 4.0 4.9 2.5 3.8 6.4 4.5
##  [613] 3.5 4.4 4.7 5.8 6.5 5.2 5.6 3.6 3.6 3.6 3.9 4.1 3.6 3.4 3.8 4.8 4.1 4.9
##  [631] 4.8 4.7 5.8 5.2 5.3 3.4 5.8 5.2 2.7 3.4 5.3 3.0 6.1 3.9 3.7 4.1 4.8 3.7
##  [649] 5.0 5.0 5.7 3.9 6.1 4.6 4.2 4.6 4.0 3.3 4.1 4.9 4.4 5.8 6.5 5.8 4.9 3.4
##  [667] 3.6 4.8 5.4 4.6 4.0 4.2 4.6 5.8 4.0 5.0 4.6 2.5 3.4 4.3 5.3 3.6 6.7 4.5
##  [685] 3.7 5.4 5.1 5.0 3.8 3.9 5.3 6.3 3.7 4.0 4.3 3.9 4.0 6.7 6.2 5.7 3.1 5.4
##  [703] 5.0 4.2 4.7 3.6 4.4 5.6 4.9 3.8 6.1 3.3 5.2 5.0 3.0 4.9 3.9 4.4 3.8 5.6
##  [721] 4.9 5.4 4.3 4.8 5.5 4.1 4.7 3.8 5.7 4.8 2.0 5.8 4.3 4.6 3.9 6.8 5.3 3.6
##  [739] 4.0 4.3 4.4 5.1 4.5 5.0 4.2 5.0 4.9 5.5 4.3 6.0 4.2 4.1 3.5 6.6 4.7 5.1
##  [757] 2.6 4.7 4.2 4.1 4.2 3.8 3.5 4.5 4.4 4.2 3.6 2.6 3.4 3.9 3.6 4.8 3.1 3.5
##  [775] 4.0 4.8 5.4 5.0 4.9 3.6 3.2 5.0 4.6 5.4 4.8 3.2 3.5 4.2 2.7 4.5 3.6 4.0
##  [793] 3.8 4.5 5.2 5.4 4.3 4.0 5.0 4.7 5.4 5.9 5.1 3.1 4.5 3.5 5.1 5.9 4.5 5.8
##  [811] 2.2 3.2 4.3 3.7 5.7 4.9 3.2 6.2 5.6 3.9 3.9 4.7 3.7 3.7 4.3 4.6 5.0 5.2
##  [829] 3.4 3.9 3.8 4.1 4.7 5.5 4.5 2.6 4.3 5.6 4.3 5.7 5.1 3.9 4.2 4.0 4.0 6.4
##  [847] 4.5 4.5 3.9 5.3 6.5 4.7 4.8 4.1 3.8 2.9 4.3 2.8 3.3 4.6 3.6 4.7 5.8 3.1
##  [865] 3.7 5.7 3.4 4.2 4.5 5.8 4.7 5.6 4.7 4.2 4.1 4.8 5.7 2.9 3.0 3.8 4.0 3.8
##  [883] 5.4 5.1 5.2 4.0 4.1 5.7 3.2 3.4 3.5 5.0 4.3 5.9 6.6 4.6 4.0 5.7 3.7 7.4
##  [901] 3.6 3.7 3.2 3.3 3.6 5.7 5.7 4.2 3.8 2.7 4.7 3.3 4.2 3.8 5.6 5.0 4.4 2.6
##  [919] 5.7 5.0 3.2 4.9 4.2 5.1 4.4 6.0 5.1 4.5 3.8 5.8 4.5 4.8 6.6 5.1 4.9 5.8
##  [937] 4.0 5.1 4.2 5.4 3.8 3.3 3.3 2.6 3.2 3.6 4.3 4.0 6.7 3.1 4.7 4.6 3.3 4.5
##  [955] 4.2 5.4 5.1 7.4 4.0 4.4 6.2 3.9 5.4 5.2 5.4 6.0 4.3 4.8 3.4 5.3 5.1 6.7
##  [973] 3.3 4.4 4.8 5.5 4.4 4.5 4.9 5.8 4.8 3.6 5.7 4.6 3.6 4.4 3.9 5.0 5.5 4.3
##  [991] 5.2 5.4 6.7 4.7 7.0 5.2 4.5 5.7 4.8 4.0
## 
## $func.thetastar
## 72% 
##   5 
## 
## $jack.boot.val
##  [1] 5.50 5.50 5.30 5.30 5.20 5.10 4.80 4.76 4.70 4.50
## 
## $jack.boot.se
## [1] 1.005801
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
## [1] 0.9775653
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
##   4.873203   9.738421 
##  (2.108997) (4.439417)
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
## [1] -0.07801040  1.04399498  1.57928064  0.09091382  0.06234088  1.18508848
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
##    [1]  1.302418655  0.741228431  1.858401471  1.472990837 -0.416888336
##    [6]  1.058900543  0.232316355  1.049893454  0.770808867  0.620203481
##   [11]  0.155443108  0.577063143  0.948097990  0.474285238  1.365749020
##   [16]  0.896934465  0.630702056  0.481257145  1.575989384  1.509286342
##   [21]  0.381895669  0.629265809  0.384299177  0.337693667  0.513976139
##   [26]  0.977565326 -0.122579950  0.510201467  0.771116536  0.349639206
##   [31]  1.635096856  0.853320987  0.501881999 -0.676611764  0.797731066
##   [36] -0.162539886  1.513942377  0.731523921  1.584749737  0.881460659
##   [41]  0.480817587  1.290207503  1.297523204  0.501332888  1.673348852
##   [46] -0.081452276  0.789588213  1.852938967  0.106876608  1.091424538
##   [51]  0.122635624  0.369902399  1.356979194  0.180836157  0.929686563
##   [56]  0.978771723  0.584002367 -0.406312353  0.442740819 -0.471010753
##   [61]  0.879050883  0.811555961  0.230834459  0.812912674  1.301350425
##   [66]  0.647229013 -0.463708073  0.297832437  0.955344597  0.528105284
##   [71]  0.623427412  0.715810412  0.146309269  0.938436266  1.162801017
##   [76]  0.638317919  0.003777479  0.752416335  1.811800825  0.536910614
##   [81]  0.323914006 -0.207523231  0.460908462  1.430295737  0.774953384
##   [86]  0.667655384  0.611848517 -0.224017249  2.058020915  1.874014786
##   [91]  1.717901416  0.793014538  0.825742757  0.473429424  0.018068892
##   [96] -0.041145150  0.624582337  0.807651502  1.633356148  0.555034523
##  [101]  0.784337572  0.282498848  0.223794639  1.100754594  0.832140110
##  [106]  0.272014948  1.057359178  1.283476692  1.445318146  0.209894981
##  [111]  0.977129183  0.730750185  0.734472063  0.780458870  1.545500045
##  [116]  0.003845626  1.473568337  1.158489149  1.112723278  1.007561089
##  [121]  1.742432537  0.398677062  0.499442233  1.291786175  0.316778375
##  [126]  0.360028688  0.417294743  0.938443117  0.859462453  0.452273552
##  [131]  0.880686010 -0.259638740  0.517220228  0.077426810  1.357849012
##  [136]  0.773642972  1.283185109  0.558687890  0.018863848  0.146229926
##  [141]  0.665271875  0.699676698  1.146359962  0.613453260  0.211119179
##  [146]  0.744712334  1.052970491  0.533770349  0.735247909  0.700528245
##  [151]  1.106493635  0.384365263  1.237651851  0.104657057  0.661554835
##  [156]  0.093245192  1.756744328 -0.353195278  0.542919107  0.970062309
##  [161]  0.646285622  1.688846617  0.022973793  0.376154235  1.219631296
##  [166]  0.879050883  0.307013244  0.987945424  1.441006511 -0.042501100
##  [171] -0.058706189  1.287038283  1.140559387  1.522658564  1.652220374
##  [176]  0.409322046  1.045666054  1.396376499  1.232800341  1.072298050
##  [181]  0.867145358  1.073148599  0.962006849  0.904139355 -0.835526011
##  [186]  0.821642232 -0.059610907  0.596857642  0.176353155  1.301095724
##  [191] -0.001698954  1.141738014  0.822384256  0.960183702  0.954202636
##  [196]  0.251460010 -0.503787530  1.081801971  1.798814156  0.071147003
##  [201]  0.232174829  0.349254081 -0.053067832 -0.158823546  0.174796607
##  [206]  1.011866291  1.031720795  0.077926118  0.622674453  0.411976684
##  [211]  0.489054773  1.357168851  0.557620921  0.875438918  0.483702786
##  [216] -0.496635002  0.797403571  0.900141950  0.730157332  0.329854528
##  [221]  0.207561771  0.219601063  0.883246802  0.907595332  1.725874577
##  [226]  0.517901318  1.019021576  1.468779669  0.908934352  0.059318936
##  [231]  0.670729435  0.690779995  0.785318660 -0.107631761  1.245190002
##  [236]  0.871816229  0.729496771  1.091538776  0.335674727  0.666213069
##  [241]  0.741202435  0.286208043  0.823550574  0.547907284  0.589442685
##  [246]  1.488103714  1.426476086  0.665131537  1.229062698  0.384909788
##  [251]  0.192257493  0.181569459  0.803620572  1.220239841  1.824842286
##  [256] -0.189220268  1.017025561  0.292104596  0.686873719  1.002917142
##  [261]  0.003707311  0.904691732  0.770614966  0.659732675  1.212131448
##  [266]  0.662330010  0.957939874  0.605189878  1.064739653  0.747285394
##  [271]  0.846574162  1.784058673  0.538383043  0.011944192  0.068849915
##  [276]  0.170835284  0.195909875  1.250254534  0.993042517  1.404858763
##  [281] -0.194759714  1.257468361  0.918875123  0.346248972  1.387706752
##  [286]  0.865995339  1.906442234  0.142279288  0.590153432 -0.336664546
##  [291]  0.228483028  0.336517926  0.341150527  0.875534806  0.614700576
##  [296]  0.574507800  1.552366478  0.217792801 -0.238312010  0.688695399
##  [301]  1.182466296  1.336808356  1.057621268  0.731544707 -0.141848089
##  [306]  0.246726678  1.044135779  1.036697885  0.669243121  0.620205337
##  [311]  0.140909927 -0.106079747  1.104310895  2.141490255  0.278422533
##  [316]  0.787848390  1.370428438  0.293531535 -0.008305345  1.444188828
##  [321]  0.919685423  1.457439155  1.756220697  0.994484629  0.729359995
##  [326]  1.373749911  0.507139466  0.585186451  0.387845180  0.097685321
##  [331]  0.129496622 -0.256052844  1.368765441  1.303990248  1.546143714
##  [336]  1.520563326  0.356126001  0.056654833  1.742432537  1.214292908
##  [341]  1.205615595  0.275038634  1.041650249  0.167112044  0.308119679
##  [346]  0.611903426  0.627236414  0.912683705  0.378552263  1.852582604
##  [351]  1.371710044  1.268535463  0.439543426  1.229435483  1.747380278
##  [356]  0.219050280 -0.608936372  0.499556039  0.747466994  0.936763279
##  [361]  0.384948470  1.364646458  0.368330566  0.684014220  0.704051610
##  [366]  0.975931302  0.794701560  0.321750969  0.635618570  1.082820403
##  [371]  0.332929997 -0.258451449  0.333954678  1.139881469  1.077041368
##  [376]  0.064301542  0.251187988  0.281787650  0.610947292  0.706502223
##  [381]  1.026769014  1.582343693  0.376424364  0.267482910  0.692499564
##  [386]  0.065713984  1.133386280  0.568696917  0.674801258  0.610265256
##  [391]  1.459971407  1.654751843  1.198638586  0.925632572  1.100537265
##  [396]  1.132401124  1.238498647  0.705351761  0.821870755  1.132523256
##  [401]  0.960967671  0.478724787  1.032918315 -0.396338277  0.313958456
##  [406]  1.008297488  0.845456769  0.834691832  1.822016384  0.392918197
##  [411]  0.532132358  1.299058263  1.271809140  0.538916427  1.576561452
##  [416]  1.283485706  0.732148760  0.579798755  0.805755485  0.422288358
##  [421]  0.479564533  0.405706694  0.443025813  1.102794756  0.534561952
##  [426]  0.707179666  0.482259541  1.771573592  0.953315875 -0.434304395
##  [431]  0.769240825  1.425544021  0.561075328  0.754489294  0.299192139
##  [436]  1.202704472  0.077652235  0.854550778  1.073893433  1.042724399
##  [441]  1.084906459  0.468632885  2.298387155  0.646104779  1.081417579
##  [446] -0.611683010  0.237502161  0.512444987  0.298999466  0.321936636
##  [451]  0.829873333  0.919382171  0.882746810  0.074403591  0.830830444
##  [456]  0.670662315  0.190249981  0.369984180  0.195697487  1.603970941
##  [461]  0.830710931  0.618285327  0.046282130  1.652337914  1.979268292
##  [466]  1.115082759  1.026534760  0.804224130  1.380487225  0.722694706
##  [471]  1.012868193 -0.915412832  0.022615643  0.433955380  0.634873857
##  [476]  1.100598886  1.477295477  0.959085306  0.330376906 -0.286570815
##  [481]  0.681833984  0.719736726  0.620874690  0.793663880 -0.118566704
##  [486]  0.021242797  0.650051908  1.119447634  0.710609495  0.590014427
##  [491]  0.873191021  0.768540223 -0.557906502  0.897635419  0.168707221
##  [496]  0.119271416  1.543204939  1.443899852  0.954343163  0.686796434
##  [501]  0.407954444  1.495521899  1.601679627  0.397419574  1.476819781
##  [506]  1.404858763  1.311164442  0.962659100  0.433099626  1.043284962
##  [511]  0.049755756  0.950117842  0.062983663  0.064994470  0.604786172
##  [516]  0.500441166  0.546031177  0.543651089  1.127546561  1.302418655
##  [521] -0.619584522  0.585772135  0.049896544  0.169219899  0.978872737
##  [526]  0.567120194  0.451167826  1.172749526 -0.228748237  0.854804403
##  [531]  2.113457013  0.034178331  0.357536242  0.194931698  1.020633755
##  [536]  0.344520472  0.713488710  0.478336174  0.380037234  0.510918643
##  [541]  0.675112302  1.781312089 -0.066265019  0.227757563  0.466338856
##  [546]  0.720660603  1.218714155  0.557077245  0.832261905  0.106556581
##  [551]  0.709620848  0.742723044  0.633033967  0.344131716 -0.071323454
##  [556] -0.714079022  1.498871161  0.604350487  0.600047810  0.923622192
##  [561]  1.146184428  0.754691313  1.391020066  0.673964118 -0.854201646
##  [566]  0.465314284  0.545468969  0.008308475  0.821708899  0.283556931
##  [571]  0.218862908  1.290354170  0.967912038  0.676759460 -0.029150206
##  [576]  1.211732121  0.930391984  1.746790509  1.093301814  1.620914332
##  [581]  0.746902612 -0.254597371  0.167769008  1.184100597  0.764944924
##  [586]  1.694505367  0.740341215  0.002722238  0.093670867  1.415358191
##  [591]  0.642242369  0.335772879  0.268807733  0.841178979  0.248788080
##  [596]  0.877685783  0.280890055  0.906313720  1.253753817  1.206284653
##  [601]  1.244561077  0.952627281  0.803670409  1.086702575  1.819541697
##  [606]  1.485922553  0.164005977  0.163701147  0.824606890  1.215917629
##  [611]  1.341005532  0.918381847  0.220195482  0.709482319  0.621900658
##  [616]  0.417476154  0.450784391  0.737829980  0.112084826  0.503103341
##  [621]  0.631332129  0.257927554  1.864872310  0.894965139  1.468330857
##  [626]  0.741884071  1.532502314  1.091435641  1.005308218  0.934283374
##  [631]  0.121650909  1.203193660  1.665619850  0.439049997  0.309825238
##  [636]  1.295036379  1.332164137  0.901548695  1.125216044  0.214768224
##  [641]  0.944516517  0.670839985  0.836449349  1.272033732  0.953189435
##  [646]  0.760393365  0.436116102  1.507442582  0.912527365  0.593192540
##  [651]  1.093283013 -1.374824760  1.097135699  0.336093320 -0.008141245
##  [656]  1.146184428 -0.181636659  1.024682988  0.660897735  0.526158064
##  [661]  1.265178775  0.287513388  0.436097050  1.235570762  0.981810600
##  [666]  0.047020327  0.036842682  0.632454994  0.226905324  0.932953508
##  [671]  1.098672664  1.401858649  0.810661025  0.239175449  0.736447524
##  [676]  0.656271798  0.569263465  1.102578025  0.891202143  0.911449334
##  [681] -0.216436943  1.083508555  1.365236521  0.428533817  0.296807405
##  [686]  1.134694625  0.286224724  1.576868886  1.328762167  1.499314301
##  [691]  0.660507552  0.740136845  1.709276881 -0.161743256  0.422336927
##  [696]  0.344288730  0.384948470  0.744181830  1.000786989  0.515346558
##  [701]  0.181715843  0.798946253  0.048702880  0.792819343  1.308230618
##  [706]  1.096893232  1.836869453  0.480079399  0.385623698 -0.059594994
##  [711]  1.329896255  1.535356136  0.771290768  1.025999851  0.308895223
##  [716]  0.430340400  0.529050643  0.149331621  1.165500871  0.321750969
##  [721]  0.269740654  1.349813019  0.557107540  1.113637535  0.123078230
##  [726]  1.738812263  1.593802545  1.147148141  0.280971855  0.884504686
##  [731]  2.084710914  0.602835702  1.089220472  0.343143742  0.975143861
##  [736]  0.425059888  1.759302469  0.704933098  0.550851753  0.843591382
##  [741] -0.210073255  1.805191102  0.938897368  1.576435093  0.224583792
##  [746] -0.580663063  0.912022087  0.337249638  0.486142619  0.712970567
##  [751]  1.250680718  1.091207423  1.211144230  0.120462629  1.091883747
##  [756]  1.094856332  0.411202540  0.524853320  0.528310277  0.771168595
##  [761]  1.164016291  0.459113933  1.019540915  0.350930847  0.321274917
##  [766]  0.754627087  0.985531991 -0.555640245  1.046663267  0.425954019
##  [771] -0.348074926  0.733627019  0.559085555  0.013537891  0.648979590
##  [776]  0.901730434  0.382482652  0.460785340  0.154650231  1.148010393
##  [781]  0.922217597  1.173554072  0.781192193  1.781537472  0.750468851
##  [786]  2.066628304  0.184056853  0.241298083  1.464230981  1.018678144
##  [791]  1.634463767  0.397611869  0.724816041  1.007432022  0.725368682
##  [796]  0.767878201  0.235005924  0.963991218  0.187743653  1.025213165
##  [801]  0.586935884  0.594709280  1.811398551  0.552445157  0.403548422
##  [806]  0.690292518 -0.028521185  1.398657240  0.634931716  0.743913657
##  [811]  0.587432860  0.914010086  1.897793660  1.125409795  0.450313954
##  [816]  0.510758721  1.156926218  0.539237427  0.314034403  0.270144054
##  [821]  0.231141590  1.624518843  0.921172631  0.501558837  0.141099553
##  [826]  0.716505346  1.056833905  1.233385342  1.030461273 -0.522364535
##  [831]  0.628105803  0.471447693  0.659415898  0.596352171  0.595407322
##  [836]  1.520232551  1.540193067  0.601058515  1.125288311  0.969759317
##  [841]  0.731161588  1.405085525  0.655326471  0.402327398  0.678544288
##  [846]  0.707178175  1.202597240 -0.900918250  1.080826520  0.652251364
##  [851]  1.092104258  0.146931830  0.515231428  0.938172316  0.951698908
##  [856]  0.405863037  0.384365263  0.360570278  0.715704159  1.148403277
##  [861]  0.836468637  0.269635977  0.767141928  1.160608681  1.379248584
##  [866]  0.437113428 -0.135226973  0.493512699  1.443899852  0.786218535
##  [871]  0.742685521  0.661554835  0.788730962  0.870779369 -0.693407130
##  [876]  0.538851556  0.335573718  1.996923370  0.425193284  0.945337495
##  [881]  0.770805792  0.716234101  1.050437015  0.590019347  0.435312059
##  [886]  1.024792202  1.037024231  1.471451015  0.720908643  1.571926946
##  [891]  0.138254523  1.154479849  0.342121153  0.131270949  1.714180365
##  [896]  0.596729938  0.463058150  1.239236229  0.830710931  1.578767648
##  [901]  1.810883869  0.673254822  0.443907560  1.551957591  0.577917497
##  [906]  0.222224289  1.301809465 -0.203513912  0.670486887  0.781061549
##  [911]  1.494079018  1.019085558  0.118084020  1.100085043  0.193821881
##  [916]  0.723188082  0.557174151  0.927451556  0.515231428 -0.804579041
##  [921]  0.417197574  0.283337472  0.828506571  1.038665825  0.601745728
##  [926]  0.577860080  0.671085946  0.746447742  0.256397818  1.210870079
##  [931]  0.585005532  1.308515540  0.502778881  1.500077820  0.857533481
##  [936]  1.178028886 -0.067737618  0.564606526  0.650355032  1.020499696
##  [941]  0.872235002 -0.109173879  0.643675087  0.297594961  0.504632734
##  [946]  0.763527222  0.902042659  0.990678019  0.654051133  0.499476788
##  [951] -0.158458049  0.228120410  0.776766023  0.610155433  0.277863927
##  [956]  0.795619315  1.107975923 -0.197438008  0.737754802  0.754523220
##  [961]  0.335267040  2.033010145  1.504737563  0.133581446  0.698922413
##  [966]  0.621920797  1.165971518  0.678017527  0.574676118  0.306076782
##  [971]  0.949913912  1.750840226  1.028182354  0.777673310  1.154081722
##  [976]  0.059061610  0.579779791  1.722421359  0.782103876  0.833784377
##  [981]  0.780261361  0.903524657  1.230173437  1.047876230  1.529746676
##  [986]  0.705204607  1.076971520  0.741054039  0.177758905  0.926700754
##  [991]  0.041773509  0.609956762  0.648390897 -0.304056378  0.536713019
##  [996]  0.229607576  1.096893232  0.905633100  0.091521941  1.686953230
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
##   0.50040829   0.23940951 
##  (0.07570793) (0.05352991)
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
## [1] -0.71657487  0.05984161  0.22919204  0.88338007  0.37881379  0.60362130
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
## [1] -0.015
```

```r
sd(results2$thetastar)  #the standard deviation of the theta stars is the SE of the statistic (in this case, the mean)
```

```
## [1] 0.8874698
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
## t1*      4.5 0.01521522   0.8993589
```

One of the main advantages to the 'boot' package over the 'bootstrap' package is the nicer formatting of the output.

Going back to our original code, lets see how we could reproduce all of these numbers:


```r
table(sample(x,size=length(x),replace=T))
```

```
## 
## 1 2 3 4 6 7 8 
## 2 2 1 2 1 1 1
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
## [1] 9e-04
```

```r
se.boot
```

```
## [1] 0.9204327
```

Why do our numbers not agree exactly with those of the boot package? This is because our estimates of bias and standard error are just estimates, and they carry with them their own uncertainties. That is one of the reasons we might bother doing jackknife-after-bootstrap.

The 'boot' package has a LOT of functionality. If we have time, we will come back to some of these more complex functions later in the semester as we cover topics like regression and glm.

