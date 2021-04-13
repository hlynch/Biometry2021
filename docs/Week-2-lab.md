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
## 0 1 2 5 6 7 8 
## 1 1 1 4 1 1 1
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
## [1] 0.0057
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
## [1] 2.662944
```

```r
UL.boot
```

```
## [1] 6.348456
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
##    [1] 3.8 5.7 4.6 5.9 4.0 4.2 3.6 2.6 4.1 5.5 4.9 3.7 4.1 5.5 3.5 5.5 6.0 5.3
##   [19] 4.3 4.5 4.6 4.7 3.7 4.4 4.8 5.0 2.9 3.1 4.4 5.7 6.7 5.3 3.9 4.3 3.8 5.9
##   [37] 5.6 6.5 4.3 4.6 4.0 4.5 4.2 5.0 6.3 4.0 4.1 5.4 5.5 3.5 3.4 3.6 4.8 4.6
##   [55] 4.6 3.8 4.1 4.7 4.3 4.3 3.3 5.4 5.9 3.3 4.3 3.9 5.1 4.0 5.3 5.3 4.1 3.7
##   [73] 5.9 3.3 4.6 5.2 4.6 3.2 4.0 4.0 4.8 5.1 5.5 5.1 5.8 4.8 2.6 4.4 3.9 5.4
##   [91] 1.8 4.8 6.0 4.4 3.5 5.7 5.5 4.2 3.7 4.9 6.1 4.0 3.2 6.4 3.5 4.8 5.2 5.2
##  [109] 4.9 4.3 4.7 4.3 4.0 4.2 3.8 4.3 4.5 4.9 3.7 5.3 5.7 4.2 3.8 5.4 5.0 4.7
##  [127] 3.7 3.9 3.9 6.0 4.7 4.4 3.7 4.5 4.6 2.9 2.9 4.5 6.1 3.8 4.5 3.2 4.0 3.6
##  [145] 5.3 5.4 5.2 4.7 5.9 5.4 4.4 4.5 4.9 4.9 4.1 4.4 5.1 3.7 3.7 5.3 5.3 5.5
##  [163] 5.0 4.9 3.9 4.0 5.1 4.4 4.7 4.8 4.7 4.4 5.3 3.9 4.6 4.2 4.1 6.5 3.1 5.6
##  [181] 3.4 5.5 5.8 4.8 3.6 4.0 3.6 5.2 6.1 5.0 3.9 3.7 4.9 3.9 5.6 5.1 3.7 5.0
##  [199] 4.3 5.6 4.0 4.3 6.0 4.8 4.8 2.3 4.4 3.4 3.2 4.9 5.8 4.2 4.2 5.0 5.0 5.6
##  [217] 4.5 4.1 4.5 4.0 3.1 5.1 3.8 3.9 5.3 4.1 4.5 4.5 5.2 3.4 4.8 4.9 4.1 4.2
##  [235] 3.5 5.4 4.1 5.1 5.0 4.2 4.6 4.8 5.4 5.4 4.8 3.6 4.2 3.6 5.3 4.9 4.7 4.2
##  [253] 5.0 4.9 3.2 5.0 3.2 5.5 3.6 4.3 4.7 4.8 3.5 3.5 3.5 3.6 5.9 5.5 3.6 3.7
##  [271] 4.1 3.2 4.9 4.6 4.8 5.3 4.4 5.4 4.3 4.6 3.9 2.3 3.6 4.9 4.6 5.6 5.9 4.5
##  [289] 4.5 3.8 5.4 4.6 4.1 5.8 5.2 3.7 4.0 1.9 2.7 4.2 4.5 3.8 5.9 4.1 2.9 4.8
##  [307] 3.9 3.7 3.5 4.4 3.6 3.3 4.2 3.1 2.8 4.2 4.8 5.3 5.2 4.2 4.7 5.2 5.2 5.3
##  [325] 5.0 4.4 5.8 5.7 5.8 4.4 2.7 6.6 5.0 6.1 5.1 3.4 3.3 4.4 4.4 2.0 4.2 4.7
##  [343] 5.8 4.1 5.4 4.9 4.0 6.1 4.5 4.2 4.8 4.9 6.2 4.0 4.1 4.2 5.5 4.2 5.4 4.9
##  [361] 4.6 5.5 4.2 4.5 3.3 4.4 5.1 4.5 3.5 3.7 4.7 4.4 4.5 4.7 5.5 4.6 4.9 4.5
##  [379] 6.7 3.4 4.0 3.9 4.7 5.2 3.1 4.6 5.0 4.3 3.0 5.1 4.3 3.9 5.0 4.6 4.6 5.0
##  [397] 4.5 5.3 5.6 6.3 5.8 3.8 4.8 4.9 5.4 4.0 4.3 6.1 4.3 3.7 5.0 4.9 3.8 4.8
##  [415] 5.1 4.4 3.5 4.1 5.1 4.6 3.9 5.9 4.8 5.9 4.9 3.1 5.6 3.6 4.2 5.1 5.5 4.9
##  [433] 4.9 5.9 4.0 4.9 5.3 4.6 4.1 4.8 6.4 4.9 4.9 3.1 3.8 5.1 5.3 5.4 3.7 2.0
##  [451] 5.0 6.5 2.8 5.1 4.0 4.4 4.7 4.5 5.5 4.0 5.5 4.0 5.2 5.3 4.7 5.3 2.8 5.1
##  [469] 6.4 5.0 4.6 3.7 5.1 4.1 4.4 5.4 3.9 2.6 5.0 4.7 5.0 5.0 3.9 2.7 5.2 4.2
##  [487] 5.5 4.7 3.9 5.2 4.4 5.3 2.4 4.1 4.0 5.0 4.1 4.1 4.1 4.3 3.7 5.4 4.2 5.4
##  [505] 4.5 4.2 6.3 5.1 4.1 5.1 2.7 4.2 4.2 4.2 4.8 4.7 4.1 4.8 4.5 5.2 4.7 4.0
##  [523] 4.5 4.3 3.6 4.2 4.0 4.3 5.7 3.7 4.5 4.7 4.5 3.7 4.2 4.5 5.4 4.1 5.6 5.5
##  [541] 4.8 4.4 2.0 3.1 4.4 5.0 5.2 4.3 5.2 3.8 5.0 5.0 3.9 5.6 3.9 4.0 3.4 4.8
##  [559] 4.8 3.9 5.4 4.0 5.8 3.8 4.9 5.2 6.4 5.4 4.0 5.3 4.7 4.4 5.1 4.9 5.2 4.8
##  [577] 4.4 5.0 3.0 4.7 4.2 5.1 4.9 4.0 5.2 4.2 5.2 4.8 3.8 6.0 3.7 4.2 4.9 3.5
##  [595] 4.1 3.5 4.7 3.6 4.3 3.7 5.1 6.2 2.9 4.3 6.1 5.4 3.5 5.6 3.1 3.5 4.0 5.4
##  [613] 4.2 3.1 4.0 4.4 4.2 5.9 5.2 2.4 3.4 3.0 5.1 3.0 6.4 4.0 5.3 6.3 5.0 4.5
##  [631] 3.3 4.3 4.7 5.3 4.5 5.3 2.2 4.9 4.3 4.4 4.5 3.2 5.2 4.6 3.8 4.0 4.5 3.7
##  [649] 4.1 4.0 5.8 5.2 5.2 2.7 4.5 3.1 4.5 4.8 5.7 4.5 4.2 3.4 5.0 4.8 4.2 4.4
##  [667] 4.3 5.0 5.6 4.3 4.9 4.6 5.1 4.3 3.8 3.7 5.1 4.7 3.7 4.2 5.2 4.8 3.9 4.5
##  [685] 3.1 5.8 3.8 4.5 5.2 3.4 4.6 4.6 4.3 4.4 3.5 4.3 4.4 3.5 4.0 5.4 4.7 6.5
##  [703] 4.5 4.6 6.0 3.6 4.3 4.2 3.5 4.5 3.8 3.1 4.2 4.2 4.0 4.5 3.7 5.0 4.9 4.2
##  [721] 3.9 6.6 4.1 3.7 3.5 4.8 3.4 5.1 4.7 5.0 4.2 4.8 5.8 3.5 3.7 3.7 4.8 3.9
##  [739] 2.4 6.0 4.6 4.8 5.0 4.2 3.9 4.0 5.6 4.3 5.7 4.6 5.7 5.3 3.3 3.2 5.6 5.4
##  [757] 4.7 4.1 4.3 3.5 4.4 5.6 5.7 3.8 5.7 4.0 4.3 5.2 4.2 3.8 4.7 4.8 4.8 3.1
##  [775] 5.1 4.6 5.1 4.9 5.6 2.4 4.4 4.2 3.6 3.0 6.1 4.7 4.3 5.8 3.9 2.6 5.9 3.4
##  [793] 5.6 4.2 5.2 3.6 4.5 5.5 5.1 5.5 4.4 5.5 6.0 4.6 6.0 4.9 4.5 1.8 3.3 3.9
##  [811] 4.7 6.7 4.0 2.8 5.0 4.2 4.5 3.0 5.0 3.9 5.2 4.0 3.5 5.6 4.7 3.8 3.6 4.3
##  [829] 5.2 4.7 5.3 3.8 5.0 5.3 4.2 3.3 5.1 4.0 4.5 4.5 4.7 3.6 3.4 6.0 6.0 4.2
##  [847] 2.9 4.6 4.1 4.3 4.9 4.8 4.6 5.4 2.7 4.8 4.5 5.8 3.5 6.3 5.1 3.9 4.9 4.9
##  [865] 4.5 5.9 4.9 5.3 4.4 4.9 4.1 4.5 4.6 4.5 3.7 5.0 4.3 6.0 5.6 7.7 5.9 4.6
##  [883] 5.5 3.7 4.8 6.1 6.5 3.6 4.7 4.7 4.1 5.2 4.5 5.5 6.4 3.2 4.5 4.1 5.2 5.1
##  [901] 6.1 5.3 4.9 5.4 5.5 3.8 5.2 2.6 3.7 6.4 5.3 7.2 4.9 5.9 3.9 4.5 5.4 4.0
##  [919] 3.0 4.8 5.1 3.9 4.0 5.7 4.6 3.4 4.8 5.6 3.9 2.8 3.5 5.2 4.9 3.9 3.5 4.4
##  [937] 3.6 3.3 3.3 4.0 4.6 4.3 4.8 4.1 3.3 3.1 5.1 4.3 5.4 5.7 2.9 4.5 4.9 3.6
##  [955] 4.6 3.9 4.4 3.9 3.5 3.9 4.3 4.4 4.2 5.7 5.4 5.7 5.1 4.1 4.5 4.3 6.7 4.1
##  [973] 6.2 5.5 4.6 4.0 4.2 4.2 4.0 4.3 4.9 5.2 3.4 2.7 3.9 3.3 3.8 5.0 5.4 4.7
##  [991] 3.9 5.8 5.2 4.2 4.9 5.1 4.4 3.0 3.7 2.1
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
## 2.7000 6.2025
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
##    [1] 3.4 1.8 4.6 3.9 5.4 4.1 4.7 4.3 4.4 6.0 5.2 5.6 4.5 5.4 4.2 4.7 3.8 6.2
##   [19] 1.9 5.3 6.8 3.8 3.5 5.0 5.5 4.3 5.1 4.0 4.5 6.0 4.3 4.5 4.9 5.1 4.6 4.5
##   [37] 3.1 4.8 2.9 5.6 5.0 4.9 3.9 4.2 2.8 5.0 3.3 4.6 2.7 4.6 4.0 3.8 4.8 3.4
##   [55] 3.7 4.4 4.8 6.1 4.2 4.6 3.9 4.0 3.1 4.0 3.6 4.3 4.5 5.8 4.3 3.4 5.6 4.9
##   [73] 3.5 4.4 4.6 2.7 3.0 5.2 4.5 6.4 4.3 4.2 4.2 5.0 4.8 6.1 4.7 4.0 4.3 5.5
##   [91] 5.5 4.7 5.2 4.8 5.2 4.4 4.2 4.2 4.2 3.8 5.3 6.2 5.2 3.9 3.9 5.5 4.2 2.9
##  [109] 2.5 4.8 3.6 5.2 5.1 3.3 4.4 5.6 4.0 3.8 4.8 5.8 4.0 4.9 4.2 3.8 3.7 4.5
##  [127] 4.0 5.3 4.8 3.4 5.1 4.2 6.4 4.4 5.2 4.8 4.7 3.8 3.8 4.5 4.2 2.9 3.5 4.4
##  [145] 6.8 6.5 6.1 4.2 4.6 4.5 5.0 5.0 4.9 5.6 4.5 5.0 5.1 4.3 3.8 3.7 5.2 4.4
##  [163] 4.0 4.9 4.2 3.9 5.1 3.3 2.7 4.7 5.5 3.2 5.7 4.6 4.0 5.5 3.5 4.1 4.5 4.7
##  [181] 5.4 4.3 5.5 4.0 5.7 5.3 4.5 5.4 5.5 5.2 6.1 4.2 4.5 3.7 5.0 3.3 4.7 4.7
##  [199] 3.1 4.7 4.6 5.0 4.6 3.8 3.8 5.3 6.0 6.6 4.0 4.9 4.1 4.5 5.3 3.8 3.3 3.8
##  [217] 5.6 3.4 4.2 3.4 5.8 4.8 4.5 3.5 4.3 5.5 5.2 3.5 4.8 4.9 5.6 3.3 4.2 6.2
##  [235] 4.0 5.3 3.1 6.3 7.0 3.3 4.4 5.1 6.1 2.9 4.0 5.4 4.3 5.6 4.5 4.7 5.0 4.0
##  [253] 2.8 3.9 5.2 5.7 4.1 4.9 2.3 4.6 3.3 4.6 6.0 4.4 6.1 4.7 3.7 6.1 5.3 4.2
##  [271] 4.3 2.5 5.8 3.7 3.4 4.2 2.7 5.0 3.3 4.2 3.8 4.1 4.5 4.0 3.6 4.9 3.6 3.8
##  [289] 4.6 5.4 3.4 3.7 6.8 3.3 5.3 4.9 6.4 3.9 3.3 5.4 5.5 4.9 4.7 4.2 3.9 3.7
##  [307] 4.0 5.2 4.7 4.6 3.3 4.4 3.9 5.1 4.7 4.4 5.9 6.5 4.2 6.1 5.0 4.0 3.7 5.2
##  [325] 4.5 3.9 4.6 6.7 5.2 3.1 4.6 4.2 4.4 4.8 5.3 4.4 4.7 3.6 2.9 2.1 4.1 4.6
##  [343] 4.7 5.0 5.0 6.7 4.7 3.3 4.5 3.5 4.6 4.3 3.8 2.6 2.6 6.0 4.8 4.8 5.1 5.0
##  [361] 5.8 4.7 4.2 4.6 4.7 4.3 5.5 2.9 4.8 4.0 4.6 4.1 5.9 5.6 4.6 4.1 5.2 5.5
##  [379] 4.1 3.2 5.1 4.3 4.3 6.1 4.1 4.8 4.5 5.5 4.1 5.4 3.5 2.8 4.4 3.8 3.8 4.5
##  [397] 6.8 4.6 4.6 6.0 4.9 6.0 5.1 3.5 3.8 5.1 6.1 5.0 5.6 5.0 4.3 6.3 4.0 4.0
##  [415] 3.8 4.8 4.8 6.5 3.9 4.4 3.6 4.4 6.3 4.1 4.8 4.3 5.3 4.8 4.4 4.2 4.6 4.7
##  [433] 5.1 5.1 4.7 5.3 6.0 4.6 3.9 3.3 4.5 3.9 4.4 5.4 4.4 3.9 3.9 5.9 2.6 4.3
##  [451] 4.7 4.6 3.7 5.0 3.6 4.6 4.8 4.0 3.9 4.3 4.6 4.5 6.0 3.3 3.0 4.8 6.4 4.9
##  [469] 4.3 4.1 4.5 3.9 6.2 4.3 4.0 4.3 3.6 5.1 4.1 3.8 4.0 4.9 4.0 5.2 3.4 4.1
##  [487] 5.5 3.5 3.3 4.4 2.5 4.4 4.9 3.4 5.8 4.6 5.4 4.0 3.1 2.6 3.0 3.5 4.9 4.6
##  [505] 4.9 5.4 2.0 5.9 4.1 4.7 4.5 5.0 4.5 4.4 2.5 4.0 3.4 6.3 5.2 4.1 2.6 4.6
##  [523] 3.1 5.4 4.6 3.6 5.0 3.1 4.9 4.1 3.8 5.0 3.7 4.5 4.3 5.4 4.3 6.0 4.7 4.9
##  [541] 4.2 5.8 4.4 3.0 4.4 3.6 3.0 3.4 4.2 5.1 3.9 4.3 4.5 4.9 4.7 4.4 3.7 3.6
##  [559] 5.3 4.1 5.1 4.6 4.4 4.9 3.8 4.2 3.8 4.8 5.0 3.9 4.3 5.3 6.3 3.5 4.5 4.4
##  [577] 4.0 3.1 2.3 4.5 3.7 4.9 3.3 3.3 6.1 4.9 3.7 4.3 4.8 3.8 4.8 4.1 5.8 5.2
##  [595] 2.4 4.7 3.7 5.8 2.4 5.1 4.5 4.3 3.2 5.4 5.2 5.2 2.6 4.9 5.8 4.5 4.3 5.7
##  [613] 4.1 4.6 2.9 3.4 6.4 3.8 5.6 4.4 5.5 4.1 4.1 4.1 3.1 5.5 5.0 3.5 5.5 5.4
##  [631] 3.3 3.5 3.9 3.4 2.1 4.6 5.2 3.5 3.3 3.2 5.6 4.4 4.0 5.7 5.5 4.7 2.6 6.8
##  [649] 4.8 4.3 5.0 4.9 3.0 6.4 5.2 5.1 4.5 4.8 5.1 3.9 3.9 5.2 5.9 5.3 4.3 4.4
##  [667] 4.9 5.1 6.9 3.8 4.6 4.5 3.7 6.6 5.2 4.0 4.8 6.5 2.0 3.5 5.3 4.8 5.0 3.6
##  [685] 4.3 6.0 3.5 5.2 3.6 3.3 5.6 4.0 4.0 5.1 5.7 3.8 4.3 5.3 5.3 5.6 4.5 4.9
##  [703] 6.3 3.8 5.6 3.0 3.4 5.8 3.2 5.2 4.0 5.2 4.9 4.3 4.6 3.3 5.7 4.2 5.6 3.4
##  [721] 3.2 4.6 5.9 5.0 4.0 4.5 4.6 7.2 3.5 2.6 4.0 4.9 6.1 3.6 4.6 5.0 6.8 5.3
##  [739] 4.7 3.1 5.1 5.8 4.1 4.7 4.5 6.3 4.0 3.9 5.2 4.4 4.9 3.6 5.8 5.8 2.8 5.2
##  [757] 3.6 5.6 4.1 5.6 3.3 4.5 4.0 3.9 4.7 4.0 2.5 3.7 3.3 4.9 5.5 3.9 5.3 3.6
##  [775] 3.2 5.6 6.0 4.9 3.8 3.7 5.2 4.2 4.6 5.0 5.0 3.2 4.5 4.2 5.4 3.8 4.9 4.6
##  [793] 4.3 3.1 5.3 5.7 2.8 4.7 3.7 3.8 4.1 7.3 4.3 4.8 3.4 3.8 4.3 4.0 5.3 4.5
##  [811] 5.8 6.2 4.3 3.3 3.2 4.0 4.1 4.6 4.4 2.6 4.2 4.5 4.9 4.2 5.5 5.0 4.2 4.9
##  [829] 6.0 4.2 4.4 4.8 6.7 4.1 5.1 5.9 3.1 5.3 6.0 2.8 4.9 4.4 3.7 5.5 4.2 3.8
##  [847] 5.5 5.7 5.0 5.2 4.6 4.7 5.1 5.5 4.2 4.9 3.3 4.3 5.8 5.6 4.4 4.9 4.9 3.8
##  [865] 6.3 5.5 4.1 5.8 5.2 6.7 4.1 6.2 4.1 5.5 3.3 3.4 3.8 3.9 2.7 3.4 2.6 4.8
##  [883] 5.6 5.9 4.8 3.7 6.6 5.6 4.0 3.9 4.4 3.5 3.7 4.2 3.2 5.8 6.2 5.7 4.4 2.7
##  [901] 4.0 5.1 6.5 3.8 5.7 4.3 5.5 4.3 4.8 5.8 5.7 4.1 3.4 4.9 4.4 2.5 4.1 5.2
##  [919] 3.8 4.4 4.6 4.7 4.6 3.6 3.6 5.1 2.3 5.8 5.7 3.9 6.2 4.3 4.1 5.5 3.3 3.2
##  [937] 6.2 5.9 5.7 5.8 5.2 4.1 4.8 4.5 5.8 2.9 5.6 4.8 3.5 4.6 3.2 5.0 4.8 5.1
##  [955] 3.4 3.6 4.1 4.2 4.9 4.5 5.1 4.4 4.4 4.1 5.5 4.9 4.5 5.7 4.6 5.2 5.6 5.7
##  [973] 4.2 5.6 4.4 3.0 5.1 5.2 3.4 4.5 5.1 4.0 3.8 5.5 5.7 3.8 4.4 5.5 4.5 3.7
##  [991] 4.4 4.1 4.5 4.7 4.9 3.9 4.9 5.1 5.0 4.3
## 
## $func.thetastar
## [1] 0.0252
## 
## $jack.boot.val
##  [1]  0.54724638  0.39060773  0.32346041  0.24971591  0.02486034 -0.03187500
##  [7] -0.16239782 -0.21971014 -0.36257310 -0.55417957
## 
## $jack.boot.se
## [1] 1.008323
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
##    [1] 4.1 3.5 4.2 4.3 5.2 6.5 5.5 5.7 4.9 6.3 3.8 5.0 4.7 5.3 5.8 4.8 5.0 4.1
##   [19] 4.6 3.7 4.0 4.8 5.0 4.8 5.0 4.4 3.9 3.8 4.4 4.9 4.4 4.9 2.8 3.9 3.2 6.4
##   [37] 4.7 5.3 3.7 4.4 4.6 3.9 4.9 5.3 4.1 3.8 3.8 5.0 3.4 3.6 5.2 3.0 3.7 4.9
##   [55] 2.8 5.4 1.6 3.6 4.1 5.6 4.0 5.2 5.9 4.5 3.9 5.5 4.2 5.7 3.1 5.9 5.1 4.5
##   [73] 5.9 3.8 5.8 5.2 4.4 3.8 4.7 3.8 4.7 3.3 4.9 5.0 4.6 3.0 2.5 4.6 3.6 3.6
##   [91] 6.0 4.3 4.0 4.6 3.9 4.5 3.6 5.0 4.4 5.5 3.8 5.0 2.8 4.4 4.6 5.4 5.2 2.3
##  [109] 4.9 4.0 4.6 5.1 3.8 5.1 2.8 5.6 4.1 4.1 3.0 5.1 2.1 3.2 4.0 4.9 5.3 4.4
##  [127] 4.6 4.6 4.6 3.3 5.0 6.3 4.7 4.6 5.7 3.5 4.8 3.9 4.3 4.1 4.2 6.2 3.9 5.4
##  [145] 5.6 4.6 6.1 4.9 5.2 3.2 3.0 5.4 5.1 2.7 4.8 4.6 4.8 6.1 4.8 4.2 5.0 6.0
##  [163] 4.4 4.3 3.0 4.4 4.4 4.0 2.8 4.2 5.8 3.7 5.3 5.2 5.4 4.4 4.6 4.1 4.8 3.7
##  [181] 2.9 6.0 3.0 5.0 3.6 4.4 5.8 3.8 4.7 5.9 4.7 4.5 5.8 5.8 3.2 4.4 4.4 5.3
##  [199] 3.4 4.2 4.9 5.3 4.6 5.3 3.2 4.1 4.8 4.6 4.1 5.1 3.3 3.8 3.7 5.6 4.5 3.4
##  [217] 6.5 3.3 3.4 5.0 5.8 4.1 4.4 5.4 2.8 2.9 4.9 4.3 5.3 5.6 3.1 5.0 4.9 4.5
##  [235] 5.1 5.5 4.0 4.2 6.4 5.6 4.5 5.3 4.1 3.8 3.8 5.1 4.0 4.6 4.3 4.6 4.7 5.2
##  [253] 5.8 4.3 4.2 6.4 5.0 4.5 4.7 5.3 4.5 4.2 5.3 6.3 3.4 5.5 5.4 4.6 4.8 3.1
##  [271] 6.2 5.3 2.7 5.1 4.1 2.6 5.2 4.4 5.5 4.2 4.2 3.0 4.6 4.4 5.4 3.4 4.9 5.0
##  [289] 3.7 4.8 4.1 4.7 3.5 3.4 4.5 4.3 4.2 5.1 3.4 5.8 5.6 5.2 5.0 4.4 4.3 4.0
##  [307] 3.7 4.4 5.1 3.3 5.5 6.1 4.2 3.2 4.3 4.3 6.0 5.1 5.1 5.7 4.7 5.6 4.1 5.6
##  [325] 4.0 4.7 4.1 4.4 4.5 3.6 4.2 6.2 4.9 4.0 6.0 5.8 5.1 4.0 5.3 3.2 5.6 5.8
##  [343] 3.7 4.7 3.1 5.1 4.4 5.7 4.9 4.8 3.8 5.7 5.4 5.5 3.6 4.1 3.3 4.6 3.4 3.7
##  [361] 4.6 3.9 4.3 5.2 4.3 4.9 3.9 4.1 5.0 3.3 4.9 2.7 3.7 4.2 4.9 3.0 3.8 5.1
##  [379] 4.4 5.7 4.4 5.2 5.6 5.4 5.7 3.5 3.1 4.2 4.8 5.5 3.4 3.3 5.0 5.1 3.4 5.4
##  [397] 4.4 4.1 5.2 3.7 5.5 4.1 4.0 4.7 3.2 4.9 4.3 4.4 5.0 6.9 4.6 5.1 6.8 3.7
##  [415] 5.7 2.7 5.9 4.4 4.0 5.2 2.2 4.1 4.2 3.3 4.5 4.2 4.7 5.7 5.2 4.4 4.0 5.1
##  [433] 3.4 3.4 4.9 6.2 5.6 4.7 3.9 5.3 3.3 5.3 5.6 3.7 5.5 4.7 4.9 4.4 4.9 3.1
##  [451] 4.0 4.7 4.7 5.4 6.0 4.0 2.6 6.3 6.3 5.3 5.3 4.8 4.7 4.6 4.1 3.3 5.4 2.8
##  [469] 5.5 5.8 3.7 4.1 5.5 4.9 5.8 4.4 3.4 5.1 4.8 4.6 5.7 3.6 4.5 5.5 3.9 4.6
##  [487] 3.8 3.1 4.0 3.2 3.9 4.5 4.8 3.7 5.3 4.4 5.5 4.4 4.8 3.6 3.9 5.3 5.9 4.9
##  [505] 3.6 4.3 4.1 4.9 4.3 3.7 5.2 4.8 5.7 5.1 6.2 4.9 4.7 4.4 3.8 2.7 5.8 4.0
##  [523] 4.5 5.2 6.4 3.5 4.5 5.1 3.2 4.8 4.8 3.4 5.5 5.0 3.1 4.7 3.5 5.9 4.6 6.1
##  [541] 5.3 5.1 5.7 4.4 3.8 3.4 4.9 4.1 1.8 5.3 4.3 4.5 6.1 5.0 3.7 4.8 3.6 3.4
##  [559] 4.6 5.9 5.6 4.4 4.3 5.6 3.2 3.6 3.6 4.4 6.2 3.9 4.4 3.2 5.1 5.6 2.6 5.2
##  [577] 3.2 4.7 4.3 5.4 3.4 3.8 4.8 4.4 3.0 4.3 3.9 4.9 5.2 3.0 4.3 4.7 4.9 4.4
##  [595] 3.3 5.2 3.8 4.8 2.6 4.7 4.1 2.0 3.8 4.4 5.2 5.9 4.4 4.2 4.5 3.7 4.0 5.4
##  [613] 5.8 3.8 2.5 3.7 5.5 3.3 5.4 4.0 3.4 5.3 4.5 5.7 5.1 4.2 5.2 3.2 7.6 4.4
##  [631] 4.4 4.4 5.2 4.3 3.4 5.0 4.6 4.5 4.5 4.1 3.9 4.8 5.2 6.0 3.6 4.2 5.1 4.3
##  [649] 2.5 4.5 3.8 4.5 4.8 3.4 4.8 4.6 4.7 4.7 3.2 3.3 5.4 5.2 4.9 4.4 5.5 4.6
##  [667] 4.3 5.3 6.2 3.5 5.8 4.3 3.7 5.2 3.6 5.6 5.7 6.6 4.7 3.9 3.4 4.1 3.8 4.6
##  [685] 4.9 3.6 5.7 2.7 5.2 3.6 3.4 4.1 3.3 5.3 5.2 4.3 7.5 5.4 4.5 3.6 4.3 3.9
##  [703] 3.9 5.3 4.5 5.2 5.4 3.8 3.9 4.4 3.8 4.8 4.1 5.4 5.1 4.4 4.2 5.1 4.1 4.9
##  [721] 4.1 3.2 4.7 4.0 5.5 4.6 6.7 5.1 4.4 6.0 4.9 5.1 4.8 5.4 3.8 4.5 3.5 4.3
##  [739] 3.7 4.6 2.9 4.2 4.0 3.6 4.3 4.0 2.9 3.9 5.0 5.3 5.6 4.4 3.4 5.0 5.1 6.1
##  [757] 4.7 5.1 4.2 3.9 3.4 4.9 4.6 2.7 3.7 4.4 3.8 3.2 5.1 3.3 5.1 4.2 4.3 1.5
##  [775] 4.9 5.5 3.6 3.8 3.6 5.1 5.0 5.7 4.3 5.9 4.7 5.1 2.9 3.0 4.4 3.7 4.3 5.0
##  [793] 4.8 5.0 4.0 3.5 5.5 3.0 4.6 6.4 1.4 4.1 3.7 3.7 5.0 3.6 5.4 5.3 3.4 5.3
##  [811] 4.9 5.3 3.5 3.8 4.6 4.6 5.3 5.6 3.9 4.4 5.2 5.1 6.2 5.2 5.4 6.1 3.8 2.6
##  [829] 3.9 5.2 5.4 6.1 5.6 5.0 4.4 4.4 4.8 5.3 4.5 4.4 4.1 5.8 4.0 4.5 6.7 3.1
##  [847] 4.2 3.9 4.4 5.1 3.3 4.7 2.5 4.4 4.7 3.9 2.2 5.2 6.0 3.7 4.1 3.6 5.9 5.4
##  [865] 4.8 6.1 5.4 5.6 6.0 5.7 5.4 5.5 4.0 3.8 3.6 3.3 4.7 4.5 5.2 3.2 3.7 6.0
##  [883] 3.9 4.4 5.2 5.6 3.5 3.7 5.0 4.3 5.9 3.2 4.4 5.8 5.9 5.1 5.4 4.9 3.9 3.6
##  [901] 4.8 5.2 4.5 3.8 5.0 4.0 3.3 3.7 4.3 4.1 4.6 3.8 2.9 4.9 5.4 4.1 4.5 3.4
##  [919] 3.0 4.1 4.0 4.0 2.7 3.8 4.3 6.1 6.4 4.0 4.0 4.6 3.8 3.2 4.9 2.0 4.4 6.0
##  [937] 3.1 5.2 4.4 5.1 4.6 4.5 4.2 4.8 4.9 5.1 4.5 3.7 4.2 4.9 3.6 5.7 5.2 4.8
##  [955] 4.8 3.9 5.0 4.0 4.3 4.1 5.3 3.7 5.1 4.7 4.4 4.5 5.3 5.9 5.0 3.5 4.8 3.8
##  [973] 5.8 3.5 4.4 3.4 4.5 2.8 4.0 4.3 4.5 6.5 4.5 4.5 4.8 3.6 3.8 4.7 3.8 5.1
##  [991] 3.3 3.8 4.8 4.8 3.3 4.3 6.2 5.5 3.8 3.5
## 
## $func.thetastar
## 72% 
## 5.1 
## 
## $jack.boot.val
##  [1] 5.5 5.4 5.4 5.3 5.2 5.0 4.9 4.7 4.7 4.5
## 
## $jack.boot.se
## [1] 0.9967949
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
## [1] 1.011073
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
##    5.229714   10.560027 
##  ( 2.268131) ( 4.807051)
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
## [1]  0.004408353 -0.184470814  0.519884693  0.557071417  0.913748487
## [6] -1.074346894
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
##    [1]  0.286119467  1.654708875  0.683439630  1.157645216  1.020740430
##    [6]  0.789717136  0.923975594  0.962083304  1.229369175  0.601888364
##   [11]  0.957831471  1.038134751 -0.697719463  0.250514704  0.199390325
##   [16]  0.230837664  0.311066963  1.013532891  0.716901698  1.210733371
##   [21]  2.396347615  1.473561387  0.516069359  0.600288131  0.588941035
##   [26]  1.011073439  0.030220595  1.357242751  0.682502435  0.850378357
##   [31]  1.009067881  1.772257922  1.627477698  1.319071866  1.226355566
##   [36]  0.909077660  0.448770998  0.460728878  0.181848632  1.872590078
##   [41]  2.438595759  0.312321582  1.652474061  1.457156861  0.960753790
##   [46]  0.910914978  0.933306209  1.179952424  0.769801188  0.540494803
##   [51]  1.787866402  0.980042681  1.325036483  0.142204673  1.205857400
##   [56]  0.058132798  1.208760190  0.531839536  1.027413775  1.131299704
##   [61]  0.485790757  0.621666943  0.745852070  0.818369056  0.587379509
##   [66]  2.097187279  1.982000568  2.020508381  1.417402850  1.850709778
##   [71]  0.364927313  1.819631056  0.779113648  1.971266283  0.818928167
##   [76]  0.783357997  0.599424016  1.081657639  2.430639188 -0.291218741
##   [81]  1.747293450  0.939544965  2.018256869  0.720995536  0.106663144
##   [86]  0.366772315  0.602996886  1.304959436  1.576652649  2.231804605
##   [91]  2.064152784  0.952952222  1.066995728  1.851895763  0.812842873
##   [96]  0.613527397  0.658498701  1.388876457  0.818240607  2.031984832
##  [101]  1.585304999  0.520342245  1.088343217 -0.005656343  1.842760433
##  [106]  1.109596139  2.540110450  2.406925693  0.663997670  0.531158883
##  [111]  1.233838824  0.591304111  0.318534166  1.040524514  1.660608152
##  [116]  0.635755142  1.670560597  0.804779177  1.284211840  1.768518825
##  [121]  1.602092252  1.094887419  0.340312754  2.288561244  1.017610692
##  [126]  0.099680279  0.728707402  0.195442587  2.235045260  1.630805731
##  [131]  1.331859619  0.378069331  0.198561604  0.780418764 -0.127534291
##  [136]  0.238407236  0.971263048 -0.114948905  1.930753788  0.922611498
##  [141]  0.378230874  1.499897685  1.023386961  1.588396934  1.921963798
##  [146]  0.589515977  1.370248414  0.400870176  0.778889134  2.361247706
##  [151]  1.118207199  1.506078733  1.252115401  1.911518846  0.625481013
##  [156]  1.219881244  2.032853086  0.789347426  0.179784392  0.674580660
##  [161]  1.949470759  0.538350217  1.287584214  0.981029500  0.680676060
##  [166]  1.264201280  0.947151678  0.128187640  1.104141171  1.267433341
##  [171]  1.919884190  0.627104856  0.239880308  0.549901529  1.094997397
##  [176]  0.375011422  0.920936560  1.410704966  1.320347220  0.336520810
##  [181]  1.657476498  1.565002957  0.060203394  1.685585185  0.905150600
##  [186]  0.160354003  0.774012377  0.635913908  0.404323343  0.536360598
##  [191]  1.307280229 -0.158628743  0.919044130  0.632287978  0.056340577
##  [196]  1.910303850  1.190696051  0.960812614  1.006678814  0.555562656
##  [201]  0.725011020  0.013397731  0.987875813  1.368426171  1.739518694
##  [206]  0.718007268  0.450290785  1.584217234  0.543949992  1.773492431
##  [211]  1.120296149  1.786024202  1.616370423  1.413915146  0.316544219
##  [216]  0.544416895  0.822395147  1.008902210  0.922474954  0.815962330
##  [221]  1.100271457  0.982324091  0.348340099  0.054662921  1.364978666
##  [226]  1.122044072  0.888715600  1.336589527  0.679203677  0.211405871
##  [231]  1.248309316  0.740502865  1.618916934  1.386027517  1.020740430
##  [236]  1.006759122  1.654622091  0.532081484  0.394744690  1.363803631
##  [241]  1.323139613  0.311623011  1.084238103  0.870103880  0.916474704
##  [246]  1.143232546  1.231227038  0.001824885  0.689978959 -0.591077337
##  [251]  0.654107815  1.048602733  0.397549384  0.624489813  0.591135707
##  [256]  1.011073439  2.255743789  1.249772239  0.523786392  0.600292779
##  [261]  0.932680243  0.332853566  0.772096761  0.537750585  1.769893949
##  [266]  0.320464588  1.083202486  0.047900968  1.145402049  0.205067525
##  [271]  1.889655788  2.381815693  0.583038472  0.033560537 -0.192571779
##  [276]  1.139618208  2.421347256  0.600769131  0.672593586  0.298097765
##  [281]  1.444018594  1.046610829  0.509744247  1.097221667  0.769754881
##  [286]  0.384788672  1.068912792  0.910861512  1.178595506  1.144558119
##  [291]  1.195480562  0.388318658  1.001972758  1.268019719  0.577811823
##  [296]  0.529512888  0.326122227  0.374433644  0.260841067  1.126621666
##  [301]  0.518036227  0.258552633  1.220210311  0.497669045  1.385981900
##  [306]  0.965622941  1.178716975  0.200144721  1.317995222  1.764396388
##  [311]  1.202379880  0.560152799  0.246918167  0.918425172  0.130427196
##  [316]  1.028315652  0.609030242 -0.418740785  1.310316682  1.232744172
##  [321]  1.147975174  1.183658732  1.338597464  0.354688532  0.617180677
##  [326]  0.460198704  0.870141327  0.504569719  1.723117635  0.531004286
##  [331]  1.141809450  0.588462565 -0.053635845  0.657966232  1.432468612
##  [336]  1.391386810  0.525744262  1.618644835  0.259427760  0.629749974
##  [341]  0.567735005 -0.020547494  0.360293330  0.252701400  0.616635496
##  [346]  1.358711252  1.551182189  0.909241626  0.517499348  0.603560485
##  [351]  1.294553304  1.463475393  1.268129344  0.805904324  0.133404134
##  [356]  0.491981654  1.193708081  1.157533760  1.658362797  1.753656944
##  [361]  1.142741845  0.614513228  1.801911815  0.404109221  0.974296765
##  [366]  1.441614737  0.778923480  0.260521288  1.221073667  1.137234901
##  [371]  0.938633577  1.259710582  1.421606227  1.861589020  1.882962957
##  [376]  0.984614070  1.476926762  0.595190992  1.602004220  1.444670420
##  [381]  0.247170997  1.469967731  0.666095485  0.999727996  0.543949992
##  [386]  0.901199141  1.329833710  0.993871106  1.103977566  0.359261377
##  [391]  0.405314836  0.241858517  0.758865042  1.098928065  0.988337780
##  [396]  1.109497591  1.616121019  0.738599375  0.742385052  0.851114448
##  [401] -0.158983112  0.543781128  0.993905198  1.739518694  0.823610453
##  [406]  0.746160977  1.340119873  1.026084994  0.584009714  0.740816515
##  [411]  0.794367646  0.254565487  1.196893882  1.913358429  1.386723022
##  [416]  1.031707975 -1.044708134  1.996437013  1.069290305  1.005404998
##  [421]  1.207661555  2.083751870 -0.285161283 -1.396006449  1.707663644
##  [426]  1.936378757  1.833692469  0.494125397  0.522854651  1.668901191
##  [431]  1.186165832  1.929121351  0.987137792  0.844222009  0.360694486
##  [436]  0.661164364  0.799667297  0.814855220  0.753088240  0.627084364
##  [441]  0.651415044  1.242847325  0.592386012  1.035932735  1.248265076
##  [446]  0.359261377  0.949901330  0.674183917  1.128683341  0.468288264
##  [451]  0.453806188  1.009962679  2.437349941  0.450290785  1.227317308
##  [456]  0.129865668  0.201350129  0.742945229  0.339894428  1.335343245
##  [461]  0.257334371  0.805000583 -0.202585358  0.932680243 -0.292413363
##  [466]  0.408466325  0.101074013  0.994859545  1.500607813 -0.146429699
##  [471]  2.182126899  0.182805097  1.490207148  1.244818170  0.194721467
##  [476]  2.249131336 -0.428374706  1.210515750  1.119034721  1.553503257
##  [481]  0.579633230  1.108871989  0.842502758  0.963364429  1.178588737
##  [486]  0.917089197  0.741937477  0.930550516  0.004278763  0.489927135
##  [491]  0.538219290  0.945972077 -0.009928413  1.155691708  0.866159325
##  [496]  1.256862793  2.050455463  1.829529872  0.486866389  0.310584286
##  [501]  0.589866485  0.884203307  0.807829607  0.641839660  1.341282468
##  [506]  0.378019648  0.564088676  0.547830738  0.542437988  1.038472588
##  [511]  1.646649446  1.668766749  1.268366672  0.575522549  0.728213025
##  [516]  1.816114295 -0.063066584  1.210178981  0.460389503  1.970880297
##  [521]  1.956539171  0.366726581  1.950018514  2.323476572  0.630581508
##  [526]  1.967936912  1.452196622 -0.088342950  1.988216223  0.217083723
##  [531]  0.911563306  0.687773800  1.100945516  0.129271463  1.127845256
##  [536]  1.203156576  1.485595487  0.549608190  0.375364692  1.883793035
##  [541]  0.849173978  1.436354749  0.284088112  0.326333714  0.920146774
##  [546]  0.904151991  1.371851247  1.230612119  2.212316898  1.046958091
##  [551]  0.655268154  0.708305118  0.523694020  1.480724626  0.446438761
##  [556]  1.279929398  0.884186221  1.171559901  1.935905962  1.532230471
##  [561]  0.245411014  0.193746721  1.007956062  0.617556684  0.929470629
##  [566]  1.355222738 -0.112542107  1.273401676  0.471472405  0.900939777
##  [571]  1.886532200  0.704461391  0.285232037  1.814826421  1.773789571
##  [576]  0.906194486  0.700212606  1.239423218  1.416435588  0.449369383
##  [581]  0.179963380  0.124302606  0.328546358  1.274731665  0.985457951
##  [586]  1.007840332  1.598363729  0.978243053  1.083801794  0.101624173
##  [591]  1.191966894  1.354427418  1.076404861  1.019674822  0.543259743
##  [596]  1.203269908  0.459981794  0.725956902  1.645586358  0.680576302
##  [601]  0.672593586  1.132192136  0.129271071  0.168317987  0.232582438
##  [606]  0.627455514  1.784208635  0.213963683  0.137198400  0.819616083
##  [611]  0.627193366  1.088092095  1.445975722  0.910636781  2.200784180
##  [616]  0.650633504  0.206973586  0.923737572  1.073642174  1.485362703
##  [621]  0.453594201  1.361217769  0.960221018  0.239425181  1.622316465
##  [626]  1.442164270  1.910439046  0.741210785  1.190850651  1.364391008
##  [631]  0.560936710  0.693818907  1.001882819  1.238771248  0.657966232
##  [636]  1.424630023  0.587587219  0.866356646  0.538293114  0.126879281
##  [641]  0.968206201 -0.757238450  2.323656376  1.218497944  0.659162846
##  [646]  0.667860610  1.575885100  1.075521513  0.379237622 -0.320057756
##  [651]  1.839825158  0.836184113  0.538217155  1.444113948  1.769857578
##  [656]  0.901659120  1.343292007  0.279280291  1.281535843  0.311126513
##  [661]  0.985833242  0.396933284  1.337938478  0.479311074  0.822239662
##  [666]  0.731826826  1.107035519  2.182728051  0.098334291  0.474554803
##  [671]  0.099680279  1.057373519  0.829073256  0.623194403  1.961902773
##  [676]  0.703233460  0.675069020  0.123859030  0.392347068  0.316432910
##  [681]  1.672944254  1.209589951  0.276311002  1.082468205  1.503622221
##  [686]  0.542532323  0.974197166  1.529776906 -0.196649227  0.774012377
##  [691]  1.754049496  0.205419597  0.670270694  0.499862422  1.133811268
##  [696]  1.165962581 -0.179939958  0.796601152  0.121893238  0.804354880
##  [701]  1.574785411  0.126850786  2.072400963  2.241382221  0.714105579
##  [706]  1.355399681  1.352710020  0.719346217  0.666634843  1.367212539
##  [711]  0.731127514  0.257315373  0.405624573  1.715193205  0.812414164
##  [716]  0.350404554  0.715571447 -0.033301273  0.896203721  0.690728014
##  [721]  0.700943355  2.403167924  1.139977956  1.367031962  0.294889467
##  [726]  0.402565442  0.232720971  1.358437491  0.538882911  1.087648090
##  [731]  2.081794355  0.955780439  0.630624882  0.273981916  0.252116017
##  [736]  1.530778916  0.686802130  0.956563472  1.725730796  1.586359348
##  [741]  0.686852385  0.278106775  1.431394976  1.181939469  1.104911405
##  [746]  1.024094293  0.236341401  0.505161753  1.324313721  0.717609320
##  [751]  1.137489136  0.897906108  1.913107054  0.745457574  2.117862500
##  [756]  0.664905135 -0.071380924  1.357930222  1.080109925  1.006333145
##  [761]  0.451160999  1.511571880  0.661874363  0.678868321  0.214418098
##  [766]  0.152856362  0.401969677  0.600618386  1.432718092  0.282975687
##  [771]  1.720486437  1.499797685  0.630739932  1.240753816  0.276582100
##  [776]  1.801510440  0.620431099  1.145298621  1.255813942  0.707206288
##  [781]  0.267489070  1.083984146  0.571570010  0.944741425  2.392389961
##  [786]  1.336048645 -0.573805408  1.059535957  1.326287232  1.350575690
##  [791] -0.083977692  0.604314477  1.674293761  1.821240618  0.102758497
##  [796]  1.019771819  0.701580552  1.094586597  0.809808825  0.382550069
##  [801]  1.968576139  0.533051388  1.066885950  1.842760433 -0.464359821
##  [806] -0.032727886  1.047611947  0.790592100  0.923708562  2.317676625
##  [811]  1.080775153  1.886719137  1.200849667  0.917764963  0.022970299
##  [816]  0.460022179  0.907582724  1.089862161  0.974685013  1.156113982
##  [821]  0.403720266  0.095846980  0.347339354  1.249631882  0.896625129
##  [826]  0.356680335  0.446254904  1.198980968  1.877941296  1.273491572
##  [831]  0.849084975  0.321143567  0.484740645  1.850653945  1.022526158
##  [836]  1.691235174 -1.052851888  1.787364859  0.806520609  2.281913764
##  [841]  0.428639722  0.690295711  0.560059339  0.240497752  0.444230711
##  [846]  0.560818067  0.889028342  1.664979777  0.334715432  0.484740645
##  [851]  1.147560449  0.722239976  0.838105891  1.513359808  1.382268364
##  [856]  0.653873378  0.627310097  0.876180260 -0.011992328  0.178312251
##  [861]  0.790456590  1.092107918 -0.128081621  0.322009530  1.099786811
##  [866]  0.626976649  0.843658986  1.825145170  1.001399546  1.577803680
##  [871]  1.198087702  0.676406576  1.200630538  0.386099705  0.591833055
##  [876]  0.820294539  0.569494986  1.011842003  1.883793035  1.948177916
##  [881] -0.575327662  0.874796353  0.883041458  0.580215898  1.095602234
##  [886]  0.767447582  0.893787928  1.229912904  1.670377006  1.190836283
##  [891]  0.850253376  1.993720088  0.366996459 -0.028629411  0.802980067
##  [896]  0.922172673  1.384106596  0.863390855  1.695437739  0.753661809
##  [901]  1.134464947  2.290100719  1.321930893  0.827431320  0.537359902
##  [906]  1.318187972  1.033222279  1.562765531  1.069036521  0.152335986
##  [911]  0.273294888  0.796606674  0.575917029  1.125515402  0.656095019
##  [916]  0.806521360  0.461484634  0.600271471  0.999330233  1.303623360
##  [921]  0.459903026  1.058849994  1.359310196  1.507315812  0.336774926
##  [926]  1.295705164  0.416123859  0.391967984  0.421402639  0.697711863
##  [931]  0.264993130  0.831782642  0.381595081  2.268809539  2.397407948
##  [936]  1.270640912  0.332553263  1.181939469  0.647834484  0.324235814
##  [941]  0.547468369  1.187403046  1.114047135  1.911606457  2.151080638
##  [946]  0.250812280  0.662539906  0.834024168  1.329649929  0.413674258
##  [951]  1.274356129  0.452175210  1.330504481  0.152453705  1.752314440
##  [956]  1.103845127  1.332429100  0.513022187  0.454635601  1.872462392
##  [961] -0.194933014  1.641772453  1.596472715  0.744168218  0.622431017
##  [966]  0.398153851  0.365891309  0.634963861  0.551192805  0.769855500
##  [971]  0.507634247  2.340500770  1.834393038  1.341420640  0.964032307
##  [976]  1.346550179  0.911885016  1.346033694  0.395211902  1.358014048
##  [981]  0.850380714  0.623778180  0.278763633  1.205226611  0.426353297
##  [986]  1.153471308  0.759744215  0.808483803  0.180175279  0.179272491
##  [991]  1.003301555  1.870022031  1.201712388  0.972485531  0.957031715
##  [996]  1.259063591  1.192256361 -0.080419334  1.212413792  0.676108507
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
##   0.49522856   0.23532540 
##  (0.07441642) (0.05261677)
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
## [1]  0.3938715 -0.3304665 -0.1060644 -0.7464052  0.2486782 -0.1202831
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
## [1] -0.0485
```

```r
sd(results2$thetastar)  #the standard deviation of the theta stars is the SE of the statistic (in this case, the mean)
```

```
## [1] 0.8849525
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
## t1*      4.5 -0.03083083   0.9044149
```

One of the main advantages to the 'boot' package over the 'bootstrap' package is the nicer formatting of the output.

Going back to our original code, lets see how we could reproduce all of these numbers:


```r
table(sample(x,size=length(x),replace=T))
```

```
## 
## 0 1 3 4 5 6 7 
## 1 1 1 1 1 4 1
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
## [1] 0.0438
```

```r
se.boot
```

```
## [1] 0.8914237
```

Why do our numbers not agree exactly with those of the boot package? This is because our estimates of bias and standard error are just estimates, and they carry with them their own uncertainties. That is one of the reasons we might bother doing jackknife-after-bootstrap.

The 'boot' package has a LOT of functionality. If we have time, we will come back to some of these more complex functions later in the semester as we cover topics like regression and glm.

