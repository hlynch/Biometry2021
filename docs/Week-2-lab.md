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
## 1 2 3 4 5 6 
## 2 2 1 1 3 1
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
## [1] -0.0632
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
## [1] 2.597461
```

```r
UL.boot
```

```
## [1] 6.276139
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
##    [1] 4.1 4.9 4.0 5.5 3.9 4.8 4.7 3.3 3.8 4.0 6.4 4.0 6.0 2.3 2.9 6.7 5.4 4.8
##   [19] 4.4 5.0 5.5 4.2 6.1 3.1 5.4 4.9 4.4 4.5 3.5 3.5 5.9 3.6 5.0 4.8 6.3 5.3
##   [37] 4.4 6.2 3.3 3.5 6.1 5.4 2.9 6.1 4.1 4.5 5.2 4.7 4.3 3.0 5.0 5.6 2.9 5.6
##   [55] 3.1 4.1 4.6 5.6 4.0 4.9 3.8 3.5 4.0 3.5 3.6 4.7 5.2 4.6 4.3 5.5 2.7 5.0
##   [73] 3.2 3.5 3.6 3.6 3.6 5.2 3.8 3.9 4.4 4.6 5.1 3.5 5.9 4.1 3.6 4.4 5.1 5.7
##   [91] 4.4 4.0 5.8 5.0 4.0 4.1 3.8 4.3 3.4 3.3 6.4 3.1 4.8 2.7 2.5 3.1 4.1 5.2
##  [109] 5.2 7.0 4.7 4.0 5.7 6.4 5.7 5.6 3.3 4.9 5.0 3.6 6.3 3.3 6.6 4.3 4.0 3.3
##  [127] 3.2 4.3 5.5 3.8 2.6 3.4 5.2 5.6 4.9 3.4 5.1 5.8 4.2 4.8 3.3 5.0 5.4 4.5
##  [145] 6.1 5.3 4.5 4.6 4.1 4.5 3.3 3.8 3.8 5.1 3.2 4.9 5.2 4.4 5.0 3.5 5.5 6.9
##  [163] 4.7 4.1 6.3 5.6 3.9 3.3 3.7 6.0 4.2 5.1 4.7 4.1 5.6 3.9 3.8 3.6 4.0 5.1
##  [181] 6.2 5.4 5.2 5.3 5.4 5.2 3.5 5.4 5.2 3.6 3.4 3.1 2.8 3.7 4.2 5.6 5.8 3.6
##  [199] 5.1 3.0 3.5 4.6 4.2 2.7 4.5 6.2 4.9 3.6 3.4 3.4 2.8 4.6 4.3 5.1 4.4 3.7
##  [217] 4.4 3.5 4.1 5.1 3.6 3.8 5.6 6.1 4.5 4.4 5.4 4.4 6.1 5.4 5.2 3.6 3.9 5.0
##  [235] 5.9 3.4 5.2 5.2 4.5 6.5 3.8 4.7 4.7 4.8 4.2 5.2 3.9 4.4 5.1 4.5 5.4 2.5
##  [253] 3.8 3.6 5.1 4.0 4.4 4.5 5.7 4.6 3.6 5.7 5.8 4.8 4.3 5.8 4.5 5.3 4.6 4.2
##  [271] 5.1 6.4 3.5 5.2 5.9 3.5 4.3 4.3 6.4 4.9 4.3 4.8 4.7 5.2 4.3 3.2 3.5 3.3
##  [289] 3.7 5.8 5.3 3.5 4.8 4.6 4.2 4.3 3.7 5.4 4.3 4.0 4.5 6.8 3.6 5.1 4.2 4.7
##  [307] 5.0 3.9 4.5 5.7 5.4 4.6 3.7 4.7 5.2 6.0 4.1 5.6 4.9 4.6 4.1 5.3 6.1 4.5
##  [325] 4.6 3.5 4.3 3.7 5.7 3.9 5.9 4.9 4.1 4.5 4.5 2.9 3.7 3.5 4.4 4.1 5.5 4.7
##  [343] 2.7 5.5 4.8 4.1 5.9 4.0 4.5 6.3 4.5 5.2 6.2 4.6 5.9 4.8 3.9 5.1 4.4 4.2
##  [361] 6.4 3.0 4.1 4.1 4.2 5.0 5.5 4.1 5.4 5.7 5.7 5.7 4.9 6.3 4.8 3.9 2.9 4.7
##  [379] 2.8 3.5 3.4 3.8 4.5 5.6 4.2 5.8 3.8 2.4 4.7 5.5 4.4 5.8 3.9 3.2 5.1 3.5
##  [397] 6.0 4.4 5.1 6.3 6.0 5.1 3.7 3.7 5.1 5.8 4.2 4.4 4.4 4.3 5.3 4.1 3.3 4.9
##  [415] 3.9 4.2 4.6 4.8 3.7 5.3 4.4 4.9 5.1 3.9 4.2 4.0 2.8 5.4 4.5 5.1 5.2 4.9
##  [433] 3.5 5.0 5.3 3.4 6.3 3.5 4.1 5.2 4.3 3.8 2.9 4.3 4.9 5.2 4.1 2.5 4.7 3.5
##  [451] 4.6 3.6 3.2 6.7 5.8 5.3 6.2 4.6 5.2 5.3 5.1 4.3 4.5 3.8 5.0 4.1 5.4 6.0
##  [469] 4.8 4.2 4.7 3.9 3.3 4.0 5.5 5.6 4.0 5.2 4.5 3.4 4.1 4.4 5.7 4.6 4.8 4.6
##  [487] 4.9 3.9 4.9 5.2 4.9 5.1 5.4 4.4 4.3 2.2 4.6 3.7 3.3 4.5 5.4 4.9 4.8 6.9
##  [505] 4.9 4.3 3.5 5.6 5.8 6.5 5.5 4.0 2.8 5.2 4.8 4.6 4.2 4.8 5.2 6.9 3.2 5.4
##  [523] 3.5 4.8 5.1 4.2 5.1 2.6 4.9 4.3 4.0 6.0 4.2 5.0 4.8 4.1 4.3 3.8 5.1 3.7
##  [541] 4.0 5.6 4.6 2.9 5.3 6.0 5.0 5.1 3.8 5.2 4.2 5.9 3.6 4.3 6.3 2.3 3.5 4.0
##  [559] 5.8 5.3 4.1 4.2 3.8 5.7 3.5 5.5 5.8 3.8 3.9 5.6 5.5 5.1 4.0 6.3 4.8 4.0
##  [577] 4.7 3.6 3.3 3.2 3.9 4.7 5.3 5.4 4.9 3.3 3.2 5.1 4.9 3.6 3.2 4.4 3.6 3.7
##  [595] 5.3 3.8 3.7 4.2 4.5 4.9 5.6 3.7 4.1 4.7 3.6 5.5 4.3 4.1 3.9 3.1 3.9 4.2
##  [613] 6.6 4.1 4.2 4.1 5.6 3.9 2.9 4.5 5.0 4.9 3.2 5.6 5.1 4.9 4.7 3.2 5.6 3.7
##  [631] 3.8 4.4 4.4 3.0 6.5 5.6 4.7 5.1 3.6 4.7 4.8 5.3 2.8 5.5 4.5 2.5 3.9 6.3
##  [649] 5.0 4.4 4.4 6.0 5.7 3.3 3.3 5.0 5.6 3.4 3.6 5.4 4.9 4.9 5.8 4.4 3.7 4.9
##  [667] 2.7 4.8 4.4 5.0 4.8 3.9 4.2 3.3 4.1 3.1 5.6 3.5 3.9 4.0 5.4 4.8 4.8 5.3
##  [685] 4.7 3.3 3.8 5.7 4.2 3.1 3.3 4.5 4.2 4.0 5.2 5.0 4.4 4.7 5.3 3.8 4.8 3.4
##  [703] 3.5 5.8 3.6 4.5 4.3 5.0 4.0 4.0 4.8 3.1 5.7 3.8 4.9 4.8 4.9 4.7 6.4 6.0
##  [721] 4.8 4.4 4.4 5.2 5.7 5.7 4.8 3.8 3.2 4.2 4.7 4.3 3.7 3.4 6.1 5.2 4.5 4.7
##  [739] 3.7 5.3 4.2 4.9 5.3 3.7 3.5 3.1 4.7 4.9 6.5 4.7 4.8 4.6 5.8 5.4 5.0 5.4
##  [757] 3.1 5.3 4.3 5.6 5.0 6.1 3.6 3.0 4.6 6.1 3.6 4.9 5.1 4.7 3.0 3.5 4.2 3.1
##  [775] 4.5 3.9 4.1 4.3 3.8 5.1 4.4 5.0 6.1 5.0 3.4 3.9 3.5 4.0 2.8 4.8 4.4 5.2
##  [793] 5.0 5.0 2.2 5.1 4.6 6.0 4.9 5.8 5.3 4.4 4.6 4.7 5.4 3.2 3.1 4.8 3.4 4.3
##  [811] 5.0 5.4 4.0 3.4 3.9 3.0 4.4 3.6 4.5 4.8 4.8 5.3 5.0 5.3 3.9 4.1 3.9 4.5
##  [829] 5.2 3.7 5.3 4.7 3.2 6.7 4.5 4.4 5.1 5.7 4.8 3.4 4.7 4.8 3.5 4.2 4.7 3.1
##  [847] 3.9 5.5 4.8 5.2 5.8 5.5 4.6 4.9 4.4 3.7 4.3 6.2 5.0 4.2 4.4 4.9 4.7 3.7
##  [865] 5.2 4.1 4.9 4.0 3.6 3.2 4.5 6.0 4.4 4.7 3.6 5.3 5.3 5.0 4.6 4.9 5.0 3.6
##  [883] 5.2 3.5 3.8 4.5 3.9 3.6 4.2 4.8 4.2 4.6 4.6 5.3 3.5 3.8 5.2 4.2 4.6 5.4
##  [901] 5.4 6.3 3.2 5.0 2.9 3.4 3.2 4.5 4.1 4.3 6.6 4.6 4.4 5.4 6.1 1.7 6.9 4.4
##  [919] 5.8 5.6 5.8 4.4 4.0 4.4 5.1 5.7 4.5 5.0 3.7 4.1 3.6 3.0 4.8 4.2 4.1 4.0
##  [937] 3.7 4.2 3.3 3.8 6.1 5.9 3.9 5.7 4.6 6.5 4.3 4.0 4.0 3.3 3.2 4.6 4.1 4.5
##  [955] 4.1 5.8 5.0 6.0 5.3 4.0 5.0 3.4 4.0 6.7 5.6 4.8 3.5 5.2 4.7 5.0 3.6 5.8
##  [973] 3.4 5.0 5.2 4.0 4.6 5.1 3.6 3.4 6.4 4.0 3.6 4.6 4.2 4.9 4.5 4.2 3.6 3.2
##  [991] 5.3 4.1 4.3 3.6 3.9 3.1 4.8 5.4 4.7 4.3
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
##   2.9   6.4
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
##    [1] 6.0 4.3 4.3 3.2 2.8 6.4 5.3 4.8 4.9 3.9 3.9 4.1 3.6 3.6 5.3 5.6 2.9 5.5
##   [19] 4.4 5.9 4.5 4.6 2.9 5.6 4.4 5.2 3.4 4.4 4.5 3.9 4.2 4.4 5.0 4.3 5.3 3.1
##   [37] 4.8 4.1 5.1 3.8 3.2 3.5 5.0 5.5 2.6 5.8 4.2 3.9 4.1 5.5 4.8 3.0 3.9 4.0
##   [55] 7.1 3.6 5.7 4.4 6.7 4.9 3.6 4.5 4.0 5.1 3.0 3.9 5.8 4.1 3.8 4.7 3.9 4.3
##   [73] 3.8 3.6 4.8 4.9 4.1 3.5 4.7 4.3 4.4 4.8 4.2 4.4 4.6 4.5 3.5 3.1 4.1 5.7
##   [91] 4.1 3.3 4.5 4.8 5.7 3.1 4.4 3.6 5.2 3.7 5.4 3.8 5.9 4.5 4.9 3.8 5.9 3.6
##  [109] 5.8 4.7 3.6 4.5 4.5 5.4 5.6 4.1 3.4 4.6 5.8 4.8 4.9 4.7 5.2 5.7 3.7 4.6
##  [127] 2.2 5.8 4.2 3.5 4.1 6.3 6.5 5.2 3.4 5.0 4.8 4.9 3.5 6.0 3.4 2.8 5.2 5.4
##  [145] 4.5 3.6 3.9 5.3 6.5 4.8 3.5 5.0 2.9 2.8 4.4 3.5 3.9 4.9 5.5 5.3 4.5 5.4
##  [163] 3.3 5.8 6.2 4.0 4.7 4.1 3.5 3.0 4.8 4.8 4.4 4.0 3.7 4.2 4.1 4.5 5.5 3.8
##  [181] 3.6 4.5 6.1 4.4 6.0 3.5 4.8 4.1 3.1 5.1 4.7 3.0 3.9 4.2 5.7 5.7 5.5 5.7
##  [199] 3.9 3.8 3.5 4.0 3.7 4.2 3.7 4.4 5.1 5.6 4.2 3.6 4.6 5.3 3.6 5.0 4.9 5.0
##  [217] 3.5 5.9 5.4 3.9 5.3 2.8 4.1 4.7 3.4 6.0 5.2 3.0 4.0 5.0 3.7 5.2 4.3 4.2
##  [235] 4.3 3.7 6.4 4.7 4.9 3.9 3.7 5.9 2.8 4.6 5.4 3.8 6.1 5.2 4.7 5.7 4.3 2.8
##  [253] 4.8 3.7 4.6 3.8 5.3 5.2 4.6 5.7 3.9 4.4 5.3 4.5 3.2 5.9 4.0 5.6 5.2 4.3
##  [271] 2.7 5.3 3.8 4.5 4.4 5.2 4.0 5.3 4.9 2.6 5.2 4.5 4.5 5.7 5.3 4.8 3.9 5.8
##  [289] 4.0 4.3 3.5 4.8 4.3 3.4 4.4 4.1 5.8 4.8 3.7 3.8 4.2 4.1 4.2 5.5 4.8 4.7
##  [307] 2.9 3.8 5.0 5.1 3.6 4.8 3.8 4.9 4.6 3.1 5.5 4.0 4.4 5.6 5.5 3.5 4.9 4.7
##  [325] 3.8 3.5 5.5 4.4 4.0 4.6 3.4 4.2 3.7 4.1 5.2 3.4 4.5 3.9 3.8 3.9 4.6 3.8
##  [343] 3.1 2.7 3.0 4.3 3.4 3.9 3.4 4.9 3.5 6.5 3.6 4.5 4.0 3.7 3.0 4.3 3.3 4.0
##  [361] 5.4 5.5 2.3 3.8 5.6 2.9 4.8 4.1 6.4 4.9 6.0 4.0 4.7 4.2 4.5 5.0 5.0 4.1
##  [379] 5.5 3.2 3.2 5.4 4.0 3.1 3.9 4.4 4.3 2.6 5.6 4.7 4.5 3.5 4.2 5.0 3.7 7.2
##  [397] 3.5 4.6 5.4 3.8 6.1 4.9 5.6 4.3 3.8 5.7 4.2 3.0 5.4 3.2 4.2 5.2 5.0 4.4
##  [415] 5.5 3.0 3.8 5.9 4.1 5.0 4.3 5.2 4.3 1.7 4.2 5.4 4.5 4.4 4.4 3.3 5.0 5.8
##  [433] 5.1 4.1 4.1 5.4 5.3 6.3 4.2 5.5 3.6 3.4 4.2 4.3 4.3 5.4 3.9 4.4 3.7 4.3
##  [451] 3.7 2.8 4.4 6.0 3.6 3.4 3.8 5.3 3.8 5.5 5.6 5.0 4.4 4.2 4.8 5.2 4.9 4.4
##  [469] 4.9 6.2 4.8 3.7 4.5 4.5 4.8 6.0 4.2 4.0 5.0 3.7 4.5 3.6 6.2 2.7 3.1 4.9
##  [487] 4.7 3.7 4.5 3.0 3.6 4.0 3.0 3.8 3.9 5.2 4.3 5.0 3.9 4.3 3.0 5.1 6.2 4.0
##  [505] 4.3 5.2 4.3 4.6 5.0 5.4 3.4 5.1 4.6 6.3 4.0 6.1 4.3 2.9 5.0 5.0 3.9 3.3
##  [523] 4.6 5.1 2.8 5.2 3.3 4.7 5.1 4.5 4.2 4.5 5.0 5.4 4.0 4.3 4.5 6.6 4.6 3.4
##  [541] 3.0 2.8 4.6 3.3 4.0 3.8 6.8 4.8 5.6 5.6 3.4 5.4 5.9 4.7 5.5 5.2 4.4 5.1
##  [559] 7.5 3.6 5.1 4.2 4.1 4.2 4.9 3.3 6.0 3.9 4.9 3.3 5.1 6.2 4.4 4.2 4.8 4.6
##  [577] 4.4 5.0 3.6 4.9 3.6 4.6 4.9 4.7 2.7 5.0 5.3 5.0 4.9 3.8 5.1 3.4 4.8 3.9
##  [595] 4.7 2.5 3.4 4.5 4.0 4.7 3.0 4.3 4.2 5.7 3.7 6.3 5.1 4.0 5.3 6.1 4.3 2.6
##  [613] 3.7 5.2 4.9 4.3 5.1 4.3 5.4 5.1 3.6 4.9 4.9 4.8 4.3 5.3 3.9 4.8 5.3 4.6
##  [631] 5.4 3.7 5.2 4.2 4.3 3.5 3.9 5.9 5.2 2.2 5.4 4.9 3.8 4.0 3.1 3.6 2.9 4.6
##  [649] 4.6 3.6 5.4 6.2 4.4 5.2 4.2 3.8 4.7 3.9 5.4 3.2 6.5 2.3 5.2 5.0 3.2 4.3
##  [667] 5.7 5.3 6.4 4.7 4.4 2.2 3.4 4.3 3.0 7.1 5.4 4.4 4.7 5.6 4.8 4.3 4.3 5.0
##  [685] 4.2 4.7 4.1 5.6 4.6 4.2 6.1 4.3 5.1 3.9 5.4 5.3 4.7 1.6 5.7 5.0 4.8 3.5
##  [703] 3.9 4.3 3.4 4.0 3.3 4.7 4.6 5.8 4.0 4.2 2.6 6.0 4.4 5.8 3.5 3.5 5.6 5.0
##  [721] 4.9 4.7 4.5 4.2 5.5 4.9 3.6 5.6 3.3 3.3 2.7 6.6 4.5 4.9 4.4 4.0 4.8 4.4
##  [739] 4.1 4.2 3.5 3.6 5.3 3.7 4.5 3.4 4.3 3.8 2.7 3.7 4.7 4.2 3.6 3.5 2.7 4.2
##  [757] 4.0 4.9 4.1 4.5 4.6 4.6 5.0 6.1 3.6 4.0 4.0 4.8 4.0 5.1 4.3 4.2 4.7 3.9
##  [775] 5.2 4.3 4.9 3.3 6.0 2.5 5.1 5.2 4.3 4.2 3.7 4.0 4.5 5.2 4.6 4.1 5.3 3.9
##  [793] 4.6 3.6 2.9 5.0 5.6 4.8 4.2 5.7 4.7 5.5 1.9 4.5 4.7 3.5 3.5 4.7 4.9 5.1
##  [811] 3.5 3.8 4.1 3.0 5.2 6.2 4.9 5.1 4.5 3.9 4.0 4.5 4.0 4.1 5.2 3.4 5.6 4.2
##  [829] 5.8 4.1 4.0 4.7 2.2 4.7 3.5 4.0 6.0 4.5 3.4 4.2 5.6 5.3 4.2 5.9 3.8 5.1
##  [847] 4.1 3.2 6.4 5.2 3.6 3.2 4.7 4.9 4.0 5.7 4.1 3.1 4.9 3.2 4.2 3.3 4.9 5.3
##  [865] 4.2 5.3 3.5 3.9 4.6 4.1 6.0 5.9 4.1 4.1 5.0 4.4 4.4 4.7 5.0 4.6 4.5 4.5
##  [883] 4.7 4.9 4.1 4.2 6.3 4.7 4.4 4.5 5.4 4.0 4.0 5.7 4.0 5.4 4.3 4.0 5.1 3.4
##  [901] 4.7 4.8 3.8 5.0 5.4 4.4 4.7 4.9 4.6 4.9 2.7 5.0 4.6 4.8 6.0 3.3 4.3 4.0
##  [919] 4.0 4.2 4.6 5.6 3.5 4.4 3.7 3.5 3.9 3.5 5.1 3.2 4.3 5.4 4.9 4.9 3.4 4.4
##  [937] 5.2 4.8 4.2 4.3 3.1 5.3 3.6 3.9 4.2 6.2 5.9 3.8 3.3 4.4 4.5 4.0 4.2 3.8
##  [955] 4.9 4.2 5.7 4.9 3.7 4.1 3.5 5.4 3.6 2.6 3.3 4.5 3.5 3.2 4.4 5.0 5.1 3.4
##  [973] 4.5 3.8 5.6 5.9 4.5 5.1 3.7 4.8 5.0 4.1 5.0 4.6 4.2 5.5 4.8 4.1 5.7 4.1
##  [991] 4.7 5.3 4.4 4.7 5.3 5.6 5.2 3.7 5.2 4.6
## 
## $func.thetastar
## [1] -0.043
## 
## $jack.boot.val
##  [1]  0.48348348  0.31759531  0.22418879  0.12429379 -0.03554217 -0.09357542
##  [7] -0.20027473 -0.34894260 -0.37292225 -0.54540390
## 
## $jack.boot.se
## [1] 0.9433413
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
##    [1] 3.1 6.1 4.7 3.4 6.5 4.6 2.4 2.9 3.2 4.7 4.8 3.9 4.0 4.3 5.4 4.3 3.3 4.6
##   [19] 4.1 6.5 4.3 5.8 4.7 3.6 6.6 3.9 4.5 4.2 5.3 6.0 6.4 5.8 4.4 3.6 2.9 4.8
##   [37] 3.5 5.0 5.5 4.9 3.7 4.0 4.7 5.6 3.0 4.3 5.0 4.2 6.7 4.6 2.2 3.6 4.7 6.0
##   [55] 2.7 4.4 4.2 4.7 3.5 4.2 3.4 6.7 4.3 5.2 4.4 3.5 4.6 3.4 3.0 3.3 6.0 4.7
##   [73] 5.6 6.2 5.3 4.0 5.9 4.9 4.0 5.6 5.6 3.9 3.3 5.0 4.4 6.3 5.4 4.7 5.2 3.4
##   [91] 4.9 2.4 3.1 5.1 4.8 4.6 4.7 5.5 5.1 5.5 2.8 3.9 5.1 4.5 5.5 4.7 6.0 6.1
##  [109] 5.3 4.7 4.9 4.8 6.3 5.9 3.8 3.9 4.5 4.8 3.6 5.1 4.0 4.7 2.7 4.5 5.0 5.0
##  [127] 2.3 4.6 5.0 4.1 4.4 4.1 3.3 4.4 4.3 4.3 5.1 3.4 4.8 5.3 4.0 5.2 3.8 5.1
##  [145] 3.0 3.4 2.5 4.9 3.7 3.4 4.2 3.8 3.5 4.3 4.0 3.1 6.4 5.2 2.6 4.5 3.7 3.4
##  [163] 5.2 6.0 5.6 3.6 6.3 4.6 4.8 4.3 3.8 4.6 5.9 3.5 4.5 5.3 3.4 3.5 4.0 5.1
##  [181] 4.9 3.9 5.8 4.3 4.9 4.2 4.1 4.5 2.8 4.5 3.5 5.7 4.7 3.0 2.6 3.3 2.3 4.6
##  [199] 4.9 5.4 4.6 3.1 4.9 4.7 3.3 4.6 3.6 4.8 4.0 4.6 5.1 2.8 3.1 3.2 3.8 5.3
##  [217] 5.0 5.4 7.8 3.7 4.3 6.3 3.3 3.6 5.6 4.7 3.5 5.9 5.7 6.3 4.1 5.6 2.9 5.1
##  [235] 3.9 4.1 6.7 4.1 4.8 4.3 4.6 4.8 3.7 4.3 4.9 3.9 5.4 3.8 3.4 3.5 4.9 6.6
##  [253] 5.4 5.3 3.2 2.6 4.0 4.9 6.2 4.1 3.8 4.6 4.4 4.8 4.4 5.7 4.3 6.0 5.1 3.7
##  [271] 4.8 5.9 4.0 6.0 5.3 5.0 3.6 3.5 3.7 5.1 4.3 5.4 6.1 5.3 4.1 4.2 4.1 6.2
##  [289] 3.4 4.1 3.8 4.5 4.1 4.6 4.1 2.5 5.1 4.3 3.9 5.0 4.7 4.0 5.9 4.6 3.5 4.9
##  [307] 3.9 2.7 2.0 4.5 5.9 4.3 3.8 4.7 4.2 4.4 4.2 4.3 2.2 5.3 4.9 4.2 4.2 3.3
##  [325] 4.1 4.0 4.7 4.6 5.3 4.9 2.7 3.9 5.1 3.0 4.6 4.9 3.2 3.1 5.7 3.1 6.1 4.5
##  [343] 4.0 5.7 4.6 4.5 3.5 3.6 3.7 4.7 5.1 4.0 4.3 3.6 4.3 4.7 5.1 3.1 4.0 4.6
##  [361] 5.0 5.2 4.4 3.9 5.4 3.4 5.4 3.2 2.5 4.7 4.7 4.8 4.9 4.2 5.0 4.6 4.9 2.1
##  [379] 6.3 4.6 3.8 4.9 3.3 3.8 4.2 3.3 4.7 5.1 3.6 5.6 3.6 4.4 4.3 3.3 4.4 6.8
##  [397] 3.1 5.2 4.1 3.9 4.1 2.6 3.4 6.1 4.4 4.3 5.9 2.9 5.3 6.3 5.9 4.7 4.8 4.3
##  [415] 6.2 3.3 4.9 3.3 3.1 3.8 4.4 3.5 3.4 3.2 3.0 6.2 3.1 3.8 4.5 5.9 4.5 4.9
##  [433] 5.0 4.8 4.9 5.5 5.4 2.5 5.1 3.9 5.9 6.4 2.9 4.4 5.9 4.7 3.4 4.2 4.0 4.5
##  [451] 5.1 5.5 2.5 3.6 4.4 5.5 2.7 4.4 3.0 3.6 4.2 4.6 2.7 2.8 4.9 4.1 2.9 4.3
##  [469] 3.6 4.1 4.1 5.3 4.5 4.8 2.9 4.4 3.9 3.9 5.3 6.1 4.6 5.0 3.6 3.3 4.6 3.9
##  [487] 3.6 4.4 6.1 4.7 4.7 4.6 4.5 3.4 2.6 5.9 4.7 4.7 4.8 4.4 3.9 3.0 5.1 4.3
##  [505] 4.8 4.2 4.3 5.4 4.6 5.1 5.0 4.7 4.3 4.1 4.6 4.5 5.7 5.7 4.7 4.8 5.7 5.2
##  [523] 4.1 4.3 4.2 2.8 3.0 5.9 2.4 3.6 3.4 5.1 3.8 4.2 6.1 4.9 4.3 4.0 3.6 4.1
##  [541] 6.0 2.8 5.1 4.8 4.8 5.5 4.4 4.4 3.2 5.3 5.1 5.0 4.1 5.5 5.2 3.8 4.6 4.2
##  [559] 4.7 3.5 4.7 5.9 5.2 5.1 4.7 4.1 3.8 6.0 4.7 4.1 4.1 4.6 3.4 4.6 5.6 4.9
##  [577] 4.6 4.8 4.2 4.3 4.8 2.9 5.2 3.5 3.4 3.7 5.2 5.2 2.8 4.3 5.5 4.9 3.8 4.5
##  [595] 5.8 5.6 4.1 3.2 5.3 3.8 4.7 4.5 5.0 4.6 5.8 4.9 5.1 3.8 3.8 3.5 4.0 4.7
##  [613] 6.3 3.7 4.5 5.2 5.5 4.1 3.6 4.0 4.9 5.3 4.6 5.0 5.2 5.0 4.8 4.4 5.5 2.7
##  [631] 5.6 5.0 4.1 4.5 4.7 3.8 5.5 3.8 3.8 4.7 4.2 4.7 3.5 3.7 3.2 2.9 3.7 4.3
##  [649] 4.1 4.2 3.4 5.2 5.1 6.0 4.9 4.1 3.8 5.7 5.8 3.6 5.0 5.2 3.0 4.9 5.6 4.9
##  [667] 4.7 5.1 4.4 3.3 6.7 4.4 4.4 4.8 3.2 5.0 4.7 5.6 3.5 2.8 5.5 3.5 3.4 2.6
##  [685] 3.1 3.6 3.3 4.1 4.9 4.4 3.8 4.1 3.2 3.6 4.0 4.4 4.1 5.8 4.1 3.9 3.8 4.5
##  [703] 4.1 5.0 4.6 4.5 4.1 3.8 5.8 5.0 4.2 2.3 3.4 4.7 3.9 4.9 5.2 5.0 4.7 4.1
##  [721] 4.2 6.1 5.6 4.0 4.6 5.1 5.2 5.1 4.5 5.4 4.7 4.6 4.9 4.6 4.4 3.6 2.7 4.4
##  [739] 5.0 4.4 4.0 4.3 4.3 4.1 3.3 5.2 2.8 4.7 2.6 4.9 4.9 4.4 4.5 5.6 3.4 5.5
##  [757] 4.9 4.7 4.5 5.1 3.9 4.0 4.1 3.6 5.0 4.4 3.7 3.8 5.6 6.2 3.0 4.5 5.7 4.2
##  [775] 6.3 4.0 6.4 6.1 3.4 4.2 6.5 5.3 2.3 4.3 3.8 5.5 4.6 5.9 3.5 4.4 4.4 3.8
##  [793] 4.1 5.4 4.7 3.9 4.7 5.6 2.3 4.4 5.1 4.3 4.9 3.6 5.9 6.0 4.9 5.2 5.1 5.3
##  [811] 3.5 5.2 5.4 4.6 4.8 5.5 4.8 3.3 5.2 3.7 5.1 5.5 3.7 4.4 5.8 5.4 3.4 6.1
##  [829] 3.5 6.2 4.7 4.3 5.0 4.1 5.0 4.7 3.7 2.8 4.9 3.9 3.7 5.7 4.3 3.8 4.8 4.8
##  [847] 4.4 6.8 4.7 5.4 4.0 4.0 3.3 3.8 4.8 4.7 2.8 3.6 4.0 3.8 2.3 5.5 4.0 4.6
##  [865] 4.8 5.4 3.3 3.8 4.1 2.9 4.7 6.3 5.1 5.5 3.4 6.5 6.6 6.1 4.8 4.3 5.5 4.3
##  [883] 6.1 4.2 6.0 4.5 5.9 4.7 4.4 3.5 3.6 3.9 5.6 5.1 4.6 5.9 4.8 3.3 4.3 5.5
##  [901] 2.5 4.3 4.3 4.3 4.1 4.2 4.0 5.3 4.3 5.2 4.6 4.0 6.0 4.1 6.3 3.0 4.3 4.3
##  [919] 4.2 3.5 3.1 6.6 4.3 5.3 5.9 5.2 4.5 4.3 3.0 5.1 4.7 5.4 4.0 5.4 4.0 4.6
##  [937] 4.1 4.4 5.3 2.9 5.9 3.9 3.6 4.8 6.4 3.5 4.0 4.9 6.0 4.2 3.9 3.7 4.4 3.2
##  [955] 4.8 6.3 6.3 4.8 5.9 5.2 4.7 3.9 6.2 4.1 4.4 2.9 5.5 4.6 5.1 6.3 3.8 5.8
##  [973] 7.1 6.0 5.5 5.6 6.1 3.7 3.6 5.5 4.1 5.1 4.7 4.2 5.9 4.9 4.2 4.8 5.6 4.1
##  [991] 4.9 3.8 5.1 4.8 4.5 4.6 4.1 3.4 4.0 4.6
## 
## $func.thetastar
## 72% 
##   5 
## 
## $jack.boot.val
##  [1] 5.500 5.500 5.300 5.400 5.100 5.000 4.900 4.756 4.500 4.500
## 
## $jack.boot.se
## [1] 1.085296
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
## [1] 0.6932996
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
##   2.440150   3.807852 
##  (1.025458) (1.776215)
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
## [1]  1.08538198 -0.05791193  0.72483073  0.72805707  0.72979404  0.50595612
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
##    [1]  0.9877085872  1.0388539168  0.1496696365  0.8147133362  0.6259092634
##    [6]  1.0386918991  1.0602018919  0.4564303181  0.9198941507  0.8618109301
##   [11]  0.7031561894  0.6633235635  0.3462557470  0.7271484892  2.3510859128
##   [16]  1.0467459350  0.6790377075  0.7193508605 -0.0436049005  0.2610077200
##   [21]  0.8249795828  0.5416498730  1.1382781139  0.8482345555  0.8081321295
##   [26] -0.3823539149  0.3817410905 -0.3719730387  0.1056533746  0.4515681182
##   [31]  0.7467058435 -0.5932022818  0.5384799077  0.8660524561  0.5376690626
##   [36]  0.4817602682 -0.2499867287  1.0814680605  0.0099377982  0.4904473038
##   [41]  1.4237509414  1.4641729812  0.9074865126  0.8885713687  0.5424047966
##   [46] -0.2641586910 -0.0400980913 -1.0279350192  0.0554249402  1.4329917943
##   [51]  0.2816617621  0.7278863872 -0.5943890725  1.2145715864  0.4736594692
##   [56]  1.1256738176  1.1514644024  0.1134690063  1.4641729812  0.4899945420
##   [61] -0.1645952217  2.4838471130  0.2367670033 -0.3405747808 -0.0340412213
##   [66]  0.2535904317  0.2366889230  0.3640531706  2.4490794882  0.6788158623
##   [71]  0.9656598589 -0.1875995327  0.7219965971 -0.6277608795 -0.1427556399
##   [76] -0.0424949166  0.5445614729  0.3516283234  0.3251395805  1.5661524063
##   [81]  1.1830478842 -0.0196263919  0.6712359511  0.1096825046  0.0747365787
##   [86] -0.2230523471 -0.5360281883  0.3964793219  0.2346366269  0.6867600878
##   [91]  1.1228821886  0.6055004644  0.2774045673  0.1220441527  0.5706560632
##   [96] -0.4652714022  1.1960242120 -0.4885646536  0.8877802883  0.7234451737
##  [101]  0.1984297608  0.7532053583 -0.1917241403  2.6341043421  0.6263635292
##  [106]  0.4898498378  0.7456879871  0.8080395348  0.7213323389 -0.1013822783
##  [111]  0.0766621526  1.5626301048  0.3008700048 -0.0672957075  0.6710298832
##  [116]  1.3271222239  0.3764212524  0.4479090775  0.7453134386  0.9532007829
##  [121]  0.5073736325  0.9544551072 -0.0003798039  0.6173883906  1.1103201748
##  [126] -0.6548132882  1.0592980446 -0.3766954520  0.6204462292  0.7161498653
##  [131]  0.2705574069  0.5036247469 -0.4811975170  0.1611068491  1.6513110027
##  [136]  1.0888439661 -0.2939950396  0.5682966401 -0.6486627892  0.8902385120
##  [141]  1.1202825373  0.4621836280 -2.6192475197  1.3892343074  0.6054864609
##  [146]  1.4997720348 -0.0102225171 -0.4230785356  0.0952519010  0.4688614510
##  [151]  0.9377612011  0.5133059127  0.6025385335  0.5551546127 -0.1186373750
##  [156]  1.4408860150  0.3205277005  0.2287457790  0.2165761541  0.5796335485
##  [161]  0.8951115622 -0.0217576002 -0.3924750389  0.4381108899  0.7059208821
##  [166]  1.8700139117  0.4784477574  1.3723855789  0.8477667400  0.6166400417
##  [171]  0.5444126390  0.4496477976  0.3849870630 -0.0947317881  0.9427527267
##  [176]  0.4674929463  0.5098835139  0.9705497393  1.5777262792  0.5318922203
##  [181]  0.4340078753  0.3156500494  1.2483675086  0.3368758828  0.8728962750
##  [186]  0.8614544341 -0.9221324661  0.3324663595 -0.3162037393  1.3731666921
##  [191]  0.6430560428  0.2207203362  0.3398391064  1.1900873874  0.5866443981
##  [196] -0.3924500579  0.4888244097  0.0823544903  0.8472402918  0.9497965186
##  [201]  0.5103475602  1.4275235538 -0.3533763950 -0.7957919442  0.5820569027
##  [206]  0.5502664415  0.4267532244  1.3652909977  0.6222313897  0.7095669525
##  [211] -0.2005868444  0.9918610431  0.4398872668  0.2969432671 -0.5579624413
##  [216]  0.5083057587  0.8636793049  0.7471202588  0.3393652076  0.5475524076
##  [221]  0.3342104818  0.1955849530  0.0301517918  0.2109798145  0.3439442348
##  [226] -0.3074815138  0.6987822303  1.4576362835 -0.2867789497  1.4871093321
##  [231]  0.1948785074 -0.3580807682  1.3771499632  1.4232884958  0.5124272430
##  [236]  2.5648567287  0.4993209526  0.1116608589  1.2570516328  0.7287582572
##  [241]  0.5679863674  0.6697633811  0.3759154790  0.1766186688  0.0831319455
##  [246]  1.4254144440  0.0289001006  0.2573923263  0.8952501408 -0.1974631281
##  [251]  0.3602430528  0.3698049204  0.8707144801  0.8071753704  1.4662461066
##  [256]  0.7183665487  0.1133794437  0.2316127119  0.8981909430  1.2688725203
##  [261]  0.7071523838  0.7430030495  1.5225608678  0.4640149280  1.1251859998
##  [266]  0.0075778928  0.7020250355  0.1620159681  1.5777310103  0.7279773061
##  [271]  1.3706774310 -0.1995181722 -0.1273637626  1.0887527900  0.4910913669
##  [276] -0.2441980894  0.8781504765 -0.1739602826  0.4308906514  1.4128421670
##  [281]  0.1845284084  0.5197367889  0.1265194027  1.8702973120 -0.3491476432
##  [286]  1.0730953249  0.8193969071  0.8547653675  1.7010282662  0.8246767482
##  [291] -0.1320544499  0.9412562887  0.3733522035  1.6367047192  0.3041179657
##  [296]  0.2906922818  0.1231234417  0.1603578743  0.4364349317  1.3892934630
##  [301]  1.4466975726 -0.3844543173  0.3781062183  0.4910913669  0.3205004733
##  [306]  0.6321438796  0.0188215572  0.3318394820  1.2740593281  0.4340078753
##  [311]  0.1544667207  0.6425791016  0.4310600919  0.5421328575  0.0321199458
##  [316] -0.2614657624  0.1648206823  0.0351964070  1.0689626536  0.1434946380
##  [321]  0.8665971398  0.2580315486  0.4399823889 -0.6128074777  0.7880037270
##  [326]  1.0834207938  0.6506871464  0.0748535805  1.9364610666 -0.0575783805
##  [331]  1.1364747715  0.1879210363  0.8755729174  1.4943798715  0.1619975357
##  [336]  0.8762151484 -0.1015846446  0.7849961705  1.4779265291  0.8549733608
##  [341]  0.3168401191  0.0849782143  0.5589017907  0.4315094763  0.4344105341
##  [346]  0.8705904460  0.5472357842  0.8530499289  0.5532709532  0.7008511378
##  [351]  0.8668143161  0.0564144682  0.1486238653  1.2045330272  0.2186066321
##  [356]  1.5657912163  0.8449755137  0.6693985332  0.3907115317  0.0873123150
##  [361]  0.5872481638  0.0025428792 -0.7371313565  1.3924282867  0.7183571353
##  [366]  1.0810741968  0.4171789560  0.6506235073  0.4157366585  0.4170616756
##  [371]  0.2316127119 -0.0493268848  0.9904571205  0.3166021032  0.2696586672
##  [376]  0.7261797982  0.3060958327  0.5844353388  0.3441941194  1.1767629013
##  [381] -0.1308570754  0.3794133138  0.5363976172  0.8547825416  0.8584811097
##  [386]  1.2140141149  0.0190570517  0.8449755137  0.4454023142  0.2394019597
##  [391]  0.6026673466  0.8360387044  0.0261350979  0.6709700421  0.2773794504
##  [396] -0.0226954559  0.1827404652  0.9890769329  0.5160488045  0.4394620983
##  [401]  0.3289514675  0.1623385108  1.4671407199  0.3827610154  0.3923701400
##  [406]  0.1853207565  0.8108850285  1.2446568027  0.8269753449  0.8826948483
##  [411]  0.2950991004  0.4470325377  0.2824582842  0.2332919738  0.9930226088
##  [416]  0.3059828128  1.3921987758  0.2120413830  0.8418265445  0.0964840186
##  [421] -0.2666762712 -0.2066421664  0.2042874239  0.3004269365  0.6441295797
##  [426]  1.5865414691  1.4681672205  1.3088653858  0.2839305238 -0.0722388962
##  [431]  0.8254085718 -0.1086459117  0.8202751982  1.5869081317  0.9799591259
##  [436]  0.4792739275  0.1127128031  0.5031953393  0.6143584647  1.2051956517
##  [441]  0.4467301800  0.1107299895 -0.7799992401  0.2107192394  0.4147026388
##  [446]  1.4901282082  0.1867043856  0.4597919803  0.9498666440  1.1641938019
##  [451]  0.2960505172  0.3443387106  0.5532524758 -0.6649047839  0.0819838700
##  [456] -0.3427063078  0.5283715362  0.6153935668  2.4797043153 -0.6864301091
##  [461]  0.4154460928  0.2786103434  1.1318668605  0.8099846274  0.6801891277
##  [466]  0.2845186083  0.5662027835  1.4004246485  1.2655291013  0.3690782032
##  [471] -0.4535544923 -0.2317559464  0.7109422249  0.0501075335  0.2416430518
##  [476]  0.5070143202 -0.2868716758  1.3524968012  0.8952893381  2.1630022645
##  [481] -0.3418985005  0.3821516491  0.4840488540  1.4638050643  0.5829976202
##  [486]  1.3793752888  0.5888384565  0.0930238271  1.5778958714  1.1529906956
##  [491]  0.0797736560  0.6003145210  0.5752424275  0.1109883114 -0.0269796365
##  [496]  0.4105391082  2.3029345051  0.6239447710  0.1741066039  0.5707037906
##  [501] -0.1510527561  1.5526022613  0.0086222855  0.4030490026  1.4866653739
##  [506]  0.3860447350  0.8039371049 -0.1253501169  0.9177513962  0.6835298115
##  [511] -1.0636279641  0.0314285530  0.3401833535  0.1416285524  0.4621836280
##  [516]  1.5145283520  0.8876605680  0.8167816971 -0.3093424412  0.0545403275
##  [521]  0.3234806730  0.8714275331  0.6971530426 -0.1452659529  0.8433515559
##  [526]  0.4891717666  1.1530012934  0.2386387841  0.6848008212  0.6400924966
##  [531] -0.3012325094  0.3379184003  0.5489686184  0.9420147648  1.2350468696
##  [536]  1.4018547571  0.3390301349  0.4899837817  0.4587237747  0.8575413515
##  [541]  1.0030610121  0.4122343950  0.9009840101  0.1017577182 -0.2545457044
##  [546]  0.6932736731 -0.1281583403  0.5904849164 -0.0415304181  1.3460871402
##  [551]  0.2431331373  0.5175562005  0.4569560016  0.8609290430 -0.0185273064
##  [556]  1.4092962134  0.1349239487  0.8014003136  0.4257684673  0.1515300347
##  [561]  0.8933174823 -0.0559073424  0.4922750319  0.5167747796  0.5622392796
##  [566] -0.0167943676  1.1881328482  1.0369757080  0.6978287397  0.5930615872
##  [571]  0.8747428735  1.6127460362  0.6550453105  0.0708840635  1.2325121065
##  [576]  0.8859311516  0.1445909331  1.3831357011  0.3884227302  1.1917622148
##  [581]  0.1004637490  0.1341541719  0.7852799985  1.0162626702  0.8301294456
##  [586]  0.6113395343  0.4339591119  0.0679719585  0.3240030194  0.0379218228
##  [591] -0.0026298242  0.0290865546  0.5021996920  0.8648634403  0.8181776796
##  [596]  0.2577817920  0.5927581672 -0.5297788488 -0.2331233824  0.9787147927
##  [601]  0.5512792925  0.8283484314  1.8422627303  0.8294737237  0.8520263079
##  [606]  0.0844851740  0.2469006969 -0.1509072570  1.3358464793  1.6270139750
##  [611]  0.2328674424  0.9136516627  0.8378882018  0.4897126959 -0.4113908578
##  [616]  0.0014218951  2.0419833569  0.2508008360  1.2647222899  1.1773342817
##  [621] -0.0807019308  0.6416047129  1.4120051929  0.1755095283  0.4578814791
##  [626]  0.7359704319 -0.1164875016  0.0202277594  0.0116924589 -0.3653872438
##  [631]  0.8964181434  0.0092168398  0.4280732989  0.3135790652  0.0988954174
##  [636]  0.1537492182  1.1387895288  2.2992531544  1.1530012934  1.2365712287
##  [641]  0.4963989940  0.9340056609  1.1246911803  0.0598652637  0.8930577181
##  [646]  0.4744779214  0.0065605675  0.8852088271  0.0809206803  0.6828350258
##  [651]  0.7424879514  0.3830643287  0.0822196976  2.1240665400  1.4406848325
##  [656]  1.0676242415  0.9312575238  0.8591016882  0.8081835684  0.4733943700
##  [661]  0.1592539311  1.2044523002  0.9890340565  0.3837647608  0.8106171929
##  [666]  0.6987440129  0.6880145746  0.1826206299  0.5318249745 -0.2330578268
##  [671]  0.8883134306  0.3214159566  0.4400060775  0.2296688206  1.0530768986
##  [676]  0.3622617291  0.9231336526  0.6134296256  0.3726958425  0.8985693054
##  [681]  1.3504756649  0.8807733707  0.7440440758  1.2157430463  2.1941814421
##  [686] -0.1149421885  0.6105552854  0.4192425775  0.7223311957  0.4613818103
##  [691]  0.8709036243  0.1446305565  0.2788091495  1.1788244537 -0.2199883350
##  [696]  0.2610830687  0.1985376822  0.3975056436  1.2795004381  0.5889758870
##  [701]  0.1879210363  0.7979391848  0.8092091249  1.7791554821  0.2615016672
##  [706]  0.3355307166  0.8184009003  0.7062126425  0.0948133275  0.2078257364
##  [711]  0.5475524076  0.9236639307  0.7910591224 -0.6450784020  1.0462214805
##  [716]  1.4089557336  0.3380014464  0.9185188474  0.0859819461  0.3245570347
##  [721]  0.4792267208  0.3317348854  0.2193160399  0.5791308211  0.3935461575
##  [726]  1.4523788292  0.3239856106  0.4696630690  0.8869451859  0.7367505731
##  [731]  0.2292805545  0.4561603332  1.1847920162  1.1512753844  0.5210158011
##  [736]  0.8697612890  1.4446719145  0.7464087449  0.8207532540  0.0473085897
##  [741]  0.4598176114  1.0213751504  0.2863594653  0.2045809937  0.8281935681
##  [746] -0.1915981868  0.4822279821  0.9052665486  0.1824263397  0.4154335755
##  [751]  0.7893441461  0.1938872406  0.5070401780  0.6913666714  1.3895565976
##  [756] -0.1659919005  0.2146082976  0.0647987197  0.5283689422  0.4306638168
##  [761]  0.8707207455  0.9564084328 -1.0158480083  0.6139998473  0.8248988468
##  [766]  1.2233842201  0.5294973319  0.1699401965  0.7910591224  1.2166085092
##  [771]  0.3937552361  0.6138771383  0.4245537232  0.4978477403  0.1658803353
##  [776]  0.3735706065  0.9777670670  0.2665511889  0.7476847673  0.6058598342
##  [781] -0.0189150609  0.0260697082  1.4681672205  1.1016626307  0.4793275539
##  [786] -0.1315414276  0.1396008890 -0.5328704999  0.1328486169  0.5270728537
##  [791]  0.7448919444  0.5431796940 -0.6172332763  0.5954218944 -0.0971821212
##  [796]  1.2867565074  0.8053301054  0.0850378016  0.4524107405  0.6811156017
##  [801]  0.3363378044  0.8178270496  0.7168719332  0.5847838521  0.3843025529
##  [806]  0.5105699293  0.5949408266  0.7805697722  0.4445866438 -0.3707670751
##  [811]  0.4245383971  0.6158282682  0.0776306884  0.9849041010  0.3803666693
##  [816]  0.7059980417  0.7065869089  0.2280493816  0.1974625015 -0.4046504074
##  [821]  0.3928307679 -0.4474717079  0.4017662865  0.7343989521  0.5925030785
##  [826]  0.5754868985  0.4710580042  0.4906767015  0.7894008915  1.1001130351
##  [831]  0.1050369225 -0.0740593566  1.5876268140  2.5871280662 -0.5154420792
##  [836] -0.4079200190 -0.2956184172  0.8433678792 -0.6124651228  0.0751367116
##  [841]  1.5003235156  0.6962239373  0.5744910230 -0.2370844078  0.3663801144
##  [846]  0.4869163450  1.3840484442  0.3982890141  0.9401375226  0.8967377260
##  [851]  1.4679585571  0.4145863144  0.2277590696  0.1892160681  0.6697989331
##  [856]  0.0170755383  0.2287457790 -0.4582721865  0.7293329890  0.0699528164
##  [861]  0.7448493748  0.3541133915  0.4191385599 -0.3459925478  0.2815469662
##  [866]  1.9124252698  1.2603302335  0.5489155374 -0.1186036785  1.0001719601
##  [871] -0.4106565215  0.6835978470  1.3727854594  0.5189490989  0.0299926517
##  [876]  0.7253720358  0.2383673320 -0.3338481328  0.6105552854  0.1306767733
##  [881]  0.9420862318  0.4081993843  0.2041869355  0.1336279792  0.0761162115
##  [886] -0.0833311771  0.0656892910  0.7981690345  2.2913868899 -0.2218008787
##  [891] -0.1215352081  0.6160904927  1.0263473234  0.7253639363  0.6912420534
##  [896]  0.8553322459  0.5962176827  0.1046993800  0.3141083960  0.7283218388
##  [901]  1.0592980446  0.6284405745  0.4846858154  0.6549653525  1.0931844835
##  [906]  1.1027137674  0.2778167423  0.2490895158  1.0535086495  0.4703337786
##  [911]  0.0247782188  1.5242759343  1.6362678696  0.2857021843  0.4088238167
##  [916]  0.8426822200  0.2708580276  0.4456140119  0.7292455858  0.8850033155
##  [921] -0.0693919750  1.0535386090  1.5525658844  1.1617255426  0.7044781904
##  [926]  0.3998563108 -0.1908911638 -0.2359933847  0.6913938405  1.3540782379
##  [931]  0.2464527319  0.4357398833 -0.1013822783  1.1692179029  1.0467375239
##  [936]  0.9299381427  0.8433352660  0.6682758226 -0.2374856068 -0.1634200109
##  [941]  0.4643797625  0.0806683079  1.9161411695  0.6173883906  0.4061840420
##  [946]  0.2482985581  0.9414542366  1.4293303960  0.9323792173  0.0025789797
##  [951]  1.6292635325  0.4959369134  0.0813389592  0.1320099933  0.9052929988
##  [956]  1.1186608588  0.3327633212  0.3226832002  0.2373003377  1.7194604794
##  [961]  0.3068507918  0.4043982704  0.1545127451  0.7979405085 -0.3571186244
##  [966]  0.1944024181  1.1458655679  0.7411502907  0.8455870805  0.1880654038
##  [971]  0.7850519343  0.2272269466  0.6573453289  0.2743359394  0.3949024135
##  [976]  2.4567332828  0.5176671146  0.9075687381  1.2352715852 -0.3416096750
##  [981]  0.4959369134  0.3153434115 -0.2190010025  0.7913422085  0.9580140478
##  [986]  0.3646451602  0.3952416027  0.8250334904  0.3868662217  0.6151287440
##  [991] -0.2434135294  1.1811238200  0.4712283283  0.7351947633  0.5661656310
##  [996]  0.5107833152  0.6553064313  0.5150485490 -0.0360557939  0.5020772804
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
##   0.6408255   0.4127446 
##  (0.1305213) (0.0922908)
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
## [1] -0.8424692  0.5026277  0.3314360  0.5640741  0.1250641  0.2290250
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
## [1] 0.0014
```

```r
sd(results2$thetastar)  #the standard deviation of the theta stars is the SE of the statistic (in this case, the mean)
```

```
## [1] 0.9043539
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
## t1*      4.5 0.02292292   0.9253852
```

One of the main advantages to the 'boot' package over the 'bootstrap' package is the nicer formatting of the output.

Going back to our original code, lets see how we could reproduce all of these numbers:


```r
table(sample(x,size=length(x),replace=T))
```

```
## 
## 3 4 5 6 7 8 9 
## 2 1 2 1 1 1 2
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
## [1] -0.0148
```

```r
se.boot
```

```
## [1] 0.9458624
```

Why do our numbers not agree exactly with those of the boot package? This is because our estimates of bias and standard error are just estimates, and they carry with them their own uncertainties. That is one of the reasons we might bother doing jackknife-after-bootstrap.

The 'boot' package has a LOT of functionality. If we have time, we will come back to some of these more complex functions later in the semester as we cover topics like regression and glm.

