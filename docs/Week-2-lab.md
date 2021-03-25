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
## 0 4 5 6 7 8 9 
## 1 1 1 1 1 3 2
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
## [1] -0.0218
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
## [1] 2.670154
```

```r
UL.boot
```

```
## [1] 6.286246
```

Method #2: Simply take the quantiles of the bootstrap statistics


```r
quantile(xmeans,c(0.025,0.975))
```

```
##  2.5% 97.5% 
##   2.7   6.2
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
##    [1] 4.4 3.4 4.1 4.3 4.2 4.1 4.9 4.8 4.1 4.7 6.0 3.3 3.8 4.1 5.3 4.3 4.0 4.1
##   [19] 6.6 4.4 4.2 4.7 4.1 6.0 6.1 4.1 5.3 4.6 5.0 3.8 3.8 5.5 4.2 5.4 4.6 4.6
##   [37] 5.0 4.7 3.3 5.9 4.3 5.5 4.9 3.2 4.0 3.4 5.3 4.1 4.6 5.3 3.6 2.7 4.3 2.4
##   [55] 2.6 6.2 4.8 4.5 4.2 3.2 4.4 3.8 3.6 4.5 4.9 3.1 5.1 3.4 5.3 5.7 5.6 3.7
##   [73] 3.4 5.3 5.7 4.8 4.9 3.8 4.0 2.7 4.2 5.6 4.3 5.7 5.0 2.5 5.5 6.2 4.7 5.0
##   [91] 5.0 4.6 6.5 5.2 6.1 3.8 5.8 4.9 3.9 4.6 5.4 4.4 3.2 5.5 5.9 3.7 3.3 2.8
##  [109] 5.2 4.4 3.8 5.2 4.1 3.8 4.0 5.8 3.4 3.7 3.5 3.9 4.8 4.9 4.5 4.5 5.6 3.7
##  [127] 6.6 4.5 4.8 4.3 3.6 3.4 5.3 5.5 3.7 4.0 5.6 3.9 5.1 4.7 4.1 5.0 4.7 4.9
##  [145] 3.1 4.3 5.9 4.1 4.8 6.2 5.3 4.1 6.1 3.2 3.8 6.0 3.5 6.7 5.5 4.6 3.2 5.3
##  [163] 4.9 5.7 4.0 4.0 4.0 4.7 4.4 5.1 5.0 2.9 3.5 5.0 4.1 6.3 3.0 4.0 4.3 5.2
##  [181] 3.8 4.4 3.8 4.1 3.6 4.2 6.1 4.5 5.0 5.1 3.7 4.8 4.5 3.2 5.2 5.1 4.4 4.6
##  [199] 4.6 3.8 3.6 5.1 4.3 6.3 3.9 4.0 5.7 3.0 3.8 4.6 3.8 5.3 4.7 3.4 4.6 4.2
##  [217] 4.7 5.1 4.9 4.6 4.7 5.0 4.5 4.7 4.3 4.4 3.7 3.4 5.3 5.2 4.0 5.7 6.5 5.0
##  [235] 2.5 3.7 4.1 5.1 4.7 3.5 6.4 5.3 5.2 5.5 5.3 5.3 4.8 3.9 5.5 3.8 5.7 5.5
##  [253] 4.7 4.1 5.0 3.5 2.9 2.8 5.1 3.0 5.4 3.8 6.2 4.9 4.5 4.4 5.5 2.9 4.2 4.8
##  [271] 5.0 3.2 5.8 3.8 3.8 4.4 5.9 3.7 4.0 3.6 3.7 4.5 4.9 4.3 4.4 3.8 5.1 4.4
##  [289] 4.0 4.1 3.3 5.5 5.3 4.7 4.9 3.9 4.9 5.7 5.1 3.3 3.6 5.2 4.7 2.6 4.2 5.2
##  [307] 3.0 6.2 3.6 5.4 4.1 4.3 5.9 5.8 3.2 2.9 5.1 4.3 2.1 5.1 4.8 4.2 4.2 4.3
##  [325] 6.2 6.0 4.9 4.3 3.3 5.7 4.7 4.6 2.8 4.7 5.8 5.2 3.3 5.5 4.8 4.3 5.8 4.5
##  [343] 3.2 4.8 4.3 4.5 4.7 4.1 5.0 5.5 5.1 3.2 4.3 4.0 4.5 4.8 4.7 4.7 3.3 5.0
##  [361] 3.1 3.3 3.6 4.3 5.0 3.9 3.6 4.6 3.5 5.4 4.7 5.5 3.1 5.5 6.9 2.8 4.3 4.1
##  [379] 4.4 3.3 4.4 4.2 4.2 4.2 3.8 5.0 5.5 4.2 4.9 4.2 5.4 5.5 5.2 4.6 5.0 2.8
##  [397] 5.0 5.2 5.7 4.1 3.8 5.2 4.3 5.9 5.6 3.5 5.9 3.9 3.7 4.4 6.3 3.7 4.9 3.3
##  [415] 5.9 3.2 5.2 5.5 4.4 3.6 3.7 5.1 5.2 4.2 5.0 5.3 3.8 3.5 5.1 5.8 4.6 4.8
##  [433] 4.6 5.1 4.1 4.2 5.6 5.8 4.7 4.3 5.3 3.8 3.8 1.9 5.4 4.8 4.2 5.4 4.9 4.2
##  [451] 5.6 3.1 5.8 3.9 3.6 4.6 4.2 5.3 4.6 3.3 3.3 3.8 4.3 4.5 4.9 4.8 4.1 5.8
##  [469] 3.9 5.8 4.4 3.3 6.0 5.3 4.5 4.5 4.8 4.7 2.5 5.2 2.8 4.4 5.2 6.4 4.3 2.7
##  [487] 4.9 4.3 3.6 3.5 4.5 4.1 5.1 4.7 4.8 4.5 4.6 2.8 4.9 4.0 4.2 2.6 3.7 5.1
##  [505] 4.8 5.3 3.2 3.9 3.9 6.4 6.0 4.9 4.3 5.8 5.4 5.0 3.6 4.2 5.2 5.6 4.3 4.9
##  [523] 3.8 3.7 3.8 3.2 7.2 3.1 4.3 3.9 3.5 5.4 3.9 2.5 4.3 4.9 4.7 5.1 4.7 2.8
##  [541] 5.1 5.0 4.6 5.2 2.2 3.9 2.9 5.5 4.0 4.4 2.4 6.1 4.2 4.4 4.2 5.3 4.7 4.3
##  [559] 4.6 5.1 4.4 4.6 6.0 4.0 5.4 3.0 4.3 4.6 5.7 4.5 3.5 5.8 4.5 5.3 3.4 5.6
##  [577] 5.9 4.0 3.9 4.7 4.2 6.0 4.8 4.6 4.3 4.3 4.5 4.6 3.3 4.5 3.5 3.4 3.7 4.5
##  [595] 4.7 3.8 4.7 3.6 4.9 4.4 5.4 4.2 4.3 4.5 3.5 5.4 3.7 4.7 5.2 4.7 3.7 4.7
##  [613] 6.2 3.9 5.7 5.5 5.7 6.1 4.9 3.4 3.8 4.4 4.9 5.5 5.0 2.4 4.6 5.0 4.1 4.4
##  [631] 4.4 2.4 3.6 4.6 4.9 5.3 4.0 5.8 3.7 5.6 3.7 4.1 4.2 4.3 4.5 5.8 5.8 5.2
##  [649] 4.1 4.3 4.1 4.2 3.7 4.2 4.2 5.8 3.4 4.0 4.9 3.7 5.3 4.5 5.7 3.4 3.5 6.3
##  [667] 6.0 5.3 3.7 3.2 1.8 3.8 4.9 4.9 4.5 5.8 3.8 3.1 4.5 4.6 3.4 4.2 3.5 3.9
##  [685] 6.3 4.2 4.5 3.7 6.1 4.1 3.5 5.5 5.4 4.4 3.7 4.7 2.4 3.9 6.1 4.2 5.6 3.7
##  [703] 5.6 4.3 4.4 4.5 4.6 3.6 3.8 5.0 4.1 3.7 4.0 4.4 4.3 2.3 5.0 5.7 2.9 4.6
##  [721] 2.7 3.1 2.9 3.5 4.7 6.6 4.9 5.9 4.3 6.0 5.4 4.1 4.3 4.3 4.5 4.1 3.5 5.1
##  [739] 3.8 4.6 3.1 3.9 4.7 5.0 4.6 2.0 4.8 3.8 3.3 4.8 3.1 4.4 4.7 4.1 3.5 2.8
##  [757] 4.5 4.9 4.4 4.2 5.1 4.7 2.8 6.4 4.3 2.9 4.8 6.5 4.2 5.7 3.7 3.4 5.0 4.5
##  [775] 4.6 4.7 5.5 4.5 4.8 3.0 4.3 4.2 4.2 5.2 5.4 4.6 5.2 5.8 4.5 3.4 4.1 4.8
##  [793] 5.6 5.0 4.0 4.8 5.3 3.5 4.9 5.7 5.5 5.9 3.7 4.9 3.7 4.3 3.9 4.6 4.3 2.1
##  [811] 4.2 3.1 4.9 4.6 3.1 4.4 3.5 5.4 4.1 5.4 5.3 4.7 5.2 4.0 5.0 3.7 6.1 3.5
##  [829] 5.1 3.6 4.5 2.0 5.1 5.5 4.5 4.5 4.2 5.1 4.6 3.4 3.6 4.0 4.3 6.0 4.2 4.7
##  [847] 5.7 3.3 3.5 5.3 5.4 3.7 4.4 3.2 4.0 5.5 5.8 4.1 6.1 3.6 3.1 5.7 3.7 6.3
##  [865] 2.9 2.9 5.0 4.0 5.4 2.5 4.7 6.3 4.9 4.6 4.1 4.2 5.6 5.5 4.3 4.7 4.0 2.8
##  [883] 4.9 4.9 5.0 3.2 4.8 3.7 4.3 4.0 4.4 3.8 3.4 5.2 4.7 5.0 4.8 5.4 6.4 3.4
##  [901] 4.4 3.7 3.7 4.7 4.8 4.6 4.1 5.3 5.9 6.2 5.4 5.5 5.2 3.3 3.0 3.3 4.5 3.3
##  [919] 4.7 4.7 6.8 5.3 4.9 3.6 4.2 4.2 4.3 4.1 3.7 6.3 5.4 5.5 5.5 5.6 3.5 5.6
##  [937] 3.2 4.4 2.7 5.0 5.2 4.2 4.8 5.4 4.7 4.2 4.3 4.8 3.4 5.1 4.8 3.0 4.9 5.4
##  [955] 4.0 2.7 3.4 4.5 5.3 5.5 4.2 5.4 3.7 5.7 3.7 3.8 3.3 4.4 4.8 4.3 5.3 4.6
##  [973] 3.9 6.2 4.6 4.1 4.9 5.6 4.1 4.6 4.8 5.1 3.4 4.2 4.0 5.3 4.6 3.9 5.2 2.5
##  [991] 5.0 5.7 3.9 3.6 5.4 5.7 5.7 4.1 4.5 3.8
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
##    [1] 5.6 4.9 4.9 4.9 3.9 4.2 4.6 5.0 5.8 3.5 3.6 3.3 6.3 3.1 4.7 3.7 4.6 5.0
##   [19] 4.3 3.4 3.8 4.0 3.6 5.6 4.1 5.4 5.5 2.8 4.8 5.4 6.1 5.3 4.5 4.2 5.1 5.9
##   [37] 4.8 4.8 3.5 3.8 5.5 5.6 5.4 4.7 3.3 3.4 6.1 5.6 4.7 4.6 3.6 3.1 3.3 3.9
##   [55] 5.3 5.1 4.5 5.3 5.5 5.5 4.5 3.6 5.6 5.8 4.6 4.9 4.9 4.6 4.2 5.9 3.0 4.3
##   [73] 3.2 2.9 3.5 5.6 5.4 4.5 5.6 3.7 6.0 4.8 5.5 4.5 6.0 5.4 3.9 1.7 4.7 2.7
##   [91] 3.2 5.2 5.0 3.6 5.9 4.9 4.0 4.6 4.3 4.4 4.8 4.7 3.2 4.9 5.5 4.5 4.2 4.7
##  [109] 4.3 5.6 5.0 4.1 7.2 4.5 5.0 4.5 5.1 5.4 4.3 5.0 4.0 1.6 4.5 5.5 4.5 5.3
##  [127] 5.4 3.6 5.9 3.1 3.7 4.4 2.9 6.1 4.1 4.1 4.4 3.5 5.5 3.2 4.5 4.6 4.1 5.3
##  [145] 4.6 3.6 5.4 5.0 3.8 5.0 4.5 4.1 4.8 3.8 6.6 3.6 4.7 7.8 5.4 5.0 5.6 3.1
##  [163] 4.3 5.3 4.3 5.8 4.1 5.9 4.7 3.7 4.6 4.8 4.0 5.5 3.8 5.3 4.7 5.0 4.3 4.2
##  [181] 5.4 3.8 3.1 2.4 3.6 3.5 5.7 3.5 3.3 4.5 3.2 6.4 5.2 4.1 4.8 3.2 3.6 5.2
##  [199] 4.4 4.0 4.9 4.9 6.2 5.2 3.1 4.8 3.9 5.9 4.5 4.1 5.1 5.6 4.4 4.0 5.1 4.2
##  [217] 3.5 4.1 5.3 4.7 4.7 6.3 4.7 5.0 4.5 3.6 4.1 4.4 2.5 5.2 5.0 4.1 5.3 4.7
##  [235] 3.3 4.8 3.7 4.4 6.5 4.9 5.1 4.0 4.6 4.2 5.5 4.5 3.8 4.4 4.9 3.7 3.6 4.8
##  [253] 4.6 5.9 4.6 2.8 5.4 4.7 4.4 4.9 5.0 4.9 3.8 6.4 5.1 4.7 5.8 4.8 5.8 4.0
##  [271] 4.1 3.7 5.7 4.2 5.2 5.3 4.1 4.0 2.8 4.7 4.0 5.8 4.0 4.1 5.5 4.1 4.0 3.9
##  [289] 4.2 5.5 5.0 2.6 6.6 2.6 5.7 4.7 4.8 4.9 4.4 3.8 2.8 3.1 3.5 3.4 4.7 3.4
##  [307] 4.9 4.0 4.1 4.8 4.9 5.1 5.4 6.0 4.5 5.5 5.6 5.1 3.9 4.7 4.9 2.1 5.1 3.9
##  [325] 4.7 4.5 5.3 5.2 3.1 2.8 5.0 4.6 4.2 3.6 2.8 4.6 4.9 5.0 5.0 3.2 5.3 4.1
##  [343] 5.1 4.4 5.2 3.1 3.5 7.1 4.8 3.6 4.0 3.3 3.7 4.7 4.6 2.1 4.4 4.4 3.9 3.3
##  [361] 4.5 4.2 3.0 4.8 4.8 2.4 5.1 5.1 5.1 2.3 5.6 3.9 3.8 4.5 4.6 5.1 4.7 4.1
##  [379] 5.2 4.9 4.6 3.4 4.5 5.1 2.6 6.0 3.8 4.9 3.8 5.6 5.2 3.5 3.6 5.1 4.7 5.8
##  [397] 4.8 5.3 3.6 5.7 4.0 4.2 3.4 4.1 5.2 4.1 3.3 4.6 4.1 3.2 5.1 5.1 5.0 4.4
##  [415] 4.0 4.1 5.7 3.0 5.9 3.2 3.6 4.2 4.8 4.7 4.4 3.0 5.4 4.6 4.2 4.1 4.1 5.5
##  [433] 3.7 3.1 5.6 5.8 4.8 4.2 3.4 5.3 2.9 4.2 4.2 4.9 6.0 4.4 4.2 2.9 3.2 5.1
##  [451] 5.7 4.2 1.9 3.8 3.8 5.3 3.9 3.2 5.6 3.3 5.1 6.0 3.8 4.6 4.6 5.5 4.4 3.0
##  [469] 3.4 4.1 4.0 4.4 5.3 3.7 4.4 5.3 6.2 4.8 4.8 4.9 4.5 5.3 3.5 5.5 5.7 4.2
##  [487] 3.1 5.2 4.1 4.0 4.7 4.1 5.0 4.7 4.9 3.0 2.7 5.5 4.3 3.6 4.4 2.5 4.6 3.6
##  [505] 3.1 3.1 4.8 4.1 5.5 6.0 5.2 4.1 2.8 4.9 5.4 3.7 3.5 5.0 4.4 4.2 5.4 4.7
##  [523] 4.6 5.3 4.5 3.5 5.4 3.6 3.0 4.6 3.4 5.0 4.9 4.0 5.4 3.9 3.5 4.1 3.0 5.2
##  [541] 5.2 5.0 4.1 4.6 4.1 4.6 4.7 4.6 2.1 3.8 5.3 3.9 3.8 4.3 5.0 4.1 4.5 5.9
##  [559] 4.8 4.9 4.7 4.3 4.9 4.8 4.3 4.9 2.8 4.7 3.8 6.0 1.7 4.5 4.3 4.8 5.3 4.0
##  [577] 6.3 3.7 3.1 4.3 3.1 4.4 3.9 4.4 5.7 4.5 4.6 3.6 5.1 4.8 4.5 4.9 6.8 4.2
##  [595] 4.5 4.3 4.9 6.3 4.8 3.5 4.8 5.1 4.8 3.9 4.7 4.0 4.6 6.3 4.2 4.7 5.2 3.7
##  [613] 4.8 3.9 5.0 5.0 5.1 3.5 3.9 4.6 4.4 6.0 6.0 5.5 3.6 2.5 4.6 5.3 4.7 4.7
##  [631] 4.8 5.5 4.3 2.9 4.2 4.6 2.9 3.7 5.1 5.3 4.0 5.3 2.9 4.5 5.4 4.5 5.7 2.8
##  [649] 3.6 3.1 4.4 4.6 4.3 5.1 3.8 3.6 4.0 5.7 3.6 4.0 4.3 4.1 4.2 3.6 2.8 4.6
##  [667] 5.2 5.9 5.2 4.7 5.4 4.4 7.0 5.4 3.2 3.5 4.5 4.5 4.0 3.2 5.5 5.3 4.1 4.2
##  [685] 4.7 5.8 4.6 5.2 4.0 4.0 4.3 4.9 3.8 4.2 3.4 4.6 3.6 3.7 3.1 5.3 3.2 6.3
##  [703] 4.2 5.4 3.5 5.0 5.1 4.5 3.5 4.1 3.3 5.2 5.4 5.9 3.1 3.5 3.2 4.0 5.1 4.8
##  [721] 3.6 4.4 4.1 5.4 3.6 3.6 6.4 5.1 3.2 5.2 4.8 4.4 4.9 4.0 5.0 4.2 4.6 4.5
##  [739] 4.6 5.6 4.9 6.5 4.1 4.7 6.4 4.8 4.5 3.4 5.8 4.0 5.3 5.9 3.5 4.1 3.1 4.6
##  [757] 3.7 4.9 5.5 4.2 4.6 4.8 3.6 2.6 4.8 3.9 4.4 5.9 5.0 4.2 5.6 6.8 3.6 5.9
##  [775] 4.7 3.7 3.2 6.1 3.8 4.7 5.0 3.8 4.7 3.2 3.6 5.7 4.5 5.1 5.7 4.4 4.8 3.7
##  [793] 4.8 4.0 4.6 4.5 5.2 4.2 5.3 3.0 4.2 5.7 4.0 3.9 4.1 3.5 5.4 4.3 5.4 4.1
##  [811] 5.1 2.4 3.0 4.0 4.6 4.7 3.9 2.7 4.1 3.4 4.0 3.9 3.0 3.9 3.6 5.9 3.7 4.0
##  [829] 3.6 2.1 4.8 4.8 4.8 3.8 3.6 5.4 3.8 4.3 4.3 5.9 5.6 5.2 5.2 6.5 4.3 4.9
##  [847] 4.9 5.1 5.1 6.9 4.5 5.7 3.7 5.4 5.2 4.3 4.4 5.0 4.7 3.3 5.2 4.6 4.7 5.0
##  [865] 3.8 3.5 3.5 3.1 3.2 4.6 4.7 4.4 4.8 3.8 3.7 5.9 3.6 6.2 4.4 5.0 2.6 4.9
##  [883] 3.3 5.0 3.8 3.9 5.2 4.8 4.9 5.3 4.2 4.4 4.8 5.0 5.0 4.5 5.1 4.7 3.8 6.0
##  [901] 4.7 4.2 4.8 4.7 3.6 6.0 3.4 4.9 5.0 4.1 4.4 4.7 5.1 4.9 5.9 4.0 4.6 5.4
##  [919] 4.0 5.8 4.3 5.7 3.2 6.0 5.0 5.6 3.9 5.2 4.1 5.3 5.6 4.8 3.5 4.1 2.5 3.5
##  [937] 2.1 2.2 3.1 3.5 5.9 4.3 4.7 4.6 4.3 3.6 5.0 4.0 4.3 4.6 4.6 5.1 4.1 5.1
##  [955] 3.7 4.4 6.0 5.4 3.8 5.9 2.6 3.2 5.1 5.5 3.7 4.6 4.8 4.0 3.6 5.2 5.3 4.5
##  [973] 6.8 2.2 4.2 4.6 3.6 4.0 6.8 3.5 5.4 2.4 5.8 4.5 5.7 2.1 5.7 6.1 4.4 3.4
##  [991] 5.9 3.9 5.1 5.7 2.9 4.7 3.6 4.8 4.3 4.3
## 
## $func.thetastar
## [1] -0.0193
## 
## $jack.boot.val
##  [1]  0.497935103  0.392063492  0.227200000  0.171764706  0.001108033
##  [6] -0.079670330 -0.164285714 -0.272849462 -0.353846154 -0.607784431
## 
## $jack.boot.se
## [1] 0.9866408
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
##    [1] 4.4 5.2 4.1 3.4 3.9 5.7 3.0 5.6 4.7 4.9 3.9 2.3 5.3 5.6 4.4 4.9 5.5 5.6
##   [19] 5.2 5.2 5.7 4.2 5.5 3.3 5.9 3.9 4.7 5.4 5.5 3.8 4.3 4.2 4.1 4.7 6.1 5.0
##   [37] 4.7 5.4 4.8 4.7 4.6 5.1 4.6 5.0 4.7 4.2 7.5 3.9 4.6 2.8 5.7 4.1 4.0 4.5
##   [55] 4.5 6.3 3.5 5.0 5.5 4.0 3.9 3.7 6.1 4.5 4.4 3.4 3.8 3.7 3.3 3.7 4.0 4.0
##   [73] 2.7 5.1 4.3 3.9 4.0 3.9 3.9 4.1 5.6 5.3 3.2 5.1 4.8 4.1 4.9 6.1 5.0 4.4
##   [91] 3.4 4.4 5.1 3.2 4.9 6.0 4.7 4.8 4.4 4.5 3.4 4.6 4.6 3.4 5.2 4.1 4.0 3.7
##  [109] 3.0 3.8 5.8 4.3 4.3 3.5 6.3 5.6 6.8 3.7 4.0 5.1 4.4 3.6 5.6 5.5 3.3 4.0
##  [127] 4.4 2.9 5.7 6.6 5.5 3.5 3.6 5.4 4.1 2.8 2.8 5.1 4.7 4.1 4.9 5.4 3.7 4.8
##  [145] 6.1 4.9 5.3 4.9 4.1 4.5 3.6 4.9 4.2 3.8 5.5 3.9 3.7 4.9 5.0 4.6 3.1 4.1
##  [163] 6.0 2.9 3.3 5.6 5.7 3.9 5.0 4.7 3.0 4.1 4.2 4.7 4.8 3.1 4.5 4.0 3.1 4.0
##  [181] 4.6 4.7 4.9 5.5 3.9 4.7 4.0 5.2 4.0 3.8 4.5 5.0 2.9 4.5 4.1 4.0 4.6 3.0
##  [199] 3.6 3.6 4.6 4.2 4.7 4.1 6.0 4.3 5.3 2.7 3.7 4.7 5.0 4.5 3.9 5.4 4.2 4.8
##  [217] 3.6 2.8 4.8 4.6 3.7 4.2 4.6 5.0 4.7 3.3 7.5 4.1 5.9 5.2 3.9 5.2 4.5 6.4
##  [235] 4.2 5.1 4.3 4.0 4.7 6.1 5.2 4.1 3.7 5.4 3.5 4.3 4.7 5.1 3.8 4.2 5.4 5.6
##  [253] 4.6 5.4 3.3 3.6 3.3 4.0 3.9 2.5 4.1 5.2 4.6 3.3 3.5 5.0 4.6 5.0 4.5 4.1
##  [271] 3.3 3.7 4.3 4.1 4.4 3.8 4.7 5.0 5.8 4.8 6.1 4.2 5.0 5.8 6.0 3.8 4.9 4.7
##  [289] 4.8 7.0 4.0 4.9 4.5 4.8 4.4 4.1 5.6 5.6 3.5 4.0 4.8 4.3 4.5 4.2 4.0 3.8
##  [307] 4.2 5.5 4.2 6.3 3.8 6.7 4.3 6.4 5.2 2.6 4.3 4.2 3.0 4.2 2.5 6.4 5.2 4.8
##  [325] 7.1 6.2 5.6 2.4 3.9 3.1 3.9 6.0 5.0 4.2 3.2 5.0 5.9 4.5 3.8 4.9 4.3 4.5
##  [343] 5.7 5.9 4.6 3.0 2.8 5.6 4.8 4.0 2.8 3.0 3.2 4.8 4.5 5.4 4.5 3.9 4.3 5.1
##  [361] 4.8 5.5 2.3 2.4 3.3 5.4 3.9 4.7 3.3 4.0 5.7 4.8 3.3 3.2 5.9 3.7 4.7 3.4
##  [379] 6.1 3.9 5.7 5.4 4.5 4.5 3.6 4.1 4.3 3.2 4.8 3.9 3.1 2.5 4.3 4.5 2.7 5.1
##  [397] 5.2 4.1 3.7 5.4 4.8 3.9 5.7 3.5 5.3 4.2 3.7 4.8 2.1 5.0 5.0 6.2 4.6 5.5
##  [415] 5.6 4.2 4.1 3.5 3.8 5.9 3.8 3.3 5.3 2.7 3.7 3.8 4.1 4.7 2.9 4.8 4.0 5.2
##  [433] 5.8 6.0 4.8 4.4 4.2 3.6 5.1 4.3 5.5 4.5 3.7 3.8 4.0 5.5 4.6 3.8 4.1 4.0
##  [451] 2.5 6.0 3.3 4.1 3.7 4.1 3.0 4.6 5.3 5.3 4.8 4.9 4.4 5.5 4.5 4.9 4.5 4.4
##  [469] 4.3 4.1 5.3 4.8 4.1 5.1 5.7 5.5 2.8 3.4 4.2 4.6 5.1 3.9 4.4 5.5 4.9 4.4
##  [487] 3.4 4.3 3.9 4.5 4.0 5.4 2.7 4.3 4.2 4.8 5.4 3.9 5.6 4.2 3.8 3.4 4.8 5.8
##  [505] 4.8 5.7 4.4 6.1 5.4 3.9 4.9 2.9 3.7 3.8 3.6 4.7 5.1 4.4 4.8 3.4 4.3 5.0
##  [523] 5.2 4.2 6.3 4.6 5.7 4.5 5.0 4.0 5.5 5.2 4.5 4.2 5.5 3.5 3.3 4.4 3.6 5.7
##  [541] 4.6 5.7 5.5 5.0 5.3 3.9 4.2 6.1 5.0 6.5 4.9 4.0 3.9 3.9 4.0 5.1 3.6 4.1
##  [559] 4.5 4.4 3.4 2.8 5.6 4.8 4.8 4.9 4.8 4.4 5.2 5.1 3.8 6.4 5.1 6.1 3.6 4.1
##  [577] 4.3 5.3 4.6 3.5 4.9 5.3 2.6 5.7 3.2 2.2 3.9 1.7 4.0 4.0 4.2 5.2 5.9 3.0
##  [595] 6.4 2.6 4.3 4.3 3.5 3.9 3.6 4.5 4.4 4.3 5.0 3.0 5.3 5.7 4.7 4.6 5.1 3.8
##  [613] 4.7 6.8 3.8 3.8 4.6 4.1 4.1 5.7 5.0 3.3 5.0 3.8 4.4 4.7 4.6 4.3 4.5 6.4
##  [631] 4.9 5.0 5.1 3.6 4.1 3.9 3.8 3.4 4.2 4.7 4.9 4.1 4.5 4.4 4.4 3.0 3.1 4.5
##  [649] 4.9 4.1 4.7 3.4 4.2 5.9 3.6 5.8 4.2 4.0 5.8 2.7 4.5 4.7 3.1 3.7 5.1 4.6
##  [667] 3.9 3.6 4.6 5.1 3.2 2.7 5.2 6.4 3.4 4.9 5.5 4.7 4.0 4.1 4.4 5.4 3.0 4.9
##  [685] 5.6 3.8 3.9 5.2 5.0 4.6 6.1 3.2 5.3 4.8 6.3 4.2 3.7 4.1 6.2 5.2 2.5 4.4
##  [703] 5.1 5.0 4.7 4.2 4.2 3.0 4.3 5.0 4.8 3.2 4.8 3.4 4.2 4.7 3.0 4.6 5.4 5.1
##  [721] 4.9 5.4 4.8 4.5 3.1 5.6 5.8 4.1 4.9 4.3 6.6 4.5 4.0 3.4 5.7 5.3 5.2 6.7
##  [739] 5.3 3.6 5.3 4.5 3.8 5.6 3.4 4.4 5.3 5.7 4.9 4.4 4.2 4.0 3.8 5.3 4.4 4.4
##  [757] 3.0 4.2 5.9 4.5 4.0 4.5 4.5 4.9 4.6 5.3 3.4 5.1 3.2 5.9 3.1 5.3 2.6 4.8
##  [775] 4.4 4.4 5.2 4.4 2.4 2.7 4.2 5.3 3.9 5.5 4.0 3.9 5.0 4.5 4.7 4.8 4.1 4.0
##  [793] 5.3 5.1 4.8 5.4 3.9 4.8 4.8 4.3 5.1 6.8 4.5 5.4 5.0 3.4 5.7 3.7 4.5 4.6
##  [811] 4.3 4.0 4.7 7.2 4.1 3.8 4.4 3.0 5.8 4.8 3.5 3.9 6.5 4.2 3.6 2.6 3.9 4.3
##  [829] 4.2 5.9 4.0 3.6 2.8 5.2 5.0 4.0 3.9 2.0 4.8 4.4 4.7 2.6 4.7 4.7 4.0 5.1
##  [847] 4.1 4.3 4.7 5.3 4.8 3.0 5.5 4.9 4.6 5.0 2.8 3.2 3.8 4.6 4.7 4.4 4.7 5.1
##  [865] 4.6 4.9 4.5 5.3 3.5 5.2 5.7 6.8 4.5 5.1 5.0 4.9 3.5 5.3 5.2 6.0 3.4 4.9
##  [883] 4.0 4.3 5.1 4.6 3.7 5.6 4.3 6.1 3.8 3.4 5.7 6.0 2.8 5.0 5.2 6.4 5.4 4.5
##  [901] 4.4 4.8 4.9 4.5 2.8 2.5 3.2 4.1 4.2 4.1 4.2 3.7 5.8 4.6 6.1 5.8 5.3 4.1
##  [919] 5.1 3.2 5.4 4.6 3.2 3.7 4.3 4.5 5.1 5.5 6.2 3.4 5.5 5.1 5.0 5.3 3.2 4.0
##  [937] 4.7 5.3 6.0 3.5 3.9 3.8 3.7 5.4 4.1 4.1 7.9 4.5 5.6 4.5 5.0 3.7 3.3 4.0
##  [955] 4.6 4.0 5.7 3.3 3.7 4.8 3.5 4.3 5.3 5.4 6.9 4.7 4.4 4.2 3.5 5.2 5.4 4.0
##  [973] 5.1 4.7 4.8 3.3 5.1 4.2 5.5 4.7 5.0 6.0 3.6 3.5 3.5 4.3 3.3 4.1 5.9 4.5
##  [991] 3.6 4.1 3.9 4.6 4.2 4.1 6.0 5.1 5.9 4.3
## 
## $func.thetastar
## 72% 
##   5 
## 
## $jack.boot.val
##  [1] 5.500 5.400 5.300 5.136 5.100 5.000 4.900 4.700 4.600 4.500
## 
## $jack.boot.se
## [1] 0.9659098
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
## [1] 0.5392488
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
##   7.180400   7.444916 
##  (3.139401) (3.371540)
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
## [1] 0.09386922 0.64453859 0.84824535 0.57676378 0.24805673 1.10973290
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
##    [1]  0.484836895  1.027814910 -0.066609957 -0.196405057  0.440966960
##    [6]  0.085213990 -0.212451724  0.251564955 -0.542107114  0.473543003
##   [11]  1.266409756 -0.351676175  1.004112669  0.035943000  0.681302777
##   [16]  0.777455821  0.089717400  0.235910748  0.078544371  0.370361291
##   [21]  0.463657394 -0.119195887 -0.672285984  0.380358483  0.681811744
##   [26]  0.602380994  0.955461635  1.002607906  0.547451973  1.062887214
##   [31] -0.233631581  1.275541620 -0.612874838 -1.181696129  0.245136499
##   [36] -0.192094040  0.254816349  0.607170591  0.372150527  0.242898907
##   [41]  1.224764696  1.552535423  0.367300909  0.317720393 -0.064855172
##   [46]  1.002088130  0.231075319  0.402084252  0.840780120 -0.073093898
##   [51]  2.013350328  0.801134603 -0.587514442  0.186485816  1.100172471
##   [56]  0.829946294  0.787283094 -0.251222724  0.442164295 -0.082918478
##   [61] -0.214427451  1.534057283  0.893013367  0.372479626  0.472445821
##   [66] -0.072389428  0.372370708  0.544111329  1.086509110  0.574557240
##   [71] -0.045355607  0.381656197  0.770518192  0.394747839  0.819927729
##   [76]  1.153920025  0.493124947  0.610021891 -0.320767728  0.993816698
##   [81]  0.761832608  1.122241434 -0.141288571  0.583238076 -0.250853250
##   [86]  0.136598754  0.578398039 -0.125865695  0.983362017 -0.374327561
##   [91]  0.116254193  0.385510357  1.359697987 -0.346711015  0.594914763
##   [96]  0.600666828  2.125311022 -0.194186866 -0.076303341  0.954693019
##  [101]  0.674843985 -0.133210880  0.280511983 -0.825118876 -0.400019819
##  [106]  1.086656556  0.267399065  0.178410985 -0.423921749  0.585375479
##  [111] -0.383392428  1.755966947  0.812364348 -0.688877382  0.609337616
##  [116]  0.377328505  0.496884372  0.095803939 -0.395727857 -0.022896727
##  [121]  0.797746264 -0.175077709  0.628371344 -0.565643266  0.751203876
##  [126]  0.372318703  1.359032417 -0.640367507 -0.548551268 -0.195611301
##  [131]  0.469149226  0.159305321 -0.803039209 -1.584775481 -1.135920381
##  [136] -0.346711015  0.236625917  0.234224087 -0.458461810  0.360646055
##  [141]  0.925805546  1.412859452  0.388513181 -0.674885673 -0.658157448
##  [146] -0.413124016  0.653118825  0.639967286  0.522777495  1.064329658
##  [151]  0.284694102 -0.598514918 -0.151140419 -0.228575246  1.107405203
##  [156]  1.387591839  0.561288045  0.440117113  0.967243874 -0.460659971
##  [161]  0.680933775  0.349331472  0.169405537  0.843236984 -0.066876875
##  [166]  0.859985265  0.074335634  1.226496637  0.362216467 -0.334100425
##  [171]  0.339772118 -0.174566723  0.335065945  1.723549949  0.624902188
##  [176]  0.568990651  0.027457947  0.181243683  0.918988059  0.921109203
##  [181]  0.547264727  1.091797164  0.664600253  0.713298215  0.106358747
##  [186]  0.057815775 -0.064339574  0.944238415 -0.461442409  1.142015288
##  [191]  0.847526531  0.298320636  0.014224994 -0.325544871 -0.105640187
##  [196] -0.422338457  0.607925452  0.882066016  0.150930452  0.884001917
##  [201]  1.197047404  0.497591414  0.485776596  0.827970316  0.147480596
##  [206]  0.163147737  0.345002351 -1.219462480 -0.856336084  0.440966960
##  [211]  0.916416607  0.677425615  0.665635923 -0.552463793  0.713851150
##  [216]  0.670989331  1.426397180  0.213471456 -0.261237803 -0.180444087
##  [221]  1.256972736 -0.296836076  0.643092624  0.058219259  1.360339677
##  [226]  0.432498296  0.357763776 -0.152390779 -0.507210417  0.876436999
##  [231]  0.878600655  0.289321014  0.594914763 -0.805239791  0.726109711
##  [236]  0.584965443  0.136031798 -0.358527586  0.179828664  0.384599096
##  [241] -0.223123961  0.386586931 -0.392051396 -0.409350334  0.533463875
##  [246]  0.759972925 -0.610861687  0.846956818 -0.353505628  0.144914763
##  [251] -0.407748614 -0.188819270  0.237558251  0.247955846  0.327356057
##  [256]  0.907178573  0.552794898  0.555027883  0.674727438 -0.629611031
##  [261]  1.156701892  0.730494008  0.715431057  0.223966361  0.586456914
##  [266]  0.772197426  0.078267903  0.733224485 -0.214533756  0.372888598
##  [271] -0.439006111  1.096933490  1.382774121  0.957573015  0.008896370
##  [276]  0.933941159  0.474311489 -0.159317558  0.728603053  0.688450418
##  [281]  0.770126897  0.595604939 -0.624184996 -0.832308496 -0.914638842
##  [286]  0.362216467  0.820116021  0.964359759  0.746954414  0.586497549
##  [291] -0.055346131 -0.047440496  0.991016417  0.493968551  0.688845157
##  [296] -0.316128091 -1.117306980  0.346958864 -1.174639446  0.041936213
##  [301] -0.943996759 -0.265891481  1.136464733  0.169428989  0.873194024
##  [306]  0.773482343  0.167216074 -0.552197478 -0.481600444  2.494965058
##  [311]  1.022553268  0.059562294 -0.817312602  0.390540383  1.001966393
##  [316]  0.442713423  0.067026343  0.844689460 -0.640540558  0.069898149
##  [321] -0.729513788  0.592384213  0.512871209  0.067371040  0.437098192
##  [326]  0.349581997  1.060721641  0.650522838  0.748452499 -0.501933026
##  [331]  0.878358136 -0.422175173  0.372150527  0.994725523  0.087419783
##  [336]  0.323338641 -0.080909516  0.329971027  0.808101402  0.176731004
##  [341] -0.516789175  0.880826754  0.472867788  0.632353185  0.545087639
##  [346]  0.767787021  1.031651943  0.748932610 -0.048919703  0.042600362
##  [351]  1.451750436  0.391812647 -0.885968234  0.683982819 -0.163016669
##  [356] -0.190472936 -0.042225751  0.278547352  0.016498282  1.021387948
##  [361] -0.317735999  1.587029086  0.348432348  0.530185244  0.129780191
##  [366]  0.513766384 -0.455734755 -0.878395313  0.097808437  1.037901710
##  [371]  0.726907683  0.952143908 -0.673391749  0.257557550 -0.171623318
##  [376]  1.070890319 -0.687053150  0.173624926  0.202768211  1.344681576
##  [381]  0.774727902  0.769078827  0.452716976  1.445134381  0.041716081
##  [386] -0.241699851  1.140112802  1.043910825 -0.426646179 -0.061440555
##  [391] -0.293032047  1.104706897 -0.252575315 -0.779945095 -0.665766915
##  [396]  1.408619356  0.415812278 -0.356629520  0.562534080  0.074789839
##  [401]  1.073141236  0.519619838 -0.911186350  0.468246743  0.162206573
##  [406] -0.262265719 -0.006283538 -1.126108254  0.288982569  0.318551108
##  [411]  0.266265848  0.089548111  1.076189767  0.880614674 -0.080206423
##  [416] -0.231812058  1.686013326 -0.630201783 -0.598209807  0.057069566
##  [421]  0.862553375  2.240583993  0.171031852  1.253957166  0.084246712
##  [426]  0.436474764 -0.839551207  0.014813393 -0.201299566 -0.297249818
##  [431] -0.720097868 -0.091650381 -0.017214212 -0.383360410  0.941573047
##  [436] -0.259022070 -0.070940821  0.445376085 -0.544633433 -0.080163407
##  [441]  0.823645176 -0.064077267 -0.429547436  0.133690194  0.630404330
##  [446]  0.462949931 -0.284404859 -0.548443658  0.395882665  0.158813855
##  [451]  1.119405680  1.170601688 -0.009312875 -0.726431286 -0.336514802
##  [456] -1.431847772  0.360369665 -0.112369660  0.673829117 -0.156159148
##  [461]  1.002088130 -0.344001553 -1.385310209  0.251936167  0.191826570
##  [466] -0.378514738  0.648536113  0.372639921  0.042364537  0.186288617
##  [471]  0.597619358  0.687263300  0.896261518 -0.422836364  0.472011548
##  [476] -0.355319286  0.915890010  0.269007225 -0.220861057  0.157573933
##  [481]  1.236910203 -0.126228168 -0.085329676 -0.166158271 -0.213598880
##  [486] -0.603256696  0.266383881  0.854352501 -0.271584684  0.817124717
##  [491]  0.214555961 -0.734230442  0.477030298  0.511423329 -0.501118169
##  [496]  0.799394262  0.400832934  0.684538426 -0.119326610 -0.088177335
##  [501] -0.017719629  0.925620145  0.458000362  0.640680907  0.052049858
##  [506]  0.412708131  0.704372503  0.536038478 -0.305515944  0.154348111
##  [511]  0.208741600  0.697572357  0.190187465  1.285962924 -0.159102633
##  [516] -1.348272032 -0.039451543  0.816818277  0.935213154 -1.288629456
##  [521]  0.968428678  1.061232717  1.010248313  1.164307560  0.119056723
##  [526]  1.232365035  1.232718710 -0.044165068  0.460064043  0.347805298
##  [531]  0.005831996  0.978624656  0.718703657  0.211085795  1.410351839
##  [536]  0.305568777  0.395394829  0.997408071  0.539074652 -0.199616646
##  [541]  1.128560563 -0.713186393  1.301630719  0.572535431  0.162621689
##  [546]  0.496468408  0.586497549  0.649521001  0.363797704 -0.403186729
##  [551] -0.262681667  0.133971685  0.501622718  0.834923334 -0.227953671
##  [556]  0.583172213  0.942066000  0.105457574  0.353244852 -0.604057727
##  [561]  0.003660851  0.840834693 -0.012822447  1.080848484 -0.875858420
##  [566] -0.174566723  0.590885776  1.790776767  0.243555348 -0.197304206
##  [571] -0.563814280 -0.836679743 -1.419725425  0.323705509  0.519694169
##  [576] -0.141268650  1.121971118  0.872667098 -0.466546545  1.771063269
##  [581]  0.602379950  0.095383832  0.748441553  0.397883708  0.673286756
##  [586]  0.711821812  1.275896535 -0.295525344  0.438774402 -0.799722663
##  [591]  0.015159632  1.253341618  0.806830603  0.328910501  0.464416685
##  [596] -0.493911854  0.078114178 -0.028405208  0.393188536 -0.446120054
##  [601] -0.020607327 -0.295583008  1.014506813 -0.485243626  0.597649969
##  [606]  1.219834747 -0.315747335 -0.429883581  1.189775473  0.163691265
##  [611]  0.677944222  0.942604058  0.769755769  0.718018850  0.269623031
##  [616] -0.345331306 -0.059432617  1.100488220  0.117064017  0.164111663
##  [621] -0.337428190  0.462008717  0.684291722  0.330328406  1.061994293
##  [626]  0.202080750 -0.940582478  0.470983056  0.295994257  0.269121641
##  [631]  1.152818850 -0.157159793  0.404646130  1.195261790 -0.214436764
##  [636] -0.726551244  0.390643865  1.194034989 -0.548055009  0.565369062
##  [641]  0.878358136 -0.578552932 -0.065128459  0.696631231  1.769175794
##  [646]  0.904599785 -0.266739128  0.267159843  1.074027163  0.266059057
##  [651]  0.854160151  0.115751189  0.866915209 -0.064077267  0.648194465
##  [656]  1.195829001  1.522014277 -0.302918646  0.501438211  1.020426547
##  [661]  0.011475527  0.548549455 -0.687840173  0.053779162  0.587618077
##  [666]  1.326605245  0.484426394  0.516007415  0.671066125  0.941764210
##  [671]  0.380169360  0.469709226  0.492383587  0.641337660  0.207529823
##  [676]  1.361471851  1.511683152 -0.690186441 -0.775203625  0.771324726
##  [681]  0.060860669 -0.424509627  0.729227444  0.274144658 -0.276910511
##  [686] -0.193788940 -0.089599282  0.059477891 -0.017025824  0.963174421
##  [691] -0.217050073  0.931450396  0.479327360  1.180197749 -1.630642738
##  [696]  0.056531750 -0.397425256  0.576951436  0.784227556 -0.931506915
##  [701] -0.661964648  0.091172378  0.310261367  0.716152420  0.845046563
##  [706] -0.171741361  0.831005715  0.658413877  1.189358017  1.435389905
##  [711]  0.930569344  1.006274356 -0.118518447 -0.440306167 -0.247480379
##  [716] -0.426737865  0.572022790  0.257042934  1.188264131 -0.500312640
##  [721] -1.093536987 -0.428358251  0.842392728  1.116049624  0.461664551
##  [726] -1.173863444 -0.424509627  0.088515054  0.479968108  0.495708638
##  [731] -0.004754842 -0.429566220  0.691374407 -0.786843220 -1.214535739
##  [736]  0.320458078 -0.126169371  1.581535271  0.265388091  0.965986334
##  [741]  0.763948963 -0.126084092  0.013860521 -0.052467632  0.463139197
##  [746]  1.227643584  0.658049436  0.709462540 -0.464891326  0.827997362
##  [751]  0.209450909 -0.384095217 -0.546235147 -0.184125641  1.174299675
##  [756]  0.799587454 -0.153944790  1.085762343  0.683982819  0.842477652
##  [761]  0.547567369  0.145717176  0.480159928 -1.000684857  0.437478134
##  [766] -0.043998736  0.400926765  0.210630895  0.185517727  0.265138730
##  [771]  0.622016847  0.420119416  2.339924381  0.115314293  1.431042018
##  [776]  0.154669209  0.601620679  0.362741426 -0.502713524 -0.945698387
##  [781]  1.686013326  0.754364265  0.739835429  0.173627633  0.287418192
##  [786] -0.729447442  0.115951020 -0.214533756  0.629990757  1.323394285
##  [791] -0.732536954  1.145997537  0.900558720  0.918091750 -0.904485512
##  [796]  1.231057425  0.433655602 -0.715917450  0.504187414  0.280795343
##  [801] -0.064279424 -0.691384834 -0.071989241  0.731025303  0.539248750
##  [806]  1.056686508  0.013394109  1.006367170  0.853187368  1.583613211
##  [811]  0.831342545  1.014601749  0.708876285  0.250309394  0.736454162
##  [816]  0.919527911 -0.498404577 -0.434885427  0.286445049 -1.131032332
##  [821]  0.599244255  0.826242352 -0.498185522  0.527578010  0.809184460
##  [826]  0.626741745  1.255768172  0.312022398 -0.190620579 -0.046320058
##  [831]  0.031889778  0.964980512 -0.031843717  0.383444242  0.919248727
##  [836] -0.186676114  0.707472942  0.462949931  0.962879084  0.385510357
##  [841]  0.384113786  0.622448983  0.964910737  0.799662492 -0.300556794
##  [846] -0.152390779  0.192618343 -0.198107192 -0.126097845  0.535912036
##  [851]  0.563061117 -0.310184090  0.665092793  0.074645501  0.677583659
##  [856]  0.266307011 -0.793150891  0.595076491 -0.499372659 -0.184006060
##  [861]  1.491721723  0.586660320 -0.274143365  0.444566203  0.146208295
##  [866]  0.935129350  0.074186506  0.908412374  0.461758689  1.082628836
##  [871]  0.560300541  1.412217067  1.741299811  0.283602455  1.076062055
##  [876]  0.831487973  0.556657917  0.250628698 -0.511213147  0.853532109
##  [881]  1.153555074  0.339507433  1.143583191  0.318144316 -0.097029036
##  [886]  0.554097388  0.403768738 -0.201844982 -0.202808766  0.825549209
##  [891] -0.067212638  0.035889785 -0.825629387  0.379549025  0.552658320
##  [896]  0.219657084  0.560530927 -0.303848027  0.769622705  0.464234367
##  [901]  0.401557702  0.903897697  1.207276389 -0.507116028 -0.080598387
##  [906]  1.285962924  0.640988206  0.626461731  0.883865217  0.264451345
##  [911]  0.773150533  0.765179547 -0.919026580 -0.180310801  0.180028068
##  [916] -0.499418955  0.306212427  0.738636861  0.242347353 -0.053610052
##  [921]  0.983728289 -0.147132138  0.614020828 -0.817497062  1.098929291
##  [926]  0.346726112  0.206435894  0.177096029  0.437938343  0.072503436
##  [931] -1.019336304  0.110290063  0.119271993  0.774546741  0.565567332
##  [936]  0.402605098  1.466616610 -0.179757742 -0.390300279  1.269708865
##  [941]  0.347276736 -0.044485992 -1.118101681  0.726059227  0.117914987
##  [946]  0.213610399 -0.485741309  1.083288021  0.882339111 -0.292273463
##  [951] -0.317216600 -0.266047442  0.465495797 -0.259330306  0.403679462
##  [956]  1.146213043  0.194048084 -0.587015200  1.720456313  1.014158300
##  [961] -0.351569633  0.181785850 -0.013437960 -0.068347731  1.037590956
##  [966]  0.637891431  1.047109961  0.993506606  0.786479958  0.036742276
##  [971]  1.272929964  0.602973563  0.249720958  0.640439127 -0.392093716
##  [976]  0.744075287  0.118318027  0.810181387 -0.549347586  0.310538531
##  [981]  0.501679164  0.433898053  0.125111986 -0.557238772  1.010661994
##  [986]  1.022435465  0.342807602 -0.456520750  0.482023098 -0.151603874
##  [991]  0.588191075  0.006773802  0.204216815  1.248080763  0.732861016
##  [996]  0.606506201 -0.785039809 -0.085889174  0.873648566 -0.478244999
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
##   0.96448187   0.35231474 
##  (0.11141170) (0.07877798)
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
## [1] -0.5472016  0.1995965  0.5910535  0.1574606  0.3470628  0.2310687
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
## [1] 0.0282
```

```r
sd(results2$thetastar)  #the standard deviation of the theta stars is the SE of the statistic (in this case, the mean)
```

```
## [1] 0.8835187
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
## t1*      4.5 0.01641642   0.9029848
```

One of the main advantages to the 'boot' package over the 'bootstrap' package is the nicer formatting of the output.

Going back to our original code, lets see how we could reproduce all of these numbers:


```r
table(sample(x,size=length(x),replace=T))
```

```
## 
## 1 2 3 4 5 6 8 9 
## 1 1 1 2 1 1 1 2
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
## [1] -0.0139
```

```r
se.boot
```

```
## [1] 0.8925993
```

Why do our numbers not agree exactly with those of the boot package? This is because our estimates of bias and standard error are just estimates, and they carry with them their own uncertainties. That is one of the reasons we might bother doing jackknife-after-bootstrap.

The 'boot' package has a LOT of functionality. If we have time, we will come back to some of these more complex functions later in the semester as we cover topics like regression and glm.

