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
## 0 1 2 3 4 5 6 
## 3 1 1 1 1 1 2
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
## [1] 0.0072
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
## [1] 2.755189
```

```r
UL.boot
```

```
## [1] 6.259211
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
##    [1] 5.2 2.3 5.9 4.8 5.4 4.1 3.5 4.9 3.7 4.4 5.2 6.1 3.5 6.0 3.6 4.0 4.5 3.5
##   [19] 4.5 4.3 3.0 5.2 5.8 3.7 6.1 4.9 4.2 3.8 4.4 3.4 4.8 3.6 4.4 3.5 3.9 5.0
##   [37] 3.2 3.2 4.7 3.6 4.7 6.4 3.8 3.6 4.0 4.8 4.4 4.7 4.8 4.0 5.6 5.5 4.7 5.4
##   [55] 3.5 5.4 4.0 4.7 5.2 3.4 6.3 5.7 4.2 4.6 3.4 5.1 4.5 4.7 3.2 3.5 3.8 5.1
##   [73] 5.4 5.3 3.3 4.5 4.6 4.9 4.9 4.5 4.6 4.2 2.7 4.1 6.4 5.1 5.2 5.7 5.5 5.3
##   [91] 4.8 4.1 4.8 3.1 4.1 4.0 3.8 2.9 4.2 5.5 4.9 3.4 3.3 4.1 4.8 4.2 4.2 3.8
##  [109] 5.7 5.9 4.5 4.1 3.4 4.1 5.5 5.2 4.9 4.7 6.1 4.0 4.6 5.6 4.5 4.5 3.9 3.6
##  [127] 4.6 4.3 3.8 4.1 4.5 5.4 5.2 3.9 4.9 4.1 5.7 3.8 3.0 2.8 3.4 5.6 4.8 4.8
##  [145] 5.4 5.1 5.7 2.2 3.3 3.1 4.6 3.7 4.6 4.8 5.6 6.7 4.4 5.8 4.2 3.3 4.9 3.8
##  [163] 4.4 6.1 5.4 5.2 4.6 3.7 5.3 4.2 4.4 4.0 6.8 3.2 5.2 5.4 3.8 5.2 3.8 5.9
##  [181] 4.3 5.9 4.8 4.7 4.7 4.0 3.9 5.5 4.6 3.4 3.3 2.9 3.7 3.5 4.0 4.6 3.1 5.4
##  [199] 5.9 4.4 4.4 3.9 5.1 4.3 4.4 3.5 4.3 5.8 6.1 3.5 4.2 4.7 4.9 5.2 3.4 3.8
##  [217] 5.8 6.0 4.1 5.3 6.1 4.6 4.6 4.4 5.4 5.3 3.7 4.1 4.4 4.4 4.1 3.2 5.2 4.9
##  [235] 4.2 4.2 4.9 3.8 4.1 4.4 4.1 4.8 3.6 4.7 4.0 4.3 6.1 3.9 3.8 4.8 2.8 3.9
##  [253] 2.2 5.6 3.5 4.0 5.3 4.3 4.1 5.0 4.7 3.5 4.1 6.2 4.2 3.9 4.8 5.0 3.7 5.3
##  [271] 2.8 5.0 3.6 5.2 3.5 4.6 4.5 4.9 5.4 3.2 3.9 3.9 3.1 4.9 5.4 5.7 2.8 4.9
##  [289] 4.3 4.3 5.6 4.7 4.0 4.4 2.4 5.1 4.8 6.2 3.9 5.3 4.6 1.9 6.5 3.2 4.0 6.1
##  [307] 3.6 5.5 4.8 5.8 5.3 3.0 5.0 5.5 4.7 7.3 4.0 4.0 4.7 6.1 6.0 4.9 5.4 4.7
##  [325] 4.2 6.2 3.5 4.0 5.9 5.6 6.1 5.7 4.1 4.6 4.4 6.1 5.0 4.7 3.7 5.0 4.5 4.1
##  [343] 4.8 3.9 6.1 5.0 5.7 3.4 6.0 4.1 4.5 5.4 5.6 4.7 5.9 6.0 5.7 5.2 5.2 4.7
##  [361] 4.6 5.6 4.8 4.0 4.8 5.9 4.4 2.6 3.7 3.6 5.0 4.4 3.4 4.4 4.0 5.8 4.4 3.8
##  [379] 5.1 5.9 2.7 4.6 3.9 5.3 5.8 2.4 3.7 3.8 6.0 4.4 5.0 5.5 2.9 5.9 5.2 4.6
##  [397] 5.3 3.4 4.2 5.0 4.4 3.6 5.5 4.5 6.7 4.4 4.9 4.3 4.7 5.3 4.2 4.7 5.6 4.8
##  [415] 4.8 4.0 5.2 5.9 4.4 3.9 4.7 4.6 5.6 4.1 4.2 3.3 4.7 4.5 5.8 4.5 4.8 4.8
##  [433] 4.5 4.5 3.9 5.1 5.5 4.3 4.0 5.4 3.9 4.6 5.7 2.7 4.5 4.8 4.5 4.7 3.7 4.6
##  [451] 5.1 3.9 4.3 3.5 5.6 3.5 6.2 3.5 3.7 6.2 4.1 5.0 4.6 3.4 5.4 5.7 5.2 2.7
##  [469] 4.5 4.4 5.6 5.7 5.2 3.1 2.7 2.7 3.6 4.4 4.2 4.3 4.3 1.6 4.0 5.4 4.6 3.0
##  [487] 5.6 2.4 5.3 4.1 3.6 3.8 6.0 4.7 4.8 5.5 6.5 4.8 2.3 6.4 5.3 5.1 5.3 4.9
##  [505] 4.7 5.0 3.9 3.4 3.2 4.7 4.7 5.9 3.9 5.8 4.1 4.8 5.8 5.2 3.9 4.1 4.2 3.6
##  [523] 4.3 4.2 6.1 5.2 3.9 4.5 2.6 4.6 3.0 4.2 5.2 5.3 3.2 4.2 4.0 4.6 5.5 5.4
##  [541] 5.7 4.7 3.4 4.8 4.9 2.7 5.7 5.5 5.6 3.8 4.1 4.1 3.2 3.4 6.1 4.7 4.4 4.3
##  [559] 3.8 3.6 2.3 4.2 5.6 3.5 3.8 2.9 4.2 5.3 3.9 2.7 4.6 4.1 3.5 4.8 3.6 4.2
##  [577] 4.4 5.2 4.4 4.8 5.4 3.9 4.6 4.6 2.8 4.6 4.6 2.3 3.5 4.3 4.7 3.6 4.2 3.8
##  [595] 4.5 4.8 4.5 3.7 5.0 4.1 3.3 4.7 4.0 4.2 5.6 5.5 3.6 6.2 3.6 5.0 4.8 4.3
##  [613] 3.2 4.7 4.6 4.2 5.2 4.7 3.4 4.4 5.2 5.1 4.3 4.5 3.2 4.9 3.7 4.0 6.0 5.4
##  [631] 3.7 4.6 5.8 5.2 4.2 4.9 5.5 3.7 4.0 2.4 6.6 5.1 4.5 4.1 4.8 4.3 3.3 4.9
##  [649] 4.3 4.5 4.3 3.2 3.4 5.0 3.8 4.2 4.0 3.6 3.6 3.0 5.2 4.3 3.1 6.0 4.2 3.2
##  [667] 3.7 4.9 5.5 5.1 5.0 5.4 3.8 5.4 3.8 4.4 4.0 3.8 3.2 2.9 4.4 5.5 3.4 4.3
##  [685] 4.9 5.2 3.7 2.0 5.1 4.7 5.0 4.2 4.6 3.2 5.6 5.3 3.0 3.3 6.2 5.8 4.4 3.7
##  [703] 5.0 7.2 4.4 5.0 3.8 6.7 4.3 3.9 3.0 4.2 4.8 5.0 4.0 4.2 4.3 5.0 4.3 4.8
##  [721] 4.4 4.8 4.8 5.1 3.9 4.0 4.9 5.2 3.4 4.7 4.7 3.6 3.2 3.9 3.4 4.0 5.5 4.1
##  [739] 3.4 4.6 4.3 3.7 3.3 4.7 4.1 3.0 5.7 5.6 6.3 6.1 4.8 5.3 4.5 3.2 3.2 3.3
##  [757] 4.6 5.7 2.8 4.3 4.3 5.1 5.0 5.1 5.3 5.9 3.6 3.6 5.3 3.7 4.6 5.2 5.5 4.2
##  [775] 4.2 5.6 5.9 6.5 4.0 4.4 4.7 6.1 3.9 5.1 5.4 4.9 5.3 4.7 4.4 5.0 3.0 3.2
##  [793] 4.8 2.8 3.6 4.5 3.9 5.0 6.6 3.5 3.9 4.2 5.3 4.8 3.8 3.8 3.5 2.0 3.2 3.6
##  [811] 5.8 6.8 3.3 2.4 4.3 3.8 3.9 3.4 5.0 5.1 4.1 3.4 4.7 5.1 4.6 4.9 3.8 5.1
##  [829] 4.0 5.5 3.4 5.6 4.9 6.0 5.5 4.8 4.4 3.3 1.7 4.3 4.2 3.1 5.5 5.3 4.3 3.1
##  [847] 5.5 4.0 4.1 5.0 4.3 4.6 5.6 5.8 4.6 6.0 4.9 3.9 6.2 4.4 4.6 3.9 4.6 4.2
##  [865] 5.0 4.1 5.8 4.8 3.3 4.0 4.1 5.2 4.5 5.2 3.9 4.1 5.0 4.5 4.4 4.9 4.2 3.6
##  [883] 4.4 4.6 4.3 4.6 3.7 4.8 2.8 4.8 4.9 5.9 5.8 3.9 4.1 4.4 5.1 4.2 3.7 6.1
##  [901] 5.5 3.7 5.5 4.5 4.9 3.9 5.0 4.4 4.3 3.7 4.6 4.7 3.7 4.0 5.2 3.3 4.4 3.2
##  [919] 3.3 4.2 4.2 5.3 5.5 4.2 4.4 4.3 5.0 5.5 3.1 5.0 4.0 2.8 4.6 3.1 4.9 4.4
##  [937] 5.6 4.5 3.9 4.2 4.0 6.1 4.5 4.0 5.3 5.4 3.4 3.8 3.2 3.2 5.9 4.1 4.0 4.5
##  [955] 4.0 4.8 4.5 5.3 4.6 3.3 5.1 2.2 5.3 3.8 4.1 5.4 1.5 4.3 4.2 4.8 5.4 4.8
##  [973] 4.5 5.9 4.3 4.2 4.8 3.5 5.5 5.6 4.4 3.0 5.1 5.3 4.0 4.5 5.2 6.1 5.5 6.0
##  [991] 3.8 4.7 4.4 5.2 5.2 4.2 4.6 4.5 4.5 2.7
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
## 2.7000 6.1025
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
##    [1] 2.8 5.4 3.7 4.3 4.0 5.0 5.8 3.9 5.6 4.6 5.5 6.4 5.2 5.2 4.9 4.5 4.9 5.5
##   [19] 6.1 5.1 4.6 3.8 4.9 6.3 7.4 3.6 5.4 3.8 3.2 4.2 4.0 4.9 4.7 4.2 4.6 3.3
##   [37] 4.9 5.5 4.7 5.7 3.9 4.8 4.6 4.5 5.5 4.7 5.3 5.0 4.9 6.0 6.0 3.2 5.8 3.4
##   [55] 4.1 4.3 4.6 4.3 4.4 4.8 4.6 4.1 4.6 4.6 2.6 4.0 5.1 3.1 4.0 4.1 4.0 3.4
##   [73] 4.3 5.0 6.5 5.2 4.5 2.7 3.8 3.6 2.8 5.6 5.0 4.6 5.1 3.2 3.9 4.8 5.0 4.3
##   [91] 4.5 3.3 4.9 3.4 5.2 4.3 3.6 5.3 2.9 5.3 5.0 3.2 4.1 4.9 5.6 4.5 5.3 4.4
##  [109] 4.9 3.8 3.1 4.1 4.0 4.7 3.7 3.1 4.0 4.9 4.4 5.8 3.2 5.2 3.8 3.0 5.0 3.7
##  [127] 4.1 3.5 5.3 5.4 5.6 4.8 3.7 3.8 4.7 2.3 4.7 4.1 3.3 3.5 4.1 6.2 4.2 3.7
##  [145] 4.7 4.1 3.8 4.3 6.3 2.9 4.1 2.2 4.1 4.9 4.9 4.4 2.6 4.5 4.8 4.2 5.9 3.5
##  [163] 4.3 3.7 4.0 3.9 3.9 4.1 4.8 4.2 2.3 4.5 4.4 4.4 4.8 4.6 5.3 4.4 4.0 4.0
##  [181] 4.1 3.8 4.6 4.3 4.6 4.1 4.8 4.8 6.3 3.8 5.1 3.6 5.0 4.6 4.6 4.8 3.3 5.1
##  [199] 5.2 4.5 3.5 4.6 4.2 4.4 4.3 5.3 5.0 5.1 4.5 5.1 5.8 3.4 2.8 4.7 3.0 2.5
##  [217] 6.4 4.6 5.2 4.4 2.7 4.0 4.7 4.4 4.1 3.6 5.3 4.4 3.1 4.0 2.1 3.4 5.2 6.0
##  [235] 6.0 2.4 4.3 3.7 4.5 4.4 3.2 4.9 3.9 5.1 4.3 5.0 5.4 2.1 4.7 3.6 4.2 3.9
##  [253] 4.1 3.8 4.7 4.5 6.0 3.4 6.2 4.2 6.0 2.2 4.2 4.2 4.5 3.4 4.0 4.8 3.9 4.2
##  [271] 6.3 4.6 6.8 3.5 4.9 7.1 3.4 6.2 3.7 4.8 4.9 5.3 3.6 5.3 4.0 5.6 4.8 3.1
##  [289] 5.3 4.2 4.2 5.7 4.1 3.6 2.5 4.5 3.7 5.9 6.1 5.0 4.4 3.8 4.5 6.0 2.3 4.8
##  [307] 6.4 4.4 6.2 3.7 3.9 3.4 5.8 5.1 5.3 3.2 3.5 2.9 4.1 5.8 5.0 5.6 5.0 4.2
##  [325] 3.9 6.4 3.8 4.5 5.4 4.8 4.9 5.0 6.0 4.8 5.5 4.5 4.2 5.1 4.4 5.2 4.8 4.0
##  [343] 3.7 5.4 3.4 4.6 5.6 3.6 5.9 2.9 4.2 4.7 4.3 5.1 4.4 5.0 5.9 4.2 5.5 4.2
##  [361] 4.2 4.4 4.0 3.2 3.9 4.4 6.7 4.3 3.7 5.5 4.2 3.9 3.8 4.9 5.5 4.4 3.4 4.5
##  [379] 4.7 5.7 3.8 6.2 3.3 4.5 4.1 3.5 5.1 4.6 4.1 4.0 5.1 2.6 2.7 3.6 3.8 6.6
##  [397] 3.7 5.9 5.0 5.3 5.2 4.8 5.3 4.9 4.3 5.3 3.7 4.6 3.3 5.0 3.8 4.3 4.9 5.9
##  [415] 5.1 4.5 5.2 3.8 2.8 3.5 5.9 2.1 3.1 3.8 3.0 5.4 4.3 4.7 3.9 6.0 5.6 3.9
##  [433] 5.6 4.8 5.0 5.8 3.2 5.5 4.1 4.8 4.1 4.1 4.4 3.5 5.3 5.1 5.1 4.3 5.0 4.6
##  [451] 5.6 5.0 4.0 4.3 4.0 3.8 4.2 4.6 5.0 4.4 4.4 4.2 4.7 4.8 5.0 4.0 5.6 5.7
##  [469] 4.9 5.7 3.4 4.5 5.1 3.9 5.7 3.8 2.6 4.8 4.8 4.3 3.4 2.5 5.4 4.6 5.7 6.0
##  [487] 5.6 3.4 4.5 4.4 3.6 4.5 3.3 4.2 5.9 3.6 5.5 3.7 4.3 4.9 2.3 5.9 3.6 3.5
##  [505] 5.1 5.4 5.2 5.9 4.6 2.9 4.5 4.2 4.0 3.4 5.7 3.5 4.9 5.7 4.0 4.8 5.3 4.9
##  [523] 4.7 5.6 5.2 4.5 4.4 4.4 5.1 4.8 5.1 3.4 5.5 5.6 6.1 4.0 3.5 4.5 5.9 3.4
##  [541] 5.5 4.0 3.5 5.2 5.9 3.4 3.5 4.8 4.8 3.9 4.5 4.2 3.1 4.2 3.0 3.7 6.6 5.6
##  [559] 6.3 4.3 4.2 5.4 4.1 3.6 4.2 4.4 3.9 6.0 4.1 3.0 3.6 3.9 4.7 3.9 5.2 4.9
##  [577] 4.1 4.7 5.9 3.3 5.1 4.4 3.6 4.5 5.0 3.7 5.3 4.7 4.9 3.0 3.4 3.6 4.6 2.5
##  [595] 5.5 5.0 4.7 3.2 5.5 2.5 3.5 5.4 4.9 4.4 4.0 4.9 6.1 5.2 4.2 4.2 4.8 5.2
##  [613] 5.5 4.0 4.9 3.4 6.2 4.9 3.7 3.7 4.0 4.7 4.6 4.6 5.1 4.3 5.0 6.0 3.3 5.5
##  [631] 3.3 4.3 5.2 4.3 4.0 5.0 4.7 5.8 3.5 5.2 2.9 3.1 5.9 6.1 5.9 4.0 3.6 2.9
##  [649] 5.7 5.2 3.5 4.0 6.4 5.5 4.2 4.4 3.2 4.5 5.0 4.8 5.1 6.0 4.5 4.8 4.7 4.4
##  [667] 4.5 3.2 5.0 5.2 4.6 3.5 4.0 3.3 4.9 6.3 4.0 4.3 2.8 4.5 2.9 5.2 2.9 5.1
##  [685] 4.4 5.3 3.3 5.3 6.7 5.3 4.7 4.2 4.0 4.9 4.7 3.3 4.1 4.7 3.4 3.9 4.8 5.1
##  [703] 4.1 4.9 5.4 4.7 5.1 3.1 4.6 5.4 4.0 4.2 5.6 5.2 4.2 5.2 4.1 4.5 5.4 3.6
##  [721] 6.1 4.3 5.1 3.3 5.3 4.1 5.8 4.7 3.5 4.9 4.7 4.8 4.1 3.5 4.8 4.5 4.5 5.1
##  [739] 4.9 6.3 3.9 5.2 4.2 4.7 5.8 4.2 4.4 2.8 5.5 4.1 3.7 2.6 4.4 3.1 4.0 4.1
##  [757] 5.2 3.8 4.7 3.8 4.1 3.6 5.5 4.1 3.9 6.5 5.0 5.0 4.4 5.7 4.1 4.4 5.7 5.2
##  [775] 3.8 5.1 3.9 4.5 4.8 3.0 3.4 4.0 5.1 3.9 3.7 4.3 3.3 3.6 4.8 5.4 4.9 4.6
##  [793] 3.7 4.4 3.1 4.1 2.7 5.8 4.8 5.3 3.4 4.6 5.4 5.2 4.4 5.1 5.5 4.2 4.7 3.6
##  [811] 5.4 6.2 4.6 5.5 6.0 5.1 5.9 6.5 4.7 3.3 5.8 4.2 4.2 4.3 4.0 5.2 4.0 6.4
##  [829] 5.7 4.0 5.8 5.6 3.5 5.2 5.5 2.8 3.8 3.8 4.3 3.2 4.8 4.2 5.2 6.4 5.2 4.9
##  [847] 4.4 4.9 4.7 4.5 5.3 5.1 3.5 4.8 6.1 3.7 4.0 4.6 4.9 5.6 4.7 5.0 4.6 5.9
##  [865] 2.4 4.4 4.8 4.3 4.4 2.3 4.7 4.4 6.0 4.9 3.8 4.5 4.3 4.4 4.3 4.3 4.1 2.7
##  [883] 4.6 5.5 5.3 2.7 5.8 5.6 3.2 4.1 4.0 5.3 4.2 4.7 3.8 3.3 4.8 4.6 3.6 3.8
##  [901] 4.4 3.1 4.9 3.9 2.5 5.5 3.5 2.2 3.8 4.8 4.4 4.6 4.4 3.5 3.0 5.1 3.4 3.9
##  [919] 4.2 6.3 2.8 4.5 5.8 5.3 4.9 4.1 3.5 3.7 4.4 4.4 4.3 3.7 3.7 4.9 6.0 4.8
##  [937] 4.9 3.3 4.7 4.5 4.2 3.2 4.1 5.6 5.7 2.6 4.9 4.4 1.9 4.4 3.5 4.0 4.2 5.3
##  [955] 4.1 2.5 4.5 4.5 4.2 4.8 3.9 4.7 5.0 3.7 4.1 5.6 5.6 5.5 4.0 4.3 4.7 5.7
##  [973] 3.9 3.7 5.0 3.1 4.8 5.1 5.7 5.5 4.1 4.5 5.1 4.6 5.5 4.2 3.2 3.1 4.5 4.8
##  [991] 5.5 5.4 4.2 4.5 3.7 3.9 4.9 4.0 5.2 4.9
## 
## $func.thetastar
## [1] -0.0145
## 
## $jack.boot.val
##  [1]  0.49194030  0.33360882  0.20563003  0.13788301  0.07500000 -0.09819277
##  [7] -0.21008646 -0.32924791 -0.39824561 -0.48870523
## 
## $jack.boot.se
## [1] 0.9352939
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
##    [1] 3.7 4.8 3.2 2.8 4.5 3.4 4.4 4.2 4.1 4.1 4.6 3.7 2.9 4.3 2.8 5.1 3.8 5.4
##   [19] 4.5 5.4 4.6 3.9 3.0 3.7 5.5 4.7 4.4 4.5 4.6 6.0 2.7 4.8 2.6 3.1 5.8 3.9
##   [37] 4.6 5.4 4.2 4.2 5.5 5.3 5.5 5.0 4.2 4.7 4.4 5.2 5.0 4.8 3.8 5.2 5.6 4.3
##   [55] 3.8 5.3 4.8 3.1 5.2 5.7 5.0 2.7 4.7 3.7 4.1 5.5 3.0 4.0 3.4 4.8 4.8 4.9
##   [73] 5.3 4.5 5.4 4.7 6.4 4.4 3.4 4.9 3.5 4.2 5.1 4.8 6.4 3.8 5.6 4.7 3.9 7.2
##   [91] 3.8 4.4 3.5 4.5 5.4 3.5 5.0 4.3 3.7 3.8 5.1 3.7 5.0 3.2 4.3 4.6 4.2 4.6
##  [109] 4.3 4.2 3.3 2.2 2.1 4.5 5.0 3.9 3.5 4.8 4.8 5.2 4.6 4.7 4.1 4.3 3.0 5.3
##  [127] 3.9 4.3 4.4 4.2 2.7 5.7 4.2 4.3 3.9 3.2 3.8 4.0 3.0 4.7 4.0 5.2 6.0 3.7
##  [145] 6.0 4.0 4.3 4.2 3.7 5.8 3.7 5.8 3.7 6.0 6.0 4.4 4.9 3.7 3.9 4.0 5.7 3.2
##  [163] 3.0 3.4 3.7 3.7 5.1 5.1 6.8 3.5 4.1 3.9 3.6 3.6 4.9 4.4 3.8 3.4 5.3 4.7
##  [181] 5.5 4.2 4.9 3.7 3.9 5.9 4.4 3.3 4.6 2.9 5.4 4.5 3.1 5.2 4.2 4.3 5.0 4.7
##  [199] 2.7 5.1 5.2 4.2 4.2 4.0 3.7 4.9 4.9 4.4 4.1 5.5 4.8 3.5 5.0 5.3 4.4 4.7
##  [217] 4.3 3.7 6.7 5.2 4.1 6.2 5.8 4.7 5.3 4.5 5.3 4.9 4.4 4.7 4.3 6.3 4.7 4.0
##  [235] 5.7 5.7 5.5 4.9 3.1 5.9 4.5 4.5 3.9 3.9 3.3 5.6 4.1 3.4 4.7 4.3 3.0 4.3
##  [253] 4.9 5.8 4.7 4.3 2.6 4.9 5.8 3.3 3.4 4.7 3.7 4.0 3.5 3.7 3.8 5.5 3.9 4.6
##  [271] 3.4 3.3 4.0 4.1 4.7 4.3 4.5 3.5 5.6 5.5 3.5 4.8 4.7 5.9 3.0 5.4 4.2 5.4
##  [289] 5.6 5.7 4.8 3.7 4.9 3.8 3.7 4.8 5.3 4.9 3.5 4.6 5.1 4.6 3.5 4.7 4.3 4.1
##  [307] 2.6 5.2 4.0 4.9 4.6 5.8 6.6 5.7 3.6 5.7 5.0 3.4 4.9 6.0 4.9 5.8 4.3 4.4
##  [325] 4.8 4.4 3.9 3.4 6.2 2.6 5.5 5.4 4.7 3.1 4.3 2.8 4.4 4.9 4.8 6.0 3.9 5.2
##  [343] 5.5 3.4 4.2 4.2 3.5 4.4 5.5 3.7 5.2 5.2 4.6 4.9 4.3 1.9 4.8 5.2 4.9 4.4
##  [361] 4.5 2.7 4.0 4.4 4.3 5.3 6.2 2.9 4.7 3.8 4.2 3.2 3.9 3.9 3.8 5.6 3.9 3.9
##  [379] 4.8 5.1 2.4 6.3 3.8 4.9 5.3 3.8 3.4 4.4 4.8 5.5 2.8 4.4 5.0 3.3 5.3 4.0
##  [397] 3.8 3.4 3.7 6.0 4.7 4.0 5.6 5.2 4.5 6.6 4.7 5.0 4.9 5.6 3.7 5.7 5.1 5.2
##  [415] 4.8 4.0 4.6 4.0 4.6 4.3 3.3 5.6 4.1 5.7 4.4 5.0 5.4 3.2 4.3 4.5 4.9 3.7
##  [433] 5.0 5.1 4.0 5.8 5.3 4.2 5.4 5.8 4.2 3.5 3.9 4.3 5.0 5.5 3.0 4.6 4.2 5.5
##  [451] 4.2 5.0 2.9 4.5 5.9 4.2 4.6 5.4 4.0 3.2 4.0 4.0 6.1 3.1 4.1 3.9 4.8 4.1
##  [469] 4.1 2.8 4.2 5.4 4.2 4.6 7.1 3.4 3.6 4.2 3.7 4.2 4.8 4.3 5.7 4.7 2.1 4.6
##  [487] 4.5 5.5 5.0 3.1 4.5 4.4 6.2 5.8 4.5 5.8 3.8 3.0 4.8 4.7 1.1 3.6 3.3 3.7
##  [505] 4.9 3.4 5.8 4.7 5.1 5.8 4.4 3.0 4.0 3.3 4.6 2.6 5.4 5.0 3.2 4.1 4.6 5.3
##  [523] 3.6 4.6 5.7 3.0 2.8 5.8 5.2 4.1 4.9 6.2 4.7 3.6 5.3 3.3 3.0 4.3 4.0 5.4
##  [541] 5.1 4.1 3.5 4.3 4.6 4.5 3.6 4.5 3.3 5.4 5.8 3.7 3.6 4.6 2.7 4.8 3.4 3.6
##  [559] 5.3 3.6 3.8 4.7 6.0 6.5 5.2 4.4 4.1 2.6 4.4 6.2 3.4 4.2 4.7 4.2 5.0 5.4
##  [577] 4.4 4.6 4.3 4.9 5.3 3.0 5.0 4.6 4.6 4.9 3.9 3.8 3.8 3.6 4.8 4.4 4.9 3.4
##  [595] 6.7 4.0 3.0 4.8 4.5 5.3 5.7 3.7 5.8 4.9 5.5 5.5 4.8 2.0 5.3 3.8 5.5 4.0
##  [613] 3.3 3.5 5.4 5.1 4.9 4.9 3.4 5.4 3.7 2.7 4.2 4.0 4.6 4.9 3.4 4.9 2.3 4.3
##  [631] 5.5 5.9 3.4 2.8 4.7 3.0 4.3 4.9 3.5 2.8 4.2 5.0 4.5 4.3 4.8 4.3 5.4 5.1
##  [649] 6.0 5.5 6.0 4.0 4.4 6.2 3.8 6.0 6.0 6.0 4.4 3.9 5.1 3.7 5.1 4.8 5.5 5.0
##  [667] 4.3 2.8 5.3 4.9 3.8 6.4 4.0 5.3 4.4 5.2 2.3 4.2 3.7 2.1 5.6 3.3 4.0 4.3
##  [685] 5.2 4.6 4.4 4.4 3.9 3.1 4.0 4.1 4.3 5.2 4.0 2.6 5.6 4.8 3.1 3.0 6.1 5.1
##  [703] 5.2 3.9 4.2 5.2 3.2 4.6 3.6 5.0 4.6 4.2 3.9 6.3 4.4 4.8 2.9 5.5 5.4 3.5
##  [721] 2.4 5.7 4.3 4.8 6.0 4.2 5.5 4.3 4.3 3.9 4.3 4.1 5.7 4.2 5.7 3.8 4.4 4.0
##  [739] 4.5 4.5 3.5 5.7 1.9 5.1 5.8 3.9 4.1 5.1 4.3 2.5 5.1 3.6 3.9 3.7 6.0 4.7
##  [757] 4.3 3.2 4.6 5.0 4.8 6.7 4.6 4.7 5.6 4.3 5.9 4.2 3.1 5.3 4.0 4.5 4.1 3.4
##  [775] 3.3 4.3 5.0 5.8 4.7 4.0 4.0 4.8 4.4 2.6 5.5 4.3 2.9 4.2 3.8 6.1 3.8 4.8
##  [793] 5.0 5.1 3.9 4.4 4.8 4.2 4.0 4.6 4.3 1.8 5.0 5.2 6.8 6.5 4.7 5.3 3.8 5.4
##  [811] 3.7 5.5 5.2 4.4 5.2 4.1 3.7 4.5 3.3 5.6 4.5 4.7 2.6 2.8 5.0 3.4 4.9 3.4
##  [829] 4.9 3.5 5.7 5.2 6.2 5.0 3.9 4.7 5.8 4.9 4.3 4.0 3.2 3.3 3.8 4.1 3.7 4.8
##  [847] 6.8 5.4 3.7 2.4 3.3 5.6 4.2 5.5 4.2 4.0 4.5 4.9 5.5 5.9 3.6 3.7 5.0 3.0
##  [865] 3.9 3.5 3.4 3.9 4.0 4.5 5.6 3.6 4.8 3.1 4.0 4.9 3.8 6.6 4.5 4.2 5.0 5.4
##  [883] 4.7 7.1 4.7 4.4 2.4 4.2 4.1 4.7 5.8 2.9 4.2 4.4 4.4 5.4 4.5 5.2 4.1 5.2
##  [901] 4.3 4.1 3.5 4.0 3.8 3.9 3.7 3.1 4.4 5.0 3.5 5.7 4.2 3.9 5.0 3.8 4.6 3.5
##  [919] 5.0 4.0 5.6 4.0 5.7 5.1 4.7 5.4 5.6 5.7 3.8 5.4 5.7 2.6 5.0 3.8 4.1 5.2
##  [937] 3.4 3.0 5.6 3.3 4.6 3.8 5.3 4.5 3.3 4.2 4.6 4.2 4.7 4.0 5.6 3.9 4.4 4.6
##  [955] 3.5 4.9 4.3 4.6 4.4 4.6 4.7 4.5 5.1 4.7 4.5 3.3 5.1 3.8 4.0 5.8 3.7 3.2
##  [973] 4.0 4.8 4.5 4.6 5.9 6.2 6.6 4.4 4.5 4.7 4.8 4.0 4.5 6.2 5.2 4.8 5.2 5.2
##  [991] 4.4 4.3 3.2 3.9 5.2 3.8 6.1 4.4 4.9 4.9
## 
## $func.thetastar
## 72% 
##   5 
## 
## $jack.boot.val
##  [1] 5.4 5.4 5.4 5.2 5.1 5.0 4.9 4.7 4.5 4.5
## 
## $jack.boot.se
## [1] 1.008018
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
## [1] 0.9089665
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
##    5.597422   10.958019 
##  ( 2.432306) ( 4.981871)
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
## [1]  0.29906639  0.52160360 -0.09717057  0.85513559 -0.75001429  0.06405796
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
##    [1]  9.179039e-01  4.347243e-01  1.355284e+00  9.408993e-01  8.217611e-01
##    [6]  5.925574e-01  2.782268e-01  1.049061e+00  5.993093e-01  2.036545e+00
##   [11]  4.809306e-01  4.533214e-01  1.021200e+00  1.288039e+00  3.729376e-01
##   [16]  1.233074e+00  8.951708e-01  1.152963e+00  3.014339e-01  1.014096e+00
##   [21] -4.133608e-01  8.615652e-01  1.276394e+00  1.358623e+00 -8.729185e-01
##   [26]  8.104672e-01  3.091149e-01 -1.458462e-01  5.161919e-01 -1.341441e-01
##   [31]  5.897757e-02  4.396905e-01  6.137496e-01  7.971813e-01  8.804276e-01
##   [36] -1.212238e-01  1.837094e-01  5.218858e-01  3.661837e-01  3.421855e-01
##   [41]  1.077673e+00  3.963616e-01  6.890303e-01  3.705559e-01  4.252795e-01
##   [46]  1.376069e+00  6.925985e-01  5.074603e-01 -1.847905e-01  8.503210e-01
##   [51]  1.023550e+00  5.773460e-01  7.943339e-01 -2.022772e-01  1.425803e-01
##   [56]  1.143655e+00  1.326660e+00  5.972537e-02  1.586853e-02  5.505672e-01
##   [61]  7.938108e-01  9.340292e-01  8.523750e-01  9.201893e-01  7.212883e-01
##   [66]  8.571981e-01  1.319512e-02  1.233184e+00  2.757919e-01  9.304039e-01
##   [71]  8.740518e-01  8.578629e-01  3.117282e-01  9.154940e-01 -6.425069e-01
##   [76]  9.189599e-02  1.406099e+00  3.914063e-01  4.404474e-01  7.438436e-01
##   [81]  5.676189e-01  3.322250e-01  8.433372e-01  7.797560e-01  9.491664e-01
##   [86]  7.704270e-01  6.211501e-01  5.255232e-01  5.196636e-01  6.812765e-01
##   [91]  1.534642e+00  2.092613e-01  1.865356e+00  1.534739e-01  5.661652e-01
##   [96]  5.467682e-01  2.217007e-01  8.531283e-01  1.895196e+00  9.333815e-01
##  [101]  4.594655e-01 -8.622931e-02  5.778032e-01  1.365225e+00  8.739746e-01
##  [106] -4.019606e-01  9.311746e-01  4.323731e-01  8.204427e-01  2.842037e-01
##  [111]  1.150446e+00  6.736799e-01  6.877898e-01  7.938108e-01  1.004516e+00
##  [116]  1.054356e+00  5.501928e-01  9.486720e-02  1.563344e+00  3.992503e-01
##  [121]  5.379076e-01  1.213234e+00  1.002442e+00  4.828313e-01  4.763799e-01
##  [126]  1.198382e+00  6.885796e-01  6.905981e-01  9.161970e-01  1.100900e+00
##  [131]  4.542958e-01  3.865473e-01  5.597376e-01  5.485850e-01  1.929965e+00
##  [136]  1.307077e+00  6.090905e-02  6.132433e-01  7.581455e-01  8.784919e-01
##  [141]  8.019751e-01  1.077080e+00  1.192746e+00  8.062463e-01  1.559219e-01
##  [146]  1.842440e+00 -1.785773e-01  6.440654e-01  6.696743e-01  8.457882e-01
##  [151] -1.084184e-01  6.754039e-01  5.465461e-01 -1.243491e-01  7.777982e-02
##  [156]  2.067113e+00  1.918481e-01  1.453358e+00  1.120596e+00  6.136166e-01
##  [161]  1.520217e+00  7.777818e-01  1.130251e+00  1.524247e+00  3.737884e-01
##  [166]  9.316481e-01  5.242807e-01  8.626503e-01  4.594909e-01  5.972370e-02
##  [171]  1.035721e+00  1.063247e+00  1.245888e+00  1.137453e+00  5.888300e-01
##  [176]  1.883012e-01  6.631341e-01  4.647535e-01  9.800935e-01  5.681755e-01
##  [181]  1.779027e-01  3.561134e-01  1.893395e+00 -2.196679e-02  9.858713e-01
##  [186]  1.490207e+00  5.530169e-01  7.977691e-01  1.094963e+00  3.304232e-01
##  [191]  1.120898e+00  2.807911e-01  4.435209e-01  5.186569e-01 -1.050749e+00
##  [196]  9.527298e-01  4.743375e-01  1.897547e+00  4.422524e-01 -2.305072e-01
##  [201]  1.814817e+00  1.109993e+00  1.599582e-01  3.976573e-01  1.438794e+00
##  [206]  4.241547e-01 -1.480895e-01  1.086631e+00  1.026594e+00  1.465881e+00
##  [211]  9.949373e-01  4.694355e-01  5.032186e-01  8.047940e-01  5.950164e-01
##  [216]  3.991341e-01  1.201391e+00  3.265019e-01  3.625002e-01  2.036644e-01
##  [221]  8.176595e-01  1.004516e+00 -3.387368e-01  1.528693e+00  3.314152e-01
##  [226]  6.895236e-01  9.526611e-01  5.653253e-01  1.097154e+00 -3.263469e-02
##  [231]  8.402141e-01  1.106524e+00  1.023074e-01  7.718036e-01  1.295159e+00
##  [236]  2.491897e-01  8.321054e-01  5.997248e-01  1.475177e+00  5.141912e-01
##  [241]  1.227223e+00  3.511065e-01  2.075506e-01  4.479523e-01  8.514967e-01
##  [246]  3.837260e-01  1.216661e-01  1.352718e+00  1.785647e+00 -1.189352e-01
##  [251]  6.700877e-01  1.297434e+00  1.227858e+00  7.086122e-01  8.189058e-01
##  [256]  1.042028e+00  7.878825e-01  8.256668e-01  2.574627e-01  1.089692e+00
##  [261]  1.346616e+00  8.821480e-03  5.593684e-01  1.842377e+00 -5.355373e-01
##  [266]  1.395191e+00  9.736027e-01  1.397345e+00  5.657987e-01  3.485696e-01
##  [271]  6.801124e-01 -7.785794e-02  7.637993e-01  5.308829e-01  2.630682e+00
##  [276]  6.348676e-01  4.775080e-01  2.731690e-01  7.863645e-01  9.126825e-01
##  [281] -4.552531e-02 -2.761052e-02  3.750011e-01  9.949373e-01  6.700877e-01
##  [286]  1.027312e+00  9.242044e-01  9.356253e-01  4.989856e-01  1.426716e+00
##  [291] -2.059031e-01  1.542619e+00  9.019678e-01  6.305397e-01  8.916811e-01
##  [296]  1.082944e+00  1.255143e+00  2.613993e-01  6.058006e-01  4.639479e-01
##  [301]  2.819611e-02  1.577774e+00  1.252361e+00  8.798308e-01  1.309839e-05
##  [306]  1.105211e+00  3.867886e-01  1.069353e+00  1.267574e-02 -4.111519e-01
##  [311] -8.692028e-02  1.238269e+00  9.530472e-01  2.146092e-01  1.091724e+00
##  [316]  1.166096e+00  1.018933e+00  8.372751e-01  6.468628e-01  1.284491e+00
##  [321]  7.146733e-01 -1.375524e+00  2.630952e+00  5.250742e-01  7.777982e-02
##  [326]  9.078727e-01  7.780036e-01  1.065573e+00  9.552726e-01  1.570704e+00
##  [331]  7.027465e-01  9.857221e-01  9.476513e-01  1.238056e+00  4.783442e-01
##  [336] -2.336627e-02  9.014885e-01  5.922625e-01  1.162305e+00  1.349401e-01
##  [341]  3.881624e-01  8.155887e-01  1.868206e+00  7.470401e-01 -1.113147e-01
##  [346]  5.394876e-01  1.824581e-01  1.772502e+00  5.119038e-01  2.745794e-01
##  [351] -2.738801e-01  8.858374e-01  4.442133e-01 -6.425825e-01  1.339140e+00
##  [356]  1.068924e+00  6.047231e-01  5.215232e-01  7.418086e-01  4.224646e-01
##  [361]  7.636336e-01  1.002559e+00  4.569444e-01  1.083922e+00  3.469997e-01
##  [366]  2.509911e-01  3.754596e-01  1.107748e+00  7.120944e-02  9.543335e-01
##  [371]  9.303530e-01  1.728747e-01  1.343827e+00 -8.329624e-02 -2.022881e-01
##  [376]  1.305979e+00  6.179483e-01  1.757586e-01  1.476220e-01  4.050240e-01
##  [381]  3.622936e-01  1.060900e+00  6.985863e-01 -5.657260e-01  5.801133e-01
##  [386]  3.387552e-01  1.456749e+00  4.605343e-01 -3.386604e-01  6.620393e-01
##  [391]  3.888438e-01  4.480889e-01  1.138523e+00  1.290272e+00  4.683372e-01
##  [396]  8.583661e-01  1.008151e+00  3.988188e-01  1.712078e+00  1.883804e-01
##  [401]  1.225009e+00  7.387806e-01  1.447322e+00  7.975806e-01  2.534226e+00
##  [406]  7.441003e-01  7.329218e-01  5.753394e-01  1.058995e+00  1.785029e+00
##  [411]  7.151866e-01  1.368977e+00  7.497092e-01  2.340676e+00 -9.416270e-01
##  [416]  1.178782e+00  1.125563e-01  9.948683e-01  2.068203e-01  1.041196e+00
##  [421]  3.014932e-01 -4.046415e-01  1.033525e+00  5.439480e-01  1.070812e+00
##  [426]  1.090063e+00  9.098033e-01  4.565500e-01  7.667904e-01  6.704375e-01
##  [431]  7.744661e-01  9.958884e-01  6.066423e-01  8.023350e-01  1.219692e+00
##  [436]  5.653501e-01  2.522827e+00  6.821011e-01  1.465424e-01  9.381209e-01
##  [441]  6.038579e-01  2.807911e-01  6.820113e-01  1.105418e-01  6.345576e-01
##  [446]  1.038044e+00  4.574428e-01  6.709974e-01 -4.327717e-02  4.239269e-01
##  [451]  4.940357e-01  9.375174e-01  1.890446e+00  7.643152e-01  4.838517e-01
##  [456]  6.072295e-01  7.679311e-03 -3.551205e-01  1.062507e+00  4.482023e-01
##  [461]  1.149307e+00  1.202521e+00 -2.173225e-04  1.082094e+00  4.702345e-01
##  [466]  3.904538e-01  1.238432e+00  2.842524e-01  1.498697e+00  1.196893e+00
##  [471]  3.220753e-01  8.612902e-02  1.405373e+00  1.565746e+00  4.426745e-01
##  [476] -1.830964e-01  8.994284e-01  1.212880e+00  4.262326e-01  4.871277e-01
##  [481]  5.380422e-01  1.320228e+00  2.420188e-01  9.482916e-01  6.407819e-01
##  [486]  2.563779e+00  1.475638e+00  5.434892e-01  1.127982e+00  6.631341e-01
##  [491] -2.317266e-01  1.094239e+00  1.674689e+00  3.700098e-01  5.280209e-01
##  [496]  1.697051e+00  1.351299e+00  7.697831e-01  7.668780e-01  1.243832e+00
##  [501]  1.822656e-01  1.844745e+00  2.619287e+00  1.101379e+00  8.573311e-01
##  [506] -3.116609e-02  5.221374e-01  1.054283e-01  9.828078e-01  1.234657e+00
##  [511]  1.193239e+00  9.289765e-01  6.713434e-01  2.090978e+00 -5.689046e-02
##  [516]  9.008053e-01  8.738696e-01  7.182800e-01  1.063254e+00  6.132785e-01
##  [521]  8.822295e-01  7.942898e-01  2.698701e-01  1.128583e-01  5.974044e-01
##  [526]  4.716760e-01  8.310937e-01  6.441603e-01  1.140292e+00  5.778002e-01
##  [531]  5.486858e-01  5.721985e-01  1.830187e-01 -5.757927e-02 -3.473650e-01
##  [536]  1.164220e-01  7.617967e-01  3.058777e-01  8.155887e-01  1.537026e+00
##  [541]  9.815721e-01 -9.119238e-02  7.313536e-01  9.421906e-01  7.062059e-01
##  [546]  3.975127e-01  2.016139e-01 -1.080742e-01  2.037192e+00  1.334352e+00
##  [551]  4.696416e-01  8.628649e-01  9.949373e-01  5.110040e-01  1.582816e+00
##  [556]  6.982819e-01  8.799692e-01  2.099524e-01  9.572886e-01  1.186481e+00
##  [561]  1.123204e+00  1.365028e+00  3.281564e-01  1.528237e+00  1.061561e+00
##  [566]  1.470233e+00 -3.971402e-03  8.574899e-01  5.528495e-01  3.765969e-01
##  [571]  8.775557e-01  2.067217e+00  8.976679e-01  2.722867e-01  5.631614e-01
##  [576]  9.114099e-01  1.650306e+00  1.468549e+00 -9.929316e-02  5.336999e-01
##  [581]  1.788728e-01  7.588868e-01  8.221856e-02  5.572749e-01  6.239387e-01
##  [586]  1.285708e+00  1.361370e+00  7.763053e-01  1.037355e-01  4.328220e-01
##  [591] -4.160273e-02  8.561885e-01 -5.258768e-01  1.398830e+00 -8.147225e-02
##  [596]  1.318277e+00  2.427624e+00  5.458258e-01  2.368609e-01  6.927973e-01
##  [601]  5.165620e-01  7.296031e-01  1.522350e+00  2.392866e-01  2.552150e+00
##  [606]  6.752608e-01 -2.979178e-01  8.131599e-01 -1.287598e-01  1.245839e+00
##  [611]  9.014885e-01  1.440089e+00  8.320355e-01  4.849788e-01  5.271807e-01
##  [616]  1.878999e-01  9.378042e-01  2.110846e-01  8.613750e-01  2.296524e+00
##  [621]  1.450831e+00  1.612915e+00  1.762661e-01  1.214884e+00  1.431276e+00
##  [626]  3.658518e-01  1.508491e+00  8.145558e-01  1.016239e+00  3.920743e-01
##  [631]  9.262157e-01  1.678312e+00  6.860230e-01  4.558159e-01  5.361607e-01
##  [636]  1.552175e+00  1.615516e+00  6.288076e-01  2.344004e-01  5.903921e-01
##  [641]  4.542592e-01  9.906643e-01  1.638773e-01  5.900841e-01 -2.183577e-01
##  [646] -4.941658e-01  8.489445e-01  1.708725e+00  1.426487e+00  4.173214e-01
##  [651]  8.130203e-01  1.409477e+00  1.713801e+00  1.150930e+00  3.621805e-01
##  [656]  1.488906e+00 -8.904906e-01  8.447155e-02  1.154576e+00 -2.783415e-03
##  [661]  8.602525e-01  4.613801e-01  3.970455e-01  1.087320e+00  2.776045e-01
##  [666]  3.010494e-01  5.258797e-01  8.084641e-01 -2.957843e-01  2.239533e-01
##  [671]  7.720418e-01  4.996707e-01  7.934049e-01  1.484236e+00  1.518291e-01
##  [676]  1.165348e+00  9.124654e-01  8.421827e-01  3.322208e-01  5.857680e-01
##  [681]  1.842669e-01  3.256773e-01  5.446105e-01  8.172159e-01  2.653549e+00
##  [686]  6.401301e-01  6.190586e-01  5.407450e-01  1.361803e+00  1.340738e+00
##  [691]  3.879008e-01  9.317089e-01  1.055140e+00  1.498016e+00  8.324300e-01
##  [696]  8.443434e-01  2.470896e+00 -1.692227e-01  6.180999e-01  8.712623e-01
##  [701]  6.662604e-01  1.382828e+00  1.019677e+00  8.324402e-01  2.081026e+00
##  [706]  9.486812e-01  1.296559e-01  3.023127e-01  1.040843e+00  2.660844e-01
##  [711]  4.776587e-01  5.778873e-01  7.538618e-01  1.387860e-01  8.145042e-01
##  [716]  7.248390e-01  2.189781e-01  1.999940e-01  1.617407e+00  2.751419e-01
##  [721]  1.759302e+00  1.826624e+00  1.260905e+00  1.062521e+00  1.772520e+00
##  [726] -1.397337e-01  9.482411e-01  8.327526e-01  9.033590e-01  2.568247e-01
##  [731]  1.890041e+00  1.480627e-01  9.565110e-01 -2.349889e-01  1.733877e+00
##  [736]  1.341032e+00  5.305652e-01  1.598188e+00  8.348523e-01  9.850681e-02
##  [741]  1.054840e+00  3.537176e-01  1.311372e+00  1.010593e+00  5.060186e-01
##  [746]  2.426130e-01  8.164607e-01  4.160017e-01  3.067570e-01  1.084616e+00
##  [751]  9.882891e-01  3.724238e-01  6.620102e-01 -7.140180e-01  7.300753e-01
##  [756]  1.577219e-01  9.054571e-01  8.672137e-01  6.003078e-01  8.881439e-02
##  [761]  1.450714e+00  7.371939e-01  8.813050e-01  1.370073e+00  1.754762e-01
##  [766]  1.516500e+00  8.524729e-01  5.664393e-01  5.195485e-01  8.296365e-01
##  [771]  4.404740e-01  4.208567e-01  1.210914e+00  3.455144e-01 -1.313571e-01
##  [776]  7.782179e-01  4.110873e-01  2.001704e-01  2.785858e-01  1.046999e+00
##  [781]  7.328851e-01  6.063280e-01  1.592812e+00  6.703492e-01  2.086670e-01
##  [786] -3.252186e-03  1.516171e+00 -2.092808e-01 -7.261289e-01  8.408057e-01
##  [791]  3.961728e-01  6.958104e-01  1.060763e+00  1.311066e+00  1.896431e+00
##  [796]  1.444470e+00  1.126292e+00  1.005812e+00  9.089665e-01  9.702208e-01
##  [801]  1.049837e+00  5.696017e-01  1.000784e+00  9.105912e-01 -1.118031e-01
##  [806]  5.614127e-01  6.406076e-01  3.517861e-01  8.449405e-01  2.333463e+00
##  [811]  9.831010e-01  3.937601e-01  9.328352e-01  1.518073e+00  6.957264e-01
##  [816]  4.509429e-01 -3.945193e-02  2.634069e+00  3.276902e-01 -2.402219e-01
##  [821]  4.553346e-01  9.297104e-01  8.211335e-02  3.149208e-01  1.241050e+00
##  [826]  3.177552e-01  1.027797e+00  1.401418e+00  4.076431e-01  3.956257e-01
##  [831]  6.198874e-01  1.320609e+00  1.194288e+00  5.863463e-01  9.269156e-01
##  [836] -2.713913e-01  1.022939e+00  2.651537e-01  9.548970e-01 -9.922029e-01
##  [841]  5.470663e-01  5.115338e-01  1.126019e+00  1.542757e+00  3.514893e-01
##  [846]  4.023875e-01  7.594120e-01  8.194001e-01  9.418733e-01  6.523458e-01
##  [851]  6.446409e-01  1.609559e+00  8.606627e-01  6.655208e-01  5.179848e-01
##  [856]  1.338429e+00  9.695294e-01  2.332978e+00  4.712152e-01  1.894127e+00
##  [861]  1.059799e-01  1.137468e+00  9.575959e-01 -3.364304e-02  7.557245e-01
##  [866] -6.056080e-01  2.508311e+00  9.494881e-01  3.937601e-01  5.258797e-01
##  [871]  4.044846e-01  7.468026e-01  7.051570e-01 -2.478348e-01  5.083349e-01
##  [876]  1.223320e+00  5.195378e-02  2.626924e-01  3.510718e-01  5.834305e-01
##  [881] -7.664764e-02  4.836246e-01  1.862288e-01  1.373280e+00  1.929592e-01
##  [886]  2.289999e-01  2.021652e-01  6.086460e-01 -3.482049e-01  8.688996e-01
##  [891]  1.159607e+00  1.139108e+00  8.582507e-01  4.878922e-01  1.811369e+00
##  [896]  4.875576e-01  1.372216e+00  3.349516e-01  8.280674e-01  4.526575e-01
##  [901]  6.954205e-01  1.380078e+00  1.444320e+00  6.133215e-01  3.853755e-01
##  [906]  6.753975e-01  5.916860e-01  2.770411e-01  1.181686e+00  9.386963e-01
##  [911]  5.768844e-01  8.278391e-01  4.601030e-01  6.678635e-01  7.682568e-01
##  [916]  4.769896e-01  3.656810e-01  8.481616e-01 -2.665546e-01  4.623093e-01
##  [921]  1.267574e-02  2.344004e-01  7.683853e-01  2.390420e-01  2.154196e-01
##  [926] -2.788799e-01 -8.257558e-02  1.549374e+00  1.362439e+00 -2.204800e-02
##  [931]  1.168323e+00  9.985270e-01  6.727565e-01  8.951708e-01  6.349693e-01
##  [936]  1.928312e-01  8.305492e-01  1.058559e+00  1.234740e+00  1.365861e+00
##  [941]  9.936070e-01  6.001136e-01  2.383888e-01  1.058228e+00  9.844288e-01
##  [946]  6.948447e-01  1.434871e+00 -1.335447e-01  9.918569e-01  6.111783e-01
##  [951]  7.573554e-01  5.564200e-01  1.210209e+00  7.272664e-01  7.057156e-01
##  [956] -6.877341e-02  5.533912e-01  1.759287e+00  6.773524e-01  3.850037e-01
##  [961]  6.361800e-01  9.252708e-01  1.189280e+00  8.877908e-01  1.324969e+00
##  [966]  1.009075e+00  4.935633e-01  4.022411e-01  2.308686e+00  9.621204e-01
##  [971]  8.783990e-01  1.696462e+00  4.287349e-01  4.697207e-01  2.123158e-01
##  [976] -3.864153e-01  2.577847e-01  7.492179e-01  1.263827e+00  8.269815e-01
##  [981]  9.756019e-01 -4.127563e-01  6.128128e-01  3.438028e-01  6.185494e-01
##  [986]  9.346322e-02  2.776670e-01  6.214917e-01  1.290671e+00  8.016832e-01
##  [991]  8.976679e-01 -4.993407e-01  1.674590e+00 -3.362231e-02  4.591534e-01
##  [996]  1.216747e+00  1.052022e+00  1.111286e-01  8.292532e-01 -1.614211e-01
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
##   0.51080696   0.22910292 
##  (0.07244870) (0.05122502)
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
## [1]  0.404599821 -0.017175603  1.040550394 -0.398570174  1.557098007
## [6] -0.003155925
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
## [1] 0.0193
```

```r
sd(results2$thetastar)  #the standard deviation of the theta stars is the SE of the statistic (in this case, the mean)
```

```
## [1] 0.9078556
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
## t1*      4.5 -0.03413413   0.8572022
```

One of the main advantages to the 'boot' package over the 'bootstrap' package is the nicer formatting of the output.

Going back to our original code, lets see how we could reproduce all of these numbers:


```r
table(sample(x,size=length(x),replace=T))
```

```
## 
## 0 2 3 4 5 6 8 9 
## 2 1 1 1 1 1 1 2
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
## [1] 0.0614
```

```r
se.boot
```

```
## [1] 0.893985
```

Why do our numbers not agree exactly with those of the boot package? This is because our estimates of bias and standard error are just estimates, and they carry with them their own uncertainties. That is one of the reasons we might bother doing jackknife-after-bootstrap.

The 'boot' package has a LOT of functionality. If we have time, we will come back to some of these more complex functions later in the semester as we cover topics like regression and glm.

