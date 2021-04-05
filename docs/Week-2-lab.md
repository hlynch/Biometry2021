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
## 1 2 4 5 6 7 8 
## 2 2 1 1 1 2 1
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
## [1] -0.0245
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
## [1] 2.710498
```

```r
UL.boot
```

```
## [1] 6.240502
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
##    [1] 5.3 2.9 4.0 4.4 5.3 4.0 4.7 4.7 3.4 4.6 4.6 3.9 4.8 4.6 6.2 3.5 4.1 4.2
##   [19] 4.5 4.6 5.0 5.1 2.4 5.3 3.6 5.5 3.9 3.8 3.3 5.5 3.7 4.2 2.8 4.6 4.4 5.4
##   [37] 3.6 5.4 3.9 6.0 4.3 5.4 4.9 3.6 4.7 5.3 6.2 5.0 3.0 2.7 4.8 2.8 4.9 3.2
##   [55] 4.5 2.6 3.8 3.7 4.4 3.7 4.3 5.1 5.7 4.8 2.9 5.9 3.0 4.1 5.9 2.7 4.6 5.0
##   [73] 5.0 5.7 3.1 3.0 6.8 4.4 4.5 3.9 4.2 3.9 4.1 5.8 4.0 5.2 5.3 3.0 4.6 3.9
##   [91] 4.9 4.0 3.9 4.0 3.1 2.8 3.2 2.8 3.7 4.9 3.4 5.5 4.3 4.8 3.8 5.9 4.4 5.1
##  [109] 4.7 5.6 4.7 3.4 5.0 3.2 6.0 4.6 4.4 4.2 4.4 4.8 3.6 5.4 4.4 4.4 4.7 5.2
##  [127] 5.0 5.6 4.4 5.0 6.1 3.5 4.0 3.6 4.3 4.3 4.7 4.1 5.1 6.0 4.0 2.6 4.9 3.4
##  [145] 3.7 4.9 6.2 4.9 5.1 5.2 4.5 3.0 3.1 2.3 4.6 4.4 3.6 4.7 5.8 2.5 5.1 5.1
##  [163] 4.0 3.8 6.6 4.3 3.5 3.0 6.0 4.0 4.9 5.0 4.3 5.2 6.2 3.8 4.4 3.3 3.3 4.2
##  [181] 6.1 4.7 3.5 4.4 3.3 3.6 5.6 5.3 4.5 2.8 4.4 6.2 3.2 4.1 4.8 4.8 5.8 5.0
##  [199] 5.3 3.3 5.5 4.1 4.7 4.1 3.9 4.6 3.4 4.7 4.6 3.9 4.7 4.7 2.8 4.7 3.7 4.0
##  [217] 3.3 4.4 4.4 4.5 4.0 5.2 5.2 3.1 5.1 5.5 4.7 4.8 3.1 4.1 5.4 2.9 2.5 4.9
##  [235] 5.5 4.9 3.2 4.9 5.0 3.3 5.6 3.2 6.6 4.7 6.1 5.8 3.9 4.4 4.7 4.8 4.8 4.6
##  [253] 4.7 5.1 3.9 4.6 3.8 3.4 6.2 2.6 5.3 4.9 3.5 5.6 4.2 6.2 6.4 5.5 3.5 5.2
##  [271] 5.0 4.6 3.7 5.6 5.5 4.2 4.4 5.0 4.6 3.3 5.2 4.4 6.4 2.8 4.1 3.7 4.3 4.9
##  [289] 3.5 3.3 2.6 4.6 2.9 4.3 4.0 3.8 4.9 4.3 5.2 4.0 4.2 5.8 3.4 4.3 3.0 5.6
##  [307] 5.8 3.2 5.6 3.7 3.9 3.6 4.8 3.4 2.8 4.8 4.7 6.3 4.3 2.7 2.6 2.8 4.5 4.9
##  [325] 4.1 4.7 5.8 5.5 4.6 4.4 4.1 4.8 3.8 4.9 3.6 3.3 4.8 5.0 3.9 2.9 4.8 5.0
##  [343] 4.9 2.9 4.1 4.8 5.2 4.5 5.2 4.7 5.7 5.1 4.6 3.4 2.8 5.9 3.8 5.4 3.1 6.4
##  [361] 5.2 7.1 4.6 3.5 6.7 5.5 3.8 5.6 5.5 5.3 4.3 4.9 5.1 4.8 6.5 3.4 2.9 2.9
##  [379] 4.2 6.4 4.9 4.6 3.0 5.0 5.2 5.6 5.1 4.4 3.5 6.4 4.3 5.3 5.0 5.5 3.7 4.7
##  [397] 3.8 3.3 5.7 4.1 6.1 5.4 6.0 3.7 3.2 4.0 4.8 3.6 3.9 6.0 3.3 4.6 4.1 5.4
##  [415] 4.8 4.0 4.0 3.4 5.0 5.5 4.5 4.9 4.3 4.2 3.1 3.6 5.0 3.2 4.8 4.9 4.9 4.6
##  [433] 4.3 4.4 5.0 5.1 5.7 6.1 4.8 4.2 4.8 5.0 4.7 3.5 5.7 5.9 4.1 3.1 5.7 3.2
##  [451] 5.2 5.3 4.9 4.6 6.3 4.3 4.8 5.3 3.3 4.7 5.1 6.5 3.3 5.1 3.2 4.2 6.2 4.8
##  [469] 5.5 3.9 3.7 4.7 3.9 4.6 3.4 6.3 4.2 5.6 4.1 5.2 3.8 5.2 5.2 4.2 3.3 4.2
##  [487] 5.7 7.2 4.8 3.9 4.0 5.7 6.6 3.7 4.6 4.8 3.8 3.2 4.3 3.9 5.6 4.7 4.0 5.9
##  [505] 3.7 5.8 4.9 4.4 4.6 2.5 4.4 5.1 5.5 5.6 4.6 3.3 4.6 4.8 3.6 4.5 6.2 2.7
##  [523] 4.1 4.1 4.6 6.2 5.4 4.2 4.2 5.7 5.7 5.0 5.7 4.7 3.2 4.2 3.4 5.1 4.2 5.4
##  [541] 4.5 4.5 6.1 4.2 4.9 5.3 4.8 4.9 4.0 4.0 5.0 4.2 4.6 5.1 4.0 4.5 4.3 3.5
##  [559] 4.4 4.0 3.4 5.9 3.1 4.3 3.7 5.8 6.2 2.8 3.7 4.4 4.6 5.7 5.5 3.9 5.2 6.3
##  [577] 3.7 5.1 5.3 4.7 5.8 5.3 5.8 3.1 5.5 5.7 4.0 6.0 4.1 5.8 3.9 4.9 4.9 5.1
##  [595] 3.8 2.6 5.0 5.7 5.6 4.4 5.6 6.0 3.8 5.5 5.0 3.1 4.8 5.3 4.0 4.6 4.8 3.6
##  [613] 5.9 4.9 4.5 4.0 5.8 4.8 3.6 5.4 4.3 3.6 4.6 3.6 3.9 5.0 5.7 4.0 4.1 4.3
##  [631] 4.1 4.8 3.8 4.8 3.9 5.0 5.2 3.8 5.2 5.7 5.3 4.9 3.3 3.9 4.3 2.6 4.8 3.4
##  [649] 4.2 5.4 4.4 3.3 3.4 3.6 5.3 4.7 4.3 4.2 5.1 4.0 4.1 3.9 3.1 5.8 3.7 4.1
##  [667] 6.5 5.2 5.3 5.6 5.7 3.9 2.6 4.7 3.5 3.6 4.8 3.7 5.1 4.1 5.0 2.6 3.4 5.4
##  [685] 5.8 4.6 4.2 4.0 3.2 3.6 4.5 2.9 4.5 4.2 3.5 5.0 5.1 4.1 4.8 2.1 4.4 3.5
##  [703] 3.4 3.9 4.4 3.0 5.2 4.7 4.7 3.9 3.5 5.2 4.4 4.1 4.5 5.6 4.5 4.7 4.6 3.8
##  [721] 3.6 5.1 2.3 3.7 4.4 2.9 4.6 5.2 4.1 5.6 5.6 4.2 5.4 4.4 4.6 4.1 4.6 4.6
##  [739] 5.5 4.6 5.4 4.5 4.9 4.9 4.9 5.0 3.7 3.2 5.7 5.0 5.2 4.0 3.8 5.4 4.2 5.7
##  [757] 4.4 5.3 4.2 4.0 4.3 4.4 4.9 2.6 3.8 3.8 5.4 4.8 5.5 3.9 4.4 4.0 5.5 4.3
##  [775] 4.2 5.9 4.6 5.5 3.5 4.6 5.3 6.1 4.4 4.0 5.6 6.0 3.6 5.3 4.8 4.2 4.7 4.1
##  [793] 5.1 3.8 2.9 4.3 3.1 4.2 3.3 4.1 2.6 3.7 5.3 3.3 4.0 5.0 4.8 5.5 2.9 6.5
##  [811] 4.4 4.2 4.2 4.4 3.9 3.4 3.8 4.8 5.4 4.8 4.0 5.4 5.2 4.8 4.0 5.3 2.5 4.5
##  [829] 5.7 3.4 3.3 3.6 5.3 3.8 5.1 4.4 4.6 3.2 5.2 4.6 4.8 5.7 4.2 5.9 5.2 5.0
##  [847] 6.0 6.8 5.9 3.6 3.7 4.3 4.9 4.0 5.4 4.6 3.3 5.4 4.1 5.5 4.0 5.3 4.7 4.5
##  [865] 5.5 5.8 5.5 4.4 4.5 4.4 3.9 5.8 4.0 3.8 4.2 3.6 5.4 5.5 4.7 7.2 5.3 4.0
##  [883] 5.1 5.6 3.7 5.7 5.2 4.0 5.8 3.4 4.0 5.3 4.7 3.8 2.2 4.9 6.0 4.0 4.8 3.5
##  [901] 5.5 4.8 3.9 4.8 2.9 5.8 4.2 2.6 5.3 4.3 4.6 4.5 5.4 3.8 6.1 5.6 4.2 4.8
##  [919] 5.4 4.4 4.8 3.7 2.9 5.9 6.0 4.4 4.5 3.8 5.3 3.2 4.0 4.0 4.4 3.6 4.4 5.3
##  [937] 6.3 4.2 4.7 4.4 3.3 4.7 2.1 4.1 5.2 2.5 3.2 3.9 3.1 5.5 4.2 4.4 4.0 5.3
##  [955] 3.0 4.0 4.7 5.3 2.9 3.6 4.1 3.0 4.8 5.0 4.8 5.4 5.5 4.3 4.1 4.6 6.2 4.8
##  [973] 4.1 4.2 5.3 4.5 5.5 3.5 5.7 4.3 4.5 5.1 3.9 2.2 5.1 4.7 4.8 4.8 4.2 3.9
##  [991] 4.8 3.9 5.0 4.1 4.8 3.4 4.7 4.2 3.8 4.2
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
##    [1] 4.3 4.7 5.3 3.2 4.3 6.6 3.7 5.0 3.9 4.7 4.0 3.5 4.7 5.2 5.2 4.8 4.6 4.6
##   [19] 5.3 4.3 4.0 5.1 3.4 2.8 4.8 5.2 5.7 5.7 5.1 3.8 4.7 5.8 3.2 4.9 4.9 3.7
##   [37] 4.9 3.8 5.0 5.0 5.2 3.3 4.0 4.4 4.9 2.9 5.2 4.4 4.5 5.0 4.2 3.6 5.2 4.7
##   [55] 4.1 2.9 3.8 5.5 5.1 3.5 5.1 5.3 4.6 3.7 4.2 4.5 6.0 3.9 4.4 5.0 4.5 3.7
##   [73] 4.5 4.1 4.4 5.7 4.4 4.1 3.8 4.2 4.0 4.6 5.2 5.3 6.7 4.9 2.9 6.7 4.1 5.2
##   [91] 3.6 4.3 5.9 4.6 3.1 4.2 4.5 4.5 4.5 4.7 4.7 4.3 4.0 3.9 3.2 4.9 2.7 5.0
##  [109] 4.7 4.3 4.6 5.5 5.3 3.7 4.3 4.4 3.7 5.5 3.9 4.5 5.3 5.8 4.7 3.8 4.0 4.7
##  [127] 6.6 2.9 3.5 4.1 4.1 3.7 5.0 5.8 3.3 4.4 3.2 2.3 5.4 4.4 5.2 5.2 1.9 2.7
##  [145] 5.0 3.2 4.4 4.5 4.1 3.6 3.9 5.3 3.2 5.7 4.1 5.1 5.8 3.0 4.4 5.8 5.4 3.8
##  [163] 5.6 2.8 4.1 3.4 4.7 5.0 4.6 4.7 4.0 4.3 2.8 3.7 4.5 3.3 4.6 4.2 6.6 4.4
##  [181] 5.4 5.6 4.8 6.2 4.6 4.8 5.3 5.0 4.4 4.7 5.6 5.3 4.5 5.2 4.5 4.3 4.9 6.0
##  [199] 5.2 4.8 4.2 5.0 6.2 4.9 4.3 5.9 7.2 4.0 4.4 4.2 2.8 5.6 3.5 4.2 4.4 3.6
##  [217] 3.2 5.2 5.7 3.7 5.0 4.6 6.1 2.7 4.9 3.7 4.9 3.9 5.1 4.9 4.0 3.6 6.9 4.5
##  [235] 5.1 5.3 3.9 4.1 3.1 5.1 4.6 3.9 2.7 3.0 3.7 6.3 4.0 2.5 3.8 3.2 3.2 5.1
##  [253] 3.5 4.0 2.4 4.1 5.4 3.4 3.0 5.7 5.9 3.4 6.6 5.2 2.7 6.0 3.5 4.8 3.6 5.5
##  [271] 4.3 4.1 4.6 5.2 5.8 6.0 5.3 3.7 4.6 3.5 4.8 3.7 5.1 3.5 4.7 5.7 3.0 4.7
##  [289] 3.7 4.2 3.0 5.3 4.2 3.0 4.3 5.5 5.3 5.0 6.5 5.2 4.3 3.7 4.7 4.3 4.9 5.0
##  [307] 5.5 4.6 3.4 4.5 4.9 4.4 4.6 4.2 4.9 4.3 5.5 2.4 4.7 4.2 5.1 2.9 4.5 5.2
##  [325] 4.8 3.8 5.1 7.3 5.2 4.6 3.4 4.2 3.8 4.8 4.2 4.3 4.7 3.5 3.5 3.8 3.7 5.3
##  [343] 4.0 5.2 5.2 5.1 4.7 4.8 3.7 5.1 5.7 3.7 4.1 4.0 5.0 4.3 2.7 4.4 3.2 4.4
##  [361] 3.5 5.1 4.4 5.4 3.4 5.1 3.9 3.7 4.1 4.3 4.5 4.4 4.3 4.8 3.8 6.2 4.3 5.2
##  [379] 2.6 5.3 4.9 3.0 3.8 4.8 5.8 4.8 3.2 5.3 6.0 4.9 3.9 4.3 4.3 5.4 4.4 3.3
##  [397] 3.1 3.9 5.1 5.8 4.7 5.8 4.6 2.6 5.1 5.6 3.9 4.1 4.5 3.6 3.0 6.2 3.2 5.4
##  [415] 5.2 5.0 6.4 3.0 4.3 4.4 5.0 5.0 3.7 2.6 4.3 4.3 4.1 3.6 4.9 4.3 4.6 3.6
##  [433] 5.1 3.4 5.7 5.0 5.3 3.7 4.2 6.2 3.6 6.2 4.3 5.4 3.1 3.8 5.1 4.7 3.2 5.5
##  [451] 5.0 5.7 2.6 4.7 4.7 2.4 5.0 5.0 4.0 5.1 4.2 4.0 3.8 4.1 4.7 4.1 3.6 3.6
##  [469] 3.3 4.6 3.3 4.4 3.5 3.0 5.6 3.2 2.6 5.6 2.6 4.2 4.6 3.2 4.6 4.8 3.0 6.5
##  [487] 5.4 4.5 5.5 4.2 2.4 3.9 3.7 4.8 4.8 3.3 3.5 2.9 5.6 4.9 4.7 4.7 5.7 5.0
##  [505] 3.6 4.0 3.3 5.3 5.8 4.7 4.3 5.9 5.0 4.2 5.6 3.7 4.9 4.0 5.1 4.0 3.2 5.1
##  [523] 5.6 4.8 3.7 4.1 4.1 4.3 5.2 4.5 4.2 3.6 5.3 2.9 4.3 6.1 5.7 3.7 4.9 4.3
##  [541] 4.4 5.0 3.5 4.1 4.9 2.3 4.0 2.9 4.5 5.4 4.2 5.4 4.8 5.9 4.3 3.9 4.5 3.3
##  [559] 4.2 5.0 4.3 4.1 4.9 6.6 4.7 5.5 4.1 3.3 3.9 3.3 4.7 3.8 5.0 4.9 5.0 4.9
##  [577] 3.8 6.1 3.7 5.2 4.6 3.5 4.6 3.7 4.7 3.9 3.8 4.1 4.1 5.8 5.2 5.5 5.0 4.9
##  [595] 4.6 4.2 4.9 4.3 3.9 2.7 3.7 5.3 4.0 5.0 4.0 4.7 5.3 3.1 2.5 5.9 3.3 3.6
##  [613] 5.3 4.9 4.4 3.8 5.2 4.5 4.5 6.7 3.2 4.3 5.9 5.1 4.7 4.2 4.8 3.9 5.9 6.1
##  [631] 4.8 6.8 3.6 4.3 4.4 5.0 4.6 4.2 4.2 3.1 4.5 3.9 4.4 3.9 3.0 4.9 5.1 6.6
##  [649] 3.0 4.0 3.8 6.0 4.5 5.1 6.3 4.4 5.3 6.1 5.6 5.7 3.7 5.2 5.7 4.5 3.2 4.6
##  [667] 2.2 4.4 4.6 4.7 4.0 3.6 4.5 5.0 4.3 4.3 5.3 3.4 4.8 3.8 3.6 5.6 4.4 3.7
##  [685] 5.9 4.1 4.5 4.5 5.4 3.9 2.6 4.0 6.0 5.8 5.1 5.2 4.7 4.4 5.0 5.7 4.9 5.2
##  [703] 3.8 5.8 5.4 4.2 5.8 3.7 4.9 4.0 4.1 3.2 4.9 4.2 4.0 4.2 3.8 4.5 4.0 3.1
##  [721] 4.4 5.9 3.5 5.2 3.4 3.7 3.9 3.9 5.1 3.2 6.1 4.3 3.7 2.8 4.4 4.2 5.6 6.9
##  [739] 6.4 5.2 4.6 4.8 4.5 4.1 3.8 2.3 4.7 4.6 4.3 5.2 3.6 4.8 3.6 5.0 3.6 3.8
##  [757] 4.4 4.9 4.5 3.4 5.1 2.7 4.6 5.8 5.8 5.2 5.3 4.1 4.5 4.5 4.0 3.2 5.8 5.4
##  [775] 4.2 4.2 5.9 4.6 2.7 3.3 4.2 6.4 4.2 5.0 4.1 3.0 5.4 5.2 4.2 4.1 2.8 6.3
##  [793] 4.0 5.0 4.0 4.5 7.0 7.0 3.2 4.7 4.8 4.6 2.6 5.5 5.0 5.1 2.7 4.4 5.1 5.1
##  [811] 3.7 4.5 3.2 4.1 3.8 4.7 4.4 4.6 2.8 4.5 6.0 4.7 5.3 5.1 5.1 5.3 5.1 5.2
##  [829] 3.7 6.0 4.3 3.5 4.8 4.2 4.5 4.0 3.8 4.9 4.4 2.3 5.2 4.1 5.1 2.9 2.9 4.3
##  [847] 2.2 5.7 3.4 5.1 5.0 5.1 5.0 3.6 3.2 4.8 4.6 3.9 4.8 3.3 3.8 5.5 4.0 5.4
##  [865] 5.1 2.4 4.8 3.6 3.8 2.7 5.3 5.5 2.8 4.0 3.9 4.9 4.9 5.8 4.1 6.1 5.1 5.8
##  [883] 5.5 4.9 3.5 4.0 4.4 4.9 4.3 4.8 4.9 4.8 4.7 4.2 5.2 2.6 4.2 5.1 4.0 4.5
##  [901] 3.2 4.4 3.8 4.4 4.0 4.5 3.9 4.8 5.1 4.7 4.1 5.0 4.2 3.8 3.5 5.1 4.4 4.8
##  [919] 3.9 3.9 2.4 4.7 4.0 5.0 5.2 5.0 3.9 5.1 4.4 4.2 3.6 5.4 3.8 3.8 4.1 3.4
##  [937] 6.7 4.1 5.4 5.3 5.6 4.8 5.1 5.1 3.0 4.2 4.9 5.3 3.4 4.8 3.7 4.3 5.9 2.3
##  [955] 4.1 5.4 4.8 5.6 3.7 3.7 4.6 3.4 4.2 6.4 3.0 5.1 4.3 4.0 5.5 3.5 3.7 5.5
##  [973] 6.1 5.3 3.4 4.4 5.6 5.8 4.0 4.5 4.3 2.9 4.0 5.3 4.8 5.2 4.0 5.1 4.5 5.4
##  [991] 4.0 5.5 4.7 5.0 3.9 6.1 5.0 4.4 2.6 4.9
## 
## $func.thetastar
## [1] -0.0289
## 
## $jack.boot.val
##  [1]  0.4390769  0.3790560  0.3090652  0.1271709  0.0572973 -0.1156923
##  [7] -0.2229462 -0.3331565 -0.4011905 -0.5896359
## 
## $jack.boot.se
## [1] 1.004816
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
##    [1] 2.3 5.8 3.9 5.1 5.6 4.4 4.0 5.4 4.5 6.2 3.9 4.8 5.3 4.0 4.4 5.5 3.8 2.9
##   [19] 4.3 2.9 3.2 3.4 3.4 4.3 3.5 5.8 3.6 4.6 4.0 4.8 5.1 4.4 5.9 2.8 3.7 5.2
##   [37] 3.9 5.9 4.4 3.1 5.2 5.3 4.0 5.5 4.9 3.4 4.2 3.8 5.2 4.9 2.8 4.3 4.4 4.8
##   [55] 5.6 4.5 5.2 5.1 4.8 4.6 4.8 4.3 5.5 4.2 5.9 3.7 5.1 4.4 5.7 6.1 5.0 4.4
##   [73] 4.7 4.2 4.8 5.8 3.9 3.7 4.1 4.7 4.7 5.4 4.0 5.5 4.3 3.3 4.7 3.7 3.6 5.5
##   [91] 3.8 3.7 5.2 4.8 5.0 4.8 4.6 2.9 4.8 4.1 4.6 4.6 4.2 5.7 5.3 4.9 3.8 4.4
##  [109] 4.9 5.1 5.7 4.7 3.7 3.1 5.1 5.7 4.1 5.0 4.0 2.6 5.8 4.2 5.3 5.3 5.1 6.7
##  [127] 5.2 3.6 5.0 5.7 3.0 4.2 4.9 5.4 6.1 4.4 5.5 5.5 2.6 4.0 4.3 3.2 3.9 5.4
##  [145] 3.4 4.9 5.3 3.7 5.0 4.2 3.9 3.2 3.9 4.9 2.7 2.6 4.7 5.3 4.8 4.2 6.4 6.1
##  [163] 3.8 3.8 4.8 3.9 4.6 4.6 3.7 5.8 5.2 4.2 5.7 3.3 6.7 3.7 3.7 4.9 5.7 6.1
##  [181] 5.3 5.4 3.1 4.7 4.6 5.3 3.7 2.9 4.8 3.8 5.4 4.9 5.6 4.8 5.3 4.4 4.6 4.3
##  [199] 4.3 4.5 5.4 5.0 2.6 2.2 5.6 4.3 3.9 4.6 4.3 2.4 5.0 3.9 5.3 4.5 4.9 4.9
##  [217] 5.3 2.6 4.3 4.4 5.1 4.6 5.2 4.4 5.4 4.8 4.6 5.2 5.5 3.4 4.0 5.5 5.5 3.8
##  [235] 4.9 5.4 4.8 3.5 4.7 4.8 3.4 4.0 2.2 4.1 2.9 4.2 6.3 4.9 4.8 5.2 4.3 4.1
##  [253] 5.5 2.9 4.9 5.1 4.1 5.1 4.0 4.9 4.3 5.4 5.6 4.7 6.4 4.6 3.8 3.8 3.3 5.9
##  [271] 4.5 4.9 6.9 4.3 4.0 5.7 4.9 4.1 4.6 4.6 4.7 5.3 4.7 5.3 4.7 6.0 2.2 3.0
##  [289] 3.3 5.2 4.7 4.5 3.7 3.4 3.8 5.1 4.5 4.1 4.7 3.2 4.1 5.0 5.0 4.8 3.4 4.5
##  [307] 6.2 3.7 6.0 4.1 4.2 2.4 4.3 5.4 4.5 5.0 4.9 4.2 7.1 5.0 4.7 3.8 3.7 3.9
##  [325] 3.4 4.2 3.5 4.7 5.0 5.7 4.0 5.5 4.7 5.0 5.8 5.8 3.3 4.0 5.0 3.5 5.0 5.2
##  [343] 4.5 2.7 4.0 4.1 3.0 4.7 3.2 4.7 4.2 5.7 3.5 5.6 3.5 3.6 5.5 4.0 5.3 3.6
##  [361] 4.6 4.0 3.9 3.9 4.6 4.1 5.9 4.5 5.2 4.1 5.4 5.5 4.3 5.5 4.5 3.9 4.4 4.9
##  [379] 5.3 2.6 2.1 4.1 4.3 3.3 4.4 6.2 4.0 3.2 5.3 4.1 4.2 3.6 3.8 2.3 3.4 3.7
##  [397] 3.5 4.4 4.6 3.8 4.2 4.6 5.3 4.7 4.2 5.7 4.7 4.1 4.7 4.6 4.8 4.8 5.4 6.2
##  [415] 5.5 4.6 4.3 5.5 5.0 5.9 5.4 5.1 4.4 3.2 4.4 3.3 4.0 5.9 5.1 3.7 4.1 5.6
##  [433] 2.9 3.7 3.7 3.5 3.5 4.9 3.7 5.5 5.5 6.3 6.2 4.2 5.3 3.6 6.1 4.4 5.4 2.3
##  [451] 3.8 4.9 6.0 3.4 4.8 5.4 4.4 4.1 4.2 3.0 4.6 4.1 2.0 5.3 5.0 3.0 4.8 5.0
##  [469] 5.2 4.6 2.8 4.9 2.7 3.5 3.7 4.2 5.2 5.2 3.9 3.9 4.7 3.4 3.7 5.7 5.0 5.9
##  [487] 3.7 3.2 5.2 5.2 4.7 4.9 4.9 4.8 4.6 5.3 4.1 5.5 3.5 1.8 4.4 4.7 6.1 4.3
##  [505] 4.3 4.5 5.2 4.9 5.2 3.1 5.2 4.9 6.7 5.4 3.7 4.6 3.7 5.1 5.0 5.0 3.9 4.4
##  [523] 5.6 4.5 5.6 4.1 5.2 5.5 6.0 3.6 5.8 2.7 3.1 6.0 3.4 4.3 4.5 4.7 3.8 5.5
##  [541] 3.6 3.8 2.3 3.6 3.4 3.7 3.9 3.6 3.8 4.2 4.4 3.8 4.0 3.6 5.0 4.3 4.5 4.9
##  [559] 4.8 4.6 3.3 3.4 5.7 4.4 4.2 4.1 5.6 5.8 4.4 4.5 3.1 6.5 5.8 3.2 5.2 6.3
##  [577] 3.9 6.1 4.5 4.8 3.7 5.5 4.4 4.1 5.2 3.6 5.8 4.4 3.2 4.7 4.9 4.2 4.1 3.5
##  [595] 4.3 4.2 3.1 4.1 3.3 4.7 3.5 5.1 4.5 4.1 3.9 4.7 4.0 4.3 5.0 4.5 5.0 3.0
##  [613] 4.0 2.9 5.6 5.9 4.7 5.1 4.0 4.7 3.7 5.8 5.6 5.4 4.0 3.3 5.0 4.1 4.3 3.8
##  [631] 5.0 4.1 3.2 4.5 3.9 6.6 2.8 4.5 3.5 3.6 3.5 2.3 4.0 5.2 5.9 4.0 4.8 4.3
##  [649] 4.2 3.0 4.0 3.3 3.4 5.3 4.3 5.1 2.6 3.1 3.0 4.3 3.5 2.9 4.5 4.5 3.0 4.8
##  [667] 4.3 4.2 2.9 4.6 4.5 5.3 4.5 3.7 4.4 6.2 5.7 3.8 5.8 3.9 5.6 5.0 3.3 5.1
##  [685] 4.0 5.0 4.6 5.0 4.1 5.4 4.8 4.8 4.6 6.1 4.4 3.8 3.3 4.0 3.3 5.2 4.2 4.5
##  [703] 2.4 4.0 4.2 4.5 4.9 4.3 3.6 7.2 5.9 3.7 5.5 3.8 4.5 5.1 2.9 5.5 3.8 4.6
##  [721] 4.7 6.0 5.4 5.3 4.3 4.9 5.3 5.6 3.5 5.4 4.0 4.4 5.2 4.4 1.5 2.6 5.1 5.5
##  [739] 2.9 5.4 5.1 4.6 2.9 4.1 3.9 4.3 3.9 4.0 2.4 5.2 5.2 6.9 2.5 5.2 3.9 5.8
##  [757] 3.9 4.8 5.1 5.6 4.7 3.3 3.2 4.0 4.2 4.9 6.9 4.3 4.1 4.4 5.8 5.7 2.6 3.7
##  [775] 4.5 3.7 3.6 4.7 3.9 3.7 4.5 4.9 4.6 6.7 5.7 3.5 3.8 4.4 3.9 5.4 5.0 5.6
##  [793] 4.9 4.2 4.6 4.4 5.6 4.6 3.5 3.6 5.8 3.2 4.5 3.8 4.1 4.3 4.6 3.3 4.9 5.6
##  [811] 5.1 3.8 2.3 3.0 4.8 5.0 2.4 3.9 3.7 4.4 4.6 5.4 5.7 5.0 2.7 5.3 3.1 3.8
##  [829] 3.8 2.6 5.9 4.6 5.4 4.2 4.8 4.6 4.1 4.7 4.7 4.8 4.5 3.6 4.2 3.2 4.3 5.3
##  [847] 4.1 2.9 3.5 3.4 5.2 4.4 4.5 5.5 4.3 4.0 5.2 2.9 2.9 4.1 3.0 6.0 5.8 4.3
##  [865] 5.4 5.1 3.0 3.4 5.0 3.0 5.5 6.0 6.2 5.9 3.5 4.1 4.1 5.3 3.4 4.1 5.0 4.0
##  [883] 5.7 4.1 4.2 4.5 3.2 5.1 3.9 4.3 5.4 3.9 3.5 4.0 6.0 4.5 3.9 4.6 4.6 3.9
##  [901] 4.3 2.9 4.0 4.2 3.7 5.0 5.8 3.8 3.7 3.8 4.2 4.6 5.5 6.6 4.4 4.7 3.2 3.4
##  [919] 3.9 4.8 4.9 3.1 4.2 3.2 3.6 4.2 5.9 4.3 4.3 5.2 4.9 6.3 3.6 3.2 4.4 5.0
##  [937] 5.1 6.2 4.0 3.3 5.5 3.1 3.4 2.9 3.9 6.0 3.1 5.0 3.6 6.7 4.6 5.4 5.4 5.3
##  [955] 5.7 5.0 5.0 4.8 4.3 3.8 5.7 4.7 6.3 4.1 5.0 3.9 4.8 4.5 3.6 5.1 5.5 4.9
##  [973] 3.3 3.1 4.8 4.5 6.0 6.4 4.6 3.1 4.7 6.6 5.1 3.8 4.2 5.3 5.7 5.9 4.3 5.6
##  [991] 4.5 4.5 4.6 2.1 3.8 3.5 5.3 6.0 5.5 3.6
## 
## $func.thetastar
## 72% 
##   5 
## 
## $jack.boot.val
##  [1] 5.500 5.400 5.324 5.200 5.100 5.000 4.900 4.700 4.500 4.600
## 
## $jack.boot.se
## [1] 0.9813066
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
## [1] 0.7015565
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
##   1.7800708   3.4479303 
##  (0.7336604) (1.6392979)
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
## [1]  1.27407396  0.47791336 -0.07083241  0.34290988  0.07048136  0.75889843
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
##    [1]  0.079550346  0.255372324  0.142246649  0.848919601  0.997702704
##    [6]  0.168900159  0.481351664  0.638393857  0.767616360  0.383316549
##   [11] -0.007068761  0.467784087  0.198745506  0.545606088  0.344775792
##   [16]  0.758462897  0.326949512  0.377378350  1.391971394  0.338367829
##   [21]  0.293649941  0.300299795 -0.360029280  0.176021312 -0.047080985
##   [26]  0.642975991  0.633157505  0.397481502  0.502691927  0.348903877
##   [31]  0.047073225  0.288898069  0.264450892  1.129997048  0.764100968
##   [36]  0.328482842  0.752312283  1.198471475  0.187217555  0.539314594
##   [41]  0.589813878  0.378615047  0.488178157  1.474467727  0.548265329
##   [46] -0.345590334  0.624088668  1.004279849  0.220245050  0.189116432
##   [51]  0.453136118  0.937343826  0.120385721  0.506647752  0.551219053
##   [56] -0.089660600  0.449768656  1.043907890  1.153516751  0.605798617
##   [61]  1.022352962 -0.037923235  0.723189639  0.376508936  1.294784741
##   [66]  1.217963384  0.934968331  1.233383479  0.284467844  0.966057379
##   [71]  0.539162895  0.299105195  0.666084095  0.005545988  0.803250588
##   [76]  0.738619192  0.223434515  0.408501016 -0.410768491  1.377892303
##   [81]  0.115081362  0.567120289  0.254800331  0.656005102 -0.265182464
##   [86] -0.318876446  0.480381624  0.774482527  1.015170256  0.648849040
##   [91] -0.065498993  0.829875582  0.404830301  1.028635898  0.176468333
##   [96]  0.087886326  0.428380712  0.774177739  0.362834267  0.659251296
##  [101]  0.259122833  1.636009199  0.385864883 -0.505615776  0.605085826
##  [106]  0.468067170  1.064017550  0.153171040  0.750244136  0.136299062
##  [111]  0.095082258  0.080740320  0.486222928  1.236010197  0.126480503
##  [116]  0.456202696  0.639675811  0.266058958 -0.003873680  1.475240292
##  [121]  0.719368586  0.113482474  0.320527741  0.342028354  0.924245579
##  [126]  1.140954711  0.657043805  0.437912566 -0.413817610 -0.335238079
##  [131]  0.772953770  0.524330115  0.237968337  0.757147389  1.030762315
##  [136] -0.262528130  0.790565823  1.431952537 -0.105775354  0.397869843
##  [141]  0.848358917  0.313035794  0.464416118  1.291447107  0.853310665
##  [146]  0.284638459  1.081409241  0.595261829  0.192662825 -0.083833194
##  [151] -0.193878152  0.446685956 -0.630634474  0.296895767  0.245287562
##  [156]  0.422420807  0.402155099  0.310637853  0.582881891  1.050228384
##  [161]  1.139261094  0.189884284  0.348798219 -0.019856653  1.211442003
##  [166] -0.297396608  0.558012627  1.295595589  0.481926748  0.170131348
##  [171]  0.394294900  0.324632692  1.051998580  0.214820563 -0.237329306
##  [176] -0.409849341  0.383980730  1.512057841  0.518799198  0.456120947
##  [181]  0.896733848  0.423853158  0.707392619 -0.156768095  0.175274488
##  [186]  0.870590408  0.706362858  0.567058580 -0.304538482  0.507151202
##  [191]  0.680710872  0.452993091  0.143923872 -0.384616306  0.262410784
##  [196]  1.004927740  0.305606197  0.120322730  0.870429796  1.063178628
##  [201]  0.932665292  0.705842344  0.652549063  0.612586855  0.846756451
##  [206]  0.730347145  0.463097021  0.359761792  0.103693380 -0.115760391
##  [211]  1.114160204  0.542606134  1.070426126  1.231187069  1.025402427
##  [216] -0.119145812  0.848541648  1.265962872  0.864260816  0.786656501
##  [221]  0.716576791  0.845202850  0.117258266  0.796134876 -0.039445965
##  [226]  1.504946773  0.997936883  0.224049431  0.845945941  0.735281972
##  [231]  0.375456791  0.317027039  0.409426710  0.697341250  0.383452601
##  [236]  0.611590391 -0.361329467  0.819052961  0.280612847  0.795128331
##  [241]  1.451451073  0.132613304  0.506491571  0.087873643 -0.030369134
##  [246]  0.555123549  0.884167877  0.823936856  0.299480376  1.589668024
##  [251]  1.071348098  0.683870160 -0.008601139  0.350222963  1.539535978
##  [256]  0.605085826  0.164449233  0.984076525  0.798709733  0.863297705
##  [261]  0.463955296  0.068973577  0.402929609  1.746584405  0.772388023
##  [266]  0.448379860  0.068270728  0.764578085 -1.239363816  1.390946548
##  [271]  0.279533401  0.716813202  0.524404689 -0.330354633  0.329010957
##  [276]  0.338807797 -0.129088519  0.678052605  0.556068725  0.036476830
##  [281] -0.056726927  1.663318083  0.149638869  1.628751477  1.029729174
##  [286]  0.256072095  0.072693236 -0.412389708  0.629710631  0.749303678
##  [291]  0.483417323  0.737502387 -0.026030808 -0.128721737  0.473815513
##  [296]  0.573457329  0.198264797 -0.146383997  0.545138556  0.515781809
##  [301]  0.661449209  0.552939907  0.068143451  1.579003791  0.896370234
##  [306]  0.699396648  0.698484355  0.571351805  0.497555532  0.571429585
##  [311] -0.585137319  0.329992300  1.212097052  0.719647828  0.935161780
##  [316] -0.589430783  0.599098884 -0.245328887  0.379174192  0.787724317
##  [321]  0.431191882  0.727636957  0.816609841  0.035385681  0.648808221
##  [326] -0.165196572  0.270110228  1.098362027  0.899181090  0.658224007
##  [331]  0.396297735 -0.940749433  0.404326800 -0.095851340  0.400423006
##  [336]  0.386671578  0.367747301  0.218471986  0.494002116  0.666561286
##  [341]  0.744448597  0.793766784  1.085933714  0.746294072  0.812427621
##  [346] -0.343787675 -0.251804201  0.158023892  1.010883045 -0.026502780
##  [351]  0.956091711  0.482175561  0.658624624  1.080519346  1.160030676
##  [356]  0.293808364 -0.246244469  0.939829573  0.684050417  1.227504285
##  [361]  0.205798382  0.516515148  0.408281953  0.409218234  1.161021241
##  [366]  0.355220206  0.522280620  0.609558552  0.844512686  0.262410784
##  [371]  0.180104565  0.978959848  0.474556096  1.186418351  1.024934751
##  [376]  0.925464585  0.314828282  0.109375625  0.446287381  0.057730320
##  [381]  0.242586761  0.463873565  0.688167975  0.574776487  1.170524200
##  [386] -0.049554107  0.983633974  0.776251561  0.686707392  0.512596047
##  [391]  0.900776129  0.669607056  0.503848041  0.699810438  0.520490495
##  [396] -0.428251945  0.095428898  0.554695967  0.684052136  0.945525685
##  [401]  0.439818621  0.531902920  0.096426445  0.524756033  0.843564965
##  [406] -0.365208181  0.085684967  1.024306043  0.673705515  0.943105533
##  [411]  0.735654439  0.445508847  0.330782170  0.624119490  0.140581814
##  [416] -0.192930933  0.166872981  0.706848527  0.493651370  0.879388514
##  [421]  0.558450353  0.137001804  0.207860286  0.482784744  1.259443663
##  [426]  0.867643885  0.785917566  0.453874971  0.220763505  0.829607697
##  [431] -0.009922793  0.921723859  1.083738960  0.479804120  0.704143642
##  [436]  0.400487517  0.626124285  0.330725745 -0.127441803  0.945209268
##  [441]  0.630321691  0.134222031  0.831817842  0.438241177  0.794126585
##  [446] -0.030473468  0.190084975  1.300132602  1.566244890  0.101468544
##  [451]  0.364929383  0.769295966  0.947306405  1.113390483  1.001604080
##  [456]  0.814416516  0.442016954 -0.081558627  0.521007504 -0.049520393
##  [461]  1.938602679  1.049720761  0.543880457  0.305946274  1.113819806
##  [466]  0.861578309  0.751631709  0.646041352  0.108411833  1.097145160
##  [471]  1.021891213  0.281885055  0.522278195  1.926600619  0.757894056
##  [476]  0.893595052  0.518691309  0.671168466  0.801291875  0.071878699
##  [481] -0.068882995  0.305189990 -0.110482760  0.637365467  0.747753727
##  [486]  0.549103540  0.668497173  1.289535604  0.304881421  0.684059043
##  [491]  1.095468653  1.414718151  0.259529029  0.585555435  0.087585761
##  [496]  0.293722999  0.217641007  0.699663504  0.354499165  0.653919248
##  [501]  0.490844438  0.642092104  0.413537069  0.810582603  0.730347145
##  [506]  0.740790477  2.005243621  0.199501491  0.748596570  0.951324295
##  [511] -0.072105933  1.331304667 -0.341953572 -0.782008880  0.735322242
##  [516]  0.401641733  0.587423897  0.742234021 -0.340056153 -0.081226654
##  [521]  0.456151848  0.145297123  0.471720182  1.154358134  2.080927246
##  [526]  1.299875777  0.388354059  0.495738458  0.181587437  1.227871205
##  [531]  0.529250522  0.221799403 -0.026054722  0.016872451  0.930432598
##  [536]  0.261358658  0.618525153  0.986256900 -0.097496814  0.979799816
##  [541]  0.352161310  0.019828805  1.158432528  0.448225469  0.737402307
##  [546]  1.220558597  0.153450075  0.373858801  0.748647994  0.460811156
##  [551]  0.988085409  1.114583934  0.554635630  0.338740850 -0.155528861
##  [556]  1.149435888  0.430176781  0.488123810  1.660054065  0.261602527
##  [561] -0.270634483  0.698430824 -0.035493010  0.633700774  0.379190307
##  [566]  0.774192156  0.442284547  0.607428474  1.117686087  0.910854977
##  [571]  0.682383019  0.027366811  0.606269915 -0.588143764  0.763629149
##  [576]  0.179188397  0.201361230 -0.385207181  0.902148854  0.295356820
##  [581]  0.445951654  0.609873786  0.774036363  0.899640323  0.522217559
##  [586] -0.194060291  0.471540464  0.842698593  0.620943699  0.418427351
##  [591]  0.852122864  0.179038736  0.587577536  0.913697606  0.932684742
##  [596] -0.063553360  0.827992756  0.232110894 -0.328037747  1.541849655
##  [601] -0.185641476  1.907218354  0.092587400  0.396760559 -0.182489802
##  [606]  0.996552131  0.137766803 -0.087873781  0.812181882 -0.265295015
##  [611]  0.491864877  0.756430121  0.351238477  1.315864726  0.968993687
##  [616] -0.791323702  0.803538870  0.679730071  0.428441569  0.625721465
##  [621]  0.172309084  0.544242965  0.048028126  0.002975588 -0.122822360
##  [626]  0.718404266  0.386940469  0.474535211  0.238751599  1.284819201
##  [631]  0.808590080  0.848243869  0.328387972  0.151764115  0.196594261
##  [636]  1.022322192  0.871622347  0.811009008 -0.348384668 -0.070877574
##  [641]  1.048250976  0.983047141  1.347653545  0.279335839  0.488583728
##  [646]  0.308884343  0.589168160  1.036622655  0.335221395  1.229991305
##  [651]  0.834615464  1.084474999  0.481856704  0.418772152  0.639373465
##  [656]  0.483262927 -0.108353306  0.129850026  0.220158827  0.444022239
##  [661]  0.524741213  0.702877177  0.593184075  0.558155936  0.391194746
##  [666]  0.760787734  0.115679842  0.659968229  0.502735553 -0.232016660
##  [671] -0.011415302  0.622103485  0.447949236  0.714997513  0.719107501
##  [676] -0.087662746  1.291560294  0.311594667  0.170428680  0.531555106
##  [681]  0.667194099  0.159543417  1.299213589  0.300997106 -0.246116470
##  [686]  0.540006167  0.576158140  1.057846553  0.670536776  1.270280043
##  [691]  0.406330599  0.386806471  0.467115353  0.116130809  0.195562495
##  [696]  1.034801475 -0.667940156  0.524444086  0.644494295  0.503166427
##  [701]  0.476127416  0.040102994  0.522839322 -0.129091292  0.563930928
##  [706]  0.198301143  0.213979744 -0.408610388  0.181389818  1.042338425
##  [711]  0.819640548  0.220793670  0.411102864  1.027401486 -0.253674407
##  [716]  0.432382108  0.499206194  0.477366305 -0.104730330 -0.086518944
##  [721]  0.558150533  0.342026143  0.186244746  0.213979744  0.796225046
##  [726] -0.163906887  0.947055272 -1.034429522 -0.153083098  0.425172334
##  [731]  0.414381652  0.314126109  0.631887257  0.561813646  0.579813973
##  [736]  0.135603868  0.107778377  0.195562495  0.161875655  0.439713875
##  [741]  0.301513259  0.353998303  1.177798944  0.945789308  0.460349437
##  [746]  0.549611394  0.269882352  0.435885873  0.584358380  1.622141258
##  [751] -0.522091169  1.423308849  0.608063570  1.262455344  1.868921281
##  [756]  0.602932209  0.642446813  0.833745040  1.104140543 -0.022049063
##  [761]  0.418993218  0.115839853  0.014816975  1.234435470  0.217062292
##  [766]  0.086623728  0.590805591 -0.445165930  0.889725972  1.154652535
##  [771]  0.863240786  0.163337594  0.537783877  0.824002037  0.593293314
##  [776]  0.938431963  1.752033023  0.272091072  0.145356770  0.215539794
##  [781]  0.276652940  0.576408735  0.769314147  0.778158279  0.233937098
##  [786] -0.276284337  0.461048888  0.193953091  0.806249765  0.608296420
##  [791]  1.272125614 -0.018777084 -0.514090229  0.916180613  0.973339138
##  [796] -1.159150958  0.664019719 -0.325119027  0.577367064  0.793381171
##  [801]  0.794196236 -0.450203189  1.354475077 -0.041262175  0.945525685
##  [806]  0.354258859  0.994900792  0.573982074  0.505499369  1.611218534
##  [811]  0.883567121 -0.670554727  0.168226756  0.842839782 -0.348963503
##  [816]  1.522867771  0.718215922  0.457350082  0.304014069  0.458280067
##  [821]  1.293736077 -0.523073562  0.260782804  0.553213373  0.692249014
##  [826]  0.742780736  1.170191760  1.104778100  0.047164712 -0.010705648
##  [831]  0.228107213  0.477802866 -0.011153438  0.801684152  0.055331155
##  [836]  1.058318251  0.925101846  0.488178157  0.247580767  0.186610303
##  [841] -0.016495103  1.016323550  0.034890391 -0.057294361  1.196180649
##  [846] -0.188458916  0.751509582  0.130912838  0.999578654 -0.145392989
##  [851]  0.080939985 -0.008825742  0.016473370  0.703064757  0.419194731
##  [856]  0.991514736 -0.445917599  0.410306396 -0.355613297  0.349045092
##  [861]  0.216382124  0.260840400  0.537842498  1.132463562  0.518238302
##  [866]  0.378939974  0.162267708  0.646925679  0.161799133  0.041381304
##  [871]  0.960639997  0.920855685  0.449754308  0.618316463  0.131742813
##  [876]  0.773310923  0.321540127 -0.050075573  1.205275636  0.721095919
##  [881] -0.161322953  0.698766986  0.067978837  0.296024333  0.189514434
##  [886] -0.162374179  0.849896051  0.337478066  0.917520567  0.426763836
##  [891]  0.189503478  0.588335395  0.651023577  1.034782230  0.948242592
##  [896]  0.948965608  0.413245218  0.119933078  0.343188369  0.360321023
##  [901]  0.790880398  0.100767530  0.969310562  0.766934910  0.247280751
##  [906] -0.113195840  0.032347185  1.485155131  0.488220095  1.111470758
##  [911]  0.777891478 -0.108679615  1.510962424  0.159615740  1.051742832
##  [916]  0.560828721  0.706299315  1.383424869  0.955226832  0.497724448
##  [921]  0.731603403  0.004854182 -0.119762665  0.090027151  0.564692391
##  [926]  0.631840556  0.700648538  1.056199201  0.133556990 -0.106611004
##  [931] -0.515280374 -0.067911327  0.079112449  0.898175326  0.474535211
##  [936]  0.467298677  0.948330637  0.634262142  0.167062244  0.369137181
##  [941]  0.497500142  0.865492795 -0.736426798  0.765494861  0.865025559
##  [946]  0.742029150  0.760556273  0.867923305  0.620855362  1.110886556
##  [951]  0.595075768  0.280772361  2.003314890  0.720257230  0.146637138
##  [956]  0.219718472  0.938714320  0.373386310 -0.040368882  0.403138333
##  [961]  0.456029829  0.883904009  0.143044392  1.620524789  0.422420807
##  [966]  1.322028629  0.307127778  1.020408559  0.819474637  0.490647953
##  [971]  0.515819888  0.449545905  0.163443995  0.875393855  0.827992756
##  [976] -0.072232179  0.423758420  0.295402798  0.420262028  0.188067311
##  [981]  0.095355662  0.161799133  0.052568539  1.214127313  0.242429339
##  [986]  0.245249269 -0.100332748 -0.097191920  0.234283920  0.869382429
##  [991] -0.057755687  0.198951901  0.767978477 -0.308573699  0.867340497
##  [996]  0.253600683  2.037125104 -0.163316368  0.423895693 -0.159412178
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
```

```r
fit2
```

```
##       mean          sd    
##   0.51627157   0.36663973 
##  (0.11594166) (0.08197964)
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
## [1]  0.83024679 -0.25346990  0.92453813  0.03904544 -0.14455324  0.01316868
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
## [1] -0.0074
```

```r
sd(results2$thetastar)  #the standard deviation of the theta stars is the SE of the statistic (in this case, the mean)
```

```
## [1] 0.9031727
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
## t1*      4.5 0.03423423   0.9006729
```

One of the main advantages to the 'boot' package over the 'bootstrap' package is the nicer formatting of the output.

Going back to our original code, lets see how we could reproduce all of these numbers:


```r
table(sample(x,size=length(x),replace=T))
```

```
## 
## 2 3 6 8 9 
## 1 2 2 2 3
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
## [1] 0.0214
```

```r
se.boot
```

```
## [1] 0.9022949
```

Why do our numbers not agree exactly with those of the boot package? This is because our estimates of bias and standard error are just estimates, and they carry with them their own uncertainties. That is one of the reasons we might bother doing jackknife-after-bootstrap.

The 'boot' package has a LOT of functionality. If we have time, we will come back to some of these more complex functions later in the semester as we cover topics like regression and glm.

