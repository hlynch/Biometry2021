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
## 0 1 4 5 6 7 9 
## 2 1 1 1 2 1 2
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
## [1] -0.0072
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
## [1] 2.711728
```

```r
UL.boot
```

```
## [1] 6.273872
```

Method #2: Simply take the quantiles of the bootstrap statistics


```r
quantile(xmeans,c(0.025,0.975))
```

```
##   2.5%  97.5% 
## 2.7975 6.2000
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
##    [1] 4.2 5.2 4.3 4.7 3.5 4.8 5.2 3.4 4.7 4.9 5.1 6.4 5.2 6.0 4.0 4.4 5.1 5.9
##   [19] 5.6 6.8 5.1 5.1 4.3 4.4 4.7 5.2 4.4 4.4 4.2 5.5 3.3 3.3 3.6 4.2 4.9 5.4
##   [37] 4.6 3.6 5.0 4.8 5.8 4.3 3.1 4.5 4.5 3.8 3.1 5.3 5.6 6.0 5.8 2.8 4.1 5.7
##   [55] 5.7 2.9 4.5 4.2 4.5 4.2 4.7 4.9 5.3 4.0 3.0 3.5 4.7 4.8 4.8 5.8 2.9 5.0
##   [73] 5.9 5.3 5.1 4.8 6.5 4.4 5.7 5.0 4.6 3.5 3.6 4.3 2.3 4.4 4.8 5.0 4.6 5.4
##   [91] 3.7 4.0 4.5 4.9 4.6 3.3 4.3 5.5 4.8 5.1 4.7 4.6 4.0 5.1 5.2 5.7 4.9 4.1
##  [109] 4.4 3.8 2.9 5.0 4.3 3.9 4.7 4.5 4.8 4.7 4.1 4.4 5.1 4.1 4.2 4.2 3.3 6.3
##  [127] 4.6 4.1 4.3 3.6 4.4 4.2 5.4 4.1 4.4 5.3 4.4 6.3 5.7 3.3 4.6 5.4 4.6 4.5
##  [145] 6.1 4.8 4.8 5.7 5.6 2.6 3.6 4.9 3.2 4.1 5.3 4.4 4.3 3.7 3.8 5.1 6.5 5.2
##  [163] 3.5 6.1 6.2 5.4 4.7 6.7 5.0 4.2 6.7 3.8 3.0 4.2 4.9 4.7 4.8 3.6 5.6 6.2
##  [181] 4.9 3.1 5.1 5.2 4.2 4.9 4.2 4.8 4.7 4.7 4.6 3.4 4.9 4.5 5.1 4.7 4.0 3.5
##  [199] 4.8 3.7 4.4 3.5 4.1 4.2 4.6 3.2 5.5 3.4 4.2 3.2 4.4 4.9 4.2 3.9 5.1 3.2
##  [217] 4.5 3.4 3.3 4.6 6.0 4.3 4.0 4.1 4.5 3.5 3.2 5.3 4.8 4.3 4.8 6.3 2.7 5.4
##  [235] 5.2 5.4 5.3 3.2 3.5 3.5 4.6 5.0 4.7 5.2 5.6 5.8 4.8 5.3 4.6 3.9 4.7 5.4
##  [253] 5.2 4.3 4.0 3.9 5.2 4.1 4.8 4.4 4.6 4.8 3.9 2.5 4.8 5.4 5.1 3.7 5.3 5.9
##  [271] 4.6 5.9 3.5 4.7 4.2 3.5 5.7 4.1 4.7 2.9 3.8 5.5 4.7 3.7 4.1 4.7 6.0 3.6
##  [289] 5.7 3.4 4.0 4.7 2.6 6.0 4.3 4.3 4.5 3.3 6.1 7.2 4.3 3.7 3.5 5.5 3.1 6.2
##  [307] 5.4 4.8 5.9 4.7 4.1 5.7 5.8 6.6 5.7 5.3 2.9 5.9 2.7 6.5 3.9 5.2 3.6 6.0
##  [325] 5.8 4.3 5.1 5.4 4.6 5.2 4.9 4.5 3.2 5.3 3.8 4.2 3.7 5.3 5.9 6.0 4.6 4.5
##  [343] 4.4 5.2 5.4 4.7 5.3 5.6 5.6 4.9 4.7 4.9 4.6 4.3 5.6 5.2 5.2 3.6 6.3 5.8
##  [361] 3.2 4.2 4.6 4.4 5.1 3.6 4.9 4.1 4.1 4.2 1.9 5.6 3.9 3.2 3.8 3.2 5.4 5.3
##  [379] 4.6 6.0 3.1 5.0 4.3 4.2 6.1 5.9 4.1 5.0 5.4 6.0 3.8 5.2 2.9 3.7 6.2 4.1
##  [397] 3.5 5.4 3.3 4.4 5.3 5.9 4.3 4.5 4.9 4.3 5.0 4.1 5.4 4.8 3.8 3.8 4.6 5.3
##  [415] 4.5 4.4 4.6 6.0 4.3 4.2 3.6 5.5 4.4 6.6 2.8 2.3 3.6 3.8 4.8 3.8 5.3 4.3
##  [433] 5.7 5.4 6.5 3.7 6.2 5.4 4.4 5.2 5.5 4.1 6.3 4.9 4.3 4.5 4.2 5.3 4.2 5.5
##  [451] 5.0 4.8 4.5 3.6 4.3 4.8 4.2 5.3 4.5 3.4 4.5 5.3 4.6 3.7 3.7 5.4 5.2 4.8
##  [469] 4.0 5.6 4.9 3.1 5.1 3.0 3.6 3.4 4.6 4.0 2.9 4.6 5.2 3.9 4.3 4.4 4.0 4.0
##  [487] 4.6 4.3 5.0 4.0 2.8 3.6 4.2 5.7 6.4 4.1 4.9 5.9 3.7 4.1 4.5 3.2 5.6 3.6
##  [505] 5.6 4.4 5.1 4.5 5.3 5.0 4.4 5.9 3.9 5.4 3.3 4.4 4.6 5.4 4.8 5.7 5.0 4.2
##  [523] 5.6 5.1 3.8 5.3 4.2 2.9 2.7 5.8 5.5 4.6 3.7 4.6 5.3 5.3 3.4 5.9 4.0 5.2
##  [541] 5.5 4.9 4.8 4.5 5.3 5.0 3.2 3.8 4.2 5.7 6.1 4.0 5.2 5.6 7.3 5.4 5.1 3.9
##  [559] 5.7 5.2 3.1 2.4 4.9 4.6 3.3 2.6 2.2 4.2 3.6 5.1 2.4 4.4 4.0 4.5 4.2 6.1
##  [577] 5.3 4.6 5.8 4.7 6.3 5.2 3.9 4.1 5.7 4.3 3.4 3.1 4.0 5.4 5.5 3.8 3.2 6.4
##  [595] 4.6 4.3 5.3 4.0 4.4 4.4 6.1 4.1 4.1 3.9 4.2 2.8 4.9 3.2 4.4 5.0 5.6 6.1
##  [613] 4.4 5.3 6.2 4.8 4.5 4.0 1.9 4.6 5.2 4.8 4.2 3.9 5.4 3.8 4.4 4.9 5.1 2.9
##  [631] 3.1 3.9 4.0 3.8 5.0 5.1 3.9 4.2 4.2 4.2 5.8 5.1 3.4 5.6 4.4 6.1 4.1 5.7
##  [649] 5.9 6.7 4.2 4.7 4.7 2.6 5.5 4.4 4.6 5.1 3.0 5.8 4.7 5.3 3.6 5.3 3.5 3.8
##  [667] 4.0 5.1 5.3 4.6 4.6 4.4 4.0 6.0 4.8 4.6 3.9 3.8 6.1 5.8 4.1 2.8 4.4 5.5
##  [685] 4.4 4.8 4.1 3.6 4.1 4.5 4.9 6.0 3.2 4.9 3.9 4.7 4.4 4.9 5.9 5.2 3.0 2.9
##  [703] 3.4 3.6 3.8 5.2 4.8 4.2 3.4 4.3 5.4 4.8 4.6 5.6 4.8 5.6 3.2 4.7 3.1 3.8
##  [721] 3.1 4.9 6.7 4.0 2.9 3.1 4.8 4.2 5.4 4.3 7.1 4.5 5.1 3.6 5.0 4.0 4.7 4.2
##  [739] 3.7 3.9 6.6 4.7 6.3 3.5 3.9 5.0 4.3 3.7 5.0 5.5 4.6 3.7 3.9 5.4 5.4 5.7
##  [757] 4.6 3.3 5.5 3.5 5.1 5.4 3.4 5.3 4.3 4.6 4.1 4.4 5.8 4.7 4.6 5.5 4.8 4.8
##  [775] 4.5 4.5 2.5 3.2 5.6 5.0 4.5 4.1 4.2 4.4 5.5 4.3 4.0 5.6 5.5 6.4 3.3 6.1
##  [793] 4.7 4.9 3.4 3.3 5.1 3.6 5.0 5.1 3.8 5.1 4.0 4.8 6.5 3.7 4.6 4.9 3.7 4.6
##  [811] 4.8 4.1 4.7 4.0 4.8 4.7 3.1 4.3 3.4 5.0 5.6 5.0 3.7 3.6 5.4 4.9 3.7 4.9
##  [829] 5.3 4.1 3.4 4.7 4.2 4.0 4.7 5.3 4.6 6.9 5.2 5.3 3.4 4.6 4.3 4.5 5.3 4.7
##  [847] 2.7 4.2 4.3 4.5 4.0 4.6 5.1 5.3 3.3 5.2 5.3 4.7 4.5 4.5 3.6 4.3 3.1 4.9
##  [865] 4.0 4.2 2.7 3.6 2.1 4.2 5.7 4.8 5.9 5.3 3.4 5.3 4.5 5.4 4.0 6.0 4.6 3.9
##  [883] 3.5 4.8 4.8 4.2 5.5 5.2 4.7 4.4 5.1 4.0 4.1 4.5 4.9 3.2 6.9 4.4 5.8 5.2
##  [901] 2.9 5.9 4.6 4.4 3.1 2.6 4.4 5.4 4.7 4.2 4.1 5.8 5.1 4.4 5.6 4.9 5.9 6.0
##  [919] 3.4 4.3 4.4 5.8 2.9 6.2 4.3 4.4 3.6 3.5 3.9 3.7 4.4 2.5 4.2 5.4 4.8 4.2
##  [937] 5.8 5.3 3.6 5.7 3.8 4.3 3.8 3.3 4.0 5.2 3.4 4.7 5.1 5.0 4.8 2.8 4.9 2.4
##  [955] 2.8 4.0 5.6 5.2 5.0 4.8 4.8 3.9 4.2 6.5 4.6 5.2 4.9 4.7 2.8 5.8 4.8 4.3
##  [973] 3.9 5.1 4.7 6.2 4.4 3.7 4.7 4.7 5.2 4.3 2.8 4.9 4.2 4.9 3.8 4.3 4.7 4.1
##  [991] 5.0 5.4 3.9 4.5 4.4 3.3 4.6 5.2 6.7 4.0
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
##    [1] 5.6 4.3 4.1 3.7 4.4 4.7 2.8 3.9 4.7 4.7 5.7 5.6 4.4 6.0 4.4 4.7 3.3 5.9
##   [19] 4.8 4.4 3.2 4.8 3.7 5.2 3.8 3.9 5.4 5.1 3.2 5.2 4.3 4.3 4.1 4.3 5.7 4.8
##   [37] 5.7 3.8 5.4 5.6 5.2 4.7 5.6 4.1 4.5 3.0 4.8 4.6 5.8 3.8 3.7 4.2 2.9 3.8
##   [55] 4.8 4.0 5.6 2.9 4.5 4.0 4.7 5.4 4.1 3.4 4.4 3.3 3.4 5.8 4.2 4.3 4.6 4.3
##   [73] 3.4 3.4 4.3 6.0 4.6 3.4 4.3 5.0 4.7 4.7 4.7 4.9 5.8 5.1 4.6 6.3 5.7 5.2
##   [91] 4.4 4.6 4.3 4.1 3.9 5.5 4.9 5.6 4.8 5.2 5.6 4.9 5.4 3.4 6.0 3.3 4.0 6.0
##  [109] 3.9 4.6 6.1 3.9 4.7 5.7 3.7 5.6 5.5 5.5 1.8 5.2 4.5 4.8 4.5 4.4 5.7 3.8
##  [127] 4.4 4.2 2.8 4.8 3.9 4.7 4.1 4.2 4.6 2.0 3.4 4.8 5.0 4.2 4.1 5.0 6.5 3.9
##  [145] 6.5 4.3 5.8 4.7 5.7 2.7 3.3 4.7 4.0 5.0 6.5 5.1 5.0 5.5 4.5 4.6 5.4 4.0
##  [163] 6.5 4.0 4.3 4.2 3.9 5.3 4.0 5.5 5.6 4.6 5.2 3.2 3.6 5.4 5.8 4.6 4.9 3.4
##  [181] 3.3 3.1 4.2 5.1 4.2 4.5 4.0 3.5 2.7 5.9 5.2 3.9 4.2 3.8 3.5 4.0 4.2 3.1
##  [199] 4.4 4.3 3.8 6.2 4.1 4.7 5.8 4.0 5.0 4.0 5.0 2.9 5.5 6.4 4.8 5.4 2.0 6.7
##  [217] 2.9 3.0 5.3 4.2 3.7 3.3 3.2 3.9 5.8 5.2 3.8 2.9 4.1 4.9 3.4 2.9 3.5 4.9
##  [235] 4.2 5.9 4.0 4.9 4.8 5.0 2.8 3.6 5.9 4.3 4.9 2.8 3.7 3.2 4.5 3.6 4.8 5.8
##  [253] 5.3 5.0 6.0 3.3 4.8 5.3 4.7 4.6 3.5 4.7 5.2 4.5 4.3 4.7 4.6 4.5 4.3 2.3
##  [271] 5.5 6.0 4.2 4.9 3.0 4.8 4.9 3.8 5.0 5.2 3.5 4.1 4.6 3.6 4.8 4.9 4.9 5.3
##  [289] 6.3 4.6 3.5 4.5 3.5 5.4 4.2 4.1 3.2 3.7 4.7 3.3 7.1 5.9 4.3 4.1 4.1 5.2
##  [307] 4.9 5.8 4.7 5.4 4.9 3.7 2.8 4.9 3.6 4.0 3.8 3.3 4.9 4.3 4.9 4.4 5.0 4.4
##  [325] 4.6 4.8 3.2 4.7 5.2 3.9 4.8 5.7 3.1 4.8 5.0 4.2 4.0 3.5 3.0 5.7 5.3 6.5
##  [343] 5.1 2.4 5.7 5.6 4.4 4.9 5.0 6.1 4.1 6.1 4.9 3.4 2.4 4.5 5.4 5.2 4.5 3.7
##  [361] 4.5 3.5 4.6 2.9 3.1 3.6 4.1 5.6 5.2 6.0 3.8 4.5 4.2 4.9 6.1 3.2 3.8 4.5
##  [379] 4.2 4.1 4.8 5.0 4.7 5.5 3.0 4.6 4.4 3.5 4.8 5.2 6.6 4.9 4.4 5.6 4.3 5.1
##  [397] 5.2 5.2 4.3 4.3 3.7 4.3 4.2 3.6 3.3 4.3 3.9 4.4 3.5 2.1 5.3 4.5 5.7 5.0
##  [415] 4.3 3.9 3.6 3.0 2.9 4.0 5.2 4.9 4.9 5.0 4.0 4.8 3.8 4.0 2.0 5.4 3.9 3.8
##  [433] 5.1 3.6 2.9 3.8 4.8 4.3 2.7 3.9 2.6 3.3 6.2 5.8 4.6 4.8 4.9 4.2 4.4 5.4
##  [451] 4.4 3.8 5.0 4.6 4.0 5.2 6.3 5.8 3.7 3.3 3.6 4.8 3.4 4.9 3.8 4.1 5.8 3.4
##  [469] 3.8 2.9 5.5 4.5 6.4 4.3 4.7 4.2 5.1 2.8 4.0 4.5 4.3 3.8 3.9 4.6 3.8 4.1
##  [487] 4.7 5.0 3.8 4.3 6.1 4.0 5.2 3.9 5.5 5.6 3.4 4.4 4.3 4.8 3.4 3.8 3.5 3.2
##  [505] 4.1 3.4 5.2 5.5 4.8 5.6 4.8 3.8 3.0 2.3 3.4 4.7 4.8 3.4 5.1 5.0 6.0 4.4
##  [523] 5.9 4.3 4.0 4.9 5.4 4.9 6.7 6.1 3.5 4.5 5.4 5.3 4.3 2.8 4.5 4.0 4.5 3.5
##  [541] 4.9 3.2 4.4 5.5 4.3 3.1 4.9 2.4 4.8 3.1 4.6 3.8 3.5 2.5 5.2 4.9 5.2 5.1
##  [559] 5.3 5.9 4.5 4.5 5.9 3.6 5.0 4.7 4.9 3.8 5.1 5.4 5.4 5.6 5.4 4.2 4.2 5.4
##  [577] 4.5 4.5 5.2 4.5 3.8 3.5 4.7 3.7 4.6 3.5 2.6 3.5 4.5 4.5 5.5 4.1 4.5 4.4
##  [595] 5.9 3.5 5.6 6.1 4.3 6.2 7.3 5.1 3.9 3.5 5.0 3.2 4.2 3.0 4.7 3.7 4.0 5.0
##  [613] 4.6 3.6 4.0 4.7 5.8 5.2 4.2 4.6 3.7 5.5 3.4 5.4 6.2 6.0 3.0 4.7 3.6 4.5
##  [631] 4.4 4.5 4.1 4.2 4.3 4.7 4.7 5.2 5.1 4.6 4.5 5.3 5.1 4.5 4.0 3.7 4.7 4.6
##  [649] 5.2 4.0 5.2 5.3 4.7 5.1 5.1 4.0 5.2 3.5 5.2 4.6 5.9 4.9 4.9 6.0 4.8 5.8
##  [667] 2.9 2.1 5.2 5.1 3.8 4.9 5.9 4.2 6.2 4.8 3.8 3.7 1.4 4.3 3.7 4.2 3.8 4.3
##  [685] 4.4 4.6 4.6 3.6 3.9 3.4 3.8 4.9 5.2 3.3 5.5 5.9 4.8 4.4 4.1 5.6 3.6 3.0
##  [703] 3.5 4.3 4.8 4.5 3.3 5.8 5.1 4.8 5.6 7.0 3.6 4.8 3.3 4.5 3.1 4.4 4.3 3.3
##  [721] 3.8 3.2 5.6 6.3 3.4 4.2 3.7 4.2 3.1 4.6 2.8 5.7 4.5 4.4 5.0 5.6 3.9 4.2
##  [739] 5.5 5.6 5.5 4.7 4.6 4.6 4.2 4.3 4.7 5.9 5.4 4.2 3.3 5.7 3.9 5.7 4.0 6.0
##  [757] 4.5 6.6 4.7 4.4 4.1 3.9 4.4 5.8 4.3 4.8 2.0 5.9 5.2 4.1 5.8 4.7 3.9 3.7
##  [775] 2.5 4.3 3.5 5.7 3.9 4.9 3.6 4.5 6.1 5.6 5.5 4.9 4.9 4.7 4.9 4.9 4.1 4.0
##  [793] 4.5 3.4 3.4 5.1 3.8 4.0 4.0 4.7 5.4 4.7 5.2 4.9 4.7 3.7 4.9 3.6 5.9 4.6
##  [811] 6.4 6.3 3.8 4.7 4.4 4.7 4.9 3.5 3.6 3.2 4.2 4.3 3.8 6.3 6.0 4.7 2.2 4.6
##  [829] 3.7 5.7 3.1 4.5 4.7 3.4 4.0 4.7 5.5 5.0 5.1 3.5 5.9 5.4 4.3 3.3 4.6 5.8
##  [847] 4.2 5.5 3.5 4.7 5.6 4.4 4.1 5.0 4.5 5.6 5.5 3.5 4.8 4.4 3.5 4.2 4.5 4.2
##  [865] 2.7 4.1 4.5 3.0 4.3 5.9 5.7 4.5 3.3 4.1 5.1 5.0 4.6 6.2 6.0 5.2 4.2 4.2
##  [883] 4.9 4.6 2.8 5.5 5.1 4.3 3.2 4.7 3.7 4.1 3.7 5.1 4.4 5.2 4.1 5.1 4.9 4.0
##  [901] 4.0 5.6 3.7 6.0 3.9 4.9 5.8 5.5 5.8 3.7 5.6 4.4 5.9 5.2 5.3 4.7 4.1 6.4
##  [919] 4.7 4.1 3.4 4.4 5.9 5.5 4.4 4.5 3.9 4.5 5.7 5.3 4.4 5.8 5.7 4.4 3.9 5.1
##  [937] 4.2 5.2 3.3 6.5 4.8 5.3 5.3 5.4 2.9 3.4 4.7 3.1 3.9 3.3 4.9 5.3 4.6 3.6
##  [955] 3.6 3.8 5.6 3.4 5.7 3.7 6.6 4.5 4.0 4.3 3.6 4.7 4.2 5.4 4.7 5.3 3.7 5.4
##  [973] 4.8 5.2 4.5 4.5 4.7 5.8 5.3 4.8 5.2 4.3 4.6 4.9 3.9 5.6 3.4 4.5 5.6 5.9
##  [991] 5.1 4.6 3.9 3.7 3.6 4.0 1.4 5.3 4.8 3.1
## 
## $func.thetastar
## [1] 0.0103
## 
## $jack.boot.val
##  [1]  0.49429429  0.41002786  0.29033233  0.18813559  0.04034091 -0.02083333
##  [7] -0.13083573 -0.24204204 -0.37019499 -0.50909091
## 
## $jack.boot.se
## [1] 0.9505149
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
##    [1] 4.0 3.3 3.0 3.1 5.2 4.5 4.7 5.3 3.3 3.5 5.3 4.6 4.3 5.3 4.8 4.3 2.4 4.3
##   [19] 4.1 6.0 5.1 5.3 6.0 5.0 3.7 5.8 4.0 5.1 4.7 4.6 4.3 4.4 5.8 4.6 2.8 4.4
##   [37] 3.9 3.3 5.4 5.1 5.2 4.9 4.2 4.0 5.0 5.0 2.9 2.7 4.3 3.1 3.4 6.4 5.3 7.1
##   [55] 4.7 4.6 3.5 6.6 5.3 5.0 4.8 3.3 6.0 3.5 3.1 3.1 4.5 5.7 4.9 3.4 5.9 4.6
##   [73] 3.8 1.8 3.1 5.2 3.7 3.8 3.1 5.3 3.5 3.0 3.4 5.1 5.9 3.9 4.6 5.4 4.7 4.7
##   [91] 5.3 5.5 3.2 5.5 4.6 5.7 3.8 2.7 3.3 6.5 6.3 4.5 4.9 4.3 6.8 4.1 3.5 5.5
##  [109] 4.2 3.8 4.3 4.8 2.7 3.1 3.6 3.3 5.7 5.4 3.8 4.3 4.3 5.2 4.9 6.1 4.8 3.6
##  [127] 6.4 4.5 4.4 4.2 3.4 4.2 5.1 4.2 4.9 5.4 3.0 4.5 5.0 6.0 3.4 2.9 3.1 4.4
##  [145] 5.0 4.7 5.1 4.6 4.1 5.1 3.7 3.3 4.8 5.4 3.1 4.4 4.2 3.9 6.1 4.6 4.2 5.0
##  [163] 5.4 3.2 5.0 4.2 4.4 2.0 4.4 4.8 4.4 4.5 3.6 5.3 3.5 4.6 3.0 3.9 4.7 5.1
##  [181] 3.9 5.3 5.7 4.9 3.5 4.5 3.9 4.6 3.4 4.8 5.5 3.6 4.3 5.6 5.4 5.3 7.3 4.4
##  [199] 4.7 5.5 4.9 4.7 5.1 3.6 4.8 5.4 6.0 3.9 4.4 5.1 5.5 4.9 5.6 5.3 5.0 4.7
##  [217] 6.3 3.2 4.4 5.9 4.5 3.9 3.0 4.5 5.1 4.0 4.4 5.2 6.1 5.4 4.1 4.9 3.9 2.3
##  [235] 3.2 4.6 3.8 4.6 4.5 2.4 3.6 5.2 4.1 3.7 5.0 4.1 4.1 3.5 4.5 5.1 5.1 4.5
##  [253] 4.1 3.9 5.1 4.8 4.1 4.6 3.9 3.0 3.4 6.0 6.0 5.5 5.6 5.2 4.9 5.0 5.3 5.5
##  [271] 5.4 4.2 4.2 4.0 4.0 3.7 4.3 4.0 3.8 4.3 6.0 4.0 4.4 4.1 5.0 4.3 4.6 4.3
##  [289] 2.3 3.7 2.4 3.5 4.9 4.6 3.6 5.7 4.5 5.4 3.7 4.4 5.2 3.8 5.0 2.4 5.5 4.5
##  [307] 4.6 5.3 5.2 4.5 4.1 4.5 4.8 3.4 4.2 4.8 5.0 4.8 4.0 4.7 4.3 4.2 3.6 5.4
##  [325] 5.8 4.3 4.0 3.3 4.5 3.4 3.5 2.8 6.0 4.4 4.6 4.4 5.6 4.8 3.3 5.6 2.8 3.3
##  [343] 3.9 5.0 3.6 4.1 5.6 4.5 4.5 3.7 3.0 5.5 3.7 3.8 5.0 5.8 4.9 4.0 5.4 4.4
##  [361] 4.1 2.8 3.5 5.1 4.2 4.1 3.2 3.9 2.7 3.3 4.9 6.2 4.5 3.7 3.8 3.8 4.8 3.2
##  [379] 5.4 4.4 4.7 4.4 4.9 6.4 5.0 3.3 5.0 4.0 4.9 4.7 4.3 4.6 4.5 4.1 4.6 3.7
##  [397] 5.7 3.1 3.6 3.6 4.0 4.3 3.8 5.0 3.1 6.1 3.7 6.0 4.9 5.8 5.7 3.5 6.8 4.2
##  [415] 5.3 4.8 5.4 5.2 3.2 2.7 3.1 4.8 3.6 3.2 2.3 4.1 4.9 4.7 4.1 4.8 3.6 3.6
##  [433] 5.4 3.6 6.1 5.2 6.0 3.9 4.2 4.3 5.4 5.5 4.3 4.8 3.0 4.6 5.4 4.5 4.3 4.4
##  [451] 3.9 4.2 4.9 4.2 4.0 3.9 4.0 4.2 6.2 4.6 3.5 3.7 4.8 3.8 4.3 3.1 4.0 3.7
##  [469] 4.9 3.3 6.1 4.9 4.1 3.9 3.8 5.7 5.7 3.7 5.1 3.9 4.2 4.5 5.4 3.8 5.9 5.8
##  [487] 4.9 2.7 5.4 5.1 3.7 3.3 5.4 3.3 5.5 4.2 5.1 4.5 5.0 6.2 3.3 6.2 5.4 5.5
##  [505] 4.0 3.5 3.9 3.9 3.4 4.8 4.0 3.7 4.6 4.3 4.0 3.7 4.0 3.6 4.1 4.7 4.3 4.9
##  [523] 3.5 4.6 4.8 4.1 5.0 4.0 5.2 5.0 5.3 4.0 4.8 4.8 2.9 4.7 3.8 5.0 4.1 5.5
##  [541] 5.2 4.1 4.9 3.3 4.7 3.8 6.1 3.7 4.3 4.6 5.7 4.8 5.1 5.5 4.8 4.8 4.3 3.2
##  [559] 5.8 5.9 3.8 4.3 4.4 4.0 4.4 3.5 5.0 4.1 4.0 4.5 5.0 4.9 6.1 3.8 5.9 4.0
##  [577] 4.8 5.0 4.0 4.2 4.3 4.5 4.4 3.7 5.3 4.5 5.1 4.5 4.8 4.2 3.8 4.0 4.2 5.2
##  [595] 4.4 5.2 4.0 5.0 3.5 5.7 3.9 4.0 5.0 4.0 6.0 4.2 3.5 4.3 5.4 4.9 4.5 5.1
##  [613] 5.2 5.9 4.8 5.0 5.1 4.5 6.6 4.2 5.1 3.4 4.5 4.3 4.9 2.8 3.2 5.1 3.8 4.0
##  [631] 4.5 4.8 3.8 6.2 5.0 4.6 4.2 5.0 5.9 4.4 4.6 4.7 4.3 3.5 5.6 5.3 5.5 4.0
##  [649] 4.7 5.1 4.8 5.8 5.9 4.5 4.8 5.6 3.7 3.2 4.4 4.6 5.4 3.5 4.0 3.6 4.2 3.6
##  [667] 4.4 5.2 3.9 4.7 4.5 4.6 3.4 4.6 5.0 3.2 4.4 5.1 4.7 4.5 4.0 4.4 4.7 2.8
##  [685] 3.9 2.9 6.0 3.9 5.2 3.9 4.1 3.6 5.3 3.5 5.0 4.7 3.8 4.7 3.3 5.3 4.8 4.7
##  [703] 5.5 3.6 6.4 4.2 4.3 4.2 4.3 4.1 4.9 5.9 4.0 5.0 5.4 2.7 4.1 3.2 4.8 2.8
##  [721] 3.9 4.3 3.9 3.1 4.5 4.0 6.4 4.8 5.8 4.5 4.0 4.9 4.2 5.2 4.3 4.9 3.1 5.1
##  [739] 2.8 2.8 5.2 3.9 4.4 5.2 6.6 4.3 4.1 6.0 4.6 4.0 3.2 4.2 4.6 3.1 2.8 6.4
##  [757] 3.9 3.1 4.0 4.9 4.2 3.7 4.5 2.9 4.4 3.7 4.6 4.4 3.3 4.0 5.0 4.5 4.5 3.4
##  [775] 5.2 3.9 4.0 6.4 3.9 3.3 5.0 4.0 3.7 4.8 4.1 3.4 5.4 5.0 5.2 5.8 5.6 2.9
##  [793] 6.1 4.2 4.0 4.7 4.7 3.8 5.0 5.8 5.0 5.7 6.7 4.6 5.1 5.3 4.9 4.9 5.9 5.6
##  [811] 4.7 3.9 3.9 6.2 6.2 5.8 3.0 5.6 4.7 5.0 5.0 4.7 6.2 5.6 5.4 5.8 4.3 4.4
##  [829] 4.0 4.5 4.8 5.9 6.8 5.5 4.3 4.8 3.4 5.3 4.9 2.5 4.3 3.7 4.6 5.3 2.8 5.2
##  [847] 5.6 5.9 3.6 6.0 5.1 4.3 5.3 4.6 4.0 6.1 3.3 4.7 3.3 6.2 4.0 5.1 4.7 5.8
##  [865] 6.2 4.8 5.1 5.7 4.1 4.5 5.1 3.8 4.6 4.6 5.1 4.9 4.2 4.1 3.7 1.3 5.1 4.5
##  [883] 4.2 4.2 3.8 4.1 3.7 4.5 4.9 4.5 2.9 2.9 4.9 4.3 3.9 2.9 5.2 2.1 4.5 6.3
##  [901] 5.1 5.8 4.6 4.5 6.1 5.4 3.6 4.1 3.2 4.7 3.7 5.3 6.5 5.2 4.8 5.7 3.6 5.5
##  [919] 6.8 4.0 5.1 4.2 5.6 5.1 3.4 4.3 3.1 4.5 6.0 4.2 5.2 4.3 2.6 4.6 4.2 5.0
##  [937] 6.7 2.8 3.5 4.7 5.3 3.4 4.2 6.5 5.1 5.2 5.2 4.5 5.6 5.5 4.2 2.8 6.0 5.6
##  [955] 4.0 4.1 3.4 5.0 4.1 5.2 6.7 3.1 5.6 4.5 4.2 4.5 5.8 4.3 4.5 5.8 3.8 4.7
##  [973] 4.6 4.3 3.4 4.3 5.0 2.8 2.9 3.7 5.7 3.6 5.9 5.5 3.8 5.5 4.0 3.3 5.3 4.4
##  [991] 3.2 4.7 3.4 5.0 5.3 6.3 6.0 4.1 4.0 4.0
## 
## $func.thetastar
## 72% 
##   5 
## 
## $jack.boot.val
##  [1] 5.4 5.5 5.3 5.2 5.1 5.1 5.0 4.8 4.6 4.5
## 
## $jack.boot.se
## [1] 0.9415413
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
## [1] 0.5299747
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
##   4.509930   6.247899 
##  (1.946896) (2.853110)
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
## [1]  0.8906622  0.7561109  0.4587840  0.7892122  0.0327753 -0.5414604
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
##    [1]  0.3026304764  0.4660222673  0.3779783227 -1.2160256978  0.7143272334
##    [6] -0.3237572854  0.9917713954  0.1711315264  1.5287681224 -0.5623383449
##   [11] -0.1225900333  0.3337376509  0.5073223652  0.7500411509  1.3626839860
##   [16]  0.1042810658  1.0967430472  1.1290786113  0.4201366461 -0.2535688590
##   [21]  0.7040814070  0.4230827630  0.2622663894 -0.2191823119 -0.5665646259
##   [26] -0.9074806965  0.8412749944  1.0384126908  0.8152143528  0.4240081679
##   [31]  0.8231825310  2.2593378355  0.2321163850  0.6587490916  0.1445136665
##   [36]  0.4449347949  0.5731919020  0.6707150357  0.3785799156  0.5934549744
##   [41]  0.6230470677 -0.0065419949 -0.9334108026  0.4611721461  0.1271929440
##   [46]  0.3278843385  0.6686875718  1.7039247081 -0.3145536363  0.9523961437
##   [51]  0.4861768325  0.9844904874 -0.2254718535  0.1082466150  0.4680298884
##   [56]  0.7449197642 -0.6984498169  0.8525844660 -0.3508493233  0.1829524068
##   [61]  1.2298418310 -0.1754310023 -0.1220693744  1.0309563959  1.4539813153
##   [66]  0.7954922802  0.5120300717  0.4278006611  1.9601361808  1.8884850978
##   [71] -0.3650854586 -0.2033027965 -0.4425669353  0.9811689788 -0.6984625994
##   [76]  0.0361961777  0.8920670346 -0.8519448588 -0.4741780748 -0.4473480056
##   [81]  0.0533929423  1.0837913513  0.4957241616  0.6182930990  0.7779838711
##   [86]  0.6943569261 -0.5160117466  0.3569116257  0.0119989898  0.3004771118
##   [91] -0.2398751484  0.5801443375  0.1367660170  0.9483543638  0.7742855622
##   [96]  0.5746644784  0.7465493240  0.5273759847  1.5766156021  1.5322104793
##  [101] -0.2087595351  0.3724951879 -0.9947720779 -0.5379033080  0.8994584926
##  [106]  1.1825147599  0.4581439984  0.0007804637  0.6832670253  0.6269883055
##  [111]  0.2705375212 -0.0587219667  0.8777084977  0.5460085808  0.2589513287
##  [116]  0.3382606092 -0.0637739890  1.7247989851  0.6920454679 -0.2191162261
##  [121] -0.0865471440 -0.2046652378 -0.2130544256  0.2421256019  1.2386121637
##  [126]  0.3389844574  0.5862520735  1.1251773127  2.0250143358  0.6622521199
##  [131]  1.4742989902  0.2635621286  0.3259710468  0.2856126849  0.1073571287
##  [136] -0.8099669217 -0.0609340794  0.1729955591  0.6686915860  0.5858427399
##  [141]  0.4465336660 -0.1001874132  1.1912338860  1.2259858793  0.9588082618
##  [146]  0.6500613893  0.6906712085  0.3597099189  0.4723147794 -0.3617524574
##  [151]  1.1950033578  1.4286395014  0.1828046147  0.3853803556  0.1954509353
##  [156]  0.1317386438  0.5647241356  0.8344299826  0.6698211182  0.0309401986
##  [161]  0.6637952178 -0.0462452059 -0.8139573068  1.0023608182 -0.0959561939
##  [166]  1.1236371283  0.4179479453 -0.3672143486  0.8775976143  0.6316997467
##  [171]  0.6272117499  0.2197577167  1.1945001550 -0.0497438152 -0.3971431871
##  [176] -0.8146348576  1.4286875126 -1.0132065169 -0.9055683268  0.2069263133
##  [181]  1.1640526564  0.3621606412  1.3108849795  0.8764285786  0.6902252787
##  [186] -0.4682266720 -0.5376941451 -0.2115761917  0.6127813549  0.6490118520
##  [191]  0.6677570764 -0.1519334963  0.8067229291  0.4606658011  0.2832456619
##  [196]  0.8107910012  1.4093626275  0.9235978952  1.1156323338  0.4203709539
##  [201] -0.1581247844  0.1897001187  0.2638150454  0.4134756508 -0.6335312628
##  [206]  0.6077857733 -1.8917905096 -0.8637583514  0.6709402066  0.8431141602
##  [211]  0.7822207708  0.2969315324  0.9216079544  0.9167103623 -0.8801290895
##  [216]  0.5906545847  0.4405365777  1.0440569830  0.7619021184  0.6301877883
##  [221]  0.3867078378  0.5035305565  0.0999107171  0.8094965566  0.3435559669
##  [226]  0.2672158037  0.1153277809  0.2404209459  0.4995026987  1.7165490900
##  [231] -0.0498390569  1.3834321384 -0.3096953278 -0.2090995636 -0.5673239091
##  [236] -0.7687429489  1.0144752167  1.3867075877 -0.3623137027  0.5289975858
##  [241] -1.8806829991 -0.0077938142  0.2788644837  0.8270449926 -0.5945731971
##  [246]  1.1879985416  0.8057791976  1.6705812318  0.4916760905  0.4712139018
##  [251] -0.0782273890 -0.4070657047  0.5652900442 -0.3233344279  0.1345157106
##  [256]  0.7882118531  0.8174968169  1.0283003271  0.7242469461  0.4643419017
##  [261]  1.0168107456  0.5726375185 -0.0417482029  0.6900936068  0.3807080585
##  [266]  0.0814815232  0.6518582282  0.6865175301  1.0412287758  0.7145968345
##  [271]  0.2670999196  0.4875482271 -0.1269487087  0.3557901174 -0.4507725783
##  [276]  0.1108081026  1.4030024219  0.2724051351  0.7788774182  0.4541410763
##  [281] -0.5681101045 -0.5142449532 -0.1930784776  1.4805206626  0.6770789711
##  [286]  0.1980614424 -0.1038097318  1.4959352964  1.2410381459  0.5501448788
##  [291] -0.1951528537  0.8479985662  0.6101984327  0.2715919268  0.1360787636
##  [296]  0.7056898932  1.0141869096  0.6403732620 -0.8568249993  0.5806338910
##  [301]  1.2061644125  0.2176644849 -0.0179825138  0.2696909378  0.8920194427
##  [306]  0.4227701472  0.5647241356  0.5046024512  0.9588082618  0.0233488628
##  [311]  1.1962168475  0.4108671529  0.8047940849  0.4841951016  0.5957854088
##  [316]  0.1404733851  0.5551227164  1.1901584451  0.5875990018 -0.3970680294
##  [321] -0.7502657915  0.6741210874 -0.4334093730  0.4794389386 -0.2214424301
##  [326]  0.3038378812 -0.5139359159  0.1613544250 -0.2364895923  0.6484562342
##  [331]  1.0288832193 -0.5467220735  0.8575642181  0.0520800699  0.5373657901
##  [336]  0.9810627185  0.6338075583  0.7239509675  1.0790594573 -0.5608144326
##  [341]  0.7222107614  0.4613340005  0.6758807493  0.1855016877  0.8705810791
##  [346] -0.6534001663  0.4864276287 -0.2146944608  0.8722301142 -0.0940134293
##  [351]  1.1826853358 -1.0669206467 -0.1149938846 -1.0768049324  0.1587855926
##  [356]  0.1394706545 -0.0988273325  1.3562436177  1.2000273831  0.0962004940
##  [361] -0.0264471889  0.0423843104  0.3618435300 -0.7642759294  0.6567081245
##  [366]  0.5252285869  0.4379950383  0.8239209290  0.4919914312  0.5681484896
##  [371]  0.4260748828  0.1975780248  0.1369096941  1.4649721861  1.0065471748
##  [376]  0.5585320654  0.0131727824  1.0855594986  0.1659948758 -0.1508349757
##  [381]  0.5225657069 -0.6500095086  0.4074393589  0.1631100620 -0.5833966855
##  [386]  0.2264941675  1.3797302321  0.4364021498  0.6127813549  0.6549614974
##  [391]  0.3927596272  2.0014878148  0.8378731708  0.4618902391  0.4676099671
##  [396] -0.0630502053  0.1020520691  0.9521054511  0.8938525147  1.1776156369
##  [401]  0.3228703138  1.6117269996  0.3123651321  0.9024297851  0.7220713860
##  [406] -0.2336245599  0.0590933485  0.9430986051  0.6416641609 -0.8345995705
##  [411] -0.3591373862  0.2837015777  0.4752656446  0.4729069193 -0.3869420880
##  [416]  0.6021049845 -1.0669206467  0.7576399811  0.4211298493 -0.0341355593
##  [421]  0.6999832078  1.1671622062 -0.7811197890 -0.2174187355 -0.7139043566
##  [426]  1.6564298460  0.6056155065 -0.1237386708 -0.3936247250 -0.8626995430
##  [431] -0.7450658026  0.7957919933  0.9725000270  0.3451529182  0.2757806695
##  [436]  0.3562938771 -0.3682334694  0.3437347210  0.4505130736 -0.9119257386
##  [441] -0.3001522746  1.3147544406  0.4517519164  0.4672320864 -0.6648970210
##  [446]  0.1808528773  0.4144962689  0.1253752254  0.4598093189  1.6437152783
##  [451]  0.5201904151  0.2544182867  0.6682670173  0.5896834357  0.3249474607
##  [456] -0.9379563834  0.0123913113 -0.1997431981  1.3100519932  0.0575362912
##  [461]  1.1383783889  1.0009951203  1.4538456441  0.0594103678  1.2190974024
##  [466]  1.1158081853  0.9425720188 -0.1566971592 -0.5406684084  0.7449104999
##  [471]  1.1940080327  0.2705375212  1.6432355036  0.4187947111 -0.7739781853
##  [476]  0.5320979811  0.8221297097  0.1965546642 -0.3166971175  1.1598581312
##  [481]  1.6439680710  0.6037646948 -0.1000263213 -0.3802737516  0.4606168394
##  [486]  0.2223041330  1.2748784373 -0.2791029313  1.2077951378  0.3184361706
##  [491] -0.3135851054  0.2754180454 -0.6354126686  0.6351309265  0.8385249267
##  [496]  0.1520549235  1.1880695142  0.2381488219  2.2116269387 -0.4456750648
##  [501]  0.1116422602  0.6932435465  1.3740120888  0.3049641578  0.6869927157
##  [506] -0.2628495555  0.4329113118 -0.3624689088  0.0424071731 -0.2500915920
##  [511]  0.1668498232  0.8263866750 -0.3674749321  0.8241866858  0.8192386304
##  [516]  1.0922382907  1.6178348808  1.3560510201  0.7842111995  0.3182129830
##  [521]  0.9401835171  0.4187977804  1.2410511580 -0.9933114618  0.8763411242
##  [526]  0.3335032071 -0.0697044619  0.2675319808  0.8675738664  0.6588405620
##  [531]  0.6914402175 -1.1212432563  0.9515042022  0.7788774182  0.4944123100
##  [536]  0.8382116867  0.9670484085  0.4483347094 -0.2631231524 -0.7793761570
##  [541]  1.0425723700  0.6889420190  0.7882441210  0.9661354478  1.2225431788
##  [546]  0.9106210982  0.5429240223  0.2087158091 -0.4988523643 -0.1269398955
##  [551] -0.5497654813 -0.2588701927  0.3015160670 -0.6166195616  0.9861354640
##  [556]  0.1925392544 -0.4283881732 -1.0176496962  0.7391446097  0.6303847455
##  [561]  2.4333344769  0.2675319808 -0.5453873481  0.2896900286  2.3783698872
##  [566] -0.8500682823  1.6322732556 -0.0741239245  0.4711701053  0.9251364197
##  [571]  0.3367583910  0.9672662529  0.5148107876 -0.1706528708 -0.3564740752
##  [576]  0.8966098502  0.2161078161  0.5328796344 -0.7218793299  1.3528133561
##  [581] -0.4295440038  0.2129215276  0.0291539823 -0.9952761953  0.4755727063
##  [586]  2.0025832185  0.2313396824  0.4187855788 -0.0799440790  0.6272117499
##  [591] -0.1550051557  1.2076501204  0.1640874986 -1.1948276577  0.7211748286
##  [596]  1.0667647790  0.5504462971 -0.1938040201  0.4678284568  1.1101192165
##  [601]  0.4872708720  0.1381244796  0.4034868000 -0.0576800506  0.3889356521
##  [606]  0.7498610012  0.3067014747 -0.0642166722  0.8507081536  0.7201127895
##  [611]  0.0031930930  1.0097118887 -0.2444766135  0.4339358604  0.0886812355
##  [616]  0.5318269248  0.1739542582  0.7190898230  0.4730953662  1.4597929753
##  [621]  0.8908657889  0.4872708720 -0.2468055321  0.1728886108  0.4740184911
##  [626]  0.2360866350  0.7068186578  0.9123887914  0.6145140682  1.3383754840
##  [631]  1.0095195507 -0.0124606717  0.9727380937  0.2441925323 -0.9889544302
##  [636] -0.0207663945  0.8120893394  0.1406913580  0.3002835557  0.0667759594
##  [641]  1.1797682250  0.8881905676 -0.0737459569  0.7684189649  0.8837029789
##  [646]  0.2273412725 -0.1190618501  0.8172589178  0.6975143558  0.5822102027
##  [651] -0.0589418585  0.6853122363  0.1559041483  0.7830671897  0.0441271276
##  [656]  0.5682217878  0.9479028617  0.6307735720  0.2398194086  0.3601583784
##  [661]  0.6667406582  0.9235800977  1.3201749793 -0.8840843535 -0.5973892274
##  [666]  0.9311132707  1.2320337847  0.1040659847  0.2266618033  0.8371665111
##  [671]  0.1853452982 -0.4146489524 -0.0753886296  0.7814967471  1.3461472091
##  [676]  0.9676268205  1.0775681992  0.1924568416  0.8557221945  0.1031914192
##  [681]  0.8384614634  0.6702899565 -0.8739411177 -0.0621360234  0.5559484302
##  [686]  0.0839457108 -0.7022475728  0.2252752586  0.0899690179  0.5884270670
##  [691]  1.1006533829  1.0024169922  0.5398252819  0.1334431533  0.2193550756
##  [696]  1.4218142832 -0.3778360424  0.5175942366  1.0654855276  0.4914858294
##  [701] -0.6139672223 -0.5748019043 -0.2144056928  0.5693524132  0.6632181236
##  [706]  0.9878992777  0.3924699332 -0.3022709379 -0.5298012861 -0.2468436717
##  [711]  0.8439942476  0.6239444736  0.8050975517 -0.2526765410  0.5500075751
##  [716]  0.8658612060  0.1735915075  0.9819148852 -0.1023618058 -0.6523828509
##  [721]  1.1014464041 -0.3623112043  1.0897802058  1.3206446457  1.1313915503
##  [726]  0.2829232439  0.2527608906  0.4581427730  0.4463163745  0.8840061084
##  [731]  0.8073568166  0.9390387983  1.3022451742  0.6419015883  0.2555779451
##  [736]  1.3711396071  0.1550889459 -0.4846227301  0.8330917137  0.5269746711
##  [741]  0.9392044187  0.2828482616  0.4552139751  0.2747380989 -0.2303516636
##  [746] -0.2983820525 -0.6156918779  0.1306430806  0.0364930829  0.4797609447
##  [751] -0.0725921749  0.4211298493 -0.7282968149  1.6691352627  1.1081703954
##  [756]  0.5156663974 -0.3737104526  0.2563459525  0.5543180504 -0.7162553113
##  [761]  1.4613834538  1.0060839660  0.3199270565  0.3924751624  1.0031166171
##  [766]  0.0469743008 -0.9730669138 -1.3877426225  0.2530068956  0.3007109259
##  [771]  0.8382116867 -0.2454754994  1.3691909378 -1.4032482268  0.0220017000
##  [776]  0.6778174180  0.5937586212 -0.3674207985  0.9779915846  0.1015795769
##  [781] -0.2655837506  0.2816153841 -1.0943482354 -0.0883118307 -0.9811443857
##  [786]  0.6310661894 -0.7312571598  0.0308682680  1.7225689323 -0.0257043283
##  [791] -0.1562871753  0.4350191260 -0.7032104426  1.3387208713 -0.5288550112
##  [796]  0.1619127842  0.4797980828 -0.3180215407  1.0591452092  0.5147619030
##  [801] -0.2404888039  0.3061642977  0.5611312093  0.9805664344  1.4579294878
##  [806]  1.3555982629  0.2791129817  0.1006712868  0.5874702907  0.3556444563
##  [811]  0.6147089719 -0.2801328790  0.7844864657 -0.0605693640 -0.2898799124
##  [816] -0.2631231524  1.1428262425  0.1337523933 -0.4416548392  0.5208339272
##  [821]  0.4055871952 -0.5382646582  1.0426157526  0.0286009946  0.1753251340
##  [826] -0.0234782887  1.3751929034  0.3104280232 -0.3081737425  0.8490512731
##  [831]  1.0000409817 -0.2376898480  1.5809302995  1.2457249660  0.5832613699
##  [836]  1.0621254792  0.5566676963  0.1812965266  0.7286363557 -0.1145056987
##  [841]  1.2707008538  0.6755599264  0.6011026788 -0.0763104293  0.4189925524
##  [846]  0.7161331746  0.6214411954  0.8737345652  0.4012218300 -0.1979445334
##  [851] -0.5094729759 -0.8767456301 -0.1326258583  1.3880983175  0.0918907497
##  [856]  0.3602606042  0.5989846826 -0.4526005917  0.2472548114  1.0982766508
##  [861]  0.9193122233  1.0658912961  0.6519627726  0.0063773528  0.5427132424
##  [866]  0.2278172202  0.2171173766  0.3720957264 -0.7705773798  1.3275519888
##  [871]  0.0842659274 -0.3021502164 -0.3415955761  0.2076474680 -0.5444333892
##  [876]  0.3086850932  1.5172064860  0.2109748440  0.4893677036  0.4175346852
##  [881] -0.2536571549  0.2980464818 -0.7048085307  0.8753278031  0.6296042736
##  [886]  0.9233584634  0.4350389904  0.6242901697  0.2335526880  0.0705974227
##  [891]  0.2142713465 -0.5359317079  1.1106787334  0.9399793150  0.1036101844
##  [896]  1.4947771887  0.4003014897  0.7103874048  1.7546272749  0.5485339399
##  [901] -0.6750566257  0.9175945542  0.9234479322  0.7075896972 -1.0139446691
##  [906]  0.5208419714  0.0058993750  0.0723846968 -0.9407380583 -0.1482433881
##  [911] -0.9905757097  0.3755554559  0.2496414395  1.5772612895  0.8738695980
##  [916]  0.6850511674 -0.1225011025 -0.1515237072  0.1904318172 -0.2461322448
##  [921]  0.3815000215 -0.6147803120  0.4661911514  0.3602606042  0.1160133716
##  [926]  0.0038915695  0.8473929254 -1.6085700273 -0.5285146681  0.1838342301
##  [931] -0.1617583669  0.1334236536  0.6189533850  0.5406046489  0.4682371123
##  [936]  1.2047271028 -0.0716846496  0.9705196867  0.7157051324 -0.6437530040
##  [941]  0.5679929587  1.0594312752  0.8230807695  0.2484736548  0.4576444954
##  [946]  0.1046377860  0.3036386455 -0.5258090608 -0.2061557422  0.1658143965
##  [951]  0.5548631716 -0.0757268927  0.4447554544  0.1311483904  1.0012310621
##  [956]  0.3155040646  0.2337266328 -0.3439419484 -0.1489656558  0.1878392135
##  [961] -1.1555754488  1.1826853358  0.2070663060 -0.1775942370  0.0290455448
##  [966] -0.5981543237 -0.3673707867  0.2261846362  0.4127315030  0.2961495248
##  [971] -0.2225154878 -0.3042660855  1.0707047671  1.2406230751  1.0161598516
##  [976]  0.5539734995  0.0665426879  0.2671201510  1.2941357590 -0.9721489364
##  [981]  0.4221801761  1.2788646541 -0.5798666375  0.0230738336  0.3785799156
##  [986]  1.3372995576  0.2415744057  1.3143690062 -0.1420607377  0.8855323902
##  [991]  1.4233707381  0.3624662160  0.1013615578  0.8031169254  1.0385787767
##  [996]  0.7642338550  0.4229759562 -0.3037301481 -0.1949418807 -1.0323328892
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
##   0.72184682   0.31728939 
##  (0.10033572) (0.07094532)
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
## [1] -0.1740450 -0.6967422  0.2809919 -0.4965089  0.1760012  0.7751387
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
## [1] 0.0221
```

```r
sd(results2$thetastar)  #the standard deviation of the theta stars is the SE of the statistic (in this case, the mean)
```

```
## [1] 0.9099833
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
##     original       bias    std. error
## t1*      4.5 -0.004904905   0.9203696
```

One of the main advantages to the 'boot' package over the 'bootstrap' package is the nicer formatting of the output.

Going back to our original code, lets see how we could reproduce all of these numbers:


```r
table(sample(x,size=length(x),replace=T))
```

```
## 
## 1 2 3 5 6 7 8 9 
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
## [1] 0.02
```

```r
se.boot
```

```
## [1] 0.9252654
```

Why do our numbers not agree exactly with those of the boot package? This is because our estimates of bias and standard error are just estimates, and they carry with them their own uncertainties. That is one of the reasons we might bother doing jackknife-after-bootstrap.

The 'boot' package has a LOT of functionality. If we have time, we will come back to some of these more complex functions later in the semester as we cover topics like regression and glm.

