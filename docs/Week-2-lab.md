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
## 0 1 2 3 4 7 9 
## 1 3 2 1 1 1 1
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
## [1] -0.0417
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
## [1] 2.698511
```

```r
UL.boot
```

```
## [1] 6.218089
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
##    [1] 5.4 4.0 5.4 4.4 3.8 4.1 4.4 4.3 3.2 4.1 3.4 4.0 4.8 4.9 4.5 4.7 5.3 4.7
##   [19] 4.5 4.6 4.8 4.3 4.1 4.7 4.2 4.2 3.2 3.0 4.6 5.4 3.9 5.9 3.6 5.4 5.7 5.2
##   [37] 4.8 4.8 4.5 2.6 4.3 5.7 4.1 4.8 4.9 3.5 4.6 4.8 5.3 5.1 2.7 2.6 4.0 4.1
##   [55] 5.0 5.6 6.3 3.1 2.1 4.3 3.7 3.1 4.5 4.5 4.4 4.5 3.6 5.1 4.7 3.9 4.5 4.4
##   [73] 5.8 3.2 4.6 4.6 4.7 6.3 4.1 3.2 3.6 4.6 3.6 3.3 4.0 6.8 3.6 4.4 4.4 5.5
##   [91] 4.5 4.2 6.2 5.3 3.5 5.4 3.7 2.6 4.6 3.0 3.4 5.2 3.3 5.9 5.1 3.9 5.0 3.9
##  [109] 4.4 4.0 5.4 4.5 3.5 4.6 3.4 4.1 5.2 5.0 5.2 4.1 4.4 5.1 4.2 4.8 4.1 4.6
##  [127] 5.3 6.7 4.5 3.4 4.3 5.4 3.6 5.0 4.4 3.2 4.1 4.2 6.1 5.9 6.1 4.0 3.6 3.8
##  [145] 4.3 4.4 6.1 4.3 4.4 4.2 4.0 4.8 5.0 4.2 4.7 4.7 3.7 4.5 4.3 3.9 4.3 4.0
##  [163] 3.3 3.3 4.4 6.3 5.0 4.6 5.7 7.2 4.7 4.1 6.2 4.3 4.4 3.7 3.4 5.5 5.4 3.8
##  [181] 5.6 3.8 4.5 3.8 4.9 5.5 5.1 4.2 5.5 2.1 4.2 4.2 5.2 4.1 4.4 2.4 4.7 3.3
##  [199] 3.2 2.9 4.7 4.8 4.7 5.4 4.4 4.9 5.9 4.8 5.1 3.0 4.0 3.9 4.6 3.6 5.3 4.4
##  [217] 4.6 3.9 4.3 4.7 4.6 3.3 3.4 4.2 2.7 4.7 4.2 4.0 5.4 3.7 2.9 4.0 3.7 4.7
##  [235] 4.1 3.5 4.5 5.0 3.5 5.2 4.0 5.0 4.1 4.9 4.0 5.4 3.7 5.6 4.3 4.7 5.0 3.9
##  [253] 4.2 5.1 3.8 3.5 2.2 4.8 3.3 2.9 4.9 4.9 3.9 4.6 3.4 3.7 4.6 3.7 4.4 4.8
##  [271] 5.5 5.0 4.5 4.5 4.0 4.2 5.0 4.7 5.1 3.5 2.9 5.6 5.4 5.6 4.0 6.1 4.3 3.4
##  [289] 3.1 5.0 6.1 3.7 4.1 6.0 3.5 5.8 3.2 5.1 4.7 4.2 3.8 4.5 4.5 4.2 5.6 5.5
##  [307] 3.7 5.4 4.9 3.7 2.8 3.8 4.6 5.1 4.8 4.7 6.3 4.3 4.6 3.7 5.7 4.5 4.2 4.2
##  [325] 4.6 3.9 3.7 5.0 3.6 4.7 5.2 4.7 3.9 5.8 4.4 6.5 6.1 3.8 4.5 4.7 4.6 4.3
##  [343] 4.6 4.0 2.7 4.5 4.6 3.2 3.7 4.4 5.5 3.8 3.9 4.7 4.6 2.7 4.9 4.4 4.5 5.1
##  [361] 4.9 5.9 5.0 4.2 3.6 3.0 4.9 3.9 3.8 4.6 3.5 5.5 5.3 4.8 2.5 3.7 4.3 5.1
##  [379] 4.8 4.8 2.7 5.7 4.8 5.1 5.5 5.9 4.8 4.8 5.6 6.0 4.4 6.1 5.2 3.3 5.4 6.0
##  [397] 4.2 5.2 5.5 3.6 4.2 3.2 3.4 3.6 3.6 3.9 4.0 5.4 4.4 4.8 3.1 4.4 3.8 4.9
##  [415] 5.0 5.3 2.9 4.2 3.5 4.8 5.5 5.3 2.7 4.3 3.5 4.2 4.4 5.0 4.4 4.4 5.5 3.5
##  [433] 3.0 5.5 3.5 5.8 3.7 5.2 5.2 4.9 4.2 4.7 5.5 3.4 4.3 4.0 3.1 2.8 4.9 5.2
##  [451] 6.8 5.3 2.9 5.5 2.2 5.0 4.4 3.6 4.6 6.1 3.3 6.6 4.6 4.2 3.1 4.3 5.1 3.9
##  [469] 4.4 3.6 5.1 6.2 3.7 5.2 3.2 3.7 4.5 4.6 4.5 3.2 4.5 5.2 5.1 3.2 5.0 3.4
##  [487] 4.7 4.7 5.9 5.1 5.4 5.8 5.0 4.1 4.5 3.7 5.9 4.0 5.5 4.4 4.8 4.2 4.5 4.0
##  [505] 2.3 5.0 5.4 3.7 5.0 4.8 4.8 4.9 4.6 4.0 5.6 5.3 6.9 5.3 3.3 6.0 4.5 3.3
##  [523] 4.3 4.5 4.5 4.2 5.3 3.8 5.1 5.3 3.7 5.3 4.0 4.2 5.9 5.3 4.7 4.4 6.0 3.8
##  [541] 2.8 5.1 4.5 5.7 3.0 5.7 6.0 4.7 4.8 3.9 5.0 4.6 3.2 4.9 5.5 4.5 3.4 4.1
##  [559] 3.3 4.5 4.0 3.6 5.0 5.4 4.1 3.1 4.1 2.9 3.1 3.7 4.5 3.7 4.3 6.1 5.3 5.0
##  [577] 4.8 4.0 3.6 5.4 4.8 3.4 5.8 7.3 3.4 4.9 3.5 6.0 5.5 3.5 4.9 4.7 5.8 4.5
##  [595] 4.5 4.9 5.0 3.8 5.8 4.0 4.4 4.0 4.0 5.3 5.4 5.7 3.0 4.3 4.5 3.4 3.0 3.0
##  [613] 3.6 3.7 5.1 4.9 5.4 5.2 2.2 3.9 6.0 5.4 2.2 3.9 3.7 3.7 5.0 4.8 3.4 4.9
##  [631] 4.2 5.4 4.8 5.7 4.0 5.4 4.8 5.0 5.0 3.7 3.5 3.1 3.5 3.5 4.2 2.6 4.0 5.5
##  [649] 5.0 4.8 4.7 4.8 4.4 4.0 4.6 2.6 5.0 4.0 2.7 4.0 3.9 4.6 4.9 4.3 3.9 4.9
##  [667] 4.8 4.2 4.8 3.8 3.7 5.3 3.9 4.5 3.4 5.6 5.9 2.8 6.0 4.7 3.1 4.3 4.8 3.4
##  [685] 3.9 3.1 5.5 4.6 4.9 5.6 3.9 4.0 5.5 4.2 5.0 4.9 7.1 4.8 4.5 1.8 4.7 6.3
##  [703] 5.7 2.7 3.6 4.4 4.1 4.8 3.6 2.5 4.4 4.5 4.0 3.9 2.4 5.7 4.1 4.5 4.3 5.3
##  [721] 3.9 5.0 6.2 3.3 3.9 4.3 6.1 4.2 4.0 3.8 4.1 3.7 4.4 5.9 6.1 5.2 3.6 4.3
##  [739] 4.7 5.6 4.9 4.2 3.9 4.0 5.1 4.3 3.9 4.3 4.4 3.0 4.2 3.6 4.9 4.6 4.1 4.0
##  [757] 3.5 3.6 5.6 4.6 3.0 5.2 4.0 5.4 4.3 3.4 3.9 5.0 4.4 5.0 4.2 4.7 3.6 6.5
##  [775] 3.6 3.3 4.4 4.3 6.6 5.0 4.7 3.5 4.0 3.9 4.1 4.5 5.8 3.4 3.3 2.8 3.6 3.9
##  [793] 3.9 4.0 4.8 4.6 3.5 4.5 4.8 4.2 3.4 5.0 5.1 4.9 4.9 3.8 5.0 2.5 3.4 4.1
##  [811] 3.4 3.3 3.7 5.4 3.9 3.1 5.6 2.7 5.9 3.5 4.6 5.0 3.1 5.7 4.2 4.5 4.6 3.7
##  [829] 5.3 5.9 4.8 6.0 4.5 2.8 3.1 5.2 5.2 6.1 4.2 5.1 4.5 4.6 3.9 2.5 7.1 4.1
##  [847] 3.6 4.3 5.5 5.2 4.4 6.0 3.3 4.3 4.0 4.5 4.4 3.8 4.1 5.7 3.8 4.0 2.6 4.8
##  [865] 5.4 3.8 4.2 7.2 4.2 4.0 4.3 5.4 5.1 4.0 3.9 5.1 5.3 4.6 5.0 3.5 4.4 4.9
##  [883] 4.3 3.6 2.7 5.6 4.2 4.7 3.9 5.5 3.6 3.5 5.7 5.5 4.3 6.4 5.6 3.7 4.9 6.2
##  [901] 3.9 3.6 5.3 4.3 6.9 5.4 5.7 5.8 2.8 3.9 3.5 6.1 5.1 4.1 4.6 5.4 6.5 3.5
##  [919] 4.0 4.6 4.1 3.7 2.8 5.2 3.2 6.8 3.6 4.9 4.5 4.8 2.9 4.1 4.0 6.0 3.1 4.5
##  [937] 3.8 4.7 5.0 5.1 5.1 6.0 3.6 5.4 4.3 4.7 4.2 7.1 3.6 5.0 5.3 3.6 4.1 5.3
##  [955] 4.0 3.4 3.6 4.4 4.3 5.0 6.1 5.0 4.3 4.0 4.1 4.4 3.0 4.9 5.3 4.5 4.6 3.5
##  [973] 6.5 6.9 3.2 4.7 4.5 3.8 6.6 5.2 3.5 3.7 4.1 4.9 2.8 6.1 4.5 3.6 5.9 4.8
##  [991] 3.2 4.6 4.5 4.2 4.9 4.3 3.4 5.3 5.2 3.7
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
##   2.7   6.3
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
##    [1] 3.9 5.0 5.1 4.4 4.7 3.5 2.6 4.3 3.3 4.9 5.3 3.9 5.0 5.6 5.0 6.2 4.3 4.8
##   [19] 3.6 4.9 5.1 3.7 7.0 2.2 4.9 4.1 6.5 4.8 6.1 4.8 4.4 4.4 4.0 3.0 4.7 5.4
##   [37] 3.6 4.3 5.1 2.9 4.3 4.1 3.3 4.2 4.4 6.3 4.9 5.4 5.4 4.5 2.9 5.1 4.5 3.0
##   [55] 4.7 4.2 4.7 5.4 3.6 5.0 3.8 5.0 4.7 4.7 3.6 2.9 6.2 5.2 4.2 4.3 4.5 4.6
##   [73] 4.2 5.7 3.8 4.5 4.9 4.1 6.8 3.4 5.0 4.3 4.9 3.8 4.7 2.9 4.8 3.4 4.8 4.8
##   [91] 5.5 3.0 4.7 3.9 3.7 5.3 4.4 4.8 4.1 3.0 5.5 3.8 5.2 6.2 4.1 5.6 4.5 4.4
##  [109] 6.5 5.5 5.7 4.9 3.3 5.3 4.6 3.8 4.6 4.2 2.9 3.6 5.6 4.0 6.0 5.8 5.9 4.8
##  [127] 4.2 4.1 5.4 3.9 2.9 4.2 3.7 4.8 3.8 4.8 5.3 5.5 3.6 3.4 4.3 5.0 4.5 4.8
##  [145] 4.7 5.0 5.7 5.6 6.1 4.1 5.5 4.1 5.3 4.6 5.5 3.4 4.5 4.3 4.7 4.9 3.7 4.6
##  [163] 5.3 3.2 5.2 4.1 5.3 6.4 3.3 4.0 2.5 5.4 5.1 3.0 2.0 3.8 4.9 3.6 4.8 4.5
##  [181] 3.5 4.5 5.2 3.9 4.7 4.5 5.1 4.8 4.1 3.6 4.3 5.4 4.0 4.5 4.3 4.9 4.4 3.5
##  [199] 6.5 4.9 5.3 5.3 3.0 3.7 5.4 3.5 5.8 4.3 3.5 3.9 3.8 3.8 5.0 5.5 4.6 4.9
##  [217] 3.5 5.2 4.7 4.0 5.0 4.9 4.6 3.0 3.1 4.2 4.2 3.7 5.2 3.9 6.0 3.1 5.4 3.9
##  [235] 5.8 5.5 5.6 3.9 5.3 4.9 3.2 7.5 6.3 4.7 3.8 3.7 4.6 5.2 4.3 3.9 5.9 4.0
##  [253] 5.8 5.1 5.2 3.6 4.7 5.0 4.8 5.0 4.7 3.6 4.9 3.4 4.7 4.6 4.5 5.1 4.6 4.0
##  [271] 4.5 4.3 5.8 4.3 2.5 4.5 4.1 5.3 6.4 6.1 4.2 5.6 4.2 3.7 5.6 6.3 4.5 4.6
##  [289] 4.1 4.1 3.8 5.0 5.1 2.4 4.7 5.0 4.5 3.4 5.4 5.9 4.5 3.5 5.9 4.0 4.8 4.0
##  [307] 4.1 4.4 5.1 3.2 4.7 4.0 3.8 4.8 3.8 2.0 5.0 4.5 5.4 4.6 4.0 4.2 6.2 4.8
##  [325] 3.7 6.1 5.8 4.1 3.4 4.9 3.7 5.9 5.6 4.9 4.3 4.3 4.5 6.5 4.6 4.2 3.6 4.5
##  [343] 3.1 5.8 4.2 5.9 5.1 4.2 3.3 5.9 3.8 3.4 3.9 5.4 5.4 4.1 5.9 5.6 6.4 3.8
##  [361] 5.1 4.4 5.0 2.5 5.5 4.2 4.9 5.4 3.6 3.9 5.2 4.8 6.6 4.3 4.4 3.6 4.7 3.6
##  [379] 3.2 4.0 4.4 4.7 3.5 5.4 4.2 3.9 3.1 3.2 4.7 5.1 4.8 5.4 4.3 3.6 4.5 5.5
##  [397] 3.9 3.3 3.3 5.2 3.4 3.8 4.1 3.1 5.1 4.7 4.7 4.7 4.0 3.5 4.8 3.2 3.3 2.6
##  [415] 4.4 5.5 4.7 5.1 4.6 3.3 4.5 3.7 4.6 4.2 6.0 3.5 4.0 4.7 4.3 4.6 4.7 5.2
##  [433] 3.0 4.6 5.0 3.9 4.0 4.0 4.2 5.3 4.9 5.9 4.7 3.3 4.4 4.3 5.7 4.9 2.7 3.2
##  [451] 4.2 3.7 4.6 2.9 4.5 4.3 5.5 6.2 3.3 5.4 2.9 5.1 3.3 5.1 4.2 4.4 5.4 3.0
##  [469] 4.4 5.2 5.1 5.9 4.5 4.2 4.6 4.8 4.1 3.5 3.5 5.0 2.6 4.6 4.9 4.4 5.0 4.3
##  [487] 2.5 5.2 5.1 3.9 5.0 3.4 3.4 2.9 4.5 3.6 3.8 2.6 7.1 5.6 5.2 6.1 5.0 3.2
##  [505] 2.7 2.2 5.1 4.6 4.7 3.9 4.3 4.9 4.2 4.1 2.3 4.0 3.6 3.9 5.5 3.3 4.4 5.0
##  [523] 4.8 4.2 5.3 3.0 3.2 4.7 5.7 5.7 4.2 4.3 2.8 5.9 4.4 4.8 6.7 4.4 4.7 5.5
##  [541] 4.2 4.2 5.0 3.8 3.2 3.6 3.6 5.2 4.8 3.3 4.2 4.7 4.1 6.0 4.1 3.7 2.7 3.9
##  [559] 3.7 3.2 4.5 4.8 5.4 5.6 4.0 4.6 5.9 4.8 3.9 3.2 3.4 4.4 4.4 5.9 5.1 3.6
##  [577] 6.0 5.0 6.2 5.6 5.7 4.5 3.9 4.6 4.4 4.7 3.5 5.7 4.5 5.8 3.8 3.7 3.3 4.5
##  [595] 3.4 4.5 4.6 4.7 3.9 5.9 4.4 4.9 5.3 5.4 4.3 3.3 5.0 4.9 4.5 4.7 4.5 5.7
##  [613] 4.9 5.3 5.1 5.3 4.7 2.7 4.4 3.0 3.6 5.2 4.1 5.2 3.5 5.3 3.7 5.8 6.1 3.9
##  [631] 3.5 6.1 2.8 5.0 4.3 6.2 2.9 5.2 5.0 4.9 4.4 4.8 4.2 1.9 4.7 3.5 4.5 3.7
##  [649] 4.8 4.4 3.6 3.9 6.3 3.5 4.0 3.2 3.5 4.5 4.9 2.2 2.8 6.1 4.8 4.4 3.5 4.7
##  [667] 4.5 3.0 5.7 5.4 4.8 4.4 4.8 5.7 4.2 4.2 4.8 4.0 4.7 5.1 3.8 4.7 2.9 3.7
##  [685] 3.7 4.8 5.1 4.6 4.3 4.4 6.0 4.7 4.7 4.3 3.7 6.3 4.9 6.1 4.3 5.7 4.8 5.2
##  [703] 5.5 5.0 3.2 3.9 5.3 4.7 3.9 3.6 3.7 4.7 2.4 5.0 3.7 4.0 3.1 4.2 3.6 4.6
##  [721] 6.0 3.5 5.1 6.3 4.8 4.8 4.9 5.3 4.1 3.7 5.6 4.0 3.2 4.9 3.9 4.8 3.7 5.6
##  [739] 3.9 4.7 5.6 4.9 5.2 5.3 5.2 5.6 4.8 5.8 4.0 4.9 5.0 3.0 3.7 3.3 2.2 5.1
##  [757] 4.6 3.1 4.4 4.3 5.9 4.4 3.7 6.7 4.0 3.9 4.1 4.2 5.7 4.5 4.7 5.4 4.1 4.9
##  [775] 3.9 5.7 4.1 2.1 3.6 7.2 4.5 5.4 4.0 3.6 4.6 4.5 5.7 3.5 5.4 4.6 3.0 4.6
##  [793] 4.4 3.3 4.2 1.4 4.2 5.3 3.7 4.8 3.8 3.7 4.5 3.8 3.8 3.8 6.1 3.9 5.2 6.7
##  [811] 3.8 3.7 6.2 3.4 5.2 4.1 4.1 2.6 5.5 5.7 4.7 4.2 5.7 4.4 4.9 3.6 4.6 4.6
##  [829] 4.4 5.0 5.4 6.6 5.6 3.6 4.0 3.7 4.5 5.7 4.4 5.6 3.4 4.9 5.4 5.9 2.8 3.1
##  [847] 3.8 4.0 2.4 3.9 3.8 5.3 5.0 4.2 4.5 5.2 3.8 4.5 5.7 5.9 3.9 4.2 5.4 4.8
##  [865] 5.0 4.6 4.1 4.7 4.4 3.7 5.3 4.2 4.9 2.2 4.4 5.1 3.2 5.4 5.4 3.5 5.2 5.9
##  [883] 4.5 4.5 4.4 4.9 6.2 4.8 4.2 2.8 5.2 4.2 3.7 4.8 3.3 4.0 4.5 4.0 3.9 5.1
##  [901] 5.0 4.3 5.3 5.6 2.6 4.1 4.6 2.8 2.8 3.5 4.2 5.4 4.7 5.4 4.3 4.0 3.7 2.7
##  [919] 3.1 5.2 3.7 5.9 4.0 4.5 5.8 5.7 4.0 4.8 3.3 2.8 4.7 4.7 5.0 5.5 4.0 4.1
##  [937] 5.6 4.2 5.0 3.1 3.7 4.8 4.6 3.8 4.0 3.7 4.5 3.5 4.5 6.3 3.6 4.0 3.2 4.4
##  [955] 4.8 4.8 4.9 5.4 5.5 5.3 4.3 3.9 3.8 4.6 4.7 4.6 5.7 6.3 3.7 5.1 6.0 4.8
##  [973] 4.5 4.0 4.4 4.6 4.1 4.9 5.1 4.2 4.8 4.1 4.0 4.7 4.8 4.2 2.2 4.8 5.4 4.9
##  [991] 3.1 3.7 5.5 3.3 4.1 4.9 3.9 5.1 4.1 5.2
## 
## $func.thetastar
## [1] -0.0208
## 
## $jack.boot.val
##  [1]  0.55193548  0.35702479  0.26666667  0.16979472  0.09785933 -0.07533156
##  [7] -0.22853107 -0.32668539 -0.43090909 -0.49378531
## 
## $jack.boot.se
## [1] 1.010468
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
##    [1] 5.3 6.8 4.5 5.7 3.9 2.3 5.0 4.4 3.3 5.6 3.9 5.5 3.7 5.2 6.8 5.3 5.4 4.6
##   [19] 4.6 2.7 5.5 3.0 5.2 5.6 4.3 4.3 4.7 5.0 4.7 4.5 4.6 4.2 6.0 5.5 3.2 2.8
##   [37] 4.5 4.9 6.0 5.4 3.9 4.3 5.6 4.6 5.4 5.6 5.8 4.3 3.6 4.5 4.7 4.1 3.4 3.5
##   [55] 5.6 4.9 5.3 4.7 3.5 4.1 3.7 5.3 5.4 4.0 5.3 4.5 5.0 4.9 2.7 4.8 5.7 5.4
##   [73] 5.9 4.2 4.7 5.8 3.9 4.1 4.7 3.9 5.8 5.5 4.7 5.1 4.8 4.7 5.0 3.3 4.7 4.6
##   [91] 3.8 5.5 4.4 6.0 3.4 5.5 2.8 3.8 4.5 4.8 3.6 5.5 3.3 6.9 4.2 5.6 4.6 4.0
##  [109] 4.6 4.1 6.0 3.3 6.6 5.0 5.0 4.4 4.8 4.8 4.6 3.9 3.8 3.8 4.1 2.9 5.1 4.7
##  [127] 3.0 4.3 4.0 4.5 5.0 3.0 4.6 3.0 3.2 6.1 4.2 2.9 2.9 4.0 3.9 4.7 3.9 3.8
##  [145] 3.8 4.8 4.5 3.5 3.0 3.9 3.5 4.7 4.5 4.6 4.8 5.5 3.9 6.5 3.5 4.3 4.6 3.9
##  [163] 5.1 3.4 2.5 3.7 4.5 4.1 3.8 5.3 5.4 4.3 3.1 4.4 4.1 4.1 3.8 4.7 4.1 4.0
##  [181] 5.3 3.2 4.7 4.4 5.8 3.6 4.3 3.4 5.9 3.9 4.4 5.3 4.7 5.8 4.7 4.8 5.0 5.7
##  [199] 6.1 5.6 4.8 4.1 3.7 3.5 4.2 5.1 4.0 5.8 4.1 3.9 3.8 4.1 4.3 5.6 4.8 5.5
##  [217] 4.1 4.0 4.6 4.0 5.4 5.5 5.8 3.8 3.4 3.5 4.1 5.2 5.0 4.6 4.0 3.1 5.8 3.6
##  [235] 5.0 4.0 3.6 3.9 4.6 4.8 4.6 2.6 2.9 5.8 4.4 5.3 3.9 2.7 4.0 4.3 5.0 5.1
##  [253] 5.1 5.0 3.7 3.5 3.3 3.4 3.5 5.7 6.2 4.2 2.5 4.1 4.2 6.0 5.5 4.8 4.8 5.4
##  [271] 4.1 4.8 4.6 4.8 4.5 6.3 4.3 5.0 4.2 3.4 4.1 3.5 5.3 4.0 4.6 2.7 5.4 4.9
##  [289] 6.0 4.9 5.1 4.3 5.1 4.7 4.2 2.8 5.5 5.2 3.0 2.9 4.5 3.6 4.4 5.0 3.4 1.8
##  [307] 5.5 4.8 5.2 3.6 6.1 4.5 4.4 4.4 5.6 3.9 4.1 5.7 5.4 4.8 3.1 4.3 4.4 3.4
##  [325] 2.2 4.5 4.9 3.8 4.9 4.5 4.3 4.3 4.7 4.3 4.4 4.8 4.5 5.4 4.8 4.1 6.4 2.3
##  [343] 4.7 6.0 5.9 4.2 4.4 4.1 4.3 4.5 4.7 2.5 3.0 4.0 3.8 6.7 4.8 4.6 6.1 4.8
##  [361] 4.5 4.6 3.8 3.9 4.7 4.1 5.5 4.3 6.2 4.3 4.5 5.7 4.3 4.0 2.6 4.2 4.5 3.6
##  [379] 5.3 3.7 4.4 3.4 4.4 5.4 4.4 2.3 4.8 4.9 4.8 6.2 4.9 4.4 4.5 3.8 4.0 4.9
##  [397] 4.5 3.5 3.1 5.0 4.8 4.9 4.4 4.7 4.5 4.2 5.2 5.6 3.7 4.6 5.0 3.8 4.4 4.9
##  [415] 4.4 4.5 4.0 4.1 3.7 4.6 3.7 4.9 4.4 4.2 5.4 4.7 2.9 3.6 5.1 6.1 4.4 4.8
##  [433] 4.2 4.0 4.1 4.5 5.6 5.3 2.7 3.1 3.8 6.2 3.3 4.3 4.2 5.2 5.0 5.2 4.4 3.7
##  [451] 4.9 4.2 5.0 5.0 5.0 5.5 3.6 4.9 3.4 6.1 2.9 4.2 5.0 6.3 2.4 5.1 4.7 5.5
##  [469] 4.6 4.6 5.9 5.4 4.6 5.1 2.6 4.4 5.2 5.2 3.3 3.9 5.7 4.3 5.4 5.5 5.2 5.6
##  [487] 4.9 4.6 5.6 3.1 5.0 3.0 4.2 4.8 4.2 4.1 3.7 5.4 6.9 4.0 4.8 4.3 5.4 3.5
##  [505] 5.4 2.8 3.3 2.6 4.5 4.4 3.1 4.6 5.9 4.8 4.6 4.1 5.4 4.3 3.3 4.5 3.7 5.1
##  [523] 4.4 4.9 3.6 3.4 4.2 4.0 4.2 5.4 6.0 4.9 5.5 3.9 4.3 5.9 4.3 5.1 4.8 4.6
##  [541] 4.2 5.8 6.0 5.9 3.7 3.2 4.7 4.8 5.0 5.6 4.5 3.4 6.7 3.2 4.8 5.6 4.4 4.3
##  [559] 4.7 4.7 5.2 4.6 2.8 4.9 5.2 5.4 3.5 4.4 2.9 4.6 4.8 4.5 3.8 4.6 4.4 5.6
##  [577] 4.9 3.6 4.7 3.8 3.9 3.6 4.2 4.1 3.8 5.5 5.7 4.9 3.6 5.2 4.1 5.3 4.5 4.2
##  [595] 4.3 4.7 3.3 5.6 5.6 5.2 4.5 6.0 4.7 4.6 3.2 2.7 4.5 3.6 4.4 4.3 5.8 6.0
##  [613] 3.7 3.1 6.9 4.9 2.2 4.5 5.7 3.8 5.0 5.8 4.5 6.3 4.4 4.9 5.9 4.3 4.5 3.9
##  [631] 5.2 4.0 3.6 4.8 2.5 4.8 5.4 4.5 4.7 5.6 4.2 3.8 3.9 4.8 5.0 3.8 4.4 7.1
##  [649] 4.4 3.4 4.8 3.4 3.1 5.7 3.3 4.0 4.9 3.0 5.3 4.1 4.3 5.1 5.3 2.4 3.5 4.5
##  [667] 5.8 3.9 4.0 3.8 2.9 5.7 5.8 5.2 4.7 4.8 3.9 4.5 4.9 3.6 3.3 5.5 5.2 5.5
##  [685] 3.4 4.3 4.1 4.9 3.5 4.8 4.6 5.3 5.8 4.5 6.4 4.7 5.8 4.9 3.9 4.2 4.7 5.0
##  [703] 4.2 4.7 5.5 4.0 4.6 4.5 4.7 4.3 5.5 4.7 5.8 3.8 5.6 5.4 5.3 5.0 4.0 4.3
##  [721] 4.2 4.2 5.1 5.0 3.8 5.6 4.3 3.3 4.0 5.1 4.1 4.1 4.2 4.2 5.2 6.0 3.8 4.9
##  [739] 5.6 3.3 4.8 4.0 6.0 3.5 4.6 4.3 2.7 3.8 4.0 4.1 4.3 6.2 4.1 5.3 2.8 3.9
##  [757] 3.7 5.4 4.2 4.4 5.0 5.7 4.2 5.4 5.6 3.9 5.8 5.5 4.7 3.8 2.6 4.5 5.8 3.3
##  [775] 4.4 5.7 3.9 4.4 3.6 4.2 5.0 4.8 2.7 4.8 2.9 4.9 3.7 5.9 5.2 4.2 4.6 3.9
##  [793] 4.6 4.2 4.8 3.8 4.7 3.6 4.0 4.8 6.0 4.3 2.9 4.1 3.1 2.6 4.6 3.8 4.7 5.4
##  [811] 3.4 3.1 4.7 6.0 3.9 4.0 4.3 4.0 4.3 5.8 5.0 5.0 3.8 4.5 4.5 5.0 4.0 4.8
##  [829] 3.7 3.7 4.2 5.6 5.6 3.7 5.0 4.1 5.2 4.8 5.2 3.8 3.5 5.6 5.7 3.7 4.8 5.5
##  [847] 2.8 5.1 4.4 3.3 4.7 5.1 3.4 2.9 4.9 4.0 4.7 4.6 3.8 4.4 3.3 4.8 5.1 5.2
##  [865] 5.6 3.1 4.2 5.1 2.7 5.7 3.8 4.9 3.4 4.5 3.8 5.0 5.2 3.4 3.4 5.9 3.9 3.6
##  [883] 5.0 5.8 5.1 4.6 4.4 4.7 4.2 2.3 4.1 5.1 3.8 5.5 5.3 3.8 4.2 3.8 5.7 3.5
##  [901] 4.2 4.9 4.2 3.8 5.2 4.7 5.2 4.5 3.3 5.9 3.6 5.8 7.2 4.9 5.0 4.5 4.9 4.5
##  [919] 6.3 4.6 5.7 4.0 4.0 2.8 3.7 5.3 4.6 5.2 6.0 5.4 5.6 5.1 5.3 5.3 5.2 3.5
##  [937] 3.8 2.7 4.4 4.6 6.3 3.3 4.3 3.5 4.5 4.9 5.6 4.1 5.7 4.7 5.3 3.3 4.3 4.0
##  [955] 2.6 3.5 5.4 5.2 4.8 4.3 3.7 3.8 3.5 3.6 4.1 4.4 3.5 4.4 4.5 5.5 5.9 4.1
##  [973] 5.8 5.3 3.5 4.1 4.3 6.1 5.1 3.5 4.9 3.9 4.3 4.3 4.4 3.3 4.5 4.8 4.0 5.5
##  [991] 4.8 4.6 2.4 3.4 3.9 5.8 4.4 4.4 3.1 3.5
## 
## $func.thetastar
## 72% 
##   5 
## 
## $jack.boot.val
##  [1] 5.500 5.400 5.300 5.200 5.000 4.968 5.000 4.700 4.600 4.500
## 
## $jack.boot.se
## [1] 0.9666341
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
## [1] 0.3988247
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
##   1.7683698   2.8925242 
##  (0.7285095) (1.3758985)
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
## [1]  1.9580027  0.8722472 -0.2346603  1.4930030  1.0735421  0.6120649
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
##    [1]  0.2894301934  0.7207663958  0.4253417666  1.1944464978  0.0334235694
##    [6]  0.3212269674  0.6571306710  0.8125813181  0.8198712990  0.5166706751
##   [11]  0.2837056195  0.6862778634  0.7173767039  0.0203800837  0.6706228116
##   [16]  0.0656868412  0.5238103000 -0.6876782444  0.5225970936  0.3290656608
##   [21]  0.2480735924  0.0159692628  1.2754997170  1.2884570361  0.0941896081
##   [26] -0.2786800820  0.2737877618  0.0394523458 -0.5935239130 -0.5364352230
##   [31]  0.9753379709  0.6736037181 -1.7489922129 -1.6668861560 -0.7313139811
##   [36]  0.1063068732  0.4075882147  0.8927077605  0.2291983907  1.4194408914
##   [41]  1.0461952868  0.7180888058  0.3769012813  0.2310880314  0.4671674544
##   [46]  0.8803209855  0.2983718768 -0.7156736527  0.4198563835  0.2630265488
##   [51]  1.0019756881 -0.2200007491  0.8254552086  0.4691203926  0.9086847463
##   [56]  0.6038290543  0.4495320940  1.2410404569  0.2626617343  0.4188428774
##   [61]  0.3725619952 -0.0085880666  0.0517488063  0.5580018174 -0.0833260729
##   [66]  0.3525356363  0.3874352741  0.8237230792 -1.0341233675 -0.6077294740
##   [71]  1.0262354352  0.5379338318 -0.7161269367  0.3853741023  0.1219119667
##   [76]  0.3780627312  0.6051396852  0.5252453006 -0.1371307138 -0.7005823642
##   [81]  1.7780333902  0.8776312641  0.9146959099  0.5386092001  0.7260425406
##   [86] -0.0216693240  0.1095766922  2.0150782210 -0.1366457382  0.5481828407
##   [91] -0.1573139922 -0.3578163492  0.8603585577 -0.1310405271  0.3466118681
##   [96]  0.7788542040  0.7186009036  0.8030526357  1.0187006438  0.0747673158
##  [101]  0.3831841851  0.3908909702  0.7038193498  0.6937978392  0.4914287932
##  [106]  0.1296841627  0.6048466977  0.5534802709 -0.3526790881  0.6215796227
##  [111]  0.1393891508  1.2639165779  0.0735312611 -0.0294932809  0.0185574363
##  [116]  0.5743493696  0.3511245670  0.9509169830  0.5204249514 -0.2759073206
##  [121]  0.5477989676  1.2481822258  1.8311635597 -0.4148288846  0.5875766030
##  [126]  0.4513553943  1.1989044331  0.0302004507  0.2083748947  1.0124194811
##  [131]  1.0407281414 -0.1305326136  1.0404039630 -0.1609257643  0.8497342736
##  [136]  0.6491438486  0.9865437959  1.5898669119  1.4771135517  0.6737653085
##  [141]  0.9241776124  0.7066351425  0.7415824649  0.5913277043  0.5391515904
##  [146]  0.4201182101  0.6537167162 -0.0918534582  0.1252014350 -0.0970937565
##  [151]  0.5384605067  0.4944947438  0.7867753070  0.1271121190  0.5889003376
##  [156]  0.1250837253  0.6381638137  0.2466807473  0.0523117988  0.6855680242
##  [161]  1.1784775227  0.2588660127  1.6650647892  0.8471729410 -0.3104097448
##  [166] -0.1293841224  0.4018996044 -0.0451679952  0.3458905916  1.2549811897
##  [171] -0.0027426001  0.8566060996  0.8739067586  0.1240535448  0.7816931601
##  [176] -0.3573768702  0.2488621323  0.6122030615  0.0608494751  0.0262978052
##  [181]  0.4996230185  1.7429073664 -0.1170413524 -0.5411478093  0.2406403828
##  [186]  0.8908728387  0.8556594062  0.3644050846  0.0510160165 -0.0294931867
##  [191] -0.1273454072  0.4125676024 -0.3233544675  1.6306922525 -1.2248410000
##  [196]  0.7762361825  0.1666345205  0.4016439125  0.6358565479  1.0716620086
##  [201]  0.7787224916  0.3252571978 -0.0862294142  0.0243832116  0.2113500715
##  [206]  1.0281760466  0.1131236986  0.0026099114  1.0186585666  1.0488206768
##  [211]  0.9415298016  1.3073775813  1.0344950834 -0.0873182274  0.1197218500
##  [216]  1.2522458086  0.3840245866 -0.0179540081  0.2979469182 -0.4465471708
##  [221]  0.3338796868  0.3587401815 -0.2414794538  1.2704578961  0.3832912390
##  [226]  0.4095684777 -0.2428386875  0.0734610132  1.3491899052  0.0001985261
##  [231] -0.0358217587  1.0482626887  1.0988035575  0.8754773326  0.2410322420
##  [236]  0.0223282670  0.2759059728 -0.7415974807  0.1188318554  0.0523739111
##  [241]  0.7296430184  0.2468455347  0.3063908702  0.6751254674  0.3424002016
##  [246]  0.8613356962  0.6594223059  0.1025570386  0.6837462423 -0.0148471821
##  [251]  0.3694801932  1.1908185819 -0.4146686242  1.2818564670  0.6720405620
##  [256]  0.4792490472  0.5316942617  1.0333040597  0.9658965343  0.2301803855
##  [261] -0.3157459930  0.3338796868 -0.0240949129  1.4049766566  0.2842202835
##  [266]  0.6822017573  0.5384707102  0.9624998500  0.1693717093  0.9914878426
##  [271]  0.9012750556  1.4144164707  1.3653336862  0.2702583920  0.5009224875
##  [276] -0.3440986712  0.1260136459  0.0202174008  1.4376658627  1.0033671071
##  [281]  1.1346262803  0.5174221963  0.9583962236  0.5825180230  0.2384426590
##  [286]  0.8927806768  0.3585984082  1.0345852563  0.2743792974  1.0039335649
##  [291]  0.3842434487  0.0238388513  0.8626501424  0.9699317402 -0.7073472273
##  [296]  0.5925529475 -1.2272969215  0.1830820142  0.5357731458  0.4686979058
##  [301]  1.2644452867  0.8693179717  0.2264078181 -0.3949250760  0.6048254408
##  [306]  0.1139488369 -0.1438215831  0.5520860479  0.7874853562  0.2899551278
##  [311] -0.4079848460  0.0651233276  0.5422769738  0.4225989130  1.4695257580
##  [316] -0.2548520105  0.6445741189  0.3988762267  0.1020342620  1.1189913824
##  [321]  1.4634301534 -0.4593819179 -0.2650822714  0.7774256602  0.5001569301
##  [326]  0.3606101018  0.7402224509  0.5685301801 -0.3130974328  0.1992954414
##  [331]  1.2724576408  1.4018676130 -0.1632456357  0.3740718976 -0.5666175035
##  [336]  0.3090997913  1.0684564376  0.7143101657 -0.0668358989 -0.0158363645
##  [341]  1.0019242571  0.5105870943 -0.0285506491  1.3124795609 -0.0627924292
##  [346]  0.8576846472  0.6346968998  0.8546068552 -0.0194205094  0.9028353823
##  [351]  0.8414675211  0.2060168832  0.2627119659 -0.0182537122 -1.1305835367
##  [356]  0.5610930456  1.3855096519  0.5557055394  1.8352757459 -0.0348746143
##  [361]  0.7159822415  0.5683231944  0.3266625915  0.5325268739  0.1192646314
##  [366]  1.1884361576  0.0001208478  0.5951319106  0.8478102699 -0.5695799517
##  [371]  0.4831086976  0.2674666258  1.4458972280 -0.0214662258 -0.3281988746
##  [376]  1.8149240522 -0.7354006117  0.3490431566  1.1716824282  0.5524227524
##  [381] -0.0311925166  0.6956803088 -0.0511520939  0.3927024224  0.2482671969
##  [386]  0.5074222370  0.3541123285  0.6779668508 -0.8761894588  0.4832044625
##  [391]  0.7054515134  0.0172422187 -0.2771914613  0.1272044274  0.1487309657
##  [396] -0.0007474468  0.2781945979  0.0022638711  1.3758595593  0.3861789411
##  [401]  0.7848184621  0.6227198036  0.3984072207  0.2827686307  0.2838940735
##  [406]  0.0734610132  0.2607542477 -1.0679371907  0.0050811777  0.1012005526
##  [411]  0.3642718202  0.6798206200  0.5285119138  0.7916020398  0.6758582347
##  [416] -0.2616807213  0.4536297213  0.4288427234  0.7968720676  1.2344965045
##  [421]  0.4021764014  0.1458735155  1.6131047499  0.8434635441  0.6612280519
##  [426]  0.6689341019  0.3879188204  0.3990967856  0.4094247703  0.2667869908
##  [431]  0.4331406277  0.5348358007 -0.0038082278  0.3108244059  1.1651685043
##  [436] -0.7004001680  0.3433905425 -0.4161680528 -0.2625623866  0.4113554832
##  [441]  0.5921697006  0.9179902672  0.4687080716  0.0286439404 -0.1492148765
##  [446] -0.6562358798  0.3110335321  0.2485423040  0.1138825184  0.2764617813
##  [451]  0.1371946926 -0.0000796112  0.0687814297  0.0984038066  0.8665225851
##  [456]  0.4159388484 -0.4067392162  0.5184997557  1.1083077143 -0.1304990036
##  [461] -0.5903802016  0.0273681803  0.0777639190  0.1625519234  0.8486315788
##  [466] -0.3431316607  0.3408789387  0.5470550662  0.6531692480  0.6836842607
##  [471]  1.0620760939  0.3992003938  0.6826448218  0.0179857252  0.8908728387
##  [476] -1.1833627008  0.7302415399  0.8582597049 -0.2910162627  0.3087619323
##  [481]  0.0279459015 -0.1162241832 -0.0625384367 -0.5954628717  0.8307719130
##  [486]  1.0096621516  1.0763081206  0.8124510838  0.0013682755 -0.2218712830
##  [491]  0.2737498940  0.7170346027 -0.0218176216  0.4106379887  0.0112472928
##  [496]  1.2228261603  0.4792490472  1.5866864157  0.3385136305  0.2291307542
##  [501]  0.8414675211  0.4116386767  0.0168213632  0.5971988213  0.5953649362
##  [506]  0.5600764446  1.0015476058 -0.1739331345 -0.3192302729  0.5559720791
##  [511]  0.3263611138  1.3578812654 -0.0420241023  0.4505718465  0.6080380214
##  [516]  0.1323093280  0.8536483077  0.0811649170  0.4104561295 -0.5825219174
##  [521]  1.0994901902  0.0628283080  0.7858820451  0.4285475721 -0.3160626655
##  [526]  0.4964289664 -0.0372989296 -0.1622311694 -0.5693771479  0.9947252695
##  [531]  0.2995201129  0.5777782258  0.5083719656  0.1397762972  0.6847829102
##  [536]  1.2822068681  0.5059692155 -0.3896159662  0.2836170963  0.9300253649
##  [541]  0.6624942826  0.8555820144  0.8079939375  1.3062775666  0.0028218611
##  [546]  0.9542736826  0.3866931081  0.6812616753  0.5238377413  0.6787503792
##  [551]  0.8778891117 -0.1714572838  2.0806849028  0.9387170804  0.3774593529
##  [556] -0.2706933492  0.5318869347  0.1713825353 -0.0161426827 -0.2661148527
##  [561]  1.0874770108  0.6661791965  0.5989910567  0.9769106397  0.2291448664
##  [566]  0.6161730668  0.3989687607  1.3014063820  0.2758874096  0.7166657928
##  [571]  0.5380146349  1.6924262948  0.3981985833  0.1312861141  0.6866590222
##  [576]  0.2745106040 -0.4083814931 -0.0508261290  0.5486407054 -0.0248222707
##  [581]  0.9798225150 -0.2091672023 -0.0148471821  0.1741095876  0.1703907900
##  [586]  1.0624483482  0.1932089054  1.9631258744 -0.5486131874  1.5711472459
##  [591]  0.5538716640  0.7222727545  1.1788650254 -0.2771527293  0.2525089970
##  [596]  0.7437384307  0.7845414953  0.2641041972  0.2549704202  0.4094247703
##  [601]  0.4986437934  0.0541689374  0.3606537897  0.5618655285  0.2683047219
##  [606]  0.3216356969  0.4191138901  0.6503264683  1.0391111175  0.5234603256
##  [611]  1.7306935533  0.9131823727  1.2256828871  0.8396920568 -0.1492212301
##  [616] -0.1298357005 -0.2007809534 -0.2415121046 -0.2544467532  0.3638667076
##  [621]  0.9972859597  1.5988825331  0.8222119310  0.2540728114  0.0925037935
##  [626] -0.0057593487  0.5930573072  0.6946808667  1.3089899764 -0.2049593553
##  [631]  0.3073568784  0.5329846667  0.2403537466 -0.0319572153  0.1275065364
##  [636]  1.6329660406  0.8860010074  1.0033274363  0.0026649515  0.3510953731
##  [641]  1.0658340630  0.7665252744  1.1161317089 -0.1166683789  1.0463859682
##  [646]  0.5348358007 -0.0419837034  1.2836505082  0.3662831921  0.9631362895
##  [651]  0.3003033827  0.5530452735  0.7080215071  0.4468385082  0.2126997870
##  [656]  1.4067527234 -0.2026343205  0.0777426994  0.7413564442  0.7656201997
##  [661]  0.7659938487  0.7994198526  0.5359904555  0.1005562598  0.1569930203
##  [666]  2.1055075531  0.8806569193  0.3454226657  0.9224392623  0.2205497401
##  [671]  0.8450899516  0.8433663999  0.1436078333  0.7173767039  0.2405141326
##  [676]  0.1216773097  0.7503737316  0.2392272414  0.2252206190 -0.3535805458
##  [681]  0.2458826734 -0.3785203315  0.2320112476  1.0760628908  0.8613760348
##  [686] -0.0079617932 -0.2845496370  0.9091337880  0.0943072839 -1.5136849722
##  [691]  0.7621612050  0.6815259289  0.9844129600  0.2826630686 -0.6237605947
##  [696]  0.8697409310  0.8220617720  0.2472129347  0.8031786700  0.4886604804
##  [701]  0.6947316435 -0.0145790883  0.5228938616  0.2609261140 -1.8622866449
##  [706]  0.0524610932  0.2737131503  0.8066148413  0.5009441123  0.4488432951
##  [711]  0.0040644418  0.2763000000  0.3730569109  0.6523998682  0.9926313157
##  [716]  0.4692436869  0.7905964827  1.1739530201  1.0371005526 -0.2907101137
##  [721]  0.3627997649 -0.0070039836 -0.2839400518 -0.0547927419  0.7739828207
##  [726] -0.1620743379  0.3509326907 -0.0359534322 -0.1172795882  0.5128946963
##  [731] -0.0028636259  0.5246541505  0.8441114834  0.3762194907  0.7143615042
##  [736]  0.9298894901 -0.1374823670 -0.1444855148 -0.2644372470 -0.0650984403
##  [741]  0.5540933229  0.0268368470 -0.2725615350  0.9390141151  0.6858971375
##  [746]  0.4439309039 -0.0995893570 -0.1790445528  0.1353437893  0.7202272487
##  [751]  1.0815098239  1.3428300246  0.4350685028  0.0919778919  1.1009517677
##  [756]  0.0848740991  0.8241129676 -0.6919815781 -0.4072418776  1.0045571382
##  [761]  0.8887886526  0.0764273080  0.3822950143  0.3974778831  0.0085140267
##  [766] -0.2825294907  0.7396141386  0.8437296070  1.0343034102  1.1037884016
##  [771]  0.7320698363  2.0446713092  0.2388317634  0.2433873729  0.6418675124
##  [776]  0.6570853701  0.7134142380  0.5099546759  1.5948244838  0.0376340706
##  [781]  0.9351550593  0.7414650082  0.7722066731  0.0748524848  0.6497636092
##  [786]  0.2219111378  0.4000688477 -0.4041986897  0.4591986085  0.1656429805
##  [791]  0.2403519200  2.0385197374  0.2638708773  2.1476735977 -0.2768119767
##  [796] -0.3625758010  0.2526307042  0.4951635065 -0.7002098183  0.6637623557
##  [801]  0.6040546868 -0.0530106069  0.6478910326  1.3253291583  0.6015185221
##  [806] -0.1611806196  0.1193814607  0.4099702073  0.2817137933  0.7739204112
##  [811]  0.9317162847 -0.7083265849 -0.0334323397  0.3081430323 -0.6367674695
##  [816]  0.5424680403  0.3948615185  0.4597586303  0.1247722754  0.7982362755
##  [821]  0.7190549736 -0.0459430566 -0.0401839937 -0.2874534028  0.8489697035
##  [826]  0.2087960887 -0.2640169089  0.5299822687  0.7411200726  0.5183268789
##  [831]  0.4445276249  0.0293810042 -0.0764919760 -0.0665479523  0.8053859881
##  [836]  0.7821435788  0.3496050253  1.0451585081  1.1760527149  0.7003372670
##  [841]  0.9825731015 -0.3117161520  0.3401568878  0.3570943561 -0.0110465604
##  [846] -0.3195630835  1.1058028002  0.6102705214  1.2229506721  0.3020130684
##  [851]  0.3576979401  0.6074555683 -0.7473504110 -0.1050624183  0.7132086768
##  [856]  0.5216267267  0.3840245866  0.9180328774  0.2438246375 -0.3779746405
##  [861]  0.1571087455  0.8464893812  1.3195646981  0.1490292182  0.3867650618
##  [866]  0.2222916202  0.6154529422  0.1683642893  0.7129012240  0.8255790141
##  [871] -0.2911449368  0.0146840731  0.3329432889  0.7405438324  0.1014032055
##  [876]  0.8427170413 -0.3731181514 -0.3583416759  0.1059610643 -0.0055214067
##  [881]  0.0272360796  0.0959590347  0.5934050334  0.1653655427  1.2712742837
##  [886] -0.3040482005  0.5850581969  0.7624304133  0.8386326828  0.0662665453
##  [891]  0.3801236153  0.0594299058  0.4626729402  0.3422628572  1.2739252701
##  [896]  0.0114107202  0.1072002307  0.4031072482  0.5494465441  0.1080024140
##  [901]  0.0025137224  0.1156775707  0.1450889854  0.5402187634  0.7188783196
##  [906]  0.6746150934  0.5025947750 -0.0065833793 -0.0007042109  0.6671341132
##  [911]  0.2472046977  0.9180578806  0.6859334987  0.3379011231  0.1630122425
##  [916]  1.4317545872  0.6833423918 -0.3828237942 -0.2675190039  0.5106704495
##  [921]  0.4212038734  0.5246781888  0.8637330071  0.1956053548 -0.4281563340
##  [926]  0.7293252609 -0.5565757745  0.2387908975  1.0363916306  0.7123801096
##  [931] -0.1539078465  0.8433663999  0.1214664705  0.5623843825  0.8434635441
##  [936]  0.1931303712 -0.0148416286 -0.2705556399  0.7934288863  0.5476631879
##  [941]  0.3505821195  1.6823431316  0.8852844742  0.6243655634  0.5812605917
##  [946] -0.7562189588  0.7376066016  0.4614373984  0.6997736373 -0.2685184617
##  [951]  0.7402287478  0.2703270064  1.0015241007  0.5363529686  0.2262664174
##  [956]  0.2641625438  0.2930438384  0.4616760713  2.4679671416 -0.2939024350
##  [961]  0.6873891246  0.5392161162  0.0030044353  0.3424002016  0.1670337218
##  [966]  0.2406668603  0.1475796577  0.2576009503  0.3400206094  0.6015234473
##  [971]  0.9898589362  1.0393189477  0.5492524403  1.0382592970 -0.2554820080
##  [976]  0.9136662520  0.4957551373 -0.1009602207  0.9835319468  0.5686297384
##  [981]  0.2729073225  0.1141056124  0.4213136072  0.3495522713  0.3948856702
##  [986]  0.3257113572  0.6696379970  0.1251316588 -0.0102915891 -0.6261702965
##  [991]  0.6619073438 -0.5375291829  0.2307682696  0.5017894798  1.3358101371
##  [996]  0.5329460214 -0.2679127610  1.2255760140 -0.1550578539  0.5966188364
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
##   0.61138025   0.42776631 
##  (0.13527158) (0.09564976)
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
## [1]  0.0003557019  0.5719738190  0.2131592587 -0.2919941454 -0.4783360867
## [6] -1.4724803608
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
## [1] -0.0229
```

```r
sd(results2$thetastar)  #the standard deviation of the theta stars is the SE of the statistic (in this case, the mean)
```

```
## [1] 0.9215177
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
## t1*      4.5 -0.001301301   0.9264799
```

One of the main advantages to the 'boot' package over the 'bootstrap' package is the nicer formatting of the output.

Going back to our original code, lets see how we could reproduce all of these numbers:


```r
table(sample(x,size=length(x),replace=T))
```

```
## 
## 0 1 5 6 9 
## 3 1 2 2 2
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
## [1] -0.0073
```

```r
se.boot
```

```
## [1] 0.8932271
```

Why do our numbers not agree exactly with those of the boot package? This is because our estimates of bias and standard error are just estimates, and they carry with them their own uncertainties. That is one of the reasons we might bother doing jackknife-after-bootstrap.

The 'boot' package has a LOT of functionality. If we have time, we will come back to some of these more complex functions later in the semester as we cover topics like regression and glm.

