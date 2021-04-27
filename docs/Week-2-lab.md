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
## 2 5 6 7 8 9 
## 2 1 1 3 1 2
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
## [1] -0.0105
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
## [1] 2.682863
```

```r
UL.boot
```

```
## [1] 6.296137
```

Method #2: Simply take the quantiles of the bootstrap statistics


```r
quantile(xmeans,c(0.025,0.975))
```

```
##  2.5% 97.5% 
##   2.8   6.2
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
##    [1] 5.7 6.1 4.4 3.7 4.6 3.0 6.5 4.9 4.9 4.6 5.2 5.3 4.4 4.4 5.8 6.0 3.7 3.5
##   [19] 5.3 4.2 5.4 4.1 4.4 4.4 5.1 6.6 6.5 4.6 4.6 4.1 4.4 4.9 5.5 3.2 3.5 4.3
##   [37] 5.1 3.3 5.0 4.2 2.8 4.0 3.0 4.8 5.2 5.4 4.2 5.5 4.4 4.9 3.8 3.8 5.6 3.1
##   [55] 5.3 5.6 3.9 3.9 4.1 4.5 5.7 3.5 4.6 4.6 3.8 4.5 5.0 4.2 3.4 3.7 4.4 4.4
##   [73] 5.3 4.7 5.3 4.8 3.8 5.1 4.8 3.6 4.1 3.4 4.7 3.4 3.4 4.4 4.7 3.9 4.7 5.3
##   [91] 2.8 5.4 5.5 3.6 5.6 5.0 4.6 3.9 5.3 3.4 5.5 4.0 5.1 5.8 5.1 5.0 3.8 3.1
##  [109] 5.9 4.6 4.5 5.3 5.9 4.0 3.2 3.9 3.1 2.4 3.6 4.0 4.8 6.1 4.9 3.9 4.7 3.8
##  [127] 4.4 5.5 5.0 3.0 3.2 4.3 4.7 4.9 4.3 4.1 5.0 3.5 3.1 4.4 5.5 3.7 5.8 5.3
##  [145] 3.7 3.2 3.8 4.8 5.4 4.1 4.9 3.4 4.6 3.5 3.9 3.3 4.3 3.6 5.0 4.3 4.3 4.3
##  [163] 5.6 3.7 2.9 4.2 4.0 5.8 3.6 2.9 5.5 4.8 4.0 3.3 5.0 4.6 4.0 5.9 2.9 6.1
##  [181] 4.6 5.6 6.2 4.4 5.0 3.8 3.6 2.9 4.1 4.5 4.8 4.2 4.5 4.2 5.4 6.1 4.1 3.4
##  [199] 4.9 5.0 2.6 3.7 5.9 4.6 2.3 3.9 4.7 5.1 4.9 4.1 4.9 6.0 4.8 4.5 4.6 5.0
##  [217] 5.6 3.1 5.9 5.9 4.3 3.5 5.8 4.5 5.5 5.3 3.9 4.9 5.5 4.7 4.9 3.4 5.4 2.9
##  [235] 3.9 5.5 2.6 3.9 7.3 5.9 6.1 5.4 5.6 5.8 4.1 4.2 3.9 5.2 4.0 4.5 4.9 4.4
##  [253] 3.8 4.2 5.9 5.0 4.1 4.2 4.7 5.1 4.1 5.2 4.9 4.3 4.4 6.3 4.6 4.6 4.5 2.3
##  [271] 3.8 6.0 4.2 5.2 4.7 5.2 4.3 5.5 4.0 5.0 3.9 4.3 3.1 4.1 4.2 3.5 3.7 3.5
##  [289] 6.0 5.3 3.7 3.6 5.2 4.5 2.9 3.9 4.4 3.9 5.1 3.5 4.0 4.9 6.4 3.3 3.5 3.5
##  [307] 4.3 5.3 4.0 3.4 3.1 3.5 4.5 3.1 5.0 4.8 4.5 5.4 3.1 2.8 5.2 3.9 3.1 3.6
##  [325] 3.4 3.7 5.3 3.5 5.0 4.8 5.5 3.0 3.6 3.1 4.9 5.7 2.5 3.4 3.9 5.1 4.0 4.6
##  [343] 4.4 5.5 5.6 3.2 5.5 5.7 3.8 4.2 5.4 4.3 4.0 4.7 3.3 4.5 4.5 2.7 4.7 4.0
##  [361] 5.1 4.1 5.5 5.1 3.5 3.8 5.4 4.7 5.8 4.8 3.7 3.8 3.1 5.1 2.5 5.9 6.2 3.9
##  [379] 4.4 5.8 6.2 5.4 4.8 4.4 5.1 5.7 4.4 5.0 4.6 5.2 3.9 4.4 4.7 4.6 4.2 2.8
##  [397] 5.3 4.6 4.5 3.7 5.3 5.2 4.8 4.2 4.6 6.0 2.7 2.5 4.4 4.5 4.3 4.5 4.2 4.1
##  [415] 3.6 3.9 5.1 4.7 5.2 5.6 5.0 3.3 4.9 5.1 5.0 6.0 5.1 3.4 4.9 3.8 3.3 3.9
##  [433] 6.0 4.2 5.6 5.4 3.8 4.4 5.1 3.9 5.3 4.7 4.3 4.1 3.5 4.7 5.9 4.3 5.4 4.5
##  [451] 4.2 6.2 5.3 5.7 4.9 4.0 4.3 4.5 4.7 3.0 4.9 5.6 3.8 2.9 4.8 5.6 4.4 4.7
##  [469] 5.9 3.5 3.1 5.5 6.3 6.3 5.5 4.2 2.6 4.6 4.2 4.9 3.4 5.6 4.4 4.2 5.2 4.0
##  [487] 2.7 4.7 4.5 4.1 4.5 5.1 3.2 4.4 4.4 2.7 5.6 4.4 5.4 4.9 3.8 5.3 3.7 3.1
##  [505] 5.2 4.9 4.4 3.0 4.3 3.3 4.4 5.0 3.0 4.0 5.3 4.3 3.9 3.3 5.6 4.9 4.4 4.1
##  [523] 2.9 5.3 4.2 4.7 5.0 4.9 7.5 4.3 4.0 4.5 4.8 3.9 5.1 5.4 5.3 3.8 4.5 4.7
##  [541] 3.2 4.4 4.2 4.7 3.8 4.4 5.6 6.0 4.8 4.4 3.3 5.4 3.5 3.2 4.2 4.9 3.7 5.3
##  [559] 4.4 2.7 3.9 5.0 5.4 3.8 4.9 2.7 4.5 5.1 6.2 3.2 4.3 3.0 5.2 5.5 4.9 5.1
##  [577] 4.5 3.4 4.4 3.6 4.5 5.5 3.0 6.5 5.2 4.8 2.8 4.6 3.3 5.8 4.9 4.1 4.1 5.7
##  [595] 3.4 5.6 4.6 4.0 4.1 3.8 4.8 4.3 3.7 4.0 4.4 3.5 6.1 4.8 4.6 3.1 4.3 5.2
##  [613] 4.2 4.2 5.2 2.7 5.2 5.9 4.9 4.7 4.4 4.1 3.1 4.1 4.3 3.2 5.8 5.7 4.1 4.5
##  [631] 3.3 5.4 3.7 4.3 4.0 3.6 5.0 3.7 3.7 4.1 4.5 5.3 4.2 4.8 2.7 5.2 4.6 3.4
##  [649] 4.3 3.8 3.6 4.4 6.3 4.6 6.0 5.8 5.5 3.7 2.9 6.1 6.5 4.6 5.7 4.2 4.4 4.0
##  [667] 5.7 5.4 4.9 3.3 4.6 4.7 4.7 4.4 4.1 5.1 4.7 3.0 4.1 2.5 4.4 4.2 6.2 4.6
##  [685] 2.4 4.8 4.6 4.4 5.1 4.8 3.9 3.8 4.6 3.3 3.6 3.8 4.5 5.4 4.5 4.5 6.6 5.1
##  [703] 1.7 4.7 2.4 2.6 3.9 3.0 5.1 3.2 4.5 6.5 4.1 5.2 5.1 5.1 4.8 5.2 4.5 4.8
##  [721] 5.6 2.2 3.6 5.6 5.9 4.1 5.0 4.5 5.0 4.1 6.4 5.9 4.1 4.3 5.0 2.9 4.4 5.2
##  [739] 2.8 4.6 4.6 5.6 4.3 3.7 4.6 3.6 3.8 3.2 4.9 3.1 5.1 3.8 3.5 3.5 4.7 4.6
##  [757] 4.5 3.5 5.3 3.8 4.9 3.3 2.4 4.5 3.3 6.1 3.6 5.5 4.1 3.2 3.7 3.6 5.0 4.0
##  [775] 5.8 6.1 4.2 5.5 4.4 3.6 5.3 5.4 5.0 3.4 5.8 4.8 5.1 4.4 5.0 3.5 4.9 3.9
##  [793] 3.9 4.0 5.5 4.1 4.3 4.3 5.2 4.4 4.3 3.8 3.1 4.3 3.1 4.3 4.9 2.7 4.7 3.7
##  [811] 3.4 5.0 4.0 5.0 5.0 3.0 4.1 4.3 3.8 4.0 5.3 5.4 4.4 4.1 5.7 4.2 3.6 4.7
##  [829] 3.3 2.0 5.3 3.9 3.6 4.9 4.0 5.0 3.2 4.6 2.9 3.8 5.3 3.1 4.6 4.3 4.6 4.7
##  [847] 3.7 4.2 5.4 5.0 4.6 5.0 4.3 4.9 4.6 4.7 4.2 6.2 5.3 4.1 3.4 5.2 4.9 3.7
##  [865] 5.0 4.1 3.9 5.3 5.4 4.9 4.4 4.6 4.6 6.3 4.1 3.8 4.2 3.5 4.2 4.6 4.9 4.5
##  [883] 4.7 5.5 5.7 4.4 4.9 3.8 3.2 4.7 3.8 2.7 4.3 5.8 4.5 3.6 3.8 4.7 5.4 3.9
##  [901] 3.8 3.1 5.0 3.0 6.2 5.4 4.2 4.2 5.2 4.5 5.1 3.2 3.3 4.5 3.9 4.2 5.9 4.5
##  [919] 5.7 3.4 3.8 2.6 5.2 4.6 5.3 4.4 4.5 5.0 4.0 4.8 3.1 4.6 3.1 5.9 3.9 4.4
##  [937] 4.4 3.4 5.0 2.9 3.8 5.5 5.5 4.3 5.3 3.8 3.5 4.5 4.2 4.5 4.0 5.0 4.2 3.9
##  [955] 5.2 3.4 4.9 5.6 5.5 4.8 3.5 6.1 5.0 5.4 4.8 4.3 4.8 4.7 4.2 4.8 5.5 5.7
##  [973] 5.0 5.4 5.4 4.1 5.3 5.5 3.7 2.7 4.7 5.2 6.0 4.2 5.0 7.0 3.1 5.3 4.9 5.0
##  [991] 6.0 4.2 3.2 4.8 5.1 6.1 5.5 6.1 4.3 3.2
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
##    [1] 4.7 3.6 4.8 4.2 5.4 5.5 4.8 3.8 3.8 4.8 5.1 4.7 4.5 4.2 4.2 5.0 3.9 4.3
##   [19] 4.7 5.1 4.6 4.9 4.2 5.0 4.2 3.6 3.5 4.0 4.0 3.3 5.0 4.1 4.4 5.2 4.7 4.6
##   [37] 3.5 4.4 5.4 4.3 3.3 4.6 4.0 5.3 5.8 4.9 3.7 5.1 3.7 2.6 4.8 4.9 5.2 5.4
##   [55] 5.9 5.3 6.4 3.7 3.4 5.4 3.6 4.9 5.1 3.8 4.1 5.6 4.5 3.4 3.7 3.6 4.9 3.5
##   [73] 4.4 3.8 3.7 4.9 4.6 5.5 4.4 4.4 4.7 3.9 4.8 4.7 4.0 3.6 4.2 3.9 4.7 4.6
##   [91] 5.1 4.8 5.1 4.1 5.1 5.6 3.3 5.5 4.3 4.3 5.8 2.4 6.2 3.5 3.6 3.2 4.3 4.9
##  [109] 3.5 4.6 3.5 4.6 4.3 5.4 6.1 5.1 4.3 5.0 4.4 3.5 4.1 4.9 5.2 4.6 5.1 3.7
##  [127] 5.6 3.5 3.5 4.8 3.8 3.7 4.1 3.8 5.3 4.7 3.6 5.1 4.9 4.5 6.0 4.5 4.6 5.3
##  [145] 5.1 4.9 4.9 2.9 5.3 4.4 5.3 4.3 3.7 5.0 3.5 4.0 4.8 5.0 4.1 4.1 5.7 7.1
##  [163] 3.3 5.6 4.3 3.6 5.5 2.8 2.8 3.0 5.2 5.5 4.5 4.4 3.5 4.9 4.0 5.2 4.4 3.5
##  [181] 5.3 5.7 4.0 4.7 4.5 4.5 4.4 3.6 5.1 5.1 4.3 3.5 3.4 4.7 5.0 3.9 3.7 4.8
##  [199] 5.5 3.1 4.7 5.2 4.7 5.1 4.5 5.0 6.0 3.3 3.0 4.3 3.7 5.5 3.4 5.1 3.6 5.4
##  [217] 4.0 4.9 5.6 4.6 4.6 4.6 4.0 4.2 4.1 4.8 4.7 3.6 6.1 3.1 4.2 5.3 4.3 3.7
##  [235] 6.0 4.0 3.6 4.4 4.0 5.7 3.8 4.7 4.1 4.7 3.7 3.2 6.1 5.5 4.2 5.0 5.0 4.4
##  [253] 4.7 3.9 2.4 4.6 4.4 3.9 5.0 5.1 2.8 5.7 3.1 4.2 4.5 5.0 5.5 3.8 4.1 5.3
##  [271] 4.2 4.6 5.9 4.4 6.4 4.7 2.9 4.8 4.5 3.9 3.2 3.9 4.5 3.9 3.9 3.5 3.7 3.7
##  [289] 5.8 4.6 3.3 5.4 4.0 6.5 4.8 3.0 5.4 2.8 6.8 4.3 5.7 4.9 3.7 5.9 2.8 4.1
##  [307] 3.5 4.3 4.4 3.2 5.4 2.9 4.4 4.7 4.0 3.3 3.5 4.8 5.3 4.4 2.9 3.9 2.6 4.0
##  [325] 2.6 4.6 5.5 4.4 3.5 5.4 4.1 5.3 5.1 6.5 4.8 3.5 5.5 5.2 5.6 6.3 5.9 5.7
##  [343] 5.0 3.1 3.7 4.3 4.7 6.1 5.3 4.2 4.4 4.5 4.0 7.1 5.9 4.4 4.7 2.7 3.9 3.5
##  [361] 3.9 4.7 4.9 2.9 4.6 4.0 5.9 3.5 3.6 2.9 4.3 3.7 4.9 4.1 3.2 4.8 4.4 3.7
##  [379] 4.4 5.3 4.9 4.9 4.4 4.4 7.0 5.4 4.6 5.2 4.7 7.3 4.0 4.2 4.1 6.6 4.1 6.8
##  [397] 2.6 5.1 4.1 3.7 4.6 5.4 4.3 4.1 4.9 5.2 3.2 5.0 6.0 3.9 4.9 4.3 3.5 4.5
##  [415] 3.4 3.9 5.2 4.5 2.6 3.9 4.6 4.5 6.1 6.1 3.4 4.3 5.4 4.0 3.9 4.8 4.9 3.9
##  [433] 4.3 5.7 3.7 4.8 3.9 5.9 5.0 5.6 4.6 3.3 4.3 4.7 3.2 4.9 6.0 5.6 5.0 3.5
##  [451] 5.5 4.2 3.8 3.4 5.1 4.5 3.6 6.5 3.6 4.8 2.3 4.2 4.6 5.4 5.2 4.5 5.0 3.7
##  [469] 4.2 5.0 5.7 3.8 4.0 5.3 4.5 3.0 4.1 5.2 4.9 6.4 3.3 4.2 3.4 4.9 5.3 5.1
##  [487] 4.6 4.4 5.1 3.2 3.4 4.2 4.0 3.8 4.6 4.5 4.1 5.1 3.4 2.4 4.7 4.6 3.8 3.4
##  [505] 5.3 3.4 3.8 5.1 5.1 5.0 3.5 3.8 5.1 2.3 5.4 5.2 4.3 5.1 4.2 5.4 5.4 4.3
##  [523] 4.9 4.5 3.6 4.6 4.7 3.6 4.8 4.8 2.7 2.7 4.1 3.8 5.3 6.1 4.4 3.9 5.4 4.2
##  [541] 4.3 5.1 4.1 4.5 3.6 5.6 6.3 5.8 4.9 5.3 4.3 4.9 3.4 4.7 4.1 3.0 4.3 4.9
##  [559] 4.0 5.7 4.5 6.3 4.3 5.1 5.7 4.0 4.6 5.2 3.4 5.0 5.7 4.2 5.7 2.7 3.5 6.0
##  [577] 2.2 5.3 5.6 3.9 5.4 2.8 4.3 4.2 4.3 3.3 4.3 4.9 3.9 3.6 4.0 3.1 4.3 4.4
##  [595] 6.4 4.5 4.7 6.0 5.8 5.1 4.5 4.6 5.2 4.1 5.9 5.1 3.3 3.8 4.0 4.2 3.9 5.6
##  [613] 4.6 5.4 4.9 4.7 4.7 5.5 4.7 2.9 5.0 5.4 4.1 3.6 3.7 4.2 4.0 3.7 3.0 4.1
##  [631] 4.9 2.6 4.2 5.3 5.3 5.8 3.5 4.6 4.7 5.4 5.6 4.5 5.1 3.8 3.6 3.9 7.1 4.5
##  [649] 4.6 6.5 3.5 5.5 5.1 5.2 4.0 4.8 5.1 4.4 4.6 3.4 5.5 4.8 4.5 3.6 4.1 5.8
##  [667] 3.4 5.2 4.8 5.5 4.1 3.6 3.6 4.1 4.5 4.7 5.0 3.7 6.1 4.5 5.0 3.6 5.4 2.8
##  [685] 4.0 3.9 3.3 4.8 4.2 4.6 3.7 5.3 5.1 3.7 3.3 4.3 3.0 4.6 4.6 3.1 5.1 4.2
##  [703] 5.1 5.3 3.0 3.7 4.5 5.6 4.9 3.8 4.8 5.5 4.7 4.8 5.0 5.4 3.4 3.9 4.9 4.8
##  [721] 5.7 3.8 4.7 4.0 5.1 4.4 3.1 4.2 6.2 5.8 5.1 4.6 3.7 3.6 4.9 4.9 3.3 2.3
##  [739] 4.8 4.7 4.4 5.2 3.6 6.1 3.9 3.6 3.1 3.7 2.3 4.1 4.5 3.9 3.0 4.1 4.1 4.1
##  [757] 5.1 5.1 5.0 5.1 4.8 5.3 3.9 4.4 4.8 2.8 4.8 5.5 5.0 3.1 5.6 5.7 4.3 4.2
##  [775] 4.5 4.7 4.0 5.7 2.3 5.7 4.4 5.7 3.1 4.6 3.4 4.3 4.7 5.4 5.9 4.1 3.0 5.5
##  [793] 2.7 5.6 5.1 6.4 3.4 2.1 3.1 2.8 4.6 3.3 4.1 4.9 5.6 5.3 6.0 5.1 3.1 5.4
##  [811] 4.7 3.6 3.0 4.8 5.0 3.1 3.9 3.6 5.3 5.2 4.6 5.7 6.0 4.4 4.0 6.3 3.8 4.4
##  [829] 3.9 5.3 3.5 3.4 6.8 5.5 4.1 4.2 3.0 6.4 4.4 3.8 5.5 3.7 6.3 6.8 3.1 3.9
##  [847] 4.1 4.0 4.7 4.9 5.1 4.9 3.9 5.3 4.2 4.7 4.5 4.5 4.3 4.1 4.8 4.1 4.3 4.0
##  [865] 4.7 5.2 4.6 3.3 4.8 5.0 6.6 4.4 5.9 4.9 6.0 5.0 5.7 5.9 4.2 5.5 5.5 5.2
##  [883] 3.7 4.0 4.2 4.6 5.1 5.1 4.6 4.1 6.4 5.5 4.6 5.4 4.3 3.4 5.3 4.3 4.7 5.7
##  [901] 2.4 3.0 4.5 5.1 5.3 6.5 3.8 4.0 2.9 4.8 4.1 2.5 4.1 4.8 4.5 4.4 4.1 3.1
##  [919] 3.7 4.8 4.3 3.9 4.9 4.5 3.9 5.0 3.9 5.3 4.0 4.2 6.3 5.1 5.5 5.2 5.6 4.0
##  [937] 5.5 5.2 4.0 4.4 6.2 4.5 4.0 4.3 3.3 4.3 3.9 4.9 4.7 5.4 5.6 4.9 5.1 3.8
##  [955] 5.6 4.7 4.9 5.2 5.4 5.9 5.5 4.6 3.7 4.4 5.1 5.4 4.8 3.9 4.4 4.1 3.8 4.8
##  [973] 5.6 3.5 6.1 4.3 4.7 5.5 5.3 4.9 5.6 4.3 3.9 5.4 3.7 3.1 5.3 4.7 4.8 4.6
##  [991] 4.4 4.3 5.7 3.5 6.4 5.6 4.2 4.6 5.2 4.7
## 
## $func.thetastar
## [1] 0.0102
## 
## $jack.boot.val
##  [1]  0.52789318  0.36214834  0.21703911  0.24000000  0.01893491 -0.11203704
##  [7] -0.17031250 -0.19035088 -0.33342776 -0.46630137
## 
## $jack.boot.se
## [1] 0.91015
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
##    [1] 5.5 3.6 5.8 3.6 4.1 5.7 5.4 4.5 4.3 4.4 4.2 3.6 5.2 3.4 5.2 3.2 4.8 4.9
##   [19] 3.8 4.0 5.1 5.4 4.6 6.1 4.8 3.5 3.6 6.4 6.0 4.1 3.6 4.9 4.3 3.1 4.4 3.8
##   [37] 4.1 3.8 3.6 5.1 3.9 3.0 4.9 4.5 3.1 4.1 5.6 4.7 3.5 4.4 3.6 4.1 4.9 5.2
##   [55] 5.0 4.4 4.9 5.0 3.0 3.7 6.1 4.9 4.7 4.9 4.0 4.7 4.6 4.3 4.3 5.3 2.9 4.6
##   [73] 5.2 4.7 4.4 4.8 4.2 4.6 4.9 3.2 5.1 5.2 4.1 6.1 4.6 4.6 5.0 4.9 4.6 4.5
##   [91] 3.8 2.8 3.2 4.9 4.7 2.8 4.3 5.2 3.5 4.8 4.9 4.6 5.1 3.3 4.3 3.9 3.9 4.2
##  [109] 4.6 5.0 2.0 4.3 5.7 6.2 5.0 4.4 4.6 5.6 5.3 5.7 4.2 4.6 5.3 3.5 5.7 4.8
##  [127] 5.0 3.8 5.1 5.2 6.3 5.4 4.2 4.8 4.5 3.8 4.6 4.7 4.5 5.7 6.3 5.1 5.6 3.8
##  [145] 5.8 3.9 5.4 4.2 5.1 3.8 5.3 5.5 5.7 4.5 3.8 3.9 4.7 5.7 5.0 4.5 4.9 3.8
##  [163] 3.4 5.1 4.5 4.4 4.6 4.5 4.1 3.4 4.3 4.7 4.2 4.5 5.4 3.0 5.2 3.9 5.5 5.4
##  [181] 2.9 3.6 3.8 5.1 4.3 5.1 2.3 4.5 3.8 4.1 6.2 4.3 3.6 4.4 4.4 2.4 4.0 4.5
##  [199] 5.4 6.4 5.7 3.8 5.6 3.2 5.2 2.6 5.4 4.9 3.0 3.8 3.2 4.2 3.9 4.8 3.5 3.5
##  [217] 6.2 4.5 5.1 5.1 5.6 3.6 4.6 5.6 6.1 5.3 3.9 6.1 4.9 4.4 3.8 4.8 6.0 3.0
##  [235] 2.4 3.3 4.6 4.7 4.7 5.8 3.4 4.4 4.6 4.5 3.9 5.8 6.1 5.0 3.6 4.7 6.0 3.7
##  [253] 4.0 3.5 3.7 5.0 3.7 4.8 4.1 5.3 2.7 3.9 4.1 5.9 4.4 3.0 6.2 4.1 3.0 2.4
##  [271] 5.0 3.1 4.3 4.7 5.0 4.4 6.0 4.8 4.9 5.5 3.0 5.1 4.2 5.2 4.3 4.5 4.5 5.1
##  [289] 4.8 5.4 4.5 5.3 4.7 4.4 3.9 4.2 4.5 5.4 4.0 3.8 5.8 5.5 4.9 3.7 4.0 5.4
##  [307] 4.6 4.4 4.4 3.0 3.9 4.2 4.9 3.4 5.6 2.1 4.3 2.3 4.6 5.4 5.3 4.3 6.4 4.6
##  [325] 5.3 4.2 5.6 4.0 4.3 4.1 4.0 3.0 5.1 4.4 4.3 5.4 3.6 3.6 4.2 4.9 5.5 4.0
##  [343] 6.7 4.0 5.3 4.5 3.6 3.9 4.0 3.3 4.4 3.4 4.5 4.6 3.9 3.7 4.7 4.2 4.9 3.7
##  [361] 5.2 3.3 3.8 4.3 4.3 6.1 5.3 5.5 4.0 2.9 4.7 4.4 2.9 5.1 3.9 4.7 5.4 5.6
##  [379] 3.1 6.9 3.6 4.1 3.6 3.8 5.1 4.3 5.3 5.3 4.6 5.1 6.1 3.8 5.8 4.7 4.2 4.1
##  [397] 4.7 4.9 4.1 4.8 3.6 5.2 4.5 5.2 4.5 4.3 2.9 5.3 2.3 5.1 3.6 4.3 3.5 4.3
##  [415] 3.5 4.9 6.1 3.9 5.0 3.1 2.6 6.5 5.2 4.6 4.2 3.9 4.2 5.1 4.6 2.7 3.8 4.4
##  [433] 5.5 5.8 4.7 3.6 3.7 5.0 4.0 4.0 4.0 5.1 3.0 4.9 3.8 6.5 3.9 2.0 3.4 4.0
##  [451] 4.3 5.3 4.7 5.2 4.3 4.5 3.5 5.5 5.2 3.3 5.1 5.4 3.3 6.7 4.4 3.6 3.5 4.0
##  [469] 5.6 4.5 4.9 5.2 4.2 5.3 1.8 3.7 4.3 3.9 4.3 3.3 3.4 5.7 5.5 3.4 4.2 5.6
##  [487] 5.0 4.2 4.8 4.4 4.4 6.3 4.0 5.1 4.3 3.3 3.5 4.2 5.2 4.0 4.7 3.9 3.4 3.0
##  [505] 2.7 4.4 5.3 4.5 7.0 4.4 3.8 4.3 4.1 3.9 3.2 4.4 3.9 3.4 5.5 4.3 5.2 2.6
##  [523] 4.7 3.9 4.3 4.4 4.1 5.1 2.5 5.4 4.3 4.9 2.8 4.8 3.2 6.2 3.8 2.8 4.8 5.2
##  [541] 5.3 4.0 4.2 4.5 3.6 2.9 5.5 3.1 4.6 4.4 5.6 3.6 5.6 4.3 4.7 3.3 3.9 5.0
##  [559] 4.4 4.5 6.5 3.6 3.9 4.3 4.0 5.3 5.7 2.7 4.7 5.5 3.0 3.7 4.7 5.0 3.8 3.8
##  [577] 4.7 4.5 5.7 4.6 3.8 5.9 3.7 4.8 4.8 4.1 4.7 3.6 4.7 4.3 4.3 4.5 4.4 5.8
##  [595] 4.3 4.8 5.4 6.3 4.5 5.5 4.5 5.1 3.5 3.9 4.4 3.6 5.9 3.9 5.7 6.1 4.6 3.7
##  [613] 4.5 4.4 4.4 5.3 4.0 4.9 3.4 4.1 4.9 5.1 6.5 3.7 4.1 5.6 4.6 4.8 5.4 5.4
##  [631] 5.4 4.6 4.0 4.2 2.7 4.3 4.3 4.2 4.4 4.2 6.0 5.4 3.7 3.5 4.1 4.5 5.6 3.9
##  [649] 3.6 3.5 3.2 5.2 4.1 4.7 5.2 5.8 5.8 4.2 4.3 4.8 3.5 5.1 5.8 4.0 4.6 3.7
##  [667] 5.6 2.9 5.5 5.1 2.8 4.3 4.5 3.9 4.8 3.4 3.6 3.9 5.6 3.4 4.4 3.1 5.7 6.1
##  [685] 4.2 4.8 4.3 4.0 4.8 4.8 3.5 5.6 4.0 5.2 5.0 4.3 4.9 3.7 4.4 3.2 3.8 4.5
##  [703] 3.5 4.7 6.4 3.2 5.8 4.5 3.9 4.7 4.0 4.0 5.2 3.3 6.2 4.9 4.7 4.9 3.7 4.7
##  [721] 3.0 4.6 5.2 4.8 3.3 5.3 4.6 4.5 4.8 5.2 4.2 4.9 4.2 4.4 3.1 4.6 6.4 5.3
##  [739] 3.8 4.2 4.6 6.3 4.2 6.0 5.3 4.3 4.5 2.7 3.9 5.8 5.4 4.8 5.6 3.8 4.5 4.7
##  [757] 4.1 4.2 3.5 3.6 5.1 3.6 3.8 3.9 4.0 4.1 6.0 3.7 4.5 3.7 5.0 3.2 3.6 4.6
##  [775] 4.5 5.2 4.4 4.0 3.2 6.3 5.0 5.8 3.7 4.1 5.9 5.3 4.2 5.2 6.8 5.2 4.1 4.1
##  [793] 3.5 3.7 5.8 5.0 4.4 5.3 4.4 3.9 5.1 3.9 3.7 5.0 4.2 5.3 3.5 4.3 4.8 5.4
##  [811] 3.9 4.5 2.8 4.7 3.3 4.5 4.9 3.3 4.4 4.8 4.1 4.0 6.4 3.8 5.3 4.6 3.5 3.6
##  [829] 4.9 3.1 3.7 5.0 4.9 5.9 3.1 3.7 6.8 4.6 3.5 6.0 4.5 3.2 2.0 4.3 5.8 6.7
##  [847] 3.8 5.4 4.5 4.7 5.2 4.6 5.3 5.7 4.5 4.3 2.9 5.4 3.6 5.2 2.7 3.8 4.1 4.0
##  [865] 4.3 5.6 6.5 4.0 6.1 4.2 5.0 6.1 3.4 6.7 3.5 5.5 3.5 4.2 5.6 4.3 3.5 3.8
##  [883] 3.9 4.1 3.7 4.2 4.7 4.2 6.5 5.3 4.8 4.2 3.1 4.9 3.4 4.1 4.7 3.2 4.1 6.3
##  [901] 6.5 4.2 3.4 4.3 4.0 4.4 3.8 3.5 3.3 4.6 3.5 4.3 4.8 6.0 5.1 4.0 6.3 3.7
##  [919] 2.9 3.4 4.1 3.7 6.2 3.9 3.5 4.9 3.9 4.7 4.9 3.5 4.1 3.3 4.8 3.8 5.5 3.8
##  [937] 6.8 3.2 4.6 3.6 4.6 4.5 3.7 6.0 5.5 4.1 4.1 5.9 3.2 4.9 4.4 4.5 5.5 3.4
##  [955] 3.5 3.4 6.2 3.4 4.8 5.6 4.1 5.2 2.3 3.1 5.4 3.9 4.2 4.9 4.5 4.5 3.6 4.2
##  [973] 3.8 5.9 3.4 6.1 2.4 5.1 3.7 6.1 3.7 4.0 3.3 3.4 3.9 4.6 4.3 3.2 3.7 2.7
##  [991] 4.5 4.5 5.8 4.3 3.6 5.6 5.1 3.8 4.0 4.1
## 
## $func.thetastar
## 72% 
##   5 
## 
## $jack.boot.val
##  [1] 5.5 5.4 5.4 5.2 5.1 4.9 4.7 4.7 4.5 4.4
## 
## $jack.boot.se
## [1] 1.128893
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
## [1] 0.1044155
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
##   4.755766   7.094562 
##  (2.056587) (3.235868)
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
## [1] 1.1742121 0.2195052 0.3579725 0.3350293 0.8370648 0.5097721
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
##    [1] -0.4258298960 -0.4753964226  0.0138890291  0.2054585221  0.0459374117
##    [6]  0.6391189000  0.1251902587  0.5446645108  0.2371353820  0.3334945946
##   [11]  0.4295700524  0.5837777166  0.6552055094 -0.5082407577  0.4064691706
##   [16]  0.5444660286 -0.3098764534  0.2866867815  0.8614899778 -0.3673447916
##   [21]  0.1247274290  0.3528321675  0.2325992536 -0.3890342899 -0.0503945240
##   [26]  0.0332391584  0.1554713014  0.7840128538 -0.1121166587 -0.3686913956
##   [31] -0.5184361432  0.5333724808 -0.0333198543 -0.6309037551 -0.0651179227
##   [36] -0.3994530687  0.4148597109 -0.6618606874 -0.2999523901  0.2644019843
##   [41] -0.0578653412  0.0231768859  0.0824565624 -0.0493400716 -0.1463137241
##   [46] -0.9470060006  0.2645816971  0.4804580318  0.0461478748 -0.2542204293
##   [51] -0.2847786231  0.3496722390 -0.4108053558  0.5554074415  0.2792855130
##   [56] -0.1571790493 -0.2731683840  0.5421091085 -0.7780678255  0.4830934185
##   [61] -0.1382388871  0.3586480208 -0.0365638914  0.1917559945 -0.3043638178
##   [66] -0.6756557553  0.2546430246  0.2371455589 -0.4090919850  0.2829884573
##   [71]  0.0527819513 -0.6558235924  0.2490090547 -0.2657907967  0.6293509690
##   [76] -0.2885944425  0.0231768859 -0.3770584221  0.0345948456  0.0731997545
##   [81] -0.1696986014  0.1325669332 -0.1452606553  0.1210616212  0.2994769135
##   [86]  0.0969970204 -0.3042395057 -0.1687997875  0.1651875719 -0.2954923356
##   [91]  0.2432597331  1.0334594929  0.6255415182  0.2704548620  0.1511311324
##   [96]  0.3468402247 -0.1276080598  0.2469739951 -0.1025779787  0.2883604589
##  [101] -1.5359465601  0.2512165608  0.2575250344  0.1305512752 -0.1539057537
##  [106]  0.0217807146  0.7383118295 -0.0894531878 -0.2548959720  0.6946337991
##  [111] -0.1466482747  0.7384367696 -0.5037633547  0.2449226053 -0.0493400716
##  [116] -0.0992197130  0.4205387516  0.0483674474 -0.8318730840 -0.0576498551
##  [121]  0.3549728101  0.2879694286  0.4714668474 -0.1104250617  0.0708062475
##  [126] -0.6010307576 -0.3397205618  0.1819365149  0.2255000884  0.8253091484
##  [131]  0.5548446342  0.0611972170 -0.0218849958 -0.1822238538  0.0031674837
##  [136] -0.4363577233 -0.0767254804  0.5278440534  0.4889795940  0.3906063621
##  [141]  0.3719015062  0.1625139558  0.1044154897 -0.0402389295 -0.3080994706
##  [146] -1.4562072579 -0.4291583264  0.0442188971 -0.3171729678 -0.5401161249
##  [151] -0.3799996948 -0.0106734408 -0.3862272935  0.4174575688  0.6531498763
##  [156] -0.5504202067  0.5951201060 -0.7147102682  0.5807683997 -0.0383340829
##  [161] -0.6334879129  0.1154919889 -0.1231558945 -0.0671253363  0.1875340482
##  [166]  0.2830819492 -0.5200599303  1.0857471077 -0.2301390301 -0.4068921666
##  [171] -0.1556926124  0.2613263239  0.1861863859  0.0966160556  0.0362033412
##  [176]  0.2202508330  0.1305806704  0.1738335609 -0.3215593499  0.1670802067
##  [181]  0.0855889512 -0.4143121861 -0.6406579461 -0.3218807361  0.1907980531
##  [186] -0.4203825265  0.4111133555 -0.1327529389  0.1521376391 -0.0080476675
##  [191]  0.0271231513  0.5979811441  0.7613174498  0.3230654768  0.0713237386
##  [196]  0.1746292483  0.0128735964 -0.2036729385 -0.1824138172 -0.3055740847
##  [201]  0.2548276243 -0.5264717774  0.4112421387 -0.7059311218 -0.4974563065
##  [206]  0.2690809462 -0.2185143949 -0.0485163632 -0.1873102800  0.5030962296
##  [211] -0.1044895677 -0.4105004210  0.1767818136  0.2984810232  0.4334829006
##  [216] -0.2548686589  0.4224641276  0.2496957365 -0.7734693685  0.1690658409
##  [221]  0.2259741095 -0.2847514414 -0.3809507029  0.5943633071 -0.0678253588
##  [226]  0.0278101650  0.7594439580 -0.5465495411  0.1377358663  0.1856691322
##  [231] -0.4865496067 -0.0740836140 -0.7652135559 -0.1776348334 -0.5957918174
##  [236]  0.4114619179 -0.1760149625  0.6336756381  0.2071146028  0.1794665884
##  [241]  0.7205130121 -0.3492115936  0.0771252475 -0.5421384816 -0.5098429397
##  [246]  0.1006334991  0.0068353020 -1.0285894599  0.1195042312 -0.3867080607
##  [251] -0.3128531270  0.3484494260 -0.1423190874  0.5650683584  0.1136161978
##  [256] -0.6135762820  0.0809540033  0.9926698317  0.0460150229 -0.3061219310
##  [261]  0.2322215260  0.0179055143 -0.1849062079  0.2764268653 -0.7220613834
##  [266]  0.2491437933 -0.3686360665  0.0439895107  0.0186634104  0.4811180004
##  [271]  0.1989567484 -0.4596835077  0.5031634675  0.4172937908  0.3065857266
##  [276]  0.3349976407 -0.6048352616 -0.1608563393  0.2117319165 -0.2363029146
##  [281]  0.0788006748 -0.1421764459 -0.0500351128 -0.3692664470  0.0716993740
##  [286]  0.8153062036 -0.5184361432  0.1321017382 -0.6273193410  0.7527345863
##  [291]  0.2372280672 -0.0863339018 -0.2421655664  0.3744705020  0.0226979931
##  [296]  0.5624440131  0.9160665668  0.3948375262  0.8063965077 -0.0701852775
##  [301]  0.2436374484 -1.0953209105  0.7770999066  0.3304081121  0.0384933035
##  [306] -0.1961313373 -0.2029660174 -0.3303570433 -0.3056076126  0.4362618076
##  [311]  0.3042337535  0.1823764988  0.2934968434  0.3335056835  0.4484298346
##  [316]  0.1568141171 -0.1287237103 -0.1521704747 -0.0213736042  0.0956961587
##  [321]  0.2829632378  0.4098034272  0.8153062036  0.4892082229  0.0761274670
##  [326] -0.1152147125  0.5324682707 -0.0748872952 -0.5934565326  0.2646554995
##  [331] -0.2744483854 -0.6796232778 -1.6343579872 -0.3332729017  0.9584575119
##  [336]  0.1753346240 -0.1388454609  0.0330038250 -0.0123796435  0.1495598613
##  [341] -0.0262126242  0.0830740532  0.1407312215 -0.7176259436 -0.2625507473
##  [346] -0.3098464325  0.6348095787  0.5489026330  0.2332451366 -0.6401970265
##  [351]  0.3560997650  0.1053225912 -0.4501307569 -0.9201434238  0.0007014611
##  [356]  0.6094684476 -0.1559082998 -0.5561778506  0.6074556394 -1.0696469669
##  [361]  0.5962124564 -0.2432830271  0.6304847450  0.2070416355 -0.4162335329
##  [366] -0.1424399961 -0.0891810174 -0.5782218563 -0.4450012459  0.1176534147
##  [371]  0.4992655636  0.3877446473 -0.1747853358 -0.3832584359  0.0140266875
##  [376] -0.0340252370  0.4574999404  0.1869710339  0.2012982684 -0.0236365008
##  [381]  0.6494100006  0.0167012135  0.4638942244 -0.1041319661  0.8236467288
##  [386] -0.0548092311 -0.5968626301 -0.1315613494 -0.9838337625  0.0468422424
##  [391]  0.3034081798  0.4285007361  0.2166529796 -0.8138858983  0.0216240836
##  [396] -0.2574106331 -0.0801334480  0.4249482124  0.7507092760 -0.2612696654
##  [401]  1.3284583435 -0.0022387786 -0.6411283527 -0.5748507893  0.1752821939
##  [406]  0.3792750423  0.4454679301  0.8503724146 -0.0017297509  0.0609470752
##  [411] -0.9535331216 -0.2446508890 -0.1636425555  0.6861807229  0.7089947966
##  [416]  0.4268860782  0.0831347669  0.0245083747  0.0172436274 -0.9150073877
##  [421] -0.3197321600 -0.5542006434 -0.2976299452  0.6520528104  0.2428898905
##  [426]  0.2746043725  0.2475054622  0.5606980651  0.2742609787  0.7844618896
##  [431] -1.8263201812  0.4351733648  0.3714203764 -0.1509593898 -0.0051424748
##  [436]  0.1296510973  0.0993818399  0.4730229892 -0.4989587227 -0.5952347729
##  [441]  0.6623798913  0.1225019544 -0.2221276431  0.4185632980  0.5316353885
##  [446]  0.1044874208  0.5488776796  0.3053604174 -0.0404453101  0.0349202048
##  [451] -0.2691153442  0.1321017382 -0.4377713862  0.0186811330  0.5777576525
##  [456]  0.2649881515  0.5738989686  0.0862494278 -0.3044287562 -0.3342042729
##  [461] -0.5975241464  0.2428985196 -0.4653875193  0.8703334778 -0.8066967074
##  [466]  0.3512334122  0.0417469180  0.5874248757  0.0256911229 -0.4612414006
##  [471] -0.5007915944 -0.7730453450  0.7326895309  0.0766105288 -0.0342507638
##  [476]  0.3430367579 -0.1538900452  0.3598202317 -0.1499696252 -0.4422670611
##  [481]  1.4004001603 -0.2889163265 -0.3372219777 -0.0937203950  0.4664172432
##  [486]  0.0353904711  1.1060940476 -0.2246102322 -0.1114329462 -0.5109157134
##  [491]  0.4080984128  0.2634334541  1.4303508379 -0.4301767132  0.4418876618
##  [496] -1.2354436367  0.3183443739  0.0564248857  0.4800492107  0.2171750673
##  [501]  0.3447231191 -0.2937292772 -0.5944218507 -0.1164827052 -0.1046329650
##  [506]  0.3140338213 -0.1855595559 -0.3801999303  0.3241705953 -0.6072479410
##  [511] -0.1271544944 -0.1169642546  0.6692444933 -0.1527059749 -0.2427510155
##  [516] -0.7471023871  0.3614554979  0.0913477797  0.6208178311  0.2174492991
##  [521]  1.3328831243  0.4777144444 -0.5435247085  0.0458371919 -0.4144801109
##  [526] -0.0016808275  0.4737201693  0.2027479265  0.0231791537  0.1330488037
##  [531]  0.0885103635  0.3000504936 -0.7481895065 -0.4821417648  0.6033741640
##  [536]  0.7088140090 -0.6726587319  0.1220758759 -0.2872014795 -0.0649695046
##  [541]  0.3752443353  0.5693753665 -0.1890388764  0.5360118206  0.5618328394
##  [546]  0.1113201878  0.1710693744 -0.3906062003 -0.1361550375  0.5183152267
##  [551]  0.3703491363 -0.2251230719  0.1302075177 -0.0509286042 -0.0892662489
##  [556]  0.0616749931 -0.1332499850  0.2208751067  0.1074237609 -0.1066627800
##  [561]  0.3773404818  0.4828026351  0.7213821823 -0.2429667765 -0.6184309850
##  [566]  0.2247260150 -0.2393802304 -0.3120293967  0.1358578820  0.0808253729
##  [571]  0.3405717206 -0.1703285291  0.6380270445 -0.3479086780  0.7811714523
##  [576] -0.0449727438  0.6113018058 -0.3634694520  0.0880334039  0.0771252475
##  [581] -0.2556547812  0.4212807288  0.0486266808  0.7713951619 -0.0772072359
##  [586]  0.0260160568 -0.0226402568  0.0119400940  0.0714640449 -0.1558075371
##  [591] -0.0642805824  0.5701324592  0.0043743975 -0.3691321465  0.9475303628
##  [596] -0.0920334232  0.3286239386  0.7112243312 -0.4115787493 -0.2473902145
##  [601]  0.4346404481 -0.0188041066 -0.2239763895 -0.1550291077 -0.0722819491
##  [606]  0.3164669028  0.0771716087 -0.6465920968 -0.1296933024  0.2414564041
##  [611] -0.0063241135  0.4163995944  0.0107185166 -0.0568516295  0.9058583868
##  [616]  1.1314058616 -0.4351102811  0.4416093396 -0.0196698613  0.2258855297
##  [621]  0.4650446297  0.3956517432 -0.2835754819 -0.5443884890  0.0580275484
##  [626] -0.0452794917 -0.4392512142  0.1452464958  0.1428795080 -0.1615509441
##  [631]  0.3194820582 -0.5145259628 -0.2022783577 -0.2109844589 -0.6884482765
##  [636]  0.2428073430  0.9111680529 -0.3119743707  0.5400713829 -0.1339533774
##  [641]  0.5761030307  0.6056103466  0.1447705071 -0.2091720099  0.2390679248
##  [646]  0.3575859262  0.2663517037  0.0204584237  0.2409356993 -0.4365637911
##  [651]  0.0060410120  0.8107313051  0.0784851193  1.1894573512  0.2214717681
##  [656] -0.0238214628  0.2213444888 -0.2217567217 -0.3743233165  0.4587636887
##  [661]  0.3751545850 -0.0926433612  0.7576082241 -0.3418052436  0.1754984637
##  [666]  0.7448801580  0.2471577275  0.0201231930 -0.2236136303  0.0521565728
##  [671]  0.1298363759  0.5747818818 -0.1221342994 -0.1101717057 -0.0213736042
##  [676]  0.4711335136 -0.2997191568  0.5221735044 -0.1302310642 -0.0206063375
##  [681]  0.0491927800  0.2663638818 -0.1202677621  0.0222826619 -0.2866804346
##  [686] -0.5129960576 -0.7311829172 -0.2173131943  0.1388285656 -0.2888486286
##  [691] -0.0676748357  0.2465068908  0.1698013842 -0.0951912115  0.1200697448
##  [696]  0.5920245343 -0.4312870339  0.0265068859 -0.4189516586 -0.1250643152
##  [701] -0.4937745176  0.2825081038  0.2289786210  0.1646976587  0.6290549279
##  [706]  0.4200781262  0.2363303867  0.1193614565 -0.4168032031 -0.3090389934
##  [711] -0.4762248806  0.7427556270  0.4500261313  0.3739363181 -0.6675541565
##  [716]  0.0337964065  0.4273658924 -0.5562975054  0.3621648978 -0.0359563270
##  [721]  0.0285912255 -0.3293263938 -0.1995698717 -0.0540416486  0.8333267068
##  [726] -0.5074234800  0.5116943301 -0.3914136253  0.8209703536  0.6379970570
##  [731]  0.0550787313 -0.1766210361 -0.1704647171  0.4945221553  0.2149580969
##  [736]  0.3587344358  0.2690809462  0.8248846613 -0.7414942510  0.2787270122
##  [741] -0.2414268887 -0.5942385946  0.2107027742 -0.0961031805  0.1781342088
##  [746]  0.0021741011  0.2753481028 -0.1126349050 -0.0638425805  1.1670516090
##  [751] -0.3018441574  0.6599805185  0.3998455186  0.4420176216  0.2451712020
##  [756]  0.0618025590  0.0298424035 -0.1715566057  0.7321144908  0.3196882798
##  [761] -0.2639222518 -0.0517770603  0.1087217778  1.1897251538 -0.6630332846
##  [766] -0.0009812504 -0.0608239008 -0.4543463678  0.3026874346 -0.1277600012
##  [771] -0.1449242473  0.0698695229 -0.3482145982  0.2881677769 -0.2019274538
##  [776]  0.1671519626 -0.4239068391  0.2068268018  0.9749166205 -0.2854444866
##  [781]  0.5531377808 -0.1139453921  0.1324610110 -0.4267243196  0.2841489123
##  [786]  0.3306294336 -0.2405660981  0.3564206349  0.3543195754 -0.2981653964
##  [791] -0.5834619412 -0.0900073699  0.2574022247  0.7991616221  0.3238622814
##  [796]  0.7378584353 -0.6553772008  0.3460278049  0.7611014773  0.4663817697
##  [801]  0.3949534014  0.3293330106 -0.5248879987 -0.2526054243  0.2423302441
##  [806] -0.2919574831 -0.0604354241  0.4819565100 -0.0210589570  0.3952737797
##  [811]  0.0635548689 -0.1236025048 -0.3438097257  0.7583505566 -0.3482145982
##  [816]  0.2754079738 -0.2472535340  0.2850602074  0.7618965657 -0.0823052213
##  [821] -0.6586056665  0.2787129160 -0.3068876858 -0.0354001262  0.3297404347
##  [826] -0.1778678324  0.2736666163  0.6367219624 -0.1713908745  1.0092517591
##  [831]  0.5666855029  0.3218621646 -0.1026923331 -0.0997520479  0.0013392379
##  [836]  0.3424088017  0.1797436541  0.5458874217 -0.0403493494 -0.0624127354
##  [841] -0.8895577830 -0.2435078152  0.3374395913  0.3152629731  0.4132154534
##  [846] -0.0382286014  0.3081808648  0.3019065714  0.2116779498 -0.0672530828
##  [851]  0.6501883749  0.5406068651 -0.0764563836  0.2049214786  0.4405602370
##  [856] -0.4110760478  0.9342128488 -0.0627511696  0.0775289034 -0.1601223385
##  [861]  0.7064032493  0.1224798801  0.3145553363  0.1301461235 -0.1038523032
##  [866] -0.0100756543  0.0077432990 -0.4766003488 -0.1692503441  0.2374698404
##  [871] -0.1076911052 -0.2625229441  0.9619957269  0.1360790698 -0.2678647559
##  [876]  0.6335462131 -0.0044288488 -0.0951912115 -0.4918601387  0.4501209189
##  [881]  0.2758184495 -0.0951912115  1.0008970941 -0.4433193776  0.3570563812
##  [886]  1.5578114808 -0.2006523848 -0.0319462877  1.0332655114  0.9812154120
##  [891]  0.0862494278  0.0871396597 -0.0821292855 -0.5439035494  0.0616749931
##  [896] -0.1158463032  0.0007401171 -0.1688133137  0.8212867230 -0.0916995388
##  [901]  0.2027923684  0.6079588039  0.2081719646 -0.0599714137 -0.7352336945
##  [906]  0.4814441873 -0.0518178585 -0.6818022651  0.4583675815 -0.1231690529
##  [911]  0.6215860382 -0.0402775770 -0.2648663750 -0.4490648406  0.8748827517
##  [916] -0.4332787942 -0.1043239170  0.3993323944 -0.0684948543  0.0548420365
##  [921] -0.2046837276  0.2090563447 -0.1888357803  0.2053997560 -0.0811091420
##  [926]  0.4448705602  0.0024886081  0.6415319763 -0.0119563121 -0.4878803144
##  [931]  0.4396980324  0.5691995394  0.6419428720 -0.3100861692  0.1021267764
##  [936] -0.3870375152 -0.9070401866 -1.1286000547  1.0255100138  0.9326177757
##  [941] -0.0801292385  0.3845123117 -0.1737831698 -0.2187887560 -0.3683821880
##  [946]  0.1635722860 -0.1529602401 -0.0635019427 -0.4651369626 -0.1208114758
##  [951]  0.1942765092  0.3465979575  0.2406773191 -0.0322457525  0.9169751826
##  [956]  0.0666587415  0.9740670862  0.6786656603 -0.2323562526  0.8529575834
##  [961]  0.2141410473  0.3869642032  0.2934439685  0.2932986531  0.3425082006
##  [966]  0.2967743443  0.0461509923  0.1341158200  0.4875981505 -0.0317996430
##  [971] -0.0030031715 -0.1291468774 -0.4328381167  0.0402002260 -0.0339629870
##  [976]  0.2371455589 -0.7412912682 -0.5772342439  0.5972000453 -0.2516105524
##  [981] -0.1540635244  0.7169926476 -0.0966204219  0.0407031273  0.2392269053
##  [986] -0.3211206587 -0.3966131969 -1.6515490415  1.0033786304 -0.5025192856
##  [991]  0.3347753909  0.8876854078 -0.1602369226 -0.3291587357  0.3519264089
##  [996] -0.6762821755  0.1499254679  0.0500135394 -0.0138074595 -0.4851360506
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
##   0.67033820   0.28818872 
##  (0.09113328) (0.06443794)
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
## [1] -0.25472406 -0.74386594  0.57680141 -0.06922391  0.28113528  0.16001440
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
## [1] 0.0101
```

```r
sd(results2$thetastar)  #the standard deviation of the theta stars is the SE of the statistic (in this case, the mean)
```

```
## [1] 0.9228324
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
## t1*      4.5 -0.01221221   0.9178212
```

One of the main advantages to the 'boot' package over the 'bootstrap' package is the nicer formatting of the output.

Going back to our original code, lets see how we could reproduce all of these numbers:


```r
table(sample(x,size=length(x),replace=T))
```

```
## 
## 1 3 6 7 8 
## 5 1 1 1 2
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
## [1] 0.0073
```

```r
se.boot
```

```
## [1] 0.9227395
```

Why do our numbers not agree exactly with those of the boot package? This is because our estimates of bias and standard error are just estimates, and they carry with them their own uncertainties. That is one of the reasons we might bother doing jackknife-after-bootstrap.

The 'boot' package has a LOT of functionality. If we have time, we will come back to some of these more complex functions later in the semester as we cover topics like regression and glm.

