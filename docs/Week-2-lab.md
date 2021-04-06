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
## 0 2 3 4 5 6 
## 1 1 2 2 1 3
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
## [1] -0.0274
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
## [1] 2.692715
```

```r
UL.boot
```

```
## [1] 6.252485
```

Method #2: Simply take the quantiles of the bootstrap statistics


```r
quantile(xmeans,c(0.025,0.975))
```

```
##  2.5% 97.5% 
##   2.7   6.4
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
##    [1] 4.6 4.2 4.8 2.9 4.7 4.2 5.2 3.8 3.4 3.5 4.2 5.0 3.4 5.2 4.3 4.9 5.7 3.8
##   [19] 4.6 4.8 3.0 3.7 5.1 3.2 4.5 4.2 4.6 3.7 5.0 5.2 5.0 4.1 4.1 6.5 5.5 5.5
##   [37] 4.6 4.9 4.1 5.3 4.9 4.9 4.3 4.3 4.0 4.0 6.1 3.7 3.8 4.5 4.2 4.0 3.9 3.1
##   [55] 3.2 5.1 3.6 3.9 4.0 5.1 6.4 4.7 3.9 4.1 3.5 3.9 4.7 4.7 2.6 2.7 2.6 4.1
##   [73] 5.7 5.1 2.7 4.1 6.0 4.0 3.4 3.6 4.6 4.3 4.6 4.5 4.3 4.6 5.0 4.0 3.7 4.8
##   [91] 4.6 4.5 4.1 4.1 3.1 5.8 5.8 5.9 3.8 3.5 3.3 2.9 3.6 4.7 3.5 5.9 3.2 4.8
##  [109] 5.2 3.5 4.8 5.9 4.6 5.0 6.3 4.6 4.7 3.4 4.2 3.4 4.4 5.8 4.7 4.3 5.8 3.1
##  [127] 4.1 4.3 3.5 5.1 4.7 5.4 4.2 4.1 6.4 5.8 4.1 5.3 5.0 4.2 5.3 2.9 2.6 6.1
##  [145] 5.2 4.9 5.2 3.8 6.2 5.5 3.3 4.4 4.3 4.5 3.6 4.7 3.7 3.8 5.7 5.1 5.4 5.8
##  [163] 3.7 5.8 4.9 2.9 3.9 5.5 4.3 4.8 4.7 4.8 3.8 3.7 4.1 5.0 4.9 5.1 5.0 5.5
##  [181] 4.7 4.7 4.6 2.0 5.5 5.0 3.5 3.7 4.2 6.6 3.5 3.9 4.8 4.0 5.8 3.3 5.5 4.9
##  [199] 4.6 5.2 5.5 3.2 2.9 4.0 2.7 3.0 3.9 5.8 4.3 4.6 4.8 4.1 3.8 3.5 4.2 4.8
##  [217] 5.5 4.8 5.4 5.2 5.2 5.0 4.7 5.6 4.5 6.0 6.7 3.6 4.4 4.5 4.5 4.9 5.0 4.1
##  [235] 5.1 3.5 3.3 2.7 3.6 4.6 4.2 5.9 5.7 5.6 4.8 3.6 4.8 4.1 3.7 4.2 5.1 2.4
##  [253] 5.4 4.3 5.4 3.5 3.7 4.3 3.7 5.8 5.0 5.1 3.5 5.4 5.4 4.0 5.5 4.7 4.1 5.4
##  [271] 3.8 5.2 4.9 4.5 5.2 3.0 2.9 4.0 3.6 3.0 4.5 3.6 5.8 4.2 5.1 4.9 4.9 3.8
##  [289] 2.9 2.0 4.6 4.2 3.8 4.0 5.1 3.8 6.8 4.2 4.7 5.4 4.5 4.5 3.7 4.8 3.6 3.6
##  [307] 5.7 4.6 4.8 2.8 3.6 4.8 4.1 2.9 3.9 4.2 5.4 3.7 3.5 3.2 3.5 4.2 4.7 5.0
##  [325] 4.8 3.5 3.4 5.0 5.7 4.6 4.7 4.2 6.0 5.2 5.0 4.2 4.8 4.0 4.7 4.2 4.6 4.5
##  [343] 5.4 5.3 5.8 6.0 4.3 4.4 3.6 4.6 4.3 4.1 4.9 6.0 5.3 4.2 4.4 6.7 3.7 4.7
##  [361] 5.3 4.3 5.4 5.5 3.4 4.3 5.2 4.7 4.2 4.8 4.9 3.3 3.5 5.9 3.8 3.8 7.4 3.7
##  [379] 5.3 3.2 5.6 4.4 2.5 4.3 5.3 2.1 4.0 4.1 5.3 5.1 5.0 5.0 3.0 3.9 4.4 4.2
##  [397] 2.7 4.8 5.8 5.3 6.3 4.9 2.9 5.3 3.8 4.5 5.2 4.2 2.3 6.6 3.6 5.2 4.7 4.0
##  [415] 4.2 3.7 2.8 4.4 3.7 4.2 4.0 3.2 3.4 5.6 3.1 4.4 3.6 4.8 3.8 2.3 3.6 4.4
##  [433] 3.9 5.3 5.0 5.3 4.2 4.2 5.5 3.1 4.9 4.8 5.2 4.4 5.0 5.6 3.8 4.7 3.8 4.4
##  [451] 4.5 4.0 4.0 5.6 4.8 4.2 3.1 2.8 5.4 4.4 4.3 4.0 4.1 4.7 5.8 5.7 4.0 4.2
##  [469] 7.1 4.2 6.4 6.5 5.2 4.7 4.0 4.5 3.4 5.7 3.7 5.4 4.4 5.2 3.0 5.1 5.4 4.2
##  [487] 4.2 4.1 5.9 5.3 4.2 3.9 4.5 4.6 3.5 5.4 6.0 5.2 4.3 3.4 5.1 5.0 3.8 2.4
##  [505] 4.4 5.1 6.8 2.9 3.7 5.8 4.6 4.0 3.8 4.9 5.5 4.2 4.3 3.6 4.6 3.5 4.8 3.8
##  [523] 5.0 3.7 4.3 5.7 3.2 4.1 3.6 5.2 3.8 5.0 4.6 3.4 5.0 3.8 4.5 4.4 5.1 4.5
##  [541] 4.3 3.0 5.0 4.0 5.3 4.3 3.6 4.9 3.2 5.4 2.9 3.4 4.4 6.6 4.7 5.3 4.3 3.5
##  [559] 4.1 6.0 4.4 3.6 5.0 2.8 3.1 5.0 4.3 5.9 4.0 5.4 4.4 5.2 4.2 5.3 3.9 3.1
##  [577] 4.9 4.3 3.5 3.8 3.8 3.8 4.7 5.1 4.7 4.9 3.9 3.8 4.8 3.7 5.3 5.9 5.6 4.7
##  [595] 5.3 4.7 4.4 3.3 3.2 4.7 4.3 4.9 4.2 5.1 3.9 3.4 4.6 3.8 3.8 5.6 5.3 5.6
##  [613] 3.3 3.3 2.8 4.7 5.9 3.6 6.1 4.9 3.3 3.7 6.1 3.8 1.8 5.3 3.7 5.2 4.7 3.9
##  [631] 5.0 5.2 3.0 5.6 2.5 5.5 4.0 2.5 4.3 3.3 5.8 3.7 5.1 3.4 4.0 3.0 5.2 4.0
##  [649] 5.6 3.7 4.8 4.9 3.4 5.3 4.4 3.8 4.8 5.0 4.0 3.8 3.3 5.3 4.2 5.5 3.8 4.4
##  [667] 5.0 2.8 3.6 4.0 5.9 4.9 3.9 3.8 3.8 4.3 6.5 5.0 4.9 3.5 4.1 5.1 3.7 4.0
##  [685] 4.3 3.8 3.6 4.7 3.6 4.6 6.3 3.2 4.5 4.3 3.7 4.1 4.0 4.7 5.3 2.6 4.2 4.3
##  [703] 5.3 4.0 5.5 4.9 5.0 3.2 4.2 5.0 3.7 3.8 5.0 3.4 5.7 4.3 3.9 5.5 5.4 4.5
##  [721] 4.7 4.0 4.2 5.4 4.1 5.9 5.0 5.4 4.3 5.1 5.2 4.2 4.2 5.2 5.1 3.0 3.4 5.1
##  [739] 5.6 5.3 4.8 4.4 3.5 4.0 4.7 4.0 4.2 4.6 5.9 4.0 4.4 4.4 5.8 5.1 3.7 3.3
##  [757] 3.3 3.8 3.8 6.1 3.9 3.9 4.7 3.3 5.3 3.1 3.8 4.9 4.9 3.4 5.0 2.6 3.8 5.1
##  [775] 4.1 2.7 4.4 4.4 4.0 3.9 3.9 4.8 5.1 4.5 4.3 2.8 5.0 3.2 3.8 3.5 5.4 5.4
##  [793] 3.8 3.4 5.0 4.5 4.9 4.6 3.6 4.4 4.8 4.8 2.5 5.1 5.3 4.9 3.7 6.0 4.8 4.7
##  [811] 5.1 4.5 5.9 3.5 2.7 5.3 3.8 5.5 5.3 4.7 5.5 3.0 3.8 3.0 4.8 3.8 5.3 5.6
##  [829] 4.7 4.3 5.8 2.6 4.0 3.8 4.4 4.8 4.7 3.8 4.4 4.6 2.7 5.7 3.8 4.7 3.8 4.1
##  [847] 4.0 4.6 3.0 5.4 3.5 4.9 4.7 5.2 5.0 4.6 4.9 3.1 6.3 4.2 4.7 4.1 3.6 5.9
##  [865] 6.2 3.2 5.5 4.5 4.4 4.1 4.8 3.7 3.9 4.6 4.0 4.6 4.7 4.9 4.0 3.2 3.0 4.7
##  [883] 3.1 4.6 4.1 5.0 4.9 5.4 5.9 5.0 4.9 5.3 4.6 3.6 4.2 4.6 3.8 1.9 4.2 5.1
##  [901] 4.9 4.9 3.9 3.7 4.1 2.9 4.7 4.6 4.4 4.5 4.4 5.4 3.4 5.7 4.2 4.6 4.7 3.1
##  [919] 5.9 5.5 5.4 5.7 3.5 2.7 3.5 6.5 3.8 4.2 4.6 4.1 4.5 4.3 6.5 5.6 5.6 5.1
##  [937] 5.4 4.7 2.7 6.0 4.8 3.8 4.9 5.0 3.3 4.2 4.9 5.6 4.1 6.3 3.4 6.8 4.4 3.9
##  [955] 6.3 3.3 5.0 4.5 4.6 4.7 5.3 2.4 2.9 3.8 6.9 5.2 4.7 4.0 5.4 5.6 5.7 4.8
##  [973] 5.7 4.1 5.1 4.9 5.0 5.6 4.9 3.1 4.5 6.0 4.4 2.9 4.2 5.0 4.3 4.6 3.3 5.3
##  [991] 4.9 3.3 4.5 4.0 6.2 5.2 3.0 5.1 3.3 4.6
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
##    [1] 3.6 3.9 4.4 5.8 4.2 3.9 3.8 3.1 5.2 4.5 4.2 5.7 4.6 4.9 4.4 4.7 4.4 4.7
##   [19] 4.0 3.7 4.2 3.1 3.4 5.6 5.6 3.7 5.1 4.0 4.0 4.6 4.1 3.7 4.2 5.3 3.8 3.8
##   [37] 4.6 4.5 6.8 3.3 5.2 3.8 6.2 3.9 5.8 5.2 3.9 3.9 6.0 5.4 4.1 5.2 3.4 4.5
##   [55] 5.2 3.8 4.4 4.6 6.4 4.6 7.0 4.4 4.1 3.0 6.5 5.6 6.1 4.4 4.0 4.2 4.4 4.3
##   [73] 4.3 5.6 5.0 3.6 4.1 5.1 2.4 5.5 5.6 4.2 4.1 4.6 2.7 3.5 3.8 5.1 5.5 2.7
##   [91] 3.8 6.3 3.8 3.9 4.5 6.2 4.4 5.4 3.0 6.2 4.5 4.5 3.9 3.9 4.8 4.7 5.2 3.0
##  [109] 6.0 3.2 3.6 4.6 5.4 5.0 4.3 4.4 5.5 4.9 4.7 5.6 5.6 4.3 5.0 3.8 4.5 5.0
##  [127] 4.5 4.1 5.2 4.1 4.6 4.0 5.5 4.9 4.1 4.1 4.0 4.6 2.6 3.8 5.4 5.8 4.4 3.2
##  [145] 4.4 5.9 2.7 4.1 5.1 5.4 5.2 4.5 5.6 4.9 3.8 4.1 4.3 3.2 3.8 3.2 5.1 5.1
##  [163] 3.4 4.6 5.7 3.5 3.6 3.7 4.8 4.5 3.6 5.0 4.8 4.0 4.9 5.3 5.0 5.1 5.5 3.0
##  [181] 3.4 6.4 6.3 4.0 4.0 3.8 4.1 2.9 5.9 3.5 2.8 4.7 4.3 6.2 4.5 4.9 4.6 5.3
##  [199] 4.7 3.0 5.0 4.3 4.4 4.2 4.6 2.6 3.0 4.9 4.0 2.0 4.1 5.0 4.2 2.3 4.3 5.4
##  [217] 3.5 6.7 5.4 4.4 3.7 3.2 4.9 3.5 5.2 4.2 3.3 4.4 4.2 5.0 4.7 4.7 4.7 6.4
##  [235] 5.1 4.8 3.5 5.5 4.6 4.6 5.7 4.3 3.2 4.1 3.6 4.2 4.6 6.0 4.2 3.7 3.6 2.0
##  [253] 6.2 5.7 5.8 5.4 5.0 2.3 5.1 4.4 2.9 3.8 4.7 4.2 4.0 3.1 5.8 5.9 3.8 4.2
##  [271] 5.2 4.8 5.0 4.8 4.4 4.4 4.2 4.3 5.2 4.6 3.9 4.9 4.5 5.0 3.8 4.2 5.9 5.0
##  [289] 4.4 3.2 3.8 4.5 5.1 5.4 5.5 4.7 3.0 5.3 5.3 4.0 4.6 4.9 4.8 3.9 5.4 4.1
##  [307] 4.1 6.1 5.1 4.6 5.2 4.2 4.2 5.6 4.6 4.5 5.4 5.8 5.1 3.8 5.3 4.9 5.6 3.9
##  [325] 5.8 3.6 4.1 3.1 4.7 6.2 4.2 3.1 5.8 3.4 5.1 4.1 4.1 4.9 4.5 3.8 3.4 4.9
##  [343] 3.0 4.7 4.9 3.6 3.6 4.8 4.7 4.3 3.7 3.4 5.2 3.9 4.1 5.3 3.5 5.0 3.7 4.6
##  [361] 4.7 4.9 5.7 5.9 4.5 3.8 4.5 5.5 4.6 4.1 2.5 6.0 5.0 3.6 3.1 3.2 5.0 5.3
##  [379] 3.9 4.4 5.3 4.1 4.8 3.7 5.3 3.5 3.1 3.9 4.7 4.1 4.1 3.2 5.1 5.5 4.2 4.4
##  [397] 4.8 5.2 5.4 5.8 4.7 5.0 3.8 3.9 4.2 4.4 5.2 5.7 5.8 4.2 5.1 5.3 4.6 5.8
##  [415] 4.3 5.7 5.0 4.2 6.3 5.2 4.0 4.9 6.3 6.7 6.8 4.8 6.1 5.4 3.7 5.1 4.1 4.4
##  [433] 4.9 6.2 5.5 7.0 5.9 5.7 4.6 4.7 4.8 5.6 5.2 2.5 4.7 5.0 4.9 5.5 3.8 4.4
##  [451] 5.3 3.9 4.6 4.4 5.9 4.1 5.0 4.5 4.5 5.2 4.6 5.0 5.5 5.1 3.8 2.7 2.6 4.2
##  [469] 4.7 3.0 5.8 4.9 3.9 4.8 5.3 4.6 3.8 3.0 3.9 3.8 4.5 5.6 5.5 4.7 4.2 4.0
##  [487] 5.3 3.4 2.9 4.0 4.1 4.4 4.9 5.0 3.1 5.0 5.6 4.0 3.3 3.9 5.0 4.0 4.5 4.7
##  [505] 4.8 3.3 4.4 4.1 5.4 4.8 5.1 4.4 4.4 3.9 5.0 4.5 4.6 5.0 4.2 4.9 4.6 4.4
##  [523] 4.0 4.0 4.0 6.3 5.5 4.7 5.5 3.6 4.6 5.9 5.8 4.2 2.9 5.9 5.8 4.1 3.3 4.9
##  [541] 3.5 5.1 5.7 4.0 3.8 5.1 4.5 4.4 2.8 4.5 5.6 3.5 6.2 4.6 4.4 4.9 3.5 5.6
##  [559] 5.5 2.8 6.0 3.1 3.7 4.8 4.6 3.2 5.1 4.9 5.3 4.0 4.4 4.8 4.8 4.5 4.6 4.9
##  [577] 4.8 4.3 3.9 3.6 5.1 3.1 3.4 4.1 4.4 5.1 5.0 4.6 2.6 3.6 3.0 4.7 5.9 2.8
##  [595] 4.4 4.8 5.6 5.8 5.3 3.0 4.8 5.6 3.7 4.6 4.5 6.7 4.5 5.1 6.6 4.5 4.7 4.3
##  [613] 4.9 5.6 4.3 4.1 3.8 5.4 4.8 3.3 5.4 5.7 4.6 3.6 5.8 4.3 4.6 4.8 4.1 5.0
##  [631] 5.2 4.8 4.6 4.3 4.3 5.4 3.9 5.1 5.1 4.1 6.1 5.1 5.2 3.2 4.6 5.1 5.3 3.8
##  [649] 4.5 3.5 6.7 4.2 4.5 3.7 5.8 6.1 3.8 3.5 6.1 4.0 5.6 4.9 4.4 3.7 4.2 4.4
##  [667] 4.2 5.4 5.6 4.5 4.9 4.0 4.9 5.5 5.1 6.2 5.0 4.6 5.4 5.6 4.1 5.9 5.5 2.8
##  [685] 4.4 4.5 4.7 5.8 4.3 4.2 2.9 4.5 4.5 4.2 4.3 6.2 5.1 6.5 4.1 4.7 3.4 2.6
##  [703] 5.9 4.6 3.4 4.7 5.5 3.6 3.3 4.9 6.7 3.1 3.5 5.0 3.8 4.0 3.5 4.0 4.2 4.3
##  [721] 2.9 3.3 4.2 3.8 5.5 5.7 5.5 4.1 4.7 3.2 3.5 4.4 4.4 5.6 5.6 4.3 5.9 4.9
##  [739] 3.9 1.8 5.1 4.4 4.5 4.1 6.2 4.8 4.1 4.1 4.7 6.6 4.3 3.2 4.1 4.3 4.4 3.2
##  [757] 4.1 4.2 3.1 3.1 4.9 3.5 5.1 4.8 5.0 4.8 4.0 5.0 4.7 4.8 4.4 4.5 7.1 4.0
##  [775] 5.6 4.6 6.3 5.2 3.2 3.7 5.4 3.4 3.9 6.6 5.7 3.8 4.6 4.3 4.2 4.9 4.4 4.9
##  [793] 5.0 4.7 5.7 6.1 3.6 4.2 4.0 5.2 3.9 4.8 3.3 5.3 5.2 4.9 3.7 2.9 4.5 3.5
##  [811] 3.6 5.1 4.4 5.4 4.1 6.4 3.7 4.2 5.0 5.0 5.4 3.9 3.8 5.9 3.9 3.8 3.5 6.0
##  [829] 5.1 3.4 3.0 6.0 5.2 2.3 4.3 3.1 5.1 4.5 3.8 2.7 5.0 2.1 3.1 5.3 5.6 3.9
##  [847] 6.3 4.5 3.3 5.4 4.8 4.7 4.9 5.1 4.6 4.1 4.1 5.1 5.0 3.1 5.3 4.6 5.0 5.0
##  [865] 5.6 5.1 4.7 3.7 5.3 4.0 4.3 4.7 5.2 2.5 3.5 2.5 5.8 4.2 3.8 4.3 5.4 3.7
##  [883] 3.7 3.6 4.6 4.3 6.0 3.5 4.4 4.8 4.9 4.3 4.9 4.2 4.9 3.8 4.7 4.8 4.5 4.0
##  [901] 5.0 3.8 5.7 6.1 3.6 4.6 4.9 1.9 2.6 5.5 3.8 5.6 4.8 4.9 4.6 3.7 4.5 5.0
##  [919] 4.2 4.2 3.2 4.5 4.2 5.7 4.6 3.7 6.6 4.7 3.9 4.3 5.2 4.3 4.8 5.3 3.7 6.0
##  [937] 5.2 4.6 3.8 5.7 4.0 4.7 2.6 4.4 4.7 4.1 4.6 5.4 5.8 5.2 4.3 5.5 4.9 3.2
##  [955] 6.0 4.3 5.6 4.0 5.5 5.4 3.6 2.8 6.5 4.5 5.4 4.2 5.4 4.6 5.0 4.7 4.1 5.6
##  [973] 2.9 4.6 6.1 5.3 3.9 6.2 3.1 4.6 3.2 4.0 5.9 4.3 4.2 4.1 4.0 4.1 4.8 4.7
##  [991] 4.7 5.9 4.9 4.9 3.2 3.9 3.0 2.4 4.1 4.1
## 
## $func.thetastar
## [1] 0.0398
## 
## $jack.boot.val
##  [1]  0.51505682  0.44782609  0.28848168  0.20830946  0.15482955  0.01642651
##  [7] -0.15304348 -0.22251462 -0.34956772 -0.46028986
## 
## $jack.boot.se
## [1] 0.9505027
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
##    [1] 5.3 5.4 4.7 4.0 6.6 5.1 4.0 4.4 5.0 3.3 4.2 4.4 4.3 3.8 4.2 3.4 4.2 5.0
##   [19] 3.4 3.0 4.4 4.4 4.3 6.4 5.5 4.8 4.1 5.8 5.2 3.6 4.3 4.1 3.5 5.2 4.2 4.7
##   [37] 4.1 5.3 4.6 4.3 2.7 3.5 4.7 3.5 4.4 4.3 4.1 4.6 5.2 3.9 5.0 4.6 5.4 3.0
##   [55] 5.6 5.0 5.7 3.8 5.2 5.1 5.0 3.9 5.6 4.7 4.0 4.5 3.0 3.2 4.7 3.4 4.6 6.1
##   [73] 5.0 3.1 4.7 4.0 4.9 5.2 4.7 4.0 4.5 4.3 4.7 4.4 3.8 4.8 5.6 4.4 4.3 3.9
##   [91] 5.3 4.7 3.8 3.7 2.9 2.5 3.8 4.5 3.8 2.6 2.8 4.2 1.7 5.2 4.6 3.6 4.0 4.5
##  [109] 4.6 5.6 3.9 2.2 5.7 5.1 4.3 6.1 4.5 6.0 4.4 3.8 2.9 5.0 5.4 5.3 4.3 4.3
##  [127] 5.0 2.9 5.2 3.9 3.5 3.8 3.3 4.0 5.3 4.5 4.8 5.7 4.8 5.2 4.5 3.6 5.2 3.8
##  [145] 4.7 5.6 2.6 4.7 3.0 2.9 4.5 3.9 5.1 3.9 4.8 3.2 3.8 3.9 5.0 5.3 4.6 4.3
##  [163] 4.1 3.8 4.2 5.8 3.7 4.0 4.0 3.4 5.3 4.0 4.7 2.9 4.9 6.3 4.5 4.9 5.9 5.7
##  [181] 5.6 4.8 3.7 3.3 5.5 2.7 4.7 2.7 6.2 3.4 5.5 4.0 6.0 3.5 3.6 4.4 5.2 4.1
##  [199] 5.2 2.9 6.3 5.3 5.4 2.7 5.4 3.5 4.2 6.8 3.4 5.0 5.7 4.4 3.7 4.3 2.8 4.7
##  [217] 6.0 5.2 5.7 4.6 4.9 4.9 3.8 2.9 4.2 4.0 4.8 5.9 4.2 4.3 4.3 5.2 4.8 3.9
##  [235] 5.5 5.5 5.0 5.0 4.0 3.3 4.7 4.5 3.8 4.5 4.7 5.1 4.9 2.5 4.0 4.0 3.8 5.1
##  [253] 4.8 3.4 2.9 3.3 4.4 4.4 3.6 4.4 2.9 6.2 4.1 3.9 5.8 5.3 4.7 5.1 4.1 4.1
##  [271] 4.9 5.9 4.7 3.1 4.4 5.8 3.9 3.2 5.3 4.2 3.2 3.4 3.9 6.0 3.6 5.1 4.6 4.7
##  [289] 4.7 4.3 2.1 5.4 2.8 6.2 4.9 3.8 3.8 4.4 3.8 4.2 4.7 3.6 3.8 5.9 3.9 5.9
##  [307] 5.0 4.1 4.3 3.5 3.3 5.6 5.2 5.3 5.4 4.6 3.5 3.7 5.9 1.7 3.4 3.0 3.7 4.6
##  [325] 4.2 3.9 4.7 4.0 4.4 4.0 5.5 3.0 5.3 5.4 3.9 5.6 6.0 3.6 4.7 4.9 5.1 3.3
##  [343] 5.5 5.1 4.2 4.8 4.6 6.0 4.7 5.1 6.0 5.3 4.4 5.1 4.5 3.2 5.6 5.0 4.2 5.1
##  [361] 5.3 3.5 6.1 4.9 2.2 5.0 5.4 4.1 3.5 5.5 5.9 4.8 4.3 5.2 4.9 4.3 3.6 5.5
##  [379] 5.3 5.6 4.3 4.7 4.6 4.7 5.5 4.1 5.2 5.0 5.0 3.7 4.7 5.7 3.7 4.4 3.9 4.8
##  [397] 5.0 2.9 4.3 3.8 3.1 3.1 4.0 4.4 4.2 4.5 5.8 4.4 4.3 4.4 3.1 3.4 5.2 4.3
##  [415] 4.6 4.4 4.4 4.4 4.1 5.0 3.2 4.5 6.5 6.4 4.2 5.4 4.5 3.2 3.9 5.8 4.3 4.6
##  [433] 4.6 3.3 3.4 6.0 4.5 4.6 3.7 4.7 3.2 4.5 4.2 5.7 4.9 5.1 5.0 4.9 5.5 5.4
##  [451] 4.0 4.4 5.3 5.5 6.0 4.8 4.8 2.6 3.8 4.0 4.6 5.2 4.3 4.8 4.5 3.7 5.1 6.1
##  [469] 4.7 4.0 2.9 4.0 3.4 4.0 4.8 4.1 5.0 4.4 4.7 6.0 5.6 4.3 5.4 5.2 2.9 2.9
##  [487] 4.8 5.2 5.5 5.4 5.5 4.0 4.7 4.0 5.4 5.1 5.4 4.9 5.5 3.4 3.4 4.5 4.6 4.6
##  [505] 5.3 4.4 3.5 4.6 6.4 3.6 5.0 5.0 2.4 3.4 2.5 2.7 3.8 3.7 6.2 4.9 4.4 3.3
##  [523] 4.8 3.5 6.2 4.7 4.6 4.6 4.3 3.6 5.1 5.3 3.3 4.9 5.6 4.1 4.2 3.5 4.2 3.9
##  [541] 2.8 4.8 3.0 3.1 4.8 4.4 4.3 3.8 3.2 3.3 3.6 4.0 4.4 4.0 4.7 2.8 4.1 5.2
##  [559] 4.3 5.8 3.6 3.5 4.2 3.4 4.4 5.4 2.9 4.0 4.2 4.2 5.1 6.4 3.4 5.1 4.5 6.2
##  [577] 4.0 5.2 4.9 4.6 5.4 3.5 5.1 3.9 4.5 4.4 5.1 5.5 4.2 4.1 5.8 3.7 5.4 6.1
##  [595] 3.7 5.0 5.1 4.4 4.4 4.0 5.3 3.3 4.7 3.9 6.3 4.7 3.8 4.6 5.0 4.5 4.8 6.4
##  [613] 2.9 5.0 4.3 3.5 3.4 4.0 3.9 3.7 3.7 3.8 5.6 2.5 4.7 5.7 3.0 4.2 3.6 4.1
##  [631] 4.7 2.4 4.8 5.6 4.4 4.8 4.4 4.4 3.5 4.0 4.2 5.8 4.2 4.8 5.0 4.8 5.5 4.9
##  [649] 5.0 4.0 4.4 4.9 5.2 5.7 3.6 4.4 5.0 3.7 5.0 3.6 5.2 4.6 4.4 5.8 4.5 5.6
##  [667] 5.6 5.4 4.9 3.4 5.7 3.7 6.5 5.2 5.5 5.1 3.6 4.5 4.6 3.7 4.8 4.2 4.7 4.3
##  [685] 3.5 4.7 6.0 4.5 4.7 4.3 5.6 4.3 5.5 4.3 5.0 4.2 5.5 5.3 3.2 5.7 4.5 5.3
##  [703] 5.3 5.2 4.8 5.1 3.3 4.1 4.8 4.1 5.3 4.9 4.1 3.8 5.0 3.3 4.8 3.8 4.8 4.2
##  [721] 6.0 6.0 4.7 4.7 5.3 4.4 5.8 4.5 5.8 5.0 5.2 3.7 3.8 4.7 4.6 4.0 4.6 4.3
##  [739] 4.8 5.2 5.8 3.9 4.9 3.5 4.4 4.8 4.7 4.3 6.0 4.5 4.1 4.8 4.6 4.3 5.3 2.6
##  [757] 4.2 4.8 5.5 4.1 4.0 6.2 5.2 6.2 4.1 3.9 2.5 5.3 5.0 1.6 4.3 5.2 6.7 4.6
##  [775] 5.3 5.1 4.6 4.4 4.9 3.9 3.0 4.4 4.7 5.2 2.5 4.6 3.9 4.8 4.5 3.3 2.8 5.4
##  [793] 3.6 5.1 3.7 4.4 4.8 4.6 4.4 4.3 3.7 4.8 4.8 4.3 3.6 4.8 4.0 4.2 5.5 4.5
##  [811] 4.5 4.3 4.3 2.9 3.4 4.4 3.4 4.6 4.5 3.5 5.0 4.6 4.3 4.8 4.7 4.1 5.3 3.9
##  [829] 3.3 3.4 4.7 5.3 4.6 5.0 4.8 2.3 6.5 4.1 4.0 3.0 4.2 4.3 2.8 5.0 6.4 4.0
##  [847] 2.9 5.3 5.2 3.6 4.4 4.6 4.0 5.4 5.8 4.7 5.1 5.1 4.3 5.9 5.5 3.9 5.1 5.9
##  [865] 3.8 4.3 3.9 4.8 6.0 4.5 5.3 4.1 5.2 5.3 5.6 4.6 5.1 4.9 4.6 4.2 4.4 4.6
##  [883] 2.5 5.2 4.8 3.1 3.8 6.0 4.9 5.5 4.4 4.7 3.1 3.4 3.6 4.4 4.7 4.9 3.7 4.6
##  [901] 2.8 3.0 5.1 5.5 4.7 3.8 4.2 4.7 5.0 5.6 3.9 3.2 4.4 4.1 5.7 6.2 4.5 3.8
##  [919] 5.3 3.0 5.3 5.4 5.4 4.6 5.5 3.6 2.0 3.4 3.1 5.1 5.0 5.8 5.5 3.2 4.1 4.4
##  [937] 3.9 5.4 4.8 2.8 5.1 5.0 3.4 4.5 4.5 4.8 3.0 4.7 6.0 3.8 3.8 5.2 5.0 5.0
##  [955] 3.8 4.3 5.4 4.2 4.5 5.3 4.4 5.6 4.7 4.3 4.7 3.7 5.3 2.8 4.2 3.8 5.9 4.6
##  [973] 5.3 4.0 4.7 4.7 3.5 4.8 3.9 5.0 3.5 5.0 4.1 4.0 4.5 4.7 5.6 4.9 4.9 4.6
##  [991] 5.5 4.3 5.1 3.6 4.9 4.6 5.8 4.1 5.5 4.4
## 
## $func.thetastar
## 72% 
##   5 
## 
## $jack.boot.val
##  [1] 5.400 5.300 5.300 5.200 5.172 5.000 4.800 4.896 4.700 4.500
## 
## $jack.boot.se
## [1] 0.8449236
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
## [1] -0.06481082
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
##   2.421062   5.511177 
##  (1.016995) (2.571753)
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
## [1] -0.37350062  0.71219171 -0.03768583  1.57256611  1.15164413  0.11535178
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
##    [1]  4.640878e-01  3.182490e-01 -5.002058e-01 -9.609852e-03 -4.220180e-01
##    [6] -4.092379e-01  6.465207e-01  4.878573e-01  4.334140e-01 -5.510888e-01
##   [11] -1.145895e+00  1.381601e+00  3.744314e-02 -3.223750e-01  3.135135e-01
##   [16] -1.152248e+00 -4.526300e-01  5.608016e-01 -3.331772e-01  1.224182e+00
##   [21]  5.705156e-01 -2.551526e-01 -3.850597e-01 -5.491488e-01 -7.796918e-01
##   [26] -1.688287e-01 -9.380490e-02 -5.309945e-01  1.028247e-01 -3.882785e-01
##   [31]  3.982834e-01 -3.370569e-01  1.466065e-04 -4.357596e-02  1.748412e-01
##   [36]  1.188155e-01  6.117144e-02  3.135135e-01  1.104187e+00  1.160638e-01
##   [41] -4.434046e-01  8.384846e-01  5.353499e-01 -3.284730e-02 -4.792607e-01
##   [46]  1.025587e-02 -1.588251e-01  3.350433e-02  2.043081e-03 -3.268163e-01
##   [51] -7.712470e-01 -2.015524e-01 -1.681770e-01 -2.461033e-01  1.027675e-01
##   [56]  3.351908e-01 -2.267298e-01 -3.412345e-01  3.642686e-02  9.425355e-01
##   [61] -6.401789e-01  8.711194e-01 -6.263104e-01 -4.150397e-01  5.780917e-02
##   [66] -2.237231e-01 -4.123714e-01 -3.707063e-02 -2.034692e-02 -2.123288e-01
##   [71] -1.176078e+00  5.980758e-01  3.236976e-01  2.776899e-01  1.068655e-01
##   [76] -4.745179e-01 -9.934147e-01 -2.066513e-01 -4.893352e-01 -5.515406e-01
##   [81] -4.634705e-01 -4.811826e-01 -3.276070e-02  6.375326e-01  4.399922e-01
##   [86] -2.080007e-01 -1.247674e+00  2.812955e-01  2.168367e-02  3.086389e-02
##   [91] -1.997637e-01  2.382097e-01  3.986027e-01  5.554399e-02 -4.504823e-02
##   [96]  7.828908e-01 -8.097089e-01 -2.134596e-01  2.972137e-01  4.225514e-01
##  [101] -2.837136e-01  3.367603e-01  1.190294e-01 -5.415248e-01  6.562051e-02
##  [106] -6.734920e-01 -6.233158e-01 -8.998919e-03 -5.134339e-01 -6.902195e-01
##  [111] -3.299066e-02 -5.578114e-01 -1.080791e+00 -3.178471e-01  8.147429e-03
##  [116]  6.346652e-01 -4.471509e-01 -2.246091e-01  1.062442e+00 -4.584922e-02
##  [121]  2.931437e-01  1.023935e+00 -3.653723e-01  7.970092e-03  9.072364e-02
##  [126] -1.213686e-02 -1.940140e-01 -1.220456e-03 -8.565278e-01  7.354013e-01
##  [131]  6.757281e-01 -5.042653e-01 -5.441450e-01 -7.058138e-01  7.935194e-01
##  [136]  4.901149e-01 -8.236158e-01 -1.860095e-01 -1.277485e-01  3.004512e-01
##  [141]  2.678608e-02 -3.377493e-02 -3.452526e-01  3.097239e-02  5.247562e-02
##  [146] -2.673171e-01  2.789061e-01 -1.819709e+00 -1.027749e+00  1.293018e+00
##  [151] -1.107275e-01  1.105386e-01 -1.024101e+00  2.979104e-01 -4.063334e-02
##  [156] -1.960464e-01  4.395632e-01 -3.146426e-01 -1.223790e+00 -1.693697e-01
##  [161]  2.840479e-01  3.599957e-01 -8.758213e-01  2.249785e-01  9.553592e-01
##  [166] -5.151761e-01 -1.126901e-01 -1.612762e-01 -2.346622e-01  5.762213e-01
##  [171]  2.896019e-01  5.862737e-01  6.970302e-02 -7.398815e-02 -2.543444e-02
##  [176] -4.355104e-01 -4.516984e-02  8.760440e-02  7.994485e-02 -7.232381e-01
##  [181] -4.663557e-01  5.404699e-01 -2.497705e-01 -2.197720e-01 -5.706156e-03
##  [186]  8.556941e-01 -3.835964e-01 -6.922485e-01 -9.007874e-02  3.033781e-01
##  [191] -1.016470e+00 -1.059036e+00 -1.231473e+00 -1.888702e-01  8.928570e-02
##  [196] -3.486839e-01 -6.120188e-01  2.429591e-01 -1.991739e-01  3.979731e-01
##  [201]  8.166257e-02  4.503364e-01 -1.447919e-01 -5.518140e-02 -8.834579e-02
##  [206]  1.802464e+00 -4.357596e-02 -6.671065e-04  1.111692e+00 -2.089472e-01
##  [211] -2.056038e-01 -4.366348e-01  5.085891e-01  7.420964e-01 -2.300121e-01
##  [216] -5.407363e-01  2.358024e-02 -4.059485e-01  1.767946e-01 -7.939997e-01
##  [221] -2.214480e-01 -5.023138e-01 -3.086653e-01 -2.731057e-01 -2.753569e-01
##  [226]  1.985459e-01  4.027106e-01 -3.941701e-01 -8.265762e-02 -4.904133e-01
##  [231] -1.418480e+00  3.387825e-01 -4.495924e-01  3.622882e-01 -1.798878e-01
##  [236]  4.687570e-02  7.882899e-02  5.386066e-01 -7.746427e-01 -3.203900e-02
##  [241] -1.077448e+00 -1.216131e+00 -1.823013e-01 -5.882815e-01  3.950604e-01
##  [246] -5.624809e-01 -5.447824e-01  7.728445e-01 -8.907226e-01 -2.615919e-01
##  [251]  7.256191e-01 -5.327077e-02 -3.108550e-01 -4.689208e-01  2.721972e-01
##  [256] -5.774658e-01 -9.722028e-01 -2.978550e-01  7.445535e-01  4.487766e-02
##  [261]  2.438594e-01 -1.397038e+00  1.010172e+00 -7.187600e-01  5.534346e-01
##  [266] -1.948606e-01 -6.450992e-01  3.488480e-01 -1.422483e-01 -5.638464e-01
##  [271] -5.512726e-01 -3.023946e-01 -4.009398e-01  5.442319e-02  7.598348e-01
##  [276] -1.120292e+00 -5.299379e-02 -4.514676e-01 -1.364509e+00  2.992624e-01
##  [281]  1.615050e-01  2.108284e-01  2.631332e-01  1.877448e-01  6.140902e-01
##  [286] -4.192981e-02 -1.020624e-01  3.560393e-01 -1.821996e-03  1.405585e+00
##  [291] -3.533662e-01  1.260124e-01 -4.948176e-01 -5.483807e-02 -2.234734e-02
##  [296] -7.361249e-02  5.861194e-01 -1.545631e-01 -7.248682e-01 -2.931238e-01
##  [301] -6.674094e-01  5.626428e-01 -3.332600e-01  1.071405e+00  5.137693e-02
##  [306] -5.776397e-01 -2.257807e-01  4.345330e-01  6.152803e-01 -3.269622e-01
##  [311] -5.478645e-02  8.445430e-01 -5.918839e-01  1.970148e-01 -3.675985e-01
##  [316] -9.088063e-02 -7.711243e-01 -2.529832e-01 -6.018456e-01 -1.278710e-02
##  [321]  1.047365e-02 -3.062147e-01 -4.604093e-01  1.992301e+00  3.644772e-01
##  [326]  1.815039e-01  5.529279e-01 -4.176599e-01  1.748412e-01  8.872831e-01
##  [331] -2.286396e-01 -2.065350e-01  9.958306e-01  2.269585e-01 -5.608806e-01
##  [336]  8.205799e-01  9.990032e-02  1.108269e-01 -2.953363e-02 -1.607905e-01
##  [341] -4.509456e-02  2.907617e-01  3.184715e-01 -1.717761e-01  1.675405e-01
##  [346] -3.872688e-01 -8.391654e-01  1.953966e-01  4.792396e-01 -4.993686e-01
##  [351] -3.775329e-01  1.773914e-02 -4.370446e-01  6.892801e-02 -2.170981e-01
##  [356]  4.652071e-01 -2.838143e-01 -1.670617e-01 -6.689412e-01 -7.136467e-01
##  [361] -1.710228e+00  1.129160e-01  2.831139e-01 -6.089252e-01  4.595186e-01
##  [366]  5.462939e-01 -1.591682e-01  2.030095e-01 -2.623234e-01  3.544056e-01
##  [371] -9.467692e-02 -9.658716e-01 -3.732414e-01  8.468440e-01  7.352229e-01
##  [376]  1.962243e-01  1.314969e+00 -3.925465e-01  2.145284e-01 -1.584982e-01
##  [381] -7.167997e-01 -3.736262e-01  4.299791e-01 -1.686403e-01  2.467130e-01
##  [386] -3.874104e-01  8.178306e-01  1.909503e-01 -9.229590e-01  1.997122e-02
##  [391]  8.344130e-01  5.068218e-01 -4.593744e-01 -6.647384e-02 -4.373127e-02
##  [396] -1.001728e-01 -9.426833e-01 -2.768987e-01 -2.482933e-01 -6.682370e-02
##  [401]  1.135186e-01 -1.134742e+00  2.302856e-01  4.228943e-01  7.866926e-01
##  [406] -1.785666e-01  8.085932e-01  1.558333e+00 -3.613990e-01 -3.169900e-02
##  [411] -5.023138e-01  5.847042e-01 -6.696165e-01  5.256299e-01  3.184320e-01
##  [416] -3.309645e-01  4.216796e-01 -2.617077e-01  7.335834e-01 -8.339387e-01
##  [421] -9.970251e-01 -1.563598e+00 -8.056613e-01  2.909565e-01  5.404765e-01
##  [426]  3.224540e-01 -7.203127e-03 -1.007004e-01 -6.033510e-01 -3.822154e-02
##  [431] -1.298126e+00  5.863163e-02 -4.884015e-01 -6.407302e-01 -1.046720e-01
##  [436] -5.604714e-01  1.988180e-01 -9.186087e-01 -5.871625e-01 -7.727239e-02
##  [441]  5.855036e-03 -2.526489e-01  1.175104e-01 -1.588444e-01 -1.797968e+00
##  [446]  6.789454e-01 -1.257036e-01  5.982510e-01  1.132920e+00  2.727447e-01
##  [451]  2.309701e-01  1.698946e-02 -1.769222e+00 -1.048385e+00  2.279859e-02
##  [456]  8.047240e-01  3.453635e-01 -1.096298e-01 -6.625246e-01  7.687671e-01
##  [461] -7.381862e-01  2.569648e-01 -4.574057e-01  2.675408e-01  3.451180e-01
##  [466]  3.863761e-01  4.980022e-01  6.577477e-01  1.524788e+00  1.224587e-01
##  [471] -2.506259e-02  4.099792e-02 -1.705325e-01 -2.126852e-03 -4.019857e-01
##  [476] -4.765618e-01 -3.042213e-02 -3.053704e-01 -1.568906e-01 -2.717850e-01
##  [481]  2.730030e-01 -6.703779e-01  6.024901e-01 -4.334174e-01 -3.173527e-01
##  [486] -7.559921e-01 -1.856344e-01 -3.792012e-01  1.306335e-01 -8.595036e-01
##  [491] -1.238513e+00 -3.055705e-01  3.464776e-02 -1.420195e+00 -5.238027e-01
##  [496] -2.067437e-01 -6.905315e-01 -4.912815e-01 -4.474852e-01  2.342748e-01
##  [501]  2.141122e-01  2.304973e-01 -1.074668e+00 -2.141779e-01 -2.926389e-02
##  [506]  1.967095e-01 -7.253774e-01  1.123741e-01 -1.437793e+00 -7.193370e-01
##  [511] -4.483287e-01 -6.687720e-01  5.287495e-02  4.331474e-01 -5.762878e-01
##  [516] -2.919424e-01 -9.419995e-01  9.820838e-02 -3.239915e-03 -1.731787e-01
##  [521] -2.192590e-01 -3.991649e-01  8.609891e-01 -9.808393e-01  2.007944e-02
##  [526]  3.897080e-01  5.616901e-02 -1.577607e-01 -3.281638e-01  6.270852e-01
##  [531]  5.020033e-01 -7.104114e-01  3.262845e-01 -9.277141e-01 -1.505685e-01
##  [536]  1.303201e-01  3.467471e-01 -1.351453e-01  4.749452e-01 -6.507967e-01
##  [541] -3.201615e-01 -1.748059e-01  2.449597e-01 -1.237471e-01  8.845072e-02
##  [546]  5.114686e-02 -1.848049e+00  6.042403e-01 -4.116854e-01 -4.083560e-01
##  [551]  4.603215e-01 -2.958050e-01 -4.613712e-01 -1.402890e-01  6.113666e-02
##  [556] -6.698909e-01  5.986714e-01  1.232377e+00 -1.239072e+00  5.955363e-01
##  [561] -1.268369e+00 -3.206183e-01 -1.049310e+00  6.557427e-02  2.484024e-01
##  [566]  4.337962e-01  1.077502e+00 -8.534795e-01  7.292215e-01 -4.316801e-01
##  [571] -2.654270e-01  8.051377e-01 -1.484666e-01 -5.227654e-01  4.966994e-01
##  [576]  9.918475e-02 -8.642478e-02 -6.783304e-02 -6.696689e-01  1.455467e-01
##  [581] -4.582645e-01 -1.614089e+00  4.159356e-01 -4.015814e-01  1.560867e-01
##  [586]  5.015627e-01  9.423475e-01 -4.370446e-01 -3.087190e-01  7.808067e-01
##  [591] -7.855686e-01  1.789036e-01 -6.228969e-01  3.108688e-02 -4.086977e-01
##  [596] -8.439226e-01  9.210831e-01 -2.208984e-02  9.022761e-01 -4.022454e-01
##  [601]  3.900064e-01 -1.492659e+00  6.436372e-01  7.566828e-01  6.130068e-01
##  [606]  2.687472e-01 -1.858573e-01 -4.957129e-01 -1.866646e-01  4.657121e-01
##  [611]  1.051516e-01 -1.161707e+00  8.195682e-02 -8.108031e-01  1.964549e-01
##  [616]  7.087393e-02 -5.049068e-01  2.812569e-01 -9.251330e-03 -6.874466e-02
##  [621] -7.087462e-01 -4.183830e-01 -1.557386e+00  5.355381e-01 -2.213891e-01
##  [626] -9.544312e-01  3.219915e-01  1.478666e-01  2.609626e-02  2.301013e-01
##  [631] -4.261517e-02  9.033735e-01 -9.120966e-01  3.627820e-01 -5.061117e-02
##  [636]  3.476334e-01 -2.126282e-01 -9.436893e-01  6.506661e-01  4.237890e-01
##  [641]  1.690927e+00  2.874946e-01 -1.380304e-01 -9.267125e-01  7.590918e-01
##  [646] -1.187834e+00 -7.045652e-01 -1.227621e-01 -9.205293e-01  4.526039e-01
##  [651] -4.892253e-01  1.759302e-01  1.569994e-01 -7.351660e-01 -1.214890e+00
##  [656] -5.379546e-01 -5.050727e-02  9.971647e-03 -6.794036e-01 -2.268584e-02
##  [661] -7.390111e-02 -1.027414e+00  2.684648e-01 -6.800058e-01  2.502621e-01
##  [666] -6.855689e-01  4.947990e-01  4.077365e-02 -1.463578e+00  1.389027e+00
##  [671] -7.769315e-01 -2.494422e-01 -4.244592e-01 -7.184167e-01 -4.133164e-01
##  [676] -3.289267e-01  4.167427e-01 -2.005830e-01 -5.014523e-01 -6.197954e-01
##  [681]  8.791071e-02  1.222809e-01  3.365760e-01 -8.110778e-01  7.169156e-02
##  [686] -1.681612e-01 -2.395632e-01 -4.442023e-01 -6.414012e-02  8.415879e-01
##  [691] -1.613407e-01  1.388934e-01 -4.218597e-01 -1.888702e-01 -8.451587e-02
##  [696] -2.134220e-01 -5.295997e-01  6.702669e-01 -2.485546e-01 -1.638978e-01
##  [701]  7.594897e-01 -2.110296e+00  5.583085e-02 -1.154980e+00  3.646464e-01
##  [706] -3.552903e-01 -1.011106e+00 -8.890926e-01 -1.624577e-01 -4.479062e-01
##  [711] -1.340749e+00 -6.474964e-01  5.226805e-02 -4.748205e-02  2.842540e-01
##  [716]  4.015799e-03 -1.455979e-01 -1.192344e-01 -7.392503e-01  1.407659e+00
##  [721] -5.589990e-01 -2.004749e-01 -4.088850e-01 -1.063565e-01 -4.924073e-01
##  [726]  3.496129e-01  5.857610e-01  9.996730e-01 -1.043272e+00 -3.717310e-01
##  [731]  7.944993e-01 -1.063242e+00  4.567316e-01  1.536868e+00 -3.461696e-01
##  [736]  1.572811e-01 -3.233983e-01  4.219376e-01  4.855284e-02 -1.241520e+00
##  [741] -1.194607e+00 -6.281025e-01  5.072262e-01  5.194379e-01  7.583751e-01
##  [746] -1.057762e+00 -1.373994e+00 -1.502117e-01  7.702257e-01  2.070393e-01
##  [751]  2.324157e-02 -8.307432e-02  3.040280e-01  2.255552e-01 -6.062713e-01
##  [756] -6.813273e-01  7.935766e-02  4.440445e-01 -5.569647e-01  3.343452e-01
##  [761] -7.973108e-01 -1.537762e-01  5.091137e-01 -3.672801e-01 -8.060986e-01
##  [766] -7.095434e-02 -1.578923e+00 -2.912264e-01  1.072760e+00 -4.758938e-01
##  [771] -1.026928e+00 -4.164317e-01  9.484729e-05  7.262141e-02  2.311871e-02
##  [776] -1.686403e-01 -6.712535e-01 -3.830013e-01 -1.140438e-01 -2.178735e-01
##  [781] -9.407543e-01  3.507579e-01  1.106557e+00 -4.052975e-01  3.125111e-01
##  [786] -2.360755e-01 -7.045338e-02  2.202870e-01 -1.711680e-01  3.952809e-01
##  [791] -3.455211e-01  2.862197e-01 -1.725068e-01  1.038425e+00 -2.112565e-01
##  [796]  3.170526e-01  3.737002e-01 -7.411451e-01 -6.136295e-01  3.190990e-01
##  [801] -1.997699e+00  2.327276e-01 -5.432943e-01  8.263324e-02 -8.665061e-01
##  [806] -9.341239e-02 -3.153502e-01  3.824832e-01 -2.921452e-01 -4.480237e-01
##  [811]  8.344340e-02  6.112248e-01  1.973141e-01  2.402976e-01  3.250787e-02
##  [816] -3.981446e-01 -3.976195e-01  1.083028e+00 -1.183751e-01 -7.649558e-01
##  [821]  4.353149e-01 -8.095162e-01  3.421275e-01 -5.892012e-02 -3.325073e-01
##  [826] -1.433913e-01  5.578916e-01 -7.241657e-01 -3.541911e-01  3.852138e-01
##  [831]  4.718431e-01 -2.645731e-01  2.370627e-01 -1.717761e-01  3.437059e-02
##  [836] -7.206647e-01 -6.220703e-01  1.024945e+00  1.442294e-02 -5.266191e-02
##  [841]  2.944056e-01  2.749598e-01 -1.187591e-01 -3.293324e-01  3.322343e-01
##  [846]  2.558509e-01 -1.238944e+00 -3.435527e-01 -5.756723e-01  1.410307e-01
##  [851] -4.428806e-01 -3.057823e-01  4.365542e-01  6.164678e-01 -6.831476e-01
##  [856]  8.357114e-01  4.007070e-01 -9.374254e-01  6.387065e-01 -5.440246e-01
##  [861] -1.333380e-01  2.950753e-01  6.916226e-01  5.768064e-02  4.760370e-01
##  [866]  1.859129e-01  1.046318e+00 -4.336245e-01  8.595774e-01  8.284144e-01
##  [871]  8.617085e-01  3.555179e-01 -4.244319e-01  2.059419e-01 -9.264250e-02
##  [876] -2.375370e-01 -1.445741e-01  2.134749e-01 -6.367354e-01 -3.786432e-01
##  [881] -6.152517e-01  2.170088e-01  9.216047e-01  6.453577e-01 -1.067983e-01
##  [886]  1.477228e-01 -2.646926e-01  7.681896e-01  4.943801e-01 -2.234734e-02
##  [891] -1.501031e-02 -3.233804e-01 -1.031155e+00  3.631236e-01  6.898054e-01
##  [896] -6.238022e-01 -1.502117e-01  5.699448e-01 -6.644205e-01  6.088528e-01
##  [901]  4.064229e-01 -3.360462e-01 -1.231011e-01 -1.125899e-01 -6.686075e-01
##  [906] -1.542626e-01 -6.492639e-02  1.027200e+00  2.498614e-01 -8.210500e-01
##  [911]  4.940023e-02 -1.014444e+00  1.867152e-01  3.746956e-01  1.194173e+00
##  [916]  1.344563e+00  1.104270e-01  2.635177e-01 -1.713163e-01 -6.511029e-01
##  [921]  4.937011e-01  1.888135e-01  4.011362e-01  2.472085e-01  6.373540e-01
##  [926] -4.261089e-01  1.410307e-01  2.762380e-01 -4.128675e-01 -4.097715e-01
##  [931]  6.353218e-01  3.254275e-01  1.525387e+00 -2.177246e-01  7.544479e-01
##  [936] -2.266198e-01 -1.153338e-01  4.277269e-01  2.991700e-01  4.774522e-01
##  [941]  6.234260e-01  3.619296e-01 -7.142975e-02 -1.420589e-01 -7.654752e-01
##  [946]  7.648565e-02 -2.797082e-01  4.412216e-01 -8.484922e-01 -8.152918e-01
##  [951]  4.901149e-01 -6.146654e-01 -5.138843e-02  3.353749e-01 -5.826030e-01
##  [956] -1.252863e-01 -1.611104e+00  2.596112e-01 -4.562557e-01  6.827901e-01
##  [961] -1.268841e-01  2.104820e-01 -2.112499e-01  1.292291e+00 -5.585912e-01
##  [966] -3.652506e-01  5.396396e-01 -1.259813e+00 -1.201291e+00 -2.668273e-01
##  [971]  2.291519e-01 -1.287407e-01 -6.954538e-01  2.897021e-01 -4.764891e-01
##  [976] -1.241787e+00 -1.637009e-01 -7.885134e-01  2.758455e-01 -1.005305e+00
##  [981] -2.628110e-01 -4.799853e-01 -5.718778e-02 -7.933295e-01  8.248954e-01
##  [986]  2.205677e-01 -5.577850e-03 -3.267816e-01  3.956176e-01 -1.006200e+00
##  [991] -7.266085e-01  1.056450e-01  5.415774e-01  1.415252e+00 -7.850433e-02
##  [996] -2.301065e-01  4.192339e-01 -1.040880e-01  1.556803e-01 -1.104232e+00
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
##   0.43930000   0.24717059 
##  (0.07816220) (0.05526607)
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
## [1]  0.67042330 -0.02995041  0.88131538 -0.73048628  0.50582837 -0.03637209
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
## [1] -0.0243
```

```r
sd(results2$thetastar)  #the standard deviation of the theta stars is the SE of the statistic (in this case, the mean)
```

```
## [1] 0.9306264
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
## t1*      4.5 0.05285285      0.9327
```

One of the main advantages to the 'boot' package over the 'bootstrap' package is the nicer formatting of the output.

Going back to our original code, lets see how we could reproduce all of these numbers:


```r
table(sample(x,size=length(x),replace=T))
```

```
## 
## 1 3 4 5 6 7 9 
## 2 1 1 2 1 2 1
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
## [1] -0.015
```

```r
se.boot
```

```
## [1] 0.9075233
```

Why do our numbers not agree exactly with those of the boot package? This is because our estimates of bias and standard error are just estimates, and they carry with them their own uncertainties. That is one of the reasons we might bother doing jackknife-after-bootstrap.

The 'boot' package has a LOT of functionality. If we have time, we will come back to some of these more complex functions later in the semester as we cover topics like regression and glm.

