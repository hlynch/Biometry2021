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
## 1 2 3 5 6 8 
## 2 1 4 1 1 1
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
## [1] 2e-04
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
## [1] 2.725496
```

```r
UL.boot
```

```
## [1] 6.274904
```

Method #2: Simply take the quantiles of the bootstrap statistics


```r
quantile(xmeans,c(0.025,0.975))
```

```
##   2.5%  97.5% 
## 2.6975 6.3000
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
##    [1] 3.8 4.9 5.1 5.9 4.8 5.4 4.7 4.4 5.5 4.5 4.6 3.7 3.3 5.4 3.6 6.1 5.5 4.4
##   [19] 2.5 5.9 3.7 6.1 4.5 3.9 4.0 2.8 4.1 3.2 3.6 4.4 5.9 5.1 4.4 4.9 3.5 4.9
##   [37] 4.8 2.9 4.0 4.3 5.0 5.0 4.8 4.9 5.0 2.3 2.9 3.0 4.2 2.6 5.6 4.8 5.5 4.2
##   [55] 3.8 3.8 6.0 4.8 3.4 4.8 5.2 2.8 6.0 4.8 2.5 7.6 4.8 4.5 3.8 4.6 3.4 3.5
##   [73] 4.3 5.0 3.7 4.7 5.9 5.0 4.8 5.6 4.4 4.8 5.0 5.9 3.9 4.1 6.3 4.6 4.0 4.5
##   [91] 4.9 3.7 4.8 3.5 7.0 5.7 4.7 3.7 3.8 3.2 5.4 4.7 4.8 3.9 5.6 4.0 3.6 6.2
##  [109] 4.0 3.9 4.4 5.7 4.3 4.5 4.2 3.5 4.7 3.9 3.5 4.6 6.2 4.4 4.7 3.2 3.0 6.6
##  [127] 5.0 4.2 5.6 5.2 4.6 5.9 4.2 5.3 4.2 4.7 3.3 3.4 7.3 2.9 5.4 3.5 4.5 3.6
##  [145] 2.7 4.5 3.5 3.3 4.6 3.1 6.1 4.7 3.7 4.5 3.2 4.1 4.3 5.8 4.4 4.3 5.8 3.9
##  [163] 4.9 3.7 5.6 3.7 5.5 4.1 3.6 6.2 4.5 4.3 5.1 6.2 5.7 2.8 5.5 3.3 2.7 6.2
##  [181] 4.2 4.8 4.1 5.8 4.2 3.6 4.2 5.8 4.6 3.7 2.5 5.3 3.8 3.6 6.5 4.9 4.3 3.2
##  [199] 3.7 4.9 6.2 3.6 5.3 3.7 4.1 5.4 4.7 4.4 5.2 4.9 5.9 5.4 4.3 3.7 4.9 4.9
##  [217] 5.8 5.7 3.2 5.1 2.8 4.1 4.7 2.2 4.3 4.6 4.5 5.5 3.9 5.0 5.5 3.8 2.2 5.0
##  [235] 3.5 4.7 5.1 5.8 5.2 4.4 6.8 4.2 6.5 4.0 4.9 5.3 4.4 5.1 4.8 4.2 4.4 5.4
##  [253] 4.4 4.0 4.9 3.3 5.0 5.1 4.6 3.6 4.3 4.2 2.7 4.4 4.3 4.1 5.4 4.9 4.5 4.2
##  [271] 4.6 4.5 4.9 4.5 3.1 5.1 4.9 4.9 4.5 5.1 3.8 3.6 3.5 3.5 3.4 5.3 5.8 6.8
##  [289] 5.5 5.3 5.4 4.7 4.6 4.9 4.0 5.1 5.0 4.5 5.6 5.9 4.1 4.3 2.8 4.1 4.8 6.1
##  [307] 4.6 4.2 4.2 5.3 3.7 4.6 4.1 3.9 4.5 6.3 3.9 2.8 4.4 5.7 4.8 5.9 5.0 3.6
##  [325] 4.3 3.2 3.4 5.7 4.5 3.9 2.4 4.2 3.6 4.7 5.2 5.1 3.6 5.5 4.4 4.4 3.2 5.3
##  [343] 4.4 4.0 5.2 4.7 6.1 4.3 5.6 3.7 5.3 5.4 4.9 4.6 4.1 3.3 5.8 3.1 4.9 4.9
##  [361] 4.3 3.0 6.4 5.1 4.8 5.5 6.9 4.2 5.9 4.8 3.2 3.2 3.6 5.1 4.1 6.1 4.4 5.1
##  [379] 1.6 5.1 4.9 4.7 6.4 5.8 4.2 4.0 3.6 6.5 5.0 4.4 5.4 3.9 3.8 4.1 4.5 5.6
##  [397] 4.3 4.5 3.9 5.2 6.0 3.8 5.5 2.2 5.0 3.3 5.0 3.9 4.3 5.6 4.5 5.1 4.6 4.0
##  [415] 2.9 4.4 3.5 4.2 4.0 2.8 3.8 3.7 5.0 5.9 3.6 3.6 3.4 5.6 5.0 4.5 5.9 5.3
##  [433] 5.1 4.0 4.5 4.7 3.9 4.8 3.2 3.5 3.5 5.7 5.0 3.4 4.0 6.6 3.5 4.4 5.0 4.8
##  [451] 4.2 4.4 4.0 3.9 4.5 3.4 5.0 4.0 4.2 3.6 4.4 5.7 5.1 5.3 4.6 4.6 3.6 4.1
##  [469] 3.4 4.3 5.5 6.0 3.2 3.4 3.9 3.4 5.2 3.9 4.4 7.0 4.4 3.0 4.3 3.4 4.3 5.8
##  [487] 4.6 6.2 3.5 4.0 5.0 3.9 3.2 5.3 3.7 5.7 4.5 5.0 5.2 4.4 4.2 6.0 5.4 5.2
##  [505] 3.2 4.0 5.3 6.5 5.3 6.1 3.5 1.6 3.1 5.7 4.3 5.5 4.3 3.5 5.0 4.5 4.1 5.1
##  [523] 3.9 4.6 4.0 5.8 4.5 5.7 5.0 3.6 3.2 4.2 3.1 3.9 4.6 5.3 4.8 4.1 3.3 4.8
##  [541] 4.3 2.7 4.0 4.2 4.0 4.4 6.0 5.5 5.7 5.5 5.3 4.4 4.0 4.3 3.0 5.4 4.3 3.7
##  [559] 3.9 4.7 4.7 4.5 3.4 4.0 3.6 4.6 5.3 5.2 3.4 6.0 5.9 3.7 4.6 5.2 4.3 5.0
##  [577] 4.5 3.8 3.7 3.5 4.4 5.1 4.9 4.2 4.6 4.9 4.3 5.1 5.0 5.2 4.2 4.3 3.6 3.7
##  [595] 4.4 5.7 4.9 4.6 4.6 4.0 6.0 6.1 5.5 5.1 4.7 5.1 4.8 4.5 3.3 3.5 6.2 4.5
##  [613] 5.1 5.1 4.9 4.6 4.1 4.4 5.6 4.8 4.7 5.1 3.6 4.0 4.2 3.5 2.6 5.7 3.5 4.9
##  [631] 4.7 6.0 3.8 3.9 3.2 4.3 4.6 5.2 3.7 6.2 5.3 4.3 4.5 5.9 4.8 4.8 4.3 4.8
##  [649] 3.8 2.8 5.7 4.7 4.9 5.0 3.8 5.2 4.3 5.4 6.0 6.7 4.8 2.9 4.0 6.1 4.9 5.6
##  [667] 5.2 3.7 3.8 4.2 3.1 3.2 5.6 4.7 3.6 5.3 4.5 4.0 6.2 5.0 3.6 4.8 4.6 3.8
##  [685] 6.3 5.1 4.8 5.6 5.3 5.0 4.0 3.1 3.3 4.7 3.2 3.8 4.8 4.6 6.3 4.7 4.8 5.0
##  [703] 5.0 5.8 3.3 4.2 4.8 6.1 4.4 4.8 3.4 3.5 4.3 4.2 3.9 5.1 4.3 3.6 4.4 4.1
##  [721] 4.2 4.1 5.0 5.9 3.4 2.8 5.1 4.9 5.8 4.3 3.3 5.2 4.3 4.2 5.8 4.3 5.0 4.1
##  [739] 5.7 6.0 5.3 4.7 4.9 5.3 4.2 4.4 5.2 3.0 4.4 3.8 4.7 3.7 3.1 4.7 4.2 3.6
##  [757] 2.9 4.2 5.1 4.0 2.8 4.2 5.1 5.0 5.2 4.6 4.4 4.3 4.3 4.0 4.2 4.3 4.4 5.0
##  [775] 4.3 4.8 4.4 1.6 3.6 5.1 4.3 5.7 4.1 5.0 4.2 3.0 3.8 4.4 3.9 4.9 3.5 6.2
##  [793] 4.5 3.5 4.4 4.6 6.4 5.2 5.0 4.6 4.2 5.7 5.3 5.2 4.3 4.9 5.6 4.0 4.8 2.6
##  [811] 4.1 5.2 3.4 4.7 3.5 6.6 5.4 4.8 3.2 2.4 5.6 5.2 4.4 3.9 3.8 4.3 4.1 3.6
##  [829] 3.8 4.2 4.0 3.8 4.0 5.5 4.2 3.5 5.1 4.5 4.2 5.2 5.8 4.9 4.9 3.0 3.6 4.9
##  [847] 4.8 4.4 3.1 5.0 4.7 5.2 4.6 5.0 5.0 5.0 4.2 4.7 3.6 4.9 3.3 4.2 4.4 6.7
##  [865] 5.9 4.4 3.9 5.9 5.4 4.3 5.5 5.1 5.3 5.1 5.8 6.1 5.3 5.7 4.1 4.5 6.1 5.5
##  [883] 4.8 3.7 5.8 3.7 5.3 5.2 4.3 3.4 3.7 5.8 5.2 2.5 5.3 4.4 4.1 2.8 4.2 4.4
##  [901] 4.5 2.9 5.9 5.6 4.1 3.7 7.3 4.6 4.4 4.4 3.0 4.4 4.0 4.3 5.2 5.3 6.0 5.6
##  [919] 6.6 4.2 5.8 3.1 5.1 5.8 4.3 4.8 4.5 4.9 6.1 5.6 4.4 4.5 4.9 4.8 4.9 3.1
##  [937] 5.3 4.9 5.0 4.0 5.2 3.3 3.0 5.2 4.2 6.1 3.3 4.8 5.3 5.5 4.5 4.0 2.9 5.3
##  [955] 3.9 5.2 4.2 4.2 2.7 5.3 5.0 4.1 3.6 5.0 4.7 4.5 4.4 4.8 3.7 4.1 4.0 4.2
##  [973] 5.7 5.4 6.0 3.9 3.1 3.9 4.7 3.8 4.4 4.0 5.7 5.9 3.8 5.9 3.7 3.9 5.9 3.5
##  [991] 4.7 4.8 3.3 5.2 4.0 6.0 4.7 3.1 4.4 4.9
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
## 2.8000 6.2025
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
##    [1] 4.9 3.2 4.1 4.0 3.6 5.0 4.9 4.1 4.9 4.0 3.4 5.4 3.5 3.6 6.3 4.1 5.3 5.9
##   [19] 4.6 4.4 5.3 2.7 4.4 5.7 5.6 4.8 4.6 4.8 5.2 3.3 6.2 4.5 3.9 4.9 4.9 4.6
##   [37] 3.8 5.4 3.7 4.6 3.8 4.2 3.5 4.3 4.6 3.8 3.2 3.0 4.1 4.7 2.8 4.3 3.8 3.6
##   [55] 4.3 6.4 5.0 3.7 3.7 7.1 5.0 3.3 2.1 4.0 4.3 4.5 4.2 5.1 4.0 3.7 4.3 4.2
##   [73] 4.3 4.1 5.5 5.2 5.2 5.6 3.7 4.9 4.7 2.8 5.1 4.8 5.9 5.0 5.1 4.2 5.0 4.7
##   [91] 6.0 3.6 4.2 5.9 4.5 4.5 5.2 3.3 3.0 4.5 4.7 4.1 4.2 4.0 4.0 4.7 4.2 5.2
##  [109] 5.4 4.9 6.3 4.7 5.5 5.0 4.3 2.7 3.7 4.4 6.1 4.2 5.2 5.1 5.0 5.0 5.1 3.1
##  [127] 5.5 3.6 3.5 4.9 4.3 5.0 4.9 4.9 4.8 5.1 4.7 5.1 5.5 3.7 4.6 4.1 5.1 4.3
##  [145] 3.0 4.0 4.7 4.5 4.0 5.2 4.3 3.1 4.2 4.8 4.9 2.8 4.3 3.4 2.7 4.1 6.7 3.7
##  [163] 4.1 4.6 4.8 4.2 5.4 4.2 4.6 5.4 4.3 4.2 3.5 4.0 2.2 4.4 3.8 4.3 4.6 4.8
##  [181] 4.1 3.2 3.6 5.3 4.4 5.0 4.7 3.9 4.0 3.9 4.4 3.9 5.7 4.5 3.1 5.1 3.5 4.0
##  [199] 4.0 4.4 4.4 3.9 4.9 5.1 3.7 4.2 5.4 4.2 4.4 4.9 5.6 3.9 3.6 4.3 4.0 3.8
##  [217] 4.5 4.6 4.8 5.5 3.8 4.6 4.5 4.2 4.5 5.1 4.2 3.7 4.8 5.6 5.7 6.1 6.5 4.9
##  [235] 3.8 5.0 4.3 4.4 3.5 4.9 5.4 3.7 5.2 4.1 5.1 3.3 3.3 6.6 4.3 4.7 4.1 5.6
##  [253] 3.2 4.1 5.7 4.0 4.0 4.5 4.8 4.2 5.0 2.8 4.4 3.7 5.1 3.7 5.2 3.7 2.9 4.7
##  [271] 4.3 3.9 6.4 3.3 6.6 2.5 5.0 5.6 4.7 2.9 4.7 3.8 5.5 4.6 5.5 3.9 2.8 6.1
##  [289] 4.5 4.3 4.9 4.1 6.3 4.1 5.3 5.3 6.3 4.9 4.1 5.7 4.7 5.1 4.7 3.7 4.5 4.0
##  [307] 3.9 3.3 4.9 4.8 5.2 4.6 5.0 5.3 3.0 5.8 4.1 4.7 4.3 4.0 5.1 3.8 2.6 5.4
##  [325] 4.9 4.9 3.1 5.6 4.8 3.9 4.0 4.6 2.4 3.3 5.2 4.8 3.8 3.4 4.6 3.8 4.4 5.5
##  [343] 4.3 5.6 4.2 4.6 3.8 3.7 5.2 4.6 6.6 4.9 3.4 3.5 3.5 3.7 3.1 2.9 5.6 6.5
##  [361] 5.4 4.2 5.0 3.6 3.5 3.7 3.3 5.6 4.5 4.0 4.5 4.9 5.2 4.5 3.7 5.4 4.2 5.4
##  [379] 3.9 4.7 3.7 4.7 4.7 3.9 3.3 5.4 5.5 3.5 5.7 4.9 4.3 5.0 4.3 5.5 4.5 6.2
##  [397] 4.1 4.5 5.4 4.8 4.7 4.3 5.4 5.8 4.7 5.8 2.8 5.2 3.5 5.3 4.7 4.7 4.2 4.2
##  [415] 3.2 4.4 4.2 4.4 5.3 4.2 5.1 4.9 4.3 4.4 5.1 3.7 4.4 4.1 5.5 4.0 5.7 1.6
##  [433] 3.0 3.9 4.9 5.0 5.2 3.6 5.0 4.5 4.9 3.3 4.2 5.2 5.7 4.9 6.0 5.2 5.4 3.0
##  [451] 2.6 5.1 3.9 3.5 5.5 4.6 5.1 6.2 4.4 5.3 3.4 5.7 4.3 4.2 3.9 5.5 5.8 3.9
##  [469] 3.3 4.1 3.8 4.9 5.4 3.9 4.9 3.9 3.4 4.6 4.4 5.9 3.4 4.5 6.0 4.0 6.3 4.2
##  [487] 4.8 3.3 2.3 4.9 4.6 3.1 5.7 3.2 5.4 1.5 4.1 4.6 4.4 4.2 4.0 2.8 3.3 4.8
##  [505] 5.0 4.2 4.4 5.3 1.9 3.2 3.3 4.2 4.0 3.4 3.7 3.7 5.1 6.1 3.8 3.7 6.3 4.0
##  [523] 3.3 4.2 2.6 4.4 2.6 5.1 4.5 4.6 3.6 4.1 5.6 2.3 4.7 4.5 5.4 3.2 4.2 3.6
##  [541] 3.6 3.9 4.9 3.6 3.3 4.1 3.8 5.2 4.4 5.0 4.9 5.9 6.5 3.0 2.7 5.2 3.8 4.4
##  [559] 5.1 4.9 4.1 4.1 2.9 4.5 4.6 4.7 6.0 4.3 5.7 3.3 4.0 2.8 5.4 4.4 4.0 4.0
##  [577] 5.6 3.8 5.5 3.4 5.0 3.7 6.0 4.3 4.6 4.7 4.7 5.2 2.8 5.6 4.0 3.8 3.7 3.7
##  [595] 4.1 3.6 4.0 3.6 3.4 6.2 3.9 4.8 5.0 5.3 4.4 4.7 5.7 3.9 4.6 5.9 5.7 3.9
##  [613] 5.3 4.0 7.0 5.9 5.7 4.3 3.4 4.0 4.9 3.9 5.7 2.7 3.8 5.4 5.3 4.1 3.1 4.2
##  [631] 3.0 3.7 5.7 6.1 5.6 4.4 4.3 6.0 4.0 6.0 2.7 4.2 4.8 5.2 4.1 4.2 5.2 5.3
##  [649] 4.3 4.4 2.1 6.0 4.5 3.5 5.0 3.6 5.8 4.8 5.3 3.9 3.8 5.6 4.3 5.8 5.5 4.4
##  [667] 5.2 4.9 5.9 4.3 4.5 4.6 3.7 3.0 4.0 4.7 3.6 4.4 4.4 4.8 4.3 4.5 5.1 2.8
##  [685] 3.9 6.6 5.2 4.2 4.3 5.2 4.2 3.5 4.0 5.2 4.0 1.7 5.5 2.8 4.0 5.9 4.2 6.1
##  [703] 4.8 4.5 5.7 5.7 3.3 5.2 4.6 4.3 2.8 5.0 4.0 4.2 4.0 4.5 3.8 2.4 3.4 4.4
##  [721] 3.5 4.3 4.3 4.2 4.8 3.5 3.7 4.5 4.6 3.5 5.1 3.8 4.6 5.1 4.7 5.6 5.0 2.6
##  [739] 4.5 4.2 4.3 5.7 5.2 6.1 4.2 4.9 4.7 5.0 5.1 3.1 5.5 5.1 4.9 2.9 3.7 5.1
##  [757] 6.6 4.8 4.3 4.0 5.0 4.4 4.6 5.3 5.5 5.4 5.2 4.8 5.7 4.4 3.9 3.0 4.6 4.1
##  [775] 5.3 5.3 4.8 4.0 4.3 5.7 4.5 5.0 4.3 3.1 4.7 4.4 3.8 3.8 4.7 5.2 3.1 4.9
##  [793] 5.2 4.2 4.3 4.4 3.3 5.8 4.9 5.0 5.3 2.4 4.6 5.3 5.4 5.9 3.1 4.6 4.7 5.9
##  [811] 4.4 6.4 4.0 5.0 4.6 3.7 3.7 4.4 2.0 4.0 2.5 5.9 4.7 3.0 5.2 5.5 5.1 5.3
##  [829] 4.4 3.4 3.9 5.8 4.5 6.1 5.7 5.6 3.7 3.7 3.8 3.2 5.2 4.6 6.6 5.2 2.7 4.9
##  [847] 4.3 2.9 2.7 4.4 4.0 5.1 4.5 3.9 4.4 5.5 4.4 4.8 5.5 3.3 4.8 4.5 4.4 5.0
##  [865] 4.8 5.2 3.8 5.5 4.4 5.1 4.8 4.3 5.3 4.8 2.9 4.5 4.9 4.4 4.7 3.5 3.9 4.6
##  [883] 4.1 4.4 3.8 3.0 3.5 3.3 4.6 4.3 3.6 4.5 4.2 4.3 4.9 2.8 3.9 5.9 4.9 4.2
##  [901] 4.2 4.1 4.4 3.5 3.6 4.0 3.7 6.5 3.5 5.6 5.3 4.7 5.0 2.7 4.0 5.0 4.2 5.3
##  [919] 5.6 3.2 4.7 4.2 3.7 6.1 3.8 3.4 5.6 5.7 4.6 3.1 4.1 3.9 4.1 5.0 3.7 6.0
##  [937] 4.9 4.8 2.8 3.5 4.6 4.2 4.2 3.9 4.0 4.3 4.0 4.0 4.7 5.3 5.9 4.6 6.3 3.9
##  [955] 4.3 2.9 5.2 5.5 5.3 5.3 3.7 5.4 5.3 4.6 3.8 2.5 4.8 4.8 5.5 3.3 5.2 4.8
##  [973] 3.0 3.9 3.8 4.9 3.6 4.0 5.3 3.6 4.2 4.7 4.0 4.8 4.0 4.6 6.7 5.8 3.8 3.8
##  [991] 4.8 5.0 3.9 6.7 4.4 4.3 6.4 3.5 4.1 4.7
## 
## $func.thetastar
## [1] -0.0391
## 
## $jack.boot.val
##  [1]  0.44000000  0.35329341  0.22296919  0.05632530  0.01930836 -0.03899721
##  [7] -0.25205882 -0.25754986 -0.47833828 -0.49525140
## 
## $jack.boot.se
## [1] 0.928852
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
##    [1] 5.2 2.5 4.6 4.2 3.8 5.2 4.9 3.0 4.4 3.2 3.8 5.1 6.7 4.2 5.1 3.5 3.6 5.3
##   [19] 4.5 4.9 3.1 5.1 3.6 2.6 4.7 3.8 5.7 2.7 3.2 5.9 2.8 4.9 6.8 4.9 4.1 4.8
##   [37] 4.6 5.2 4.8 5.2 6.2 4.0 3.6 4.4 5.1 4.0 4.9 6.3 5.9 5.7 5.2 3.8 5.5 5.9
##   [55] 4.2 3.1 4.8 5.4 6.2 5.9 3.4 6.1 5.0 3.4 4.4 2.8 4.7 4.1 4.2 6.1 5.2 6.1
##   [73] 2.9 4.4 5.5 4.4 4.7 4.0 2.9 4.2 3.4 5.5 5.6 5.9 4.0 3.5 4.4 3.5 2.9 4.0
##   [91] 4.8 4.1 4.6 5.2 4.8 4.5 4.4 5.9 3.2 5.0 5.8 4.9 3.9 3.7 3.7 4.4 3.5 4.5
##  [109] 3.3 3.4 6.2 4.7 5.3 3.2 5.3 4.2 5.6 4.6 5.0 5.1 6.2 2.4 6.4 4.4 4.5 6.0
##  [127] 4.7 4.9 4.7 5.0 4.7 3.8 3.0 5.1 4.2 6.3 4.3 3.3 4.4 3.6 6.1 4.6 4.2 5.7
##  [145] 5.1 4.5 5.4 3.9 5.0 4.7 4.6 4.7 3.6 5.5 4.0 6.0 3.6 3.6 4.5 2.7 3.5 4.2
##  [163] 4.5 3.5 6.0 5.5 5.2 4.7 4.6 4.7 4.4 3.6 6.4 4.0 4.5 4.8 5.0 3.5 3.8 3.0
##  [181] 4.4 3.7 6.1 5.6 3.9 5.0 5.2 5.7 5.1 4.3 4.0 4.3 3.6 4.9 4.5 5.1 5.6 4.2
##  [199] 3.5 3.9 4.1 6.8 5.3 5.6 4.4 3.6 4.7 4.8 4.0 5.2 4.0 5.1 4.1 4.6 4.7 4.8
##  [217] 3.2 4.2 5.3 5.2 5.2 3.8 5.4 3.0 4.4 6.5 3.2 4.0 3.7 5.1 5.6 3.7 3.8 3.1
##  [235] 3.8 5.5 5.5 3.0 2.7 4.3 4.2 3.8 5.2 4.5 4.7 3.5 4.8 3.9 4.2 6.3 4.4 4.0
##  [253] 4.3 5.3 3.6 3.9 6.4 5.8 5.2 5.1 3.7 5.0 4.4 7.0 2.4 4.1 4.1 4.3 5.8 5.9
##  [271] 3.0 3.7 2.5 3.6 4.4 4.0 6.1 4.5 3.3 4.1 4.8 5.6 4.1 4.6 3.2 5.1 5.0 4.2
##  [289] 4.8 4.9 5.8 3.5 4.4 3.3 4.5 3.9 3.5 3.9 4.3 3.6 4.0 4.6 5.5 4.5 5.4 5.1
##  [307] 5.1 3.5 3.5 4.9 3.8 5.7 3.8 5.5 4.7 3.4 4.6 5.7 4.3 4.6 4.0 4.3 4.7 3.8
##  [325] 3.0 5.7 4.9 6.1 5.6 5.0 5.4 5.0 2.6 4.7 3.4 5.7 3.8 3.7 2.9 4.9 5.8 5.1
##  [343] 3.6 4.7 4.3 4.9 4.2 4.3 5.1 4.5 4.4 4.7 4.9 3.4 5.2 3.2 4.3 4.6 3.9 3.1
##  [361] 5.0 4.0 3.4 4.6 5.4 3.8 4.6 4.4 3.9 5.0 4.2 4.0 4.4 5.8 4.8 4.1 4.4 5.6
##  [379] 4.5 6.1 4.9 4.4 4.2 5.6 3.8 4.5 4.8 5.3 4.2 3.7 5.3 1.8 3.1 4.2 4.3 4.5
##  [397] 5.8 6.6 5.3 4.8 5.6 3.9 4.0 3.0 5.2 3.0 4.6 3.5 4.1 2.8 5.3 5.4 5.1 5.1
##  [415] 5.8 5.4 4.0 2.1 4.6 3.6 4.3 3.5 4.6 3.8 4.5 4.8 3.5 3.0 4.3 3.7 5.0 3.9
##  [433] 5.5 4.9 4.5 6.9 4.1 6.5 6.0 5.3 4.9 3.3 4.6 4.6 4.5 3.2 3.4 4.0 4.2 3.8
##  [451] 5.4 5.9 4.5 5.4 4.2 4.7 5.6 5.8 4.8 3.0 4.8 5.2 3.6 4.9 5.2 4.9 3.7 4.2
##  [469] 5.5 3.5 6.0 5.1 4.2 4.2 5.4 4.4 3.3 4.3 2.7 5.2 3.8 4.9 5.2 3.8 3.6 4.2
##  [487] 5.6 4.8 4.5 4.8 6.2 4.7 4.7 5.2 4.4 4.5 4.8 4.8 4.1 4.3 5.7 3.5 5.0 4.7
##  [505] 5.5 4.2 5.5 6.7 3.4 5.2 3.6 4.1 4.5 3.7 3.8 5.9 3.7 3.2 4.7 5.4 4.6 4.3
##  [523] 5.0 5.1 2.8 5.9 4.8 6.0 5.6 3.6 5.2 4.8 3.7 4.8 3.7 5.4 4.2 4.0 3.8 3.2
##  [541] 3.9 4.5 3.4 5.5 5.7 5.1 4.8 4.2 4.3 4.1 5.8 3.9 5.5 4.3 4.0 6.0 5.6 4.0
##  [559] 5.3 4.1 4.4 5.2 4.3 4.7 4.4 4.4 4.6 6.5 4.2 4.0 2.8 6.0 4.5 3.8 4.6 3.7
##  [577] 4.6 4.9 3.3 5.1 4.3 3.7 4.8 3.6 3.4 5.7 4.5 5.2 4.1 5.2 4.8 4.4 4.3 3.6
##  [595] 6.9 4.6 2.7 2.9 6.2 4.0 3.8 3.8 3.7 3.1 4.9 6.1 5.5 3.7 5.6 4.1 4.4 5.9
##  [613] 3.4 4.5 2.4 4.5 4.6 4.7 4.0 4.8 4.1 5.2 5.4 4.9 4.6 4.0 3.3 6.0 5.0 4.1
##  [631] 3.6 5.5 6.0 5.5 4.7 3.8 4.3 2.6 6.0 3.8 4.5 5.4 4.7 5.8 4.8 4.0 5.8 4.6
##  [649] 4.1 6.1 4.0 4.3 5.2 4.0 5.1 3.3 4.9 5.5 5.3 4.2 4.1 5.5 4.1 5.6 4.7 4.1
##  [667] 4.0 4.2 3.8 3.7 4.8 4.9 3.5 6.2 3.1 5.1 3.2 3.6 4.8 6.5 3.4 3.9 3.9 4.2
##  [685] 4.1 3.9 5.3 4.5 4.0 3.2 6.4 4.1 4.6 4.6 4.8 2.6 4.5 4.4 4.5 4.2 3.2 4.8
##  [703] 3.0 4.8 3.4 4.4 3.0 3.1 3.3 4.8 4.6 6.2 4.3 5.2 3.5 6.5 4.5 4.2 4.6 6.0
##  [721] 2.7 3.6 3.7 5.3 4.5 4.4 5.0 5.1 3.4 3.4 3.6 4.7 4.9 3.7 4.5 4.6 4.6 4.5
##  [739] 6.3 2.6 5.1 3.7 4.1 3.1 4.8 4.3 4.3 4.4 4.8 4.3 5.2 3.5 4.2 3.6 5.1 4.7
##  [757] 2.1 3.7 4.4 5.0 3.2 4.8 3.6 7.2 3.1 3.8 3.9 4.9 2.8 5.3 5.2 4.9 5.6 3.0
##  [775] 3.7 4.1 4.4 4.1 5.9 4.6 6.2 3.1 3.7 4.9 4.3 3.6 3.5 3.6 6.0 4.3 4.8 5.3
##  [793] 4.3 4.5 4.9 3.8 4.3 4.3 4.6 3.9 5.9 4.9 3.2 4.3 4.1 4.6 4.9 3.6 3.7 2.8
##  [811] 3.4 5.3 3.7 5.0 4.4 3.8 2.7 4.5 5.4 3.8 4.1 3.3 3.9 5.0 3.5 3.2 4.8 5.1
##  [829] 4.7 2.5 3.4 3.5 5.2 4.9 3.9 4.2 6.0 3.5 4.1 6.4 4.7 4.8 5.5 3.5 3.5 4.0
##  [847] 5.1 6.1 4.3 6.7 4.8 4.1 5.1 3.3 3.9 4.7 4.3 5.6 5.2 2.7 3.8 4.7 5.5 2.9
##  [865] 5.0 2.9 6.3 4.7 4.1 2.9 4.4 5.4 4.7 5.2 5.2 3.4 5.0 3.9 4.7 5.4 5.4 3.5
##  [883] 3.0 5.1 5.0 5.7 3.8 4.7 3.9 4.3 4.3 4.7 5.4 6.1 5.3 4.6 4.9 3.8 4.3 5.0
##  [901] 3.5 3.5 4.2 4.3 5.3 3.4 4.1 4.1 6.2 4.5 5.6 6.0 4.8 4.6 4.5 4.5 3.7 4.1
##  [919] 5.5 5.4 4.6 6.0 2.8 4.0 5.0 5.7 5.5 3.4 4.5 5.1 5.4 4.8 4.5 4.1 4.2 3.2
##  [937] 3.3 5.3 3.7 4.5 3.5 4.4 4.6 5.3 5.7 5.2 3.1 3.9 5.8 5.5 3.8 7.2 3.2 5.0
##  [955] 5.4 5.6 5.8 4.4 5.0 5.0 6.1 3.4 4.7 3.3 4.5 4.4 3.7 4.3 3.5 3.6 5.2 3.5
##  [973] 5.2 6.0 4.8 3.6 1.9 4.2 3.0 4.5 3.5 3.1 4.1 3.4 4.5 4.1 5.2 5.1 3.6 4.9
##  [991] 4.1 3.2 5.6 2.9 3.6 5.5 6.7 4.2 4.3 6.4
## 
## $func.thetastar
## 72% 
##   5 
## 
## $jack.boot.val
##  [1] 5.500 5.400 5.308 5.200 5.200 5.016 4.800 4.800 4.600 4.428
## 
## $jack.boot.se
## [1] 1.015558
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
## [1] 0.8813886
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
##   4.689010   7.410484 
##  (2.026798) (3.381014)
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
## [1] -0.1975504  0.5102819  0.4082862 -0.1996824  0.9943933  0.6484109
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
##    [1]  0.830923148  0.342289221  1.878894647  0.122308537  1.325965405
##    [6]  1.112189658  0.693770391  1.322846895  0.793441253  0.312543767
##   [11]  0.586014602  0.946799716  1.416777962  0.548549258  0.209340501
##   [16]  1.056354500  0.888348132  0.207379274  0.259303556  1.181670893
##   [21]  0.403407782  1.293280873  0.858943098  1.103909961  0.889571257
##   [26]  0.880922889  0.855956540  1.032726537  1.219640462  0.251644217
##   [31]  0.865243664  0.169820434  1.303944697  1.472352982  0.660409095
##   [36] -0.221107710  0.359349607  1.268360173  1.184669587  0.352885076
##   [41]  0.170133176  1.397514620  1.268469258  1.359626276  1.314145339
##   [46] -0.112056676  0.980896557  1.589840140  1.567201533  1.377919850
##   [51] -0.165265718  0.482759450  0.882237456  0.120003151  1.896956759
##   [56]  1.188992288  1.288754602  0.557789949  0.503171937  2.339005153
##   [61]  1.164637164  1.545035343  1.608320089  0.998908585  1.443147898
##   [66] -0.096453762  0.859607844  2.167190410  0.971894224  0.141805107
##   [71]  2.493918926  0.611888774 -0.200333041  0.910333732 -0.250869947
##   [76]  0.962098375  0.582995848  1.109353484  1.591626668  0.942297010
##   [81]  0.301290351  1.026066428  1.169645196  0.108743686  2.137122733
##   [86]  0.369730242  1.449275862  1.619286418  1.123246865 -0.213097119
##   [91]  0.440850456 -0.022741056  0.437425631  0.915666180  1.016744466
##   [96]  2.277221657  0.793368188 -0.906731433  0.624096158  0.772584365
##  [101]  2.472598861  1.652143700  0.153045925  0.395502192  1.068901619
##  [106]  1.123224786  0.709483888  1.081569660  1.120007868  0.151416563
##  [111]  0.981677682  0.930526221  1.106013214  2.337371222  0.970310984
##  [116]  0.366938620  0.081507535  0.407030006  0.964545391  0.468568968
##  [121]  0.968800313  2.263529449  0.030020771  0.862844865  0.139738591
##  [126]  0.873495944  0.955752488  0.821601011  0.335701967  1.229869456
##  [131]  0.460317869  0.351428283  1.412022255  1.676999545  0.205692155
##  [136]  1.350564320  0.547245208  0.507688102 -0.228342869  0.200832650
##  [141] -0.301309627  1.401271776  0.864578327 -0.518651453  0.089622580
##  [146]  0.453895446 -0.433041559  0.219056766  0.203254862  0.584999235
##  [151]  0.882749717  1.008078186  0.358969333  1.399987032  0.699189865
##  [156]  0.561553015  1.445215431  1.525685304  0.756101186  0.186391471
##  [161]  1.399707113  1.449377430  0.864888685  0.722976771  0.442464818
##  [166]  0.513158650  0.966871555  1.418711837  1.828871914  0.633307512
##  [171]  1.623309688  0.891130891  0.946440898  0.938171657  0.772824096
##  [176]  1.891756238  0.494364249  0.909496899  0.670105075  0.762660540
##  [181]  0.471924762  0.982220493  1.731035488 -0.225014451  1.425550159
##  [186]  0.108298506  1.052993265  1.394420801  1.567196343 -0.080519667
##  [191]  0.075375922  0.894526053 -0.421391417  0.212389639  1.123786713
##  [196]  0.971894224  0.974634307  0.746487283  1.741024672 -0.310732809
##  [201]  0.520759022  1.145013075  1.387108762  1.663286439  1.501864665
##  [206]  1.088831737  0.793352203  1.612038092  1.749226642  1.010931267
##  [211]  0.631315175  0.137073070  0.830625562  1.000307869  2.074320676
##  [216]  0.524910389  0.195406940  0.410946468 -0.151010725  0.879012121
##  [221]  1.088909089  0.003495621  0.571453885  2.078982829  1.073037992
##  [226]  0.407065357  1.567859171  0.494548574  1.413231750  2.138807608
##  [231]  0.954897349 -0.297245453 -0.170327754  1.537920599  0.091071668
##  [236]  1.993546047  1.620894378  1.869101585  0.639770750  1.309637590
##  [241]  0.958461419  0.830155340 -0.054164157  0.376697746  0.998139978
##  [246]  0.914370695  0.850016918  1.620046398 -0.607543599  2.457246107
##  [251]  0.759192994  1.326101093  2.303878766  0.508663693  0.216791471
##  [256]  0.494364249  0.400902730  1.016918887  0.214877237  1.785548512
##  [261]  1.500725443  0.552358384  0.552289741  0.843815816  1.972820269
##  [266]  1.159170909  0.472844424  0.549319310  0.511623127 -0.053031060
##  [271]  1.313350168  0.048944020  0.524051807  0.101441481  1.300542099
##  [276]  0.706371553  0.750386550  1.539701175  0.770830610  1.044433265
##  [281]  0.349432752  1.401828310  2.216123610 -0.233457210  0.582197318
##  [286]  0.620530316  1.600543831  0.527375034  1.412038113  1.334917468
##  [291]  1.391449261  1.211589985  0.106133957  0.991638785  2.421801815
##  [296]  0.081511803 -0.198335914  1.167781908  1.298928550 -0.155375908
##  [301]  0.897093183  0.369682191  0.968383545  0.873980785  0.888551711
##  [306]  0.914624053 -0.538812634  0.272699796  0.991340869  1.555644676
##  [311] -0.405109372  0.447633692  0.350487998  1.491939452 -0.932578021
##  [316]  0.582411722  1.802486711  0.853706141  1.087197093  0.984074742
##  [321]  0.581417712  0.348823848  0.942079477  1.438537154  0.313784083
##  [326]  0.839171207  0.984704015  0.079915648  0.354794998  0.172998059
##  [331]  0.852624858  0.907699302  0.423550022 -0.582429946 -0.591761232
##  [336]  0.743373596 -0.595551599  1.401238732  1.708396184  0.070359311
##  [341]  0.026681110  1.261333510  1.320184012  0.491671760  0.329532790
##  [346]  1.014379914  0.332965261  1.623774100 -0.175250889  1.007997136
##  [351]  0.469049547  0.090331641  1.308259496  0.377412828  0.136810176
##  [356]  0.585177072  0.467969079  0.356087782  0.327433122  1.554628857
##  [361]  0.926439669  0.631979021  1.409461574 -0.177119744  0.850636476
##  [366]  0.888526317  0.429570092  0.185152796  1.296630414  0.781292306
##  [371]  0.962612352  0.197756105  1.491902835  1.554576881  0.375739003
##  [376]  0.629985976  0.953237511  0.382713096  0.184135585  1.076428146
##  [381]  1.084164788  2.361887780  0.866727478  0.302081013  0.059991978
##  [386]  0.192567367  0.798052232 -0.093387583  0.100748051  1.408906180
##  [391]  0.632692607  0.268721420 -0.418182069  1.233999745  1.103906732
##  [396]  0.859170012  0.202698404  0.197189117  0.643231304  1.374831093
##  [401]  1.054018240  0.787267859 -0.139423837 -1.495148480  0.544313077
##  [406]  0.553111278  0.316563618  0.489414357  2.077287814  1.808616600
##  [411]  1.438074065  0.684741435  0.951079750  0.713747210  0.538041315
##  [416]  1.364082415  2.180628818  0.834017580  0.156611783  1.180429991
##  [421]  0.245062653  1.099965087  0.496897209  1.359627075  0.384931344
##  [426]  0.885147035  2.428255541  1.062951096  1.412731431  0.075980138
##  [431]  0.189258957  2.173211373  0.998139978  0.706000178  0.621342952
##  [436]  1.928769975  1.299449755  2.501697903  0.455313529  0.239793909
##  [441]  1.764700043  1.310415436  1.445673555 -0.553860991  0.759166042
##  [446]  0.794570711 -0.436275742  0.864368946  0.425408382  0.925128761
##  [451]  0.824534160  1.574356626  0.090776941  0.548356729  0.735092748
##  [456]  1.391826488 -0.309220777  0.848477986  0.233039691  0.787768078
##  [461]  0.438229275  0.512397386  0.127402018  1.327026845  0.520176037
##  [466]  2.000392944  0.203674483  1.341689714 -0.146250471  1.070540023
##  [471] -0.154795161  0.896870753  1.104293921  0.761571778 -0.233622694
##  [476]  0.966939548  0.170839176  0.095682745  1.282565561  0.724211951
##  [481]  1.066234615  0.148614293  0.354015221  0.753417497  2.113093788
##  [486]  0.563157471  0.349189737  1.243465735  0.567237703  0.061065344
##  [491]  0.875016470  0.944498123  1.200230800  0.774591014  1.286529027
##  [496]  0.548562396  1.295553975  0.584650665  0.865383219  1.144743291
##  [501]  1.589254794  1.470706162  1.457457985  1.051087966  1.515653521
##  [506]  1.457017394  0.455196535  0.268201722  1.147518957  1.143928081
##  [511]  0.003706970  0.135145232  2.276904699  0.549937789  0.950656022
##  [516]  0.926052536  0.840159508  1.133611328  1.763070622  2.578131001
##  [521]  0.852870159  1.010551171 -0.225844825  1.960288431  0.215244845
##  [526]  1.762881829  1.479782376  0.545831491  0.564642835  0.732889804
##  [531]  0.830140523  0.566839649  0.297045201  0.646902676  0.601184298
##  [536]  1.163239916  0.764374294 -0.599250436  0.822641025  2.277195549
##  [541]  0.528390324  2.219185406 -0.155075184  1.013254585  0.570160287
##  [546]  0.861683033  1.024853021  0.692028780  1.343040158  0.669005633
##  [551]  0.588075552  0.494212624 -0.272364237  0.411454361  0.510284905
##  [556]  0.875896158  1.549593757  0.852387849 -0.322577219  0.348425171
##  [561]  0.730987629  0.785768610  1.519328178  0.847024509  1.001308418
##  [566] -0.150991649  0.823858395  0.071265811  0.912032389  1.788002585
##  [571]  0.067633398  2.305565672  1.591207987  0.597516602  0.563538496
##  [576]  1.066573650  0.147756595  0.752791552  0.848409366  0.079877057
##  [581]  0.411540327  0.088800967  2.014114917  0.599172440  0.734645145
##  [586]  0.415111283  0.774078194  0.456410269  0.828867446  0.335798788
##  [591]  0.953719602  1.190686778  2.196083003  1.399742314  0.957880187
##  [596]  0.876741935  0.876146672  0.575860358  0.599833598  0.060776943
##  [601]  1.118671920  0.992609703  2.421966917  1.416777962  0.451308752
##  [606]  1.423148038  0.799362623  0.878798387  2.245544143  1.380933831
##  [611]  0.547245208 -0.038039635  0.177158905  1.436942821  0.499269557
##  [616]  0.014668589 -0.417291097  1.821165317  0.179239460  1.387093511
##  [621]  0.676118776  1.699379151  2.544952942  0.946617815  0.462237088
##  [626]  1.387108762  1.775195969  1.347487072 -0.715914251  0.810754013
##  [631]  0.584636403  0.842023582  1.641707232  2.153242692  0.318701469
##  [636]  1.098872116  1.877829576  0.902034819  0.841663109  0.597706246
##  [641]  0.462791824  1.447450420  0.859490721  1.939419218  2.420685919
##  [646]  0.860023034  0.421122475  0.106121782  0.162915266  1.303312468
##  [651]  0.316801579  0.762668180  0.564798748  0.942148462  0.404435847
##  [656]  0.819654950  0.484112292  0.261296032 -0.151849485  1.008593864
##  [661]  0.132707396 -0.078535817  0.203454373 -0.628502839  0.998521869
##  [666]  0.488148944  0.828456779  0.356374406  1.419073549  1.514889204
##  [671]  0.569526888  1.389720272  0.564728913  1.317996341  1.558279914
##  [676] -0.430514084  1.278607861  0.892638103  0.460472473  0.810185315
##  [681]  0.841513188  2.463416375  0.447662374  1.419185299  0.175865819
##  [686]  0.560569819  1.537796277 -0.234605333  0.148479457  0.793441253
##  [691]  0.711864576  0.501617582  1.587809033  2.184223891  0.786808860
##  [696]  0.464995476  0.476218387  1.155125544  0.649376842 -0.489273280
##  [701]  0.553558388  0.645718599  0.918580876  1.473425511  0.899776231
##  [706]  1.902707854  0.994421395  1.553306208  0.447144587  0.545014792
##  [711]  1.833518847  0.812343203  0.998139978  0.849947946  0.059388835
##  [716]  0.471504921  0.694496673  1.071100028  0.501096027  1.445271233
##  [721]  0.543444426  0.089115841  0.358676553  0.342289221  0.975192505
##  [726]  0.869223127  0.749054245  1.024884015  1.972803340  0.163164639
##  [731]  1.008636394  0.875954845 -0.075475748  0.610202621  1.090304177
##  [736]  0.568951079  0.379247517  1.019098070  0.478500406  0.168532165
##  [741]  2.265074672  0.949836466  0.219722351  0.821086937  1.261398963
##  [746]  0.221732080 -0.178541201  0.925973680  1.332275505  2.337465232
##  [751]  1.329259817  0.546454941  0.380637677  1.442622973  1.479591355
##  [756]  0.949270384  1.419264268  1.121515027  0.908029183  1.519344449
##  [761]  1.421017653  0.937212277  0.438180266  0.248406342  1.521316343
##  [766]  0.154006277  0.006518253  1.400830342  0.756440875  0.495689184
##  [771]  0.598143272  1.320914750  0.817578078 -0.305855000  0.566705939
##  [776]  0.009724059  0.734244233  0.788818860  0.465299686  0.498095078
##  [781]  0.693819547  1.522038228  1.145845658  1.264181890  0.709069459
##  [786]  1.582167804  0.934136698  0.074826621  1.005261580  1.429583118
##  [791]  1.324020732  0.356567911  0.665320046  0.370392234  0.738654836
##  [796]  0.147446569  2.167190410  0.834934389  1.261353903  1.494137843
##  [801]  0.542960991  0.481627647  0.469911532  2.021132234  0.957195247
##  [806]  0.423415097  1.167781908  0.130234728  0.567696161  0.190672858
##  [811]  1.000467652  0.815116311  2.600849561  1.051436392  0.174524363
##  [816]  0.885147035  1.223773048  0.471592125  1.000480520  1.082778324
##  [821]  1.336515793  0.749020141  1.373317464  0.495709445  0.944498123
##  [826]  1.150994019  0.959209084  1.185821469  1.005061906  0.037816017
##  [831]  0.555460898 -0.146262251 -0.200652247  0.215524978  1.333402473
##  [836]  0.568478545  0.891130891  1.402713548  0.940357091  0.675659117
##  [841]  0.768870614  0.347000133  0.268201722  0.309419793  0.566330309
##  [846] -0.283545431  0.954481848  0.196497804  0.888330126  0.997583299
##  [851]  1.000211408 -0.233692485  0.894528770  0.392147689 -0.176765133
##  [856]  1.517349029  1.154959644  1.370301608  0.921954676  0.789618762
##  [861]  0.279421174  1.279006464  2.046525594  0.304314600  1.405267196
##  [866]  0.986177191  1.252857031  0.586107301  1.283316701  0.751429165
##  [871] -1.021848071  0.812739306  1.458762676  0.576726376  1.023299977
##  [876]  0.938398212  1.463389293  0.881536678  0.995398572 -0.244502251
##  [881]  0.510241047  0.870290918  1.164068859  2.016971964  0.565998481
##  [886] -0.137853668  2.244857574  0.931914155  1.401238732  0.858282127
##  [891]  1.561690522  0.982923979  1.203182737  1.049299094  0.554372172
##  [896] -0.984079554  1.615749842  0.591177825  1.445191024  1.010741794
##  [901]  1.043961122  0.442955504  1.513371050  1.091433170  0.478917553
##  [906]  1.023172340  0.917247027  1.018909506  0.539200361  0.633388160
##  [911]  0.943676662  1.598093677  0.348616156  0.720683319  1.829202751
##  [916] -0.393739031  0.455557655  0.934870683  0.508255279  0.317230505
##  [921]  0.883599591  0.458000179 -0.023047587  1.823116870  0.491174208
##  [926]  0.343336050  2.127483272 -0.199786808  0.602304776  0.950468686
##  [931] -0.228608453  1.188507567  1.308186486  0.195215856  1.218962543
##  [936]  0.624439807 -0.141927764  0.542975327  0.227932024  0.977242675
##  [941]  0.390208148  1.028407209  0.095682745 -0.551661779  0.586089524
##  [946]  1.782262287  1.358711007  0.140861043  0.789876607  0.873564900
##  [951]  1.555208603  0.227840303  2.241459098  0.855253341  2.056644048
##  [956]  1.551096486  1.253178624 -0.018305495  1.052940086  1.049410129
##  [961]  0.937858702  0.464227457  1.146940839  0.375384271  0.566124295
##  [966]  1.341584488  0.618243554  1.347518875  0.481054871  0.514648251
##  [971]  0.909496899  1.231747831  0.578078416  2.293258232  1.757258291
##  [976]  0.851672705  0.463335298  0.443580943  1.493858973  1.126187273
##  [981]  1.463389293  0.593597725  0.483214780  1.014392655  0.687359095
##  [986]  0.755690884  0.554153070  1.445268207  2.571963266  0.212389639
##  [991]  0.359875229  0.991444511  0.415462394  0.733358015  0.446320870
##  [996]  0.918878722  2.331312025  0.135100158  0.601365356  1.332969288
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
##   0.63274890   0.31374730 
##  (0.09921561) (0.07015036)
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
## [1]  0.4540686 -0.8844809 -0.2051522 -0.2963754  0.5240416  0.5563924
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
## [1] 0.007
```

```r
sd(results2$thetastar)  #the standard deviation of the theta stars is the SE of the statistic (in this case, the mean)
```

```
## [1] 0.8968921
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
## t1*      4.5 0.01091091   0.8810355
```

One of the main advantages to the 'boot' package over the 'bootstrap' package is the nicer formatting of the output.

Going back to our original code, lets see how we could reproduce all of these numbers:


```r
table(sample(x,size=length(x),replace=T))
```

```
## 
## 1 3 6 7 8 9 
## 1 1 2 1 4 1
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
## [1] -0.0052
```

```r
se.boot
```

```
## [1] 0.9239733
```

Why do our numbers not agree exactly with those of the boot package? This is because our estimates of bias and standard error are just estimates, and they carry with them their own uncertainties. That is one of the reasons we might bother doing jackknife-after-bootstrap.

The 'boot' package has a LOT of functionality. If we have time, we will come back to some of these more complex functions later in the semester as we cover topics like regression and glm.

