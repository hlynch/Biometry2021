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
## 1 2 3 4 5 8 9 
## 1 2 2 1 2 1 1
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
## [1] -0.0516
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
## [1] 2.681582
```

```r
UL.boot
```

```
## [1] 6.215218
```

Method #2: Simply take the quantiles of the bootstrap statistics


```r
quantile(xmeans,c(0.025,0.975))
```

```
##  2.5% 97.5% 
##   2.7   6.1
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
##    [1] 4.7 4.3 2.0 3.5 3.9 1.7 4.5 4.9 5.5 4.4 4.1 2.8 4.7 5.4 3.8 3.7 4.0 5.5
##   [19] 4.2 4.8 5.7 3.7 5.1 6.1 3.5 4.2 4.5 5.6 3.9 4.7 3.7 5.1 4.3 5.1 5.4 5.1
##   [37] 5.1 3.6 4.7 3.6 4.5 4.7 3.5 5.9 5.4 6.3 3.8 3.8 4.9 4.8 2.9 4.4 4.2 5.0
##   [55] 4.2 4.0 4.8 2.9 4.5 4.7 3.7 4.1 4.0 4.2 7.5 4.9 3.5 3.8 3.1 4.2 5.7 4.1
##   [73] 3.7 4.0 4.1 3.4 5.5 3.9 4.1 3.2 4.1 5.3 5.0 3.8 3.6 5.1 4.5 3.9 4.6 4.4
##   [91] 5.4 5.2 4.8 4.2 5.4 4.6 5.4 4.7 5.0 4.6 4.1 4.0 4.0 3.1 4.9 4.0 3.6 3.3
##  [109] 5.8 4.6 3.8 3.4 4.2 4.9 3.6 4.4 3.8 3.9 4.1 5.5 6.9 4.3 3.9 5.8 3.5 4.1
##  [127] 3.1 6.0 3.9 5.0 4.4 5.0 4.9 2.9 4.0 3.4 4.2 4.9 4.1 4.4 5.1 5.2 5.1 5.3
##  [145] 4.2 5.3 4.6 5.8 5.3 3.2 4.9 3.5 4.6 4.7 5.0 4.0 4.9 4.2 5.1 4.2 5.0 4.3
##  [163] 3.8 4.3 4.4 5.6 3.9 4.1 4.6 4.8 4.4 3.1 3.7 5.5 5.0 4.4 5.2 5.1 5.4 3.8
##  [181] 3.9 4.4 3.7 4.4 4.3 4.3 6.3 4.7 4.2 4.5 3.8 5.8 5.6 5.3 4.0 5.4 4.5 5.0
##  [199] 5.9 2.7 4.3 4.4 4.8 5.1 4.3 5.8 5.3 5.5 4.2 4.3 4.0 3.4 4.9 5.0 4.6 5.0
##  [217] 3.9 4.5 2.4 3.9 5.0 3.1 5.6 5.0 4.2 4.7 4.4 5.9 3.3 5.9 5.7 4.5 5.6 6.3
##  [235] 3.9 4.0 4.2 6.3 5.0 4.7 3.7 5.6 4.7 6.0 5.9 4.6 3.5 4.3 2.6 6.1 4.4 3.3
##  [253] 6.1 4.0 4.1 4.5 5.1 4.8 6.2 2.6 4.3 5.2 5.5 4.6 3.9 6.0 5.6 5.5 4.9 3.4
##  [271] 3.8 5.2 5.4 3.3 4.6 5.2 4.2 3.5 4.2 4.5 3.4 3.9 4.7 4.7 6.1 5.7 3.9 3.7
##  [289] 4.4 6.2 4.9 4.8 6.6 6.1 4.1 6.4 3.7 3.9 4.7 3.5 4.5 6.8 6.1 4.3 4.4 5.5
##  [307] 5.0 4.7 3.9 5.0 2.6 4.8 4.8 3.0 4.5 6.2 6.3 4.0 5.2 4.5 6.1 3.8 4.3 2.9
##  [325] 5.8 5.4 4.6 5.8 5.0 4.6 4.0 5.5 4.5 6.1 6.2 5.6 3.5 5.2 4.8 4.9 3.7 5.9
##  [343] 4.1 3.4 3.4 5.2 4.0 3.6 6.3 3.4 2.1 2.9 4.8 4.1 3.8 5.1 4.6 3.7 3.7 5.4
##  [361] 4.4 4.1 5.6 4.2 3.4 3.6 3.3 5.5 4.4 2.2 4.4 3.8 4.1 4.5 5.8 4.8 5.1 4.3
##  [379] 5.1 3.0 2.6 4.6 4.1 3.3 4.0 4.1 4.5 4.7 4.1 4.3 3.4 6.0 4.0 4.3 4.9 2.8
##  [397] 4.1 5.8 4.8 5.5 5.3 4.6 3.6 4.2 3.3 4.3 5.7 4.5 4.2 5.3 4.8 5.7 4.7 4.3
##  [415] 3.4 3.8 4.7 5.1 5.0 4.0 4.8 2.7 4.3 4.2 4.5 4.6 4.2 4.4 5.2 3.0 3.5 5.6
##  [433] 3.7 4.0 5.9 3.2 4.5 3.1 6.6 4.8 3.7 4.0 4.5 6.2 4.7 4.4 3.7 2.7 5.1 5.5
##  [451] 4.4 5.2 4.4 3.0 3.2 4.9 5.3 4.0 6.2 5.8 4.5 6.0 6.0 3.5 3.8 5.5 6.2 4.4
##  [469] 3.6 5.1 4.7 3.3 5.4 3.4 4.2 4.5 3.7 4.0 5.0 3.9 3.5 6.2 5.2 5.2 4.2 7.0
##  [487] 5.3 5.0 4.6 4.5 3.9 3.2 3.5 4.3 3.8 4.1 5.0 4.9 3.9 4.1 3.9 6.1 4.6 3.9
##  [505] 4.5 6.1 4.3 3.6 5.7 2.7 3.8 4.9 5.3 4.9 4.4 6.6 5.4 6.7 5.8 4.2 4.9 4.1
##  [523] 4.7 5.3 4.9 3.7 5.3 4.4 3.0 4.8 4.5 2.0 4.1 5.5 5.0 5.4 5.5 3.2 4.2 5.1
##  [541] 5.3 2.3 4.6 2.8 5.4 4.0 4.4 5.0 3.4 4.7 4.8 4.9 3.6 6.1 2.3 3.8 5.5 4.8
##  [559] 3.4 5.6 3.4 2.2 4.3 4.9 4.2 2.8 4.4 5.1 4.8 4.5 3.5 3.8 2.0 5.4 4.6 3.8
##  [577] 3.7 5.7 4.5 3.8 4.9 4.3 6.0 6.8 4.3 3.3 4.6 4.7 3.0 2.8 3.5 3.5 5.4 4.8
##  [595] 5.2 4.6 4.3 4.8 5.3 4.5 2.3 5.1 5.3 3.3 5.2 5.4 3.9 5.9 5.4 4.3 4.3 3.5
##  [613] 6.2 5.7 6.2 4.0 3.9 4.6 2.9 3.9 4.3 5.4 4.3 4.6 4.8 2.8 5.0 4.8 4.0 4.1
##  [631] 4.4 3.3 3.5 3.7 4.5 2.8 4.2 4.9 3.8 4.1 4.6 3.0 4.0 4.5 5.7 3.5 5.1 5.1
##  [649] 5.4 6.4 5.6 5.6 3.8 5.6 4.3 4.8 4.3 2.7 4.1 4.1 3.6 2.9 3.6 3.5 4.8 3.8
##  [667] 5.0 4.5 4.0 2.5 4.2 5.8 6.0 4.4 3.9 4.4 4.6 5.4 3.1 4.1 4.5 4.1 3.4 3.4
##  [685] 6.4 3.2 5.0 2.7 5.1 5.7 4.7 5.3 2.7 4.8 5.0 4.5 4.6 4.1 5.7 5.5 6.3 4.1
##  [703] 3.9 5.6 4.0 3.9 3.6 4.6 5.4 4.1 4.7 3.7 4.0 5.0 2.6 3.7 4.4 5.5 4.3 4.8
##  [721] 4.2 3.7 4.3 3.5 4.9 3.9 4.8 3.0 5.8 5.4 4.0 5.0 3.1 4.3 5.3 5.0 4.6 4.2
##  [739] 5.1 5.0 4.2 3.4 5.5 2.6 4.1 5.1 3.0 4.4 4.3 3.8 4.6 4.8 4.4 4.6 5.5 3.9
##  [757] 4.6 4.7 6.3 3.8 5.5 3.3 4.1 4.2 3.6 3.4 4.4 5.6 4.5 5.3 5.2 3.8 3.8 4.3
##  [775] 5.6 4.2 4.1 3.1 4.4 3.9 4.9 5.5 4.6 3.4 4.3 4.5 4.8 4.2 4.4 5.0 5.0 4.2
##  [793] 6.4 4.1 5.2 4.6 5.0 4.4 3.0 3.6 4.4 5.7 5.8 4.6 5.4 5.4 3.9 5.2 5.8 5.3
##  [811] 5.1 4.4 3.9 3.3 6.1 5.4 5.9 6.5 5.6 5.5 4.4 4.3 5.5 1.9 4.2 6.1 5.5 3.2
##  [829] 4.3 4.5 4.1 2.8 5.9 4.7 5.2 3.6 4.3 3.4 4.1 3.8 3.5 3.4 4.4 3.4 3.5 4.7
##  [847] 5.2 4.7 3.0 4.3 4.9 5.2 3.8 5.4 5.0 3.1 3.7 4.9 5.8 4.6 4.7 4.5 2.9 4.5
##  [865] 5.1 4.8 4.8 3.2 5.2 5.4 3.1 5.0 5.4 4.4 4.7 4.7 5.5 3.7 4.6 5.9 5.3 5.2
##  [883] 4.2 3.0 4.8 3.4 4.5 4.1 5.4 3.4 3.8 4.3 4.0 4.5 3.3 4.6 4.7 4.8 4.6 4.6
##  [901] 5.3 2.9 3.9 2.6 4.8 4.2 5.1 3.9 5.6 4.6 4.8 5.0 4.2 4.1 4.0 4.9 4.5 4.4
##  [919] 2.3 3.6 4.4 4.6 5.4 4.5 4.3 5.6 4.9 5.5 4.9 5.6 5.2 4.6 4.1 6.1 3.8 3.6
##  [937] 5.1 3.5 6.1 3.4 5.4 3.5 4.3 4.2 2.7 4.1 4.0 4.2 3.9 4.9 4.6 6.6 5.0 5.5
##  [955] 4.4 5.1 6.3 5.4 5.8 5.0 4.9 5.2 5.3 2.9 3.4 5.3 3.5 3.2 3.7 2.9 5.2 4.9
##  [973] 5.4 3.3 5.0 4.9 4.3 4.4 4.3 3.5 4.9 4.8 3.4 4.1 4.9 4.2 1.5 3.5 3.3 2.7
##  [991] 4.2 3.9 5.4 4.2 5.0 3.6 5.1 5.1 4.4 5.1
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
##    [1] 3.6 5.3 4.9 4.5 5.4 3.8 5.0 3.6 4.6 4.7 4.2 5.5 3.2 5.2 3.6 3.9 6.2 4.2
##   [19] 4.1 4.0 4.1 5.5 2.9 3.1 3.8 5.7 5.1 3.4 5.1 4.8 5.1 4.1 3.9 4.5 5.2 3.7
##   [37] 5.3 3.1 3.7 4.5 5.3 4.7 4.6 3.9 4.0 3.2 3.8 4.3 4.0 3.5 4.7 5.4 4.9 5.1
##   [55] 4.4 4.6 4.3 5.7 3.7 3.4 4.6 2.7 4.6 3.8 3.7 4.0 4.8 4.7 5.2 4.6 4.8 3.2
##   [73] 6.1 6.6 4.0 4.5 4.7 5.5 4.6 7.3 5.3 4.1 5.0 4.2 4.1 5.4 5.7 6.4 4.2 4.4
##   [91] 3.2 4.2 4.0 4.3 3.1 5.8 5.2 4.9 4.3 5.0 5.7 4.7 4.0 5.0 5.2 4.3 4.7 5.1
##  [109] 3.7 3.9 3.4 4.3 3.0 5.0 3.5 4.1 3.7 5.5 3.8 4.3 3.9 4.6 5.0 4.3 5.5 5.4
##  [127] 4.1 5.7 3.7 4.2 4.3 5.0 5.9 4.5 4.8 4.3 3.8 5.0 4.5 5.1 4.4 4.9 3.7 2.5
##  [145] 5.1 5.1 5.2 5.9 4.6 4.8 5.9 4.6 4.0 5.3 5.4 2.0 4.6 5.3 3.5 4.9 4.2 6.2
##  [163] 4.9 4.9 3.9 4.6 4.4 4.5 5.0 2.2 3.3 4.1 4.7 4.3 4.9 4.6 4.4 5.1 3.9 4.1
##  [181] 5.2 6.4 4.1 6.4 5.1 6.1 4.4 4.0 5.2 5.7 4.8 5.3 5.0 6.0 4.1 3.3 3.8 5.2
##  [199] 4.1 4.8 3.4 3.5 3.2 5.5 3.3 4.7 5.8 5.4 4.7 4.4 4.3 4.8 4.3 4.9 3.0 4.4
##  [217] 4.5 4.2 3.3 5.5 5.4 5.6 4.1 6.5 4.3 5.4 3.2 2.6 6.5 5.9 3.9 4.3 2.7 5.0
##  [235] 4.2 3.8 5.6 3.6 5.2 5.1 3.4 5.8 5.8 2.8 4.3 4.0 6.0 4.2 5.6 4.3 4.8 4.7
##  [253] 2.6 5.4 5.1 4.6 4.2 5.0 3.9 3.9 4.4 3.8 4.5 4.6 4.4 4.5 4.9 4.1 4.2 4.8
##  [271] 3.9 4.9 4.8 5.4 3.8 5.7 3.1 2.8 3.8 3.5 5.3 3.8 5.4 3.3 5.1 4.3 3.2 3.2
##  [289] 2.3 4.7 5.3 4.6 4.6 4.7 5.2 3.5 5.9 5.4 5.4 3.6 4.6 3.8 5.8 5.4 3.5 5.5
##  [307] 5.3 4.9 4.5 3.8 4.2 3.7 3.4 5.9 4.6 4.2 5.4 4.7 5.7 3.0 4.3 4.7 3.5 3.6
##  [325] 3.4 4.5 3.8 3.9 4.7 4.9 3.8 5.7 5.5 3.9 5.2 4.2 5.7 5.8 4.5 6.8 4.2 4.1
##  [343] 3.9 4.3 5.2 4.6 5.8 4.4 4.7 4.1 4.7 2.8 4.3 4.7 4.0 5.0 5.5 5.3 4.5 5.1
##  [361] 5.4 4.0 4.0 5.4 4.7 5.5 6.2 5.8 4.2 5.1 5.5 3.6 4.9 4.5 4.0 6.5 3.2 4.2
##  [379] 4.2 5.1 4.8 5.0 5.2 5.5 4.3 5.8 6.3 5.6 5.0 5.5 4.0 4.5 5.4 4.0 4.4 4.8
##  [397] 5.9 3.1 3.7 3.7 2.9 5.7 5.5 3.7 4.9 4.5 4.6 5.1 3.9 6.6 4.1 5.1 4.8 4.4
##  [415] 4.5 4.5 5.3 4.4 2.4 6.5 2.9 2.5 4.7 4.3 3.3 6.0 6.2 5.0 3.0 6.0 5.4 5.0
##  [433] 3.4 4.1 4.3 4.6 5.5 3.8 3.6 4.0 4.3 3.0 4.4 4.7 6.3 5.4 5.3 4.2 3.3 5.9
##  [451] 2.7 4.3 3.8 4.1 5.4 4.6 3.7 4.4 3.9 4.2 5.6 3.8 5.6 4.1 6.1 5.0 3.9 3.4
##  [469] 4.2 4.8 5.2 3.4 3.2 4.0 4.9 5.7 4.6 3.8 3.7 4.1 3.6 4.2 4.7 3.7 5.1 3.9
##  [487] 4.7 4.6 3.7 4.6 4.7 3.1 4.2 6.3 5.0 5.7 5.3 4.0 5.2 4.3 3.9 4.0 5.8 6.3
##  [505] 3.8 4.7 3.2 4.6 5.1 5.7 3.1 4.1 7.1 5.6 4.6 5.2 4.1 4.8 4.9 5.9 4.8 4.3
##  [523] 4.7 3.4 2.6 3.5 4.7 3.9 4.6 3.3 3.8 6.5 4.3 5.5 4.1 4.1 3.4 5.2 5.7 5.9
##  [541] 3.7 4.7 4.2 3.4 3.7 4.8 6.3 4.5 4.5 3.8 2.7 5.0 5.6 4.8 5.7 4.6 3.4 4.6
##  [559] 5.0 6.8 6.1 3.9 4.3 4.9 6.0 5.0 6.3 4.3 5.0 4.3 4.8 5.4 6.7 5.3 3.7 4.2
##  [577] 4.8 4.7 3.5 3.6 5.0 4.8 4.3 5.0 4.7 3.7 3.7 4.9 6.0 3.8 4.1 6.3 4.8 4.7
##  [595] 3.0 3.9 6.0 3.5 4.7 2.8 4.8 4.5 4.0 6.4 5.3 2.5 4.6 5.3 5.0 4.2 4.4 3.3
##  [613] 3.4 4.3 4.2 5.4 2.3 4.9 4.1 5.9 5.7 4.5 6.1 5.5 4.1 2.5 4.7 6.2 4.8 4.6
##  [631] 4.2 5.7 3.7 4.7 3.4 4.9 4.2 4.9 3.0 4.8 4.2 5.1 4.9 3.8 4.0 4.4 4.8 6.0
##  [649] 2.1 4.1 4.3 4.8 3.8 4.9 3.5 4.7 2.7 4.9 3.2 4.8 3.5 3.8 4.3 4.5 4.9 4.7
##  [667] 4.7 3.8 3.4 4.7 4.3 4.4 4.1 4.5 3.8 3.8 4.9 3.8 4.2 3.2 3.7 6.3 4.4 6.3
##  [685] 3.2 3.6 3.2 4.0 3.1 5.0 4.4 3.6 4.7 3.3 3.8 4.4 3.3 4.7 4.3 4.9 2.9 4.1
##  [703] 2.6 4.0 5.2 5.1 4.1 4.8 4.1 4.8 3.2 5.4 6.0 5.4 2.7 3.0 4.5 3.5 6.5 5.5
##  [721] 4.8 4.1 4.6 6.0 5.4 4.5 3.7 3.7 5.3 5.6 5.1 4.7 3.4 5.8 4.1 3.3 4.1 4.0
##  [739] 4.0 4.5 2.9 5.1 4.1 5.0 3.2 4.0 4.1 6.0 3.8 4.2 3.6 4.7 3.4 4.9 5.2 5.2
##  [757] 4.2 3.9 4.3 6.8 3.1 4.7 4.3 3.2 4.1 4.4 4.6 4.2 4.3 5.5 4.5 3.6 4.5 4.8
##  [775] 4.9 5.1 4.8 4.5 3.9 5.4 4.1 6.0 6.2 4.4 5.9 4.0 4.3 3.7 4.1 4.0 3.4 4.4
##  [793] 3.2 4.9 3.3 4.3 3.6 3.8 3.9 4.4 3.8 5.0 5.1 5.7 6.1 4.9 5.3 4.3 4.6 5.0
##  [811] 3.7 3.6 3.9 4.2 3.1 3.1 4.2 5.5 4.6 3.9 4.0 5.0 4.5 4.2 4.4 4.3 4.9 4.6
##  [829] 4.8 5.8 5.1 2.7 4.0 4.5 5.2 5.2 4.2 5.1 5.0 2.9 3.9 4.1 4.2 4.4 4.3 5.5
##  [847] 4.1 5.0 2.9 5.1 3.3 4.1 3.9 5.1 4.4 4.9 4.1 3.3 4.5 5.2 3.4 3.2 4.6 4.1
##  [865] 4.6 3.7 3.9 3.2 4.1 5.1 4.4 5.0 2.4 4.3 3.6 6.1 5.4 4.5 4.6 4.5 4.8 3.8
##  [883] 4.9 4.1 3.5 3.7 5.9 2.1 4.2 3.4 3.8 5.9 3.8 3.5 3.0 5.8 3.8 4.7 4.1 6.4
##  [901] 6.3 4.3 5.1 4.8 4.4 3.8 3.6 5.4 4.7 4.4 3.8 6.2 4.5 4.6 4.5 5.4 3.6 2.3
##  [919] 5.0 4.9 5.0 4.8 3.3 4.9 5.7 4.1 5.1 5.2 5.3 3.2 4.2 3.9 5.4 5.1 5.4 5.6
##  [937] 5.7 4.2 5.4 3.7 3.2 4.0 5.4 5.0 4.3 5.6 4.2 5.9 5.2 3.4 3.5 5.4 4.5 4.8
##  [955] 5.0 5.2 5.6 4.5 5.0 4.3 5.1 5.4 4.8 4.4 3.7 4.8 5.7 6.0 5.9 3.5 3.9 4.7
##  [973] 4.8 4.4 4.1 3.7 5.2 6.2 2.4 4.3 5.4 5.5 4.7 5.1 4.6 5.8 5.4 5.8 4.9 5.0
##  [991] 3.6 4.6 4.3 3.7 4.0 4.8 5.4 5.5 6.3 2.9
## 
## $func.thetastar
## [1] 0.0164
## 
## $jack.boot.val
##  [1]  0.49065156  0.42543860  0.32634921  0.15950413  0.09402174 -0.06434316
##  [7] -0.10924370 -0.22550725 -0.35412088 -0.52611465
## 
## $jack.boot.se
## [1] 0.9613843
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
##    [1] 5.4 4.5 4.4 5.8 2.8 5.3 4.5 4.2 5.4 5.3 4.0 4.3 5.6 3.5 3.9 3.8 4.0 5.2
##   [19] 4.9 5.2 5.9 5.8 6.4 5.1 5.7 3.8 2.9 4.2 3.5 4.3 3.3 5.7 3.6 4.3 5.8 4.9
##   [37] 4.0 5.3 5.4 6.6 4.8 4.0 4.5 5.0 6.1 4.7 4.3 5.1 3.8 6.0 4.6 4.4 5.1 3.6
##   [55] 4.8 3.6 5.2 5.4 5.5 6.7 3.8 6.7 4.9 4.7 4.3 2.4 4.2 4.0 4.8 5.2 3.3 4.1
##   [73] 5.8 4.2 5.9 3.8 1.9 5.7 4.8 4.1 5.1 4.3 5.2 4.7 3.8 5.3 5.7 5.4 4.0 4.0
##   [91] 3.6 5.8 4.9 4.9 5.7 4.0 4.8 3.9 2.8 5.6 4.8 5.3 4.1 2.1 3.4 4.9 2.4 5.2
##  [109] 4.8 5.1 4.7 4.2 6.3 4.5 6.2 4.6 4.2 5.0 3.1 4.1 3.6 4.3 4.5 3.9 5.8 5.5
##  [127] 6.0 6.0 3.7 4.5 4.2 4.7 4.8 5.0 4.9 4.7 4.2 3.6 4.3 5.1 4.7 4.0 3.9 2.2
##  [145] 4.7 5.7 2.9 4.0 4.3 4.1 2.6 4.3 4.7 4.3 3.2 3.6 4.3 3.4 3.6 5.1 3.7 5.1
##  [163] 4.7 3.7 3.7 5.8 5.0 4.3 5.0 4.9 4.4 4.5 4.2 3.5 3.4 4.8 4.6 5.0 5.9 3.1
##  [181] 3.8 6.0 5.7 4.9 4.6 5.4 3.9 4.1 3.8 4.0 2.5 4.6 5.5 5.8 3.9 4.8 6.3 4.6
##  [199] 5.5 5.3 5.5 5.3 4.6 3.7 4.1 4.1 5.0 5.0 4.2 4.3 5.1 4.7 5.2 5.5 5.5 4.7
##  [217] 4.7 5.4 2.9 4.2 5.1 5.1 2.8 4.7 5.1 5.0 4.7 3.1 3.4 5.3 3.1 4.1 3.7 2.2
##  [235] 2.4 4.3 4.8 4.0 3.4 5.7 4.2 3.1 4.2 5.0 5.4 4.3 4.5 3.8 6.3 4.8 3.6 4.0
##  [253] 4.1 4.1 6.4 4.2 5.8 4.7 3.6 5.1 4.1 5.5 3.0 3.9 3.6 2.7 3.8 4.3 5.9 4.6
##  [271] 4.9 4.7 4.0 6.0 5.1 4.0 5.2 4.2 4.6 4.6 5.2 3.1 4.2 3.8 5.1 4.2 3.7 4.2
##  [289] 4.0 3.1 4.6 5.7 6.0 2.9 2.5 3.8 3.5 3.9 6.6 4.6 4.1 3.6 5.0 5.4 3.9 5.4
##  [307] 5.6 3.8 5.6 4.4 4.9 4.4 3.4 5.0 4.3 5.3 4.2 5.1 4.0 6.0 4.1 4.9 3.8 3.7
##  [325] 5.6 4.1 4.1 5.7 4.3 4.3 6.0 3.5 5.4 5.7 4.3 6.1 4.4 5.7 4.8 4.1 3.4 2.9
##  [343] 3.3 4.9 5.2 6.7 4.8 4.7 4.5 2.3 6.3 3.8 4.2 3.5 4.5 4.0 5.1 4.4 4.6 4.2
##  [361] 4.3 4.4 4.1 3.9 5.0 4.8 5.2 4.5 5.1 3.6 6.0 5.9 3.0 3.6 5.3 5.0 4.4 5.0
##  [379] 5.4 5.5 5.3 4.7 2.5 4.1 3.4 4.5 4.4 4.5 6.5 3.9 5.2 5.6 3.4 4.3 4.9 3.4
##  [397] 4.9 5.4 4.3 7.3 5.6 3.5 4.3 4.8 5.4 5.3 3.8 3.8 6.0 4.3 5.0 4.0 5.6 5.0
##  [415] 5.7 4.5 3.5 5.0 4.4 3.3 6.0 3.1 5.5 3.7 3.2 4.9 5.5 3.9 2.7 5.1 3.9 4.3
##  [433] 4.3 5.9 3.3 6.4 4.8 4.8 5.4 4.7 5.8 3.8 5.2 5.1 3.2 2.7 4.5 3.6 4.0 4.5
##  [451] 5.5 4.3 4.2 3.1 3.2 4.8 4.1 4.2 5.4 4.7 4.6 3.9 4.1 5.8 4.2 4.0 6.2 5.5
##  [469] 3.5 4.4 4.8 4.8 4.3 4.6 4.4 4.7 4.8 4.2 4.6 4.1 5.3 5.3 6.4 5.3 5.2 2.2
##  [487] 5.1 4.0 3.5 5.0 5.2 4.9 4.4 4.6 4.4 5.0 4.9 2.8 3.4 5.0 5.9 4.2 5.1 4.6
##  [505] 4.4 5.1 4.7 4.1 4.4 5.4 3.7 4.8 4.8 4.9 5.6 3.5 6.0 4.2 3.5 4.6 4.9 4.8
##  [523] 4.3 4.3 3.8 4.2 4.4 5.4 3.4 6.6 3.6 6.3 4.3 3.2 3.0 5.3 3.6 4.4 2.8 5.2
##  [541] 4.7 4.0 4.0 3.5 4.3 6.1 5.4 4.0 3.6 5.9 3.1 3.2 5.2 4.6 4.1 3.6 4.7 3.7
##  [559] 3.7 5.8 4.5 2.9 5.0 4.5 5.3 3.8 6.6 5.4 4.2 4.0 4.3 4.9 4.8 4.9 3.9 5.6
##  [577] 4.9 3.1 2.6 4.6 3.5 4.6 3.6 3.2 5.3 5.1 4.6 3.5 4.6 4.8 6.1 4.6 4.5 5.5
##  [595] 3.4 4.8 5.2 4.7 4.0 6.1 5.1 5.4 4.7 3.3 4.1 3.6 3.4 3.7 3.7 3.5 4.9 4.1
##  [613] 3.9 4.1 5.7 6.3 4.3 3.4 4.6 4.8 5.6 4.9 5.1 3.7 4.5 4.0 4.3 6.0 5.1 6.3
##  [631] 5.0 4.1 3.0 4.4 4.4 5.3 4.8 3.7 3.4 2.4 4.9 5.9 4.0 5.8 4.5 6.2 4.8 4.7
##  [649] 7.1 4.3 5.6 3.5 5.3 4.2 4.1 3.2 3.8 4.4 3.4 6.2 4.4 4.0 4.5 3.9 4.6 3.7
##  [667] 4.6 3.8 4.5 5.0 3.7 4.3 3.2 4.0 4.7 5.5 5.2 3.9 5.8 5.4 6.1 6.3 3.8 4.2
##  [685] 5.1 3.7 5.1 4.4 5.3 4.0 4.1 3.5 3.6 3.9 4.6 4.0 4.9 2.7 3.2 4.7 5.3 5.1
##  [703] 5.0 4.7 5.3 6.3 4.9 5.0 3.9 3.5 3.6 2.4 6.1 5.3 2.4 4.2 3.5 4.3 4.5 5.5
##  [721] 5.1 4.8 3.8 5.3 4.8 3.4 6.8 4.9 4.7 4.5 4.0 4.2 5.4 3.7 5.0 4.1 3.7 2.6
##  [739] 4.8 5.4 3.6 6.0 4.1 3.7 4.8 4.7 5.7 4.2 4.1 4.8 4.5 4.9 2.9 4.2 4.5 5.5
##  [757] 4.3 4.1 5.2 5.6 4.7 3.5 3.8 5.4 3.1 2.8 4.5 5.0 4.4 6.4 4.2 4.4 5.5 3.4
##  [775] 3.4 5.9 4.5 4.1 3.3 5.5 4.2 3.4 4.1 4.6 6.1 5.3 3.6 3.3 1.7 3.6 3.2 5.0
##  [793] 3.4 5.3 6.4 4.0 2.3 4.1 5.0 3.5 2.8 5.0 5.2 3.7 2.9 3.2 3.6 4.4 4.0 4.4
##  [811] 5.8 4.9 3.0 4.6 4.4 4.7 4.2 6.1 5.7 5.5 4.1 4.1 4.3 3.9 3.3 5.7 4.8 4.4
##  [829] 3.5 4.0 5.3 3.0 5.0 5.4 2.5 4.3 3.4 6.4 4.1 4.3 4.9 2.6 5.2 3.1 3.6 2.1
##  [847] 4.0 3.5 4.4 3.9 4.0 4.1 3.7 4.6 4.6 3.3 4.9 4.0 2.8 4.6 4.6 4.3 3.1 6.0
##  [865] 5.3 3.6 4.7 3.4 4.3 4.1 3.4 2.7 5.8 5.1 5.1 3.6 5.2 4.1 4.7 4.5 4.3 4.4
##  [883] 3.8 4.2 5.2 4.1 5.5 3.2 3.5 5.6 3.7 3.8 4.5 4.6 3.8 4.3 5.2 6.0 5.8 5.3
##  [901] 3.8 5.0 4.8 3.1 2.2 5.2 6.9 5.1 4.1 3.9 4.7 6.1 4.3 5.0 3.7 5.5 4.0 5.2
##  [919] 6.2 4.5 5.0 3.4 4.4 4.5 4.7 4.6 3.4 6.0 3.7 3.9 6.4 5.7 3.0 5.0 6.0 4.4
##  [937] 4.9 5.6 4.7 5.5 4.7 3.5 5.1 5.3 4.4 4.8 5.4 4.4 6.3 3.4 4.2 4.5 3.5 4.1
##  [955] 4.0 3.7 5.9 4.3 4.7 3.6 4.8 3.8 4.0 4.0 4.7 4.9 2.5 4.8 2.8 4.2 4.1 4.3
##  [973] 2.3 5.9 5.3 3.6 3.5 5.3 4.8 3.4 4.6 5.4 4.0 5.2 5.7 3.3 4.6 2.3 3.2 3.6
##  [991] 6.3 3.2 6.2 4.4 5.4 4.8 4.3 4.7 4.8 5.3
## 
## $func.thetastar
## 72% 
##   5 
## 
## $jack.boot.val
##  [1] 5.600 5.400 5.300 5.300 5.116 5.000 4.800 4.700 4.600 4.400
## 
## $jack.boot.se
## [1] 1.105853
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
## [1] 0.3031338
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
##   2.593861   3.949881 
##  (1.093653) (1.837121)
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
## [1]  1.07613142  1.20904610  0.50700261  0.07998862 -0.15891548  0.66632459
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
##    [1]  0.857984666  0.001929169  0.136659750  1.248462619 -0.024422615
##    [6]  1.487462815  0.359328402  0.833745515  0.663106516  0.588344240
##   [11]  0.106233289 -0.228591329  0.336243347  0.118543310  0.788816723
##   [16] -0.534042044  0.625521524  1.438483902  0.319506841 -0.213263468
##   [21] -0.187049148  1.181238327 -0.256475234  0.583561598  0.184665351
##   [26] -0.075010691  0.219111066  0.425998243 -0.746574751 -0.252722577
##   [31] -0.382152804  0.596160701  0.452766524 -0.284943029  1.073333500
##   [36]  0.215155807 -0.025623605 -0.261508970  0.107619289  0.333092909
##   [41] -0.173691404  0.294277127  0.307050924  0.700037781  0.346380143
##   [46]  0.554074860 -0.158345927  0.522683688  0.266205256  0.273634669
##   [51]  0.660921101 -0.210559576  0.447284538 -0.308956824  0.012940414
##   [56] -0.083554441  0.206110576  0.970865276  0.425153167  0.333792879
##   [61]  0.062417031 -0.368110896  0.879051619  0.880395672  0.387489409
##   [66]  0.221008250  0.272465608  0.317467259  0.443967523  0.772220533
##   [71] -0.735231579  0.251146904  0.281457358 -1.276535536  1.164293104
##   [76]  0.666053467  0.058443156 -0.414323666 -0.203220367 -0.324799305
##   [81]  0.573301944 -0.417061798  0.436101807  0.297519552 -0.254022031
##   [86]  1.007611449  1.146168975  0.621547838  0.145970056 -0.114912537
##   [91] -0.139138139  0.309368420  0.231571097  0.169678645  0.806841817
##   [96]  0.510909664 -0.703505776  0.337718782  0.460240731  0.198521537
##  [101]  0.038248311  1.045873115 -0.382391212 -0.496994403 -0.854153371
##  [106] -0.194257539  0.402082593  0.013118432  0.988382680  0.520432078
##  [111]  0.084559213 -0.050346340  0.102924787  0.318935760  0.991278079
##  [116] -0.048340442 -0.166437606  0.355288055 -0.172285070  0.298827998
##  [121]  0.606810192  0.090897514 -0.002323009 -0.040391097  0.002877035
##  [126]  0.699925123 -0.005449377  0.913004067  0.535956514 -0.056593594
##  [131]  0.917548838  1.170008583  0.744830918 -0.126623354 -0.041372169
##  [136]  0.118680728 -0.059525794  0.822856817  0.452133210  0.257587291
##  [141]  0.019162002 -0.664091602  0.350687353  0.618709831  0.086936468
##  [146]  0.690110855  0.556807216 -0.273784769  0.703348744  0.303124503
##  [151] -0.165210190  0.847514045  0.300277110  0.911478618  0.412734832
##  [156]  0.064420061 -0.424428156  0.223413660 -0.246476101  0.552558101
##  [161] -0.516540399  0.486760830  0.116107305 -0.142208780  0.443264539
##  [166]  1.069432739  0.147708298  0.309054545  0.163542706 -0.450191821
##  [171] -0.395807686  0.242831695 -0.149271442  1.183085839  0.860755790
##  [176]  0.200728094 -0.144012986  0.402776152  0.679758837  0.464883831
##  [181] -0.250111112 -0.066454303  0.960692231 -0.301598465  0.593375756
##  [186]  0.678146180  0.507763603 -0.274021543 -0.014892849 -0.419818300
##  [191]  0.757816774  0.595483816  0.821220892  0.058134691  1.604024731
##  [196] -0.831637764  0.618471278  0.264774378  0.717907200  1.111292290
##  [201]  0.312601991  0.207461752  0.371346234 -0.083815078 -0.178274348
##  [206]  0.086882016  0.299850625  0.672644796  0.209837308  0.289040766
##  [211]  0.595081253  0.173637219 -0.390744806  0.955276946  0.941784926
##  [216]  0.137604351  0.437560516  0.593604690  0.754724231  0.065866474
##  [221] -0.057705223 -0.413509739  0.856695350  0.002887747  0.285110722
##  [226]  0.351753884 -0.101220444  0.508883803 -0.450715560  0.050194521
##  [231]  0.056023291  0.085646240 -0.399598563  1.777750206 -0.012915415
##  [236]  0.535949700  0.186406878  0.274651448  0.257780160  0.450943554
##  [241]  0.493724940 -0.058899179  0.129887543 -0.074863474  0.133611809
##  [246]  0.486915181  0.880157883 -0.816145628  0.093920448  0.314127912
##  [251] -0.252968786 -0.130238337  0.480345779 -0.605968818 -0.042596334
##  [256]  0.091302947 -0.548753637  0.635764717  0.337306820  0.433916247
##  [261]  0.641042342  0.781491248  0.142535238  0.190356834  0.608253406
##  [266]  0.763788270  0.663490080 -0.048975101  0.396953637  0.067747162
##  [271]  0.588249121  0.198187700  0.159258960  0.304708890 -0.753268271
##  [276]  0.020689891  0.326673156 -0.076319983  0.406729277  0.265209845
##  [281] -0.669855134 -0.008593992 -0.359756095 -0.018700594  0.409568252
##  [286] -0.329928528 -0.078781331  0.140660245  0.214697408 -0.235316683
##  [291]  0.614353660  0.476004902 -0.007899366 -1.005001639  1.103912372
##  [296]  0.200173274  0.233097707  0.641009922  0.593531654  1.374533111
##  [301]  0.271962342  0.588428541 -0.246805104  0.186719314  0.632936487
##  [306]  0.256798748  0.302731359  0.053026280 -0.279144222  0.604843920
##  [311]  0.606010832  0.530677447  0.271090488  0.170629623  0.225720451
##  [316]  0.329741427  0.688864029 -0.067283824  0.898176897  0.194570534
##  [321]  0.547295711  0.306183817 -0.126443528  0.545941480 -0.266615186
##  [326] -0.435404961  0.725834022  0.452288954  0.172916989  0.046916927
##  [331]  0.647708633  0.964921137  0.946920380 -0.045268381 -0.865955564
##  [336]  0.382168807  0.182218437  0.079026971 -0.183805733  0.399091916
##  [341]  0.258110015  0.930006304  0.925942841 -0.413559407 -0.305763605
##  [346]  0.223349855  0.271584668  0.637117435  0.347105568  0.398996950
##  [351]  0.656011810  0.802035662  0.121749196 -0.874891883 -0.843042656
##  [356]  0.406041556  0.352644950  0.196706803  0.922658189 -0.250958777
##  [361]  0.629756403 -0.442915059 -0.137010710  0.526238962  0.296534890
##  [366]  0.630260539  0.075607536 -0.457107273 -0.056121817 -0.541415309
##  [371] -0.095056164  0.362068145  0.557041249  0.021216679 -1.294412207
##  [376]  1.538400016  0.638204087  0.097348779  0.615144636  0.160030812
##  [381]  0.145574375 -0.132118630  0.123093133 -0.182532025  0.337544063
##  [386]  0.874359619  0.585690751  0.719105249  0.449965289 -0.028470496
##  [391]  0.219853308  0.504932725 -0.423399297  0.161921336  0.199824935
##  [396]  0.388573400  0.195848327  0.348822544  0.999605925  0.766166950
##  [401] -0.652956229  0.175599269  0.441270027 -0.259267921  0.308032387
##  [406]  0.323997507  0.544660452  0.066334032 -0.198221024 -0.280277886
##  [411]  0.179805958  0.225950902 -0.077752355 -0.492327499 -0.198750527
##  [416]  0.336012037  0.584098397 -0.557711974  0.771078731  0.866419491
##  [421]  0.321010494  0.712360862  0.574950419 -0.189733029  0.820080790
##  [426]  1.129956574  0.557803324  0.145434061  0.010065241  0.163608556
##  [431] -0.266594606  0.325861555  0.286981073  0.087606006  0.554182029
##  [436]  0.521740295  0.487612021  0.397479341  0.466202497  0.153086820
##  [441] -0.232584895 -0.234532829  0.370658428  0.156293132  0.065846721
##  [446]  0.652345253  0.683170807 -0.123854006  0.988927459  0.547377637
##  [451]  0.308240912 -0.069650206  0.551064176  0.878824902  0.614717849
##  [456]  0.849614140  0.207951053  1.153039581 -0.162688220 -0.183918258
##  [461]  0.851858012 -0.461165309  0.433313422 -0.091316759  0.511138872
##  [466]  0.236333375  0.133761251 -0.281184716 -0.274889713  1.006065323
##  [471]  0.744257343 -0.224849578 -0.158462939  0.382052865  0.086263927
##  [476]  0.100002360 -0.036956200  0.900015881  0.491246568 -0.261364759
##  [481] -0.921284278  0.017810083  0.320082786  0.067908927 -0.577639028
##  [486]  0.187639836  0.418834030 -0.115993367  1.146049930  0.098365372
##  [491]  0.036499209  0.411156195  1.807574342  0.659711540  0.714717074
##  [496]  0.845809361  0.082364296  0.493256562  0.289649483  0.179535751
##  [501]  0.534168371 -0.246049531  1.035333316  0.294547593  0.370093637
##  [506]  0.855057072  0.674412029 -0.501552862  0.454355744  0.110305363
##  [511]  0.178199023  0.282076027  0.283903331  0.892080198  0.022712849
##  [516]  0.299547721  0.303950209  0.389076065  0.540349559  0.415092156
##  [521]  0.765070054  0.100845754 -0.247212633 -0.102355300  0.062730170
##  [526]  0.056799975 -0.099309401 -0.042539492  1.123173673 -0.252667565
##  [531] -0.604578749  0.089402409  0.848334435 -0.083815078 -0.010981894
##  [536]  0.212523912  0.050188768  0.373587269 -0.052104119  0.220254931
##  [541]  0.228199243  0.418178359  0.070444716  0.304393368  0.932879830
##  [546]  1.120588424  0.271332793 -0.915929686  0.299926279  1.249145782
##  [551] -0.408130989  0.546085912 -0.302200149  0.287993478  1.107324453
##  [556]  0.219111066  0.062791783  0.516294107 -0.060022929  0.190971158
##  [561]  0.955591822  0.779204778 -0.268488494 -0.322171836  0.199824935
##  [566] -0.473523976  0.348120728  0.423522628 -0.415957478  0.357289286
##  [571]  0.482651738  0.328671475  0.202100433 -0.292212305  0.553742756
##  [576]  0.583563947 -0.543145638  0.291739369  0.513087342  0.639338257
##  [581]  0.651815239  0.332526041  0.570996315  0.796042648  0.594339350
##  [586]  0.526740696  0.108569868  0.625251045  1.598547592  0.918298021
##  [591]  0.692414089  0.841310176  0.573697682  0.867883231 -0.260698418
##  [596]  0.400394332  1.581486100  0.677114568  0.308644998  0.461631302
##  [601]  0.246131477  0.281813897 -0.256185995  0.508419825  0.540996807
##  [606]  0.111868408  0.980233119  0.839032288  0.063434482  0.863737476
##  [611]  0.032869380  1.173707587 -0.030933594  0.270200559  0.272932156
##  [616] -0.663138284 -0.095422345  0.102205405 -0.094780589  0.425914699
##  [621]  0.398402484  0.130653927  0.031737037  0.010599424  0.050669529
##  [626]  0.126031563  0.488071591  0.326320362  0.256664232  0.666188019
##  [631] -0.288930559  0.266895198  0.037082332  0.328786140 -0.074340465
##  [636]  0.044340137  0.104836801  0.076973956  0.155982327  0.455812737
##  [641]  0.069398926  0.693937334  0.489065559 -0.083068334  0.285040257
##  [646]  0.692059904  0.525990840  0.672753679  0.414677456  0.550762418
##  [651]  0.875889948 -0.102743632 -0.293969191  0.160563626  0.605199390
##  [656]  0.145024866 -0.166875071 -0.062060246  0.662880300  0.744936482
##  [661]  0.687443810  0.258146735 -0.032926450 -0.114432737 -0.043712885
##  [666]  1.077202431 -0.102121077  0.089896974  0.182744610  0.122872423
##  [671]  0.662706566  0.364709062  0.163655767 -0.373582031  0.684650450
##  [676]  0.080811201 -0.074526414  0.489035293  0.851369747 -0.136848779
##  [681]  0.465810632  0.550401400  0.225089363  0.251995319  0.313612484
##  [686] -0.063472254  0.740042119  0.650039168  0.435016492 -0.111131044
##  [691] -0.383508075 -0.628492516 -0.052740845  0.073023593  0.041490713
##  [696]  0.279946782 -0.421820993  1.086278316  0.326384537  1.189952806
##  [701]  1.353916724  2.022024117  0.433256727  0.530293877 -0.250409240
##  [706] -0.132395635  0.306319173  0.781787165  0.295860854  0.221993481
##  [711]  0.293933556  0.670654615  0.496120012  0.557434581 -0.074081644
##  [716] -0.216480814  0.429085310  0.157937415 -0.586168653 -0.474395332
##  [721]  0.175975905  0.477831744 -0.200241024  0.028679090 -0.193491821
##  [726]  0.702510235 -0.809984754 -0.309104295  0.203795357  1.277051988
##  [731]  0.128329176  0.159983049  0.126187514  0.203497693  0.138988035
##  [736]  0.179178189  0.727366563  0.062284704  0.772182706 -0.237332148
##  [741]  0.645189831  0.426939023  0.158530814 -0.267805142  0.210170122
##  [746]  1.028502549  0.188320570 -0.033683360  0.149409515 -0.119507009
##  [751]  0.647248107  0.316354156  0.148896999 -0.149632417 -0.203054205
##  [756] -0.006940123  0.325593842  1.407244178  1.010628444  0.218100719
##  [761]  0.627436370  0.486488322  0.201267234  0.620268290 -0.037326631
##  [766]  0.550449293  0.450658750  1.043195205  0.569144732 -0.268142599
##  [771] -0.106998800  0.676117652 -0.082667664 -0.233223315 -0.229888856
##  [776]  0.577347322 -0.135771686  0.873842983  0.211361118  0.185736477
##  [781] -0.252504015 -0.178350872 -0.167232867  0.219111066  0.259066150
##  [786]  1.041754836 -0.241228989  0.336531045  0.263801197  0.184370432
##  [791]  0.161757372  0.228719120  0.714056225  1.094616481  1.166365053
##  [796]  0.414729627  1.013991830 -0.347071105 -0.309576591  0.812870192
##  [801]  0.146912786  0.273280359  1.000122222  0.421670804  0.181768021
##  [806]  0.300028138 -0.045692271 -0.185444744  0.337510931  0.330638107
##  [811]  0.027265545 -0.429702375  0.229510756  0.533778260  0.372047325
##  [816] -0.486658371  0.776432835  0.315233566  1.212761798  0.188177481
##  [821] -0.212718354  0.234981203 -0.501156812 -0.009400819  1.225807866
##  [826]  0.618016188 -0.457576541  0.182524269  0.839439700  0.045751540
##  [831] -0.296876575  0.327974247  0.268334246  0.740766912 -0.252868854
##  [836]  0.402793892  0.542160789  0.610976855  0.293022490 -0.032304146
##  [841] -0.730090812  0.917536973 -0.141615408  0.001748127  1.274886319
##  [846] -0.008592132 -0.102712075  0.098970589  0.431649405  0.434216054
##  [851]  1.150628985  0.542203987  0.366118342 -0.301147796  0.493724940
##  [856]  0.143888994  0.556432221  0.014218573 -0.773301102 -0.037241434
##  [861]  0.063402688  0.480292564  0.875889948  0.108387295 -0.359856388
##  [866] -0.003201746  0.018712741  0.007890949 -0.032521991  0.210170122
##  [871] -1.000421147  0.152909422 -0.463875344  0.120115108  0.797010993
##  [876]  0.198306138 -0.668227171  0.046057292 -0.465545789  0.401669163
##  [881]  0.153606331 -0.128972493 -0.113537416  0.061789082  0.578537531
##  [886]  0.115948568  0.827217109 -0.218919909 -0.883951380 -0.798404568
##  [891]  1.349770283 -0.715914888  1.074821731  0.554355788  0.021719055
##  [896]  0.070004772  0.333554551  0.210981545  1.074429171  0.026313699
##  [901] -0.298671217 -0.038923073  0.540128184  0.159508808  0.303200067
##  [906]  0.403888237 -0.310339294  0.705287873 -0.304027926  0.825917659
##  [911] -0.655145946  0.278102640  0.836886164  0.090601462 -0.171551409
##  [916]  1.180317986  0.167039560 -0.252165420  0.365035417  0.750142737
##  [921] -0.035407380  0.038040510  1.095717310  0.042763156  0.676770030
##  [926]  0.551810891  0.058583404 -0.334074172  0.144302413 -0.573528515
##  [931]  0.033618044  1.102114463  1.205226190 -0.325481631 -0.087845876
##  [936] -0.437586857  1.168010484  0.400226982  0.473038492  0.050924973
##  [941] -0.243468246  0.581190801  0.125175434  0.723974760  0.160755266
##  [946]  0.091103033  0.693082112  0.054739632 -0.122292350  0.714195472
##  [951]  0.469255509  0.526534805  0.846987905  0.234224172  0.591818789
##  [956]  0.502822504  0.155041931 -0.396507426  0.365541986  0.442532402
##  [961]  0.032282670  0.055298286  0.250057317  0.396138479  0.136892614
##  [966] -0.681947005  0.122017943  0.013093831  0.844074066 -0.195367531
##  [971] -0.865982803  0.432881908  0.076953713  0.151275626  0.050785426
##  [976] -0.240101020 -0.490136180  0.023825933  1.216151879  0.873188098
##  [981] -0.510277558  0.275113888  1.132252935 -0.186998750  0.368199127
##  [986]  0.320007473  0.375231483  0.521344670  0.334789964  0.493078404
##  [991]  0.774928424 -0.158462939  0.131501666  0.006938090  0.307744488
##  [996] -0.151535496 -0.032064596  0.568191999  0.659468278 -0.441445974
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
##   0.65669266   0.37906607 
##  (0.11987122) (0.08475928)
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
## [1]  1.19989683  0.08160386 -0.96361691 -0.73060347  0.65310180  0.25283015
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
## [1] -0.0169
```

```r
sd(results2$thetastar)  #the standard deviation of the theta stars is the SE of the statistic (in this case, the mean)
```

```
## [1] 0.9085427
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
## t1*      4.5 0.004104104   0.9334645
```

One of the main advantages to the 'boot' package over the 'bootstrap' package is the nicer formatting of the output.

Going back to our original code, lets see how we could reproduce all of these numbers:


```r
table(sample(x,size=length(x),replace=T))
```

```
## 
## 0 1 2 4 5 6 8 
## 1 2 1 1 1 1 3
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
## [1] 0.0012
```

```r
se.boot
```

```
## [1] 0.9043763
```

Why do our numbers not agree exactly with those of the boot package? This is because our estimates of bias and standard error are just estimates, and they carry with them their own uncertainties. That is one of the reasons we might bother doing jackknife-after-bootstrap.

The 'boot' package has a LOT of functionality. If we have time, we will come back to some of these more complex functions later in the semester as we cover topics like regression and glm.

