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
## 0 1 4 5 6 
## 2 3 2 1 2
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
## [1] -0.0113
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
## [1] 2.675305
```

```r
UL.boot
```

```
## [1] 6.302095
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
##    [1] 5.6 4.7 4.3 3.1 5.0 5.3 5.8 3.5 3.1 3.9 4.4 3.2 4.1 4.9 4.5 5.2 4.0 6.3
##   [19] 4.2 3.4 4.0 3.3 2.0 4.5 3.9 5.0 5.0 4.0 4.4 5.5 4.9 5.5 4.5 3.4 5.1 5.1
##   [37] 4.6 5.2 5.1 4.0 4.6 4.9 4.9 3.9 4.2 7.2 4.4 6.4 4.4 4.7 4.1 5.4 4.2 4.4
##   [55] 5.2 5.7 2.5 4.3 5.1 3.9 5.6 3.9 3.3 5.2 4.1 3.5 4.7 3.7 4.1 5.7 5.5 3.8
##   [73] 6.4 3.5 5.7 5.7 6.0 5.0 3.0 5.4 6.0 4.7 5.0 3.4 5.5 2.8 5.7 3.0 5.3 5.3
##   [91] 4.9 3.5 4.4 3.7 5.9 5.5 3.8 4.1 3.7 2.9 3.6 4.5 4.4 3.3 5.1 4.1 2.9 4.4
##  [109] 5.7 3.9 4.6 4.9 1.7 5.9 3.5 3.0 2.7 2.3 4.6 2.9 5.2 4.4 2.4 5.1 4.6 5.9
##  [127] 5.8 3.9 3.5 4.1 3.5 5.3 5.6 4.4 4.7 4.3 3.6 4.1 3.8 6.7 4.7 4.1 4.1 5.3
##  [145] 2.9 4.4 5.2 5.2 4.7 5.3 5.6 4.0 3.5 3.3 3.8 5.5 3.8 4.2 6.2 3.5 4.0 4.7
##  [163] 5.4 4.5 3.4 4.7 5.0 4.2 6.0 4.7 4.7 6.7 3.8 2.9 5.4 4.4 4.4 5.2 4.3 3.7
##  [181] 4.8 3.8 3.2 5.2 4.0 4.9 5.0 4.2 4.1 4.0 4.2 4.2 4.0 4.4 5.6 4.8 3.5 3.4
##  [199] 5.9 4.0 3.2 4.7 4.6 4.8 4.7 3.9 4.3 6.2 5.6 6.4 3.5 4.4 3.6 4.6 4.3 6.1
##  [217] 3.6 5.2 3.7 4.7 5.0 3.2 4.8 3.6 3.0 4.5 4.6 5.3 4.0 5.4 3.4 3.4 2.4 4.4
##  [235] 5.7 3.2 3.7 4.1 4.2 5.1 4.2 3.7 3.7 4.9 3.5 3.6 5.8 4.4 5.0 5.8 5.0 4.4
##  [253] 4.2 4.1 3.9 4.4 2.2 2.5 4.8 5.4 4.4 3.6 5.6 4.9 3.0 4.5 6.3 5.0 4.5 5.2
##  [271] 5.4 4.8 4.8 3.2 4.1 3.9 4.4 4.4 4.9 4.5 3.9 3.1 4.1 5.0 3.5 4.9 4.5 4.1
##  [289] 4.7 5.5 3.9 4.9 4.1 3.8 4.2 5.6 4.1 3.3 4.9 6.9 3.2 5.5 5.0 3.5 4.5 4.4
##  [307] 3.7 4.9 4.5 4.1 5.6 2.2 5.9 6.0 4.6 4.8 6.3 4.1 5.1 6.0 3.6 4.7 5.5 2.9
##  [325] 5.4 4.8 1.7 3.5 5.4 4.8 3.3 3.0 2.4 4.5 4.5 5.0 5.0 6.0 4.6 4.1 5.2 4.1
##  [343] 5.3 5.3 5.5 5.5 2.4 6.8 3.5 4.8 4.7 3.3 4.0 3.9 5.0 4.5 3.9 4.0 3.2 4.4
##  [361] 6.1 4.9 4.7 4.3 5.5 4.2 3.5 5.2 3.8 6.0 5.2 2.0 4.8 5.1 5.8 4.2 5.2 4.9
##  [379] 5.5 4.4 5.4 4.0 3.9 3.3 4.3 3.4 4.6 4.5 3.8 4.5 3.6 5.5 3.4 4.6 5.7 6.1
##  [397] 2.8 4.4 4.1 4.0 3.1 3.4 4.8 4.3 4.3 4.9 5.4 3.8 4.0 3.8 4.7 3.6 3.6 3.5
##  [415] 3.9 3.7 3.0 4.7 5.1 4.4 4.2 3.6 3.8 2.8 4.1 4.9 4.3 6.2 4.6 4.7 2.6 5.8
##  [433] 4.0 4.5 3.0 4.6 5.2 4.9 4.7 5.1 5.0 3.6 4.1 4.9 3.9 3.9 5.2 4.3 4.2 4.1
##  [451] 3.6 3.9 3.6 3.3 4.9 3.6 5.2 3.8 3.8 4.9 4.6 3.8 3.7 4.0 3.6 4.4 5.5 4.6
##  [469] 5.0 3.9 3.7 5.2 5.3 3.8 5.0 4.2 3.6 5.0 5.0 4.2 4.8 4.4 4.8 4.5 4.1 4.3
##  [487] 5.1 4.0 3.8 4.0 4.9 3.0 3.5 5.1 3.6 4.9 4.3 4.6 5.5 4.5 4.2 4.9 5.2 4.2
##  [505] 6.6 4.0 5.5 4.0 5.3 3.3 4.3 4.8 3.7 3.0 5.2 3.5 5.8 4.4 4.9 4.6 5.3 3.6
##  [523] 3.6 4.1 5.1 4.6 5.3 3.0 5.0 4.4 4.8 5.1 3.1 4.5 5.2 3.8 5.3 4.4 3.0 5.3
##  [541] 5.6 5.0 3.6 4.1 5.0 4.6 4.5 4.0 4.5 4.8 4.4 4.1 4.4 3.2 3.8 3.6 4.0 5.9
##  [559] 4.7 3.9 3.2 4.4 5.5 5.5 5.6 4.8 6.2 4.5 5.0 4.9 4.1 3.6 2.7 3.6 4.2 4.9
##  [577] 5.0 3.2 3.0 3.5 4.6 4.3 4.7 5.9 6.5 6.8 4.7 4.3 4.7 5.2 5.0 5.4 4.5 4.8
##  [595] 4.4 3.2 3.1 4.1 6.3 3.7 4.8 3.1 3.6 4.4 3.8 4.6 6.1 5.1 3.6 4.9 4.4 5.2
##  [613] 2.7 3.5 3.4 3.5 4.5 4.5 2.2 3.2 3.0 4.1 3.8 3.4 3.9 5.2 4.3 3.7 5.1 5.0
##  [631] 6.3 4.4 3.7 3.2 3.4 4.4 5.9 4.3 6.6 5.4 3.1 5.8 5.0 3.8 2.6 4.5 6.0 4.9
##  [649] 4.3 3.2 5.1 4.4 5.5 3.8 5.3 2.1 3.9 4.3 4.8 4.4 4.1 4.5 2.2 3.3 3.9 4.5
##  [667] 5.8 5.4 4.8 4.5 4.2 5.9 2.2 5.2 4.5 5.1 4.9 4.0 4.0 3.7 4.9 4.5 4.6 3.7
##  [685] 3.8 5.4 5.4 3.8 4.2 3.8 4.8 3.9 3.5 5.5 4.4 5.1 5.7 5.4 4.4 2.8 5.3 4.9
##  [703] 5.5 4.5 4.6 4.0 4.6 4.3 4.3 2.6 4.4 5.0 4.6 4.8 5.2 4.5 5.3 3.9 5.0 4.4
##  [721] 4.4 4.2 5.8 4.6 5.5 4.8 5.4 4.4 5.6 4.2 4.5 2.8 4.1 4.2 5.9 5.3 2.6 4.4
##  [739] 4.8 4.5 3.2 4.6 5.9 6.3 6.5 4.4 4.9 5.5 4.6 4.0 5.0 4.9 6.2 5.6 4.2 4.6
##  [757] 5.1 3.4 3.1 4.2 4.5 5.0 4.7 5.0 3.4 5.1 3.6 3.3 5.9 3.7 4.0 5.6 3.6 5.1
##  [775] 3.4 4.0 4.3 4.5 4.1 4.5 5.7 5.1 2.8 4.9 5.5 3.2 2.9 5.5 5.1 4.4 3.0 3.4
##  [793] 4.2 6.4 5.0 4.3 4.7 6.2 3.4 4.5 4.7 4.9 5.4 3.3 4.8 2.9 4.9 5.0 3.7 5.9
##  [811] 4.1 4.3 4.3 4.6 4.4 5.2 5.2 4.9 4.9 4.9 3.4 5.3 5.5 5.0 4.6 3.7 4.7 3.4
##  [829] 3.9 5.1 5.7 4.3 4.9 6.5 4.6 5.9 4.6 4.4 4.7 3.5 4.9 2.7 5.3 5.5 4.5 4.6
##  [847] 4.5 3.4 3.9 3.2 2.9 2.8 3.9 6.8 3.4 5.1 5.7 4.6 4.4 3.9 4.5 6.0 4.4 4.3
##  [865] 5.1 5.5 3.4 6.3 5.1 5.4 4.2 5.0 3.7 4.3 3.2 5.1 5.0 4.6 3.6 3.5 4.1 3.0
##  [883] 5.7 4.3 2.4 4.5 4.6 4.3 4.4 4.4 5.4 5.0 4.4 3.5 5.3 4.7 5.4 5.0 3.6 5.4
##  [901] 4.6 2.4 3.8 5.7 5.7 4.5 4.0 4.9 4.9 3.2 3.5 2.9 3.2 4.6 3.7 4.6 5.6 3.2
##  [919] 5.6 3.8 3.9 3.3 4.3 4.8 3.3 5.1 5.2 5.5 5.4 4.6 3.1 1.9 5.2 4.2 3.0 5.0
##  [937] 5.0 3.5 4.8 5.2 5.9 4.2 5.2 6.2 5.4 2.8 4.1 5.2 3.8 3.4 5.0 3.6 3.8 3.2
##  [955] 2.7 5.1 5.0 3.3 4.6 4.8 4.8 6.6 6.0 5.1 4.4 4.6 4.6 4.1 6.2 4.0 4.2 3.9
##  [973] 4.0 3.3 3.0 5.4 3.2 4.0 4.6 5.7 5.9 4.3 4.6 5.3 4.8 5.1 5.8 3.7 4.5 5.0
##  [991] 2.9 4.1 4.6 2.7 5.2 3.3 3.7 5.8 3.9 5.0
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
##    [1] 4.0 5.7 3.8 4.3 5.0 4.2 4.7 3.8 5.2 3.9 4.0 3.8 2.9 6.2 4.4 5.5 4.4 4.6
##   [19] 5.1 5.0 3.3 5.1 2.5 4.7 5.2 3.7 4.1 3.6 3.7 5.6 3.4 5.5 4.0 4.3 4.6 4.8
##   [37] 3.8 5.2 3.1 4.2 5.4 4.7 4.5 4.2 4.7 4.3 5.0 4.9 4.2 4.3 4.6 4.7 5.5 4.4
##   [55] 3.5 3.2 3.8 2.4 4.8 4.4 5.8 6.2 5.1 2.3 4.2 4.4 4.4 4.8 4.0 6.1 4.0 4.3
##   [73] 4.9 4.1 3.8 4.8 5.1 3.8 4.4 4.9 4.7 5.1 4.7 3.8 5.7 3.8 5.4 4.3 2.8 5.2
##   [91] 3.9 5.1 5.4 5.3 5.5 4.0 4.0 4.2 4.6 3.7 5.9 4.3 4.4 5.4 2.8 5.5 4.9 5.4
##  [109] 2.5 4.6 4.7 4.8 6.4 5.1 3.7 4.0 2.8 3.9 3.9 3.9 5.5 3.7 4.0 6.0 3.4 5.3
##  [127] 4.1 3.9 5.0 3.2 3.3 3.8 4.7 5.2 2.0 5.2 4.0 5.4 4.6 5.4 4.3 4.2 4.4 5.1
##  [145] 3.2 5.2 5.7 5.6 5.7 3.4 4.7 3.2 5.0 3.3 4.4 5.8 6.0 5.3 4.1 3.9 4.7 4.3
##  [163] 3.9 4.6 4.7 4.7 3.3 2.9 3.3 5.4 3.6 5.8 4.7 3.4 4.7 5.4 6.1 4.8 4.7 3.0
##  [181] 3.6 4.2 4.3 4.9 4.8 5.5 4.1 3.2 4.8 4.8 4.8 5.2 4.1 3.1 3.6 4.1 4.6 3.3
##  [199] 4.3 4.5 3.4 5.2 4.9 4.5 4.2 3.4 3.8 4.8 4.4 3.3 4.5 3.3 4.3 5.9 5.5 5.8
##  [217] 6.5 4.7 4.2 4.6 3.0 5.5 4.8 4.9 4.9 4.7 5.6 4.0 3.0 4.0 3.7 5.1 5.0 3.3
##  [235] 3.9 3.3 4.0 4.1 2.8 3.8 3.7 3.4 6.4 4.2 4.4 4.8 4.9 5.6 4.6 4.9 4.5 4.4
##  [253] 2.3 5.3 5.1 4.5 3.8 2.6 3.1 4.2 4.4 5.9 4.2 3.8 7.3 4.9 5.4 3.2 4.0 5.6
##  [271] 4.2 4.3 5.0 5.4 4.2 4.6 3.7 3.2 4.5 3.7 5.3 5.8 4.1 4.6 3.8 6.2 4.0 4.9
##  [289] 4.4 5.7 4.6 4.1 3.6 4.8 5.4 4.5 2.9 5.8 3.0 4.9 5.6 5.6 4.5 3.8 4.9 5.9
##  [307] 3.4 5.2 5.0 4.3 5.2 4.6 3.9 4.5 5.3 4.6 5.7 5.1 4.6 4.8 3.1 4.0 5.3 6.0
##  [325] 3.9 6.7 4.4 2.8 4.6 4.6 4.4 5.9 5.3 4.5 4.0 4.3 6.1 4.2 4.0 5.7 4.0 4.5
##  [343] 5.4 5.5 3.8 5.4 4.3 3.6 4.6 5.4 4.2 5.0 4.0 3.5 4.6 4.8 3.8 4.6 5.3 5.2
##  [361] 4.1 4.2 6.4 3.6 5.1 3.8 4.5 5.7 4.6 3.5 4.5 5.6 4.5 4.8 3.9 3.6 5.0 6.1
##  [379] 3.8 5.0 4.0 5.7 3.3 5.8 2.9 5.0 3.4 6.1 3.6 3.4 5.5 3.0 3.9 6.5 5.1 3.3
##  [397] 6.3 4.7 5.2 4.4 3.0 4.2 5.5 4.8 6.6 3.6 4.7 4.8 4.3 4.6 5.4 3.5 4.5 5.2
##  [415] 5.8 3.9 5.0 4.0 3.3 5.5 5.8 4.5 3.7 3.6 7.7 5.2 3.9 3.3 3.4 3.6 4.8 3.1
##  [433] 5.6 3.0 5.6 3.3 3.7 5.3 5.3 6.3 3.7 5.3 4.8 5.6 3.6 3.3 4.4 5.5 5.9 3.7
##  [451] 5.3 4.7 4.7 4.1 6.1 6.1 5.5 4.8 6.0 4.8 4.7 3.5 4.9 6.6 5.3 5.2 4.5 4.7
##  [469] 3.5 5.0 5.1 4.1 5.3 4.9 3.4 3.8 4.2 4.7 6.6 4.3 4.1 4.9 5.7 6.9 5.3 7.1
##  [487] 4.8 4.0 5.0 5.1 3.5 4.2 6.5 3.6 4.3 3.1 4.6 4.2 4.7 5.0 3.6 4.8 4.1 2.6
##  [505] 6.6 5.7 4.6 3.6 4.3 4.2 4.6 4.6 2.8 5.5 4.0 3.0 4.6 5.3 2.8 4.5 4.9 4.0
##  [523] 5.2 5.2 4.9 5.6 4.9 5.2 5.4 4.8 5.2 3.2 3.2 3.8 4.4 4.3 5.2 5.5 5.2 5.0
##  [541] 4.7 4.1 5.3 5.9 4.6 5.1 3.9 3.9 5.6 4.8 4.8 5.1 5.5 5.3 5.3 4.1 3.0 4.4
##  [559] 4.3 5.7 5.0 2.6 3.3 6.8 4.1 4.8 5.4 5.7 3.9 5.1 3.5 5.8 5.3 3.7 3.8 4.6
##  [577] 6.0 4.8 4.5 3.9 4.7 5.6 4.2 5.2 5.0 3.7 3.8 4.9 3.0 3.6 3.2 5.0 4.7 4.8
##  [595] 4.5 4.9 3.4 4.5 4.4 5.0 6.6 4.1 4.6 5.1 5.6 4.1 3.1 5.5 5.1 3.5 4.8 4.9
##  [613] 5.0 5.7 4.7 4.2 5.2 4.5 4.0 4.5 5.1 5.8 3.9 5.0 5.0 4.8 3.3 5.5 2.8 4.8
##  [631] 3.3 4.4 3.7 4.5 6.0 4.8 3.5 4.6 4.4 5.2 4.4 4.5 4.1 4.5 3.8 3.8 4.0 4.9
##  [649] 5.2 6.0 3.1 4.0 5.3 5.6 4.9 3.5 4.5 4.5 4.5 6.4 4.5 4.7 2.5 4.3 4.5 4.8
##  [667] 3.6 4.8 3.3 4.4 2.8 4.8 4.4 3.0 4.5 4.4 5.0 3.8 4.5 6.3 4.8 3.6 4.7 5.4
##  [685] 2.6 4.8 5.6 3.9 3.9 3.7 4.5 5.8 3.8 2.4 3.8 4.9 2.7 3.5 4.7 4.6 4.2 4.9
##  [703] 3.9 3.7 5.7 3.9 5.3 5.0 3.3 5.4 4.5 5.5 2.9 2.1 6.8 5.6 4.3 4.0 6.3 5.4
##  [721] 4.6 3.9 5.3 2.7 3.5 4.6 4.5 4.2 5.9 5.9 4.5 4.3 3.8 5.8 6.2 4.1 4.2 5.2
##  [739] 4.3 4.5 5.2 5.2 3.2 4.7 4.2 4.6 4.0 5.0 3.9 3.6 4.7 4.1 4.3 3.0 2.6 3.7
##  [757] 4.1 3.5 4.7 2.5 4.6 4.3 4.4 5.2 3.7 6.1 5.0 4.2 3.5 3.8 4.8 3.8 5.6 3.9
##  [775] 4.1 3.8 6.1 5.2 4.7 4.5 5.0 4.0 2.2 4.5 5.1 4.8 4.8 4.9 6.3 4.5 4.4 5.7
##  [793] 5.0 4.6 4.2 4.4 4.8 4.4 5.3 5.4 4.7 3.6 4.9 3.8 3.2 4.4 5.2 6.2 3.6 4.2
##  [811] 4.7 5.6 5.6 2.9 5.2 3.9 4.7 2.5 4.6 4.3 4.2 3.6 6.3 4.7 4.6 3.1 4.7 4.1
##  [829] 4.2 3.2 5.3 5.7 3.4 6.3 6.3 5.3 4.9 4.4 3.1 5.6 2.7 5.1 5.8 4.8 3.2 3.3
##  [847] 4.9 5.8 4.9 4.7 6.4 6.3 5.6 4.9 2.6 6.1 4.5 4.1 5.7 4.2 5.7 3.3 3.4 5.6
##  [865] 4.5 5.4 4.9 5.1 4.3 3.1 6.4 5.5 3.3 5.0 4.9 5.3 4.2 2.8 5.0 5.2 4.9 3.6
##  [883] 5.6 4.5 4.8 4.9 3.9 3.7 6.6 5.4 4.1 6.4 5.4 4.8 2.0 5.7 4.8 5.0 4.4 4.3
##  [901] 3.8 3.6 5.0 4.7 3.7 5.4 3.1 4.5 6.2 3.9 5.4 4.1 5.1 4.4 4.0 4.2 5.0 3.9
##  [919] 4.9 4.4 4.1 6.3 3.9 4.5 4.0 4.4 5.0 3.2 5.0 3.8 3.9 3.6 4.7 5.0 4.5 3.3
##  [937] 4.0 4.7 6.1 4.7 4.0 4.7 5.4 3.2 5.7 5.3 5.2 4.3 4.5 4.6 3.9 3.4 5.3 2.7
##  [955] 4.3 3.6 2.9 4.5 5.4 2.9 5.2 3.6 4.9 3.7 5.1 4.5 4.1 3.7 4.7 4.6 4.1 4.1
##  [973] 5.2 4.1 3.5 5.3 5.2 5.5 3.2 4.1 4.5 4.6 3.8 4.4 4.6 4.7 4.4 5.3 3.7 4.6
##  [991] 4.2 4.2 3.1 3.2 3.7 5.0 5.9 4.5 5.5 4.8
## 
## $func.thetastar
## [1] 0.0284
## 
## $jack.boot.val
##  [1]  0.5048969  0.4949853  0.2833811  0.1435028  0.1415902 -0.1123944
##  [7] -0.1714286 -0.1798271 -0.3048851 -0.4639769
## 
## $jack.boot.se
## [1] 0.9449692
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
##    [1] 5.9 3.6 3.6 3.9 4.3 5.0 3.8 3.8 3.6 4.1 4.7 4.2 5.9 3.9 5.7 3.6 4.1 4.4
##   [19] 5.0 3.5 3.4 3.9 4.5 4.4 3.0 4.2 5.4 6.0 3.7 4.1 5.1 3.9 5.8 3.3 4.2 3.2
##   [37] 3.1 4.1 5.4 5.5 4.5 4.2 5.0 6.5 4.1 4.9 4.9 5.5 3.1 6.4 3.8 5.2 5.0 5.4
##   [55] 4.1 5.3 3.5 4.3 3.7 3.6 4.9 4.8 4.0 5.1 4.0 5.2 7.8 4.9 4.7 5.5 2.3 3.3
##   [73] 5.3 2.6 5.2 4.4 5.4 3.9 4.6 4.5 5.5 5.2 4.7 3.3 4.7 4.4 4.4 4.8 3.8 6.0
##   [91] 4.6 4.1 4.7 3.9 3.5 6.0 4.1 6.2 5.2 5.5 5.2 4.6 4.9 3.6 5.0 4.4 3.4 4.1
##  [109] 4.9 4.7 4.7 4.3 4.3 3.4 5.2 4.5 5.4 4.2 3.0 3.7 4.0 4.3 4.1 2.3 4.9 5.1
##  [127] 5.3 3.8 4.4 4.5 4.8 3.5 3.9 4.2 4.1 3.5 6.2 4.7 3.7 5.3 5.2 3.6 5.7 3.6
##  [145] 3.8 3.7 4.1 4.0 4.7 4.9 5.5 4.0 4.9 5.5 4.2 5.7 4.2 4.6 3.8 5.2 5.3 3.3
##  [163] 5.4 5.0 3.8 4.5 3.4 4.7 4.4 5.1 3.7 3.3 4.5 3.2 5.6 3.5 4.2 5.6 4.5 3.9
##  [181] 5.6 4.6 5.4 5.5 4.1 4.1 3.9 4.1 4.7 4.4 4.0 4.4 3.8 4.7 4.3 5.5 4.3 4.6
##  [199] 3.1 4.1 4.9 3.2 3.1 4.2 4.4 5.3 3.9 4.8 4.0 5.1 4.6 4.1 4.3 3.5 5.4 4.3
##  [217] 3.1 4.5 7.2 6.3 4.4 5.2 5.5 4.8 5.0 5.0 4.8 3.6 5.2 3.5 6.0 4.9 5.1 3.3
##  [235] 4.3 7.5 3.8 5.0 4.6 4.1 4.7 4.8 3.5 5.0 2.9 3.8 4.2 4.4 4.3 2.3 6.0 3.5
##  [253] 4.3 5.3 4.9 4.7 6.3 5.7 4.3 4.1 4.9 3.7 4.8 3.7 5.7 4.1 4.8 5.8 5.1 4.7
##  [271] 6.0 4.9 5.5 4.4 2.3 3.0 4.7 4.7 5.1 3.3 4.6 4.8 3.8 3.8 6.1 3.3 5.0 5.5
##  [289] 4.4 3.7 2.7 4.4 5.0 6.1 4.4 4.4 4.8 3.7 4.0 5.4 3.4 3.5 3.4 4.3 6.6 3.1
##  [307] 4.6 5.1 4.8 5.8 4.3 4.8 4.8 5.6 4.5 4.9 4.1 3.0 4.8 5.1 4.9 3.4 4.1 5.0
##  [325] 4.0 4.5 5.5 4.9 6.3 5.2 3.3 4.1 3.5 4.2 5.0 4.1 5.4 3.7 3.5 4.9 6.0 7.2
##  [343] 3.6 2.7 5.3 5.3 4.6 3.2 5.5 5.4 3.9 5.2 6.1 4.2 6.3 3.7 5.4 4.9 5.8 3.2
##  [361] 3.4 5.2 3.5 4.7 3.3 1.8 4.0 5.2 3.6 3.6 4.5 3.0 3.8 5.8 6.4 5.4 4.6 3.8
##  [379] 4.2 3.6 7.2 5.4 4.4 4.6 5.2 5.5 4.3 3.7 4.8 5.6 4.3 4.3 5.8 4.3 5.7 4.3
##  [397] 4.0 4.0 4.9 5.2 5.9 2.8 3.7 2.9 4.8 3.6 4.4 4.2 3.2 4.9 4.7 4.8 3.0 4.8
##  [415] 6.2 4.0 4.5 5.7 5.0 5.8 6.3 4.5 2.9 4.8 5.3 3.1 4.4 3.5 5.5 3.6 4.1 3.4
##  [433] 5.6 4.4 5.7 4.3 3.8 5.1 4.7 2.9 3.9 5.6 4.2 4.1 3.9 2.8 4.3 5.0 3.8 4.7
##  [451] 4.7 4.2 5.3 3.8 3.5 4.3 4.0 3.2 3.6 3.4 6.4 3.7 5.8 3.3 5.4 3.5 4.6 3.7
##  [469] 4.9 5.7 4.4 3.4 6.1 4.1 4.7 3.3 4.0 3.8 3.9 5.0 2.3 4.1 5.2 6.3 5.7 4.4
##  [487] 5.7 3.6 6.2 5.4 5.4 3.0 4.6 4.0 4.4 4.0 3.8 4.6 4.6 3.7 5.9 4.6 4.1 4.0
##  [505] 6.8 5.3 5.0 4.1 5.6 5.9 5.4 5.0 4.9 4.4 4.3 3.6 4.8 3.6 4.4 6.1 4.1 4.2
##  [523] 4.6 3.3 6.0 4.2 4.1 4.8 3.9 5.0 4.8 6.1 6.4 4.8 4.7 5.5 5.9 3.7 4.2 4.1
##  [541] 4.7 3.1 4.5 3.9 3.9 3.9 5.5 4.8 4.1 4.8 4.2 3.6 5.2 4.6 4.1 4.0 3.3 3.8
##  [559] 5.1 3.5 5.5 3.1 5.0 3.9 4.2 4.9 4.5 3.9 5.1 4.4 5.0 5.7 4.6 4.1 4.3 4.0
##  [577] 4.1 4.4 3.8 4.9 4.1 5.2 4.1 4.8 5.4 3.1 5.0 3.7 4.6 5.4 5.0 6.1 4.4 5.3
##  [595] 3.1 4.9 5.0 5.0 4.5 2.2 2.8 3.4 3.4 5.2 3.9 3.9 3.9 3.8 3.2 4.8 3.3 3.4
##  [613] 4.2 4.6 4.8 3.4 4.7 6.1 3.7 4.4 5.4 4.9 4.1 3.9 5.7 3.6 4.0 3.9 5.2 3.7
##  [631] 3.6 4.5 3.7 4.6 4.7 4.8 5.1 6.8 4.5 6.5 4.2 6.5 5.2 2.9 3.3 5.6 5.7 4.8
##  [649] 5.6 5.0 6.1 4.2 5.4 3.5 2.7 4.4 3.2 4.6 3.6 4.4 4.1 3.3 5.0 5.1 5.4 7.0
##  [667] 4.5 5.0 4.9 4.7 4.7 3.4 4.8 2.8 4.4 4.4 4.2 4.6 5.7 5.7 3.8 4.2 4.8 4.9
##  [685] 5.3 5.3 4.1 3.5 4.6 4.0 5.3 5.3 5.5 5.6 5.6 3.0 3.4 4.9 4.1 3.3 5.1 4.6
##  [703] 4.2 4.9 5.1 4.6 4.8 4.4 3.7 5.1 4.4 3.5 4.9 4.7 4.5 4.7 3.8 3.7 3.9 2.7
##  [721] 4.1 6.2 4.3 4.6 3.6 4.7 4.4 2.6 4.0 5.5 4.4 4.9 5.4 3.6 4.5 5.6 4.6 5.4
##  [739] 4.7 3.5 4.6 6.0 5.7 5.4 5.9 4.3 3.9 6.0 5.1 4.6 3.6 4.7 5.0 3.6 4.3 5.5
##  [757] 5.0 4.2 3.8 5.0 5.4 4.4 4.9 3.1 4.4 4.2 5.2 3.2 5.2 3.8 4.5 5.2 4.4 5.6
##  [775] 3.0 5.4 5.0 5.1 3.9 4.3 3.7 5.3 6.9 4.6 5.6 2.8 4.2 2.5 5.9 4.5 4.2 6.7
##  [793] 6.7 5.0 3.5 4.2 3.5 4.8 5.0 3.9 5.5 4.3 5.1 3.3 3.9 4.9 4.0 4.1 5.8 4.4
##  [811] 5.1 3.4 4.1 3.7 3.9 3.7 4.3 3.4 2.5 5.3 4.0 2.8 5.2 4.4 4.3 4.2 5.0 3.9
##  [829] 4.4 5.9 3.0 3.1 4.0 6.6 5.0 4.3 4.2 4.4 4.7 4.9 4.6 3.8 4.0 6.0 4.3 6.1
##  [847] 4.7 6.1 3.6 4.6 3.1 4.8 6.0 2.5 5.1 4.2 3.5 3.8 4.5 3.6 5.0 3.6 7.5 3.8
##  [865] 5.0 3.7 5.4 4.8 5.5 5.4 5.1 3.9 4.8 4.8 5.1 5.3 3.2 7.5 4.6 5.4 5.8 4.5
##  [883] 4.9 4.8 6.4 4.6 4.1 3.9 4.1 4.7 5.7 4.6 4.5 3.1 2.9 5.5 5.4 5.2 4.3 6.5
##  [901] 5.0 4.6 5.4 5.4 4.6 4.9 4.5 4.5 4.8 3.9 2.5 4.4 5.2 5.5 3.3 6.2 3.1 5.6
##  [919] 4.4 4.3 4.4 4.7 4.0 3.4 5.0 5.3 5.4 3.9 3.7 5.5 3.3 4.4 4.5 3.1 4.6 5.2
##  [937] 3.9 5.7 4.0 6.2 4.4 3.8 4.0 6.6 4.5 4.6 5.9 3.7 4.1 5.6 4.3 3.5 6.2 4.9
##  [955] 4.3 3.6 3.8 5.2 5.2 4.0 4.9 5.1 4.8 5.2 4.0 4.2 4.0 5.7 5.0 3.4 3.8 3.4
##  [973] 4.0 5.0 3.0 3.3 4.2 5.0 5.6 3.9 5.5 3.7 3.9 4.5 4.9 4.9 4.8 5.3 5.1 4.3
##  [991] 5.0 5.4 3.7 4.2 2.8 4.8 4.7 3.2 4.0 5.3
## 
## $func.thetastar
## 72% 
##   5 
## 
## $jack.boot.val
##  [1] 5.4 5.4 5.4 5.2 5.1 5.0 4.9 4.8 4.8 4.5
## 
## $jack.boot.se
## [1] 0.8720665
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
## [1] 0.4073946
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
##      shape       rate   
##   11.359124   18.028595 
##  ( 5.007131) ( 8.125085)
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
## [1] -0.03509314 -0.28290134 -0.60369595  0.36192685 -0.16010030  0.18308065
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
##    [1] -0.3506254577  0.9690395623  0.7965817913  1.1364168035  1.0037144732
##    [6]  0.4391540201 -0.1829434626 -0.2437150167 -0.3277983629 -0.3297570727
##   [11]  1.1384704204 -0.6181303759  0.3228680382  0.2628489573 -0.2967079061
##   [16]  0.2779179810 -1.2634143020  0.9786514775  1.0380635983  0.8786298837
##   [21] -0.8550850137  0.4491504931 -0.3711632386 -0.0487462093 -0.1465372717
##   [26]  0.4041430159  1.5225476659  1.0973125294 -1.2920448021 -0.3102229784
##   [31]  1.1328991726  0.2876485952 -0.1835719023  0.5776616348  0.7990034091
##   [36]  0.2890915175  1.6306747591 -0.0343848521  0.4649222537  0.3087171157
##   [41] -0.2264037872 -0.8568388997 -1.3970908911  1.1251817571  1.0772947478
##   [46]  1.0577160617 -0.5190722441  0.9990067279 -0.0063676818  0.2484300365
##   [51]  0.2269072950  0.5207868775  0.2170453875  0.4110061590  0.7628777289
##   [56]  0.4267475385  0.2242839264  0.7721363436 -0.9894490059  1.1840012208
##   [61]  0.1835474308  0.8543340758  1.0383201382 -0.0257475834  0.5083206484
##   [66]  0.4040303839  0.8325736342  0.2708592560  0.0606805774 -1.1131100189
##   [71] -0.1743238178 -0.4289429085 -1.2430551587  0.2056939715 -0.0152439511
##   [76]  0.9242493891 -0.7755012453  0.8574092753 -0.1522395991  1.1924650643
##   [81] -0.2068583971  0.4875509298  0.9217691085  0.5957716241 -0.2009354045
##   [86]  0.0136283914  1.2763357634  0.7281510325  0.2342268113  0.1210416990
##   [91]  0.1479543445  0.4271541901  0.1783826015  1.5294898459 -0.3842192658
##   [96] -0.0724588291 -0.3593426331  1.2745136234 -1.5348932713  0.5069935722
##  [101] -0.1104423554  0.3412859766  0.9735157163  1.1191581407 -0.0712630096
##  [106]  0.1166042293 -0.2534412758  0.6025539063  1.4962297987  0.0005874198
##  [111]  0.3350594374  0.6948965051  1.2769654692  0.4036876182  0.6699478141
##  [116]  0.1136341767  0.8098514867 -0.4970946815  0.6921956227  0.2247163424
##  [121]  1.0151940503 -0.3368156237  0.6135097411 -0.3947910295  0.5390581312
##  [126]  0.0022556678  0.8223973522 -0.4866020946  0.5889264049  0.7143145721
##  [131] -0.3070521350  0.0029819921  0.6377931452  0.1211150916  0.0397893944
##  [136] -0.1423054342 -0.5183415020  0.7689596995  0.0882534309 -0.1189954659
##  [141]  0.5166479974 -0.1323113076  0.3769986123 -0.2352897174  0.5148844725
##  [146]  0.1569411955  1.1328991726 -0.7321595823  0.0400498115  1.3451565750
##  [151]  0.5384211267  0.1937224186 -0.5502441263  0.7266605109  0.2973401784
##  [156]  0.2389466885  0.3110628483  1.0907244188  0.8357335552 -0.3290114906
##  [161]  0.5042392141  1.0439855281  0.6833279784  0.7910492163  0.7241326025
##  [166]  0.4874122972  0.7254954482  1.2344689172  0.3918688326  0.8771082074
##  [171]  0.1631453994  0.5287169985 -0.4172295178  0.2586310637  0.5984149691
##  [176]  0.7717162941  0.7358661008  0.5853498613 -0.0626603573  0.4135665253
##  [181]  0.5806062003 -0.0316795622 -0.0327306784  0.9292180595  0.4936117254
##  [186] -0.9305480799  0.0124965197  1.0268227243  0.4105587167  0.4285125976
##  [191] -0.7808270254  0.4401173147 -0.3237862098  1.4725923823  0.1070768639
##  [196]  0.8214807084  1.0959227943 -0.4613670752  0.5966833269  0.2004934385
##  [201]  0.5012072669  0.6344163725 -0.2120963772  0.3633576021  0.9743992942
##  [206]  0.3202321682 -1.2363929145  0.1031656480 -0.6645485738 -0.5740896188
##  [211] -0.1258498097  0.5631400599  0.1458261029  0.3708746508 -0.1496601395
##  [216] -0.4795794358 -0.1017105412  0.5971487923  1.8368903028  0.7879621692
##  [221]  0.8591965875  0.2705560698  0.0929400173  0.5881070154 -0.6491913494
##  [226]  0.1439355316  0.0436083651  0.6974972564 -0.5046998822 -0.3343477446
##  [231] -0.5552058352 -0.1992472223  0.1285670625  0.1140509692  0.0795194201
##  [236]  1.1878005500 -0.6378891577  0.4908000855  0.6967900670  0.6070462714
##  [241]  0.5036651257  0.8761507184  1.2502674354 -1.0196653374  0.3382479425
##  [246] -0.5399700077  0.3171433128 -0.1192782987  0.4514772775 -0.2924677197
##  [251]  1.0563328028  0.4161411032 -0.3218180323  0.0282823599 -0.3758872327
##  [256]  1.2262854118  0.2342524928 -0.0821528609 -0.0299282101 -1.0535562274
##  [261] -0.1005861507  0.4911173347  0.6221845291  0.6451670399  0.4928393113
##  [266]  0.2440091948  0.2749867636  0.5847233523  1.1794496257  0.8764595815
##  [271]  0.5036992080  0.6100233416  0.2715628176  0.9467615957  0.4247477697
##  [276]  0.3988417833 -0.3374362475  0.5912261956 -0.2743798582 -0.0620789289
##  [281] -0.2590792947 -0.3360171354  1.2687614079  0.4361916354  0.5912261956
##  [286]  0.7614602503  0.0205838618  1.6426811574  0.6919461808  0.2693143128
##  [291]  0.4344876227  0.4028691682 -0.5203515026 -0.0321235163  0.1373646256
##  [296]  0.6174329908 -0.0694303315  1.1705605511  0.2614440857  0.9346012022
##  [301]  0.4420630643  0.7279702767  0.7864135663 -0.1896676814  0.0883401584
##  [306]  1.0202987672  0.0071421481  2.1345908388  0.5279750527 -1.6152378830
##  [311]  0.6600406845 -0.3591599615 -0.1585446726 -0.0213664433  0.5597454569
##  [316]  0.8274630245 -0.4602236954  0.5384211267  0.5604282788  0.3568374237
##  [321]  0.2551323765  0.2872865425 -0.2085661050  1.4094349528  0.3980763774
##  [326]  0.4349597446  0.8174084243  1.3478995746  0.1365109083  0.5390581312
##  [331]  0.3842119522 -0.0446149959 -0.5017507584 -0.0030180419  0.8965899954
##  [336] -0.1023515696 -0.1161562333 -0.4141268633  1.3016275807  2.1309582076
##  [341]  0.0232285606  0.4362585555  0.3130644464 -0.1753431292  0.8738746629
##  [346]  1.1040815008  0.0207808518  0.0252047029  0.8282146739  1.2116786812
##  [351]  0.3363107730  0.2656971445 -0.0929090103  0.5186158903 -0.4639541845
##  [356] -0.8573960533  1.3222111820  0.3622867763  0.8909220853 -0.3890975154
##  [361]  0.0153358209  0.2673333352  1.3311959258  0.4001804366  1.1260942488
##  [366] -0.3340253953  0.5126793119  0.7608236836  0.6131634661  0.2716046331
##  [371] -0.1134319595  0.3781821931  0.5851273894  0.6282558560  0.4487297376
##  [376]  0.5840871989  0.6853309268 -0.2820221163  0.2577021366  0.6479598169
##  [381]  0.6745043241 -1.1270789064  0.6882928862  1.3272561435  0.5577598704
##  [386] -0.4137276287  0.2892120892  0.7932702404 -0.4136849476  0.4903441736
##  [391]  1.1530974522  0.1058299340  0.2670429707  0.4633333900  0.3271593688
##  [396] -0.1923993823 -0.5180062303 -0.1220170059 -0.3487602452 -0.2388298381
##  [401]  0.6945057740 -1.1298848011 -0.4894871149 -0.0477866947  0.3782274790
##  [406]  0.2084951957  0.7449615022  0.4092186807  0.1638432830  0.0699322621
##  [411]  0.0475780908  0.1488204986 -0.1215283462  0.0101136969  0.7239523453
##  [416]  0.1679682474 -0.4458199937 -0.4343032199 -0.5830474742 -0.0094941667
##  [421]  0.2382244071  0.2567075472  2.2381325889 -0.0579023236  0.2305659736
##  [426] -0.1816536190  0.6407135170  1.2647510238 -0.2292949331  0.1475040710
##  [431] -0.8787158203 -0.7094370230  0.4778090357  0.1267263818  0.9043856922
##  [436] -0.3233471365  0.7934193298 -0.4879468841  0.2919754968  0.3151287394
##  [441]  1.1925742842  0.9679628412  0.8212880984  1.2066647939 -0.3453005470
##  [446] -0.4224227392 -0.4976834832 -0.1555941941  0.0519431483 -0.1154695260
##  [451]  0.9022042818 -0.0985845308  0.4186606549  0.6260627465  0.4072765518
##  [456]  0.4476209713 -0.4802135462  0.6199319019 -0.4053786890  0.4427120556
##  [461] -0.4297317391 -0.8548360516  0.2155016252  0.0393185521  0.6304560895
##  [466]  1.1935848007 -0.1450039976  0.9370526339  0.8020246056  0.1289572477
##  [471] -0.5276451278  0.6905969139  0.2286913387  0.7825737368  0.4968423747
##  [476]  0.3211930842  0.1973207617  0.5098067859  0.2201271440  0.6191912069
##  [481]  0.6659143871  0.6736102970  1.4802395648  0.0144368599  0.7440715578
##  [486]  0.0411858045  0.6283371244  0.5964192548  0.6189971107 -0.2071186397
##  [491]  0.2031369847 -0.2927802928  0.5086902033  0.5021021926  0.6461259230
##  [496] -0.6611294300  0.5439570451  0.5562007275 -0.5394273630  0.0633128420
##  [501] -0.0278314571  0.7106780721  0.9707670018 -0.0261707087  0.8502521581
##  [506] -0.4487555205  0.8408229210  0.4076348728 -0.5721148583  0.4058401243
##  [511]  0.0509757761  1.0040452218 -0.1459520530 -0.4519306713  0.5000649594
##  [516]  0.3114002863  0.1458899687  0.0698469817  0.4784081466  0.1409926822
##  [521]  0.3968366040  0.2622440408 -1.2402374843  0.6990353597  1.0245855065
##  [526]  0.1954015850 -0.5824569512  0.0185166999 -0.1118671951  0.6127361403
##  [531]  0.7255365270 -0.3136366217  0.2515759810  0.6047038999  0.7629276179
##  [536]  0.2204183019  1.0374651584  0.4546675213  0.0389211935 -0.1522713090
##  [541]  0.1659781342 -0.1972674446  1.3734455976  0.6333053332 -0.4520090811
##  [546]  0.9918107847  0.3317592066  0.1792564687  0.5203384549  0.7380788141
##  [551]  0.6245230050 -0.0560594270 -0.1033839685  0.6400467648  0.3756281055
##  [556]  2.0436316766  1.1967790512  0.3114288827  0.0488470233  0.0936704730
##  [561]  1.3190433387 -1.5608786377 -0.3773628175  1.1553884563  0.3441877873
##  [566]  1.3483564393 -0.1096115537 -0.4143702103 -0.0894568823  0.3513829153
##  [571]  0.9953925930  0.5822503119  0.5700542649  1.5044483119  0.0337381565
##  [576] -0.5100417410  0.8588596320  1.1829540944 -0.0641620258  0.5491744123
##  [581]  0.1557243611  0.1409926822  0.1477366652 -0.3814411217  0.1361577004
##  [586]  0.7805001201 -0.9618109684 -1.3248065888  0.2815941732 -0.2268689582
##  [591]  0.1698668765  0.4868367981  0.9627964126  0.2444983896  0.7793121949
##  [596]  0.4388680253  0.3322473263  0.5570064039 -0.1086449573 -0.0492177555
##  [601]  0.0914272936  0.3287789762 -0.4726471764  0.5599821662  1.0802188141
##  [606] -0.1168858837 -1.1582303305  0.2903232710  0.3026988491 -0.1522323511
##  [611]  0.8639357755 -0.9388778266 -0.0646262416 -0.4708047371  0.5381805210
##  [616]  0.2827708476  0.6220084371  2.0500669685 -0.6597742872  0.8436808404
##  [621]  1.0114666803  0.0128058542  0.4965799935  0.3407561419 -0.6212656433
##  [626]  0.1997639522  0.7520283703  0.9417396896  0.0721108731  0.1012948235
##  [631] -0.3379083586  0.3435674821  0.3654062076  0.4610776611  0.3664958092
##  [636]  0.8285879933  0.4359992015  0.4899324903  0.5094406880  1.4802395648
##  [641]  0.9407315100 -0.3339163905 -0.8127727914  0.5779188152 -0.1438035405
##  [646]  0.0138967848  0.1930117325  0.6448051916  0.7751118603  1.2543616517
##  [651] -0.1781315512  0.7727685145  1.3783932233 -0.0766725882 -0.3700044935
##  [656] -0.6181303759  0.0435569946  0.3184674182 -0.7765702164  0.1263466673
##  [661] -0.0636421567  0.4900977568 -1.6672118653  1.9736458120  0.3958173242
##  [666]  0.4029957649 -0.3273385463  0.3437190274  0.5681599569 -0.0307196897
##  [671]  0.0613340419  0.0635741336 -0.2896745601 -0.5224617359 -0.0115827849
##  [676]  0.8708457866 -0.4770132637  0.8976254277 -0.7051541250 -0.6702170921
##  [681] -0.2514490606 -0.6865827033  0.5088964766 -0.0723970784 -0.6840248606
##  [686]  0.3282625227  0.4215816953  0.7296451701  0.1838735418 -0.1009680121
##  [691] -1.3046157855  0.8105776030 -0.0691813619  0.8975799671  0.6807517870
##  [696]  0.2655878620  0.3060632088 -0.4325974826 -0.6062316063 -0.1140175141
##  [701]  0.7712805025  0.3701483137  0.4470797723  1.6734268496 -0.0541116585
##  [706] -0.3270751081 -0.4033123008  0.5325889855  1.4263949559  0.0791902666
##  [711]  0.3147481405  0.5443220834  0.1423908810  0.3527155629 -0.2596386896
##  [716]  0.0076079012  1.2090913229 -0.1947025888  0.3183321102 -0.1252906699
##  [721]  0.8998787227  1.1307142619 -0.3117239940  0.0970645623  0.6709741138
##  [726]  0.0477260022  0.6263348894  0.3386882328  0.5931425576  0.2745905680
##  [731] -0.0360448606 -0.1355990964 -0.1214068417  0.3889523360  0.7811084232
##  [736]  1.1027951205  1.3696277486  0.0401322523  0.0783208152  0.7747325168
##  [741]  0.5383531547 -0.0595059763  0.2612733288  0.0421617178  0.4256350003
##  [746]  0.1408388142 -0.1722084833 -0.5362921832  1.0891628282  1.0274545957
##  [751]  0.0804574898  1.1580980528 -1.2230734469  0.0468719188  0.1528172417
##  [756]  1.1896738756 -0.0315026513  0.0810387947  0.6986433958  0.4479098030
##  [761] -1.5103865053  0.1813431443 -0.3170361462  0.3288254253  0.0695709103
##  [766] -0.1017930563  0.2930030917  0.0787054284 -0.5616596800  0.2710127715
##  [771] -0.0461516100  1.3907370413  0.1402672433 -0.6994507811  1.1184056871
##  [776]  0.2430138704  0.1790981038  1.0828451233  0.3283392725 -0.0846941486
##  [781]  0.0622492603  0.5770162035  0.0207573244  1.1776372984  0.1556094316
##  [786] -1.0340120662  1.5355146375  0.0230881506 -0.3076424836  0.8140652294
##  [791]  0.2884943321 -0.6758340587  0.5228190958 -0.5141345679 -0.3451108021
##  [796]  0.2206054629  1.8085487497 -0.0942124159  0.5823799185  0.9042573229
##  [801]  0.0821747914 -0.2914085933  0.8121696642  0.1142806722  0.1376480336
##  [806]  0.6102443253  0.5395392874 -0.1876246868 -1.0694221433 -0.0086426577
##  [811] -0.0380540557  0.7425146984  1.0331112167  0.5287169985  0.5434096850
##  [816]  0.3389151898  0.6771039235  0.4888528277  0.4379364036  1.7683146484
##  [821]  1.2086729257  0.0434581130  0.2080106070  0.0565070977  0.7547039929
##  [826] -0.1061059188  0.4989737921 -0.1122444012  0.2996212369  0.8156662064
##  [831]  0.3329254628  0.9197240613 -0.4618594195  0.5468661591  0.4641900712
##  [836]  0.4617862756  0.2789813601  0.8342315646 -0.6349178481 -0.3626035899
##  [841]  0.4519917031 -0.7252401111 -0.1384295137 -0.0210166285 -0.8245296148
##  [846]  0.8907832869  0.6912463797  0.6706265410 -0.3165485616  0.0972341660
##  [851]  0.1000960330  0.4085987835  0.6850881912 -0.1519514886 -0.3939958384
##  [856]  0.6329158386 -0.6175494800  0.4394972211  0.2980215905  0.9690277886
##  [861]  0.8748448127  1.1260942488 -0.7269847434  1.0199657701  1.1163624291
##  [866] -0.4277703105 -0.3334597680  0.6407135170  0.2203111830  0.6847567549
##  [871]  0.5804993257  0.2250604716 -0.8609372423  0.8562375532  0.6471482845
##  [876] -0.6394003812  0.5600588127 -0.7011534337 -0.3947189314  1.9883329302
##  [881]  0.0618974712  0.0015217984  0.6372775600 -0.4308206983  0.2410792115
##  [886] -0.0467840002  1.3647340132  2.1777413819  0.2608418645  0.3679610838
##  [891]  0.2719987347  0.1444923710  0.5018375117 -0.5027769763  0.3005384658
##  [896]  0.0803127532 -0.0994435939 -0.3126443938 -0.0142961683 -0.0046969515
##  [901] -0.7306808593  0.1978800834  0.2382672490 -0.0442682059  0.1528172417
##  [906]  1.7518534934 -0.5912349126 -0.1121784618 -0.7464555128  0.7406907767
##  [911]  0.9155385542  0.4885490458  0.9273571213  1.3041163355  1.0422239626
##  [916] -0.0613661041  0.8988441622  0.0877144522  1.7739846997 -0.5456732794
##  [921]  0.0584948956  1.3227088473 -0.2491262789  0.3067036855  1.0936476743
##  [926] -0.0076631661  0.0132112277  0.6026435686 -1.0330013255 -0.7446841840
##  [931] -0.7479483275  0.2410088733  0.7317116846  0.2656382041  1.5810684319
##  [936]  0.8288401558  0.5326727430 -0.1176094506  0.1611307070 -0.0474514038
##  [941]  0.7770438994  0.8838604237  0.0368328914  0.3367389011  0.1973860353
##  [946]  0.0107970864 -0.6959811233  0.7403637054 -0.0321559500  0.3406748324
##  [951] -0.7069267699 -0.6645485738 -0.0105119296 -0.1618916906  0.4342758330
##  [956]  0.8059774678 -0.0574060545  0.0089794998  0.8342340671 -0.2663105550
##  [961]  0.5742610893  0.0337322017  0.1350939521  0.8967223603 -0.3500851907
##  [966]  1.1126961381  0.5239894513 -0.8344426761  0.3926022914 -0.0133409012
##  [971]  0.3649754257 -0.1968034472 -0.0173757693  0.0435390137  0.5576910018
##  [976] -0.2931553163  1.1519695297  0.7900099811 -0.1487047510 -0.3266040738
##  [981]  0.5070994815  1.4463740222  0.0763916297 -1.2741834534 -0.0740210730
##  [986]  0.5228993283 -0.1120724244  0.1493089005 -0.4915075954  0.6313590807
##  [991] -0.4576867371  0.9094388187  0.3381333157 -0.1406917190 -0.5181313059
##  [996] -0.6563575523  0.1945707497  0.2194275963  0.4309327333  0.8781765155
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
##   0.63005352   0.18450959 
##  (0.05834705) (0.04125255)
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
## [1]  0.06328287  0.18342966 -0.60645073  0.17776459  0.03615058 -0.05754394
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
## [1] -0.0277
```

```r
sd(results2$thetastar)  #the standard deviation of the theta stars is the SE of the statistic (in this case, the mean)
```

```
## [1] 0.9480198
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
## t1*      4.5 0.001101101   0.8973624
```

One of the main advantages to the 'boot' package over the 'bootstrap' package is the nicer formatting of the output.

Going back to our original code, lets see how we could reproduce all of these numbers:


```r
table(sample(x,size=length(x),replace=T))
```

```
## 
## 0 2 3 4 6 9 
## 1 2 1 2 2 2
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
## [1] 0.0158
```

```r
se.boot
```

```
## [1] 0.8922031
```

Why do our numbers not agree exactly with those of the boot package? This is because our estimates of bias and standard error are just estimates, and they carry with them their own uncertainties. That is one of the reasons we might bother doing jackknife-after-bootstrap.

The 'boot' package has a LOT of functionality. If we have time, we will come back to some of these more complex functions later in the semester as we cover topics like regression and glm.

