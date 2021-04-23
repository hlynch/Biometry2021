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
## 0 1 2 3 5 6 7 9 
## 2 1 1 1 1 1 2 1
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
## [1] 0.021
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
## [1] 2.75662
```

```r
UL.boot
```

```
## [1] 6.28538
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
##    [1] 4.7 4.1 4.8 3.8 5.0 4.8 3.7 4.1 3.1 3.7 3.8 4.9 4.9 3.7 5.0 4.1 3.3 4.3
##   [19] 3.6 4.4 5.2 6.3 5.2 4.3 5.4 4.0 4.4 5.7 5.6 5.0 5.9 4.5 6.0 3.4 5.8 4.6
##   [37] 4.0 5.2 5.1 5.1 4.9 4.2 3.6 5.0 6.7 6.7 4.0 4.0 4.1 3.7 4.5 5.5 6.1 4.8
##   [55] 3.8 5.3 4.7 4.6 4.4 3.3 5.7 5.6 4.5 4.0 5.9 4.7 5.7 4.4 5.2 4.1 4.8 3.8
##   [73] 4.3 6.1 4.2 5.9 4.9 4.2 4.8 3.8 3.4 5.7 3.2 3.6 5.0 4.5 5.6 5.7 4.9 4.6
##   [91] 3.8 5.2 4.0 2.6 4.4 4.6 4.8 3.3 5.9 6.4 4.0 4.8 5.8 4.3 5.0 3.3 4.0 4.6
##  [109] 3.5 4.5 4.9 3.8 4.5 4.0 4.9 4.0 3.2 4.8 4.5 2.4 4.6 5.4 3.2 4.6 6.1 4.0
##  [127] 5.2 5.1 4.6 4.3 5.1 5.5 4.5 4.2 3.9 3.7 5.2 4.0 3.7 6.5 4.6 3.6 5.1 4.6
##  [145] 3.5 4.8 4.9 5.3 4.6 4.6 6.0 6.7 4.8 5.8 3.1 5.5 5.4 2.6 4.0 5.4 4.1 4.4
##  [163] 3.2 4.8 4.0 4.6 4.6 4.7 4.7 3.2 4.2 4.4 4.5 4.4 4.7 4.3 3.2 4.2 3.9 5.4
##  [181] 4.6 5.5 6.2 3.4 5.1 5.3 2.9 4.6 5.3 4.7 5.3 5.1 4.2 4.0 4.7 5.0 6.1 5.1
##  [199] 4.4 4.2 4.1 2.9 3.8 2.8 5.0 4.4 4.0 4.5 4.6 5.9 4.0 5.5 3.8 4.4 5.3 5.4
##  [217] 2.8 5.2 5.8 4.6 4.3 3.9 5.9 5.1 5.6 3.7 4.1 3.5 4.1 5.7 5.6 3.2 6.5 5.3
##  [235] 4.6 4.9 4.7 3.0 5.3 1.6 5.1 4.1 4.2 4.4 4.1 4.7 4.5 4.2 2.9 4.3 6.1 3.1
##  [253] 4.9 3.0 5.7 4.1 4.8 1.9 4.2 4.9 4.8 4.4 4.2 3.2 4.9 3.4 2.3 6.1 2.4 3.9
##  [271] 5.7 5.9 4.0 5.9 4.8 2.9 4.7 3.9 4.6 4.7 4.9 5.8 4.3 5.2 6.2 6.2 4.0 4.3
##  [289] 4.0 4.2 4.3 5.2 4.3 5.4 5.0 4.9 5.0 4.7 4.0 4.3 4.1 4.3 4.0 3.9 3.1 5.6
##  [307] 4.9 6.0 4.9 3.5 4.5 4.0 4.7 4.9 4.6 4.9 4.5 4.9 5.9 4.9 3.5 5.0 4.2 3.3
##  [325] 4.3 4.1 4.3 4.3 4.0 5.5 3.8 4.5 4.7 4.6 5.0 6.3 3.0 5.0 4.6 4.3 4.4 6.5
##  [343] 4.4 3.9 3.3 4.8 4.5 5.9 4.1 6.1 3.5 3.2 5.3 5.9 4.5 4.1 4.7 4.5 4.5 3.4
##  [361] 3.2 4.3 3.0 3.9 5.7 3.7 4.2 4.7 3.5 2.7 5.7 3.4 5.0 4.3 4.4 3.5 4.5 2.7
##  [379] 3.3 4.9 4.4 3.4 3.3 2.6 3.7 2.8 3.6 3.4 4.1 4.0 3.9 5.0 3.9 5.7 4.9 4.7
##  [397] 4.8 3.6 4.7 6.0 5.5 3.3 3.4 4.4 6.9 4.2 3.8 3.1 6.0 3.6 4.7 6.0 3.6 2.8
##  [415] 5.8 4.7 3.0 4.5 5.6 5.0 3.9 5.5 4.9 4.8 5.5 5.6 4.8 4.4 3.9 3.8 4.0 3.9
##  [433] 3.4 4.5 3.4 4.6 3.1 5.3 4.9 5.9 3.7 3.7 4.4 6.1 3.4 4.1 4.8 3.7 5.3 4.5
##  [451] 4.6 2.8 4.6 3.3 5.5 5.7 5.7 4.7 4.2 4.3 5.2 5.8 5.9 3.3 5.6 5.0 5.5 6.5
##  [469] 3.3 5.1 4.7 4.8 5.9 5.2 4.8 4.4 3.2 5.5 6.4 3.9 3.6 5.4 4.7 5.6 4.6 3.8
##  [487] 5.0 4.2 3.1 3.7 5.0 4.5 5.5 4.7 4.3 4.9 2.9 5.0 4.7 5.2 5.3 5.6 4.1 5.1
##  [505] 3.7 5.8 4.5 3.2 4.6 3.7 4.4 5.7 4.4 4.3 4.5 3.7 2.5 4.7 5.1 3.2 6.1 5.3
##  [523] 5.0 3.7 2.9 4.5 2.6 5.6 4.1 5.6 6.1 4.6 3.6 3.9 3.7 4.6 5.5 5.5 4.6 3.8
##  [541] 3.8 3.4 4.3 4.5 6.1 3.4 5.1 4.6 3.1 4.2 5.3 2.4 4.7 5.7 5.3 4.1 4.2 4.0
##  [559] 3.9 4.4 3.6 5.2 4.6 5.4 5.4 6.2 3.6 3.5 4.8 5.3 4.4 4.0 4.7 4.8 4.8 6.2
##  [577] 4.9 5.3 3.9 3.8 4.3 5.6 4.1 6.1 3.0 4.2 4.6 3.4 6.3 4.2 3.3 3.5 3.6 4.4
##  [595] 3.7 3.2 4.9 4.1 5.0 2.9 4.7 5.0 4.7 3.5 4.5 6.2 5.0 5.1 3.7 4.4 4.3 5.5
##  [613] 4.7 3.6 3.4 4.1 4.8 4.9 5.2 4.2 4.6 4.1 4.8 3.4 5.1 5.1 3.0 4.0 4.6 3.7
##  [631] 3.9 5.3 4.7 5.3 4.2 5.5 4.7 4.6 4.3 5.3 3.6 4.5 2.7 4.1 3.8 6.1 3.9 4.5
##  [649] 5.0 4.7 3.0 4.7 4.7 3.6 3.4 4.3 4.8 4.1 4.7 4.8 4.5 3.9 4.5 5.5 4.6 4.1
##  [667] 5.2 4.4 4.5 4.6 4.7 4.3 6.1 5.5 5.5 2.7 3.7 5.0 4.4 4.2 3.8 3.5 3.5 5.0
##  [685] 4.6 5.0 3.4 4.1 5.3 6.4 5.4 5.9 5.2 3.1 3.9 3.8 4.5 4.0 4.0 4.1 3.2 4.3
##  [703] 3.6 3.9 5.7 5.7 4.0 4.0 3.5 3.8 4.3 4.3 5.0 4.8 3.8 3.7 4.7 5.4 5.1 5.1
##  [721] 5.5 4.0 5.4 5.1 5.1 5.2 3.5 4.9 2.8 4.3 4.9 3.2 5.9 4.5 4.8 5.6 4.8 5.0
##  [739] 6.6 5.1 3.3 4.8 4.1 6.5 4.2 5.5 3.7 6.0 2.6 4.4 4.5 5.7 6.1 3.5 6.1 3.5
##  [757] 3.5 4.7 3.8 4.6 4.9 5.6 6.0 4.3 4.0 4.1 4.1 3.8 3.6 3.5 3.4 3.0 4.4 4.8
##  [775] 4.2 5.0 5.6 2.8 4.2 5.5 3.6 3.8 6.5 5.1 4.5 3.7 5.5 4.9 5.2 5.4 4.3 4.0
##  [793] 4.3 3.9 3.5 5.2 4.5 5.8 4.6 4.0 4.8 3.8 3.7 5.2 5.5 4.7 4.3 3.3 3.9 2.4
##  [811] 3.1 3.3 3.6 4.2 4.6 5.5 6.0 5.0 5.6 4.8 3.0 3.9 6.3 3.6 4.7 4.1 5.8 2.8
##  [829] 4.4 4.2 3.7 4.1 3.4 5.0 3.8 6.1 4.2 6.8 4.2 4.8 3.6 4.7 3.9 3.1 4.4 5.3
##  [847] 4.4 4.5 5.1 4.8 4.7 5.1 5.2 4.7 4.0 4.1 5.1 5.0 4.4 5.3 5.8 4.9 4.7 5.9
##  [865] 3.8 4.1 5.2 5.4 2.7 4.9 4.6 4.1 4.8 3.0 4.7 5.1 4.5 5.7 4.2 5.2 4.3 4.8
##  [883] 4.6 4.5 4.0 5.0 3.9 4.4 5.7 4.2 4.4 4.8 3.5 5.3 5.5 4.2 5.5 4.9 5.8 5.6
##  [901] 3.8 4.1 4.2 4.5 4.0 4.7 4.3 5.5 5.4 4.2 6.0 3.8 5.6 5.9 4.4 5.8 4.7 3.9
##  [919] 5.8 5.0 5.0 4.8 4.6 4.3 5.5 3.9 4.6 5.6 4.0 4.2 4.4 5.1 5.4 5.9 4.8 3.9
##  [937] 3.6 3.7 5.2 3.2 5.6 4.5 5.4 4.9 5.0 5.7 5.7 3.4 3.7 5.8 4.5 5.2 3.8 3.7
##  [955] 2.9 5.4 3.1 6.0 4.2 5.8 5.8 3.6 5.0 5.2 4.5 4.4 4.6 4.3 4.2 3.0 5.2 4.8
##  [973] 3.7 4.4 5.2 4.9 4.0 5.7 6.4 4.2 5.7 5.6 4.6 3.0 4.2 5.4 3.7 4.3 5.7 5.0
##  [991] 5.1 4.7 5.5 5.0 4.3 5.2 6.6 5.2 3.5 4.8
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
##   2.8   6.2
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
##    [1] 5.2 5.2 5.6 5.2 3.3 4.4 4.7 4.9 4.1 4.3 4.4 3.5 4.4 5.4 6.7 5.6 6.3 4.3
##   [19] 4.2 4.0 2.6 4.6 6.0 4.3 3.2 4.8 3.4 5.5 3.2 3.6 5.1 5.2 4.3 3.8 5.1 3.8
##   [37] 5.3 6.7 3.9 3.5 3.1 5.0 4.5 4.5 5.7 3.1 4.0 3.8 4.6 4.2 5.5 4.8 3.4 4.5
##   [55] 4.5 5.2 6.4 4.8 3.3 5.5 3.3 4.6 4.4 5.0 5.2 4.3 4.7 3.3 4.8 2.9 5.1 4.7
##   [73] 5.5 5.4 4.2 4.6 5.8 4.7 4.4 4.9 5.8 4.9 5.4 4.4 4.4 5.1 5.2 3.1 3.8 4.2
##   [91] 4.9 5.5 5.2 4.8 4.4 5.7 4.9 3.8 4.2 4.1 4.3 4.6 4.1 3.3 5.2 4.0 4.2 4.1
##  [109] 4.7 4.8 4.4 2.9 4.3 3.0 4.8 5.4 2.9 3.8 6.2 4.4 3.6 4.6 3.6 4.9 5.1 4.0
##  [127] 2.5 2.9 5.6 4.0 3.8 4.6 5.4 5.3 4.7 4.3 3.9 5.4 4.4 5.2 5.0 4.3 4.3 6.6
##  [145] 5.2 4.1 5.4 5.7 4.5 5.8 2.9 3.6 5.7 4.3 4.6 3.8 5.3 3.5 5.0 4.4 5.6 5.1
##  [163] 4.7 4.2 5.7 3.7 4.8 5.4 3.2 2.5 3.8 4.8 6.5 5.3 4.8 4.0 6.7 2.6 4.7 3.9
##  [181] 2.9 4.8 4.7 5.1 5.2 6.0 3.8 3.6 4.9 4.8 4.9 4.3 3.6 3.2 5.4 4.2 3.9 3.2
##  [199] 5.5 4.8 2.3 3.6 5.2 5.4 5.5 4.1 3.1 3.1 4.1 3.7 3.7 4.2 4.3 4.1 4.3 4.5
##  [217] 5.3 5.6 3.4 3.7 5.2 5.2 5.5 6.0 4.0 5.5 4.4 3.8 4.0 4.0 3.5 3.4 3.8 4.5
##  [235] 3.8 4.0 3.9 2.7 3.7 4.7 3.8 4.7 5.1 5.0 6.7 5.6 6.2 2.9 5.0 4.6 4.8 5.6
##  [253] 4.8 4.0 3.6 5.6 3.8 2.9 3.7 5.5 3.5 5.3 5.8 3.0 6.0 4.4 2.4 4.6 5.5 5.5
##  [271] 3.0 4.7 4.5 5.3 4.3 4.1 5.3 5.5 4.8 6.3 5.4 3.9 4.2 4.9 2.6 5.0 2.7 5.3
##  [289] 5.0 3.4 3.7 4.3 4.6 6.0 4.7 3.4 4.9 3.2 6.0 5.1 4.0 5.2 4.7 3.8 4.7 4.5
##  [307] 4.5 3.3 5.5 5.5 5.3 2.5 5.2 5.5 4.4 4.5 3.5 3.6 4.7 5.1 5.4 4.7 4.2 3.6
##  [325] 3.7 6.4 3.6 4.8 4.5 4.5 4.4 6.2 3.9 3.4 3.5 3.9 3.9 5.4 5.2 3.5 3.0 5.4
##  [343] 3.5 5.7 5.0 5.2 3.4 3.9 4.3 3.8 4.7 5.8 4.2 4.7 4.1 2.7 4.6 4.8 3.2 4.6
##  [361] 3.3 5.2 5.1 5.0 6.3 4.2 5.1 4.2 4.1 4.1 3.3 4.5 5.4 2.5 5.0 4.3 4.4 3.9
##  [379] 3.9 4.7 3.7 5.2 5.7 6.2 5.1 4.1 4.9 4.0 4.2 4.1 4.8 5.6 5.4 5.2 4.0 4.1
##  [397] 4.7 5.9 3.5 3.6 4.7 4.2 2.8 6.3 5.3 4.6 4.3 3.8 5.7 5.2 4.4 3.9 4.4 5.5
##  [415] 4.0 4.0 4.8 3.0 4.1 3.1 2.5 4.9 5.3 5.8 3.8 5.5 3.5 3.8 4.9 4.5 4.8 3.2
##  [433] 4.7 3.3 5.2 4.8 2.9 4.8 2.8 4.8 3.6 3.3 6.0 5.0 4.6 4.6 3.3 3.2 4.4 4.9
##  [451] 4.0 3.5 5.5 4.1 4.0 5.1 2.7 4.0 4.4 4.9 6.1 4.5 4.6 4.7 3.6 2.9 4.8 4.4
##  [469] 3.9 2.7 5.5 3.1 3.9 4.3 3.9 4.3 4.8 5.4 4.5 3.9 3.8 5.6 6.0 5.0 4.0 5.0
##  [487] 4.8 3.5 4.6 5.6 4.6 4.8 4.4 4.5 4.3 4.8 5.1 3.8 5.1 4.5 3.1 4.6 5.0 5.0
##  [505] 5.8 6.2 5.5 5.1 4.9 3.7 5.4 5.3 4.6 5.0 5.7 3.2 3.6 4.2 4.5 5.8 4.1 4.9
##  [523] 3.7 3.0 3.5 4.0 5.4 6.2 4.6 4.1 5.1 4.1 4.6 4.7 5.3 4.9 4.5 4.1 3.3 3.4
##  [541] 6.5 5.5 4.1 5.2 4.0 3.1 5.9 4.6 5.5 6.4 5.2 4.0 4.9 3.9 3.8 4.9 4.1 4.2
##  [559] 3.6 5.0 3.8 4.1 4.5 7.0 4.3 3.6 4.5 4.6 6.8 4.4 5.5 4.7 4.0 5.2 6.0 6.1
##  [577] 4.0 5.9 5.9 5.1 3.5 4.7 5.0 3.7 5.6 5.4 6.5 4.7 4.3 4.7 4.1 5.0 4.1 4.2
##  [595] 4.1 3.4 3.0 4.9 3.1 4.7 4.7 5.0 4.2 4.0 4.0 5.2 3.9 4.7 2.9 5.8 3.1 3.5
##  [613] 4.2 3.9 5.0 4.0 3.9 2.5 6.3 3.7 4.5 4.0 3.1 5.3 3.6 4.3 4.0 3.4 4.9 3.8
##  [631] 4.3 3.7 3.7 4.3 4.2 6.3 3.7 3.1 3.5 3.7 5.0 3.0 4.2 3.5 5.3 5.7 4.4 2.6
##  [649] 4.9 4.0 3.5 5.0 5.0 2.1 3.6 4.6 5.8 5.3 4.7 4.8 4.9 3.5 4.6 5.5 4.1 6.1
##  [667] 4.3 4.4 5.6 4.8 4.6 4.5 4.5 4.1 4.9 4.8 4.2 4.1 3.1 4.4 4.5 5.9 4.1 5.6
##  [685] 5.7 4.7 5.9 4.8 3.1 3.6 5.4 3.9 5.1 5.3 5.5 4.6 4.3 5.2 2.9 5.2 4.4 5.5
##  [703] 2.7 5.3 4.4 5.3 4.1 4.3 4.0 5.8 3.5 6.3 3.4 3.3 5.9 2.5 3.9 5.3 5.9 4.6
##  [721] 4.0 6.6 5.6 4.7 3.3 5.2 5.1 5.5 3.8 5.1 4.8 6.3 1.8 5.5 3.2 4.8 7.1 4.5
##  [739] 4.0 4.2 4.0 4.9 5.3 2.8 3.1 4.8 4.2 3.3 4.2 2.9 4.2 3.9 3.4 5.0 5.3 5.1
##  [757] 4.6 5.3 4.9 4.5 5.2 5.3 3.7 3.3 4.6 5.2 2.7 5.7 5.1 5.2 3.4 5.7 2.5 4.5
##  [775] 4.1 3.1 4.6 4.0 4.6 3.2 4.1 3.8 3.9 4.6 3.6 4.0 4.6 3.8 6.4 2.9 4.9 5.6
##  [793] 5.1 5.8 3.4 5.6 4.8 4.1 5.2 6.7 5.0 4.8 5.1 5.6 4.4 4.2 4.7 3.9 5.2 3.6
##  [811] 3.9 5.3 4.0 5.8 5.4 4.9 4.8 5.5 4.5 3.9 4.2 3.5 5.1 4.6 4.3 5.2 3.7 4.8
##  [829] 3.0 4.1 5.3 5.0 4.5 4.4 3.7 3.7 5.3 4.4 5.1 5.6 4.5 5.2 4.2 5.8 4.1 4.6
##  [847] 3.2 5.0 5.1 4.6 4.0 2.9 4.0 3.8 6.8 4.7 5.1 3.5 3.3 3.2 5.0 4.7 4.8 3.4
##  [865] 3.3 3.2 6.5 3.3 3.9 5.0 4.5 3.8 3.0 4.1 3.5 4.5 3.5 5.0 5.7 4.5 4.4 4.3
##  [883] 3.6 4.5 5.8 5.9 4.7 5.2 3.0 4.6 5.9 6.5 4.0 4.8 5.0 5.5 5.3 2.7 4.0 6.0
##  [901] 5.5 5.3 3.7 5.7 4.3 4.1 4.0 4.4 5.1 5.1 5.5 6.2 4.0 3.0 5.5 4.4 5.0 5.7
##  [919] 3.7 3.7 5.2 4.0 4.5 5.6 5.4 4.6 5.4 5.4 6.2 5.7 4.5 3.2 4.3 3.9 3.3 6.1
##  [937] 4.1 3.9 3.9 5.4 3.8 4.8 5.6 3.5 3.8 2.8 5.8 4.1 2.8 5.4 3.7 5.3 4.1 2.7
##  [955] 4.6 5.1 3.1 2.7 4.1 4.4 4.2 5.0 3.2 5.3 4.3 3.9 4.5 4.6 3.7 6.0 4.8 4.0
##  [973] 5.9 3.1 5.4 6.0 6.3 5.3 4.8 5.3 4.5 5.2 3.5 5.9 4.7 5.1 4.3 3.2 4.1 4.4
##  [991] 5.5 5.3 5.0 5.1 3.8 4.2 3.3 6.3 5.3 3.9
## 
## $func.thetastar
## [1] 0.006
## 
## $jack.boot.val
##  [1]  0.507530120  0.434972678  0.255752212  0.207826087  0.020343840
##  [6]  0.002162162 -0.139117647 -0.300000000 -0.368555241 -0.487106017
## 
## $jack.boot.se
## [1] 0.9658304
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
##    [1] 4.8 5.7 3.8 4.9 4.2 4.2 4.3 4.8 4.8 4.3 4.5 2.4 5.6 4.4 4.4 5.4 4.7 4.4
##   [19] 5.0 4.2 3.8 5.1 3.5 5.4 4.1 5.6 5.7 4.9 2.7 5.6 5.0 4.1 5.1 4.4 4.5 4.3
##   [37] 4.2 5.2 3.2 4.1 3.5 5.6 5.0 3.1 6.1 3.4 5.4 5.7 3.6 3.2 3.1 4.8 4.1 6.0
##   [55] 4.7 4.2 4.6 4.0 5.1 4.0 5.9 6.4 3.6 2.8 3.1 5.1 6.2 3.6 3.7 4.4 4.3 3.9
##   [73] 5.0 5.9 3.9 3.0 4.1 4.4 4.7 3.7 4.5 4.0 6.9 5.0 3.6 3.6 4.6 4.6 4.6 5.1
##   [91] 3.8 4.3 4.4 3.7 5.3 4.1 5.1 3.6 4.2 5.4 5.5 5.3 5.1 4.6 3.5 5.7 4.2 2.7
##  [109] 5.7 4.6 4.3 5.8 4.2 4.8 5.0 5.0 4.2 4.1 4.5 4.4 5.7 5.7 6.0 5.2 4.5 5.3
##  [127] 5.8 4.8 4.1 3.0 5.9 3.1 4.5 4.6 3.5 5.7 4.7 5.9 6.2 4.0 5.0 4.4 5.1 4.9
##  [145] 3.3 5.1 5.0 4.8 4.8 5.5 3.3 3.4 4.5 4.3 4.6 4.2 4.2 3.7 4.8 4.8 3.8 3.9
##  [163] 4.1 5.5 4.1 4.7 6.3 4.6 3.6 3.8 4.6 4.6 4.9 4.0 4.5 5.5 4.0 4.3 5.4 4.2
##  [181] 4.6 5.5 4.0 4.0 4.3 2.4 3.6 5.1 4.9 4.6 5.5 4.8 4.5 5.9 5.6 5.9 5.3 4.0
##  [199] 5.2 4.6 3.6 5.2 4.3 3.7 4.4 5.7 5.0 3.1 3.7 5.6 5.1 3.7 4.7 5.1 3.8 3.4
##  [217] 5.5 5.6 2.8 4.0 4.5 5.5 3.4 4.1 4.1 6.0 4.7 4.6 6.3 4.4 4.3 3.9 4.3 6.0
##  [235] 4.0 3.6 4.5 4.5 4.5 3.9 2.8 5.2 3.9 2.5 5.6 2.7 4.0 4.4 3.5 5.3 2.6 3.7
##  [253] 6.0 3.9 5.3 5.0 4.2 5.5 6.0 6.0 4.5 6.5 3.3 3.8 3.7 5.1 5.7 3.6 3.8 4.1
##  [271] 5.5 5.1 3.7 3.4 5.6 4.5 2.7 5.0 6.8 4.4 4.1 4.3 4.9 4.0 4.0 4.9 4.2 4.8
##  [289] 4.0 4.1 4.4 4.8 4.1 4.7 5.1 4.1 3.6 4.3 4.7 5.0 3.2 5.6 3.5 3.4 3.8 5.9
##  [307] 4.4 5.3 3.7 4.7 2.5 4.8 5.0 4.1 4.0 4.1 5.0 5.0 2.7 5.0 5.4 5.4 4.4 6.2
##  [325] 4.2 3.6 2.6 4.6 4.7 3.8 3.7 3.1 4.2 3.6 3.5 3.7 3.7 5.0 4.3 3.1 3.8 3.6
##  [343] 4.3 5.6 4.0 4.2 3.7 4.6 5.4 5.0 5.3 5.0 4.5 3.6 4.6 3.9 4.2 3.4 3.7 5.4
##  [361] 5.1 4.5 6.1 3.1 4.0 3.8 4.5 5.0 5.6 3.7 3.4 5.3 5.8 5.4 5.4 4.8 4.7 4.3
##  [379] 4.1 5.6 3.7 3.9 4.3 4.3 4.8 5.0 6.1 3.7 3.0 4.5 4.2 3.0 6.5 6.2 4.9 4.9
##  [397] 5.1 4.6 4.1 4.2 3.9 5.1 4.2 4.0 6.2 3.3 3.7 5.3 5.0 3.3 5.5 2.6 5.5 4.0
##  [415] 6.2 5.9 6.4 3.9 3.3 4.5 4.6 4.3 4.6 5.4 3.8 5.4 4.6 4.5 4.6 3.0 5.9 3.9
##  [433] 4.3 2.7 5.1 4.0 3.4 5.2 4.9 4.0 4.2 3.7 5.2 5.9 4.0 5.3 4.5 5.6 3.8 4.6
##  [451] 5.5 5.1 4.8 5.3 5.3 5.8 4.8 3.6 4.0 5.1 6.7 3.6 5.2 4.4 5.9 5.1 4.3 2.5
##  [469] 2.4 5.0 4.6 3.7 4.6 4.7 4.0 5.7 4.4 4.0 4.2 4.8 4.5 4.9 5.1 5.5 6.0 5.4
##  [487] 5.2 5.5 4.8 6.2 4.3 4.1 4.0 6.9 6.7 3.9 3.9 4.9 5.2 4.4 4.5 3.7 2.8 3.6
##  [505] 3.9 4.3 4.5 4.9 4.0 4.4 2.3 4.5 3.9 3.7 4.6 3.4 4.5 4.8 5.8 3.2 4.5 3.9
##  [523] 4.2 6.3 4.9 4.7 3.0 3.3 4.7 3.3 2.9 4.3 4.6 4.2 4.4 6.4 4.3 4.7 4.4 4.7
##  [541] 3.3 4.9 5.6 6.2 4.0 4.3 5.0 3.9 3.9 5.0 2.6 5.1 2.4 3.7 5.3 5.0 5.7 4.1
##  [559] 4.5 3.5 3.7 3.8 3.5 5.2 4.7 5.5 4.8 4.8 6.2 3.5 2.6 4.9 4.8 4.6 5.3 3.8
##  [577] 5.2 4.9 3.1 4.5 3.2 3.9 5.5 6.1 5.9 6.2 3.3 4.4 4.2 4.0 4.6 4.1 4.6 5.1
##  [595] 4.6 3.5 5.6 5.1 3.3 4.2 5.9 4.1 4.8 4.0 6.3 4.4 4.9 5.4 5.1 3.5 5.1 3.6
##  [613] 4.6 4.4 3.8 5.6 3.6 3.8 5.1 3.8 6.0 5.3 3.6 6.4 5.5 4.1 4.0 4.5 4.0 4.5
##  [631] 5.3 4.9 3.1 4.0 4.8 3.9 3.6 5.9 4.2 3.1 3.8 3.9 3.5 5.9 4.6 4.6 3.2 3.0
##  [649] 4.3 4.5 3.3 4.1 3.7 6.0 5.5 5.0 4.3 4.2 5.3 4.1 4.7 3.6 4.5 3.7 5.9 2.8
##  [667] 4.8 4.0 5.1 3.6 3.7 4.6 4.1 3.5 4.1 3.9 4.4 3.3 5.7 2.3 4.7 5.2 5.2 7.4
##  [685] 5.3 5.0 6.6 4.7 4.5 3.9 5.5 4.7 4.6 4.8 5.9 4.5 4.3 4.9 4.8 3.6 5.8 4.3
##  [703] 4.3 5.3 5.2 3.4 5.8 5.0 4.5 5.4 5.6 4.4 5.6 3.1 4.6 4.2 6.0 5.4 5.3 4.9
##  [721] 3.8 4.6 5.4 3.1 2.4 5.1 4.1 4.5 4.4 4.1 3.5 5.9 3.8 4.0 4.1 3.9 4.6 4.0
##  [739] 5.3 4.3 4.4 3.6 3.8 3.4 5.4 3.2 3.5 5.0 5.0 3.8 5.6 6.2 3.0 4.1 5.5 2.8
##  [757] 3.8 4.0 3.4 4.9 4.7 5.2 3.8 6.3 5.9 4.6 4.7 5.5 3.8 4.6 6.0 3.7 3.0 4.6
##  [775] 5.1 5.4 4.0 5.0 4.6 2.5 3.8 4.3 5.3 6.1 4.4 5.5 4.4 5.1 5.3 3.9 4.1 4.7
##  [793] 3.7 3.4 4.9 3.2 4.0 3.1 3.7 4.7 3.8 5.1 2.2 4.3 2.9 5.7 4.3 3.0 6.5 5.9
##  [811] 5.2 5.6 4.2 4.0 2.9 3.3 2.9 3.4 3.7 5.8 5.3 2.8 3.2 4.7 4.0 4.1 4.7 3.9
##  [829] 4.2 4.1 4.2 5.6 5.1 3.4 4.0 5.4 4.9 3.3 4.3 3.7 4.7 6.0 4.4 4.7 4.3 4.4
##  [847] 4.8 4.8 4.8 4.7 3.9 5.9 4.7 4.7 2.6 3.4 4.7 2.9 6.9 7.8 5.5 3.7 3.3 3.0
##  [865] 5.2 5.2 4.7 6.4 5.6 3.6 6.8 4.2 3.2 3.1 5.2 4.2 3.0 3.0 3.6 2.8 3.4 4.2
##  [883] 4.0 4.4 5.1 5.2 5.9 3.9 4.5 5.9 4.6 4.0 4.2 4.1 4.5 4.3 4.1 5.3 3.0 5.7
##  [901] 4.7 5.9 4.2 4.9 4.2 5.1 3.8 4.7 6.1 4.9 6.5 4.4 5.5 3.0 4.7 5.4 2.8 3.0
##  [919] 4.7 4.7 3.4 2.8 6.2 3.6 4.6 3.1 4.1 3.5 6.9 3.0 4.2 3.8 4.2 2.6 5.2 4.9
##  [937] 5.7 3.7 3.4 2.5 4.1 7.3 4.3 6.0 4.2 4.9 4.3 4.2 6.3 4.4 3.2 5.4 4.2 5.0
##  [955] 3.9 5.2 4.6 5.7 5.7 3.9 5.0 1.9 4.0 4.6 3.6 4.2 4.7 2.9 5.5 3.6 4.8 6.1
##  [973] 4.2 3.8 4.9 4.4 5.8 4.1 4.2 5.4 3.7 5.0 3.6 5.6 4.8 4.7 2.6 4.1 6.1 3.2
##  [991] 3.4 5.9 2.8 4.3 4.5 4.2 5.2 2.6 4.8 4.5
## 
## $func.thetastar
## 72% 
##   5 
## 
## $jack.boot.val
##  [1] 5.600 5.500 5.400 5.400 5.100 5.000 4.800 4.764 4.500 4.400
## 
## $jack.boot.se
## [1] 1.209028
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
## [1] 1.600079
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
##   2.1170741   4.0512627 
##  (0.8823811) (1.9043139)
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
## [1]  1.2319765  1.4002592 -0.3828810  0.1684133  1.0566434  0.9545307
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
##    [1] -7.085676e-01  1.466920e+00  6.916472e-01  5.916579e-01  6.510162e-01
##    [6]  1.915284e+00  1.537394e+00  1.844980e+00  3.227182e-01  1.698078e+00
##   [11]  1.292882e+00 -1.295609e-01  1.612871e+00  8.149594e-02  1.547856e+00
##   [16]  1.698078e+00  1.701716e+00  3.347039e-01  4.679723e-02  1.044190e+00
##   [21]  1.954229e-01  1.019115e+00  4.001746e-02  1.763317e+00  9.929086e-01
##   [26] -6.791566e-01  1.037045e+00  1.487446e+00  2.066526e+00  2.068412e+00
##   [31] -4.706087e-01  1.715242e+00  1.893024e+00  5.212656e-02  4.628001e-01
##   [36]  1.025281e+00  1.719910e+00  1.726905e+00  1.850149e+00  9.123330e-01
##   [41]  9.136422e-01  2.467015e-01  6.509780e-01  1.171388e+00  1.749455e+00
##   [46]  1.579457e+00  2.179907e+00 -7.612147e-01  1.846079e+00  6.825760e-01
##   [51]  1.895720e+00 -3.517989e-01 -7.600955e-01  1.566038e+00 -5.901391e-02
##   [56]  1.676391e+00  1.460882e+00  1.530287e+00  6.375178e-01 -3.123720e-01
##   [61]  1.501395e+00  9.620036e-01  1.053815e+00 -1.695857e-01  5.652549e-01
##   [66]  2.120916e+00 -5.093361e-01 -2.183935e-01  1.038353e+00 -1.134645e-01
##   [71] -1.365257e+00  1.684830e+00  1.460882e+00 -9.625889e-02 -3.208774e-01
##   [76]  1.049432e+00  1.615014e+00  1.844724e+00  1.897747e+00  8.371808e-01
##   [81]  1.333708e+00  1.898745e-01  1.158483e+00  1.814963e+00  1.067934e+00
##   [86]  3.583690e-01  1.622505e-01  1.590701e+00  7.774163e-02  1.754317e+00
##   [91]  1.973907e+00  7.030878e-01 -2.620927e-01  1.646837e+00  1.650983e+00
##   [96]  1.113094e+00  1.163132e+00  1.082903e+00  9.315826e-01 -1.913937e-02
##  [101]  1.796259e-01  1.681132e+00  1.147134e+00  6.380552e-01 -2.465988e-01
##  [106]  1.397935e+00  1.637103e+00 -2.450586e-01  9.104782e-01  1.083512e+00
##  [111]  1.184649e+00  1.465588e+00  5.517552e-02  8.371809e-01  1.059695e+00
##  [116]  8.447266e-01 -9.708669e-01  1.211832e+00  1.846843e+00  3.793670e-01
##  [121]  2.161556e+00 -8.475357e-01  5.440548e-01  1.269987e+00  1.236476e+00
##  [126] -1.923827e-01  3.785579e-01  1.208887e+00  1.606456e+00 -2.430573e-01
##  [131]  6.737377e-01  1.932641e+00  1.444950e+00 -6.452466e-04 -2.704895e-01
##  [136]  1.113065e+00  1.933845e+00  1.594534e+00  1.310276e-01  6.843866e-02
##  [141]  1.725167e+00  4.022221e-01 -2.347023e-01  1.400074e+00  1.087321e+00
##  [146]  1.675094e+00  1.747440e+00  8.974930e-01  5.501534e-01  1.722758e-01
##  [151]  1.084511e+00  2.043011e+00  1.770480e+00  1.302825e+00  2.241968e-01
##  [156]  1.717184e+00  1.414781e+00  1.934141e-01  4.712910e-01  9.152603e-01
##  [161]  1.715376e+00 -3.063096e-01  7.347511e-01  4.283067e-01 -1.381043e+00
##  [166] -8.330135e-01  1.469253e+00 -2.022313e-02  1.609760e+00  1.423693e+00
##  [171]  1.559588e+00  1.799310e-01  2.119013e-01  1.011389e+00 -4.133121e-02
##  [176]  9.270913e-01  3.762521e-01  1.605977e-01 -1.842027e-01 -1.331607e+00
##  [181]  6.535764e-01  1.429967e+00  2.077257e+00  9.105871e-03  1.298036e+00
##  [186]  1.151131e+00  2.148154e+00  1.702365e+00 -2.003539e-02  2.048939e-01
##  [191]  1.753421e+00 -7.415866e-01  1.321446e+00  1.253989e+00  1.501194e+00
##  [196]  1.311712e+00  1.049194e+00  6.721010e-01  5.049139e-05  1.920130e+00
##  [201]  5.694229e-02  4.652327e-01  1.704817e+00  8.273031e-02  1.349232e+00
##  [206]  1.347355e+00  1.262907e+00  1.787329e+00  1.737800e+00  1.859939e+00
##  [211]  1.438798e+00  1.834813e+00  1.004449e+00  6.565735e-01  1.512511e+00
##  [216]  1.617644e+00  1.391659e+00  1.663685e+00  6.630938e-01  1.935067e+00
##  [221]  1.027292e+00  2.029083e+00  1.061110e+00  1.692911e+00  5.418951e-01
##  [226]  1.280234e+00  9.481922e-01  4.321606e-01  3.179494e-01  1.183707e+00
##  [231]  1.832558e+00  1.830910e+00  7.202938e-01  1.566707e+00  5.032630e-01
##  [236]  8.361143e-01  2.010524e+00  1.072352e-01  4.896830e-01  1.009059e+00
##  [241]  1.815343e+00  1.632799e+00  9.116628e-01  1.068322e+00  1.498959e+00
##  [246]  1.098557e+00  1.406009e+00 -3.151670e-03 -1.589701e+00  1.923544e+00
##  [251]  1.794671e+00  1.808969e+00  1.719681e+00  1.213921e+00  1.858834e+00
##  [256] -2.793273e-01  8.966351e-01  1.227222e+00  1.007948e+00  8.583456e-01
##  [261]  5.319415e-02  3.815235e-01 -1.059818e+00  6.682733e-02  1.866051e+00
##  [266]  1.039012e+00  1.094886e+00 -3.775124e-01  1.217969e-01 -2.037776e-01
##  [271]  1.116600e-01  1.294578e+00  8.161973e-01 -2.338258e-01  1.535198e+00
##  [276]  1.578614e+00  1.482552e+00  1.910873e+00  9.480922e-01 -4.211334e-01
##  [281]  4.216686e-01  9.174691e-01  3.579738e-01  1.140470e+00  1.845821e+00
##  [286]  2.481719e-02 -5.701979e-01  1.032470e+00  1.729273e+00  1.313812e+00
##  [291]  1.314783e-01  1.972472e+00  7.689671e-01 -5.766739e-01  4.798237e-01
##  [296]  1.801513e+00  1.127915e-01  1.079745e+00  1.094746e-01  1.576635e+00
##  [301]  5.629498e-01  8.957980e-01  1.584073e+00 -1.203829e+00  1.097298e+00
##  [306] -6.197476e-01  1.812376e+00  1.729740e+00  1.628478e+00  1.634660e+00
##  [311]  1.938467e-03  1.128341e+00  2.071410e+00  1.716246e-02  1.710768e+00
##  [316]  1.696410e+00  1.328644e+00  1.088930e+00  1.715370e+00  1.803252e+00
##  [321]  3.670039e-01  4.662850e-01  5.977569e-01  5.854306e-01  7.936137e-01
##  [326] -2.275375e-01  5.395859e-01 -3.766589e-01  1.054814e+00  1.099842e+00
##  [331] -1.046341e+00  1.172112e+00  7.592217e-01  1.298544e+00  1.343932e+00
##  [336] -6.337182e-01  2.242386e+00  1.612871e+00  9.488186e-01  5.433320e-01
##  [341]  1.825622e+00  1.298265e+00  9.721588e-01  1.995476e+00  1.686943e+00
##  [346] -2.893129e-02  1.658269e+00  1.340220e+00  6.044605e-02  1.724790e+00
##  [351]  1.932266e+00  1.380108e+00  7.085057e-01  1.715489e+00 -3.082543e-01
##  [356]  1.340605e+00  1.045301e+00  1.632373e+00  1.182948e+00  4.393500e-01
##  [361] -1.765677e-01  4.209825e-02  1.629277e+00  1.036924e+00  1.593163e+00
##  [366]  1.782656e+00  1.175618e+00 -1.126889e+00  1.136131e+00  1.639290e+00
##  [371] -2.244997e-01  8.572140e-01  1.768386e+00  1.506791e+00  3.391226e-01
##  [376]  5.838897e-01 -1.619532e-01  1.023462e+00 -1.716589e-01  1.778530e+00
##  [381]  2.097551e-01  1.793577e+00  3.617911e-01  9.926745e-01  1.338543e+00
##  [386] -1.136141e-01  1.104485e+00  1.205514e+00  2.055582e-01  8.528588e-01
##  [391]  1.975674e+00  1.139225e+00  1.154341e+00  1.118243e+00 -9.854628e-01
##  [396] -1.566827e-02  1.814450e+00  1.730486e+00 -2.147144e-02  7.106522e-01
##  [401]  1.679930e+00  1.648206e+00  1.609799e-01 -2.700610e-01  2.022284e+00
##  [406]  8.137062e-03  6.912827e-01 -4.008236e-01  1.033951e+00  2.296229e-01
##  [411]  1.085156e+00  1.501624e+00 -1.312295e-01  9.342254e-01  4.718229e-01
##  [416]  1.164185e+00  1.173195e+00  1.635693e-01  1.941433e+00  1.617193e+00
##  [421]  1.969756e+00  1.395185e+00 -1.929224e-01 -1.459230e+00 -1.183297e-01
##  [426]  4.231550e-01 -2.506619e-01  2.108578e-01  8.468125e-01  1.861532e+00
##  [431]  1.005187e+00  5.933694e-01  6.207141e-01 -2.998301e-01 -8.238737e-01
##  [436]  4.359055e-01  1.448818e+00  1.887496e+00  9.875585e-01  1.964678e+00
##  [441]  9.292759e-01  1.137520e+00  7.876108e-01  1.701069e+00  8.595213e-01
##  [446]  1.608687e-01  1.286061e+00  7.283958e-01  1.948705e+00  1.125129e+00
##  [451] -2.336000e-01  1.984284e+00  1.951138e+00  1.179910e+00 -2.171299e-01
##  [456]  1.508663e+00 -3.824739e-01  1.756316e+00 -5.823967e-01  1.039051e+00
##  [461]  7.982125e-02  2.010690e+00  1.143430e+00  4.008286e-01 -8.463844e-02
##  [466]  1.064124e+00  2.183154e-01  1.193205e+00  1.766346e+00  1.579019e+00
##  [471]  1.181669e+00  1.153318e+00 -9.367796e-01  1.019307e+00  1.586397e+00
##  [476]  1.894067e-02  1.610097e+00 -1.012916e-01  1.144711e+00  1.540221e+00
##  [481]  1.751370e+00  9.887605e-01  1.863052e+00  1.135102e+00  1.064792e+00
##  [486]  1.198113e+00  1.198113e+00  8.554079e-01  6.601041e-01  1.884394e+00
##  [491]  1.904625e+00  8.271742e-01  1.609034e+00  1.157904e+00  6.080717e-01
##  [496]  1.264647e-02  1.047513e+00  1.407803e+00  1.821043e+00 -1.605832e-02
##  [501]  1.621674e+00  8.271033e-01  2.343899e-01  5.291154e-01  1.635186e+00
##  [506]  1.812451e+00  1.654635e+00  3.168526e-01  8.484132e-01  7.092693e-01
##  [511]  1.645011e+00  1.660613e+00  1.738614e+00  1.409732e-01 -2.647192e-01
##  [516]  1.448179e+00 -5.430765e-01  1.803208e+00  8.475497e-01  1.434988e+00
##  [521]  1.341528e+00  1.610548e+00  4.326865e-02  1.983436e+00 -7.979960e-01
##  [526]  1.615781e+00  4.804180e-01  1.181085e+00 -3.299779e-01  3.819565e-01
##  [531]  1.812640e+00  2.656240e-01  5.285137e-01 -3.761259e-01  8.162372e-01
##  [536] -6.432862e-01  2.658027e-01  1.819524e+00  1.562377e+00  1.662481e+00
##  [541]  1.718970e+00  1.602921e+00  1.027058e+00  3.704573e-01  8.236345e-02
##  [546]  8.586915e-01  1.833590e+00  1.879305e+00  6.113639e-01  1.379495e+00
##  [551]  2.171980e+00 -5.884699e-01  1.386414e+00  2.128605e+00 -2.067009e-01
##  [556] -5.985487e-01  1.143648e+00  1.385447e+00  1.052089e+00  4.116521e-01
##  [561]  1.728900e+00 -7.068893e-01  8.325954e-01  6.095399e-01  1.536111e+00
##  [566]  1.751210e+00  1.289271e+00 -3.356287e-01 -8.027305e-01  2.171980e+00
##  [571]  8.082044e-01  5.399731e-01  1.631045e+00  4.880228e-01  1.318860e-01
##  [576]  1.167947e+00  1.841630e+00  1.636244e+00  1.656566e+00  8.973473e-01
##  [581]  1.569254e+00  7.387838e-02  8.046983e-02  8.755464e-01 -8.555449e-01
##  [586]  8.255129e-01  1.279645e+00  9.633748e-01  9.687847e-01  1.078371e+00
##  [591]  9.811667e-01  1.400122e+00  1.142312e+00  1.239681e+00  6.093658e-01
##  [596] -1.114703e+00  1.647221e+00  1.486777e+00  7.774506e-01  1.746929e+00
##  [601]  2.129455e+00  4.234974e-01  1.477624e-02  7.909046e-01  1.065689e+00
##  [606]  1.704672e-01  1.717136e+00  3.643891e-01  1.007381e+00  4.558496e-01
##  [611]  4.803167e-02  6.591338e-01 -6.491565e-01  5.838214e-01  1.910873e+00
##  [616] -4.147871e-02  8.674626e-01  1.453855e+00  1.787401e+00  4.833893e-01
##  [621]  1.015935e+00 -7.681547e-01  1.854080e+00  1.042275e+00  1.782181e+00
##  [626]  1.077147e+00  7.708620e-03  1.510717e+00  1.021552e+00  1.888058e+00
##  [631] -5.825354e-01  1.613075e+00 -4.669810e-01  1.009146e+00  1.213603e+00
##  [636]  9.593068e-01  9.571360e-01  4.393665e-01  4.556794e-02 -2.070693e-01
##  [641]  1.349991e+00  7.857192e-02  1.793191e+00  7.414526e-02  1.417403e+00
##  [646]  1.706262e+00  1.832690e+00 -1.013572e+00  3.021465e-01  1.366642e+00
##  [651]  1.045149e+00  1.136803e+00  1.300668e+00  8.494922e-01 -1.637254e-01
##  [656]  8.331148e-01  9.350176e-01  1.205275e+00  4.760986e-01  1.516366e+00
##  [661]  1.027874e+00  1.399609e+00  2.505853e-02  1.089235e+00  1.164277e+00
##  [666]  1.787510e+00 -1.933234e-01  2.201733e+00 -7.906074e-01 -1.843465e-01
##  [671]  4.228314e-01  4.516804e-01  1.455778e+00 -2.866921e-01  6.531776e-01
##  [676]  1.340447e+00  2.022514e-01  1.097473e+00  9.871079e-01  1.151761e+00
##  [681]  5.416380e-01  2.077663e+00  1.477965e+00  1.392460e-01  1.190890e+00
##  [686]  1.672074e+00  1.416261e-01  1.617788e+00  1.503190e+00  1.180288e+00
##  [691] -2.238674e-01  6.063308e-01 -2.532478e-01  1.603412e+00  9.013557e-01
##  [696] -2.139239e-01  1.455688e+00  1.409710e+00  1.457031e+00  4.305888e-01
##  [701]  3.962955e-01  1.980425e+00  1.612871e+00  1.077149e+00  5.966048e-01
##  [706] -1.439517e+00  1.088970e+00  1.409260e+00  1.765631e+00  1.744125e+00
##  [711]  1.538164e+00  1.235146e+00 -2.972198e-01  4.065907e-01  5.278156e-01
##  [716]  2.709900e-01  1.243849e+00  1.327609e+00  1.740826e+00 -4.604341e-01
##  [721]  1.281538e+00  1.565049e+00  1.192966e+00  9.727569e-02 -2.034867e-01
##  [726]  1.130144e+00  1.415417e+00 -5.650651e-01  1.810221e+00  1.314783e-01
##  [731]  1.681331e+00  1.348679e+00  1.540061e-01  6.938214e-01  1.481930e+00
##  [736]  1.313597e+00  1.755423e+00  2.263114e+00 -3.329419e-01  1.004743e+00
##  [741]  1.975983e+00  1.323755e+00  7.086281e-01  1.679758e+00  1.633077e+00
##  [746] -6.714224e-03  1.646444e+00  1.123350e+00  6.987094e-01 -1.764507e-01
##  [751]  1.555322e-01  2.065711e+00  1.283221e+00  9.194466e-01  4.113788e-01
##  [756] -1.023597e-01  4.633389e-02  1.759464e+00  1.794448e+00  1.885060e+00
##  [761] -5.417032e-01  1.037571e+00 -1.059418e+00 -6.088251e-01  1.555975e+00
##  [766]  1.601620e+00  3.779003e-02  6.352055e-01  3.221557e-01 -5.466382e-01
##  [771]  6.289391e-01  8.169998e-01  1.642982e+00  1.895794e+00 -7.488971e-01
##  [776]  5.541494e-01  1.358862e+00 -5.135670e-01  6.115343e-01  1.796802e+00
##  [781]  8.842941e-02 -3.903048e-01 -3.086056e-02  1.885599e+00  2.149431e+00
##  [786]  1.177241e+00  1.181291e+00  1.062895e+00  8.106873e-02 -8.942028e-02
##  [791]  1.244536e+00  9.962563e-02  1.403039e+00  2.310079e+00  4.283067e-01
##  [796]  1.775074e+00  1.267219e+00  1.948705e+00  9.832296e-01  1.195211e+00
##  [801]  1.964678e+00  5.829682e-01 -1.490693e-02  9.664358e-01  1.417185e+00
##  [806]  1.372503e+00  1.905998e+00 -2.204100e+00  2.147553e+00  1.645172e+00
##  [811]  1.477569e+00  1.622533e+00  1.298438e+00  1.150801e+00  2.892533e-01
##  [816] -6.334600e-01  4.913557e-01  1.358313e+00 -2.250374e-01 -2.699442e-01
##  [821] -2.264862e-01  1.507114e+00  7.082902e-01 -1.789107e-01  1.231944e+00
##  [826]  8.668884e-02 -3.495376e-01  7.519248e-01  1.859855e-01  9.263572e-01
##  [831]  4.785260e-01  2.216121e-02  2.368244e-01  1.964088e+00  3.617292e-01
##  [836]  5.574844e-01 -2.702405e-01 -2.188216e-01  1.764009e+00  1.062033e+00
##  [841]  1.621639e+00  1.766858e+00  1.149533e+00  3.730043e-01  1.265231e+00
##  [846]  1.055124e+00  1.774088e+00  1.572043e+00  1.137305e+00 -2.441452e-01
##  [851]  1.199623e+00  1.876836e+00  6.365438e-01  8.206128e-02  1.671865e+00
##  [856]  1.365540e+00  1.814747e+00  1.172196e+00 -1.182556e+00  1.400122e+00
##  [861]  1.611672e+00  1.431096e+00 -5.074762e-01  1.591600e+00  1.494209e+00
##  [866] -6.946820e-01  1.039094e+00  1.159197e+00  1.606305e+00  1.758304e+00
##  [871]  1.605442e+00  1.819480e+00  1.621702e+00  1.340426e+00  1.618900e+00
##  [876]  1.609312e+00  1.233171e+00  1.610874e+00  1.268559e+00  1.943826e+00
##  [881]  1.477430e+00  7.000649e-01  1.918282e+00  1.280635e-01 -1.606418e-01
##  [886]  1.419860e+00  1.776719e+00  9.251224e-01  1.117519e+00  1.095538e+00
##  [891]  4.481344e-02  3.481187e-01  1.572550e+00  1.187688e+00  1.965598e-01
##  [896]  3.863184e-01  4.692375e-01  5.084818e-01 -5.682760e-01  2.077891e+00
##  [901]  1.500603e+00  7.781974e-01  1.157382e+00 -6.407133e-01  1.754062e+00
##  [906] -2.591483e-01  1.566302e+00  1.262907e+00  1.568365e+00  1.832588e+00
##  [911]  1.724195e+00  5.711931e-01  8.161973e-01 -1.771460e-01 -2.039906e-03
##  [916]  1.045785e+00  6.549156e-01  1.193123e-01  1.487279e+00  2.019245e+00
##  [921]  1.026788e-01  1.666001e+00  5.515050e-01  1.181291e+00  1.503190e+00
##  [926]  6.019399e-01  1.522608e+00 -4.555701e-01  2.137924e+00  1.670994e+00
##  [931]  1.396779e+00 -2.039906e-03  1.135117e+00  1.498920e+00  2.340860e+00
##  [936]  1.640673e+00  1.197221e-01  1.296080e+00  1.824403e+00  1.226217e+00
##  [941]  2.056387e+00  1.760886e+00  3.307333e-01  4.016758e-01  1.712827e+00
##  [946] -3.788409e-01 -1.489643e+00  1.606151e+00  6.171803e-01 -6.533252e-01
##  [951]  1.660613e+00  1.623944e+00  1.950225e+00  1.028076e+00  1.479930e+00
##  [956]  1.561594e+00  1.107196e+00  1.011740e+00  2.019098e+00  1.268914e+00
##  [961]  1.065837e-02  7.249272e-01  1.727043e+00 -1.333010e-01  1.672234e+00
##  [966]  2.126016e+00  1.542087e+00  1.526623e+00 -2.156540e-01  1.395185e+00
##  [971]  8.306521e-01  7.904108e-02  9.207986e-01  2.029923e+00  1.038618e+00
##  [976]  1.222842e+00  1.119974e+00  9.985999e-01  1.387758e+00  4.087461e-01
##  [981]  2.026271e+00 -8.104323e-01  1.114249e+00  1.386758e+00  1.009542e+00
##  [986]  1.702190e+00  8.974930e-01  1.577027e+00 -3.415412e-01  1.555322e-01
##  [991]  1.819879e+00  1.660287e+00  6.390466e-01  1.122225e+00  1.627202e+00
##  [996]  1.852074e+00 -2.257468e-01  1.052500e+00  6.034519e-01  9.827079e-01
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
##   0.52257140   0.39106997 
##  (0.12366718) (0.08744357)
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
## [1] -0.003082878 -0.390124409 -0.228022912 -0.166548859  0.088470077
## [6] -0.168511083
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
## [1] -0.0579
```

```r
sd(results2$thetastar)  #the standard deviation of the theta stars is the SE of the statistic (in this case, the mean)
```

```
## [1] 0.8966055
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
## t1*      4.5 0.03963964   0.9320341
```

One of the main advantages to the 'boot' package over the 'bootstrap' package is the nicer formatting of the output.

Going back to our original code, lets see how we could reproduce all of these numbers:


```r
table(sample(x,size=length(x),replace=T))
```

```
## 
## 2 3 4 7 8 9 
## 3 1 3 1 1 1
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
## [1] 0.0293
```

```r
se.boot
```

```
## [1] 0.8795369
```

Why do our numbers not agree exactly with those of the boot package? This is because our estimates of bias and standard error are just estimates, and they carry with them their own uncertainties. That is one of the reasons we might bother doing jackknife-after-bootstrap.

The 'boot' package has a LOT of functionality. If we have time, we will come back to some of these more complex functions later in the semester as we cover topics like regression and glm.

