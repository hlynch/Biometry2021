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
## 0 1 2 3 4 7 8 9 
## 1 2 1 1 1 1 2 1
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
## [1] -0.0024
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
## [1] 2.81308
```

```r
UL.boot
```

```
## [1] 6.18212
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
##    [1] 4.9 5.6 5.7 5.1 4.5 3.7 3.9 4.8 4.4 4.6 5.1 5.6 3.6 4.4 4.6 4.6 4.4 4.1
##   [19] 4.8 4.6 6.6 4.9 4.2 3.8 3.9 4.7 4.9 4.8 4.1 3.8 4.0 4.2 4.3 2.5 3.6 4.1
##   [37] 4.1 4.1 3.8 3.5 3.4 4.1 4.1 1.5 4.3 4.2 6.0 4.7 4.0 3.6 4.2 3.0 5.1 4.1
##   [55] 3.7 5.8 4.7 5.8 6.3 4.8 3.9 6.0 5.8 3.1 5.9 4.8 3.9 5.5 3.6 5.3 6.0 3.7
##   [73] 4.2 5.1 5.0 4.5 5.1 3.6 5.4 3.7 4.7 5.2 3.7 5.0 4.6 4.8 3.3 5.6 3.3 4.0
##   [91] 5.5 4.5 4.0 5.0 4.1 4.9 3.4 4.7 4.8 3.3 3.9 5.0 2.9 4.5 5.0 4.2 4.1 5.7
##  [109] 3.1 4.2 4.7 4.6 6.5 4.8 4.4 3.5 4.3 2.8 4.8 4.8 4.7 5.2 5.3 4.5 4.2 4.2
##  [127] 4.1 5.2 3.8 4.5 5.4 4.2 4.0 4.8 3.8 4.9 4.9 5.0 2.7 4.8 5.9 3.8 4.6 5.8
##  [145] 4.3 5.3 4.6 5.9 4.3 5.6 4.5 4.4 4.7 5.2 4.6 3.9 4.1 4.3 5.4 6.3 4.9 4.4
##  [163] 4.3 4.5 4.2 6.3 4.4 3.6 4.8 4.8 6.0 4.9 3.7 4.5 5.5 2.5 5.0 5.4 4.4 3.9
##  [181] 6.0 3.2 4.3 4.2 5.2 3.8 4.0 3.8 4.7 4.1 5.0 5.9 6.2 4.5 4.9 4.8 4.9 4.4
##  [199] 3.7 5.5 4.0 3.7 2.9 4.6 5.0 5.9 4.4 5.3 4.4 4.6 3.9 2.6 4.0 5.0 4.1 3.2
##  [217] 4.6 4.5 5.9 4.6 3.7 4.6 4.9 5.0 5.8 5.1 4.2 3.4 3.3 4.2 3.9 5.0 4.5 4.0
##  [235] 6.6 3.4 4.4 5.5 5.3 3.5 6.6 4.8 5.7 5.0 4.3 4.0 4.0 2.8 3.1 4.6 4.7 3.6
##  [253] 4.5 5.6 5.3 4.6 4.4 4.2 3.9 4.8 3.9 3.0 5.8 5.5 5.0 4.9 3.7 5.1 3.3 3.9
##  [271] 4.7 4.8 4.3 5.6 3.4 5.5 4.3 4.0 5.4 6.3 5.3 4.1 4.7 5.1 3.1 5.4 5.2 3.4
##  [289] 3.6 3.2 3.7 4.3 4.4 5.2 3.2 5.4 4.8 5.0 5.0 4.6 3.5 3.7 3.1 3.7 5.1 5.6
##  [307] 3.8 4.0 6.1 4.8 3.0 5.7 6.3 4.9 5.3 5.4 4.9 5.0 3.3 3.8 5.5 3.2 4.7 3.6
##  [325] 3.7 6.7 4.6 1.9 3.8 4.4 4.2 5.1 6.1 4.4 4.7 5.6 5.5 4.7 3.7 3.6 4.5 4.8
##  [343] 6.3 4.0 3.8 3.5 4.5 4.0 4.4 3.4 5.5 4.7 5.3 4.0 5.6 2.8 5.7 5.8 5.0 4.3
##  [361] 3.5 3.9 4.8 3.6 5.0 5.0 2.1 6.7 4.1 3.6 5.2 4.2 3.9 4.8 4.8 4.3 4.1 3.6
##  [379] 5.5 3.8 3.6 3.2 4.6 3.3 5.7 4.3 3.9 5.5 4.6 5.3 4.3 5.0 4.0 4.0 4.2 4.6
##  [397] 4.4 4.4 3.5 5.6 5.4 4.2 4.4 3.2 3.8 4.3 3.5 3.4 5.7 3.9 4.3 5.7 3.2 4.3
##  [415] 5.8 6.1 2.9 5.0 4.3 6.4 3.6 4.4 3.4 4.7 6.1 6.2 5.0 4.1 3.8 4.4 3.9 5.6
##  [433] 4.0 4.7 3.9 4.6 5.0 3.8 3.8 3.7 6.1 3.1 5.1 5.0 3.8 5.3 3.6 4.0 3.8 6.3
##  [451] 5.0 3.4 5.0 3.4 6.5 3.8 3.3 4.1 4.9 5.2 3.1 4.5 4.3 5.2 2.6 4.5 6.3 4.2
##  [469] 3.6 3.9 4.9 4.9 4.9 5.2 4.4 5.1 3.3 4.0 5.3 5.4 4.1 2.8 3.7 4.2 3.7 3.5
##  [487] 4.4 3.9 5.3 4.3 3.6 4.3 3.9 4.8 4.0 3.0 4.6 5.2 4.2 3.4 6.4 5.2 4.2 3.9
##  [505] 4.5 3.6 5.0 4.7 2.8 4.4 6.3 2.7 3.1 4.8 2.7 6.4 4.3 4.4 4.8 5.9 6.0 5.0
##  [523] 3.4 4.6 5.2 5.4 4.1 3.3 5.0 5.3 5.4 6.0 4.7 3.9 4.1 4.8 5.7 5.4 4.2 4.8
##  [541] 3.6 2.8 4.3 3.1 4.6 3.2 4.6 5.4 4.5 4.6 3.6 5.3 4.1 5.0 5.8 3.8 4.4 3.9
##  [559] 4.7 4.1 3.8 7.3 4.6 4.1 5.1 4.2 4.7 5.3 4.4 4.6 4.3 3.9 6.1 5.8 4.8 2.6
##  [577] 4.1 3.8 4.7 4.5 6.0 4.5 3.6 3.0 4.1 4.9 3.8 4.7 4.2 4.2 4.9 3.8 5.1 3.8
##  [595] 5.1 5.2 5.9 3.7 4.8 5.4 5.5 5.6 3.6 4.3 3.7 4.8 4.8 6.2 4.8 4.3 3.8 4.5
##  [613] 3.5 5.4 6.1 3.4 4.4 3.5 2.4 2.7 3.6 4.2 4.1 3.5 3.5 4.7 4.3 4.6 3.3 5.5
##  [631] 5.8 5.3 3.6 4.1 5.8 5.0 6.4 6.2 4.1 6.3 4.1 5.6 4.1 6.1 6.6 4.7 2.8 5.8
##  [649] 5.3 4.3 4.0 4.6 5.1 4.6 3.5 2.5 3.9 4.4 4.7 5.5 4.3 3.9 3.3 6.1 5.0 4.3
##  [667] 5.2 4.9 4.5 3.7 4.0 5.2 3.4 4.4 4.2 4.1 5.1 5.4 4.7 2.8 4.0 5.7 3.6 4.7
##  [685] 3.1 3.5 5.0 4.3 5.6 4.4 4.3 4.0 3.7 5.2 5.5 4.0 5.4 5.3 4.5 3.6 5.2 5.5
##  [703] 4.0 4.0 3.6 6.3 4.8 4.6 6.0 5.7 5.6 4.9 5.4 3.1 4.4 2.5 5.0 2.8 2.5 5.9
##  [721] 4.1 5.2 4.8 5.9 5.6 4.0 4.4 4.7 4.7 3.1 3.5 4.4 3.9 4.1 3.8 5.0 4.0 5.1
##  [739] 5.2 4.4 4.6 5.0 4.2 5.2 5.4 4.6 5.6 3.1 4.7 4.5 4.3 5.2 3.7 5.3 4.9 6.1
##  [757] 2.3 3.8 6.1 4.4 4.9 3.6 4.4 3.4 2.9 4.3 5.3 3.4 3.2 4.2 4.9 4.1 3.8 3.4
##  [775] 4.9 4.2 3.2 3.5 4.9 3.4 5.7 5.0 3.8 4.1 5.5 4.9 5.2 5.0 5.0 3.9 3.4 5.2
##  [793] 4.9 3.7 4.5 4.1 5.5 4.6 5.1 5.9 3.1 4.9 3.8 3.4 3.8 3.3 3.8 3.2 3.9 4.4
##  [811] 3.5 5.2 5.4 5.8 3.8 2.6 6.3 4.5 4.5 4.0 3.4 4.1 4.6 3.4 4.9 4.5 6.0 4.3
##  [829] 4.8 3.1 5.7 4.1 5.8 4.7 5.4 6.4 4.4 4.5 3.4 3.7 3.4 4.7 5.6 3.2 5.7 3.6
##  [847] 3.6 4.0 4.0 3.4 3.6 4.3 3.9 3.3 4.6 5.1 3.5 4.7 2.9 4.5 5.2 4.7 5.0 4.5
##  [865] 4.2 3.8 4.5 4.5 4.1 4.7 5.3 4.2 3.5 4.8 6.0 4.4 3.6 5.0 5.1 3.7 4.5 6.4
##  [883] 4.9 4.7 4.7 5.2 5.9 4.4 3.3 3.3 4.7 4.2 4.7 3.9 4.2 5.0 3.6 3.6 4.8 5.8
##  [901] 3.5 4.5 5.1 4.7 5.7 4.0 5.2 4.6 3.9 4.2 5.0 4.9 5.9 6.1 4.6 3.7 5.1 3.9
##  [919] 4.9 5.4 5.7 4.1 2.9 4.0 4.7 3.3 4.7 4.2 5.4 4.4 4.5 4.9 4.2 6.4 4.5 6.1
##  [937] 4.1 2.2 4.4 3.8 3.6 5.2 5.3 4.1 3.5 3.0 6.1 4.7 2.9 4.6 5.6 4.2 4.8 4.2
##  [955] 5.4 5.2 2.9 5.2 3.3 4.7 4.4 5.0 4.2 4.8 4.9 4.9 3.8 4.4 4.7 3.9 4.5 4.2
##  [973] 4.0 4.9 4.5 4.5 5.0 4.9 6.2 4.7 4.8 4.9 4.0 3.9 6.1 5.2 4.8 4.1 4.8 5.2
##  [991] 3.6 3.5 4.4 5.2 4.8 5.4 5.5 4.7 3.9 4.5
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
##    [1] 5.4 5.4 3.9 4.5 2.1 5.2 3.5 5.5 3.8 3.8 6.2 2.9 6.7 5.7 4.1 3.1 7.1 5.3
##   [19] 4.7 3.9 5.1 5.2 5.0 4.0 4.2 6.0 3.4 5.7 3.3 5.3 4.2 4.1 4.4 4.5 4.9 3.4
##   [37] 4.7 5.9 5.3 3.6 3.9 4.4 5.3 5.4 4.7 3.2 4.2 4.9 3.6 4.0 4.4 5.2 4.4 5.2
##   [55] 3.5 5.9 4.5 5.0 3.4 4.5 3.1 5.3 5.2 3.9 6.4 4.1 2.6 4.4 4.4 4.8 5.1 5.3
##   [73] 3.8 4.3 5.3 4.8 3.3 5.8 3.8 3.4 4.4 3.1 4.0 3.9 4.1 5.0 4.1 6.2 4.4 4.0
##   [91] 3.7 4.0 5.8 4.5 4.8 3.8 4.9 4.2 5.7 3.9 4.5 5.9 5.3 5.5 5.9 5.4 3.9 4.1
##  [109] 4.2 4.3 6.1 4.3 4.5 3.6 4.9 4.8 5.3 5.1 3.8 4.0 3.6 4.2 4.1 4.2 6.1 4.1
##  [127] 5.0 4.3 3.5 4.5 7.3 4.4 2.9 3.8 5.0 4.4 4.6 3.9 3.6 4.5 3.6 4.7 5.8 4.6
##  [145] 3.7 4.9 4.6 5.4 3.9 3.0 5.1 5.6 3.7 6.0 4.5 5.8 4.3 6.1 3.6 4.5 5.7 5.3
##  [163] 2.9 3.2 3.6 4.3 4.2 4.2 5.6 3.3 4.8 4.4 5.2 3.7 3.9 4.1 6.3 2.8 6.2 3.9
##  [181] 4.0 5.1 5.2 2.9 3.4 3.8 5.2 4.6 3.3 5.7 5.2 5.0 4.8 5.2 4.7 4.0 3.2 4.3
##  [199] 4.5 5.4 5.2 5.7 3.5 3.4 3.7 4.5 5.9 3.5 3.9 5.8 3.5 3.5 5.6 3.6 4.7 4.7
##  [217] 3.7 5.0 4.5 4.1 6.4 4.0 5.0 3.3 4.5 3.6 4.6 5.2 3.4 6.2 3.6 5.3 4.0 3.4
##  [235] 3.7 4.9 3.4 4.2 4.4 5.5 4.9 5.1 4.0 5.0 3.4 7.1 4.8 4.3 3.7 4.7 4.4 3.6
##  [253] 3.7 4.4 4.3 3.3 4.0 5.2 4.5 5.0 3.9 4.1 3.4 4.5 3.8 5.0 5.8 3.8 4.3 4.2
##  [271] 4.0 3.7 4.1 3.8 4.7 5.0 5.4 3.6 4.5 3.9 4.7 2.5 5.2 6.7 5.7 6.0 4.7 4.8
##  [289] 5.8 5.3 2.7 3.8 4.5 6.0 4.3 4.0 4.7 5.1 2.2 4.9 4.7 7.0 3.8 6.3 4.4 4.1
##  [307] 2.9 4.9 5.8 3.2 4.4 5.7 5.1 4.2 6.1 5.7 2.5 5.3 4.0 5.0 5.3 3.9 4.9 3.9
##  [325] 4.5 3.7 3.0 4.7 4.6 4.8 5.1 4.6 4.8 6.8 5.0 2.9 4.6 5.0 4.9 4.6 4.1 3.7
##  [343] 5.3 4.0 5.2 5.9 4.2 4.0 4.8 5.8 4.3 4.2 4.7 4.6 4.8 4.9 5.0 5.0 4.4 6.6
##  [361] 4.0 4.1 4.6 3.9 5.6 5.5 5.6 3.4 3.9 4.8 4.0 6.6 4.6 4.8 4.1 3.1 4.2 3.0
##  [379] 3.0 6.2 4.5 4.7 2.3 6.8 5.3 4.3 4.1 5.0 4.4 3.5 5.6 4.2 5.2 4.2 4.1 5.4
##  [397] 4.0 4.3 3.7 3.4 4.2 5.2 3.6 4.7 5.6 4.6 4.6 3.9 5.4 5.3 5.6 4.9 4.6 4.6
##  [415] 4.9 3.1 4.9 3.2 3.7 3.8 4.5 4.4 3.8 6.0 2.7 2.1 5.3 4.2 5.3 4.9 5.2 5.0
##  [433] 3.9 4.0 5.0 4.8 5.5 5.1 4.2 4.4 5.1 5.3 2.9 4.5 4.0 3.0 5.4 4.7 2.6 3.5
##  [451] 4.2 4.4 4.4 6.1 3.7 5.1 3.0 4.5 4.6 4.8 2.8 6.0 5.2 5.2 3.7 4.6 4.8 5.0
##  [469] 4.4 3.7 5.3 4.7 4.5 6.0 5.2 3.7 5.0 5.0 3.8 4.7 5.2 5.8 3.9 5.6 4.0 5.2
##  [487] 3.2 4.1 6.3 4.3 3.4 4.6 4.7 3.9 4.2 4.4 4.7 5.0 4.3 5.0 6.0 5.1 3.7 4.9
##  [505] 4.5 4.3 3.5 5.0 4.1 4.4 3.6 3.6 4.5 5.1 3.0 5.2 3.8 4.7 4.7 4.0 5.4 4.2
##  [523] 3.9 5.5 3.8 5.3 5.3 5.6 5.5 3.3 3.3 3.5 4.7 4.2 4.9 5.7 4.9 2.7 3.5 5.2
##  [541] 3.6 5.7 6.6 4.2 4.3 4.0 3.6 3.9 4.1 3.0 2.4 4.6 5.9 4.2 5.3 4.8 5.7 4.7
##  [559] 4.2 6.3 3.3 6.6 5.1 5.1 3.8 4.5 3.0 3.2 6.7 4.7 5.4 3.3 4.2 4.2 4.7 4.9
##  [577] 5.8 3.9 5.9 3.5 4.7 5.8 4.5 4.4 5.1 4.4 4.5 3.4 3.4 4.8 3.8 4.7 3.8 4.1
##  [595] 4.4 4.9 4.0 2.5 3.5 4.6 4.2 7.0 4.1 3.9 4.4 3.6 4.2 3.5 4.9 5.7 3.3 3.9
##  [613] 5.0 4.4 4.6 5.3 2.8 4.1 4.6 4.5 4.5 5.8 4.9 5.3 4.6 3.3 4.9 3.8 4.3 5.3
##  [631] 5.6 3.7 4.9 4.7 5.2 3.5 4.1 3.6 3.9 4.2 3.3 4.7 4.4 4.0 4.1 4.4 5.4 4.0
##  [649] 4.4 4.3 4.2 4.3 4.5 3.4 3.5 4.5 6.5 5.0 5.3 4.6 4.8 4.5 4.1 5.2 5.2 5.9
##  [667] 2.7 5.8 5.2 4.0 5.1 3.6 5.9 4.4 3.5 3.1 4.8 5.3 3.8 6.2 3.2 5.9 4.6 5.2
##  [685] 4.3 3.6 4.6 4.7 5.1 2.6 4.0 3.2 5.3 4.0 4.6 5.2 4.0 4.7 2.5 3.8 4.2 4.0
##  [703] 5.7 2.2 6.2 3.7 4.6 4.3 4.6 4.2 4.1 4.2 4.7 4.1 4.5 5.8 4.0 3.7 6.4 4.3
##  [721] 3.8 4.0 3.4 3.3 4.8 3.6 5.6 4.8 4.5 3.1 4.9 6.2 5.7 3.8 6.3 5.2 3.5 4.4
##  [739] 2.9 3.9 4.4 6.1 4.3 5.2 5.5 5.7 5.1 5.8 2.5 3.7 4.7 4.4 5.2 5.1 6.2 5.3
##  [757] 3.5 4.5 4.6 4.6 5.4 5.7 4.8 3.4 4.6 5.5 5.1 5.5 4.8 5.8 6.1 3.7 5.3 3.4
##  [775] 3.5 3.8 4.3 4.0 4.5 4.8 4.1 3.5 3.8 4.3 4.6 6.3 3.2 5.2 4.6 5.9 4.5 5.8
##  [793] 4.5 4.5 4.5 4.5 4.2 3.5 4.0 5.4 3.4 4.1 4.6 5.5 5.6 4.8 4.5 4.9 5.5 3.2
##  [811] 4.1 4.8 4.7 4.4 5.4 6.6 5.4 3.4 4.0 4.1 3.6 4.4 4.2 4.8 3.9 3.7 3.2 3.6
##  [829] 5.1 4.8 7.0 6.1 4.7 5.3 5.4 4.3 3.6 4.2 3.6 3.7 4.7 6.0 4.3 3.9 3.7 4.3
##  [847] 5.7 5.6 5.5 4.3 4.7 4.2 3.6 3.5 4.1 4.6 4.1 4.7 5.2 3.0 3.2 4.6 3.2 5.3
##  [865] 4.2 5.4 4.3 4.0 3.7 5.5 4.5 4.3 4.1 3.8 4.8 4.9 5.0 5.5 5.0 5.0 5.3 6.1
##  [883] 4.2 4.8 4.2 6.3 3.9 4.5 3.9 4.2 5.1 4.6 4.3 3.2 3.5 3.6 2.3 3.5 5.4 2.4
##  [901] 5.3 4.1 4.2 4.4 5.0 4.0 5.2 4.4 5.4 5.9 4.5 6.0 4.4 4.0 4.4 6.2 3.4 3.9
##  [919] 5.2 4.7 3.7 6.0 3.4 5.2 3.8 5.1 4.9 3.6 4.5 3.8 4.8 3.7 4.5 4.4 3.9 5.3
##  [937] 5.1 4.4 3.2 4.1 5.2 4.1 3.7 5.0 3.8 3.0 5.1 6.3 5.2 3.4 4.4 4.4 3.1 4.3
##  [955] 3.1 4.2 5.1 5.4 4.8 4.7 5.3 4.7 3.4 5.8 4.8 6.0 5.5 5.8 5.9 4.8 6.4 3.3
##  [973] 3.1 5.3 4.1 5.2 3.7 3.7 2.8 5.6 4.0 5.0 5.2 4.1 4.1 2.8 4.5 5.0 5.7 5.1
##  [991] 5.0 3.8 4.4 3.5 4.1 4.6 5.3 4.5 4.0 4.2
## 
## $func.thetastar
## [1] 0.0148
## 
## $jack.boot.val
##  [1]  0.472375691  0.362921348  0.318786127  0.199719101  0.070679012
##  [6]  0.003216374 -0.158244681 -0.250867052 -0.317451524 -0.460112360
## 
## $jack.boot.se
## [1] 0.8980049
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
##    [1] 4.6 3.6 3.7 6.1 5.1 3.8 4.4 6.0 5.0 4.8 3.7 3.3 3.0 4.0 5.9 3.7 2.9 4.9
##   [19] 4.4 6.0 5.6 4.2 4.9 5.0 3.6 4.5 4.2 4.4 4.4 2.7 3.6 4.1 6.0 4.6 2.2 6.5
##   [37] 6.2 4.6 3.9 3.2 3.5 3.1 3.8 4.0 4.7 5.7 5.6 5.0 6.7 5.2 6.0 4.0 4.6 4.9
##   [55] 3.7 7.0 6.2 5.0 5.2 4.6 6.2 5.2 5.6 5.8 4.6 5.4 3.5 5.9 3.9 5.5 5.3 5.6
##   [73] 4.2 5.9 5.5 6.2 5.3 4.7 3.6 3.7 4.9 5.7 4.0 4.1 4.9 4.2 4.4 4.5 4.6 3.5
##   [91] 5.9 4.6 4.1 4.9 2.8 4.8 2.8 3.9 6.1 5.4 3.2 5.2 3.1 3.3 6.4 4.9 4.7 5.5
##  [109] 3.7 3.7 5.2 3.9 4.1 5.6 5.1 5.3 4.5 6.8 5.4 3.9 3.9 3.6 3.2 3.3 4.6 4.7
##  [127] 4.8 3.9 4.4 2.4 3.3 5.1 5.4 4.3 4.5 4.8 5.1 3.1 4.8 5.8 3.1 5.1 4.1 4.9
##  [145] 4.2 4.4 4.1 5.2 3.8 5.2 4.9 4.4 4.9 4.3 3.3 6.5 5.9 4.9 4.4 4.0 4.2 4.3
##  [163] 4.4 2.9 4.7 5.8 4.8 4.8 4.0 4.1 3.5 3.1 3.7 3.8 3.7 3.3 4.0 4.1 4.9 5.1
##  [181] 3.8 4.2 4.0 4.8 4.3 5.6 4.3 5.3 5.8 4.5 3.5 5.3 5.1 6.0 4.6 4.5 3.6 6.1
##  [199] 4.5 5.7 4.2 3.3 4.3 3.7 4.8 6.5 4.2 4.5 4.5 4.4 3.2 4.3 4.7 4.7 6.3 4.1
##  [217] 4.2 5.2 5.1 3.2 5.3 4.6 5.3 4.9 5.5 5.0 5.6 2.7 4.4 3.9 3.2 4.8 5.5 3.6
##  [235] 4.1 5.3 5.0 5.3 3.9 4.5 3.4 4.0 4.6 4.7 4.2 3.8 5.6 3.4 6.6 2.9 5.1 2.6
##  [253] 5.0 4.7 6.1 3.4 5.6 4.3 4.6 5.6 4.1 5.0 4.6 3.2 4.3 6.6 4.6 3.8 5.3 6.0
##  [271] 5.3 4.4 4.0 5.6 3.4 4.5 3.4 5.6 4.3 3.1 5.1 6.9 4.3 3.5 5.2 5.4 4.2 4.5
##  [289] 4.5 5.7 3.8 3.7 4.8 4.1 4.3 3.0 4.7 3.4 5.3 5.7 3.6 4.4 4.2 5.1 4.7 3.7
##  [307] 4.5 4.4 3.8 4.3 3.8 4.5 5.0 4.5 4.0 4.1 3.6 4.3 4.7 4.8 3.5 3.8 4.5 4.9
##  [325] 2.9 5.0 5.6 4.1 5.9 5.0 5.3 4.7 3.6 5.1 5.2 2.6 5.7 4.3 3.8 3.7 4.7 5.0
##  [343] 4.1 3.3 5.6 2.6 4.9 3.5 3.6 5.3 4.5 4.3 5.0 3.2 4.6 3.7 3.6 3.6 4.9 5.9
##  [361] 3.2 4.9 4.3 4.9 3.6 4.8 4.9 5.5 5.9 4.1 6.4 3.8 4.0 4.1 5.5 6.4 5.1 5.3
##  [379] 4.3 4.0 5.5 4.2 3.9 3.2 4.2 4.3 4.9 4.0 5.0 5.4 3.9 4.7 4.9 5.9 4.1 3.6
##  [397] 4.9 4.5 5.2 2.9 1.5 2.7 4.3 3.7 4.9 4.2 3.9 4.4 3.1 4.5 5.4 5.3 4.4 3.1
##  [415] 5.6 5.9 4.7 4.5 3.9 5.0 5.6 3.3 4.1 5.0 3.8 3.9 3.8 3.5 5.6 3.0 4.5 5.2
##  [433] 4.6 3.1 4.5 5.4 7.2 4.0 4.7 5.2 5.6 4.4 3.9 4.1 5.3 4.5 4.9 4.2 5.1 4.3
##  [451] 3.4 4.2 3.9 5.5 6.5 4.4 5.5 5.3 5.1 4.9 4.1 4.2 5.5 4.5 5.8 3.9 3.9 5.6
##  [469] 4.8 3.1 3.9 3.8 5.3 4.4 3.9 4.4 3.7 4.8 3.5 3.9 3.4 3.6 4.8 4.0 4.4 3.9
##  [487] 4.8 5.4 2.8 5.0 5.7 6.3 4.1 5.7 4.6 4.9 5.6 4.4 7.5 5.5 3.4 4.9 6.5 3.8
##  [505] 3.5 5.1 5.1 6.0 4.2 4.7 3.4 5.3 3.8 4.2 4.1 3.8 3.1 5.7 4.4 5.5 6.3 6.3
##  [523] 4.5 4.8 3.9 5.6 4.2 4.7 4.9 5.7 4.5 5.2 3.8 4.6 4.3 3.2 4.3 3.3 4.5 3.8
##  [541] 4.6 4.5 4.5 4.2 4.2 4.7 3.6 5.5 5.0 2.4 3.2 4.0 6.6 5.2 4.0 3.7 2.8 4.8
##  [559] 4.3 4.1 3.7 3.7 4.4 2.7 4.1 5.4 4.9 6.7 5.7 5.1 3.7 5.4 4.8 4.2 4.8 5.6
##  [577] 4.1 2.8 4.5 3.6 5.5 4.3 4.9 5.3 4.1 4.8 4.2 3.9 5.1 5.0 4.0 4.8 5.1 5.0
##  [595] 3.9 4.4 5.0 3.3 4.2 4.5 4.7 5.2 3.9 3.8 3.5 5.7 3.4 4.7 4.1 4.3 4.7 2.7
##  [613] 3.4 5.5 5.2 4.0 3.7 5.8 5.1 3.3 4.3 4.2 3.7 5.0 5.2 5.0 4.4 4.9 5.6 3.6
##  [631] 3.9 4.2 4.8 4.7 5.6 4.6 6.0 3.8 4.4 4.8 4.8 2.8 4.3 4.1 4.2 6.1 2.8 5.0
##  [649] 3.2 4.6 5.5 3.7 3.4 3.4 2.3 3.3 6.6 5.6 6.5 4.6 5.3 4.0 3.9 3.3 4.3 4.5
##  [667] 3.9 5.8 3.5 4.5 3.2 5.5 4.1 4.3 4.2 3.6 5.1 3.3 4.1 3.4 3.7 4.1 5.4 4.1
##  [685] 5.1 4.5 4.4 3.7 5.6 4.4 4.1 5.1 6.8 5.1 4.3 4.3 5.7 3.3 4.2 5.7 3.9 4.8
##  [703] 3.1 4.3 3.1 4.8 6.0 3.8 4.5 5.3 4.3 3.9 4.2 5.5 5.2 7.3 5.1 5.2 4.3 5.2
##  [721] 3.4 4.6 4.8 3.2 3.8 5.0 5.3 5.0 3.6 4.1 3.7 4.7 5.0 4.7 3.4 3.9 3.4 5.9
##  [739] 5.9 5.8 3.5 4.1 5.2 4.4 5.3 3.0 3.7 4.8 6.1 5.7 4.3 5.5 3.7 4.5 4.8 4.9
##  [757] 3.8 2.6 5.0 4.8 5.0 4.1 5.7 3.1 3.3 4.3 3.8 4.3 4.2 4.4 4.3 4.3 3.7 5.0
##  [775] 4.9 4.1 4.9 2.3 3.7 4.4 3.6 5.1 5.4 5.1 5.2 2.2 5.8 5.5 5.5 3.1 4.7 3.4
##  [793] 5.1 4.9 4.1 5.7 5.5 5.2 4.2 3.6 5.6 4.1 5.5 4.3 4.4 4.3 4.2 4.0 5.5 3.0
##  [811] 5.1 4.1 5.4 4.4 5.6 3.3 3.6 4.8 3.7 4.3 4.7 3.9 4.5 5.1 5.3 5.2 4.4 4.1
##  [829] 4.3 5.1 4.6 5.3 5.3 4.8 4.6 5.0 3.1 4.0 6.3 4.9 4.6 5.6 4.7 3.1 3.7 4.8
##  [847] 4.2 4.7 4.0 3.8 4.7 3.2 3.7 3.1 4.9 2.9 4.2 2.9 4.5 5.6 5.5 4.6 3.8 4.6
##  [865] 5.9 5.4 6.1 3.4 3.9 4.7 5.2 5.4 3.8 5.2 4.7 3.5 6.1 4.2 4.5 4.3 4.8 2.9
##  [883] 4.3 3.2 4.9 5.7 4.6 5.1 4.9 4.3 2.6 5.0 5.5 3.9 5.9 5.8 4.7 5.5 5.9 2.8
##  [901] 4.3 3.9 5.3 3.0 4.5 5.0 3.7 5.5 5.7 5.6 5.7 5.9 4.8 4.6 4.3 4.8 4.8 5.3
##  [919] 3.9 4.4 4.2 6.0 4.8 3.9 3.5 3.1 5.8 4.6 4.3 6.1 4.4 2.7 4.8 5.6 3.8 5.9
##  [937] 5.1 3.2 4.8 4.6 4.2 5.0 5.2 5.2 3.6 4.2 4.6 4.6 5.3 3.6 3.8 3.8 4.6 3.8
##  [955] 3.3 3.9 6.5 5.0 4.2 4.4 7.6 4.7 3.1 5.1 4.7 3.3 3.4 3.5 4.8 6.5 3.2 4.6
##  [973] 3.6 4.0 5.3 4.6 4.5 3.5 6.0 5.5 3.1 4.8 4.4 3.0 2.7 4.9 4.5 4.7 4.5 4.6
##  [991] 4.1 2.8 4.2 3.6 5.3 4.9 5.3 3.2 3.9 4.6
## 
## $func.thetastar
## 72% 
##   5 
## 
## $jack.boot.val
##  [1] 5.500 5.500 5.372 5.200 5.100 5.000 4.900 4.800 4.700 4.500
## 
## $jack.boot.se
## [1] 0.9701026
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
## [1] 0.4488522
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
##   2.841484   3.754352 
##  (1.203637) (1.739303)
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
## [1]  0.31258029 -0.05280945  1.04425242  1.56438871 -0.32712011  0.67465511
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
##    [1]  0.8068041390  0.1452394885  0.3623013545  0.6154221182  0.2057111872
##    [6]  0.5179057362  1.3103945027 -0.1993308181  0.4451628127 -0.5777107699
##   [11]  0.1265383651  0.6406593555  0.2858034991  1.5007790567  0.6844891604
##   [16]  0.5033145296  0.4837972532  0.8203509234  0.7923416069  0.9747980351
##   [21] -0.1706603419  0.3804952609  0.3185851055  0.1199930710  0.4879924203
##   [26] -0.5320416088 -0.0331932406 -0.8206124351  0.7895551536  0.3128921825
##   [31]  1.2405031849  0.1867565240 -0.1661948511  0.7478619319  0.2272175661
##   [36]  0.6124011872 -0.1589351785  0.1969381911 -0.2947251024  1.3491152236
##   [41]  0.2241959729  0.7823759332  0.2317323496  1.3859122639  0.5315291350
##   [46]  0.6973419914  0.8317779187  1.8235466683  0.1181117844  0.9709647998
##   [51]  0.6791256594 -0.5786318389  0.3680496811  0.2453287847  1.7048081272
##   [56]  0.6142035765 -0.3818891737  0.9171153789  0.3992422474 -0.2805588713
##   [61]  0.8877339619 -0.1957188227 -0.5183705545  0.2169813301  0.2166658521
##   [66]  0.0244307419  0.8236091316 -0.1243780971  0.4729663746  0.2844829089
##   [71]  0.3166107682 -0.3865953842  1.8939951063  0.2909969366  1.4576351707
##   [76] -0.4145853399  0.0051527979  0.8486064181  0.1661433776  0.4918836249
##   [81]  0.2779940094  0.0958642500  0.3107949533  1.0250946070  1.3387969851
##   [86]  0.3489737011  1.0565322012  0.1102613451  0.3129307985  1.1124263124
##   [91]  0.2387041368  0.2508214395  0.8796337755  0.2154049411  0.1490890124
##   [96]  1.3088492000  0.1324444764  1.3025620870  0.1854056224  0.3927124000
##  [101]  0.3484170512  0.7658949649 -0.1530054418  0.3809497893  0.6139567643
##  [106]  0.7886036038  0.5017093669  0.8460965446  0.2553595660  0.3648122554
##  [111] -0.2209487051  0.8136384633  0.4762247622 -0.5688355861  1.2812828761
##  [116]  0.6031307294  0.6251772875  0.1322479284  1.2522908169  0.3490126338
##  [121] -0.1458893951  1.1246808964  0.4188641884  0.1297210122  0.5005939764
##  [126]  0.0134930076 -1.0569509747  0.0065769452  0.1579047318  1.1806868622
##  [131]  1.4259913710  0.5768077541  0.1815604115 -0.1822845185  0.0749483777
##  [136] -0.0802346360  0.0277500572  0.5307219416  0.3198051066  0.7015670046
##  [141] -0.3130207255  0.4978045245 -0.1583187888  1.3145383221 -0.4165967315
##  [146]  0.1783859615  0.5966050993  0.0590285494  1.0574691205  1.2287667083
##  [151]  0.4434437989  0.0198174011 -0.3814804128 -0.5777107699  0.2751108205
##  [156]  0.4935426857  0.0593206935  0.2070109926 -0.1887885041  0.0343925766
##  [161]  0.3921458944  0.7760035886  0.1587558637  0.0743866785 -0.0734043896
##  [166]  0.9721858360  0.4522761981 -0.0931260377  0.0589440708 -0.0940719888
##  [171]  0.9648354029 -0.3055312552  0.4782020407  0.2091922828 -0.0849575031
##  [176]  0.8959458183  0.3706089482  0.3341664244  0.6590047258  0.3721325065
##  [181]  2.1244906048  0.6456436184 -0.2575707368  1.0639082013 -0.1795087224
##  [186]  0.0051527979  0.4806122354  0.4798618633  0.0331836041 -0.4613478732
##  [191]  0.4154819270 -0.0540053635 -0.0115127282  0.4759246069  1.1687767404
##  [196]  0.2397575363 -0.8176470416  0.7948256981  0.2058629889  0.3319135431
##  [201]  0.4602554962  0.4326404006  0.4159446727  0.3615234684 -0.2054589111
##  [206] -0.4133845868  2.1181952482  0.3761830660 -0.2107076469 -0.1141404722
##  [211]  0.1223153392  0.7383249225  0.6003054762  0.7822150060 -0.7074190962
##  [216]  0.1864633793  1.0789034969  1.5201697307  1.2741357656  0.6196506596
##  [221]  0.0952722383 -0.6765120667  0.3228198012 -0.5807440781  1.6880919718
##  [226]  1.0201353729  1.3547938679  0.0533974483  0.8610723676  1.3415606316
##  [231]  0.2507326124  0.2607055684  0.7168761501 -0.0124965616 -0.5398318661
##  [236]  0.9048892640  0.1156595950  0.5735082375  0.2502293123 -0.2428251467
##  [241]  1.2035805284 -0.2768328302  1.3017763437 -0.1800032844  1.1622934839
##  [246]  0.3860540165  0.8709245411  0.0815262515  0.1974026524  1.9892326676
##  [251]  1.2671746097 -0.3282439565  0.7401373971  0.5070859324  0.7830692775
##  [256]  0.6562859576  0.4668922420  0.3582600479 -0.2710356614  0.7149615112
##  [261] -0.3289915719  0.3921458944  0.9989396071  0.5186475496  0.8519766916
##  [266]  0.3565559412  0.9841581138 -0.4030737039  0.4782020407 -0.0477679903
##  [271]  0.8930014978  0.0698019341  0.6871902664  0.1558354321  0.8199646641
##  [276]  1.0594729441  0.1473799731 -0.0048021412  0.5325269850 -0.1912020553
##  [281]  0.6287878424  0.9542615977  0.1904165404  0.4084655066  0.4273048085
##  [286] -0.0598781103  1.3352586305  0.6285598431  0.3627768041  1.2522908169
##  [291]  0.4458823320 -0.6496875239 -0.8247788406  1.0226311329  0.2294661673
##  [296]  1.0128075726  0.8884154765  0.4835138146  0.6369593192  1.0863101754
##  [301]  0.0597887552  0.4521794290 -0.5391082764  0.8820415142 -0.2018646925
##  [306]  0.8720792001  0.3162436855  0.0091350780 -0.5533131074  1.8242611207
##  [311]  0.5102514227 -0.2077800673  0.3475454925  0.0213763795  0.9041375518
##  [316]  0.7228946854 -0.2429319384 -0.5290210236  0.2239140839  0.9324927561
##  [321]  0.3251466734  0.1306182821 -0.2649544117  0.1660016809  0.2273813243
##  [326] -0.7126629267 -0.4805003646  0.3732358028  1.4052913242 -0.5527178178
##  [331] -0.9255574812 -0.1808292955 -0.8358286801 -0.0778800465  0.5342513385
##  [336] -0.6831609753  1.0852319812  0.1140314406  1.1976472227  0.2497444808
##  [341] -0.4532744063  0.5189500708 -0.0436977342  0.1425286971  0.9463306314
##  [346]  1.0098052453  0.4727178090  0.3703319554  0.3387086302  0.2904852550
##  [351]  0.1622228396  0.1646369613  0.4110078002  0.2123339491  0.3840805234
##  [356]  0.8475937516 -0.1875765284 -0.5298124864  0.0403241168  0.6018144869
##  [361]  0.3764695544  1.1518085955  0.0599047900 -0.1936503126  0.1777343866
##  [366]  1.0537691967 -0.4291568385  0.5051560735  0.2608607177  0.1289183601
##  [371]  0.7033475414  0.7205304270  1.3216229766  1.5341998495  0.1439792319
##  [376]  0.0322570145 -0.5620706226  0.3943704363  0.0872444711  0.8615407853
##  [381]  0.2814546715  0.1983070050  0.1617921683  0.0380713445  0.7906197696
##  [386]  0.1144044956  0.4756936792  0.1016307279  0.6285598431 -0.5025178736
##  [391]  0.9845484184  0.7744901285  0.9225273091  0.1850176651  0.1722621023
##  [396]  0.4681104986 -0.6382003409 -0.0441963660  0.5852807694  0.8979214342
##  [401]  1.0241261557  0.7639731236  1.0738066159  0.9243249339  0.6047571151
##  [406]  1.1748183951 -0.2808174722  1.0802980440  0.5571596196  0.5524476162
##  [411] -0.4066540227  0.7339066855  0.8526512146  0.3425749264  0.3953339379
##  [416] -0.0720993304 -0.5305122056  0.1175595499  1.0276905923  0.9231035734
##  [421]  0.9427319165  0.6471326931  0.7344489668  0.3742553461  0.8998457453
##  [426]  0.2993323576 -0.2808170825  0.0347352796  0.7547414661  0.2795290958
##  [431]  0.8884154765  0.3521474508 -0.6951695503  0.0131348964  0.4751114491
##  [436] -0.2208356408  0.1183671545  0.0674807334  0.0276545740 -0.8055917391
##  [441]  1.2425408498 -0.7906079686  0.7362850225  0.4885136894  0.5854892668
##  [446]  0.0457018761  0.4375655052  0.2015477057  1.5092032742  0.5296557722
##  [451]  0.0498505426 -0.4553615684  1.6318309994 -0.1768095473  0.2892723054
##  [456] -0.0755329713  0.2880289426  0.0423639626 -0.0009374136  0.5182149493
##  [461]  0.5025354986  1.1537545745  0.7914440217  0.5425785680  1.6316215216
##  [466] -0.5698824365  0.9316575700  0.3927897672  0.0389236729  0.1181841362
##  [471]  0.3718445177  0.7434512882  0.6919502028  0.8833068401  0.2403337843
##  [476]  0.3598246386  0.4402179553 -0.4194468466  0.8463528586 -0.3084232332
##  [481]  1.0820897694  0.0092821517  0.9948750965  1.7359263993 -0.2080324780
##  [486]  0.8180510841  0.4573785078  0.7449261578  0.8976689587  1.1817092938
##  [491]  0.8141042829  0.0792801164  1.3891170491  0.4397610456 -0.8469316446
##  [496]  0.2581107120  0.9631960616  0.2332704997  0.2886600187  1.3669492479
##  [501] -0.5152878132  0.5389356927  0.7482511286  0.3726496055 -0.0074829176
##  [506]  0.6489237964  0.2302613599 -0.1615951939 -0.2608895640  0.4765310079
##  [511]  0.7484595375  1.2339136715  0.2725401799  1.3836654103  0.0541947809
##  [516]  0.1848153021 -0.1418269407  0.4698266894  0.0615649377  1.1776286217
##  [521]  0.5699330950  0.4707677655 -0.3817676492  0.4748870081 -0.1045247496
##  [526]  0.5276626261  0.8670597802  1.1689270777  0.0721406723  0.7318716681
##  [531]  1.3342532932  1.3663327687  0.4737217997  1.2625945953  1.2754783033
##  [536]  1.9994673199 -0.1396696432  0.1983070050  0.0885087074 -0.3367676408
##  [541]  0.1308963442  0.0506361783  0.5877103238  1.2362114308 -0.0121543453
##  [546]  0.0197215751 -0.1005307585  0.8613170717  0.4962415412  0.3392480155
##  [551]  0.4628899444  0.4266306765  2.0073697893  0.4687622391 -0.4928138248
##  [556]  1.3150508311  0.5424891840  0.4137291355  0.7476075840  0.1826369882
##  [561]  0.2058902428 -0.2205359019 -0.8538766820 -0.1882352169  0.5084568876
##  [566]  0.9713083936  0.4650630975  0.6340118054  0.1591761963  1.0095874501
##  [571]  0.1525765096  0.2516678722  0.4675018948  0.2748485025  0.3418590116
##  [576]  0.1585948912  0.3965317641  0.3718423435  0.4889078212 -0.1840624508
##  [581] -0.2522990256  0.0486590424 -0.1126535769  0.5853848983  0.4530950613
##  [586]  0.2817100177  0.0542268090  0.5966050993 -0.3639805529  0.0461117118
##  [591]  1.2680247851 -0.0552850808  0.1425688949  1.3816805969 -0.2564237895
##  [596] -0.5493604846  1.0448629117  0.1935541997  0.7452663240  0.2960471745
##  [601]  0.8573311659 -0.1729776193  0.4409988555  0.4238382105 -0.4620777172
##  [606]  0.0711136634  0.0443696869  0.4145634174 -0.4090592488  0.4026632202
##  [611]  0.0680370181 -0.0367234804  0.4090520020  0.5995914129 -0.3384873705
##  [616]  1.4598829266  0.8108155969  0.7824219644  0.2912759229  0.2075940630
##  [621]  0.7157443189 -0.0231452631 -0.0090487675  0.4281748838  0.1414024990
##  [626]  1.2037634256  0.4516964005  1.7197978253  0.7904083215  0.5042312018
##  [631] -0.1129916675 -0.1619848399  1.2606549532  0.3795259870 -0.0059195137
##  [636] -0.0309616168 -0.5557171717  0.7185408780 -0.0023119775  0.7568306071
##  [641]  0.7615994798  0.0658889980  0.3746538870 -0.0751112733  0.0915416828
##  [646]  1.2583595411  0.5919095544  0.2144170192  0.7819963983  0.9490491548
##  [651]  1.2466557324  0.5356758756  1.5264660980  0.2067589011 -0.0118933930
##  [656]  0.1834839362  0.7472143989  0.6212338927  0.4977574575 -0.0874864610
##  [661] -0.2937860724  0.1698418583 -0.1196323810  0.4128426355  1.4866382835
##  [666]  0.7331969062 -0.3374770160  0.7878298318 -0.1305043854 -0.1958242198
##  [671]  1.4429341572  0.9266201479  0.7045971478  0.1983070050  0.1770158370
##  [676]  1.0900642576 -0.5226323483  0.7694741644 -0.0211532658  1.3414480266
##  [681]  0.8844789920 -0.0239527488  1.2585480498  0.5887293419  0.8050516299
##  [686]  0.6518481351  0.7854790961  0.2923397631  0.1175786844  0.5481513114
##  [691]  0.1747221129  0.0200663365  0.1530373339  0.9034296043  0.8159458207
##  [696]  0.0148448460  0.4362529763  0.2301839623  0.4581904022  0.9895118684
##  [701]  0.4467410581 -0.1770718534  0.4567732474 -0.6590287095 -1.2396903574
##  [706]  0.2665764284  0.1430408932  0.0602078348  0.8057576589 -0.1888963755
##  [711]  0.3058360323 -0.1962718592  0.5461389738  0.6821394753  0.4109004580
##  [716] -0.1064204350  0.1857602360  0.1384209384  0.1396634189  0.6057176373
##  [721]  1.3312718111  0.1265767040  0.7321337774  0.0519078120  0.0534984471
##  [726]  1.3090531944  0.2048138023  0.8779923445  0.4209677549  0.3400318344
##  [731]  0.7945918132  0.3607379356 -0.8184529760 -0.8581854041 -0.1358394085
##  [736] -0.1346819766 -0.1724862924  0.1447378247 -0.0667784823 -0.1796483281
##  [741] -0.1492183323  0.4578838934  0.8096491256 -0.1940433557  0.6730117119
##  [746] -0.4719489724  1.2583003236  0.9172639321 -0.5380717090  0.0277752126
##  [751]  0.7350051385  0.3081253271  0.6724998889  0.3445373134  0.6812144750
##  [756]  0.5047718407  0.4443042683 -0.1675922831  0.6827005031  0.4616552587
##  [761]  1.2583595411  0.3445787938  0.5816523802  0.4300282765  0.8521683571
##  [766]  0.9738792497  0.1826522603  0.2895952250  0.4424156475  0.0918895708
##  [771]  0.5251733855  1.9093560848  0.8076936814  0.2190934701  0.7648557319
##  [776]  0.2917077110 -1.0209060125  0.1411439492  0.0359501293 -0.2379623748
##  [781]  0.0726879129  0.5229240734 -0.1671653371  0.2682684718  0.6060485385
##  [786]  0.9720332572 -0.1073103232  0.5895771647  0.9102813096  1.0528052518
##  [791] -0.7364809709  0.7003627542  0.5571596196 -0.1560227838  1.4675898551
##  [796]  0.6057969149  1.1100313396  0.7595424730  0.4032184468  0.8322994829
##  [801]  0.9230439214  0.5567470262  0.5034309613  0.5329944820  1.1605753929
##  [806]  0.1471665728  1.3606514685  0.4762247622  0.7866112379 -0.0473899129
##  [811]  0.5928241024  0.8920801379  1.1332269645  0.7298178622  0.1847727768
##  [816]  0.3574709316 -0.0975506905  0.0880218381  0.5471277409  0.1864589238
##  [821] -0.2369295345  0.6327221994  0.1592500606 -0.0938713605  0.7742372233
##  [826]  0.3472058878  0.2429825708  0.1797369527 -0.0012783526 -0.0741494919
##  [831]  0.4314903933  0.3830120059 -0.2283526107  0.3169473472 -0.4295916490
##  [836]  1.0433480677  0.5047121571 -0.2931357262 -0.0270075817  1.1665988841
##  [841]  2.1781909076  0.3184108242 -0.8541395685  0.2088557270  0.1769733909
##  [846]  0.1299189709 -0.2489959039  0.1672894883 -0.6950507589 -0.1752740366
##  [851]  0.1048984867  0.0727901767  0.1252312776  0.3711762200  1.4527748468
##  [856]  0.2701439963  0.7657924601  0.9972258243  0.7917799943  0.4659995837
##  [861]  0.2702625721  1.2178008905  0.4034957092  0.7821102890  0.0123613926
##  [866] -0.4216602748  0.9285562981  0.4847996385  0.2592072461  0.0514998650
##  [871] -1.1706096890  0.5131096050  0.4616552587  0.6369593192  0.4440167905
##  [876] -0.4247224627 -0.1482575159  0.5908757118  1.3540666030  0.1125641890
##  [881]  0.0939000548  0.7618814367  0.3245614099  0.7313194283  0.1076254307
##  [886]  0.3868855627 -0.1039797293  0.8429555741  0.0970233341  0.5256836081
##  [891]  0.4491258582 -0.2261444344  0.3592436212  0.1115966427  0.0603251752
##  [896]  0.7923416069  0.1479125194 -0.7477679725 -0.2917266711  0.4191150105
##  [901]  0.2830815392  1.2625429988  0.9789329842  0.6079471803  0.2570626309
##  [906]  1.3581804141  0.1149841959  0.3579273740  0.4995303049  0.8649143838
##  [911]  1.0065387023  0.0080856980 -0.4230314743 -0.1428501139 -0.3373888469
##  [916]  0.9809968083  1.0758460400  1.3971309452  0.7896172795  0.3149501246
##  [921] -0.0969010692  0.0916676113 -0.5623762832  0.2591883679 -0.3335863234
##  [926]  2.3251823668 -0.1510524261  0.8744466025  0.8777875267  1.0051714657
##  [931]  0.2302613599  1.1498927184  0.7579389963 -0.3471652241  0.4330022598
##  [936]  0.7268395756 -0.3501641904  1.2159675272  0.5165140386  0.5293775953
##  [941]  1.5197147205  1.2387303220 -0.4479183201  1.3473982082  0.4277616058
##  [946]  0.7823782823  0.3525023692  0.3872450227  0.2166948299  1.3336331979
##  [951]  0.6011241479  0.0623652468  0.1018324038  0.1789711411  0.6114843859
##  [956]  0.8747616171  0.6800753392  0.1178971611 -0.5782028634  0.5808939231
##  [961]  1.2419316719 -1.2324546550  0.3823378716  0.1059889418  0.2315104525
##  [966] -0.0690060967  0.3934681652 -0.5639096749  0.8257066293  0.1945892241
##  [971]  0.7300330552  0.5265188292  1.4306837834  0.5525404038  0.5717571462
##  [976]  0.5763850486  0.5003997841 -0.1222703743  0.3605547408  0.2306307450
##  [981]  1.5935680062  1.1130568629  0.3989042538  0.3749090238  0.7488972581
##  [986]  0.9042401524  0.6124725170  0.6902530590  0.5715593821  0.8364587347
##  [991] -0.1142452659  0.1174342816  0.8155486697  0.3601408790  0.0645667361
##  [996] -0.1827948628  0.6563364080  0.3820303250  1.0708862914  0.0019459211
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
##   0.75685047   0.43681814 
##  (0.13813402) (0.09767261)
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
## [1]  0.5845668 -0.1619391 -0.2621568 -0.7591021 -1.0878519 -0.4113954
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
## [1] 0.0086
```

```r
sd(results2$thetastar)  #the standard deviation of the theta stars is the SE of the statistic (in this case, the mean)
```

```
## [1] 0.9093145
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
## t1*      4.5 -0.03553554   0.9121362
```

One of the main advantages to the 'boot' package over the 'bootstrap' package is the nicer formatting of the output.

Going back to our original code, lets see how we could reproduce all of these numbers:


```r
table(sample(x,size=length(x),replace=T))
```

```
## 
## 0 2 3 7 8 9 
## 1 1 3 1 2 2
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
## [1] 0.0211
```

```r
se.boot
```

```
## [1] 0.8869676
```

Why do our numbers not agree exactly with those of the boot package? This is because our estimates of bias and standard error are just estimates, and they carry with them their own uncertainties. That is one of the reasons we might bother doing jackknife-after-bootstrap.

The 'boot' package has a LOT of functionality. If we have time, we will come back to some of these more complex functions later in the semester as we cover topics like regression and glm.

