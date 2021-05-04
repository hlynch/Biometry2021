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
## 0 1 3 5 7 9 
## 2 1 3 1 2 1
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
## [1] -0.0523
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
## [1] 2.682932
```

```r
UL.boot
```

```
## [1] 6.212468
```

Method #2: Simply take the quantiles of the bootstrap statistics


```r
quantile(xmeans,c(0.025,0.975))
```

```
##  2.5% 97.5% 
##   2.6   6.1
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
##    [1] 4.3 4.8 5.2 5.4 4.7 3.7 5.4 5.7 5.3 4.2 2.8 4.4 3.9 4.7 3.6 4.0 4.4 4.3
##   [19] 4.8 4.3 4.4 4.2 4.1 4.8 3.4 4.7 4.7 4.8 3.2 4.3 4.6 5.6 5.8 4.5 5.5 6.0
##   [37] 4.2 4.5 4.1 4.6 4.0 4.2 3.9 6.1 5.1 3.9 3.2 3.3 5.0 5.9 2.3 4.7 7.4 3.4
##   [55] 4.7 5.6 4.7 5.5 5.3 4.3 6.3 2.8 4.4 3.6 5.4 4.8 4.8 4.4 4.3 3.4 4.4 4.4
##   [73] 4.5 3.8 4.0 5.4 5.3 4.8 4.6 5.1 3.8 5.4 4.2 4.6 3.9 4.6 5.3 5.1 3.7 5.3
##   [91] 5.4 4.8 3.6 2.0 3.7 5.7 3.9 4.4 4.3 5.1 6.7 2.5 5.4 3.0 4.4 4.9 4.1 4.6
##  [109] 5.6 4.9 4.3 2.7 5.5 5.8 3.7 3.1 4.5 4.4 3.9 5.0 4.7 4.6 5.1 5.3 5.1 5.3
##  [127] 4.7 4.9 4.0 5.2 4.7 4.7 2.9 4.0 2.9 4.1 4.6 4.9 4.5 3.9 4.5 4.9 5.0 5.0
##  [145] 4.2 4.2 2.6 5.2 3.9 3.3 4.2 4.2 4.9 5.6 4.6 5.4 4.7 4.8 4.8 2.8 3.0 3.5
##  [163] 4.8 5.7 4.6 5.3 4.1 5.9 4.8 2.5 4.8 5.5 3.8 3.8 6.0 4.5 5.1 4.0 5.2 5.9
##  [181] 5.4 4.7 4.4 4.7 4.2 4.1 4.9 2.5 4.3 4.3 3.9 4.3 3.7 5.8 3.6 4.1 4.5 5.9
##  [199] 3.6 5.8 3.7 5.1 2.9 3.7 5.5 3.1 3.2 4.6 3.2 4.1 5.7 4.8 5.0 5.3 4.3 3.6
##  [217] 3.7 3.6 4.4 3.4 4.9 3.8 4.5 5.2 4.2 4.0 5.5 5.3 3.4 5.5 4.9 4.6 4.0 3.7
##  [235] 5.4 5.9 5.1 5.0 2.8 3.5 6.4 3.7 5.4 3.9 5.0 4.8 4.6 4.6 5.2 4.6 5.8 4.1
##  [253] 5.6 5.3 4.0 3.9 4.7 4.1 4.3 4.5 4.5 3.9 4.5 4.3 6.5 4.6 4.1 3.3 5.2 5.8
##  [271] 3.3 3.3 4.0 4.2 3.9 5.4 5.0 3.7 4.3 6.2 5.6 3.8 4.4 4.3 4.5 5.2 4.4 2.8
##  [289] 4.4 5.1 5.1 3.9 5.0 5.9 3.7 6.4 2.8 4.3 5.6 4.4 3.2 2.4 3.7 5.1 5.0 5.7
##  [307] 4.9 6.5 5.1 4.9 6.0 5.3 5.6 2.8 4.4 4.9 4.8 6.1 4.5 2.8 4.1 4.0 3.3 6.5
##  [325] 5.6 3.7 2.7 2.9 3.9 3.8 5.3 4.6 4.5 4.2 6.2 3.2 5.9 4.2 4.6 5.0 4.2 4.8
##  [343] 5.0 3.9 5.4 5.7 6.5 3.3 4.1 4.3 4.3 3.4 4.5 5.5 4.7 5.2 3.7 3.7 4.8 4.4
##  [361] 6.1 5.6 4.0 4.4 3.9 4.1 5.4 4.9 5.0 4.0 4.2 4.0 5.1 3.3 4.8 3.3 4.4 5.6
##  [379] 5.7 4.5 4.3 2.9 4.5 5.0 4.3 2.9 5.3 4.9 4.5 4.2 5.2 3.6 5.0 5.4 4.8 5.4
##  [397] 4.1 3.0 4.1 2.9 6.4 5.9 4.3 5.7 3.1 4.7 3.7 4.1 4.7 3.4 4.7 4.1 3.8 4.4
##  [415] 5.1 5.3 4.5 4.2 5.2 5.3 4.7 3.8 4.8 4.9 4.9 4.7 4.8 2.7 4.8 3.1 4.1 4.8
##  [433] 5.5 5.4 4.5 4.9 4.3 3.0 4.0 3.0 5.6 5.1 4.1 5.1 3.9 2.1 2.9 3.1 4.1 6.3
##  [451] 4.7 3.9 4.9 4.3 5.2 4.7 3.6 5.9 3.7 4.5 4.2 3.8 3.1 5.6 5.9 3.5 5.0 3.1
##  [469] 4.1 5.2 4.3 7.0 4.8 3.8 6.1 2.5 5.3 5.0 3.0 5.1 3.6 4.1 4.7 4.1 5.9 5.4
##  [487] 4.5 3.9 5.2 5.3 4.9 5.0 4.9 5.9 4.3 3.9 5.0 3.6 3.9 4.0 4.8 3.9 3.6 4.2
##  [505] 5.7 4.0 4.4 4.3 5.2 3.4 5.0 4.6 2.4 6.1 6.1 5.1 3.7 3.9 5.1 5.1 3.6 3.9
##  [523] 4.5 5.7 5.0 2.8 5.4 3.7 4.5 3.5 5.0 4.3 4.3 3.2 4.5 5.2 4.9 4.1 3.4 5.2
##  [541] 4.7 5.2 4.1 4.5 4.4 5.9 4.3 4.1 3.8 4.6 4.5 5.1 5.4 5.5 5.2 3.9 5.3 3.8
##  [559] 4.9 2.2 3.9 5.5 4.6 2.4 5.3 4.4 4.7 4.5 4.3 3.5 3.2 3.3 3.5 5.9 4.1 5.4
##  [577] 5.1 5.3 3.4 5.0 4.4 6.1 2.3 4.3 4.1 3.6 4.4 3.4 6.4 4.7 4.1 4.8 5.1 2.9
##  [595] 3.9 6.0 4.0 4.8 6.1 3.9 3.7 4.5 4.8 3.8 5.5 2.9 5.8 5.4 5.8 5.7 2.2 4.9
##  [613] 5.9 4.2 5.6 3.0 4.8 3.9 6.3 4.7 5.5 5.1 5.2 4.5 4.5 3.1 4.9 6.0 4.2 4.1
##  [631] 3.4 3.8 2.9 5.3 4.1 5.4 5.7 3.4 2.0 4.5 4.2 5.7 4.0 6.1 5.0 5.0 5.2 5.0
##  [649] 3.8 4.6 4.0 2.5 5.7 5.6 4.7 5.5 4.7 4.4 3.4 3.8 4.1 5.1 4.8 6.6 4.8 4.0
##  [667] 4.7 3.8 4.6 5.7 3.3 3.9 4.4 4.4 6.3 5.0 4.6 5.7 4.5 4.5 5.0 5.9 5.4 4.5
##  [685] 4.5 6.3 5.6 5.1 4.5 5.1 3.8 4.9 3.5 3.2 3.7 5.1 3.6 4.5 4.2 4.9 3.9 3.2
##  [703] 5.5 4.4 4.1 4.7 4.9 4.4 3.8 4.8 4.9 4.5 5.0 5.5 5.4 5.6 4.8 5.2 4.7 4.4
##  [721] 4.2 3.7 5.5 7.4 5.8 4.6 4.8 5.0 5.9 5.7 4.5 3.9 4.7 3.7 4.3 4.8 3.4 4.1
##  [739] 4.3 6.0 4.0 3.0 6.4 4.4 4.4 4.2 4.9 2.1 6.6 3.6 5.3 5.2 6.5 5.3 4.1 5.2
##  [757] 4.1 4.3 4.0 4.7 6.2 4.1 3.3 3.1 3.8 2.5 5.4 4.4 4.3 4.8 3.8 6.4 4.4 4.1
##  [775] 6.1 5.2 4.9 3.5 5.4 3.5 4.6 5.0 3.1 3.3 5.5 4.9 4.3 2.4 3.8 3.9 5.5 5.7
##  [793] 5.1 2.8 4.5 5.2 5.2 6.2 2.8 5.3 4.6 4.0 4.0 5.5 4.9 4.5 4.3 4.1 4.4 4.8
##  [811] 3.4 4.2 4.5 3.9 4.3 5.7 5.3 6.9 3.8 3.9 3.7 4.5 4.6 4.8 5.1 4.4 3.7 3.6
##  [829] 4.9 4.6 4.9 2.5 5.0 5.0 4.4 4.5 4.2 2.7 6.1 4.0 4.2 4.5 6.4 4.5 4.8 4.3
##  [847] 5.7 4.1 6.1 6.3 5.0 5.4 4.5 4.8 6.2 5.0 5.3 3.8 5.1 2.9 4.8 3.5 3.4 5.6
##  [865] 4.5 4.5 4.5 5.1 3.6 3.5 4.3 4.9 5.9 5.9 4.9 5.3 4.3 4.1 2.8 5.8 3.8 4.0
##  [883] 3.6 4.9 5.1 3.0 3.4 5.4 4.6 4.5 5.5 4.3 4.9 5.3 4.4 4.4 3.7 4.3 4.2 4.8
##  [901] 3.8 4.4 6.0 4.3 4.5 4.5 5.9 4.4 4.2 4.0 3.4 5.3 4.1 5.3 4.4 4.3 3.8 5.6
##  [919] 2.8 5.4 3.3 4.2 5.8 4.4 4.2 5.1 4.8 5.0 3.6 4.7 3.7 5.9 4.2 2.9 3.5 5.1
##  [937] 4.5 4.9 5.7 3.3 5.3 2.7 2.1 3.2 4.6 2.7 4.4 4.9 4.8 5.0 5.7 3.1 4.0 4.0
##  [955] 2.4 4.2 3.7 4.7 2.9 6.1 4.1 5.0 4.7 5.5 4.6 4.4 4.9 2.7 4.9 3.3 4.4 4.1
##  [973] 2.9 5.2 4.9 5.0 4.5 4.9 4.6 5.2 2.7 4.9 5.7 4.4 3.0 5.9 4.7 5.3 4.9 3.9
##  [991] 3.8 5.1 3.5 2.3 3.9 5.2 4.4 4.3 4.8 4.9
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
##    [1] 4.0 4.6 3.6 5.3 4.1 4.4 3.0 4.8 6.6 5.5 4.3 3.7 4.4 5.1 5.5 5.3 3.8 4.3
##   [19] 3.9 3.9 3.6 3.9 4.2 4.2 4.5 2.4 4.8 4.3 3.2 5.4 4.5 4.0 5.9 4.2 4.4 5.7
##   [37] 4.4 2.5 4.0 4.9 3.5 4.3 2.9 4.4 3.5 2.5 5.6 3.9 4.3 2.6 5.9 4.1 5.2 4.4
##   [55] 4.7 3.3 4.4 5.6 4.3 4.5 5.8 4.6 5.3 4.9 5.1 4.2 5.2 3.7 3.0 6.5 5.2 3.5
##   [73] 4.8 3.8 4.9 4.8 3.4 5.2 4.4 4.7 4.1 3.8 4.5 4.6 6.2 3.0 5.5 5.3 3.9 4.8
##   [91] 5.1 1.9 3.8 6.5 6.0 4.5 4.2 4.0 4.2 4.4 4.4 3.3 4.2 2.7 6.7 5.1 4.3 3.6
##  [109] 3.9 3.8 4.4 5.6 5.0 3.6 3.7 2.7 3.6 3.4 5.5 5.4 4.6 4.8 5.8 4.4 5.1 5.2
##  [127] 3.0 4.4 5.1 4.1 4.3 5.3 4.9 5.4 5.3 4.2 3.0 3.7 3.5 4.7 4.6 6.0 5.4 5.2
##  [145] 4.1 3.1 4.2 3.1 4.4 3.6 2.2 4.6 4.4 5.4 5.0 5.8 3.6 5.9 4.4 4.8 4.4 4.2
##  [163] 5.0 4.1 3.4 3.1 5.5 4.8 5.1 4.5 4.1 7.5 4.2 5.9 6.2 5.2 3.1 4.2 4.9 4.7
##  [181] 4.4 3.8 4.9 3.9 2.8 3.2 2.9 3.6 5.3 3.2 3.8 5.0 3.6 5.2 4.3 4.6 3.4 5.8
##  [199] 4.4 4.8 5.3 5.5 4.1 3.3 3.1 2.6 5.4 3.9 4.6 3.6 4.1 4.0 5.7 5.5 5.9 5.8
##  [217] 5.2 5.5 3.6 4.9 4.4 3.6 4.8 5.2 3.2 4.2 5.3 4.4 5.0 4.9 4.4 4.0 5.3 6.3
##  [235] 5.6 3.5 5.1 4.4 4.2 3.8 4.3 5.1 6.0 5.3 3.9 4.9 3.8 5.5 5.3 4.8 5.3 4.1
##  [253] 3.8 3.2 4.7 4.3 5.7 2.8 4.5 6.1 4.4 5.1 6.8 4.9 3.2 4.4 4.6 4.8 4.7 2.4
##  [271] 5.2 5.6 6.4 3.8 4.8 4.5 4.2 4.9 4.3 5.6 4.0 4.9 6.1 4.0 4.2 4.9 4.0 3.9
##  [289] 3.9 4.4 4.4 5.3 3.4 3.8 6.2 5.7 4.1 3.5 3.6 6.0 3.3 4.7 5.1 4.4 3.9 3.3
##  [307] 3.8 4.1 5.7 5.1 4.7 4.3 4.3 6.0 4.9 4.0 2.4 3.5 4.3 5.4 4.8 5.0 4.0 4.9
##  [325] 4.5 3.5 4.2 4.3 3.8 4.3 6.8 5.5 4.4 3.7 3.4 5.7 5.0 3.9 5.2 4.6 5.1 4.1
##  [343] 4.4 6.1 5.1 4.2 4.2 5.6 3.8 4.5 4.4 3.7 5.9 6.7 4.2 4.2 3.3 5.5 6.1 3.9
##  [361] 4.8 5.9 5.0 6.5 4.6 4.8 3.4 5.5 3.4 4.0 3.7 3.5 4.5 4.1 4.9 5.1 3.9 3.7
##  [379] 5.1 3.8 5.7 3.7 5.3 5.8 3.4 5.0 7.4 4.6 5.4 4.5 3.8 5.3 4.3 3.8 6.2 4.0
##  [397] 4.9 4.4 4.4 5.3 3.8 3.4 4.5 2.5 4.8 5.6 5.1 5.4 3.7 5.5 4.4 6.1 3.2 4.5
##  [415] 3.5 5.2 3.8 6.1 3.8 5.1 3.3 3.7 5.2 6.2 4.4 4.7 3.6 5.9 4.0 4.1 4.0 4.7
##  [433] 6.0 5.0 4.6 3.7 3.2 5.0 5.6 5.2 3.7 3.9 3.6 5.5 4.4 3.8 4.1 4.3 4.5 5.4
##  [451] 3.6 3.4 5.1 3.4 4.9 4.4 2.4 3.6 4.1 4.4 5.1 4.1 4.5 5.4 4.9 5.5 5.7 6.4
##  [469] 4.9 5.0 2.9 3.6 4.0 5.7 5.5 4.4 3.8 4.2 6.6 3.9 4.3 3.8 3.1 4.3 4.6 5.1
##  [487] 4.6 4.3 3.3 4.7 3.9 3.2 4.0 5.4 5.2 2.9 5.0 5.0 2.2 5.1 3.7 3.5 5.0 4.5
##  [505] 5.8 5.8 4.2 3.2 5.0 3.8 3.8 3.8 5.2 4.7 3.8 3.5 6.9 5.9 4.0 4.8 4.4 4.3
##  [523] 3.7 3.8 2.9 5.4 5.4 4.6 5.6 5.3 4.3 4.9 5.8 3.4 4.4 5.8 4.8 4.8 4.1 3.0
##  [541] 4.9 6.2 3.2 4.0 5.2 4.4 3.7 2.7 3.9 4.1 3.0 4.1 3.1 3.6 4.7 3.9 4.8 3.7
##  [559] 2.5 5.1 5.8 3.6 4.2 4.7 5.1 6.6 4.1 3.3 4.1 5.4 5.7 4.8 4.1 4.1 4.9 4.0
##  [577] 3.4 4.5 3.5 4.2 3.4 3.1 5.3 3.7 3.8 6.8 5.0 4.7 5.3 5.0 6.8 4.2 3.9 3.8
##  [595] 4.0 4.1 4.8 4.3 4.0 5.5 4.3 4.8 5.0 4.8 5.7 4.8 4.8 2.8 5.6 4.4 4.0 5.7
##  [613] 4.6 5.0 5.7 4.4 4.0 6.2 4.4 3.6 5.1 2.6 5.0 5.3 3.1 5.2 5.2 5.4 3.7 5.2
##  [631] 3.7 5.1 5.5 3.6 5.4 2.9 4.7 4.6 4.9 4.6 4.7 5.1 5.1 4.6 4.7 4.2 4.8 3.5
##  [649] 5.2 5.1 5.6 3.4 5.0 5.0 6.3 3.6 5.8 5.8 4.2 5.0 3.5 4.7 4.4 6.0 5.7 3.1
##  [667] 3.6 5.8 5.9 3.4 5.2 4.4 4.6 4.1 3.2 3.1 3.3 3.9 3.7 4.7 2.2 4.7 4.4 3.9
##  [685] 6.6 6.0 5.5 4.8 5.1 5.3 4.8 5.3 3.3 3.5 3.5 5.1 4.3 4.6 4.0 5.7 4.9 3.4
##  [703] 3.5 4.6 4.5 5.2 5.1 4.3 5.0 4.3 5.5 5.8 6.1 5.1 4.7 5.3 5.1 1.9 5.4 3.8
##  [721] 4.9 5.2 2.9 4.1 5.5 5.8 5.3 5.2 5.2 3.8 5.2 5.8 6.1 4.0 5.4 4.5 4.1 3.7
##  [739] 3.8 4.0 4.7 5.1 3.7 4.4 4.1 5.4 3.8 4.7 4.3 5.4 4.2 5.8 4.9 4.3 4.3 6.0
##  [757] 4.0 3.4 5.9 4.7 6.9 4.6 5.1 4.9 4.9 3.6 5.4 3.5 3.1 4.6 3.6 5.1 4.7 3.8
##  [775] 3.8 5.1 5.8 5.3 5.3 3.4 3.9 3.9 2.9 4.4 5.9 3.5 4.2 4.3 4.8 4.5 4.3 4.4
##  [793] 3.7 6.0 3.4 4.7 4.4 2.5 5.6 5.6 5.0 3.9 4.8 5.8 6.5 4.9 4.6 4.8 3.6 4.5
##  [811] 4.0 3.1 5.8 4.4 3.6 3.6 4.5 7.3 5.8 4.2 6.3 6.8 4.5 2.0 5.2 6.6 5.7 4.5
##  [829] 4.6 4.3 5.3 4.1 5.0 3.4 4.7 5.8 4.2 3.2 6.6 3.3 6.4 4.4 5.7 5.8 3.6 4.7
##  [847] 4.8 5.2 3.3 5.1 3.5 5.5 6.2 6.1 4.2 4.5 4.1 4.4 4.2 5.2 4.9 4.7 5.1 3.7
##  [865] 4.4 3.4 4.9 5.1 5.5 3.8 5.2 4.1 3.6 4.8 4.5 4.5 5.2 3.9 5.3 5.2 3.5 2.2
##  [883] 4.0 5.4 4.0 3.8 6.2 4.2 3.3 4.7 3.6 5.3 3.4 3.0 4.0 3.9 4.3 3.8 3.4 4.3
##  [901] 5.2 4.1 5.2 7.0 3.5 5.2 3.8 4.1 3.1 4.5 5.9 5.9 2.9 3.9 5.7 3.9 4.7 3.0
##  [919] 2.4 4.3 4.1 5.7 5.1 4.7 2.8 4.6 2.4 2.8 2.6 4.7 5.3 5.4 3.8 2.7 4.0 5.2
##  [937] 3.9 5.5 4.9 5.4 4.8 5.1 5.2 4.8 4.5 4.1 4.2 5.5 3.5 4.5 3.8 4.2 4.4 4.1
##  [955] 5.8 5.3 2.7 4.7 4.9 4.8 4.4 4.5 4.0 5.9 5.4 3.9 3.8 4.8 4.3 3.9 6.0 3.9
##  [973] 6.0 4.5 4.6 3.4 4.5 4.7 4.9 4.0 3.5 4.5 5.0 4.4 2.0 3.8 4.2 4.5 4.2 4.1
##  [991] 4.7 4.9 3.7 3.3 4.9 6.8 5.4 3.3 4.2 3.5
## 
## $func.thetastar
## [1] 0.0185
## 
## $jack.boot.val
##  [1]  0.54424779  0.43784530  0.31590214  0.24986595  0.09485714 -0.07939394
##  [7] -0.20434783 -0.23693182 -0.37848837 -0.53414634
## 
## $jack.boot.se
## [1] 1.033811
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
##    [1] 3.9 3.0 4.5 3.3 3.9 3.1 3.3 2.3 4.5 4.7 5.4 5.5 4.5 4.0 3.8 4.7 3.3 3.6
##   [19] 2.2 4.4 4.3 3.5 5.7 5.5 4.6 4.3 4.3 4.8 4.1 3.2 4.5 4.4 4.1 5.0 3.7 6.5
##   [37] 4.9 4.9 3.9 4.4 5.6 3.1 5.8 3.8 3.4 4.0 1.9 4.3 5.4 5.4 5.3 6.4 4.2 4.8
##   [55] 4.6 6.3 5.4 4.8 3.5 4.1 4.5 5.2 4.1 4.8 4.3 3.6 3.4 4.3 6.1 5.1 5.0 5.9
##   [73] 5.1 5.9 4.2 3.4 5.4 3.7 4.4 3.9 2.3 4.5 5.1 2.1 3.9 4.3 6.2 5.7 5.2 4.4
##   [91] 5.2 3.5 3.5 3.8 3.8 3.8 4.5 3.7 3.3 4.6 5.0 4.8 4.1 5.1 5.2 4.9 5.0 4.9
##  [109] 3.2 6.3 3.9 5.4 4.5 4.5 4.3 5.2 5.1 4.0 2.7 6.1 5.4 5.2 3.7 4.3 4.6 5.1
##  [127] 5.7 3.7 4.4 4.3 5.1 4.1 2.8 4.5 6.0 5.2 2.3 4.5 4.9 4.3 5.0 3.7 3.9 3.7
##  [145] 4.5 4.7 4.3 4.1 3.9 4.4 5.1 4.0 4.4 6.1 4.0 3.4 4.7 3.9 3.9 4.4 4.5 4.3
##  [163] 4.7 3.4 5.3 4.3 4.7 4.5 4.2 4.3 5.2 5.3 3.0 4.2 4.9 4.1 3.4 4.6 5.7 2.5
##  [181] 4.4 6.0 4.1 2.9 4.2 6.4 4.7 4.6 4.6 5.1 4.3 4.4 3.7 4.2 4.3 4.6 4.1 3.4
##  [199] 4.5 3.5 3.2 4.3 4.2 5.6 3.2 4.2 4.2 5.2 4.6 4.6 4.5 4.0 4.9 4.0 5.0 3.8
##  [217] 3.9 5.0 4.3 3.7 5.2 5.4 5.8 3.5 5.2 4.3 5.0 4.9 3.3 4.0 2.7 3.4 2.2 4.4
##  [235] 5.2 5.7 4.6 5.5 4.5 6.0 3.3 4.1 4.0 5.0 5.5 4.1 3.9 4.9 4.8 6.8 4.8 2.6
##  [253] 4.5 3.1 4.5 3.3 5.1 4.4 4.6 4.6 3.5 5.6 2.9 6.1 6.3 4.8 6.8 3.4 3.5 5.1
##  [271] 7.1 5.5 4.8 4.9 5.6 4.5 2.9 4.5 4.7 3.5 4.8 5.8 5.9 2.8 3.4 4.3 4.3 3.9
##  [289] 4.4 3.4 4.1 3.8 5.0 4.9 4.3 4.3 4.7 4.6 2.5 4.5 4.3 3.5 2.9 4.0 4.4 4.8
##  [307] 4.8 3.7 4.5 3.5 6.4 4.3 4.8 5.7 4.7 4.3 5.3 5.1 4.5 4.3 4.4 3.0 5.6 3.9
##  [325] 4.1 4.9 4.1 4.0 3.5 4.5 5.6 5.9 5.0 2.9 3.3 4.9 2.8 4.5 4.8 4.8 5.6 4.3
##  [343] 5.3 4.8 2.7 5.1 4.6 4.4 3.7 6.6 5.4 4.2 3.4 4.9 4.4 3.9 4.4 3.7 3.8 3.9
##  [361] 4.1 5.1 5.3 5.9 5.1 3.8 5.0 4.7 5.0 5.4 4.2 4.4 5.3 2.9 3.0 2.5 4.9 4.3
##  [379] 4.0 4.9 4.5 5.1 4.1 4.2 6.1 4.3 4.3 4.3 4.5 4.5 4.7 3.9 3.3 5.0 6.9 4.9
##  [397] 2.5 4.3 6.3 4.5 5.1 3.3 4.2 4.3 4.7 3.7 4.8 3.7 3.9 3.7 5.1 4.3 5.1 4.2
##  [415] 3.7 4.2 5.6 5.2 5.1 5.3 4.8 4.1 5.4 5.2 5.4 3.6 3.2 3.8 4.6 4.5 4.3 5.0
##  [433] 4.0 3.4 4.2 5.7 2.4 2.3 5.2 5.0 3.4 3.7 5.1 3.7 5.9 3.5 4.2 4.9 4.6 3.6
##  [451] 4.1 4.5 5.4 4.7 4.6 5.2 3.7 4.5 5.3 6.2 3.6 5.3 3.4 3.5 3.5 5.3 5.4 4.9
##  [469] 3.2 4.3 4.9 5.0 4.9 4.7 3.2 3.9 4.6 3.4 5.4 5.0 3.3 4.8 5.5 4.1 5.8 3.4
##  [487] 4.5 6.2 4.5 4.4 4.2 5.4 3.3 5.2 3.3 6.7 4.4 5.3 4.2 4.3 3.3 7.1 3.1 4.4
##  [505] 5.4 4.7 5.4 4.0 4.1 5.4 5.2 4.4 4.5 3.4 4.7 5.3 4.7 4.4 5.0 3.5 2.8 6.6
##  [523] 4.5 3.7 4.1 5.0 5.3 4.9 4.5 4.2 5.4 6.0 3.0 4.1 2.9 3.8 4.4 4.7 4.0 4.9
##  [541] 5.0 4.2 3.7 4.6 5.5 4.2 4.2 3.4 4.8 4.5 4.1 4.5 4.3 3.4 4.0 2.7 4.2 3.0
##  [559] 5.4 4.6 4.7 5.0 4.7 3.6 4.3 4.3 3.5 4.2 5.0 4.5 4.0 6.6 3.4 5.4 3.0 3.2
##  [577] 3.3 4.3 4.8 5.8 4.1 3.7 4.1 3.3 3.5 5.9 5.7 3.2 4.5 7.1 3.6 4.4 4.7 4.3
##  [595] 6.4 4.4 3.6 3.3 3.4 4.8 5.7 4.1 5.7 4.6 4.3 4.5 5.0 4.8 4.0 4.1 4.1 5.0
##  [613] 3.5 6.7 4.9 4.3 4.5 3.7 4.6 4.8 5.0 4.3 4.3 3.9 3.6 4.8 3.2 3.5 5.4 3.7
##  [631] 4.1 4.6 6.0 3.6 4.0 5.1 4.2 2.4 4.6 2.8 3.9 4.0 4.7 4.3 3.5 3.8 3.4 3.6
##  [649] 5.7 4.5 4.7 5.6 4.7 4.0 4.8 3.3 5.9 3.6 3.9 7.0 4.4 3.2 3.2 5.5 4.6 4.2
##  [667] 4.0 5.4 4.1 4.1 5.6 5.1 4.2 4.6 4.1 4.1 4.4 4.0 5.2 3.5 3.1 3.7 4.0 3.9
##  [685] 4.7 3.7 4.1 4.0 4.3 4.9 3.9 3.6 5.5 4.8 5.8 5.2 3.4 3.3 3.9 6.3 4.6 6.6
##  [703] 4.4 3.5 3.2 4.1 4.5 6.1 4.1 5.1 4.8 4.3 4.4 4.5 4.8 4.4 3.6 5.3 4.6 5.7
##  [721] 4.8 4.7 4.1 4.6 4.9 4.5 5.7 4.7 3.9 3.5 4.4 3.4 4.8 3.3 4.7 2.9 3.9 3.5
##  [739] 5.0 5.3 4.0 3.0 2.7 3.9 5.0 4.4 3.5 4.1 6.1 4.3 3.2 5.8 4.0 5.4 3.0 3.5
##  [757] 4.7 6.0 3.6 4.1 4.5 4.8 5.3 4.0 4.3 2.9 4.5 4.6 3.2 5.8 4.1 3.4 4.3 4.6
##  [775] 5.9 4.6 4.2 5.7 4.1 4.1 4.6 3.4 6.1 4.1 5.5 3.9 5.5 4.7 3.0 5.6 2.6 4.3
##  [793] 4.7 5.5 5.8 4.0 4.8 5.0 3.7 5.9 3.8 2.9 5.3 5.4 4.3 4.6 4.4 4.6 3.1 4.7
##  [811] 5.8 4.9 3.3 3.6 4.3 4.5 4.1 3.1 6.0 3.8 4.3 4.9 3.4 3.6 5.2 6.5 4.5 3.8
##  [829] 5.6 3.4 5.2 4.9 3.2 4.7 5.4 4.2 3.8 4.2 4.5 4.4 6.7 7.0 4.3 5.6 4.5 5.1
##  [847] 3.7 4.8 4.6 4.3 5.4 2.8 3.1 4.9 4.7 7.1 6.1 5.8 5.8 4.3 4.3 3.9 3.8 4.7
##  [865] 4.5 4.4 5.2 4.8 5.8 4.3 4.7 4.4 5.0 3.5 4.9 4.3 4.1 3.2 4.7 5.9 5.2 5.3
##  [883] 4.2 4.0 6.0 4.5 4.2 4.7 4.9 4.2 3.6 4.2 3.3 5.5 5.4 4.4 1.9 4.6 4.5 4.7
##  [901] 5.0 2.3 4.6 4.3 4.8 3.9 6.4 4.8 4.6 5.6 3.6 2.2 4.5 4.5 5.6 3.6 4.7 2.6
##  [919] 4.7 5.1 4.5 4.8 5.9 6.2 4.7 5.6 6.1 5.3 5.3 5.6 3.5 6.3 4.0 3.9 4.5 2.4
##  [937] 5.1 4.9 4.7 3.9 5.0 4.2 3.7 5.8 4.9 5.6 4.9 4.6 2.6 4.5 4.5 4.2 2.7 4.7
##  [955] 6.4 2.8 4.3 4.1 4.2 4.2 3.0 3.9 5.3 5.4 4.8 3.5 4.3 3.5 4.1 3.8 4.7 5.1
##  [973] 3.7 5.5 4.3 5.1 4.5 5.3 5.1 3.6 5.0 4.5 4.7 4.7 3.3 6.2 5.0 4.6 4.2 3.6
##  [991] 3.4 4.2 5.1 3.3 5.3 4.1 4.7 4.4 4.1 4.1
## 
## $func.thetastar
## 72% 
## 4.9 
## 
## $jack.boot.val
##  [1] 5.400 5.300 5.200 5.200 4.900 4.916 4.700 4.700 4.500 4.400
## 
## $jack.boot.se
## [1] 0.9838858
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
## [1] 1.31495
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
##   3.377903   5.243804 
##  (1.442279) (2.413924)
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
## [1]  0.3545775  0.5423127  0.5272525 -0.5293040  0.9249700  1.5014014
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
##    [1]  0.217463765  0.073747016  1.348751313  0.510836842  0.959290952
##    [6]  0.754371818  0.494889860 -0.227096484  0.892783045  0.729264421
##   [11] -0.278231532  1.340674089  0.103515529  0.422277694  0.762495722
##   [16]  0.349348005  0.503180420  1.541596897 -0.068536852  1.048048844
##   [21]  1.900235611  1.512391708  0.237860663 -0.486743278  1.655273799
##   [26]  0.178572825  1.248610928  1.272540316 -0.204943507  0.491731698
##   [31] -0.248747538  1.033404270  1.662017977  1.468565105  1.435807916
##   [36]  1.089824775  1.004088749  0.246538015  1.319283713  1.548709486
##   [41]  0.113119588 -0.177906927  1.146023103  1.446904544  0.511108922
##   [46]  1.249291795  1.344968425  1.226390531  1.037177977  1.492931913
##   [51]  1.514501839  0.436290824  0.030230617  1.435590284  1.128478741
##   [56]  1.308048741  1.769292973  2.171418441  0.766308921  1.104697051
##   [61]  0.671556948  0.993133899  1.202345540  1.444388735  2.220546349
##   [66]  1.484355570 -0.088683826  0.511108922  1.459715801 -0.156190034
##   [71]  0.280677145  0.902248301  1.339117716  0.919675395 -0.485717744
##   [76] -0.050590019  1.594179528 -0.712258370 -0.201400533  1.236501838
##   [81]  1.753150294  0.863387931  0.199870184  0.281177078  0.525131115
##   [86] -0.649976118  0.190539720  1.749496955  1.210747501  0.601919140
##   [91] -0.021175571  0.583963032  1.662686269  1.518076522  0.147411651
##   [96]  1.627084151  0.534720911  1.980915484 -0.005839520  0.693506278
##  [101]  0.504843793 -0.101060881  1.399324865  1.235299391  0.074230979
##  [106]  1.627940040  0.338305514  0.742392991  0.691586529  0.567892211
##  [111]  0.096056462  1.206175554  1.162030031  0.135200523  1.211579017
##  [116]  1.364941989  0.920978690  0.945581919  1.787747623  1.452409272
##  [121]  0.458771629  0.766768638 -0.029831136  1.913010219  0.977890519
##  [126]  1.286090722  1.909656624  0.520003984  1.409081876 -0.251744046
##  [131]  2.272668873  0.836639904  1.943724720  1.046999637  0.824651938
##  [136] -0.770557888  0.655737268 -0.082442498  0.869556270  0.501547576
##  [141]  1.699845934  1.393387998  0.364723150  0.997394703  0.991171398
##  [146] -0.120219305  0.957854504  1.237908540  1.157058058  1.743598742
##  [151]  0.564564547  1.572419817  0.745560032  0.599300264  1.225750456
##  [156]  0.172856916  0.304050671  1.090167324  0.359291311  0.017326504
##  [161] -0.115662017  0.033991626  0.695578320  1.933785169  1.261349207
##  [166]  0.841261424  1.181479688  1.878267157  0.138780595  0.822400232
##  [171] -0.330317678  1.962842689  0.561694953  1.069387910  1.487331543
##  [176]  0.414633989  1.282906198  1.948995308  0.729270371  2.059812223
##  [181]  1.516282010  0.032493489  0.262865776  0.753554530  1.347565851
##  [186]  1.261580727  0.546134392 -0.192216324  0.569430484  1.054566294
##  [191]  1.851313018  1.118270461  1.074094307  0.681855035  1.275121176
##  [196]  1.415274271  1.383521343  1.630779977  1.505603177  1.210839862
##  [201]  0.083148998  1.274177426  0.595350660  1.712956645  1.217452874
##  [206]  0.292112185  1.616443681  1.055231950  0.610432958  0.536620926
##  [211]  1.045449763  1.437044942  1.952091226  1.413970153  1.111210771
##  [216]  0.022398878  0.098765691  1.720961167  0.465652606  1.493052984
##  [221]  1.781229311  0.208867108  0.347666270  1.659879941  0.997612620
##  [226]  0.354029200 -0.102729371  0.956094813  1.013312693 -0.394987686
##  [231]  0.278377565  1.980334429  0.785015056  0.485550305  0.598925358
##  [236]  1.284437549  1.600736930  0.857993259  1.189627512  0.627162623
##  [241]  1.176973485  1.455116031  0.689390502  1.455469668 -0.096861248
##  [246]  1.154887948  0.182589725  0.443430486  1.837511302  1.194147909
##  [251]  1.269899861  1.080115671 -0.240203877  0.024324274  0.569377037
##  [256] -0.352660618  0.751805337  1.764854092  0.718752575  0.502161086
##  [261]  0.910624250  1.459702238  1.954283023  1.342843222  0.914482229
##  [266] -0.582003458  1.675835958 -0.559461128  0.616245332  0.330746499
##  [271]  0.810586753 -0.241557146  0.407545380  2.085290259  1.056239533
##  [276] -0.187109962  1.358669119  0.967689898  0.267692703  1.554689501
##  [281]  0.093224536  0.134702328  1.027166618  1.089357038 -0.314525951
##  [286] -0.004758864  1.464943103  0.650280867  0.225500436  0.706453907
##  [291]  0.587542624  1.197356383  1.864598586  0.669783269  1.595855481
##  [296] -0.172622274  0.139443453  0.056659360  1.838017194  1.732193251
##  [301]  0.977176495 -0.612144887  0.777382028  1.137677739  0.572110001
##  [306]  0.496105719  0.958150164  0.782844226  0.392818493  1.331068011
##  [311]  0.735954003  1.505570804  0.253172934 -0.272022268  0.104210290
##  [316]  1.592857532  1.089824775  0.553266535  0.690445347  1.685750297
##  [321]  0.997372351  0.398905456  0.441202508  0.284116603  1.251505372
##  [326]  2.173872012  1.945384738  0.552419686 -0.142464836  1.137677739
##  [331]  1.549504803  0.383921327 -0.223998571  0.378441087  0.117796770
##  [336]  0.959369476  0.625443990  0.922629148  0.533478049  1.328775296
##  [341] -0.301584762  1.915191350  1.414354569 -0.240278938  1.478853963
##  [346]  1.206602381  0.140739268 -0.125392697  1.420834823  0.425790518
##  [351]  0.128391492  1.277228333  1.236775505  1.962103513  1.444081397
##  [356]  0.649981357  1.551017119  1.099012025  0.889352390  1.813746441
##  [361]  0.937711520 -0.599077946  0.725247857  0.237357490 -0.012476500
##  [366]  2.046892920  0.473849758  1.102674001  0.733231425  1.701563016
##  [371]  1.096996084  1.589987348  1.109192364  1.580692749  2.084149202
##  [376]  0.807972911  0.127282085  0.175377643  0.498782172  1.094762105
##  [381]  1.529800929  1.277900032  1.141788637  0.675586986  0.800171784
##  [386]  1.699443732 -0.084429904  1.405645842  0.066319261  1.129767824
##  [391]  1.413221260 -0.264002656 -0.676658506  1.496828469  1.078886660
##  [396]  0.316595537  1.326480786  0.211301236  1.493052984  1.384293794
##  [401]  0.660263040  0.890523087  1.614419730 -0.089685723  0.023013372
##  [406]  1.844195264 -0.135848934 -0.261654751  0.819013993  1.254039263
##  [411]  2.170363030  1.190646323  2.037142054 -0.028147013  1.122645550
##  [416]  0.694994786  0.954157058 -0.496445985  0.274082915  1.175997056
##  [421]  0.989712940  0.902489632  0.708127342  1.146260821  0.657115503
##  [426]  1.098568501  2.173190981  0.274814043  0.352878135  0.577212261
##  [431]  1.155283349  1.554160199  0.287794983  1.158457593  0.870947065
##  [436]  0.802884739  0.940313982  1.311484497  0.568011459  0.885327537
##  [441]  1.154673851  0.108184276  1.912933128  2.083724857  0.850711839
##  [446]  1.015925552  1.269298621  0.816091884  0.103916859 -0.476467217
##  [451]  0.678166852  0.038196679  1.688872062  1.282184508  0.222964547
##  [456]  0.360278114  1.860505621  0.365148714  0.775474670  0.526792867
##  [461]  0.758137586 -0.209648592  1.154548781  1.691331566  0.817686782
##  [466]  0.572757951  0.238838765  1.333331705  1.042835035  1.105865625
##  [471] -0.803655255  1.015802909  1.039281001 -0.329560855 -0.072122281
##  [476]  1.586839827  1.076373171  1.068414376  0.300231874  2.142023370
##  [481]  0.355684546  1.586279557  0.867925213  0.964440682  1.116131714
##  [486]  0.512113111  0.036125567  1.133549115  0.622202414  0.685546469
##  [491]  1.574236561  1.018955028  0.203593182 -0.509217869  1.017330774
##  [496]  1.771657667  0.859178357 -0.005893735  0.498782172  0.048903290
##  [501]  1.461476412  0.740725508  1.778306751 -0.373673574  0.917166521
##  [506]  1.967254190  0.630202204  1.075556032 -0.112410005  0.532298165
##  [511]  1.119217480 -0.360367286  0.435200598  0.110614693  0.962570380
##  [516]  1.578891790  1.987337827  1.058894245 -0.272781589  0.879982153
##  [521] -0.466836836 -0.142413818  1.639510928  0.869382372  0.104265823
##  [526]  1.882138507  1.491481465  0.980669037  0.509643253 -0.586274197
##  [531]  1.571719321  0.723953933  0.566030833  1.526399955  0.990687394
##  [536] -0.214852292 -0.283879582 -0.147764894  1.531773609 -0.074412426
##  [541]  1.264613671  0.173436800  0.407638157  1.759349444  1.519423325
##  [546]  0.922962712  1.050362457  1.738899098  1.267388634  0.941167431
##  [551]  1.259104590 -0.673085918  1.736339470  1.608975099  0.446215674
##  [556]  1.310487066  1.409257786  1.838527931  0.578821056 -0.203142115
##  [561]  1.449965976  1.191169788  1.615435121  1.543325915  1.646312717
##  [566]  1.354570077  0.941167431  1.224262300  0.900645099  0.284019119
##  [571]  0.909711220  0.626515887  0.256667545  1.192845476  1.307054083
##  [576]  1.408908247 -0.115678094  0.316868413  1.313324615  0.146974607
##  [581]  0.538217267  0.677929523  0.031468829  0.542740851  0.562586823
##  [586]  0.310483375  0.388593668  0.741031089  0.508229293  1.085578488
##  [591]  0.156589128  1.479379799  1.584508092  1.266860662  0.723479901
##  [596] -0.787738380  0.601495429  1.878858912  1.296619897 -0.607245116
##  [601]  0.097789547  1.574808201 -0.726452780  1.105344338  1.072140016
##  [606]  0.589390087  1.001748653  1.030871726  1.165741364 -0.539436979
##  [611]  1.053830135  0.751945474  1.076140567  1.051943611  1.631881690
##  [616]  0.349316776  0.268687832 -0.420287801 -0.005622150  0.223649350
##  [621]  1.302713644  1.612178175  0.654015493  0.796425896  0.779650077
##  [626]  0.229389647  1.121510577  1.002078546  1.284437549  1.051504688
##  [631]  1.079957027  1.160207741  0.444412689  0.336372683  1.572343296
##  [636]  0.591672834  1.114972200  1.759477969  1.141277466 -0.232415877
##  [641]  0.436388071  1.248538795 -0.102894185  0.742160230  1.730067173
##  [646]  1.285456223 -0.082803934  0.711126077 -0.172965913 -0.416612627
##  [651] -0.387279051  0.850886270 -0.272273519  1.793643874 -0.339446092
##  [656]  0.153454394  0.534798379  0.216691375  1.001406726  1.559486617
##  [661] -0.118953873 -0.177562572  0.918223165  0.131524444  1.428905097
##  [666] -0.003542928  0.315711274  1.205207254  0.917275455  0.888606257
##  [671]  0.843748477  1.477382810  1.187957600  1.223356761  0.488606475
##  [676]  0.162863666 -0.106200578  1.398455669  1.462979573  0.957896518
##  [681]  0.083524516  1.128478741  1.377472186  1.879077600  1.383663192
##  [686] -0.301030304  1.595855481 -0.162525287  0.336493476  1.527205330
##  [691]  0.856335663  1.079957027  1.063751807  1.061670380  0.354858026
##  [696]  0.482347522  0.454112036  0.298652701 -0.602846888 -0.111561334
##  [701]  0.855487239  0.801140957 -0.337853679  1.603410042  0.074293385
##  [706]  1.113116928  1.198464365  1.082598924  0.106187492  0.083027141
##  [711] -0.024357748  0.698022761  1.255863102  0.045049268  0.015942576
##  [716]  0.481602018 -0.187725859  0.335824091  1.088686536  1.118901392
##  [721]  0.799385689 -0.500161418  0.747081725  0.670655362 -0.425455687
##  [726]  1.083977567  1.013603501  1.038526245  1.380734089  1.193877902
##  [731]  0.340337768  0.719368914  0.911640134 -0.325383809  0.695374647
##  [736]  1.290327822  0.737040161 -0.045011785  0.634516098  0.617671371
##  [741]  1.014630921  2.172435670  0.382167055 -0.011483765  1.490684840
##  [746]  0.343166607  1.020243040  1.072142531  0.101216764 -0.205913204
##  [751]  1.593603586  1.630855228  1.146446786  1.977845698  0.723305264
##  [756]  0.675806101 -0.965501307  0.806958959  0.415098660  1.464944512
##  [761]  0.202586257  0.655612522  0.061512851  1.660984348  1.051451306
##  [766]  0.907444278 -0.087707193  0.168034522  1.537366322  1.723945490
##  [771]  1.085015056 -1.703053414  2.028362655  1.732193251  0.402679244
##  [776]  0.842927884  0.566346973  1.158501654 -0.432258176  0.912322239
##  [781]  0.256312016  0.526741843 -0.432316051  0.446567525  1.055969373
##  [786]  1.418223735  1.406283718  1.525781781  1.373225126  1.292310401
##  [791]  0.299800648  1.484834102  0.791705842  0.880089207  1.562870980
##  [796] -0.458003819  1.787543864  1.390695936  1.218439053  0.394415659
##  [801]  0.114376813  1.309048253  0.300442355  0.476759549  1.663957674
##  [806]  0.490297999  1.387176042  0.959993361  1.593941481  0.107877429
##  [811]  0.624534971  1.678000942 -0.139002594  1.820752713  1.329775400
##  [816]  0.802174478  1.241448697  0.710481504  0.136903416  0.160985540
##  [821]  1.586419512  0.362144387  0.978990241  0.961689735  1.401811580
##  [826] -0.139099667 -0.028558198  1.606514512  1.588506539 -0.451650636
##  [831]  0.266210524 -0.721950327  1.257010348 -0.024731720 -0.369418291
##  [836]  0.354150637  0.964440682  0.951973477  0.208042497 -0.373759781
##  [841]  0.285494874  0.695720284  0.898699152  0.937716409  0.467032045
##  [846]  0.889976104  1.419476963  0.423593379 -0.014121780  0.183541542
##  [851]  0.550963632  1.549222260  0.138110505  0.709673457 -0.184465798
##  [856]  1.102556315  1.202016959  1.519346835  0.114609416  1.730104338
##  [861]  1.022926010  0.962581576  1.101871427  1.511267860  1.082283614
##  [866]  1.023296544  0.774137876  0.991710940  1.198993484  1.191534617
##  [871]  1.586888866 -0.543564945  1.201463233 -0.129978765  1.409178804
##  [876] -0.053925600  2.019119514  0.772806078  0.978838687  0.890910526
##  [881]  0.493134726  0.459459354  1.122757753  1.207382902  0.825724891
##  [886] -0.384853833  0.407646444  0.369648038  1.531773609  0.902064573
##  [891]  1.547512306  1.360314485  0.812514442  0.882567772  1.557442098
##  [896]  0.348762805  1.081303633  0.659400278  0.038026859  1.619644444
##  [901]  1.115813500  1.473886431  0.949531966  1.322356252  1.003072222
##  [906]  1.269936326  1.050497930  0.968852146  1.122505635  0.067176884
##  [911]  0.882335029 -0.550333407  1.642358987  0.416212840  0.920530378
##  [916]  0.563238191  0.766308921  0.919629720  1.014178402  1.372933651
##  [921]  1.093733826  0.957180386  0.751587485 -0.141063333  1.647708425
##  [926]  2.033035264  0.965470927  0.545859246  0.886786260 -0.241369893
##  [931]  0.431076952 -0.082670268  0.206229947 -0.157847062  0.593659350
##  [936]  1.214995233  0.967492203  0.457851759  1.139459113  1.435079712
##  [941] -0.178261598  1.240874928  0.692336745  1.736686091  0.897246753
##  [946] -0.657056375  0.473154001  0.177912952  1.590167632  0.806454012
##  [951]  1.621422744  1.577620270  0.350305046  1.214420806  0.471820643
##  [956]  0.764193222  1.823586013  1.675934595  0.697867510  0.611906872
##  [961]  1.343017381 -0.069050889  0.844435583  1.116286551  2.330354740
##  [966]  0.179688906  2.000448234  0.349634667  0.688322396  0.821738009
##  [971]  1.138946199  0.759303185 -0.349913239  0.544246336  0.356945776
##  [976]  1.446846868  1.358369344  0.796156773  0.731980617 -0.921060632
##  [981]  0.849970531  1.867215500 -0.893669727  1.360192069  1.164699649
##  [986]  0.562192068  0.446076847  1.320805859  0.564428761  0.265174671
##  [991]  0.486069510  1.142176568  0.586007489  0.858473177  1.057551509
##  [996]  1.358417385  0.740135589  0.962098023  0.829947238  1.656737109
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
##   0.64416796   0.37639301 
##  (0.11902592) (0.08416171)
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
## [1] -0.21094250 -1.11083543 -0.96858920 -0.08273692 -0.90933057 -0.01743893
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
## [1] 0.0473
```

```r
sd(results2$thetastar)  #the standard deviation of the theta stars is the SE of the statistic (in this case, the mean)
```

```
## [1] 0.8663505
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
## t1*      4.5 -0.04134134   0.9098011
```

One of the main advantages to the 'boot' package over the 'bootstrap' package is the nicer formatting of the output.

Going back to our original code, lets see how we could reproduce all of these numbers:


```r
table(sample(x,size=length(x),replace=T))
```

```
## 
## 0 2 3 4 7 9 
## 3 1 1 2 1 2
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
## [1] 0.0051
```

```r
se.boot
```

```
## [1] 0.9432358
```

Why do our numbers not agree exactly with those of the boot package? This is because our estimates of bias and standard error are just estimates, and they carry with them their own uncertainties. That is one of the reasons we might bother doing jackknife-after-bootstrap.

The 'boot' package has a LOT of functionality. If we have time, we will come back to some of these more complex functions later in the semester as we cover topics like regression and glm.

