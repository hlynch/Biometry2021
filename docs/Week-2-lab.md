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
## 1 2 5 6 7 8 9 
## 1 1 2 1 1 3 1
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
## [1] 0.0227
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
## [1] 2.753036
```

```r
UL.boot
```

```
## [1] 6.292364
```

Method #2: Simply take the quantiles of the bootstrap statistics


```r
quantile(xmeans,c(0.025,0.975))
```

```
##   2.5%  97.5% 
## 2.7975 6.4000
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
##    [1] 5.2 5.5 5.3 5.0 4.9 5.1 4.2 5.9 5.3 4.2 4.3 5.1 5.2 5.5 4.8 5.1 6.5 3.8
##   [19] 4.3 5.1 4.8 6.0 5.4 5.9 2.7 3.4 5.1 4.0 4.1 4.3 2.9 4.1 3.3 4.2 5.6 3.8
##   [37] 4.9 5.2 3.4 4.0 3.7 3.7 5.0 4.7 3.9 4.7 4.6 4.9 5.8 3.6 5.9 4.7 4.2 4.5
##   [55] 4.1 4.2 6.2 4.7 4.0 4.1 3.4 3.2 3.6 4.7 3.2 5.5 2.9 5.5 5.5 3.9 5.3 4.8
##   [73] 4.6 3.8 4.3 5.1 5.2 5.7 5.9 3.9 4.8 3.8 3.2 5.1 6.6 3.9 6.3 3.4 5.6 5.3
##   [91] 2.9 3.7 4.9 5.5 4.4 4.5 4.3 4.5 6.3 4.0 6.0 3.9 5.3 3.9 2.5 6.6 3.8 4.7
##  [109] 4.7 3.0 3.4 4.2 3.9 4.9 5.4 5.7 6.1 3.7 5.1 4.5 2.4 4.4 4.2 4.9 3.3 5.3
##  [127] 4.8 5.3 4.9 5.6 3.3 3.8 2.8 4.1 4.6 3.2 6.0 4.0 3.5 4.4 5.8 3.1 4.3 4.3
##  [145] 4.1 3.4 3.9 4.3 5.4 3.7 5.7 2.5 5.4 4.2 3.1 4.9 4.3 5.8 4.4 4.2 4.3 5.2
##  [163] 4.5 5.4 4.5 4.7 3.2 4.6 5.9 4.0 3.1 5.6 4.2 3.0 5.0 4.1 4.2 3.1 5.1 5.2
##  [181] 3.9 5.1 5.1 3.5 6.2 5.3 6.2 4.3 4.7 5.1 4.0 3.0 4.0 4.3 5.0 3.7 3.8 4.3
##  [199] 4.1 2.6 4.8 3.4 2.9 4.4 4.4 4.8 3.6 4.0 5.6 4.0 4.6 5.3 3.1 3.5 5.7 3.6
##  [217] 4.1 4.5 3.9 5.7 5.6 4.6 5.7 5.5 3.3 3.9 2.6 5.0 4.4 5.0 3.4 3.6 4.1 4.0
##  [235] 3.8 5.0 4.9 4.0 2.8 4.9 5.0 5.3 3.7 3.9 4.1 2.9 4.6 4.8 5.2 5.2 4.9 3.0
##  [253] 3.4 6.1 4.1 3.4 6.5 3.9 4.5 5.0 5.8 2.8 4.1 5.4 5.0 5.7 5.9 3.8 4.6 5.7
##  [271] 3.7 4.3 3.5 5.0 4.6 4.4 4.1 3.9 4.1 3.6 5.7 3.2 3.9 3.6 3.4 4.7 3.5 5.3
##  [289] 4.4 5.7 4.3 4.1 4.4 4.5 5.1 2.2 6.2 5.7 5.1 4.1 4.4 4.6 4.4 5.2 6.5 4.3
##  [307] 4.2 3.6 3.6 5.5 3.5 6.8 3.7 4.5 5.1 4.6 6.2 4.4 4.7 4.6 4.2 4.8 4.4 3.1
##  [325] 3.8 4.6 4.0 4.4 4.3 5.3 5.6 3.2 4.6 4.0 5.8 4.4 4.1 5.2 4.2 5.4 3.6 5.6
##  [343] 4.6 4.9 5.7 6.0 4.9 4.1 4.0 6.0 4.3 4.9 4.7 4.6 3.2 3.8 4.8 3.3 5.3 5.2
##  [361] 4.9 4.6 5.7 4.7 4.8 4.1 4.5 3.4 4.2 5.0 4.0 6.2 5.6 4.4 3.6 4.7 4.4 6.0
##  [379] 3.7 3.7 4.5 3.2 5.2 6.3 5.8 3.8 4.1 2.5 3.9 4.2 4.3 3.0 3.5 3.8 3.4 4.6
##  [397] 3.7 5.1 5.3 4.0 3.8 2.1 4.8 4.9 3.3 4.9 4.7 5.5 5.0 4.3 4.3 5.1 4.2 5.1
##  [415] 2.7 6.1 2.7 5.4 4.0 3.6 5.4 5.6 2.7 3.3 4.1 4.8 4.2 4.4 5.7 3.6 4.8 4.4
##  [433] 4.3 5.9 4.1 3.5 3.4 3.5 5.5 5.6 3.6 3.9 2.8 3.2 4.8 5.8 5.8 4.5 4.8 4.6
##  [451] 3.3 3.5 4.2 3.7 5.0 4.3 5.2 3.9 4.7 6.3 4.9 5.3 3.5 4.6 5.5 3.8 5.5 3.5
##  [469] 4.2 4.6 4.3 4.6 4.5 5.5 5.0 3.7 4.9 4.2 5.0 3.3 5.3 3.7 5.6 5.4 3.9 5.3
##  [487] 3.5 4.4 5.6 4.8 3.8 4.5 4.8 5.3 3.3 4.9 3.5 4.7 3.7 4.2 3.2 3.1 4.0 4.6
##  [505] 5.5 7.2 4.5 5.2 4.3 4.3 3.9 5.8 3.7 6.0 5.5 4.2 3.1 4.1 3.7 4.3 4.7 5.0
##  [523] 5.2 3.1 4.5 5.4 3.5 5.0 4.6 2.9 4.3 3.9 6.0 3.9 4.2 4.3 4.5 4.4 3.9 3.6
##  [541] 4.7 6.3 6.2 5.6 3.5 4.6 3.9 3.8 5.3 4.3 4.0 3.9 6.1 4.8 3.0 6.3 4.5 4.9
##  [559] 5.0 5.0 4.5 4.6 6.1 6.0 5.0 3.7 4.4 4.7 3.5 4.7 3.9 4.8 3.7 3.8 4.2 4.5
##  [577] 3.8 3.6 6.3 4.0 3.5 4.9 3.3 4.4 5.0 4.8 2.6 5.7 5.0 5.8 4.7 3.4 5.9 5.9
##  [595] 4.1 5.5 4.2 4.7 4.2 5.2 3.2 3.9 5.4 3.7 4.0 5.0 4.9 6.3 3.8 4.9 4.5 3.5
##  [613] 4.2 4.9 5.0 4.7 4.1 5.8 4.8 3.3 4.0 3.1 4.8 4.7 3.7 3.1 3.8 3.4 5.5 4.1
##  [631] 3.9 5.7 3.1 6.2 5.7 3.8 3.8 3.7 4.7 4.7 4.4 5.3 4.5 4.7 4.7 5.4 4.9 4.7
##  [649] 7.2 5.7 4.5 4.2 4.1 5.2 4.6 3.6 3.7 2.9 4.9 3.9 5.7 3.2 4.5 4.0 3.4 5.3
##  [667] 5.2 5.6 4.8 4.9 3.0 3.7 5.1 4.4 5.0 1.9 4.0 5.3 3.8 5.1 4.7 3.5 4.0 3.2
##  [685] 4.3 4.8 5.5 5.5 5.4 5.4 6.1 4.7 4.8 3.5 3.6 4.1 5.1 3.1 3.6 4.2 4.6 3.7
##  [703] 3.8 3.5 4.6 5.0 3.6 4.9 4.1 4.7 3.5 4.3 5.1 6.0 4.0 5.2 2.5 3.9 4.2 4.7
##  [721] 6.2 4.1 3.7 4.3 4.6 4.7 4.7 6.1 4.8 3.9 4.2 4.1 4.7 6.1 5.4 5.1 5.8 2.9
##  [739] 3.6 5.4 3.9 5.4 2.4 3.1 4.8 6.0 4.4 3.6 5.9 4.2 6.6 5.1 4.9 5.2 4.5 5.0
##  [757] 4.7 5.4 5.5 3.4 3.2 4.0 5.0 4.6 4.7 4.7 2.1 3.4 4.8 6.0 4.0 4.3 5.8 4.6
##  [775] 5.1 5.8 5.4 5.8 4.2 6.2 2.7 3.6 5.6 5.6 6.4 4.2 4.5 3.6 5.3 3.5 4.2 3.8
##  [793] 5.7 3.9 4.5 3.3 5.5 3.8 3.1 4.0 5.2 4.5 5.0 5.0 4.6 6.6 3.9 4.5 4.6 5.0
##  [811] 4.1 4.3 4.0 5.8 6.3 3.7 3.0 4.1 4.8 3.3 5.4 4.8 5.9 4.9 3.6 4.9 4.5 2.5
##  [829] 4.8 4.2 4.8 4.0 3.9 4.8 6.1 4.9 4.4 4.5 3.3 5.1 4.1 3.8 2.3 5.4 3.7 3.7
##  [847] 4.2 4.1 5.4 6.0 5.4 4.5 5.4 6.0 5.0 4.2 4.6 4.1 3.6 4.8 3.8 4.1 4.0 4.6
##  [865] 4.3 3.6 3.9 3.9 3.4 4.9 5.3 3.6 3.2 4.1 4.2 3.9 5.3 3.8 6.6 4.4 4.9 5.8
##  [883] 5.8 4.7 4.6 5.0 5.8 2.9 4.3 5.0 3.6 3.0 4.7 3.6 4.7 5.0 5.2 5.4 5.6 5.9
##  [901] 4.7 4.4 6.3 2.7 4.0 5.5 4.7 4.9 5.1 4.4 3.7 5.2 4.7 4.8 5.9 4.3 2.6 5.3
##  [919] 4.5 3.4 5.1 4.1 4.2 4.4 6.5 3.2 5.6 3.7 3.3 4.5 1.7 3.9 4.6 4.9 5.1 5.5
##  [937] 3.4 4.0 4.9 4.4 3.6 5.1 4.6 5.4 2.9 4.3 5.1 3.8 7.0 4.5 3.1 4.9 4.4 5.0
##  [955] 6.3 4.1 4.1 4.6 3.6 3.2 5.2 5.1 4.4 4.2 2.8 4.1 5.1 5.2 4.1 5.1 3.1 4.4
##  [973] 3.9 4.5 5.6 2.0 4.8 3.3 4.0 3.1 4.9 3.5 3.1 3.4 3.9 4.7 6.0 5.4 4.7 5.3
##  [991] 4.4 4.3 3.5 4.6 4.9 5.7 3.5 4.1 6.0 4.2
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
##    [1] 3.4 4.0 5.3 6.1 2.9 5.3 4.9 6.2 5.6 4.3 4.7 5.3 5.5 4.8 4.3 2.8 4.2 4.8
##   [19] 5.0 4.0 4.1 4.6 5.1 5.0 5.3 5.4 6.2 5.3 4.7 4.4 5.9 7.0 4.3 3.7 4.7 3.3
##   [37] 6.0 3.9 5.3 5.3 5.0 4.7 2.9 4.4 3.9 5.7 3.7 5.3 4.1 5.3 5.8 3.2 5.0 4.4
##   [55] 3.5 5.4 3.3 4.7 3.1 4.9 4.8 3.0 4.1 3.9 4.0 4.9 5.9 4.8 4.8 4.2 5.6 3.5
##   [73] 3.5 4.5 4.3 4.0 3.5 4.9 3.9 4.3 5.1 5.8 4.4 5.4 3.7 4.4 2.7 6.2 6.4 5.2
##   [91] 3.7 3.8 2.5 4.1 4.2 3.7 4.6 4.5 4.8 3.9 2.7 5.8 6.4 3.7 3.6 4.3 5.9 5.2
##  [109] 4.4 4.9 3.5 4.1 3.3 4.9 4.3 4.6 6.1 3.8 5.7 3.5 4.0 2.6 3.0 5.1 3.0 3.2
##  [127] 4.3 4.4 2.7 5.8 4.0 2.7 5.0 4.3 4.1 3.9 4.7 3.2 4.9 3.6 5.6 3.4 4.2 3.4
##  [145] 3.4 4.7 4.0 3.8 3.8 3.2 2.8 6.1 4.2 4.8 4.9 5.9 5.1 5.3 3.6 5.5 3.4 5.0
##  [163] 2.1 4.8 3.6 3.4 5.8 4.8 4.0 3.3 5.1 4.4 2.5 4.4 5.5 2.6 4.5 4.4 4.6 2.1
##  [181] 3.6 3.7 3.9 5.3 5.4 6.5 4.2 2.7 4.1 5.0 2.6 5.1 5.7 3.7 4.0 5.3 4.0 5.5
##  [199] 5.6 3.9 3.7 4.8 5.8 4.9 5.5 3.2 3.1 4.1 4.1 4.1 3.0 3.7 4.3 4.0 4.2 4.8
##  [217] 4.8 3.7 5.0 4.8 3.9 5.6 5.0 4.5 5.3 4.7 3.1 4.4 2.3 6.0 6.3 4.4 4.0 4.5
##  [235] 4.6 4.5 5.2 4.0 3.0 4.7 6.5 4.8 5.3 5.0 4.5 4.8 5.6 4.4 3.4 4.7 4.4 5.9
##  [253] 5.1 4.1 3.9 5.2 4.1 3.0 4.9 4.3 4.7 3.9 4.5 4.8 4.0 3.7 4.3 3.7 6.1 2.3
##  [271] 4.7 4.3 4.2 3.8 5.9 5.2 5.4 6.5 3.2 4.9 6.1 3.2 5.6 4.3 4.0 3.8 5.4 3.5
##  [289] 4.6 4.7 6.0 3.5 4.4 5.9 4.4 5.4 4.7 4.2 6.0 5.0 4.0 4.3 5.8 4.0 3.2 4.0
##  [307] 4.8 5.0 2.7 6.4 5.4 5.6 4.3 4.1 4.7 3.8 3.9 4.0 5.7 5.3 5.4 4.0 5.2 5.2
##  [325] 4.0 3.5 5.0 2.9 6.1 4.5 4.8 3.2 4.7 4.4 2.4 4.8 4.5 3.0 4.4 5.0 5.8 4.3
##  [343] 3.3 4.8 4.6 3.3 4.7 4.3 4.5 4.6 4.2 3.8 2.7 3.9 5.1 4.1 4.3 5.0 3.9 3.4
##  [361] 6.4 4.5 3.8 3.7 5.1 4.3 4.7 4.6 4.3 2.9 5.1 5.3 4.5 4.6 4.6 5.6 5.5 4.7
##  [379] 5.3 4.5 5.4 4.5 5.9 6.2 3.6 2.9 3.1 4.9 5.4 5.7 4.8 6.3 5.2 4.6 4.4 4.2
##  [397] 5.1 4.0 3.7 3.3 4.7 3.3 5.4 4.2 4.4 3.7 4.4 3.9 4.4 5.6 4.8 5.1 6.4 4.0
##  [415] 4.3 2.5 5.5 4.1 5.5 4.0 5.4 5.9 4.5 4.5 5.9 6.6 6.0 3.4 4.7 3.5 5.5 5.6
##  [433] 2.3 5.7 4.0 5.4 4.0 5.9 4.6 4.3 4.0 5.7 3.2 6.4 4.8 5.1 4.8 4.2 3.7 3.2
##  [451] 6.0 4.4 6.0 5.7 3.7 3.4 5.6 4.4 4.0 5.3 3.3 5.6 3.6 4.6 4.7 5.1 3.8 4.0
##  [469] 4.8 5.8 6.7 5.5 5.1 2.9 4.3 3.1 4.0 4.8 4.1 3.7 4.8 2.4 3.2 5.4 5.2 6.4
##  [487] 3.6 4.4 5.4 3.9 4.4 5.6 5.4 5.5 5.0 6.0 4.7 5.0 3.9 4.4 4.2 3.8 2.9 5.0
##  [505] 4.7 5.8 5.6 4.1 4.0 4.7 3.4 3.5 4.6 5.4 4.7 4.8 5.1 5.8 3.9 5.6 4.2 4.0
##  [523] 4.9 5.0 5.9 4.3 5.8 3.1 3.2 5.3 3.9 3.2 4.5 4.9 4.2 3.8 6.1 3.1 5.5 3.9
##  [541] 3.5 4.5 7.0 3.6 4.1 4.5 4.7 4.8 3.3 6.5 4.2 2.7 5.2 5.1 3.6 5.0 4.6 5.0
##  [559] 6.1 6.0 6.9 2.5 4.6 4.4 4.4 4.1 3.5 3.1 6.1 4.0 5.1 4.5 5.3 4.4 4.4 3.2
##  [577] 3.0 4.8 5.3 4.2 4.4 3.2 5.9 4.9 5.3 4.3 6.0 5.1 4.3 5.0 4.8 3.5 4.4 4.0
##  [595] 5.4 4.6 4.7 3.7 4.9 4.9 5.3 5.4 4.1 3.8 3.8 3.3 4.8 4.7 3.9 3.2 3.6 3.8
##  [613] 3.9 3.2 2.6 5.5 5.2 4.7 3.9 5.2 4.0 5.0 3.3 5.7 3.7 5.3 2.8 4.7 4.5 3.8
##  [631] 5.5 5.9 4.1 5.8 4.1 4.7 3.3 3.1 3.9 6.0 4.1 3.8 2.9 3.8 5.2 4.6 5.2 5.1
##  [649] 4.5 4.5 4.6 3.4 4.0 4.1 3.5 4.2 4.8 3.8 5.4 5.7 6.6 5.3 5.7 5.0 5.4 3.1
##  [667] 4.3 4.1 2.9 3.4 4.0 4.0 4.5 4.1 4.2 2.9 3.5 3.3 5.1 4.1 6.2 4.4 3.5 3.4
##  [685] 4.9 3.7 3.0 5.0 4.0 4.4 3.7 5.5 5.5 1.9 6.3 5.3 4.6 3.5 5.4 4.3 5.0 4.0
##  [703] 5.5 3.6 2.2 4.4 5.8 3.6 4.1 3.8 3.7 3.3 4.3 5.3 5.4 3.0 5.2 2.4 3.8 4.6
##  [721] 4.1 5.0 3.9 5.8 4.4 4.3 5.1 4.2 3.8 5.2 5.5 5.8 3.5 4.4 4.4 3.0 5.0 5.2
##  [739] 4.4 4.0 5.3 6.3 5.4 5.3 4.5 5.6 4.3 5.4 3.8 4.1 5.2 5.1 4.6 4.7 5.4 2.8
##  [757] 4.1 4.2 3.2 3.9 4.4 3.4 4.1 5.7 5.1 6.1 3.2 3.1 4.9 3.2 4.1 5.1 5.2 4.0
##  [775] 3.4 5.0 5.2 4.1 3.3 3.7 4.3 2.5 3.9 5.0 4.3 3.9 4.9 4.0 3.5 3.9 3.9 4.2
##  [793] 1.7 6.2 5.0 4.0 4.8 2.7 4.6 4.6 3.6 4.1 4.9 3.8 3.3 5.0 4.8 3.3 4.8 3.7
##  [811] 3.6 3.5 6.0 4.8 3.7 3.9 4.4 5.1 2.4 4.6 3.7 4.0 5.5 5.5 6.8 3.9 3.9 2.7
##  [829] 5.8 5.2 4.1 3.7 5.8 4.5 4.7 5.3 4.4 4.0 3.8 4.1 4.9 5.2 5.2 3.5 3.9 4.1
##  [847] 3.6 3.6 3.8 3.9 4.4 4.7 4.5 3.6 3.9 4.3 3.5 4.7 5.3 3.2 3.5 3.7 3.9 3.8
##  [865] 5.0 3.5 4.5 5.3 4.1 5.3 5.8 3.6 4.5 4.6 3.7 5.5 5.1 4.1 5.1 3.8 4.0 3.3
##  [883] 4.6 5.2 5.6 4.6 4.4 5.1 6.0 5.5 4.0 4.9 4.4 4.6 3.0 5.8 3.3 4.5 4.3 5.9
##  [901] 5.8 4.1 3.5 3.4 4.5 2.5 4.1 4.8 5.7 4.2 5.7 5.4 4.3 5.3 3.6 5.9 4.5 5.2
##  [919] 2.7 5.2 4.2 2.5 4.3 4.1 5.9 4.8 5.3 3.3 4.1 3.8 3.4 3.4 6.4 4.5 5.0 4.3
##  [937] 3.9 4.9 3.6 4.0 5.1 4.1 6.5 5.1 4.5 4.5 4.9 5.4 4.2 6.3 4.9 4.6 4.5 4.9
##  [955] 6.0 4.8 3.9 5.1 6.1 3.8 4.1 5.4 3.6 4.6 5.8 6.6 4.1 5.1 3.5 5.6 4.0 5.2
##  [973] 4.7 3.6 3.8 4.1 4.8 5.9 4.1 4.2 5.6 4.2 4.3 4.3 4.6 5.4 4.7 4.0 3.9 4.2
##  [991] 3.7 4.5 4.4 4.7 4.3 4.0 4.5 3.6 4.9 4.6
## 
## $func.thetastar
## [1] -0.025
## 
## $jack.boot.val
##  [1]  0.45527066  0.35891239  0.28103448  0.20245232  0.09109792 -0.14447674
##  [7] -0.20523256 -0.33863636 -0.43352436 -0.48525469
## 
## $jack.boot.se
## [1] 0.9770386
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
##    [1] 3.5 5.5 3.5 3.5 4.7 3.0 3.4 5.9 5.4 4.9 3.9 4.9 4.0 4.2 4.3 5.7 5.3 2.9
##   [19] 4.6 4.1 4.9 4.7 4.8 3.9 4.6 4.4 3.0 4.7 3.6 4.5 3.7 4.8 4.8 3.9 3.6 5.2
##   [37] 4.1 5.3 4.3 3.3 4.0 5.4 5.4 4.0 5.0 2.7 5.3 4.0 5.3 3.3 5.7 4.7 4.7 4.9
##   [55] 5.3 4.7 4.2 4.6 4.8 2.9 4.9 4.3 5.3 4.5 4.3 5.6 4.6 3.2 3.8 2.7 4.3 3.6
##   [73] 5.1 4.5 4.1 3.9 5.2 4.1 4.7 5.4 5.2 5.0 4.9 5.7 6.3 6.2 4.5 4.4 4.8 3.8
##   [91] 4.7 4.4 3.3 3.8 4.5 4.6 4.8 4.0 5.2 4.9 4.4 4.3 5.6 3.9 4.2 6.0 4.6 4.8
##  [109] 3.8 4.5 4.6 5.0 4.1 5.7 3.1 6.3 3.7 3.8 3.8 4.9 4.0 5.1 4.9 3.3 4.5 3.1
##  [127] 3.6 4.6 4.5 5.6 3.9 3.7 4.2 4.1 4.9 4.2 5.0 4.2 4.4 3.2 4.2 2.1 3.2 4.2
##  [145] 4.3 3.2 3.6 4.9 5.9 5.9 3.9 4.5 3.9 3.6 4.4 5.0 3.7 4.3 3.0 3.7 3.4 4.3
##  [163] 4.7 4.9 5.0 3.7 3.4 4.0 4.0 4.6 4.8 4.2 4.9 4.9 5.2 2.4 4.4 3.9 4.6 4.4
##  [181] 5.9 4.4 4.1 5.1 4.8 3.6 5.8 4.8 5.0 4.4 4.7 4.9 6.3 4.7 4.5 3.2 4.5 4.2
##  [199] 4.9 5.3 4.2 5.3 4.2 5.1 4.7 3.7 5.4 5.3 5.7 4.4 4.4 4.2 5.4 4.5 4.0 4.4
##  [217] 4.2 4.8 3.7 4.1 3.4 5.7 5.0 5.7 4.8 2.7 4.0 4.0 4.7 5.4 5.7 3.7 4.4 3.5
##  [235] 3.6 3.0 5.9 5.7 3.5 5.8 4.5 4.2 3.7 3.3 4.8 4.7 5.6 4.7 5.5 4.9 4.9 3.5
##  [253] 4.7 4.7 4.2 6.2 5.0 5.3 4.5 5.5 7.2 2.9 4.1 4.9 7.1 4.0 5.7 3.1 4.8 4.8
##  [271] 7.2 4.8 4.4 4.0 4.8 4.5 4.0 5.2 5.4 5.5 5.0 2.6 5.4 3.7 3.1 5.6 4.4 4.3
##  [289] 5.8 6.0 5.5 6.8 3.8 3.5 6.4 4.0 4.7 6.6 4.1 3.9 3.0 3.0 4.3 4.6 5.2 3.5
##  [307] 5.7 3.7 4.4 4.1 5.3 4.4 2.5 4.5 5.2 3.9 4.6 4.2 3.1 5.0 5.4 4.2 5.6 4.1
##  [325] 4.4 4.6 4.0 3.8 2.9 5.3 6.3 4.4 4.7 5.2 4.9 6.3 4.1 6.2 4.5 4.5 5.3 3.4
##  [343] 4.0 3.4 4.5 4.4 5.7 3.5 3.8 3.9 5.3 3.8 4.1 3.9 4.7 4.9 4.1 3.0 4.8 3.6
##  [361] 4.5 3.4 3.6 3.4 5.4 5.5 3.1 5.1 4.4 3.0 5.3 5.5 5.7 3.7 4.4 4.9 3.4 5.6
##  [379] 3.8 4.7 5.0 3.9 3.2 2.8 4.9 4.9 5.4 5.5 5.2 4.1 4.1 5.9 4.4 4.9 4.4 2.9
##  [397] 3.9 5.9 4.5 5.1 4.2 4.1 4.5 4.6 3.0 4.4 3.7 2.8 4.8 3.5 3.4 5.5 5.5 4.8
##  [415] 3.1 4.2 3.6 4.5 5.4 3.5 3.4 4.8 5.1 5.0 5.2 3.5 6.0 4.2 4.0 5.5 4.1 4.9
##  [433] 4.5 5.4 5.5 4.6 4.4 6.2 4.7 4.9 5.6 4.8 4.7 3.7 3.6 5.5 3.2 4.2 4.9 4.4
##  [451] 4.6 5.2 5.1 4.4 3.1 3.6 3.6 2.7 2.6 2.9 5.3 5.3 5.1 5.0 4.8 5.7 4.0 4.8
##  [469] 5.9 2.9 5.4 6.4 4.5 4.3 3.7 3.0 6.0 4.9 5.4 3.9 4.2 6.0 3.8 3.9 4.6 4.6
##  [487] 5.1 4.3 4.6 5.2 4.7 5.1 4.8 5.2 4.6 5.5 4.7 4.3 4.8 4.1 5.3 5.0 4.5 4.8
##  [505] 4.8 4.6 4.1 3.2 5.2 4.6 4.6 4.9 4.6 4.4 4.5 4.6 4.3 5.0 4.1 3.3 2.8 4.4
##  [523] 5.0 2.4 6.0 3.6 5.1 5.1 2.7 6.3 3.4 4.9 5.3 4.5 4.1 3.1 5.3 5.6 4.6 4.4
##  [541] 3.8 3.2 5.6 5.9 5.1 5.1 3.5 5.3 4.7 5.1 5.1 4.6 4.1 4.9 3.9 5.7 3.9 3.3
##  [559] 4.6 5.8 3.9 4.0 5.0 4.8 4.6 4.2 4.0 4.6 5.3 5.6 6.2 4.1 5.3 5.2 4.2 5.6
##  [577] 4.3 4.7 4.3 3.3 4.6 5.9 4.9 4.0 3.9 4.0 4.5 4.9 4.7 2.7 4.0 3.5 4.0 4.1
##  [595] 3.4 5.3 5.2 4.2 2.7 3.7 3.7 6.9 3.4 4.5 4.1 5.1 4.3 4.6 5.1 5.5 5.2 3.9
##  [613] 3.5 2.9 6.7 4.4 4.5 6.0 2.4 5.0 4.1 4.2 4.7 4.8 3.7 4.5 6.9 4.5 4.7 3.7
##  [631] 4.9 4.2 4.4 4.0 6.1 4.3 4.2 6.1 4.5 6.7 3.8 4.5 3.0 3.9 5.0 4.0 4.4 5.3
##  [649] 4.8 6.2 4.5 4.7 5.4 3.5 3.3 2.7 5.3 2.8 5.2 5.4 4.9 5.2 5.0 3.5 4.0 2.9
##  [667] 4.1 4.0 4.3 5.0 5.9 5.9 4.4 2.9 5.3 4.9 4.7 4.5 4.5 4.5 3.5 5.2 6.5 4.9
##  [685] 4.9 4.5 5.4 3.0 5.0 5.3 4.9 5.9 4.4 6.6 4.6 5.1 5.4 4.2 4.9 3.9 4.7 4.6
##  [703] 2.8 5.1 4.5 4.1 4.9 4.1 4.2 4.5 2.9 4.9 5.4 5.1 5.0 4.6 6.1 4.5 4.0 5.0
##  [721] 4.5 5.6 5.6 4.2 4.5 5.8 5.1 5.2 3.7 4.0 5.4 3.5 5.2 3.7 5.6 4.2 3.7 4.3
##  [739] 4.7 3.9 4.3 3.6 4.7 6.4 6.2 4.9 4.9 5.0 4.0 2.1 3.9 3.5 4.5 6.2 5.1 4.6
##  [757] 4.1 4.7 4.8 4.7 4.5 5.2 4.3 4.5 4.6 4.0 5.3 3.0 2.8 2.9 5.2 5.4 5.7 4.0
##  [775] 4.8 5.2 4.0 4.9 4.5 6.0 3.7 5.5 3.5 4.5 5.5 3.4 3.8 3.2 5.9 6.2 4.2 5.7
##  [793] 6.6 4.1 4.9 4.8 4.1 4.7 4.2 4.6 4.5 4.3 3.9 5.6 4.4 3.5 5.7 5.9 5.1 4.9
##  [811] 4.5 3.5 5.4 5.6 4.4 5.9 6.2 4.6 4.2 5.1 4.3 4.2 3.0 4.6 6.6 5.2 5.0 4.2
##  [829] 5.6 3.9 4.1 5.5 3.3 4.7 5.1 4.0 4.9 3.3 5.2 4.0 5.8 4.3 2.8 4.0 3.2 3.6
##  [847] 4.5 4.4 4.0 5.2 4.5 4.5 4.1 4.7 5.3 4.8 4.8 4.5 3.6 5.0 3.5 3.8 5.5 3.8
##  [865] 5.4 4.2 3.5 4.3 4.8 5.1 4.1 4.1 4.9 5.2 4.7 5.4 2.9 5.2 5.9 4.0 5.4 4.5
##  [883] 4.7 4.0 3.8 3.9 5.7 5.3 4.2 4.8 4.5 3.1 5.3 3.5 3.3 5.2 3.7 7.0 6.3 3.8
##  [901] 4.5 4.5 7.0 6.3 3.5 3.9 4.1 2.3 3.0 5.1 3.6 4.4 4.5 3.9 3.7 4.6 3.4 3.0
##  [919] 4.8 4.0 4.2 3.8 5.1 5.0 3.5 4.4 4.8 4.5 4.5 5.0 4.6 4.6 3.6 5.8 5.2 4.8
##  [937] 3.2 4.6 3.7 3.3 4.8 4.0 6.5 4.0 5.9 4.2 5.1 2.8 4.3 5.2 5.2 4.2 4.5 4.0
##  [955] 5.2 4.2 4.3 3.1 5.1 4.6 4.3 5.7 3.6 5.5 3.8 4.7 5.3 5.5 5.3 3.6 4.7 5.9
##  [973] 5.3 3.9 4.6 4.0 4.3 5.6 4.9 4.2 5.1 3.6 2.8 2.8 4.2 5.5 5.3 5.4 6.1 3.4
##  [991] 2.8 4.1 3.9 4.9 4.7 3.5 5.4 4.6 4.6 4.2
## 
## $func.thetastar
## 72% 
##   5 
## 
## $jack.boot.val
##  [1] 5.4 5.3 5.3 5.2 5.1 5.0 4.9 4.7 4.7 4.5
## 
## $jack.boot.se
## [1] 0.8637708
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
## [1] 1.364979
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
##   3.246828   5.421209 
##  (1.383927) (2.498934)
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
## [1] 0.6848393 0.1269251 0.9007643 2.1262062 0.9217344 0.9962905
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
##    [1]  1.9311390429  2.4859165030  1.0312632355  1.8875935875  0.6558621341
##    [6]  0.2787830175  1.4900182206  1.1241260352  1.3978841198  1.4257941840
##   [11]  1.2700713923  0.7721477875 -1.7033639858 -1.0771162322 -1.7716661542
##   [16]  1.9263005218  0.8374729357  0.6793825813  0.3585233040  0.1724706133
##   [21]  1.7365285409  2.4366901026  0.3372538826  1.9212378927  0.0368448975
##   [26]  2.5132822706  0.8728209316  0.6680393215 -0.7209503658  0.1151652601
##   [31]  0.7331266234 -1.6197803101  1.9168071913  1.2922380007  1.2286915115
##   [36]  0.3921974743  1.4987252558  0.6249708937  0.8746055773  0.4958253440
##   [41]  0.9529762108 -0.3429967222  1.8831309001  1.8681377223  0.7125961479
##   [46]  0.9472739719 -1.7617764749  0.6420599923  1.4493863812  0.2572380436
##   [51]  1.3696482306  1.1581418135 -1.9991848517  0.0197567479  1.1267861883
##   [56]  1.7244886758  1.4692246689  1.1545776109  1.3937550403  1.7107146784
##   [61]  0.0971084890  0.4394784649 -0.2433925815  1.6290147888  2.0008592017
##   [66] -0.8364480983  0.6379180537  1.8743883974  1.2238846716  0.2434110863
##   [71]  0.1998849518  1.9305363247 -0.0538855645  0.4580984881  1.3122112483
##   [76]  0.7624898313  1.0257003211  1.3972784876 -0.5656667914  1.9128730491
##   [81]  1.0406584778  0.3586028503 -1.9005206861  0.9239176862  2.0001860364
##   [86]  2.5154249696  0.1354152705  2.0149584902  1.1225279726  1.5531941148
##   [91]  0.1735027494  1.0722544742  1.9469510081  1.9304842270  1.4141056852
##   [96]  0.0499667711  0.8027901557 -0.1466340274  1.8831622741  1.2550082194
##  [101]  1.1388203741  1.0588745159  1.7979184198  1.1079794314 -0.1077842095
##  [106]  1.0157416238  0.0687106246  0.3166114112  0.7451074168  1.1512835195
##  [111]  1.7849370030  1.2154807436  1.1070155387  2.3324041077  0.8432400684
##  [116]  0.6203140094  0.1023646819 -1.8238514139  0.0697002054  1.4597099557
##  [121]  1.1645460167  0.2255052292  1.6674244974  1.1778483546  1.9126846843
##  [126]  1.1877690594 -1.4838967284  0.6399858770 -2.0395670146  2.5678589217
##  [131]  0.9331509824  0.6078255585  1.0279788926  1.1675714129 -0.0859186610
##  [136]  0.6703384610  1.2837833146  2.0561315259  2.4728807075  0.6465171138
##  [141] -2.0249953248  0.9917047339 -1.8932104582  0.5259594714  1.9867780212
##  [146]  0.7739101632  1.6865539907  0.0496684998  1.4496744156  0.7277014308
##  [151]  0.8721380262  0.8837536171  0.1862134959  0.2286310440 -1.9875394272
##  [156]  0.2278255686  1.1452767816  1.9843538891  0.2506045116  0.7271292180
##  [161]  1.1185978962  1.1194003766  0.8007861049  1.9543807232  1.4662390414
##  [166]  1.4727795953 -0.0895207553  0.7438287619  1.1783940578  1.2897165331
##  [171]  1.0056954097  2.2013751172  1.0049166564  2.5789115222  0.8755103552
##  [176]  0.0616243851  0.2768172214 -0.1569560402  0.8275350157  1.4796858737
##  [181] -1.9265115053 -1.2151093176  0.1277995808  0.2445251052  0.5206989438
##  [186]  0.8451317423  0.2986431365  1.0095768342  1.5020983445  0.6394472092
##  [191]  2.5204532759  1.2476733565  1.9224992637  1.1457660994  0.1992362095
##  [196]  0.6804463911  0.8354848207  0.9456633847 -0.2297786692  0.2930971370
##  [201]  1.4611552577  1.9279193498 -0.1584383045  0.1443068115  1.8275344020
##  [206]  0.3648074657  0.8386305904  0.2352536939 -1.0694929881  1.8829999811
##  [211]  1.1760438317  0.9168771444  2.5912381604  1.9036495347  0.4178191807
##  [216]  0.6501224069  1.8671675961  1.1283256983 -1.8843821381  1.8654011557
##  [221] -0.0670849731 -0.1502062427  1.3763151564  1.9580180891 -1.1501607847
##  [226]  1.8402047410  1.4437880782  1.2666457487 -2.1048607409  0.4039315674
##  [231]  0.7781341276  1.4363450075  1.4601884296  0.8996229063  0.8411788333
##  [236]  0.3211678649 -1.0779953401  1.1483385109  0.8682703770  1.1785184590
##  [241] -0.0006724696  1.3183600799  2.5513026508  0.2626818071  0.5913997492
##  [246]  2.5401447231  0.1297829203  1.8977403398 -1.8476538590  0.2166460855
##  [251]  1.8898037177  2.1991446404  1.9252335183 -0.0589694079  0.8142424719
##  [256]  0.9282889401  0.1505155112  0.2617318697  1.5144946954  1.9030610735
##  [261]  1.2007092688  0.2421203312  0.6599002086  1.4569177188  1.5760027421
##  [266]  0.7881753709  1.0755631675  1.9625782352  0.7534302435  1.5054267869
##  [271]  1.4520858117  0.8673951966  1.7716323506  1.4958845382  0.1533720375
##  [276]  2.5319262876  2.6081312985 -0.4101851962  1.3349719745  1.1771808831
##  [281] -0.6403453374  1.9192294450  1.9036495347  0.4430997451  2.5419894375
##  [286]  0.2207577233 -1.1717863256  0.3945345390  1.4510332689  1.5949178549
##  [291]  0.4171589566  2.6268171970  1.4413810582  1.1538677273  1.8687755079
##  [296]  0.1854562834  0.8740460898  1.6942869157 -1.8442117396 -2.1174196818
##  [301] -1.6959157988  0.0362906488 -0.0101040617  1.5164718465  1.4510332689
##  [306]  0.2783413517  2.4903360334  0.8881393799  2.4232555406  1.3599595929
##  [311]  1.0218921960 -1.0887445741  1.0237603218  0.1938300196  0.1226253227
##  [316]  1.4900258231 -1.9868241567  2.1047481075  0.6179532456  0.8185974650
##  [321]  0.6661870371  0.7056865528  1.6078646729  1.4270290737  0.6602443826
##  [326]  0.1504973654  1.1486207750  0.9159406415  0.7734443201  0.8202641421
##  [331]  0.5361237130  0.7550132508  1.4494506417  1.2176903000  0.5904420676
##  [336]  0.2103188189 -0.4014363678  0.4186120275  1.2818628455  0.7176775090
##  [341]  0.7080683567  1.8627273281  1.4144369066  1.9062391706  0.0235643464
##  [346]  0.3901608870  2.2333692669  2.5029106603  1.1947589286  1.4043094065
##  [351]  1.1048888744  1.0470271464  1.5980625922  1.5392807201  0.3774599719
##  [356]  1.4154157574  1.2956177588  1.4424359438  2.4439478528 -1.8915638126
##  [361]  0.6883747397  0.6568814965  1.1815259128  0.6298392759  1.5568849662
##  [366]  1.4722504896  1.0956815537  1.1517374961 -1.1985356336  0.3181494730
##  [371]  0.6654074898  0.6300429741  0.8318329091  0.6075014891  1.1583331734
##  [376]  0.0836161643  0.3313104266  1.8821985937  0.0637175751  1.5159759637
##  [381]  1.1449841231  0.0497723533  1.9874106278  0.9987235084  1.4361639408
##  [386]  1.1757627035  0.2975801832  1.3915918641  1.4778554104  1.6385253272
##  [391]  0.4015196147  1.1515013625  0.4903434544  0.3837481279  1.2951301314
##  [396]  0.4399545501  1.3009640458  0.8647756806  1.3990605782  0.6033916344
##  [401]  1.2787096781  1.0553296880  1.1203211054  1.4027834662  1.1751131126
##  [406] -0.3605777849  0.0694530120 -0.5742697245  0.7498844342  1.3291997554
##  [411]  0.3770538031  0.7320715746  1.6065579969  1.8933054714  0.6405001057
##  [416] -1.9190260577  2.6159563588 -0.0438069523  0.5897531290  1.1228092696
##  [421]  1.6320921461  1.8765149504  1.9861379085  1.2364781260  1.7183405944
##  [426]  0.8472777368  0.9405520246  0.9740319929  0.3514298393  1.3670596775
##  [431]  2.5048768236  0.4283274179  2.5437877999 -0.3142895097 -1.9770136567
##  [436]  1.7317981449  1.9817405227  0.8364648350  0.6754343457 -2.0159511227
##  [441] -0.3068423593  2.1043292632  1.2195385769  1.4506052908  1.4430282252
##  [446]  0.9441776097  0.4388509076  1.4395866688  0.7712730057  0.0127783165
##  [451] -0.1223691737  0.2710376679  0.1447928991  1.8993722490 -0.1132634040
##  [456]  0.1526068647  0.8340691943  1.5265119351  2.2391797960  0.1902316127
##  [461]  0.9813736206  1.9043445696  0.4252549483  1.6487786170  0.3708398602
##  [466]  0.7889026714  1.9192502025 -0.0925477697  1.8160461957  0.3266735075
##  [471]  1.3178034934  1.4489519897  0.4076371956  1.8855249478  0.9011862858
##  [476] -0.6596079332  0.4818858770  0.6838478223  0.5753565207  1.8998285189
##  [481]  1.8639132184  0.2436523871  0.9035563125  1.4256913204 -0.0822268776
##  [486] -1.4721308924  1.1432536661  1.4272469285  1.0169671958  0.9014061457
##  [491]  1.9430677738  0.9088719196  1.9238996362 -0.1125435467 -1.2700123015
##  [496]  1.4437880782  0.8519295974  2.1155860362  1.8130766778  1.3990444647
##  [501]  0.8941244051  0.5718812924  0.7441260440  1.2720262832  1.0993836982
##  [506]  1.4317564120  1.4179223051  0.3862562943  0.6677333198  0.4323245679
##  [511]  0.6831655321  0.6492798211  0.1266369005  1.1140780357  0.0057250565
##  [516]  0.6471587946  1.0540973514  0.9367234589  1.9203419146  0.6028947505
##  [521]  0.6641461600  1.9161518818  0.4148900586  0.0805934765  1.1363104855
##  [526]  1.1694406706  0.7757893279  0.6736684974  1.9435463211  0.9132484790
##  [531]  1.1701694478  0.3774071979  0.9592911498 -1.1673673367  0.8326992829
##  [536]  1.9525834250  2.4498604168 -1.9373193838 -1.2169912931  0.1194746459
##  [541]  1.3220389605  0.1144617592  1.2800412126  0.7402123135  0.7449799939
##  [546]  1.2292493963  1.5650820674 -0.1558091532  0.9446003126  0.3360164885
##  [551] -0.1962135649  0.8516912304  0.9921951213  0.1759750528  2.5008210444
##  [556] -1.6940235254  2.1240805122  2.2377827666  1.9440867040  1.9641313456
##  [561]  1.1667355200  2.2349217804  1.9964864911  0.8873334806  0.8933505407
##  [566]  0.1135049687  0.6835651106  0.8638185638  1.7431142771  0.9432000935
##  [571]  1.0611913621  1.1565703029  1.1697380501 -0.2084013930 -0.0783788057
##  [576]  0.9667703539  0.3902939098  2.1576289402 -1.8690867768  1.1610491353
##  [581]  0.3601832927  1.1289068726  0.9896602124  1.6737332196  1.1067315855
##  [586]  1.6861984128  1.8116958014  1.8059011269  1.2107793791  0.4458209506
##  [591] -1.8897255537  1.4543708472 -0.1333060182  0.1976358266  0.9479959699
##  [596]  0.1897879232  0.8479330456 -0.3678827926  0.1285321504 -0.7113474709
##  [601] -1.0210995558 -2.2283662214 -0.0472686466  0.2118979741  1.4092177793
##  [606]  0.7444134944 -0.6693383098  0.4449645479  1.4412205090  1.9005911795
##  [611]  1.1339269761  1.8965137966  0.5833604663  1.3378905438  1.1365941554
##  [616]  1.9278253153  0.8254563118  0.8052587689  1.1038329114  2.4561596477
##  [621]  2.3989773706  2.4421113288  0.9294755771  0.8876808927  2.5974456757
##  [626]  0.8194962605  1.8152171865  2.4488760747  1.8657548642  2.1342607112
##  [631]  1.9425515331  1.3337932826 -0.1158328195  0.1572865717  1.0680402238
##  [636] -1.1841377403  0.6398527699  0.9186150183  0.3670573098  0.0097149402
##  [641]  1.4607065698  1.8620800798  2.5543768841  0.1708794715  0.6334731822
##  [646]  0.5945135405  0.2118979741 -1.1524234950  1.5472095897  1.9304842270
##  [651]  0.8625147055  1.4682197501  1.0658214589  1.1780157232  0.0420480591
##  [656]  0.4192543781  0.2514849985 -0.0023411092  1.8594435184  1.9590796653
##  [661]  1.9236004018  1.1568877814 -1.8418392689  0.5713485397  2.4904673302
##  [666]  0.9294755771  0.7101516523  2.1047481075  0.5849872489  1.8445598828
##  [671]  0.3602476496  1.0140347123  1.9066142272  0.4363809978  0.6588802862
##  [676]  1.1674904037  2.5447423755  0.8383289635  1.1341251126  1.1785184590
##  [681]  0.6703277212  1.4553007978  0.9189371278  0.7445607038  1.3542034107
##  [686]  2.0675775424  0.4227639284  1.0410286159  1.7149215792 -0.0753595231
##  [691]  1.1328810883  1.3863643930  1.1460051081  0.6223036492  1.9038279680
##  [696]  1.9144671928  0.5904420676  2.4592656008  1.0447764530  0.0908416710
##  [701]  0.5843894490  1.0675363287  1.1323820059 -0.5988220237  0.2297496771
##  [706]  1.5046477723  1.4526147321  0.1591817319  1.1406871539  0.2474833115
##  [711]  1.8770853092  1.1624802437  2.0660644782  1.0553803046 -0.0582915200
##  [716] -0.6405405437  0.1258509107  0.6938611792  0.6541849939  1.1116977967
##  [721]  1.5172367250  0.8921233605  1.8544857389  0.6532088906  1.0119955019
##  [726]  0.1718293428  1.4002093585  0.2676980355  1.8899025358  1.4287456014
##  [731] -1.4854852526  0.1015007155  0.6319012761  1.0103110488  1.1016018143
##  [736]  0.9083983584  0.8167796805 -1.0957292069  1.2145408879  2.2979683722
##  [741] -1.4760933550  1.9486145200  0.6698458825  0.7311844058  1.4659518130
##  [746]  0.0093508673  0.0396547908 -1.7816365234  0.7842255397  1.4441814316
##  [751]  1.9082018527  0.9018032667  0.4109738808  1.8857749818  0.5438974748
##  [756]  0.2121566185  0.5153497650  0.3494847153  2.5734924365  1.0633622524
##  [761]  0.8466628182  0.4226988046  1.3942940213  1.1343805853  1.1631812956
##  [766]  1.9645550927  1.9373459959  1.4118984385  1.6809272475  1.0174809210
##  [771] -0.0695609313  0.0628781068  1.9511309720  0.6100287635 -0.5658006079
##  [776]  1.3798429158  1.1545003990  0.6227856028 -0.0150776970  1.4628814927
##  [781]  1.0611265825  1.1425091768  0.2994830041  0.3710597322  1.9063786029
##  [786]  0.6475568117  1.4717105355  1.2462603428  1.9148260987  1.1515458892
##  [791]  1.7006805839  2.0615852038  1.1399134101  0.8468573536  1.2105049236
##  [796]  2.6073332840  0.5417161374  0.0909683345  0.2106293570 -1.7797666730
##  [801]  0.4559731119  0.4065106951  1.1669680501  1.2919491227  0.6964834847
##  [806]  1.3695060392  1.0754067654  1.2083151357  1.4529627848  1.4400100669
##  [811]  0.7176899902  1.3390204374  0.4789050967  2.2486921740  0.0742224191
##  [816]  1.0580673055  0.6737892737  1.3577988526  0.3965630368  1.1466933635
##  [821]  0.6616042308 -1.8907632215  1.4411789757  1.1334603165  1.7227786920
##  [826]  0.9462881085  1.1770756106  1.8279473977  0.2910667094  1.0848228831
##  [831]  0.8922030164  2.4807562768  0.6472679922  0.8533926746  0.6894740993
##  [836]  1.2893883553  0.8177322749  1.9205271647  0.6357471483  0.9089344904
##  [841]  0.8057829883  1.5179254060  1.2176114928  0.6477503720  0.8720669308
##  [846]  1.2110235459 -1.8442331847  1.7431142771  1.1515729477  2.4899218181
##  [851]  1.4476312005  0.6538803153  1.9330933876  0.5215161298  1.8612443839
##  [856] -0.6033669448  0.5776415645  1.1713869342 -1.3440324359  1.3994406759
##  [861]  2.5528851275  0.1944807737  0.3377279640  1.3710379924  1.1485671532
##  [866]  0.3081096622  0.4905900280  1.8274514185  1.1326159799 -2.1339938525
##  [871]  0.5932132764  1.9575283185  0.3187781678  0.1886540923  1.1266589012
##  [876]  1.1396758990  0.8420384010  0.2976921640  0.9518751576  1.9967012167
##  [881]  0.5899997455 -0.8346536890 -1.1707847729  1.2620625730  1.2329960822
##  [886]  0.3076764582  0.1994315128  1.9727356134  0.8664933481  1.9478695144
##  [891]  1.1429345096  1.1395869724 -0.1323177024 -1.3429702495  0.6078255585
##  [896]  1.4624517150  1.9790047978  1.4067504810  0.6946903005  0.8914995537
##  [901]  1.9231246617  1.1591973422  0.7880067736  0.5376643657 -0.1159088644
##  [906]  2.1036146849 -1.1124951898  1.7382829740  0.9354361391  0.6067563924
##  [911]  1.1592242236  1.9734284578  1.5604496025  2.4586381915  0.6047279069
##  [916]  0.9174750629  2.4756148674  1.8595520167  0.8523679824  2.5181237990
##  [921]  1.4760585122  1.2928308085  2.4788730020  1.0896240318 -0.1765026291
##  [926] -1.1051356209 -0.0503457974  0.8943862523  0.9051609109  1.1893467409
##  [931] -0.0636933464 -1.5462817742  1.0520093548  1.3467715405  1.9760328598
##  [936] -2.2604065241  1.4637647323  0.7445233050  1.6612564636 -0.0153758656
##  [941]  0.5759487115  1.1699813756  0.8432500611 -2.0306655416  0.6037530699
##  [946]  1.3264130969  2.5155395605  0.9479959699  1.4843369072  1.1690785669
##  [951]  2.3655853467  0.3124377314  0.6317344947  2.0870201053  0.6146724742
##  [956]  0.0796828888  1.8601882462  1.4666300474  0.3829318161  0.5455161047
##  [961]  1.1051782853  0.6428004175  1.9068296123  1.1840997389  1.2912231493
##  [966]  0.9433525540  1.2191131669  1.0174931334  0.9916444831 -0.7834410572
##  [971]  1.9844836007  1.1676640384  0.2076210014  1.1292539844  1.4439485568
##  [976]  1.5289396478 -1.1338721190  0.3773360189  1.5424170915 -0.0750200854
##  [981]  2.5012746807 -0.6533740501  1.3306734898  1.1304616453  0.2784362965
##  [986]  0.8180578237 -0.5517225993  1.8401352113  1.8839482154 -0.1008998732
##  [991]  1.4846032528  1.5636345750  1.0496853744  0.1332506988  1.9266841069
##  [996]  2.4653949746  1.9047574351  0.9869257460  0.9352237056  0.9047261094
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
##   0.59890863   0.34112233 
##  (0.10787235) (0.07627257)
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
## [1]  0.002498169  0.722497635  0.124881072 -0.327823470  0.669293445
## [6] -0.294282441
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
## [1] 0.0348
```

```r
sd(results2$thetastar)  #the standard deviation of the theta stars is the SE of the statistic (in this case, the mean)
```

```
## [1] 0.8934916
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
## t1*      4.5 -0.01741742   0.8530019
```

One of the main advantages to the 'boot' package over the 'bootstrap' package is the nicer formatting of the output.

Going back to our original code, lets see how we could reproduce all of these numbers:


```r
table(sample(x,size=length(x),replace=T))
```

```
## 
## 1 2 4 5 6 8 9 
## 1 2 1 1 3 1 1
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
## [1] -0.0171
```

```r
se.boot
```

```
## [1] 0.9188917
```

Why do our numbers not agree exactly with those of the boot package? This is because our estimates of bias and standard error are just estimates, and they carry with them their own uncertainties. That is one of the reasons we might bother doing jackknife-after-bootstrap.

The 'boot' package has a LOT of functionality. If we have time, we will come back to some of these more complex functions later in the semester as we cover topics like regression and glm.

