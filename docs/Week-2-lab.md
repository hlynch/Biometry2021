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
## 0 1 2 4 6 9 
## 3 1 2 1 1 2
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
## [1] 0.006
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
## [1] 2.705029
```

```r
UL.boot
```

```
## [1] 6.306971
```

Method #2: Simply take the quantiles of the bootstrap statistics


```r
quantile(xmeans,c(0.025,0.975))
```

```
##  2.5% 97.5% 
##   2.6   6.3
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
##    [1] 5.4 6.7 4.8 3.9 4.4 3.8 4.9 5.3 4.3 5.1 4.9 3.3 4.4 4.9 5.6 3.4 4.7 4.4
##   [19] 5.2 4.2 4.3 5.2 4.3 4.1 5.1 3.0 3.3 3.2 4.9 3.8 4.0 6.3 4.3 6.0 4.3 2.8
##   [37] 6.2 2.8 4.9 3.7 4.2 4.6 5.6 4.8 4.3 2.5 3.9 7.0 4.1 4.5 3.6 5.8 4.7 3.9
##   [55] 3.7 4.7 4.3 5.4 4.5 3.7 4.7 5.6 5.4 5.2 4.2 5.4 5.7 5.5 4.3 5.3 4.1 4.7
##   [73] 6.3 4.4 4.9 3.3 4.0 4.0 4.5 6.1 4.8 5.4 6.6 4.1 5.0 4.1 4.2 5.7 3.3 2.6
##   [91] 5.0 4.2 4.6 3.5 3.8 4.1 4.0 5.3 4.2 5.8 5.3 5.5 4.4 4.1 2.8 3.9 4.8 4.0
##  [109] 4.1 4.9 5.0 4.6 5.6 5.5 5.6 5.2 5.6 4.1 4.4 2.9 4.5 4.9 4.0 3.8 5.6 4.5
##  [127] 4.0 4.2 3.9 4.5 4.5 4.4 4.2 7.1 4.1 5.6 3.6 4.1 4.0 5.9 4.7 4.8 5.0 3.9
##  [145] 4.9 4.0 5.1 4.4 4.2 5.9 4.0 5.4 5.6 4.5 3.1 5.3 4.8 5.3 4.3 4.4 5.9 5.8
##  [163] 4.0 5.2 4.9 3.5 5.4 4.0 5.6 3.6 4.9 2.2 4.3 3.2 5.0 4.7 6.5 5.2 5.3 4.9
##  [181] 3.7 5.1 6.0 5.1 4.4 4.3 4.3 4.5 4.2 4.5 3.2 4.1 4.8 3.5 4.6 4.1 4.7 4.2
##  [199] 2.9 3.5 4.7 4.9 3.7 4.6 3.6 2.9 3.9 3.9 5.4 4.5 3.6 5.3 4.3 5.9 5.2 4.5
##  [217] 6.3 4.4 3.3 3.7 5.9 4.4 5.8 4.9 2.1 5.8 5.1 4.6 3.4 3.9 1.8 4.2 3.9 4.7
##  [235] 5.0 4.8 5.2 5.5 5.0 5.3 3.7 4.7 2.7 3.6 3.9 4.6 2.8 3.3 5.4 4.4 3.8 2.1
##  [253] 4.1 4.5 4.0 2.8 5.7 4.5 6.1 4.3 4.1 2.3 4.7 4.7 6.3 4.5 4.5 5.0 4.8 4.8
##  [271] 6.1 4.9 4.4 2.9 4.3 4.5 4.0 5.9 4.3 5.4 3.7 3.8 4.1 4.7 6.4 5.0 3.9 2.7
##  [289] 3.4 6.6 4.3 5.0 3.8 5.2 4.8 5.3 3.6 4.8 5.1 4.9 6.5 2.6 3.6 4.0 5.8 4.8
##  [307] 4.1 4.1 3.6 5.6 5.4 5.8 4.0 5.5 5.8 4.8 5.3 4.3 4.1 3.0 5.0 5.1 4.7 4.5
##  [325] 4.8 4.7 3.2 3.4 4.3 3.3 6.0 4.9 3.9 4.9 2.4 5.5 5.1 5.9 4.6 4.4 4.6 4.2
##  [343] 5.3 4.9 5.7 4.9 4.0 5.1 4.6 2.8 4.0 5.2 3.7 5.5 4.3 3.0 6.2 5.5 5.3 3.5
##  [361] 4.7 5.5 4.2 4.4 4.6 4.4 4.4 3.3 3.0 4.0 4.5 5.8 3.7 5.0 5.4 3.2 4.7 4.8
##  [379] 4.0 6.0 5.0 4.7 3.5 4.8 2.8 5.1 4.0 5.0 5.3 5.3 4.3 3.1 5.0 5.5 5.0 6.1
##  [397] 4.7 3.5 4.7 4.5 3.5 4.7 6.3 3.6 3.7 3.9 3.4 4.9 4.0 4.3 3.0 5.8 5.4 3.7
##  [415] 5.0 4.6 3.5 4.7 4.4 5.2 4.3 4.6 6.2 5.7 4.4 4.3 3.7 4.1 5.5 3.6 3.9 5.1
##  [433] 4.7 4.1 4.7 3.7 4.1 5.5 3.2 4.8 3.5 4.5 3.1 4.6 3.5 5.6 5.7 4.0 3.0 3.9
##  [451] 4.2 4.8 2.9 4.9 4.4 2.8 5.1 6.2 4.7 4.4 4.2 4.2 4.8 4.3 4.2 5.1 5.1 3.6
##  [469] 3.1 4.4 4.4 5.9 5.0 5.1 6.0 4.1 3.9 3.1 4.1 3.8 5.3 5.6 3.5 3.9 3.3 4.4
##  [487] 3.6 5.0 4.9 5.2 3.0 5.4 4.5 4.5 4.6 4.5 5.5 4.3 4.6 5.2 4.3 5.0 4.5 4.3
##  [505] 6.6 3.3 3.5 4.4 5.0 3.9 4.2 5.5 5.0 5.1 6.6 5.2 5.5 4.1 5.5 4.4 5.3 4.9
##  [523] 3.5 3.2 5.3 3.8 4.0 5.5 4.8 3.7 3.9 4.8 4.4 3.8 3.4 4.2 4.9 3.8 4.4 4.7
##  [541] 3.9 4.4 3.8 5.4 5.3 5.2 4.5 3.1 3.9 6.2 5.2 5.4 3.9 4.8 5.9 6.7 5.2 3.2
##  [559] 4.9 4.1 2.5 5.8 3.4 5.8 3.5 3.6 4.7 5.9 5.1 4.8 5.3 5.4 5.1 3.0 6.2 5.2
##  [577] 6.1 4.0 5.1 5.1 2.6 3.7 4.7 4.2 4.0 4.6 5.0 5.0 4.9 3.9 4.3 4.5 4.8 5.6
##  [595] 3.9 4.5 4.3 5.3 4.8 5.5 3.5 3.5 3.6 4.9 3.4 5.0 4.7 5.3 5.0 4.0 4.6 5.6
##  [613] 4.9 4.6 4.3 5.4 4.0 3.5 5.4 6.0 4.6 5.8 5.5 5.0 4.9 5.0 3.7 4.7 4.0 3.8
##  [631] 6.8 5.8 6.0 4.5 3.8 5.2 2.8 4.1 4.2 5.3 5.2 5.1 4.5 6.9 4.9 4.8 5.5 4.4
##  [649] 4.4 5.9 5.3 4.4 3.0 4.9 4.2 4.6 3.0 5.2 4.4 5.2 4.4 4.6 3.8 5.4 5.0 5.1
##  [667] 3.0 4.4 4.6 4.2 3.7 3.6 3.4 6.6 4.7 4.5 5.9 4.1 4.2 3.2 5.2 3.1 4.3 5.2
##  [685] 5.0 4.0 5.2 4.7 5.1 4.0 5.5 4.3 5.6 3.4 4.3 4.5 5.4 5.2 4.7 4.3 4.7 4.7
##  [703] 4.9 3.8 4.7 5.2 4.6 5.8 3.9 3.4 3.3 4.4 4.9 6.6 5.7 3.1 3.7 4.0 5.0 5.0
##  [721] 4.6 3.1 5.6 5.5 3.3 4.8 5.0 4.2 5.3 4.3 3.7 5.8 4.5 4.3 3.9 4.1 3.7 5.7
##  [739] 4.9 4.3 3.6 2.6 4.0 6.1 3.6 4.1 5.5 4.1 4.6 5.1 4.2 3.2 4.5 5.1 4.5 4.0
##  [757] 5.7 3.9 4.5 5.4 3.3 4.5 4.9 4.7 3.9 4.7 5.0 4.7 5.3 4.2 4.3 4.2 6.1 4.2
##  [775] 2.9 6.0 4.4 4.3 5.2 4.8 5.0 5.0 4.6 4.0 4.5 4.3 4.6 5.9 5.8 4.4 5.7 4.4
##  [793] 5.4 2.3 4.2 4.3 4.5 3.3 5.5 5.3 4.6 5.4 3.6 5.9 3.0 3.6 4.4 4.3 6.0 4.2
##  [811] 4.3 4.8 3.8 3.7 3.6 4.4 4.0 5.5 3.2 5.4 5.9 5.7 5.3 5.7 3.9 6.8 4.7 5.5
##  [829] 6.1 5.0 3.3 4.0 4.5 5.9 3.4 3.5 4.8 4.5 4.2 4.6 4.5 5.2 4.5 5.0 4.2 5.6
##  [847] 4.5 4.1 5.5 5.1 1.9 4.2 3.6 3.2 5.3 4.9 3.2 6.5 4.3 4.1 3.7 5.3 3.6 3.1
##  [865] 3.9 4.5 5.1 2.6 4.5 4.1 4.2 4.5 5.4 3.7 4.6 3.7 3.8 6.4 3.8 5.2 5.2 4.1
##  [883] 3.7 4.3 4.8 3.9 3.9 2.7 3.9 4.3 3.0 3.5 5.7 4.9 4.4 4.4 5.0 2.7 5.2 4.6
##  [901] 5.3 3.3 3.6 4.1 5.1 4.1 5.7 6.0 5.9 4.8 4.3 4.6 3.7 5.6 4.2 5.1 6.2 4.5
##  [919] 3.9 4.3 4.5 3.4 4.5 3.0 3.8 4.3 5.7 3.5 3.6 3.8 4.6 4.5 2.9 4.4 3.9 3.1
##  [937] 4.4 4.8 5.3 3.6 4.7 4.0 5.8 4.8 4.3 2.9 3.6 2.8 6.3 4.2 6.4 3.1 3.7 5.6
##  [955] 2.8 5.5 5.4 5.4 4.5 5.3 4.7 4.7 4.4 4.6 4.5 3.2 3.6 5.2 5.5 4.5 3.4 4.6
##  [973] 4.5 4.0 5.0 5.0 4.8 3.2 5.1 5.1 4.1 4.0 5.2 4.1 4.2 5.4 4.4 5.2 4.4 5.4
##  [991] 4.2 4.4 4.2 6.2 4.4 2.7 4.4 3.6 5.2 1.9
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
##    [1] 4.0 4.9 4.8 4.4 4.5 2.5 1.7 6.7 4.5 4.0 4.3 5.4 4.1 5.7 4.1 6.7 4.0 5.4
##   [19] 4.4 5.3 5.2 4.8 3.3 4.9 5.9 4.8 4.6 4.3 5.6 4.6 3.3 5.0 2.3 4.5 3.6 5.8
##   [37] 4.9 5.2 4.9 4.3 3.4 4.7 5.8 3.3 5.5 5.3 4.4 3.8 4.9 4.7 4.3 4.7 5.7 3.7
##   [55] 4.3 3.9 4.7 4.4 5.5 3.8 3.7 3.4 3.5 5.0 3.2 6.1 5.3 5.6 4.6 5.6 4.3 5.0
##   [73] 3.8 5.3 5.0 4.5 5.6 3.8 2.9 4.8 4.4 4.9 4.1 4.0 3.6 3.2 5.0 5.8 5.3 4.6
##   [91] 3.0 4.5 4.5 4.2 5.6 2.3 4.7 5.9 4.2 3.8 3.7 4.5 3.7 3.9 3.2 5.5 5.2 3.1
##  [109] 5.0 4.2 3.1 4.1 2.8 2.8 6.0 2.9 3.5 5.6 5.4 5.5 5.5 3.4 5.2 5.1 5.2 4.1
##  [127] 3.8 3.0 4.6 4.6 4.3 6.6 5.1 5.1 5.4 3.8 4.3 5.8 3.8 5.1 3.8 4.4 5.1 4.2
##  [145] 4.4 2.7 3.5 5.7 6.0 5.7 5.1 4.2 3.5 3.4 5.5 5.0 3.7 4.6 5.2 5.4 3.0 4.8
##  [163] 6.1 3.4 4.9 5.2 4.6 4.0 5.4 5.1 4.4 4.1 2.9 4.7 3.7 3.4 4.0 6.1 4.8 4.7
##  [181] 5.8 3.1 6.0 3.1 5.7 4.8 4.1 4.6 3.5 3.8 3.6 5.1 4.2 4.3 5.4 4.6 3.0 4.4
##  [199] 4.5 5.1 5.4 3.2 4.6 4.9 4.4 3.5 4.9 4.3 4.9 3.3 3.6 3.9 2.9 5.4 5.6 5.1
##  [217] 5.1 4.6 2.7 6.0 4.5 5.5 5.4 5.4 4.3 2.2 3.8 5.1 4.3 4.7 4.2 5.6 4.5 4.9
##  [235] 7.0 4.1 5.1 5.8 6.4 5.4 5.4 6.4 3.5 5.0 5.8 5.1 4.9 6.0 3.3 4.9 3.7 3.5
##  [253] 4.5 4.2 6.4 2.7 3.6 5.7 2.7 3.5 3.5 4.0 2.9 4.1 5.2 6.1 3.5 3.8 3.4 3.2
##  [271] 4.5 4.3 4.5 3.4 4.0 4.0 5.7 4.3 4.1 3.6 3.3 3.2 4.6 4.0 4.6 3.0 4.6 5.4
##  [289] 3.6 4.5 4.1 5.0 4.3 5.5 5.1 2.6 4.5 5.0 2.0 4.7 5.3 4.7 4.9 4.2 3.6 5.9
##  [307] 4.3 4.1 4.1 5.5 3.8 5.5 5.3 6.2 4.3 5.1 3.5 5.0 5.7 4.0 5.1 5.5 5.3 3.2
##  [325] 3.3 6.3 6.1 5.9 5.1 3.9 3.3 3.7 4.6 4.3 5.1 4.0 4.3 5.3 6.8 4.4 5.2 5.1
##  [343] 1.7 3.6 5.7 5.9 5.4 5.4 5.3 3.6 3.4 4.5 3.9 6.0 4.0 3.9 5.6 6.2 6.4 4.2
##  [361] 3.8 3.7 4.2 5.7 4.4 3.3 2.7 4.5 4.4 4.1 6.1 4.4 4.1 4.6 4.3 4.6 4.0 6.1
##  [379] 3.1 4.2 5.5 4.6 5.3 4.2 5.8 5.3 5.6 3.7 2.8 4.0 5.1 6.5 5.2 4.5 4.4 4.2
##  [397] 2.8 3.3 4.9 4.1 5.6 4.0 3.8 3.0 4.9 3.7 4.2 4.8 3.4 4.7 5.8 4.0 3.4 4.8
##  [415] 4.5 3.9 4.8 4.7 5.2 3.5 6.0 4.1 7.0 5.6 5.1 4.2 3.9 3.8 3.7 4.3 4.3 3.9
##  [433] 3.6 4.8 5.5 3.3 5.1 3.5 5.5 5.2 4.1 4.2 4.9 4.5 3.9 4.4 4.9 4.0 5.6 3.6
##  [451] 4.3 4.8 3.1 3.9 4.8 3.8 5.5 3.6 3.8 4.4 4.1 4.0 4.0 4.6 4.6 6.2 4.3 5.5
##  [469] 5.4 5.4 5.3 4.2 4.3 2.5 5.5 2.4 3.2 4.7 5.8 4.5 3.5 5.0 4.4 4.4 4.3 4.9
##  [487] 3.5 4.8 2.9 3.4 4.7 5.5 3.8 3.8 4.2 5.4 4.9 4.6 5.3 4.2 7.2 6.3 3.7 2.8
##  [505] 3.3 4.8 3.4 5.4 2.8 3.2 4.5 5.2 3.1 5.1 6.4 3.5 3.4 6.1 5.1 3.0 5.3 4.0
##  [523] 4.9 3.0 4.4 3.8 4.1 3.7 3.6 4.5 5.6 3.2 3.5 4.4 5.9 5.1 5.0 4.5 4.0 3.9
##  [541] 3.0 5.2 3.0 5.6 5.0 4.9 4.8 4.7 4.5 4.9 2.7 4.4 4.0 5.7 4.6 4.7 4.8 4.3
##  [559] 2.6 4.7 5.1 5.0 5.1 4.3 2.8 3.4 4.4 4.8 6.0 4.4 6.3 5.8 4.6 5.2 4.9 5.0
##  [577] 3.6 4.0 5.0 4.6 4.8 4.5 4.9 3.9 5.9 3.6 5.5 4.5 4.4 3.1 3.4 4.6 5.1 4.2
##  [595] 4.2 3.4 3.0 5.5 3.6 3.8 3.6 4.1 5.8 5.0 6.0 4.2 5.4 5.3 4.8 3.6 4.2 4.7
##  [613] 2.9 5.7 5.8 3.8 4.2 5.3 3.9 4.2 4.7 4.6 5.3 3.9 5.4 4.7 3.8 4.5 4.7 4.8
##  [631] 3.4 3.9 3.1 4.3 6.4 3.8 2.6 4.0 3.4 4.8 4.8 4.6 4.7 3.4 3.7 4.3 5.5 6.4
##  [649] 3.3 3.3 3.5 4.8 4.6 6.4 3.6 5.5 5.6 5.4 5.2 5.5 4.3 4.9 5.5 5.0 6.1 6.2
##  [667] 4.1 5.5 5.1 5.8 4.3 4.3 4.2 4.5 5.9 4.4 4.2 4.2 5.2 5.3 5.1 4.3 4.4 4.5
##  [685] 5.8 3.2 2.4 3.8 5.2 4.4 4.8 4.6 5.8 4.0 3.6 4.7 4.9 3.6 3.4 5.0 5.2 6.4
##  [703] 4.0 5.6 5.2 3.3 3.9 4.1 5.3 4.0 3.3 4.5 4.3 5.2 3.1 7.1 4.2 5.1 3.4 5.3
##  [721] 3.6 4.1 5.3 4.6 3.9 5.3 3.6 4.2 5.5 3.2 4.4 5.0 3.1 5.5 3.8 4.2 3.5 5.8
##  [739] 4.5 5.7 4.1 5.2 4.0 5.6 4.9 4.1 3.8 5.0 4.6 5.2 3.9 3.4 4.3 4.6 4.4 3.7
##  [757] 4.3 4.8 6.0 3.9 4.5 6.1 2.6 4.5 3.9 6.3 5.5 4.2 4.4 4.5 4.6 3.9 3.8 4.8
##  [775] 3.6 4.4 4.0 5.1 4.6 4.6 3.5 4.2 4.0 4.9 4.1 4.7 3.8 3.8 5.8 4.3 4.3 5.7
##  [793] 5.1 3.9 6.4 5.3 5.5 3.3 6.0 4.8 3.7 5.8 4.4 3.0 4.9 4.7 4.9 5.7 5.8 4.5
##  [811] 4.0 6.1 4.6 4.7 4.3 4.9 5.0 4.1 5.6 2.9 5.0 3.3 5.0 4.4 4.5 2.9 4.3 3.4
##  [829] 3.2 3.7 6.9 4.8 5.1 4.2 3.6 3.4 4.9 5.4 3.8 4.2 2.2 2.7 4.3 5.0 5.2 4.3
##  [847] 5.3 3.0 5.7 4.7 6.8 5.9 3.7 6.0 5.5 4.2 3.3 2.9 4.6 6.1 2.7 3.8 4.7 2.9
##  [865] 4.7 3.4 3.8 1.7 4.1 3.3 4.4 4.3 6.5 4.9 3.6 3.2 3.3 5.7 3.1 4.0 4.8 4.7
##  [883] 4.2 4.2 4.6 5.4 4.3 5.9 4.4 3.2 2.7 5.6 4.7 5.0 5.5 4.1 3.0 5.5 4.9 5.6
##  [901] 3.6 3.3 3.4 3.3 3.4 4.3 4.0 1.9 5.5 4.0 5.0 5.1 4.0 3.6 4.8 4.3 3.3 4.9
##  [919] 4.8 3.5 3.4 4.9 6.1 4.3 4.7 6.1 4.9 6.7 6.3 5.7 3.9 3.9 4.1 5.5 4.0 4.8
##  [937] 4.2 3.8 3.1 2.7 4.2 3.4 3.1 5.0 4.8 3.4 3.1 4.1 5.7 4.7 4.8 4.6 5.1 4.0
##  [955] 4.4 5.3 5.4 6.5 4.2 3.3 3.4 3.6 3.7 5.3 3.5 5.5 4.2 3.9 4.0 4.3 3.2 4.1
##  [973] 5.2 3.2 4.2 4.1 4.7 5.4 5.6 3.7 2.9 5.0 4.9 6.0 4.5 6.2 3.9 4.0 3.8 4.6
##  [991] 4.7 4.9 6.1 4.0 5.1 4.9 3.9 3.7 4.8 5.6
## 
## $func.thetastar
## [1] -0.0125
## 
## $jack.boot.val
##  [1]  0.543055556  0.448459384  0.258571429  0.209523810 -0.003804348
##  [6] -0.113314448 -0.242028986 -0.233423913 -0.402005731 -0.497109827
## 
## $jack.boot.se
## [1] 1.013506
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
##    [1] 3.2 6.4 4.5 5.3 4.2 4.1 3.8 3.1 4.8 3.1 3.9 5.3 6.5 5.6 3.1 4.8 5.0 4.5
##   [19] 4.2 4.3 5.5 2.8 5.6 4.8 4.2 3.3 4.5 5.0 6.2 5.1 5.8 3.8 4.4 5.0 5.6 5.0
##   [37] 4.1 4.3 4.8 5.2 6.6 2.2 4.9 4.0 5.1 5.0 2.8 4.9 5.5 4.5 5.3 4.9 4.8 5.4
##   [55] 5.0 3.7 3.1 4.3 3.1 6.6 6.3 3.7 3.2 3.1 5.4 4.5 5.1 6.4 3.5 3.7 3.8 4.6
##   [73] 4.4 5.0 3.0 4.3 6.0 6.0 5.4 4.5 3.6 6.2 4.3 3.9 3.7 4.0 5.5 4.4 4.5 3.2
##   [91] 4.9 4.9 3.6 4.5 4.8 4.2 3.7 2.6 4.0 3.8 4.2 5.9 4.1 7.0 5.1 3.6 5.3 2.9
##  [109] 3.0 4.4 3.5 5.4 4.5 3.6 5.5 4.3 5.6 4.7 4.5 4.5 5.0 5.2 4.1 2.6 6.1 3.6
##  [127] 3.5 2.6 4.5 5.8 5.3 4.9 4.3 3.5 3.8 3.2 5.9 3.3 3.1 6.4 6.1 3.3 5.0 4.3
##  [145] 3.9 5.0 5.8 6.5 5.6 4.5 5.6 4.2 5.2 3.1 3.6 3.7 4.3 4.1 3.6 5.8 4.9 2.9
##  [163] 3.3 4.1 4.6 6.5 4.2 5.6 4.0 5.7 5.6 4.9 5.3 4.1 4.3 4.0 2.8 5.2 6.0 4.5
##  [181] 3.9 4.5 3.2 3.9 5.9 4.4 3.4 3.1 5.8 7.0 4.0 4.7 4.3 3.7 5.4 2.0 3.4 5.6
##  [199] 5.3 5.3 5.4 3.2 6.7 4.3 3.7 4.2 4.7 3.0 3.5 2.7 5.0 4.8 6.2 4.8 3.9 4.5
##  [217] 4.8 5.2 5.8 5.0 4.1 5.1 6.0 2.6 3.0 5.0 4.8 3.1 4.1 3.2 4.3 6.0 3.4 5.2
##  [235] 3.3 2.7 5.1 3.9 4.3 4.3 5.1 4.3 5.3 4.6 3.2 2.1 3.7 5.8 5.2 4.7 4.0 3.1
##  [253] 5.2 5.7 4.8 4.7 5.5 3.9 4.6 6.0 4.0 4.3 3.7 5.4 4.4 5.3 4.9 4.2 3.4 2.6
##  [271] 5.2 5.1 2.9 4.6 4.3 3.0 4.7 4.6 4.0 5.5 5.5 4.5 3.9 4.7 4.7 4.5 6.0 4.6
##  [289] 4.5 2.0 4.0 4.2 4.8 4.4 3.7 3.2 4.7 5.8 4.4 4.1 5.0 4.3 4.5 5.6 3.8 3.1
##  [307] 3.6 6.1 6.2 4.7 4.8 3.9 6.0 3.5 5.3 4.5 2.9 4.1 3.0 6.0 4.1 6.5 3.3 4.5
##  [325] 3.2 3.6 3.1 4.1 4.2 3.8 5.2 5.5 3.4 5.2 5.4 4.9 5.1 2.6 4.4 5.6 5.0 3.3
##  [343] 4.5 3.8 5.1 5.3 5.0 4.1 5.3 3.3 3.6 4.4 5.2 4.8 3.0 4.5 4.6 4.8 5.7 4.4
##  [361] 6.7 5.3 6.3 3.6 4.0 3.2 4.6 2.6 4.7 2.8 4.7 4.2 2.6 5.3 5.7 4.4 4.6 2.8
##  [379] 2.7 3.6 4.5 4.3 4.3 3.8 5.1 4.9 3.5 7.1 4.5 4.1 5.2 6.0 4.5 4.0 5.8 6.5
##  [397] 4.1 5.4 3.2 6.0 3.6 3.8 5.7 3.3 4.5 4.8 4.8 5.4 3.7 1.7 3.6 4.8 4.5 3.8
##  [415] 4.3 5.8 3.6 6.2 4.9 4.3 5.1 4.0 5.1 4.4 5.1 3.6 4.4 3.4 5.9 4.6 4.5 5.0
##  [433] 4.2 4.0 4.6 2.9 2.7 5.4 2.9 5.4 4.4 6.1 4.5 5.0 5.1 2.5 4.6 2.5 3.6 5.7
##  [451] 5.0 4.8 5.3 4.5 4.0 3.9 5.4 5.5 3.5 4.5 3.6 4.0 4.0 5.6 6.0 4.8 5.6 4.0
##  [469] 3.3 3.4 4.5 4.3 3.7 5.7 3.5 5.6 4.0 4.1 4.3 4.4 6.2 3.4 3.9 4.8 4.7 3.4
##  [487] 4.3 4.5 4.1 4.9 3.9 2.5 4.3 4.5 4.4 4.9 4.2 5.1 5.2 4.5 4.9 3.5 3.6 3.5
##  [505] 4.9 4.5 5.0 3.7 5.4 4.0 4.8 4.2 4.2 2.9 4.6 5.7 4.4 3.2 5.7 3.3 4.6 4.8
##  [523] 4.7 4.4 5.1 3.9 3.3 4.1 3.4 4.3 4.2 4.0 2.8 4.1 5.2 4.5 3.6 5.0 3.5 4.0
##  [541] 5.6 4.0 4.0 5.1 3.9 4.1 3.6 4.0 5.8 4.9 3.9 4.8 5.2 5.7 5.3 4.4 3.1 3.9
##  [559] 4.3 4.6 3.3 5.0 5.5 5.7 4.5 4.2 4.7 4.9 4.5 3.6 5.6 3.8 4.4 3.2 4.6 4.0
##  [577] 3.8 4.7 4.7 5.0 5.5 4.4 5.3 4.8 4.5 3.3 4.4 5.9 4.3 4.6 3.0 4.3 4.6 4.2
##  [595] 5.1 4.4 4.7 4.3 5.6 2.7 2.2 3.8 5.1 3.6 3.7 4.9 4.1 4.5 4.0 6.3 5.2 3.9
##  [613] 4.7 5.2 5.5 5.0 5.1 3.0 3.7 4.1 5.8 5.2 6.1 4.4 4.7 5.6 2.6 4.6 7.8 4.1
##  [631] 5.2 6.2 5.8 3.0 5.0 4.7 4.2 5.1 3.2 5.3 4.2 5.0 2.8 4.3 4.7 4.2 4.0 2.6
##  [649] 4.4 3.6 3.7 3.4 4.2 3.9 4.2 4.1 6.4 4.5 5.0 4.5 5.1 2.9 4.9 2.8 3.7 6.2
##  [667] 6.4 4.8 3.8 5.1 5.3 4.9 5.4 4.4 4.4 4.6 4.0 4.9 4.9 3.8 4.3 3.7 3.0 3.5
##  [685] 4.7 3.4 5.5 4.2 4.9 3.5 4.4 4.2 3.5 5.0 3.4 4.3 4.8 5.0 3.4 6.1 5.5 3.2
##  [703] 4.6 3.7 4.8 4.9 3.9 4.7 5.2 4.0 5.6 3.2 4.3 3.3 4.5 3.2 4.9 5.5 3.5 4.7
##  [721] 5.6 4.0 4.5 5.8 5.1 3.6 4.0 3.9 5.0 4.8 4.6 6.4 4.0 6.1 3.7 3.7 5.0 2.9
##  [739] 3.3 4.3 4.8 2.9 3.9 5.2 3.6 5.3 4.6 3.6 4.5 4.9 5.9 6.5 4.8 4.9 4.2 5.1
##  [757] 5.3 2.7 3.3 4.4 4.2 7.0 5.2 4.7 5.0 3.8 3.6 3.7 4.2 3.8 3.6 3.4 4.7 4.4
##  [775] 3.8 5.8 3.6 5.0 3.5 4.7 4.3 3.9 4.1 3.9 3.8 3.4 4.2 6.3 4.0 4.9 6.3 3.9
##  [793] 4.9 4.3 5.2 3.8 5.0 4.5 4.4 4.4 5.6 3.0 4.8 4.5 3.6 4.4 4.5 4.3 3.6 4.9
##  [811] 4.0 3.9 3.8 5.3 3.7 3.7 5.0 4.9 4.3 5.1 3.9 4.4 6.4 4.3 3.5 4.4 5.7 4.1
##  [829] 3.2 2.3 4.5 3.5 4.2 5.3 5.3 5.5 4.5 5.3 5.1 6.4 4.2 5.1 3.0 4.0 3.3 4.8
##  [847] 4.1 5.4 4.6 5.1 5.8 5.1 4.0 4.5 5.1 3.9 5.2 5.1 4.8 4.3 5.1 4.2 2.9 3.9
##  [865] 5.5 4.5 4.6 4.4 4.7 4.6 4.3 4.3 5.9 4.3 3.8 5.5 3.5 4.9 4.7 4.8 2.8 4.0
##  [883] 4.6 4.9 4.2 2.0 3.7 3.3 4.0 4.6 5.4 5.0 3.0 5.6 3.7 6.0 4.5 3.7 4.2 3.5
##  [901] 4.4 5.4 5.3 4.7 4.9 4.1 4.7 5.7 5.7 4.0 6.5 5.6 4.4 4.9 4.1 5.8 4.2 5.4
##  [919] 5.1 4.2 4.3 3.4 4.2 4.6 3.5 3.9 3.4 3.9 4.1 5.5 2.6 5.0 4.2 3.2 4.8 4.0
##  [937] 4.5 4.9 4.3 4.3 4.2 4.5 6.1 5.8 4.7 3.9 3.9 5.9 3.0 5.0 4.1 5.3 5.1 4.7
##  [955] 5.6 4.5 3.9 4.7 3.9 3.1 4.3 6.1 3.5 3.9 3.4 4.5 5.9 3.8 2.6 3.5 4.3 2.9
##  [973] 4.2 5.0 3.0 4.5 5.0 3.7 4.1 4.9 4.1 5.5 5.4 4.2 5.8 3.3 6.7 6.4 2.0 3.4
##  [991] 4.8 4.2 4.8 2.2 4.7 5.7 3.9 4.5 4.1 5.6
## 
## $func.thetastar
## 72% 
##   5 
## 
## $jack.boot.val
##  [1] 5.5 5.4 5.4 5.1 5.1 4.9 4.9 4.7 4.6 4.5
## 
## $jack.boot.se
## [1] 0.9990495
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
## [1] 1.577085
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
##   3.813398   7.967904 
##  (1.636287) (3.654313)
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
## [1] 0.9648322 0.2978322 0.1600926 0.4677651 1.6222204 0.1543230
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
##    [1]  1.0342443114  0.9729321902 -0.5183619391  1.5390465548 -0.6774228415
##    [6]  0.4684911087  0.9154597032  0.9365352076  0.9649895487  1.2162238795
##   [11]  1.0396155619  2.1030144162  0.9234338536  1.2861792566  0.9329876533
##   [16]  1.0438023127  1.3563510966  1.2722481737  1.0832092589  0.6247653217
##   [21]  1.7612247286  1.2951881458  1.4131706931  1.4188741663 -0.1997077591
##   [26]  0.6988608229  1.9121781521  1.1925685525  1.8813775429  1.4832660162
##   [31]  0.3987393361 -0.1478871894  0.9447535496  0.8885090384  0.7767629829
##   [36]  2.2411262297  1.3777518609  1.2980230120  2.0360705944  0.6863343307
##   [41] -0.4683549274  0.5032760908  1.2837340845  1.6631016796  1.3967868132
##   [46]  2.1317954746  1.1928889177  1.1519552554  2.0840394101  1.7182416623
##   [51] -0.5845505110  1.0152445424  1.1406440418  1.6701616125 -0.1006233842
##   [56]  1.1358167744  2.1541516567  0.9163312434  1.1821753116  0.5850277582
##   [61]  0.5575148592  1.0141559427  0.3769456409  0.7882821723  0.6788621278
##   [66]  0.3864057078  1.8963259553  1.7715599568  0.6606841956  0.4327808900
##   [71]  0.1186027501  1.1206165360  0.9286663824  0.6976899998  1.4130907729
##   [76] -0.5251544058  1.5644252047  0.0591453486  1.1659951574  1.1938997241
##   [81]  1.0376307681  0.9819587338  0.8298970997 -0.1495245557  0.9370678505
##   [86]  1.1804086949  1.2027870205  2.1949923690  0.4812868744  1.2549139351
##   [91]  0.7672797062  1.5871049971  1.1475004706  0.6581356781  2.0964748289
##   [96]  1.2906359988  1.8595653607  0.4385554608  0.9048473993  1.0557680593
##  [101]  0.5624576645  1.3415156524  0.8740026101  1.0483912621 -0.2341399660
##  [106] -0.7366114619  0.4894435012  0.3426878402  1.0559287176  1.7844298445
##  [111]  0.3962831751  1.0036849809  0.3429952547  1.9993361780  1.4334940819
##  [116]  1.7569228692  1.2710193604  0.3094523918  1.0430495210  0.5072858069
##  [121]  1.1816254886  2.2021146572  1.4601576515 -0.6784440334  0.6977282522
##  [126]  1.3865791891  1.6144694545  0.9601193977  2.0515751266 -0.6170348141
##  [131]  0.6860392468  0.7529737277  0.6327288306  2.2057202589  0.3235355000
##  [136]  0.3464183953  1.0664528752  0.9374996387  2.0450480341  0.0816670075
##  [141]  0.7931660009  1.6197353086  1.4334116348  0.6897302226 -0.3181886321
##  [146]  0.2079000936  1.7816728667  0.9692142208  0.9555726268  2.0856187703
##  [151]  1.5007301286  0.8262065551 -0.2559098868  0.6846378561  1.0411216625
##  [156]  2.1543941875  1.7132701149  1.0947069256  1.3751842095  2.0960000171
##  [161]  0.3693077209  1.4137854448  0.6882133482  0.2374845132  1.0353117956
##  [166]  1.7633714352  1.2312031373  2.3486647342  1.2137064565  0.7076104079
##  [171]  0.5853912861  1.1571470143  1.6446631312  0.7175536984  0.1764485883
##  [176]  0.8568098256  1.3648151196  0.6615999157  1.2547917503  1.6425302146
##  [181]  1.0938377574  1.8279213350  0.7416088218  0.9717847497  0.6853237154
##  [186]  2.1537886300  1.0533519498  2.3173348314  0.2667197725  1.2741824362
##  [191]  1.4358651613  0.8050293015 -0.2517812629  0.2124283120  1.2913268778
##  [196]  1.0154652455  0.3287139148  1.4306455230  1.3314971502  1.9651086224
##  [201]  1.0525258412  0.3068957936  0.8327507283  1.2496630665  1.1192138829
##  [206] -0.2678704920  0.5907462638  2.1071178034  1.7071117298  0.5383229966
##  [211]  2.2468863499  1.9578938147  0.4783179292  0.8586939575  1.1136313495
##  [216]  1.5781638422  1.6048974240  1.5554012676 -0.3415995953  1.4993369568
##  [221]  1.7241528125  0.5426938950  0.3944097067  1.5469117652  0.0770206484
##  [226]  2.0964748289  0.5112761614  0.6739847074  0.0935741243  1.4857220009
##  [231]  0.7357772613  1.1160823238  1.5677616513  1.4752231798  0.1136769106
##  [236]  0.2261536683  1.0398528670 -0.0643073581  1.2746687942  1.0782640637
##  [241]  1.6305829431  0.8088467401  0.4137468502  1.3119135479  0.2044783228
##  [246]  0.4899134338  0.9348596075  0.7625633395  1.6905803101  0.0947574217
##  [251]  1.2978166528  1.1582950333  0.9631525833  2.2163830153  1.5392655730
##  [256] -1.1320377070  1.0786083546  0.9779577723 -0.3954092929  1.3479560534
##  [261]  0.4043854730  0.2954320862  0.7322881981  1.2308129041 -0.0348283681
##  [266]  1.1204590925  1.5441915809  1.1505783916  2.3292280778  1.7067607792
##  [271]  1.2925797881  0.3878448605  1.2247286051  0.3985359227  0.2046065252
##  [276]  0.7882821723  1.6730331769  1.4909858035  0.6004096341  0.1976561817
##  [281]  1.1772624652  1.1679214412 -0.5469196663  1.5312705275  0.4038596087
##  [286]  1.0402274912  1.2156358662  0.3786062291  1.1980857473  1.3691361213
##  [291]  0.6387566425  1.4160342646  1.1667383976  0.4495826137  1.8681307806
##  [296]  0.5774113402  2.2765868595  1.1349413658  0.4129207068  1.3221633644
##  [301] -0.0750482183  1.4486866128 -0.0353623989  1.0843812398  0.7603300152
##  [306]  2.4146031182  0.4788426828  2.0541860205  1.0262654804  1.7603601788
##  [311] -0.0405273731  0.5428379987  0.8431834866  0.7139424679  2.0315090324
##  [316]  1.2930203511  1.5160473757  2.2838217535  2.2278044387  0.1976393345
##  [321]  1.6979023442  1.4694322574  0.2631658083  1.4577633316  0.6582056274
##  [326]  0.5993823780  1.2126522902  0.9553128827  1.8941016444  1.4054890335
##  [331] -0.2646100804  1.1352662307  1.2515983730 -0.1570912393  1.1108618527
##  [336]  1.1861510468  0.7131798105  1.1041520738  0.6434589041  0.8345777129
##  [341]  0.7607696580  0.3847817030  1.7045643808  2.1488349949  1.2123987192
##  [346]  0.1404163986  0.3116003690  0.6980831577  2.1678673525  0.0925825401
##  [351] -0.7450699297  1.3945598864  1.1459615525  2.4185647336  1.4015341955
##  [356]  1.9153541239  0.7315109849  1.2364426490  2.3035834202  1.1384785529
##  [361]  0.0826244923  0.3921075384  1.3123775974  1.1393401748  1.2835111454
##  [366]  0.7112789649 -0.1077221624  0.7457681493  0.7173174165  2.0943143883
##  [371] -1.0791401002  1.2423736193  0.8409102131  1.1134004430  0.7231754603
##  [376]  0.8061304367 -0.0868192614  0.5626132245  1.5328467052  1.2202460138
##  [381] -0.1749501101  1.6787279535  2.2761171388  0.6705931973  0.8722776247
##  [386]  0.7448853159  0.4029754849  1.4009125222  1.4015341955  0.9080311809
##  [391]  0.6231571867  0.5661881287  0.5555324367  0.6078153375  0.1513518677
##  [396] -0.3327777207  0.5595212881  1.2101654500  1.3152129229  1.6964843946
##  [401]  1.3010941568  1.8923027817  0.1124601531  1.6578921653  1.6271272897
##  [406]  2.0824662143  1.1476427083  1.0725635862  2.0652331616  0.8010699895
##  [411]  1.0140237324  1.3938352869  1.7285895389 -0.0687576989  2.2161112265
##  [416]  0.8627378496  1.1269227650  0.4523022130  1.2410785626  0.4008223890
##  [421]  1.5663080129  1.5688009946  2.0225182887  1.4104589274 -0.1406026659
##  [426]  1.0238673595  1.9308321589  0.4684141329  0.4299913171  1.2037100031
##  [431]  1.7015972448  1.3193349030  0.1855592855  1.9840245719  1.5384436018
##  [436]  0.7899439705  0.3749074847 -0.0522973317  1.7668485184  1.2258616321
##  [441]  1.6536194640  1.3090267109  0.6124279916  1.2663146135  0.9512702931
##  [446] -0.3603540484  1.3666017552  0.4701763162  0.6922091791  0.7897221164
##  [451]  0.9418851272  0.7457236823  1.7728324728  0.5739143109  0.6040566283
##  [456]  0.5738718505  1.0150378373  1.9686414938  0.6492805564  2.2455661402
##  [461]  0.6795550311  0.8053501729  0.1344031387  0.8915034439  0.8102413842
##  [466]  1.1967065496  0.7215135169  1.5122472587  0.3389909758  0.3966189095
##  [471]  0.8689036419  0.5477319791  0.7544786517  1.1607410532 -0.1236338393
##  [476]  0.8119620504  0.1121635916  0.3068957936  0.2687228035  1.0091180364
##  [481]  0.6582056274  0.9012004930  1.6256008062  1.5120603856  0.4659701875
##  [486]  0.9243814719  1.4220275577  2.0612725049  1.9077546203  0.9217351910
##  [491]  1.2790892339  1.2541621085  0.9118452745  1.0712183288 -0.6589486213
##  [496]  1.3448371188  0.9351362458  1.1629544395  1.6796224131  2.0570814923
##  [501]  0.0116977186  0.3853134841  0.7584483255  0.9972126819  1.2389442844
##  [506]  0.4855252890  1.2491863325 -0.3775088372  0.7695618880  1.5558090041
##  [511]  0.2336613240  0.9002909592  1.7647704389  1.2828335462  0.8117562915
##  [516]  0.3148294741  1.7265387664  0.8933915098  1.1677640730  1.3273581808
##  [521]  1.2693626092  1.9325408608  0.6346059247  1.5806552219  1.5176143572
##  [526]  1.2871198446  0.5532445492  1.7491544207  1.2595746309 -0.5482793866
##  [531]  1.8376551300  0.9860479890  1.4691040810  0.3689057719  1.3259200226
##  [536]  1.2299476854  1.6112110444  1.4830878535  0.6283672624  0.9755502304
##  [541]  0.9373833087  1.3035105870  1.7223069636  2.2264541605 -0.0756157625
##  [546]  1.3324586652  0.5148811243  2.1743576154  0.5734506626  0.8015223932
##  [551]  1.2008586804  0.2867458423  0.9447574273  0.9938105272  1.7075452576
##  [556]  1.7081366137  1.2608633987  1.3516814889  1.2266065471  1.5542417828
##  [561]  1.1666714207  0.8092772438 -0.5936398135  1.2603316060  0.7696396593
##  [566]  1.0369816927 -0.0823231675  1.1290206094  1.4641602185  0.5314147533
##  [571]  0.4379554344  1.3458467169  1.1674488510  0.5049556688  1.7485784675
##  [576] -1.2052751756  1.8162273491  1.0940350234  1.9252794854 -0.1216475066
##  [581]  0.5896125151  1.0434519936  0.7644031427  1.5390813357  0.1254094947
##  [586]  0.5335045416  1.0964487878  1.6264794410 -0.4611385977  2.1025679675
##  [591]  0.4497410888  0.5289561638 -0.0532685321  0.9581558424  1.9696027367
##  [596]  1.0026185543  1.3866790433  0.7420174928 -0.0291214634 -0.0873914599
##  [601]  0.2189498752  1.1162568872  1.3674566562  1.6322124195  2.0202786696
##  [606]  0.1759984529  2.1554025138  0.8671829160  0.7940419543 -0.1928887874
##  [611]  0.2478766949  1.1959449440  1.3327589787  2.1445052535  0.5264275146
##  [616]  2.0960000171  0.8272746972  0.2200796073  1.6821855037  0.6842333242
##  [621]  1.3782328639  0.6697627286  2.0669997530 -0.3001343004  0.7936690756
##  [626]  1.6469075808  1.0203565803  0.8491985362  1.7523412521  1.0523689685
##  [631]  1.7846138221  1.0092569348  1.8277067229  0.8055261650  1.2629563534
##  [636]  0.6183400363  1.5940600254  0.7458549210  0.9712558682  0.0893814336
##  [641]  1.7640711217  0.6180961685  1.1068809850  0.2032873173  1.7879495510
##  [646]  1.9636545480  1.5366764928  0.7206027816  2.0936010585  0.2662301573
##  [651]  0.7371602133  1.8130067515  0.9427613436  2.2148529429  0.1710767161
##  [656]  1.6453707995  2.3071858485  1.2126522902  0.1612619580  1.8481800381
##  [661]  1.1219576435  1.1337684970 -0.5064955619  0.1117645577  1.3304158020
##  [666]  1.0894811513  0.5150983418  0.7684360255  0.5907462638  0.6754102306
##  [671]  2.2963641715  0.6863558340  0.9394654599  0.9939381619  1.2741374128
##  [676]  0.8594690422  0.6069640553  0.5789657474  0.4357335728  1.7032667934
##  [681]  1.1313923126  1.5782921622  2.1025380894  1.0367116819  1.5332237951
##  [686]  1.9560523061  0.7036568531  0.9550905956  2.1103296403  1.1759231702
##  [691]  1.0043632908  2.2433079758  0.7682692055  1.3189162517 -0.0582530599
##  [696]  2.3159200349  0.7140615318  0.7025916434  0.9830368296  0.7617853376
##  [701]  1.5773977571  0.7928172879  1.6884490135  0.0809523553 -0.0677754698
##  [706]  2.0543655134  1.2151658704  0.2722215909  1.6019265993  1.7731270124
##  [711]  1.0951881046  1.8487437270  1.8139115026  0.9435769188  1.5590843473
##  [716]  1.6796653667  2.1251262912  0.6883648036  0.4310456670  2.1502082496
##  [721]  0.9160714609  1.3825697862  2.0851114297  1.6707007070  1.2196617019
##  [726]  0.4939339762  1.7439861819  1.5561298660  1.7890646317  0.2379150516
##  [731]  0.4803087279  0.5093892404  0.7534260078  0.4690071409  0.1939741615
##  [736]  0.4427896679  1.3123755517  1.7211273287  0.8228528818  0.9231646332
##  [741]  1.5877611621  2.1373747245  0.0584453579  0.5862654587  0.5478353574
##  [746]  1.2476827404  1.2641109234  0.0844055089  2.0765128750 -0.1870128005
##  [751]  0.8178691816  0.0938647481 -0.3272606903  1.3088144234  0.0968742801
##  [756]  1.1559257197  0.7231754603  1.7198079094  1.1276702273  1.3836932196
##  [761]  0.4161317098  0.2452618748  0.8689503057  1.4292441517  1.3676920627
##  [766]  0.6181744724  0.0404508560  2.1188752179  1.9428717229  2.0972570548
##  [771]  1.7012422833  1.4673745946  0.9384859562  1.7695755414  1.0316015941
##  [776]  1.2794097147  0.8748572717  0.4062489939 -0.1821959446  1.7600999918
##  [781]  0.5840222751  0.9028467603  2.0610190657  0.9065601856  1.6762161808
##  [786]  0.6582056274  0.5378913416  0.8093075301 -0.6439839195  1.3008927060
##  [791] -0.4436469433 -0.1066551222  1.6216126176  2.0069234870  1.7211273287
##  [796] -0.0927070556  1.5585315245  1.4338591605  1.4518113070  0.7010529851
##  [801]  0.6700707259  0.6316973727  0.5952346329 -0.5017813037  1.4863341874
##  [806]  1.6843549491  2.5280847004  0.7548544998  1.5850920192  0.8008658193
##  [811] -1.5460456502  0.9677229161  0.9305122164  1.1629795485  0.9724221715
##  [816]  1.6315697728  1.2736025672 -0.1465314910  0.6140222623 -0.2709166377
##  [821]  1.3147046848  1.6049622563  1.6161243951  0.9660614190  1.8621528119
##  [826]  1.4230972173  1.3116419488  0.6219044485  1.7491544207  2.1051129331
##  [831]  2.1412445424  0.7556050950  0.2323931348  1.3750343872  0.7159412788
##  [836]  0.9891650442  0.7833033865  0.4770648220  0.0392195873  0.7753324120
##  [841]  0.3146907498  0.1309236033  1.1633601399  0.5249917813  0.4582836693
##  [846]  1.2910275957  1.4204236596  0.6540065941  0.6877002651  1.1947800255
##  [851]  1.1077105271  0.5145885741  1.3235874694  0.6402710583  0.8956831029
##  [856]  1.1793156640  0.4119406354  0.7717924922  1.3525581146  0.6370509415
##  [861]  0.2836457695  0.8323731871 -0.1321232844  1.5874858163  2.2115763595
##  [866]  0.6110128336  1.6506900885  0.8001201347  1.0721123058  2.1858366122
##  [871]  2.1164913108  1.3697759863  1.9969589374  0.9426011827  0.8038050246
##  [876] -0.2501890579  0.5552160305  1.3000738078  0.6246568101 -0.0378713105
##  [881] -0.2095174297  1.5081497900  1.0553870943  1.0849262534  0.2865278977
##  [886]  0.0299701295  1.0634929070  0.2520218161  0.7574270507  1.5406549600
##  [891]  1.7539069065  2.2031957672  1.0678782832  0.7532519811  0.6995364789
##  [896]  1.1166228940  2.1082033105  1.7594575454  1.8704995983  0.7065287930
##  [901] -0.7585954168  1.7034594075  1.1906058618  1.3100438094  0.4500354630
##  [906]  1.0257920795  1.4622958216 -0.4870525217  1.5160178750  0.7902429816
##  [911]  0.4010253500  1.1511333000  1.1335560543  0.7365623290  1.1969848484
##  [916]  0.7935884119  0.7220899606  1.8237408218  1.5995962098  0.5357886894
##  [921]  1.6621936036 -0.3738457028  1.0468428609  0.7291856152  1.0901015751
##  [926] -0.0896158427  1.2572546356  1.0657782570  1.1437242863  0.0002437956
##  [931]  0.7022646684  2.3126404833  1.3199590102  0.1798088010  1.0805637059
##  [936] -0.0692780069  1.3785430104  1.9999037124  1.1617929438  0.8236372507
##  [941]  1.3769988552  0.2101283100  1.0689354441  0.9396390769  0.9128223335
##  [946]  1.1069817073  0.9342779075  1.3376314559  1.8244573585  1.2407074069
##  [951]  0.0431105497  1.3499105716 -0.8735810666  1.0297383637  0.8054427258
##  [956]  1.1017250197  1.6571562289  0.5616665286  1.2643712247  2.1562508592
##  [961]  1.8865683532  1.2218880895  0.7220095909 -0.5245682950  1.2328222298
##  [966]  1.4869095936  0.1450156474 -0.2300008486  1.2118079748  1.3458467169
##  [971]  0.3944613215  1.1263032074  1.3360771061  0.2686814423 -0.5356645536
##  [976]  1.1958086073  0.0009124546  0.9555726268  0.2530243817  1.4578669335
##  [981]  0.8413869817 -0.0923465347  1.8162273491  1.6735445488  1.5877784946
##  [986]  2.0131800262  0.7106626367  0.7095025165  0.4016732715  0.7183568887
##  [991]  1.5096360424  2.2254251274  1.2561311422 -0.2387177472  0.1756588140
##  [996]  1.2250914277  1.8643347649  1.3476688740  1.0178851244  1.3695551431
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
##   0.47859478   0.27939942 
##  (0.08835385) (0.06247298)
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
## [1] -0.3356062 -0.2439235 -0.3330007 -0.1547865  0.3333690  0.3681591
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
## [1] -0.0029
```

```r
sd(results2$thetastar)  #the standard deviation of the theta stars is the SE of the statistic (in this case, the mean)
```

```
## [1] 0.8720792
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
## t1*      4.5 -0.06106106   0.8924731
```

One of the main advantages to the 'boot' package over the 'bootstrap' package is the nicer formatting of the output.

Going back to our original code, lets see how we could reproduce all of these numbers:


```r
table(sample(x,size=length(x),replace=T))
```

```
## 
## 1 2 3 5 6 7 8 9 
## 1 1 2 1 1 1 2 1
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
## [1] 0.0342
```

```r
se.boot
```

```
## [1] 0.8617616
```

Why do our numbers not agree exactly with those of the boot package? This is because our estimates of bias and standard error are just estimates, and they carry with them their own uncertainties. That is one of the reasons we might bother doing jackknife-after-bootstrap.

The 'boot' package has a LOT of functionality. If we have time, we will come back to some of these more complex functions later in the semester as we cover topics like regression and glm.

