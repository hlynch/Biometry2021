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
## 0 1 2 3 4 6 7 8 
## 2 1 1 1 1 2 1 1
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
## [1] -0.0039
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
## [1] 2.761818
```

```r
UL.boot
```

```
## [1] 6.230382
```

Method #2: Simply take the quantiles of the bootstrap statistics


```r
quantile(xmeans,c(0.025,0.975))
```

```
##   2.5%  97.5% 
## 2.8000 6.2025
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
##    [1] 4.5 3.9 3.9 3.3 4.2 4.9 4.3 4.5 4.5 4.5 4.5 6.6 4.0 4.6 4.2 5.5 4.8 5.2
##   [19] 4.9 6.1 4.1 3.9 3.9 5.5 6.7 3.8 4.7 4.2 5.5 4.4 5.0 2.1 5.6 2.4 3.9 3.9
##   [37] 3.5 4.1 3.7 5.5 5.6 5.1 3.1 4.8 3.5 5.2 5.3 4.7 3.2 4.9 5.2 5.5 5.1 3.7
##   [55] 5.1 4.3 5.2 4.1 4.3 5.4 4.7 3.7 4.8 5.0 2.1 3.9 6.2 5.7 3.5 4.5 4.9 4.7
##   [73] 4.1 3.7 2.8 4.1 5.1 4.2 4.9 6.4 5.5 2.6 5.5 3.6 4.9 3.4 4.7 3.8 4.5 6.1
##   [91] 5.9 6.3 4.2 5.9 4.2 3.9 3.4 5.6 4.6 5.5 5.0 4.2 5.6 3.5 4.9 4.5 3.8 4.6
##  [109] 5.2 5.3 3.9 3.1 4.6 4.1 2.5 4.3 5.5 4.4 4.1 3.8 5.3 3.9 5.7 5.4 4.1 4.4
##  [127] 4.6 4.1 2.5 3.4 3.9 4.8 4.4 4.5 4.2 3.8 4.9 4.2 3.8 3.4 4.1 3.9 4.7 4.5
##  [145] 4.6 3.4 5.0 6.2 5.5 5.5 3.9 4.2 5.2 4.5 4.1 4.1 2.9 4.2 4.8 4.1 6.0 4.2
##  [163] 4.7 4.0 3.6 3.3 5.0 5.5 2.8 2.9 4.8 3.6 5.7 5.4 4.2 3.1 4.6 5.4 4.7 4.5
##  [181] 5.2 4.3 4.9 2.9 4.0 4.8 4.6 3.7 5.6 4.5 4.5 5.2 4.0 4.1 5.1 4.7 6.2 5.9
##  [199] 4.9 5.6 3.3 3.2 3.1 5.1 4.8 5.0 4.9 4.0 5.5 5.2 3.6 4.4 4.6 3.3 2.5 3.7
##  [217] 2.8 4.5 4.0 5.1 4.1 4.1 6.0 3.6 3.5 4.3 4.1 4.0 3.6 5.0 5.6 5.4 4.3 4.5
##  [235] 4.5 5.1 5.6 3.7 5.0 3.0 4.4 4.9 6.2 4.5 3.8 4.0 3.8 4.1 2.9 4.8 5.0 4.4
##  [253] 3.9 4.6 4.5 3.2 4.6 3.6 7.3 4.2 5.8 3.4 5.3 5.4 4.4 4.2 3.9 3.7 4.3 5.2
##  [271] 4.4 3.8 3.1 4.6 4.6 4.9 4.6 2.6 5.4 5.5 6.8 3.6 5.4 3.4 4.4 3.4 4.0 5.2
##  [289] 4.5 4.8 2.7 6.4 3.8 3.8 3.2 4.5 5.4 4.4 3.5 4.9 3.4 2.9 4.3 5.4 3.1 3.8
##  [307] 4.1 5.2 4.7 3.8 3.6 4.2 5.0 5.3 3.6 4.9 3.2 3.4 4.7 4.3 3.4 2.5 3.3 3.0
##  [325] 5.2 4.8 3.4 5.0 6.0 5.5 2.6 5.0 3.3 5.2 4.5 5.0 4.1 5.5 4.4 5.1 5.1 6.8
##  [343] 5.8 3.5 4.4 2.5 3.9 4.9 4.3 4.9 5.5 4.9 5.5 5.0 4.6 5.7 3.2 5.2 5.0 3.4
##  [361] 4.4 4.0 4.1 4.5 5.5 5.0 3.9 2.6 5.2 4.0 3.9 5.3 3.7 4.6 3.4 4.9 5.2 4.1
##  [379] 3.7 6.3 5.0 4.3 4.2 4.7 4.2 4.7 5.5 5.5 4.1 3.1 3.9 4.5 3.4 4.8 4.8 5.4
##  [397] 4.2 2.2 3.3 4.2 4.4 6.0 4.9 3.8 4.8 4.8 4.9 4.2 6.4 3.1 4.9 6.1 3.9 6.1
##  [415] 3.5 5.2 4.3 4.7 3.9 4.2 4.0 4.7 4.2 4.8 3.7 2.6 4.9 5.4 4.0 5.7 3.0 3.8
##  [433] 3.3 3.2 2.8 4.8 5.3 5.7 4.0 5.2 3.6 4.6 3.6 5.5 4.0 3.5 6.2 4.3 5.1 3.7
##  [451] 4.8 3.6 6.4 3.6 5.4 3.6 5.8 5.0 3.2 4.7 5.5 4.3 4.9 4.6 4.6 3.7 4.7 3.6
##  [469] 4.5 4.4 4.2 3.2 3.6 5.4 3.8 4.1 3.9 4.2 4.3 5.7 3.7 3.6 5.1 4.7 6.5 5.7
##  [487] 4.9 3.5 3.6 3.7 4.3 5.3 4.6 4.2 4.0 4.2 4.0 4.4 3.4 4.2 4.9 4.2 3.1 4.9
##  [505] 4.0 6.2 5.3 3.9 5.7 4.7 3.3 4.6 3.5 5.4 4.9 4.1 4.8 3.5 3.6 4.4 4.1 7.1
##  [523] 4.1 5.0 3.5 5.1 5.9 4.8 4.6 3.8 4.3 3.4 4.9 4.3 3.4 4.6 4.2 5.7 4.2 4.9
##  [541] 4.2 3.9 4.5 4.0 4.7 5.5 4.0 4.9 4.6 4.1 4.7 3.3 3.8 5.1 4.6 4.0 4.2 5.5
##  [559] 4.3 4.8 3.5 3.4 3.6 4.0 3.5 5.9 6.1 4.3 4.9 5.3 4.9 5.9 4.2 4.3 5.2 4.0
##  [577] 2.9 4.0 3.5 4.8 4.7 4.8 3.4 3.2 4.8 2.9 4.1 4.9 4.7 4.9 3.7 4.6 5.1 5.5
##  [595] 4.5 4.7 5.3 4.3 3.6 3.6 5.3 3.6 5.2 3.2 5.0 2.9 3.6 3.7 6.3 3.7 4.7 5.4
##  [613] 4.1 4.2 3.9 5.1 5.2 4.8 5.4 7.0 3.3 4.2 4.7 5.1 3.8 2.7 4.1 5.1 6.0 2.5
##  [631] 3.5 4.4 4.4 5.1 3.4 4.4 4.7 4.3 4.3 4.9 5.6 4.6 4.2 5.1 4.8 3.2 3.8 6.3
##  [649] 4.0 4.5 3.9 4.5 4.1 3.8 2.5 5.6 2.6 4.7 5.6 5.3 4.6 4.2 2.9 4.6 4.7 4.9
##  [667] 4.1 4.0 3.4 3.2 5.0 4.0 4.1 5.3 5.0 5.0 3.2 5.1 6.4 5.5 3.2 4.2 4.0 5.0
##  [685] 4.7 5.9 5.0 4.6 3.3 4.1 4.9 4.7 5.3 5.4 5.8 5.3 2.1 4.8 4.6 4.4 4.9 5.3
##  [703] 4.5 4.4 6.4 4.4 5.6 4.0 4.8 5.1 4.4 3.8 4.6 5.6 5.2 5.0 4.7 4.9 5.3 4.1
##  [721] 6.7 3.8 4.3 4.0 3.8 3.1 4.4 4.9 5.3 4.7 5.7 5.7 6.2 5.2 4.9 5.4 5.2 4.2
##  [739] 4.4 6.2 3.3 4.6 6.0 4.6 5.2 4.2 3.6 5.3 5.0 4.3 3.9 4.9 4.8 3.8 6.1 4.3
##  [757] 4.7 3.0 3.7 4.2 3.9 3.7 5.2 4.9 3.4 4.7 5.3 4.1 2.2 3.7 2.3 3.5 4.9 3.1
##  [775] 4.7 4.0 4.6 2.8 5.5 5.2 5.6 4.3 4.7 4.2 4.7 5.4 4.6 5.1 4.3 3.2 2.8 5.6
##  [793] 5.2 3.3 6.4 4.5 5.0 3.9 5.6 4.6 4.9 5.5 4.7 4.3 4.5 4.0 3.6 3.9 4.9 6.5
##  [811] 3.4 4.9 4.1 5.1 4.3 5.0 2.9 3.1 3.8 5.6 4.5 4.5 4.5 4.3 5.0 2.3 2.7 4.3
##  [829] 5.3 2.8 2.7 5.3 4.7 4.5 5.8 4.2 4.7 3.6 3.1 5.3 3.2 3.2 6.0 3.0 4.6 5.3
##  [847] 4.9 4.7 3.9 4.3 5.7 4.9 2.9 3.8 4.5 4.6 5.7 3.1 4.5 4.7 4.4 4.2 2.4 3.5
##  [865] 5.5 4.8 5.5 4.6 4.2 3.0 4.4 3.8 6.1 3.5 5.0 4.5 4.3 5.8 5.8 4.0 4.5 5.8
##  [883] 4.7 5.4 4.2 5.0 3.4 4.3 6.5 4.1 4.4 5.1 3.9 3.5 4.1 5.4 6.3 5.8 2.2 5.8
##  [901] 5.1 4.0 3.1 5.0 4.8 4.3 4.9 3.7 2.6 4.3 4.8 4.2 5.8 5.1 4.5 5.0 5.0 4.3
##  [919] 4.1 4.2 4.2 4.6 5.0 4.9 3.7 4.6 2.5 2.5 4.4 4.5 4.9 4.8 5.1 5.3 4.6 4.8
##  [937] 4.1 5.3 5.2 5.4 4.0 3.4 4.5 4.5 4.4 5.1 5.3 3.3 4.2 4.4 4.2 6.0 6.1 6.4
##  [955] 4.8 4.7 3.6 3.7 3.4 5.0 6.4 7.5 5.0 4.5 3.6 4.6 4.5 3.7 5.2 5.1 4.8 4.8
##  [973] 4.6 3.7 4.1 3.2 3.4 3.3 5.7 5.0 5.1 4.3 4.4 5.3 5.4 5.3 4.6 3.8 5.0 4.5
##  [991] 5.1 4.2 3.3 4.8 4.3 5.3 5.6 5.8 3.5 3.8
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
##   2.6   6.3
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
##    [1] 4.2 5.5 4.2 3.6 4.0 4.4 3.2 5.1 3.8 5.0 3.2 3.3 3.6 5.0 5.3 4.6 5.4 3.8
##   [19] 3.9 6.0 5.4 5.4 5.3 4.0 4.8 5.3 2.2 4.8 3.6 5.3 5.7 6.9 4.4 5.8 4.0 4.8
##   [37] 4.8 3.8 2.8 5.6 3.5 3.6 1.7 3.5 5.1 2.8 2.5 4.0 4.2 5.4 4.5 4.9 6.3 3.4
##   [55] 5.4 4.3 4.1 3.4 4.4 3.6 4.1 5.0 4.7 5.2 5.1 4.7 4.2 4.4 4.5 5.0 3.7 4.9
##   [73] 5.0 4.3 3.6 4.6 5.2 4.8 2.8 4.7 3.9 4.9 4.1 4.1 4.5 4.5 4.7 5.0 2.7 3.6
##   [91] 4.1 3.9 4.1 5.5 5.2 4.1 3.5 5.5 5.3 5.6 3.8 4.2 4.8 4.2 5.1 4.4 4.4 4.1
##  [109] 1.7 4.3 6.0 4.4 4.4 5.0 4.9 3.0 3.7 3.4 3.9 4.1 4.4 3.2 4.6 3.7 5.4 4.4
##  [127] 3.7 5.6 3.8 4.5 4.8 4.0 4.0 3.2 3.8 3.3 4.8 4.4 4.5 4.7 4.1 6.0 3.8 5.1
##  [145] 4.7 4.8 4.0 4.4 4.1 3.3 4.4 5.1 5.1 3.2 4.5 5.7 5.9 6.1 5.0 4.7 5.3 4.8
##  [163] 4.1 3.7 4.4 3.7 4.6 3.7 4.7 4.8 4.7 3.8 6.6 3.2 5.0 3.9 4.3 4.5 5.5 5.1
##  [181] 4.4 4.8 4.2 3.4 3.8 4.4 3.8 5.1 6.5 3.7 4.3 3.6 6.6 5.4 4.9 3.3 3.9 4.5
##  [199] 3.7 4.3 4.8 3.2 4.3 5.2 4.0 3.1 4.0 3.9 4.7 4.8 4.2 4.7 5.9 3.3 3.6 5.2
##  [217] 5.0 5.6 4.2 5.8 3.5 4.5 4.8 4.1 5.1 4.5 3.4 5.3 3.4 4.4 6.5 4.4 5.3 5.4
##  [235] 2.7 3.7 5.1 5.1 3.4 5.3 6.4 4.4 5.0 4.9 4.7 4.9 3.8 5.0 3.3 4.4 5.6 4.6
##  [253] 5.5 4.7 3.7 5.5 4.4 4.1 3.2 4.3 5.2 6.6 5.6 2.5 3.8 5.0 2.4 5.2 5.8 4.4
##  [271] 4.5 4.3 6.6 4.8 3.6 3.0 5.4 4.8 4.3 4.7 3.7 3.7 3.5 5.1 5.5 5.8 6.3 4.3
##  [289] 4.1 6.5 4.2 5.4 4.2 4.6 3.0 3.9 4.1 3.4 3.2 4.9 3.1 5.2 3.7 5.4 4.0 4.3
##  [307] 4.9 6.1 5.0 4.5 3.8 4.0 5.2 4.8 3.8 3.3 3.3 3.4 4.5 4.5 4.9 5.1 4.2 4.2
##  [325] 4.0 3.8 4.4 3.8 2.5 5.9 3.5 3.5 4.7 2.8 4.2 4.6 4.6 3.6 4.4 4.9 2.5 4.8
##  [343] 3.6 5.2 6.3 4.4 3.9 4.4 5.0 5.5 4.7 5.7 3.7 3.7 3.3 6.4 5.1 5.3 3.6 4.0
##  [361] 3.4 5.7 4.8 5.6 3.8 5.0 4.4 4.6 3.9 4.5 3.7 4.9 4.9 6.3 4.0 4.0 4.2 4.9
##  [379] 4.5 4.1 2.7 4.3 4.8 2.8 5.3 4.3 3.8 2.7 4.6 5.9 3.9 4.1 5.3 4.6 3.5 3.2
##  [397] 2.7 3.2 5.4 5.2 4.5 4.1 3.7 4.3 3.4 4.9 3.4 5.0 5.5 4.6 5.5 5.2 4.7 5.1
##  [415] 4.4 4.6 5.9 4.1 5.9 5.4 3.9 5.9 2.7 3.7 2.8 5.3 4.0 4.1 4.2 3.8 3.7 6.0
##  [433] 3.0 4.4 5.1 5.0 2.4 5.7 2.2 4.1 3.7 4.1 3.8 5.4 4.8 5.0 4.2 4.8 5.3 3.8
##  [451] 5.3 4.3 3.9 3.0 4.8 4.8 5.2 4.3 2.8 4.7 5.5 4.0 5.9 5.8 5.1 4.5 5.6 4.0
##  [469] 2.2 5.0 4.9 5.0 4.6 4.5 5.1 3.2 3.1 3.3 2.7 5.6 4.5 6.3 4.3 5.4 5.2 6.3
##  [487] 5.5 4.8 3.8 3.4 5.9 5.3 4.9 3.9 4.6 4.8 4.4 4.1 6.1 4.0 4.9 4.0 5.0 4.1
##  [505] 4.6 4.6 4.6 4.8 5.4 3.6 3.5 3.8 3.5 5.7 3.9 5.9 4.0 4.6 3.7 3.1 6.2 5.5
##  [523] 3.4 4.6 4.8 5.2 3.3 4.8 5.0 4.7 2.8 5.1 4.1 5.3 4.4 4.3 2.5 4.6 5.0 4.5
##  [541] 6.6 4.4 5.2 4.6 5.5 5.0 4.8 5.0 2.6 3.6 5.0 5.4 4.3 5.7 4.7 3.5 5.0 4.2
##  [559] 5.0 3.6 4.6 4.1 5.4 2.9 5.0 4.9 4.1 3.7 5.1 4.8 5.2 5.3 6.6 4.7 2.5 5.9
##  [577] 4.1 4.6 5.0 5.3 4.4 4.2 5.6 4.1 4.7 3.8 4.7 4.9 3.3 4.5 5.5 4.5 5.6 6.1
##  [595] 4.3 4.5 4.8 3.4 4.1 6.1 4.2 5.4 4.1 3.7 5.9 5.0 4.4 3.7 5.0 4.1 5.2 2.4
##  [613] 4.1 4.2 5.6 3.0 5.1 4.9 3.2 6.0 5.8 3.8 6.4 5.0 4.4 2.8 3.2 4.5 5.8 6.2
##  [631] 5.4 5.7 3.0 3.3 5.7 4.3 3.2 3.7 5.6 6.0 5.0 4.0 5.0 5.0 4.2 3.7 5.5 5.7
##  [649] 5.2 5.0 6.0 3.3 2.8 4.6 4.1 5.0 4.2 5.4 4.3 4.6 3.1 3.6 6.2 4.9 4.0 3.9
##  [667] 4.6 5.9 5.1 3.9 3.4 3.7 3.1 5.2 4.1 4.2 4.4 2.6 5.1 3.2 5.5 4.1 5.6 5.1
##  [685] 5.5 2.8 4.6 5.2 4.2 4.6 3.9 3.1 5.2 5.5 3.4 3.9 4.7 4.0 3.8 4.6 4.4 4.8
##  [703] 5.4 3.9 5.7 3.4 3.2 4.5 2.4 3.2 5.6 4.1 2.2 2.9 4.1 4.8 4.3 3.8 3.9 5.0
##  [721] 3.9 4.0 4.9 5.3 6.0 3.3 3.2 4.1 5.3 4.9 5.4 4.7 5.7 3.9 4.1 4.0 4.6 3.7
##  [739] 4.3 3.9 5.5 3.2 4.1 3.9 5.8 4.2 3.9 4.9 4.8 4.2 5.0 4.4 3.2 3.3 3.8 3.0
##  [757] 5.7 3.3 4.9 4.5 3.9 4.5 4.6 4.8 3.3 3.6 4.1 4.8 3.5 3.6 5.6 5.2 4.6 3.5
##  [775] 4.1 4.2 4.7 5.0 4.0 4.0 4.6 5.2 3.3 4.5 4.3 5.0 4.2 6.2 5.2 5.1 3.3 3.9
##  [793] 3.4 5.9 4.7 5.8 4.4 4.0 4.5 4.1 5.4 5.0 4.0 4.5 3.6 4.4 4.1 4.5 5.4 4.8
##  [811] 4.8 4.2 4.4 3.5 4.9 4.3 3.9 4.6 5.3 4.4 3.4 4.5 4.7 5.3 5.1 3.9 4.7 5.0
##  [829] 5.0 4.1 3.1 4.4 3.3 4.0 5.1 4.6 3.5 4.8 3.9 3.9 4.0 3.6 4.8 3.1 6.5 2.8
##  [847] 4.0 2.9 5.9 3.4 4.5 4.4 3.6 4.1 2.9 3.2 4.9 5.8 3.8 4.4 2.6 4.7 4.2 5.0
##  [865] 5.6 4.9 3.3 4.7 3.9 4.0 4.1 5.6 4.9 5.8 4.3 5.6 3.8 3.4 4.8 4.1 5.9 5.1
##  [883] 5.8 4.6 3.6 4.0 5.1 3.6 4.1 5.0 3.5 4.6 4.6 4.3 3.5 5.4 4.8 3.9 4.3 2.0
##  [901] 4.8 3.4 5.7 4.4 5.5 2.7 5.1 5.2 4.7 5.1 5.6 4.7 4.6 3.7 4.3 4.4 3.9 4.7
##  [919] 4.4 3.5 3.8 5.2 3.3 5.2 5.0 3.8 3.6 4.5 4.2 5.8 3.8 4.7 5.2 5.6 3.2 4.1
##  [937] 4.7 4.4 3.5 5.8 5.2 4.0 4.0 2.9 5.1 4.0 5.2 3.1 2.5 3.8 6.8 3.3 4.6 4.9
##  [955] 6.0 4.4 6.4 4.1 6.2 4.8 6.3 4.4 5.3 3.8 4.1 3.5 5.2 4.6 4.3 3.8 3.3 5.8
##  [973] 5.9 4.9 5.0 3.4 4.7 5.1 4.7 3.8 4.7 4.3 5.5 5.5 5.2 6.6 4.0 5.3 4.4 4.4
##  [991] 4.6 4.1 4.6 5.5 4.2 4.4 3.6 3.2 4.8 5.7
## 
## $func.thetastar
## [1] -0.0425
## 
## $jack.boot.val
##  [1]  0.417441860  0.290988372  0.212534060  0.152760736 -0.008108108
##  [6] -0.093993994 -0.203133903 -0.238418079 -0.410000000 -0.527200000
## 
## $jack.boot.se
## [1] 0.8817704
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
##    [1] 4.7 5.8 5.6 5.0 5.1 4.8 6.0 5.0 5.3 5.0 5.7 3.4 5.1 5.9 2.9 5.3 5.6 3.9
##   [19] 4.9 4.3 5.3 4.3 4.0 5.6 4.4 4.7 5.8 5.3 4.1 3.9 3.8 3.7 3.0 3.0 4.3 3.7
##   [37] 5.1 3.5 3.8 3.1 4.0 5.5 3.9 5.1 3.1 3.8 4.8 3.8 3.8 2.9 2.5 4.7 2.9 3.5
##   [55] 4.1 4.7 3.8 3.4 4.3 4.8 3.6 6.1 4.6 4.3 6.7 5.2 3.5 5.8 5.1 5.4 4.6 4.3
##   [73] 3.8 4.0 5.0 4.9 4.6 4.2 5.4 5.7 5.0 4.3 4.4 4.5 4.0 5.4 3.5 3.1 5.0 4.1
##   [91] 5.4 5.8 4.2 4.1 4.7 4.1 3.2 4.3 5.2 6.0 4.8 4.5 5.8 4.7 4.4 3.6 2.4 4.3
##  [109] 4.5 3.8 5.5 3.6 4.8 4.1 3.7 5.3 3.3 4.6 4.2 5.0 5.6 5.2 5.9 5.5 3.3 4.8
##  [127] 3.7 4.1 3.1 5.3 4.1 4.2 3.3 5.5 3.5 3.8 4.6 4.0 5.8 4.5 6.3 4.8 3.6 4.0
##  [145] 4.7 4.7 4.4 4.5 4.7 6.0 5.8 3.9 3.9 5.4 3.6 4.1 4.5 4.1 3.6 5.5 6.4 4.5
##  [163] 3.0 3.5 6.9 5.9 3.2 4.9 5.4 5.1 5.9 4.6 4.2 4.2 4.5 3.4 4.0 4.5 3.4 4.1
##  [181] 4.3 4.2 2.9 3.6 3.2 5.5 5.8 7.2 5.1 5.9 4.0 4.7 4.7 4.1 2.7 6.2 4.8 4.9
##  [199] 3.8 3.6 5.3 4.2 4.3 4.8 5.5 4.0 4.9 4.2 2.8 5.0 4.5 4.5 5.1 5.9 4.5 3.8
##  [217] 3.3 4.7 4.2 5.5 3.8 5.3 4.2 4.3 4.4 3.8 4.9 5.1 5.1 4.3 5.6 4.6 5.9 4.0
##  [235] 5.5 6.4 5.1 2.8 4.3 4.7 5.1 5.4 4.3 5.6 4.2 4.3 3.7 5.1 3.4 4.2 4.1 3.8
##  [253] 5.3 2.8 3.5 4.9 4.6 5.0 3.1 3.8 5.2 3.4 3.9 3.8 5.5 5.2 4.5 4.1 4.6 5.0
##  [271] 5.0 4.3 5.2 3.9 5.2 3.7 4.3 4.6 4.5 5.0 3.4 4.1 4.5 5.5 5.3 4.0 3.8 5.4
##  [289] 4.5 4.7 3.2 5.1 3.8 4.8 4.9 3.0 3.7 4.3 5.4 5.3 5.8 4.6 5.1 3.9 4.7 4.5
##  [307] 4.0 5.3 3.3 5.7 5.3 5.0 4.4 3.9 4.5 4.4 4.6 4.2 4.3 4.9 5.5 5.4 3.2 3.8
##  [325] 3.9 2.8 5.1 7.9 3.0 4.3 5.7 4.8 6.2 4.6 6.0 4.1 3.1 2.8 4.2 4.9 5.2 2.9
##  [343] 4.2 5.0 4.3 5.8 5.6 4.5 4.9 3.2 4.7 5.1 4.9 4.5 5.3 3.9 3.9 5.3 4.9 3.9
##  [361] 5.3 4.9 3.4 4.0 5.0 3.9 5.5 4.5 5.3 4.5 2.5 5.0 4.5 5.4 4.9 4.5 3.3 2.3
##  [379] 4.2 4.9 5.5 5.1 5.5 5.3 4.6 4.0 5.7 3.3 5.2 4.3 4.9 5.4 3.3 5.3 4.0 5.4
##  [397] 4.3 4.1 3.0 4.4 5.1 6.3 2.2 5.3 4.4 4.0 4.0 4.4 5.0 4.0 5.4 3.2 2.8 5.8
##  [415] 5.1 4.5 5.2 5.6 4.8 3.8 3.9 3.0 6.6 5.0 3.5 4.6 3.1 4.9 5.2 4.3 4.6 3.7
##  [433] 3.7 3.9 3.8 2.1 5.7 2.5 4.8 4.5 4.6 3.7 3.9 4.5 4.2 3.8 4.6 4.3 4.2 4.6
##  [451] 3.3 4.3 4.4 5.0 6.0 5.7 5.3 5.5 5.1 4.6 4.7 4.0 4.8 5.9 3.6 4.8 6.6 3.7
##  [469] 5.8 5.2 4.1 3.7 5.9 4.9 4.2 5.1 4.7 6.5 4.5 6.3 3.4 4.1 4.4 4.4 3.9 4.9
##  [487] 5.0 4.3 4.4 2.4 4.9 5.3 4.6 3.2 3.2 5.1 5.5 5.0 4.9 4.7 4.0 5.5 4.7 5.7
##  [505] 2.1 4.4 4.4 5.6 5.3 4.2 4.8 4.0 4.4 4.0 4.7 5.1 5.0 4.0 4.8 5.1 3.5 4.9
##  [523] 4.7 4.5 6.1 4.5 5.0 5.3 4.4 2.8 4.0 4.7 5.3 5.5 4.9 4.6 3.1 3.6 2.7 4.3
##  [541] 3.7 3.7 3.5 4.9 5.1 3.9 4.9 3.9 4.4 4.5 3.8 5.3 5.2 4.1 5.6 4.8 2.9 4.4
##  [559] 4.3 5.0 4.3 5.8 5.5 5.3 4.6 3.8 4.3 4.2 4.5 3.9 5.2 2.8 2.8 6.7 4.4 5.5
##  [577] 3.4 4.2 5.9 4.4 5.8 4.0 4.7 5.2 5.2 5.1 4.4 5.9 4.7 2.9 4.9 3.8 5.8 3.4
##  [595] 4.1 4.7 3.5 3.7 5.3 2.9 4.5 5.3 5.7 4.5 5.3 4.9 5.8 5.5 4.5 5.4 3.1 5.0
##  [613] 4.9 6.1 5.6 4.9 3.6 5.2 4.2 4.5 5.5 3.9 3.8 2.6 3.8 3.8 3.9 4.4 2.8 4.7
##  [631] 3.3 5.1 5.4 5.0 3.2 4.4 4.2 3.4 4.2 5.1 4.5 5.7 4.8 5.2 5.8 3.1 5.9 2.9
##  [649] 5.6 3.6 5.0 5.5 5.4 4.5 2.9 4.7 4.2 4.9 5.9 5.3 4.7 4.5 3.7 5.8 5.5 4.0
##  [667] 3.2 2.7 3.9 5.2 5.0 2.7 3.2 4.3 4.5 3.8 5.9 4.7 4.5 4.9 3.1 5.1 5.6 4.9
##  [685] 4.8 4.1 5.0 4.2 4.0 3.9 5.0 6.4 4.5 4.0 3.6 4.9 4.3 4.8 4.0 3.6 5.1 2.6
##  [703] 4.7 5.6 5.5 5.8 4.1 4.8 3.6 6.8 4.2 4.4 3.0 3.9 2.8 4.5 3.9 4.3 4.2 5.9
##  [721] 4.1 4.7 3.1 3.5 3.5 5.7 5.2 4.3 5.0 5.4 4.1 4.7 4.9 4.8 4.0 4.7 4.3 3.7
##  [739] 3.4 5.2 3.6 5.5 5.3 5.7 5.4 4.0 4.7 3.8 4.0 4.4 4.8 5.7 4.2 3.5 4.9 4.1
##  [757] 5.8 4.5 5.0 4.2 4.6 6.3 5.1 5.3 4.6 6.1 5.7 3.8 4.2 4.4 4.3 3.4 3.5 5.2
##  [775] 5.6 3.9 4.8 4.7 5.3 4.5 6.5 3.2 4.7 3.8 4.8 5.0 3.0 5.4 4.4 6.2 4.5 4.8
##  [793] 5.7 4.2 3.4 4.1 3.1 5.5 4.7 5.1 3.9 3.9 4.3 5.0 5.1 6.5 4.2 5.3 3.4 4.5
##  [811] 5.3 6.1 5.3 5.7 4.9 4.5 2.9 4.2 6.2 4.0 2.8 3.7 5.0 4.6 5.6 5.1 4.1 3.9
##  [829] 3.8 4.1 4.2 5.9 4.1 2.8 5.3 6.3 4.2 4.2 4.9 3.9 4.0 4.3 4.3 4.5 3.6 4.4
##  [847] 6.0 4.4 3.9 5.0 4.6 5.7 3.2 4.8 2.9 4.4 3.6 3.1 3.9 5.2 4.7 5.2 3.7 3.9
##  [865] 5.7 5.2 3.5 5.2 5.0 4.1 4.6 5.1 4.6 3.3 4.4 3.6 4.3 5.1 4.7 5.5 5.2 4.1
##  [883] 4.2 2.7 5.2 4.1 4.9 3.8 5.0 4.0 3.1 4.6 3.4 4.2 3.7 4.0 4.7 5.8 6.2 3.4
##  [901] 4.6 4.9 5.4 4.2 3.9 3.5 2.9 4.4 5.6 4.3 3.5 3.7 5.3 6.0 2.9 5.7 6.9 3.5
##  [919] 3.8 5.0 5.1 4.7 3.7 5.3 4.8 5.1 4.5 5.6 3.3 5.2 5.0 4.3 4.3 4.1 4.3 4.4
##  [937] 4.7 3.0 4.9 3.8 6.0 4.5 6.5 4.1 4.9 3.8 5.7 3.8 5.7 5.9 5.6 4.5 4.8 2.2
##  [955] 3.7 4.8 4.3 5.3 4.7 4.0 3.7 2.9 5.7 6.3 4.0 5.4 2.7 4.3 3.7 3.3 3.9 4.2
##  [973] 5.4 3.6 3.3 4.5 4.3 5.3 3.8 4.9 2.7 2.7 4.8 5.7 4.7 4.7 5.5 6.3 4.4 3.4
##  [991] 6.3 4.4 3.7 3.8 4.2 3.4 3.7 5.9 3.3 5.1
## 
## $func.thetastar
## 72% 
## 5.1 
## 
## $jack.boot.val
##  [1] 5.500 5.500 5.300 5.300 5.100 5.000 4.900 4.900 4.636 4.500
## 
## $jack.boot.se
## [1] 0.9730579
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
## [1] -0.5707962
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
##   1.7144144   3.0302708 
##  (0.7047716) (1.4447768)
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
## [1] 0.4862155 0.8592702 0.7935357 1.0363249 0.6138119 0.8555230
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
##    [1]  0.234653841 -0.005404201 -0.935070742  0.321158884 -1.739221732
##    [6] -1.080589947 -1.023410569 -0.223881770 -0.551895086  0.024231523
##   [11] -0.568276364 -0.556071958 -1.077133096 -0.164222108 -0.467485540
##   [16] -0.556503671 -1.090117519 -0.581142344 -0.212168005  0.028098496
##   [21] -0.961211452 -0.481964619 -0.458716635  0.434160135 -2.470214265
##   [26] -0.813070471 -1.216275767  0.370383157 -1.366876244 -0.201588346
##   [31] -0.274933362  0.303817680 -0.111686953  0.241456820  0.246091631
##   [36] -1.160171429 -0.556686309  0.637809456 -0.991289712 -1.098600285
##   [41] -0.543881422 -0.268163220 -0.982048748 -0.532921322 -0.916765171
##   [46]  0.219850594 -0.030471942 -0.039223139 -0.659319032 -0.923321065
##   [51] -0.601432456 -0.115571708 -0.670903602  0.057554427 -1.058073957
##   [56] -1.229202045 -2.173655885  0.024634595 -1.023743567 -1.602432583
##   [61] -0.515891433 -1.037572283 -1.621692371 -1.055394331 -0.988116262
##   [66]  0.157710967 -0.562604822 -0.555212883 -1.351309935 -0.555573747
##   [71] -0.998786380  0.285548992 -0.307987494 -0.372602490 -0.057808935
##   [76] -0.665202302 -0.839340490 -0.912063576 -1.664120543 -0.517267708
##   [81] -0.601816324  0.106584867 -0.169524061 -0.183981655 -0.238485797
##   [86] -0.229929463 -0.326906316 -0.411343477 -0.906130921 -0.473606605
##   [91] -1.306642095 -0.829623639  0.522582384 -1.169406535 -0.944959785
##   [96] -2.074222671 -0.609235238 -1.420048512 -1.496712589 -1.023606745
##  [101] -0.148481460 -0.569646392 -0.099546307 -0.276849479 -0.101297904
##  [106] -0.974947299 -0.224814285 -1.450131125 -0.618311292 -0.499989645
##  [111] -0.447577780  0.087413909 -1.677963856  0.488573037 -0.918136798
##  [116] -0.954288414 -1.945288341 -0.370763287 -0.314161530 -2.501289700
##  [121]  0.077215632 -0.067043208  0.416555292 -0.077899488 -0.584191890
##  [126]  0.608807624 -1.185630078 -1.538471834 -0.601209281  0.502866740
##  [131] -1.374669041 -0.711696834 -1.138192521 -0.104424376  0.146160291
##  [136]  0.299219745 -0.997350110 -0.534635359 -0.562604822 -0.085825084
##  [141] -0.196021480 -0.613495365  0.705341958 -1.167968094 -0.170872249
##  [146] -2.300516399 -1.395742609  0.270874263  0.494756909  0.484582480
##  [151] -1.033764669 -1.040389663 -1.008778411 -0.200487825 -0.147249191
##  [156] -0.879139633 -0.956770412 -0.991287474 -0.990624397 -0.959882340
##  [161]  0.188937954 -1.329338027 -0.855369850 -0.522857949 -1.160171429
##  [166] -0.668560231 -0.258807199  0.318541192 -0.511738004 -0.914486840
##  [171] -1.134151402 -1.429755276 -0.533068353 -1.232393350 -0.963511108
##  [176] -0.937351818  0.455615882 -2.169871669 -1.355432926  0.155676133
##  [181] -0.378351326 -0.467907394 -1.450869124 -0.379602552 -0.901566491
##  [186] -0.448071968 -1.502492100 -1.064497818 -0.290518235 -0.494391127
##  [191] -0.290056305 -0.089388221 -0.576274325  0.313631001 -0.390623065
##  [196] -0.164567204 -0.912740436  0.761354903 -0.887026602  0.483802548
##  [201] -0.403834220  0.418599016 -0.891707899 -0.271358813 -1.536828421
##  [206] -0.462523056 -0.599076540  0.160845199 -0.524070903 -1.363105040
##  [211] -1.040389663 -0.822912806 -0.572355191 -2.102115016 -0.094560561
##  [216]  0.154632900 -0.382330797  0.408085278  0.551124304 -0.571982977
##  [221] -0.522440311 -0.861907562 -0.857912007 -0.126014441 -0.465807517
##  [226] -0.556422073 -0.301768053  0.180024605  0.897787624 -1.632493501
##  [231] -0.473059614 -0.843449223  0.170177432 -0.289047319 -0.504689908
##  [236]  0.291656369 -0.844334773 -0.171183149 -0.460060115 -0.364604874
##  [241] -1.541239178 -1.489470956 -0.321451835 -0.804036352 -0.494564303
##  [246] -0.149134861 -0.455575398 -0.111686953 -0.986897372  0.212768798
##  [251] -0.996003193 -1.779484731 -1.502492100 -1.608764261 -0.349084774
##  [256] -0.179090838  0.145224293 -0.522839580 -0.147531668 -0.871194440
##  [261]  0.266784805 -0.499626433 -0.945400208 -0.147705539 -0.628385882
##  [266] -0.293515592  0.557926340 -0.670129301 -0.485616698 -1.051153365
##  [271] -0.886610866 -0.548960480 -0.124584120 -2.366976602 -0.090750699
##  [276] -0.284205337 -0.971560190 -0.101297904 -0.805998990  0.135255022
##  [281] -1.057132795 -0.463605022  0.396368868 -0.136953727 -1.561695167
##  [286] -0.482538722 -0.990006586 -1.629416588 -0.311413965  0.321486625
##  [291] -0.532908271 -0.456934636 -0.497314308 -1.140113976 -0.849615625
##  [296] -0.437064960 -0.923617759 -0.549985557 -0.109876663 -0.980982817
##  [301] -0.422233664  0.289805640 -0.919920562  0.821428538 -0.181818506
##  [306] -0.209799676  0.706653320 -1.476023345 -0.243352699 -0.159566395
##  [311] -0.682889102 -0.883545802 -1.416581871  0.688028187 -1.287475242
##  [316] -0.033383123 -0.290182005  1.289183559  0.302703764 -0.135501184
##  [321] -0.801261238 -0.456552047 -1.583009313  0.262888991 -1.319791252
##  [326] -0.647925944  0.096136635  0.537036071 -0.795481554 -0.519039794
##  [331] -1.422577076 -0.108666395  0.688733239 -0.403046720 -0.704358827
##  [336]  0.321525375  0.265218990 -0.594791722 -0.494557590 -0.914585138
##  [341] -1.384212207 -0.172155241 -0.161435573 -1.057136762  0.039401776
##  [346] -0.393317275 -0.106351899 -0.825901792 -0.012813239 -0.388740799
##  [351] -0.055635740 -1.448626869 -1.937021626 -1.079017799 -0.566147154
##  [356] -0.124897674  0.272196285 -1.355904937 -0.479015056 -0.382550079
##  [361] -0.513772699 -1.050612582 -0.903537556 -0.855338148 -1.421575779
##  [366] -0.966003322 -0.402290727 -1.114950000 -2.114374418 -0.779475708
##  [371] -0.084642058 -1.509225098 -0.290034513 -1.145349807 -0.473500851
##  [376] -0.001835868 -1.133776810 -0.556523720 -0.277751154 -0.232632235
##  [381] -1.637216294 -0.896486013 -0.120201065  0.289842301 -1.548380369
##  [386] -0.710338477 -0.261089391 -0.993057707 -1.049513120 -1.072048991
##  [391] -0.755467334 -0.095395694 -0.337782920  0.585975619 -0.790345564
##  [396] -0.035004688 -0.460793632  0.291656369 -0.931684710 -0.748502107
##  [401] -0.516309979  0.100346378 -1.009488360 -0.684967121 -1.449701968
##  [406] -2.338745760 -0.941350954 -0.207084670  0.233995472 -2.327865595
##  [411] -0.786288215 -0.212059331 -1.763908031 -0.702800133 -2.235795008
##  [416] -1.393149180 -1.629590329 -0.517368449 -0.910164536 -0.140703136
##  [421] -0.408728405 -0.686523426 -0.397814357 -1.069866842  0.723217669
##  [426] -0.970483346  0.190061215 -0.950085476 -0.954859596  0.132078048
##  [431] -0.112159971 -0.728322337 -0.199333046 -0.882958531  0.034483141
##  [436] -1.435248150 -0.444029000 -0.967627617  0.127760295 -0.565680090
##  [441] -1.435179610 -0.162576244 -1.519183832 -0.220987528  0.900405413
##  [446] -1.478638127 -0.297007366 -0.841023890 -1.222323908 -0.478487154
##  [451] -0.618508435 -0.813287761 -0.265163892 -0.401971386 -1.030110138
##  [456] -0.049838995 -1.369358537 -1.499896364  0.283220035 -0.212674142
##  [461]  0.344048656 -0.121791508 -0.321336152 -1.023410569 -1.322365051
##  [466] -1.829850616 -2.513367552 -0.525079248 -1.293667347 -1.166837547
##  [471] -0.814980971 -1.399934297  0.072599407  0.183358505 -0.188889714
##  [476] -0.947978409 -1.038190077 -0.921106008 -0.904033645 -0.325820249
##  [481] -0.220082124 -1.366591674 -0.992505600 -0.109133115  0.368458707
##  [486] -0.524403285  0.023450897 -0.629388098 -0.057258122  0.289829403
##  [491] -1.499896364  0.150181736 -0.150959080  0.316334116  1.268545614
##  [496] -1.165088775 -1.414609182 -0.913476001 -0.206578655 -0.514968840
##  [501] -0.422189311 -0.767575313 -0.083535398 -1.291464392 -0.182447206
##  [506] -0.642355890  0.349354736 -0.166274732 -1.224880545 -0.568190284
##  [511] -1.527819812 -1.001783126 -1.474791872 -0.688018364 -0.215761080
##  [516] -0.886880616 -0.484958944  0.121821887 -1.026682387 -0.523991005
##  [521] -0.269483665 -1.483915126 -0.003597180 -0.164204366 -0.786131610
##  [526] -1.752197688  0.051820996 -0.439170764  0.718407456  0.169823727
##  [531] -0.411224179 -0.928699315 -0.536076398 -0.623211778 -0.009220203
##  [536] -0.665651936 -0.194082177 -0.104554619 -0.766968491 -0.842983409
##  [541] -0.160162143  0.434928919 -1.052251658 -0.944625976  0.015331079
##  [546] -1.564973996 -0.237006713 -0.928103343 -2.338533184  0.143549231
##  [551] -0.487114413  0.099926828 -0.525083955 -0.624918471 -1.123925620
##  [556] -2.511876411 -0.958260154 -0.516309979 -0.555816250 -0.160903734
##  [561] -2.097755888 -0.998326411 -0.374175484 -0.532189511 -0.616223076
##  [566] -0.826578140  0.236505866 -0.538331544 -0.094795990  0.465262681
##  [571] -1.616946990  0.098451669 -0.259672388  0.057882073 -0.894907943
##  [576] -0.422730493 -0.855227003 -0.134837124 -1.128891139 -0.997350110
##  [581]  0.709371376 -0.206509603 -1.393070156 -0.652135508  0.080319023
##  [586] -0.451188763 -0.084592356 -0.812023468  0.176513139 -1.311628691
##  [591]  0.008820376 -0.486569490  0.081816700 -0.338820949 -0.973877114
##  [596] -2.276889244 -1.660934966 -0.988252562 -0.498306433 -0.553214413
##  [601] -1.353670064 -0.853776063 -1.450943617 -1.087132570 -0.997279489
##  [606] -0.745913217 -0.126057264 -0.582438518 -0.089640418 -0.433845921
##  [611]  0.033678785 -0.031361618 -0.957527780 -0.804048959  0.020147438
##  [616] -0.863146111 -0.921966763  0.102880811 -0.457689164  0.099802473
##  [621] -0.448718881  0.126293018 -0.531540019 -0.993135049 -1.000299150
##  [626] -0.062674522  0.793167992 -0.554723298  0.896481367 -0.184268990
##  [631] -1.458648167  0.157917966 -0.401237828 -1.387755880 -0.257340451
##  [636] -1.391125975 -1.394596406 -0.911772757 -0.518707868 -0.138272420
##  [641] -0.833985714 -1.380780458  0.164502146 -2.379162636 -0.182447206
##  [646] -0.910656244 -0.252594998 -0.702609762 -0.427434788 -0.608786132
##  [651] -0.656969588 -1.069739988 -0.502425188 -0.149265704 -0.053838146
##  [656] -1.752197688 -0.561553093 -1.220393984  0.119861173 -0.419245104
##  [661] -0.101735918 -2.073162747 -0.855855815 -0.094012289 -0.009520767
##  [666] -0.499082727 -0.824048153 -1.412651696 -1.452068795 -0.495949299
##  [671] -1.387005515 -0.385738319 -0.697999146  0.045399289 -1.455683537
##  [676]  0.492452417  0.003594380 -0.612952410 -0.447314800 -0.874037740
##  [681] -0.835193567 -0.214843571 -0.291352526 -0.270063025 -0.734986692
##  [686] -2.201422022 -0.516433692 -0.067681517  0.919542060 -0.263668854
##  [691] -0.023437535 -0.269728490 -1.227884135 -0.863413477 -0.026908587
##  [696] -0.510170317 -0.024723133 -0.578199521 -0.214413815 -2.421824022
##  [701] -2.433987089 -0.104644139 -1.054646277 -0.175313022  0.078808889
##  [706] -0.258264231 -0.219325633  0.330980044  0.139116055 -0.428939792
##  [711]  0.054717167 -0.903105687 -1.526958719 -0.237627602 -0.522712455
##  [716] -0.851362165 -0.050739415 -0.630935525 -1.050313877 -0.547896519
##  [721] -1.116490827  0.110852305 -0.349331825  0.020937181 -0.602665428
##  [726] -0.943330144 -1.471827669  0.708771134 -0.825345180 -0.996356260
##  [731] -0.503800223 -0.570817942 -0.578033562 -0.987262929 -0.570608606
##  [736] -0.524096322 -0.463655195 -0.606320020 -0.141409384  0.015514493
##  [741] -2.394144003 -0.554608288  0.865943345 -2.095928624 -1.416549887
##  [746]  0.275992975 -0.823443603  0.081908590  0.191720162 -0.639637089
##  [751] -0.892995831 -0.218587260 -0.633729671  0.571902105 -2.099342824
##  [756] -0.534828844 -0.484529931 -0.861529450 -0.524403285 -0.920966336
##  [761] -2.417541964 -0.511401261 -0.528395406 -2.235783487  0.302610330
##  [766] -0.475572110 -0.471326483 -1.163086716 -1.238639234 -1.412869699
##  [771] -0.912307470  0.585800451 -1.078511243 -1.232393350  0.345569849
##  [776] -0.174698000 -0.442155166 -0.454877703 -0.440825579 -0.796395670
##  [781] -0.460382535 -0.240968998 -0.524190722  0.123242875 -1.166999106
##  [786] -2.527841904 -0.563248553  0.869741840 -0.960346659  0.533273875
##  [791] -0.372543605 -0.526971515 -0.367975338 -0.101752482 -0.255919151
##  [796] -1.668695108 -0.207169383 -0.653655792  0.070747686 -2.171355074
##  [801] -0.049936718 -1.247244520 -0.513597331 -0.239300383 -0.599330536
##  [806] -0.760122031  0.139493674 -0.627206795 -0.619184924 -0.653161846
##  [811] -0.252059703 -0.597749600 -0.905609367 -2.526758168  0.096743375
##  [816] -0.420037928 -0.902810709 -1.053448253 -0.092755963 -0.501435571
##  [821]  0.079764616 -1.541425386 -0.325571044  0.538820795 -0.475575634
##  [826] -0.978731928 -0.487059430 -0.302367165  0.074871804 -1.766277048
##  [831] -1.543363867 -0.857399532 -1.063165334 -0.459400552 -0.272766634
##  [836] -0.306751755  0.083173321 -0.135590675  0.684425881  0.126733274
##  [841] -0.257061587 -1.001249116 -0.448710677 -0.329100773 -0.246126732
##  [846] -0.674074263 -1.485621281 -0.942384858 -0.576911388 -0.893530644
##  [851]  0.122574658 -0.937878170 -0.167277733 -0.049486582 -0.223419918
##  [856] -1.439876301 -1.491402047 -0.157421796 -1.955324203 -0.789390603
##  [861] -0.525220094 -0.991847810  0.192927450 -1.067237700 -0.496420456
##  [866] -0.230032908 -0.584684833 -0.468920519 -0.911584873 -0.475434526
##  [871]  0.115437679 -0.755805564 -0.567808880  0.132036552 -0.249885481
##  [876] -0.204178832  0.237514680 -2.112968246 -0.513816020  0.074392128
##  [881] -0.627666216 -0.051594417 -0.385845989 -1.360414222 -0.498538660
##  [886] -1.632629626 -0.138252373 -0.049668581 -0.257329368 -0.110323631
##  [891] -1.154645793 -0.324627691 -2.047623067 -1.442119100 -1.577553051
##  [896] -0.040562775  0.376215129 -2.486771689 -0.142380260 -0.091891926
##  [901] -0.228042565 -0.847844191 -0.523383746 -0.613010666  0.760508043
##  [906] -0.765628828 -0.197536990  0.084092095 -0.230205939 -2.332193782
##  [911] -0.416630743 -0.514340598 -0.775855536  0.121035413 -0.629052253
##  [916] -0.820550617 -0.521897104 -1.095358715 -0.930242155 -1.161980011
##  [921] -1.496383959 -0.600058769 -0.925342019 -0.241688534 -0.215299161
##  [926] -2.520506375 -0.968256365  0.008657759 -0.122899607 -0.777750153
##  [931] -1.257859098 -0.891597150 -0.486569490 -0.675723338 -0.311064784
##  [936] -0.576989058 -0.560128246 -1.456078554 -2.327177501 -1.086160253
##  [941] -1.315463216 -0.100527671 -0.981670525  0.158259070 -0.883320401
##  [946] -0.555030052 -0.412915653  0.076112717 -0.485071444 -0.448215266
##  [951]  0.486487448 -1.006877991  0.112697927 -0.987822753 -0.147364769
##  [956] -1.417625886 -0.963127966  0.299004263  0.242147176 -0.080228659
##  [961] -0.533798751  0.184069467  0.366772580 -1.452335783 -0.524391879
##  [966] -0.380001018 -1.507946346 -0.671345616 -1.272735283 -0.507959122
##  [971] -0.114365409 -0.824514742  0.100829367  0.218328303 -0.128559258
##  [976] -0.166468216 -0.202162038 -1.453899115 -0.185805832 -0.485577838
##  [981]  0.237910412 -0.613521745 -1.330526706 -0.274599329 -0.144319316
##  [986] -1.292541509 -0.136135039 -0.154755927 -1.453313413 -0.442575010
##  [991] -0.405746855 -0.354595645 -0.116504321 -0.774249410 -0.443933411
##  [996] -1.679960980 -0.921780850 -0.558759494 -0.441334107  0.143754521
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
```

```r
fit2
```

```
##       mean          sd    
##   0.56576911   0.30412899 
##  (0.09617403) (0.06800328)
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
## [1]  1.1952908  0.5265331  0.2361982 -0.8943204 -0.5469568  0.1211282
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
## [1] 0.0565
```

```r
sd(results2$thetastar)  #the standard deviation of the theta stars is the SE of the statistic (in this case, the mean)
```

```
## [1] 0.9058025
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
## t1*      4.5 -0.01111111   0.9176878
```

One of the main advantages to the 'boot' package over the 'bootstrap' package is the nicer formatting of the output.

Going back to our original code, lets see how we could reproduce all of these numbers:


```r
table(sample(x,size=length(x),replace=T))
```

```
## 
## 2 3 4 6 8 9 
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
## [1] -0.0376
```

```r
se.boot
```

```
## [1] 0.9249767
```

Why do our numbers not agree exactly with those of the boot package? This is because our estimates of bias and standard error are just estimates, and they carry with them their own uncertainties. That is one of the reasons we might bother doing jackknife-after-bootstrap.

The 'boot' package has a LOT of functionality. If we have time, we will come back to some of these more complex functions later in the semester as we cover topics like regression and glm.

