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
## 1 3 5 8 9 
## 2 2 3 1 2
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
## [1] -0.0459
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
## [1] 2.679588
```

```r
UL.boot
```

```
## [1] 6.228612
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
##    [1] 3.3 4.5 4.9 5.8 3.0 5.1 4.6 5.3 3.8 3.8 4.4 2.8 4.2 5.0 5.5 3.5 5.0 4.1
##   [19] 2.7 4.6 4.9 5.5 4.2 5.9 4.5 4.2 4.0 3.8 5.6 2.9 4.1 4.0 3.4 3.7 5.0 4.3
##   [37] 3.3 5.4 4.2 3.1 4.2 4.3 5.6 4.5 4.7 5.0 6.3 5.2 4.6 5.1 3.8 5.4 3.6 4.3
##   [55] 3.3 6.5 5.9 5.5 4.7 5.2 5.2 5.1 3.5 4.1 5.8 3.0 4.3 5.8 4.9 3.0 4.6 5.5
##   [73] 4.8 4.0 4.7 4.8 5.2 5.1 4.0 4.5 3.9 3.1 4.3 4.9 5.2 5.9 3.0 3.6 3.5 5.6
##   [91] 5.4 6.5 4.1 3.1 3.3 4.4 5.5 5.1 4.4 4.6 5.5 4.9 5.5 3.4 3.4 3.1 5.5 4.5
##  [109] 4.1 4.4 5.2 3.3 4.8 5.1 5.0 3.1 3.6 5.5 5.7 5.0 4.9 5.0 6.0 3.5 3.9 5.3
##  [127] 4.5 4.5 3.7 5.2 4.4 3.8 3.8 3.5 5.7 4.1 4.1 5.4 4.7 4.8 4.7 3.8 5.0 6.1
##  [145] 4.8 6.3 4.1 4.0 4.3 5.4 3.4 4.2 3.6 3.3 4.7 6.3 5.1 4.8 4.7 2.6 5.2 3.2
##  [163] 2.4 5.5 5.2 3.2 4.0 6.4 5.2 5.1 5.0 3.3 4.6 4.9 5.1 6.5 4.6 4.9 5.7 3.2
##  [181] 5.6 5.3 3.3 5.1 4.0 5.2 3.1 3.9 3.9 5.6 3.4 5.2 5.7 4.4 4.4 5.5 2.9 4.6
##  [199] 4.4 4.2 4.7 3.9 4.1 5.6 6.0 4.2 2.7 5.9 3.7 3.5 3.2 4.8 3.3 4.7 4.7 4.7
##  [217] 5.1 4.7 3.6 3.6 5.6 4.9 2.8 3.2 3.2 6.1 5.6 3.9 3.8 2.9 4.7 2.7 4.9 4.3
##  [235] 4.6 4.4 5.6 3.4 4.5 3.5 3.9 4.7 4.1 3.8 3.7 3.8 5.1 5.1 4.0 3.1 5.3 3.3
##  [253] 2.9 4.2 4.1 5.0 5.3 4.4 4.9 4.6 4.1 5.3 4.4 3.6 4.8 4.5 3.8 4.3 6.2 4.2
##  [271] 6.4 3.6 5.6 3.3 3.8 5.2 2.8 6.1 4.7 4.7 3.6 5.9 5.5 4.6 4.6 3.5 6.5 5.5
##  [289] 3.3 4.8 5.6 5.1 3.6 3.5 3.6 4.8 3.8 4.3 5.7 4.4 3.9 6.2 4.3 5.3 4.6 2.4
##  [307] 3.0 3.6 4.4 4.1 3.6 2.5 3.6 3.5 4.2 4.4 5.7 5.0 4.3 5.1 5.1 4.8 3.9 4.4
##  [325] 5.0 5.2 5.1 3.9 4.1 5.3 4.2 4.2 4.4 3.4 4.4 4.7 3.8 4.8 3.8 4.8 4.0 5.2
##  [343] 4.7 4.8 3.0 6.0 4.9 4.3 3.5 5.6 5.5 5.0 5.0 4.6 5.0 5.8 3.3 4.8 5.2 5.8
##  [361] 5.3 3.9 5.1 4.8 3.9 6.4 5.9 3.8 3.5 5.1 4.9 4.7 4.5 4.9 4.7 4.0 3.5 4.4
##  [379] 5.4 4.5 4.5 3.7 3.4 4.5 6.2 3.0 4.5 3.6 6.0 4.1 4.7 3.2 5.0 5.1 2.0 4.3
##  [397] 3.7 4.5 4.3 5.4 5.6 4.2 4.0 3.7 4.1 3.8 5.0 5.6 5.4 5.6 3.7 3.9 4.9 5.6
##  [415] 4.0 7.0 5.8 3.2 6.9 4.8 7.2 4.5 4.8 4.1 3.5 5.7 5.3 3.8 2.6 4.7 3.9 4.4
##  [433] 5.2 5.0 4.6 4.0 5.0 4.6 3.1 5.0 3.7 4.8 3.2 4.7 4.3 4.6 5.0 4.5 4.0 3.7
##  [451] 3.8 3.9 3.9 4.3 5.5 3.5 4.6 3.4 4.3 4.9 5.1 5.2 6.2 5.0 3.1 3.1 4.9 4.2
##  [469] 3.5 4.6 4.2 4.1 4.8 4.0 4.5 3.6 2.5 6.8 4.5 5.0 3.9 4.3 5.0 3.9 3.8 5.3
##  [487] 4.0 4.0 4.8 3.8 5.6 6.3 3.8 4.3 4.1 5.6 5.1 5.6 3.7 4.9 3.3 6.4 5.1 4.7
##  [505] 4.6 5.4 5.0 4.2 5.0 3.2 4.2 6.0 3.5 5.4 4.3 3.4 3.8 3.8 4.7 5.1 3.9 2.1
##  [523] 4.3 4.7 5.2 3.1 3.0 5.1 5.1 4.3 3.5 5.4 4.8 3.7 4.4 3.3 6.0 5.3 4.2 5.7
##  [541] 4.1 3.5 4.3 3.4 5.3 5.2 5.8 3.5 3.9 6.3 5.6 4.5 4.9 3.4 3.6 5.1 4.4 4.5
##  [559] 6.4 5.5 3.2 6.0 3.7 4.2 4.8 4.8 4.2 4.1 4.4 3.5 3.7 3.2 5.1 5.1 5.0 5.0
##  [577] 4.2 3.5 5.1 3.1 2.7 4.7 5.3 5.8 5.4 3.4 4.7 6.1 6.1 5.4 4.6 6.0 6.7 4.2
##  [595] 5.2 4.7 3.9 4.1 5.0 4.6 5.1 4.9 4.6 4.7 4.6 2.3 4.9 6.0 4.2 4.7 3.5 7.1
##  [613] 5.4 4.6 6.2 2.8 4.6 2.7 5.6 3.7 4.2 3.2 5.6 5.6 4.7 3.7 3.3 6.4 5.2 4.6
##  [631] 4.3 2.8 3.3 3.8 4.7 4.8 4.1 4.4 4.0 4.3 4.1 5.0 4.1 3.5 4.2 3.3 5.9 4.0
##  [649] 4.4 4.8 4.1 3.9 4.1 4.3 5.3 4.7 7.2 4.2 5.6 3.3 4.9 4.4 5.3 3.6 3.7 4.7
##  [667] 3.4 5.2 4.7 4.7 3.9 4.1 4.0 3.0 2.9 3.7 4.9 4.6 4.5 5.0 3.4 5.2 4.7 4.8
##  [685] 4.7 4.6 5.3 3.7 4.8 4.1 3.5 4.6 4.5 5.8 3.6 4.9 4.9 5.9 5.5 5.8 3.9 3.4
##  [703] 5.0 5.1 3.4 2.8 4.6 5.1 4.7 3.0 3.4 3.7 4.8 3.6 4.9 5.5 3.5 4.6 5.5 4.3
##  [721] 2.3 2.0 3.0 5.2 4.8 4.9 6.1 6.1 4.4 1.8 5.6 5.1 6.5 5.6 3.8 4.7 6.1 3.5
##  [739] 5.2 3.5 5.2 3.5 4.9 3.0 5.5 6.3 5.8 5.2 4.7 5.1 4.5 5.4 4.5 3.9 4.0 3.0
##  [757] 5.1 3.2 4.9 5.5 3.5 5.9 6.0 5.5 5.1 4.2 5.7 5.2 4.6 5.8 4.0 5.1 5.6 4.0
##  [775] 4.1 4.3 4.8 4.9 3.4 5.8 5.1 4.0 4.1 3.7 3.6 5.6 5.0 5.9 4.4 5.8 3.6 5.2
##  [793] 4.1 4.5 5.7 5.1 3.4 5.5 3.3 4.4 4.5 5.4 2.8 5.0 5.0 4.6 4.3 4.2 3.6 4.5
##  [811] 5.6 5.0 3.7 3.5 5.1 4.0 5.1 6.3 4.9 4.5 5.9 3.8 5.3 5.1 5.5 3.8 5.3 4.7
##  [829] 4.0 3.4 4.1 3.8 4.0 6.1 5.3 5.0 4.3 4.5 5.0 5.1 3.8 4.9 5.1 4.4 5.3 4.5
##  [847] 5.1 4.4 5.9 5.5 4.5 4.0 5.6 5.0 5.1 4.6 4.1 4.4 5.7 5.4 3.2 4.5 5.3 5.4
##  [865] 3.2 5.0 4.9 7.0 3.9 3.9 4.7 4.9 5.1 4.9 3.5 5.6 4.6 3.2 5.7 5.1 4.6 4.5
##  [883] 4.5 3.6 5.4 4.4 3.6 5.4 4.1 3.7 2.8 5.2 4.4 5.0 5.6 3.7 4.8 4.8 4.9 4.0
##  [901] 4.8 4.9 6.3 4.9 5.6 3.6 4.1 3.4 3.9 4.4 4.0 3.3 4.0 4.8 5.8 4.0 4.9 4.6
##  [919] 4.9 4.2 3.9 4.8 5.3 5.7 5.7 5.1 3.9 5.0 3.8 3.9 2.8 5.4 5.4 4.5 6.1 5.8
##  [937] 3.3 3.4 3.3 4.2 4.6 5.7 3.2 1.9 4.4 4.3 5.4 3.5 5.9 5.9 3.9 4.3 2.6 4.5
##  [955] 3.2 4.6 6.2 3.0 4.9 5.3 5.7 3.7 5.7 5.2 5.0 3.3 6.6 3.9 4.7 3.8 5.3 4.2
##  [973] 4.2 4.6 3.7 5.1 4.3 6.5 3.4 4.0 3.9 5.4 4.1 5.4 3.5 5.1 5.3 4.7 5.0 4.7
##  [991] 5.0 5.2 3.7 5.8 4.7 5.8 3.0 3.9 5.2 4.9
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
##    [1] 2.4 4.1 4.9 3.8 4.2 4.2 3.7 4.1 6.0 4.2 5.2 2.0 4.0 4.3 4.5 4.5 5.6 4.8
##   [19] 3.5 3.4 4.1 5.5 4.2 3.7 3.4 4.6 4.3 5.7 3.6 3.7 5.9 4.7 4.1 4.3 5.6 4.3
##   [37] 5.0 4.8 5.7 3.3 4.9 6.5 4.7 4.2 5.0 4.4 4.6 4.1 4.4 6.8 4.7 4.8 3.3 4.0
##   [55] 4.9 4.7 4.2 4.2 5.0 4.8 4.6 5.0 1.9 5.4 4.4 5.8 3.2 5.9 3.9 6.5 5.3 4.9
##   [73] 4.2 2.5 4.6 4.6 3.6 4.7 5.2 2.1 5.8 3.5 4.8 5.8 3.3 6.8 3.9 5.0 5.8 3.3
##   [91] 5.7 4.4 4.7 3.0 5.2 4.9 3.1 4.7 5.1 4.7 6.3 2.9 2.9 4.3 6.5 2.8 4.8 5.8
##  [109] 5.1 4.1 4.2 6.0 5.9 4.3 3.7 5.0 4.3 6.0 4.1 4.2 4.8 5.1 5.6 5.5 4.0 4.5
##  [127] 4.8 3.1 3.9 4.1 3.3 4.6 4.8 4.0 4.4 4.6 6.1 3.9 4.5 5.6 4.3 5.1 2.7 2.1
##  [145] 2.7 3.2 5.7 4.8 4.3 3.4 4.8 3.9 4.8 4.0 5.0 3.8 5.1 4.3 3.9 5.1 4.8 6.0
##  [163] 5.0 5.5 4.9 4.6 3.3 5.3 4.1 2.6 4.8 4.8 3.8 6.0 4.8 5.0 4.8 5.4 4.6 4.3
##  [181] 4.7 4.1 3.4 3.7 4.3 4.1 5.4 4.8 2.2 4.9 3.5 4.5 3.6 5.7 3.7 3.4 4.1 4.5
##  [199] 3.9 5.3 4.7 5.6 4.7 3.6 4.8 3.8 4.1 3.1 4.5 3.9 6.0 4.7 4.7 6.0 4.2 3.6
##  [217] 5.8 3.0 5.4 5.4 6.1 5.6 5.6 5.4 4.5 3.8 4.1 4.2 3.7 4.1 4.8 4.3 3.4 4.0
##  [235] 4.8 4.8 4.3 2.5 2.8 4.9 4.3 5.3 4.7 4.2 3.5 2.1 5.0 3.1 4.8 3.9 2.6 5.1
##  [253] 5.3 4.6 5.3 4.3 4.1 4.3 4.3 3.3 6.4 5.0 3.3 5.1 4.6 3.9 4.4 5.1 4.4 4.3
##  [271] 3.5 4.5 4.9 5.7 3.6 4.6 5.5 3.7 4.3 4.0 4.3 3.3 5.1 3.2 2.2 4.0 5.1 2.9
##  [289] 5.1 2.9 5.0 5.6 5.2 4.4 4.7 3.8 7.4 5.7 5.2 4.3 3.6 5.9 4.5 4.0 4.2 4.6
##  [307] 3.3 2.9 5.5 5.2 4.2 5.4 5.2 4.5 3.6 5.9 5.0 6.1 4.2 4.1 3.6 5.4 5.1 4.6
##  [325] 5.3 5.0 3.6 5.5 5.1 6.2 4.1 3.8 3.8 3.6 3.9 5.7 3.5 3.1 4.8 5.4 4.5 3.2
##  [343] 5.8 3.7 5.6 5.3 4.2 5.2 3.8 3.8 5.0 4.4 3.9 3.9 5.1 3.9 4.2 4.4 4.4 5.3
##  [361] 4.7 3.7 4.9 3.5 4.2 5.7 5.6 4.9 5.2 5.7 4.9 4.3 4.8 4.4 4.3 5.5 4.4 4.3
##  [379] 5.1 4.2 5.8 5.3 4.2 5.1 4.5 4.4 3.3 5.3 5.9 4.1 5.7 4.6 3.7 2.3 5.2 4.7
##  [397] 3.6 4.2 2.7 5.6 4.7 5.9 4.6 5.7 5.8 5.3 5.0 5.1 4.2 4.4 3.5 2.5 4.5 4.2
##  [415] 6.3 3.6 5.1 2.8 5.0 4.6 4.7 4.5 4.8 4.0 5.7 3.4 4.6 4.9 6.0 4.7 4.1 4.2
##  [433] 3.3 4.6 4.4 5.4 4.3 5.2 4.9 5.9 5.9 3.7 3.3 4.1 5.0 2.7 3.9 6.5 5.2 5.8
##  [451] 4.0 4.0 4.4 4.6 4.9 5.8 4.3 3.4 3.7 4.5 5.1 5.9 5.8 5.3 3.4 4.9 3.4 5.3
##  [469] 5.3 4.8 5.1 5.7 3.1 4.8 5.5 4.9 4.8 4.0 4.6 4.2 4.1 4.7 5.4 4.5 5.3 4.1
##  [487] 3.5 4.3 5.3 6.1 4.2 4.7 4.4 5.0 3.4 5.4 5.3 5.0 4.7 6.1 3.5 1.7 6.1 4.4
##  [505] 4.3 3.8 4.7 3.2 4.8 4.3 3.9 5.2 4.6 4.8 6.2 4.9 4.0 6.0 4.1 4.3 4.3 4.4
##  [523] 3.3 4.2 3.2 4.4 3.2 3.7 3.7 5.5 4.5 4.3 6.2 4.0 4.2 3.1 3.6 4.2 4.2 3.9
##  [541] 5.5 5.3 5.6 5.5 3.9 5.1 4.7 4.6 5.9 2.8 5.2 3.5 3.6 5.5 3.4 5.9 4.2 3.3
##  [559] 4.5 3.6 4.2 4.2 3.5 5.1 2.8 4.5 3.2 4.4 5.4 4.8 5.6 3.6 5.0 4.9 4.4 6.0
##  [577] 5.9 4.9 4.5 2.3 5.8 4.2 3.9 5.0 3.5 5.0 5.4 4.8 5.7 4.9 2.5 4.8 5.3 6.0
##  [595] 5.0 4.3 4.0 3.7 4.2 5.6 3.8 3.8 4.8 5.9 5.2 3.8 5.0 5.1 4.9 5.7 5.1 4.4
##  [613] 4.9 4.6 4.1 4.0 4.3 4.9 4.1 4.9 4.8 5.5 4.9 5.5 2.9 5.3 3.2 4.5 3.9 5.0
##  [631] 5.9 4.4 5.7 3.6 6.6 3.2 5.8 3.9 4.3 3.0 3.5 1.4 5.2 4.0 4.2 5.2 4.6 3.0
##  [649] 5.9 3.7 4.6 4.7 5.5 3.9 3.7 5.3 4.6 5.3 5.2 4.3 7.0 3.7 6.0 4.4 4.1 5.0
##  [667] 4.4 4.1 3.6 4.3 4.5 2.8 5.2 5.7 4.3 4.7 4.0 3.5 4.8 4.7 5.9 2.8 6.6 3.4
##  [685] 3.5 6.0 3.8 3.1 4.6 6.5 3.1 4.3 5.3 5.3 5.5 5.2 4.3 4.6 4.0 5.9 4.0 4.6
##  [703] 3.6 4.8 5.3 4.3 6.6 4.6 5.9 4.9 5.2 5.8 3.6 4.6 4.1 5.1 4.1 4.5 4.6 4.4
##  [721] 4.9 3.1 3.4 4.2 4.3 5.0 3.9 3.6 3.6 3.9 4.3 5.2 4.7 5.5 4.1 6.4 4.3 4.9
##  [739] 3.8 4.7 3.2 4.8 6.1 4.3 5.7 5.9 4.2 5.6 3.4 4.8 4.8 4.0 4.6 4.7 4.8 2.8
##  [757] 5.2 3.6 4.3 4.2 4.3 3.4 3.5 3.2 3.9 5.9 5.3 5.3 4.4 4.0 3.3 5.9 5.1 3.8
##  [775] 4.4 5.7 2.9 3.1 4.4 2.9 3.8 3.8 3.8 4.3 4.9 6.1 3.6 2.3 5.1 4.1 3.9 4.8
##  [793] 4.0 6.3 5.2 3.2 4.4 3.0 4.9 4.9 3.8 6.0 3.7 4.3 4.9 5.6 4.5 3.7 2.3 3.7
##  [811] 4.7 3.3 2.8 3.7 6.0 4.8 3.5 5.9 4.5 4.4 4.6 4.0 4.1 3.8 4.0 2.8 5.5 3.3
##  [829] 4.1 1.9 4.1 3.0 3.8 5.9 4.3 5.4 5.9 4.6 3.9 5.1 5.5 4.5 4.1 5.3 5.9 3.8
##  [847] 5.7 3.9 3.5 4.3 3.4 5.5 5.6 5.8 4.0 2.7 3.8 5.0 3.5 4.1 2.6 4.1 6.2 3.6
##  [865] 3.8 3.8 4.2 5.1 4.7 4.8 3.5 4.1 3.6 4.9 4.4 5.5 2.9 4.1 3.7 2.8 5.5 5.1
##  [883] 3.4 4.3 2.9 3.5 5.7 2.9 4.5 5.6 3.6 3.9 5.0 4.4 3.1 4.7 3.4 2.8 2.9 4.5
##  [901] 4.4 5.5 4.6 4.6 4.1 4.0 3.7 3.3 2.4 4.9 5.1 5.0 3.5 5.4 5.0 4.5 5.3 4.0
##  [919] 5.7 4.6 3.1 2.6 5.0 5.1 4.4 4.4 5.6 4.2 6.1 3.9 5.2 4.4 3.5 4.7 5.5 3.3
##  [937] 4.8 3.3 4.4 3.8 3.7 3.7 3.1 5.4 3.5 5.1 5.2 3.8 3.9 4.3 4.3 2.7 4.8 1.8
##  [955] 4.7 2.8 3.9 4.5 4.2 3.9 4.9 4.1 4.6 5.1 5.7 5.2 3.1 5.2 3.8 4.5 3.2 3.8
##  [973] 4.7 4.0 5.3 4.5 4.4 3.8 4.3 6.0 4.3 4.0 4.6 5.8 4.0 3.4 3.0 5.3 5.5 4.0
##  [991] 4.9 4.2 4.8 4.4 2.1 4.0 3.7 5.0 4.4 4.1
## 
## $func.thetastar
## [1] -0.0329
## 
## $jack.boot.val
##  [1]  0.43126844  0.38166189  0.28937500  0.08597884  0.01671309 -0.06428571
##  [7] -0.23897436 -0.32478632 -0.40332410 -0.58338279
## 
## $jack.boot.se
## [1] 0.9840783
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
##    [1] 5.2 5.3 4.9 4.5 5.7 4.0 4.5 4.3 4.1 5.1 4.0 5.0 4.9 4.0 4.2 3.3 5.2 3.6
##   [19] 2.7 6.2 4.9 3.2 5.3 3.2 4.6 3.2 5.8 4.7 2.5 4.5 4.4 5.7 4.0 4.0 3.8 2.9
##   [37] 3.5 4.7 5.8 4.6 4.4 4.0 4.6 3.1 5.0 4.4 5.2 4.6 5.2 5.7 4.6 5.3 4.5 3.4
##   [55] 3.5 5.3 3.8 6.4 5.3 4.2 4.2 3.6 4.1 5.4 3.1 6.1 4.4 4.3 4.3 5.4 4.0 3.8
##   [73] 4.5 1.9 5.9 4.8 5.6 2.9 5.3 4.8 4.4 5.1 5.4 5.6 5.7 3.1 5.4 2.7 3.7 3.3
##   [91] 4.0 4.0 3.3 4.8 5.6 6.7 4.1 3.7 3.4 4.1 4.0 3.7 4.5 5.2 3.3 3.7 4.0 3.3
##  [109] 3.3 3.4 3.5 5.1 3.5 6.1 2.4 4.5 5.3 4.6 4.0 5.8 3.2 5.4 5.0 5.5 3.6 2.9
##  [127] 6.2 5.0 4.0 6.4 3.5 3.7 4.8 4.7 5.6 6.2 4.6 4.6 5.6 4.5 4.6 4.1 4.7 5.0
##  [145] 4.6 5.0 4.1 3.8 4.8 5.2 4.2 4.7 4.9 4.5 4.6 4.0 5.3 3.8 5.3 5.2 5.2 4.7
##  [163] 5.1 6.3 5.0 5.1 2.7 5.7 6.5 4.5 4.0 3.3 3.3 3.8 5.4 3.5 5.2 4.5 3.4 3.3
##  [181] 3.5 5.1 3.9 3.1 6.0 3.3 4.1 4.2 4.5 4.9 3.8 4.0 3.9 5.1 6.2 3.9 4.2 5.7
##  [199] 4.4 4.1 3.8 4.6 5.3 2.7 4.5 5.8 2.9 5.0 4.5 4.1 5.7 4.5 4.6 3.5 4.6 5.6
##  [217] 3.3 3.5 5.5 3.7 4.8 4.7 3.7 4.2 4.7 2.8 5.1 4.3 3.9 4.0 5.3 4.1 3.8 3.9
##  [235] 4.5 4.2 4.5 4.1 5.0 5.1 2.9 5.2 4.9 3.6 4.0 4.9 5.5 5.2 4.5 3.9 5.5 3.8
##  [253] 4.7 5.5 5.0 4.3 4.0 3.9 6.2 4.5 4.4 4.7 4.9 3.7 3.2 4.1 6.4 3.8 4.1 3.8
##  [271] 3.7 3.5 4.4 3.0 3.8 3.8 4.7 4.3 5.0 5.4 4.9 4.3 3.3 3.5 5.9 6.8 5.1 4.8
##  [289] 4.5 3.7 5.2 4.8 3.7 4.3 3.4 5.8 3.5 4.6 4.3 3.8 5.0 6.2 4.5 5.3 4.5 3.6
##  [307] 2.5 5.3 5.4 5.2 5.3 5.5 5.1 5.2 5.4 5.8 5.5 4.7 3.8 4.9 4.4 5.3 3.8 3.6
##  [325] 3.3 4.2 3.4 5.0 4.5 3.3 5.5 4.2 5.7 4.6 4.4 3.5 4.0 3.0 4.2 4.8 5.9 4.5
##  [343] 5.4 6.3 5.0 4.1 6.1 5.2 3.4 3.7 3.8 6.4 3.9 3.2 5.5 2.9 5.7 4.2 4.5 4.5
##  [361] 4.0 5.0 5.2 4.3 4.7 4.8 2.9 3.5 4.9 5.9 4.3 3.8 5.4 4.1 5.1 5.7 3.8 4.5
##  [379] 4.9 5.1 4.5 5.0 4.7 5.3 4.6 5.2 5.6 4.7 2.4 5.0 3.8 4.6 4.9 3.2 4.8 4.6
##  [397] 4.0 4.1 5.4 3.2 4.2 3.9 5.2 3.2 4.3 3.9 3.6 5.1 3.4 5.8 4.7 3.7 4.4 5.9
##  [415] 4.8 4.6 3.5 5.1 5.4 4.0 4.9 2.4 4.6 4.8 2.4 4.4 4.2 4.1 6.0 4.5 3.7 5.4
##  [433] 4.7 5.7 4.9 4.2 3.5 4.5 5.5 4.3 6.1 4.5 3.8 5.2 6.5 3.4 3.6 4.2 5.4 6.3
##  [451] 3.2 6.7 4.0 4.4 6.1 4.3 3.8 3.6 3.9 4.2 3.3 4.2 4.6 5.1 5.5 4.7 5.9 5.0
##  [469] 5.7 4.1 5.3 4.7 5.3 5.6 4.3 5.0 3.5 5.4 5.3 2.5 4.1 5.2 5.6 4.5 5.2 4.2
##  [487] 4.8 6.3 4.2 4.7 4.6 5.4 4.6 3.6 4.7 4.9 5.8 4.0 4.0 3.6 5.9 4.2 5.7 5.3
##  [505] 4.7 3.8 4.6 5.2 3.9 4.3 1.8 5.0 5.5 5.3 4.0 3.7 5.0 6.0 4.5 4.3 4.1 4.4
##  [523] 4.8 5.9 4.0 6.0 4.5 6.3 5.2 4.3 4.4 4.8 4.6 5.4 3.8 4.3 3.5 4.8 4.1 4.2
##  [541] 5.4 5.3 4.9 4.8 3.5 4.9 4.3 5.5 4.8 4.6 4.1 5.4 5.3 4.6 3.3 3.0 4.2 4.0
##  [559] 5.9 3.7 4.3 3.5 4.5 3.8 3.6 4.4 4.9 3.9 3.7 4.1 5.5 4.9 5.1 3.3 5.3 4.5
##  [577] 3.4 5.0 4.0 5.5 5.4 6.0 4.7 2.8 4.4 4.6 4.2 4.9 3.8 3.7 4.8 4.4 4.2 5.1
##  [595] 4.5 4.9 3.7 4.9 5.9 4.4 6.2 5.1 4.2 5.4 3.4 5.0 4.9 3.4 4.4 4.4 3.8 4.9
##  [613] 5.6 4.1 4.6 5.4 4.8 2.2 6.2 3.4 4.0 5.6 4.7 4.0 3.5 4.3 4.8 4.1 4.4 4.1
##  [631] 5.4 4.1 4.5 5.6 4.7 4.0 5.2 5.2 4.1 4.6 3.9 4.3 4.6 4.9 3.2 2.5 5.0 4.9
##  [649] 3.4 5.3 3.6 3.9 4.8 5.1 5.5 5.1 2.6 3.9 5.3 3.8 4.4 5.9 5.1 4.1 4.8 6.7
##  [667] 4.4 3.2 4.3 5.0 4.5 4.5 5.1 2.6 5.1 4.0 4.4 4.4 5.5 4.7 5.7 4.5 3.1 4.6
##  [685] 4.5 3.4 5.1 5.1 4.1 4.1 5.1 4.2 4.9 4.5 4.0 4.0 5.2 3.6 4.4 3.1 4.6 5.6
##  [703] 4.4 5.9 5.5 4.1 3.0 3.5 4.9 3.8 2.7 3.4 4.1 5.0 5.1 4.0 4.6 5.2 4.3 5.0
##  [721] 5.3 5.6 5.4 5.0 5.8 4.2 5.2 2.6 6.1 3.6 3.3 3.7 4.6 4.0 4.7 4.1 4.8 4.0
##  [739] 4.7 6.3 5.2 3.5 4.6 4.5 4.8 4.7 5.1 3.9 5.9 5.2 3.9 3.4 3.8 5.0 2.8 4.4
##  [757] 2.6 6.2 2.9 4.0 6.0 4.3 4.5 3.4 5.0 3.7 5.5 4.9 3.5 4.1 4.8 5.0 4.3 3.7
##  [775] 4.9 4.9 4.4 6.2 4.3 3.7 5.7 4.8 6.2 4.2 3.4 5.1 5.1 3.5 4.4 5.2 4.7 4.2
##  [793] 5.5 5.2 4.3 4.7 4.7 5.3 4.2 4.0 5.1 5.1 4.6 3.7 4.3 3.7 3.6 4.0 3.3 3.9
##  [811] 3.1 5.0 3.6 5.2 4.4 4.7 4.2 4.4 4.2 3.9 5.7 6.3 4.7 3.5 4.2 2.5 4.2 4.2
##  [829] 5.0 4.5 5.7 4.4 5.0 4.5 4.2 5.5 4.5 4.2 4.1 5.3 3.2 3.5 3.4 3.4 5.2 4.3
##  [847] 5.5 3.3 3.9 4.1 4.6 5.9 4.1 5.2 4.3 4.4 6.5 4.0 5.4 5.7 6.1 3.7 4.5 4.0
##  [865] 3.9 5.0 5.9 2.4 5.3 4.1 3.3 3.9 4.9 3.3 5.3 3.7 5.6 5.1 4.7 4.8 4.6 3.5
##  [883] 3.1 5.3 4.8 3.3 4.4 4.1 4.4 4.2 4.5 5.5 4.9 5.7 5.5 3.8 4.9 6.4 4.9 4.9
##  [901] 5.7 6.1 5.5 3.9 3.9 3.7 4.0 4.8 4.5 4.2 4.8 7.3 4.8 3.8 4.6 4.7 5.4 5.2
##  [919] 3.8 4.6 3.4 3.3 3.9 4.8 4.3 5.3 5.4 6.0 4.3 4.4 4.9 5.6 5.8 4.7 5.8 5.8
##  [937] 4.4 6.1 5.3 4.0 3.6 3.9 4.9 3.6 3.2 6.1 4.5 5.4 2.9 3.3 4.7 4.4 5.1 5.1
##  [955] 3.7 4.9 5.1 4.0 3.4 4.9 4.1 5.3 5.3 4.1 4.5 5.5 4.9 4.8 5.2 3.7 4.7 3.5
##  [973] 5.4 5.9 5.3 5.8 5.2 5.2 3.9 3.4 4.2 3.8 5.7 5.1 2.3 4.4 2.6 3.2 4.9 2.6
##  [991] 5.3 3.0 5.6 3.1 4.8 4.4 4.1 4.6 5.4 4.8
## 
## $func.thetastar
## 72% 
## 5.1 
## 
## $jack.boot.val
##  [1] 5.4 5.4 5.3 5.1 5.2 5.0 4.9 4.7 4.6 4.5
## 
## $jack.boot.se
## [1] 0.9338629
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
## [1] 0.8397673
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
##   5.115702   9.033740 
##  (2.217236) (4.114058)
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
## [1]  0.2032886  0.6515886  0.2340799  0.8174289 -0.3968843 -0.4806364
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
##    [1]  0.8004320458 -1.1271507447 -0.1548992823 -0.5315124225  0.2321403580
##    [6]  0.4077962157  0.4840803810  1.6280954521 -0.0597140721  0.7064598201
##   [11]  0.7163010519  0.6315380755 -0.3983391550 -0.5070841214 -0.3237960190
##   [16]  1.0025270839  0.9446852730  0.9329762875 -0.7754532692  0.6566465127
##   [21]  0.1276649495  0.8298039556 -0.0297412522  0.4237616385  1.1210623204
##   [26]  1.2230552208  1.2824991223  1.2982877612 -0.0512181187  0.0027795753
##   [31]  0.5742914118  1.0161130551 -0.2631244967 -0.4204209526  1.3284533165
##   [36]  0.7289023711  0.9343756275  1.0064894682  0.6681587369  0.1925339110
##   [41]  0.8393328472  1.1805356223  1.1690549946 -0.0218172489  0.8478475441
##   [46]  0.3887709025  0.5416427023  0.5824221278  0.9860493732  0.9560898865
##   [51]  0.6948819278  0.5813779643  0.7016627265  0.6720693502  0.7931261588
##   [56] -0.2371745345  0.3344898116  1.0382821931  0.3160825773  0.5959246185
##   [61] -0.0697582406 -0.3527495283  0.3429212765  0.4567321229  0.6903589338
##   [66]  1.0223045829 -0.2005097760  0.8364478704  0.9208279547  0.2152071898
##   [71] -0.0003605232  1.4215788971  0.2669531692 -0.2922634584  0.0732309619
##   [76]  0.3073387008  0.5488904334  0.6833259062  0.5344500504  1.7245927736
##   [81]  0.9574041890  1.1500484211  1.0282656823 -0.2461500466  0.0898902775
##   [86]  0.6797962212  0.5269902131  0.2283333011  0.6636799679  0.8312234974
##   [91]  1.0119471648  0.0132267326  0.4393405936  1.2192127056  1.2253788650
##   [96] -0.2637046949  0.8834178285  0.6674324808  1.0167735736  0.6929512074
##  [101]  1.4669447696  0.8763639497  0.8561340691  0.9796952594  1.1014319724
##  [106] -0.3680325202  1.0762040007  0.3887709025  0.8012466982  0.9988066050
##  [111]  1.3340204833  0.8464588360 -0.7561838388  0.0687901367  0.3720293684
##  [116] -0.5784906020  0.2967418949  1.3275800023  0.7063305095  0.4050099450
##  [121]  1.1196820271  0.5829738278 -0.8512505269  1.6504186778  1.5210108055
##  [126]  1.9162724179  0.1756718200 -0.2669897136 -0.6472144204  0.1827224893
##  [131]  0.3419953625  1.1178901993 -0.1759871538  2.0478191564  1.1840179472
##  [136]  1.0576671269 -0.4467835103  0.2570879303  1.1538940959  0.7494132946
##  [141]  0.2694551038 -1.0055716891  0.6296890013 -0.4516860204  0.4629273263
##  [146]  0.9686968098 -0.0803935973  1.1997482472  1.2784975924  0.0715842048
##  [151] -0.1365741694  1.7130727809  0.4633613942  0.7129928197 -0.3219956899
##  [156]  0.7080084448  0.0100215524  1.0346640038  0.5916901963 -0.1409856012
##  [161]  0.9980428304  0.1554329299  0.6538200649  0.3283716457 -0.4342272994
##  [166]  0.2370927701  1.1934931167  0.4846781497  1.0440273879  1.1528985211
##  [171]  1.2186203178  0.2504496217 -0.0031741025  0.4338581545  0.5281104149
##  [176]  0.0852841581 -0.1234256693  0.5577381390 -1.2633613406  0.7049154278
##  [181]  0.7284729632  1.6313970266  0.9739780642  0.2127783000  0.4040100318
##  [186]  0.8929433230  0.8526952489  1.1600956773  1.9377819131  0.3538547621
##  [191]  0.9126865203  1.0228602261  0.4380438243  1.7206169333 -0.4263763622
##  [196]  0.8075210409 -1.9519974756  0.2477353157  1.5583486803 -0.6646727979
##  [201]  0.8445905787 -1.4010358974  0.9678919396  0.4819386620  0.7758197118
##  [206]  0.5778409145 -0.5649091914  1.3611266897 -0.4930881061 -0.0338258730
##  [211]  0.7247748874  0.6665445905  0.4237242107 -0.1231294454  0.8378968270
##  [216]  0.7237304985  0.4348628384  0.0671065988  0.0382587062  0.4845298106
##  [221]  0.8388653857  0.3863072708  0.6504772586  1.0209098474  0.1761190837
##  [226]  0.8009260453  0.6406656954  0.2864529183  1.7506588093 -0.4825559119
##  [231]  0.3006013795  0.3200953372  0.0279653396  0.1329339054  0.7269740311
##  [236]  0.8552859529  0.2495531579 -0.4212661297 -0.0785399278  0.3605183403
##  [241]  0.2712895378  0.6028113768  1.0806228525  0.2919301434  0.7368251407
##  [246] -0.6489680936  0.5829738278  0.5051300826  0.1655814746  0.9598416382
##  [251]  1.1477809182  0.3989856080  0.3573315358  0.8890138528  0.5497606780
##  [256]  0.3390392944  0.6627479378  0.6076134431 -0.0551182700 -0.8847367737
##  [261]  0.8493667854  1.6775035588  0.6149034752  0.5880602543  1.1495484951
##  [266]  0.7515956195  0.6375903710  1.4916908396  0.5810932746  0.3571450210
##  [271]  0.5376104269 -0.0883240255  0.4857031932  0.5528628353 -0.0236592721
##  [276]  1.1452675360  1.1524272724  0.5149744474  0.6837204844  0.1761827207
##  [281]  0.9521002590 -0.2114737961  0.2980201816  1.7075968156  0.4649351350
##  [286]  0.4002526209  1.0389472235  0.7855581570 -0.2453448752 -0.1118456424
##  [291]  0.5452916441  0.9509011448  0.1083223667  1.0847523397  1.9741136622
##  [296]  0.0852147417  0.6519719785  0.5594415428  1.1244183600  0.2540134188
##  [301]  1.2317797209  0.9512098046  0.3794274964  1.0919349000 -0.2587235687
##  [306]  0.9248768657 -0.3144758170 -0.2395242245  0.6731421045  1.2877950015
##  [311]  0.4280434814  0.9327864194  0.7624595396  0.3630619839  0.7645619273
##  [316]  1.4533675938 -0.4370632298  0.9688413446  1.0816217452  0.7392849420
##  [321]  1.6981707968 -0.3622151212 -0.0494301322 -0.3124831865 -0.0075424504
##  [326]  1.9498067262 -0.2950799631  0.1670028949  0.7626611071 -0.5272296623
##  [331] -1.0940772670  0.9597915248  0.1307761503  0.6794363877  1.4259695963
##  [336]  0.3434305624  1.1154692979  0.2944013269  0.1216862404  0.7943870805
##  [341]  0.5308411430  0.6369600067  1.7408496046  1.3905912911 -0.0273971010
##  [346]  1.4523714297  0.7046762344  1.2489349265  1.2207027243  1.4602500112
##  [351]  1.4955412145  0.8393328472  0.3913267854  0.7285526358  0.3991856539
##  [356]  1.0027083950  1.0112064503  1.1161935357  1.7660829810 -0.2050099450
##  [361]  1.1089620346 -0.1625369265  0.1965995074  0.3053436126  0.4751195367
##  [366]  1.2283028939  0.9439324204  0.4932290559 -0.6371914441  0.2866510047
##  [371]  0.6089714202  0.1026988048  0.1920070219  0.7760777649  0.2761701618
##  [376]  0.7801942775  0.2708714829  0.7845343047  1.0764438025  1.3351621726
##  [381]  0.3043731309  1.2801728218  0.8885946381  1.3314728981  0.5943441838
##  [386]  1.1039106644  0.1611389224  0.2405839572 -0.6090409654  1.6242696910
##  [391]  0.1271356044 -0.9399624175  0.4935489570 -0.5100220924  1.0550827401
##  [396] -0.0204387470 -0.1773085663  1.0544384226  0.7436888317 -1.0379618704
##  [401]  0.9935215775  0.6894790376  0.0131844229 -0.2243525927  0.9931532398
##  [406]  0.8843112910 -0.9147917585 -0.1648205188  0.9638355292  0.6108903347
##  [411] -0.1905931126 -0.1017807203  1.0839537348  0.5042252943  0.7892584391
##  [416]  0.2402734965  0.9791402194  0.3202699131 -0.0315976580 -0.9360418530
##  [421]  0.5702574806  1.0522093251  0.6137661316  0.1521161246  0.6782100154
##  [426]  0.5824221278  0.4828370635  0.4287378491 -0.6261255919  0.6939078120
##  [431]  1.3593103480  0.9312263759 -0.0284478569  0.4838301000  0.2357096841
##  [436]  1.6001789333  0.4303745894  0.9501207058 -0.2088029675  1.1419373804
##  [441]  0.9937019857  0.7760777649  0.0927353645  0.6230124405  0.5577931837
##  [446]  0.4525818414  0.7920790909  1.3727414989  0.9744463503  0.0098735473
##  [451] -0.4193552296 -0.1468762813  0.3313787860  1.3508571909  0.7196095806
##  [456]  0.3012484594  0.1691042675  1.0058060320  1.1298330115  0.3411942593
##  [461]  1.2025843314  1.1356680179  0.8528853134  0.3312975525  0.5804774891
##  [466]  0.5661882635 -0.4449358589  0.4438861203  0.9024434128 -1.0631534558
##  [471]  0.6198459398  1.4771776915  0.7427732283  0.2737897517  0.9264568011
##  [476]  1.0121877498  1.3185669424  0.8321764292  1.6833859306  1.0554342736
##  [481]  1.0113961927  0.2757788202  0.9444146886  0.9713013953 -0.3804950003
##  [486] -0.5692132512  1.0053579517  0.4405845521  0.9271005423 -0.2305769345
##  [491]  1.1616239838  0.7489296752  0.8477127358 -0.3316690077  0.6193898862
##  [496]  0.1545103811  0.3567191199  0.7747855857  0.1462646710 -0.3781145714
##  [501]  0.4736288856  0.9989852542  1.9491582695  0.9489407212  0.5689155585
##  [506]  0.9428912232 -0.5532385192  0.8312752134  1.6658864965  0.7638263767
##  [511] -0.2243525927  0.0535620410  0.3654748971  0.5102665548  0.8600731185
##  [516] -0.1540342417  0.2683544701  0.0773157451  0.8262207097  0.4707174823
##  [521]  0.2215142120  0.4937372076  0.4455589134  1.3118527711  0.1396061924
##  [526]  0.3284451435  0.0080510886  0.3398424644  0.6782100154  0.4676164265
##  [531]  0.6395755886  0.5431889470  1.3057459616  0.9519942594  0.0011496611
##  [536]  0.1854612064  0.7452600330  0.9217171571  1.3249714358  1.2722244453
##  [541] -0.5460914676  0.3427087047 -0.0871840241  0.4130325807  1.3903034569
##  [546]  0.6039377891  0.8746732045  0.0281560631  0.0910021674  0.9914756367
##  [551]  0.4062710913  0.5460138984  0.2684916458  0.9580826680  0.9201688502
##  [556] -0.7302315059 -0.1149788822  0.2455900174  0.8691733835  0.8217741402
##  [561]  0.4530913217  1.0190610334  0.6366087039  1.1123859320  0.5784622484
##  [566]  0.3836289349  1.5460950280  0.8338798150  0.4769996167  0.4737750898
##  [571] -0.1772375542  0.5603111532  1.1726363490  0.6088827432  0.5454469792
##  [576]  0.5368734849  0.9911269397  0.8318461673  0.7559394525  0.6727861025
##  [581] -0.2711074677  0.0686343142  0.9239047081  0.8199948197  1.1014319724
##  [586]  0.7314629551  0.2714669868  0.3293962910 -0.2763129254  0.4252996947
##  [591] -0.9234970643  0.7796993105  1.4166843200  0.8522002851 -0.0573468380
##  [596]  0.4908076465  1.1856986844  1.0539752266  1.1525192602 -0.3527495283
##  [601]  0.7068191005 -0.5277302124  1.3651103705  0.0909714525  1.1974239224
##  [606]  0.8538576561  0.8145606811  0.7645619273 -0.6463157296  0.0026116363
##  [611]  0.2883228625  0.6873929850  1.3429079141  0.7694021849  0.7674846738
##  [616] -0.0306826158  0.1264734839  1.7762813848  1.4845209914  0.3266395983
##  [621]  0.6978286380 -0.7580347100  0.7372669902  1.0470141081  0.2584802830
##  [626]  1.1206550650  0.0177664178  0.7198241263 -0.2906688504  1.6091746007
##  [631]  1.5736203177  0.4870904013  1.6620602104  1.1329256435 -0.4917168149
##  [636] -0.1375918361  0.8262207097  0.9281425999 -0.0734716359  0.6282711016
##  [641] -0.0325127635  0.5041231677  0.7164455318  1.4560470682  0.7747855857
##  [646]  0.9059621935  0.8818548390  2.0999627265  0.2397167314  1.2142547152
##  [651] -0.0647641886  0.6531578187 -0.0671116623  0.6689291460  0.3110614101
##  [656]  0.0318216742  0.1329173869 -1.1844722694  1.2468808772  0.9639114982
##  [661]  1.3368941939  0.2854405601 -0.3771030385  0.2819788596  0.1281584482
##  [666] -0.1084375168  1.2280822174  0.5084200135  1.5516012406  0.1985120491
##  [671]  1.1496841758  0.5988171747  1.6398935738 -0.2001630388  0.0135114131
##  [676]  0.6385118931  0.3328404657 -0.0003757918 -0.1894179267  1.2938592884
##  [681]  1.0408894533  0.0761982332  1.1982195447  0.2874248193 -0.3881039810
##  [686]  0.2423928831  0.8821679493  1.0898786251 -0.2608520419  0.0379261811
##  [691]  0.8106489831  1.4440821751  0.8094144289 -0.2624199390  0.7825980622
##  [696] -0.4768695649  0.3500169853 -0.1490890257  1.7502689754  0.5168197057
##  [701]  0.0774746193  0.7252686852  0.9593370783  0.1118926898  0.2269434387
##  [706]  0.9006440010  0.7445342362  0.9990202496  0.0359943026  1.3925063123
##  [711] -0.2939529864  0.0665209873  0.6026965758  0.5773523107  1.9491582695
##  [716] -0.6259582797  1.2056600090  0.9708500841  0.7980257833  0.3684799952
##  [721]  1.2676261109  0.3147690776  1.0485772620  0.8164087062  0.2911301463
##  [726]  0.1808452211  1.1997398042  1.0763025449 -0.5621163047 -0.1287433787
##  [731]  0.4889923649  1.1106618442  1.0052414870  0.1864000907  0.7014832534
##  [736]  0.5247898090  0.1269307760 -0.4224866353  0.9726899689  0.6618259432
##  [741]  1.0378642943  0.9498160966  1.6161540140  0.8022220695  0.3208631766
##  [746]  0.3384551966 -0.2257244637  0.9426368396  0.6928126706  0.1864000907
##  [751]  0.5360379887 -0.1585116199  0.6939043985  0.2993914037  1.1816831038
##  [756]  0.6248433813  0.1307321534 -1.0142551235  0.7050159371  0.0104914619
##  [761] -0.2609926778  0.8373430879  0.8308821379 -0.1249613133  1.0832973991
##  [766]  0.1116239377  0.6877630400  1.0411894585  0.4655232330  1.7223783709
##  [771]  0.0610807951  0.7382879030  0.7873837915  0.0633531823  0.6949059357
##  [776]  0.0711215754  0.1920716804  0.6263887957  0.4186207832  0.5766950823
##  [781]  0.9681643001  1.3521474630 -0.9752239007  1.0133701619  0.1000622591
##  [786]  1.5299222384  0.1702510937  0.3856591483  0.7932345479  0.9437803743
##  [791] -0.5709971977  0.2138764817  0.9537691905  1.0470141081  0.2678760106
##  [796] -1.0916497085  0.6370067575  0.8331234811 -0.2123190524  0.1978756553
##  [801]  0.4716618112  0.8713178991  0.2951537751  1.5448838391  0.7095355963
##  [806]  0.7672568878  0.6555854056  0.2490052188  0.8753041445  0.9581640517
##  [811]  0.0007537174 -0.7177064729 -0.0230615202  1.3439793796  0.3039429250
##  [816]  0.9479272259  1.4569866878  0.9465180846  0.5296038112  0.2807742250
##  [821]  0.5813997460  1.6865459423  0.4846917693 -0.0583311287  1.0587537541
##  [826]  0.3879126051 -0.3976954972  1.6820403919  1.8592308451  0.5378322810
##  [831] -0.6299115962 -0.2507188710  0.6725875719 -0.0010712354  0.3604704408
##  [836]  0.4567665847  1.4461428611  0.1564403100  0.7921709318  0.6336174821
##  [841]  0.3094251773 -0.8010297707  1.0190679606  0.6088616041  0.4109233596
##  [846]  0.9905141884 -0.3187139719  0.5847333404  0.5965618288  0.5144080842
##  [851] -0.3034362885  0.5619408527  0.1771188473  1.2288674550  1.2375103483
##  [856]  0.5557675856  0.0009765767  0.9482675044  0.9899992447  1.2360448566
##  [861]  1.0201985271  0.5402130853 -0.4479814407  1.3616662373  0.9326467649
##  [866]  0.9018667631  1.1263225690  1.9839590720  0.8777818878  1.0645465426
##  [871] -0.1861494113  0.9207508403  1.1815318495 -0.2793391732  0.8418709079
##  [876] -1.2600930256  0.7841456991  1.0108287541 -0.0407874410  1.5976138123
##  [881]  0.2806501583  1.0240094556 -0.3090729861  0.4015588566 -0.0894006415
##  [886]  0.5001819888  0.6368637180  0.3460722861  0.0536306332  0.6122374423
##  [891]  0.3737082171 -0.3542068491  0.7530306721  0.1270384189 -0.4944714447
##  [896]  0.2135005843  1.1473695747  0.1082013280  0.8122429207  1.2982877612
##  [901]  0.0376889650  2.1261652982  0.9065172482  0.3288127373  1.0190009501
##  [906]  1.2371614961  1.2235945526  0.9504001477  1.8541583326  0.2148125835
##  [911]  1.3406070222  0.2988339116  0.8590736417  0.0653323410  0.7021880977
##  [916]  0.6741775836  0.4945660149  0.8130799497  0.2310057304  1.0986369962
##  [921] -0.4784583533  1.4773648515  0.6839003509  0.6151538275  0.6886923283
##  [926] -0.0878424735  0.8190255983  0.2692025332  1.1894692005  1.0028812702
##  [931] -0.5460914676 -0.3504219023  0.0161639057  0.9707390932  0.6707117273
##  [936]  0.9760223103  0.4965854442 -0.3159212240  0.1554635158  0.6144025085
##  [941]  0.8295388149 -0.7415595413  1.3376625622  0.4455599327  0.3878994756
##  [946]  1.0296952055  0.5725370205  0.6487352902  0.1309882929  0.4568973394
##  [951]  0.9808427311  0.8482775911  1.3368941939 -0.4824065816  0.1643207569
##  [956]  0.4326829089  0.4024216496  0.7935466711  1.5465811992 -0.9476643942
##  [961]  0.7020790169 -0.8090770879  0.6789494337  0.7374241091  0.6606003427
##  [966]  0.3077527832  0.5146782591  0.4292472374 -0.0430232218  0.5722827168
##  [971]  0.8630380572  1.3350192146  1.4302699215  0.1590429003  0.2621758851
##  [976] -0.6077152806  0.3463416872  0.6116488658  0.8383351249 -0.7184238879
##  [981]  0.1759041735  0.5437195672  0.2982872745  0.6761749110 -0.1790357903
##  [986]  0.7452948697 -0.0899645826  1.2914930146  0.8385831344  0.6576695362
##  [991]  1.0390894396  0.8346227525  1.0534091475  0.1171765795 -0.1206602863
##  [996] -0.1823119109  0.4986468068 -0.5272296623  0.5698991124 -0.1664396495
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

## Warning in densfun(x, parm[1], parm[2], ...): NaNs produced
```

```r
fit2
```

```
##       mean          sd    
##   0.56628537   0.25551852 
##  (0.08080205) (0.05713235)
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
## [1]  0.8034671  0.3461327 -0.3538218  0.4567202 -0.9043123 -0.7085273
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
## [1] -0.0462
```

```r
sd(results2$thetastar)  #the standard deviation of the theta stars is the SE of the statistic (in this case, the mean)
```

```
## [1] 0.8731368
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
## t1*      4.5 0.01971972   0.9282914
```

One of the main advantages to the 'boot' package over the 'bootstrap' package is the nicer formatting of the output.

Going back to our original code, lets see how we could reproduce all of these numbers:


```r
table(sample(x,size=length(x),replace=T))
```

```
## 
## 0 1 2 3 4 5 
## 2 2 1 1 2 2
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
## [1] -0.0074
```

```r
se.boot
```

```
## [1] 0.9049332
```

Why do our numbers not agree exactly with those of the boot package? This is because our estimates of bias and standard error are just estimates, and they carry with them their own uncertainties. That is one of the reasons we might bother doing jackknife-after-bootstrap.

The 'boot' package has a LOT of functionality. If we have time, we will come back to some of these more complex functions later in the semester as we cover topics like regression and glm.

