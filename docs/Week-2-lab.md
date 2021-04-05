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
## 0 1 2 3 4 6 8 9 
## 1 1 1 2 1 1 1 2
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
## [1] 2.781933
```

```r
UL.boot
```

```
## [1] 6.213267
```

Method #2: Simply take the quantiles of the bootstrap statistics


```r
quantile(xmeans,c(0.025,0.975))
```

```
##   2.5%  97.5% 
## 2.8975 6.1025
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
##    [1] 4.0 5.9 3.6 5.8 3.5 5.2 4.8 5.7 3.7 3.4 3.9 3.6 4.0 4.4 5.5 6.0 5.1 3.6
##   [19] 3.3 4.1 5.5 5.6 4.1 5.5 4.3 4.0 4.9 3.0 5.0 5.5 4.0 4.5 4.8 5.0 4.3 4.6
##   [37] 4.2 5.7 4.4 4.0 4.2 5.4 4.6 5.4 3.9 6.1 4.3 5.6 4.6 3.7 3.5 5.3 5.5 3.1
##   [55] 3.3 4.8 6.6 4.2 4.0 3.8 5.5 5.6 6.5 3.2 5.1 4.5 3.6 5.5 4.1 4.2 4.2 3.6
##   [73] 4.1 3.9 5.1 4.4 4.6 4.4 5.5 4.4 4.5 3.8 4.7 3.0 4.8 5.1 3.7 5.9 3.1 4.3
##   [91] 3.0 2.9 3.5 4.0 4.1 4.4 4.2 4.7 4.3 4.8 3.8 3.3 3.3 4.1 3.0 3.8 3.6 2.4
##  [109] 4.4 3.2 4.6 5.1 4.4 4.8 4.4 4.2 4.5 3.6 3.7 6.2 3.9 4.8 4.4 4.6 4.4 3.7
##  [127] 5.6 5.3 4.9 3.9 1.9 5.1 5.8 5.1 4.8 4.7 3.8 4.8 3.9 3.3 3.6 4.6 5.4 3.3
##  [145] 5.2 4.5 4.7 5.3 4.2 3.1 5.0 4.4 4.8 4.0 5.0 5.4 3.6 5.4 2.9 4.3 4.8 5.3
##  [163] 2.6 4.7 3.7 5.5 5.2 3.0 4.9 4.5 5.2 3.6 3.2 3.4 4.6 4.6 4.4 4.9 3.4 4.8
##  [181] 3.8 3.4 4.5 5.1 5.4 4.6 4.2 4.7 3.9 5.7 2.6 5.4 2.8 4.2 4.3 2.8 2.9 7.0
##  [199] 5.4 4.9 2.5 4.5 3.4 4.2 4.9 5.4 4.2 3.2 3.7 5.1 5.8 6.5 5.0 4.5 4.7 5.8
##  [217] 4.3 5.6 5.4 3.1 2.9 3.7 5.1 6.0 5.7 3.8 3.8 4.5 3.4 4.5 3.8 5.6 6.8 4.0
##  [235] 6.4 4.2 5.1 6.5 5.3 4.1 4.5 4.4 2.7 3.9 4.0 4.5 3.5 4.2 4.0 6.0 4.4 4.4
##  [253] 4.8 4.7 4.4 3.9 4.5 5.6 4.2 4.4 4.0 3.9 3.4 4.4 4.3 4.2 3.7 5.6 4.9 4.5
##  [271] 4.7 5.1 4.3 4.9 5.6 5.5 4.5 3.8 6.2 3.1 4.3 4.1 4.4 4.4 6.0 3.7 3.9 4.4
##  [289] 3.6 4.7 4.3 4.3 3.8 4.3 5.0 3.8 6.1 5.2 5.5 4.7 6.3 2.6 6.5 4.9 4.8 3.4
##  [307] 4.2 5.3 6.3 3.9 5.3 3.6 3.0 5.6 4.0 7.0 4.1 4.3 4.9 4.9 5.3 6.6 4.6 3.3
##  [325] 3.8 5.0 5.9 6.5 2.3 6.3 4.6 4.2 3.1 4.0 5.4 3.9 2.3 3.8 5.9 4.4 4.1 4.1
##  [343] 4.5 3.7 4.8 5.5 6.5 5.4 4.0 4.3 5.1 4.3 4.6 5.1 5.1 4.3 5.5 4.8 4.4 4.4
##  [361] 5.1 4.2 5.0 5.1 5.4 4.5 3.2 5.1 5.1 3.9 5.5 4.6 5.1 5.0 3.9 4.9 3.4 4.9
##  [379] 5.1 4.5 3.5 5.0 4.0 4.2 2.7 4.0 3.9 5.1 5.1 4.4 3.2 4.9 4.0 3.5 5.6 3.2
##  [397] 5.4 4.8 4.5 5.1 5.4 4.2 4.3 3.3 4.4 5.1 3.6 4.7 5.9 6.0 4.2 5.6 4.9 3.3
##  [415] 4.4 4.5 4.5 4.8 4.1 3.7 4.6 5.3 3.8 5.2 5.0 5.8 5.2 4.8 5.3 6.1 6.3 5.1
##  [433] 5.1 3.9 5.5 3.6 4.6 6.0 4.5 6.1 4.3 3.8 5.6 4.8 3.4 4.9 4.2 5.1 3.7 3.4
##  [451] 3.8 5.5 4.4 4.2 5.6 4.0 3.9 3.2 3.3 3.4 4.0 4.5 2.5 5.1 2.5 4.9 5.2 5.4
##  [469] 2.6 5.9 4.1 5.7 4.1 5.1 4.9 4.0 3.3 4.5 5.4 5.6 3.6 4.2 3.2 4.3 4.2 3.3
##  [487] 6.3 4.4 4.3 4.9 3.9 5.9 4.4 5.4 3.4 5.3 5.2 3.8 5.1 6.0 3.8 4.0 5.3 4.6
##  [505] 6.0 3.5 2.6 6.1 4.5 5.0 4.6 4.5 2.9 4.7 3.1 3.9 5.0 4.1 4.3 4.9 3.9 3.8
##  [523] 3.9 4.9 2.6 4.8 3.7 3.3 2.2 3.3 2.0 4.4 4.9 5.0 4.8 3.3 5.0 3.5 4.5 3.9
##  [541] 5.0 5.4 5.0 3.0 5.4 4.7 4.1 5.6 4.9 4.0 5.0 4.2 3.0 4.4 5.6 5.1 3.7 3.9
##  [559] 7.1 4.6 3.8 5.0 2.9 4.9 5.6 5.9 3.5 4.5 5.3 4.2 5.3 6.2 5.5 4.7 5.5 4.8
##  [577] 3.7 4.3 5.3 5.4 4.2 3.6 5.0 4.4 4.0 4.4 5.0 4.4 5.2 3.7 4.3 4.5 6.4 5.9
##  [595] 6.1 2.6 5.0 5.4 4.8 5.8 4.1 4.7 2.8 4.1 4.3 5.2 4.3 5.0 4.9 3.6 5.5 3.8
##  [613] 5.2 5.9 4.8 5.0 6.7 4.4 3.9 4.6 5.6 3.8 4.7 3.7 4.7 3.1 2.9 5.3 6.3 3.4
##  [631] 4.8 5.0 5.4 5.7 4.5 3.4 4.9 5.1 5.1 5.5 5.2 6.1 5.4 3.4 2.9 4.0 3.4 4.5
##  [649] 4.3 4.2 3.3 4.8 4.7 4.4 4.6 4.9 2.8 4.5 5.3 4.9 3.7 4.6 3.5 3.4 4.1 4.9
##  [667] 3.5 4.3 4.9 3.9 3.9 4.7 4.5 4.7 4.1 4.6 5.0 4.6 4.8 4.1 4.8 3.5 4.6 4.8
##  [685] 5.9 3.3 2.5 4.6 6.3 4.2 3.8 5.1 4.5 4.2 5.4 3.7 4.5 4.8 4.4 3.9 6.0 4.3
##  [703] 4.0 4.1 3.3 2.5 3.8 5.4 4.1 6.8 3.0 5.5 5.3 4.8 4.0 4.9 4.5 3.4 3.6 4.2
##  [721] 4.2 4.8 4.4 3.5 4.3 6.5 3.5 5.1 6.4 4.9 4.9 6.1 4.9 4.1 4.1 4.4 4.8 4.9
##  [739] 5.1 4.8 4.1 3.7 5.8 5.7 4.4 4.6 5.2 4.5 4.4 4.3 5.1 4.3 4.5 5.9 6.2 4.5
##  [757] 4.4 4.7 2.2 4.5 4.8 4.9 3.2 4.2 3.2 5.0 4.7 4.9 3.7 3.7 5.1 4.7 3.9 5.5
##  [775] 3.5 5.2 4.1 4.5 4.5 4.9 4.6 5.8 4.5 6.1 4.5 3.9 4.8 4.4 2.3 3.5 5.5 4.4
##  [793] 5.1 3.8 2.5 5.4 3.9 6.0 5.8 3.9 3.7 4.1 5.8 4.2 5.1 4.9 5.5 4.6 4.4 4.2
##  [811] 4.7 4.2 3.2 3.9 5.2 3.2 4.8 2.0 4.8 3.9 3.7 5.8 4.5 3.4 2.6 5.3 3.9 4.5
##  [829] 5.9 5.3 5.2 4.6 4.8 5.7 5.7 5.5 3.0 3.9 4.0 5.6 6.1 3.8 4.4 3.4 5.6 2.9
##  [847] 5.0 3.9 3.8 5.1 4.1 4.7 5.3 3.1 3.6 4.7 5.5 4.5 3.7 5.0 5.8 4.5 5.3 4.7
##  [865] 4.3 4.9 5.9 4.4 4.4 4.8 5.4 3.9 3.1 4.2 4.2 5.6 5.1 5.0 5.7 3.8 4.3 4.8
##  [883] 6.4 4.6 5.2 3.1 5.8 5.1 4.8 5.7 3.5 3.7 4.2 6.0 6.9 3.1 3.7 4.5 4.6 4.7
##  [901] 4.9 4.9 4.7 4.9 2.5 4.1 6.2 4.0 3.8 4.2 4.0 3.9 4.5 4.8 4.7 3.9 5.1 4.1
##  [919] 5.3 5.7 5.2 5.0 4.2 4.9 3.6 3.2 4.5 4.9 5.0 4.1 3.8 4.8 4.5 4.9 6.1 3.3
##  [937] 3.1 4.9 5.4 3.9 4.9 4.2 6.0 6.5 3.6 4.4 4.5 4.2 2.7 4.8 4.1 4.1 4.5 4.2
##  [955] 3.3 3.7 4.0 5.5 4.5 5.1 3.6 3.9 4.8 5.9 4.8 4.6 5.3 2.2 3.6 4.9 4.8 4.3
##  [973] 3.4 5.4 5.0 4.8 4.5 2.8 4.8 4.2 4.7 3.0 3.4 4.6 4.0 4.8 4.6 4.0 4.8 4.4
##  [991] 6.3 3.7 2.9 4.2 4.0 5.2 4.5 5.2 4.3 4.0
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
## 2.6975 6.3000
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
##    [1] 4.1 3.1 5.0 5.8 5.2 4.2 3.0 5.2 4.9 4.5 3.9 4.8 6.1 5.0 4.0 3.9 4.0 4.8
##   [19] 2.8 4.4 5.0 5.3 5.0 4.6 3.4 4.4 4.8 3.6 4.2 4.2 4.9 3.4 5.4 3.9 3.9 4.8
##   [37] 4.6 3.3 4.9 4.1 3.7 4.3 3.9 2.8 4.5 4.3 5.7 4.0 6.2 6.0 5.1 3.4 4.9 4.7
##   [55] 6.1 5.2 4.9 3.7 4.7 6.3 6.0 4.5 4.6 4.5 5.4 4.0 4.9 2.9 4.2 4.8 5.4 5.3
##   [73] 5.1 5.3 4.7 4.7 3.9 4.0 3.7 3.5 3.7 3.6 4.9 4.5 4.2 2.2 2.9 3.5 3.1 4.8
##   [91] 6.3 4.2 5.2 5.9 4.8 3.8 1.7 4.1 5.1 4.8 5.5 4.2 4.1 4.7 4.3 5.0 3.7 4.7
##  [109] 6.0 4.5 5.2 4.9 4.3 4.9 3.7 3.5 2.7 3.0 4.8 3.2 3.3 5.4 4.7 4.6 4.7 3.9
##  [127] 4.6 5.1 3.9 4.3 4.5 3.0 3.9 5.5 5.5 6.3 5.1 3.9 4.9 4.3 4.0 4.2 5.2 4.3
##  [145] 4.7 6.1 4.3 4.9 4.5 4.2 3.9 3.7 6.2 3.1 4.7 2.8 4.3 4.6 4.9 5.0 4.4 5.4
##  [163] 2.6 6.8 6.5 5.5 4.0 3.7 4.1 4.9 6.3 4.8 3.1 6.0 5.2 5.6 4.7 5.2 4.3 5.6
##  [181] 5.4 3.7 6.0 4.8 4.6 4.5 5.2 4.6 4.0 4.9 4.0 4.9 4.7 3.8 4.5 4.3 2.4 4.9
##  [199] 6.0 4.7 5.4 4.6 5.1 5.2 6.3 5.6 5.6 3.2 4.2 4.1 3.7 6.3 6.5 5.5 3.2 4.3
##  [217] 3.2 4.4 3.5 4.5 2.9 4.8 5.2 6.1 4.3 6.4 5.4 4.1 5.8 4.5 4.3 4.6 3.5 3.9
##  [235] 3.6 4.8 3.3 4.4 3.4 4.4 4.0 5.9 3.8 4.2 5.2 4.1 4.7 2.4 5.6 6.5 3.6 2.6
##  [253] 3.5 5.0 6.1 4.5 4.3 4.6 3.6 4.3 3.7 4.5 3.7 6.3 5.5 4.3 6.0 4.4 5.7 4.8
##  [271] 3.1 3.3 3.5 3.8 5.3 5.0 3.6 3.1 2.9 4.8 4.0 5.8 3.5 5.8 4.9 5.6 4.6 4.5
##  [289] 5.3 3.9 3.4 3.6 4.1 4.6 4.9 5.3 3.5 3.9 5.8 5.8 5.3 3.9 3.3 4.9 3.7 5.1
##  [307] 5.9 5.8 4.5 4.8 3.1 4.6 4.2 3.9 5.8 5.8 4.9 3.5 4.5 4.1 2.8 5.7 3.3 5.1
##  [325] 4.0 4.5 3.9 4.0 4.4 4.6 3.9 3.7 4.5 4.8 2.9 4.9 3.9 4.2 4.7 4.8 3.9 4.1
##  [343] 5.0 5.0 4.7 4.7 3.9 5.4 4.0 4.6 4.1 4.4 5.2 4.9 4.9 4.0 4.1 4.5 4.0 5.1
##  [361] 4.0 4.4 4.7 3.1 3.5 5.5 5.6 4.0 3.5 5.4 5.8 4.4 5.5 3.5 3.9 4.7 5.2 5.3
##  [379] 4.0 5.7 5.3 4.8 6.3 7.3 5.4 6.5 4.4 4.6 3.8 3.7 5.7 4.8 5.7 3.9 5.5 3.9
##  [397] 3.3 4.3 4.3 4.3 4.2 5.3 4.3 3.6 3.8 3.4 4.7 4.8 4.5 3.7 4.2 3.8 3.6 3.3
##  [415] 2.9 4.2 4.3 5.2 4.1 4.7 3.4 3.9 2.6 5.2 4.1 3.6 6.0 5.6 4.5 5.3 3.1 4.4
##  [433] 4.8 3.2 5.4 3.7 5.5 4.0 4.2 5.1 5.8 6.5 5.8 2.9 3.2 4.3 4.4 4.9 5.5 3.0
##  [451] 5.2 4.2 3.7 2.3 5.1 4.4 4.3 3.9 3.4 3.2 4.3 6.4 3.8 3.8 4.7 5.1 4.6 5.9
##  [469] 5.7 3.4 3.8 4.8 3.0 5.2 4.3 5.5 3.9 5.6 4.0 4.4 5.5 5.8 5.4 3.2 4.4 4.3
##  [487] 5.3 3.8 5.5 2.9 3.3 3.8 4.8 5.2 2.7 5.1 4.3 4.5 3.9 2.0 4.4 5.2 4.0 3.7
##  [505] 5.6 4.1 4.3 3.9 4.0 4.8 4.8 4.5 4.0 4.8 4.7 5.9 4.2 5.5 3.8 4.1 5.8 5.9
##  [523] 3.6 4.5 5.5 3.3 4.5 3.2 3.6 3.3 4.7 3.9 3.1 5.3 4.9 4.9 3.4 6.7 5.1 4.3
##  [541] 4.6 4.2 6.0 5.0 4.0 5.0 3.4 5.3 4.0 5.0 4.1 5.0 3.7 4.0 3.4 6.3 5.1 4.6
##  [559] 4.3 3.1 4.3 4.0 5.9 4.3 3.6 5.8 5.3 5.8 3.9 4.3 4.1 5.0 4.1 3.9 4.7 4.0
##  [577] 5.7 3.6 4.9 5.1 4.9 5.4 3.9 4.1 3.5 3.3 3.8 3.4 4.1 3.3 2.5 6.1 5.4 2.9
##  [595] 3.6 4.8 4.5 3.8 4.9 6.1 4.7 3.0 4.1 5.3 3.4 6.4 4.6 5.7 5.0 4.2 5.6 5.4
##  [613] 2.4 4.5 4.8 4.8 4.9 3.6 4.5 5.2 4.6 4.4 5.5 3.6 4.1 3.6 5.8 4.8 4.5 3.9
##  [631] 5.7 3.5 2.9 5.8 2.7 2.7 6.5 4.9 3.4 5.9 3.5 4.3 4.8 2.3 4.9 4.7 4.5 3.3
##  [649] 3.7 6.4 3.3 3.9 5.4 4.4 6.0 4.3 5.3 2.9 4.7 6.1 5.7 2.4 5.0 5.9 5.2 4.8
##  [667] 3.3 2.8 5.0 4.6 4.1 4.7 4.4 3.8 4.5 3.5 4.8 3.5 5.3 3.9 3.2 3.7 4.5 3.6
##  [685] 5.5 5.3 5.0 5.6 3.9 4.6 4.3 4.7 4.6 4.3 4.6 5.2 5.4 2.6 3.1 3.6 5.3 3.0
##  [703] 3.6 5.6 4.0 5.8 3.6 5.6 4.3 4.9 4.8 4.7 4.5 4.4 5.2 6.0 3.8 6.6 7.0 3.9
##  [721] 4.2 6.6 3.3 5.2 5.0 3.5 3.0 4.3 4.3 5.2 4.5 4.9 5.3 5.6 3.6 3.0 3.5 3.4
##  [739] 5.5 3.5 5.0 5.3 3.1 4.4 4.7 5.6 5.1 4.6 2.7 4.4 6.8 4.3 2.7 5.4 4.5 1.4
##  [757] 5.1 5.0 4.9 5.6 5.2 3.8 4.7 5.2 4.0 4.9 3.9 5.3 4.4 4.2 4.0 4.3 5.5 4.2
##  [775] 4.5 3.8 3.8 3.4 3.9 3.6 4.5 3.8 5.6 5.8 5.0 4.1 3.0 4.2 3.5 3.9 4.8 4.4
##  [793] 4.2 6.5 5.3 6.0 4.6 6.2 3.8 3.3 5.8 2.7 4.5 5.7 4.1 6.9 5.0 5.0 4.1 3.2
##  [811] 3.0 6.2 3.9 2.8 6.3 4.7 6.6 5.4 4.4 4.5 5.0 4.7 5.9 3.6 4.8 4.9 3.8 6.3
##  [829] 5.2 3.0 4.4 5.2 3.9 4.4 4.1 6.2 3.6 6.3 3.9 4.2 5.5 4.8 6.6 3.3 4.6 4.9
##  [847] 3.4 3.7 5.5 4.8 4.2 3.6 2.4 4.3 5.4 4.4 2.6 3.9 3.3 3.3 5.0 2.7 4.6 3.9
##  [865] 3.2 4.3 4.8 4.0 6.0 3.9 4.2 4.1 5.0 4.1 3.7 2.7 4.5 4.0 5.4 3.3 5.2 4.5
##  [883] 4.5 3.9 4.6 4.1 4.4 2.9 5.8 5.6 4.6 4.8 4.9 4.8 5.9 4.4 3.9 4.3 4.1 4.3
##  [901] 4.1 5.8 4.8 4.7 4.7 5.2 4.0 4.2 4.0 4.2 3.8 4.6 2.8 5.4 5.5 6.5 5.4 5.0
##  [919] 5.2 4.3 4.4 3.9 3.6 7.2 2.7 5.5 5.0 3.3 4.0 5.3 4.8 4.3 4.3 5.4 4.8 3.5
##  [937] 4.1 5.3 3.6 4.2 4.3 3.1 4.2 5.0 4.8 4.8 5.4 4.9 4.7 4.6 3.4 3.5 6.5 3.7
##  [955] 5.9 5.1 5.2 5.9 4.8 4.4 5.4 4.5 4.1 4.9 5.1 4.1 5.2 6.0 4.4 6.3 4.7 5.6
##  [973] 3.3 5.0 4.8 4.0 5.4 5.1 4.8 3.7 4.0 3.8 4.2 4.1 4.4 6.5 5.8 5.0 3.5 4.7
##  [991] 5.2 4.6 3.4 4.6 5.3 4.6 2.5 6.3 6.3 4.7
## 
## $func.thetastar
## [1] 0.0064
## 
## $jack.boot.val
##  [1]  0.53274336  0.42112676  0.26873156  0.14893617  0.03229462 -0.01853933
##  [7] -0.15615142 -0.26957746 -0.38159341 -0.51104816
## 
## $jack.boot.se
## [1] 0.9768857
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
##    [1] 3.4 5.8 3.6 5.6 4.4 3.5 5.2 4.7 3.7 4.1 4.4 4.0 3.7 5.5 4.2 7.7 5.0 3.4
##   [19] 4.2 2.6 3.4 5.2 3.6 3.1 4.5 3.1 6.2 6.1 4.3 3.0 5.3 4.2 5.4 3.5 5.5 4.9
##   [37] 3.9 4.7 3.4 4.4 4.7 4.3 4.0 5.0 4.4 4.5 4.0 5.4 6.3 4.5 6.3 3.7 5.2 3.1
##   [55] 4.2 5.4 5.6 4.0 2.9 4.7 2.8 5.6 3.5 4.0 3.7 4.7 5.4 6.1 4.2 2.7 4.7 5.7
##   [73] 4.2 4.7 4.4 4.6 1.5 4.9 5.1 5.2 3.6 5.9 3.0 4.5 3.6 5.2 3.6 6.2 3.5 4.5
##   [91] 3.4 3.0 4.2 3.0 3.9 4.9 2.7 3.9 3.9 4.0 4.2 3.8 3.9 5.4 5.0 3.7 5.4 6.5
##  [109] 4.1 3.4 6.1 4.2 4.5 2.1 4.5 4.3 2.8 3.1 5.0 4.0 5.0 3.1 4.7 5.2 4.4 5.3
##  [127] 3.3 6.2 4.8 3.9 5.3 6.1 4.5 4.1 4.3 4.4 5.1 5.4 5.5 5.5 3.8 3.9 3.0 6.1
##  [145] 3.6 4.6 4.2 5.4 4.0 5.4 4.4 4.1 5.4 4.5 5.1 2.2 4.8 6.0 4.2 4.8 5.2 4.9
##  [163] 5.2 4.0 5.8 5.3 4.7 4.7 3.4 4.2 4.1 4.9 3.5 5.2 3.2 5.7 3.3 5.3 5.9 3.3
##  [181] 3.4 5.1 4.3 4.9 3.8 4.4 4.8 3.0 1.9 3.6 3.7 4.4 4.3 4.1 3.5 3.8 4.7 4.2
##  [199] 3.9 4.7 2.7 5.5 4.6 3.9 5.5 5.4 4.8 4.9 3.1 4.9 5.0 4.6 5.4 4.2 5.0 3.4
##  [217] 4.0 5.0 4.9 4.7 4.2 5.2 5.0 5.9 4.1 3.7 5.8 7.0 5.0 6.6 3.1 4.9 4.2 4.5
##  [235] 3.5 3.2 3.7 5.7 5.7 4.3 5.6 4.9 3.2 4.1 4.4 5.9 4.8 4.2 6.0 4.6 2.3 4.0
##  [253] 4.9 3.5 3.0 3.9 4.4 4.6 5.7 3.8 5.8 4.8 4.1 5.3 3.3 3.6 5.3 5.0 4.4 4.6
##  [271] 3.4 4.2 3.7 3.3 4.3 5.2 4.5 3.1 5.1 4.8 3.3 3.9 5.1 4.8 4.0 4.5 4.9 5.4
##  [289] 5.3 3.7 4.2 5.4 5.4 5.5 5.5 4.7 4.7 4.8 4.7 3.9 5.8 4.1 5.1 2.6 6.5 4.2
##  [307] 5.5 3.7 4.9 4.5 5.5 5.3 4.2 4.0 2.9 2.9 5.7 5.3 6.4 5.0 4.3 5.0 4.8 3.4
##  [325] 3.8 3.5 6.1 3.6 5.6 4.0 5.1 5.7 4.3 4.0 2.8 3.9 5.8 4.7 5.5 3.3 3.1 4.6
##  [343] 4.5 3.2 6.0 4.6 5.4 4.5 5.4 3.9 4.6 6.0 3.5 4.9 3.8 4.4 4.9 4.6 4.7 5.6
##  [361] 4.4 4.6 4.2 4.7 3.6 5.4 5.0 3.7 3.0 4.2 4.6 4.3 3.9 2.2 3.8 4.0 4.0 5.0
##  [379] 3.8 4.0 2.8 5.1 5.8 5.9 4.0 3.3 3.6 4.6 4.1 2.7 5.1 5.8 3.6 4.4 5.0 4.0
##  [397] 5.0 4.2 3.9 5.7 4.5 4.5 5.4 3.7 4.3 4.3 2.6 3.4 4.8 3.8 4.4 5.9 4.2 2.6
##  [415] 4.7 5.1 5.2 4.8 3.1 4.2 4.6 3.0 5.0 3.8 4.4 3.8 4.0 4.2 4.8 5.0 5.1 5.5
##  [433] 5.0 3.7 4.2 5.2 3.7 4.6 4.1 3.2 4.7 4.1 4.5 4.8 6.2 4.4 3.1 4.3 4.0 5.0
##  [451] 5.2 4.2 5.7 3.9 3.6 4.2 4.5 4.2 4.4 3.0 4.8 5.0 4.4 5.9 4.5 4.1 4.9 4.5
##  [469] 5.5 4.2 4.2 6.3 2.3 3.2 4.8 5.1 4.3 4.1 3.0 3.5 5.5 4.8 3.9 4.0 5.4 3.8
##  [487] 5.1 6.5 3.9 4.7 4.0 4.2 4.6 5.0 3.7 5.5 5.5 4.9 3.9 4.2 3.8 3.5 4.6 6.3
##  [505] 3.3 5.5 4.7 5.2 4.4 4.7 4.7 3.7 5.9 4.2 4.8 5.9 5.0 5.9 4.9 5.8 2.8 5.5
##  [523] 5.6 5.6 5.1 4.6 4.1 2.8 5.2 4.1 4.2 4.3 4.5 3.8 5.5 4.7 3.3 3.8 4.4 3.8
##  [541] 5.2 5.0 5.1 4.6 4.3 5.5 6.5 4.2 4.6 5.0 5.8 4.5 4.7 3.8 4.1 3.8 3.0 5.1
##  [559] 5.2 4.1 3.8 4.5 3.8 5.2 6.1 4.1 6.4 3.3 2.4 3.6 4.4 5.5 4.6 5.8 4.3 5.2
##  [577] 4.0 6.1 2.6 5.4 4.0 5.7 5.0 5.6 4.3 5.4 4.4 5.3 4.8 4.6 4.5 6.2 4.2 4.2
##  [595] 3.6 3.9 4.9 4.2 4.9 3.7 2.9 5.4 2.7 4.6 5.5 3.1 4.3 3.3 3.6 4.2 4.9 4.2
##  [613] 4.7 3.3 3.9 3.6 4.1 3.6 4.4 3.3 4.8 4.9 5.4 4.4 3.6 5.5 3.9 4.9 4.9 2.7
##  [631] 5.9 5.4 5.8 4.6 4.2 4.9 3.6 3.6 4.8 5.6 4.5 3.7 3.7 3.9 4.5 4.8 3.9 3.2
##  [649] 5.1 5.1 3.4 4.5 3.1 4.6 3.4 6.7 4.8 4.4 4.8 4.7 5.9 3.6 2.5 3.8 5.1 4.2
##  [667] 4.3 5.9 4.0 4.1 5.3 6.3 5.1 4.8 4.7 4.0 6.3 5.2 3.4 5.9 4.7 5.5 4.0 3.6
##  [685] 4.8 5.4 5.4 4.8 3.4 3.8 4.2 3.0 5.3 6.0 3.4 3.5 5.1 4.4 4.6 5.2 4.9 5.6
##  [703] 4.7 4.3 4.3 3.1 4.3 4.6 2.6 3.5 4.8 5.1 5.2 2.6 4.2 4.5 4.7 3.9 2.8 4.1
##  [721] 3.9 5.1 4.4 4.0 4.9 3.6 4.6 4.4 5.7 3.6 4.2 4.3 3.9 3.7 4.5 5.5 5.6 4.6
##  [739] 3.9 3.1 4.6 4.5 5.9 5.0 3.3 4.9 3.1 4.6 4.9 3.7 4.6 3.3 4.4 4.5 4.2 5.3
##  [757] 4.2 4.8 5.4 4.4 3.4 4.2 5.0 3.4 3.3 6.2 4.9 4.6 6.7 4.8 5.5 5.5 4.6 5.3
##  [775] 4.6 5.4 4.6 5.0 3.3 5.1 2.9 2.5 3.6 4.4 4.7 4.7 3.9 4.8 5.7 4.1 5.1 4.6
##  [793] 5.0 3.3 4.3 4.2 3.4 4.2 5.0 3.8 5.0 5.1 4.6 3.6 5.4 5.6 3.8 3.7 5.0 3.8
##  [811] 7.1 3.5 6.9 5.2 4.9 4.2 4.6 3.9 4.1 6.0 4.4 6.0 3.0 3.7 4.6 4.8 6.1 4.7
##  [829] 4.8 4.9 5.5 4.5 4.6 5.0 4.2 4.7 5.7 4.4 4.3 4.2 4.6 5.4 3.3 4.0 4.0 5.0
##  [847] 4.4 4.0 5.3 5.5 3.4 3.7 5.4 4.7 3.6 4.2 3.7 5.3 5.5 3.0 3.3 3.6 5.6 3.8
##  [865] 3.1 6.7 3.5 5.5 3.3 4.2 4.0 3.5 5.8 4.5 4.4 4.3 4.4 4.4 5.8 4.8 5.7 3.8
##  [883] 4.0 3.5 3.9 3.1 3.2 7.3 4.2 6.3 4.7 3.9 4.6 3.3 5.5 3.9 4.7 4.7 5.7 4.8
##  [901] 3.3 6.5 4.9 4.5 4.2 5.1 4.4 4.2 2.9 5.7 3.4 4.7 4.0 6.8 4.8 6.1 5.3 4.5
##  [919] 5.8 4.5 4.2 4.4 4.3 5.2 5.0 2.6 4.6 3.7 3.5 6.0 5.2 4.2 3.5 5.3 4.9 6.0
##  [937] 4.8 4.8 5.1 3.9 2.7 4.6 2.1 5.5 4.3 3.0 4.2 3.1 3.8 4.9 4.8 3.3 4.0 4.7
##  [955] 3.8 3.7 4.0 5.7 2.2 4.1 4.6 4.2 4.3 3.6 3.7 4.1 4.3 4.0 5.2 5.4 4.2 5.7
##  [973] 5.7 5.5 5.1 3.7 4.9 5.0 5.2 6.3 2.4 5.3 4.1 2.1 4.5 4.1 5.4 4.3 4.7 4.1
##  [991] 4.6 5.4 6.2 4.2 3.9 5.7 4.0 3.7 5.0 6.0
## 
## $func.thetastar
## 72% 
##   5 
## 
## $jack.boot.val
##  [1] 5.500 5.400 5.400 5.200 5.100 4.932 4.828 4.800 4.600 4.552
## 
## $jack.boot.se
## [1] 0.9698762
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
## [1] 0.5818867
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
##   10.860544   18.597582 
##  ( 4.784240) ( 8.384637)
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
## [1]  0.1485816  1.0182068  1.2314207  1.6080656 -0.3431003 -0.6646901
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
##    [1]  0.9566974300  0.3727930811  0.0460015212  0.9477503779  0.0513179753
##    [6]  0.7249526525  0.9139392114  0.4097900537  0.2247556978  0.4201032054
##   [11]  0.1463539843 -0.1404799225  0.3832759663  1.5461516100  0.1514535604
##   [16]  0.9225753846 -0.1583125693  0.9601359244  0.3985680419 -0.2846143244
##   [21]  0.6582965617  0.2229657206  0.6679728657  0.3351110244  0.4051963632
##   [26]  0.5151910190 -0.0795150476  0.0725603026  0.7885494854 -0.3736543712
##   [31]  0.7032826498  0.3355025999  1.2926811061 -0.1249775437  0.9789060630
##   [36]  1.0782634076  0.1984494480 -0.2377996292  0.0261831917  0.4789478277
##   [41]  0.4580503442  0.2113886146  0.2955693053  0.1365094262  0.6571652765
##   [46]  0.7386996901  0.5995885386  1.2562463933  1.2722916852  0.8660232519
##   [51]  0.8148847325  0.4557171949  0.2829986779  0.5180865556  0.4507590514
##   [56]  0.1080393458  0.5281438263  1.1059599657  0.6406647825  0.4439299485
##   [61]  0.7074224877  0.1449658111  0.5556892314  1.1256637839  0.5128821915
##   [66]  0.2275714342  0.9831882264  0.9886560281  0.0456083605  1.1490388208
##   [71]  1.3696289255  1.1195695520 -0.3735101369 -0.1263757572 -0.2532503498
##   [76]  0.7296844647  0.6370088049  0.1492118868  0.4357351388  0.5213980221
##   [81]  0.2324052155  0.7556262531  0.5152200270  0.3495061904  0.1864256824
##   [86]  0.8777955252 -0.0040315711  0.8314338495  0.4923746413  0.3023305874
##   [91]  0.3237860241  0.3729596032  0.5678340406 -0.3262698345  0.4276398304
##   [96]  0.7289177064  0.1742460384 -0.0685748886 -1.1156753287 -0.0539103355
##  [101]  0.1597279418  1.5510995001  0.9512222115  0.1264587283  0.5177639747
##  [106]  0.7786030656 -0.1621186852  1.0210692400  0.8699512279  0.9341532212
##  [111]  0.7667700662  0.4013656373  0.2637595521 -0.4040703349  0.9755355907
##  [116] -0.3328318869  0.7168011981  1.4825213711  0.2737748948  0.0493213893
##  [121]  0.2248771238  0.0299078887  1.2121370453 -0.3925034718  1.5232375727
##  [126]  0.8864868513  0.7873350850  0.3663314166  0.7648984151  0.7223780724
##  [131]  0.6164312953  1.1951079140  0.7359591447  0.9290055260  0.3483082971
##  [136]  0.1456089033  0.7733254867  0.9640397938  0.3502089857  0.0232568409
##  [141]  0.4874971384  1.6500201461 -0.0649744350  0.3705183860  0.3881050519
##  [146]  0.2946486453  0.5364837441  0.3075863926  0.6382386540  0.5208655459
##  [151]  0.2906308057  0.5651996952  0.4159929842  1.2538932064  0.6712214948
##  [156]  0.6406586159 -0.0462027618  0.6753276305  0.2532222857  1.0819716583
##  [161]  0.8862406169  0.5194514128  0.4252124685  0.0333252606  0.2155184131
##  [166]  0.3539759593  0.0137716502  0.1792164340  0.7898753005  0.6625097393
##  [171]  0.2211077725  0.1416928625  0.1418906224  0.2116364217  0.1999530624
##  [176]  0.1613488534  1.2462423555  0.2663206283  0.9815986188  0.4893487697
##  [181]  0.2967413349  0.3387672958 -0.0606387224  0.9750060244  0.2307339761
##  [186]  0.9665548535 -0.6601618733  0.2859654356  0.3295850070  0.8127966757
##  [191]  0.7329236714  0.4856153909  0.8344643816 -0.1341290209  0.8095411540
##  [196]  0.9205076935  0.4571287424  0.6544928989  0.1365325185  0.5123773748
##  [201]  1.2881145190  0.9536550208  0.5977898537  0.6849362775  0.2811913175
##  [206]  0.3258950668  0.3507742709  0.3146910526 -0.0224218119  0.4084873991
##  [211]  1.6649326130  0.9455207926  1.0019155384 -0.2508842736  0.2012470762
##  [216]  0.5067281566  0.6575463964  0.5486064420  0.4232638029  0.5195637518
##  [221]  0.5843127422  0.3826945756  0.5013956527  0.1948911193  0.6025688453
##  [226]  1.1602847856  0.9311638034  0.2235768759  0.0404760646  0.9438301774
##  [231] -0.4027766337 -0.2054447004  1.3212262362  0.6687939440  0.3092759591
##  [236] -0.1638595521  0.5060875792  0.6943368165  0.3595530400  0.3673989206
##  [241]  0.2672112454  0.7441399409  0.3745060655  0.6811616506  0.4126920390
##  [246] -0.5505044245  0.4541170418  0.5326005268 -0.6684704379  0.5737569734
##  [251]  0.5462843216 -0.2747522907  0.8341483129 -0.8246439420  0.3998443161
##  [256] -0.1286461392  0.0850854648  0.7620275822 -0.7212092567  1.0908567293
##  [261]  0.9221593660  0.8253219924  0.2139111092  1.0214978636  0.2824532878
##  [266] -0.0754068575  0.2455882081  0.1033157135  0.9785939082  0.8035775905
##  [271]  1.9240671575  0.2137125658  0.5487417678  0.3838808330  0.1474151823
##  [276]  0.2376274597  0.7209240847  0.1058445358  0.4077547089  0.4026769316
##  [281] -0.1052672175  0.8044043730  0.7084523567  0.0823781887  0.0455256548
##  [286]  0.3283486236  0.3120490799  0.7143966294  2.1478783806  0.5661179287
##  [291]  0.7938902460 -0.2835693994  0.7471901620 -0.2883086265  0.5793304448
##  [296]  0.7701053219  1.3923641025  0.2990814990  0.0753531674  0.4304835481
##  [301]  0.8436702287  0.0317583821  0.8129211914  0.3467281412  0.1446567336
##  [306]  0.7405092312  0.6690632995  0.2362119465  0.7307110663 -0.1370178945
##  [311]  0.4596355150  0.9167936767 -0.3030565317  0.9380544308  0.8578793853
##  [316]  1.2283053927  0.0834503077  1.0186360845 -0.3701084969 -0.8359137902
##  [321]  0.1683146255  0.7687544597  0.3454229396 -0.0303218606 -0.8653916046
##  [326]  0.5493543233  0.7181040647  0.8657074056  0.1499989940  0.6855165210
##  [331] -0.2321947902 -0.3237880520  0.6339189815  0.2404111210  0.5710034144
##  [336]  0.3896905924  0.6133605349  1.4054869028  1.0889137830  0.1688112871
##  [341]  0.1102365332  0.8235272642  0.4624250097  1.0916001467  0.2009172780
##  [346]  0.6197133989  0.7735076108  0.5373876372  0.2197612738 -0.8653916046
##  [351]  0.1068011404  0.6646152024  0.4853587965  0.4119581290  0.6765052232
##  [356]  0.0754269553  0.2817223373 -0.7628768344 -0.5895749154  1.2914806309
##  [361]  0.4656749661  0.7326133664  0.9816190192  0.4913261305 -0.6921774354
##  [366]  0.7173269535  0.8555085832  0.6340620275  0.2833546296  1.1086147594
##  [371]  0.1818396922  0.5088733843  0.3486156316  0.0972170236  0.4962604543
##  [376]  0.2666981401  0.7873060257  0.8891496848 -1.4399931728  0.3307564621
##  [381]  1.4227486055 -0.8624871142  0.7613423124 -0.0388106155  0.7065392336
##  [386]  0.2699157087  0.3196870689  0.9811491816  0.6199092115  0.4943638330
##  [391]  0.9975937244  0.8600465541  1.0306262016  0.5097009280 -0.0618068718
##  [396]  0.9011915193  0.2462622872  0.0297040154  0.5891616962  0.2777785422
##  [401]  2.0050613213  0.1160506424 -0.3170200643  0.7358340553  0.4391254656
##  [406]  0.7425005432 -0.8331444918  0.5599467821  0.4837429712 -0.0594215022
##  [411]  0.8386697942 -0.0304197743  1.0368664904  1.1438123018  0.0565021309
##  [416] -0.0011485342  1.0410364346  0.5964615784  0.2524323064  1.4408851804
##  [421]  0.0493689700  0.7563898546  1.0521814675  0.6509513256  0.1570801639
##  [426]  0.6809612411 -0.3109140669  0.7711579884  0.2137125658  0.4694313908
##  [431]  0.5184790276  0.3782404049  0.3387803756  0.1791917981  0.3418314956
##  [436] -0.7841944103 -0.1061683141  1.0232937635  0.6042350862  1.3136092739
##  [441]  0.8958217802  0.4827024081  0.3025597030  0.4106930130  0.7106096738
##  [446]  1.3996601490 -0.7532040246  0.1127280336  0.3664641097  0.6705335732
##  [451]  0.1069801811 -0.5747630477  0.4500317844  0.6205672104  0.5905753088
##  [456]  0.5339039262  1.4822487257 -0.0033836478  0.6302863561 -0.7415026991
##  [461] -0.0115996300  0.7098962988  0.0420069716  0.9654519749 -0.0860223199
##  [466]  0.4742360116 -0.1168020025  0.8386489511  0.9324476404  1.0008421427
##  [471]  0.5427759828  0.3356065277  0.7802805641  0.9292823922 -0.5743001230
##  [476]  0.5533985378  0.6237248937  0.5075918106  0.8766655964  0.0151284668
##  [481]  1.1196436458  0.3234560554  0.8570685042  0.1033741579  0.6416466430
##  [486]  0.6335977837  0.7016424275  0.7239160208  1.0908567293  0.4517891235
##  [491]  0.7376092719  0.9820519181  0.6003071764  0.0776429494  0.6760711700
##  [496]  0.8596598679  0.9645317912  0.2491969815  1.0103821518  0.9353021200
##  [501]  0.3913239311  0.4433164252  0.3935208109  0.9216853493 -0.2872670921
##  [506]  0.0408441857  0.5156991117  0.6138588160  0.1160213192  0.7728966265
##  [511]  0.7695589425 -0.0044359363  0.2854773551  0.1392131260  0.5442279484
##  [516]  0.7931437536 -0.0011313024 -0.1152471874  0.6697308647  0.3544727434
##  [521]  0.4647579579  0.4693487918  0.4367856311  0.9531710190  0.0886076652
##  [526]  0.3529929538  0.2853130662  0.6313416660  0.6860977033  0.5522889873
##  [531]  0.5889947550  1.6276052998  0.9420824955  0.5373260757  0.3601423677
##  [536]  0.9495221918  0.8171845597 -0.0873306713 -0.8241960476  0.2745151157
##  [541]  0.2363145781 -0.1422587534  0.9265581612  0.8684597676  0.2476951657
##  [546]  0.2531155389  0.4144010421  0.8264182414  0.7337959728  0.4872086983
##  [551]  0.4068912674  0.7475802675  0.3858952971  0.8038558151  0.4751258818
##  [556]  0.1252104717  0.0614294697 -0.0982169205 -0.1795581874  0.6625308229
##  [561]  0.4693487918  0.2817768681  0.6804209811  0.1456089033  0.7263976459
##  [566] -0.3422796051  0.8572037066  0.5524853854  0.7928885845  0.4727296832
##  [571] -0.5442373540  0.3700314073 -0.8246231316  0.6061619169  0.5820252877
##  [576]  0.4970790881  0.6923397495  0.9992172758  0.7089009043  0.9536787960
##  [581]  0.9399821728  0.2386603811  0.5221231291  0.0158078656 -0.1333245167
##  [586]  0.3216612348 -0.5111734445  0.0662178033  0.4870246921  0.5963434578
##  [591]  0.8532162334  0.4399611127  0.4937045924  0.3169396835  0.1312985099
##  [596]  1.1183351349 -0.1460828034  0.1211454922 -0.3637736720  0.7485270540
##  [601]  0.0620457125  0.9898247206  0.2445567682  0.5717502404  0.4843393122
##  [606]  0.6638421910  0.5566071562  0.1057499666  0.1170258370  0.0606633931
##  [611]  0.5312185609  0.4856177957  0.4545054755  0.6322319446  0.3950564849
##  [616]  0.4836634257  0.7322003319  0.7422118890  0.4874903813  0.4067937895
##  [621]  0.2332752437  0.7130212617  0.0488817976  0.3352963523 -0.1323229480
##  [626]  0.9999058994  0.0172965539  0.2196025682  0.2540261808  0.3463831400
##  [631]  1.2366589217  0.0845321038  1.0029285101  0.4529214674  0.7083937450
##  [636]  0.1133198009  0.5893286377  0.4326164940  0.7656086069 -0.0087886095
##  [641] -0.1571595891  1.3033916583 -0.0511722029  0.3511910344  0.1392131260
##  [646]  0.0951503678  1.3386661623  1.8516945241  0.4905110344  0.8924930779
##  [651]  1.1268654390 -0.3883161923  0.2763964014 -0.6759047807  1.5737576572
##  [656]  0.9963970910  0.0907306123  0.3951440798  0.8396234204  1.2914806309
##  [661] -0.1503076316  0.5179811986  0.3164276687  0.3295526808 -0.2253131759
##  [666] -0.0109190270  0.0030046101  0.6071977104  0.5171042534  0.8788647637
##  [671] -0.0006594056  0.2523558306  0.6193765351  0.1738764467 -0.0753814987
##  [676]  0.6356864032 -0.0200038639  1.1887645421  0.3884442356  0.0927110495
##  [681]  0.4885588057 -0.3472039362  0.2406384572  0.5499357176  0.4840796365
##  [686] -0.5829450527  0.4760121326  0.6643354685 -0.0500090761  0.2803674627
##  [691]  0.3277478598 -1.4396524244  0.0542709753  0.7511318845  0.6436269889
##  [696]  0.5806440700  0.3767850939  0.6127735243  0.9717383122 -0.3402612330
##  [701]  0.3590501133  1.5892800550  0.4439297707  0.4178924489 -0.0351454645
##  [706] -0.1035149425  0.2040958643  0.4639047625  0.4241807355  0.8709583805
##  [711]  0.1471896868 -0.4112860613  0.1784263916  0.1969163977  1.1243678103
##  [716]  0.2185970273  0.1266592605  0.4351997542  0.2952884186  0.5314586035
##  [721] -0.0717010237  0.6072217954  0.6215112937  1.0329462562  0.2506956293
##  [726]  0.1421777737  1.3636077241  1.7759421256  0.4103579745  0.8655610660
##  [731]  0.1723286052  1.1573415754  0.4978420289 -0.1496593944 -0.2924496523
##  [736]  0.8793364150  0.4763740335  0.3544328561  0.6560668022  1.0357964970
##  [741]  0.4098865673  1.0482173098 -1.0938104831  1.4916407416  0.8675367839
##  [746] -0.5829450527  1.0706789522  0.1324949602  0.0231379293  1.2099598140
##  [751]  0.5139781346  0.5013651896  0.9240330110  0.3394543917  0.4176254964
##  [756] -0.0528345835  0.6724402847 -0.6568086867  0.6364313309  0.2604835355
##  [761]  1.8440134968  0.1255189465  1.0567887636  0.5999343296  0.2893348798
##  [766]  0.9045656572  0.3044332527  0.5120006125  0.5527674179  0.0588218322
##  [771]  0.6296167516  0.5899083816  0.8129290076  0.2481947419  0.6804613194
##  [776]  0.2299524747  0.4856463114  0.7547495804  0.3862267516  0.6021223572
##  [781]  0.6197735473 -0.8246507810  0.5620364435  0.9047850857 -0.3913122694
##  [786]  0.3610595672  0.5142332324  0.5295429637 -0.1505872957  0.8573888838
##  [791]  0.7420346087  1.3310384173  0.1815015158  0.3605365034  0.8740084258
##  [796]  0.3950950688  0.6805874912  0.4199608019  0.9754824525  0.9351364845
##  [801]  0.4199306742  0.3449443145  0.7225166501  0.6877551312 -0.0136647952
##  [806]  0.8666305109  1.2668191449  0.9266231783  1.0259430103 -0.4677439946
##  [811]  1.7440497830  0.4866956856  0.0773782880  0.5416222114  1.5330411627
##  [816]  0.7981473109  0.6530917735  0.0150037057  0.0837240727  0.7584509879
##  [821]  0.6369835896  0.2088220761  0.4319585269  1.0720825947  0.2119182016
##  [826]  0.5893961255  0.9762180490  0.4120046564 -0.5185790349  0.3780855808
##  [831]  0.7244495527  0.3036203598  0.2560794481  0.3936525515 -0.0383316476
##  [836] -0.0457993172  0.2960427542  0.6561356045  0.7453606456  0.3758920643
##  [841]  1.2281279381  0.5664329732  0.4021818333  0.3104907173  0.5098130439
##  [846]  0.1059195975 -0.0775263994  1.3385704609 -0.1261408931  0.0102274663
##  [851]  0.8449060107  0.9047871814  0.4073689530  0.2198364597 -0.3021948895
##  [856]  0.4989136932  1.0577193085  1.0277332225 -0.5951170697  0.5856583596
##  [861]  0.7891282260 -0.2712030127  1.0902714503  0.3758663496  0.8618478188
##  [866]  0.2524694612  0.3777148923 -0.1797834892  1.7155332002  1.0201997311
##  [871]  0.1579195053  0.0046311066  1.0010517406  0.4171992945  0.2954865299
##  [876]  0.5475128290  0.5165613128  0.3575353235  0.4992512113  1.0374628578
##  [881]  0.0344460315  0.7838970807  0.1736379015 -0.1420879442  0.5552811543
##  [886]  0.7830786959  0.5851801516  0.2950334272  1.1754247725  0.4168084923
##  [891]  0.4439297707 -0.1330940620  0.1202331341 -0.3296797937  0.9245607848
##  [896]  0.4754513582  0.2326279610  0.4066190216  0.5049514402  0.7257162930
##  [901]  0.5863151443  0.7461215118  0.6656626618  1.1197877272  0.7453319196
##  [906] -0.0094971456  0.0592968922  0.8391874569  0.4734529425  0.6901711011
##  [911]  0.8341651389  0.7976799774  0.0486003592  1.8472134405 -0.0653719718
##  [916]  0.5678013942  0.3058789512  0.2972337689  0.6943092554  0.7840463257
##  [921]  1.2538932064  0.3333902153 -0.1261081346  0.4241191066  0.3056088273
##  [926]  0.6306751558  0.6124219982 -0.0279739312 -0.2939325502  0.9871843651
##  [931] -0.5638829455  0.1761988407  0.4061152767  0.3427740847  0.5776094226
##  [936] -0.2620604479  0.8507553488  0.2047821140  1.5817124125  0.3782404049
##  [941]  0.5196387891  0.9828683997  0.2734289188  0.6358243392 -0.5331113018
##  [946]  0.2575817083  0.4785368026  0.0465094051  0.3155780016  0.9978835395
##  [951]  1.0022644458  0.4741228580 -0.2905407467  0.3133714718  0.9149537979
##  [956]  0.3428272976  0.6274762860  0.2204782837  1.0471234509  0.7200747177
##  [961]  0.3998708493  1.0241468833  0.5566423914  0.0307842003  0.3960704912
##  [966]  0.3792415755  0.7177267649  0.3078045847  0.9518220675  0.8525307825
##  [971]  0.9553291590  0.1875870497  0.0507781520  0.2487788324  0.0898053423
##  [976]  0.8122594334  0.2766189143  0.6242967033  0.5900970603  0.2614707811
##  [981]  0.3982964503  0.6703273420  0.6321154326 -1.0341478335  0.5079916172
##  [986]  1.0192259104  0.9119606631  0.0442495659  0.8350466289  0.3269377040
##  [991]  0.2009394806  0.3657739422  0.1283230701  0.0628053700  0.4078560367
##  [996]  0.2874338509  0.1621152388 -0.0764526551  0.7280681194  0.2424480547
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
##   0.58396293   0.18099013 
##  (0.05723410) (0.04046648)
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
## [1] -0.277003097 -0.488992274  0.537874083  0.002969472 -0.168089511
## [6] -0.614349752
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
## [1] 0.0235
```

```r
sd(results2$thetastar)  #the standard deviation of the theta stars is the SE of the statistic (in this case, the mean)
```

```
## [1] 0.9124419
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
## t1*      4.5 0.02042042   0.9012924
```

One of the main advantages to the 'boot' package over the 'bootstrap' package is the nicer formatting of the output.

Going back to our original code, lets see how we could reproduce all of these numbers:


```r
table(sample(x,size=length(x),replace=T))
```

```
## 
## 0 1 4 6 7 8 9 
## 1 1 1 1 2 3 1
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
## [1] 7e-04
```

```r
se.boot
```

```
## [1] 0.9181571
```

Why do our numbers not agree exactly with those of the boot package? This is because our estimates of bias and standard error are just estimates, and they carry with them their own uncertainties. That is one of the reasons we might bother doing jackknife-after-bootstrap.

The 'boot' package has a LOT of functionality. If we have time, we will come back to some of these more complex functions later in the semester as we cover topics like regression and glm.

