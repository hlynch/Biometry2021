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
## 0 1 4 5 6 7 
## 1 1 3 1 1 3
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
## [1] -0.0502
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
## [1] 2.638317
```

```r
UL.boot
```

```
## [1] 6.261283
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
##    [1] 5.9 6.4 3.6 5.0 2.9 4.5 3.5 6.1 3.7 4.8 5.6 4.4 4.1 5.5 4.9 4.8 4.1 3.4
##   [19] 4.1 3.6 3.2 5.0 2.9 3.9 5.2 4.0 3.3 3.0 5.1 6.1 5.6 5.9 5.4 3.1 4.3 4.1
##   [37] 3.2 7.6 3.2 5.0 4.1 3.6 4.1 3.7 3.8 4.6 6.3 4.1 2.4 6.2 4.3 4.1 2.8 4.5
##   [55] 3.8 4.3 4.7 6.2 4.1 5.4 3.9 4.3 4.1 5.0 4.7 4.7 3.6 4.3 4.2 4.6 4.0 4.2
##   [73] 3.6 3.5 3.0 3.8 4.6 5.2 4.8 5.2 4.5 5.7 3.7 5.5 4.9 5.0 4.0 5.2 4.7 3.7
##   [91] 5.1 4.4 4.3 4.6 4.1 3.3 3.7 4.6 5.2 4.3 5.0 5.3 3.7 4.2 5.1 3.8 3.2 4.9
##  [109] 5.7 3.4 3.9 5.3 3.1 6.9 5.2 3.7 4.8 4.0 4.8 2.3 2.8 5.6 4.7 3.9 4.4 5.1
##  [127] 4.1 3.3 3.2 3.6 4.6 5.0 5.4 3.7 5.2 3.4 5.3 4.8 4.3 3.6 5.4 3.9 4.2 3.9
##  [145] 4.6 4.3 4.7 3.7 6.0 5.1 4.6 4.1 4.0 3.8 5.6 4.5 3.9 4.1 6.2 3.9 5.4 4.5
##  [163] 4.8 4.3 3.5 4.5 4.7 4.0 4.5 3.8 4.8 2.9 3.8 5.2 3.6 4.4 4.6 3.2 5.1 6.2
##  [181] 4.6 5.4 4.3 3.6 5.5 4.9 3.6 5.4 2.9 3.8 3.9 5.2 3.2 5.2 5.3 5.4 4.0 3.3
##  [199] 4.2 4.6 2.0 2.4 4.8 3.8 4.9 3.7 5.2 4.3 4.5 4.5 4.8 5.2 4.4 3.1 3.4 4.5
##  [217] 3.4 4.9 4.2 3.9 5.5 4.4 5.8 5.2 4.0 4.9 4.5 4.3 3.8 5.3 4.0 3.8 5.5 5.0
##  [235] 5.3 4.3 7.0 5.6 4.1 4.3 3.5 4.2 4.5 4.3 5.4 3.0 5.4 6.6 3.7 4.7 3.7 4.5
##  [253] 4.0 4.8 3.1 4.1 2.7 3.5 4.1 5.0 4.8 3.3 3.6 2.9 4.5 4.0 5.0 5.1 4.4 4.9
##  [271] 4.4 3.5 4.4 5.3 3.9 3.8 4.7 4.3 4.3 5.5 3.7 4.2 4.3 5.4 4.5 4.0 2.9 4.6
##  [289] 3.7 4.6 3.5 5.3 3.0 4.8 6.2 4.7 4.1 4.5 5.1 3.8 3.9 5.7 4.2 4.3 4.0 4.4
##  [307] 3.7 4.4 5.0 4.4 5.8 5.5 6.2 3.5 5.7 4.4 3.8 4.5 5.0 5.3 5.8 4.3 3.9 6.3
##  [325] 2.8 6.2 4.1 6.1 4.9 3.7 5.6 5.2 4.0 5.2 5.5 4.3 5.0 5.2 3.5 4.2 3.6 5.5
##  [343] 3.8 4.0 6.1 6.0 4.0 5.5 4.6 3.7 3.6 3.8 5.7 5.3 5.8 4.5 4.1 4.9 3.6 4.4
##  [361] 4.2 4.2 2.9 5.0 4.3 4.3 4.6 2.8 5.2 4.4 4.1 3.6 3.9 4.3 5.1 5.9 4.2 4.5
##  [379] 3.9 4.2 5.5 4.1 3.4 5.4 3.9 3.4 5.2 3.9 4.9 5.3 6.1 3.4 5.2 4.2 6.0 5.4
##  [397] 3.1 4.5 4.0 4.9 3.1 3.9 4.7 6.3 4.1 5.6 4.2 3.9 3.4 3.5 2.4 2.7 2.7 4.0
##  [415] 4.9 7.0 4.9 3.1 5.8 5.4 3.8 5.9 3.1 2.7 3.9 5.7 5.4 5.4 4.6 4.8 3.6 3.6
##  [433] 5.8 3.9 1.9 2.2 4.9 5.6 3.3 6.3 3.2 2.7 5.2 3.2 4.9 3.8 6.1 4.2 5.4 3.7
##  [451] 4.5 3.8 4.5 5.4 4.7 5.6 5.3 5.4 5.0 3.2 4.0 5.5 3.3 5.3 4.1 4.8 3.7 4.3
##  [469] 3.5 4.5 4.8 4.7 4.7 4.3 3.5 3.9 5.0 5.2 3.9 5.6 4.5 4.0 5.9 5.4 5.2 4.7
##  [487] 5.2 6.1 4.3 4.7 4.1 6.0 4.6 5.9 4.1 5.3 4.1 2.4 5.4 4.5 4.1 4.3 5.6 4.8
##  [505] 4.7 4.9 4.7 4.6 6.0 2.9 5.4 4.4 3.8 3.5 4.4 3.4 3.5 4.4 5.5 6.4 4.9 3.4
##  [523] 6.7 5.2 3.1 4.0 5.2 4.0 4.0 5.4 4.8 4.1 4.9 5.8 5.5 6.0 3.3 3.5 3.4 3.6
##  [541] 4.1 5.0 3.8 4.5 2.9 4.7 3.3 7.0 3.5 5.9 4.9 5.0 4.6 3.3 5.4 4.9 4.0 3.6
##  [559] 5.0 2.5 4.0 5.2 3.0 4.1 3.6 2.8 5.6 4.7 4.5 4.9 4.5 3.7 6.3 5.2 6.0 4.4
##  [577] 4.4 3.1 4.8 4.4 3.4 4.5 3.7 4.8 5.2 4.3 4.6 3.8 3.8 3.6 5.1 6.3 4.6 3.9
##  [595] 4.7 5.4 4.2 5.9 4.7 3.5 5.2 3.9 4.1 4.6 5.9 4.2 2.5 5.0 4.8 5.1 5.1 4.8
##  [613] 5.2 4.6 2.9 6.2 4.0 5.7 5.7 5.8 3.6 4.0 5.5 4.8 4.3 4.5 5.1 3.5 5.4 3.8
##  [631] 4.4 4.6 4.8 5.3 6.1 5.2 5.8 4.0 4.1 3.7 4.1 5.4 3.6 3.9 3.8 5.1 4.0 3.4
##  [649] 3.2 5.0 4.5 2.9 4.4 3.7 5.3 3.5 4.3 4.8 4.7 3.5 3.9 5.0 5.2 4.8 4.6 4.1
##  [667] 3.9 3.1 3.6 4.6 3.9 4.7 7.0 4.4 4.0 5.8 5.7 3.7 4.6 5.2 5.7 3.2 2.8 2.5
##  [685] 4.3 6.0 5.1 5.8 3.7 3.9 4.1 4.7 4.9 4.9 5.2 2.8 3.3 4.6 4.2 4.7 4.2 4.2
##  [703] 6.4 4.9 4.8 4.5 4.3 3.1 4.6 6.1 3.6 5.2 5.3 4.1 5.4 3.4 4.8 3.7 6.2 4.8
##  [721] 4.3 4.3 1.9 4.7 4.9 4.1 4.8 3.8 5.8 4.1 4.4 4.4 5.1 3.8 3.6 6.9 5.7 3.6
##  [739] 2.8 3.8 5.0 4.5 4.0 4.2 5.4 3.5 4.3 5.1 2.6 5.9 4.7 5.5 4.2 5.0 5.2 3.6
##  [757] 5.3 4.8 2.1 5.4 4.7 3.5 5.8 4.9 3.0 4.1 3.5 3.5 3.8 4.3 5.1 5.1 5.4 3.0
##  [775] 4.4 1.8 4.3 5.6 3.6 5.2 6.2 4.5 4.1 5.7 5.2 4.1 5.5 4.0 4.5 5.8 4.0 3.5
##  [793] 5.6 5.5 3.8 4.8 1.6 4.1 4.7 1.6 4.8 5.4 2.3 5.3 3.9 3.4 5.4 5.8 5.4 3.7
##  [811] 5.5 3.9 3.9 3.1 3.6 4.5 5.3 3.6 5.6 4.0 4.7 5.0 6.3 4.4 4.1 4.3 6.5 4.2
##  [829] 5.2 3.5 3.8 5.8 6.0 3.6 6.7 4.5 5.7 4.4 3.8 5.0 3.9 4.3 4.6 3.1 5.5 4.9
##  [847] 3.5 5.5 3.6 3.7 3.2 7.0 5.2 3.2 5.4 3.6 5.5 4.3 6.0 3.7 4.4 3.6 4.6 5.8
##  [865] 5.4 4.9 5.0 4.7 4.7 4.2 3.5 3.9 5.1 6.3 4.5 6.3 3.5 6.1 6.1 4.6 5.5 4.5
##  [883] 3.1 5.4 2.8 3.3 4.0 5.2 4.9 6.4 4.1 5.0 4.7 3.6 5.9 4.4 5.4 5.7 5.8 4.9
##  [901] 3.8 5.8 5.5 4.0 3.0 3.1 3.8 3.6 4.6 4.4 4.5 5.2 5.7 4.6 5.5 3.6 5.5 3.8
##  [919] 3.6 4.8 5.5 4.6 5.4 4.1 4.3 3.8 3.7 4.8 3.6 4.9 4.4 4.7 4.1 4.8 4.1 3.1
##  [937] 4.1 3.5 4.6 3.0 6.4 4.4 5.1 4.1 3.8 5.0 4.1 3.3 4.3 2.3 3.6 4.3 5.1 4.9
##  [955] 4.2 5.2 5.1 4.3 4.2 5.7 3.3 5.0 4.4 4.4 4.2 3.0 5.2 4.7 4.1 5.0 7.3 4.9
##  [973] 6.2 3.9 3.3 4.4 2.8 5.9 3.4 4.4 4.1 3.8 5.6 3.9 4.4 3.3 4.4 6.1 4.0 4.2
##  [991] 3.5 5.3 5.8 4.0 4.2 5.4 3.2 6.3 3.9 5.4
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
##    [1] 4.3 4.4 2.9 2.9 5.4 5.0 3.9 5.0 2.9 4.7 6.6 5.5 5.0 4.5 4.0 2.8 4.2 4.3
##   [19] 3.7 5.1 5.4 4.9 3.7 6.2 3.0 4.0 2.8 3.7 5.0 5.1 3.9 4.3 5.4 5.6 4.8 3.6
##   [37] 4.1 5.2 4.5 1.7 3.7 4.9 4.8 3.3 4.0 4.6 5.1 6.0 4.4 3.5 5.6 5.7 4.7 3.2
##   [55] 4.7 5.7 5.1 4.1 5.1 5.4 3.7 4.3 3.0 4.7 5.9 4.5 5.3 4.6 3.0 3.3 2.2 4.5
##   [73] 3.6 2.5 4.7 5.6 4.8 4.4 4.0 4.4 4.3 4.3 3.7 3.9 4.6 5.7 3.9 3.9 4.5 3.0
##   [91] 4.1 3.4 5.5 5.6 4.2 4.2 5.3 5.5 4.1 5.6 3.7 4.7 4.8 4.3 4.7 4.0 3.1 3.1
##  [109] 6.4 5.5 4.8 5.4 4.7 5.3 4.9 5.5 2.9 4.7 4.6 4.7 4.8 4.4 4.8 5.2 5.5 5.7
##  [127] 5.2 5.7 4.7 4.4 3.0 2.4 5.6 4.4 3.6 5.6 6.0 5.0 3.6 2.6 5.1 4.8 5.2 4.5
##  [145] 2.9 3.7 4.8 5.6 4.6 3.4 3.4 3.9 5.9 5.5 4.4 4.6 5.0 6.1 4.7 4.7 6.2 3.1
##  [163] 2.7 5.5 4.8 4.1 4.8 2.6 4.8 5.8 4.5 5.3 3.7 4.0 2.9 3.9 6.1 4.5 3.4 3.5
##  [181] 5.0 3.9 2.6 2.8 3.2 5.7 3.7 4.5 5.5 4.9 5.7 6.1 6.0 4.6 4.3 4.4 3.4 5.0
##  [199] 5.2 3.9 5.5 3.8 3.6 3.1 2.8 4.7 3.2 4.4 5.0 5.0 3.8 3.1 3.0 5.7 5.1 3.1
##  [217] 4.4 3.6 6.2 2.8 3.3 3.0 4.6 4.7 6.4 4.7 4.9 4.1 4.8 2.2 3.9 5.0 5.0 4.6
##  [235] 5.8 5.2 3.3 4.4 6.0 5.2 4.2 4.1 3.6 3.9 4.0 4.7 2.7 3.3 4.4 3.6 3.9 4.6
##  [253] 6.0 4.6 3.7 4.0 4.6 3.7 5.9 3.2 4.2 4.4 3.8 4.9 3.9 4.8 4.9 5.2 6.0 4.1
##  [271] 4.3 5.5 4.4 5.4 4.0 2.4 4.9 4.1 4.0 3.5 3.3 5.2 4.2 3.9 4.3 2.6 4.5 3.8
##  [289] 5.4 4.0 5.3 2.6 5.3 3.9 5.8 6.5 3.9 5.6 6.5 5.6 3.7 3.1 3.8 4.8 5.1 3.4
##  [307] 4.2 4.9 5.5 4.4 4.2 4.2 4.5 5.6 5.8 4.0 4.4 4.7 4.2 4.9 5.3 5.6 4.0 3.9
##  [325] 6.5 3.6 4.9 4.8 3.1 3.3 4.2 3.8 3.4 3.5 5.2 3.7 4.3 5.7 5.9 5.5 2.9 5.2
##  [343] 5.3 3.1 4.5 4.1 3.9 4.6 4.8 4.7 3.5 3.3 5.3 5.7 4.8 5.1 3.4 5.3 5.3 3.0
##  [361] 5.0 2.8 3.7 3.6 4.5 6.0 4.4 4.8 5.3 4.7 5.2 5.4 6.6 4.7 3.9 4.4 3.1 2.9
##  [379] 5.0 5.0 4.8 4.9 3.9 5.3 7.3 4.7 3.6 4.5 3.8 3.4 3.6 3.6 5.3 5.0 5.9 4.7
##  [397] 2.8 3.8 5.7 4.3 4.3 5.2 4.9 4.2 5.8 4.7 2.9 3.9 5.0 3.5 5.3 4.4 5.8 3.8
##  [415] 5.4 5.7 4.1 3.1 3.7 5.6 3.5 4.9 3.4 5.6 5.4 4.1 4.9 3.7 4.6 5.2 4.8 6.4
##  [433] 6.1 3.7 2.9 4.0 5.5 4.6 3.7 4.4 5.3 4.6 4.6 3.0 5.3 3.8 2.6 4.3 5.9 4.7
##  [451] 4.8 3.8 3.4 4.6 3.8 3.3 5.1 2.7 4.5 5.5 3.1 3.8 6.0 5.0 1.6 4.1 5.3 3.9
##  [469] 3.9 4.5 4.1 4.1 5.9 3.2 5.6 4.3 4.8 4.0 5.1 6.4 4.4 2.7 5.5 4.0 4.8 4.8
##  [487] 3.8 4.9 4.1 4.5 5.9 4.1 5.3 3.9 4.2 3.9 5.4 5.8 5.3 3.9 3.3 3.6 4.6 5.3
##  [505] 6.2 3.8 4.3 4.5 4.4 4.6 4.9 4.1 4.0 3.0 4.7 4.7 4.3 3.5 4.6 4.6 6.1 4.2
##  [523] 5.7 5.4 4.5 4.3 5.3 4.4 5.5 4.1 5.1 4.8 4.9 4.0 5.0 5.1 4.5 4.7 5.6 4.1
##  [541] 3.8 5.0 4.6 4.2 5.4 3.9 5.3 4.5 4.5 3.5 3.9 4.8 3.9 5.8 3.2 5.1 4.8 4.1
##  [559] 4.4 3.3 2.6 3.1 4.9 4.0 4.0 5.2 4.0 5.6 3.5 6.1 4.7 3.4 4.4 5.3 6.1 3.5
##  [577] 5.0 5.2 4.8 5.7 3.4 5.4 2.8 6.5 4.1 2.1 3.8 3.5 3.3 5.0 4.4 4.1 6.3 3.2
##  [595] 5.7 3.2 5.0 5.3 6.1 4.4 3.5 3.7 5.4 5.6 5.4 5.1 3.7 4.9 4.6 3.8 3.7 6.4
##  [613] 4.8 3.2 5.5 3.6 4.9 6.9 4.1 3.7 4.8 3.4 4.6 5.3 3.9 5.3 5.5 5.4 4.8 4.9
##  [631] 5.0 2.8 5.4 6.5 4.7 5.3 5.0 3.0 4.7 4.6 4.6 5.2 5.3 4.9 4.2 3.6 5.3 3.6
##  [649] 5.4 3.8 5.6 5.3 6.2 5.2 4.9 3.9 6.4 3.4 3.4 4.3 4.3 5.4 3.3 6.0 4.6 5.2
##  [667] 4.9 3.6 3.1 4.0 4.1 6.3 4.5 5.2 3.7 4.3 5.0 5.4 3.4 5.4 6.3 4.1 6.4 4.2
##  [685] 4.7 5.8 4.9 5.1 3.1 4.4 4.4 6.3 3.9 6.4 5.2 4.0 5.4 3.0 3.5 2.6 4.1 4.8
##  [703] 5.1 3.9 4.4 4.2 6.1 5.1 5.0 3.3 4.5 5.3 5.8 5.2 4.1 4.9 3.4 3.7 2.2 3.8
##  [721] 3.8 5.0 3.5 3.6 5.9 4.9 5.2 4.1 5.5 4.9 5.2 3.8 3.8 5.2 5.5 4.1 4.4 4.3
##  [739] 5.0 5.8 5.4 4.7 5.2 3.5 4.9 3.1 5.8 5.2 3.4 4.2 4.5 4.5 4.5 3.0 4.3 5.2
##  [757] 3.4 3.8 4.5 5.6 5.3 4.7 5.7 4.9 5.6 6.2 3.9 4.8 5.3 3.8 5.3 4.9 4.7 3.5
##  [775] 4.5 5.3 4.0 4.0 5.2 4.0 4.5 4.3 5.7 3.7 3.6 4.7 4.5 4.5 4.6 5.9 3.8 3.9
##  [793] 4.5 5.4 4.8 6.3 4.4 4.3 4.1 6.1 4.5 4.7 4.8 4.5 4.9 5.0 4.3 5.4 4.2 3.0
##  [811] 3.8 4.0 4.1 2.8 4.0 4.7 4.8 4.0 4.1 5.0 3.9 5.1 4.7 4.0 3.7 4.3 6.1 4.4
##  [829] 4.9 4.7 3.4 5.2 3.9 3.8 4.4 5.1 6.5 6.1 5.0 4.6 4.3 4.2 3.9 4.9 4.9 5.7
##  [847] 3.5 5.0 5.3 4.8 4.9 6.0 3.0 5.1 4.7 3.8 4.1 3.2 5.1 3.6 2.5 5.0 4.3 2.9
##  [865] 5.5 4.3 5.3 4.6 4.5 3.4 3.4 5.3 5.2 4.0 4.0 3.9 4.5 4.1 4.3 5.3 4.5 4.2
##  [883] 4.7 3.8 3.3 4.9 5.5 3.8 2.7 4.7 3.1 6.3 4.9 4.3 4.0 5.9 3.5 5.6 4.7 4.6
##  [901] 4.5 6.3 4.2 4.5 3.9 4.3 3.4 5.8 4.9 4.9 3.2 6.1 3.8 5.2 6.0 4.7 5.8 3.3
##  [919] 4.5 4.7 6.0 3.8 3.6 4.7 3.8 5.6 4.8 3.0 4.3 4.8 4.3 3.7 5.8 5.3 4.9 5.5
##  [937] 3.2 4.8 5.3 5.7 4.4 5.6 4.5 2.8 5.1 4.6 3.3 5.8 4.1 4.0 3.3 5.0 6.2 3.6
##  [955] 4.0 5.0 5.3 4.4 3.9 3.8 4.9 5.2 4.0 3.8 4.7 5.2 4.5 5.2 4.3 4.4 5.0 2.7
##  [973] 4.7 4.3 3.6 3.6 1.7 5.0 4.4 4.9 6.8 4.5 3.2 2.7 4.8 4.2 5.2 5.1 2.8 6.4
##  [991] 5.4 4.0 4.9 5.5 4.6 4.1 5.4 3.2 4.2 5.4
## 
## $func.thetastar
## [1] 0.0058
## 
## $jack.boot.val
##  [1]  0.509631728  0.431073446  0.285507246  0.185161290  0.026210826
##  [6] -0.008108108 -0.164265130 -0.304255319 -0.379829545 -0.563501484
## 
## $jack.boot.se
## [1] 1.014439
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
##    [1] 4.1 6.1 5.7 3.9 4.6 4.3 5.6 6.1 4.4 5.4 2.8 5.6 3.7 4.3 4.9 5.3 4.3 5.7
##   [19] 5.0 3.9 4.7 4.0 3.7 4.6 3.6 4.3 4.6 5.1 5.0 6.2 4.4 5.0 4.7 5.6 4.7 5.2
##   [37] 3.5 4.1 4.0 4.0 3.6 4.1 5.7 3.3 4.5 3.1 3.5 3.4 4.3 4.4 6.1 2.7 4.2 5.2
##   [55] 3.8 4.5 5.7 6.1 4.9 5.9 3.8 3.2 3.4 4.3 3.5 3.9 4.1 4.3 4.0 3.4 4.9 3.4
##   [73] 5.1 3.9 3.9 4.9 4.2 5.5 4.1 4.2 5.5 5.4 4.2 5.8 5.9 5.0 3.3 5.0 4.3 5.0
##   [91] 3.8 4.0 3.8 5.0 5.7 5.5 4.8 2.5 4.6 3.1 4.8 4.6 4.2 5.3 4.3 3.6 3.3 3.6
##  [109] 5.5 4.3 4.2 5.3 4.2 3.6 5.4 3.3 4.7 6.0 5.1 3.8 3.2 3.8 2.9 4.4 4.0 4.0
##  [127] 2.9 3.5 4.2 3.5 5.1 2.8 5.4 4.2 4.5 4.2 4.4 5.6 5.0 3.3 3.3 4.2 4.8 3.8
##  [145] 5.0 6.5 4.9 5.0 4.2 5.6 3.2 5.0 5.7 4.1 6.2 3.6 3.4 2.2 3.5 5.1 2.6 3.9
##  [163] 3.0 4.3 5.1 4.2 4.4 4.5 4.5 3.4 5.1 5.7 4.0 4.6 6.4 4.3 4.3 5.6 5.1 4.1
##  [181] 5.7 3.6 3.7 4.2 4.8 3.4 5.4 4.6 3.0 6.7 4.2 4.1 3.8 6.7 4.4 3.3 4.1 4.3
##  [199] 4.9 6.0 4.6 4.7 5.3 5.1 5.6 5.7 3.2 3.4 3.6 3.0 6.1 4.2 4.6 4.8 4.2 4.4
##  [217] 4.9 4.1 3.6 5.6 3.9 4.3 5.3 4.4 2.8 5.6 3.3 5.4 4.7 3.5 3.8 4.4 4.3 3.1
##  [235] 4.0 4.4 4.0 6.1 4.0 5.6 5.1 3.7 4.7 4.1 4.9 4.7 5.4 4.2 3.8 3.6 3.4 3.2
##  [253] 2.4 4.1 4.3 5.2 4.6 4.1 6.2 3.3 4.2 5.3 4.9 3.8 4.9 3.7 4.6 4.0 4.4 3.4
##  [271] 5.3 3.7 4.9 3.8 5.0 4.7 5.1 4.7 3.2 5.1 4.0 3.4 3.7 4.7 5.1 4.8 3.8 4.1
##  [289] 4.1 4.9 5.7 4.5 3.8 4.6 4.9 2.7 3.3 4.7 3.6 4.7 5.1 5.4 4.3 5.0 4.4 4.3
##  [307] 3.8 4.6 5.5 3.7 3.7 2.6 5.5 4.5 4.6 4.2 6.2 2.3 5.3 5.1 4.9 4.7 4.8 3.9
##  [325] 4.5 4.3 4.8 4.4 4.8 4.4 4.6 3.6 3.8 5.3 3.3 6.2 5.7 3.1 4.4 4.3 4.9 6.1
##  [343] 6.3 3.4 3.8 4.3 4.1 4.5 2.2 3.8 2.7 6.7 5.6 4.0 5.1 3.7 4.3 6.4 5.2 3.1
##  [361] 3.4 5.0 2.9 5.6 3.9 3.6 2.5 4.5 5.3 4.5 5.3 4.2 3.3 4.2 5.7 5.3 2.9 3.6
##  [379] 4.4 3.7 4.6 4.4 4.3 4.7 3.5 4.2 4.8 5.3 4.5 4.7 4.3 5.8 4.6 3.6 3.4 4.0
##  [397] 4.9 5.4 2.8 5.2 4.5 5.0 4.1 3.6 5.1 4.0 4.5 4.3 4.9 4.6 4.4 3.8 3.7 3.8
##  [415] 3.9 6.5 4.4 4.0 3.9 3.1 4.3 6.9 3.0 4.3 4.0 3.8 3.3 5.0 5.1 4.7 4.1 4.7
##  [433] 2.0 2.9 4.2 4.4 5.4 4.3 5.7 4.8 4.0 5.6 4.7 3.3 5.4 5.8 4.7 6.0 4.5 4.5
##  [451] 4.2 3.8 4.4 5.4 6.0 4.9 2.8 5.5 4.8 5.6 5.2 2.8 4.8 2.8 5.6 5.2 4.0 3.8
##  [469] 4.0 4.0 5.3 4.9 4.8 3.2 5.6 5.2 4.3 4.6 5.2 5.0 5.3 5.3 4.4 4.5 6.0 5.7
##  [487] 5.1 5.5 4.0 4.7 3.9 4.1 4.4 4.5 3.8 3.8 4.3 4.9 4.7 4.1 4.3 3.5 4.0 3.9
##  [505] 3.0 4.1 2.9 6.7 4.3 3.5 4.9 5.1 4.4 2.5 5.3 4.6 3.5 3.5 3.4 4.9 4.4 5.3
##  [523] 5.6 5.1 5.1 2.9 5.6 5.2 3.5 3.6 4.7 5.6 5.2 3.6 4.2 5.1 5.5 4.4 2.7 4.4
##  [541] 5.0 5.2 5.7 4.3 2.1 2.3 5.0 4.2 3.2 4.2 5.8 4.3 4.2 5.4 4.4 4.5 5.0 5.2
##  [559] 6.3 3.5 4.0 5.0 6.0 4.1 4.8 3.5 4.1 5.6 4.6 4.3 4.5 3.8 6.5 4.3 3.8 4.6
##  [577] 3.6 4.9 4.2 4.0 5.0 5.0 5.2 5.0 5.3 4.6 5.4 5.2 4.1 3.3 4.3 3.5 4.5 3.5
##  [595] 4.2 4.8 4.5 4.8 5.3 4.0 3.5 3.2 4.3 6.2 3.4 5.3 6.5 3.8 4.2 4.6 6.9 5.9
##  [613] 3.7 3.8 3.1 5.5 5.1 6.4 5.7 3.4 5.8 5.2 5.5 4.4 5.0 4.3 5.2 4.9 4.1 4.7
##  [631] 3.9 4.8 3.5 3.4 5.4 5.2 4.7 6.5 3.9 6.1 5.6 4.3 4.8 3.3 3.6 3.8 5.8 4.4
##  [649] 5.7 4.1 4.0 5.3 4.9 7.0 2.8 4.2 3.5 5.7 5.0 4.2 3.6 4.0 4.7 4.5 3.8 4.2
##  [667] 5.3 6.2 4.7 4.3 3.5 3.0 4.4 4.7 3.4 4.3 5.4 4.5 5.4 4.3 5.6 4.7 4.5 5.7
##  [685] 4.2 5.3 4.8 3.4 5.1 5.6 4.1 4.6 5.7 4.6 3.7 6.1 3.5 4.3 5.8 2.4 3.9 3.7
##  [703] 4.4 2.3 5.9 5.5 2.8 4.1 4.2 4.7 5.6 3.4 3.4 5.5 4.2 3.4 3.4 4.1 3.4 5.8
##  [721] 3.3 5.3 4.1 3.9 5.2 4.1 5.9 4.2 2.9 3.9 4.7 4.3 4.4 5.7 4.8 4.0 4.3 2.9
##  [739] 2.8 5.3 4.5 4.2 3.9 5.1 2.3 4.4 3.6 6.4 3.3 3.8 4.4 4.0 3.4 3.0 3.1 4.1
##  [757] 4.6 4.2 5.7 6.6 4.1 4.4 4.0 4.6 5.3 4.9 5.4 4.5 5.0 4.3 2.8 4.2 4.0 4.4
##  [775] 6.4 5.7 5.6 3.6 5.8 5.3 4.6 4.6 3.8 3.2 5.8 3.0 3.1 4.0 4.8 3.3 4.4 5.0
##  [793] 5.3 2.8 5.9 4.3 4.0 4.0 4.9 5.1 3.8 4.6 3.1 5.4 4.8 5.0 4.7 4.6 6.1 4.1
##  [811] 5.3 3.7 4.7 4.5 3.7 4.4 2.2 5.4 4.1 5.2 4.9 4.7 4.9 3.5 3.2 5.1 3.9 4.0
##  [829] 4.0 4.5 2.3 5.1 5.4 4.8 3.0 5.0 3.8 3.8 3.7 5.1 5.3 5.4 3.9 5.5 5.4 5.1
##  [847] 3.7 2.8 6.0 3.8 3.7 4.6 3.6 4.0 2.3 4.8 5.5 4.4 5.0 3.8 3.0 5.3 3.0 4.6
##  [865] 4.3 5.1 3.9 5.2 1.9 4.4 3.8 3.5 6.0 4.6 4.5 6.7 3.0 4.8 4.4 3.3 3.4 3.9
##  [883] 4.9 3.7 3.5 3.0 6.4 4.4 4.9 3.2 5.2 3.6 4.4 4.1 4.0 5.0 6.8 5.3 5.9 4.5
##  [901] 4.1 3.8 2.5 4.7 4.4 4.8 3.7 3.9 3.7 4.2 5.3 3.9 5.8 5.9 6.4 3.9 5.9 3.0
##  [919] 4.1 3.8 4.8 3.1 4.2 4.6 4.7 4.9 4.9 4.7 5.1 5.5 5.0 4.4 4.7 3.4 2.5 5.5
##  [937] 5.2 5.7 4.5 4.9 4.7 5.6 4.5 4.3 6.0 5.8 4.1 4.3 4.1 4.8 3.2 5.2 4.1 3.0
##  [955] 3.9 2.6 3.4 3.6 3.6 3.6 5.6 5.3 4.3 5.5 4.2 5.6 4.6 5.3 3.3 4.6 4.4 4.3
##  [973] 5.6 5.6 4.1 3.9 3.5 3.8 6.3 5.5 3.6 4.3 5.2 3.4 5.0 6.2 4.9 3.2 5.6 4.7
##  [991] 5.2 7.4 4.6 4.7 5.8 4.8 4.8 1.8 5.1 4.0
## 
## $func.thetastar
## 72% 
##   5 
## 
## $jack.boot.val
##  [1] 5.50 5.40 5.30 5.20 5.10 5.00 4.80 4.60 4.50 4.44
## 
## $jack.boot.se
## [1] 1.086249
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
## [1] 0.3557013
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
##    7.976570   11.215429 
##  ( 3.495170) ( 5.072277)
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
## [1]  0.4944246  0.3639078  0.1585412 -0.5219437  0.1466399  0.2998710
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
##    [1]  0.296716732  0.584212866  0.395954881  0.716782597  0.444555989
##    [6]  0.495907223 -0.273432577  0.541733368 -0.363747796 -0.265736001
##   [11]  0.481772892 -0.511429086 -0.956923614  0.060503676  0.182451386
##   [16]  0.155909932 -0.283214565  0.028687447 -0.065259211  1.269223109
##   [21] -0.549628240  0.254237180  0.389821501  0.510492668 -1.000096542
##   [26]  0.044145452 -0.279337140  0.389395584 -0.317778437  0.259027766
##   [31]  1.032417918  0.094774262  0.229529211  0.607526284  1.242174906
##   [36] -0.531818038 -0.159488613  0.730947534 -0.817182006 -0.203585598
##   [41]  0.174821661  0.523720041 -0.463424280  1.104547200  0.296804615
##   [46]  0.312883114  0.738789738  0.040193690 -0.546630741  0.458249318
##   [51] -0.163078879 -0.422854411 -0.203661744 -0.169473798 -0.546970686
##   [56]  0.188853085  0.204002037  0.182811729  0.008922207 -0.299319441
##   [61]  2.159200720  0.560432245 -0.094558159 -0.141072805  0.019264576
##   [66]  0.392932919 -0.217949931  0.082043675  0.487061515 -0.027748661
##   [71]  0.226294293  0.297143847 -0.280449913  0.803868839  0.366697737
##   [76]  0.981869898  0.041502108  0.082574580 -0.897380963  0.546433553
##   [81]  1.553893066  2.147323790  0.137637219  0.127577542  0.068810181
##   [86]  0.804127660 -0.339427230  0.943648992  0.042602563  1.523563463
##   [91]  0.488390949  0.440015182  0.744666436  0.418058599  1.081228897
##   [96]  0.261514605  0.453003908  0.391329600 -0.220924883 -0.629583875
##  [101]  0.791285487  1.192729610 -0.078404970 -0.240030347 -0.297661891
##  [106]  0.054930801 -0.185051102  0.169807067 -0.088703805  0.892758542
##  [111]  0.618069755  1.897047712  1.273696399  0.815790683  0.522532519
##  [116] -0.069574130  1.372168486 -0.675708380  0.440251237  2.136135500
##  [121] -0.244462504  0.805282456 -0.190382794 -0.032507113  0.575478942
##  [126] -0.592657530  0.393622421  0.911379808  0.577388312  0.566749108
##  [131] -0.749044785  1.383795023  0.307713625  0.778923524  0.893897057
##  [136]  1.064498394  0.057108579 -0.062146281  0.366266360 -0.200422166
##  [141]  0.534302570 -0.783845911  1.273860154 -0.561077989  0.832833086
##  [146]  0.352569474  0.867943871 -0.648104753 -0.880497952  0.643248164
##  [151]  0.451387671  0.427044085  0.708964594 -0.242629819  0.785229290
##  [156] -0.050430551  0.403203767  0.373324088  1.209797618 -0.124428299
##  [161] -0.259543160  0.027156376  1.385579131  1.514361117  0.569607589
##  [166]  0.506839956  0.560736019 -0.481385579  0.847430877 -0.213011049
##  [171]  0.503840118 -0.221257776 -1.139177343  0.465187889 -0.207635023
##  [176]  0.025881889  0.201030364  0.433719920  0.041765692  1.288482340
##  [181]  0.587035307 -0.407438057  0.056246918 -0.088367990 -0.411217749
##  [186]  1.327246201  0.570141853  0.109404029 -0.086496110  0.168032814
##  [191]  0.991553602  0.569093374 -0.615708215  1.019571771  0.039305127
##  [196]  0.351523497  0.762026999  0.596648989 -0.325605742  0.120613585
##  [201]  0.689977479  0.307395969  0.944559540  0.872505176  0.749859351
##  [206]  0.983449417  0.513110019 -0.296402304 -0.548359668  0.416102634
##  [211]  0.208981529  0.073285231  0.606830072  0.778327554  0.777728125
##  [216] -0.301805838 -0.819108139 -0.431091564  0.225364193  0.359265195
##  [221] -0.182991710 -0.170465373  0.988462631  0.387749848 -0.612151139
##  [226]  1.513175735 -0.183801104  0.046435889  1.218926396  1.561166906
##  [231]  0.967590311 -0.140441958  0.195456834  0.991630848  0.349785063
##  [236] -0.125980441  0.070511829 -0.039065264  1.268140293  0.076583082
##  [241]  0.573343759  0.549075059  0.037899696  0.080053802  2.308780965
##  [246]  0.297397314 -0.136213715  0.733745486  0.522532519 -0.054427047
##  [251]  0.215313286  0.462659214 -0.350133734 -0.329761284  0.611496004
##  [256]  0.380583116  0.979907712 -0.074933158 -0.111870071  0.226555049
##  [261]  0.163795614 -0.734291321  0.352708068  0.069779937  0.154621234
##  [266]  0.219088302  0.823645855  1.397087164 -0.077381420  0.507977257
##  [271]  0.882220847  0.356200345  0.226723653 -0.352297386  1.607989487
##  [276]  1.246699871  0.093792505  0.383625900  0.488889899  0.980229851
##  [281]  0.522886226  0.549709477 -0.003447864  0.555039836 -0.022442170
##  [286] -0.074928350  0.388086130  1.275845654 -0.032507113 -0.210634997
##  [291]  0.776699520 -1.071004390  0.573566759  0.065199350  0.594538184
##  [296] -0.040395941 -1.249202934  1.622415078  0.148997661  0.741348258
##  [301]  0.886648921 -0.231451695  0.778096030 -0.169526495  0.556462498
##  [306]  0.285404018  0.544785920  0.335760574 -0.116512627  0.517811409
##  [311] -0.344631327  0.518804748  1.610277751  0.251818391  0.488381293
##  [316] -0.097704457  0.334453664  0.310534188 -0.155831795  0.712925167
##  [321]  0.078470621 -0.192327798  0.257250751  0.021644015 -0.361830031
##  [326]  0.342727702  0.414823135  0.013393997 -0.297442154  0.394268375
##  [331]  0.430112394  0.342481696  0.973830190 -0.150727856 -0.125771028
##  [336]  0.204332692  0.639617223  0.790140739  0.721854991  1.027123763
##  [341]  0.732768764  0.155348316  0.758843979 -0.015030461  0.399979955
##  [346] -0.025081714  0.350334408  1.379955933  0.802969968 -0.343683563
##  [351] -0.578184208  0.377156125 -0.054092810  0.150382640  0.419815794
##  [356]  0.581742661  0.343535230  0.553573222 -0.044878709 -0.319241161
##  [361]  0.349398619  0.323048885 -0.109129597  0.742834900  0.344682443
##  [366] -0.733625989  0.832833086 -0.334755004  0.547796084  0.015168161
##  [371] -0.404362253  0.392836778 -0.268670082 -0.097847752  0.615690241
##  [376] -0.041543510 -0.143182906  0.250578218  0.303897647  0.756359184
##  [381]  1.437386133  0.370479863 -0.191555443  0.128581977  0.332223234
##  [386]  0.046355243  1.214339828  0.735216650  0.964669136  0.333123850
##  [391] -0.252232218  0.317096384  0.594179689 -0.866876393  0.061314324
##  [396] -0.347818576  0.757146187  0.454788239 -0.079989939 -0.789099773
##  [401]  0.737136502  0.478659292 -0.730379911 -0.499378682 -0.747547312
##  [406]  0.190867213  0.416314346  1.310442492 -0.032202230  0.006919948
##  [411] -0.203459538 -0.147121345  0.328135677 -0.407735223  0.967157823
##  [416] -0.033325336  0.300772038  0.189468999 -0.138434434  0.166077929
##  [421]  0.060513046  0.040862167 -0.110037591  0.835349490 -0.765111822
##  [426] -1.493936985  0.515578429 -0.226917423  0.754748104  0.747491832
##  [431]  0.204917468  0.884564833 -0.217016864  0.891103337 -0.322126143
##  [436]  0.376376015 -1.383395287 -0.342277820  0.347025260 -0.160128375
##  [441] -0.222118780  0.405607768 -0.019958354  0.730824181  0.569855014
##  [446]  1.261630048  1.898189810  0.490322476  0.041479299 -0.012251654
##  [451] -0.495901493  1.222360350 -0.121386969  0.781767414  0.415645550
##  [456]  0.741841640  0.373426461  1.380142913 -0.115488486 -0.303391534
##  [461] -0.112390024  0.573334827  0.969268979  0.043254786  0.385061870
##  [466]  0.355650779 -0.396659754  0.489381024  1.233995970  0.601594399
##  [471]  0.597969085  0.926529341  0.588939649  0.714860287  0.419913012
##  [476] -0.282448382  0.259004494 -0.010554912  0.024706563  0.062035107
##  [481]  0.333920405 -0.605941381  0.574372494  0.302941847  0.716224759
##  [486]  2.101110475  0.158130183 -0.160516135 -0.459551455 -1.559589863
##  [491]  0.004070161  0.601594399  0.779173028  0.738191540 -0.202909487
##  [496]  0.104242887  0.352437544  0.079749127  0.219208509 -0.375707068
##  [501]  0.312611773 -0.132488283  0.163074792  0.367021974  1.979297106
##  [506]  0.466730394  0.550641286 -0.133652804 -0.536432320  1.218712106
##  [511]  0.603451127  1.012601405 -0.150395941 -0.088956164 -0.704018823
##  [516]  1.213857038  0.327756722 -0.623294659  0.202520383 -0.255022405
##  [521]  0.502185859  0.528736208  1.487557023  0.133806703 -0.385451911
##  [526]  0.173149148 -0.292585139  0.431548778 -0.245763339 -0.195132252
##  [531]  0.085836644  0.187431269  0.506905418 -0.363100721 -0.098740231
##  [536] -0.778357420  0.205420051  1.767309763  0.756619379  0.029342789
##  [541]  0.719562493  0.386295566  0.011136274  1.228080778  1.432905545
##  [546] -0.153924342 -0.076595914 -0.145813955  0.182575566  1.232587507
##  [551]  0.479216699  0.593534962 -0.327026504 -0.503995279  0.072150739
##  [556]  0.059733974  0.542468455 -0.096779323  0.374332178 -0.017767837
##  [561]  0.078155401  0.170780244  0.044317151  0.426805258  0.183022176
##  [566]  0.070152555 -0.013958625  0.034010877  2.104136471 -1.062538956
##  [571]  0.568689780  0.330588087 -1.267561774  1.267539539  0.628468216
##  [576]  0.397477294  0.324279483  0.566842605  0.632448712 -0.170152538
##  [581]  1.301298344  0.319217252  0.772439166 -0.027369715  0.597908732
##  [586]  0.213486395  0.557726860 -0.152661437  1.070910436  0.162799078
##  [591]  0.568253502 -0.071864032  0.537108985  0.265918963  0.554230626
##  [596] -0.011269540 -0.051143215  0.061674517  0.428893089  0.520532653
##  [601] -0.402335685  0.778337645 -0.232626487  0.205886161  0.355510777
##  [606]  0.046598698  0.100988729  0.767276699  0.523767290  0.405600432
##  [611]  0.123914367  0.257663568 -0.465899025 -0.361272827  0.757307471
##  [616] -0.285310572  0.544784347  0.322975613 -0.342387213 -0.534497650
##  [621]  0.545737204  0.566429784 -0.447090685  0.291722918  0.755881183
##  [626]  1.293255012  0.319697687  0.591598995  1.186539437 -0.068813663
##  [631]  0.731894331  1.025570386 -0.719245348 -0.144630151  1.437419388
##  [636]  0.017087438  0.712107965 -0.085219033  1.269045943 -0.372657333
##  [641] -0.090276285 -0.155279139 -0.610056839  0.013241375  0.847478024
##  [646]  0.227121264  0.090229257  0.784308923 -0.109848287 -0.001911470
##  [651]  0.952007215  0.236524894  1.035650998  0.441791887 -0.757996538
##  [656]  0.312669165  0.550490526  0.871332708  0.747175116  0.184366103
##  [661] -0.480438885  0.268594769  0.446599875 -0.398482690  0.564675659
##  [666] -0.465177728  0.545883786  0.561312200  0.732301279  0.556680123
##  [671]  0.077120660  0.477994293  0.724668869 -0.428878662  0.064827320
##  [676] -0.194583299  0.153918780  0.990137558 -0.190052905  0.498540476
##  [681]  0.575573641  0.915013714  0.991454678  0.778301247  0.115330917
##  [686] -0.450098088  0.460091619 -0.705824735 -0.190538622  0.751076831
##  [691]  1.002657393  0.907904236  0.167931397  0.542439687  1.364514174
##  [696] -0.011711062 -0.064106321  0.750703525  0.977615761  0.184916164
##  [701]  0.929614832 -0.348717496  1.581581213  1.220505571  0.719469281
##  [706]  1.018743396  0.378769987 -0.236749018  0.481382661  0.373191242
##  [711]  1.122117732 -0.015372818  0.542396134 -0.163078879  0.085222179
##  [716] -0.677164001  0.554597265 -0.069412495 -0.284836502  1.175881822
##  [721] -0.198651223  0.304958157  0.076407757 -0.385571460  0.747127064
##  [726] -0.057955197  1.249451584 -0.106659896  0.726521333 -0.364741000
##  [731] -0.330667749  0.060529774  1.525421194  0.259010830  0.244161457
##  [736] -1.282838562 -0.197036019 -0.352483748  0.729364943  1.337375238
##  [741]  0.393936687  0.576893402 -0.032507113  0.289297754  0.765687151
##  [746] -0.510707579 -0.127327998  0.200599905  0.941024157  0.471138545
##  [751] -0.174231478  0.068180528  0.517811732 -0.453185311 -0.389479636
##  [756]  0.011281168  0.809270354  0.086202927  0.158130183  0.618169362
##  [761] -0.834304685  0.905715230 -0.123989768  0.850645086  0.181680105
##  [766]  0.222355264  1.303575569 -0.805490327  0.061677996  0.266400742
##  [771] -0.036652810 -0.330201417  0.358212026  0.336256538  0.359023471
##  [776]  0.135236067  0.496309322  1.217796728 -0.038064197  0.319697687
##  [781]  0.278375700 -0.748781694  0.146925255  0.002783500  0.575038523
##  [786]  0.050279846  0.449985626  0.916493387 -0.106491436  0.137232759
##  [791]  0.205614327  0.786335455  0.345821476  0.343142670 -0.829538986
##  [796]  0.462831103  0.178599059 -0.329637393  1.911952941  2.517202674
##  [801]  0.937287243  2.470379376  0.058745186  0.456127024 -0.182455082
##  [806]  0.465049396  0.266059596  0.065393703 -0.189327146  0.788972446
##  [811]  2.513129264  0.355701286  0.230581220  0.542693311  1.625538314
##  [816]  0.562742659 -0.150395941  0.537864757  0.001049221  0.063198645
##  [821] -0.517226208  0.604979581  0.510364574  0.268377016  1.229000106
##  [826]  0.761792971  0.484224638  1.269916692 -0.307182818  0.207640812
##  [831] -0.124428299  0.362592086  2.043757169  0.960966406 -0.575564429
##  [836] -0.147049964 -0.180127099  0.221908180 -0.285518650  0.733869470
##  [841]  0.229868258  0.191207170 -0.553022456  0.386749299  0.695641975
##  [846]  0.778682588  0.489095631  0.294346724  0.229663957 -0.048617775
##  [851]  1.242174906 -0.895078467 -0.378092895 -0.306895380  0.369784963
##  [856]  0.629517062 -0.049372009  0.013183112  0.594117419  0.530785135
##  [861]  0.342618461  0.337913320  0.159307478  1.255678688 -0.631114172
##  [866] -0.128191641 -0.094367148 -0.032984973 -0.028011144  0.522243855
##  [871]  1.368790551 -0.113989693  0.229075689  0.757531063  1.006747182
##  [876] -0.615708215  0.964534697 -0.328181662  0.234708126  0.358539527
##  [881]  0.541806979 -0.052196732  0.797813790  0.065231127  0.750479995
##  [886]  0.753244249 -0.888184293  1.014330529 -0.107062761  0.158066120
##  [891] -0.216117797  1.973627863  1.004675173  0.587596535 -0.242154110
##  [896]  0.240264876 -0.446825655 -0.317412538 -0.248689314  0.156083112
##  [901] -0.173941421 -0.490000283  0.036217526  2.186206728  1.401973274
##  [906]  0.580627337  0.245731456  0.133238286 -0.312604211  0.226586947
##  [911] -0.206528329  0.795892320  0.581354548  0.501517235 -0.318573299
##  [916] -0.717995662 -0.083289372  2.155447424  0.355853491  0.913608076
##  [921]  0.982819160  0.773176363 -0.055248580  0.603961028 -0.491782498
##  [926]  0.995751851  0.989042557 -0.144606364 -0.166944460 -0.087917046
##  [931]  0.564670886  0.123145783  0.441791887 -0.109573373  1.036342379
##  [936]  0.234094329  0.183130500  1.973723112 -0.050202897  0.789087828
##  [941]  1.231216458  0.065454899  2.499691188  0.516300039  0.058909701
##  [946]  0.447749899  0.861206934  0.570603913 -0.394012347  0.545603644
##  [951]  0.354765571  0.703946153 -0.225081097  0.560754575 -0.034860768
##  [956]  0.095318718  0.590700033 -0.134315904 -0.008780392  1.281549322
##  [961]  0.294539900  0.979079773  0.174138268  0.477425353  0.440656625
##  [966]  2.287680450 -0.095586822  0.420447625  0.032542493  1.257427275
##  [971]  0.491605973 -0.142330553  0.575489548  1.272728937  0.560799094
##  [976] -0.536041319  0.123145783  1.696313679  0.985897280 -0.086760456
##  [981]  0.731953460  0.637152772  1.521956833 -0.184720209 -0.685792959
##  [986]  0.754203381  0.027669524  0.766466050  0.734253298 -0.341200366
##  [991]  0.553723177  0.602996052  0.380166593 -0.402056831 -0.008643585
##  [996]  0.019567229  0.137129916  0.535001744  0.368169833  0.994346724
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
##   0.71121431   0.25287757 
##  (0.07996691) (0.05654016)
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
## [1]  0.7729648 -0.3150362  0.2699170  0.9282798  0.2281797 -0.7490065
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
## [1] -0.0376
```

```r
sd(results2$thetastar)  #the standard deviation of the theta stars is the SE of the statistic (in this case, the mean)
```

```
## [1] 0.9075957
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
## t1*      4.5 -0.03853854    0.896512
```

One of the main advantages to the 'boot' package over the 'bootstrap' package is the nicer formatting of the output.

Going back to our original code, lets see how we could reproduce all of these numbers:


```r
table(sample(x,size=length(x),replace=T))
```

```
## 
## 3 5 6 7 8 9 
## 2 1 1 4 1 1
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
## [1] -0.0619
```

```r
se.boot
```

```
## [1] 0.9332788
```

Why do our numbers not agree exactly with those of the boot package? This is because our estimates of bias and standard error are just estimates, and they carry with them their own uncertainties. That is one of the reasons we might bother doing jackknife-after-bootstrap.

The 'boot' package has a LOT of functionality. If we have time, we will come back to some of these more complex functions later in the semester as we cover topics like regression and glm.

