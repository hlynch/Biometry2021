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
## 0 3 4 5 6 7 8 9 
## 1 1 1 3 1 1 1 1
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
## [1] 0.0191
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
## [1] 2.759339
```

```r
UL.boot
```

```
## [1] 6.278861
```

Method #2: Simply take the quantiles of the bootstrap statistics


```r
quantile(xmeans,c(0.025,0.975))
```

```
##  2.5% 97.5% 
##   2.9   6.4
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
##    [1] 3.6 4.7 4.5 5.0 4.6 6.7 3.5 3.8 2.4 4.2 3.9 3.1 4.1 3.7 5.6 4.2 2.3 5.2
##   [19] 3.6 4.0 3.2 4.0 4.4 3.7 4.4 3.8 4.2 3.9 4.5 4.8 4.2 3.0 4.1 4.1 4.0 5.3
##   [37] 4.8 3.4 4.2 5.2 4.9 4.9 3.0 3.3 5.7 5.0 4.5 4.6 3.6 3.9 3.3 2.5 5.7 3.2
##   [55] 4.4 4.8 5.2 5.2 4.1 4.4 5.3 4.5 4.0 4.4 5.2 4.5 5.1 5.1 5.2 4.1 4.6 3.6
##   [73] 4.1 5.2 3.4 4.4 4.0 2.7 3.2 3.9 4.1 4.6 3.4 3.4 3.4 3.8 4.0 4.8 5.1 2.6
##   [91] 4.8 5.2 2.7 4.4 3.7 3.7 4.2 4.7 4.3 2.8 5.4 5.4 4.5 4.8 4.5 4.8 5.1 5.8
##  [109] 4.8 4.8 3.9 3.3 4.3 4.4 4.2 4.0 4.5 3.3 2.7 5.1 3.1 4.6 4.0 6.2 3.8 4.9
##  [127] 5.2 5.6 4.9 4.2 3.2 6.5 2.1 4.3 4.8 4.4 3.3 4.8 4.9 4.5 4.2 2.4 4.9 4.7
##  [145] 5.9 3.5 5.3 3.5 4.3 4.5 4.7 4.1 4.8 2.5 5.6 5.1 4.4 4.4 4.0 5.6 6.8 4.0
##  [163] 3.5 4.5 3.9 4.5 6.4 3.7 5.1 4.4 5.1 3.8 5.4 3.7 5.0 4.8 3.4 4.4 4.8 4.6
##  [181] 4.6 5.0 3.8 4.1 5.1 4.4 5.5 3.7 5.1 4.5 6.4 4.9 4.5 5.3 4.1 4.2 5.2 3.1
##  [199] 4.3 6.7 5.9 5.5 4.2 2.7 4.4 4.7 5.0 4.6 5.3 4.3 3.9 3.5 6.1 4.0 5.7 5.3
##  [217] 4.2 4.0 3.3 5.7 4.9 4.8 4.5 3.9 3.5 4.2 3.4 2.5 5.1 3.6 4.0 5.1 4.1 4.0
##  [235] 4.6 6.7 5.3 4.4 5.5 5.0 5.7 5.1 7.2 5.9 4.2 5.0 3.0 5.4 4.4 4.4 5.1 4.4
##  [253] 6.8 3.8 3.5 5.2 5.3 4.4 7.2 5.4 4.6 5.9 3.9 4.6 4.2 5.0 3.4 5.1 4.6 4.4
##  [271] 3.5 3.0 5.0 4.3 3.5 4.2 4.0 4.1 4.8 4.0 5.3 4.4 5.4 3.9 4.3 3.7 3.6 4.8
##  [289] 4.9 4.5 3.1 5.1 5.1 4.1 4.7 5.8 3.7 5.2 4.9 4.5 5.6 5.3 3.9 5.0 3.7 4.6
##  [307] 3.7 4.5 5.0 6.0 4.6 2.5 3.5 4.6 4.5 4.4 3.9 3.9 4.8 5.2 5.6 5.6 4.1 4.5
##  [325] 6.2 4.5 5.0 3.6 4.1 2.6 3.5 5.1 3.7 4.9 4.2 4.6 4.8 3.1 3.4 3.6 5.2 4.3
##  [343] 4.0 4.6 4.8 4.9 5.6 4.2 4.0 5.1 5.3 3.4 4.1 6.4 7.0 4.1 5.5 3.5 2.4 5.8
##  [361] 4.6 3.7 4.6 3.8 3.3 5.7 4.7 3.0 4.2 6.0 4.4 5.2 3.6 5.5 4.9 3.8 3.7 5.1
##  [379] 5.2 3.7 4.2 6.6 5.1 4.7 4.7 4.8 4.4 5.4 5.8 4.9 5.2 4.6 2.9 3.5 4.8 4.2
##  [397] 4.3 5.8 4.9 5.1 4.0 4.0 5.8 5.0 6.1 3.9 3.1 5.5 4.3 4.3 3.5 5.1 4.2 4.7
##  [415] 5.5 5.0 3.6 3.9 4.8 3.1 3.7 4.5 3.5 5.4 3.3 4.0 4.9 3.7 5.6 3.8 5.1 3.0
##  [433] 5.2 3.5 2.7 3.3 3.1 3.1 4.7 5.2 4.1 4.8 4.1 3.2 2.5 2.7 5.6 4.8 5.0 3.4
##  [451] 4.0 5.6 5.3 3.3 3.6 6.1 5.2 4.6 3.8 5.7 5.6 4.4 6.2 5.8 4.4 4.8 5.5 3.5
##  [469] 5.9 4.3 3.1 5.0 3.3 6.1 4.8 6.4 5.2 4.6 5.6 4.5 5.0 4.2 4.7 3.7 4.8 3.2
##  [487] 4.4 4.2 4.0 3.5 3.3 3.4 3.8 5.6 4.4 4.3 4.5 4.4 5.5 5.4 5.5 3.2 4.6 4.3
##  [505] 3.9 3.0 3.5 5.2 4.9 5.2 4.1 4.4 4.1 3.5 3.7 5.0 4.6 3.8 4.8 3.5 4.2 5.3
##  [523] 4.2 4.5 4.3 5.6 4.2 5.3 3.3 5.0 3.8 2.2 5.3 5.0 3.3 4.9 5.1 5.0 3.5 4.0
##  [541] 4.3 4.3 4.4 4.9 6.5 5.8 4.4 4.4 4.3 4.1 4.0 5.0 4.0 4.4 1.8 5.4 3.8 3.9
##  [559] 4.6 4.7 4.5 3.1 4.1 3.0 2.7 5.5 4.8 3.4 4.8 4.2 5.5 4.2 2.9 2.4 4.7 4.9
##  [577] 3.2 4.1 5.1 5.6 4.2 3.5 5.5 4.8 5.2 6.6 3.9 6.1 5.3 3.0 4.6 5.5 3.7 4.0
##  [595] 5.0 5.6 4.1 4.3 5.2 2.6 5.9 4.2 2.3 3.9 5.1 3.8 4.1 2.8 5.7 4.4 4.3 3.6
##  [613] 4.0 6.0 3.9 2.9 4.7 5.4 5.2 4.4 4.3 5.3 5.0 3.5 4.7 3.9 4.3 3.1 5.6 5.8
##  [631] 4.2 4.6 4.1 2.8 5.6 5.5 4.6 3.5 5.3 6.1 4.2 5.4 4.5 5.1 5.4 5.2 6.1 5.0
##  [649] 6.1 4.7 4.2 4.1 3.5 4.8 4.9 5.6 5.1 4.9 3.3 3.9 4.3 2.7 4.5 5.0 4.2 3.2
##  [667] 5.6 4.5 4.8 4.5 4.5 5.9 4.5 4.3 4.7 5.1 3.7 5.8 6.0 3.5 4.9 5.9 4.1 4.5
##  [685] 5.8 4.4 6.9 5.0 4.2 4.0 3.8 4.2 6.2 5.8 4.3 5.0 3.8 3.9 4.6 5.3 4.6 4.4
##  [703] 5.1 4.5 5.6 4.9 5.3 3.9 3.4 4.9 5.5 4.4 5.1 4.0 4.6 5.9 6.3 5.4 4.8 5.5
##  [721] 4.6 3.9 4.4 4.5 3.8 3.3 4.7 3.6 3.8 3.8 4.4 3.6 4.2 5.7 5.1 3.6 3.7 4.5
##  [739] 3.5 5.3 5.0 3.7 5.6 4.9 5.2 3.5 3.6 4.4 3.2 4.5 5.1 3.5 5.2 4.2 3.4 3.5
##  [757] 4.3 4.0 4.3 4.3 3.4 4.9 5.8 4.1 4.9 3.0 4.1 5.8 4.8 5.3 5.6 4.6 4.6 5.4
##  [775] 5.9 3.4 4.3 5.5 6.3 4.0 4.5 4.3 5.0 4.1 5.2 4.2 5.5 3.7 4.2 5.4 4.4 4.7
##  [793] 4.4 5.3 4.4 4.5 4.5 5.3 5.2 5.1 4.3 3.5 3.9 4.5 4.5 3.5 4.7 5.3 5.1 3.4
##  [811] 3.5 3.7 5.3 4.6 4.7 5.3 3.8 6.2 5.0 4.3 2.9 4.3 4.7 5.0 5.9 3.7 4.1 4.4
##  [829] 5.3 3.8 2.8 4.0 4.7 4.8 3.9 3.3 3.3 4.2 5.1 4.8 3.7 4.2 4.0 4.9 4.6 2.2
##  [847] 4.4 5.3 5.0 3.7 4.9 5.1 3.5 3.5 3.9 4.2 4.2 4.3 5.9 6.5 4.6 4.7 4.5 3.5
##  [865] 5.2 4.3 5.1 4.2 5.1 5.3 3.8 4.4 4.8 3.9 2.7 2.5 3.3 3.8 4.1 2.5 4.7 5.4
##  [883] 3.6 3.8 4.9 4.9 3.4 4.4 2.9 6.3 4.5 2.9 5.5 4.4 4.9 5.5 6.6 4.7 5.8 5.1
##  [901] 5.0 2.8 4.9 4.5 3.5 6.5 4.2 3.0 5.0 3.6 4.5 4.3 5.4 3.5 4.5 5.4 4.0 3.7
##  [919] 3.8 6.0 4.0 5.4 6.4 5.5 4.5 4.2 4.8 4.2 4.0 5.7 3.6 5.4 5.0 5.2 3.2 3.7
##  [937] 2.6 5.1 5.2 5.1 4.0 4.8 4.2 4.9 4.1 4.5 4.8 4.0 4.0 6.8 4.9 5.9 2.5 4.9
##  [955] 3.2 4.5 5.7 4.9 6.6 3.8 3.6 4.5 4.8 3.2 5.3 5.4 4.8 4.4 5.6 3.6 4.6 3.8
##  [973] 4.0 4.7 4.7 5.6 3.0 4.0 5.7 5.9 4.2 3.7 4.6 4.5 5.2 4.1 4.9 3.2 4.4 4.6
##  [991] 4.9 3.9 3.1 4.3 4.3 4.3 5.2 2.8 4.0 4.0
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
##   2.7   6.3
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
##    [1] 4.5 3.8 4.5 5.0 4.7 6.0 4.9 4.3 5.0 4.3 6.6 5.0 3.6 5.2 3.4 3.7 4.5 3.8
##   [19] 3.5 3.4 2.4 3.8 5.5 4.8 2.8 4.6 3.4 5.4 6.9 3.7 6.0 2.4 5.2 3.7 4.5 3.1
##   [37] 2.6 5.8 4.8 4.5 3.1 3.4 4.2 3.6 5.4 5.2 4.4 6.0 5.2 5.3 5.2 3.8 4.9 3.8
##   [55] 5.5 5.0 4.6 4.3 3.8 4.8 5.4 4.6 4.9 4.6 4.5 5.5 4.4 3.7 4.6 5.1 5.9 4.3
##   [73] 4.1 6.6 5.8 4.9 4.0 3.5 5.6 4.3 4.5 5.1 4.7 4.5 6.3 3.1 6.8 3.7 5.9 4.5
##   [91] 4.2 5.4 3.6 4.9 4.5 2.9 5.2 5.0 4.5 5.6 5.2 4.3 4.4 4.0 5.4 4.9 3.2 3.9
##  [109] 5.3 4.0 6.0 4.6 3.7 5.3 5.8 4.2 4.9 4.7 4.8 4.1 5.0 4.4 4.6 2.8 5.1 4.9
##  [127] 5.6 6.5 2.1 3.8 2.2 3.4 4.2 4.3 5.3 5.9 5.6 4.7 5.3 5.3 4.8 2.8 5.1 4.4
##  [145] 3.6 3.6 3.2 4.4 5.1 3.1 6.0 5.6 4.2 5.8 3.8 4.0 4.6 4.6 3.7 4.3 3.1 3.7
##  [163] 5.5 5.3 5.1 5.0 4.7 4.5 4.6 5.0 5.1 3.8 4.4 4.6 4.8 4.6 4.9 5.1 5.6 3.8
##  [181] 3.6 4.8 3.3 4.6 7.5 4.5 3.4 3.1 5.7 4.5 4.5 5.8 4.1 3.8 2.6 4.6 5.0 3.7
##  [199] 4.4 6.0 4.0 5.1 2.6 4.7 3.9 5.1 5.8 3.2 4.7 4.2 5.8 6.0 4.8 3.9 3.3 4.1
##  [217] 4.0 4.0 2.0 4.7 3.8 5.5 4.6 3.7 4.8 5.3 6.7 5.2 4.7 3.7 3.7 4.6 3.3 4.2
##  [235] 4.2 3.8 4.0 4.5 4.2 3.0 5.3 6.2 4.5 2.1 4.3 3.7 5.8 4.6 5.7 4.4 2.7 5.5
##  [253] 3.7 3.8 4.8 5.6 5.9 3.0 4.7 4.9 4.6 4.4 5.1 4.2 4.7 5.4 4.1 6.3 4.8 4.5
##  [271] 3.7 3.0 4.2 5.2 4.8 4.1 4.1 4.5 4.2 4.6 4.0 5.4 4.6 4.0 3.1 4.7 3.2 2.8
##  [289] 4.0 3.7 4.8 2.6 4.7 3.6 5.2 5.5 4.9 3.9 5.2 3.9 4.8 4.4 4.9 5.7 4.7 4.2
##  [307] 3.5 5.2 5.1 4.9 4.9 4.0 4.3 4.8 4.9 5.1 4.5 5.5 5.6 5.1 3.8 2.6 3.6 3.8
##  [325] 3.9 5.6 5.2 5.0 5.5 4.4 4.2 5.3 3.4 3.0 5.0 2.9 5.5 4.7 6.2 6.2 2.5 5.4
##  [343] 3.8 3.2 3.2 5.5 4.5 4.9 4.6 5.6 5.0 4.3 4.6 3.7 4.7 5.1 3.7 3.7 5.5 5.7
##  [361] 4.4 4.9 4.6 3.2 4.0 4.3 3.7 6.1 5.6 4.6 5.5 4.9 4.7 2.8 5.2 4.5 4.4 4.7
##  [379] 5.9 4.8 2.5 4.8 4.1 4.2 5.1 4.3 4.6 4.8 5.8 4.7 5.5 4.4 2.9 4.5 3.5 4.4
##  [397] 3.0 4.1 5.2 6.0 5.0 5.0 4.0 6.0 4.3 6.0 3.3 3.8 4.3 5.7 3.9 5.5 4.5 6.3
##  [415] 4.0 3.3 4.3 4.7 5.7 4.6 5.1 2.7 5.0 5.2 5.6 2.4 2.8 5.6 5.0 4.2 3.6 3.9
##  [433] 4.1 4.3 3.2 3.5 4.2 5.7 5.7 4.6 5.6 3.3 4.2 4.2 4.0 3.6 4.7 3.8 4.7 5.4
##  [451] 4.1 5.6 2.1 5.0 4.9 5.3 5.3 3.9 5.7 6.4 4.2 4.1 4.7 4.4 4.0 4.5 5.4 6.3
##  [469] 4.5 5.4 4.5 4.5 4.9 3.1 4.0 4.4 5.0 4.5 5.2 3.0 4.1 5.0 4.9 6.3 5.1 4.4
##  [487] 4.7 5.5 4.7 3.4 4.9 5.2 4.4 3.6 4.4 5.1 5.2 4.1 4.9 4.1 4.3 5.0 3.6 3.1
##  [505] 5.3 5.1 4.5 4.8 3.2 5.4 4.6 3.1 4.0 4.3 4.4 2.9 5.6 4.2 3.1 4.5 5.6 5.2
##  [523] 3.8 4.6 5.9 4.3 5.0 3.9 3.1 3.7 4.5 5.0 3.1 3.3 4.7 3.3 5.1 4.0 5.2 5.3
##  [541] 4.5 5.4 4.4 4.2 4.3 3.1 5.3 5.8 5.0 3.8 4.7 4.8 3.1 5.0 4.4 4.7 4.1 4.0
##  [559] 4.6 4.0 5.8 5.5 4.7 3.9 5.1 4.9 4.6 3.8 4.5 2.4 5.3 4.3 5.5 3.5 4.9 4.3
##  [577] 3.5 4.3 2.9 3.6 5.1 4.8 3.3 3.5 3.7 3.4 3.3 3.7 5.2 3.4 5.8 4.5 4.6 6.0
##  [595] 4.1 4.4 5.4 4.1 3.5 3.5 3.0 3.5 3.0 4.0 2.8 4.8 3.8 4.3 5.1 4.1 3.7 4.9
##  [613] 5.0 3.5 6.4 3.8 4.7 6.0 3.8 5.0 3.1 3.7 4.4 4.8 3.9 4.8 5.9 2.6 4.9 5.2
##  [631] 5.2 3.0 5.3 4.6 3.4 3.5 5.9 4.1 5.0 4.0 5.0 4.2 3.8 4.2 4.5 5.6 4.0 6.6
##  [649] 4.2 3.8 3.6 5.5 5.1 6.4 3.3 3.5 4.5 5.2 3.5 2.4 6.7 3.7 4.7 3.6 4.7 2.9
##  [667] 3.2 5.9 2.8 5.4 4.1 5.4 4.1 5.6 5.1 3.2 3.5 6.3 5.3 2.8 3.5 4.7 4.9 5.0
##  [685] 3.4 5.1 4.5 4.3 5.3 2.8 4.2 5.6 3.7 3.8 2.3 5.3 5.1 5.1 3.6 4.6 4.1 3.5
##  [703] 5.1 5.2 5.0 5.1 6.6 4.8 4.4 4.6 3.8 3.6 4.3 5.1 2.9 4.1 4.7 5.2 4.8 6.1
##  [721] 5.8 3.6 3.6 6.5 4.2 4.9 2.9 3.3 6.5 3.4 4.1 3.9 5.2 4.7 3.9 4.3 3.6 4.6
##  [739] 4.2 5.2 5.6 3.7 4.9 5.4 3.9 4.4 5.8 4.2 4.8 3.6 4.8 4.0 4.8 4.2 4.4 5.1
##  [757] 4.6 3.8 5.8 6.1 4.5 4.0 4.6 3.6 3.7 3.3 4.1 5.3 3.6 4.4 4.5 5.0 6.1 3.2
##  [775] 4.8 4.0 4.1 4.5 4.1 4.4 4.2 5.5 6.6 4.5 6.7 4.8 6.2 4.5 5.0 3.9 4.7 5.1
##  [793] 4.2 3.5 3.9 5.2 4.2 3.0 3.7 3.9 3.0 4.1 3.6 3.3 5.4 3.4 5.7 5.2 3.3 4.0
##  [811] 5.0 5.3 6.0 4.8 5.4 3.4 5.8 3.1 3.8 5.0 6.9 4.8 4.7 4.3 5.0 3.4 3.4 4.9
##  [829] 4.6 4.1 3.4 4.6 6.2 4.1 5.4 3.2 4.6 6.0 4.6 3.5 3.8 3.5 4.4 2.9 4.7 3.3
##  [847] 5.4 4.1 3.4 4.2 4.2 3.6 5.5 3.6 5.6 5.7 5.0 4.4 3.7 5.2 4.0 3.8 6.1 5.0
##  [865] 5.1 2.1 5.4 5.9 4.7 4.2 5.1 3.5 3.7 4.3 4.4 5.4 4.9 4.8 6.0 3.5 4.7 4.0
##  [883] 4.2 6.1 4.6 4.0 5.4 3.2 5.0 4.7 3.9 4.6 5.5 6.3 6.5 3.7 4.6 4.7 3.4 4.9
##  [901] 6.2 5.2 4.4 4.7 3.8 4.4 3.2 3.4 5.2 2.4 5.5 6.0 3.4 3.4 3.4 4.6 5.5 3.3
##  [919] 4.8 4.5 6.5 4.6 6.0 4.8 4.3 5.9 5.3 4.1 5.9 4.6 4.3 5.5 4.6 4.4 5.4 5.1
##  [937] 6.0 5.5 5.9 4.5 4.5 4.1 5.7 4.4 3.3 2.5 4.2 4.9 5.5 4.1 5.7 4.2 2.9 4.4
##  [955] 5.3 3.6 5.2 5.2 3.2 4.9 5.9 4.6 5.3 4.4 4.7 5.2 5.1 4.9 4.1 4.2 3.2 6.3
##  [973] 4.0 4.9 4.3 4.3 3.8 2.7 4.5 3.9 3.4 3.2 6.6 3.4 5.3 4.8 5.7 3.4 4.5 2.8
##  [991] 4.0 2.7 4.6 5.0 4.7 4.5 5.3 5.0 5.9 4.2
## 
## $func.thetastar
## [1] 6e-04
## 
## $jack.boot.val
##  [1]  0.51306818  0.35915493  0.30189189  0.22824207  0.11719745 -0.08978979
##  [7] -0.17843137 -0.31895604 -0.39517045 -0.50797721
## 
## $jack.boot.se
## [1] 0.9972226
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
##    [1] 4.8 6.0 4.3 3.1 4.0 4.7 4.3 4.1 6.1 5.3 5.0 3.7 4.9 3.7 4.3 5.2 5.6 4.9
##   [19] 3.7 5.2 3.7 4.4 4.3 3.4 2.8 4.8 2.8 4.2 4.7 3.1 3.9 3.5 4.0 4.1 5.9 4.5
##   [37] 5.5 4.0 3.3 3.1 4.3 4.1 4.3 4.1 3.8 3.9 4.1 5.3 5.0 4.3 4.7 4.4 5.9 4.7
##   [55] 3.9 5.1 3.2 4.0 5.0 3.6 4.9 1.9 5.7 3.7 3.2 6.0 6.3 3.2 4.4 4.0 4.5 3.6
##   [73] 5.3 3.3 5.5 3.5 4.0 4.3 3.3 5.5 4.5 4.1 4.2 3.1 6.2 5.5 3.9 3.7 4.2 4.9
##   [91] 4.0 6.3 5.9 5.8 3.5 2.8 5.2 4.4 3.5 4.9 3.2 4.7 5.4 3.8 3.9 4.3 4.7 3.6
##  [109] 5.6 5.1 5.9 4.0 3.2 5.7 6.2 3.8 4.1 3.8 4.3 4.7 5.4 5.5 4.0 4.0 4.4 4.3
##  [127] 4.4 4.3 3.7 4.7 2.6 6.1 3.7 4.3 4.3 4.4 4.6 3.9 4.5 3.8 5.6 5.3 4.7 3.4
##  [145] 4.5 4.0 4.2 5.8 4.8 4.6 6.9 3.4 5.0 4.2 4.7 4.4 4.1 4.0 4.5 3.2 3.4 5.4
##  [163] 5.0 6.3 5.5 4.7 4.1 3.0 4.6 6.3 4.3 4.4 6.3 3.6 3.8 4.8 5.7 2.6 4.9 4.0
##  [181] 4.8 5.4 5.6 3.9 3.6 4.0 5.6 3.7 5.6 2.6 2.9 4.3 3.9 4.5 4.0 3.5 5.0 4.5
##  [199] 5.2 4.3 4.4 5.6 3.8 3.5 4.2 6.0 3.9 5.0 4.2 5.7 5.5 3.9 3.9 3.6 3.9 5.1
##  [217] 5.0 3.8 5.4 2.5 4.0 4.9 5.2 3.8 4.7 4.2 4.6 5.5 5.6 4.2 4.8 5.5 5.1 5.1
##  [235] 3.4 5.8 4.6 5.9 4.2 6.8 4.9 4.9 4.9 3.0 4.4 4.6 3.3 5.1 4.4 5.0 2.9 5.9
##  [253] 5.1 5.7 4.9 5.6 5.1 2.9 3.7 6.5 4.8 2.9 3.9 4.5 2.9 5.9 3.4 4.8 4.2 2.9
##  [271] 3.9 4.6 5.2 3.3 4.6 4.6 5.5 4.5 5.8 5.6 4.0 3.9 4.4 5.3 3.5 4.9 4.7 4.6
##  [289] 5.2 3.5 3.9 3.9 6.9 5.0 4.8 4.4 4.2 3.9 5.0 2.9 3.9 3.5 4.2 4.6 4.4 5.0
##  [307] 4.2 2.6 3.4 6.2 3.2 4.2 5.3 4.6 4.0 4.9 5.0 4.6 4.8 5.2 5.2 4.3 3.9 6.3
##  [325] 5.7 6.2 3.6 5.5 3.6 4.3 4.1 4.3 5.2 5.5 4.9 5.2 4.0 5.4 3.8 3.8 7.5 2.9
##  [343] 4.2 3.9 3.4 3.1 5.0 4.3 4.6 4.6 4.6 5.0 5.5 3.8 5.1 5.0 5.3 4.3 4.5 5.5
##  [361] 3.6 5.2 5.2 3.5 5.5 4.5 4.5 3.9 3.2 3.6 6.1 4.3 3.0 5.5 3.9 5.4 5.1 5.5
##  [379] 5.0 5.8 3.8 3.8 3.6 4.6 4.2 6.8 3.3 5.6 6.0 3.3 5.1 4.0 3.0 4.4 4.3 3.0
##  [397] 4.9 4.0 4.8 3.9 3.6 5.0 2.3 5.0 3.7 3.9 4.1 4.7 4.9 4.3 4.2 2.9 5.0 3.7
##  [415] 3.7 3.8 4.2 3.1 4.3 4.1 3.7 3.0 5.1 4.4 3.6 4.7 2.5 3.8 4.5 3.0 4.3 3.3
##  [433] 5.1 4.6 5.3 4.3 4.3 4.3 4.0 4.4 5.5 7.1 5.0 4.6 5.0 6.7 3.6 3.8 5.2 5.6
##  [451] 6.1 4.2 4.6 3.0 5.1 5.0 4.0 4.2 3.6 2.5 4.1 4.4 4.6 4.4 5.7 5.9 4.5 4.3
##  [469] 5.0 4.7 5.3 4.2 3.6 4.8 4.2 2.9 3.9 3.9 3.8 3.8 5.2 6.1 4.6 4.9 3.8 4.7
##  [487] 4.2 3.9 4.8 3.5 4.8 4.9 4.2 3.2 4.1 4.0 5.4 4.3 4.7 5.3 4.2 5.6 3.0 4.6
##  [505] 6.2 4.9 4.8 3.7 4.8 3.1 5.6 5.6 4.9 5.8 4.9 5.2 4.6 6.3 4.5 3.7 5.2 4.2
##  [523] 5.5 3.5 5.9 4.0 5.3 4.4 4.8 4.3 5.1 4.7 5.8 5.5 3.3 5.4 6.5 4.4 5.1 5.2
##  [541] 5.5 5.5 4.2 4.6 6.7 3.3 5.4 4.4 3.9 5.3 5.9 4.7 6.0 5.3 5.5 5.9 3.5 6.4
##  [559] 3.5 3.5 2.3 2.4 3.8 4.5 2.9 4.4 4.2 5.0 7.3 3.3 5.2 5.4 5.3 3.9 5.8 4.4
##  [577] 3.2 5.4 2.7 4.8 6.0 5.4 5.7 4.3 5.1 5.6 3.5 3.9 4.8 4.7 5.1 4.4 3.9 5.0
##  [595] 5.5 4.1 4.4 3.2 3.8 3.8 4.2 5.3 3.5 3.5 5.0 4.9 5.1 4.9 5.3 4.4 3.7 4.3
##  [613] 4.7 3.6 5.8 4.3 4.5 5.3 4.1 5.2 4.5 4.4 5.1 2.7 5.3 4.4 4.8 5.4 4.3 5.5
##  [631] 3.4 3.9 6.8 2.7 4.6 4.5 4.4 5.3 3.5 4.5 5.5 5.5 4.3 5.4 4.5 4.4 5.3 5.0
##  [649] 3.5 4.3 6.2 3.7 3.6 5.0 4.0 5.9 4.7 5.7 4.4 3.0 4.3 3.8 4.3 2.9 4.0 3.3
##  [667] 3.4 4.2 6.5 4.6 4.9 6.5 5.9 3.4 3.9 4.4 3.4 5.7 5.2 4.9 4.2 4.1 5.4 3.0
##  [685] 5.0 3.8 3.6 5.4 4.8 5.8 4.6 3.8 4.7 4.9 4.0 3.8 5.0 4.1 5.1 4.5 3.0 5.2
##  [703] 4.7 4.1 4.4 4.8 4.3 3.7 4.3 5.9 4.5 4.6 4.9 5.4 6.3 4.5 6.3 4.5 4.7 5.8
##  [721] 4.0 4.6 5.5 3.9 4.8 5.6 6.9 5.3 3.8 4.4 3.9 4.9 4.0 4.1 4.1 3.8 4.6 4.1
##  [739] 4.8 4.5 4.3 3.9 4.1 5.4 3.4 5.3 5.6 3.1 4.6 4.8 5.1 5.8 3.9 5.3 4.5 4.1
##  [757] 3.5 4.1 2.0 3.6 3.5 4.2 4.2 6.0 4.6 2.5 3.4 3.9 5.1 3.8 3.7 4.4 5.0 4.1
##  [775] 5.2 5.4 5.6 4.9 6.0 5.6 5.2 4.5 3.8 3.9 2.3 6.3 3.3 4.2 4.9 4.6 5.5 4.6
##  [793] 5.9 5.5 3.6 4.2 3.7 5.2 4.7 5.6 5.6 5.2 4.8 6.2 3.7 4.3 3.3 2.8 4.8 6.6
##  [811] 4.0 7.1 5.6 4.2 4.1 5.1 4.4 4.3 3.7 4.4 4.2 5.1 5.8 5.2 5.6 4.4 4.5 4.3
##  [829] 3.8 2.5 4.2 4.3 3.3 1.8 2.0 4.2 4.0 5.0 5.2 3.6 4.0 4.6 5.1 5.4 5.1 5.5
##  [847] 4.2 3.9 4.6 5.1 5.4 5.3 4.2 6.1 4.1 5.7 5.0 4.1 4.7 4.1 4.3 6.6 5.5 4.4
##  [865] 2.6 3.6 2.5 3.8 4.8 4.8 5.7 4.0 5.7 5.2 5.6 3.3 4.9 5.6 5.0 4.6 3.8 4.5
##  [883] 3.8 5.0 2.3 6.3 3.9 7.0 4.1 3.0 3.7 5.9 5.0 4.3 4.5 4.0 3.9 4.3 6.1 6.3
##  [901] 4.7 3.9 5.2 5.5 4.4 3.5 3.6 4.0 4.3 4.3 5.4 3.8 5.9 5.9 5.0 4.8 4.6 5.6
##  [919] 4.6 5.5 3.9 6.4 3.6 4.2 4.1 3.7 5.1 4.9 5.2 5.0 3.7 4.7 4.9 4.4 3.6 3.8
##  [937] 4.6 3.6 5.2 5.1 4.3 4.2 4.6 6.8 4.0 5.3 5.2 5.8 6.4 6.2 5.1 4.9 3.6 5.6
##  [955] 4.8 3.4 4.0 2.7 5.7 5.6 4.7 4.1 4.3 5.7 5.1 6.6 3.7 3.0 3.6 5.4 4.3 5.4
##  [973] 3.9 4.7 5.2 4.4 4.9 5.5 4.4 4.5 5.0 3.4 4.9 4.2 4.0 5.4 4.2 2.2 5.1 3.7
##  [991] 5.7 4.0 5.7 4.5 4.3 5.0 4.2 4.8 3.9 4.5
## 
## $func.thetastar
## 72% 
## 5.1 
## 
## $jack.boot.val
##  [1] 5.600 5.500 5.500 5.300 5.100 5.008 5.000 4.900 4.600 4.500
## 
## $jack.boot.se
## [1] 1.072666
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
## [1] 1.766606
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
##   2.982242   4.537271 
##  (1.266213) (2.097984)
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
## [1]  1.6883233 -0.1970964 -0.6205984  0.3270388 -0.1364516 -1.0392873
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
##    [1]  0.102690564  0.575000877  1.087195297  1.166558652  2.189715441
##    [6]  1.963059939  1.946161204  1.121806123  0.271697239  1.737292271
##   [11]  1.266356600 -0.022339860  1.162417935  1.151952753 -1.406058997
##   [16]  1.047499968  1.123791073  1.184826791  1.276278270 -0.274806923
##   [21]  1.787411952  1.584037450  2.124466096  1.813548750  1.815489077
##   [26] -0.228805115  0.711244052  0.774629416  1.535538425 -0.185797095
##   [31]  1.092547317  1.823685678 -1.207033182 -0.286336570  1.236118132
##   [36]  0.867751389 -0.768335786  1.951627275  0.090828149 -0.286770659
##   [41]  1.861831090  1.352758317 -1.215905673  0.059313218 -0.381703586
##   [46]  2.202974109  1.879251779  1.147298429  1.055053113  1.179530271
##   [51]  2.074906272  2.202301312 -0.807067987  1.813910751 -0.288392902
##   [56]  1.161339882 -0.161632852  1.366370670 -0.204532927  0.788382353
##   [61]  1.868378045  1.529735262 -0.306049989  1.201578992  1.223722608
##   [66]  1.166774221  1.132405901  0.066123864  0.074286356  1.648739617
##   [71]  1.775401233  1.634358772  0.632319088  1.931288554  1.937575078
##   [76] -0.410981452  1.124360563 -0.487140308  0.023802296  1.169501335
##   [81]  2.037344162 -0.145612160  0.892090046  0.123472168  1.740477670
##   [86]  1.989182858  0.419553136  0.272304407  1.085283516 -0.626452686
##   [91]  1.817896505  1.130508442  1.783407198  0.023360053 -0.337236020
##   [96]  0.394826823  1.988429552  0.546539232  0.053330020  1.696735952
##  [101]  1.210644895  1.620417249  1.589667804  1.774626462  1.926873920
##  [106] -0.861790907  1.148931919  1.800925334  1.816408242  2.075013381
##  [111]  1.138049915 -0.010434217  0.110728623  0.105146243 -0.581969975
##  [116]  1.854014927  2.174715058 -0.210268773  1.060188232  1.132234701
##  [121]  1.741940248  1.315699109  0.253953491  2.003671392 -0.273292174
##  [126]  1.077772228  0.069868530  1.974304014  1.058394975  0.088517030
##  [131]  1.884225510  1.955781059  1.829631513  1.895027997  0.261448548
##  [136]  0.464615492  1.637042591  0.242868232  1.602826235  1.210854630
##  [141]  0.555836865  1.729962201  1.190452431  1.712786754  0.647310348
##  [146]  2.414070571  1.714956650 -0.784921827 -0.769761567 -1.119480595
##  [151]  1.947404874  1.147040938  1.160006778  1.119849902  1.856915897
##  [156]  1.073522979  1.632365073 -0.293174427  1.730786786  1.621747000
##  [161]  1.444331600 -1.265063769  0.546539232  0.170781959  1.874649288
##  [166] -0.617426784  1.869091267  1.890946774 -0.224316257 -0.280237563
##  [171]  1.145810367  1.760407329 -0.285960336  1.780713395  1.704205416
##  [176]  1.188016259  1.930368187  1.722549453 -1.002198068 -0.660931958
##  [181]  1.918929814  1.902057157  1.863963529  0.127996334  1.743154170
##  [186]  0.498830492  1.990126002  1.093850716  0.032820829  1.176585099
##  [191] -0.531170169  0.382469572  0.365399085  1.939692367  1.627505080
##  [196] -0.706207139  1.822168606  1.852945154 -0.643923403 -1.476092793
##  [201]  1.475528529  1.687803957  1.513521861  0.578175406  2.045483038
##  [206]  1.530013625  1.160618774 -0.685875893  0.448823598  1.247933440
##  [211]  0.394712303  0.605104163  1.744697361  0.021073476 -0.229003295
##  [216]  0.121781037  0.790094622  1.113079453 -1.115748352  0.149035178
##  [221] -0.044104988  1.109381580  1.509155628  1.868563555  1.052969500
##  [226]  0.433009302  0.581069345 -0.686748168  1.719326723  1.693009939
##  [231]  1.786607970 -0.577623237  0.105146243  0.268922763  1.036509911
##  [236]  1.882364209 -0.294047995  1.745803988 -1.339502662  1.878437468
##  [241]  1.219445810  1.770248693  1.770103085  1.743570603  0.066995357
##  [246]  0.072477690 -0.210700785 -0.693429873 -1.238009700  1.173693601
##  [251]  1.696533715  1.579757728  0.591482014  0.376756605 -0.400663762
##  [256]  0.182667608  2.007649497  1.874908704  2.175193223  1.683419701
##  [261]  1.327661242 -0.280362790  0.690591498  1.631840797  1.759814757
##  [266]  2.080119949  1.756692949  1.998186787  0.115328292  0.783774820
##  [271]  1.060050984  0.111859356  1.221215198  1.955031303  2.198546134
##  [276]  1.099828642  0.406340587  2.090186409  1.078292232 -0.313295176
##  [281]  1.948140194  1.806708994  1.191328065  1.994928517  1.845919566
##  [286]  1.154249056  1.126135410  1.717784689  1.271258657  0.121781037
##  [291]  1.267639941  1.626124332  0.706403228  1.748733837  1.219811324
##  [296]  2.027020352  2.403289620 -0.653281339 -0.711432388 -0.177042429
##  [301]  0.836522458  1.838555949  0.660182042 -0.862465496  2.166685576
##  [306]  0.358520391  1.885586045  2.116175378  0.688678739  1.217463700
##  [311]  0.994718569 -0.001365253  1.964019917  0.375186943 -0.624548589
##  [316] -1.094721179  0.516567644  0.716374607 -0.392504832 -0.317734181
##  [321]  1.923738940 -0.334612868  2.010158030  1.888136934  1.419060697
##  [326] -0.251256926  1.593312170  0.524755161  0.107991797  1.800567613
##  [331]  1.207392705  0.349002705 -0.621823777  1.227445197  1.532689678
##  [336]  1.821515287  1.850916790 -0.207466905  1.163196147  1.247793437
##  [341]  1.296726541  1.701674026  1.698226099  1.416616004 -0.623800150
##  [346]  0.376468781  0.677530419  1.194988280 -0.317864176  1.209459650
##  [351]  1.856215164 -0.738274182 -0.308278615 -0.421610962 -1.250942567
##  [356]  0.749397704  1.268132666 -0.258519857  1.791002002  1.047982771
##  [361]  1.325481592  0.126021506  1.554880925 -0.285032374 -0.033021530
##  [366]  1.910751633  1.053830039  1.508845520  1.827752538  1.868716075
##  [371]  0.520768126 -0.680004388  2.150312719  2.014676234  0.482214499
##  [376]  1.768356025 -0.332573293  1.859965393  2.062340351  0.715169854
##  [381]  1.028170378 -0.625699974  1.926377726  0.972523455  1.134573791
##  [386]  1.087927829  1.872506187 -0.057113335  0.523008224  0.615126922
##  [391]  1.992219157  1.474840347  1.123610493  1.933170978  1.771155880
##  [396]  0.599282509  0.922586240  1.907307220  1.808189830  1.913217507
##  [401]  1.494551458  1.672969732  0.063100623 -0.024161840 -0.806389911
##  [406]  1.895606271  0.016561963  1.187375889  1.143132608  1.575773086
##  [411]  1.157954033  1.008703772 -1.135427752  0.219075428  0.599876468
##  [416]  1.243494725 -0.317088766  1.727316163  2.023394686 -0.379349207
##  [421]  1.942087555  1.215726234  0.101316783  0.421236187 -0.335761646
##  [426]  0.863279518  1.923863428  1.858575299  2.082221187  0.869504197
##  [431]  0.378304938 -0.030507112  0.027253756  1.162675191  1.728429575
##  [436]  1.038420917  0.725770089  1.131436517  1.204099460  1.131651201
##  [441]  1.895791282  0.855073680 -0.276833346  1.722986657  1.921304967
##  [446]  1.617765172  1.067562569  1.816262311  1.134573791  1.866816143
##  [451]  1.233104383  1.943133228  0.722142939  0.668149262  0.353713539
##  [456] -0.565874526  1.664640273  1.119753666  1.787647315  0.072341725
##  [461]  0.592724077  0.609659187  1.194113861 -0.073704598  1.880993752
##  [466]  0.448175040 -0.368709560 -0.234791730  0.545581768  1.888232899
##  [471]  0.071812332  1.679303622  1.171126942  0.076070917  1.756005354
##  [476] -0.379416318  2.053716655  1.275251750  1.190167710  1.333741052
##  [481]  1.117442417  0.610703119  1.777197129 -0.773277066  1.655585496
##  [486]  1.893937872  1.069383488  0.479393799  1.113334771  1.180134463
##  [491] -0.296362063  1.589154311  1.619629950  1.741431177  1.317258190
##  [496] -1.072945968  0.376756605  0.182754177  1.304479690 -0.519202257
##  [501]  0.010599350  1.558899366  1.795344780  1.727539034  1.758985888
##  [506]  1.225911857  0.408035934  1.120492084  1.798882377  2.082130286
##  [511]  1.821062070  1.124046184  0.212628720  1.210644895  1.141014627
##  [516]  1.142870415  2.035578467  2.090354291  0.643045854  2.024444742
##  [521] -1.112106924  1.658184225  0.052273671  1.186311012  1.587305791
##  [526]  1.156660468  0.685276210  1.224631682  0.432768133  1.903348466
##  [531]  2.065041600  1.224859483  1.924838260 -0.624585614 -0.188531642
##  [536]  1.335520843  1.126087048  2.053745081  0.105786074  0.896616529
##  [541]  1.177406547  1.700595945  1.715574244  1.686311026  1.986488158
##  [546]  1.450909563  0.696280980 -0.364042378 -0.133226771  1.531981415
##  [551]  1.269523871  1.025580435  0.015864137  1.211844919  1.094867520
##  [556]  1.007117813  1.766089826 -1.213920502  0.695050083 -0.695918782
##  [561]  1.939358484 -1.147104919  1.819462773  0.029529723  0.434559279
##  [566]  0.994718569 -0.534640167 -0.866148321 -0.643192626  0.447883749
##  [571]  1.894109721  0.852980462  0.865827266  0.675770437  1.747733319
##  [576]  0.788587725  2.108283538  1.780877380  1.799714630  0.313431582
##  [581]  0.046565820  1.545750634  1.773318807  1.735465067  2.178382216
##  [586] -0.125961414  1.238917237  0.782015331  0.421098070 -0.197882482
##  [591]  1.715916101  1.594433648 -1.080935307 -1.139251664  0.131371543
##  [596]  1.243160226  1.911065769  1.767217379 -0.155741489  1.645756770
##  [601] -0.715651613  1.509322625  1.780482432 -1.107472053  0.151695040
##  [606]  1.080030004  1.789717056  2.285818664  0.429771353 -0.173642855
##  [611]  1.112940534  1.165366886  1.156135924  0.078078089  1.284609721
##  [616]  1.938213031  0.249911169  0.267720114  0.390404694  1.905541197
##  [621] -0.785544307  0.103864226 -0.698785647 -1.040072025  0.889501858
##  [626]  1.143589555  1.246531282  2.000560783  1.117960328  1.690612941
##  [631]  0.914909056  1.769966589  1.194659291  1.689824635  1.489647932
##  [636]  1.469582279  1.030153516  1.824944610 -0.253073401  0.477671174
##  [641]  0.544121321 -0.207660140  1.991383196  2.318463268  1.195566353
##  [646]  2.061924737  1.857451147  0.830573026  2.029588208  1.426884245
##  [651]  1.290469417  1.192974387  1.912342162 -0.795980235  1.162049466
##  [656] -1.350037514  1.851087227 -1.071784158  1.932668796  2.035796813
##  [661]  0.957463860  1.094034341  1.163340488  1.215364755  1.624515503
##  [666]  1.259367608  0.016544673  1.141351467  1.849156991 -0.250000091
##  [671]  1.830495759  1.734672063  1.864223029 -0.493706647  1.606765564
##  [676]  1.126054998  1.877037690  0.677700674  2.055410314  1.057211212
##  [681]  1.164155606  1.818909472  1.318784172  0.411788899  0.084608718
##  [686]  1.654800635  0.853175080  1.867998031  0.623194998  0.060283685
##  [691]  1.047310793  1.060086101  1.847292559  0.751970808  1.787371554
##  [696]  0.938240496  2.060152991  1.138062274  1.244917676 -0.354026167
##  [701]  2.110219361  1.913537532  0.064002101  2.268559334 -1.323355252
##  [706]  1.638058120  2.020424578 -0.976541994  0.450031654  1.717760754
##  [711]  1.646208170  1.110305026 -0.376751783 -0.546806339  1.805716020
##  [716]  1.807137703  1.865172466  1.964019917  1.193537326  0.576870333
##  [721] -0.046899596  0.116573442  2.113017068  0.706131601  1.576175602
##  [726]  1.419060697  1.809289232  1.196787837  1.950634919  1.390200924
##  [731] -0.665806399 -0.271274101  1.925796216 -0.288324502  1.947458941
##  [736]  1.802211530  0.146271101  1.820175791  1.144221805  1.102781510
##  [741]  0.782939306  2.094803773  1.775595004 -0.319146709  0.786441878
##  [746]  2.081402362 -0.255718906  1.890697540  0.825403306  1.836000968
##  [751]  1.086461775  1.828624523  1.978150184 -0.506296381  1.865172466
##  [756]  0.427727104  1.264588209  1.839287016  0.323847022  1.519247839
##  [761]  0.564558393  2.048159116  1.225911857  1.790411537  0.661848275
##  [766]  0.635623793  1.792404492  1.740960375 -1.382728080  0.522501015
##  [771]  2.014676234  1.942801345  0.734806791  1.834149985  1.122649031
##  [776] -0.629132292  0.064149600 -0.277086754 -0.380533044 -0.340904186
##  [781] -0.440032762  1.885030490  0.639212231  1.164747821 -0.646202669
##  [786]  0.549913119  1.622812847  0.323914062  1.724504001  2.082342009
##  [791]  0.695148018  0.440884679  1.780988385  0.579232234  0.031372665
##  [796]  1.847292559  1.623679246 -0.626737664  0.998316407  0.940521647
##  [801] -0.331180248  1.192404079  0.344683408  1.734303181  1.822439295
##  [806]  1.176492770  1.031860261  1.707583086  0.765316665 -0.190122070
##  [811]  0.644631904  0.045773897  1.167000276  1.815815444 -0.690959095
##  [816]  1.736033041  0.248187182  1.694226907 -1.053336391  2.000662324
##  [821]  1.964621686  2.057650469  1.896040799 -0.622084358 -0.017515250
##  [826]  1.278758541  1.928048705  1.256877559 -0.328624875  0.085450274
##  [831]  1.413877580  1.186761785 -1.273750481  0.545341539 -0.293541421
##  [836]  1.810631823  1.919983138  1.552523377  1.709797476  0.389752227
##  [841]  1.072094325 -0.596462599  1.744631781 -0.683421230 -0.430115562
##  [846] -0.348342141 -0.417161747 -1.180573443  1.852942982  0.669685260
##  [851] -0.857973743  1.842047363  1.790808273  1.981097951  1.909959810
##  [856]  0.495761154  1.153505660  0.437047062  1.046792089  1.236995624
##  [861]  1.069174850  0.744022452  2.078624347  1.758182490  2.140365412
##  [866]  1.833941085 -0.201217859  1.885030490  2.536597527  1.340335122
##  [871] -0.296019690 -1.058011266  0.268972327  0.430633493 -0.226490219
##  [876]  0.316080072  1.908839090  2.060872005  0.432190572  1.686372531
##  [881] -0.335921356  0.030715500 -0.132340835  1.148810271  1.648739617
##  [886]  1.213044318  1.304745556  2.063964207  1.233985648  0.071812332
##  [891]  1.761722417  1.970966572  1.945389212 -0.160489136  1.752021339
##  [896]  2.022456373  1.097899058  1.151345213  0.081165099  1.726853073
##  [901]  1.761738672 -0.478707707  1.868464138  1.739925526 -0.666327200
##  [906]  1.144523176  0.320310038  1.639190086  0.084786495  1.144633408
##  [911]  0.058405868  0.604900316  1.178046381  1.104607568 -0.286833698
##  [916] -0.155870441  1.172311813  0.099610689 -0.721808160  0.384802432
##  [921]  1.972035871  1.055607387  0.897278694  0.374099994  0.532858428
##  [926]  0.686748369  1.771807581  1.092244439  1.229209898  0.463940102
##  [931] -0.596462599  0.628203505  0.149691561 -0.124432219  1.867472376
##  [936]  0.726206469  2.185932148  1.712699008  0.012108313 -0.334612868
##  [941] -0.270061032  1.906051098  1.237767355 -0.272090000  1.820274441
##  [946]  0.843870837  1.817326583  1.748068496 -0.224434529 -0.989046376
##  [951] -0.521697761 -0.258000882 -0.310612343  1.937575078  1.809210978
##  [956] -0.711003362  0.642712524  1.260525393 -0.039288598  1.510573627
##  [961]  1.600140243  2.473641889 -1.595520622  2.057596575  1.482214088
##  [966]  0.120100385  0.144170694  1.756297130  1.046774699  1.299155878
##  [971] -1.388456720 -0.628534712  1.866144306  2.017666328 -1.541276042
##  [976]  1.918609868 -1.320596496  1.597626152  1.752138179  2.191335629
##  [981]  1.319561807  1.058126663  2.142205616  1.094490902  1.807688380
##  [986]  0.398584773  1.630400506  0.075026481  1.778265617 -0.369078236
##  [991] -0.347641756  0.360438032 -0.313566670  0.461590505  1.782958224
##  [996]  1.064651232  1.626071704  1.907467003 -0.179436721  1.113788747
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
##   0.65724384   0.44231678 
##  (0.13987285) (0.09890593)
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
## [1] -0.9381802 -0.2158264 -0.2559062 -0.4677836 -0.3545510 -0.3621433
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
## [1] 0.0167
```

```r
sd(results2$thetastar)  #the standard deviation of the theta stars is the SE of the statistic (in this case, the mean)
```

```
## [1] 0.9113955
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
## t1*      4.5 -0.04144144   0.8932077
```

One of the main advantages to the 'boot' package over the 'bootstrap' package is the nicer formatting of the output.

Going back to our original code, lets see how we could reproduce all of these numbers:


```r
table(sample(x,size=length(x),replace=T))
```

```
## 
## 1 2 5 8 9 
## 1 1 3 3 2
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
## [1] 0.0174
```

```r
se.boot
```

```
## [1] 0.8988688
```

Why do our numbers not agree exactly with those of the boot package? This is because our estimates of bias and standard error are just estimates, and they carry with them their own uncertainties. That is one of the reasons we might bother doing jackknife-after-bootstrap.

The 'boot' package has a LOT of functionality. If we have time, we will come back to some of these more complex functions later in the semester as we cover topics like regression and glm.

