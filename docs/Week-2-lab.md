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
## 0 2 3 5 7 8 9 
## 1 1 1 1 2 2 2
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
## [1] -0.0534
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
## [1] 2.698148
```

```r
UL.boot
```

```
## [1] 6.195052
```

Method #2: Simply take the quantiles of the bootstrap statistics


```r
quantile(xmeans,c(0.025,0.975))
```

```
##  2.5% 97.5% 
##   2.7   6.2
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
##    [1] 4.5 4.3 4.4 6.4 4.2 4.4 4.0 3.1 4.8 5.0 4.3 4.3 6.6 4.4 6.2 3.2 2.6 5.1
##   [19] 5.8 5.4 5.0 5.3 3.7 5.2 4.8 3.5 4.2 4.5 5.2 4.3 4.3 4.4 4.7 4.8 4.6 4.8
##   [37] 4.5 3.3 4.0 3.4 4.2 5.3 2.5 5.1 3.9 4.9 4.9 5.2 4.0 5.6 5.1 5.0 3.4 4.3
##   [55] 3.2 4.5 7.0 4.5 4.3 4.2 3.6 4.4 3.3 4.8 3.3 4.1 4.9 2.9 4.6 4.2 4.7 3.5
##   [73] 3.3 3.6 6.0 3.7 5.8 4.2 4.9 5.5 6.2 3.4 4.8 3.2 4.2 3.9 2.8 3.4 4.5 6.6
##   [91] 3.4 4.3 5.6 4.8 3.6 5.2 5.7 4.0 4.0 4.6 4.2 5.9 4.4 2.8 5.0 4.9 3.7 4.9
##  [109] 4.1 6.3 4.9 4.6 5.1 5.0 4.2 5.3 1.8 4.7 5.6 4.1 5.2 3.8 5.6 5.1 5.3 3.5
##  [127] 4.4 5.1 5.1 4.7 5.1 3.8 3.4 3.2 3.9 4.6 4.3 3.6 4.8 3.8 3.5 5.3 5.2 3.0
##  [145] 3.9 4.6 5.2 4.8 4.4 3.9 4.6 2.9 3.6 3.2 4.8 3.6 2.5 4.7 3.9 4.4 5.6 3.6
##  [163] 4.3 4.1 4.8 4.3 2.6 5.4 3.1 5.1 3.4 6.2 4.9 4.8 5.5 5.2 4.2 3.9 4.6 5.3
##  [181] 3.8 3.8 3.6 4.6 3.4 3.5 5.0 5.6 3.4 5.2 4.9 3.8 3.6 4.6 4.7 3.5 5.5 2.4
##  [199] 5.8 5.7 3.9 5.5 5.0 2.8 2.9 3.8 4.4 4.1 3.2 3.4 6.2 2.9 5.1 4.8 4.6 5.4
##  [217] 5.0 3.9 4.0 4.4 5.8 3.8 5.0 5.2 4.0 3.4 5.7 3.2 5.7 5.0 4.6 3.5 5.0 3.3
##  [235] 4.6 4.6 3.4 5.4 4.3 4.6 6.2 4.0 3.2 4.8 5.9 7.0 5.4 4.3 5.2 5.6 4.2 4.2
##  [253] 4.4 4.1 6.5 4.4 5.4 4.9 5.8 4.3 6.1 5.6 5.1 2.5 4.7 4.8 3.8 5.4 5.4 4.7
##  [271] 5.4 2.5 4.2 5.1 4.8 4.3 5.2 4.1 5.1 4.8 2.5 3.5 3.7 3.6 2.9 3.8 4.1 5.9
##  [289] 3.4 2.9 5.6 3.8 5.7 2.9 3.9 3.6 4.3 4.6 4.2 5.3 5.1 4.6 5.2 3.5 4.2 5.0
##  [307] 4.6 2.5 5.2 3.8 4.8 4.7 5.2 4.3 6.4 3.6 3.5 5.4 4.9 4.7 5.4 3.7 2.9 4.6
##  [325] 3.7 5.0 3.6 2.7 5.5 4.1 5.0 5.1 5.5 3.7 4.8 3.6 3.6 5.4 4.6 5.3 3.4 4.5
##  [343] 5.1 5.1 4.2 5.3 4.2 2.1 3.8 2.1 6.0 5.1 3.6 4.6 4.6 6.5 5.0 3.4 4.4 2.6
##  [361] 4.1 2.4 6.1 4.6 4.5 4.6 4.4 5.8 5.3 3.2 4.0 4.1 5.0 5.3 4.7 5.0 4.1 5.4
##  [379] 4.7 3.4 3.8 4.4 5.1 3.3 4.3 3.7 3.6 5.2 3.7 4.1 3.8 5.0 3.9 4.5 4.1 5.2
##  [397] 5.4 4.1 4.7 4.8 5.1 4.6 4.5 6.0 4.3 4.8 4.5 4.5 5.1 5.1 3.5 5.3 5.2 5.5
##  [415] 5.4 3.3 5.2 2.9 2.8 3.8 3.4 3.5 4.4 4.4 3.6 5.2 4.0 5.1 5.2 5.0 4.0 4.4
##  [433] 5.3 5.9 5.1 4.9 4.9 4.5 5.0 4.4 3.8 2.1 5.0 4.7 4.3 5.5 5.4 4.7 2.9 3.2
##  [451] 4.1 4.4 6.4 4.0 4.6 4.4 3.9 3.2 3.6 4.0 5.6 4.4 4.5 5.2 3.0 5.0 5.1 5.2
##  [469] 2.9 3.8 3.5 4.7 5.0 4.6 4.6 4.9 4.5 4.5 4.7 4.9 6.4 4.4 5.1 4.8 4.5 4.2
##  [487] 5.8 5.0 3.5 5.4 4.2 4.4 4.0 4.7 5.6 5.0 3.9 6.5 2.9 3.3 3.7 3.4 6.2 5.5
##  [505] 3.3 3.9 5.4 3.7 4.2 4.4 5.9 4.5 6.8 3.8 6.0 3.8 4.5 6.5 4.9 2.2 3.1 5.7
##  [523] 4.5 4.2 3.3 5.4 4.5 4.0 4.5 4.9 5.0 4.7 4.0 5.0 3.5 4.4 4.1 5.3 5.1 4.8
##  [541] 5.4 3.9 6.0 2.9 4.9 4.8 3.8 3.4 3.3 4.8 5.2 3.3 5.5 3.1 4.5 5.1 4.2 5.1
##  [559] 2.4 4.8 4.3 3.5 1.9 4.7 4.0 6.1 5.9 2.8 4.3 3.0 4.1 4.0 4.9 4.6 5.4 5.3
##  [577] 6.6 2.5 4.6 3.6 6.6 3.2 3.6 3.5 5.7 3.5 4.8 4.6 5.1 3.7 5.2 3.4 5.5 4.5
##  [595] 5.3 4.9 3.8 3.6 4.2 3.7 3.6 4.9 5.3 4.1 4.5 3.7 5.4 5.2 5.0 3.2 5.0 4.9
##  [613] 5.4 5.3 4.8 3.6 4.0 4.0 5.0 5.3 4.9 5.2 5.5 4.2 4.0 4.1 3.6 5.5 3.3 4.5
##  [631] 5.4 3.7 5.0 4.2 6.1 4.9 4.8 3.1 5.2 5.5 4.5 5.6 4.3 4.0 5.0 5.0 5.9 4.2
##  [649] 2.9 3.9 6.4 3.8 4.9 4.1 5.0 3.0 5.3 4.4 2.7 3.4 3.9 3.1 2.1 4.0 5.7 5.1
##  [667] 3.8 3.7 3.7 3.3 3.4 4.9 4.3 4.4 3.5 4.3 3.9 6.3 5.3 3.0 4.2 2.7 4.2 4.1
##  [685] 3.4 4.7 3.5 4.4 4.0 4.0 3.4 6.3 5.3 4.0 6.7 4.7 4.5 3.4 4.9 4.6 4.0 5.4
##  [703] 4.2 5.4 4.9 4.3 5.0 3.9 2.9 4.0 4.3 4.5 3.0 3.2 3.6 3.2 4.6 4.9 3.4 4.4
##  [721] 5.4 3.6 5.4 3.2 3.5 4.8 3.7 5.1 4.4 6.5 4.2 4.6 4.1 6.5 4.1 5.6 3.5 4.9
##  [739] 5.2 5.4 6.3 6.7 6.4 6.4 3.6 6.2 3.1 3.3 3.9 4.7 4.4 4.8 4.2 4.9 3.4 6.5
##  [757] 3.6 4.3 4.1 3.9 5.1 4.2 3.5 4.6 2.1 4.2 4.0 4.7 4.7 5.9 4.5 5.3 6.1 4.3
##  [775] 6.0 4.7 6.3 4.3 3.4 3.4 3.4 2.6 5.0 4.1 4.4 4.6 3.6 3.2 4.4 4.8 4.5 4.3
##  [793] 5.3 5.0 4.3 3.6 3.8 6.9 5.1 4.5 2.7 2.3 3.7 4.0 5.7 4.4 4.4 4.1 3.7 4.4
##  [811] 5.4 3.0 4.8 6.0 6.4 4.9 2.3 3.6 3.9 5.8 4.2 4.2 2.3 4.8 4.2 4.7 5.7 2.8
##  [829] 5.2 4.5 4.7 5.8 4.7 5.4 5.1 3.8 3.6 5.7 4.5 4.6 4.2 5.4 4.6 4.5 4.7 3.3
##  [847] 5.1 4.6 3.6 3.1 5.4 5.7 5.1 4.7 3.6 4.3 4.5 2.9 6.5 4.4 4.6 5.1 5.3 4.7
##  [865] 4.4 4.5 4.1 4.2 5.9 4.9 5.8 5.9 3.7 4.9 5.1 4.8 4.2 4.8 6.8 3.3 2.4 5.6
##  [883] 2.4 3.3 5.3 4.1 5.5 4.2 4.3 5.1 3.3 5.3 4.6 3.8 4.9 4.6 6.4 3.3 3.0 5.1
##  [901] 5.4 5.2 5.0 6.2 3.6 4.5 4.0 3.2 4.9 4.9 3.7 3.3 4.0 3.3 4.2 4.5 7.9 4.3
##  [919] 5.4 3.8 4.7 5.1 4.7 3.9 4.1 5.2 5.9 4.0 5.1 4.9 4.3 4.1 3.5 4.4 4.2 4.6
##  [937] 3.7 5.0 4.1 6.4 4.4 3.8 3.9 4.9 2.8 5.5 5.9 5.3 4.8 3.0 3.9 3.7 4.9 3.9
##  [955] 4.2 4.5 4.5 5.9 5.6 5.6 4.6 6.9 5.8 5.5 5.8 4.6 5.1 3.9 5.8 4.0 4.8 5.7
##  [973] 4.2 4.9 4.0 5.6 4.9 4.7 1.2 3.9 3.5 5.6 4.3 3.6 4.6 4.1 5.0 5.4 4.3 5.5
##  [991] 4.2 3.9 4.0 3.5 3.7 3.4 3.8 3.5 6.3 5.8
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
##   2.6   6.4
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
##    [1] 4.1 4.9 5.5 3.5 3.7 4.7 3.4 5.6 3.2 5.3 3.9 4.1 5.4 3.8 4.7 3.5 4.9 5.1
##   [19] 5.5 2.8 2.5 5.4 3.5 3.5 4.9 4.8 5.0 3.7 3.5 4.8 4.2 5.9 5.8 4.3 6.0 3.6
##   [37] 3.5 5.5 3.7 5.8 2.3 5.2 6.1 4.4 5.5 4.1 4.9 3.8 4.0 5.8 5.4 3.4 4.9 6.0
##   [55] 5.8 4.3 4.7 5.4 4.5 3.4 5.5 3.6 2.7 6.2 4.3 5.4 4.5 4.9 4.6 5.4 4.7 4.9
##   [73] 4.3 3.8 5.4 4.4 3.8 4.8 5.1 4.5 4.9 4.1 3.3 5.3 4.9 5.5 3.0 4.1 4.7 3.9
##   [91] 4.4 4.8 4.7 5.2 5.7 4.6 4.8 5.8 3.0 3.9 4.1 4.4 4.1 5.7 4.3 4.1 4.0 3.3
##  [109] 5.0 5.3 4.5 5.2 4.3 3.4 4.5 5.3 4.5 4.6 5.7 3.2 4.1 3.7 4.4 4.0 6.2 3.3
##  [127] 2.5 3.2 6.5 4.9 5.3 2.5 3.4 2.9 4.9 4.7 4.0 5.0 4.6 4.2 5.6 4.8 5.2 3.8
##  [145] 5.5 3.4 4.9 4.2 5.0 4.8 1.7 3.3 4.6 4.4 2.9 3.6 4.5 4.6 5.9 4.7 1.0 5.2
##  [163] 4.5 4.1 3.5 3.8 3.8 5.5 4.6 4.9 3.5 4.2 5.0 4.6 5.6 4.4 3.5 4.7 4.7 4.7
##  [181] 3.7 3.1 4.4 4.0 4.1 4.6 5.1 5.5 3.7 2.9 5.9 4.0 4.4 3.6 4.2 4.9 5.3 4.7
##  [199] 4.8 4.6 3.8 3.3 5.1 3.3 5.1 2.0 3.8 5.1 4.4 4.1 4.9 5.7 3.3 3.8 5.1 4.7
##  [217] 5.2 5.5 5.6 5.3 3.7 4.7 2.2 5.0 3.3 6.6 3.7 4.0 4.2 4.8 6.3 5.5 5.5 3.8
##  [235] 5.7 4.7 4.2 5.0 3.8 4.5 3.3 5.5 4.5 5.0 4.4 5.3 5.9 4.4 4.3 3.1 5.1 3.2
##  [253] 4.2 3.7 4.8 4.3 4.6 3.8 3.2 5.7 3.4 7.0 6.3 4.4 4.8 4.8 4.5 5.2 3.5 5.0
##  [271] 3.2 4.1 4.9 3.7 3.6 4.1 3.9 5.1 6.6 3.6 4.3 4.5 4.2 5.8 5.8 4.4 5.8 3.5
##  [289] 5.9 5.5 3.8 6.0 5.3 5.2 5.2 3.0 5.2 5.0 4.4 4.0 3.9 6.6 3.3 3.7 6.1 5.0
##  [307] 5.8 6.1 5.0 4.8 3.3 3.9 4.3 3.0 6.7 3.7 4.7 3.0 4.8 3.8 4.3 5.8 4.2 3.7
##  [325] 5.0 5.1 5.2 3.7 4.3 3.6 5.5 5.1 4.6 6.3 5.9 4.9 4.8 4.4 4.6 3.8 4.7 4.7
##  [343] 4.7 3.6 3.6 5.5 2.4 4.2 4.2 4.0 3.3 3.7 5.3 5.2 3.7 4.8 4.9 4.7 4.6 5.2
##  [361] 5.5 4.2 5.5 5.4 3.8 5.3 3.6 4.8 2.9 5.3 4.6 4.8 4.0 4.1 5.5 4.7 2.8 4.1
##  [379] 4.3 5.4 3.4 4.6 5.6 5.0 3.6 3.7 4.0 3.7 3.9 5.3 5.7 5.3 5.2 3.2 4.2 5.0
##  [397] 3.8 3.9 5.0 3.9 5.2 3.7 4.7 2.7 2.8 6.0 4.4 6.3 4.6 4.3 5.7 3.1 5.1 2.8
##  [415] 3.4 6.8 4.0 3.9 3.6 4.6 3.5 5.0 4.0 4.7 4.9 4.2 5.2 4.4 4.2 3.9 4.0 2.7
##  [433] 4.3 4.2 4.2 3.5 5.1 4.5 3.1 5.4 5.0 4.2 4.2 5.8 5.0 4.9 6.0 4.8 3.2 4.6
##  [451] 5.0 6.6 5.0 5.1 5.1 3.3 3.3 5.9 3.6 3.0 4.0 4.7 5.2 5.9 4.5 5.5 3.3 5.7
##  [469] 5.3 3.3 4.8 4.6 3.4 4.7 2.8 4.6 4.0 3.1 6.3 3.8 4.3 3.0 4.1 5.4 3.8 4.6
##  [487] 3.7 3.9 5.0 4.1 3.8 4.8 3.9 5.4 4.9 4.7 2.9 3.2 5.0 6.5 5.1 4.5 4.8 5.9
##  [505] 3.8 5.2 4.3 4.1 3.8 4.8 3.9 4.8 5.1 6.5 6.7 5.1 4.1 5.3 2.8 5.3 5.3 4.5
##  [523] 2.9 5.0 4.8 4.4 4.9 3.4 4.0 4.4 4.3 3.5 3.9 3.5 5.9 4.5 5.1 2.4 4.3 5.3
##  [541] 2.8 4.8 2.0 5.3 5.1 4.4 3.6 4.3 4.4 4.2 3.3 5.8 3.0 5.7 3.4 4.0 4.0 4.1
##  [559] 5.0 4.3 3.8 4.3 5.1 4.9 4.0 5.7 3.1 5.0 4.4 5.1 3.4 3.6 3.2 5.6 3.9 4.4
##  [577] 3.7 3.7 3.7 4.1 5.5 3.4 3.4 4.8 4.4 5.1 4.1 6.1 3.5 5.2 4.8 5.2 5.3 5.0
##  [595] 4.6 4.8 5.0 4.6 4.4 5.0 4.2 4.1 4.5 4.9 4.6 4.7 6.0 4.7 5.0 5.1 4.0 5.2
##  [613] 5.5 4.9 3.4 3.4 5.3 4.0 4.5 6.7 5.5 4.2 3.3 4.7 3.5 4.2 4.1 5.3 4.1 4.9
##  [631] 5.0 6.6 4.1 5.2 5.0 4.4 5.7 4.6 2.1 3.2 4.2 4.2 3.6 4.9 2.9 4.4 4.8 5.1
##  [649] 4.5 4.7 4.2 4.7 3.8 5.5 6.4 4.3 2.8 5.5 3.2 2.9 3.0 5.3 4.4 4.7 4.7 3.8
##  [667] 3.5 4.8 4.0 6.1 3.1 3.2 5.4 3.8 3.5 6.8 6.8 3.9 3.2 5.2 4.3 3.9 3.5 3.4
##  [685] 6.0 4.4 5.7 6.6 5.4 3.7 4.2 3.6 4.2 3.4 5.8 5.0 3.6 5.6 5.1 6.0 5.1 5.0
##  [703] 4.6 5.4 3.7 4.6 4.9 4.9 4.7 5.1 5.1 3.4 2.0 3.6 3.3 3.8 3.6 4.5 5.6 4.4
##  [721] 3.9 3.9 4.7 5.2 4.2 3.2 5.9 4.4 4.0 5.7 4.6 5.4 4.1 4.2 6.1 6.2 3.7 5.0
##  [739] 4.2 5.5 5.8 5.4 2.9 3.2 5.3 3.6 3.6 4.7 4.3 4.3 4.5 3.1 4.3 3.2 4.5 3.7
##  [757] 3.7 5.1 3.9 4.9 5.5 5.1 3.7 5.6 3.2 6.0 5.0 5.5 5.4 5.6 5.3 4.0 7.2 4.9
##  [775] 3.2 3.7 5.9 4.5 5.0 3.6 4.7 4.7 6.4 4.9 6.2 3.5 6.0 4.5 5.4 5.6 4.2 4.3
##  [793] 3.9 5.8 3.0 5.4 5.8 2.8 4.4 5.1 6.6 5.0 5.4 4.2 4.9 4.1 5.1 5.4 5.2 5.6
##  [811] 2.0 3.4 4.6 5.9 4.8 5.9 4.5 4.5 4.9 5.1 5.2 4.7 3.7 3.9 4.4 4.8 3.3 4.4
##  [829] 5.8 3.8 5.8 3.8 5.1 5.0 3.4 3.8 5.4 1.9 4.7 5.0 3.9 4.0 3.0 4.0 5.3 3.6
##  [847] 5.2 4.6 5.7 4.1 4.5 3.9 5.8 4.3 4.0 5.8 5.0 4.0 4.7 3.3 4.5 3.4 3.9 5.2
##  [865] 6.2 6.3 3.0 4.0 4.1 5.0 5.3 5.5 6.2 4.0 5.1 3.6 3.0 3.4 5.9 4.6 4.9 4.9
##  [883] 3.6 4.0 4.8 5.1 4.1 5.3 4.3 4.7 5.6 4.2 3.9 5.4 4.6 2.8 5.5 4.9 5.5 2.5
##  [901] 5.6 3.8 4.4 2.8 2.4 4.1 5.5 6.0 5.1 4.6 4.8 4.5 4.3 3.5 4.6 3.8 5.6 4.3
##  [919] 3.2 4.1 4.1 4.5 5.2 4.4 4.3 4.2 4.0 4.4 5.4 3.2 5.1 5.0 4.3 5.9 4.1 2.5
##  [937] 4.6 4.2 5.3 6.8 4.1 3.1 5.4 5.3 5.7 5.1 3.9 5.1 5.7 5.2 4.2 4.8 4.5 5.0
##  [955] 2.4 4.0 3.8 5.3 6.6 2.7 5.8 5.2 4.9 2.5 5.0 4.9 4.5 4.6 4.1 4.9 3.5 6.4
##  [973] 2.6 3.7 4.9 4.0 5.8 4.5 4.3 4.1 3.4 3.9 6.2 5.7 4.6 5.0 4.9 3.9 4.0 4.5
##  [991] 4.9 4.5 4.5 5.2 6.0 4.6 3.5 4.8 4.8 5.2
## 
## $func.thetastar
## [1] 0.0104
## 
## $jack.boot.val
##  [1]  0.54186747  0.38879056  0.30921788  0.19387755  0.11749271 -0.01137026
##  [7] -0.13901734 -0.43873874 -0.33994778 -0.47758621
## 
## $jack.boot.se
## [1] 1.015267
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
##    [1] 5.3 3.9 3.8 4.1 5.6 4.9 4.9 4.9 4.8 3.2 3.5 5.3 6.9 5.7 5.5 4.7 4.8 5.4
##   [19] 3.9 3.8 4.2 5.0 3.8 4.1 4.2 4.9 3.8 5.0 3.9 4.4 5.4 3.9 4.1 4.8 4.1 4.6
##   [37] 3.1 4.3 4.0 5.4 4.9 4.1 4.3 4.5 3.3 2.4 5.1 2.5 4.1 4.4 4.6 5.1 5.9 5.6
##   [55] 4.8 3.4 7.6 4.7 2.8 4.8 3.5 2.3 4.9 3.4 5.7 5.4 6.0 5.2 5.8 5.0 4.6 3.7
##   [73] 3.7 5.8 3.9 4.7 5.4 5.9 5.9 3.4 4.8 5.5 5.1 5.9 3.1 4.8 4.7 3.9 3.9 3.6
##   [91] 4.5 5.4 4.1 4.2 5.3 3.2 4.2 6.2 4.3 3.9 5.4 4.1 4.4 5.0 5.2 7.3 4.6 3.9
##  [109] 4.1 3.1 3.8 4.7 3.9 4.5 5.3 4.8 3.7 3.9 4.7 4.5 4.6 3.8 4.5 4.9 4.1 3.4
##  [127] 4.5 4.5 3.1 5.1 6.6 4.5 4.2 4.1 4.3 4.5 3.6 5.7 5.7 3.9 4.3 3.2 4.4 5.8
##  [145] 5.4 3.1 6.1 4.0 3.7 4.3 4.5 4.2 5.5 3.6 4.2 5.4 5.2 5.2 3.4 3.7 5.8 3.8
##  [163] 4.0 5.0 4.9 4.8 5.4 4.5 6.0 3.8 4.6 4.8 5.0 3.9 5.0 4.9 4.1 4.5 6.0 4.2
##  [181] 4.6 4.7 3.9 3.7 5.4 5.8 5.5 5.5 3.6 5.2 5.4 4.7 3.6 5.5 3.7 5.3 3.1 4.2
##  [199] 3.5 4.8 7.1 3.5 5.1 5.0 6.3 3.6 5.5 4.3 5.0 4.7 4.0 5.4 4.7 4.9 5.7 4.9
##  [217] 3.2 5.0 3.7 4.0 5.0 4.4 6.1 4.8 5.2 4.1 4.1 3.8 4.7 3.7 4.6 4.2 4.4 4.2
##  [235] 5.3 4.6 3.7 3.7 4.9 4.8 4.6 5.3 3.7 4.2 4.4 4.5 5.7 5.9 3.4 3.4 2.8 2.9
##  [253] 4.2 5.5 5.7 3.7 4.6 4.8 3.8 4.3 3.5 4.1 3.2 5.8 3.5 5.0 4.0 4.4 3.9 4.6
##  [271] 3.6 5.1 4.2 2.7 4.1 5.5 3.7 5.1 5.8 4.7 4.7 2.6 4.5 5.2 3.7 5.6 2.7 5.9
##  [289] 4.8 3.4 6.6 4.8 4.2 4.5 3.8 4.0 6.1 5.2 4.7 5.0 5.1 3.8 5.2 3.0 3.6 3.6
##  [307] 4.7 5.1 4.3 3.6 4.4 3.9 2.8 4.5 4.0 3.0 3.0 3.6 5.7 3.4 5.1 3.4 5.7 3.9
##  [325] 4.3 4.7 3.1 5.6 3.2 4.7 3.7 5.3 5.7 4.1 5.5 4.9 3.8 5.7 3.3 4.5 5.4 4.4
##  [343] 5.1 5.2 5.0 3.7 4.9 6.6 5.0 3.9 4.7 4.2 5.6 5.3 5.2 4.9 5.4 2.5 3.5 3.5
##  [361] 5.5 4.6 5.0 5.3 5.6 3.9 3.8 4.1 2.8 4.1 3.9 3.9 4.3 3.6 4.3 3.7 2.8 4.1
##  [379] 5.2 2.4 5.0 3.4 4.9 4.3 5.3 4.2 5.0 4.9 3.7 2.8 4.1 3.8 4.0 5.7 5.7 4.7
##  [397] 3.2 4.2 4.7 3.6 3.7 4.0 3.9 4.6 4.1 3.3 5.0 4.0 4.2 5.4 3.6 3.1 3.6 4.4
##  [415] 3.7 5.9 6.4 3.9 4.6 4.6 5.0 4.7 4.9 4.0 5.0 2.6 5.0 4.5 5.9 4.8 5.8 4.5
##  [433] 4.9 5.0 3.0 3.7 4.5 4.2 4.2 4.2 4.8 4.0 4.4 5.6 5.9 4.5 4.6 4.5 4.5 3.5
##  [451] 5.2 4.1 4.1 4.4 2.5 5.7 7.4 4.9 2.4 4.6 5.2 3.6 4.6 6.1 2.7 3.2 3.1 4.5
##  [469] 4.0 5.4 3.8 4.6 3.7 4.0 6.0 3.9 5.1 5.8 4.7 5.3 4.1 3.0 5.1 3.4 5.1 3.7
##  [487] 3.5 3.4 4.5 4.2 3.1 3.2 4.6 4.0 3.5 4.1 4.1 4.5 4.3 4.5 4.0 4.3 4.7 5.4
##  [505] 5.4 2.6 5.1 4.1 3.5 5.2 5.1 4.2 5.2 3.0 4.8 4.5 3.9 4.2 4.4 5.0 5.1 5.4
##  [523] 3.5 5.5 5.1 5.0 4.8 5.2 3.0 4.7 2.8 5.3 2.9 5.3 4.6 3.5 5.5 5.0 4.2 5.3
##  [541] 3.7 4.5 4.4 3.3 3.7 3.3 3.2 5.7 2.9 4.2 5.0 5.7 4.3 4.7 3.7 5.7 5.5 4.3
##  [559] 3.9 4.9 7.1 4.9 4.2 3.3 5.1 4.2 3.8 3.1 4.9 5.1 4.4 4.7 5.1 3.9 4.4 4.5
##  [577] 4.3 4.1 4.5 5.1 4.6 2.9 3.6 5.0 4.0 5.0 5.5 4.1 4.9 3.8 6.0 4.7 6.6 3.5
##  [595] 4.3 5.8 5.3 4.7 5.9 4.3 3.5 3.8 3.5 4.3 2.6 5.1 4.4 3.1 5.8 3.6 4.2 4.1
##  [613] 6.0 3.5 4.1 3.9 3.7 3.7 4.9 5.2 4.8 4.2 3.4 5.2 4.7 5.8 4.0 4.2 5.0 4.2
##  [631] 3.5 6.4 6.0 4.1 5.0 3.0 4.4 4.7 4.8 5.6 3.7 5.8 4.9 4.7 4.7 5.3 5.8 2.6
##  [649] 4.2 5.3 2.7 5.1 3.6 3.7 4.0 4.6 3.1 6.0 5.2 5.2 3.1 4.3 5.1 4.1 4.9 5.5
##  [667] 4.4 3.7 3.1 5.6 3.7 4.7 4.4 4.0 6.7 4.2 3.6 5.3 3.7 5.4 4.4 4.4 4.7 4.3
##  [685] 3.7 3.8 6.5 7.0 4.8 4.3 5.0 5.8 4.6 4.6 3.5 3.9 4.9 3.6 4.0 5.1 3.5 4.3
##  [703] 5.0 3.2 5.0 3.9 4.7 4.4 5.5 4.0 3.6 7.0 4.5 5.1 4.9 4.9 2.9 4.4 4.6 5.4
##  [721] 3.7 5.2 3.5 3.5 4.9 5.6 3.2 4.7 5.0 3.1 5.0 4.9 4.9 5.8 4.0 4.4 3.9 4.2
##  [739] 4.8 3.2 5.3 2.6 5.6 4.3 5.7 4.1 3.6 5.7 4.8 4.1 5.0 4.8 5.5 3.8 3.4 2.5
##  [757] 7.0 4.7 3.1 4.8 3.3 5.7 4.2 3.9 5.7 4.6 3.4 3.2 4.4 4.3 4.8 3.6 6.0 3.7
##  [775] 3.9 3.9 3.4 5.6 4.7 4.3 4.3 2.8 4.7 3.4 4.9 4.6 6.1 5.7 2.8 4.7 4.3 4.5
##  [793] 1.8 5.5 4.5 3.4 4.8 6.0 4.2 4.5 3.1 4.2 4.7 5.8 4.0 4.6 3.9 4.4 4.9 4.6
##  [811] 4.6 4.1 5.6 3.4 3.8 4.8 3.0 3.8 4.5 4.5 3.3 5.0 5.3 5.2 4.1 3.0 4.6 3.6
##  [829] 4.4 5.8 5.3 4.3 3.9 4.5 3.8 4.7 5.2 3.7 4.3 4.9 3.9 5.1 3.7 4.0 4.8 4.2
##  [847] 5.1 6.4 3.5 4.2 3.8 5.3 3.4 4.2 5.7 3.4 4.0 4.0 3.8 5.8 4.3 6.4 3.8 2.1
##  [865] 4.6 4.8 4.0 3.7 6.9 4.0 5.3 3.2 2.9 6.4 5.0 5.2 4.4 2.5 4.3 3.6 5.2 4.7
##  [883] 4.6 5.5 3.8 3.7 6.1 4.8 4.1 3.8 4.6 3.4 4.2 5.3 4.4 5.4 5.2 4.5 5.3 2.7
##  [901] 5.7 3.8 4.4 5.4 3.5 5.7 5.3 4.4 4.6 4.5 5.4 3.0 4.4 4.4 4.4 3.7 5.1 4.9
##  [919] 4.2 5.9 4.9 4.8 3.9 4.6 4.8 5.9 3.7 4.3 4.8 4.9 4.3 5.1 4.7 3.6 5.5 4.7
##  [937] 4.8 4.7 4.0 4.6 4.7 3.9 5.5 4.0 3.1 3.4 5.3 5.4 4.3 6.1 5.5 4.7 5.6 5.4
##  [955] 4.8 6.0 3.8 4.7 4.4 4.7 4.0 5.7 3.5 4.2 4.2 4.6 4.9 5.4 5.6 3.6 5.3 3.8
##  [973] 3.9 4.0 2.6 4.9 5.5 4.2 3.6 4.7 4.2 3.1 3.7 3.1 5.2 5.1 4.4 4.3 2.7 5.3
##  [991] 2.4 5.8 4.5 4.1 5.0 5.2 3.9 6.1 4.2 5.6
## 
## $func.thetastar
## 72% 
##   5 
## 
## $jack.boot.val
##  [1] 5.4 5.3 5.3 5.2 5.1 5.0 4.9 4.7 4.7 4.3
## 
## $jack.boot.se
## [1] 0.980867
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
## [1] -0.1576571
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
##   3.984512   5.762069 
##  (1.712563) (2.639439)
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
## [1]  0.6711694  0.7630019 -0.1654172  0.1035048  0.3574922  0.7024930
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
##    [1] -1.481671e-01 -1.641908e-01 -1.079753e+00 -8.103328e-01 -4.399774e-01
##    [6]  7.115045e-01  2.069458e-01 -2.276363e-01  1.419732e-01 -3.691873e-01
##   [11] -2.935550e-01 -5.708948e-01  3.281413e-02 -3.370228e-01 -1.100515e-01
##   [16]  6.808226e-01  7.140763e-01 -1.625997e+00 -4.203156e-01  4.232130e-01
##   [21] -2.373755e-02 -4.964126e-01  6.665069e-01 -1.220203e-01  2.153168e-01
##   [26] -1.638218e-01  1.981237e-01 -7.664671e-01  3.564732e-01 -9.620447e-01
##   [31] -6.638773e-01  7.219646e-01 -2.404561e-01 -8.721025e-02 -2.002364e-01
##   [36] -3.253210e-01 -6.948233e-01  1.894741e-01  2.112558e-01 -5.984338e-01
##   [41] -7.322577e-01 -4.942656e-01 -1.683091e-01 -2.615861e-01  7.513026e-01
##   [46] -1.623183e-01 -8.873957e-02  2.968655e-01  8.294247e-01  1.700448e-01
##   [51] -1.703979e-01 -1.285342e-01 -6.255065e-01 -2.696575e-01  2.048962e-01
##   [56]  1.744786e-01 -2.815950e-01 -7.400765e-01 -2.098543e-02 -1.047183e-01
##   [61] -8.372750e-02  1.445341e-01 -1.180934e-01 -7.879105e-02 -8.392305e-01
##   [66]  2.780891e-01 -7.489207e-01 -4.807086e-01  5.780636e-01  3.400233e-01
##   [71] -1.576571e-01 -3.540455e-01  9.714831e-01 -4.902350e-01 -8.908946e-02
##   [76] -2.540507e-01  6.011232e-01 -2.574136e-01 -5.926892e-02 -8.250068e-01
##   [81]  5.293433e-01 -5.153459e-01  8.965616e-02  2.937916e-01 -6.222821e-01
##   [86]  7.277545e-01  1.319513e-01 -3.437290e-01 -3.160773e-01  6.815748e-01
##   [91] -4.605043e-01 -3.165174e-01 -3.077494e-01 -2.373755e-02  3.859460e-01
##   [96] -2.501475e-01 -3.084043e-01  1.987495e-01 -2.895264e-01  2.766971e-01
##  [101]  6.399169e-01  2.853928e-01 -1.102942e-01 -7.831823e-01  6.850030e-01
##  [106] -3.265915e-01  7.903350e-02  2.648937e-01  3.532234e-01 -3.937567e-01
##  [111]  2.846110e-01  2.002003e-01 -2.323818e-01 -2.681647e-01 -2.173588e-01
##  [116] -1.954210e-01  7.057674e-02  1.564028e-01 -3.750500e-01  5.951940e-01
##  [121]  1.046795e-01 -2.799637e-01  1.351993e+00 -9.422423e-02 -2.618931e-02
##  [126]  4.522857e-01  4.276429e-01  1.280980e-01  1.420752e-01 -9.364118e-01
##  [131]  2.598971e-01 -5.301296e-01 -2.295826e-01 -5.028632e-01 -3.176787e-01
##  [136]  1.714529e-01 -1.246893e+00 -5.517451e-01 -2.613040e-01  5.110460e-01
##  [141] -1.485552e+00 -3.991789e-01 -3.622334e-01  6.010733e-01 -1.441749e-02
##  [146]  8.323812e-02 -6.046446e-01 -4.085549e-01 -6.470603e-01 -2.684006e-01
##  [151]  3.016597e-01  1.322350e-02  8.569700e-02 -1.234417e-01 -5.923465e-01
##  [156]  1.181541e-01 -2.566887e-01 -1.444661e-01  1.172774e-01 -2.291916e-01
##  [161] -1.823459e-01  5.962869e-01  2.002562e-01 -3.712038e-01  1.509673e-01
##  [166] -9.032191e-01 -4.496003e-01 -2.586073e-01 -4.385218e-02 -4.177884e-01
##  [171] -1.912279e-01 -9.321713e-01 -6.145898e-01  3.453328e-01  2.447272e-01
##  [176] -2.819519e-01  3.107272e-01  2.487398e-01  5.807350e-01 -4.823885e-01
##  [181]  5.449980e-01 -2.769306e-01 -6.222821e-01 -4.085572e-01  3.543251e-01
##  [186]  2.924374e-01 -4.469897e-01  9.749022e-02  1.539159e-01 -6.855356e-01
##  [191]  8.754536e-01  5.039430e-01  4.473386e-01  9.951855e-02 -4.715449e-02
##  [196] -5.229806e-01 -8.860195e-02 -5.776458e-01  9.998345e-02  8.220922e-02
##  [201] -4.402540e-01  2.963556e-01 -3.411525e-03 -1.634328e-01 -1.876791e-01
##  [206] -2.267174e-03  2.777248e-01 -3.235813e-01  5.110148e-01 -6.320295e-01
##  [211]  2.423978e-01 -4.599410e-02 -3.466451e-01  1.485174e-02 -5.836741e-02
##  [216] -5.266010e-01  8.906142e-01  2.380611e-01 -7.167064e-02  1.509673e-01
##  [221]  2.384574e-01 -1.146112e+00  3.696521e-01  6.532445e-01 -7.122893e-01
##  [226] -2.203719e-02 -3.956640e-01 -1.023800e-01 -2.487180e-01 -3.599387e-01
##  [231]  1.008928e+00 -4.010446e-01 -1.384793e-01 -3.270468e-01  1.671141e-02
##  [236] -4.178386e-02  1.307614e-01  2.770939e-02  7.247042e-01 -4.154415e-01
##  [241]  1.127770e-01 -8.713787e-02 -1.489653e+00 -3.493428e-01 -3.907382e-02
##  [246]  5.252915e-01  9.784863e-02 -5.739833e-01  2.653802e-01  4.012293e-01
##  [251] -6.990581e-01  1.510875e-01 -4.242775e-02 -1.070683e-01  2.358427e-02
##  [256] -7.220045e-02 -1.031744e-01 -5.356053e-01 -1.353229e-01 -7.938213e-01
##  [261] -4.942046e-01 -4.108055e-01 -9.277011e-01 -1.757132e-01 -4.213231e-01
##  [266]  1.405077e-01 -1.195191e+00  4.443722e-01 -9.352861e-02  5.626709e-02
##  [271]  4.839457e-01 -3.863381e-01  3.922694e-01 -2.050124e-01 -3.523421e-01
##  [276]  7.602301e-02 -7.813331e-01 -4.442372e-01  4.783078e-01 -1.118288e-01
##  [281]  1.409042e-01 -4.862528e-01 -5.348325e-01 -4.906384e-01 -3.900322e-01
##  [286] -4.256349e-01 -8.884147e-01  1.949002e-01 -2.555296e-01 -4.808283e-01
##  [291] -3.504644e-01 -2.645432e-01  3.263147e-01 -2.576426e-01  5.264315e-01
##  [296]  1.566503e-01 -6.176896e-01 -7.659190e-01  8.769313e-01 -5.365630e-01
##  [301]  2.763175e-01 -6.097573e-01  2.268423e-01  3.468967e-01 -3.427722e-01
##  [306]  3.624325e-01 -6.604649e-01 -2.081806e-01  2.937916e-01  3.210250e-01
##  [311]  5.390501e-01 -5.236756e-01 -3.074636e-01  6.979226e-01  2.910422e-01
##  [316] -1.038216e+00  5.152583e-01 -2.171749e-01  4.193833e-01 -3.973079e-01
##  [321]  5.380146e-02  5.621287e-01  7.343997e-02 -6.445592e-02 -3.143005e-01
##  [326]  4.294459e-01 -9.498927e-02  2.663988e-01  1.955965e-01 -2.915584e-01
##  [331] -7.273046e-01  6.008063e-01  4.872142e-01 -7.957431e-01 -1.390814e-01
##  [336] -2.065664e-01  8.683034e-01  4.120592e-01 -6.322271e-01  3.410703e-01
##  [341]  2.177557e-02 -4.406888e-02 -5.437935e-01  4.645924e-01  3.643499e-01
##  [346] -4.961647e-01 -3.723909e-01 -1.771416e-01 -9.575396e-02  2.213211e-01
##  [351]  1.240460e+00  3.510928e-01  3.418043e-01 -1.320247e+00  5.725384e-02
##  [356]  3.185064e-01  6.518682e-01  3.260035e-01  4.532160e-01 -1.223256e-01
##  [361]  8.388123e-02 -5.520994e-01 -4.960580e-01  3.800922e-01  3.907350e-01
##  [366] -1.178089e+00 -2.518408e-01 -1.732018e-01  6.035194e-02 -1.714006e-02
##  [371]  1.692413e-01 -5.035157e-01 -5.510505e-02  5.488873e-01 -5.025650e-01
##  [376] -4.240999e-01 -6.102155e-01  3.163391e-01  1.916471e-01 -5.277690e-02
##  [381]  6.862199e-01 -9.129074e-01 -4.342708e-01 -4.495970e-01 -6.406503e-02
##  [386] -1.117910e+00  1.638136e-02  2.459735e-01  2.763175e-01 -2.967190e-01
##  [391] -5.556373e-01 -1.464958e-01  2.398958e-01 -9.778281e-01 -5.005421e-01
##  [396]  5.568909e-02 -3.528967e-01 -2.986021e-01 -3.078500e-01  2.206239e-02
##  [401]  2.529727e-01 -1.094334e+00 -2.464259e-01  4.158774e-01 -5.083175e-01
##  [406]  1.711285e-01  4.012746e-01 -4.611104e-01 -9.629635e-01 -3.971627e-01
##  [411]  4.726291e-01 -2.449706e-02 -3.120088e-01 -1.019217e+00 -9.956126e-01
##  [416] -5.928574e-01 -3.486108e-02 -1.997024e-01 -1.701513e-01 -1.455255e+00
##  [421] -6.039099e-02 -3.869573e-01  2.439510e-03 -4.343284e-01 -1.851796e-01
##  [426]  7.878837e-02 -2.783231e-01 -3.612502e-01  7.054500e-01 -5.069238e-01
##  [431] -1.152145e+00  5.932229e-02  2.338375e-01  9.014228e-03  2.165077e-01
##  [436]  5.068594e-01 -8.520423e-01  9.598496e-02 -9.258327e-01  1.582965e-01
##  [441] -2.036266e-02 -5.494922e-01  5.461283e-01  1.548599e-01  3.920071e-01
##  [446]  4.352473e-01  3.694384e-01 -5.081828e-01 -1.235317e-01 -6.061765e-01
##  [451] -1.001073e-01 -5.407796e-01 -1.609299e-01 -1.314717e+00  2.737251e-02
##  [456] -6.237965e-02 -2.461930e-01 -5.242575e-02  7.472758e-01  1.267390e-01
##  [461]  2.710962e-01  1.702640e-01  3.502227e-01 -7.993985e-01 -1.000807e+00
##  [466]  7.755389e-01 -3.080626e-01 -1.172860e-01 -1.331286e+00 -2.510701e-01
##  [471] -9.392699e-01  8.971433e-01 -5.437935e-01  3.770421e-01 -3.609931e-01
##  [476]  4.344240e-01  1.042510e+00 -3.823840e-01 -6.604405e-01  3.329402e-02
##  [481] -4.451696e-01 -5.006659e-01  2.632305e-01  7.235193e-03 -9.032191e-01
##  [486] -9.500250e-02 -6.922399e-02 -5.823852e-01 -6.058354e-01 -4.524996e-01
##  [491]  3.212731e-01 -5.573416e-01 -4.592916e-01 -5.586690e-01  2.386637e-01
##  [496]  7.638418e-01 -1.371905e-01 -1.021017e-03  4.385487e-01 -2.753667e-01
##  [501]  7.726429e-02  8.354071e-01  1.722818e-01  1.423476e-01 -5.291835e-01
##  [506] -4.984545e-01 -1.459252e-01 -3.979028e-01  1.331093e+00 -3.760679e-02
##  [511]  3.285106e-02 -3.604008e-01  5.740026e-01 -8.679611e-01 -2.008613e-01
##  [516] -4.621863e-01 -8.924495e-01  3.565193e-01  9.014228e-03 -3.088932e-02
##  [521]  1.075301e-01 -1.102934e+00 -8.330506e-01 -3.721761e-01  6.867162e-01
##  [526]  5.464882e-02 -6.654586e-01  1.707152e-01 -3.455463e-01  6.292066e-01
##  [531]  7.639067e-01 -4.608343e-02  3.798931e-01  2.783417e-01  2.610619e-02
##  [536] -3.375238e-02 -4.774759e-01  1.884676e-01 -1.354336e-01  5.804395e-01
##  [541] -2.373755e-02  3.216618e-03 -2.680915e-01 -1.620129e-01  5.743347e-01
##  [546]  5.121099e-01 -3.533460e-01 -6.293423e-02 -6.511759e-01  1.253833e+00
##  [551] -6.430576e-01 -1.973602e-01  1.186725e-01  2.304014e-01 -3.730290e-01
##  [556] -2.660601e-01  1.187537e+00 -1.710310e-01 -5.049554e-01  8.934116e-01
##  [561]  1.058566e-01  1.237571e-02  2.466649e-01  1.059987e-01 -3.172380e-01
##  [566]  2.571371e-01 -2.879924e-01 -3.360303e-01  2.997846e-01 -9.443452e-01
##  [571]  4.629437e-01 -8.089185e-01 -3.135796e-01 -2.398033e-01 -1.359493e-01
##  [576]  8.451121e-02 -3.616370e-01 -1.356264e-01 -4.055720e-01 -8.302982e-01
##  [581]  1.007059e-01  7.103716e-01 -9.497462e-01  2.275937e-01  6.420734e-02
##  [586] -7.322267e-02  3.445290e-01 -3.956048e-01 -3.620013e-02 -2.753172e-01
##  [591] -7.344015e-01  2.732391e-01 -3.441684e-01 -1.383721e-01  3.936217e-01
##  [596] -5.745595e-01 -3.676579e-01 -1.801165e-02 -1.004090e+00  3.923773e-01
##  [601]  5.715630e-01 -8.745696e-01  4.158171e-01 -6.498437e-01  2.808695e-01
##  [606] -2.293410e-01 -5.179892e-01  2.663988e-01 -4.969280e-01 -4.239339e-01
##  [611] -1.196760e-01  2.654360e-01 -1.116197e-01 -2.373992e-01 -2.182070e-01
##  [616] -6.162814e-01  3.122226e-02 -4.528433e-03  9.423968e-01  2.680739e-01
##  [621] -6.709341e-01 -7.806732e-01 -3.277526e-01 -6.496896e-01  6.540561e-01
##  [626] -5.682307e-01 -4.671738e-01  1.195524e-01 -7.164154e-01 -2.966010e-01
##  [631]  2.301351e-01  4.447301e-01 -1.012965e+00 -5.578013e-01  8.919976e-02
##  [636] -4.154775e-02 -2.535251e-01  8.170933e-01  6.731800e-01 -6.410893e-01
##  [641] -7.932637e-01  4.378349e-01 -3.836547e-01 -6.157489e-02  4.610860e-01
##  [646]  2.905419e-01 -4.404814e-01 -9.330414e-01 -6.405468e-01  2.485494e-01
##  [651] -3.619771e-01  3.935127e-01 -5.100589e-01 -5.632154e-01 -3.584803e-01
##  [656] -7.768683e-02 -2.556086e-01  2.736185e-02 -2.838991e-01  5.768116e-01
##  [661]  1.737668e-01 -1.003621e-01 -6.782214e-02 -5.437935e-01 -3.748336e-01
##  [666] -5.207672e-01  1.738643e-01  5.558709e-02  3.938885e-01 -1.396593e+00
##  [671] -3.106409e-01  1.268682e-01 -2.336272e-01 -5.777997e-01  2.548756e-01
##  [676] -1.473288e+00 -3.152297e-02 -3.468371e-01  4.742944e-01 -1.986049e-01
##  [681] -2.992958e-01  1.045852e-01  1.538537e-01  2.991395e-01 -2.870091e-01
##  [686] -1.501853e-01 -6.575971e-01  3.302966e-01  5.441273e-02 -4.259736e-01
##  [691] -3.756968e-01 -1.940516e-01  6.608807e-02 -2.563489e-01  7.394303e-01
##  [696] -1.984379e-01 -3.538986e-01 -1.132094e-01  4.066393e-01 -2.206323e-01
##  [701]  6.918659e-01  2.563661e-01 -1.173432e+00 -3.294024e-01  1.234232e-01
##  [706]  3.696521e-01  2.305427e-01 -4.996020e-02  2.812888e-02  1.929463e-01
##  [711]  3.672332e-01 -9.660327e-01  1.324961e-01  2.539695e-01 -2.010016e-01
##  [716] -1.680892e-01  5.541960e-01 -7.942263e-01 -9.164587e-01 -1.384843e-01
##  [721]  1.818426e-01  3.094003e-01 -4.556611e-01 -1.779892e-01  3.835447e-01
##  [726] -7.387047e-01  1.260889e-01 -4.214406e-01  3.450771e-01 -2.463806e-01
##  [731] -5.979524e-01  8.225993e-03 -2.108700e-01  2.875754e-01  6.103563e-02
##  [736] -1.130456e+00  3.027255e-02 -6.948481e-01  7.052829e-01 -5.848068e-01
##  [741] -1.518906e-01 -1.359517e-01 -1.178642e-02 -8.238676e-01 -4.852835e-02
##  [746] -1.279977e+00 -3.778752e-01 -7.128792e-01 -2.784479e-01  3.674502e-02
##  [751]  1.727353e-01  6.648662e-01 -4.305631e-02  9.052302e-01  5.828031e-01
##  [756] -1.287316e+00  2.831845e-01  8.604889e-01  3.649851e-01  1.517908e+00
##  [761] -1.215638e-01  5.669360e-01  1.302167e-01 -2.958970e-01 -3.093459e-01
##  [766] -1.215150e+00  3.978384e-01 -6.315372e-02  3.144164e-01 -4.515684e-01
##  [771]  6.035194e-02 -1.652372e-01 -4.714920e-01  2.783690e-01 -5.993758e-01
##  [776] -6.361796e-01 -9.128090e-01  4.217030e-01 -2.381134e-01  2.886776e-01
##  [781] -5.691035e-01  1.823342e-01 -7.735038e-01  6.829672e-02 -4.417371e-02
##  [786] -5.002370e-01  7.944845e-01 -3.044445e-01  8.477732e-01  4.912417e-01
##  [791]  1.719828e-02 -1.289834e-01 -3.492006e-01  3.448209e-01  1.776092e-01
##  [796] -4.209334e-01  2.381914e-01  5.081844e-01 -6.169928e-01 -2.379938e-01
##  [801] -2.407259e-01 -1.068405e-01  1.134806e-05 -3.246108e-01 -1.649831e-01
##  [806] -4.178067e-01  8.015319e-02 -2.599418e-01  3.770342e-01 -6.039099e-02
##  [811]  3.631962e-01  4.104052e-01 -3.118830e-01 -7.295920e-01  2.977894e-01
##  [816]  8.026838e-01  3.162686e-01  5.294924e-01 -4.628767e-02 -7.352809e-01
##  [821]  1.604867e-01  5.195160e-01 -1.515122e-01  2.055166e-01 -1.241747e-02
##  [826]  6.010631e-02  4.731953e-02 -1.842271e-04 -1.416565e+00 -3.687124e-02
##  [831] -4.414772e-02  8.451121e-02 -3.883250e-01  5.819893e-02 -7.017436e-01
##  [836]  2.019206e-01 -1.125524e-01 -1.926514e-01 -1.669776e-01  1.326842e+00
##  [841]  5.438868e-03 -7.906603e-01 -6.274431e-01  9.246534e-01 -9.272689e-02
##  [846]  2.318031e-01 -2.058150e-01  4.126876e-02 -4.514004e-01 -2.091854e-01
##  [851]  1.947824e-01 -2.488894e-01 -8.095076e-01 -1.074716e-01 -5.592731e-01
##  [856]  4.590108e-01 -3.727298e-01 -7.253627e-01 -3.758151e-01  8.283479e-01
##  [861] -1.130864e-01 -4.274187e-01  1.247790e-01 -4.602185e-01  8.370494e-01
##  [866] -2.340965e-01  2.877656e-02 -4.062928e-01 -1.691308e-01 -4.151207e-01
##  [871] -1.158454e+00  6.196860e-01  4.346906e-01  9.000693e-01 -7.488245e-01
##  [876] -3.579197e-01  4.647028e-03 -2.031874e+00 -3.647165e-01  1.465616e-01
##  [881]  4.293171e-01  1.876571e-02 -1.465952e-01  4.676402e-01 -8.320777e-01
##  [886]  1.009630e+00 -1.666431e-01  3.181711e-01  1.392054e+00 -3.599755e-01
##  [891] -1.003621e-01 -1.562796e-01  3.823159e-01  2.151822e-02 -8.820406e-01
##  [896] -5.872855e-01  4.012384e-01  7.570869e-01 -1.560365e+00 -1.364814e-01
##  [901] -7.404791e-01 -3.162907e-02  2.039909e-01  2.701355e-01 -1.349824e+00
##  [906] -5.709007e-01 -2.617726e-01 -5.121278e-01 -8.862225e-01  1.745621e-02
##  [911]  2.557415e-01 -2.523770e-01 -7.262455e-02 -1.199366e-01  6.993126e-01
##  [916] -4.336576e-01 -2.279470e-01 -1.700519e-01 -6.006052e-01 -1.355176e-01
##  [921] -3.732170e-02  2.424470e-02  1.685151e-01 -4.153754e-01 -1.197023e-02
##  [926]  8.962091e-03  5.697123e-01  1.728744e-01 -3.952046e-01 -7.115587e-02
##  [931] -1.314262e+00  6.203638e-01  1.110593e-01  1.133321e+00 -1.467556e-01
##  [936] -1.554291e-01 -8.479998e-01 -2.105238e-01 -2.566467e-01 -1.168758e-01
##  [941] -1.837559e-02 -2.306450e-01 -4.346687e-02  6.668996e-01 -2.581454e-01
##  [946] -7.055966e-02 -5.656775e-01  2.238634e+00 -3.688869e-01 -5.285993e-01
##  [951] -6.960047e-02 -7.361498e-01 -3.673026e-01  5.796158e-02 -4.292214e-01
##  [956]  1.219596e+00 -1.272991e-01  3.592794e-01 -9.814901e-01 -3.570451e-01
##  [961] -9.165966e-01  3.139017e-02 -2.502727e-01 -8.914466e-01  2.431875e-01
##  [966] -9.360241e-01 -7.347207e-01  1.039800e-01  5.779533e-01  1.457383e-01
##  [971] -5.631653e-01  2.245332e-02  2.278225e-01  4.278669e-01 -5.615401e-01
##  [976]  3.630150e-01 -2.082338e-01  7.573933e-01 -9.502162e-02  1.848253e-01
##  [981]  6.339914e-01  2.509990e-01 -1.654958e-02  1.200651e-01 -2.144307e-01
##  [986] -8.713626e-01  3.809615e-02  3.225183e-01 -6.049539e-01  4.115521e-01
##  [991] -1.680892e-01 -4.598960e-01 -9.422515e-01 -4.099210e-01  3.267827e-01
##  [996]  6.379050e-02 -1.073208e+00 -3.535445e-01  1.268441e-01 -5.446842e-01
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
##   0.69150694   0.29614282 
##  (0.09364858) (0.06621637)
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
## [1]  0.70029521 -0.31021014  0.37646897  0.22254572  0.67190115 -0.08515691
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
## [1] -0.0157
```

```r
sd(results2$thetastar)  #the standard deviation of the theta stars is the SE of the statistic (in this case, the mean)
```

```
## [1] 0.9248128
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
## t1*      4.5 0.01621622   0.9276953
```

One of the main advantages to the 'boot' package over the 'bootstrap' package is the nicer formatting of the output.

Going back to our original code, lets see how we could reproduce all of these numbers:


```r
table(sample(x,size=length(x),replace=T))
```

```
## 
## 0 1 3 4 5 7 8 
## 1 1 2 1 2 1 2
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
## [1] -0.0039
```

```r
se.boot
```

```
## [1] 0.8975524
```

Why do our numbers not agree exactly with those of the boot package? This is because our estimates of bias and standard error are just estimates, and they carry with them their own uncertainties. That is one of the reasons we might bother doing jackknife-after-bootstrap.

The 'boot' package has a LOT of functionality. If we have time, we will come back to some of these more complex functions later in the semester as we cover topics like regression and glm.

