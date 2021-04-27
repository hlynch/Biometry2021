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
## 0 1 4 6 7 9 
## 1 3 1 3 1 1
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
## [1] 0.026
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
## [1] 2.833566
```

```r
UL.boot
```

```
## [1] 6.218434
```

Method #2: Simply take the quantiles of the bootstrap statistics


```r
quantile(xmeans,c(0.025,0.975))
```

```
##  2.5% 97.5% 
##   2.9   6.3
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
##    [1] 3.6 3.5 3.7 4.9 3.6 3.3 4.6 5.5 5.0 4.8 4.5 4.2 3.8 3.8 5.1 5.3 6.0 3.8
##   [19] 4.1 4.7 5.2 4.7 5.9 3.7 3.9 3.6 3.3 6.1 4.9 4.0 5.8 5.5 4.2 4.0 3.8 4.2
##   [37] 5.3 4.4 4.6 5.0 5.3 6.3 6.3 4.5 6.9 4.8 3.7 4.1 3.2 4.9 3.5 4.2 4.5 5.4
##   [55] 4.2 3.2 3.8 4.1 4.4 5.8 4.2 4.6 2.5 5.1 3.4 5.9 4.3 3.9 5.3 4.6 4.1 5.0
##   [73] 4.1 4.2 3.4 5.1 3.4 5.7 3.6 2.8 3.8 3.1 5.2 5.4 4.6 4.3 4.6 3.1 2.8 4.7
##   [91] 6.5 4.8 3.7 3.9 3.1 5.6 5.3 5.9 4.6 2.8 4.9 4.7 4.0 5.2 4.3 4.2 2.7 5.4
##  [109] 4.8 5.0 5.9 3.3 4.7 5.3 5.0 4.8 4.8 3.6 5.1 3.7 5.0 3.5 4.6 5.9 5.4 4.0
##  [127] 4.9 4.0 3.7 3.5 5.4 5.2 4.8 4.5 3.2 4.9 4.8 4.4 6.1 4.4 4.3 4.2 3.6 5.6
##  [145] 5.2 4.1 3.1 4.4 4.4 7.1 3.4 5.3 2.7 4.2 5.6 4.4 3.6 4.1 5.4 4.3 5.0 4.7
##  [163] 5.1 4.7 3.4 4.4 5.0 5.1 4.8 5.3 4.7 3.4 4.6 4.8 3.6 4.7 5.8 3.0 3.9 5.0
##  [181] 5.1 2.7 5.2 5.1 4.7 5.5 5.6 6.1 4.4 5.2 4.4 3.5 5.1 4.8 4.2 5.8 4.2 4.6
##  [199] 4.3 4.7 3.7 5.6 5.6 5.2 2.4 5.9 4.5 3.1 4.3 4.5 3.7 3.5 5.6 3.9 4.5 4.0
##  [217] 4.0 3.0 3.8 5.7 4.4 5.6 4.1 3.7 3.4 5.3 3.6 5.3 4.4 4.5 5.5 3.3 4.4 5.1
##  [235] 2.8 4.9 5.1 4.5 3.2 3.9 4.8 3.8 4.2 4.9 3.6 4.5 2.4 5.5 4.2 3.4 4.0 5.8
##  [253] 5.3 5.5 4.6 4.6 4.4 5.6 4.4 5.4 5.2 4.9 3.9 5.6 4.3 3.0 2.8 3.9 6.1 4.2
##  [271] 5.5 4.9 4.2 4.4 5.3 3.9 4.3 5.9 4.9 5.4 4.1 5.8 4.6 3.5 5.6 4.2 4.4 5.0
##  [289] 4.1 4.7 4.1 5.5 4.5 2.9 5.1 4.6 6.1 4.6 3.7 7.1 4.1 5.3 6.1 4.0 4.6 5.3
##  [307] 3.7 5.4 4.0 3.6 5.0 4.5 3.3 5.4 5.9 5.1 4.5 3.0 4.0 4.4 5.3 5.8 3.3 2.9
##  [325] 4.0 5.3 4.7 5.8 4.6 4.3 5.1 3.4 4.6 5.0 5.0 4.5 4.0 6.3 4.1 3.4 3.6 6.4
##  [343] 3.4 5.0 5.4 3.7 3.8 3.9 2.7 3.7 4.7 4.5 4.4 3.2 5.0 5.4 5.1 4.4 4.9 5.3
##  [361] 4.4 5.0 4.2 5.0 2.3 4.7 5.3 4.0 3.9 3.6 5.1 3.7 4.4 3.9 3.9 5.8 4.8 6.3
##  [379] 4.8 5.5 5.1 4.8 4.5 4.6 4.9 3.9 4.0 4.4 5.5 4.5 3.8 5.1 4.9 5.1 5.7 4.4
##  [397] 3.0 4.3 3.8 5.0 5.4 6.1 4.8 4.7 3.1 4.4 5.5 4.3 3.9 5.1 2.6 5.1 3.5 3.7
##  [415] 5.9 6.0 6.0 4.4 5.1 3.9 3.4 5.8 5.5 3.1 4.0 4.8 4.2 5.6 4.4 5.0 5.6 4.6
##  [433] 5.4 4.3 6.0 4.1 5.3 5.3 5.8 4.5 5.4 3.7 3.5 4.7 5.4 5.4 4.4 4.7 5.7 4.2
##  [451] 6.2 3.0 4.9 5.0 3.0 3.8 3.7 3.9 5.7 4.6 5.3 5.1 6.5 3.1 4.3 5.6 2.5 4.3
##  [469] 2.7 4.4 6.6 4.5 3.4 4.9 5.5 4.7 6.2 3.8 4.0 3.3 3.4 3.5 4.4 4.2 5.4 5.2
##  [487] 3.8 5.0 3.3 2.9 4.8 3.9 4.9 5.4 5.7 3.4 3.2 5.0 4.5 6.2 6.0 7.2 6.5 2.7
##  [505] 5.9 2.7 5.3 5.5 5.6 2.8 5.0 4.3 4.4 4.5 3.4 5.5 4.0 3.3 5.2 4.7 4.4 5.2
##  [523] 3.7 4.7 5.2 5.5 5.5 2.9 5.5 5.4 5.4 5.1 3.1 4.8 4.6 6.2 6.1 5.5 4.2 4.9
##  [541] 5.9 6.0 4.8 3.7 4.0 4.2 4.5 5.3 3.9 3.2 3.9 4.6 5.2 4.0 2.7 4.1 3.7 3.8
##  [559] 5.6 3.0 4.2 4.6 3.9 3.4 4.3 3.5 3.4 4.8 3.4 5.5 5.1 4.1 3.7 5.4 4.6 3.2
##  [577] 5.7 4.5 5.0 5.0 3.7 3.8 3.8 5.5 4.6 4.4 3.0 5.5 4.5 4.7 3.2 4.6 4.0 5.6
##  [595] 3.2 2.9 6.8 3.1 4.7 4.7 4.5 5.3 5.1 5.7 5.1 4.8 5.1 6.1 3.9 2.8 4.5 3.3
##  [613] 4.4 4.5 4.8 3.7 4.2 3.4 3.6 4.5 5.5 4.9 5.1 5.1 4.6 3.0 4.5 3.9 5.8 3.7
##  [631] 4.8 5.3 3.9 5.1 6.3 5.5 4.7 4.2 3.4 5.1 4.4 4.1 3.8 4.9 5.5 4.5 4.9 5.4
##  [649] 3.9 3.3 6.5 5.0 4.6 4.6 4.9 4.3 4.0 5.0 3.9 5.0 6.7 4.8 4.0 2.6 3.8 4.9
##  [667] 5.3 4.3 3.5 4.8 5.0 4.8 5.4 4.4 4.3 1.8 4.5 5.4 4.6 3.6 3.6 5.8 3.4 5.5
##  [685] 4.8 4.0 3.2 4.5 4.8 5.3 3.8 5.0 2.7 3.3 5.6 5.5 4.8 4.8 3.9 4.7 3.5 5.7
##  [703] 5.2 4.5 3.5 3.7 5.1 4.3 4.6 4.7 4.1 4.6 3.7 3.6 5.3 5.8 5.8 4.4 4.2 4.7
##  [721] 4.8 3.8 5.4 3.2 3.6 4.4 4.6 5.0 4.0 4.2 5.6 3.8 3.3 4.0 5.5 5.7 5.0 4.5
##  [739] 5.3 5.3 4.2 4.2 5.3 2.9 5.1 4.5 3.5 4.5 3.6 6.5 2.3 2.9 5.4 3.4 5.1 1.8
##  [757] 2.5 4.0 3.7 3.9 3.5 4.0 4.2 5.9 5.2 4.6 4.6 5.0 3.9 6.5 5.2 5.5 4.3 5.1
##  [775] 4.9 5.9 4.1 4.2 5.1 2.9 3.7 5.5 5.4 4.4 3.9 3.8 3.0 4.5 3.5 4.3 6.1 3.0
##  [793] 4.3 4.7 2.9 4.0 4.9 5.5 3.9 3.7 4.5 5.3 3.7 5.3 5.0 3.3 5.1 4.9 4.4 5.2
##  [811] 4.7 5.4 4.5 4.7 3.6 7.3 5.2 4.6 5.5 5.1 4.8 4.5 3.8 3.1 5.1 4.4 3.8 6.8
##  [829] 5.7 4.2 4.8 5.2 4.2 4.6 4.7 4.3 4.0 6.0 3.9 2.7 4.4 4.3 4.7 3.7 4.8 4.8
##  [847] 3.5 4.2 6.1 4.8 5.0 5.1 4.4 3.2 3.9 3.9 2.3 3.6 4.9 2.8 5.8 3.7 4.8 5.3
##  [865] 4.3 4.7 4.9 3.3 4.7 3.9 5.5 4.3 6.0 2.9 3.9 5.8 2.5 5.0 3.9 3.2 3.5 4.8
##  [883] 4.0 4.0 3.9 5.4 4.2 4.1 4.3 4.5 4.8 4.5 4.4 3.4 4.2 5.5 5.7 4.2 3.7 5.6
##  [901] 3.9 5.7 3.8 5.8 4.8 3.7 4.3 3.4 5.2 5.2 5.4 4.3 4.7 3.5 4.0 4.8 4.2 5.7
##  [919] 4.2 4.1 5.7 2.6 5.0 3.6 3.3 2.8 4.2 4.2 5.2 2.1 5.0 4.3 2.8 4.8 5.1 3.8
##  [937] 3.3 3.4 6.2 3.6 4.1 5.0 4.9 3.6 3.8 4.9 6.4 3.3 4.3 5.3 4.4 5.7 5.0 6.0
##  [955] 5.8 4.5 4.0 4.4 4.6 4.9 5.5 2.1 2.9 4.8 4.4 3.7 4.9 3.6 4.4 5.8 4.9 4.3
##  [973] 6.6 3.6 3.8 4.0 5.1 4.4 2.8 3.5 5.0 3.0 4.3 4.7 3.5 4.6 4.6 4.3 3.4 4.7
##  [991] 5.6 5.0 5.1 3.6 5.8 5.3 6.1 5.0 5.4 6.1
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
##   2.7   6.2
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
##    [1] 6.0 4.2 3.9 4.3 4.1 3.3 5.7 6.2 4.3 3.2 5.0 5.8 5.4 5.5 4.2 4.5 4.8 3.3
##   [19] 5.5 5.2 4.1 5.0 4.8 3.5 3.6 5.6 3.7 4.4 4.7 4.9 3.8 4.2 3.2 4.1 2.9 4.5
##   [37] 5.6 4.6 4.6 5.9 4.0 2.8 3.9 3.7 4.3 4.7 2.9 4.9 4.5 5.0 5.2 4.8 3.8 5.5
##   [55] 6.5 4.3 4.2 4.1 4.8 4.7 3.1 3.7 5.4 4.3 3.3 4.3 5.4 4.3 3.6 3.3 6.0 4.7
##   [73] 4.6 5.4 5.0 4.9 1.9 4.2 4.9 4.0 4.1 3.9 6.4 4.5 5.4 5.1 2.1 5.9 3.8 4.3
##   [91] 3.8 4.2 4.4 3.0 3.9 3.3 4.9 4.0 4.2 4.9 5.5 4.0 5.0 4.9 4.2 4.4 4.6 2.7
##  [109] 4.3 3.3 3.2 4.5 5.2 3.9 4.9 5.3 3.8 3.5 2.9 4.0 4.5 4.4 5.1 4.5 4.3 3.8
##  [127] 5.9 5.0 3.4 3.8 5.0 4.6 5.0 4.0 4.3 5.7 4.7 2.8 4.3 4.6 6.3 6.4 3.6 4.1
##  [145] 4.2 4.9 3.0 3.6 4.7 3.8 4.2 4.0 4.7 4.5 4.2 6.1 4.4 3.1 2.9 3.7 5.6 5.2
##  [163] 6.7 4.7 4.8 5.1 3.7 4.7 5.5 4.1 4.0 3.8 5.3 3.3 4.7 4.4 4.3 5.0 4.0 5.4
##  [181] 4.9 4.3 4.0 4.4 5.4 6.0 5.6 4.0 5.3 4.8 3.8 4.7 3.9 6.1 5.0 4.2 5.3 5.1
##  [199] 3.2 2.4 4.4 3.6 4.9 6.4 5.6 5.2 4.1 3.6 5.5 4.3 4.9 5.1 4.3 2.5 3.2 3.7
##  [217] 4.6 4.5 3.6 3.2 6.5 4.8 4.3 4.2 4.7 3.8 6.0 4.3 4.4 2.6 5.6 4.9 4.1 4.7
##  [235] 3.3 4.8 4.4 2.2 5.1 3.6 5.1 5.5 3.9 3.9 5.3 5.2 5.4 3.2 5.2 5.2 3.1 5.5
##  [253] 6.1 4.5 5.8 4.6 3.2 5.0 5.0 3.8 4.1 3.3 5.1 3.6 3.1 4.9 5.0 5.1 4.3 4.8
##  [271] 2.4 2.9 4.9 3.0 4.0 2.8 5.5 4.0 3.9 2.1 4.0 3.4 2.7 4.9 6.6 3.8 4.6 4.8
##  [289] 6.5 4.2 6.1 4.1 6.5 5.2 5.4 5.2 4.6 3.6 4.5 4.1 6.5 5.0 4.6 4.9 4.2 5.3
##  [307] 5.1 4.8 4.8 5.0 5.5 4.5 4.3 4.1 3.4 5.2 4.4 4.1 3.8 3.7 5.0 3.6 4.0 3.9
##  [325] 2.5 4.2 3.5 4.6 2.3 5.6 4.6 3.4 5.2 3.9 4.3 4.7 5.3 5.0 5.2 3.9 4.9 3.9
##  [343] 3.6 4.7 5.5 5.2 4.1 3.0 6.3 4.1 5.0 4.4 4.3 4.2 6.6 4.0 4.5 5.6 5.8 4.0
##  [361] 5.3 4.9 2.6 5.8 3.8 5.7 4.2 6.1 5.1 2.5 3.3 6.3 5.5 4.4 5.2 4.3 6.6 5.5
##  [379] 5.7 3.7 4.4 6.1 4.2 4.2 4.5 4.2 5.0 5.2 4.4 4.7 4.5 5.3 3.4 3.2 4.8 5.3
##  [397] 5.3 5.3 4.8 4.1 2.2 2.8 2.2 4.1 2.7 4.1 3.9 4.0 4.2 3.7 4.1 4.1 2.8 4.1
##  [415] 4.1 3.5 6.4 6.2 5.1 5.7 5.1 4.9 5.3 3.2 3.6 5.0 3.9 4.4 4.4 5.0 4.3 4.3
##  [433] 2.6 5.3 4.6 4.4 3.4 5.1 3.8 5.0 4.1 3.7 4.0 4.6 4.7 5.2 5.3 4.8 5.3 4.6
##  [451] 5.1 3.6 5.7 3.5 2.8 4.2 5.3 4.2 5.5 5.0 4.5 5.2 5.3 5.4 5.8 5.8 6.0 5.0
##  [469] 5.5 2.3 3.8 4.8 2.9 4.0 5.2 5.1 5.0 5.0 4.1 4.0 4.9 5.5 2.9 4.1 5.2 4.3
##  [487] 4.8 3.4 5.8 3.7 3.7 4.1 2.8 3.7 3.6 4.0 4.0 3.4 4.4 4.5 3.4 4.3 5.6 3.9
##  [505] 5.5 5.1 5.6 3.5 6.1 5.0 4.2 4.4 3.7 5.1 3.8 4.9 6.0 4.4 4.9 4.2 4.7 6.3
##  [523] 3.1 5.1 4.0 2.4 3.4 5.3 4.3 3.3 4.7 4.0 3.4 4.4 2.6 5.3 3.7 4.4 3.6 4.2
##  [541] 4.5 3.3 4.0 5.0 4.6 4.9 3.1 6.5 3.0 5.3 4.7 5.2 5.0 4.2 5.8 4.3 4.9 3.9
##  [559] 5.5 3.4 4.6 5.8 5.3 4.4 3.4 2.7 4.3 3.9 3.9 5.6 4.2 5.1 4.0 5.8 4.4 7.0
##  [577] 4.7 3.9 3.4 3.4 3.6 5.3 6.4 4.9 4.4 4.9 4.0 5.2 4.8 4.1 5.0 3.9 5.0 4.6
##  [595] 4.4 3.4 5.0 3.8 5.7 5.2 2.6 7.2 4.0 3.4 5.2 5.0 5.1 3.7 5.6 3.8 4.4 5.2
##  [613] 4.2 5.0 3.3 4.8 3.3 3.8 4.2 6.2 4.0 4.7 5.4 5.5 6.5 4.0 2.5 4.6 6.2 4.2
##  [631] 6.0 5.6 4.1 3.8 4.7 4.7 3.0 4.7 5.7 4.7 5.1 3.4 5.2 2.9 4.9 5.2 4.0 4.6
##  [649] 3.2 4.5 3.2 6.1 4.4 1.9 5.0 4.1 7.2 4.1 4.0 3.7 4.2 6.3 4.6 6.4 4.3 5.4
##  [667] 5.0 4.1 5.0 4.1 4.2 6.5 5.5 5.1 5.7 4.1 6.1 4.7 3.1 4.0 5.5 3.4 5.3 3.6
##  [685] 3.4 2.7 4.6 5.4 4.4 4.7 5.2 3.7 3.4 5.8 4.8 3.5 5.4 3.8 4.3 5.1 2.9 5.9
##  [703] 4.5 4.2 4.6 4.1 5.5 7.1 3.4 5.7 4.4 3.8 3.2 2.9 3.8 4.4 4.5 5.3 4.3 4.9
##  [721] 4.8 3.4 5.2 4.5 5.2 5.1 7.0 5.2 5.2 4.1 2.7 3.4 5.1 5.9 4.0 5.2 4.1 4.6
##  [739] 4.8 3.7 4.1 4.0 4.6 3.6 6.5 4.0 4.6 4.9 4.3 4.9 4.0 5.2 5.1 3.2 5.4 4.7
##  [757] 4.1 3.3 5.3 3.9 4.3 5.5 5.4 4.5 4.8 4.6 5.3 4.6 5.3 3.2 6.1 5.0 4.9 3.4
##  [775] 4.7 4.0 5.1 3.1 4.5 4.5 4.7 2.5 4.9 4.1 4.8 4.6 3.9 5.1 5.4 4.8 5.3 4.2
##  [793] 5.6 3.5 3.9 3.1 5.2 3.9 3.6 3.9 5.1 6.0 4.0 5.7 5.4 5.2 6.6 5.9 6.2 5.3
##  [811] 4.7 5.5 4.1 5.3 3.7 3.7 5.5 4.4 5.6 2.3 3.2 3.1 5.1 3.7 5.1 4.2 3.9 5.2
##  [829] 4.9 4.8 4.1 3.7 4.1 3.5 4.0 4.2 3.5 5.3 3.3 4.4 3.9 4.2 3.2 4.9 3.2 4.7
##  [847] 5.7 3.1 5.1 3.3 3.0 4.6 4.3 3.6 3.6 4.3 7.1 4.6 5.0 5.0 3.9 6.8 4.2 4.7
##  [865] 3.0 4.6 4.5 5.2 4.0 4.5 5.3 6.1 5.1 3.9 3.0 4.9 4.9 5.0 3.4 6.2 2.6 5.4
##  [883] 4.0 5.0 4.1 3.3 4.3 3.8 3.4 5.0 5.2 5.2 3.7 5.1 5.4 5.6 5.5 4.4 3.5 4.2
##  [901] 4.8 3.6 4.5 5.1 5.1 6.0 4.4 4.3 3.1 5.9 4.6 3.8 4.5 4.8 4.3 5.2 3.6 4.0
##  [919] 4.9 3.1 3.9 6.0 3.9 4.0 4.2 4.3 4.8 5.6 3.1 5.0 4.5 4.2 5.4 4.9 5.1 5.0
##  [937] 4.7 2.0 4.5 4.0 3.9 4.5 5.3 4.4 5.2 3.9 4.5 4.1 4.0 4.3 4.7 4.1 4.4 6.0
##  [955] 4.9 3.9 3.9 3.5 5.2 4.4 6.8 5.5 4.6 3.7 3.5 5.6 5.4 5.3 4.6 3.7 4.5 3.0
##  [973] 4.9 5.5 3.2 3.0 3.9 4.2 1.7 3.9 3.1 5.2 4.5 4.3 6.1 5.0 5.4 3.3 6.9 5.4
##  [991] 3.8 5.7 6.0 4.7 3.9 3.9 4.9 5.2 4.0 2.5
## 
## $func.thetastar
## [1] -0.0091
## 
## $jack.boot.val
##  [1]  0.51563422  0.37203390  0.29596774  0.17119205  0.06126374 -0.07905028
##  [7] -0.13422619 -0.31005917 -0.40527704 -0.54153846
## 
## $jack.boot.se
## [1] 0.9961483
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
##    [1] 5.3 5.4 3.2 3.8 4.6 5.1 3.4 2.8 2.6 4.4 3.6 4.1 3.5 5.6 4.5 5.8 5.0 3.6
##   [19] 5.3 4.2 3.0 4.2 6.3 3.8 4.3 3.1 3.7 5.1 4.0 3.9 4.5 4.6 5.6 4.9 5.8 3.8
##   [37] 3.0 2.4 4.5 5.9 4.5 3.0 5.2 2.7 2.7 3.8 4.7 3.7 3.7 5.4 3.7 2.7 4.0 4.2
##   [55] 4.1 4.1 6.2 3.7 4.8 5.2 4.2 4.9 3.4 3.7 3.8 2.6 4.4 5.0 5.9 4.6 5.7 3.8
##   [73] 3.5 4.0 6.2 3.5 4.8 3.7 5.4 4.8 2.5 3.8 3.9 4.8 5.0 3.3 5.6 5.6 4.6 3.3
##   [91] 4.2 4.3 6.3 2.9 5.1 3.2 3.5 4.8 4.8 3.8 3.7 4.9 3.4 3.8 3.8 4.6 3.5 4.7
##  [109] 3.6 4.6 5.8 3.7 5.0 4.1 4.5 5.3 5.6 5.5 4.6 3.8 4.9 4.7 2.7 3.3 4.7 3.3
##  [127] 4.7 3.7 5.3 5.9 7.0 3.0 3.4 3.6 3.7 4.6 5.3 3.9 5.0 2.8 5.8 5.4 5.4 4.8
##  [145] 5.9 5.5 5.1 5.4 6.0 5.7 6.5 3.6 5.3 4.0 5.0 4.8 2.9 4.3 3.6 4.0 3.0 4.9
##  [163] 5.3 3.7 3.7 4.9 6.3 3.8 5.6 4.1 3.4 3.0 4.6 5.0 4.9 6.0 3.7 4.1 5.6 5.8
##  [181] 5.2 4.8 5.0 4.7 5.8 4.5 3.8 6.8 3.5 3.2 4.8 5.0 5.4 4.6 5.9 4.9 4.2 3.5
##  [199] 4.1 5.5 5.2 4.3 4.3 5.6 4.4 4.4 3.7 3.0 4.9 5.5 4.6 3.3 3.9 5.4 3.8 2.6
##  [217] 4.8 4.0 4.4 4.5 5.0 4.8 4.1 4.3 4.1 3.7 4.4 4.3 3.2 3.4 3.9 3.7 4.0 4.0
##  [235] 4.3 4.8 5.3 3.3 2.9 4.6 3.8 5.2 4.1 5.8 4.0 4.5 5.3 4.1 4.6 3.9 4.8 5.0
##  [253] 4.4 5.0 3.4 6.0 6.0 4.1 5.7 5.8 5.3 5.3 5.2 5.0 3.5 3.4 3.9 4.5 2.6 4.9
##  [271] 5.1 4.6 2.4 4.0 3.9 3.5 4.3 4.0 3.6 3.6 2.7 1.8 3.8 6.5 3.9 3.5 4.0 4.6
##  [289] 4.7 5.0 2.8 4.6 3.6 4.2 4.4 3.5 3.8 6.7 4.0 5.7 3.8 5.9 4.1 3.0 5.0 5.3
##  [307] 6.0 3.7 3.4 4.2 3.9 5.2 5.1 4.9 3.1 4.7 3.9 4.1 3.9 4.6 4.6 2.8 3.4 3.7
##  [325] 5.0 4.1 5.1 3.8 3.9 3.3 4.3 4.2 4.3 4.4 3.3 5.3 4.9 5.7 4.6 3.6 5.5 4.2
##  [343] 4.0 5.3 5.2 5.7 4.5 3.8 5.5 4.2 4.8 4.1 3.7 4.1 4.5 3.9 4.1 4.4 4.3 5.6
##  [361] 5.7 4.9 3.4 5.9 3.9 1.8 3.6 3.3 4.5 2.9 4.9 6.2 3.8 7.3 4.3 4.5 3.5 4.3
##  [379] 6.1 5.2 3.6 3.4 5.3 5.7 2.8 4.5 4.7 3.9 3.5 4.2 5.8 4.5 4.5 3.0 4.8 5.5
##  [397] 5.4 3.5 3.9 3.4 5.1 4.8 5.6 4.9 4.0 5.8 4.6 2.7 4.5 5.0 5.2 4.6 3.4 3.6
##  [415] 2.9 4.1 5.3 3.6 4.1 4.4 4.6 4.9 5.4 4.0 5.4 4.1 3.4 4.0 6.2 5.0 5.9 4.0
##  [433] 4.1 5.1 3.5 4.8 3.3 3.3 4.9 4.6 4.1 2.9 4.0 4.9 3.6 6.0 4.4 2.4 3.1 6.1
##  [451] 2.9 4.6 1.6 3.1 4.8 5.8 4.6 4.0 4.5 3.2 3.8 5.2 3.5 4.9 4.6 4.3 4.3 4.3
##  [469] 4.0 2.7 4.5 3.6 6.0 4.2 5.4 4.0 6.0 2.8 4.4 3.7 3.4 5.7 3.8 4.2 4.0 5.5
##  [487] 5.3 4.1 3.2 3.8 3.4 3.6 5.6 4.7 2.9 5.0 5.7 2.7 4.7 4.9 3.9 3.4 5.4 6.1
##  [505] 6.1 3.6 5.3 3.6 5.9 5.2 4.2 4.0 3.9 5.1 4.8 4.0 4.7 4.5 5.4 4.8 4.0 5.4
##  [523] 4.5 5.9 4.7 5.3 2.7 4.8 6.1 4.8 4.4 5.6 4.8 3.4 5.1 3.7 6.1 5.2 4.3 6.0
##  [541] 3.3 4.7 3.9 3.9 5.0 5.5 4.9 4.7 5.7 3.8 4.4 5.4 5.7 4.9 5.0 3.8 4.9 3.6
##  [559] 3.7 4.6 3.7 4.6 5.0 4.9 4.3 2.7 5.4 4.5 4.3 5.8 5.1 3.5 4.6 4.6 2.6 4.0
##  [577] 5.5 2.4 4.1 4.6 5.7 5.8 5.0 5.5 3.4 3.7 5.8 4.0 5.3 5.4 4.3 4.4 2.7 3.8
##  [595] 5.7 4.4 2.8 5.2 5.2 6.4 4.2 5.0 5.1 5.2 5.2 5.9 2.8 4.7 5.6 6.2 4.9 4.3
##  [613] 3.1 3.8 5.2 4.2 3.6 4.6 4.6 4.9 5.0 6.2 5.3 4.4 4.4 4.5 6.1 4.3 4.5 4.4
##  [631] 4.6 5.2 3.1 4.2 4.2 5.1 6.1 4.3 5.3 5.1 4.1 2.9 5.3 5.1 3.5 4.1 4.2 4.2
##  [649] 5.2 2.5 5.7 4.9 4.5 5.5 5.7 4.2 4.0 5.6 3.5 5.2 3.9 5.1 2.9 4.2 4.4 4.7
##  [667] 3.5 6.0 5.6 4.2 3.7 5.3 4.8 4.2 5.4 3.4 3.3 4.5 5.0 6.7 3.3 4.1 3.6 3.4
##  [685] 4.2 3.7 3.8 6.4 4.9 4.6 4.4 3.2 5.2 5.5 2.9 4.0 4.3 3.3 6.4 4.6 3.6 4.7
##  [703] 3.8 4.7 4.6 3.9 5.6 5.7 3.0 3.9 5.0 3.6 4.8 2.5 5.8 5.7 3.8 3.0 4.1 3.1
##  [721] 6.3 3.8 4.4 3.9 6.0 3.4 4.0 4.5 4.9 3.7 4.8 5.3 3.8 4.6 5.9 4.9 2.3 4.2
##  [739] 4.0 3.9 2.8 4.2 4.4 4.5 5.2 4.7 5.4 4.4 4.8 5.1 4.5 4.9 4.2 3.3 5.1 5.5
##  [757] 3.0 6.0 3.6 4.4 3.0 6.3 5.2 4.6 2.6 3.1 5.6 4.9 3.8 5.4 5.6 5.0 5.2 4.5
##  [775] 3.2 3.6 5.8 3.7 4.8 4.3 5.6 5.5 3.0 3.5 3.0 5.2 5.0 2.7 3.6 4.7 4.3 4.7
##  [793] 5.7 5.8 5.8 4.8 4.7 5.2 3.9 2.0 4.9 4.1 4.2 4.9 3.8 4.4 5.0 4.8 4.2 4.4
##  [811] 4.4 5.8 4.9 3.3 3.7 3.6 3.4 5.0 3.5 4.4 3.0 6.1 3.2 4.5 5.2 4.6 5.4 6.3
##  [829] 5.4 4.4 4.3 4.3 4.2 4.6 6.0 3.9 4.0 5.0 4.3 5.1 3.6 3.2 3.9 5.1 4.4 4.0
##  [847] 4.1 4.1 4.9 3.1 5.2 5.6 5.9 4.7 3.7 4.9 5.2 4.5 3.1 5.2 5.0 5.5 4.0 4.6
##  [865] 3.7 5.6 5.4 4.8 4.5 4.6 5.0 4.4 4.9 6.7 4.2 3.6 4.0 3.4 5.2 2.8 4.5 5.2
##  [883] 5.6 4.2 3.1 4.4 4.3 4.8 3.8 4.1 5.2 4.0 4.8 3.2 5.7 2.9 3.7 4.6 3.9 4.1
##  [901] 3.6 5.0 4.9 5.2 5.4 4.2 5.5 3.4 4.2 3.1 4.7 2.2 4.0 5.8 4.4 4.8 4.0 3.9
##  [919] 4.1 4.8 5.2 5.6 3.7 5.5 5.6 4.8 4.3 4.3 4.2 3.0 3.3 4.4 5.4 4.7 4.0 4.7
##  [937] 4.2 3.4 5.1 4.9 5.9 4.4 3.0 3.7 5.0 3.7 4.9 4.3 4.4 3.7 3.3 2.6 5.1 4.4
##  [955] 5.1 5.8 4.6 3.1 6.0 4.0 3.7 2.6 6.3 4.1 5.0 2.7 3.7 3.8 3.8 4.8 4.8 2.2
##  [973] 3.8 3.6 4.4 3.6 4.4 3.1 4.6 4.5 2.5 3.5 3.3 3.2 5.2 4.0 6.8 4.7 3.5 4.4
##  [991] 5.6 4.8 5.6 4.1 6.3 3.1 5.7 4.4 5.9 4.6
## 
## $func.thetastar
## 72% 
##   5 
## 
## $jack.boot.val
##  [1] 5.58 5.40 5.30 5.20 5.02 5.00 4.80 4.80 4.50 4.40
## 
## $jack.boot.se
## [1] 1.088632
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
## [1] 0.6754559
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
##   4.102150   8.382083 
##  (1.765014) (3.836631)
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
## [1] 0.3785263 0.4193977 1.0925741 0.2471757 1.4523682 0.9546382
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
##    [1]  1.110005348  0.138846676  1.600911458  0.465070892  0.694162637
##    [6]  0.679750976  0.689214540  0.817771552  0.367027468  0.166222363
##   [11]  0.683258600  0.570815421  1.873417762  0.572070072  0.039974497
##   [16]  0.892474922 -0.231543868  1.138647069 -0.130099070  0.802194211
##   [21]  0.554139516  0.728352836 -0.294183491  0.528476292  0.925950199
##   [26]  0.048171763  0.425393093  0.270523843  1.577701996  0.127340799
##   [31]  0.671639939  0.923775320  0.485745261  0.297537355  0.272265054
##   [36]  0.841275707  0.367758389  0.631968654  1.888261398  0.859825013
##   [41]  0.446149248  0.942846084  1.084822646 -1.123372747  1.331483381
##   [46]  0.077517729  0.704865282 -0.264210500  0.539114497  0.543727777
##   [51]  0.668798286  0.003481200  1.066896286 -0.017836102  0.270234138
##   [56]  0.805427916  0.023382459  0.973520975  0.650029140  0.692940573
##   [61]  0.450824940 -0.185241580  1.251710803  0.300017634 -0.030270955
##   [66] -0.011161987  1.274625394  0.842009853  0.583278332  1.301508684
##   [71]  0.440246954  0.312825211  1.221239252  1.012489737  0.044091240
##   [76]  0.891600249  0.832344180 -0.111870408  1.067430713  0.629239373
##   [81] -0.304360220  2.119404203  0.408978933  0.554135578  0.657364098
##   [86]  0.007694652  1.048365226  1.026911788  0.445871587  0.048421321
##   [91]  0.767331112  0.441536077 -0.081366308  1.031154867  1.075949048
##   [96]  0.730478566 -0.003275071  0.392698268  0.666815501  1.190118936
##  [101] -0.382396032  0.386852005  0.280329208  0.604176403  0.789896025
##  [106]  0.045721245 -0.036181668  1.551009510 -0.402855982  0.628312117
##  [111]  0.951201626  0.551286997  1.312991968  1.235219266  0.749069744
##  [116]  1.087564230  1.211982431  0.408460955  0.541305733  0.746187174
##  [121] -1.022830776  0.038916944  0.833417132  0.034434333  0.780918788
##  [126]  0.701953813  0.636619450  0.847334503 -0.404501728  0.786801045
##  [131]  0.730330609  1.139933768 -0.019115654  0.235634614  0.297503293
##  [136] -1.252775230  0.633803145 -0.032635867  0.389597810 -0.788401580
##  [141]  1.147412618  0.592789518 -0.512603290  0.726275009  1.291613367
##  [146]  0.367355431  0.637220671  0.769487693  0.038495539  0.856112515
##  [151]  1.048338749  1.160312989  0.223737308 -0.350781828  1.599671197
##  [156]  1.860258208  0.455279049  0.789770568  2.043705056  0.455327791
##  [161]  1.105307862  0.437413604  1.534746768  0.451654382  0.295309724
##  [166]  1.242543148  0.379369725  1.268812671  0.861771467  0.797052177
##  [171]  0.731285325  0.758317748  0.635827267  1.118225622  0.839113052
##  [176]  0.254884740  0.327554721 -0.813081339  1.225875181  0.259570476
##  [181]  0.822293738  1.201379015  1.046593130  0.838515920  0.895130708
##  [186]  1.014280817  0.410875535  0.782937577  0.247288941  0.644787926
##  [191]  0.328840060  0.584373008  1.272061074  1.533307168  0.758982923
##  [196]  0.687596695  0.573999462  1.116669030  0.340704510  0.867442881
##  [201]  1.117748238  0.298501257  0.571184799  0.735496912  0.408456544
##  [206]  0.675455940  0.413306962  0.436789664 -0.276543653  0.714325624
##  [211]  1.603957053  1.262600438  0.274120380  0.531257777  0.634806007
##  [216]  0.414481398  0.014883330  1.798492941  0.770158075  0.553390106
##  [221]  0.604692919  1.063652073 -0.800587550  0.495857770  0.638512223
##  [226]  1.567506900 -0.343231644 -0.180399485  0.398032013  0.943183840
##  [231]  0.635318923  0.302108735 -0.069026453  0.513510116  1.264942181
##  [236]  1.067147102  0.003827363  1.337840338  0.354754964  1.184360252
##  [241]  1.022645161  0.586581714  1.215439239 -0.409543162  0.245470967
##  [246]  0.693350467  0.652060193  0.732101844  0.856927665  0.344499861
##  [251]  0.760893076  0.045672061  1.611716282  0.283415882  0.783520352
##  [256]  0.554272470  0.561999397  1.381788756  0.076059779  0.068220618
##  [261]  0.413796194  0.711042738  0.017881834  0.449311970  1.269270232
##  [266]  1.155314178  0.558457873  1.244968298 -0.087578085  0.988106790
##  [271]  1.146583222  0.726748909  0.731648860  0.959467631  0.610839674
##  [276]  0.876634371  0.978724456  1.303889215  1.399824787  0.688868769
##  [281]  0.669001159  0.408323668  0.977146692  0.589437603  2.033187870
##  [286] -0.065941889  1.418926654  0.166857586  1.210542262 -0.679838815
##  [291]  0.588671224  1.269286525  0.004144966  0.441022870  1.470629261
##  [296]  0.822945776 -0.332986785  0.292539443  1.020470219  0.544636656
##  [301]  0.458460666  1.012699501  0.289149948  0.391880211  0.112271151
##  [306]  1.120162477  0.580736240  0.033102086 -0.211473629  0.762402405
##  [311]  1.143841423  0.654157464  0.493791231  0.387774991  0.736486895
##  [316]  0.170112252  0.401905125  0.093161550  1.916687089  0.166222363
##  [321]  0.602567429  1.783274875  0.826345116  0.903367032  1.460289654
##  [326] -0.275634546  0.054962364  0.847045130 -0.306035333  0.111683239
##  [331]  0.340512205  0.569614114  1.886494118  0.035123621  0.814951869
##  [336]  0.889885298  1.627726676  0.689489043  1.137079836  1.824362113
##  [341]  0.446567723  0.348826029  1.934550400  0.010171329  0.371528832
##  [346]  0.722017382  0.682592113  0.026486821  0.442533926  0.411363125
##  [351] -0.471889466  1.248809784  0.903433028  0.244699688  0.592931972
##  [356]  0.671165048 -0.003461782  1.061775242  0.553746697  0.327606584
##  [361]  0.768470185  1.768123577  0.774090021  1.502550005  1.243426783
##  [366]  0.578827697 -0.006654591  0.829753594  0.722504896  0.410018494
##  [371]  0.394628762  1.257875886  0.772929020  0.033556298  0.716659108
##  [376]  0.580206887 -0.022865017  1.116218177  0.600078770  0.745743741
##  [381]  0.653028962  0.720094400  1.105409554  0.538359178  0.377412659
##  [386]  1.161584712 -0.387477036  1.725412700  0.258364563  0.868863596
##  [391]  0.607363537  0.733427993  0.357318037  1.086251994  0.633697845
##  [396] -0.357497667  0.843526817  0.505607538  0.624497553  0.359414670
##  [401]  0.373601947  1.276687522  1.121824032  0.963076496  0.404545158
##  [406]  1.639436406  1.069365127  0.292063164  0.715645650  0.732745830
##  [411]  0.750233434  0.342773966  0.740219551  0.453549131  0.954337788
##  [416]  0.392754459  1.236994234  1.288514246  0.360167704  0.020587239
##  [421]  1.221309990  1.877985766  1.020934104  0.313771150 -0.288301640
##  [426] -0.492140483  0.736912337  1.391528617  0.705656949  0.374210602
##  [431]  0.328774977  0.781029308  0.238721863  0.357747942 -0.193453575
##  [436]  1.172663795 -0.396119978  0.990874302  0.407254207  1.117425309
##  [441]  1.147356337  1.688157556  0.371281479  1.151055163  0.677493262
##  [446]  0.648151480 -0.034744093 -0.279551620 -0.017762298  0.970794176
##  [451]  0.346468680  1.435080347  1.522546731  0.488297848  1.323865726
##  [456]  1.536261140  0.355873345 -0.319191166  1.302608682  0.576665232
##  [461]  0.650831790  1.390402891  0.668257867 -0.940994578  0.099065056
##  [466]  0.647906623  0.454354246  0.519246653  0.226119055  1.145002058
##  [471]  0.297555783  0.721288270  1.035954837  0.378763947  1.329032875
##  [476]  1.320889242  0.494373093  0.587790404  0.774167455  0.643354366
##  [481]  0.815779235  0.666086359  0.628762270  0.451959600 -0.078983140
##  [486]  1.209501810  0.678108784 -0.704589695  0.406725711  0.980518581
##  [491]  0.349811341  0.418011876  0.699732504  0.551096018  0.111486080
##  [496]  0.732657618  1.173593271  1.135907737 -0.025246463  0.511245502
##  [501] -0.627235110  1.148323860  0.998776683  0.684018579  0.748842371
##  [506]  0.006428235  0.719020731 -0.284178872  0.331518220  0.593882180
##  [511]  0.736766258  1.400319783 -0.687412356  0.784580866 -0.070234163
##  [516]  0.996486355  0.756744069  0.353751170 -0.618796015 -0.060281436
##  [521] -0.300206895  0.710513164 -0.393640722  1.216612941  0.346601967
##  [526]  1.416492192  0.604224107  0.782694354  1.147975143  0.844196814
##  [531]  0.841117758  0.488996882  0.359928942  0.307059017  1.595614305
##  [536]  1.357786615  1.293042120  0.628169053  0.898179829  0.592789518
##  [541]  0.019751603  0.646801235 -0.717641037 -0.377323003 -0.039056138
##  [546] -0.369864069  0.588501612  0.698552822  0.381855606  0.441084062
##  [551]  1.583565269  0.711394843  0.657728225  1.264464629  1.001469709
##  [556]  0.996499063  0.749991292  0.852187713  0.647355216  0.379994476
##  [561]  1.270958439  1.618372300  1.309198728  1.085145783  2.213709440
##  [566]  1.455680220  1.622318537  0.393629555  0.163851179 -0.264249877
##  [571]  1.988318127  0.229758091  0.023334216  0.866345784  1.053329939
##  [576]  1.334931141 -0.683228463  1.058476785  0.715064713  0.320559018
##  [581]  0.681511636  1.253127068  1.142142195  0.084444252  0.008049599
##  [586]  0.692873923  0.672535018 -0.255079180  1.030395601 -0.077517713
##  [591]  0.672614090  0.653223061  0.351864580  0.690528310  1.004022417
##  [596]  1.028758640  0.577890583  1.177873635  0.015794849  0.961170589
##  [601]  0.443594288  1.376683547  0.886259888  0.378763947  0.660250875
##  [606]  0.981448145  0.403456255  0.461414002 -0.178323877  0.384276344
##  [611]  2.055109549  0.292808043  0.600730342  1.247315514  0.796093633
##  [616]  1.334264270  0.678569141  0.299825378  0.585538772  1.796422309
##  [621]  1.731532060  0.614225326  1.397303721  0.581492409  0.863749492
##  [626]  0.255663302  0.794836404  0.717671304 -0.529541617  0.770639945
##  [631]  1.731776352  0.797305954  0.953257375  0.269673035  0.979083107
##  [636] -0.337976561  0.768601217  0.685663870 -0.406759640  0.046929171
##  [641]  1.622043237  0.632082695  0.487280237  0.963394921  0.580072883
##  [646]  0.984363018 -0.104173104  1.173467799  1.246508165  1.244456093
##  [651]  0.145822409  0.057367684  0.246837132  0.048106672  0.190481565
##  [656]  0.056529931 -0.198511431  2.380168637 -0.125016249  0.732596211
##  [661] -0.056745202 -0.107251131  0.429131628  1.205088267 -0.048456192
##  [666]  0.129723142  0.841401816  0.337024688  0.001436213  0.692946480
##  [671] -0.594187794  0.771135462  0.635196451  1.357786615  0.763107168
##  [676]  0.943857389  1.153049347  0.736408713  0.674485210  0.561400064
##  [681] -0.089400008  0.287049103  0.731961850  0.260686679  0.743818020
##  [686]  0.534691834  0.539543359  0.623719117  0.395112044  1.703047293
##  [691]  0.480528175 -0.333848711  0.881353547  1.018081720  0.347691271
##  [696]  1.168337358  1.238211068  1.107622904  0.413796194  0.348198661
##  [701]  0.343053338  1.261129039  0.652150026  0.799929845  0.038804255
##  [706]  1.070061662  0.460822068  1.263662329  0.031164949  0.787261874
##  [711] -0.343388477  0.302168821  0.481925369 -0.012430629 -0.038366054
##  [716]  0.417553343 -1.210794889  0.259081068  0.693613230  1.101859811
##  [721]  0.020422832  0.266378333  1.602254335  0.587292965  0.008071140
##  [726]  0.762454439  0.841955491  0.403423878 -0.387407575  1.145815917
##  [731]  0.806446267  0.423997413  0.712719663  0.850773908  1.157008724
##  [736]  0.805343486 -0.407788874  0.622790880 -0.510659000  1.232530890
##  [741]  1.393825804  1.351933852  0.079147009  0.768350248  0.513686555
##  [746]  0.419496855  1.118626732  0.745570434  0.577883618  1.874724393
##  [751]  0.844755125  0.821511492  0.471300308  0.884512316  1.626696293
##  [756]  0.704562882  0.974638804  0.017802755  0.318193022  1.223100207
##  [761]  0.360020134  0.393973469  0.391598587  1.168590533  0.327263337
##  [766]  0.344763099  0.848069763  0.080960984  1.413828107  0.090239864
##  [771]  0.454221424 -0.011787554  1.156342657  0.490560392  1.114667537
##  [776] -0.746601321  1.146602568  0.598144304 -0.300755274  0.282040275
##  [781] -0.708960201  1.088431105 -0.050816899  1.615798564  0.190904755
##  [786]  1.253127068  0.168427634  0.997211037  0.870694546  0.303087677
##  [791]  1.156829901  1.159086384  0.334683805 -0.033764915  0.360595595
##  [796]  0.038558407  0.651914469 -0.125950723  0.948456401  0.694587737
##  [801]  0.604962654  0.121132748  0.767537361  1.759842507  0.715907474
##  [806]  1.306556524  0.363360445  0.927188174  0.193010608 -0.034130486
##  [811]  1.589147474  0.338298859  0.874225306  0.260688375  1.510647490
##  [816]  1.038767928  0.058975643  0.018386559  0.369332245  0.363351097
##  [821]  1.292523977  0.699302231  1.700667132  0.044170893  0.695956488
##  [826]  0.519473042 -0.036724320  0.503962659  0.649662715  0.445203119
##  [831]  0.826712712  0.130228295  1.351228891  1.359241079  0.671317836
##  [836]  1.161850308  0.282507788  0.463608913  1.177587336  0.294052614
##  [841]  1.094063552 -0.407720895  0.769530569 -0.185521545  1.290446429
##  [846]  0.493158683  0.533239305  0.347983379  0.289529558  0.644336256
##  [851]  1.223950353  0.750811714  0.542847566  0.252815757  1.337207274
##  [856]  1.546695951  0.383096709  1.161579959  0.762228668  1.368634198
##  [861]  0.848942631  0.421097000  0.825743227 -0.735220417  0.249158758
##  [866]  0.666574797  0.708582649 -1.095812690  0.387450004  1.056678093
##  [871]  0.605089521  0.930551667  0.669132350  0.596627149  1.088818956
##  [876]  0.926563589  1.222478247  1.290867395 -0.043883076  1.104058476
##  [881]  0.343398033 -0.127497548 -0.353461409  1.177656551  0.698855908
##  [886]  1.172173541  1.831876941  0.672646761  1.110206054  0.851674442
##  [891]  1.754318628  0.602817994  0.379656761  0.158856527  1.134737026
##  [896]  0.777787531  0.718769249  0.802334667  1.488229402  0.064355061
##  [901]  0.061334653  0.179774898  0.825168927  0.801679540 -0.320244669
##  [906]  0.737721182  0.440646334  1.089088091  0.412643152  1.476298641
##  [911]  0.650590210  1.638197566  0.700153350  0.260719252  0.409116347
##  [916]  0.664503204  0.670366009  1.084451504  0.997211037  0.369012287
##  [921]  0.802023069  1.469332824  0.406910279  0.507069607 -0.056757693
##  [926]  1.119060441  1.766526497  0.455473760 -0.308876908  0.711042738
##  [931]  0.646536315  0.669096801 -0.371058162 -0.700884544 -0.389939260
##  [936]  0.344948882  0.226921690  0.299861524  0.421680435  0.636004818
##  [941]  0.175558693  0.693350467  0.831974584  0.176742861  0.242994200
##  [946]  0.236514445 -0.006229784  0.314680629 -0.446736502  0.056333550
##  [951] -0.044557485 -0.018903245  1.877985766  0.446197144  1.492079235
##  [956] -0.409482906  1.217397591  0.672104658  0.063431485 -0.576210279
##  [961]  1.651077612  1.460298455  0.870925354  1.351881617 -0.061667725
##  [966] -0.371425088  0.740870400  0.886874869  0.255424049  0.395837375
##  [971]  0.602844090  1.170432878  0.769435796  0.094391324  0.536552531
##  [976]  0.712459571  0.359954182  1.859157657  0.662087311  2.026439039
##  [981]  1.535925366  0.413306962  0.663761809 -0.321737300  0.021138145
##  [986]  1.440900201  1.306614377  0.172048872  0.009617261  0.166994212
##  [991]  0.149909721  1.061443782  1.359241079  1.045380022  0.749400319
##  [996]  1.646376231  0.708457849  0.420496306 -0.023612523  0.707289049
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
##   0.48939371   0.24880845 
##  (0.07868014) (0.05563191)
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
## [1] -0.69401482 -0.33573881  0.58349006  0.46668012 -0.79406862 -0.03912669
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
## [1] 0.007
```

```r
sd(results2$thetastar)  #the standard deviation of the theta stars is the SE of the statistic (in this case, the mean)
```

```
## [1] 0.8855254
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
## t1*      4.5 0.04364364   0.9285049
```

One of the main advantages to the 'boot' package over the 'bootstrap' package is the nicer formatting of the output.

Going back to our original code, lets see how we could reproduce all of these numbers:


```r
table(sample(x,size=length(x),replace=T))
```

```
## 
## 0 1 2 5 7 
## 1 4 1 3 1
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
## [1] -0.0557
```

```r
se.boot
```

```
## [1] 0.8902809
```

Why do our numbers not agree exactly with those of the boot package? This is because our estimates of bias and standard error are just estimates, and they carry with them their own uncertainties. That is one of the reasons we might bother doing jackknife-after-bootstrap.

The 'boot' package has a LOT of functionality. If we have time, we will come back to some of these more complex functions later in the semester as we cover topics like regression and glm.

