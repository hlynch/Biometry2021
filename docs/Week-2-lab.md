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
## 1 3 5 6 7 8 
## 2 1 2 2 2 1
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
## [1] -0.0211
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
## [1] 2.794664
```

```r
UL.boot
```

```
## [1] 6.163136
```

Method #2: Simply take the quantiles of the bootstrap statistics


```r
quantile(xmeans,c(0.025,0.975))
```

```
##  2.5% 97.5% 
##   2.8   6.1
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
##    [1] 5.3 3.2 5.4 4.4 4.5 4.4 5.1 4.0 3.4 5.6 3.9 3.9 5.1 4.9 4.0 4.9 4.1 2.6
##   [19] 5.1 4.7 4.3 3.7 4.4 6.4 4.6 4.3 5.8 5.7 7.4 4.0 5.3 5.0 3.9 3.8 7.2 4.2
##   [37] 3.7 4.1 4.1 6.3 5.2 4.9 4.2 4.1 3.5 3.5 4.4 4.6 4.7 4.0 4.7 3.2 3.6 4.3
##   [55] 5.2 6.4 3.1 5.3 2.8 4.1 5.3 5.2 3.8 4.9 4.8 4.5 4.7 4.3 3.6 4.9 4.3 4.4
##   [73] 4.4 5.3 4.5 4.2 6.5 4.4 2.7 5.5 2.8 5.5 4.8 4.7 5.3 4.8 4.3 5.0 5.5 4.8
##   [91] 4.0 4.1 4.5 4.7 5.0 5.0 5.4 4.3 5.2 3.9 6.2 4.2 5.3 5.1 6.6 3.8 5.4 3.7
##  [109] 3.7 3.8 4.1 4.9 5.5 4.2 5.4 3.0 5.0 5.2 5.1 3.5 3.3 5.5 4.0 3.6 5.8 3.5
##  [127] 4.0 4.6 4.2 4.4 4.1 5.0 5.3 3.8 4.0 2.7 5.4 4.6 5.9 5.0 6.0 6.4 6.0 4.9
##  [145] 4.9 3.5 4.6 5.3 4.7 3.2 4.5 4.7 3.8 5.3 2.9 4.9 4.9 5.1 3.8 4.7 3.8 3.8
##  [163] 5.7 5.0 6.8 5.1 4.6 3.8 5.9 2.2 4.1 3.1 5.0 6.2 4.6 4.9 4.0 4.0 3.1 4.3
##  [181] 4.4 5.9 1.9 3.4 5.2 6.0 4.6 3.5 3.8 3.8 3.4 4.8 5.8 5.0 4.6 3.6 4.8 4.0
##  [199] 5.7 3.9 4.6 4.0 4.6 4.2 6.5 4.4 1.8 3.9 3.6 5.1 5.2 5.8 3.4 4.1 5.8 2.9
##  [217] 3.6 5.3 3.8 3.2 4.6 4.0 5.4 4.8 4.1 6.0 3.5 4.3 3.3 3.3 4.7 4.8 3.0 5.8
##  [235] 4.7 5.0 3.4 4.1 3.6 4.9 4.5 2.4 7.3 7.2 4.5 3.9 3.9 3.4 4.0 5.2 4.9 4.4
##  [253] 4.1 5.2 5.7 5.5 5.7 4.9 4.2 3.9 3.8 4.3 4.2 4.5 4.3 4.8 4.8 4.4 5.0 5.0
##  [271] 5.6 3.8 4.3 3.4 3.8 4.1 5.5 4.5 5.1 6.1 4.6 5.6 4.2 5.5 5.6 4.8 3.4 2.5
##  [289] 4.0 5.0 4.0 3.8 5.0 4.9 5.3 3.7 4.7 5.0 3.6 4.6 3.1 4.4 3.5 4.1 5.3 4.8
##  [307] 4.1 6.7 4.3 6.3 3.6 5.6 3.4 4.2 5.4 4.4 3.6 5.1 6.3 3.6 4.9 5.1 4.4 5.8
##  [325] 4.3 5.4 3.3 3.3 4.4 4.1 3.5 4.5 4.5 4.4 5.9 3.1 4.7 4.9 4.6 4.0 3.8 3.5
##  [343] 5.1 4.7 4.1 2.9 4.0 3.7 4.3 5.2 4.9 4.7 5.6 4.5 5.7 4.5 3.2 3.3 6.0 5.3
##  [361] 4.4 3.0 2.9 6.1 5.2 3.9 4.2 4.8 4.9 5.8 2.8 6.7 1.8 3.6 4.6 3.8 2.8 4.8
##  [379] 2.3 5.3 4.1 5.6 5.0 4.2 3.1 4.6 3.9 4.8 4.7 4.4 4.9 5.8 4.4 5.1 4.5 4.3
##  [397] 5.8 5.0 4.4 3.5 3.1 3.9 5.3 3.6 4.6 4.4 3.6 3.8 3.6 5.4 2.7 5.1 4.7 4.4
##  [415] 3.5 5.0 6.9 3.2 3.9 3.6 3.3 4.1 4.1 4.3 4.3 4.6 4.2 4.7 5.1 4.7 6.8 5.1
##  [433] 4.4 3.3 3.9 5.2 4.0 3.2 5.3 5.1 4.7 5.7 4.2 4.5 5.7 5.3 4.3 3.6 5.0 3.4
##  [451] 4.4 4.3 3.3 3.3 5.4 3.6 3.4 4.6 4.8 3.9 6.3 6.7 2.9 5.0 4.1 4.9 4.4 4.8
##  [469] 4.6 4.4 4.2 4.1 4.6 4.9 2.7 5.1 4.2 3.3 4.7 5.5 4.6 5.1 6.3 4.3 4.5 5.5
##  [487] 5.7 5.0 4.9 2.6 4.3 5.7 5.0 3.2 4.9 3.1 4.3 6.1 3.9 5.2 4.1 4.9 3.3 5.8
##  [505] 2.8 3.3 4.3 5.6 3.6 6.1 3.5 4.9 3.8 4.1 5.7 3.3 4.2 3.8 6.3 4.5 4.8 4.8
##  [523] 4.9 4.0 3.6 4.8 5.2 4.3 4.0 4.0 2.8 3.7 5.2 2.5 5.2 5.5 5.9 3.7 5.3 3.3
##  [541] 4.0 4.9 4.6 5.1 3.9 6.1 5.1 3.9 6.7 5.2 4.0 5.6 4.3 6.0 3.8 5.2 3.9 4.8
##  [559] 5.4 6.0 4.3 5.6 4.5 5.1 5.3 4.9 5.1 4.2 5.0 3.4 4.4 3.1 4.3 5.6 5.1 5.4
##  [577] 5.0 4.8 4.3 5.2 5.7 5.4 4.9 3.9 5.2 2.8 4.9 5.4 3.0 3.8 4.5 6.4 3.3 3.7
##  [595] 4.4 4.0 3.8 5.5 1.9 4.8 4.5 5.2 4.8 3.5 5.2 5.3 5.4 3.9 4.0 5.6 4.6 3.3
##  [613] 5.5 4.5 2.8 4.0 3.4 5.4 4.9 4.2 4.4 3.1 3.8 4.3 5.4 2.8 3.9 3.5 3.8 4.4
##  [631] 4.8 3.7 5.8 6.3 4.4 6.2 3.6 5.0 5.3 5.4 4.3 4.5 3.5 4.0 5.1 4.6 5.8 5.8
##  [649] 4.3 3.2 3.7 5.1 5.1 4.4 4.9 5.7 5.1 5.1 5.5 4.7 5.9 4.4 5.2 4.2 5.1 3.9
##  [667] 2.9 4.7 4.8 4.8 3.5 3.1 4.6 4.8 4.8 5.5 4.2 3.9 3.9 5.1 4.3 4.0 4.1 5.9
##  [685] 3.5 3.9 4.6 4.5 3.6 3.3 6.7 3.9 4.7 5.4 4.7 4.4 4.9 5.9 2.9 6.2 3.8 5.9
##  [703] 4.0 4.0 4.1 4.1 3.7 3.0 4.2 3.4 5.0 6.1 3.4 4.8 3.2 5.5 5.3 4.2 4.8 5.0
##  [721] 5.2 4.8 4.0 4.3 4.8 3.9 4.7 4.0 3.1 4.3 3.9 5.2 4.0 4.2 4.3 5.1 3.3 5.2
##  [739] 3.5 5.6 4.9 4.4 3.9 5.4 4.3 3.5 2.9 3.9 4.4 3.2 5.2 3.0 4.8 4.9 4.4 2.4
##  [757] 5.5 4.9 4.5 1.5 4.2 4.7 4.9 4.2 4.2 5.5 5.7 6.7 5.2 4.7 4.2 4.1 6.1 4.0
##  [775] 5.0 4.8 4.6 4.8 3.7 3.4 5.1 3.9 6.1 4.5 6.6 5.6 4.2 4.0 4.5 4.6 3.8 5.2
##  [793] 4.3 3.3 4.5 3.2 4.4 4.4 4.1 3.0 3.8 4.9 5.2 5.6 5.0 3.8 4.5 5.5 2.7 4.7
##  [811] 4.5 4.5 5.7 3.7 4.7 4.3 4.7 2.9 5.2 4.4 3.9 4.7 5.9 5.6 5.0 3.2 4.6 4.6
##  [829] 4.8 5.1 5.2 4.4 4.5 4.7 4.6 4.2 4.3 3.4 3.8 5.0 3.8 4.6 3.0 3.5 4.7 6.3
##  [847] 6.4 4.4 4.0 5.0 3.4 5.4 2.8 5.0 4.3 5.7 3.1 5.1 5.1 5.9 4.6 4.5 5.4 3.9
##  [865] 3.3 4.7 4.4 5.3 5.8 3.9 4.2 5.2 3.9 4.7 4.8 4.8 3.9 5.2 5.5 5.4 5.5 3.2
##  [883] 3.8 4.1 5.4 3.4 5.4 4.3 4.2 4.7 5.0 6.1 4.4 4.2 5.6 4.4 5.1 4.6 4.4 2.9
##  [901] 5.3 5.6 2.7 4.1 3.7 4.2 3.4 6.7 4.6 3.8 5.1 4.4 3.6 5.6 5.9 3.6 4.3 5.6
##  [919] 3.6 5.3 5.4 4.0 5.3 3.7 5.9 5.8 2.9 5.8 4.9 2.6 5.2 5.5 4.8 5.4 4.9 4.5
##  [937] 4.2 4.7 3.5 2.4 4.1 5.1 4.5 4.9 3.3 2.2 4.6 4.2 6.4 4.9 4.8 3.2 4.6 6.7
##  [955] 4.8 5.5 5.5 5.4 5.7 3.1 3.0 4.2 5.5 3.2 5.1 4.9 4.6 5.9 5.1 5.6 3.9 4.6
##  [973] 3.9 4.0 4.1 4.3 5.1 3.6 4.9 4.3 5.4 5.1 3.8 4.5 4.0 3.7 3.5 3.8 4.9 5.2
##  [991] 4.7 4.4 5.4 4.4 3.7 5.0 5.8 5.3 4.6 5.1
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
## 2.8000 6.3025
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
##    [1] 4.8 4.6 1.5 6.6 4.3 4.5 5.3 4.7 4.0 4.3 4.1 4.8 4.5 5.1 5.1 5.0 4.1 5.1
##   [19] 5.6 4.3 5.7 4.9 5.3 5.4 3.8 4.5 4.9 5.6 4.3 2.8 5.0 4.7 4.8 4.3 4.8 3.9
##   [37] 4.5 3.7 5.9 4.1 3.7 5.1 3.0 4.0 3.0 4.9 5.6 4.7 4.8 5.6 5.1 3.4 4.4 4.2
##   [55] 5.2 4.8 4.0 5.5 4.7 5.1 5.0 4.8 4.2 5.0 4.4 4.1 3.6 4.3 4.2 5.2 4.8 4.9
##   [73] 3.6 4.2 4.4 5.3 3.7 4.9 4.7 4.1 5.0 5.1 5.1 4.0 4.4 5.3 4.4 4.8 3.4 3.3
##   [91] 4.6 5.6 3.9 3.3 3.5 4.3 4.8 4.1 4.6 5.7 4.0 4.7 3.7 4.3 3.8 4.3 5.8 4.3
##  [109] 4.9 4.7 5.1 3.6 4.2 3.1 4.8 4.1 5.1 3.8 4.1 4.6 5.0 4.1 4.0 4.8 3.1 4.4
##  [127] 3.8 3.5 4.8 4.8 3.2 5.1 4.9 2.7 4.7 4.6 4.7 4.7 5.9 4.4 3.9 3.1 4.4 5.3
##  [145] 4.3 4.0 4.3 3.5 3.4 4.5 4.5 3.9 5.1 4.3 3.9 4.5 4.7 3.8 1.8 4.4 4.0 3.7
##  [163] 4.5 4.9 4.6 3.9 3.9 4.2 4.9 6.2 4.8 5.6 4.1 5.5 5.0 4.7 4.1 4.5 3.2 3.2
##  [181] 6.0 5.0 4.5 5.3 4.9 4.3 4.1 3.1 4.3 3.7 4.4 4.3 4.9 4.4 5.4 3.8 4.5 4.7
##  [199] 3.9 5.1 5.8 4.7 3.2 4.8 4.5 4.2 3.6 5.5 5.0 4.2 3.8 5.0 5.1 3.4 6.7 3.2
##  [217] 3.3 4.9 6.7 5.4 6.2 4.5 3.4 4.4 3.9 3.8 3.2 4.0 3.9 4.3 3.9 3.9 4.9 4.6
##  [235] 4.5 4.2 4.8 4.6 5.2 4.2 5.0 6.3 2.7 5.2 4.8 4.9 3.2 4.3 4.4 5.3 4.7 5.6
##  [253] 3.8 2.9 3.4 3.4 4.5 4.3 3.5 3.9 5.9 6.1 5.2 4.1 3.6 2.0 4.8 3.9 4.7 6.8
##  [271] 4.1 5.1 4.8 4.3 4.9 4.8 3.5 4.0 4.9 5.1 3.0 3.9 3.8 4.7 5.2 3.2 5.6 5.0
##  [289] 4.7 4.0 3.8 3.4 4.8 5.0 3.8 4.1 3.8 5.0 3.8 4.1 3.6 4.6 6.7 4.3 6.1 4.9
##  [307] 4.4 4.9 5.8 3.9 4.2 4.2 5.2 4.8 5.1 5.5 4.4 5.2 5.5 2.9 3.2 5.1 4.0 4.7
##  [325] 5.5 3.5 4.4 4.8 3.9 5.0 5.2 6.5 4.3 4.2 4.2 5.6 4.4 5.0 5.2 3.5 4.4 4.3
##  [343] 5.0 5.3 3.1 5.2 4.4 5.9 6.0 3.7 4.2 5.0 5.1 3.8 5.1 3.6 4.1 5.3 5.5 5.2
##  [361] 5.1 4.3 3.5 4.7 3.8 4.9 4.3 3.6 5.1 3.5 3.6 5.9 3.9 3.9 4.2 4.8 5.2 4.6
##  [379] 6.2 3.0 2.4 4.5 4.4 4.5 4.4 4.6 3.8 5.3 3.5 4.7 4.1 4.1 5.1 5.6 5.1 4.6
##  [397] 4.1 3.9 5.3 5.1 5.2 4.6 5.2 3.9 5.3 4.6 4.2 4.0 5.8 4.0 4.1 4.5 3.6 5.1
##  [415] 6.2 6.2 3.6 4.6 3.1 3.8 2.6 5.3 4.3 4.5 4.6 3.4 3.7 4.1 4.1 6.4 5.0 3.2
##  [433] 2.7 3.8 4.3 5.2 5.2 5.0 4.8 4.3 3.8 4.6 4.0 5.0 3.8 3.3 3.8 4.6 4.7 4.4
##  [451] 5.6 4.3 5.5 4.2 4.9 5.5 4.1 3.7 4.9 3.6 4.3 4.0 4.6 3.7 5.0 4.2 2.5 5.1
##  [469] 6.2 4.9 4.0 4.5 4.4 6.7 4.6 5.1 6.6 5.1 4.5 4.4 4.5 3.3 4.2 4.3 3.4 5.6
##  [487] 4.3 3.2 4.2 4.5 4.7 4.4 6.3 3.9 5.5 2.9 4.3 5.6 3.8 4.1 3.9 5.8 3.2 3.4
##  [505] 5.8 5.5 4.5 4.8 4.5 3.1 2.2 2.7 3.5 3.9 6.1 5.7 4.5 4.5 6.0 3.4 4.0 4.9
##  [523] 3.1 4.1 3.5 3.8 6.6 4.3 5.3 4.6 4.2 5.0 2.7 4.3 4.9 4.5 6.0 3.8 4.6 6.8
##  [541] 5.7 4.8 6.0 4.8 5.5 3.5 4.7 3.1 3.6 3.8 3.0 5.3 5.2 3.7 3.9 4.5 4.2 4.6
##  [559] 2.9 3.5 5.5 4.2 4.4 4.5 4.1 4.8 4.7 3.5 4.8 3.0 4.5 5.2 5.3 3.7 5.0 4.0
##  [577] 4.1 4.2 4.5 4.4 2.6 5.2 6.5 4.3 4.7 4.1 5.1 4.3 4.3 2.7 4.1 4.9 3.1 5.2
##  [595] 5.3 6.0 5.5 4.0 4.8 4.1 4.0 4.4 4.4 4.4 2.8 5.1 2.7 3.9 3.9 4.3 3.3 4.8
##  [613] 5.3 6.6 4.1 4.2 4.3 4.4 4.4 4.1 4.8 5.8 4.9 5.5 5.8 4.6 4.0 5.0 4.1 5.0
##  [631] 5.5 3.9 4.0 3.9 4.4 4.0 4.6 4.4 6.0 5.9 4.4 4.3 5.5 3.9 4.0 3.4 5.6 3.9
##  [649] 4.6 5.8 2.1 3.8 6.3 4.0 3.4 3.4 2.0 6.0 4.4 4.5 4.0 2.8 5.5 5.0 4.3 4.9
##  [667] 6.7 4.5 3.7 5.6 5.9 3.6 3.8 4.0 4.9 4.6 5.0 6.0 2.3 5.1 4.3 5.0 3.7 3.6
##  [685] 3.7 4.1 4.7 4.9 4.8 6.2 4.1 5.5 5.9 5.7 4.1 3.8 4.4 5.5 5.2 3.5 3.3 4.1
##  [703] 5.8 3.6 5.3 2.5 4.8 3.5 6.8 3.9 4.7 3.2 4.1 4.1 4.4 5.8 4.7 4.0 7.5 4.9
##  [721] 4.6 4.3 6.0 5.7 4.5 5.3 4.2 2.8 3.7 5.2 4.1 5.3 3.6 4.4 4.6 5.1 4.6 4.6
##  [739] 3.7 5.5 3.6 4.1 5.8 3.3 4.4 4.1 3.0 3.2 6.1 3.8 5.8 3.8 5.0 6.5 4.8 5.9
##  [757] 6.8 5.0 5.7 3.4 4.4 3.1 4.0 6.9 4.6 3.8 4.7 3.5 5.3 4.1 4.6 5.1 4.6 3.7
##  [775] 4.3 4.6 3.8 5.8 4.1 5.4 4.2 3.9 5.2 4.6 5.4 3.7 4.4 2.8 3.5 3.2 3.5 4.1
##  [793] 5.1 3.5 3.2 4.5 3.2 3.8 5.2 3.8 3.7 4.9 6.5 4.2 5.8 4.2 5.4 5.7 4.9 4.6
##  [811] 5.3 5.4 5.4 6.0 3.6 4.1 5.1 4.5 3.8 5.0 5.4 5.8 4.2 5.5 2.7 5.2 5.1 5.8
##  [829] 6.3 4.3 6.5 4.4 4.0 4.5 4.3 4.8 3.5 5.2 4.8 6.3 6.0 3.3 4.6 4.9 3.7 4.0
##  [847] 6.5 4.1 4.1 4.9 6.5 3.4 3.3 3.0 3.9 5.4 2.4 5.8 4.4 4.2 5.7 5.2 3.2 4.9
##  [865] 4.8 4.0 5.1 3.8 2.3 4.5 5.6 3.2 4.0 4.8 4.4 3.5 3.8 5.0 2.9 5.5 3.3 4.4
##  [883] 4.8 3.2 4.3 5.2 5.2 5.5 4.7 5.6 4.1 3.6 4.6 4.2 3.5 4.5 4.1 5.8 4.7 5.0
##  [901] 4.1 5.8 5.3 5.3 4.2 4.3 5.6 4.2 3.5 4.6 3.9 5.0 5.0 3.6 6.0 4.7 5.2 4.1
##  [919] 5.9 5.2 3.5 5.6 3.9 5.8 5.6 4.9 3.7 4.3 1.2 2.4 4.5 5.5 5.9 3.9 6.2 1.8
##  [937] 4.7 4.8 4.6 5.2 5.2 2.9 4.4 5.5 4.9 5.9 5.1 3.4 2.1 4.0 4.3 4.2 3.9 3.2
##  [955] 4.5 5.7 4.5 2.9 4.0 5.4 4.4 4.8 5.1 3.3 5.2 4.6 4.5 5.9 4.9 3.4 3.2 5.6
##  [973] 4.5 3.8 4.7 2.8 3.8 4.1 5.4 3.3 4.0 3.9 5.5 5.7 3.8 2.6 4.6 5.4 5.8 6.4
##  [991] 5.2 3.7 6.0 3.8 4.7 5.2 4.3 4.4 3.9 5.3
## 
## $func.thetastar
## [1] -0.0032
## 
## $jack.boot.val
##  [1]  0.52076023  0.41810089  0.25129683  0.15114943  0.08678679 -0.06289855
##  [7] -0.20845070 -0.23249300 -0.31869688 -0.46685237
## 
## $jack.boot.se
## [1] 0.9288301
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
##    [1] 3.7 5.1 5.3 3.7 3.0 5.8 2.8 3.8 5.5 4.8 4.0 3.4 3.9 4.8 5.2 5.8 4.8 3.7
##   [19] 4.4 4.9 3.0 5.2 5.9 3.6 3.5 5.4 5.4 5.0 5.0 5.1 4.8 4.2 4.3 4.0 6.1 4.5
##   [37] 5.9 4.8 4.0 5.5 4.6 4.3 3.0 5.2 3.9 3.1 5.0 4.2 4.4 5.2 4.6 5.9 4.4 4.3
##   [55] 4.7 5.8 3.9 4.3 4.8 4.7 3.1 5.0 3.2 4.2 5.9 5.0 4.9 5.0 3.7 4.0 6.3 4.0
##   [73] 4.9 3.2 4.0 4.6 4.5 5.3 4.2 2.6 3.7 3.2 3.7 5.0 5.0 5.3 4.3 4.4 4.7 4.3
##   [91] 5.5 4.4 4.9 3.5 4.0 5.1 3.3 4.0 5.5 4.4 5.2 4.9 2.8 4.0 5.3 2.6 4.5 3.7
##  [109] 4.8 4.3 5.3 5.2 4.4 4.4 4.4 4.5 3.4 4.1 4.5 5.3 4.6 5.1 4.3 3.9 3.6 4.1
##  [127] 4.6 3.8 3.6 6.3 5.0 5.0 3.7 3.0 4.4 3.7 7.1 3.8 3.6 5.7 5.0 5.8 4.8 5.3
##  [145] 4.9 3.4 4.3 5.4 6.0 4.9 4.7 5.6 4.6 4.7 3.9 4.5 4.7 4.6 5.2 5.2 4.5 3.8
##  [163] 4.0 3.2 5.3 4.1 5.0 3.7 6.1 4.7 4.1 5.2 3.9 4.7 5.0 4.8 4.1 4.6 4.4 3.4
##  [181] 4.9 5.3 5.0 4.7 4.2 5.3 3.5 3.1 6.2 3.5 5.0 4.2 4.4 3.8 5.1 4.1 6.3 4.5
##  [199] 6.2 4.4 5.1 4.9 3.3 4.5 3.8 4.8 5.4 3.1 4.2 6.1 4.2 3.4 4.1 5.3 6.5 3.7
##  [217] 4.4 5.3 4.0 3.5 4.4 5.4 4.3 5.1 4.8 3.6 4.2 6.7 4.1 5.5 4.7 4.8 6.1 6.9
##  [235] 3.9 4.7 5.3 2.7 5.2 4.1 4.2 3.6 5.2 3.8 6.3 3.6 4.8 4.1 3.7 4.3 4.4 4.0
##  [253] 3.9 3.8 3.5 3.5 3.8 3.5 4.7 4.4 4.3 4.0 4.3 4.5 4.8 4.5 3.5 3.7 2.8 5.3
##  [271] 4.1 3.9 3.8 6.3 4.5 3.0 5.6 4.4 3.7 4.5 2.5 4.6 5.0 4.7 2.8 5.3 6.0 5.1
##  [289] 6.9 4.9 4.6 4.3 4.0 5.0 5.1 3.4 4.3 4.8 3.6 4.9 3.7 4.8 3.9 4.0 4.8 6.3
##  [307] 4.1 4.9 4.3 4.6 4.4 3.7 5.0 4.8 4.2 5.8 5.1 4.1 4.4 4.4 6.1 4.3 5.3 5.7
##  [325] 5.6 4.3 4.2 1.8 4.9 4.0 3.9 4.6 6.4 4.2 3.4 4.8 4.5 3.5 5.4 3.6 4.1 5.1
##  [343] 5.2 4.8 5.6 4.3 3.8 2.6 5.3 3.5 4.0 4.7 5.8 4.9 4.1 3.4 4.7 4.9 5.5 3.7
##  [361] 4.6 4.7 3.6 3.4 5.1 6.0 4.2 3.4 3.9 4.6 3.6 4.2 5.0 6.3 4.8 3.8 4.7 3.6
##  [379] 5.1 3.8 5.6 3.7 5.5 4.0 6.1 5.1 4.1 3.6 5.2 4.5 3.4 4.5 4.4 4.0 2.9 4.3
##  [397] 2.2 3.9 4.5 4.0 6.1 5.1 6.7 4.7 3.4 2.7 5.4 4.4 3.1 4.4 4.0 4.0 6.3 5.4
##  [415] 4.6 4.3 4.8 5.5 6.3 5.3 5.0 5.4 4.4 5.3 5.8 3.9 3.5 4.4 5.2 4.6 5.6 2.8
##  [433] 4.2 4.3 4.0 4.5 4.7 3.5 4.3 4.1 3.3 4.6 3.1 5.1 4.8 4.6 5.7 2.7 3.8 4.8
##  [451] 4.4 5.2 5.1 6.0 5.3 3.7 4.6 4.9 3.9 5.7 3.9 5.0 4.6 5.8 4.7 5.1 5.8 5.5
##  [469] 5.9 4.6 4.9 5.7 4.9 3.1 5.4 5.3 4.7 5.2 2.6 3.6 2.8 5.5 5.0 4.5 5.4 4.3
##  [487] 6.2 4.0 2.8 4.0 6.0 4.0 4.1 5.8 4.5 5.2 5.6 6.0 4.8 2.8 3.3 3.8 4.9 5.5
##  [505] 6.5 4.8 4.5 3.5 3.4 4.1 3.9 3.8 5.0 5.3 4.7 5.9 4.1 3.8 5.4 4.8 3.6 5.5
##  [523] 6.7 3.8 6.0 5.9 3.5 4.3 4.7 3.8 4.5 3.3 4.7 4.1 3.7 4.6 5.3 5.5 5.0 5.4
##  [541] 5.1 4.6 5.0 4.9 4.8 5.5 4.6 4.6 5.8 5.0 5.2 4.4 5.8 3.8 4.9 4.6 4.1 4.8
##  [559] 3.7 4.6 4.6 4.2 4.2 3.9 2.7 4.2 4.4 5.2 4.4 3.6 5.1 5.3 2.1 5.2 4.0 4.2
##  [577] 3.6 2.8 6.0 2.3 3.7 3.6 2.7 3.6 4.1 4.6 4.1 3.8 3.6 6.8 4.8 3.8 3.9 4.9
##  [595] 2.5 3.8 5.5 4.7 4.2 5.1 5.2 4.7 4.8 5.4 5.6 5.2 5.2 5.5 4.3 3.4 5.7 4.0
##  [613] 3.1 5.2 6.4 5.9 4.1 4.6 5.3 5.6 4.6 2.7 5.5 4.6 5.3 4.5 4.6 4.3 4.3 4.3
##  [631] 3.7 4.8 5.2 5.4 6.0 4.8 5.6 6.5 4.5 6.8 3.3 4.3 4.8 4.0 4.6 4.1 4.0 4.3
##  [649] 3.5 5.8 3.8 3.7 3.9 5.3 7.1 5.6 3.7 4.5 4.7 4.9 5.0 6.3 5.5 3.3 4.4 5.2
##  [667] 3.3 5.8 4.1 5.5 3.0 4.8 4.2 2.5 5.6 3.6 4.7 3.4 3.5 4.3 5.0 4.7 3.9 4.9
##  [685] 5.3 5.1 4.4 4.3 3.8 4.9 4.7 5.4 3.5 4.0 4.9 4.4 4.8 4.2 5.7 2.4 3.4 2.8
##  [703] 4.2 5.0 5.0 4.9 3.2 5.0 4.8 2.6 3.6 6.0 5.5 4.1 4.2 4.1 4.1 5.7 5.1 3.9
##  [721] 3.8 4.9 5.0 4.4 3.3 3.7 4.2 3.7 4.1 4.6 4.3 5.1 4.2 4.7 4.9 4.1 6.2 4.0
##  [739] 3.7 3.1 5.3 4.2 5.9 4.5 4.4 4.2 4.9 4.3 5.2 4.9 5.6 4.9 2.5 4.5 6.2 4.2
##  [757] 5.7 3.7 3.2 5.0 5.2 5.9 5.1 4.5 3.8 6.3 4.4 4.3 4.6 4.2 5.5 3.5 4.9 3.6
##  [775] 5.3 2.3 4.9 4.5 5.5 5.5 3.1 4.3 4.5 4.5 4.9 3.8 4.4 4.3 4.0 5.7 4.5 4.0
##  [793] 4.9 5.0 4.9 3.6 3.8 5.9 4.0 5.5 3.0 5.4 4.1 5.9 5.3 4.3 5.1 5.4 4.6 3.4
##  [811] 5.6 4.0 4.7 4.7 4.7 5.1 4.4 5.6 3.2 2.6 3.3 4.2 4.8 3.0 4.6 4.5 4.1 4.7
##  [829] 4.9 4.4 5.9 6.1 4.0 4.7 4.8 3.5 3.7 3.5 4.5 5.6 4.3 5.2 3.7 5.8 4.9 3.4
##  [847] 3.9 4.8 6.3 6.4 5.3 5.7 3.7 5.0 3.5 4.6 3.8 3.6 6.0 7.0 3.3 3.5 4.4 4.2
##  [865] 4.4 3.3 4.0 4.8 4.5 6.7 4.6 4.5 5.5 5.8 4.8 3.7 4.9 3.2 4.9 2.1 2.9 3.8
##  [883] 5.4 4.4 4.5 3.8 6.3 6.2 6.1 5.5 5.0 5.3 4.3 3.3 3.6 4.3 4.9 3.8 2.8 5.2
##  [901] 5.9 4.4 5.2 3.8 4.8 4.4 5.5 5.3 5.0 3.9 4.7 4.9 3.2 7.0 4.9 5.8 4.4 4.1
##  [919] 3.3 6.3 3.3 4.9 3.9 3.6 4.4 3.6 4.7 5.0 6.3 3.7 5.9 4.8 6.2 5.3 4.5 3.8
##  [937] 5.2 5.2 4.0 3.5 4.7 4.4 3.5 6.2 6.4 4.6 3.8 5.6 4.4 5.7 4.8 3.7 4.8 5.7
##  [955] 4.5 6.0 3.3 5.3 4.2 5.4 3.9 4.3 4.4 5.1 5.3 5.5 3.1 4.8 5.7 5.0 4.6 2.1
##  [973] 4.0 4.0 6.1 3.7 4.1 5.7 3.7 4.7 5.1 3.8 4.6 3.3 4.8 3.5 3.4 4.7 4.2 3.8
##  [991] 4.4 6.3 6.0 4.7 4.0 4.1 2.6 3.1 4.4 4.6
## 
## $func.thetastar
## 72% 
##   5 
## 
## $jack.boot.val
##  [1] 5.500 5.500 5.300 5.300 5.156 5.100 4.900 4.800 4.700 4.500
## 
## $jack.boot.se
## [1] 0.9708059
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
## [1] -0.1101931
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
##    9.422240   15.567774 
##  ( 4.141297) ( 7.027849)
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
## [1] -0.15406970 -0.04331664  1.22270643 -0.18529799  0.45016743  1.50530603
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
##    [1] -0.2871561632  0.2292945339  0.6051068640 -0.6446550476 -0.0245023234
##    [6]  0.4320705878 -1.2517619746 -0.5455189735 -0.3428517087 -0.3220903496
##   [11] -0.0267159331 -0.5006705054 -0.6705817809  0.2796478161  0.9944420280
##   [16]  1.0341457505 -0.3920192660 -0.6223897990  0.4156227429 -0.3362137827
##   [21] -0.7744830174 -1.5900607532 -0.6264599089  0.6770061053  0.1450372958
##   [26] -0.3922552102 -0.0651728842  0.2452212442 -0.1028398013  0.3479350130
##   [31]  0.5434034695  0.0808795301  0.1055365844  0.0750507106 -0.3079605038
##   [36] -0.0800548073  0.0155268730 -0.2238991495 -0.0563518573 -0.4343046455
##   [41] -0.7084117743 -0.3040453023  0.1273704827 -0.0477535137 -0.1095433959
##   [46]  0.0661963840 -0.0281395083 -0.6355870378 -0.5385799120 -0.0778742230
##   [51] -0.2180274941 -0.4821472541 -0.0518874552  0.2143715772 -0.8540166958
##   [56] -0.5038081905 -0.0915659903 -0.7579048365 -0.0055877977  0.1993017028
##   [61] -0.8345610233  0.2442815359 -0.2376167437  0.0642589558  0.3172692840
##   [66] -0.0052788833  0.1484249370  0.1333735935 -0.5937348658 -0.1092840857
##   [71]  0.1190458900 -0.2271827943 -0.3529658990  0.4229904389 -0.0259943634
##   [76] -0.2319617731 -0.2603801795  0.1292347673  0.0569008112 -0.5093485567
##   [81] -0.6523918041 -0.2792251005 -0.0240669776 -0.8806294161  0.3820373709
##   [86]  0.4269563181 -0.0377034574  0.3638211729 -0.3492815599 -0.1140786021
##   [91] -0.2094938692 -0.7622659515 -0.3831037276  0.2604691535 -0.7512889097
##   [96]  0.4019599040 -0.2942121029 -0.1091381008  0.1935174224 -0.5265600516
##  [101]  0.1130220496  0.3993615791 -0.2988455859  0.1127286321 -0.4296527192
##  [106] -0.3177107037 -0.4892358948  0.3326234214 -0.5777117647 -0.2905882938
##  [111] -0.8272065164 -0.3785588138  0.1672837815  0.0265455975  0.1506283746
##  [116] -0.1820876535 -0.0445669913 -0.2026830953  0.3239343776 -0.0626678773
##  [121] -0.9986815387 -0.3548764397  0.4972391488 -0.0209921697 -0.5848647773
##  [126]  0.1167670719 -0.1924378496  0.0963578324 -0.3736742964  0.3833454412
##  [131] -0.1183873818 -0.4743113647  0.4648109284 -0.2768651621  0.0245081894
##  [136] -0.0319213428  0.8636385976 -0.3614471486  0.3279298801 -0.1492437885
##  [141] -0.0272371841  0.1139099747  0.1215596038 -0.0757648766  0.0581496902
##  [146]  0.2297759449 -0.2025387872  0.6738030054  0.3495123601  0.0540154100
##  [151] -0.0334364626  0.2042905286 -0.4433726253  0.3115232903 -0.6271648261
##  [156]  0.5464872259  0.1897887182  0.0065302985 -0.5676973063  0.5800174344
##  [161] -0.6442228321  0.0164937169 -0.2800597035  0.0708808547 -0.3459435202
##  [166] -0.3158169082 -0.1371600848  0.0801246862  0.5209583449 -0.8596813745
##  [171]  0.7888053594 -0.3717001353  0.2033594589 -0.0519506874  0.6835804043
##  [176]  0.2653872687 -0.2080390029  0.2219651215 -0.1718162702  0.0875003010
##  [181]  0.1478142567  0.3338188407 -0.1396710056  0.2933346860  0.6570677396
##  [186] -0.0215754596  0.2118057322 -0.4890486094  0.5489172076  0.2136037775
##  [191] -0.9997074148 -0.0334498737 -0.2007719868  0.2217073707 -0.6004109599
##  [196] -0.3702023974  0.1750293514  0.3419111574 -0.1154620624  0.3159732707
##  [201]  0.3519014243 -0.4921491202  0.0316635267  0.1536573749 -0.4117204696
##  [206]  0.0851752120 -0.1135669612  0.5723205792 -0.5338765206  0.3731378999
##  [211]  0.4219413278  0.2695341794 -0.2658691390  0.5501580072  0.0535114286
##  [216] -0.2226125250  0.1773839081 -0.5093142473  0.0683676369 -0.2560381241
##  [221] -0.4321599016  0.1489149413  0.2301628905  0.2636324342  0.7588200999
##  [226] -0.1786058906 -0.0730814320 -0.5800198406 -0.5167598239 -0.1075989190
##  [231]  0.1328221922 -0.0231787994  0.2013270335 -0.6912520269 -1.3302579614
##  [236] -0.1088388738 -0.1770285218 -0.4347060253 -0.0482217112 -0.1698419507
##  [241]  0.2082616262  0.0706241655 -0.0002684380 -0.3856691238  0.6428592881
##  [246] -0.3534857840  0.3780581050  0.2049275809 -0.0616416459 -0.4515679006
##  [251] -0.9830207162 -0.1657235706 -0.1152906901  0.2554865711 -0.0978062861
##  [256] -0.2512686276 -0.0674804928 -0.1950247983 -0.0283505125  0.7451860905
##  [261] -0.5900172036 -0.1410190921 -0.1267709076  0.2027278649  0.1691441133
##  [266] -0.7356835247  0.1348638714  0.1640977574  0.0989208634 -0.0966118441
##  [271]  0.3272302959  0.1610421379  0.1385388824 -0.0107947838 -0.0715270930
##  [276] -0.3977147648 -0.9107218690  0.2769137242 -0.0745205617  0.5768349314
##  [281] -0.4210433058  0.0430444213  0.3249358956  0.1148120246  0.1659236045
##  [286] -0.8073081335 -0.1165838931 -0.3417919224  0.0823618101  0.4349028944
##  [291] -0.2478071181 -0.1932986678  0.0634973282  0.7489379920 -0.2018193768
##  [296] -0.2516964248  0.1915449761  0.4644584473  0.6938943438  0.0447844546
##  [301]  0.1534462034  0.4316805413 -0.0205434869 -0.3465531492 -0.8648811772
##  [306] -0.3177879531  0.2130255928  0.3387844603  0.6014736395 -0.6080638663
##  [311]  0.4466021496 -0.2094936067 -0.4544887957 -0.1632112699 -0.0131975959
##  [316] -0.1138905155 -0.1009467461  0.2599980637 -0.9934602409  0.2041555587
##  [321]  0.4176809727 -0.2130363239  0.6069275871  0.0850480808 -0.2013844489
##  [326]  0.2441808088  0.2389600975 -0.2610258745 -0.6033890379 -0.2586957058
##  [331]  0.5414919456 -0.0255275319  0.1485084776  0.0890676983  0.3082601527
##  [336] -0.5601815018 -0.0372987418 -0.1108704801 -0.4001111115 -0.3157722102
##  [341] -0.3573760127 -0.1431765240 -0.2966935921 -0.2986096020  0.4943238910
##  [346]  0.1235245967 -0.0280540802  0.1526269405  0.1375140572 -0.1901375392
##  [351]  0.7551011813 -0.2303109650  0.2769721538 -0.0656120746 -0.1700780530
##  [356] -0.6351710937 -0.7015131554  0.8701032741 -0.5240644547 -0.2351328262
##  [361] -0.1457122136 -0.0847705854 -0.4249988932  0.1569207442  0.0060040681
##  [366] -0.3536836974  0.5674212836  0.4058819843  0.3724768085 -0.5758538364
##  [371] -0.5597185315 -0.2252717731 -0.3610359354  1.2892127229  0.0504256238
##  [376]  0.3479350130  0.2789112026 -0.1279826884 -0.0618746836  0.0574798410
##  [381] -1.1099712501 -0.5392896449 -0.1390244702 -0.0265384348 -0.4940795325
##  [386] -0.0857373807  0.0071032050 -0.1072210067  0.4245577580 -0.5526663436
##  [391]  0.2297839394 -0.3529658990  0.0218941224  0.4824972717 -0.0965956417
##  [396]  0.0297465677 -0.0811155176  0.5256789256 -0.0769514637 -0.4205031321
##  [401]  0.1021064258  0.4041910959 -0.4309421412  0.0844464077  0.0967056737
##  [406] -0.0419809888  0.5740313100 -0.4433492024 -0.6705817809 -0.1546478462
##  [411]  0.2605130526  0.4730571857 -0.3506608035  0.4290824829  0.3773337242
##  [416] -0.3943979005  0.2782270097  0.0287226254  0.1277527909 -0.0395434440
##  [421] -0.6510385514 -0.0986947779  0.5804874300  0.3198746093 -0.2792981982
##  [426]  0.5738036830  0.2252409450  0.3700579694  0.2980254853 -0.0206871126
##  [431]  0.5703274850 -0.6210347017  0.4174019162  0.1650272610 -0.5647081357
##  [436]  0.1610469653 -0.1557415945 -0.0547376727 -0.8121524817 -0.4094566820
##  [441] -0.1737884553 -0.5063761189  0.0045115368 -0.0945864096 -0.3761607108
##  [446] -0.1070017865  0.4003726851 -0.2043848136  0.1254221077 -0.7775357208
##  [451]  0.1358649841  0.6807520707  0.5787538733 -0.8131594447  0.2049107233
##  [456] -0.4191895137  0.1005749601 -0.0705150273 -0.8432548832  0.2848950891
##  [461]  0.7592016349 -0.0264634931 -0.7504475125 -0.5397578601  0.4277434973
##  [466]  0.1772185499 -0.5888077912 -0.5028179993 -0.0605991821 -0.6850819982
##  [471] -0.7071476117 -0.2902911329  0.1945945236 -0.7423254320  0.5212219749
##  [476] -0.3162323710 -0.3272537527 -0.1389562856 -0.4746317801  0.4122105273
##  [481]  0.1573455451 -0.2261167328 -0.5529202087 -0.3732991536 -0.1982849121
##  [486]  0.0618541930 -0.5040748474  0.0112838854  0.2764144522  0.5565954067
##  [491] -0.4529188142  0.1926395437 -1.3388223101 -0.4627294155 -0.2994870793
##  [496] -0.5313721584 -0.1721181598 -0.5204671862 -0.1856731626 -0.2469373745
##  [501] -1.2817446665 -0.1086808424 -1.0712351766 -0.4973173201 -0.3647267354
##  [506] -0.3488504366 -0.2197674401 -0.6301257242 -1.0598521026 -0.4613522859
##  [511]  0.6750003501 -0.1332783603 -0.0586995097  0.5744111069  0.2091819231
##  [516] -0.1530615771 -0.2985758749  1.3316749105  0.2960834256  0.3789612403
##  [521] -0.2078419336  0.1556512875 -0.3485632583  0.0770589509 -0.2841847290
##  [526] -0.1104782568  0.1590137524  0.4992390998  0.5427047139 -0.5232805274
##  [531] -0.3907649164 -0.2628511169  0.0110654430 -0.4710594529  0.5397926592
##  [536] -0.0009140441  0.1317425593  0.1156128655  0.6069924349 -0.9072452520
##  [541] -0.3912599115 -0.0680350931  0.2307406663  0.6001680689  0.0892875181
##  [546]  0.7016532581 -0.3490346031  0.6415916716 -0.4285191010 -0.4086008830
##  [551]  1.0842551726 -0.2305196211  0.4906845927 -0.4395588899 -0.0248357350
##  [556]  0.2599201094  0.7304823613  0.0843413128  0.3925571935  0.3419111574
##  [561]  0.3168747348 -0.5004866260  0.4507826632  0.8842389059  1.0070295632
##  [566]  0.5305536648  0.4290364071  0.6054351701 -0.7912806871  0.5808647218
##  [571] -0.7562587327 -0.5923417356 -0.0199629859  0.3060703900 -0.2701607182
##  [576] -0.0733177245 -0.8977487704  0.2178434537  0.0640966398  0.6004201994
##  [581]  0.3051643162 -0.0467421526 -0.7169005338 -0.0735711199 -0.6213570030
##  [586] -0.0624201366  0.3524768178 -0.3660480778 -0.3836192078 -1.1335270305
##  [591] -1.3792933777 -1.0858047101  0.0249207744  0.0940327003 -0.5038081905
##  [596] -0.1445137313  0.5453007163 -0.2540943102 -1.3637032865 -0.1327796977
##  [601] -0.8546148837 -0.4433950742 -0.3340113254 -0.1217994889 -0.4467683121
##  [606]  0.0033141551  0.9869264717  0.9446821576  0.6511168608  0.3338188407
##  [611]  0.5353543631 -0.6973887545  0.1205282605 -0.5352367222  0.0390395031
##  [616] -0.7153446745 -0.0690037539  0.2422830361  0.1760798822 -0.7259486802
##  [621] -0.1539020242  0.7999304769 -0.1077557158 -0.0799162288 -0.2668055531
##  [626] -0.0947282562  0.0997941432 -0.5433172848 -0.1704580008  0.0043355115
##  [631] -0.0840377701 -0.0227804854 -0.1487265276  0.1129072694 -0.2088087819
##  [636] -0.3087415029 -0.3367356962  0.1667320740  0.0095012251  0.0372845905
##  [641]  1.3695425093  0.1654826115 -0.0226093537  0.2140920612  0.2379396264
##  [646] -0.4646099507  0.2219651215  0.1305850335  0.1637816227 -0.1054665005
##  [651]  1.0531965606  0.0337406391  0.1965264078 -0.0444885764 -0.5595461236
##  [656] -0.0684230675 -0.2928432403 -0.0370929412  0.1195153636 -0.1130247821
##  [661] -0.4784189988  0.2546713171 -0.8645175199 -0.0412377279  0.7766247845
##  [666]  0.8424368077 -0.2194196576 -0.1330228020 -0.8011416122  1.2264529317
##  [671] -0.3923327094  0.0270745935 -0.1337824150  0.0249724734 -0.6101969744
##  [676] -0.0061322803  0.5735785502 -0.1438862682 -0.9223145757 -0.0716488586
##  [681]  0.1820737389  0.2206916750  0.7052278968 -0.2951970611 -0.3789072288
##  [686] -0.1551726314  0.2136355275 -0.2897755505 -0.2283224807  0.4307434088
##  [691] -0.0189526451 -0.0271546294 -0.5782873783  0.1817718922 -0.3157959854
##  [696] -0.9085104030 -0.2817048024  0.3401200724 -0.6440762527 -0.0152430406
##  [701] -0.1105692449  0.0469613515  0.2212904675 -0.0795643380 -0.3593641268
##  [706]  0.8600497225 -0.2004998291 -0.3952074939 -0.2107037313  0.1196044055
##  [711] -0.1337712161  0.3546988774 -0.7973812527  0.2466480836 -0.7155756061
##  [716] -0.3156063977 -0.3114313761 -0.0618854142 -0.2005747501  0.0598603180
##  [721] -0.2548457691  0.3013846320  0.5837798320 -0.0138670719 -0.5712170311
##  [726]  0.4453400126 -0.8132538995 -0.1039587200 -0.2595144467 -0.0252202516
##  [731] -0.2183085421  0.4368905892 -0.3490346031 -0.0624201366  0.6685688181
##  [736]  0.2522084724 -0.3030923668  0.5341038930 -0.0109683045 -0.0309557116
##  [741] -0.5027457381 -0.6038337125 -0.1539020242 -0.7331533493 -0.2709408582
##  [746]  0.1043685282 -0.3295339036 -0.5211018811  0.0884766416 -0.0917517559
##  [751] -0.2361474337 -0.1002267654 -0.0703488436 -0.7004570399 -0.8970156721
##  [756]  0.2215999744 -0.1140573016  0.1700261981 -0.8011416122 -0.5352367222
##  [761]  0.4277582919  0.0252992449  0.3280676503 -0.0196795847 -0.3748361496
##  [766] -0.2294593630 -0.3711939118 -0.3214790614 -0.0286924689 -0.4640798932
##  [771] -0.4957884071  0.1377951334 -0.0774465277 -0.0635501781 -1.0252102744
##  [776]  0.0810386791 -0.4876575918 -0.5965564464  0.1304995662  0.1010192139
##  [781] -0.6635362812  0.2204175342  0.1597776905  0.5409539929  0.1315575221
##  [786]  0.2068104146 -0.3689575933 -0.6507323079 -0.1910374410  0.2193677057
##  [791]  0.3343074676  0.1382322515 -0.1162871222  0.1386188145 -0.3524707910
##  [796] -1.1681767612 -0.3698245202 -0.1488729634 -0.1830120246 -0.1351595847
##  [801]  0.1835788272  0.1675890442 -0.3766396318  0.2628350665 -0.6275299350
##  [806] -0.8574818562  0.6643503304  0.6035877387 -1.2904718560 -0.0610793127
##  [811] -0.3211695203 -0.3611300539 -0.2710315182  0.1710190869 -0.5045362263
##  [816] -0.3781272966 -0.0002775544 -0.2692484486  0.2444427125  0.2300415232
##  [821] -0.0890368598  0.1989788799  0.0482483825  0.4006185025 -0.9946066363
##  [826] -0.4355007773 -0.2749192714 -0.5152405405 -0.5617366636 -0.1584712308
##  [831] -0.4917005746 -0.1186638727 -0.4627622359 -0.0063495364 -0.5340619694
##  [836] -0.1704488814 -0.1379179887 -0.0247733866  0.0456285729 -0.3716711202
##  [841]  0.4333264650  0.3767114676 -0.2449970464 -0.5545926857 -0.3697676569
##  [846] -0.3776572477 -0.4106865038 -0.0128588697  0.1600789949 -0.1300999773
##  [851] -0.3495858150  0.1492379244  0.0706299323 -0.3507458130  0.5118432421
##  [856]  0.5260491942 -0.1650122553  0.4480232898  0.2072374270 -0.6003320203
##  [861] -1.1582320039  0.6495287719  0.4728544499 -0.4674257348  0.7094331227
##  [866]  0.2821715064 -0.1005214495  0.1408005971  0.1749766270  0.6259279779
##  [871] -0.1211545978  0.3379044487 -0.6239560524 -0.3477933321 -0.5772526801
##  [876] -0.1038637179  0.1285343345 -0.2412994452 -0.9260342484 -0.5456433260
##  [881] -0.2748551575 -0.7625807140 -0.6147817317  0.3708118839 -0.1560381971
##  [886]  0.0824635837 -0.1036837453  0.2088544598 -0.4926547140 -0.0520014864
##  [891]  0.5530437340  0.2406849215 -1.0028651231 -0.0912721951  0.3216674999
##  [896] -0.9043617892 -0.1288592062  0.2727187652  0.1332685242  0.2505826587
##  [901] -0.4236474759 -0.7092061498 -0.1646768085 -0.0975993305  0.2221091619
##  [906] -0.8396866797  0.0811352988 -0.0106008923 -0.3781272966  0.1077355446
##  [911]  0.2104556553  0.0099101112 -1.4636589862  0.5705262340 -0.0687238643
##  [916]  0.0489976989 -0.6610707599 -0.9155091391  0.6108967827 -0.0525660742
##  [921] -0.2037235275  0.2391656799 -0.0729668880  0.4466021496 -0.7813026921
##  [926] -0.6020209221 -0.3230799447  0.0505836691  0.1704425314 -0.5203499783
##  [931]  0.3882324903 -0.7765162727 -0.3610841046 -0.6909882239  0.0642556938
##  [936]  0.1467374459  0.2041568545 -0.2373445476 -0.1717790859 -0.4004669626
##  [941] -0.3411964003 -0.4882306627 -0.5139102595 -0.0577004158  0.0509347975
##  [946] -0.2192705017 -0.2931604173 -0.0857373807 -0.1840738778 -0.6020252598
##  [951]  0.0515068241 -0.2841474567 -0.1382374161 -0.3537829726 -0.0799501885
##  [956]  0.3133369262  0.4526463304  0.0390635074 -0.0537645812 -0.0866974381
##  [961]  0.2924867579  0.1069818407  0.4866186794  0.2832205010  0.4962884159
##  [966] -0.2646325300  0.2219126145 -0.4122780958 -0.2784724800  0.0684209369
##  [971] -0.1932525549  0.3835044654 -0.3863308356  0.3709616448 -0.9726558753
##  [976]  0.7412566459 -0.2709304169 -0.3918591140 -0.4669050195 -0.0100945163
##  [981] -0.8560328302 -0.3403873388 -0.2775148322  0.3847248562  0.2888139734
##  [986]  0.7185929039  0.9793046240  0.2032481968  0.0922944879 -0.6930316330
##  [991]  0.5841491137 -0.2758643322 -0.4965739308  0.4919965808  0.0252097447
##  [996] -0.4627148599 -0.5894171243 -0.2757624124  0.0585960778 -0.1919592215
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
##   0.60524318   0.18648772 
##  (0.05897259) (0.04169618)
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
## [1] -0.4886037  0.8702811 -0.7735795  0.9407835 -0.0838929 -1.5086505
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
## [1] -0.012
```

```r
sd(results2$thetastar)  #the standard deviation of the theta stars is the SE of the statistic (in this case, the mean)
```

```
## [1] 0.9012259
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
## t1*      4.5 0.007507508   0.9025757
```

One of the main advantages to the 'boot' package over the 'bootstrap' package is the nicer formatting of the output.

Going back to our original code, lets see how we could reproduce all of these numbers:


```r
table(sample(x,size=length(x),replace=T))
```

```
## 
## 0 1 3 5 6 8 9 
## 1 1 2 1 3 1 1
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
## [1] 0.0726
```

```r
se.boot
```

```
## [1] 0.8963887
```

Why do our numbers not agree exactly with those of the boot package? This is because our estimates of bias and standard error are just estimates, and they carry with them their own uncertainties. That is one of the reasons we might bother doing jackknife-after-bootstrap.

The 'boot' package has a LOT of functionality. If we have time, we will come back to some of these more complex functions later in the semester as we cover topics like regression and glm.

