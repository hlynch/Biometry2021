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
## 0 3 4 5 7 8 
## 1 1 2 2 3 1
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
## [1] 0.0081
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
## [1] 2.75114
```

```r
UL.boot
```

```
## [1] 6.26506
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
##    [1] 4.8 5.6 3.9 2.7 4.8 2.6 4.2 5.5 4.6 3.8 5.1 4.0 4.5 3.8 3.7 4.0 3.4 4.6
##   [19] 4.9 3.7 3.3 5.3 3.6 3.5 4.0 3.2 4.4 6.1 4.6 4.4 5.1 4.7 2.2 4.4 4.0 4.6
##   [37] 4.7 5.6 5.3 5.2 6.0 6.2 4.4 5.4 5.3 4.6 3.7 4.0 5.3 3.2 6.7 4.1 3.5 4.5
##   [55] 3.7 4.4 3.5 3.4 4.1 2.9 3.8 3.7 4.1 3.0 4.1 4.0 3.5 3.6 5.0 5.2 3.5 5.9
##   [73] 1.8 5.1 4.8 2.9 5.5 3.6 4.6 3.9 3.1 2.9 5.2 5.1 4.1 4.5 3.6 3.4 5.4 5.3
##   [91] 4.3 4.7 5.7 4.1 4.1 2.1 5.4 4.5 4.6 4.5 4.4 3.5 5.4 4.1 3.7 3.6 4.6 3.1
##  [109] 3.6 4.4 3.6 4.4 3.2 3.5 4.2 5.3 1.9 3.5 5.1 4.6 4.0 4.3 5.9 4.3 6.0 4.7
##  [127] 5.1 4.1 5.3 6.6 4.8 5.7 4.1 5.1 2.9 4.2 3.5 4.5 5.2 5.3 5.2 4.2 3.0 4.0
##  [145] 3.3 4.4 4.3 5.4 4.9 2.6 6.7 5.1 4.3 2.8 3.7 3.9 4.0 3.4 3.3 3.3 5.5 3.9
##  [163] 3.1 5.9 5.2 4.3 5.2 4.7 2.4 5.2 3.8 5.0 6.2 4.3 4.0 6.0 5.1 3.7 4.0 4.9
##  [181] 4.5 4.4 3.1 4.3 5.7 3.3 3.8 4.0 4.5 5.0 3.9 4.5 5.3 4.5 5.1 3.8 4.5 4.2
##  [199] 3.9 4.9 3.8 3.5 5.6 5.2 4.0 5.5 5.2 5.0 5.3 2.9 5.9 5.2 5.0 4.8 3.5 5.0
##  [217] 3.8 6.1 3.1 3.8 3.0 6.2 5.1 4.3 5.7 3.3 3.7 3.3 6.1 3.7 3.5 5.6 5.2 6.3
##  [235] 5.2 4.4 4.8 3.9 4.2 4.9 3.4 6.1 5.4 5.0 4.0 5.8 3.8 5.7 5.1 5.4 4.6 4.5
##  [253] 4.0 3.8 4.9 4.8 5.8 3.6 3.7 4.9 5.0 4.5 4.6 4.5 5.4 5.4 5.8 4.7 3.2 4.0
##  [271] 5.0 5.6 4.6 3.9 4.1 4.6 3.7 2.9 4.4 4.8 4.8 3.5 5.7 4.0 3.5 4.3 3.0 3.6
##  [289] 3.6 5.1 4.7 4.7 3.5 5.2 4.2 5.3 4.8 4.1 3.6 3.7 4.8 4.7 3.2 4.4 4.4 4.0
##  [307] 4.6 4.1 5.6 4.8 6.9 5.1 6.1 4.2 4.1 4.8 3.6 4.4 2.9 4.6 4.2 5.2 4.9 3.5
##  [325] 4.3 3.8 3.6 5.4 4.7 5.2 3.6 5.0 3.2 5.0 5.3 3.0 5.9 5.4 5.1 6.5 4.8 5.1
##  [343] 5.6 5.5 6.2 2.0 3.5 5.4 5.1 4.7 5.1 3.4 3.8 4.6 3.6 4.0 5.0 4.2 5.3 5.4
##  [361] 5.0 3.1 3.8 5.1 4.6 3.6 3.5 5.6 5.4 3.7 3.6 2.2 4.4 4.6 5.3 3.8 5.5 4.9
##  [379] 3.6 4.5 5.5 4.0 3.5 2.1 3.9 4.2 4.6 4.8 3.9 5.2 5.1 3.8 5.1 4.1 5.0 4.9
##  [397] 6.3 4.3 5.7 3.6 6.2 4.6 4.3 5.0 5.1 3.4 4.1 3.7 2.8 5.0 5.7 4.0 5.1 5.7
##  [415] 4.1 4.7 4.3 4.7 3.9 6.6 3.7 2.7 5.8 4.7 5.2 5.4 5.5 5.4 4.3 4.2 3.8 4.2
##  [433] 3.9 6.0 5.0 5.5 4.3 5.1 3.5 4.9 4.7 4.1 5.8 4.0 4.5 3.7 5.1 5.1 3.3 4.5
##  [451] 3.9 4.5 6.2 5.2 5.6 4.6 4.8 2.6 4.7 4.3 3.8 5.0 3.7 5.1 3.9 5.3 3.9 4.0
##  [469] 3.2 4.2 3.9 6.3 3.9 4.8 3.5 4.9 4.2 4.8 4.2 4.1 3.1 4.6 5.4 4.5 5.8 2.8
##  [487] 4.5 5.1 4.4 4.0 5.1 4.7 4.3 4.3 6.0 4.8 3.4 5.2 5.4 4.8 5.2 4.5 3.6 5.1
##  [505] 4.2 4.6 5.4 4.7 4.9 4.3 4.2 2.7 5.8 5.1 4.9 2.1 4.4 5.6 3.8 4.4 5.5 4.6
##  [523] 4.7 3.9 5.4 3.2 5.6 4.0 4.4 5.1 4.1 4.1 7.2 4.5 4.7 4.1 5.3 5.4 5.1 5.2
##  [541] 5.1 5.0 4.0 3.6 4.9 4.8 4.6 4.0 4.8 4.4 2.8 3.9 3.3 4.2 5.6 4.4 3.7 3.1
##  [559] 3.4 5.0 3.8 5.6 4.7 5.3 4.3 5.2 3.2 2.7 3.4 5.9 4.8 3.1 3.8 3.4 4.5 4.4
##  [577] 3.3 3.8 3.3 4.8 6.6 4.0 4.3 4.5 5.7 5.0 5.0 4.4 4.0 4.4 3.2 5.5 5.3 4.3
##  [595] 4.1 4.4 5.7 5.4 4.5 3.6 4.7 5.5 4.5 5.3 4.0 5.2 5.6 3.9 5.2 6.1 5.6 3.6
##  [613] 5.4 5.9 5.2 4.6 5.8 4.0 5.9 4.6 5.2 3.5 4.7 4.0 3.7 4.1 3.6 4.8 5.1 2.7
##  [631] 4.9 3.1 3.8 3.6 5.4 5.3 3.4 4.3 4.4 4.7 4.2 4.7 3.9 4.4 4.4 3.9 5.0 3.9
##  [649] 6.6 5.3 2.4 5.2 5.3 3.4 3.1 5.0 4.8 3.8 4.2 5.7 3.5 4.9 5.5 5.1 4.1 5.0
##  [667] 5.0 4.7 5.7 4.2 6.8 5.1 2.8 3.0 4.1 3.3 4.8 5.1 5.0 5.3 5.3 3.7 4.0 3.0
##  [685] 6.0 5.8 4.7 3.4 2.6 5.3 6.3 3.4 3.7 3.6 4.6 5.7 6.0 4.3 5.0 5.0 4.5 5.0
##  [703] 3.9 5.5 3.8 4.5 3.3 4.0 4.8 5.8 5.0 3.8 2.6 4.6 4.5 5.1 5.6 3.4 4.2 5.9
##  [721] 3.5 5.8 2.9 4.8 3.6 4.4 4.4 5.5 5.2 3.0 3.6 6.7 4.2 4.1 6.4 4.2 3.7 4.9
##  [739] 4.5 3.8 4.9 3.0 4.2 3.8 6.3 5.1 4.4 5.4 4.6 2.8 5.3 5.3 5.8 5.4 4.3 3.3
##  [757] 4.8 3.5 3.1 4.0 6.4 4.0 4.2 3.0 5.5 4.8 5.2 4.3 5.0 4.2 4.8 2.4 3.9 5.6
##  [775] 5.6 1.9 5.0 5.1 3.8 5.5 4.3 3.7 3.1 3.5 4.0 3.3 6.0 5.4 5.4 3.9 3.9 4.1
##  [793] 5.1 4.0 4.0 5.3 4.1 3.9 3.3 4.9 3.5 5.5 3.5 3.6 4.9 4.3 4.9 3.3 4.9 4.4
##  [811] 4.1 5.9 2.3 4.1 4.2 5.5 5.1 3.7 5.8 3.8 6.0 4.5 3.1 4.0 4.0 4.2 4.9 5.5
##  [829] 3.8 4.0 4.5 3.3 5.5 5.3 6.2 5.6 4.0 2.3 3.9 4.2 4.8 3.6 3.8 6.0 5.2 4.7
##  [847] 4.6 4.5 3.7 5.1 4.1 4.2 3.5 4.0 5.6 4.7 5.9 4.2 4.9 4.4 4.3 2.6 3.7 5.9
##  [865] 4.7 4.1 2.7 2.8 5.2 3.7 5.5 5.3 4.2 3.1 4.8 4.3 2.7 4.0 3.7 5.3 5.4 5.6
##  [883] 3.3 5.5 4.8 6.4 5.4 4.3 4.2 4.4 5.3 5.0 4.7 5.9 4.4 3.7 6.4 5.8 3.8 2.4
##  [901] 4.3 3.4 5.9 4.6 4.8 5.1 4.7 3.9 5.3 4.2 3.6 3.3 2.6 3.6 5.2 4.0 4.3 6.3
##  [919] 4.2 5.0 3.9 6.0 5.1 4.5 3.0 4.7 4.1 3.5 3.6 5.0 4.8 3.7 5.2 4.8 3.6 3.8
##  [937] 3.8 4.5 5.7 4.5 5.5 5.5 4.3 4.5 5.2 3.0 4.7 5.5 4.4 3.1 4.1 4.4 5.2 5.6
##  [955] 6.3 4.3 3.5 5.3 6.0 5.5 4.7 5.1 4.4 5.2 4.0 4.0 4.3 4.2 4.3 2.9 4.4 5.0
##  [973] 3.6 2.8 6.0 4.5 4.0 5.8 4.3 5.6 5.3 3.8 4.0 4.6 4.7 4.3 4.0 5.0 3.8 4.1
##  [991] 2.3 5.2 3.4 4.9 4.7 3.8 4.5 4.6 2.8 4.6
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
##    [1] 4.8 5.2 4.9 4.6 4.7 5.9 2.4 5.3 5.2 6.0 2.5 5.5 5.0 4.1 5.9 5.2 4.8 5.4
##   [19] 3.7 3.4 5.1 4.7 5.2 4.3 4.6 4.5 4.7 4.1 5.3 3.4 5.5 3.1 2.9 3.7 3.6 4.3
##   [37] 3.5 4.8 3.3 3.5 4.7 4.9 4.1 5.6 4.4 4.1 7.2 3.2 4.3 3.2 3.5 5.1 4.6 3.4
##   [55] 3.3 3.9 6.0 3.6 4.9 5.0 2.9 4.7 4.1 5.5 4.9 3.8 4.5 4.6 2.5 4.8 4.7 2.5
##   [73] 5.9 4.2 4.9 5.1 4.1 4.6 3.9 4.8 4.4 3.7 4.7 5.5 4.8 4.9 4.9 4.2 4.7 4.9
##   [91] 4.9 6.0 3.8 5.4 5.3 5.0 2.4 3.9 3.2 4.8 4.2 4.2 3.5 6.0 5.1 5.4 3.0 5.7
##  [109] 5.5 4.1 5.4 3.4 4.7 4.7 3.8 5.0 4.9 3.5 5.5 3.0 5.4 4.6 5.7 3.5 4.5 4.4
##  [127] 4.4 3.8 4.3 5.3 4.4 3.6 5.9 5.2 5.1 5.4 4.7 3.4 2.8 4.6 5.1 3.9 4.1 4.0
##  [145] 5.1 5.9 4.1 4.5 3.6 4.5 6.8 4.1 5.1 4.3 5.1 3.7 3.1 5.9 4.4 4.3 3.9 4.6
##  [163] 4.0 4.9 4.5 4.9 4.7 5.5 6.2 5.1 5.1 4.9 4.5 4.7 4.0 4.8 3.9 4.8 5.2 5.0
##  [181] 4.7 5.5 3.8 4.9 6.5 4.7 5.2 5.0 3.7 3.6 6.5 5.1 5.5 4.7 3.6 4.1 4.8 5.2
##  [199] 3.9 3.7 3.3 5.5 7.1 4.5 4.5 4.3 4.7 3.6 4.4 4.3 3.9 5.3 3.9 3.1 3.9 3.9
##  [217] 5.2 3.5 4.2 5.1 4.3 3.7 6.2 4.7 5.0 3.6 6.1 4.1 4.4 4.2 4.0 4.3 5.1 5.9
##  [235] 3.5 4.7 4.4 5.0 4.1 3.3 5.7 3.2 4.9 3.8 4.1 6.1 4.5 4.1 3.4 3.9 4.9 3.9
##  [253] 4.2 4.0 5.6 5.0 3.4 4.5 3.6 5.1 5.1 4.7 6.3 4.2 4.1 4.4 5.4 4.2 4.5 4.8
##  [271] 4.1 4.7 5.1 5.3 4.1 5.6 4.4 5.4 5.6 4.0 4.4 3.7 5.9 3.4 5.5 3.9 1.9 3.0
##  [289] 4.2 4.6 3.5 6.1 5.7 5.3 3.2 4.1 4.6 4.7 4.7 4.7 5.8 5.2 3.3 5.0 4.7 2.7
##  [307] 3.3 3.8 3.4 3.9 4.8 3.8 3.9 4.9 4.0 2.9 4.7 4.2 4.9 3.8 4.4 5.4 7.0 4.4
##  [325] 4.4 5.2 4.7 4.6 4.3 3.3 6.2 4.2 4.3 4.6 2.9 3.3 4.8 4.5 3.6 5.4 5.8 4.3
##  [343] 4.1 5.2 5.6 3.7 3.1 3.6 5.1 3.3 5.9 3.9 4.2 5.3 4.4 3.4 3.3 3.5 6.0 4.2
##  [361] 4.1 4.1 5.1 4.6 4.8 5.1 3.8 5.7 4.5 4.7 3.8 4.9 5.5 3.6 2.3 3.8 3.3 4.7
##  [379] 3.8 3.9 4.2 5.6 4.7 4.5 5.1 3.2 6.1 4.7 3.7 5.2 5.5 3.9 4.3 3.6 5.0 3.4
##  [397] 4.8 6.2 5.5 5.7 4.4 5.5 4.0 5.1 4.4 4.1 5.2 4.9 3.9 5.0 5.2 3.6 4.6 3.9
##  [415] 3.1 3.5 5.6 6.7 3.7 5.5 4.3 3.1 4.6 3.8 5.0 4.8 5.5 5.0 4.5 5.8 4.5 3.9
##  [433] 4.5 2.8 5.9 6.0 4.3 5.5 4.1 5.9 5.3 4.5 4.0 5.2 4.9 5.2 4.9 5.6 4.7 3.9
##  [451] 4.0 4.5 4.4 3.7 5.1 5.1 4.9 4.2 4.8 4.4 5.5 5.2 5.1 3.9 5.6 5.0 5.5 5.4
##  [469] 3.8 4.7 4.5 5.9 4.4 5.3 5.5 3.8 2.4 3.7 4.0 4.7 4.5 4.5 3.5 3.7 4.1 5.4
##  [487] 4.8 4.3 4.2 3.0 3.4 6.0 4.6 5.1 5.0 4.7 6.9 4.7 4.7 5.1 3.0 4.5 3.8 4.7
##  [505] 4.3 3.6 4.0 4.4 3.6 3.4 3.7 4.1 4.2 5.2 4.1 3.8 4.8 6.1 3.2 4.8 4.7 3.8
##  [523] 4.0 4.6 4.9 3.9 5.1 5.5 4.1 3.5 5.3 4.8 3.5 4.1 5.3 4.0 4.0 4.6 4.0 6.1
##  [541] 4.5 3.8 4.5 3.3 5.3 5.0 5.3 4.9 5.0 5.0 4.2 5.8 4.8 3.9 4.0 4.6 4.3 6.2
##  [559] 3.3 4.6 4.5 5.5 5.7 4.2 3.8 4.0 6.2 4.7 4.7 3.8 5.2 4.7 4.2 4.7 5.0 3.5
##  [577] 5.0 2.9 5.2 4.7 4.7 6.7 4.4 4.3 4.4 5.0 4.4 3.3 5.7 5.1 3.4 3.7 4.2 4.4
##  [595] 6.2 4.6 4.8 4.8 5.3 4.3 5.5 4.6 5.3 3.9 5.7 4.8 3.5 3.7 4.1 6.4 4.6 5.0
##  [613] 5.0 3.2 5.4 5.8 4.1 3.7 4.5 5.4 4.9 4.3 4.6 6.6 3.4 4.1 4.9 3.5 6.0 3.2
##  [631] 4.7 4.7 3.3 3.0 4.7 4.5 5.1 5.8 5.2 3.4 5.2 5.2 5.2 3.0 2.8 4.8 5.1 5.1
##  [649] 3.4 5.4 3.6 3.8 3.9 5.5 6.0 3.4 4.4 4.5 4.7 4.6 3.0 4.3 5.7 5.4 3.9 3.9
##  [667] 5.1 3.1 3.5 4.6 3.9 5.0 4.9 5.9 4.3 2.6 3.9 4.5 6.1 2.8 3.8 3.7 4.1 5.6
##  [685] 4.6 4.8 4.1 3.7 4.8 6.2 4.6 4.0 5.2 4.9 3.7 4.7 3.3 4.8 5.3 3.9 3.6 2.8
##  [703] 3.9 5.4 5.7 3.7 3.6 3.4 3.6 4.2 5.1 4.3 6.7 3.5 5.3 4.6 3.5 4.9 5.7 3.3
##  [721] 5.4 3.9 4.3 3.4 5.1 5.5 5.1 3.3 2.2 3.1 4.1 5.2 4.8 5.5 5.0 4.6 5.4 3.8
##  [739] 4.6 4.0 3.5 4.4 5.1 4.4 4.0 3.6 3.7 4.3 4.8 4.9 3.9 3.4 4.8 4.4 2.9 3.2
##  [757] 2.3 4.5 4.9 4.8 2.7 4.8 5.5 3.7 4.8 7.0 5.1 5.8 5.5 6.0 5.4 3.9 4.4 3.3
##  [775] 4.4 4.0 3.5 3.3 4.2 6.6 4.6 3.2 4.9 4.3 4.2 4.1 5.3 4.0 6.0 3.8 2.4 4.2
##  [793] 5.0 4.8 3.4 5.4 2.9 3.6 4.9 4.0 4.4 4.1 3.6 4.4 4.7 5.8 2.7 6.0 5.6 4.5
##  [811] 5.4 5.1 3.8 3.5 3.4 4.4 4.8 4.1 4.3 3.1 5.3 4.3 3.0 4.5 5.0 5.3 5.8 4.2
##  [829] 4.2 6.7 4.5 4.6 4.3 5.5 5.9 5.8 3.7 3.6 4.2 4.3 5.8 3.7 3.6 3.6 4.0 5.2
##  [847] 3.0 3.1 4.9 4.5 7.2 4.5 5.0 5.8 3.6 5.1 5.5 5.5 4.6 6.1 3.8 4.5 4.1 2.9
##  [865] 5.4 4.8 5.1 2.8 5.5 3.6 3.2 4.7 6.1 3.6 4.6 5.0 5.9 5.4 3.5 5.8 3.9 5.6
##  [883] 2.6 4.4 4.9 3.7 6.5 2.5 3.3 4.5 2.9 3.0 4.0 4.3 4.2 5.0 3.5 3.6 3.6 5.5
##  [901] 5.4 4.5 5.4 5.6 6.3 3.5 4.2 6.4 4.9 4.0 3.6 4.2 4.2 3.1 4.8 5.1 4.9 6.6
##  [919] 4.0 5.2 4.9 6.5 2.8 4.9 4.7 4.7 4.0 4.9 5.2 3.7 3.8 3.9 3.5 3.2 4.6 4.2
##  [937] 2.4 7.3 3.5 4.5 3.3 6.2 5.7 3.6 6.4 3.8 4.2 3.1 6.1 4.5 4.0 4.2 4.2 5.3
##  [955] 4.4 5.4 4.0 4.6 4.7 5.1 4.2 4.5 5.0 3.3 5.4 4.2 5.3 4.8 5.7 5.9 5.0 5.3
##  [973] 5.2 4.8 3.3 2.9 4.0 4.1 3.2 3.3 4.4 4.3 5.0 3.5 5.2 4.3 5.4 5.1 6.2 5.2
##  [991] 4.7 4.3 5.0 5.3 5.9 4.2 5.0 5.9 3.9 3.9
## 
## $func.thetastar
## [1] 0.0149
## 
## $jack.boot.val
##  [1]  0.49198895  0.37620321  0.31428571  0.17794562  0.15125000 -0.01984127
##  [7] -0.22042683 -0.19310345 -0.36914601 -0.45816619
## 
## $jack.boot.se
## [1] 0.9309937
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
##    [1] 5.6 3.1 5.9 3.4 5.2 4.6 4.2 4.6 5.3 5.2 6.2 4.2 4.6 3.0 5.4 4.6 4.8 4.1
##   [19] 4.9 6.0 5.3 4.2 5.0 6.1 5.4 3.6 4.0 4.1 4.8 3.7 3.7 4.9 5.4 5.9 3.4 4.9
##   [37] 3.8 4.5 4.2 4.7 2.8 4.6 4.5 4.6 5.0 6.3 3.9 5.1 5.6 3.6 5.5 4.1 5.2 5.2
##   [55] 4.1 4.6 5.1 4.6 4.5 3.4 4.3 5.0 4.5 3.8 3.5 2.9 5.5 4.7 5.5 4.3 5.0 4.8
##   [73] 4.6 5.1 3.9 4.3 3.5 3.9 4.9 4.5 4.4 4.6 4.6 4.3 4.6 6.5 4.7 3.7 3.8 5.9
##   [91] 3.6 4.3 3.1 2.2 4.8 4.7 3.1 3.2 4.3 6.5 4.6 3.4 2.9 4.1 5.4 5.2 4.7 4.7
##  [109] 6.4 5.5 5.1 4.7 5.2 5.2 5.3 2.7 4.4 3.6 5.7 3.5 3.6 6.4 5.1 4.5 5.0 4.0
##  [127] 4.9 5.0 2.4 5.0 4.9 3.8 4.9 5.7 6.1 4.2 4.5 3.5 5.2 5.2 4.3 4.7 4.7 4.1
##  [145] 4.9 4.3 4.3 5.5 5.0 3.0 5.1 4.3 3.5 2.8 5.1 4.9 5.2 3.0 4.7 4.2 3.6 6.0
##  [163] 3.3 3.7 4.9 4.9 3.4 3.9 5.8 4.3 2.7 6.5 4.0 4.4 3.7 4.9 4.7 4.6 3.3 2.8
##  [181] 3.5 5.4 4.1 4.0 4.1 5.7 5.2 5.1 6.0 3.5 5.3 3.7 3.9 4.7 5.0 3.3 5.0 3.7
##  [199] 4.4 3.8 4.9 4.5 4.8 5.7 6.1 5.1 3.6 3.7 6.2 5.4 5.7 5.7 4.7 3.2 4.7 6.5
##  [217] 4.2 4.2 4.8 3.5 4.2 3.2 5.7 5.4 4.3 5.4 5.7 2.7 5.7 3.3 4.3 4.9 4.1 5.0
##  [235] 4.8 5.2 3.8 5.0 4.5 4.7 4.1 4.0 4.3 2.7 5.9 4.9 5.4 4.5 4.5 3.4 3.0 4.7
##  [253] 3.3 4.5 3.9 4.4 3.8 3.1 5.4 4.5 5.3 4.0 5.1 4.2 5.8 5.4 3.8 6.1 6.5 4.4
##  [271] 5.9 5.9 4.8 4.9 4.1 6.0 5.0 4.6 6.5 2.8 4.4 4.7 4.8 3.2 4.6 3.1 4.4 4.0
##  [289] 4.6 5.3 6.2 3.1 4.6 4.8 4.0 3.1 5.0 4.0 5.4 3.8 5.4 4.5 4.5 4.9 3.5 5.0
##  [307] 5.3 5.1 4.6 5.1 3.8 4.9 4.6 5.7 2.3 4.1 4.7 5.2 4.9 5.0 3.5 4.4 5.0 5.4
##  [325] 5.1 2.8 3.9 6.5 3.1 4.2 5.1 4.0 4.6 5.8 5.0 3.4 5.5 4.8 6.4 4.4 2.8 5.3
##  [343] 4.4 4.6 4.6 5.1 6.1 3.0 3.8 2.7 4.5 2.9 3.6 3.5 4.0 3.5 5.3 4.8 2.2 5.0
##  [361] 3.3 4.4 4.8 4.1 5.3 2.8 4.7 4.1 5.5 5.2 5.1 4.1 3.4 5.9 3.8 4.8 4.1 4.2
##  [379] 3.7 4.3 4.5 5.0 6.3 4.3 6.1 4.3 3.3 3.9 3.3 3.5 3.9 3.3 3.3 3.5 4.3 2.4
##  [397] 2.9 4.7 4.6 4.2 3.6 4.9 5.5 4.2 4.5 3.9 3.4 5.2 4.5 5.7 3.5 3.8 3.8 4.3
##  [415] 5.3 3.3 3.9 4.3 5.8 3.9 4.6 5.0 4.3 5.5 4.8 4.9 3.3 5.2 2.3 3.6 5.0 4.2
##  [433] 4.0 4.1 5.0 3.5 5.5 4.3 5.3 6.0 4.8 4.4 3.5 4.8 3.2 2.8 5.9 3.8 4.7 3.3
##  [451] 3.8 4.2 4.6 3.1 4.3 3.3 3.3 4.1 5.0 3.9 4.1 5.4 4.5 4.3 4.4 3.4 2.6 5.1
##  [469] 5.0 5.4 3.7 2.2 4.8 4.8 5.6 3.9 5.0 3.8 3.0 5.4 3.7 4.2 3.0 3.5 4.7 5.5
##  [487] 4.4 5.1 6.2 5.4 4.1 4.4 6.0 4.7 3.6 5.8 3.8 4.3 4.6 4.8 5.0 4.1 4.5 3.8
##  [505] 3.8 4.5 4.1 6.1 3.0 5.1 5.1 6.0 5.2 4.3 4.7 4.6 5.6 4.4 4.2 4.1 2.9 6.1
##  [523] 6.0 3.0 4.7 5.8 3.8 4.7 4.1 3.8 5.5 4.1 4.6 5.8 3.6 4.3 4.4 3.9 4.8 3.7
##  [541] 4.4 3.2 4.1 5.7 3.7 4.8 4.9 5.0 3.3 4.9 4.3 6.1 6.5 4.8 4.5 4.0 4.5 5.0
##  [559] 3.3 5.8 6.1 3.7 4.7 3.4 5.3 5.5 3.8 3.8 4.2 5.3 4.8 3.4 4.3 4.6 5.5 4.7
##  [577] 4.0 6.0 4.0 3.3 7.3 3.8 3.3 3.5 4.0 3.9 4.6 5.9 5.1 3.0 3.0 3.8 2.5 4.0
##  [595] 4.3 4.8 4.6 5.2 5.6 3.7 4.2 3.2 4.2 4.2 4.5 6.4 4.0 4.0 5.8 4.8 2.9 7.0
##  [613] 4.5 4.3 3.8 6.6 3.7 3.7 3.8 5.2 5.0 4.2 5.4 4.1 6.3 2.4 4.5 3.9 3.7 3.8
##  [631] 4.8 2.7 4.6 4.1 4.3 4.6 4.5 5.1 5.7 4.5 3.1 6.2 6.1 4.1 3.7 5.2 3.5 3.8
##  [649] 4.6 5.5 5.6 3.9 5.8 4.5 6.2 6.2 4.1 4.3 5.2 5.1 6.1 3.6 4.8 4.2 3.3 4.0
##  [667] 4.8 4.4 4.2 4.4 3.1 4.8 4.8 6.1 6.2 7.1 5.7 3.3 3.3 4.2 3.9 4.0 4.5 4.6
##  [685] 4.2 5.3 5.8 3.8 3.9 4.2 5.0 4.5 4.2 4.2 3.4 4.6 4.0 4.4 5.3 3.0 5.4 5.4
##  [703] 3.9 4.2 4.3 4.9 5.0 5.0 4.0 4.7 4.1 4.1 5.3 3.4 3.8 6.1 5.0 3.0 3.7 3.2
##  [721] 4.6 5.3 4.5 4.0 3.9 3.1 3.9 3.6 4.7 3.6 5.5 2.3 3.7 5.8 4.4 4.2 3.8 3.4
##  [739] 4.7 4.6 5.8 2.7 5.1 3.5 4.2 4.6 3.5 2.6 5.3 5.3 4.6 5.2 3.9 3.9 2.1 4.4
##  [757] 5.1 4.4 4.8 3.4 5.0 5.0 4.2 5.5 3.3 6.7 4.7 4.2 4.4 4.9 4.7 6.0 4.1 5.6
##  [775] 4.0 3.8 4.2 4.4 5.7 5.0 5.9 5.3 6.8 4.4 3.9 4.0 4.3 4.6 4.2 4.8 4.9 4.4
##  [793] 3.1 5.4 3.7 5.0 5.2 5.7 3.7 5.4 4.3 3.3 5.8 4.3 4.9 3.5 4.8 3.9 4.9 4.4
##  [811] 4.2 5.2 4.1 4.6 4.2 5.0 5.2 3.3 5.7 5.8 6.1 3.4 3.9 4.5 6.5 4.1 3.8 6.3
##  [829] 5.6 5.6 4.2 4.8 3.8 3.9 4.7 4.7 3.9 5.7 5.7 4.3 4.3 4.7 4.2 3.3 5.7 3.9
##  [847] 3.7 4.6 5.7 5.3 3.8 4.9 4.9 4.8 6.1 5.7 4.6 4.0 5.7 4.5 5.7 4.5 3.2 4.7
##  [865] 4.8 4.1 4.1 4.4 3.8 4.5 4.6 3.6 3.6 3.8 4.7 4.9 4.5 3.4 4.9 4.7 3.7 4.3
##  [883] 3.3 4.1 5.5 4.1 4.2 5.1 4.6 5.8 4.6 3.6 4.2 3.9 5.0 6.1 4.9 4.1 5.0 3.8
##  [901] 3.2 4.6 4.1 3.5 3.7 5.9 3.4 6.4 4.4 4.0 6.6 3.9 2.9 3.4 6.2 3.8 3.8 5.8
##  [919] 5.1 3.0 3.3 5.2 4.7 4.2 5.2 5.0 4.1 4.4 4.4 5.1 5.1 3.9 3.8 5.7 5.8 4.2
##  [937] 4.7 5.0 4.4 5.7 3.1 5.3 4.4 4.8 5.3 5.6 4.3 4.4 3.3 4.8 4.2 4.3 3.7 4.9
##  [955] 5.0 4.3 3.8 5.3 3.7 3.9 3.1 5.1 5.1 4.0 5.1 2.9 3.0 5.0 6.2 2.2 3.6 4.7
##  [973] 3.5 6.0 4.0 3.4 4.3 3.8 2.5 5.3 4.1 4.3 4.8 3.2 4.2 4.4 4.3 5.5 4.6 4.1
##  [991] 5.6 5.2 4.4 4.7 4.2 3.0 4.1 3.6 4.4 4.9
## 
## $func.thetastar
## 72% 
##   5 
## 
## $jack.boot.val
##  [1] 5.5 5.4 5.3 5.2 5.0 5.0 4.8 4.7 4.7 4.4
## 
## $jack.boot.se
## [1] 1.003992
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
## [1] -0.04076961
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
##   2.523521   4.960962 
##  (1.062438) (2.310339)
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
## [1]  0.42880221  0.92734039  1.48563511 -0.07634401 -0.30635233  0.47677905
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
##    [1] -1.220528700 -0.463804783  0.394846902 -0.298198568  1.454697761
##    [6]  0.735818763  0.793747312  0.775652067 -0.007765045 -0.836094524
##   [11] -0.787243507 -0.083957193  0.310008679  0.663988031  0.006635307
##   [16] -0.099406931  0.034281731 -0.326790350 -0.046357633 -0.765856609
##   [21] -0.746228030 -0.044139325 -0.020527792 -0.397303782  0.743379822
##   [26]  0.298429670  0.413098844 -0.332213798  0.341505390  0.769235068
##   [31]  0.773780740  0.330606240  0.352359257 -0.323465650 -0.365283659
##   [36] -0.218028305 -0.333502037 -0.118948506 -0.355914837  0.809472867
##   [41] -0.361068280 -1.300726119  0.769592359 -0.137540682 -0.459041639
##   [46]  2.280026122 -0.110477523  0.241937183 -1.260446547 -0.317847598
##   [51]  0.754592869  0.314072453  0.412501697  0.623345836 -0.368240282
##   [56]  0.376186268  0.679261549  0.284439069 -0.025803483 -0.022052211
##   [61] -0.475259101 -0.326790350  0.291449855 -0.513316406 -0.074340806
##   [66]  0.017626586  0.031818503  0.255089822 -0.081536731  1.418246449
##   [71] -1.320921060  0.320983536 -0.029851453  0.453479204 -0.467089929
##   [76] -0.286291573 -1.104946448  0.438649174 -1.294478389  1.953032691
##   [81] -0.343404970 -0.417510627  0.636152141 -0.428204637 -0.360962244
##   [86]  0.221873436 -0.034985248  0.036010400 -0.768368413  0.667033162
##   [91] -0.070384147 -0.886712822 -1.236520984 -0.068624167 -0.785088504
##   [96] -0.468222770  0.276719595 -0.732633567  0.021079064 -1.108249942
##  [101]  0.731269706  0.848636575 -0.293150015  0.251085421 -0.902132135
##  [106] -0.351495488 -0.850419092 -0.031934671 -0.105603190 -0.910998642
##  [111]  0.368035357  0.401757577  0.039156077  0.363919486 -0.371086743
##  [116] -0.927800560 -0.442207116 -0.867131171 -0.050981558 -0.088371948
##  [121] -0.299202605 -0.656719390  0.406969392  0.793446954 -0.377511910
##  [126]  0.614194915 -0.057310687 -0.620208577 -0.016935698 -0.009141326
##  [131]  0.339100814 -0.074148795 -0.396829314 -0.363820405  0.725222072
##  [136] -0.275288093 -0.056331856  0.339115824 -0.117597897 -1.357258133
##  [141]  0.326100552 -0.408827650 -0.485321608  0.639704104  0.286713515
##  [146]  1.252799553 -0.869840233  0.045853056  0.015981006 -1.362909886
##  [151] -0.365917108 -0.704108120  0.757121171  0.304348733  0.326107050
##  [156]  0.688071534 -0.869743596  0.042568108 -0.363956255 -0.498569881
##  [161] -0.334554120  0.343947820  0.768873970  0.086440415 -0.400940732
##  [166] -0.862031829 -0.470362397 -0.401097420  0.253496330 -1.362289959
##  [171] -0.124199295 -0.352086781  0.821758117  0.390977302 -0.818022742
##  [176] -0.870160125  0.291097670  1.408138431  0.676236804  0.284994857
##  [181] -0.017279270 -0.119746881  0.301116062  0.042213193  0.688612325
##  [186]  0.022802244 -0.686060207  0.296815333 -0.698440218 -0.061668284
##  [191] -0.549944971  0.720942363 -0.288349860 -0.072972561 -0.827223611
##  [196]  0.261142257 -0.452115034 -0.039594944  1.923981768  0.329741683
##  [201] -0.777721765  0.691184567  0.348904083  1.281873367 -0.464488718
##  [206] -0.417842898 -0.071287684  0.277011714 -0.046340680  0.847237628
##  [211]  0.327193160  0.762300565 -0.358759553 -0.805386064  0.303782955
##  [216] -0.429060123  0.727556423 -0.384009051  0.708128790  0.567912220
##  [221]  0.808424255 -0.451831877  0.421440068 -0.402367319  0.288839022
##  [226] -0.492852299 -1.891304606 -0.361094797 -0.006248588  0.649693823
##  [231] -0.841286788  0.339000404 -0.623962470  0.196275385  0.756516174
##  [236] -0.442975844 -0.837372232 -0.111517456 -0.020931803 -0.053126631
##  [241] -0.395250184  1.413349123 -0.138289485  0.331256725  0.363919486
##  [246]  0.908450734 -0.100513007  1.145676222  0.800835700  0.740390930
##  [251] -0.016125911  0.343399384  0.205309338  0.288140430 -1.982360435
##  [256] -0.281946160  0.024364679 -0.356305848 -0.916395778 -0.836278663
##  [261] -0.865708751 -1.159969647  0.051476918 -0.827020159  0.635145610
##  [266] -0.337714148  0.241371944  0.019825671  0.775871732 -0.479078538
##  [271]  0.373603134 -0.344834203 -0.348233867  0.049375449 -0.418090741
##  [276] -0.020842002  0.309999138  0.361316496 -0.085239742 -0.416869511
##  [281] -0.417828961 -0.932199631 -0.009697209 -0.022219178 -0.386752689
##  [286] -0.805088966 -0.298808828  0.594265857 -0.051711994 -0.772154256
##  [291]  0.316122844  0.768015295  0.280692750 -0.775696657 -1.196088518
##  [296]  0.629004809 -0.377965084 -0.522853281 -0.457083572 -0.455881782
##  [301] -0.409160082  0.246223930  0.003674913  0.304397162  0.369115932
##  [306] -0.054590503  0.301791185 -0.400721247 -0.071903799  0.073283719
##  [311] -0.129321234  1.950284728 -0.479983071  0.864534333  0.783344699
##  [316] -0.423372595  0.751344026 -0.347355937 -0.112098757 -0.879935570
##  [321] -0.016405686  0.048808496  0.554339868  0.734153589  0.422295695
##  [326] -0.833312172 -0.057017135 -0.071711403  1.061473445 -0.853039055
##  [331] -0.062732433  0.677237901 -0.078602797 -0.126786930  0.374845430
##  [336] -0.101324395 -0.053141963  0.040601602  0.012437195  0.329349085
##  [341]  0.611497319 -0.726421128 -0.369927230 -1.358740409  1.044631389
##  [346]  1.195848030  0.418255240  0.320102356 -0.354371736  0.276719595
##  [351] -0.059432344 -0.512543794  0.720709462  0.398830898  0.410809343
##  [356] -0.494964585 -0.074501475  0.237265999  0.405555424  0.343969289
##  [361]  0.287142945  0.825888510  0.468515172  0.418783531 -0.769718728
##  [366] -0.290950318 -0.256213301  1.137948800  1.161642816  0.017810362
##  [371] -0.474229453 -0.501593056  0.771500843 -0.357511875 -0.093056262
##  [376] -0.311970289 -0.529450138  0.025300786 -0.002777765 -0.023931805
##  [381] -0.957027568  0.771225499 -0.137536515  0.653484718  0.351451510
##  [386] -0.835136529 -0.501980936 -0.400509060  1.417410188 -0.345094108
##  [391] -0.384569727 -0.343240622 -0.378762961 -0.404700816  0.256170700
##  [396]  0.059423738 -0.484000055  0.759383088  0.292520708  0.319760041
##  [401] -0.016947355  1.274868383  0.282735319 -0.665730277  0.680705520
##  [406]  0.733027415  1.776169693  0.236138307  0.392806610  0.110802400
##  [411] -0.001011744 -0.116334867  0.858111805  0.720373892 -0.450819677
##  [416] -0.361993044  0.380502135  0.252968344 -0.422227382  0.795106123
##  [421] -0.857335155 -0.097646247 -0.939478480 -0.047485752 -0.508297105
##  [426] -0.322816638  1.349585353 -0.001978845 -0.854696047 -0.029851453
##  [431] -0.430231959 -0.478040343 -0.078191369  0.364531117  0.244672918
##  [436]  1.319899764 -0.469646213  1.100505820  0.363310282  0.340358900
##  [441] -0.018855410  0.434456227 -0.445464915 -0.031370876 -0.460963865
##  [446]  0.989410044  0.245683948 -0.049275403 -0.726421128  0.613051979
##  [451] -0.349330647 -0.024257175  0.298429670  0.990519237 -0.852756054
##  [456]  0.249377361  0.839538091 -0.409184093 -0.020456994  1.404661548
##  [461]  0.700457287  0.646039928  0.801462184  0.054251615 -0.129321234
##  [466] -1.276586485 -0.756887523  0.051036303 -0.052888467 -0.788971208
##  [471] -0.325029767  0.580348456 -0.395406809 -0.086482876  0.361915582
##  [476] -0.118133696 -0.418071697  0.342244568 -0.780375137  0.366930606
##  [481]  0.848980088 -0.356328301 -0.060439665  0.778149582  0.421307334
##  [486]  0.374369536 -0.842890019 -0.942727013 -0.779877396  0.045284719
##  [491]  0.700962097  1.276957237  0.033120402 -0.063804274 -0.391965946
##  [496] -0.314943645  0.367864417  0.730656158  1.391234361  0.332243222
##  [501] -0.025860854 -0.455965974 -0.912928184 -0.034063909  1.149140238
##  [506] -0.009014134  1.081123215  0.564454822 -0.367639389  0.647567475
##  [511]  0.294636706 -0.347972116 -0.025391710 -0.085724903 -0.405976003
##  [516]  0.396857734  2.252011204 -0.129694282  0.109748382  0.048720508
##  [521]  0.341252522  0.328630704  0.751071992  0.717951723 -0.034501591
##  [526]  0.016265634  0.386691750  0.106018994 -1.049916103  1.321804410
##  [531]  0.292327923 -1.287936431  0.331311256  0.341977220 -0.743378989
##  [536]  0.321186912 -0.858172394 -0.027325852  0.833461111 -0.125364551
##  [541]  0.056413874  0.358747819  0.601405440 -0.811469253  1.263770815
##  [546] -0.107902631 -0.472502307 -0.392308785 -0.810884803 -0.447303692
##  [551]  0.444514718 -0.317229127 -0.002169633  0.398169317  0.031722649
##  [556] -1.362291567 -1.436925718  0.274078033 -0.800843428  0.395083187
##  [561]  0.683706074  0.868861881  0.408104675  0.273523537 -0.393641290
##  [566] -0.124086389  0.332619397  0.813744499 -0.764314869 -0.376464649
##  [571] -0.073676399 -0.417406876 -0.494637052  0.379360481 -0.010737986
##  [576]  0.273572448 -0.048322781  0.010898069  0.382672976  0.817674397
##  [581] -2.398552604  0.725059577  0.669379462 -0.007008701  0.694840999
##  [586] -0.353189015  0.333694174  0.334746237 -0.496997027  1.193364839
##  [591] -0.401283333  0.738504832  0.679314165 -0.404431047 -0.314393770
##  [596] -0.300245345 -0.051938766  0.381495330 -0.841986059 -0.939073154
##  [601]  0.301564847  0.261522915  0.284095508 -0.294468067 -0.356118802
##  [606]  0.367544501  0.771071251  0.357886806 -0.774267662 -0.778753492
##  [611]  0.374889290  0.847776403  0.390523846 -0.398457338  0.815605587
##  [616]  0.054118196  0.432990549 -0.408514368  0.031742741 -1.258654630
##  [621] -0.437768018  0.665718155  0.034664901 -0.840555812 -0.481517107
##  [626]  0.350092120 -0.356222296  0.332428590 -2.094116395  0.330606240
##  [631]  0.391868212 -0.457515993 -0.405831385 -0.783364009 -0.897373458
##  [636] -0.479114477 -1.318080814  0.797080206  0.645391149 -0.081233300
##  [641]  0.488263474 -1.374159805  1.189535488 -0.371579308  0.241960109
##  [646] -0.006542046  0.619411380  0.360011660 -0.348621757 -0.121949472
##  [651] -1.320852354 -0.059662735 -0.108637674 -0.818355340  0.719923448
##  [656]  0.014268055  1.435375611 -0.345780887  0.287284566 -0.719586935
##  [661]  0.325830772  0.080141386 -0.870614856 -1.291289323  1.235613858
##  [666] -0.068766120 -0.013434389 -0.038114963 -0.026913569 -0.381914890
##  [671]  1.102258594  0.674047889 -0.389019636 -1.341680358 -0.061216995
##  [676] -0.848518249  0.715094565 -0.427359205  0.330102708  0.308659158
##  [681]  0.706250021 -0.299908510 -0.016564056 -0.910027103  0.009869268
##  [686]  0.058181392 -1.382724784 -0.045375766  0.264687634 -0.422450971
##  [691] -0.004474834 -0.714385023  0.408493289 -0.862659902 -0.358759553
##  [696] -0.859600232 -0.357577848  0.312362037 -0.416061559 -0.769599267
##  [701]  0.302672199 -0.080276985  1.075158259 -0.275288093  0.658226490
##  [706]  0.339939502 -0.976075754  0.066459722  0.249797930  0.759053920
##  [711]  0.008728136  0.411412631 -0.005966205 -2.594883992  0.684872819
##  [716] -1.318473302  0.276960013  0.348617225 -0.874291214 -1.519043634
##  [721] -0.312649332 -0.080118927 -0.092559377 -1.054847255 -1.466843553
##  [726]  0.708234988 -0.803419557 -0.450334472  0.391996583  0.388606985
##  [731] -0.005107751  1.103553669 -0.366965451  0.434569188  0.024364679
##  [736] -0.415496661  0.326660733 -0.121259104  0.037548408 -0.094923299
##  [741] -0.472125398 -0.257121561  0.239215361 -0.333678708  0.787517711
##  [746] -0.301851713 -0.106905884  0.608189715 -0.403418162 -0.831606680
##  [751]  0.010898069 -1.416328190 -0.394008859 -0.027583005 -0.053799342
##  [756] -1.463728093 -0.030633630  0.035957065 -0.742321242 -0.359734926
##  [761]  0.720219758 -0.351255537 -0.819237702  0.318144251  0.298613701
##  [766]  0.699006931 -0.493740906  0.662171698  0.696807269 -0.053980736
##  [771]  0.727991025  0.475075882  0.640875535 -0.766628312 -0.807992273
##  [776]  0.583936124 -1.232973549 -0.058080503 -0.827775813 -0.055789573
##  [781]  0.795704385  0.321262909 -0.762058055 -0.426218466  0.042340455
##  [786] -0.418370767  0.048399398  0.666147160 -0.517893090 -0.358353853
##  [791]  0.445877232 -0.764165450  0.030199267  0.623536180  0.676900421
##  [796] -0.299457774 -0.730766669  0.368512729 -0.852931345 -0.336769761
##  [801] -0.312972874  0.417204548 -0.353686459 -0.056566604 -0.072163746
##  [806]  0.681831721 -0.322669576  0.027059253  0.331867918  0.722598106
##  [811]  0.300231039 -0.689179270 -0.947624481 -0.361957952 -0.394012062
##  [816]  0.076098653 -1.234709270 -0.060497404 -0.299864305  0.453876402
##  [821] -0.324046869  1.328201407 -0.004676831  0.052936791 -0.478657721
##  [826] -0.332213798 -0.071104963 -0.288017489 -0.730954897  0.788595581
##  [831] -0.830462421  0.384225908 -0.027626984 -0.117697977  0.322794833
##  [836]  0.337960620 -0.891090565 -0.015136697 -0.359450491  0.622507999
##  [841]  0.487760277 -0.423421121 -0.481517107 -1.109719446  0.039313869
##  [846] -0.739472162 -0.811285263 -0.102650868 -0.007521865  0.388017388
##  [851] -0.306557911  0.048136287  0.288083897  0.312644558 -0.404358739
##  [856] -0.061355632 -0.793462524  0.724392508 -0.037972447 -0.777620188
##  [861]  0.369882319 -0.742273496 -0.789190964  0.092164584 -0.083526591
##  [866]  0.808899942 -0.914891159 -1.361138578 -0.371298620 -0.414043780
##  [871] -0.037084181 -0.006680162 -0.009141326  0.282660110 -0.458532142
##  [876]  1.109671048 -0.393641290 -0.357242345  0.396383047  0.774537844
##  [881] -0.366482351  0.876804799  0.677106368 -0.804029541 -0.353745111
##  [886] -0.762249271 -0.927972661 -0.882471066  0.295926406  0.412263932
##  [891] -0.047690004  1.107851555 -0.017211918 -0.734279405  0.023292554
##  [896] -0.881626560 -0.016935698  1.269513936 -0.117652836  0.306345395
##  [901]  0.046459302 -0.471325177 -0.094141596  0.796748607 -0.389638386
##  [906]  0.232658962 -0.468776548 -0.343344870  0.046803115 -0.862690813
##  [911] -0.066719391 -0.502595543 -0.046739484  0.402137187  0.813899091
##  [916]  0.900484255  0.040723652 -0.382327049 -0.363381339  0.328345208
##  [921] -0.018552406  0.001635772 -0.437233370  0.430001194  0.325838274
##  [926]  0.585478383 -0.852087033 -0.845327157 -0.396829314 -0.320480583
##  [931] -0.896401568  0.417058018  0.351235254 -1.344249339 -0.903974898
##  [936] -0.117652836  0.464668341  1.193975361 -0.493450932  0.430722806
##  [941]  0.343082232 -0.048308233  1.186258297  0.044551509  0.329466618
##  [946] -0.129694282  0.364254431  0.330969059  0.391713583 -0.375943897
##  [951] -0.254764342 -0.402161648  0.339616603 -0.790859856 -0.343240622
##  [956]  0.775025844  0.057190014 -0.386827449 -0.127798205 -0.070699608
##  [961]  0.406867583  0.253427057 -0.374430150  0.754307898 -0.033515233
##  [966] -0.800835666 -0.716378050  1.330333625 -0.908854274  0.614029812
##  [971]  0.285374879 -1.400420783  0.072896123  0.816344484 -0.656240892
##  [976]  0.301564847 -0.338513365 -0.064383774 -0.030219396 -0.365546695
##  [981] -0.134243399 -0.901593186  0.471836038  0.319531937  0.441687527
##  [986] -0.822045292 -0.367171453 -0.040769609  0.279244027  1.212586178
##  [991]  0.018058390 -0.941262321 -0.035307399  0.344464824 -0.379612432
##  [996] -0.422317515  0.007663025  0.366565736  0.689026368 -1.159516573
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

## Warning in densfun(x, parm[1], parm[2], ...): NaNs produced

## Warning in densfun(x, parm[1], parm[2], ...): NaNs produced
```

```r
fit2
```

```
##       mean          sd    
##   0.50867319   0.28543984 
##  (0.09026400) (0.06382002)
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
## [1]  0.27618711 -0.27468452 -0.69388813 -0.13217103  0.15284288 -0.08496204
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
## [1] 0.0449
```

```r
sd(results2$thetastar)  #the standard deviation of the theta stars is the SE of the statistic (in this case, the mean)
```

```
## [1] 0.877145
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
## t1*      4.5 -0.03113113   0.8778628
```

One of the main advantages to the 'boot' package over the 'bootstrap' package is the nicer formatting of the output.

Going back to our original code, lets see how we could reproduce all of these numbers:


```r
table(sample(x,size=length(x),replace=T))
```

```
## 
## 2 3 5 6 7 8 9 
## 3 1 1 1 1 2 1
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
## [1] -0.0299
```

```r
se.boot
```

```
## [1] 0.8859802
```

Why do our numbers not agree exactly with those of the boot package? This is because our estimates of bias and standard error are just estimates, and they carry with them their own uncertainties. That is one of the reasons we might bother doing jackknife-after-bootstrap.

The 'boot' package has a LOT of functionality. If we have time, we will come back to some of these more complex functions later in the semester as we cover topics like regression and glm.

