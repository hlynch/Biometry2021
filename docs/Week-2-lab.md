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
## 0 1 2 3 4 6 7 
## 2 1 1 1 1 3 1
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
## [1] 0.0453
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
## [1] 2.736134
```

```r
UL.boot
```

```
## [1] 6.354466
```

Method #2: Simply take the quantiles of the bootstrap statistics


```r
quantile(xmeans,c(0.025,0.975))
```

```
##   2.5%  97.5% 
## 2.6975 6.3000
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
##    [1] 3.2 5.4 2.9 4.1 5.4 3.9 3.7 4.8 3.2 5.3 4.3 4.2 4.7 4.8 5.7 4.3 4.2 2.2
##   [19] 5.6 4.3 4.1 4.4 3.2 5.3 5.1 4.0 4.2 4.3 4.7 6.0 3.5 4.4 4.8 4.3 5.1 2.5
##   [37] 4.8 5.4 2.3 2.9 4.8 4.8 3.8 2.8 3.2 4.4 4.8 4.1 6.5 5.3 4.0 4.3 5.2 4.3
##   [55] 5.0 4.0 4.2 5.4 4.1 3.6 3.5 4.4 5.6 4.6 5.6 5.5 4.3 3.4 4.1 3.4 4.6 4.2
##   [73] 6.0 5.3 5.4 3.6 4.8 2.6 3.8 4.5 4.6 3.5 3.7 3.3 5.5 3.9 4.4 4.8 4.7 4.9
##   [91] 4.3 4.8 3.1 3.1 5.5 4.8 4.9 2.5 5.3 4.8 5.0 3.9 5.5 4.3 4.4 5.5 4.1 5.0
##  [109] 3.2 5.9 3.1 4.5 4.7 5.3 4.7 3.0 3.0 4.8 4.3 4.0 4.6 3.5 2.8 5.6 4.8 3.9
##  [127] 5.5 4.7 3.8 5.1 4.3 4.5 4.4 3.4 3.1 5.2 4.6 6.1 3.9 4.2 3.9 3.4 5.2 5.5
##  [145] 4.6 3.7 4.6 4.2 5.1 3.5 6.0 2.5 4.5 4.6 5.8 5.2 4.9 4.8 5.0 6.5 6.4 3.9
##  [163] 4.5 4.0 3.8 5.6 4.9 3.1 3.5 2.9 4.4 5.7 3.0 4.4 4.7 4.8 4.4 3.6 4.3 4.2
##  [181] 4.0 4.2 4.8 3.1 3.7 3.6 4.5 3.0 4.1 5.8 4.5 5.8 3.5 4.2 5.2 5.2 3.9 4.4
##  [199] 5.8 3.4 5.3 4.4 3.7 4.8 4.3 5.3 4.4 3.3 5.7 6.5 5.3 5.9 3.2 5.1 5.0 3.9
##  [217] 3.4 5.4 2.8 4.7 3.5 6.3 2.4 7.0 3.8 4.1 3.9 4.2 3.5 2.8 5.6 5.9 4.6 4.3
##  [235] 4.4 4.3 5.4 4.6 5.7 4.9 6.4 3.4 4.0 5.3 3.3 6.1 5.3 3.9 5.2 3.4 4.6 3.9
##  [253] 4.3 4.9 4.0 3.4 4.3 5.1 5.7 4.8 5.9 4.0 4.9 4.8 7.1 3.0 6.7 3.5 5.3 5.0
##  [271] 4.0 4.8 4.4 4.5 5.0 5.1 3.7 5.8 3.5 3.7 2.9 4.9 3.8 4.0 5.3 5.4 3.4 3.7
##  [289] 4.1 4.6 4.1 3.2 3.3 5.0 4.9 5.6 3.7 5.5 4.2 5.3 3.9 4.9 3.0 3.9 3.1 4.3
##  [307] 4.8 4.4 3.7 3.6 4.4 5.0 4.5 4.6 6.1 3.7 4.3 4.7 5.5 3.5 4.2 4.5 4.7 5.8
##  [325] 3.7 4.9 5.8 4.2 3.8 4.4 2.4 4.5 5.0 3.9 4.7 4.6 1.9 2.9 5.1 4.4 4.1 6.0
##  [343] 4.2 5.0 4.4 3.4 4.5 3.9 2.4 4.7 4.7 5.6 5.8 3.6 4.6 3.5 4.6 6.6 6.0 3.8
##  [361] 5.5 4.4 5.2 3.4 5.5 6.1 5.9 4.1 3.7 3.7 3.2 3.9 3.7 4.1 4.5 4.5 4.3 4.2
##  [379] 5.7 6.3 3.2 2.5 4.1 6.9 5.9 4.8 3.2 4.3 4.5 3.2 5.0 2.1 5.4 3.3 3.4 6.0
##  [397] 4.7 6.5 3.2 4.8 5.7 4.0 4.5 5.1 4.1 2.7 4.5 5.9 4.7 6.5 4.5 4.5 5.8 5.7
##  [415] 5.4 4.9 5.3 3.2 4.1 3.9 4.5 4.7 4.8 5.4 4.2 4.9 4.6 3.9 5.4 5.0 4.6 4.4
##  [433] 4.1 3.9 4.4 2.5 4.7 2.8 3.6 3.8 4.6 5.7 4.8 4.7 1.9 3.5 4.5 4.5 4.4 4.6
##  [451] 6.5 5.8 2.6 5.2 4.0 6.1 2.9 5.7 5.3 3.0 4.4 4.7 4.6 5.1 6.6 5.1 4.9 5.7
##  [469] 6.2 5.0 4.4 4.1 5.1 3.1 5.3 3.0 5.5 3.9 4.7 4.6 4.4 4.3 5.4 3.9 5.6 4.1
##  [487] 4.1 5.6 4.5 4.8 5.8 4.2 4.7 4.5 4.1 4.6 3.3 3.9 6.0 4.9 4.8 3.5 5.6 5.3
##  [505] 4.4 4.9 5.5 4.8 4.4 5.1 4.5 4.5 4.0 4.6 2.0 5.7 3.8 3.7 3.7 4.6 4.1 3.7
##  [523] 4.0 3.5 5.5 4.1 6.3 5.5 3.7 4.1 4.8 4.5 4.6 6.2 4.3 5.2 4.8 5.9 3.8 4.4
##  [541] 5.7 5.0 5.1 4.4 4.0 4.5 3.3 4.5 4.6 4.6 4.9 4.3 5.1 4.3 6.1 4.9 4.1 3.8
##  [559] 5.6 4.6 6.1 2.7 5.2 4.5 3.5 3.6 4.8 4.4 3.9 4.8 4.5 6.4 5.0 5.6 4.3 4.3
##  [577] 3.2 2.9 4.1 6.0 4.6 5.0 4.8 5.6 6.1 5.2 5.2 5.0 5.4 4.6 5.5 4.8 5.8 4.9
##  [595] 4.0 4.4 3.9 3.1 4.5 4.6 5.5 3.4 5.2 3.4 1.1 2.7 4.5 5.7 4.8 3.1 4.0 4.1
##  [613] 3.4 5.2 3.6 4.0 3.1 3.6 3.8 3.8 4.6 5.0 3.5 5.1 4.9 4.6 4.3 6.0 4.0 4.0
##  [631] 5.7 3.1 4.9 4.1 4.1 4.2 3.7 4.3 4.7 5.2 4.9 4.1 2.8 4.7 4.8 4.8 3.5 4.4
##  [649] 4.7 5.2 5.1 5.4 4.9 3.7 4.8 3.7 4.0 4.9 2.4 3.8 4.7 6.8 3.8 2.8 3.9 4.4
##  [667] 3.9 2.8 4.2 3.5 2.6 2.5 4.4 5.1 4.1 4.7 3.2 3.1 4.3 3.2 4.8 3.5 5.4 5.8
##  [685] 5.2 4.1 4.6 2.9 5.0 4.7 4.4 6.0 5.7 3.8 4.3 4.4 3.9 6.8 3.6 3.5 5.1 4.8
##  [703] 3.9 3.7 5.5 6.8 4.4 5.2 5.1 4.6 5.2 3.7 4.0 5.7 4.6 5.7 4.8 3.8 4.7 5.7
##  [721] 4.1 4.7 5.2 4.5 4.5 5.0 5.6 6.0 5.3 5.3 5.9 3.4 5.6 4.3 5.6 6.0 5.1 4.2
##  [739] 4.5 3.0 4.2 4.6 3.6 4.7 4.6 4.5 4.0 4.8 5.0 4.5 5.0 5.2 5.3 4.6 4.5 3.1
##  [757] 4.6 5.0 3.0 4.8 5.4 3.8 4.2 4.2 4.9 4.1 3.6 4.8 3.2 5.9 6.0 5.1 3.2 2.1
##  [775] 6.5 3.0 4.5 4.8 5.0 5.1 4.7 4.8 4.9 5.2 5.5 5.4 4.0 5.9 4.4 4.6 3.5 6.2
##  [793] 5.2 5.6 4.0 4.5 4.5 3.6 4.5 4.5 5.9 4.2 3.6 4.6 3.7 4.8 3.9 3.4 3.9 3.7
##  [811] 5.1 4.8 4.2 5.4 4.3 4.7 2.4 3.9 3.1 5.0 3.7 4.6 4.0 5.0 5.8 2.9 5.0 4.7
##  [829] 4.4 4.5 3.2 4.7 4.5 4.4 4.8 3.7 5.0 4.6 4.5 4.0 3.7 4.4 4.4 5.0 3.9 4.7
##  [847] 5.2 3.7 4.9 5.5 5.0 4.6 5.0 4.6 4.8 4.3 4.3 2.3 4.7 4.6 6.3 5.2 4.2 3.5
##  [865] 4.7 4.5 5.1 3.3 4.7 3.0 4.7 5.1 4.3 6.1 5.6 4.8 5.0 3.3 4.2 5.2 4.1 3.7
##  [883] 4.6 5.1 3.3 3.5 3.7 5.4 4.9 6.4 4.6 4.8 4.1 5.2 3.5 3.6 4.2 4.1 3.9 3.4
##  [901] 3.8 4.8 3.3 4.3 2.4 5.1 2.8 4.2 5.9 5.5 4.0 4.8 5.1 4.8 4.3 5.0 5.2 4.1
##  [919] 3.0 4.4 5.6 4.7 3.8 3.0 4.5 4.9 4.2 5.4 6.1 5.3 4.8 4.8 4.6 4.1 4.3 4.8
##  [937] 4.8 4.3 5.4 2.0 4.8 5.1 4.3 4.9 4.6 4.1 4.5 3.3 4.2 3.3 5.7 3.6 4.9 5.2
##  [955] 6.0 3.5 4.6 4.6 4.8 5.8 3.4 4.9 4.5 3.3 4.2 5.7 4.6 5.4 5.8 5.1 4.5 6.0
##  [973] 4.3 4.7 5.4 4.4 5.3 5.2 3.4 5.1 2.7 4.3 4.3 5.3 4.8 5.6 4.0 5.2 3.7 4.9
##  [991] 5.9 5.0 5.3 3.4 2.9 6.8 4.7 2.8 6.2 3.9
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
## 2.6975 6.2025
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
##    [1] 5.3 5.7 4.4 5.0 5.2 4.7 5.6 6.3 5.5 4.7 2.7 3.8 5.7 5.0 4.2 5.4 2.4 3.7
##   [19] 5.6 3.9 5.2 5.6 5.7 4.1 4.0 6.3 6.1 4.6 4.7 6.3 3.7 3.8 4.4 4.9 4.5 4.1
##   [37] 5.1 5.8 3.1 3.5 4.7 4.3 3.7 5.1 3.8 4.9 4.0 4.4 4.4 3.9 6.2 4.7 6.0 4.1
##   [55] 5.9 3.1 3.3 4.1 4.8 5.1 3.0 4.0 4.9 5.2 4.5 3.2 3.7 4.3 4.3 6.1 3.0 4.5
##   [73] 5.7 4.5 3.3 6.6 4.9 4.0 5.0 0.9 6.0 4.2 3.3 5.7 5.2 4.1 5.5 4.7 2.5 3.4
##   [91] 5.7 5.3 3.1 6.1 3.9 4.6 5.9 4.5 2.9 4.4 3.6 5.6 5.5 4.1 5.1 5.1 5.0 5.3
##  [109] 4.9 4.2 4.7 5.6 5.2 4.2 3.8 4.8 5.5 4.9 5.3 3.1 4.5 4.3 5.2 4.4 5.0 4.7
##  [127] 3.9 4.3 4.7 5.3 5.1 3.9 5.3 4.8 5.2 4.4 4.7 6.1 4.5 4.9 5.1 4.2 4.2 3.0
##  [145] 4.2 4.2 3.3 4.3 3.4 4.1 3.8 5.9 5.4 5.2 5.9 3.8 3.4 3.6 5.5 3.5 3.4 4.4
##  [163] 5.6 4.1 5.0 4.6 4.0 4.6 4.8 5.2 5.3 3.9 6.0 4.6 4.3 5.0 4.4 4.3 3.5 4.6
##  [181] 5.7 4.1 4.1 3.9 3.7 3.8 5.1 6.0 4.1 3.5 5.3 3.0 3.8 5.5 3.6 2.0 4.5 5.0
##  [199] 3.4 5.6 4.1 5.9 4.5 3.8 5.8 3.2 2.6 4.6 4.4 5.1 5.0 6.1 5.3 3.4 5.0 4.4
##  [217] 4.5 3.9 3.5 5.1 4.1 4.3 3.6 5.4 5.7 5.3 4.8 4.5 4.7 4.5 5.3 5.5 4.2 4.5
##  [235] 6.1 3.0 3.9 3.4 3.5 5.5 3.7 5.2 3.5 3.1 4.1 5.0 4.8 3.9 4.3 3.2 6.5 4.3
##  [253] 3.4 5.3 5.6 5.6 4.1 4.3 3.2 3.3 5.1 6.4 4.2 5.4 4.4 5.3 4.4 5.5 4.6 3.6
##  [271] 5.3 5.1 3.2 2.3 5.2 5.4 3.9 3.8 5.0 2.9 3.1 4.5 3.0 4.4 3.8 5.1 4.9 3.8
##  [289] 4.4 4.7 3.1 3.4 4.1 4.9 5.6 4.4 5.0 4.2 5.4 3.7 5.4 5.9 6.4 3.9 4.1 5.2
##  [307] 2.6 1.8 4.7 4.9 6.1 5.3 4.5 4.9 5.2 5.3 3.1 4.8 3.0 4.1 2.9 5.8 4.4 5.7
##  [325] 3.4 4.1 3.5 3.8 4.8 3.1 3.2 4.0 6.2 5.9 3.9 3.4 3.9 4.4 4.2 3.6 4.5 4.8
##  [343] 2.5 5.3 4.9 4.2 4.1 4.8 4.2 2.8 4.8 5.0 6.2 5.3 5.8 4.7 2.9 4.2 4.4 4.0
##  [361] 4.9 4.2 3.4 6.8 5.5 3.4 4.7 6.3 3.2 5.5 4.9 5.0 4.9 4.4 3.5 4.9 4.4 4.4
##  [379] 4.6 3.6 4.3 5.0 3.6 4.1 6.1 3.0 4.5 5.2 4.0 5.4 5.9 5.6 3.4 5.0 5.4 4.6
##  [397] 4.2 4.5 6.5 3.5 5.2 3.9 4.1 4.3 3.6 5.2 4.4 4.9 4.6 4.7 3.8 4.8 4.2 5.2
##  [415] 4.6 3.9 3.5 3.6 4.4 5.2 4.6 4.6 4.7 2.5 2.7 4.8 4.1 4.6 5.3 4.3 3.4 6.6
##  [433] 4.4 4.2 3.3 4.7 4.6 5.2 4.2 3.5 4.3 5.4 3.7 5.6 4.3 4.5 3.6 5.3 5.1 5.7
##  [451] 3.3 4.0 5.2 4.5 3.7 3.0 4.7 5.0 4.4 4.8 4.0 3.8 5.0 4.1 4.6 5.6 4.4 4.6
##  [469] 4.5 4.2 3.1 6.1 4.8 4.9 4.8 2.3 5.4 4.4 4.1 4.4 4.4 6.6 4.5 2.7 4.0 3.4
##  [487] 5.5 4.8 4.2 5.6 4.8 5.6 4.0 5.4 4.1 3.2 4.4 4.6 5.4 5.7 6.3 3.5 4.5 5.9
##  [505] 6.5 5.3 3.6 4.2 5.8 4.8 3.6 3.7 4.6 4.4 4.3 2.6 2.9 3.8 3.5 4.5 4.4 5.4
##  [523] 4.6 3.9 5.2 3.1 6.0 4.9 3.2 3.8 4.4 3.8 2.9 3.6 3.1 3.9 5.7 4.5 4.3 4.9
##  [541] 5.1 6.3 5.2 5.7 4.3 3.9 3.5 4.1 4.5 4.7 3.4 4.1 4.3 6.6 3.8 2.9 4.8 4.3
##  [559] 4.3 5.4 2.8 2.9 6.3 3.9 4.6 3.8 2.7 4.5 5.0 5.1 3.5 5.9 4.2 4.4 3.6 5.0
##  [577] 5.2 4.6 4.6 5.2 4.5 5.4 3.5 3.1 3.4 3.4 4.6 2.4 3.2 4.9 4.2 5.3 4.6 4.5
##  [595] 4.6 4.7 4.0 4.2 5.6 3.3 2.9 6.8 4.5 6.0 4.0 4.8 4.3 3.3 3.4 4.4 5.6 4.0
##  [613] 5.1 4.3 4.5 4.4 3.9 4.9 4.5 3.2 3.5 5.3 4.0 4.3 6.6 5.0 3.9 5.8 3.3 5.0
##  [631] 5.0 3.6 4.1 3.4 4.3 4.3 4.1 5.3 4.2 4.6 4.4 3.3 5.8 4.8 5.6 4.9 4.8 5.3
##  [649] 2.3 3.8 3.4 4.6 4.9 6.8 4.8 5.4 4.8 5.2 4.0 4.9 4.8 2.5 5.0 4.3 3.7 4.0
##  [667] 4.2 5.6 4.5 4.1 6.6 4.7 4.0 4.6 5.5 5.5 2.1 3.6 2.8 4.7 5.2 5.6 5.3 3.6
##  [685] 3.9 3.1 5.6 5.9 4.0 5.1 3.9 4.0 5.3 3.7 6.1 3.8 3.9 3.7 5.1 5.8 4.6 3.4
##  [703] 2.9 3.9 5.4 5.2 4.7 4.2 6.1 4.9 7.2 3.7 3.9 4.6 4.4 5.6 4.2 4.8 5.6 4.8
##  [721] 4.6 3.5 4.6 4.2 5.0 6.1 4.5 3.9 4.0 5.3 4.7 5.7 4.1 4.6 5.8 4.4 5.6 5.0
##  [739] 4.2 4.8 4.2 4.1 4.6 4.9 5.0 4.2 4.4 6.5 6.5 5.0 6.4 3.3 4.8 6.6 3.3 4.7
##  [757] 3.5 3.7 5.4 4.0 4.2 4.2 3.7 3.9 3.6 4.3 5.8 5.3 4.8 4.7 2.7 4.7 4.7 3.9
##  [775] 5.9 3.9 4.7 4.6 5.9 3.2 4.6 4.4 4.4 5.6 3.9 5.3 4.6 3.7 5.0 3.9 4.7 5.0
##  [793] 5.4 5.1 4.8 4.1 4.3 6.1 5.1 4.1 4.8 4.8 2.9 4.1 5.5 6.9 3.6 4.3 6.2 4.0
##  [811] 4.9 3.0 4.0 7.3 5.2 6.2 5.1 4.3 2.5 4.9 4.9 3.7 5.3 3.4 4.2 4.6 5.0 5.3
##  [829] 4.2 5.1 3.6 4.4 4.2 3.8 4.7 4.3 4.2 3.6 2.7 2.8 4.4 6.4 7.2 4.0 4.6 5.1
##  [847] 3.6 5.1 4.2 5.6 4.1 5.9 4.5 7.4 4.5 4.5 5.7 5.4 6.9 4.9 4.0 4.9 4.6 3.8
##  [865] 5.5 3.1 4.1 5.3 6.1 4.5 4.7 5.3 5.4 3.0 4.6 5.4 5.1 4.7 6.3 4.4 5.7 4.8
##  [883] 3.2 2.9 5.3 4.3 4.7 3.5 6.0 5.1 4.9 3.7 4.7 4.9 4.5 4.2 5.1 3.1 5.8 4.8
##  [901] 7.4 3.3 5.1 5.2 4.1 5.0 4.8 5.3 4.9 4.5 5.7 5.5 3.5 5.0 2.8 5.4 4.1 5.0
##  [919] 4.5 3.9 3.8 4.0 4.3 5.9 4.2 5.7 4.3 4.0 4.6 4.9 4.9 4.7 5.6 4.5 5.5 5.9
##  [937] 4.8 4.7 5.7 5.0 4.6 3.0 5.2 3.6 5.4 5.0 4.7 3.6 3.3 5.4 4.5 2.7 5.8 5.1
##  [955] 4.6 5.2 3.5 4.1 5.3 5.8 5.0 5.1 5.6 5.7 4.9 5.6 2.6 4.7 4.1 4.7 4.7 4.1
##  [973] 4.9 6.0 4.6 3.3 4.0 4.4 4.9 4.1 5.7 5.4 2.3 4.1 4.5 3.7 2.9 4.8 3.7 3.4
##  [991] 3.4 4.2 4.8 4.2 4.2 4.0 4.6 3.1 5.8 5.2
## 
## $func.thetastar
## [1] 0.0363
## 
## $jack.boot.val
##  [1]  0.53871866  0.42390110  0.26581197  0.24057143  0.08795181 -0.04036939
##  [7] -0.20087977 -0.19874608 -0.38490566 -0.46064690
## 
## $jack.boot.se
## [1] 0.9677687
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
##    [1] 5.3 5.1 3.2 4.7 5.1 2.8 5.7 5.1 6.6 3.5 5.4 4.4 5.4 5.3 4.8 4.6 5.7 3.7
##   [19] 5.5 5.2 4.1 3.2 4.5 6.3 4.7 6.8 4.8 4.5 3.8 6.2 4.6 4.4 4.4 4.5 4.6 3.6
##   [37] 5.5 4.3 4.7 3.1 4.8 3.1 4.3 4.0 4.6 5.0 4.7 6.0 4.7 4.8 4.7 5.3 4.5 4.1
##   [55] 4.4 6.3 5.6 4.1 6.8 5.2 3.7 5.8 5.8 5.0 4.8 4.1 4.2 3.6 3.2 5.6 5.1 5.0
##   [73] 3.4 5.3 5.2 4.0 5.3 3.2 5.5 5.0 3.9 4.3 4.1 4.1 5.6 3.9 7.1 4.8 4.5 5.0
##   [91] 3.4 4.0 5.1 4.0 5.6 4.7 5.2 5.3 3.7 5.6 4.1 5.1 5.7 3.5 4.7 4.5 5.9 3.3
##  [109] 4.5 4.3 4.2 5.2 3.5 5.4 5.9 4.4 4.4 4.3 4.6 5.6 4.5 2.9 4.4 4.3 5.8 4.6
##  [127] 4.1 4.6 3.5 4.7 4.2 2.8 4.3 4.8 4.6 6.3 3.9 6.6 4.5 4.5 5.4 3.1 3.8 5.5
##  [145] 4.4 5.2 4.5 3.9 5.1 5.0 2.6 4.8 4.3 4.2 4.4 3.5 4.6 4.4 4.7 3.9 3.7 5.2
##  [163] 3.1 4.9 2.7 5.2 3.7 4.5 4.6 5.5 5.3 5.0 5.5 3.6 4.8 4.0 3.6 5.1 3.4 3.3
##  [181] 5.4 4.0 4.4 4.2 4.5 5.2 3.6 6.2 6.4 5.9 2.9 4.8 3.5 4.3 4.4 3.7 2.9 5.3
##  [199] 4.0 5.1 5.0 5.4 5.0 3.8 4.3 6.1 5.4 2.5 3.4 3.4 4.1 4.6 5.6 6.0 5.1 3.5
##  [217] 4.1 3.3 4.5 2.4 3.1 5.2 2.8 3.8 3.6 4.9 4.8 3.0 4.1 4.4 4.3 6.3 3.7 4.8
##  [235] 4.8 5.6 5.1 4.5 5.5 3.6 4.4 5.3 2.5 5.3 5.7 5.4 3.5 4.1 5.9 5.5 5.6 4.7
##  [253] 3.3 2.6 5.7 4.3 5.4 3.7 4.6 5.6 5.0 5.3 5.2 3.1 5.1 4.2 3.9 2.7 4.1 3.9
##  [271] 3.4 2.7 3.8 4.5 4.2 3.5 3.0 4.1 6.0 5.6 3.3 4.2 3.0 5.1 4.3 4.0 2.8 4.5
##  [289] 5.9 5.5 5.1 3.7 3.6 4.1 3.6 3.6 3.6 4.1 3.2 3.8 3.2 4.1 5.4 6.0 6.1 5.1
##  [307] 4.1 3.9 4.5 4.2 4.7 5.1 4.6 3.5 3.9 5.2 4.9 5.4 3.8 5.2 4.5 3.4 6.7 6.1
##  [325] 3.7 4.1 4.2 5.4 3.9 5.9 5.5 2.8 4.2 4.3 4.0 3.9 5.0 4.0 2.8 5.1 4.3 4.1
##  [343] 5.5 4.5 4.3 3.7 3.7 4.9 3.9 4.1 4.5 6.2 4.6 5.0 5.8 4.1 5.0 4.6 3.7 4.2
##  [361] 4.1 5.0 4.5 4.0 3.8 4.2 3.5 5.0 3.3 3.7 5.6 2.8 5.1 3.4 4.2 3.5 2.1 4.1
##  [379] 3.2 6.6 5.2 4.5 4.5 6.4 4.8 4.9 2.6 4.6 4.4 3.6 4.8 5.6 5.4 4.0 3.0 5.2
##  [397] 5.0 3.0 3.6 2.4 4.3 5.0 6.3 5.4 5.2 3.9 5.3 5.3 5.2 3.9 6.1 5.5 3.7 4.2
##  [415] 2.9 3.6 3.3 6.1 3.0 5.7 5.6 6.0 4.0 3.8 4.8 3.0 4.5 4.2 3.1 4.3 5.6 5.1
##  [433] 4.6 5.1 3.6 4.2 5.2 4.9 5.3 5.5 4.4 3.7 5.0 4.4 4.1 5.0 5.2 4.3 2.8 4.0
##  [451] 4.6 4.1 4.1 4.5 4.8 4.5 4.5 4.5 3.5 5.5 4.9 3.7 5.9 3.4 4.0 4.7 4.8 4.9
##  [469] 4.2 4.0 5.9 5.4 4.6 7.1 5.0 5.2 3.9 3.9 3.2 5.6 4.5 4.0 5.6 3.6 4.2 4.7
##  [487] 5.3 6.4 4.8 4.9 5.1 4.0 5.9 5.1 5.2 4.7 4.1 4.7 2.5 4.1 3.5 4.4 4.3 3.9
##  [505] 5.2 4.1 5.0 6.1 4.8 4.5 4.1 3.7 3.0 3.9 4.3 5.4 4.2 5.0 4.4 6.5 4.9 4.1
##  [523] 4.2 5.6 4.5 5.1 6.0 4.8 4.4 4.9 3.9 5.4 5.6 4.0 4.3 5.6 3.8 3.9 4.5 3.7
##  [541] 4.7 3.3 3.8 4.4 4.6 3.2 4.3 4.8 4.4 3.9 4.5 5.6 4.3 4.4 4.3 5.3 5.2 4.6
##  [559] 6.2 3.7 3.4 5.8 3.6 4.8 3.0 6.7 4.7 4.7 5.1 3.4 5.5 4.6 3.3 5.3 3.3 4.3
##  [577] 3.3 3.1 4.0 4.7 4.4 5.2 4.3 2.7 5.8 4.6 5.0 5.5 4.4 5.1 6.3 4.7 4.6 2.8
##  [595] 5.6 4.0 6.4 4.7 3.9 4.3 4.9 5.8 5.1 5.3 4.7 3.7 5.3 4.3 4.5 5.3 3.7 5.0
##  [613] 4.0 6.6 6.0 3.8 5.3 3.2 4.7 4.8 5.1 4.1 4.8 4.3 5.3 6.1 5.5 4.7 2.9 5.7
##  [631] 5.1 5.7 3.9 5.4 4.0 5.1 4.0 3.7 3.1 4.8 5.0 4.3 3.2 3.3 4.5 3.5 3.2 5.2
##  [649] 5.2 4.3 4.0 5.9 4.3 5.4 5.1 5.0 6.2 5.9 5.1 4.2 4.5 5.5 2.2 3.9 2.9 4.8
##  [667] 4.0 3.6 4.7 3.5 5.0 4.3 4.9 4.4 4.5 4.0 3.0 5.3 5.4 3.1 3.0 5.4 4.1 4.2
##  [685] 4.1 3.2 4.6 4.4 5.1 4.8 4.2 4.5 3.5 4.4 5.1 4.8 4.7 4.3 5.5 5.1 4.5 4.5
##  [703] 3.9 6.0 5.3 5.3 4.6 3.9 7.0 3.1 5.0 5.0 5.0 3.8 4.3 4.2 4.4 4.5 5.2 6.3
##  [721] 4.8 5.0 4.6 6.2 3.8 4.8 4.2 5.6 3.5 3.2 5.4 4.7 4.7 3.7 4.4 3.5 5.4 5.0
##  [739] 4.5 2.7 4.8 6.1 4.9 1.9 5.6 5.2 6.3 5.6 3.9 4.5 4.1 4.7 4.7 2.7 5.8 4.0
##  [757] 5.5 4.4 3.6 4.2 4.8 4.3 5.1 4.7 4.2 4.4 4.7 4.3 3.7 5.5 4.6 4.0 3.3 4.6
##  [775] 4.7 5.3 3.1 4.6 4.0 4.8 5.2 4.9 3.8 6.8 3.9 5.8 3.4 4.0 4.2 2.7 3.8 4.8
##  [793] 5.0 2.5 5.6 5.5 6.6 4.2 3.9 4.8 4.4 4.3 3.9 3.1 4.6 3.1 5.0 4.1 6.8 4.0
##  [811] 5.1 5.0 2.0 5.4 3.6 3.8 3.3 3.6 5.2 5.1 5.4 4.9 4.3 5.3 3.7 4.6 4.4 5.0
##  [829] 3.9 3.6 5.4 5.7 3.8 3.5 4.1 4.1 2.9 3.3 5.0 4.7 3.9 3.5 5.4 5.6 4.3 6.1
##  [847] 4.8 5.4 4.5 3.6 5.0 4.2 4.7 4.7 6.5 5.5 5.8 4.7 3.5 4.8 4.9 7.9 6.0 4.0
##  [865] 3.2 4.4 7.3 5.3 5.2 3.5 5.8 3.3 5.7 5.0 6.1 2.8 4.3 5.8 4.5 4.9 4.5 4.7
##  [883] 3.9 4.9 5.7 4.6 4.8 5.4 5.9 4.9 4.2 3.2 3.0 3.7 4.3 5.4 4.4 4.7 5.5 4.1
##  [901] 5.4 3.4 5.7 4.9 5.3 3.3 5.6 5.2 4.0 6.0 5.3 3.7 4.3 6.7 4.4 4.4 5.6 3.6
##  [919] 4.5 4.2 2.5 4.2 3.2 5.0 3.6 5.3 3.9 4.0 5.0 3.3 5.0 6.1 4.2 5.3 4.1 3.1
##  [937] 5.3 5.2 3.3 3.6 3.9 3.6 4.2 2.7 4.7 4.7 4.5 4.4 4.3 4.7 3.5 3.0 4.1 4.8
##  [955] 3.6 2.6 3.2 3.8 4.8 5.1 4.5 4.5 4.2 5.3 3.7 3.7 5.5 4.0 5.4 3.9 6.5 4.3
##  [973] 5.1 4.5 4.9 5.0 4.8 4.8 3.6 5.9 4.4 4.4 4.3 4.1 3.4 3.9 5.4 4.2 5.2 5.4
##  [991] 3.3 4.0 4.6 5.8 2.8 5.1 3.4 5.7 3.4 4.5
## 
## $func.thetastar
## 72% 
## 5.1 
## 
## $jack.boot.val
##  [1] 5.500 5.400 5.400 5.300 5.200 5.000 4.832 4.800 4.700 4.500
## 
## $jack.boot.se
## [1] 0.9805373
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
## [1] -0.297959
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
##    8.360803   14.060058 
##  ( 3.666885) ( 6.355292)
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
## [1] -0.0004897419  0.7472413865 -0.5483833245  0.6829583314 -0.4094210101
## [6] -0.1696467951
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
##    [1] -0.0210310203  0.8211144933 -0.6010071464 -1.2517772782 -1.1370255377
##    [6]  1.8866531296  0.7781950731 -1.1386359028 -0.5777771260  1.4901939449
##   [11]  0.3600027531  1.0018034638 -0.8685144641 -0.6702157268 -0.1374101440
##   [16] -0.9200095451  1.3208429720  1.0993999337  0.4631205971 -0.6432012352
##   [21] -0.1432216406 -0.1164766324 -0.1877174973 -0.9561006036  0.5405920342
##   [26]  0.8982628917  0.0068977353  0.5717090200 -0.7616653167 -0.3538424628
##   [31] -0.2688974512 -0.3767447923  0.0057873499 -0.3410978780  0.2607784863
##   [36]  0.6542924336 -0.0084637686 -0.2766289489 -0.6918328858 -0.1315722992
##   [41]  0.8972315126 -1.2669349851 -1.3348628112 -0.2035507978 -1.1872378977
##   [46] -0.0172806428 -0.3103274746 -1.0522841579 -0.3977661890 -0.4555244299
##   [51]  0.4821576709 -1.4404690805 -1.1663369651  0.4168205273  0.5712334097
##   [56] -0.7982354612 -0.2606318545 -1.8664164239  0.4637483779 -1.2924928532
##   [61]  0.3268323847 -0.2693403291 -0.1186566478 -1.1650674054  1.5655994648
##   [66] -0.4438592276  0.3552799282 -0.2848864252 -1.3134751592 -0.9147122546
##   [71]  0.7252663965  1.0282495300 -0.5689692907 -0.6692633258  0.6841278976
##   [76] -1.2481216248  0.2386543590 -0.3004950729  0.8257760238 -0.3865707004
##   [81] -0.4634322433 -0.1107650978  1.0347983041 -1.1878482579  1.3049041096
##   [86] -0.4725331891  1.8429925318  0.3022533833 -0.7629822016 -0.1864751216
##   [91] -0.6898702483  0.2632560243  0.8179298343  0.3513819786 -0.0032313175
##   [96]  0.9338389010  1.0709471948  0.0638086650  0.8608789549 -0.5947917948
##  [101]  0.6572498159  0.4029269890 -1.0950957504 -0.3785896630 -0.3073047119
##  [106]  1.3086461910 -1.0209393835 -0.3275314988 -0.4178642109 -1.1759674351
##  [111]  0.8233499021 -0.2960233881  0.0724792665  0.9521861026  0.0128509913
##  [116] -0.9137729925  0.7695726703 -1.8717880553 -0.6677096538 -1.6718242468
##  [121]  0.6376941350  0.2619766118 -0.6248410916  0.1558453659  0.5182863192
##  [126] -1.4051828182  0.7496041671 -1.2213055526  1.1600723039  0.4298752034
##  [131] -1.1139674181 -0.8149548764  1.1002302088  0.0882082799 -1.0016720601
##  [136]  1.3833236166  0.0004560372 -1.4900022051 -0.8224799034  0.1953106392
##  [141] -0.1596198827 -0.3524080501 -0.2852007560 -0.8956339866 -0.0600726855
##  [146] -0.8532586246 -0.2880002668  0.0224389361  1.1906603139 -0.8736663794
##  [151] -0.2322885222 -0.6437411916 -0.5576791602 -0.6049649498 -0.2829197806
##  [156] -1.1927981965  1.2784999309 -0.7807395600  0.1235552567 -2.2589956963
##  [161] -0.2960756371 -0.0515696495 -0.3278008982  0.1259061505 -0.2728310812
##  [166]  1.0763271054 -0.7517155224  0.8628604903  0.2732182548 -0.3571394838
##  [171]  0.0172943307 -0.5721937975 -0.0390330962 -0.2204386374 -0.0539443996
##  [176] -0.9863893755  0.0545035087 -0.5744513307 -0.6617721803  0.5177192428
##  [181] -0.8342439925  0.2689359093 -0.1940265308  0.2477573652 -2.1080627948
##  [186]  0.5104778647  0.8284390626 -0.1925726100 -0.4429712992 -0.8354518214
##  [191] -0.6205875868 -0.4518084493  1.3659876603  0.1873741507  0.2140423789
##  [196]  0.3818454334 -0.9343515668 -0.0404282230 -0.4555244299 -0.9225036565
##  [201] -0.2860038481  0.6290833407 -1.0192001889 -0.0232125116  0.0617444327
##  [206]  0.2009715943 -0.2794911558  0.2084227709  0.0326972388 -1.0644416437
##  [211] -0.0527214189 -0.1734899296 -0.4743646359 -0.0381745707 -0.1188171111
##  [216] -0.3388323113 -0.3330299796  0.3970749945 -0.8132952932 -0.1776575555
##  [221]  1.8437228614  0.8871766625 -0.3262911287 -1.5551788209 -0.6003295343
##  [226] -0.6380398593  0.5634811672 -0.4557257323 -0.0370185765 -0.3211487469
##  [231]  0.0846740794  1.4105547769 -0.2622648402 -0.4740109598 -0.0101122849
##  [236] -0.8160368509  0.1015147736 -1.4239456396 -0.5091310680 -0.4915013824
##  [241] -0.2965489206  0.1352572297 -0.4930556070 -0.8027332248 -0.7926825193
##  [246] -0.4197072685  0.5232746759 -0.2866265624  0.7003059558  0.7634402686
##  [251] -0.2501075955  0.0875295844  0.2302134164 -0.7956383384 -0.8620345601
##  [256]  0.3879437626 -1.2538314367  0.2286161087  0.4702311479  0.8155948715
##  [261]  0.4728937375 -0.2854045755 -0.7789109686  0.5872544895 -0.0258933668
##  [266] -0.2069778245 -0.2271443321  1.6121372605  0.3258295742 -0.2629486778
##  [271] -1.2694484879  0.1934089781 -0.5946605998  1.2691831364  0.0114882174
##  [276]  0.2594557932  0.3108147850 -0.1689607387  0.2174865413 -1.2677040093
##  [281] -1.1778578149 -1.9511889159 -0.4446461631 -0.5155930020  0.4266963670
##  [286]  0.0059043895  0.7070223866  0.6594922673 -0.3995436048 -0.4920219064
##  [291] -0.4648329386 -0.3585021506  0.1290751707 -0.7673817373  0.1214905474
##  [296] -0.8656013876 -0.0479573547  1.4239990612 -1.4623814616  0.8335325828
##  [301] -1.0655337739 -0.1856777729  1.1803741641 -0.2370545752 -0.5638328197
##  [306]  0.4847530735  1.1858363322 -1.6521282298 -0.1384303014  0.0174498179
##  [311] -0.5503096818  0.0299971750 -1.2497880089 -0.1220008957 -1.6237543753
##  [316]  0.1322316828  1.3750613513 -0.1796942588  0.3850356669 -0.0260811766
##  [321] -0.5539454182 -1.1840368098 -0.9075084509 -0.8881963192  0.0631987232
##  [326]  0.1921079325  0.1173345997 -0.2337375716  0.4189383216  0.0737194029
##  [331]  0.0847137803 -1.0275049379  1.0739494632  0.1509810360 -0.9812369154
##  [336] -1.1764877441 -0.7199361699 -0.0934866971 -0.2094502635 -0.4604404585
##  [341] -0.0745998307 -0.2532156660  0.1025967472 -0.4339222070 -1.1140877666
##  [346] -1.4892278922 -0.0953826781  0.4125345448 -0.2268961178 -1.1510485968
##  [351]  0.1144549826 -1.2045071864 -0.1061273005  1.7030621775 -0.0016623962
##  [356] -0.6875616704 -0.3986257284 -1.0631875789 -0.3361563662 -0.0079118535
##  [361] -0.2233593300 -0.0507547907 -1.0871822343  0.8489205305  0.6641194988
##  [366]  0.3848684635 -2.0887669493  0.4898374405 -0.1490562985 -0.1296949857
##  [371]  0.7194374404 -0.1890611276  0.9228724019  0.4364598229 -1.1322759621
##  [376] -0.4749391433  0.0300831712  1.3166885690  0.4354723533 -0.3000238205
##  [381]  0.5868169259  1.0054632178  0.0995938612 -0.4491516643  0.9516886924
##  [386]  0.5974701358  0.0505823924 -0.6402602035 -0.4558473664  1.1510636945
##  [391] -0.8391214074 -1.3237379559  0.3455652859 -0.4924686467  1.1341476122
##  [396] -0.2738212876 -1.0648862097  0.3131677027 -0.4535993621  0.3641023404
##  [401] -1.5174082641 -0.0328187855 -0.4520354607  0.5615312007 -0.6980830293
##  [406] -0.2335202316 -0.9393661109  1.1300709763 -0.4433256468 -0.7899318855
##  [411] -0.8563487172 -0.7924605991  0.9947256594 -0.2559758372  1.1751641469
##  [416]  1.2959821972  0.9082191388 -0.7672647196 -0.4470162318  0.7628368177
##  [421] -1.5538620415 -0.3235898247  0.6123717499  0.0707741610  0.3222441718
##  [426] -0.4707857898 -1.4471152295 -0.7006946551  0.5964781402 -0.8236110702
##  [431]  0.6294660855  0.3877901535 -0.0987379012 -0.6965477155 -0.2954728066
##  [436]  0.8078003944  1.2350947062 -1.8138314126 -0.0803793609 -0.0695938363
##  [441]  1.0692535297  0.1408992511 -0.5772494534 -1.0968454687 -0.3480991360
##  [446] -0.4362388904  0.9915582336 -0.2031292475  2.5518456811  1.2034845591
##  [451]  0.1917777641  0.7471144513 -0.7714650798 -0.3488107265 -0.6095419292
##  [456]  1.0042631355  1.7336342767 -0.1725573341 -0.0477855699 -0.9001478214
##  [461]  0.3233309124 -0.6636869885 -0.7956009302  0.2198537457  0.0374116676
##  [466] -1.1987587119 -0.4855969360 -0.1164766612 -1.3334040690 -0.3907759287
##  [471] -1.2104918733 -0.3239174641 -1.6297323185  1.8916968639 -0.1695662273
##  [476] -0.2214447868  0.0693573357 -0.9374060615  1.1118125409  0.7752982215
##  [481]  2.4912061264 -0.2072782024 -1.1737377614 -0.8870141264  0.0464755277
##  [486] -0.3857543886 -0.3331425914 -2.4930039119 -1.5393116324 -0.5770083585
##  [491] -1.6276318927 -0.3354740338  1.1727773608 -1.1259051583 -2.3440272635
##  [496] -0.1761340014 -0.3474341038 -2.4299943412  2.0457283823 -0.3977661890
##  [501]  1.2181393784 -0.4415933138 -0.4948836472  0.2764832222 -0.5263080567
##  [506] -0.3961986212 -0.2670479083  0.2587602725 -0.1113705567 -1.2367450353
##  [511] -0.9337772754 -0.5770968560  0.8024105149 -0.1322748133 -0.2740944741
##  [516]  0.1022015825 -0.5787970257  0.1052037744  0.2096386850 -0.2779826292
##  [521] -0.1266268711 -0.3508185510 -1.0838063004 -0.1687800170 -0.3377647671
##  [526]  1.0242351736 -0.5383139290 -0.3213325739  0.0928291732 -1.6276641279
##  [531] -0.2830470904 -0.7264954250 -0.0049502644  0.2319470686 -0.2248154352
##  [536] -0.0059214297 -0.2545371029 -0.5335594640 -1.4231826929 -0.3492932396
##  [541] -0.5867232248 -2.1985191207  0.8915971380  1.2866890893 -0.3329608872
##  [546] -1.2179826286 -0.2705402843 -1.1793624054  1.2826434848 -0.2378551815
##  [551] -0.1592181855 -0.3419543233  0.0575031292 -1.3114733950  0.4913054824
##  [556] -0.0412916937  0.0898189926  0.3605913761 -0.4067069061  0.2878308530
##  [561]  0.2322914746 -0.3872857981 -0.3253139178 -0.2511176484 -0.4650042745
##  [566] -0.3448389076 -0.0848159669 -1.1274508121  0.0014094746 -0.3010297957
##  [571]  0.1389014085 -0.5203993846  0.1765289729 -0.2471184111 -0.9438588168
##  [576]  0.5406597289 -0.3309025650 -0.1993263622 -0.0913080619 -0.6765493182
##  [581]  0.3127212604 -0.7017441856 -1.2549525948  0.2629016900 -0.4406435817
##  [586]  0.0546215920 -0.2638061575 -0.4232577633 -0.4751366219  0.3380936346
##  [591]  0.3099979975 -0.8125734772 -0.5314150373  0.9927781651 -1.3107036038
##  [596] -0.9314257489  0.4071507765 -0.9232193345 -0.3440202989 -1.2583413872
##  [601] -0.0294695201 -1.5153728714 -0.6174356611 -0.3525861591 -0.2390854334
##  [606]  0.5904577418 -0.1025223687 -0.1072149796  1.3895111245  0.2817465371
##  [611] -1.0801463963 -0.0895978719  0.9848717846  0.2368593503  1.4392814357
##  [616]  0.5267617832  0.1217524535 -0.8103285938 -0.4277073310 -0.5409913686
##  [621]  0.9710593107 -0.3055087814  0.6578540782 -0.1052838397  0.3620481204
##  [626]  1.5255933480  0.8812763505 -1.3465690663 -1.1951651695  0.1480707976
##  [631] -1.7645562692 -1.9414344751  0.0695031787  1.1552949553 -1.1518548611
##  [636] -0.3406585666 -0.0594063093 -1.2080464756  1.3728059709  0.9491958537
##  [641] -0.3020708755  1.8182175441  0.5175678228 -0.1620354189 -1.5157976061
##  [646] -1.3476887278 -0.0289352105  0.4553342917  0.6077791238  0.2509544246
##  [651]  0.7579461790  1.3957810157 -1.8663780268  1.2951330484 -0.0512274770
##  [656] -0.0748672887  0.2602259921  0.1579783893 -0.0195438581 -0.3861737553
##  [661] -1.2868594345 -0.2280935852  0.1867642629  0.7903609953 -0.6799639824
##  [666]  1.0709287226 -1.5756596710  0.6017605548  0.3096311092 -1.6059481582
##  [671] -0.3668603317 -0.0517436655  0.8082474833 -0.3202103405  0.6004059057
##  [676] -0.4983983323  1.3129742195  0.8440509351 -0.4528515557 -0.3826658085
##  [681]  0.9737588262  0.3619593440 -0.8039209265  0.8616883758  0.1532485596
##  [686] -0.3697339932  0.2795472960  0.0770121461  0.2126266679  0.0410570263
##  [691] -0.3349171587 -0.7804651369  0.5076441271 -0.4933370331 -1.3857001216
##  [696]  1.2011562962 -0.6657793104 -0.8995787936  0.7037370964  0.5679701020
##  [701] -1.6497677061  1.6886516274  0.1900963539 -1.5960073205  0.0631987232
##  [706] -0.4261621651 -1.1008124518 -0.7924373147  0.4816803871 -1.2317353874
##  [711] -0.1178910364  0.8502644399  0.8999906536 -0.9559076955  1.7897609911
##  [716] -0.2226336638  0.8893039130  1.0616789477  0.1652110612  0.9991094406
##  [721]  0.8220583546  0.2207938943 -0.2357393678 -2.6114316017  0.3249609786
##  [726]  1.4231949410 -0.2051760056 -1.6682160821 -0.1167286691 -0.3606915674
##  [731] -0.4531623767 -1.6778616601 -0.6512383568 -1.3947580062 -0.0947075462
##  [736] -1.0127074815  0.0545377952 -0.6125144175 -0.5788333218 -0.0543789803
##  [741]  0.9045668421  0.5609911772 -0.8842922877  0.4320167427  0.0811743446
##  [746] -0.7557156878 -0.6695769147  1.1235743072 -0.3103274746  1.4115437649
##  [751] -1.6452329569  0.1185214369 -0.2618577355 -0.4439344239 -0.1889664539
##  [756] -0.5182699294 -0.2124082515  0.5026382636 -0.5638211003 -0.1304700091
##  [761]  0.2367242596  1.8829356162 -0.6178722312  0.8693819582  0.4331638136
##  [766] -1.3066416217 -0.5471960832  1.1596532492 -0.8456462996 -0.9687664006
##  [771] -0.3188491883  0.0582870346 -0.2237838890  0.3583231925  0.9850768240
##  [776] -1.3110276734 -0.1402470375  0.9859992985 -0.1542244587 -1.2823517008
##  [781]  0.9341057457  0.0354026923 -1.2477078011 -1.1229198714 -0.3263776699
##  [786]  1.0101366136 -1.3525813000 -0.9733714076  0.8152643859  0.2056687519
##  [791] -0.2415278855 -0.8106210786 -1.3349888240  1.0990667913 -0.8481533197
##  [796] -1.2256944981  0.5299263998 -0.0201101958  0.0416662370  0.3701337303
##  [801] -0.1594655589  0.0977699365 -1.1326535412 -0.2889750779  0.8079286527
##  [806]  0.1612402232  0.3575959563 -0.1062504625 -0.7291513375  0.9935672642
##  [811] -0.2110950602  0.1065301805  0.6470673841  0.3766772247  0.5509285301
##  [816] -1.4708243857  0.2174283701  0.8652239599  0.2955510451  1.5936771067
##  [821] -1.1466845889  0.2437853233 -0.4758801282 -0.5001828179  1.2171563769
##  [826]  2.3946200025 -0.8669325997  0.1472085211  1.2116288501  0.4002496888
##  [831] -0.2062015108 -0.6491997548  0.7019558560  0.0706475971  0.3213865014
##  [836] -0.8184684491 -0.5759154701  0.7107750143  0.9756441908 -0.4817247266
##  [841] -0.6954655884 -0.0385649871  0.1877521144  0.0827227968 -1.8034395060
##  [846]  0.0841852643 -0.3624087432  0.0862072169 -0.4495385478 -0.7681423613
##  [851] -0.4098190856 -0.3656941245 -0.2254734193  1.5074143361 -0.1611760777
##  [856] -0.3475538150 -1.1024241158 -0.3133880207 -0.6870613429  0.7875145759
##  [861] -1.1805108954 -0.9418238243  0.2633540513 -1.5044345243 -0.3944205128
##  [866]  0.4275176918 -0.7410437196  0.7854604175 -0.3762365549  0.9907041268
##  [871] -0.0089815353 -0.0337016737 -0.6875501322  0.4858597145 -0.4948031113
##  [876] -0.0939071230  0.5744138556 -0.8005269256 -0.3522189230 -0.4422094371
##  [881]  0.7003059558  0.2995669957 -0.3247549828 -1.3357149807 -0.2960756371
##  [886]  0.1192792705 -1.7532077998 -0.8900715405  0.0173396140 -1.6558394951
##  [891]  0.9720509258 -0.7095960697 -0.2451409116 -0.3363330939 -0.4650042745
##  [896] -0.5322507216 -1.6725562678  0.2131936509  0.0037228885 -1.3801025746
##  [901] -1.7218127548 -1.4150682521  0.2375060975 -1.8940504491 -1.0142997279
##  [906] -0.3602670598  1.2190958826 -0.1973882194 -0.5222960003 -2.1991694180
##  [911] -2.0230434565 -1.3424173827  0.9257704884  0.2045649152 -0.4244758892
##  [916]  0.8536546882 -1.1028207307 -0.4208099139 -0.9337772754 -0.8303861217
##  [921]  0.1198690687 -2.5047499670 -0.1013677609  0.5191244052 -0.3421843611
##  [926]  0.3405849793 -0.3801626767 -0.1506956935 -0.1284367893 -0.5473703973
##  [931] -1.0541308685  0.3967237039 -0.4341411019 -0.5886198570 -0.5093198367
##  [936] -1.2751052450  1.5053464414 -1.0645611957 -0.9980034171  2.4668718232
##  [941]  0.6456944774 -0.4326209401 -0.3591038611 -0.7498104085  1.7121165480
##  [946] -0.7369154808  0.5564431368  0.8307038774 -0.9668207166  0.0842489971
##  [951]  0.6277040452 -0.4541504608  0.3642730979 -0.0260185684 -1.0274290313
##  [956] -1.6239310336 -0.1119671948 -0.7118211278  0.3820679890 -0.2045558506
##  [961] -0.1043734698 -0.6156572482  1.2920688929 -0.5496778984 -1.4022158548
##  [966] -0.1405130043 -0.2239921163 -1.8050288948 -0.3569918844  0.4890817947
##  [971] -0.5458705433 -1.1601201333 -0.3648202337 -0.4228869318 -0.8694776123
##  [976] -1.1466845889  0.3579915600  0.7742417552 -0.5540344417  0.0375415270
##  [981]  0.2464786272 -0.2471282859  1.0340516578  0.3176911986 -0.4915013824
##  [986]  1.4915042815 -1.6534101645 -0.1050267313  1.0414081862 -0.7588294286
##  [991]  1.9242431321 -0.3751530076 -0.0781568503 -0.5983349458 -0.0833975490
##  [996]  0.0473013031  1.1353153517 -0.0054577769 -0.2210476508 -0.5024604768
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

## Warning in densfun(x, parm[1], parm[2], ...): NaNs produced
```

```r
fit2
```

```
##       mean          sd    
##   0.59464267   0.18143140 
##  (0.05737365) (0.04056349)
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
## [1] -0.921337996  0.009119563  0.059434194  0.069114009 -0.264210061
## [6] -1.775097690
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
## [1] -0.0298
```

```r
sd(results2$thetastar)  #the standard deviation of the theta stars is the SE of the statistic (in this case, the mean)
```

```
## [1] 0.9149148
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
## t1*      4.5 0.05265265   0.9177882
```

One of the main advantages to the 'boot' package over the 'bootstrap' package is the nicer formatting of the output.

Going back to our original code, lets see how we could reproduce all of these numbers:


```r
table(sample(x,size=length(x),replace=T))
```

```
## 
## 2 3 5 6 7 8 9 
## 1 1 1 3 1 2 1
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
## [1] 0.0212
```

```r
se.boot
```

```
## [1] 0.9235494
```

Why do our numbers not agree exactly with those of the boot package? This is because our estimates of bias and standard error are just estimates, and they carry with them their own uncertainties. That is one of the reasons we might bother doing jackknife-after-bootstrap.

The 'boot' package has a LOT of functionality. If we have time, we will come back to some of these more complex functions later in the semester as we cover topics like regression and glm.

