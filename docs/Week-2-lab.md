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
## 1 2 4 6 7 8 
## 2 1 1 4 1 1
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
## [1] -0.0287
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
## [1] 2.639418
```

```r
UL.boot
```

```
## [1] 6.303182
```

Method #2: Simply take the quantiles of the bootstrap statistics


```r
quantile(xmeans,c(0.025,0.975))
```

```
##   2.5%  97.5% 
## 2.6975 6.2000
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
##    [1] 3.8 3.2 4.1 4.3 3.8 6.2 5.3 4.6 4.6 4.2 4.8 6.4 4.0 4.7 4.6 5.1 6.7 5.4
##   [19] 4.9 4.6 4.0 5.6 5.2 3.7 3.5 5.5 5.2 5.9 4.7 4.0 3.9 4.8 2.8 4.3 5.7 3.9
##   [37] 3.3 4.3 3.7 5.7 4.0 5.5 3.6 5.8 5.2 4.2 3.4 3.5 3.6 4.3 5.8 3.9 3.0 4.2
##   [55] 3.7 6.4 4.7 3.9 4.3 3.7 4.0 4.8 2.6 3.8 4.4 4.7 4.5 4.5 4.3 4.5 3.4 3.6
##   [73] 4.0 6.3 3.8 3.5 4.7 4.6 4.2 3.5 4.0 4.6 5.0 5.5 5.1 2.8 3.7 7.0 4.1 6.3
##   [91] 3.0 5.1 4.5 5.6 4.0 5.8 3.8 5.4 2.2 4.9 3.7 4.8 3.7 4.8 4.2 4.1 4.1 4.9
##  [109] 3.8 4.2 2.0 4.9 6.2 5.4 4.6 5.1 4.0 3.3 4.8 5.5 4.0 4.8 3.0 4.2 5.7 5.7
##  [127] 4.1 5.2 3.7 5.6 4.2 4.5 3.5 3.1 4.2 5.6 5.2 4.3 4.3 6.2 3.7 4.4 3.2 5.1
##  [145] 6.0 4.8 4.6 4.3 6.3 4.7 5.9 7.1 3.4 4.1 3.8 4.7 4.9 6.2 4.0 2.8 4.4 2.6
##  [163] 5.5 4.5 5.1 2.9 5.8 2.9 4.5 4.4 4.8 4.1 3.9 5.5 5.0 3.8 3.8 4.0 3.9 4.3
##  [181] 4.9 4.3 3.0 3.2 4.8 4.2 3.6 5.2 3.5 5.3 3.9 3.4 6.0 5.0 4.3 3.5 4.0 4.1
##  [199] 3.7 5.0 3.4 4.1 3.8 3.7 4.3 3.3 5.2 3.9 3.2 4.1 2.8 4.5 4.7 5.2 3.4 4.5
##  [217] 3.5 3.8 4.1 3.4 6.2 4.5 4.5 4.8 5.5 5.4 4.6 4.5 5.1 5.9 4.0 5.9 4.0 6.1
##  [235] 3.9 4.6 6.2 3.7 5.3 2.0 5.0 3.9 3.6 3.2 6.4 3.7 4.3 5.1 3.9 4.2 4.6 3.6
##  [253] 3.9 5.1 3.5 3.5 3.3 4.1 5.9 5.4 3.6 3.0 5.4 4.3 3.0 4.9 3.4 4.4 5.4 3.9
##  [271] 4.9 4.8 6.7 4.3 5.1 4.1 4.8 6.1 4.8 3.4 6.6 4.2 3.9 4.1 6.0 5.4 5.3 3.8
##  [289] 4.5 5.0 4.7 4.0 5.3 5.8 4.2 4.3 4.9 2.3 4.1 5.6 4.1 3.5 6.4 5.2 3.8 5.4
##  [307] 3.5 4.8 4.2 5.9 6.4 6.0 4.9 3.6 4.3 4.9 4.3 5.4 4.5 4.0 3.8 4.7 5.7 5.1
##  [325] 4.5 5.0 5.4 4.6 6.1 3.9 4.3 5.3 5.5 3.3 2.1 4.8 5.3 3.8 5.9 3.4 3.4 5.1
##  [343] 3.5 4.9 3.4 2.5 5.4 3.5 5.0 5.3 4.2 3.2 4.0 3.7 5.9 2.9 2.9 4.2 4.3 2.6
##  [361] 4.3 5.1 4.9 3.8 4.3 4.9 6.6 4.9 5.2 4.5 4.5 5.1 6.1 4.3 5.6 4.4 4.8 5.2
##  [379] 4.0 2.3 5.1 4.0 4.5 4.8 5.0 5.3 4.8 4.5 5.4 4.7 5.2 4.2 4.1 5.5 5.6 4.9
##  [397] 2.1 4.5 2.3 3.8 3.9 5.2 5.4 4.8 2.8 3.0 4.3 4.5 3.0 3.6 2.9 5.5 4.5 4.3
##  [415] 3.3 4.2 4.0 3.8 6.4 4.0 4.2 4.7 3.1 3.3 3.6 6.7 5.2 2.9 3.3 5.0 5.4 4.7
##  [433] 5.9 4.5 5.8 4.9 4.7 3.4 4.1 5.8 3.6 3.8 6.1 5.2 3.8 5.0 5.8 4.9 3.6 3.9
##  [451] 4.1 5.1 3.6 5.5 5.1 4.8 4.5 4.1 4.0 4.3 4.6 3.5 3.5 5.0 3.2 4.9 5.8 5.9
##  [469] 5.3 4.3 7.2 4.1 5.7 3.8 5.0 4.3 4.3 6.8 5.7 3.8 3.3 4.1 4.4 3.7 3.8 4.2
##  [487] 3.9 5.8 4.8 3.7 4.7 3.8 4.7 6.7 4.9 3.9 4.0 4.3 4.9 4.6 3.5 3.7 3.2 5.2
##  [505] 4.0 3.8 4.5 5.4 4.9 5.8 3.4 3.4 4.4 5.3 4.3 4.9 5.2 3.3 5.8 4.4 4.6 5.6
##  [523] 5.0 3.3 4.9 4.8 4.5 3.5 5.1 4.8 4.7 2.9 5.2 4.3 6.0 4.1 4.4 5.4 6.5 5.8
##  [541] 4.9 5.1 7.0 5.4 4.1 3.2 5.0 3.3 5.1 4.7 6.0 4.4 5.0 4.4 5.7 5.4 5.2 6.7
##  [559] 4.2 2.2 3.0 6.3 5.2 4.2 3.9 5.8 5.3 4.2 4.2 4.8 5.5 4.2 4.3 6.0 6.0 4.2
##  [577] 4.1 4.1 5.8 5.6 6.2 5.1 4.9 3.8 5.0 6.0 5.0 3.3 3.9 3.4 4.6 3.6 4.3 3.5
##  [595] 4.4 5.3 4.2 4.0 3.6 4.8 3.9 4.2 4.3 5.6 3.9 3.8 5.5 4.6 5.1 5.4 4.6 5.9
##  [613] 2.8 5.9 5.5 5.2 5.9 4.1 4.7 1.9 5.3 4.4 5.6 3.5 3.4 4.6 5.1 4.2 5.0 4.7
##  [631] 3.6 6.5 4.0 5.2 5.7 4.0 4.6 3.9 4.9 4.4 3.8 3.8 5.7 4.4 4.4 4.6 4.4 5.5
##  [649] 4.2 3.9 4.8 4.0 3.8 6.9 5.0 4.5 4.6 4.8 4.8 4.4 5.1 3.9 5.3 6.5 5.1 4.2
##  [667] 5.1 6.0 5.3 4.1 6.1 4.6 5.7 4.8 6.3 4.7 4.9 4.4 4.0 3.8 4.4 4.3 4.0 4.1
##  [685] 4.3 5.1 6.5 4.5 4.6 4.1 4.9 3.5 5.3 4.6 4.9 3.7 4.8 4.5 5.1 5.3 3.5 5.5
##  [703] 4.5 5.5 3.3 4.6 5.2 5.6 6.6 4.3 4.3 6.3 6.0 4.6 6.8 4.2 2.5 4.7 5.1 3.9
##  [721] 5.5 3.7 4.9 4.4 4.7 4.0 5.4 5.6 4.0 4.5 4.0 4.5 4.9 4.8 5.3 4.0 4.1 3.9
##  [739] 4.0 4.5 4.0 5.2 5.2 3.9 4.0 4.2 3.9 4.9 4.8 3.3 4.2 5.6 3.2 4.6 4.5 4.0
##  [757] 3.9 4.5 5.6 3.6 5.7 5.2 4.5 4.7 5.2 4.5 4.3 4.6 4.1 4.6 4.0 4.0 3.9 4.4
##  [775] 4.2 5.6 4.6 4.8 4.5 3.9 3.2 4.7 3.8 4.6 3.2 3.7 3.5 4.0 4.1 4.0 4.2 4.3
##  [793] 4.6 5.4 5.8 4.3 5.3 4.6 4.0 3.1 5.0 5.3 3.7 4.4 4.2 4.7 3.2 5.4 4.6 3.8
##  [811] 5.5 6.3 4.3 3.6 5.2 5.5 5.3 4.9 4.4 4.3 4.6 3.7 6.1 4.4 5.3 5.3 4.6 3.3
##  [829] 3.7 7.4 5.8 4.2 4.1 4.9 5.4 5.8 5.8 4.5 3.3 5.9 3.6 4.1 4.7 4.1 3.7 4.0
##  [847] 3.3 3.6 2.9 4.2 4.6 5.6 4.8 4.0 3.3 5.5 3.5 5.8 4.3 4.7 4.9 4.7 2.8 4.0
##  [865] 5.4 5.0 3.2 4.2 3.4 4.8 5.0 5.4 2.1 3.1 3.9 3.4 5.3 5.2 5.3 4.0 5.1 4.1
##  [883] 4.1 5.2 5.5 5.5 5.1 4.5 4.4 4.5 5.4 6.3 4.3 3.7 4.5 5.9 4.3 5.4 4.9 5.9
##  [901] 4.5 5.2 5.3 4.3 4.3 4.5 3.7 3.9 3.4 4.1 5.0 5.0 3.7 6.0 2.9 4.5 5.7 4.0
##  [919] 5.0 4.0 4.8 5.4 4.4 4.4 3.5 5.7 4.0 5.0 3.8 5.0 4.2 3.8 4.9 5.2 4.8 5.0
##  [937] 5.5 5.2 3.6 3.9 5.3 4.8 5.3 3.6 5.6 5.1 5.1 3.7 4.8 3.5 3.7 5.1 4.5 5.5
##  [955] 3.6 4.2 3.7 2.9 4.2 4.2 3.9 4.8 3.8 4.0 5.1 3.6 5.5 3.3 4.6 4.8 4.5 6.0
##  [973] 5.2 4.1 4.1 5.6 2.4 5.4 3.5 4.6 5.3 3.9 4.0 5.0 5.2 5.4 5.3 3.3 5.0 5.3
##  [991] 4.8 5.4 3.9 3.3 2.6 5.0 4.3 4.1 2.9 4.2
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
## 2.8975 6.4000
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
##    [1] 3.6 3.0 4.2 5.1 3.9 5.9 3.5 5.0 4.9 4.0 5.2 4.1 4.8 4.5 5.4 3.1 4.4 4.3
##   [19] 5.2 3.8 4.4 4.4 5.6 4.7 4.3 3.4 3.9 4.2 4.7 3.5 4.5 3.9 3.1 5.0 3.3 3.6
##   [37] 5.1 5.2 5.8 4.6 3.3 3.7 4.1 5.1 5.3 5.3 3.9 6.0 4.6 4.5 4.5 5.4 4.9 5.7
##   [55] 2.4 2.1 4.5 4.1 1.9 5.9 4.8 4.6 4.3 4.1 4.0 5.0 3.2 4.9 4.6 3.6 4.3 3.9
##   [73] 3.7 5.5 3.2 5.8 4.6 4.1 5.1 3.1 4.4 6.1 4.1 3.9 5.1 3.5 3.5 5.1 5.1 4.2
##   [91] 3.9 4.8 5.6 4.3 4.1 5.2 4.9 5.9 5.1 3.2 5.4 5.3 3.5 6.4 5.1 4.6 3.6 3.7
##  [109] 3.5 4.1 2.9 5.3 4.8 4.7 3.8 5.2 5.4 6.1 4.8 5.9 3.9 3.5 3.7 5.5 3.3 6.0
##  [127] 3.2 3.0 5.0 5.0 3.6 5.4 5.8 4.6 4.6 4.6 4.5 4.9 5.6 4.8 3.8 6.2 4.5 3.9
##  [145] 3.8 3.3 5.7 4.4 4.0 6.2 5.9 2.4 3.9 3.6 5.9 4.8 5.8 5.6 3.6 3.9 4.6 6.6
##  [163] 5.5 3.6 4.0 4.8 3.9 4.5 6.2 4.7 5.1 4.1 3.4 4.7 4.0 4.9 5.2 3.7 5.7 4.6
##  [181] 5.5 3.6 4.7 5.2 4.0 3.3 3.4 4.7 3.9 4.7 4.0 2.7 4.7 4.2 4.3 6.2 2.9 3.6
##  [199] 4.3 6.5 4.6 4.7 3.9 3.6 5.9 4.5 4.7 4.8 4.1 5.1 4.0 4.0 4.2 5.1 6.2 3.4
##  [217] 4.0 2.7 3.9 3.5 4.5 4.2 4.2 3.6 3.5 4.3 4.0 5.1 4.0 3.5 5.2 4.5 4.4 5.4
##  [235] 5.2 5.3 5.2 3.8 3.5 3.5 4.7 4.0 3.8 5.1 4.1 4.2 3.8 4.9 4.3 5.3 4.8 5.4
##  [253] 4.3 4.8 4.6 4.6 3.6 3.6 5.9 1.9 3.5 5.0 2.4 5.6 4.1 4.0 3.2 5.2 6.0 6.0
##  [271] 5.0 4.7 4.9 3.3 3.6 5.0 3.7 3.9 4.1 4.7 6.8 4.4 5.7 4.9 2.9 5.7 5.0 3.6
##  [289] 5.8 5.8 4.6 5.1 5.9 5.2 6.7 6.8 6.4 4.0 4.7 4.5 4.0 3.3 4.6 5.8 5.3 5.3
##  [307] 4.4 2.9 5.6 4.8 3.0 4.6 3.5 5.1 4.6 3.3 5.6 4.5 4.2 3.5 3.5 4.1 5.1 5.8
##  [325] 5.0 5.1 3.3 5.6 3.6 4.9 4.1 3.3 4.4 4.6 2.9 4.8 4.6 5.0 3.9 4.7 4.6 2.6
##  [343] 4.2 6.1 3.7 4.6 5.8 5.1 5.5 3.2 5.6 2.6 5.9 4.6 3.4 4.6 4.5 6.2 6.4 4.9
##  [361] 5.1 6.2 3.9 5.1 4.0 6.5 3.9 4.6 3.8 3.3 3.7 5.7 3.1 6.1 3.8 5.1 4.9 5.3
##  [379] 4.2 5.5 5.3 3.6 3.9 4.6 3.9 4.1 6.6 4.5 4.0 4.4 6.5 4.9 3.3 4.0 4.4 3.8
##  [397] 4.4 5.1 3.2 3.9 3.3 4.2 3.0 2.6 4.7 3.8 5.2 4.5 3.7 4.3 4.1 4.2 4.0 5.1
##  [415] 3.6 4.8 5.7 4.1 4.5 5.1 4.9 7.1 5.0 3.5 3.6 7.0 5.5 5.2 3.2 3.9 3.3 4.8
##  [433] 4.1 5.6 5.5 5.0 3.3 3.3 5.1 4.6 4.3 4.9 3.8 3.8 5.0 4.1 4.6 5.2 5.3 4.6
##  [451] 5.0 4.7 5.4 4.6 3.6 4.6 6.1 3.4 3.7 4.1 4.4 5.9 4.3 4.6 6.4 3.4 4.4 5.5
##  [469] 4.1 3.8 3.7 5.2 4.0 4.4 4.8 4.4 5.9 4.0 4.4 4.8 5.0 3.8 4.6 6.3 3.3 4.4
##  [487] 4.1 5.9 3.5 4.0 4.6 4.9 3.7 3.4 4.4 2.4 3.9 6.7 4.7 4.1 4.4 4.7 5.1 2.8
##  [505] 5.3 4.3 3.4 5.1 5.3 4.1 6.3 4.7 3.0 4.5 4.8 4.4 5.2 5.4 5.6 5.7 6.6 3.4
##  [523] 4.5 4.2 4.0 3.5 3.2 3.4 5.5 5.8 4.7 3.7 3.4 3.2 2.1 4.3 4.7 5.1 3.3 5.1
##  [541] 3.6 4.0 4.2 3.4 3.0 5.8 5.6 5.7 5.3 4.4 3.5 4.7 5.9 4.4 4.7 2.9 5.2 3.1
##  [559] 6.2 4.2 5.4 5.4 5.4 4.2 5.9 3.7 4.2 2.6 4.1 4.5 5.3 5.9 5.7 4.6 4.5 5.2
##  [577] 6.7 4.7 4.7 4.7 6.1 3.1 3.7 4.5 3.5 4.0 3.2 5.3 4.9 3.3 5.8 4.1 4.0 4.2
##  [595] 4.6 4.0 4.4 4.6 5.1 3.4 6.2 4.4 4.8 3.5 4.9 5.0 5.1 4.1 5.1 3.7 3.8 4.7
##  [613] 3.3 4.2 6.0 5.1 4.7 5.3 4.8 4.4 3.6 3.4 4.8 3.4 4.7 3.8 5.9 2.7 4.8 5.9
##  [631] 3.6 4.6 5.5 5.0 5.5 3.0 4.4 3.2 3.3 4.5 4.4 5.1 5.7 5.0 5.3 6.2 3.4 3.9
##  [649] 3.7 4.6 5.1 3.1 4.4 3.5 4.0 3.4 4.1 3.9 4.9 5.5 5.9 5.0 2.9 4.5 3.7 4.5
##  [667] 5.6 4.5 4.5 3.3 5.5 4.5 4.8 4.9 5.3 4.8 5.9 4.2 4.2 4.2 4.0 3.9 5.5 5.6
##  [685] 4.8 3.8 6.0 4.2 5.2 4.9 3.9 5.2 4.4 3.7 3.6 3.6 4.6 5.9 5.0 6.1 5.8 3.9
##  [703] 5.3 5.0 5.0 5.9 4.8 3.8 3.7 4.7 3.0 4.1 2.8 6.1 4.2 1.9 4.7 6.2 5.4 5.5
##  [721] 4.3 4.8 5.8 4.1 5.1 3.2 4.6 5.5 3.4 4.8 4.3 5.0 5.0 2.9 5.1 5.5 5.2 4.4
##  [739] 5.2 2.3 2.8 3.6 4.0 5.1 4.6 4.7 5.3 4.5 5.0 4.9 4.5 4.2 5.0 5.0 5.6 3.8
##  [757] 4.7 5.9 4.6 2.7 4.9 4.3 6.3 4.5 5.3 5.5 6.8 4.7 5.1 4.2 6.1 3.1 4.8 5.0
##  [775] 5.5 3.7 4.7 3.9 4.8 4.5 4.4 5.0 3.0 4.4 5.5 4.8 4.0 4.7 5.4 4.5 2.8 2.5
##  [793] 4.5 4.1 4.6 5.0 4.5 4.1 4.4 3.7 3.1 4.7 5.2 3.5 4.9 4.8 4.4 3.4 4.5 5.5
##  [811] 4.4 4.4 4.0 3.9 3.9 5.3 6.4 5.0 5.0 3.9 5.1 3.7 4.8 5.1 5.2 2.6 3.8 4.3
##  [829] 4.3 4.6 4.9 5.1 3.7 4.6 5.8 3.8 4.0 5.1 4.4 5.6 3.9 5.0 3.5 2.4 3.2 6.3
##  [847] 2.3 4.9 4.0 4.4 5.0 6.1 3.8 4.1 4.1 4.8 4.0 4.4 3.8 2.8 4.1 4.2 5.4 5.3
##  [865] 5.2 4.6 4.1 5.6 4.4 4.7 5.8 6.1 4.5 5.6 5.4 4.4 5.2 4.1 4.8 5.3 5.4 4.5
##  [883] 3.8 5.9 4.9 4.5 3.4 4.9 3.3 5.4 4.9 4.7 4.4 5.3 4.8 2.7 3.8 3.9 4.5 4.0
##  [901] 4.0 3.2 4.3 6.2 5.7 3.9 4.6 3.0 4.7 5.6 3.1 4.6 3.5 5.1 3.9 5.9 5.9 4.7
##  [919] 3.8 4.5 4.4 4.9 4.6 4.8 3.4 2.7 3.4 3.3 4.4 5.2 6.4 5.3 5.0 4.6 4.6 4.2
##  [937] 5.0 5.4 3.5 4.3 4.9 5.9 2.7 5.1 6.2 4.5 3.5 4.4 2.3 4.4 4.1 5.2 5.1 6.3
##  [955] 3.9 3.1 4.6 4.8 3.8 4.9 3.3 4.5 3.2 2.5 4.0 2.8 5.5 2.8 5.1 5.4 4.4 3.4
##  [973] 3.5 3.4 3.0 4.0 6.3 3.7 5.7 4.7 5.0 3.2 4.7 3.9 3.4 5.0 3.0 3.9 4.8 3.4
##  [991] 4.5 3.5 6.1 3.4 4.7 5.3 5.5 5.0 3.7 5.9
## 
## $func.thetastar
## [1] 0.0034
## 
## $jack.boot.val
##  [1]  0.52034384  0.37891738  0.29098837  0.15058140  0.08422619 -0.04934037
##  [7] -0.12761628 -0.21463415 -0.36889535 -0.58696884
## 
## $jack.boot.se
## [1] 0.9829665
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
##    [1] 5.3 4.8 5.3 4.3 3.8 2.7 3.9 5.3 5.5 7.1 4.2 2.9 3.9 4.2 4.5 5.9 4.7 4.4
##   [19] 4.9 4.4 5.0 4.3 4.2 3.6 4.6 4.5 4.5 4.9 5.5 3.7 4.4 5.3 3.8 4.5 5.2 3.4
##   [37] 5.2 3.7 4.2 4.1 6.6 4.6 3.4 4.1 5.7 3.5 2.9 5.8 5.1 5.0 3.4 4.3 4.6 5.8
##   [55] 3.2 3.7 4.2 4.6 4.7 4.4 5.8 3.8 4.9 2.6 3.9 5.2 4.3 4.0 3.7 5.7 4.5 5.2
##   [73] 3.2 5.8 5.3 4.9 5.8 4.4 2.5 3.4 4.6 5.3 6.2 5.3 3.7 4.9 4.6 4.0 4.0 3.6
##   [91] 4.1 5.1 3.6 4.2 4.6 4.0 4.5 3.7 3.7 5.1 4.8 4.5 4.6 4.0 4.0 3.4 2.8 5.0
##  [109] 5.8 4.0 3.7 4.6 6.4 5.2 3.9 3.6 4.6 2.9 6.5 5.6 4.1 3.6 3.6 3.6 4.5 3.5
##  [127] 5.0 4.0 4.4 5.4 4.4 4.5 5.3 3.7 4.0 3.7 6.1 3.3 4.7 3.3 5.2 4.9 5.3 3.4
##  [145] 4.2 3.6 3.5 4.4 3.1 5.2 5.1 4.7 5.4 4.3 4.5 4.8 3.5 4.6 4.6 4.6 4.4 4.8
##  [163] 4.4 5.8 5.0 3.2 4.2 5.2 3.2 3.8 3.7 5.7 4.3 3.8 3.5 5.0 2.8 5.8 5.1 5.7
##  [181] 6.9 4.0 4.7 3.1 2.7 4.7 3.1 4.9 2.9 3.6 7.0 6.1 6.3 5.0 4.7 5.1 4.6 3.1
##  [199] 5.2 5.1 3.8 5.1 3.9 3.7 2.4 4.9 6.3 3.9 5.3 5.2 5.6 4.1 4.7 4.9 3.3 4.0
##  [217] 5.3 5.0 4.8 3.5 4.7 4.3 3.9 3.8 3.1 3.1 2.9 2.9 5.1 4.8 5.1 5.7 5.6 5.1
##  [235] 5.9 4.8 3.5 3.9 4.6 3.7 5.0 6.7 4.0 3.9 5.1 4.5 4.7 4.6 5.8 5.4 3.8 4.8
##  [253] 4.8 3.9 5.5 4.0 5.3 5.7 4.6 4.5 4.6 4.1 5.0 6.5 5.4 4.8 4.2 3.7 5.2 5.6
##  [271] 3.2 5.3 5.6 2.9 3.9 2.6 3.9 5.8 4.9 5.2 6.1 4.0 4.9 4.2 4.2 2.1 4.7 4.8
##  [289] 3.9 3.5 3.5 4.3 3.1 5.3 4.2 4.5 5.3 3.4 6.5 4.5 5.3 5.9 4.8 4.2 4.6 4.7
##  [307] 5.9 5.7 5.7 5.3 2.0 4.4 5.6 3.4 3.2 3.1 4.4 3.6 4.6 5.3 5.2 4.8 4.2 2.4
##  [325] 4.5 5.1 5.4 5.2 4.3 2.6 4.9 4.9 4.8 4.9 3.9 4.7 5.6 5.4 3.3 5.8 3.7 3.9
##  [343] 3.5 2.5 4.2 3.8 6.3 5.7 2.7 4.2 5.3 4.5 4.8 4.7 4.4 3.2 4.8 2.8 3.9 3.4
##  [361] 4.6 4.5 4.2 5.0 3.8 4.7 5.5 4.4 3.8 4.5 3.9 3.7 4.4 3.8 2.4 4.9 3.1 3.8
##  [379] 2.5 5.7 4.5 4.4 4.3 4.3 4.1 4.7 4.1 4.9 4.4 4.3 4.4 4.4 4.4 6.2 4.0 5.5
##  [397] 4.7 5.9 4.9 7.0 5.0 6.2 3.7 3.7 6.8 5.8 5.9 3.8 2.8 4.7 5.2 5.7 3.0 5.0
##  [415] 5.3 3.0 6.0 5.6 1.9 4.5 3.8 4.1 5.1 5.4 5.4 3.0 5.3 5.1 3.9 4.7 4.6 4.0
##  [433] 4.5 3.7 4.8 3.5 5.2 4.6 2.1 4.5 4.6 6.7 5.2 5.0 4.6 4.6 5.8 2.6 5.0 6.4
##  [451] 4.8 3.7 5.8 3.0 5.5 5.1 5.9 4.3 4.8 4.2 3.9 3.8 4.1 5.2 3.2 4.3 4.1 3.9
##  [469] 3.6 3.5 4.5 5.0 4.0 4.3 5.8 2.3 3.1 4.9 4.8 7.0 4.1 4.0 5.7 3.6 3.0 4.3
##  [487] 2.1 4.5 4.2 4.8 4.1 2.5 5.8 4.9 4.5 4.2 5.6 3.7 3.2 4.3 3.8 4.3 4.2 6.2
##  [505] 3.2 3.3 4.0 4.5 4.8 4.1 5.3 5.4 6.8 5.5 3.9 2.7 2.8 5.1 3.8 4.9 4.7 2.8
##  [523] 4.3 4.8 4.7 3.5 3.9 5.0 6.0 5.8 4.9 5.0 4.9 4.4 5.2 4.6 3.9 6.3 4.9 2.0
##  [541] 2.8 3.8 4.1 4.6 5.2 4.8 3.7 5.8 6.6 3.9 3.8 5.5 4.2 6.0 3.5 3.9 3.8 5.5
##  [559] 4.8 3.4 4.3 4.2 5.7 3.9 4.0 5.9 3.7 3.9 3.5 4.0 4.1 5.9 3.3 4.0 4.0 4.5
##  [577] 5.1 5.2 4.8 3.5 6.7 4.9 4.4 6.9 5.4 4.8 4.8 6.0 3.4 4.5 4.4 4.5 5.1 3.9
##  [595] 3.6 3.6 4.0 5.6 3.7 3.7 4.8 4.0 4.0 5.2 2.7 3.3 5.1 4.6 4.5 4.7 3.6 4.7
##  [613] 3.2 5.6 6.6 3.9 3.4 3.9 4.3 3.7 4.6 4.4 4.8 3.4 4.0 5.9 4.3 4.7 4.1 5.9
##  [631] 5.6 4.0 4.1 4.8 4.1 4.0 4.6 4.2 6.4 5.5 2.2 4.0 2.8 4.1 3.5 5.7 5.4 4.7
##  [649] 4.6 3.4 3.7 5.1 4.5 3.6 4.9 5.1 4.5 4.4 5.8 6.4 4.5 4.3 3.9 4.3 4.9 5.2
##  [667] 5.6 6.6 3.2 4.5 5.0 3.1 4.5 5.1 4.1 4.6 4.1 3.6 4.4 5.2 3.8 5.5 4.9 4.3
##  [685] 6.5 5.1 4.4 2.9 5.8 5.6 4.0 3.3 3.5 4.8 5.3 2.5 5.0 3.5 3.8 5.9 5.5 3.2
##  [703] 4.8 2.8 5.7 4.9 4.2 6.0 3.9 3.6 3.9 3.4 5.8 3.5 4.6 5.7 3.9 5.2 4.0 3.3
##  [721] 3.1 5.3 4.3 5.7 5.9 3.4 4.2 5.1 4.4 3.9 4.0 4.2 4.4 3.8 3.4 3.5 3.5 4.4
##  [739] 5.2 3.6 3.5 4.8 6.3 6.7 2.3 4.8 3.9 3.9 5.1 3.5 6.7 6.9 5.5 4.6 4.5 5.4
##  [757] 5.6 6.6 4.8 5.6 3.6 4.2 6.6 4.6 4.3 5.3 4.2 4.0 4.4 5.0 4.2 5.1 4.7 6.0
##  [775] 4.6 4.1 4.5 3.8 3.1 4.1 5.4 5.5 3.9 5.0 2.9 5.8 5.4 4.4 5.6 6.0 3.3 3.8
##  [793] 4.8 3.5 4.2 4.1 5.4 4.1 4.3 4.9 4.8 4.1 4.4 4.6 4.1 3.6 4.5 3.3 5.2 4.9
##  [811] 3.7 4.8 2.8 4.6 6.5 4.1 6.4 4.2 3.3 3.7 5.5 5.8 5.5 5.0 3.7 3.4 3.1 6.0
##  [829] 6.0 3.3 3.5 2.3 5.0 5.9 4.2 3.5 6.6 2.8 6.6 2.4 3.1 6.2 4.8 2.5 4.2 4.1
##  [847] 4.6 6.6 3.9 5.5 2.9 3.4 4.4 4.5 4.6 4.3 5.5 3.0 4.9 2.8 4.2 3.1 5.0 4.5
##  [865] 4.9 3.9 3.2 4.7 4.6 3.1 4.9 5.8 3.7 3.9 2.9 3.9 4.4 6.1 4.6 5.4 4.8 3.2
##  [883] 4.3 3.8 3.9 4.2 5.6 3.8 5.8 5.6 5.6 3.9 5.2 3.2 5.2 5.8 5.2 5.3 3.9 3.9
##  [901] 4.3 6.4 4.6 4.8 4.2 6.1 4.4 4.8 3.8 5.5 6.8 4.4 5.2 5.5 3.4 6.2 4.7 6.2
##  [919] 5.9 4.8 4.3 4.5 3.3 4.2 5.8 5.5 4.8 3.3 4.1 5.9 3.8 5.2 3.3 4.8 3.7 3.6
##  [937] 3.5 4.6 3.6 4.6 3.9 3.4 6.4 3.0 5.0 3.3 2.8 4.8 5.1 4.6 5.5 5.8 4.4 5.2
##  [955] 3.2 4.8 5.6 3.8 4.0 3.5 3.5 4.6 5.3 3.6 3.2 4.4 5.7 5.7 5.1 5.0 3.8 3.4
##  [973] 3.9 5.8 3.9 3.4 3.4 5.7 5.2 4.4 3.7 3.9 4.4 4.3 4.8 4.3 3.3 4.8 5.2 5.6
##  [991] 5.1 3.3 4.4 2.7 3.9 6.5 4.7 4.6 4.4 4.9
## 
## $func.thetastar
##   72% 
## 5.028 
## 
## $jack.boot.val
##  [1] 5.7 5.6 5.5 5.3 5.0 4.9 4.9 4.7 4.6 4.4
## 
## $jack.boot.se
## [1] 1.267123
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
## [1] 1.324073
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
##   2.413059   3.862105 
##  (1.013447) (1.802523)
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
## [1]  0.83357950  0.23619872  1.22844453  0.36167511 -0.04910956  0.08779250
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
##    [1]  1.70618903  0.42638925 -0.08474916  1.69587065  1.07443517  1.67143056
##    [7]  0.80076427  0.99861028  0.65684982  1.33173069  1.67048432  1.64371688
##   [13]  1.32906141  1.67412919  0.76566748  2.50322271  0.17480587  2.46207109
##   [19]  1.07755102  1.96720121  1.45944539  1.13103062  0.98587947  2.34516366
##   [25]  1.44581783  0.84176709  0.24059658  0.65141792  0.32035869  2.34917987
##   [31]  1.67630705  1.00529758  2.59740268  2.23448531  1.98040382  2.05516122
##   [37]  0.38235788  0.30559355  1.51790488  0.36177135  2.43493521 -0.04002771
##   [43]  0.48087709  1.03391695  1.09537906  2.03942029  1.63273931  0.60493407
##   [49]  0.52739013  0.64255499  1.82021822  1.43799048  2.06795866  1.20328993
##   [55]  0.97085622  1.66237112  0.53132648  0.80607725  0.20550304  1.33590119
##   [61]  0.50984806  1.29024412  0.85361112  0.61060211  0.29782142  1.45618025
##   [67]  2.43876158  1.67432253  1.02789892  1.61506943  1.02377284  0.09903251
##   [73]  1.49743071  1.11224783  1.36053151  0.67892430  1.32980670  0.15232480
##   [79]  0.50384366  0.69110145  1.39858101  0.37955931  0.80737462  1.66774183
##   [85]  0.35881035  1.99061864  1.35305410  0.52486090  1.00984626  0.85500958
##   [91]  1.47552780  1.46745966  1.12674629  1.27001524  0.91084231  0.50248502
##   [97]  1.11660052  1.01604778  0.52392643  1.67900572  0.68214137  2.63933381
##  [103]  0.74775342  1.45697875  1.68869205  1.53751435  1.31738768  2.38354892
##  [109]  1.25313844  1.59860003  1.11797261  1.35217435  0.43217046  0.20689228
##  [115]  0.22153364  0.99697949  2.54902152  1.47357342  1.32227679  0.62310844
##  [121]  1.28383998  1.67089123  0.41059034  0.67229278  0.47548566 -0.16709331
##  [127]  0.23690363  1.31585932  1.06617653  0.79165211  0.85319452  0.53611790
##  [133]  1.94347948  0.86525774  1.10132495  0.63451930  2.30557103  1.35639729
##  [139]  1.29425365  1.66280685  1.69708778  1.28065461  0.87001971  0.62716049
##  [145]  0.99421795  0.85606347  1.32507491 -0.07076204  2.04144169  0.78982122
##  [151]  0.44586750  2.41713640  1.22436000  1.52485893  1.02867215  1.00663417
##  [157]  2.40347135  2.48740218  0.50421039  1.32155529  1.47651089  1.60071709
##  [163]  1.29669671  1.11960276  1.00364615  1.34954623  1.70828510  1.00957630
##  [169]  1.70781501  0.61022148  0.85875710  0.45452810  0.46461793  1.10332110
##  [175]  1.99591512  1.10582668  2.07076981  2.38030304  0.67263568  0.85320696
##  [181]  2.54483273  1.18620933  1.45022247  1.94828974  1.33570611  0.60190488
##  [187]  1.41850789  1.20880227  1.46663930  0.85620638  0.78893353  0.20799779
##  [193]  1.11376613  0.96714604  1.32616925  0.22300708  1.56375994  0.86693220
##  [199]  1.08423978  1.49354589  0.66245297  1.10790550  1.62455397  1.26070048
##  [205]  0.83721596  0.11615619  1.60747378  1.31565035  0.62375820  1.00882049
##  [211]  0.84894143  0.60981226  0.45194120  0.51379475  0.85774031  1.12775785
##  [217]  1.41993664  1.05429803  0.96029498  1.20894995  1.20479520  2.01063968
##  [223]  1.32604039  0.40480011  1.68591087  1.45730780  0.89408657  2.32330790
##  [229]  1.10829057  0.89590044  1.88172869  0.53121679  1.60249296  0.15852744
##  [235]  0.41332226  1.70222756  1.49552838  1.10152736  1.18751986  1.38138039
##  [241]  2.56574171  1.72193483  1.36299777  1.00927462  0.83612411  2.03684541
##  [247]  0.80076786  0.99343001  2.57360889  1.37057012  2.27433916  0.14282126
##  [253]  0.69137659  1.32616925  0.80236064  0.83842651  1.33568290  0.99097006
##  [259]  2.37809594  0.76100125  1.91794692  1.25846985  0.62673559  0.48989866
##  [265]  1.05284532  0.88660357  1.66886135  0.85798176  0.66632608  1.45908216
##  [271]  1.06409845  1.01999809  1.20901585  0.28454483  1.34216270  1.20169874
##  [277]  2.29420110  0.77756742  0.83250884  0.61795846  0.36615002  0.99590432
##  [283]  2.33183297  1.27527046  0.21664057  0.66693064  1.33270210  0.13731857
##  [289]  0.97974851  1.02036615  0.81820551  0.59758972  0.56049182  2.39318091
##  [295]  0.81587714  0.79492346  1.00624461  2.39255691  1.58101063  1.20851041
##  [301]  2.40098814  2.28273406  1.01874779  1.31713510  0.99858392  0.16724466
##  [307]  1.01685687  0.78799676  2.54607093  1.37718761  0.66826832  1.59050165
##  [313]  1.63856364  0.80924343  2.26905662  0.13125534  0.79767584  1.01008934
##  [319]  0.50679687  0.66145042  2.00016009  0.55718786  0.60177246  1.18168466
##  [325]  1.97638054  1.00517851  0.24835420  2.28845793  1.26213127  1.53786613
##  [331]  2.20679495  1.62105715  0.83898956  0.53560007  0.81647135  1.65420148
##  [337]  1.07098506  1.02789892  1.33269656  1.03014192  0.80741898  1.40186831
##  [343]  0.84573295  0.54660147  1.31259606  0.62086686  0.84985393  1.38537673
##  [349]  0.67242818  2.52093830  1.47683813  0.79341157  0.14848115  1.19926989
##  [355]  1.00500415  1.97056759  0.79925344  2.02446572  1.39264384  1.21421414
##  [361]  0.83970641  0.28947515  0.38057912  0.87993255 -0.44661621  1.34152010
##  [367]  0.90186004  1.01334572  0.47969513  0.62927796 -0.04002771  0.77379943
##  [373]  1.67000805  1.32273556  2.49798163  1.33330637  0.30613171  1.07287523
##  [379]  1.33628552  0.86889203  0.45303199  2.01636236  2.37874174 -0.07582513
##  [385]  1.45019055  0.77832487  0.57094961  1.32997299  0.41780407  1.01547025
##  [391]  1.45303280  1.01934672  1.67671675  1.36161622 -0.22252355  2.49145910
##  [397]  2.50000951  1.64653636  1.37197833  0.66146293 -0.19920806  0.11667528
##  [403]  1.01067324  0.91506594  1.59502356  1.21264335  0.99180574  1.36109156
##  [409]  0.49932171  2.54848473  0.81818528  1.53739311  0.53551921  1.34961924
##  [415]  1.00529758  0.55730839  2.61773863  2.03389331  1.63360461  1.03667775
##  [421]  1.64891780  1.46597927  1.22804560  1.16947120  0.48002693  0.31340274
##  [427]  1.66569583  1.49024290  0.76533260  0.93033550  1.31235854  0.63397923
##  [433]  2.39302660  1.19746583  0.52279122  1.43684911 -0.12810706  1.31827821
##  [439]  1.02475351  0.52081597  0.79713968  0.97643920  2.55732285  0.97428295
##  [445]  2.52249412  0.80467550  2.00681104  1.60471494  0.67353948  1.67442680
##  [451]  1.44310499  1.63149551  2.56158545  1.08863108  0.66935280  0.56845114
##  [457]  0.37350635  0.46107859  0.66389665  0.43769155  0.35472963  0.99086424
##  [463]  2.38589110  2.03210716  1.04942279  0.87752913  2.46228265  1.21461983
##  [469]  1.31648215  1.34961924  1.61388870  1.31552435  2.33580817  1.01469147
##  [475]  1.46246070  1.55789419  1.50080733  1.68891751  1.21725946  0.90732118
##  [481]  0.81142037  1.04636718  0.84425228  2.35239510  0.86833293  0.39876117
##  [487]  1.10872121  1.62918518  2.02643105  1.36401038  2.54506300  0.83592218
##  [493]  0.66887429  2.03792275  0.58065989  0.38345893  1.64529860  0.22424138
##  [499]  1.84248701  0.98560658  1.41452185  1.53720134  1.40762692  0.48796977
##  [505]  0.37198883  2.02522101  0.94559888  1.32395095  1.02972818  1.92497518
##  [511]  0.24396614  2.00075246  1.66770396  2.33912679  1.93170899  1.33021614
##  [517]  1.46093170  0.81943024  0.87606834  0.54757493  1.57791671  0.66635050
##  [523]  1.68611726  0.53178841  1.32547367  2.56699459  1.57580987  1.55012162
##  [529]  0.35659417  1.59486113  1.29470959  0.28967809  0.62603673  1.63068319
##  [535]  2.53062367  0.38415925  0.93833143  1.70025307  0.78551023  0.45398958
##  [541]  2.65181327  0.94384364  1.97769430  1.55214244  1.44976242  0.64082983
##  [547]  1.30826650  0.17317734  1.01043015  0.80723246  0.63747930  2.56642694
##  [553]  0.92939880  2.52126815  1.11551599  0.62603673  1.45473282  1.18000581
##  [559]  0.51272130  0.81789264  1.36944132  1.01243851  0.76115913  1.64256049
##  [565]  2.05689332  1.63938101  0.79988665  0.64038499  0.58270537  0.45877671
##  [571]  1.21304068  2.01310521  1.12393682  0.52194889  0.62582546  0.37860095
##  [577]  2.48658754  0.98579232  2.58175322  1.98681662  1.88729881  1.20659961
##  [583]  1.25311164  1.00006306  1.66829355  1.89184900  2.11968830  2.48855255
##  [589]  0.90015441  1.43075761  2.04151216  0.26185102 -0.46944838  0.37350635
##  [595]  2.29548072  2.14753038  1.02574165  1.33037417  0.55190231  1.00504619
##  [601]  1.00223328  2.34289043  1.33625190  0.81067526  0.11177024  2.01896520
##  [607]  1.66390582  2.25655450  1.48767359  1.95559864 -0.43410099  0.31636112
##  [613]  0.96972483  0.52780649  1.37699873  1.56305318  1.58647410  1.00416323
##  [619]  0.61457700  0.97375294  0.94579373  1.20548621  2.06508277  0.98725768
##  [625]  0.48942878  1.33459637  0.96813549  1.43981169  1.06011307  0.88057939
##  [631]  1.64323304  1.69767038 -0.28104982  1.51071600  2.53217700  2.45123055
##  [637]  2.06006623  1.25813151  0.35187820  0.49861073  0.90544154  0.66124072
##  [643]  1.08879768  1.41518351  1.03519727  0.41488661  1.01979390  1.57146833
##  [649]  1.62134878  1.00420721  2.03156307  0.57537883  1.01290674  1.34343922
##  [655]  0.86506902  0.23031029  0.29533851  2.26654013  2.54796650  0.98804463
##  [661]  0.90198315  1.44168869  0.84262691  1.44277609  0.54099641  1.20356906
##  [667]  0.79980680  1.21220771  1.35498703  1.77998821  2.25949970  2.26142532
##  [673]  0.79755634  0.94243477  1.62701007  2.49054766  0.08832694  1.52632919
##  [679]  2.49560703  1.59573775  1.46940547  1.68051851  1.34225583  0.54513438
##  [685]  1.09583817  1.26458654  1.29708906  2.52845340  1.56695111  2.45509267
##  [691]  1.25287074  0.22565029  1.40004808  2.51046710  0.29824143  0.23930066
##  [697]  0.50692210  2.12707225  0.55190583  0.62264755  0.30339089  1.20616698
##  [703]  1.46776333  1.06206135  0.95093610  2.19599953  2.00505976  0.07601874
##  [709]  0.56286852  0.17923847  0.59058037  0.90099468  0.67176681 -0.06991486
##  [715]  1.66569583  1.44546736  0.97586390  0.69945899  2.17893267  0.57478225
##  [721]  2.45426640  1.01772540  0.52324361  1.63723450  1.51186741  1.68806117
##  [727]  1.00350108  1.37056323  2.02094573  1.66989315  0.68428016  0.53224692
##  [733]  0.83278115  0.83182533  0.12681429  1.38106870  0.90486265  1.97302437
##  [739]  2.44824584  1.34377761  1.38119697  1.42253023  0.75981412  1.29727121
##  [745]  1.63752951  0.84613044  1.30094561  0.50958797  1.19588667  2.21256282
##  [751]  0.97890527  0.13134551  0.61951738  0.52194889  1.00024574  1.11224715
##  [757]  1.26829744  2.41424192  1.41263123  0.85450650  1.46583512  0.65052402
##  [763]  0.79295049  2.55779517  2.02312815  0.37921559  1.98874803  2.21111257
##  [769]  0.43986329  0.37431642  0.95863656  0.56913806  0.78062597  1.71792983
##  [775]  1.71382282  2.48974015  1.61486769  1.32642280  0.53113405  2.44098707
##  [781]  1.10832460  1.34525757  1.44746598  1.56843748  2.62022909  1.71600125
##  [787]  2.06137589  0.77832487  0.55025593  0.77890893  1.13968145  0.85397006
##  [793]  1.21372891  0.51640930  0.27367915  0.91409780  1.00106084  0.42573748
##  [799]  0.90822375  1.20616698  1.33081004  1.43916716  1.08109616  0.53517928
##  [805]  0.84843250  1.68757769  1.63287476  1.31664838  1.24462135  0.80288542
##  [811]  2.56755001  0.83046775  1.00724764  0.79746527  1.50154776  0.82241099
##  [817]  0.81939858  1.21220757  1.15393617  1.32346155  0.61606407  1.32769205
##  [823]  0.97617617  1.31938268  1.51626916  0.51048086  0.64838141  0.60676053
##  [829]  1.18620933  1.01553819  0.63258885  0.64105175  2.03259809  1.69684848
##  [835]  1.54973980  0.44905563  1.00671943  1.46624291  1.05976975  1.53382268
##  [841]  1.31713510  0.75920883  1.52395575  2.26794434  0.63166243  1.63217302
##  [847]  0.89655827  1.00390680  0.98726610  1.31217485  1.66280685  1.67478471
##  [853]  1.67111868  0.65673068  1.43280482  0.79114231  1.65423463  0.97586781
##  [859]  1.18379782  1.31447301  0.34592680  0.66416784  0.59489024  1.30815350
##  [865]  1.44518116  1.30216692  0.57478173  0.95545963  0.80392186  0.39466670
##  [871]  0.88469770  1.35681029  0.24744358  2.56555066  0.90879847  1.90250969
##  [877]  1.64231972  0.79149743  1.36937575  0.97072344  1.47720627  0.63833943
##  [883]  0.98890702  1.02546914 -0.18084194  2.53630504  0.99271946  2.00373612
##  [889]  1.33651886  1.58139387  1.32471872  1.59860003  1.58606419  1.39265032
##  [895]  1.20880227  1.99810197  1.03206483  0.95618506  1.65659516  1.31821814
##  [901]  0.99358884  2.30153503  2.00926616  1.67861483  2.38610409  0.87199323
##  [907]  0.92472291  1.57621135  1.01163926  0.85697504  1.51945502  1.17033737
##  [913]  1.46741604  2.09499938  0.39831789  1.19476440  0.52189225  1.11488513
##  [919]  0.56786457  2.51597526  0.80441641  1.70828510  2.64252810  1.02320898
##  [925]  1.07198315  0.86007414  0.80497031  2.13619099  0.63396478  0.39728711
##  [931]  1.68807058  1.08255678  1.26315011  1.20122881  2.38405373  0.79755634
##  [937]  0.80197000  0.07633758  1.44312390  0.85797614  1.67361232  0.44236815
##  [943]  0.96516105  2.04995805  2.55783494  2.50230334  1.20358744  1.21898025
##  [949]  0.68223011  0.17975459  2.02312815  2.49526693  1.32511846  0.99858779
##  [955]  1.32565355  0.53447123  0.81119474  0.44807469  0.48087709  0.83929704
##  [961]  1.00671943  0.23467165  1.07973291  1.42233822  1.20205464  1.50715738
##  [967]  0.57672511  1.44018094  1.19479681  2.60701100  1.00524014  0.53900253
##  [973]  0.92380563  0.10852691  2.05945545  0.90685708  1.47091530  0.79608033
##  [979]  0.13692438  0.93426329  1.03104767  2.42996947  0.66145042  1.27565692
##  [985]  0.60190488  1.50005083  1.31782543  0.44607999  0.44561938  0.79962676
##  [991]  1.20837983  0.89231884  0.56936013  1.12802794  2.14288643  0.42998353
##  [997]  1.64620150  2.53548461  0.67084701  0.79603215
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
##      mean         sd    
##   0.6248006   0.4653732 
##  (0.1471639) (0.1040588)
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
## [1] -0.07264063 -0.83621635  0.92922592 -0.28194563  0.45650290  0.48412033
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
## [1] -0.0294
```

```r
sd(results2$thetastar)  #the standard deviation of the theta stars is the SE of the statistic (in this case, the mean)
```

```
## [1] 0.900459
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
## t1*      4.5 -0.02752753   0.9061209
```

One of the main advantages to the 'boot' package over the 'bootstrap' package is the nicer formatting of the output.

Going back to our original code, lets see how we could reproduce all of these numbers:


```r
table(sample(x,size=length(x),replace=T))
```

```
## 
## 0 1 2 3 6 8 9 
## 1 1 1 2 2 1 2
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
## [1] -0.008
```

```r
se.boot
```

```
## [1] 0.9119472
```

Why do our numbers not agree exactly with those of the boot package? This is because our estimates of bias and standard error are just estimates, and they carry with them their own uncertainties. That is one of the reasons we might bother doing jackknife-after-bootstrap.

The 'boot' package has a LOT of functionality. If we have time, we will come back to some of these more complex functions later in the semester as we cover topics like regression and glm.

