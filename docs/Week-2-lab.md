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
## 1 4 5 7 9 
## 1 2 3 1 3
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
## [1] -0.0015
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
## [1] 2.72537
```

```r
UL.boot
```

```
## [1] 6.27163
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
##    [1] 3.2 5.3 4.7 3.1 3.6 5.3 5.8 3.5 5.6 3.6 5.1 4.4 3.7 5.1 4.6 3.7 2.6 3.8
##   [19] 4.4 5.4 5.2 3.1 4.3 5.6 5.6 3.6 4.8 5.6 5.6 5.8 4.8 6.0 5.0 3.6 6.3 5.3
##   [37] 3.0 5.9 4.4 5.5 7.0 4.9 4.7 5.7 3.6 5.2 4.1 4.4 4.1 5.3 5.8 3.9 4.4 5.8
##   [55] 6.1 3.9 4.2 3.3 3.4 3.9 2.6 6.3 4.6 4.6 2.5 4.0 4.8 3.4 3.6 3.9 3.2 5.5
##   [73] 4.1 4.9 4.5 5.9 4.0 3.9 5.5 3.9 3.1 6.0 4.6 4.0 5.2 5.5 4.5 4.6 5.1 3.3
##   [91] 5.4 5.9 5.1 4.4 4.3 7.1 4.6 6.1 3.7 3.0 5.3 4.8 4.7 4.2 4.6 5.2 3.7 4.5
##  [109] 5.5 4.4 4.6 5.8 3.1 3.0 4.7 5.4 3.7 4.8 2.7 4.0 5.0 4.8 4.3 3.6 6.5 5.1
##  [127] 5.5 5.7 4.5 5.4 5.8 3.9 4.1 3.0 4.9 4.6 4.7 4.9 4.1 3.7 3.1 3.7 5.6 5.6
##  [145] 4.4 4.2 3.8 2.5 4.0 6.3 4.0 3.9 3.4 4.1 5.1 3.6 4.1 3.9 5.2 5.3 6.1 4.1
##  [163] 4.4 3.3 4.1 3.6 4.4 5.4 2.8 3.0 4.0 5.0 3.5 3.2 4.2 4.9 3.5 5.0 4.9 5.0
##  [181] 6.0 3.8 5.7 3.1 4.4 5.3 2.1 6.0 4.2 4.9 4.4 4.8 5.8 5.0 4.8 5.1 3.0 5.0
##  [199] 3.8 4.9 3.0 5.1 5.7 4.3 3.4 6.3 4.7 5.1 4.4 4.1 5.0 3.5 3.7 5.8 2.8 2.2
##  [217] 3.6 5.1 3.7 5.0 4.1 3.7 5.1 3.8 5.1 5.6 4.1 5.7 4.7 4.1 5.2 5.6 4.2 3.8
##  [235] 6.1 4.8 5.1 2.8 4.0 5.4 4.5 4.0 4.0 3.4 4.4 4.4 5.8 3.4 4.9 4.7 3.0 4.4
##  [253] 3.3 3.6 4.5 5.1 4.5 4.6 3.6 4.3 3.7 3.4 5.0 5.0 3.1 4.5 3.8 5.0 3.5 5.1
##  [271] 6.1 5.4 5.9 3.4 3.4 3.1 2.6 4.0 4.0 4.8 4.0 4.8 4.9 3.8 5.5 4.5 3.5 5.1
##  [289] 3.5 4.2 2.9 4.0 4.6 3.8 4.4 5.2 5.1 3.8 3.7 4.6 2.5 3.5 4.3 3.4 3.4 5.7
##  [307] 5.0 4.3 5.2 5.4 4.0 5.4 4.6 4.3 3.9 3.4 3.2 3.7 4.9 5.1 6.0 4.9 5.6 4.3
##  [325] 4.1 4.0 5.2 4.6 6.4 6.0 5.0 4.1 5.6 3.3 4.1 3.5 4.3 3.1 3.2 4.0 4.4 4.4
##  [343] 4.4 5.2 5.1 4.2 5.2 3.7 4.2 2.5 5.0 3.6 2.8 4.7 4.7 5.5 4.3 4.9 3.4 3.2
##  [361] 5.0 2.9 6.2 5.4 2.6 5.0 3.8 5.3 4.5 5.8 5.5 4.3 4.7 7.1 5.5 4.6 5.4 3.5
##  [379] 2.9 4.7 5.6 5.0 5.2 3.5 3.6 4.5 5.5 6.2 3.4 5.4 4.8 3.9 4.5 4.3 4.2 4.7
##  [397] 4.4 4.0 5.3 4.5 4.8 4.8 4.8 5.4 4.3 4.6 4.4 4.6 4.8 6.0 3.6 3.4 4.2 2.7
##  [415] 5.1 3.2 3.3 3.5 3.7 5.2 4.5 5.0 5.3 3.7 3.0 4.7 3.8 5.0 3.8 3.8 3.7 3.5
##  [433] 5.8 3.2 5.3 3.4 3.5 3.7 3.4 4.4 4.2 4.6 4.6 6.1 4.3 4.3 4.4 3.5 3.9 6.9
##  [451] 4.1 3.8 5.3 3.7 4.5 4.2 4.0 4.9 4.9 2.4 5.1 5.5 4.7 4.7 4.1 5.7 4.0 2.7
##  [469] 2.6 5.3 3.8 5.4 5.5 3.1 3.5 4.2 5.4 4.1 2.8 6.7 4.7 4.4 5.6 4.8 4.3 4.4
##  [487] 3.9 6.4 4.7 4.8 4.9 3.3 4.2 4.7 4.2 3.7 4.8 4.5 4.8 2.8 5.6 5.0 5.3 4.6
##  [505] 3.3 4.3 3.5 4.1 5.0 5.8 4.3 3.7 5.3 4.3 3.2 6.2 4.5 4.0 3.4 4.2 3.8 4.2
##  [523] 5.7 4.8 4.1 4.4 4.6 3.6 3.4 3.8 5.0 3.7 3.8 5.6 4.4 4.9 5.0 5.4 5.4 3.8
##  [541] 4.2 4.4 4.0 4.9 4.3 3.9 4.2 3.8 3.5 4.1 3.8 5.4 6.3 4.0 6.0 4.5 4.2 4.0
##  [559] 5.2 5.0 6.1 5.3 4.0 4.8 3.8 5.6 3.3 3.6 4.1 4.7 3.6 4.8 5.1 5.2 5.8 4.2
##  [577] 3.6 5.7 3.8 3.1 5.3 4.2 5.3 6.2 2.9 4.7 5.3 5.0 4.9 6.1 4.5 4.6 4.9 3.8
##  [595] 4.1 5.5 5.4 5.7 4.6 4.9 2.9 4.2 4.4 5.6 4.4 4.0 5.5 4.2 4.5 5.4 5.0 4.8
##  [613] 6.3 4.7 5.4 5.7 3.0 5.9 5.5 3.2 3.0 4.1 3.5 3.8 5.1 5.6 4.5 6.2 3.7 4.1
##  [631] 4.4 5.8 4.8 3.7 4.5 5.0 3.8 5.2 4.8 5.9 5.4 3.9 4.4 4.3 3.0 2.8 4.6 5.1
##  [649] 5.0 6.2 2.6 3.7 5.9 3.2 5.6 5.2 5.9 3.9 4.2 4.3 4.7 4.5 4.7 4.0 3.9 5.8
##  [667] 4.2 5.1 6.4 4.9 4.9 4.6 5.0 3.0 4.2 4.7 4.2 3.1 4.1 5.1 4.0 2.1 4.7 3.5
##  [685] 4.4 3.4 4.7 4.0 5.9 3.4 4.6 5.3 3.8 5.2 4.5 4.7 4.7 5.4 3.6 6.0 3.3 6.1
##  [703] 6.0 5.1 6.6 3.9 6.7 4.4 6.0 6.1 3.8 3.6 5.0 4.9 4.5 5.4 5.3 4.6 4.8 4.7
##  [721] 4.3 3.9 4.2 4.1 5.4 5.6 5.6 3.8 4.4 5.7 4.4 4.8 5.5 5.6 4.6 4.2 3.0 5.9
##  [739] 4.9 5.2 4.6 4.3 4.9 4.1 4.3 4.0 3.4 3.7 3.9 4.9 3.9 5.2 5.3 4.1 4.3 5.6
##  [757] 5.6 5.1 4.7 4.4 4.2 4.0 3.8 4.9 4.2 4.7 3.5 2.7 5.9 3.6 4.0 6.8 3.9 5.4
##  [775] 5.0 4.4 3.9 4.7 6.1 5.3 4.5 5.9 5.8 5.9 5.6 6.0 4.6 5.1 3.9 2.8 3.4 4.9
##  [793] 4.6 4.5 5.3 4.4 4.1 4.4 6.7 4.5 6.4 4.1 4.5 4.2 5.4 6.0 5.5 4.2 3.9 4.7
##  [811] 3.3 4.3 3.6 3.2 5.7 4.9 3.9 4.6 5.2 5.3 4.3 3.9 4.3 3.5 5.1 5.0 1.9 2.8
##  [829] 3.2 5.0 5.4 5.7 4.4 3.6 3.7 5.1 3.9 3.6 5.6 4.3 4.4 5.2 4.0 4.6 4.0 4.4
##  [847] 6.3 4.1 5.2 2.9 4.3 3.9 4.7 6.1 4.2 4.0 3.1 4.7 3.8 4.2 3.0 4.1 5.2 4.8
##  [865] 5.0 5.0 3.6 4.2 5.0 4.7 1.9 4.7 3.7 4.5 3.8 4.9 4.0 4.0 4.2 4.0 3.8 6.4
##  [883] 3.8 4.6 5.5 4.9 5.4 4.2 2.3 3.9 3.2 4.9 4.7 4.8 4.3 5.4 6.6 4.7 3.7 5.6
##  [901] 4.8 5.0 3.7 3.8 4.5 5.2 4.8 4.4 5.3 4.2 5.9 4.7 2.8 4.4 3.1 3.0 3.9 4.2
##  [919] 4.3 4.3 4.2 4.3 4.6 4.3 4.9 5.4 3.6 4.3 5.2 5.3 4.1 4.4 3.9 4.2 4.3 4.1
##  [937] 5.2 3.6 4.8 4.3 5.3 5.3 4.6 4.2 2.0 3.4 3.2 4.6 5.1 4.3 5.0 5.7 5.1 3.6
##  [955] 4.9 4.7 4.8 4.9 4.4 4.7 5.1 4.6 4.3 7.0 4.4 3.4 5.3 3.9 4.9 4.5 5.1 3.1
##  [973] 5.2 4.7 4.4 3.3 4.2 2.2 5.6 4.2 4.0 5.4 4.3 4.2 4.9 4.3 2.4 5.3 6.2 3.0
##  [991] 5.3 4.0 4.6 5.2 4.9 3.6 4.0 4.7 5.1 4.3
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
##   2.8   6.2
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
##    [1] 6.1 6.1 5.3 3.6 5.7 5.1 4.9 5.1 4.1 4.1 4.4 5.4 4.3 4.1 3.5 4.9 3.6 5.9
##   [19] 4.3 2.9 3.9 3.7 3.9 5.1 3.7 4.6 3.0 4.3 5.5 4.4 4.4 5.0 4.4 4.9 5.2 5.0
##   [37] 4.5 4.3 5.7 4.5 3.6 5.0 4.0 2.7 5.4 5.0 5.8 4.2 4.4 6.9 4.2 4.2 3.3 3.8
##   [55] 6.2 4.4 4.0 5.6 5.2 3.2 5.4 4.8 5.2 4.5 4.1 4.3 4.3 6.1 5.3 5.6 3.3 6.5
##   [73] 3.5 3.9 4.5 3.2 3.7 4.6 5.6 5.8 4.4 3.6 3.3 5.3 3.9 5.0 3.8 4.7 3.8 4.0
##   [91] 4.3 3.3 5.3 4.3 4.4 4.4 4.4 6.1 4.6 3.0 4.8 4.0 4.3 4.3 4.1 3.5 4.4 4.1
##  [109] 3.1 3.1 4.7 3.7 3.7 4.2 5.7 5.2 4.2 2.9 4.3 6.9 5.7 4.4 5.0 3.5 3.2 5.4
##  [127] 6.0 6.0 4.6 4.2 5.4 4.1 5.0 3.3 4.5 3.9 5.1 4.4 2.1 5.2 5.3 3.2 4.4 4.5
##  [145] 5.1 3.7 3.0 5.2 5.3 4.4 4.8 5.3 4.9 4.9 4.0 5.2 4.9 3.7 3.8 4.1 5.1 4.2
##  [163] 4.7 4.0 5.1 4.0 4.2 3.4 6.0 3.6 3.0 3.8 5.8 3.7 4.4 4.8 5.0 4.3 3.7 4.9
##  [181] 5.3 4.1 3.2 5.7 4.6 6.3 3.7 3.7 6.0 3.3 4.6 3.3 5.7 5.4 5.0 5.3 3.1 4.4
##  [199] 5.1 5.7 4.1 4.8 4.4 4.3 2.2 3.9 4.8 4.5 6.0 4.3 4.6 5.0 4.7 5.6 5.0 3.7
##  [217] 4.7 5.0 5.2 3.7 4.1 3.9 2.8 4.0 5.1 3.7 2.5 6.1 3.7 3.2 4.5 4.5 3.3 5.0
##  [235] 4.5 3.9 5.1 4.0 4.6 6.2 3.7 4.8 5.2 4.8 2.8 5.7 4.4 5.3 4.9 4.1 4.9 3.8
##  [253] 6.1 5.4 3.4 4.1 3.2 4.4 3.9 5.6 4.6 4.1 4.8 4.0 3.2 3.5 4.2 5.4 5.8 4.3
##  [271] 3.4 5.3 5.0 4.5 3.1 4.5 5.5 4.4 3.9 4.2 3.4 5.7 3.6 5.5 4.9 3.3 6.3 5.4
##  [289] 2.9 4.9 4.7 4.6 3.9 3.5 4.4 5.0 3.7 2.6 5.3 3.8 4.6 3.9 3.6 5.7 4.8 4.2
##  [307] 4.6 4.8 4.7 5.8 3.2 2.6 4.5 5.1 3.2 4.1 3.0 5.2 3.1 5.9 4.4 5.7 4.5 4.3
##  [325] 5.9 4.9 5.9 3.5 4.4 5.1 4.1 3.9 6.3 3.6 3.5 4.4 4.0 6.5 4.7 4.1 5.3 5.7
##  [343] 5.1 3.6 2.7 4.4 6.2 4.7 5.9 4.3 3.9 5.5 5.7 4.5 5.9 4.9 5.7 3.7 3.7 4.0
##  [361] 3.0 5.0 7.3 4.4 6.0 4.7 4.1 5.5 5.4 3.8 3.7 4.3 4.4 5.0 3.6 3.8 4.7 3.9
##  [379] 5.6 5.6 4.0 4.9 4.4 4.2 5.9 3.7 3.6 5.3 4.0 5.4 4.8 5.5 3.6 1.8 4.7 2.9
##  [397] 3.7 7.1 4.5 4.6 4.9 2.3 5.4 4.0 3.9 4.2 3.0 3.2 5.3 4.9 3.7 3.3 5.0 2.3
##  [415] 3.7 5.6 3.6 3.9 5.1 3.8 5.6 4.9 4.1 4.2 4.3 4.7 5.7 4.1 5.4 5.6 3.5 3.9
##  [433] 5.3 3.2 4.6 3.8 4.2 5.2 3.8 5.0 3.5 3.4 5.0 4.3 5.9 4.5 3.6 5.4 5.6 5.9
##  [451] 4.7 5.0 4.3 3.5 4.3 2.6 4.7 5.2 4.3 6.4 4.3 4.5 5.7 4.4 3.7 3.6 3.4 3.0
##  [469] 5.7 2.8 4.5 5.6 3.7 5.5 3.6 4.7 2.8 4.7 5.4 6.8 4.5 3.2 4.4 5.1 5.0 5.0
##  [487] 4.7 4.6 3.6 5.3 3.8 3.9 2.8 4.9 4.0 4.0 4.4 4.4 4.8 4.2 3.4 3.5 3.6 4.2
##  [505] 5.9 5.0 4.3 3.7 4.3 3.4 3.9 4.0 5.1 4.0 4.0 4.5 6.2 4.7 4.4 4.5 3.9 4.3
##  [523] 5.7 5.1 2.8 4.7 2.9 4.2 5.4 4.3 4.0 3.4 4.3 2.2 4.7 5.0 5.4 3.7 3.4 3.9
##  [541] 4.9 3.7 3.9 5.9 5.6 5.2 5.1 5.3 5.0 3.4 5.0 3.5 4.5 4.6 5.3 6.1 3.2 5.7
##  [559] 4.1 4.9 5.2 5.8 4.8 5.6 4.8 4.5 3.4 5.2 4.9 3.3 4.3 6.6 2.9 4.8 5.0 4.8
##  [577] 2.9 3.8 4.0 4.6 5.4 2.7 5.5 3.2 5.2 3.7 3.5 4.5 4.7 4.4 2.7 4.0 4.0 3.9
##  [595] 4.4 4.3 3.9 6.0 4.5 5.8 6.2 4.9 4.2 5.3 3.1 5.2 4.9 4.3 4.3 5.2 5.8 4.8
##  [613] 6.0 4.5 4.8 2.6 3.8 5.4 4.6 4.9 1.9 5.1 4.9 4.8 4.7 4.4 4.0 3.2 3.6 3.9
##  [631] 6.0 3.5 4.2 5.1 2.8 5.0 3.6 3.6 3.3 4.7 4.4 3.7 3.5 4.9 5.3 4.1 5.6 3.7
##  [649] 4.4 5.4 3.3 4.8 4.5 5.7 4.0 5.4 4.0 4.1 3.4 3.6 4.6 5.0 6.2 3.0 3.3 3.7
##  [667] 5.3 4.9 4.5 4.4 3.9 5.0 4.2 6.0 3.9 3.6 3.9 4.1 4.5 5.6 4.6 5.0 4.2 5.1
##  [685] 3.8 4.9 5.0 5.1 4.9 5.0 5.3 5.5 4.8 4.4 4.3 4.3 5.3 4.3 4.4 4.8 5.9 4.6
##  [703] 4.1 3.8 6.6 3.2 4.4 5.5 5.3 4.7 4.9 4.7 4.9 3.1 4.2 4.8 3.2 5.2 3.7 4.4
##  [721] 3.6 3.2 4.5 5.6 5.0 4.8 6.4 4.1 3.4 6.6 4.8 5.2 5.2 4.9 4.2 5.9 3.7 3.1
##  [739] 4.4 3.5 2.5 6.4 5.7 4.2 5.4 5.1 4.2 4.8 4.1 5.9 5.5 4.7 5.5 2.5 3.9 4.2
##  [757] 4.0 3.9 4.0 3.5 3.8 4.4 4.5 3.2 5.7 6.1 5.1 3.9 3.8 3.6 3.9 4.6 4.4 4.6
##  [775] 5.4 4.2 6.0 4.4 2.8 3.8 5.2 4.2 5.4 5.0 2.5 4.9 4.4 6.9 4.5 3.7 4.8 4.9
##  [793] 5.6 2.7 5.0 3.1 5.1 4.3 2.7 4.4 4.6 4.6 5.0 5.3 5.1 5.2 4.2 6.0 5.1 4.5
##  [811] 5.1 4.9 3.3 3.7 3.7 4.2 3.7 2.5 4.6 4.4 4.4 4.8 4.8 4.7 3.6 5.5 5.3 2.7
##  [829] 3.7 4.8 4.2 4.2 2.9 5.1 6.5 5.1 4.6 4.1 5.3 3.4 5.0 3.5 4.3 5.0 4.7 5.4
##  [847] 3.9 5.9 3.8 4.8 4.3 5.4 5.4 3.1 4.9 4.7 3.8 4.5 4.1 4.6 5.0 4.9 4.8 4.4
##  [865] 3.4 6.7 4.7 3.9 4.8 4.8 4.0 5.4 4.0 4.2 7.6 4.4 4.9 5.0 5.0 2.9 4.8 4.8
##  [883] 5.6 4.8 2.9 2.8 4.2 4.6 3.7 5.0 4.5 7.1 3.9 6.6 2.3 5.0 3.1 4.5 3.9 4.6
##  [901] 5.1 4.0 4.1 4.8 4.2 5.2 4.5 5.9 2.8 4.1 5.4 5.1 4.8 5.0 4.7 4.6 4.4 3.7
##  [919] 5.4 4.4 5.4 5.8 4.0 4.7 3.8 5.8 4.5 4.5 5.6 4.0 5.0 3.2 4.7 4.6 3.3 5.6
##  [937] 5.2 2.6 4.7 3.7 4.7 5.0 5.4 5.0 3.6 4.8 3.5 5.1 3.7 4.7 6.2 5.5 5.3 5.0
##  [955] 4.7 6.0 4.4 4.5 5.9 5.2 6.3 4.9 3.8 4.1 2.7 3.0 4.9 2.0 3.4 6.3 5.0 3.3
##  [973] 4.0 4.3 5.5 3.9 3.3 5.3 3.0 4.1 5.0 4.5 3.6 5.9 5.0 4.2 4.4 4.5 5.0 4.4
##  [991] 4.8 3.9 5.4 4.5 5.0 4.6 2.9 4.7 4.7 4.5
## 
## $func.thetastar
## [1] -0.0078
## 
## $jack.boot.val
##  [1]  0.50277008  0.42478134  0.25131195  0.06695157  0.04156627 -0.06257143
##  [7] -0.17674419 -0.30492308 -0.43267606 -0.51908832
## 
## $jack.boot.se
## [1] 0.9872806
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
##    [1] 3.7 5.6 6.5 3.0 4.2 5.9 3.3 4.3 4.8 5.1 4.1 5.2 5.5 3.8 3.5 5.5 5.0 4.6
##   [19] 5.6 3.4 4.5 4.4 5.3 4.2 4.8 6.2 3.9 5.3 4.2 5.4 5.7 5.7 4.2 3.9 3.8 5.0
##   [37] 4.9 4.4 5.3 4.7 3.8 5.4 6.0 3.9 5.3 4.7 4.9 3.4 4.8 4.1 6.2 5.5 4.3 5.4
##   [55] 3.1 3.4 3.4 6.0 4.2 3.7 4.1 4.9 3.1 4.6 5.3 4.5 4.8 4.4 3.8 5.0 5.0 5.2
##   [73] 3.7 4.0 3.1 4.5 4.1 5.4 5.9 4.4 6.3 4.0 6.9 4.3 4.9 2.9 3.6 5.2 4.5 5.3
##   [91] 3.4 4.0 2.4 5.1 5.4 3.9 3.9 4.1 6.1 4.6 3.0 3.4 3.2 3.8 4.0 3.0 4.9 5.7
##  [109] 5.1 4.4 5.4 3.4 3.0 6.5 6.0 4.5 5.4 3.4 5.0 3.3 5.3 5.1 4.3 5.2 5.8 3.8
##  [127] 5.0 4.4 3.9 5.2 6.0 3.9 3.8 4.3 5.6 3.8 4.8 5.5 5.9 4.4 6.0 3.5 5.4 6.1
##  [145] 5.7 3.5 3.6 4.5 4.2 2.8 6.2 2.7 3.1 4.4 5.8 4.7 5.8 3.5 4.3 5.4 5.8 3.1
##  [163] 2.9 3.7 5.3 4.7 4.8 5.1 3.6 4.1 4.0 4.1 3.5 4.2 5.2 4.9 4.8 3.1 4.9 4.8
##  [181] 3.2 3.7 3.9 3.2 4.2 4.4 4.4 2.3 5.0 4.9 5.4 4.2 3.2 4.9 5.7 4.6 5.3 5.5
##  [199] 3.0 4.1 4.4 4.8 3.7 5.2 4.2 4.8 4.4 4.4 3.8 4.5 5.1 4.6 2.7 5.2 4.6 3.4
##  [217] 5.0 4.4 4.8 5.8 4.8 4.4 4.0 5.1 5.4 3.5 4.9 3.7 3.3 6.2 4.1 4.5 3.4 5.3
##  [235] 6.2 4.7 3.2 3.9 4.7 4.2 3.8 2.7 4.2 4.2 6.2 5.2 4.2 5.2 4.5 5.8 4.8 5.2
##  [253] 4.4 4.2 5.9 4.0 4.0 4.3 5.0 1.7 4.7 6.3 3.7 4.2 4.5 3.7 4.8 4.9 4.2 5.1
##  [271] 3.2 4.8 6.1 4.7 4.2 4.0 5.4 4.7 5.8 5.8 2.9 4.5 4.7 4.4 4.6 6.1 5.7 3.0
##  [289] 4.0 6.1 5.7 4.0 5.2 3.6 4.2 4.8 3.4 4.4 4.4 4.4 5.3 5.2 4.1 4.8 4.1 2.8
##  [307] 5.1 4.5 4.7 4.3 4.3 3.3 5.0 4.1 5.0 4.1 5.3 4.4 2.5 3.5 4.0 4.3 4.8 5.1
##  [325] 5.2 3.0 4.9 3.8 4.5 3.6 5.4 3.0 4.1 5.1 5.1 4.4 4.4 4.0 4.2 4.6 5.5 5.2
##  [343] 3.2 6.0 6.5 2.3 5.6 4.9 3.7 4.9 4.0 6.6 5.2 3.2 4.5 7.5 5.4 4.8 5.8 5.0
##  [361] 5.2 3.9 4.0 5.2 5.4 5.1 3.9 3.7 4.7 6.4 3.9 4.1 2.9 3.5 5.0 3.6 5.1 3.9
##  [379] 6.1 4.8 5.5 3.8 5.1 4.9 3.1 3.7 4.9 5.7 4.0 2.9 3.7 4.6 3.6 4.9 4.9 5.3
##  [397] 3.3 5.0 4.8 6.5 4.4 5.6 4.2 3.5 4.5 4.2 2.3 3.3 3.6 3.4 2.5 4.4 2.6 3.9
##  [415] 5.3 3.5 3.5 4.7 3.5 4.9 3.3 3.9 5.2 4.4 5.1 4.0 5.5 3.6 3.9 5.2 3.8 4.1
##  [433] 5.6 5.8 4.4 4.1 5.1 4.7 4.8 6.4 4.0 4.4 4.5 3.3 3.9 4.4 5.2 5.2 4.2 5.1
##  [451] 3.1 4.7 4.4 4.4 2.8 4.4 4.6 3.4 5.8 4.6 4.7 2.4 5.7 4.6 4.7 6.3 4.6 4.1
##  [469] 3.7 2.3 4.8 3.7 4.3 4.6 5.2 4.1 4.3 3.3 3.2 3.2 3.0 5.6 4.1 3.3 3.9 4.9
##  [487] 6.2 3.9 4.5 4.0 6.2 4.2 5.3 4.8 3.1 5.6 5.3 3.3 4.8 2.4 5.8 6.0 4.5 3.8
##  [505] 4.9 3.5 5.1 5.8 5.0 5.2 4.3 7.1 3.9 5.1 2.9 3.6 3.7 5.3 5.8 3.0 3.7 5.1
##  [523] 6.0 3.6 4.5 3.5 3.5 2.6 4.5 3.9 3.8 3.0 4.2 4.7 3.6 5.0 2.9 4.1 5.3 4.5
##  [541] 3.3 2.8 4.4 3.4 2.2 3.8 5.0 5.3 3.6 5.1 5.2 5.1 4.4 3.8 5.0 3.5 5.3 5.6
##  [559] 4.7 5.1 4.2 4.4 3.9 2.0 4.5 4.7 4.7 4.3 4.4 5.0 4.1 4.7 4.6 4.9 5.4 5.3
##  [577] 4.2 5.1 4.8 4.6 3.3 4.3 4.0 3.5 4.6 4.9 4.2 5.6 3.9 3.6 4.0 3.8 3.1 6.0
##  [595] 6.5 5.6 4.8 5.4 5.9 5.2 4.6 5.0 4.9 5.0 2.8 4.8 4.2 3.9 2.7 4.2 4.7 5.0
##  [613] 3.5 4.5 6.1 4.8 4.2 4.8 4.3 4.0 5.8 4.2 4.6 4.9 4.9 4.2 4.4 5.1 3.2 3.1
##  [631] 4.8 5.2 3.9 3.5 4.2 5.2 4.9 5.6 4.5 3.5 4.1 5.7 4.2 5.7 5.1 4.3 3.2 4.1
##  [649] 3.2 4.4 4.4 4.0 3.8 4.7 5.4 6.2 3.9 3.9 3.9 4.0 5.0 4.6 4.5 4.0 4.0 5.3
##  [667] 4.3 4.4 5.4 5.6 2.8 4.6 6.2 4.6 4.3 4.5 5.5 4.2 5.6 5.2 4.5 4.5 4.9 3.9
##  [685] 5.5 3.9 4.0 4.7 5.0 3.9 3.2 4.8 5.5 2.9 3.7 4.5 5.1 4.2 4.4 4.8 4.1 3.9
##  [703] 3.2 2.4 3.9 4.7 5.0 5.3 3.4 5.7 4.4 3.7 4.6 4.4 4.8 5.0 4.2 4.1 5.1 4.4
##  [721] 4.6 6.0 3.5 4.1 3.4 5.5 4.3 4.1 4.8 4.6 3.1 5.1 5.6 4.9 3.7 4.3 5.6 4.4
##  [739] 2.7 3.9 4.3 2.7 7.5 5.4 4.9 3.7 2.7 4.1 4.7 4.3 3.3 4.4 3.9 4.7 4.4 5.1
##  [757] 6.1 5.2 4.5 4.0 6.4 3.9 4.7 5.1 4.7 5.1 4.0 3.2 5.7 4.7 5.2 4.9 5.1 3.5
##  [775] 3.5 3.3 4.8 4.6 4.6 4.1 4.5 5.6 5.6 5.8 4.3 5.8 5.9 6.2 4.1 5.5 4.5 5.2
##  [793] 6.1 4.1 3.4 3.5 3.5 3.7 4.8 3.9 4.4 6.4 2.7 5.2 4.8 3.5 5.7 3.9 4.6 5.4
##  [811] 5.7 4.2 4.9 5.8 6.1 6.6 4.2 4.3 4.9 5.2 4.6 4.8 4.6 5.9 6.3 5.3 3.8 5.6
##  [829] 4.3 4.3 3.3 3.4 5.4 2.9 5.4 3.7 5.8 5.1 3.6 3.2 5.5 3.0 3.6 4.5 4.7 5.2
##  [847] 3.9 4.5 4.4 5.8 3.0 6.4 3.7 4.3 2.7 3.8 5.0 5.3 3.5 4.0 5.0 3.0 3.9 5.2
##  [865] 5.0 6.5 3.9 3.9 4.9 3.2 4.5 4.6 4.7 4.5 5.2 3.8 4.5 4.2 4.8 5.0 4.9 3.7
##  [883] 3.5 4.8 4.7 4.8 3.9 5.0 4.1 4.0 5.0 2.0 5.9 4.5 4.8 6.0 4.8 4.5 5.4 4.3
##  [901] 5.0 3.9 5.2 5.2 3.9 3.2 6.3 5.2 5.0 5.4 3.7 5.6 3.1 6.4 5.6 4.4 3.9 4.2
##  [919] 4.3 3.2 4.3 4.3 3.8 5.1 5.1 3.2 5.5 3.4 3.6 4.8 4.2 3.0 5.6 4.8 5.1 5.3
##  [937] 3.4 4.1 5.9 5.5 4.9 5.2 2.7 6.1 4.3 4.4 5.8 5.0 5.7 3.3 4.1 4.4 2.6 4.6
##  [955] 4.8 2.7 4.2 4.7 4.3 3.8 5.5 5.7 3.4 3.2 5.0 4.1 4.6 6.1 4.2 5.2 4.0 4.7
##  [973] 4.6 5.6 3.3 5.7 5.1 3.0 3.8 6.1 3.7 4.4 4.3 5.3 3.2 5.8 3.3 4.3 5.4 3.8
##  [991] 4.3 2.7 3.4 6.0 5.7 5.7 3.5 4.8 3.2 3.5
## 
## $func.thetastar
## 72% 
## 5.1 
## 
## $jack.boot.val
##  [1] 5.4 5.4 5.4 5.2 5.1 5.1 4.9 4.8 4.6 4.5
## 
## $jack.boot.se
## [1] 0.9410632
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
## [1] -0.3818607
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
##   3.134803   4.644475 
##  (1.334074) (2.143558)
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
## [1] -0.1878791  1.2456753 -0.2588819  0.2110482  0.5525791  1.8643518
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
##    [1] -0.372166152  1.179379955 -1.123521438 -0.452578012  0.020822591
##    [6] -0.307060848  0.621812598 -0.446480171 -0.240853499  0.345357683
##   [11]  0.030501661 -0.261093673 -0.455938883 -1.084614486 -1.765921461
##   [16] -0.580281696 -0.099755488 -0.869992696 -0.636541096 -0.738832496
##   [21] -0.098840237 -0.103153567  0.135318003 -0.961698646  0.054525575
##   [26] -0.201835454 -0.896238471  0.739875178 -0.949000684 -0.319959313
##   [31] -0.548959486 -0.514925439 -0.914646827 -0.269864135 -1.000295009
##   [36] -0.676936504 -1.000715689 -0.594386423 -1.286131710 -0.940415226
##   [41]  0.518260098 -0.381694953 -0.458443338 -1.358708012  0.053239104
##   [46] -0.136547068  0.882373327 -1.230623144 -1.404601073 -0.581736883
##   [51]  0.482149339  0.097410424 -1.160390286  0.004889623 -1.253499879
##   [56] -0.435297843  0.147849356 -1.392089348 -0.066830689 -0.220605792
##   [61] -0.162393987 -1.230623144 -0.958458241  0.206738832 -0.729895315
##   [66] -0.949854485 -0.135570845 -1.393655885  0.131852833 -2.338661053
##   [71] -0.355437613 -0.221803361 -0.751228148 -0.637839191 -0.762821605
##   [76]  0.217904279 -0.222733359 -1.443041904 -0.803536371 -0.444893085
##   [81] -0.179793763 -1.331217218  0.021702831 -0.826739104 -0.761988558
##   [86] -0.092582053 -0.358896435 -0.399533245 -0.479685269  1.232741011
##   [91]  0.223696905 -0.313804736 -0.170493227  0.090354216 -0.704848264
##   [96] -1.150533212 -0.676936504 -0.509276890 -0.627044166 -0.607903258
##  [101] -1.131120240 -0.559270711 -0.859711761  0.618696362  0.018181339
##  [106] -0.252278825  0.393734283 -0.100931597  0.018343049 -0.559237817
##  [111] -0.662892676 -0.584565129 -0.128621718 -0.512788952 -0.067188599
##  [116] -0.529898893 -0.008324578 -0.949882955 -0.090844075 -0.102795633
##  [121]  0.050225282 -0.427375224  0.212448971 -0.428394322 -0.978128261
##  [126] -0.572136454 -1.373529367 -1.114040385  0.377419744 -0.288523851
##  [131] -0.719907806  0.033294572 -0.721276785 -1.043786436 -0.140148985
##  [136] -0.146491808 -0.272608435  0.262818433  0.154081280 -0.906753323
##  [141]  0.295838346 -0.219232800 -0.933358839 -0.077206328 -0.974138261
##  [146] -0.064939535 -0.954085140 -0.341550225  0.016775563 -0.296908482
##  [151] -0.831085141 -1.037315248 -0.901795250 -0.471921371 -0.847252891
##  [156]  0.251826459 -0.809188497 -0.135936221  0.539728060 -0.601302122
##  [161] -0.592977470 -0.227225034 -1.718953069 -0.170180618  0.452542563
##  [166]  0.030458530 -0.835533077 -1.151224984 -1.103723255 -0.198224130
##  [171] -0.965946722 -0.173584874 -0.130572555 -0.511212775 -0.867810178
##  [176] -0.179939955 -0.194172580 -0.550898747 -0.592524351 -0.747235081
##  [181] -0.320817582 -0.036563519 -0.383019196  0.521606778 -0.196712284
##  [186]  0.525875513  0.050980289  0.614308256  0.208128978 -0.450856141
##  [191]  0.237186531 -0.713175857 -0.310372432  0.405374158 -0.336063651
##  [196] -0.197138780 -0.939094371 -0.725786391  0.197810085 -0.046080291
##  [201] -0.635060375  0.007360053 -1.052237401  0.233449121 -0.932739346
##  [206] -0.411536432  0.232124616  0.387430878 -0.240094045 -0.306735643
##  [211] -1.031205316 -0.695713301 -1.012594756 -0.338964615 -0.399631702
##  [216] -0.749454993 -0.634727376 -0.241270018 -0.095344646 -0.363772247
##  [221] -0.317758540 -1.018413459 -0.521182564  0.200853530  0.151343769
##  [226] -0.024595972 -0.696885343 -1.025189141 -0.503269241 -0.070549725
##  [231] -0.544682559 -0.768161642 -1.161248206 -0.569198375 -0.643913426
##  [236] -0.713528008 -1.158813411  0.509938106 -0.048969980 -0.516233722
##  [241] -0.832737599  0.149995972 -0.435301295 -1.746295361 -0.181581425
##  [246] -0.859764742 -0.670105495 -0.185548911  0.261807651 -0.653630847
##  [251]  0.418875918 -0.637680737 -0.655697225 -0.405241039 -0.097703310
##  [256] -0.334814020 -0.714603300  0.387715394 -0.117153419 -0.949895648
##  [261] -0.146565069  0.428076177 -0.356137532  1.244004872  0.493502687
##  [266] -0.061187495 -0.856190668 -0.512953873 -0.370662772 -0.310747951
##  [271] -0.734840395  0.569345381 -0.907485695 -0.047684355 -0.647531183
##  [276]  0.082293945 -0.806493008 -0.138035083 -0.245283481  1.450351120
##  [281]  0.343087016 -0.587732736  0.015698485  0.736904726 -0.317441144
##  [286] -1.361301443 -0.282100574 -0.598482797 -0.384179787  0.444697871
##  [291] -0.968347035  0.386093725  0.313810018 -0.740553816  1.029666480
##  [296] -0.454266731 -1.016852190 -0.332707213 -0.171644897 -0.117459771
##  [301] -1.310254367 -0.119453636 -0.088841437 -0.853860868 -0.358316189
##  [306]  0.105295750 -0.011244481 -0.839572087 -0.370399437 -0.150201843
##  [311] -0.416503507 -0.685785853 -0.143136813 -0.955813907 -1.277581486
##  [316]  0.432165995 -0.027613417 -0.238220632 -1.001335880 -0.206083249
##  [321]  0.348470508 -1.299876503 -0.550228467  0.359067264 -0.205751864
##  [326] -0.640757123 -0.059707033 -0.974507340 -0.503581538 -0.276073928
##  [331] -0.280942990 -0.327956652 -0.506993773 -0.305902159 -0.257801257
##  [336] -0.211637081 -1.191003160 -0.980200179 -0.378384084  0.266944632
##  [341] -1.010833451 -1.029229205  0.234008229 -0.505846470 -0.242791566
##  [346] -1.076713017 -0.855762239 -0.238774618 -1.360615299 -0.252137779
##  [351]  0.021206831 -0.128504420 -0.246755036 -0.842206660 -0.319906376
##  [356]  0.086689092  0.448619461 -0.401499779 -1.188116490  0.299977444
##  [361]  0.061682626  0.539664146 -0.242748949 -0.052184900  0.201040540
##  [366] -0.769465795 -0.506993773 -0.453875420 -0.154084276  0.133695620
##  [371] -0.801009521  0.773551029 -0.196892967  0.087006348 -1.187904256
##  [376] -0.759465385 -0.501996892  0.219987729  0.239909038 -0.454268904
##  [381] -0.500015698 -0.518377787 -0.676745763 -0.378405669 -0.677431578
##  [386]  0.180045617  0.088485414 -0.074758384 -0.759465385 -0.672441020
##  [391] -0.184454714 -0.078996740 -0.639236221 -0.267160126 -0.217905952
##  [396]  0.204298459 -0.628110735  0.145170084 -0.016598128 -0.547564327
##  [401] -0.992873881 -0.272650255 -0.379792350 -0.312752215 -0.471969909
##  [406] -0.655892475 -0.091770975  0.060857736  0.384134000 -1.410870100
##  [411] -0.528996027  0.245897926 -0.384805666 -0.887101501  0.156844916
##  [416] -0.711202474 -0.252433832 -1.049211337 -0.049260026 -0.737764236
##  [421]  0.238535075  0.041693107 -0.018824895 -1.085156033 -1.260035895
##  [426] -1.216880509  0.140401306 -0.871834137 -0.674884385 -0.279923589
##  [431] -0.554434555 -1.050864164  0.189103473 -0.337812481 -0.607903258
##  [436] -0.285322380  0.903687387  0.390079766 -0.706227187 -0.490516629
##  [441]  0.343085556 -0.850190843 -1.463906179 -0.122155795 -1.362761363
##  [446] -1.009066336 -0.401236674 -0.634832601 -0.121590047 -0.336321014
##  [451] -0.931491198 -0.710589011 -0.850141329 -0.884099417  0.157172783
##  [456]  0.014902411 -1.007814538 -0.949000684 -0.992185275 -0.157983120
##  [461] -0.290676516  0.253323330 -0.401589322 -0.032353586 -0.511001029
##  [466] -0.762168735 -0.583491640 -0.162355319 -0.486870263 -1.405794115
##  [471]  0.072043345  0.236219791 -0.357263420 -0.350014562  0.015385452
##  [476] -0.284755171  0.011580125 -0.053411134 -1.052392200 -0.403472954
##  [481] -1.619973306 -0.993023592  0.536302284 -0.144079575  0.601981701
##  [486] -0.367370589 -0.595306020 -0.947270182 -0.286389676 -0.364374838
##  [491] -0.444164049 -1.564100464 -0.252520633 -1.484512306  0.616176025
##  [496] -0.223156059 -0.105266910 -0.303486045 -0.190396274 -0.640687768
##  [501] -0.719389094 -1.037850915 -0.351833734 -0.936787275 -0.500165452
##  [506]  0.095601408  0.444549207 -0.183097060 -0.839999571 -0.998689419
##  [511] -0.226803790 -1.406546333 -0.452272317 -0.064744066  0.214224296
##  [516] -0.771240420 -0.598039724 -0.886875545 -0.015088306 -1.320133825
##  [521] -0.687861186  0.384000608 -0.798756721 -0.332785659  0.141074067
##  [526] -0.533278787  0.265775909 -0.878725027 -0.180942421 -0.476531129
##  [531]  0.338965012  0.091103711  0.125022455 -0.194389453  0.067721780
##  [536] -0.807601748  0.233479546 -0.474704474 -0.307718044 -1.043385745
##  [541] -1.281550281 -1.231596131 -0.252881230 -0.329915819 -0.627381324
##  [546] -2.490144731 -0.556487104 -0.512730724 -0.304158810 -0.243947788
##  [551] -0.319986480  0.239174930  0.407703093 -0.536103049 -1.034078206
##  [556] -1.151964090 -0.496773686 -0.588261644  0.670241402 -0.689201907
##  [561] -0.338095599 -0.376182730 -1.104361154 -2.648939477  0.001540554
##  [566] -0.320247345 -1.486355734 -0.044219826 -0.474838668 -0.121515268
##  [571] -0.589842114  0.057328733  0.401344761 -0.147291416 -0.600913731
##  [576] -0.089320008 -0.158529012 -0.305866130 -0.688836306 -0.738779013
##  [581] -0.683734168 -0.454024162 -0.811917982 -0.645191605 -0.635771828
##  [586] -1.111373289 -0.712528158 -0.353432702 -0.643774852 -1.089128161
##  [591]  0.141234273 -1.069401594 -0.936251155 -0.395756976 -1.399466950
##  [596] -1.364925923 -0.387408970 -0.338964615 -1.412295056 -1.029064700
##  [601] -0.688150653 -0.377662915 -0.243528639 -1.424083303 -0.889971115
##  [606]  0.338355242  0.057717463 -1.166100394 -0.124260639  0.010330628
##  [611] -0.725381145 -0.303213631 -0.711465591 -0.875225534 -1.330948591
##  [616]  0.305317679 -0.654242395 -0.024978506 -0.641147854 -0.284306610
##  [621] -0.142661828 -0.059930053  0.023010447 -0.366452786 -1.778818728
##  [626] -0.398520653 -0.240853499 -1.332368571 -0.694368599 -0.182818216
##  [631] -0.305909644 -0.833143549 -0.871300017  0.092426949  0.455699382
##  [636] -0.307068266  0.571098335 -0.232571328 -0.902486803 -0.205630496
##  [641] -0.336647271 -0.883132666 -0.627371489 -0.213555340 -1.170738575
##  [646] -0.608656400 -0.386335009  0.401521066  0.376853399 -0.794104735
##  [651] -1.125615695 -0.307466679 -1.161126344 -1.076293968 -1.007814538
##  [656]  0.355679491 -0.205616777 -0.741590586 -0.038684795  0.741732702
##  [661] -0.655606057 -1.207382613 -0.878814391  0.557994637  0.058682728
##  [666]  0.525976861  0.728604328 -0.247126564 -0.378336272  0.200546967
##  [671]  0.064685800 -0.438560149 -0.350984294 -0.159712474 -0.459256649
##  [676] -0.378522401 -0.960463179 -0.129365478 -0.327491674 -0.275345286
##  [681] -0.063309270 -0.561399568  0.408317267 -0.328799128  0.018313227
##  [686] -0.814965961  0.373533234 -0.829925788 -1.282868008 -0.824757683
##  [691]  0.173284687 -0.009217245 -0.504177216  0.500080870  0.066281946
##  [696] -0.366883081 -0.330628702 -0.289643368 -0.623861648 -0.512666861
##  [701] -0.336063651  0.508890993  0.247903069 -1.378529794 -0.840963587
##  [706] -1.120353239  0.319393508  0.011141584 -0.019010732 -0.590701510
##  [711] -0.115061402 -0.587038495  0.489055284  0.102608884  0.267893022
##  [716] -0.378873149 -1.130160618 -0.961142466 -0.995733660 -1.116768943
##  [721]  0.498247098 -0.614435232 -0.548911850 -0.388554798 -0.586523416
##  [726] -1.171167376  0.068328218 -0.009804434 -0.811612754 -0.024138055
##  [731] -0.569706982 -0.529008414 -0.011736482  0.283450408  0.787648897
##  [736]  0.604391776 -0.840057451 -0.676396055  0.198407128 -0.094199841
##  [741] -0.164323259  0.432758503  0.342598185 -0.284667255 -0.593608882
##  [746] -0.672719182 -0.641262269 -0.769363942  0.194679671  0.026822349
##  [751]  0.469829247 -1.419134367  0.220273299 -0.367343745  0.330209126
##  [756] -0.170493227  1.163381128 -0.189936242  0.159099819  0.122398668
##  [761] -0.549707661 -0.535855585 -0.799501526 -0.854143021 -0.616815843
##  [766] -0.274819481 -1.269026626  0.789299560 -0.601809129 -0.074874403
##  [771] -0.381084985 -0.418830706 -0.485507384 -0.518046366  0.147849356
##  [776]  0.181483624  0.522626323 -0.561392023 -0.544230372 -0.754148161
##  [781] -1.054447165 -0.599792460 -0.408874846 -0.169051776 -0.876457229
##  [786] -0.156839013 -0.175137720  0.762960348 -0.697637928 -0.485775006
##  [791] -0.859324924 -0.909391716  0.765737367 -0.306326260 -0.677437579
##  [796] -1.311244310 -0.919133489 -0.797117507 -0.638008166 -1.012293197
##  [801]  0.298143284 -0.899115024 -1.393470382 -1.040680550 -1.859388830
##  [806]  0.326518178  0.101778506 -0.599963944  0.501502783  0.131780626
##  [811] -0.524343450 -0.562250429 -0.028711783 -0.808634094  0.297340785
##  [816]  0.049372916 -0.921969718 -0.982764178  0.041790691 -0.544699983
##  [821]  0.024141543 -0.783574833  0.731688734 -1.589101160 -0.875968867
##  [826] -0.735505064 -1.293037157 -0.253456063  0.715280459  0.458587733
##  [831]  0.038911180 -0.209187688  0.490245419 -0.124820835 -0.377541757
##  [836] -0.495311005 -0.467162105  0.270742764 -0.342757087 -0.491811004
##  [841] -0.088558812 -0.063171196 -0.524130130  0.254706323  0.678791562
##  [846] -0.367113882 -0.676330593  0.095306381 -0.215679087 -1.291436895
##  [851] -0.172977388  0.037096026 -0.697060761  0.582335355 -0.155857107
##  [856] -1.177236871 -0.289088796  0.014620919  0.433937207 -1.581232354
##  [861]  0.493502687 -0.453156085  0.043635169  0.206659281 -0.794936450
##  [866]  0.024417315 -0.294891729 -0.482645514 -0.276674376 -0.034638950
##  [871]  0.068553550 -1.012778267 -0.842312078 -0.958458241 -1.398919444
##  [876]  0.314234982 -0.581625452  0.094908189 -0.691466845 -0.873919156
##  [881] -0.181275390  0.816411444  0.273729155  0.164567175  0.002321740
##  [886] -0.400244183  0.164821565  0.405374158 -0.210068450 -0.750443479
##  [891] -1.041588825  0.407690628 -0.553937572 -0.153488309 -0.577418864
##  [896] -0.394124689 -0.211165633  0.414412945 -1.354495266 -0.649416134
##  [901] -0.550778138  0.022110195 -0.118759962 -0.306195188 -0.368485159
##  [906] -0.427476206 -0.597412372 -0.696885343 -1.448974036 -0.901425728
##  [911] -1.625784523 -0.181671598  0.170229708 -0.574219555  0.079712928
##  [916]  0.055705319 -0.359130069 -0.555541007 -0.177974657 -0.539050927
##  [921] -0.428140464 -0.509392501 -0.893411285  0.219291034 -0.366745603
##  [926] -0.213434938  0.113906922 -0.249310123 -0.123613211 -0.030994727
##  [931] -0.389125094  0.182741911 -0.003362660 -0.146805140 -0.337949573
##  [936] -0.127987901 -0.438771925 -0.416388435 -0.286072488 -0.214109390
##  [941]  0.237952130 -0.718661245 -0.088841437 -0.382090258 -0.224071588
##  [946]  0.136536368 -0.709046714 -0.317935665 -1.131116211  0.088779358
##  [951] -0.397543357 -0.148016558 -0.661223828  0.233159042 -0.414131844
##  [956]  0.044500500  0.341449063 -0.339176992 -0.400325029 -1.535414918
##  [961] -0.044786218  0.234879521  0.242219052 -0.745684369 -0.528086824
##  [966]  0.031871102 -0.239554308 -0.023697373 -0.346446412 -0.283286027
##  [971] -0.381923489 -1.191361712 -0.572399863 -0.503243838 -0.343156128
##  [976] -0.003085813  0.386433755 -0.418954347 -0.830710596 -0.140078824
##  [981]  0.193342085 -0.212065894  0.192975354 -0.451269346 -0.307118552
##  [986] -0.674768182  0.070411460  0.127496840 -1.109440528 -0.385467495
##  [991]  0.004889623 -0.366379426 -0.408789349  0.352505403  0.232713238
##  [996] -0.317854065  0.144655322 -0.059630403 -0.165975335 -0.600854847
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
##   0.67494163   0.31667891 
##  (0.10014266) (0.07080919)
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
## [1] -0.01628463 -0.91534406  1.36366857 -0.57888497 -0.82306194 -0.29458382
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
## [1] -0.0345
```

```r
sd(results2$thetastar)  #the standard deviation of the theta stars is the SE of the statistic (in this case, the mean)
```

```
## [1] 0.9176174
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
## t1*      4.5 -0.05135135    0.917387
```

One of the main advantages to the 'boot' package over the 'bootstrap' package is the nicer formatting of the output.

Going back to our original code, lets see how we could reproduce all of these numbers:


```r
table(sample(x,size=length(x),replace=T))
```

```
## 
## 0 2 4 5 6 7 8 9 
## 1 1 1 2 2 1 1 1
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
## [1] 0.0077
```

```r
se.boot
```

```
## [1] 0.9054339
```

Why do our numbers not agree exactly with those of the boot package? This is because our estimates of bias and standard error are just estimates, and they carry with them their own uncertainties. That is one of the reasons we might bother doing jackknife-after-bootstrap.

The 'boot' package has a LOT of functionality. If we have time, we will come back to some of these more complex functions later in the semester as we cover topics like regression and glm.

