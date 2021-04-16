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
## 0 1 2 3 5 6 9 
## 1 1 1 1 1 1 4
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
## [1] -0.0241
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
## [1] 2.673985
```

```r
UL.boot
```

```
## [1] 6.277815
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
##    [1] 3.0 4.8 3.2 3.1 4.7 4.7 3.1 4.5 4.0 4.1 6.2 4.7 3.2 3.9 4.2 6.7 3.3 4.6
##   [19] 5.0 5.5 6.0 4.7 5.3 4.3 5.3 4.5 3.5 4.8 2.8 4.5 3.4 3.2 5.2 4.1 5.9 3.7
##   [37] 5.3 6.2 4.4 4.4 5.4 5.2 3.7 5.4 3.9 4.4 5.3 4.8 5.0 4.4 3.0 5.4 4.0 5.2
##   [55] 5.3 5.4 3.5 5.3 4.6 3.7 4.2 3.1 5.0 4.3 4.5 5.0 5.3 4.1 3.9 2.0 5.0 5.4
##   [73] 5.1 3.4 5.1 3.3 4.0 3.6 3.3 4.6 5.6 6.0 5.9 4.9 5.8 4.1 3.2 3.9 3.5 4.7
##   [91] 5.8 5.2 4.5 6.1 4.5 4.5 6.9 3.8 5.8 5.1 3.5 4.2 6.0 4.8 5.2 3.7 5.2 3.3
##  [109] 4.7 6.9 4.9 5.0 5.7 5.2 4.4 4.0 5.1 6.1 4.1 4.6 2.5 6.2 4.7 4.2 4.2 5.0
##  [127] 4.1 3.4 4.3 5.4 3.4 4.4 4.4 6.3 5.7 3.9 5.2 5.4 4.2 4.7 5.1 5.0 3.8 3.7
##  [145] 3.6 4.3 3.6 5.1 5.8 4.8 4.2 4.8 3.9 3.6 4.3 4.5 6.3 4.1 4.7 3.7 4.6 5.1
##  [163] 4.6 5.2 4.3 5.0 3.7 3.8 5.4 4.1 4.5 2.5 7.4 3.4 5.2 5.2 3.3 4.9 5.1 2.6
##  [181] 5.0 4.8 4.5 3.8 3.5 4.5 3.6 3.9 4.4 4.4 4.4 4.7 4.8 4.7 1.9 2.6 3.7 4.0
##  [199] 6.1 3.9 4.5 4.3 4.8 5.8 4.8 5.4 4.3 4.7 5.6 5.3 4.2 5.0 4.5 2.9 3.4 4.9
##  [217] 3.8 5.7 3.9 6.0 5.8 6.6 5.2 4.6 3.9 3.9 3.5 5.2 5.2 5.3 4.4 2.7 4.3 6.0
##  [235] 5.6 5.0 5.1 5.6 5.2 6.2 6.1 3.4 3.9 4.5 6.0 5.2 4.2 4.9 3.9 3.7 4.4 3.5
##  [253] 3.6 5.5 3.9 3.3 4.3 5.6 4.5 4.1 4.0 5.2 3.6 4.4 3.8 4.2 5.6 4.0 4.3 5.3
##  [271] 5.2 5.2 5.8 4.2 3.9 4.3 4.6 5.0 3.1 3.7 5.7 4.0 4.0 3.8 3.3 4.6 4.7 5.5
##  [289] 3.3 5.2 4.5 5.0 4.5 3.3 6.5 4.9 4.5 4.9 3.9 5.7 4.0 3.8 3.7 4.0 4.5 5.4
##  [307] 5.2 4.9 3.9 5.8 3.2 3.7 4.8 5.9 4.1 4.4 3.0 4.4 6.5 5.7 5.4 4.7 3.8 3.9
##  [325] 4.3 5.2 5.3 5.0 3.7 5.2 6.0 5.5 3.2 4.6 3.3 4.8 5.5 4.5 3.9 4.3 3.3 2.5
##  [343] 5.0 3.2 5.2 2.9 5.8 3.4 4.5 5.9 3.7 4.3 3.8 5.3 4.7 4.8 3.9 4.2 4.3 5.3
##  [361] 5.7 5.5 5.3 5.4 3.4 3.4 4.9 5.1 4.9 3.5 5.6 4.4 4.4 3.7 4.0 7.0 5.4 3.1
##  [379] 4.5 5.7 5.0 5.0 4.0 4.0 5.6 5.1 4.8 4.6 4.0 4.2 4.5 4.8 5.8 5.6 4.1 4.5
##  [397] 6.5 3.9 4.8 3.2 4.6 3.4 4.9 5.1 5.8 3.9 5.2 3.8 6.7 3.0 4.5 2.7 5.3 6.0
##  [415] 3.3 3.8 5.0 4.6 4.3 5.1 5.5 3.4 3.7 4.5 5.3 4.3 5.1 2.8 5.6 4.6 4.3 5.4
##  [433] 3.3 4.5 2.8 6.4 3.4 3.9 4.1 4.9 4.5 3.8 5.1 4.2 5.5 3.2 6.0 4.5 7.0 5.6
##  [451] 3.8 3.1 3.7 4.5 4.9 5.6 3.1 4.2 3.3 4.4 5.4 5.6 3.7 4.5 5.7 3.2 3.2 4.0
##  [469] 4.6 4.7 2.9 5.4 3.5 3.8 5.4 3.8 6.1 4.1 5.4 3.9 6.0 4.4 2.7 3.9 4.3 2.9
##  [487] 5.3 5.6 4.1 4.4 5.6 6.0 5.0 5.6 5.8 4.6 5.4 5.8 4.9 3.2 4.2 3.3 4.1 5.8
##  [505] 3.9 4.9 3.2 5.1 5.9 4.6 2.8 5.8 5.6 5.0 4.4 5.6 4.0 7.0 4.3 4.7 3.9 5.4
##  [523] 3.8 5.4 4.1 5.4 5.3 4.1 4.9 3.2 4.7 6.9 4.1 5.2 5.8 5.0 4.7 4.1 3.8 4.4
##  [541] 6.4 4.2 4.6 2.0 4.0 3.8 5.9 5.2 4.4 3.7 4.5 4.3 4.3 5.3 4.7 5.8 4.5 3.9
##  [559] 5.9 4.4 5.2 4.7 2.5 4.5 4.1 2.1 4.7 3.7 3.3 5.4 4.4 6.7 5.7 4.1 3.6 6.1
##  [577] 4.3 5.4 3.9 3.9 3.5 3.8 5.2 3.2 6.1 5.3 5.0 7.1 3.5 4.1 5.2 5.3 3.8 5.6
##  [595] 4.0 5.4 4.5 3.8 4.9 6.7 4.4 6.3 6.5 4.0 6.1 3.2 2.5 5.0 3.8 4.4 4.9 3.7
##  [613] 4.7 4.2 3.8 3.7 3.9 3.7 3.9 4.1 2.7 4.4 2.9 4.2 3.8 4.5 4.3 2.7 5.8 4.1
##  [631] 3.3 4.2 2.5 3.6 4.4 5.8 4.3 3.5 5.1 4.3 5.6 4.2 5.7 3.9 4.3 4.5 5.2 3.1
##  [649] 3.0 3.9 5.6 3.5 5.3 5.1 3.9 4.5 6.2 4.7 4.4 4.6 4.0 3.1 5.4 3.7 4.6 3.5
##  [667] 3.6 4.8 3.3 6.0 4.0 4.0 3.9 4.3 5.1 4.8 4.8 3.3 2.7 5.7 3.8 4.0 3.8 4.1
##  [685] 5.7 3.6 5.5 3.3 3.5 3.6 5.0 4.1 4.9 2.9 5.1 2.9 5.7 5.2 5.1 5.1 3.5 4.4
##  [703] 5.8 5.4 5.2 5.7 5.1 4.7 5.1 5.1 5.0 5.4 5.7 3.5 4.1 4.1 6.1 3.7 6.2 5.4
##  [721] 4.7 4.1 4.1 6.0 6.0 5.3 4.3 4.3 2.4 4.3 4.2 5.0 5.0 5.3 2.7 3.8 4.7 4.4
##  [739] 5.0 3.7 4.1 4.7 5.0 4.0 3.9 3.9 4.7 4.0 2.9 3.9 4.8 3.5 2.9 3.1 4.1 6.2
##  [757] 4.0 5.6 4.6 3.9 4.3 3.9 5.1 4.0 3.0 4.7 4.8 4.4 3.9 5.3 5.1 4.6 4.0 4.0
##  [775] 5.2 3.6 4.2 4.8 4.7 4.8 6.0 4.5 3.1 4.2 4.4 4.9 5.4 3.6 5.2 3.3 5.7 5.7
##  [793] 4.8 5.7 4.9 4.4 5.2 4.9 4.3 3.6 5.6 6.0 4.9 3.9 3.3 4.1 4.4 4.4 4.6 3.1
##  [811] 4.7 4.0 4.8 5.3 6.4 4.8 4.4 4.9 3.9 4.0 3.8 5.4 5.8 3.7 3.2 3.7 4.9 4.4
##  [829] 5.1 4.8 4.3 3.6 2.9 5.5 4.3 2.9 3.9 4.7 3.8 3.4 4.6 5.3 5.2 4.3 4.6 6.0
##  [847] 4.6 3.1 5.0 3.8 3.9 4.4 3.9 4.2 5.1 4.9 5.2 5.1 5.1 4.3 5.3 3.2 3.6 3.7
##  [865] 3.9 4.4 3.8 3.6 4.7 5.0 5.3 5.1 5.7 4.3 5.4 5.2 3.9 5.3 3.7 4.3 4.7 4.0
##  [883] 4.3 4.4 5.9 3.8 3.3 4.0 3.5 3.9 3.5 3.3 4.8 4.8 5.1 4.1 4.7 3.6 5.8 4.3
##  [901] 4.7 4.4 4.8 4.9 3.0 4.7 4.6 4.1 4.1 5.5 6.0 5.0 5.4 3.4 5.9 5.1 3.7 3.8
##  [919] 3.2 4.3 3.7 4.5 4.2 4.4 3.0 4.8 3.6 5.8 3.6 3.9 4.7 4.2 4.8 6.9 4.2 4.8
##  [937] 4.6 4.6 5.1 5.1 5.4 5.8 5.5 4.0 2.2 5.3 5.1 4.5 4.2 3.8 3.8 6.2 5.4 4.5
##  [955] 5.2 4.7 3.4 4.4 3.8 4.8 3.9 4.3 3.1 2.6 3.6 3.6 4.4 3.3 4.6 5.5 4.8 4.7
##  [973] 5.0 4.8 3.3 4.5 4.7 5.4 2.7 3.6 3.2 4.7 4.5 5.1 3.7 3.3 3.5 6.2 4.1 4.1
##  [991] 5.2 4.5 3.8 4.5 3.0 3.3 5.7 4.6 3.1 2.0
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
##    [1] 3.8 5.2 4.3 5.1 4.8 5.2 4.2 4.9 5.1 4.3 5.4 3.3 4.3 4.3 3.1 4.2 5.6 2.5
##   [19] 5.4 5.2 3.9 4.0 4.7 4.2 4.6 5.0 2.7 3.2 4.9 3.3 4.1 4.9 5.1 3.9 5.9 3.8
##   [37] 2.9 4.7 4.6 2.9 5.3 3.8 3.7 3.0 3.0 5.4 5.3 3.0 4.4 4.6 5.2 4.6 3.6 3.8
##   [55] 4.7 6.7 4.3 5.7 4.9 5.0 2.6 4.1 4.7 4.4 5.6 4.0 4.6 5.1 4.4 4.9 5.1 4.5
##   [73] 3.4 5.2 4.3 4.6 5.1 3.9 4.4 4.2 3.2 3.2 3.0 5.2 5.8 3.4 5.6 3.3 5.0 4.5
##   [91] 4.1 4.8 2.8 5.2 5.8 3.1 5.3 4.3 6.5 4.8 6.6 4.8 3.6 3.4 5.0 4.5 4.6 6.4
##  [109] 4.9 5.4 4.8 2.6 4.3 4.6 4.3 3.2 4.7 4.3 4.7 5.2 4.5 5.8 5.3 3.2 3.4 5.0
##  [127] 4.8 5.6 3.7 4.7 5.0 3.8 4.2 2.0 3.8 2.5 5.9 5.8 3.3 3.5 4.8 2.7 3.8 6.2
##  [145] 5.2 4.4 6.0 4.1 5.3 5.7 5.0 5.4 4.9 4.5 5.0 3.4 5.0 4.6 5.1 4.3 5.0 3.1
##  [163] 3.4 4.0 4.6 5.5 5.1 6.0 4.8 4.5 5.4 3.9 5.7 3.0 4.9 3.9 4.9 4.5 4.3 5.6
##  [181] 4.9 3.8 5.3 4.0 4.6 6.1 5.9 3.2 5.1 5.0 4.4 4.3 3.9 5.9 5.6 4.2 5.3 4.4
##  [199] 5.6 4.8 4.1 4.6 3.8 3.4 3.9 3.8 3.4 5.5 4.1 4.9 3.2 5.1 4.2 3.8 2.7 4.8
##  [217] 3.9 4.6 4.3 2.7 4.2 4.3 5.0 4.9 4.7 5.8 5.7 4.7 4.2 3.4 5.8 2.8 4.3 5.5
##  [235] 4.3 3.6 4.2 4.4 4.6 4.7 2.4 5.6 4.7 4.6 5.8 5.5 4.5 5.3 4.4 4.0 4.5 4.6
##  [253] 4.9 5.8 3.0 4.4 3.5 3.2 5.2 5.5 5.1 5.0 4.9 4.2 4.3 4.0 5.2 4.1 5.0 3.8
##  [271] 5.5 3.5 2.9 4.2 5.9 4.7 5.5 5.2 3.9 4.5 5.2 4.1 5.9 3.8 4.0 5.3 3.3 3.7
##  [289] 3.9 5.4 4.7 5.1 6.6 5.7 3.4 5.1 4.7 3.3 4.8 6.0 3.9 4.6 4.6 4.5 4.0 3.5
##  [307] 4.4 4.3 5.1 4.9 3.8 4.5 4.7 5.8 4.6 2.8 3.9 4.6 4.0 3.9 6.0 4.1 4.2 3.1
##  [325] 4.5 5.5 5.2 4.5 5.3 3.7 5.2 3.7 5.0 3.7 5.0 5.1 2.6 4.6 3.7 6.4 3.5 4.3
##  [343] 3.7 5.1 6.2 5.7 3.5 2.4 3.7 3.9 3.6 3.2 4.8 5.2 2.4 3.5 4.8 5.2 4.1 5.5
##  [361] 4.7 4.3 5.5 3.5 6.3 4.7 4.3 5.0 4.9 4.1 4.6 3.7 4.7 4.2 3.2 3.1 4.0 4.5
##  [379] 5.6 3.5 3.9 3.2 5.0 4.4 4.5 5.6 3.9 5.1 5.0 3.7 5.8 4.2 3.1 5.0 4.6 4.9
##  [397] 6.0 3.1 4.5 4.8 4.8 3.3 4.8 3.3 4.3 3.9 4.6 3.5 3.9 5.7 5.4 3.6 3.9 4.8
##  [415] 3.9 5.0 5.1 3.5 5.5 5.4 5.9 5.4 5.7 2.8 4.9 4.6 4.4 4.5 5.1 3.9 4.0 4.6
##  [433] 3.2 2.5 4.3 3.4 5.4 4.4 5.0 3.0 3.0 5.0 5.1 4.7 3.7 4.1 3.6 3.9 4.6 4.7
##  [451] 4.2 4.8 4.3 6.0 5.5 5.1 3.8 3.7 4.3 5.2 4.1 4.5 5.3 6.1 5.1 3.6 3.7 3.9
##  [469] 3.4 5.9 3.2 3.8 5.8 5.2 6.0 5.6 3.6 5.9 3.8 3.5 4.9 3.7 3.9 5.4 4.2 4.5
##  [487] 4.5 4.8 6.4 5.7 4.7 5.2 3.9 4.1 3.6 3.9 5.5 3.8 5.7 3.5 4.1 4.4 4.2 5.3
##  [505] 5.0 4.3 5.0 5.2 4.7 3.9 6.6 6.0 3.9 4.9 3.3 4.7 3.0 3.3 3.5 3.4 4.2 3.6
##  [523] 2.9 3.9 4.4 3.0 4.2 4.3 3.5 3.8 5.7 2.8 5.0 5.3 3.8 5.4 3.7 3.4 4.3 5.8
##  [541] 4.5 4.4 4.3 4.5 4.1 4.7 4.0 3.8 4.5 5.0 5.5 4.4 4.4 3.9 4.8 3.9 3.3 3.3
##  [559] 6.2 4.3 4.8 3.2 6.3 4.2 4.4 3.1 4.6 4.5 3.1 4.4 4.3 5.4 5.0 3.4 2.8 4.4
##  [577] 4.5 6.0 4.6 4.8 4.3 4.7 3.8 4.2 4.1 3.9 4.6 3.7 4.9 5.2 5.4 4.7 4.8 4.2
##  [595] 4.7 3.1 4.0 4.0 5.2 6.2 3.8 4.5 5.6 4.1 4.5 2.9 5.3 4.2 3.9 5.3 5.8 3.9
##  [613] 4.9 3.9 3.8 3.4 4.0 3.4 2.5 4.1 4.9 4.8 4.7 5.3 5.9 4.9 3.6 4.0 6.0 4.3
##  [631] 4.5 4.3 4.5 4.9 3.3 4.3 4.9 5.1 5.2 6.2 3.8 3.5 5.4 5.6 5.5 1.7 3.7 4.7
##  [649] 4.6 4.5 4.3 4.5 3.4 5.2 3.9 6.3 5.0 4.4 6.7 3.3 4.2 3.9 3.5 4.0 4.4 5.2
##  [667] 3.0 4.5 4.3 5.6 3.8 6.3 6.8 4.7 4.8 4.0 4.9 4.8 4.2 4.9 5.1 6.0 4.1 4.6
##  [685] 2.8 3.6 3.4 5.1 4.6 4.3 4.6 3.8 4.4 3.8 4.8 2.3 5.3 4.7 4.6 5.1 5.0 5.0
##  [703] 4.7 6.8 4.9 3.7 4.7 4.6 4.5 5.2 5.4 3.0 4.2 5.6 4.3 2.8 4.5 3.5 3.6 3.7
##  [721] 4.9 4.9 3.3 4.0 3.8 5.9 4.8 5.0 3.9 4.1 4.3 4.5 4.1 5.2 4.7 4.6 5.9 4.6
##  [739] 3.2 3.9 4.0 3.8 4.5 4.6 5.7 3.6 4.6 4.7 5.3 4.2 3.5 3.2 4.0 3.0 4.7 5.5
##  [757] 3.4 4.5 4.3 3.7 5.5 4.8 4.1 5.0 5.7 4.0 3.7 2.3 4.3 3.8 4.3 5.2 3.8 4.0
##  [775] 4.5 4.3 4.2 4.5 2.4 5.4 4.2 4.4 5.2 3.7 4.9 5.0 5.3 5.4 5.0 4.7 5.6 5.2
##  [793] 4.6 2.9 4.9 3.8 4.1 3.7 3.6 5.0 5.2 5.0 3.8 2.4 4.3 4.3 3.8 4.4 4.5 2.7
##  [811] 4.8 3.9 5.2 5.3 3.8 4.0 3.8 5.8 4.9 4.4 5.0 4.5 4.7 5.5 4.4 4.3 5.9 3.1
##  [829] 4.7 5.3 4.7 2.9 5.6 5.1 4.5 5.6 3.5 5.3 4.1 3.5 3.2 3.2 4.2 5.3 3.7 4.8
##  [847] 4.8 5.9 3.7 5.5 2.8 4.3 4.4 5.3 4.7 5.0 5.0 4.2 3.3 4.4 6.2 4.3 3.9 3.4
##  [865] 5.2 4.6 4.9 4.1 5.5 3.3 3.6 3.9 4.3 3.8 4.7 5.3 4.2 3.7 3.7 5.5 5.0 4.3
##  [883] 3.5 5.4 4.8 4.3 4.8 3.7 3.7 3.9 5.6 3.7 4.6 5.1 4.4 3.8 4.9 4.3 5.2 4.1
##  [901] 5.5 4.4 3.6 2.8 1.9 3.8 4.4 3.4 4.7 3.2 4.2 4.9 4.5 3.2 4.5 4.3 5.0 5.9
##  [919] 4.4 4.7 5.1 3.6 2.0 5.0 3.5 4.5 5.3 3.6 4.7 4.8 4.8 4.8 5.7 4.1 5.1 3.3
##  [937] 4.8 5.9 3.1 4.6 4.2 3.9 4.5 4.8 3.9 5.1 4.6 3.7 4.1 5.3 5.8 5.1 3.6 3.5
##  [955] 6.7 4.9 5.3 6.3 5.3 5.5 5.0 5.0 2.9 4.6 4.9 4.1 3.7 4.6 5.8 4.3 4.9 4.4
##  [973] 4.5 3.8 2.9 5.4 3.4 5.1 6.2 4.0 5.0 4.4 4.9 4.5 4.4 4.4 4.0 5.0 4.2 2.5
##  [991] 7.1 5.3 4.9 4.1 5.7 3.1 4.3 5.6 5.3 4.7
## 
## $func.thetastar
## [1] -0.0391
## 
## $jack.boot.val
##  [1]  0.48205882  0.36666667  0.26160991  0.15335366  0.03198847 -0.11540616
##  [7] -0.19108635 -0.36886228 -0.38051576 -0.49411765
## 
## $jack.boot.se
## [1] 0.9619548
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
##    [1] 4.2 5.1 5.4 4.5 4.3 4.0 6.0 5.2 3.7 5.1 2.3 5.3 3.2 5.6 3.8 4.3 5.5 4.1
##   [19] 4.6 3.2 3.0 3.6 4.3 5.5 3.5 5.7 4.8 6.1 6.3 3.3 3.9 3.4 3.5 4.2 5.3 4.5
##   [37] 3.3 6.0 4.5 4.5 5.7 5.0 5.2 4.6 4.4 5.3 4.7 5.9 6.2 4.5 3.7 2.9 5.5 4.8
##   [55] 5.0 3.6 5.1 5.2 3.5 3.0 6.2 6.0 4.8 3.6 4.1 4.1 3.7 6.8 3.7 4.6 4.4 4.8
##   [73] 3.1 4.2 5.6 5.0 3.6 3.1 3.7 5.3 3.7 4.0 4.0 3.3 6.3 5.6 6.0 4.4 4.4 5.1
##   [91] 4.1 5.4 5.8 4.7 5.1 3.5 5.1 4.3 4.3 4.4 5.4 5.5 2.5 3.5 4.2 4.6 5.0 2.9
##  [109] 4.9 4.0 4.7 4.3 4.4 4.6 4.3 5.3 4.6 2.9 5.3 4.1 3.3 5.2 5.1 5.2 5.8 4.7
##  [127] 4.0 6.0 5.9 6.5 3.4 5.0 4.0 5.5 4.2 4.0 4.4 4.2 4.5 5.0 5.1 3.5 5.1 4.7
##  [145] 5.2 3.9 4.5 3.7 4.4 3.8 3.3 3.8 5.6 3.6 5.1 4.0 5.8 4.0 5.7 3.6 3.9 6.2
##  [163] 3.6 5.6 4.6 5.3 4.2 5.1 4.0 4.8 4.6 4.9 6.1 4.3 3.6 4.1 4.1 5.3 4.9 5.0
##  [181] 4.8 5.7 2.7 4.7 4.8 3.8 5.2 5.7 4.9 4.2 5.2 3.4 5.0 3.3 3.4 4.1 4.7 6.1
##  [199] 2.3 4.6 3.5 5.4 4.8 3.4 4.1 4.7 2.8 4.8 3.9 3.0 4.6 5.7 6.4 4.1 3.4 3.5
##  [217] 4.0 5.0 4.2 5.3 4.5 4.3 3.3 3.7 2.3 5.4 4.5 3.9 5.0 2.8 4.8 5.2 5.2 4.9
##  [235] 5.0 4.4 3.1 5.4 5.7 4.7 4.9 5.0 3.0 3.7 3.4 4.6 5.2 4.8 3.5 4.0 4.1 4.1
##  [253] 5.8 4.2 4.7 3.8 3.5 4.5 5.3 4.8 2.9 2.8 4.1 5.2 4.1 3.3 5.2 3.5 5.0 4.3
##  [271] 2.9 6.4 5.1 3.6 3.0 6.0 3.2 4.5 3.5 5.0 4.4 5.3 5.4 5.2 4.8 4.3 3.5 6.2
##  [289] 4.5 5.4 4.7 4.4 4.1 4.6 3.3 5.6 4.4 4.3 5.0 3.4 3.9 3.1 4.1 3.9 5.2 4.9
##  [307] 3.7 2.8 4.3 3.7 4.7 4.8 4.5 5.3 4.1 4.9 2.9 4.5 4.8 6.1 6.0 4.2 4.1 5.0
##  [325] 5.1 6.7 5.1 6.0 5.1 5.3 4.3 4.5 5.5 3.9 3.8 4.8 5.2 5.3 5.5 3.5 3.6 5.4
##  [343] 4.0 2.7 4.2 4.1 5.4 4.1 4.3 4.5 4.4 4.4 3.6 5.1 3.7 4.1 5.6 3.3 3.4 5.3
##  [361] 4.7 4.6 5.4 4.2 3.5 5.4 6.4 4.8 2.7 4.1 4.2 5.9 5.5 1.3 5.3 4.9 4.3 5.7
##  [379] 5.1 5.5 3.9 3.7 2.7 4.0 4.4 4.8 4.0 4.5 4.1 4.8 4.0 4.9 6.2 4.0 5.3 3.3
##  [397] 4.4 3.4 3.9 4.9 3.9 5.1 4.9 5.3 3.7 5.6 4.3 2.9 5.5 5.1 4.4 5.4 2.8 5.8
##  [415] 4.2 4.1 3.2 3.7 4.5 4.7 5.1 1.9 5.0 5.2 5.1 5.2 4.2 4.8 4.2 5.3 5.2 3.5
##  [433] 6.2 4.3 3.6 5.2 3.0 4.6 4.6 6.1 3.8 3.2 4.9 3.2 3.6 4.3 4.6 6.4 5.0 6.2
##  [451] 4.5 5.3 4.1 6.5 5.6 3.9 5.6 5.0 5.0 5.1 5.9 3.8 5.8 4.0 3.5 5.6 6.4 4.9
##  [469] 5.3 3.8 4.5 5.6 5.2 4.1 3.7 5.0 3.3 5.1 4.5 4.2 3.3 4.0 5.1 3.8 4.5 5.8
##  [487] 4.5 5.1 4.8 3.9 3.7 3.5 3.7 3.7 4.5 4.3 3.4 4.6 2.5 5.4 5.9 4.9 3.4 2.9
##  [505] 4.9 4.6 4.7 4.1 4.3 4.7 5.2 4.1 5.0 4.4 4.9 4.6 2.2 3.2 4.1 5.8 4.6 4.1
##  [523] 5.6 3.8 5.4 4.2 3.7 4.6 5.6 4.9 5.6 4.8 5.5 5.8 4.0 4.6 4.3 5.3 5.7 5.0
##  [541] 4.1 4.8 3.6 4.2 3.8 4.6 4.5 3.9 3.9 6.7 5.0 4.5 6.4 5.5 3.8 4.4 5.3 2.9
##  [559] 5.0 3.2 4.6 3.3 4.9 4.3 4.8 4.5 4.7 5.2 4.3 3.5 5.2 3.9 3.8 4.8 3.8 3.9
##  [577] 3.3 2.3 3.8 5.1 5.5 4.2 5.7 4.2 3.7 3.9 5.1 5.9 4.7 3.6 3.7 5.3 4.6 3.6
##  [595] 4.9 5.2 3.5 4.7 4.4 5.1 2.6 4.0 5.5 6.1 6.3 4.3 5.0 6.2 5.2 5.0 5.8 4.9
##  [613] 4.5 4.4 4.9 5.9 4.5 4.2 4.4 5.5 3.4 4.4 6.3 5.1 6.4 4.5 5.3 5.2 4.5 4.3
##  [631] 4.0 6.2 4.8 3.2 5.9 5.6 5.0 4.6 6.0 2.6 4.1 4.8 4.0 5.7 3.3 3.6 4.8 4.0
##  [649] 2.5 5.5 4.7 4.3 4.5 4.8 5.4 3.7 5.3 4.3 3.6 5.3 4.8 4.5 5.7 6.6 5.6 4.9
##  [667] 4.0 6.4 4.8 3.8 3.4 4.8 4.2 3.2 4.8 3.7 4.9 5.4 5.4 5.6 2.0 3.9 4.5 4.8
##  [685] 3.9 6.4 5.0 4.5 3.3 5.1 3.5 3.5 4.6 4.5 3.5 4.9 5.6 4.5 5.5 4.4 5.3 5.3
##  [703] 4.4 5.0 3.1 3.9 2.2 4.4 4.5 4.5 4.8 3.6 3.4 5.7 5.6 6.4 5.0 4.3 2.6 6.5
##  [721] 5.8 3.2 5.4 4.4 4.3 5.4 4.8 4.1 3.9 3.5 4.2 4.2 3.8 5.1 4.5 5.3 5.0 4.7
##  [739] 4.6 4.2 4.7 3.8 3.8 4.6 3.8 5.1 5.9 5.0 4.5 4.5 5.8 5.4 4.7 5.1 3.0 7.4
##  [757] 3.5 4.5 4.5 5.6 3.2 4.3 5.8 5.3 5.0 3.6 5.7 4.5 4.5 4.3 4.1 4.3 4.2 3.9
##  [775] 3.9 3.8 3.7 5.5 3.9 4.4 5.6 5.3 5.1 4.9 3.5 3.0 5.2 3.9 2.5 3.4 4.3 4.3
##  [793] 3.7 3.4 3.6 4.4 3.3 3.9 2.4 3.9 5.9 5.3 3.4 2.6 5.3 4.0 3.6 3.9 5.2 4.1
##  [811] 4.8 4.9 3.9 6.0 4.3 3.9 4.1 4.5 3.6 6.0 3.2 3.4 4.3 2.9 5.5 4.8 4.4 3.1
##  [829] 5.2 5.7 5.0 4.4 5.5 3.4 4.4 3.3 3.8 3.6 4.1 6.1 4.4 5.5 4.9 3.9 4.7 3.9
##  [847] 5.0 4.0 4.6 4.7 5.1 5.5 4.7 4.5 3.7 5.5 4.2 5.4 4.0 4.3 5.1 4.6 4.7 4.6
##  [865] 4.7 5.2 4.1 4.0 4.9 4.1 5.1 5.1 4.7 4.5 4.0 4.7 3.8 5.2 5.3 3.4 4.9 2.7
##  [883] 4.4 4.5 1.4 4.1 4.6 2.6 3.7 3.5 5.1 3.6 3.5 3.8 3.0 5.1 4.3 4.2 5.0 4.9
##  [901] 4.6 6.7 5.3 4.2 4.0 4.2 3.3 5.7 3.5 5.1 6.1 3.4 4.2 4.7 3.2 5.2 5.1 5.1
##  [919] 4.1 5.1 6.0 4.7 4.1 4.4 4.9 6.0 5.9 5.0 2.2 5.2 4.9 4.4 2.8 4.9 5.0 4.3
##  [937] 5.3 1.6 5.1 4.5 3.6 6.0 4.5 4.7 4.4 4.0 6.3 5.8 4.3 4.3 3.2 3.9 5.0 4.7
##  [955] 5.5 5.0 3.3 6.2 5.1 5.3 4.9 3.1 4.1 4.4 4.1 4.3 5.4 5.1 5.9 5.9 4.7 2.6
##  [973] 5.2 3.2 3.7 3.2 2.4 3.5 4.1 3.7 3.9 3.8 5.2 3.8 4.3 3.8 5.1 3.0 4.5 4.8
##  [991] 6.6 3.4 3.5 5.0 3.7 3.3 4.1 5.6 5.1 4.8
## 
## $func.thetastar
## 72% 
## 5.1 
## 
## $jack.boot.val
##  [1] 5.500 5.352 5.300 5.300 5.000 5.000 5.000 4.800 4.600 4.500
## 
## $jack.boot.se
## [1] 0.9446493
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
## [1] 1.56168
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
##   3.593683   6.675358 
##  (1.538383) (3.066881)
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
## [1]  0.8463876  1.0424363 -0.4814369  0.7722495  0.6864303  0.3470842
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
##    [1] -0.028764481  1.344513262  1.142158279  0.950184010  0.932854687
##    [6]  2.132682626  1.394799866  1.476346208  1.721812214  1.480007462
##   [11]  1.912012282  0.064955339 -0.512150003  0.696909268  0.712514812
##   [16]  0.737975315  2.126432393  1.191052603  0.779384984  1.902087765
##   [21] -0.688061767  0.795299951  1.258848728  0.518134219  1.718836420
##   [26]  1.353862396  2.315013019  0.854748124  0.901522571  1.284298061
##   [31]  1.341724501  1.579480258  0.790613014  1.368921762  1.488024763
##   [36]  0.966178551  1.092060869  0.937501929 -0.155204060  0.131322288
##   [41]  0.783801396  1.646016405 -0.182261341  1.283146308  0.864329164
##   [46]  1.162442946  1.494384308  1.721714994  0.606173340  1.710608792
##   [51] -0.580167862  1.546595158  1.559787207  1.557117857 -0.099795075
##   [56]  0.579773685  1.558067670  2.168346642 -0.059287634  2.109305233
##   [61]  1.941264202  1.727732174  0.278111281  0.563764212  0.974959028
##   [66]  0.994639014 -0.334281346 -0.680651850  0.416173771 -0.382152047
##   [71]  0.328860105  1.072661303  0.126407743 -0.485528932 -0.045802464
##   [76] -0.047307601 -0.512508457  0.933385385  2.296016494  1.239874983
##   [81]  2.285870407  0.956248976  1.119771208  0.757621096  1.534199497
##   [86]  1.732555044 -0.058878699  2.379043716  1.774538463  1.033061813
##   [91]  1.033838111  1.963392690  1.959820388  1.249816832  1.805913103
##   [96]  2.213751652  0.627376287  1.570471619  1.410177328  0.678300529
##  [101]  1.559480366  1.947305555  0.485287937  1.192369164  0.059759527
##  [106]  1.327682432  0.600861596  1.553090932  2.267588112  1.758522073
##  [111]  2.031540747  0.526347020  0.210098405  1.597158682  2.344595148
##  [116]  0.839095755  1.815367076  0.098529590  0.077250679  1.386667069
##  [121]  0.640636736  1.147226758 -0.183458460  1.082919190  0.856884051
##  [126] -0.120508511  0.696270162  1.242445844  2.033445781  1.804135958
##  [131] -0.178466006  0.414502914  1.207441890  1.997912818  0.704349577
##  [136]  1.208346168  1.065155326  0.972857456  1.684917834  2.509246621
##  [141]  0.558408178  1.747595088  2.339474251  0.147503737  0.658510211
##  [146]  1.872097843  1.855955324  0.690832145 -0.574807622  0.003568954
##  [151] -0.012194826  0.792641332 -1.241869969 -0.224112227  0.649572649
##  [156]  1.343434948  1.716847894  0.296806281  1.109469390  1.424207357
##  [161]  1.270734668  1.046294950  1.536289250  1.236013573  1.688636090
##  [166] -0.243255686  1.032037103  1.169953285 -0.306770946  1.368196772
##  [171]  0.059714122  1.159484324 -0.499870017  1.785639354  1.417551468
##  [176] -0.535141665  1.520731638  0.155945584  0.710899534  0.950168349
##  [181]  1.815192501  1.717032430  1.129219053 -0.149782986  1.382751370
##  [186] -0.366729808  0.315455763  2.434955238  0.692249263  1.152005631
##  [191]  2.237996716  0.962670743  1.243180759  0.990131051  1.081278804
##  [196]  2.065043363  0.656832225  1.174143306  0.403154284  0.583624371
##  [201]  2.088538433  1.176970897  2.301193674 -0.319207231  0.642151939
##  [206]  1.711521197  1.923286443  1.401273457  1.686012768  1.012564608
##  [211]  1.572396202  1.163650960 -0.095763632  1.583753056  0.602974950
##  [216]  0.505172112  0.404877660  1.057476987  0.742155895  1.656101786
##  [221] -0.166768103  0.828306255  0.136012290  2.293907198  1.036775164
##  [226]  1.462523414  1.044176266  0.252180775  0.642781152  2.035503221
##  [231]  1.259606670  1.916790204  1.036761324  1.499354220  1.792767237
##  [236]  0.779188949  0.273792090  1.320071187  0.643860002  1.686793471
##  [241]  0.360921017  1.713928129  1.363486895  0.300360513  1.041382782
##  [246]  1.514398312 -0.203094659  1.204977705 -0.051886203 -2.023185618
##  [251]  1.724552271 -0.104464183  0.522509769  1.758365602  1.117614262
##  [256]  1.923686017  2.123098875  0.593336975  0.195240416  1.322891220
##  [261]  1.416612135 -0.223792479  2.301645628  2.157151070  1.902779746
##  [266]  0.273743927  0.552030031  0.582247946 -1.854395954  2.262828181
##  [271]  2.321003566  1.703309045  0.907807944 -0.443843910  0.582479337
##  [276]  0.962873475  0.197064685  1.757545141  1.528590046  1.727313339
##  [281] -0.122127011  1.624875765 -0.751474690  0.961888995  0.639113268
##  [286]  0.093754671  1.615301412  1.655502779  0.982932287  1.946932800
##  [291]  1.495695980 -0.279406957  1.689803381  1.554840754  1.810424488
##  [296]  0.583376255 -0.100185887  2.026386476  1.182459045  1.327642634
##  [301]  1.945226271 -0.127700870  2.056624519  1.128763905  1.098120255
##  [306]  1.208497290  0.220470274  1.567106238 -0.153871218  0.943827711
##  [311]  1.268809187  0.289905987  0.163271539  1.967542438 -1.935321617
##  [316]  1.321067447  1.363886177  1.776949504  1.271544247  1.033780570
##  [321]  0.410582612  0.953819377  0.969186247  2.106216161  0.443253697
##  [326]  0.582145648  0.892573259  1.031467472  1.812423829  1.257880244
##  [331]  1.719579719  1.638689921  2.186304642 -0.352118581  1.414838154
##  [336]  1.165478754  0.666825466  2.116590583  1.095879628  1.301171760
##  [341]  1.662106086  0.281537240  1.266952075  1.734109395  1.301555703
##  [346] -0.194305247 -0.511212546  2.322552548 -0.047307601  0.246516941
##  [351]  1.207092866 -0.197459943  0.676678222  0.525867634  1.250074862
##  [356]  0.417215965  0.181870518  1.210674434 -0.327940249  0.891175943
##  [361]  2.158055043  0.306446415  1.336773107 -0.424696064  1.868059224
##  [366]  0.729520815  1.156845029  1.028224189  0.718961707  0.344688514
##  [371]  0.632254799  0.653536208  1.789999219  2.006059543  0.390744856
##  [376]  1.743289270  1.368177979  0.358485270  1.542475054 -0.099518256
##  [381]  1.407285086  0.682325172  1.052988657  1.685235951  1.010578078
##  [386]  1.362996045  2.442595531 -0.627070976  1.541362511  0.101701338
##  [391]  1.854845466  0.922323083  1.578495359  1.056454720  1.076266567
##  [396]  0.916547823  1.170685285  1.163346863  1.728865121  1.006781233
##  [401]  0.043155115  0.849736991  0.563718562  1.430424177  0.115453574
##  [406]  1.071880388  0.001157264 -0.067540604  1.439215083  2.149883259
##  [411]  1.955259062  1.224634604  1.384413518  1.109000347  0.602311590
##  [416]  0.402848037  1.600817686  1.907917009  1.199252732 -0.114179888
##  [421]  0.640475166  2.034184673  1.192567826 -0.609498768  1.933041478
##  [426]  0.043243565  1.142638822  1.258972765 -0.109508407  1.255565605
##  [431]  0.914705217 -0.057814909 -0.001215050 -0.472964091  0.222350673
##  [436]  1.142230986  0.902155875  1.361041798  1.569665433 -0.219329765
##  [441]  0.232334126 -0.346938719  0.021465771  1.190412027  1.706990705
##  [446] -0.327531566  0.562256588 -1.217721866  0.792562251 -0.355389774
##  [451] -0.171248733  1.802321020 -0.192239885 -0.616752636 -0.352118581
##  [456]  1.055268703 -0.549755959  0.569839721  1.225148136  0.239240835
##  [461]  0.143453192  1.921035136  0.478765204  1.814233834  1.287925491
##  [466]  0.691150372  0.565029427  0.338297064  1.911073693  1.036634212
##  [471]  1.438191378  0.233993250  1.383459927  0.505633053  1.214621303
##  [476]  1.280885684  0.830508692  1.937641588  0.953231334  1.072242243
##  [481]  2.281647827  0.258020248  0.740351543 -0.587747288 -0.255002514
##  [486] -0.008540066  1.322040951  0.358392339  1.569071108  2.591122533
##  [491] -0.214127376  1.110227705 -0.194361825 -0.050033130  1.035935333
##  [496]  1.743145683  0.117677401  0.880190653  1.787687081  1.492409432
##  [501]  0.625652743  1.824304389  1.085355881  1.791684013  0.703754856
##  [506] -0.320193511 -0.400567415  1.727163085  0.520941436 -0.534099554
##  [511]  0.143888266  1.505927787  1.962183336 -0.637737401  0.259748676
##  [516]  1.722726571 -0.283984490  0.274795973  1.511859736 -0.529104556
##  [521]  1.030160356  1.709373389  0.252743325  2.612665339  1.836590550
##  [526]  0.762816116  1.208430915 -0.249207967  0.932026197  0.153563589
##  [531]  2.128037758  2.204920462  1.366673026  1.391417815  1.892793238
##  [536]  1.724313982  1.719078878  0.962670743  0.623804495 -0.421625663
##  [541]  0.504946454  1.565861676  0.247930803  1.559476140  0.260831015
##  [546]  1.260225778  0.268122429  1.956765909  1.370698755 -0.607309036
##  [551]  1.060487481  1.373322605  1.850331087 -0.599745972  1.317060970
##  [556]  1.700263404  2.080301590 -0.539225470  1.041509838 -0.576454681
##  [561]  0.570658239  1.362750121  1.803004044  1.806747204 -0.244798741
##  [566]  1.517194633  1.641246023  2.389270594  1.463004241  1.731418218
##  [571]  0.488075496  0.711936439  2.266697640  1.162726258 -0.265359909
##  [576]  1.382390138  1.366673026  1.943424530  0.091443192  1.444016362
##  [581] -0.433409297 -1.016602275  1.868674268  1.197050752  0.728777573
##  [586]  2.091495277  0.671794529  0.588143116  0.241053472  0.856829455
##  [591]  1.595298541  0.122560181  1.101122036  1.506633051  0.014566546
##  [596]  1.202623222 -0.117935764 -0.468193448  0.682978953  1.244945204
##  [601]  0.939340985 -0.090959608 -2.130793568 -0.529540373  1.954340658
##  [606]  0.307565976  1.771339167  1.973286386 -0.063967472  0.156769862
##  [611]  0.946530951  1.062610952  0.959127516  1.884237512  1.166806236
##  [616]  1.120880191 -0.098300206  1.293267315  1.650755244  1.500736747
##  [621]  1.449253243 -0.580478708  1.764165164  0.942256144  1.155379827
##  [626]  1.189552971 -0.395683672 -0.673366218  1.595237092  1.225126301
##  [631]  2.582352188  0.328467401  2.154817966  1.509912345  1.581262501
##  [636] -0.520624872  1.350531107  0.904401900  0.082928263  1.545757672
##  [641]  1.685235951  0.513455916  0.303603652  0.790325103  1.356199540
##  [646]  0.482014542 -0.576425013  1.028496418 -0.081552170  1.409521731
##  [651]  1.732270484  1.052796795  0.591840902  1.010890985  1.387924300
##  [656]  1.932604562  1.355783094  0.340603671  1.258390238  0.207965784
##  [661]  1.409274478  0.490669073  1.933365658  1.521991309  1.242319275
##  [666] -0.035884524  1.384728288  1.190052796 -0.783504547  1.570644392
##  [671]  1.425436777  1.259095373  2.149601944  2.039305242  0.074372344
##  [676]  1.797047507  0.690832145  1.558992575  0.547059597 -0.121181852
##  [681] -0.582545515  1.387720147  1.591443201  0.744059516 -0.124645146
##  [686]  2.267589448  1.303653485 -0.317877521  0.176997718  1.544790259
##  [691]  2.185443514 -0.578852836  1.489528693  0.570029092  1.448458666
##  [696]  1.082330843  1.943353504 -0.107838290  1.054542682  1.012336799
##  [701]  0.612046311 -2.562833406 -0.542943999 -0.892887910  1.069990049
##  [706] -0.155504350  1.195612671  2.449131867 -0.561281231  1.879806614
##  [711]  0.162554259  0.543220706  0.204553368  0.997042217  1.247993497
##  [716]  0.959987713  1.435586666  1.898098516  1.127467774  1.145118243
##  [721]  0.534727495  1.724552271  1.867820228 -0.029464514  1.893404082
##  [726] -0.010679651  1.681204961  1.232583852  0.920718392  1.714721714
##  [731] -0.061198870  0.629475461  2.136969644  1.238969358  0.271752213
##  [736]  2.373819421  0.294301228  0.825652892  0.225621476 -0.111877520
##  [741] -1.042236273  1.760361601 -0.278479922  0.970693040  1.830374889
##  [746]  1.422564888  1.168939072  0.673737662  2.205624646  1.345605303
##  [751]  1.034161047  1.975055862  1.251887522  1.756719579 -0.066776732
##  [756]  1.025858503  1.354963317  0.842149999  2.388598277  0.753840204
##  [761] -0.172128107  0.571493799  1.615897223 -0.039941296  0.299359682
##  [766] -0.439763686  1.728966994 -0.970527833 -0.352064387  1.431451770
##  [771]  0.653513313  2.062662033 -0.543389977  1.163116034  1.387821118
##  [776]  0.410731023  1.448043713  1.139230127 -0.365036565  1.272881546
##  [781]  1.601364337  2.068180728  0.384806530  1.864949893  1.933159590
##  [786] -0.242426170  1.557425736  1.125303067  0.614773990  1.214193105
##  [791]  1.808094336  1.294857953  0.431368442  0.169885699  0.630230556
##  [796]  1.832236390  2.010839374  1.846252846  2.306956734  0.928128790
##  [801]  0.698773441  0.486309271  0.234577404  0.051167954  2.376622674
##  [806]  2.322488278  1.177987312  0.797869023  1.884017197  1.647064376
##  [811]  0.590803814  1.321932732  0.847595728  0.270770590  0.634183110
##  [816]  2.242638293  0.780050502  0.462989770  1.190052796  0.623804495
##  [821]  1.246207651 -0.567481586 -0.301611000  1.650936487  0.003246364
##  [826]  1.742739635  0.584997153  0.228879187  1.030301401  0.087920038
##  [831]  1.215827456  1.029041762  1.182570489 -0.587277067  2.099415794
##  [836]  1.080757293  0.957885055  0.123165195  1.291109401  1.342649391
##  [841] -0.178413024  0.571652424  1.053679551  2.190277906  1.095623329
##  [846]  1.544705519 -0.660979059  0.643085873 -0.164945235  1.276303422
##  [851]  1.750805021  1.122407918  1.448979676 -0.244662539  1.212313694
##  [856]  1.854280013  0.078482042  2.230883232  1.694384507  1.592257849
##  [861]  0.983808902  1.056845427  1.343521991  1.510811660  1.280741932
##  [866]  0.757192458  1.370230052 -0.804415048  0.684174571  1.697111774
##  [871]  1.000550670 -0.221251104  1.375185239  1.391941202  1.031561557
##  [876]  0.299663049 -0.940845010  1.628086343  1.897444624  1.946054839
##  [881]  1.414263258  1.616708166  0.385459658 -0.472797115  2.279703748
##  [886]  0.905824045  2.273653191  0.359967629  0.435172636  1.918540000
##  [891]  1.124627200  1.892957252  1.565741036  1.453876716  0.343482663
##  [896]  0.279574918  0.957403830  0.444284040  0.927881532  1.077476959
##  [901]  1.593756019  1.858416059  2.184888512  1.744690253  0.817724046
##  [906]  1.529466672  1.897941990  1.367819422  1.596138110  0.527381747
##  [911]  0.345674066  0.140838220  2.073990434  1.689503133  1.362423274
##  [916]  2.215714449  1.531947681  1.709377815 -0.015439161  0.844057361
##  [921]  1.044710161  1.953540348  1.429793320  0.947001251  1.524406470
##  [926]  0.532412757  1.753906105  2.073211739  1.508145574  1.322497926
##  [931] -0.278942931  1.956765909  1.940162488  0.258530382  1.383410613
##  [936]  1.874479726 -0.254989769 -0.167533852 -0.421584179  1.929562387
##  [941]  1.222139943  1.953925139  0.194346571 -0.410916396 -1.308448482
##  [946]  0.162342112  0.939398365  1.965539448 -0.415259108  1.195081424
##  [951] -0.531500924  0.678945689  1.698879088  1.287937406  1.429437684
##  [956]  1.164340689  1.043824696  1.553555473  1.508236937  0.004055616
##  [961]  1.605794121  1.935496928  2.119607370 -0.157579243  1.891864759
##  [966]  1.974312766  0.549951100 -0.497785754  0.969818837  0.246641883
##  [971]  1.418148874  1.362882968  1.610168846  2.071119756 -0.047313163
##  [976]  0.440699022  1.301234886  1.180066893  1.532570713  0.166110491
##  [981]  2.196114726  1.735573467  0.621752915  1.376621586  1.699155223
##  [986]  1.385454820  0.341419852  1.031467472  0.524463907  0.314650765
##  [991] -0.284606762 -0.621047152  0.559285212  1.990156711 -0.497747362
##  [996] -1.147641858  1.719204780  1.165993844  0.391415556  0.934565237
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
##   0.53834964   0.31110825 
##  (0.09838107) (0.06956204)
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
## [1]  0.09497126  0.20319431  1.02034598  0.86188373 -0.60397928  0.06414232
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
## [1] 0.0047
```

```r
sd(results2$thetastar)  #the standard deviation of the theta stars is the SE of the statistic (in this case, the mean)
```

```
## [1] 0.8927289
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
## t1*      4.5 -0.00960961   0.9204467
```

One of the main advantages to the 'boot' package over the 'bootstrap' package is the nicer formatting of the output.

Going back to our original code, lets see how we could reproduce all of these numbers:


```r
table(sample(x,size=length(x),replace=T))
```

```
## 
## 0 1 2 5 6 7 
## 1 2 2 1 1 3
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
## [1] 0.0146
```

```r
se.boot
```

```
## [1] 0.891943
```

Why do our numbers not agree exactly with those of the boot package? This is because our estimates of bias and standard error are just estimates, and they carry with them their own uncertainties. That is one of the reasons we might bother doing jackknife-after-bootstrap.

The 'boot' package has a LOT of functionality. If we have time, we will come back to some of these more complex functions later in the semester as we cover topics like regression and glm.

