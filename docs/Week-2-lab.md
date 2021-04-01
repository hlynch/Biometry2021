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
## 1 2 3 5 6 7 8 9 
## 1 3 1 1 1 1 1 1
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
## [1] 0.0112
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
## [1] 2.714986
```

```r
UL.boot
```

```
## [1] 6.307414
```

Method #2: Simply take the quantiles of the bootstrap statistics


```r
quantile(xmeans,c(0.025,0.975))
```

```
##  2.5% 97.5% 
##   2.8   6.2
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
##    [1] 2.2 3.9 4.7 6.5 5.4 3.8 4.2 4.7 3.4 5.1 5.4 4.7 4.8 3.9 6.0 4.9 4.6 3.4
##   [19] 4.7 4.6 3.4 5.4 4.2 4.5 4.7 4.7 6.8 4.2 5.9 6.7 4.3 4.1 5.9 3.4 4.9 4.1
##   [37] 5.3 3.6 5.9 5.0 5.0 4.5 5.5 3.4 4.5 5.7 4.5 4.2 4.9 4.4 6.2 7.2 5.5 4.4
##   [55] 5.6 4.0 3.0 5.7 4.0 4.4 3.4 4.7 4.8 3.6 3.4 3.6 5.1 5.0 4.4 4.3 5.8 4.7
##   [73] 2.6 5.2 2.9 3.4 6.5 4.6 4.5 4.1 5.7 4.8 2.9 4.0 3.5 4.5 3.9 6.0 3.9 3.5
##   [91] 3.7 3.6 4.6 5.3 2.4 3.9 2.6 5.4 5.4 6.1 3.6 5.3 4.0 4.4 5.6 5.3 4.0 3.2
##  [109] 3.4 4.2 4.1 3.4 4.1 3.3 5.7 2.9 5.1 5.4 4.4 3.8 5.0 4.6 2.3 4.1 3.5 4.0
##  [127] 3.5 4.6 4.1 4.7 5.4 5.7 4.7 4.7 4.3 5.1 2.4 5.1 4.6 4.9 4.5 4.6 4.2 4.5
##  [145] 5.4 6.3 3.9 3.7 5.0 4.5 4.9 4.4 3.5 5.2 5.0 4.6 4.3 5.3 3.4 3.0 4.1 4.5
##  [163] 4.9 3.1 2.7 3.8 5.6 3.3 4.4 6.3 2.6 4.5 5.9 5.1 3.8 4.6 4.4 3.9 4.8 2.3
##  [181] 4.5 4.0 6.3 3.1 6.2 2.7 4.2 3.7 4.9 6.7 3.7 6.5 4.1 4.4 5.1 5.8 4.4 4.1
##  [199] 3.5 4.4 4.5 3.9 3.2 4.0 4.3 4.9 3.5 5.1 5.8 2.8 4.5 5.0 4.7 4.7 3.8 5.8
##  [217] 4.1 4.3 3.9 4.9 4.0 2.9 7.0 4.7 2.6 4.4 4.8 2.7 4.1 4.3 5.8 5.9 5.8 4.4
##  [235] 4.2 4.8 3.4 4.5 5.3 4.6 3.5 4.1 3.7 3.6 3.7 3.0 4.0 3.7 5.1 3.5 3.3 5.1
##  [253] 5.5 4.6 5.3 2.5 3.7 3.8 5.2 4.8 4.9 3.7 4.0 3.7 4.3 4.2 4.2 5.1 5.5 2.9
##  [271] 5.1 4.4 6.0 5.3 4.0 4.6 4.4 5.3 2.9 3.2 6.4 3.4 6.0 5.5 4.6 4.2 5.4 3.3
##  [289] 4.9 3.7 4.0 3.4 4.3 4.6 5.6 4.6 5.4 4.7 4.9 4.1 4.3 4.5 4.6 5.9 4.6 5.0
##  [307] 4.0 4.4 4.9 3.2 4.3 5.0 4.7 4.7 2.6 3.0 4.3 3.9 4.6 6.4 6.0 2.5 4.1 4.6
##  [325] 3.1 3.1 4.2 2.9 4.9 3.9 2.7 4.2 3.6 4.3 5.1 4.8 2.7 5.6 3.0 4.0 4.1 4.7
##  [343] 4.8 4.9 4.2 4.4 5.0 5.1 6.1 3.6 5.6 3.3 5.4 5.2 4.8 3.6 5.6 4.1 4.4 3.6
##  [361] 4.7 5.4 5.5 4.6 4.6 3.7 3.7 3.9 3.9 5.1 2.8 4.1 3.6 5.8 5.5 5.1 5.0 4.5
##  [379] 5.5 5.7 4.2 4.6 5.8 2.8 3.3 4.1 4.2 5.5 4.6 5.2 4.5 5.4 3.9 4.9 4.7 3.2
##  [397] 5.3 3.8 3.5 3.7 4.5 2.5 6.1 4.3 4.0 3.3 3.5 6.0 4.9 4.2 4.5 3.8 5.7 3.9
##  [415] 4.2 5.1 2.3 5.4 5.6 4.9 3.5 5.5 5.1 4.5 4.4 4.2 5.2 6.0 4.4 5.3 5.7 5.5
##  [433] 3.7 4.5 5.2 4.6 4.1 3.5 5.2 3.2 3.6 6.7 4.0 4.0 3.7 3.9 6.2 4.1 4.2 4.8
##  [451] 4.0 3.4 3.6 4.0 4.0 4.7 4.7 5.1 2.9 4.6 4.4 3.4 6.0 4.7 3.7 5.9 4.9 2.3
##  [469] 4.6 3.8 5.8 5.0 6.3 6.2 4.1 4.6 3.8 4.7 3.9 5.4 5.4 4.4 2.9 5.6 3.6 4.8
##  [487] 5.8 6.5 4.4 4.2 4.2 2.9 7.2 4.9 4.0 5.0 5.3 4.1 3.9 4.3 6.4 3.9 5.7 6.3
##  [505] 3.3 5.1 4.9 5.1 4.7 4.3 4.5 4.9 4.6 3.3 4.6 5.7 4.8 5.0 4.7 4.0 5.2 4.9
##  [523] 3.5 3.3 3.3 4.3 3.6 5.0 4.7 4.6 2.3 4.9 4.5 3.2 5.1 3.9 3.6 4.7 2.8 4.2
##  [541] 6.6 4.3 5.0 6.1 5.0 4.1 4.0 3.2 3.9 3.7 4.3 4.9 5.4 3.8 4.3 6.7 5.7 5.8
##  [559] 4.3 2.9 3.9 5.5 5.5 4.1 4.3 4.4 3.1 3.7 5.1 6.0 4.4 2.4 5.2 4.3 4.9 4.3
##  [577] 5.1 5.2 4.3 5.0 4.7 5.8 4.0 5.5 3.9 4.7 4.6 4.4 4.3 5.1 6.0 3.2 5.8 3.2
##  [595] 3.5 5.4 3.1 4.5 5.8 3.9 4.2 3.4 4.8 4.1 3.2 4.3 4.0 3.9 6.4 4.8 5.0 5.0
##  [613] 4.5 3.7 5.0 3.5 4.5 4.7 5.1 3.6 5.9 3.6 5.4 5.1 4.6 5.3 3.9 4.0 3.4 3.7
##  [631] 4.2 4.8 4.7 3.4 4.6 4.3 4.1 5.5 3.5 3.7 4.4 3.4 4.5 5.4 3.8 4.2 5.5 4.3
##  [649] 5.3 3.7 5.2 4.2 4.4 1.8 4.8 3.7 5.9 4.4 4.5 2.9 3.9 4.5 4.1 5.4 4.8 4.0
##  [667] 5.1 5.3 6.3 4.0 4.9 4.2 5.0 6.8 4.2 5.9 5.3 3.4 4.9 4.2 4.3 3.5 4.7 6.5
##  [685] 4.3 6.0 4.1 4.7 6.5 3.7 4.1 4.7 5.4 4.9 3.3 2.4 5.2 2.0 4.3 4.0 4.7 4.2
##  [703] 6.9 5.4 4.6 2.9 4.6 4.5 6.0 4.5 4.1 4.4 3.5 5.2 4.3 5.5 3.1 4.3 3.0 4.7
##  [721] 4.0 4.0 5.3 3.5 4.5 4.2 4.3 4.0 5.7 5.4 6.3 3.0 5.3 5.7 4.6 5.2 4.5 5.2
##  [739] 3.7 4.3 4.0 4.4 5.2 5.0 4.2 4.2 6.3 4.6 4.2 4.7 4.8 4.2 3.0 4.5 4.6 4.7
##  [757] 5.9 4.9 3.8 4.1 5.9 5.1 3.7 3.7 4.0 2.8 4.5 4.4 3.9 4.7 2.0 5.5 4.3 5.7
##  [775] 5.8 5.9 3.3 5.1 4.8 6.4 4.3 4.6 4.6 5.2 3.2 3.6 4.0 4.6 3.4 3.8 5.0 6.0
##  [793] 5.3 5.0 4.1 3.8 4.6 4.3 4.7 3.9 4.8 4.6 3.0 4.4 4.5 3.9 4.8 5.1 4.8 3.3
##  [811] 5.7 3.6 5.1 3.4 4.6 7.0 2.8 5.7 5.7 4.3 5.6 4.6 5.4 5.7 4.6 3.9 4.3 4.4
##  [829] 4.8 5.7 3.5 3.7 5.2 5.5 6.6 5.3 4.5 5.1 4.9 4.0 4.1 3.4 5.0 4.4 5.5 3.4
##  [847] 2.8 4.6 4.7 4.1 4.4 3.7 4.0 3.3 4.0 3.7 4.4 4.9 4.5 5.0 4.0 3.5 5.2 4.4
##  [865] 3.7 3.1 4.3 3.8 3.9 3.3 4.8 3.2 6.6 4.9 5.1 4.5 4.1 5.2 4.7 3.9 3.4 3.9
##  [883] 4.9 4.2 2.6 5.1 3.8 3.4 5.4 5.1 4.7 5.2 4.2 5.5 4.1 3.7 3.7 5.5 4.8 3.7
##  [901] 4.7 5.4 3.5 5.2 4.9 5.0 2.5 4.5 3.2 4.6 6.4 3.8 4.7 3.7 5.5 6.0 4.4 4.3
##  [919] 3.7 2.9 4.1 5.3 3.2 3.9 5.2 4.9 4.9 5.2 4.1 4.9 5.0 3.4 5.3 3.6 3.7 3.8
##  [937] 2.3 3.6 4.0 2.8 3.6 4.0 3.9 4.4 3.3 5.1 4.3 5.0 4.6 4.6 4.0 4.7 6.0 3.5
##  [955] 5.5 3.9 5.6 5.5 5.7 2.8 5.0 5.3 5.7 3.4 4.7 4.1 3.7 5.0 5.4 4.6 3.0 3.4
##  [973] 3.2 4.4 2.8 4.0 2.6 4.7 5.5 4.6 4.8 4.7 4.2 4.9 3.8 4.2 4.5 2.6 4.5 4.5
##  [991] 5.0 5.7 4.7 3.9 4.5 3.3 3.8 4.8 3.9 4.6
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
##    [1] 3.7 4.3 4.1 4.2 4.8 4.6 3.6 5.6 6.3 5.6 4.1 3.5 5.0 5.2 3.4 3.7 4.0 4.3
##   [19] 5.4 5.9 4.1 5.2 4.7 4.1 3.5 5.1 4.4 2.8 5.8 5.0 4.3 3.5 4.2 5.0 4.4 2.5
##   [37] 4.3 4.3 6.4 4.9 4.3 3.7 5.0 6.0 4.9 5.9 4.5 3.7 7.1 4.9 5.7 3.2 5.0 4.1
##   [55] 5.2 4.0 3.2 4.6 4.3 4.3 4.9 3.3 3.3 4.8 6.9 4.5 4.9 4.4 3.6 4.5 4.7 4.4
##   [73] 3.8 4.8 4.1 5.1 5.7 4.7 5.4 3.6 6.0 4.6 4.3 3.9 5.0 4.0 4.0 6.8 5.1 3.9
##   [91] 5.0 3.4 5.9 3.4 4.4 4.5 5.8 4.2 5.5 4.3 5.0 4.5 5.2 3.7 3.2 3.8 5.0 5.9
##  [109] 6.0 4.5 5.1 3.8 4.7 4.9 4.6 2.9 3.3 4.0 5.4 4.2 3.8 4.8 5.6 5.4 4.0 5.3
##  [127] 4.1 4.3 5.0 4.3 6.1 3.9 4.5 4.0 4.2 4.2 5.7 4.0 3.1 3.6 4.5 5.8 3.6 5.5
##  [145] 5.3 4.3 5.3 3.5 4.8 3.7 6.8 4.9 3.3 5.6 2.8 4.5 4.9 5.1 3.4 4.9 4.8 4.2
##  [163] 5.0 4.7 5.6 5.3 4.8 3.0 5.5 3.4 2.6 3.6 4.8 4.6 5.1 4.6 5.1 5.4 5.5 3.7
##  [181] 5.2 4.3 4.0 4.3 5.1 5.4 4.4 6.6 5.3 5.1 3.9 3.5 6.0 4.5 2.6 5.4 3.9 5.3
##  [199] 4.4 4.9 4.4 6.0 4.6 4.9 4.1 3.3 4.6 4.4 4.4 3.1 3.9 3.9 4.4 4.3 4.3 3.1
##  [217] 4.6 5.5 5.2 2.8 3.8 3.0 4.6 4.2 3.2 3.3 5.8 3.8 2.7 3.8 4.5 4.9 3.9 2.9
##  [235] 5.6 4.3 3.7 4.9 3.9 4.4 4.8 5.0 5.6 4.6 4.6 4.2 4.8 4.6 4.2 4.2 2.3 5.6
##  [253] 4.1 4.4 2.8 4.1 4.7 4.5 4.2 3.5 4.1 3.0 4.1 3.0 5.3 3.6 4.4 4.8 4.3 5.0
##  [271] 3.9 5.2 5.9 3.9 5.4 5.0 5.2 5.3 4.7 3.0 3.2 4.7 4.0 4.7 3.3 5.8 5.0 4.2
##  [289] 7.0 4.8 4.4 4.6 5.1 3.5 4.4 4.3 4.7 4.2 4.6 5.8 4.5 4.4 5.0 5.7 5.8 4.4
##  [307] 3.9 4.1 4.7 6.6 4.6 4.7 5.4 2.6 4.1 5.9 5.2 4.9 5.4 5.2 4.1 4.2 4.6 5.5
##  [325] 5.4 4.0 5.6 5.1 4.3 5.4 4.1 3.5 2.7 4.8 3.2 4.5 3.7 5.4 3.4 4.8 3.5 2.9
##  [343] 5.4 5.0 4.0 4.7 3.8 5.2 3.8 5.0 5.3 3.7 3.8 4.1 3.5 6.0 3.1 4.2 4.9 3.4
##  [361] 4.2 4.1 3.6 4.0 4.9 4.7 4.2 3.8 5.3 3.1 5.9 4.4 4.8 4.0 4.9 4.9 4.8 3.1
##  [379] 6.8 5.1 4.2 5.3 3.8 4.8 4.8 5.0 4.2 5.1 4.5 4.6 4.6 5.5 3.9 3.4 3.7 4.9
##  [397] 4.7 4.9 2.9 4.7 6.4 6.0 4.7 4.4 3.4 3.3 4.6 3.4 4.8 4.7 5.1 5.4 3.2 5.1
##  [415] 5.9 4.5 5.1 3.8 5.9 5.0 4.5 5.4 3.7 5.4 4.3 5.5 4.8 3.5 5.0 4.5 4.9 5.8
##  [433] 4.0 2.8 5.4 4.6 6.4 3.8 7.0 4.5 4.3 4.4 4.0 4.6 4.1 3.1 4.1 3.7 4.0 4.3
##  [451] 2.8 3.5 2.7 4.4 3.9 3.2 3.1 6.0 3.6 5.3 5.1 5.3 5.9 4.6 5.0 3.2 4.7 5.4
##  [469] 2.6 4.5 3.6 5.1 4.3 4.8 4.7 3.6 4.2 4.7 4.4 4.6 5.5 5.3 4.4 4.1 4.8 6.0
##  [487] 4.0 4.0 4.0 3.7 6.0 3.8 5.0 5.6 4.9 6.9 5.3 4.3 3.8 2.3 3.5 4.1 3.7 4.4
##  [505] 3.9 3.7 3.3 4.3 4.1 4.2 3.5 4.7 4.2 5.2 4.6 3.1 5.0 4.9 6.5 2.8 4.7 4.5
##  [523] 6.1 4.6 5.0 5.2 3.2 4.8 5.9 2.4 2.4 5.1 4.1 5.0 6.0 3.6 3.9 4.4 4.2 5.8
##  [541] 4.2 4.5 4.8 5.1 5.7 5.6 3.9 4.1 2.9 4.8 5.1 5.3 4.0 3.4 5.1 4.8 4.3 4.3
##  [559] 3.6 3.8 3.3 5.1 6.1 3.0 3.2 3.8 5.2 4.7 5.1 5.2 2.7 3.8 4.6 5.6 5.4 5.2
##  [577] 4.0 5.8 3.3 5.0 5.1 4.8 5.8 4.3 4.7 5.2 6.6 6.9 2.9 3.7 4.3 5.6 5.0 3.6
##  [595] 4.4 4.0 4.5 3.2 3.6 4.0 4.5 5.8 5.4 5.0 4.0 6.0 4.8 5.4 5.3 5.9 4.2 4.8
##  [613] 4.9 4.6 4.5 4.3 4.8 4.0 4.7 4.8 7.2 5.0 6.9 6.2 3.7 4.7 6.6 5.5 6.4 3.8
##  [631] 3.7 3.0 5.0 4.6 5.0 5.7 4.4 4.0 4.1 4.3 5.1 5.1 5.8 4.7 2.4 4.6 4.9 4.4
##  [649] 5.9 6.2 4.0 4.6 4.4 4.5 5.4 4.8 5.3 3.9 4.4 5.6 4.5 4.7 3.2 4.7 4.6 5.1
##  [667] 3.5 3.7 3.9 4.2 2.7 5.9 4.0 4.3 4.7 4.5 4.4 4.4 4.6 6.3 4.0 5.3 3.8 3.4
##  [685] 5.4 5.1 3.1 4.5 3.6 4.6 5.3 4.1 3.7 3.4 5.7 4.7 4.6 5.0 3.2 5.4 3.7 3.9
##  [703] 4.3 2.9 4.6 2.9 3.8 4.6 4.8 5.9 4.2 4.4 5.7 4.2 3.3 3.9 5.5 4.0 4.5 3.9
##  [721] 5.1 4.3 3.8 4.4 5.0 5.4 6.1 5.0 3.9 4.0 4.2 3.2 4.3 5.9 6.9 3.1 3.6 4.5
##  [739] 3.4 4.9 2.7 4.4 6.9 3.9 4.3 3.0 4.8 3.7 3.4 4.6 6.2 4.7 6.6 3.8 4.9 5.7
##  [757] 5.6 3.8 2.5 5.8 4.8 5.4 4.1 5.2 4.8 5.3 2.8 4.2 2.3 4.6 6.1 5.9 5.4 4.6
##  [775] 6.1 5.2 5.4 5.1 3.3 4.1 4.7 5.3 3.7 3.2 4.0 3.7 4.6 4.4 4.3 5.9 5.6 4.7
##  [793] 5.8 5.5 2.2 3.9 4.9 3.8 2.6 2.7 4.2 4.9 4.5 5.1 4.4 4.4 6.1 6.2 4.0 4.7
##  [811] 5.0 2.7 5.2 4.7 4.2 3.1 5.3 3.8 3.7 4.7 6.2 3.2 3.1 4.3 3.9 4.6 4.8 5.8
##  [829] 4.7 4.4 4.1 3.5 4.3 4.6 3.9 2.5 4.8 4.1 3.0 3.7 3.7 2.7 4.8 3.0 5.0 3.5
##  [847] 5.0 4.4 4.9 3.8 4.8 2.9 4.9 5.0 2.7 4.7 4.5 3.3 5.1 4.8 4.9 5.9 5.9 4.5
##  [865] 5.0 5.8 4.4 4.6 4.6 4.0 5.9 3.5 6.5 5.6 4.0 4.8 3.0 3.4 5.5 4.2 5.6 4.2
##  [883] 3.6 4.8 4.5 4.0 4.5 3.4 4.9 3.3 3.6 5.6 3.2 4.0 4.4 4.0 4.6 3.9 2.7 4.4
##  [901] 5.2 4.1 3.6 4.8 4.2 4.9 4.3 3.1 5.5 3.8 4.7 5.7 3.3 3.8 5.2 4.2 3.4 4.3
##  [919] 2.9 4.5 4.3 4.6 5.0 4.3 4.6 4.8 4.5 3.5 3.0 4.4 6.2 3.0 4.3 2.7 5.4 3.9
##  [937] 3.7 6.3 4.9 4.4 3.4 5.1 5.8 4.4 5.0 3.5 5.4 4.7 4.5 6.2 3.7 2.4 3.5 5.7
##  [955] 6.3 4.7 5.4 3.5 4.9 4.1 5.4 4.1 4.5 5.7 6.1 4.4 4.1 4.5 4.8 4.4 5.4 4.4
##  [973] 4.8 3.6 4.7 3.9 3.5 4.1 5.0 4.8 4.1 4.5 3.0 4.7 4.2 3.1 4.7 5.2 3.0 5.8
##  [991] 3.9 4.1 4.5 5.4 4.3 3.3 5.7 4.2 4.6 3.6
## 
## $func.thetastar
## [1] -4e-04
## 
## $jack.boot.val
##  [1]  0.51029412  0.37942857  0.30919220  0.20206490  0.08176292 -0.10028409
##  [7] -0.19756098 -0.33876404 -0.36069519 -0.44466292
## 
## $jack.boot.se
## [1] 0.9670008
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
##    [1] 5.5 4.8 3.7 5.4 3.4 3.3 5.4 5.0 4.3 2.4 6.2 4.4 4.2 4.9 5.2 3.6 4.8 4.6
##   [19] 2.2 4.4 4.8 4.8 6.2 4.8 3.9 4.4 4.6 4.1 4.0 4.0 4.5 4.6 4.3 3.4 4.1 4.1
##   [37] 5.2 5.8 4.4 4.2 3.7 5.6 4.0 4.1 3.2 3.7 5.2 4.4 5.4 5.1 6.1 4.9 3.7 4.3
##   [55] 4.4 4.9 6.2 4.3 5.0 5.0 4.6 4.7 4.1 4.4 4.5 6.1 5.2 3.9 3.3 4.9 4.0 4.5
##   [73] 4.7 4.9 5.0 5.4 3.9 3.8 3.7 4.1 4.3 4.1 4.0 5.4 5.1 5.3 4.7 4.1 3.9 3.4
##   [91] 4.8 3.6 4.2 6.2 5.7 4.2 3.3 2.3 3.8 4.2 4.7 4.3 3.3 4.2 6.1 3.9 4.5 4.1
##  [109] 4.5 3.2 3.7 5.1 5.4 4.8 4.5 3.7 4.8 5.4 5.8 5.2 3.9 2.9 5.4 2.7 4.3 4.6
##  [127] 4.4 5.4 3.4 5.2 4.1 6.4 4.9 5.9 4.6 3.6 3.5 4.3 4.5 4.0 3.8 5.5 2.6 4.1
##  [145] 6.0 5.5 5.3 4.0 4.4 4.5 4.4 4.1 4.9 5.7 5.7 4.4 4.1 4.0 4.8 4.1 3.0 4.9
##  [163] 6.0 3.1 6.2 4.2 4.4 3.4 3.2 4.8 6.0 3.0 5.1 3.9 6.6 3.2 5.5 4.8 4.8 5.5
##  [181] 4.5 4.7 3.3 4.9 2.6 5.6 4.3 4.9 4.2 3.3 4.0 2.6 3.9 5.4 5.9 3.5 5.4 6.1
##  [199] 5.1 4.4 2.8 4.4 3.6 4.3 5.4 6.1 5.5 4.7 5.1 5.0 3.5 5.8 6.4 4.5 4.8 4.4
##  [217] 4.6 5.4 5.2 4.4 3.5 3.6 3.3 4.0 5.9 2.9 5.3 6.1 5.9 4.9 5.7 5.1 4.3 4.1
##  [235] 4.1 5.4 4.3 4.8 3.1 6.4 4.6 4.4 5.3 5.4 4.9 5.1 2.9 6.1 4.8 6.9 5.2 3.4
##  [253] 4.6 5.1 5.7 4.8 4.5 5.4 4.5 4.8 4.1 5.1 3.7 3.9 4.5 3.2 5.3 4.6 5.3 2.9
##  [271] 5.0 5.3 4.8 3.3 4.5 4.2 5.4 3.8 4.8 4.1 5.5 3.3 4.7 5.7 4.5 4.1 4.5 5.7
##  [289] 4.9 3.2 5.9 4.6 6.0 4.7 3.9 3.4 5.6 4.9 3.3 4.3 4.3 5.7 6.3 4.1 3.7 3.7
##  [307] 3.7 3.1 5.3 4.6 5.0 4.6 4.0 4.9 5.6 4.8 4.1 4.1 4.8 4.7 3.6 4.6 4.8 4.1
##  [325] 4.6 3.9 3.9 4.6 6.0 5.5 5.1 3.5 5.2 5.5 5.8 4.1 5.0 3.9 3.9 3.5 4.6 4.7
##  [343] 7.1 4.0 5.3 4.1 3.5 1.5 3.8 5.2 5.4 4.3 5.2 3.6 4.0 4.1 4.9 4.3 4.8 4.7
##  [361] 4.2 4.0 4.7 4.4 5.0 4.1 5.4 4.3 4.1 3.0 4.0 5.2 5.1 4.5 3.6 5.0 5.0 3.5
##  [379] 5.3 3.1 5.5 4.7 3.3 4.6 3.5 4.3 3.9 5.7 3.5 6.2 4.0 4.2 4.1 4.9 4.8 2.1
##  [397] 4.2 4.4 2.8 2.7 3.1 4.0 5.5 3.4 5.0 4.2 4.1 4.2 3.0 4.8 2.3 2.0 4.3 6.0
##  [415] 3.6 4.2 5.8 4.8 3.7 3.6 4.2 2.9 5.7 4.3 5.3 5.8 4.4 5.3 5.6 5.3 6.1 3.7
##  [433] 4.7 4.4 4.3 5.3 3.6 5.8 4.6 4.7 6.2 4.1 4.4 4.9 3.7 3.4 4.7 3.7 4.8 3.2
##  [451] 3.7 4.2 6.5 3.8 5.8 4.1 5.4 4.0 5.2 4.1 4.7 5.0 5.1 4.0 4.2 4.1 3.9 4.2
##  [469] 3.9 4.6 5.1 3.1 6.9 5.7 3.1 3.2 4.2 5.4 4.7 4.0 5.0 3.4 3.9 5.4 3.8 3.2
##  [487] 4.4 4.3 3.4 6.7 5.2 3.8 3.2 5.5 5.0 5.9 5.1 5.3 2.5 3.7 6.1 4.1 4.8 5.5
##  [505] 3.6 4.1 5.0 4.3 2.7 5.0 4.7 3.5 4.6 5.0 6.5 4.6 3.2 5.5 5.6 2.0 3.7 3.7
##  [523] 4.5 5.7 4.2 4.0 4.2 5.1 5.3 5.0 5.0 3.7 4.2 5.1 5.6 5.6 4.3 2.9 3.4 5.8
##  [541] 4.5 3.5 4.9 4.9 5.4 4.4 4.6 3.9 4.2 5.8 3.7 5.7 6.1 5.2 4.0 3.4 3.3 5.2
##  [559] 5.4 4.4 5.1 5.0 1.3 5.8 4.7 5.6 4.1 5.2 3.4 5.6 4.5 3.1 3.9 4.2 2.9 6.7
##  [577] 5.3 4.1 5.4 4.7 5.1 4.9 4.2 3.2 3.5 5.7 4.7 4.2 4.9 3.7 5.5 2.7 5.4 5.3
##  [595] 6.1 3.5 2.8 4.1 3.6 4.3 3.5 4.2 4.8 3.1 5.2 3.3 3.7 3.4 4.8 5.3 4.0 5.1
##  [613] 4.2 4.6 4.4 3.8 4.0 3.9 5.4 5.9 4.5 4.6 3.6 3.1 5.4 5.4 3.9 4.6 4.6 5.7
##  [631] 5.0 4.1 3.3 4.5 2.1 5.7 4.0 3.9 6.0 2.7 5.7 3.1 4.4 4.2 4.6 3.9 5.3 4.0
##  [649] 4.5 3.1 4.5 4.4 3.4 4.3 4.4 4.3 3.5 5.5 4.0 4.8 6.2 2.0 4.2 4.6 3.7 5.0
##  [667] 3.6 5.1 4.7 4.7 3.7 4.5 4.0 4.0 5.3 4.4 5.1 3.6 5.6 5.5 3.5 5.2 5.6 5.0
##  [685] 4.5 5.2 3.7 5.6 4.3 3.9 4.4 2.2 5.1 4.0 3.2 4.2 3.8 6.5 3.5 4.5 4.3 4.6
##  [703] 5.0 4.2 3.6 4.0 5.5 3.0 3.7 3.1 5.2 3.3 5.3 4.5 5.5 4.3 3.5 6.2 5.5 5.7
##  [721] 5.4 5.5 4.3 4.1 3.4 5.7 3.2 3.0 3.6 3.5 4.9 3.4 2.8 5.4 2.8 5.2 5.1 3.3
##  [739] 3.9 4.7 3.7 5.5 3.5 4.8 3.9 4.2 4.7 4.2 3.7 4.7 4.0 4.3 2.3 5.2 2.7 3.8
##  [757] 4.7 5.3 5.4 3.6 4.3 5.0 6.4 5.5 4.2 5.3 5.1 4.3 4.8 3.8 5.7 3.4 4.3 4.4
##  [775] 2.7 2.8 4.0 5.7 5.4 4.2 4.2 4.9 3.7 4.8 4.9 4.7 4.6 4.3 5.2 4.0 3.7 3.6
##  [793] 4.4 4.1 3.9 4.9 4.5 2.6 4.8 4.6 4.4 4.2 4.2 4.9 4.8 5.3 4.5 4.8 3.3 4.2
##  [811] 5.0 3.2 5.1 2.9 5.6 4.3 3.0 3.6 4.5 3.9 3.8 5.7 5.2 4.4 4.7 6.0 5.6 5.8
##  [829] 3.1 4.9 5.3 4.1 4.2 4.4 5.4 3.7 4.8 4.9 3.8 6.3 5.0 3.8 5.3 3.8 3.1 4.8
##  [847] 4.3 3.4 3.4 4.0 4.0 4.0 5.8 5.0 6.7 4.4 5.0 6.4 6.9 4.6 6.6 5.5 4.3 4.0
##  [865] 4.3 5.5 4.4 4.0 5.7 3.8 3.5 3.9 4.5 5.0 4.1 5.4 5.6 4.9 3.2 6.3 5.7 5.0
##  [883] 4.7 4.3 4.3 3.6 2.3 5.7 4.8 4.4 4.1 4.0 4.8 4.0 4.7 5.0 4.4 6.5 2.5 4.0
##  [901] 4.4 6.3 4.2 2.9 4.4 5.4 3.9 4.3 3.9 3.6 3.1 3.8 3.7 4.0 4.4 5.3 3.5 4.4
##  [919] 2.4 5.4 4.8 4.9 5.0 5.4 2.9 4.5 2.0 3.3 3.6 3.9 4.1 3.8 4.6 5.0 5.4 4.9
##  [937] 4.1 4.2 5.3 5.9 3.9 3.6 4.7 3.6 4.1 3.7 5.7 5.1 3.6 4.8 4.4 4.1 3.5 4.6
##  [955] 4.1 4.6 4.2 4.3 4.6 4.9 2.5 4.9 5.4 4.5 5.9 4.8 4.4 3.6 3.6 3.5 5.5 4.2
##  [973] 4.0 2.0 4.8 3.3 4.8 3.7 4.2 4.0 4.5 4.4 4.8 4.3 5.7 3.4 4.0 2.6 4.2 4.6
##  [991] 4.1 5.2 3.7 3.8 5.3 2.6 5.7 4.2 3.2 5.2
## 
## $func.thetastar
## 72% 
##   5 
## 
## $jack.boot.val
##  [1] 5.500 5.400 5.300 5.232 5.100 5.000 4.900 4.700 4.544 4.400
## 
## $jack.boot.se
## [1] 1.052646
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
## [1] 1.808106
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
##   2.885895   3.216770 
##  (1.223376) (1.489318)
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
## [1]  0.91085231 -0.41348082  0.97036706 -0.03730647  0.72896664  1.47288477
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
##    [1]  0.820964222  0.592082201  2.133534592  1.431652893  0.381179484
##    [6] -0.565943162  1.966697218  1.097484161  2.090918890  0.939076605
##   [11]  1.747864532  0.010903201  1.294523007  1.150704280  1.717013758
##   [16] -0.883666661 -0.015865966  2.446783943  1.359275646  0.730082263
##   [21]  2.196767593  1.751764092  0.828936710  1.255356490 -1.338488186
##   [26]  0.282072425 -0.606257259 -0.637960529  0.769717355  0.616456777
##   [31]  1.808106221  2.491935028  2.125865242  1.504572751 -0.076589891
##   [36]  1.754415757  0.569421813  2.086191421  1.260058792 -0.657106383
##   [41] -0.842560165  1.100972960  1.490782241  0.589221959  2.150543428
##   [46]  2.067952023  2.511617587  2.434974950  0.356271112  2.455156344
##   [51]  1.084160597 -0.144340553  2.473813249  0.236554059 -0.853788214
##   [56] -1.361673394  1.619666420  0.916971120 -0.476893707 -0.752831770
##   [61]  0.380840073  2.113936650  1.491827576  1.284987371  0.136423769
##   [66]  1.554865857  1.234405122  1.266177564  1.829800853  0.284668560
##   [71]  2.175320770  1.409904562  1.213819360  1.836575017  0.609851676
##   [76]  0.646872657  0.921061719 -1.651576758  1.869564767  1.782751262
##   [81]  2.052281708  0.531497139  1.806011205  0.949653576  1.083030049
##   [86]  1.282979385  2.049647075  2.221244323  1.351791000  1.319073294
##   [91]  0.924959125  1.489370744  1.785564256  1.410307822  1.096068381
##   [96]  1.865992412  1.239641384  1.385427076  1.913047280  1.602277126
##  [101] -1.126051103  1.531908252  0.171125569 -0.044352660  0.467974567
##  [106]  1.780513425  2.184937075  0.979058338  0.686044980  1.230588913
##  [111] -0.650906871  1.450287864  2.236284900  1.319472647  0.938597284
##  [116]  0.705538965  1.415202237  1.871241917  1.096573164  0.160164165
##  [121]  1.294653483  1.150677497  1.074915359  1.985151756  1.234729234
##  [126]  1.226249932  0.448780495  1.431329722 -0.324145751  2.063938009
##  [131]  0.188263256 -0.682461732  1.820766746  0.154155127  2.458270024
##  [136]  0.553806647  0.791824800  1.868428259  1.045740599  1.539975185
##  [141]  1.428198176  2.097212969 -0.025778413  1.632024559  2.198678349
##  [146]  2.208764780  2.448635223  1.668496243  1.468272938  2.157488358
##  [151]  0.820799019  1.265456845  1.335330510  1.325080038  1.122993529
##  [156]  0.727283334  0.558816761  0.159041320  1.205024792  0.126393770
##  [161]  1.776977022  2.117727121  1.276905167  1.889256990  2.495197425
##  [166]  1.747839917  1.789749370  0.439009518  2.421614931  1.496846820
##  [171]  0.703482284  2.186234494  1.329788986  1.887550884  0.710294059
##  [176]  0.014652892  1.417462376  0.832019133  2.146455683  1.229522700
##  [181]  2.126916549 -0.519123032  0.678254229 -0.237909370  1.509036746
##  [186] -1.022789134  1.249909754  0.469070075  0.832354538  0.004550760
##  [191]  1.770859762  1.289919857  1.840456523  0.439163135  1.422352116
##  [196]  1.060283582  0.507392434  1.775408290  2.507178924  0.398419595
##  [201]  1.775695780  1.249811023  1.211346645  2.457164362  2.023328823
##  [206] -1.014914782  2.066004904  1.279958752  2.037953885  1.614336102
##  [211]  1.813099916 -0.523661838  0.252911375  1.301653879  1.785043339
##  [216]  2.431396744  2.150714088 -0.580649761  1.602214197  2.085519290
##  [221]  1.170632749  0.966652413 -0.454265690  0.859288029  1.563171443
##  [226]  0.259477896  2.010548202 -1.443551658 -1.185512979  1.795308852
##  [231]  0.013723828  2.045984506  0.582975599 -0.441729285  0.672563651
##  [236]  1.260385603  2.490241802 -0.460076937  1.221969298  2.172224542
##  [241]  1.411242935  1.756858599  1.325554658  1.046766915 -0.552396709
##  [246]  1.414949337  2.100258399  1.627113740  2.536010923  2.172004110
##  [251]  1.118720742  1.896947198  1.160248089  0.916749606  0.528077302
##  [256]  1.409359576  2.187819916  2.474900123  1.716130964  0.128405448
##  [261]  1.792469630  1.542392944  2.480829432  0.304242148  1.995383462
##  [266]  0.449801471  0.371085383 -1.453924579  1.946081107  0.069300515
##  [271]  1.783742188  2.147080267  1.821776624  2.244009545  1.423984996
##  [276]  0.006744387  1.629939905 -0.311178032  1.872211882  1.541565369
##  [281]  2.277407867  0.708346108  1.737382198  1.494995354  0.852687067
##  [286] -0.552396709  1.289972347  2.155334996  1.134436670 -0.153836417
##  [291]  1.092401617  2.127536181  0.159041320  1.445227402  0.337147029
##  [296]  0.140581961  0.589637994  2.189825919 -0.291696175  0.707386190
##  [301]  1.100739649  2.097830952 -0.296180219  1.362303451 -1.098793277
##  [306]  1.165718661  1.773820306  0.376232700  0.470073069  0.400453391
##  [311] -0.985983339  2.090444810  0.708226961  0.462196478  0.585376886
##  [316]  0.312461840  1.494617165  1.570540508  1.286726424 -0.770859484
##  [321]  0.211504470  0.183295860 -0.933704133  1.807071198  2.478275518
##  [326]  1.297014765  0.259199169  1.710438795  1.246299172 -0.313540606
##  [331]  0.287946765 -0.129643121  1.691451952  1.593918637  1.222957430
##  [336] -1.846201158  2.113753001  1.770366305  2.546421678  1.800510704
##  [341]  1.739582968  1.708227333  1.259842780 -0.421781186  1.024462934
##  [346]  0.965588686  1.090303213  1.666819350  0.780435854  0.996689419
##  [351]  0.756932205  2.127864873  0.708531759  1.984633969  1.973296684
##  [356]  1.875425229  0.635068833  0.087944360  1.587642505  1.302874133
##  [361]  0.637229818  1.174618301 -0.296807804  0.054059223  1.264324922
##  [366] -0.423799747 -0.056380503  1.859895948  1.851327339 -0.429100199
##  [371]  1.933970570  1.807077517  1.273917878  2.152471265 -0.877361413
##  [376]  2.126916549  1.528311952  1.307619501  0.528600934  0.239561987
##  [381] -0.929663239  0.358614501  1.062688090  0.560306806  1.347156704
##  [386]  0.712482757  1.869704001  0.324939047 -0.837983413 -1.436310191
##  [391]  1.238772222  1.949113143 -1.262103035  2.001659801  1.645793856
##  [396] -0.440865725  2.563677388  1.496189964  1.771492306  2.449658488
##  [401]  1.264473632  2.185131863  1.000424361  1.472893690  1.246714495
##  [406]  2.217331011  0.211332883 -0.605239788  0.254043459  1.955627314
##  [411]  2.510492730  0.719803428  1.351791000  2.126507459  2.489284250
##  [416]  2.267071443  0.941521580  1.424998765  0.835382676  1.617316818
##  [421]  1.748010441  2.440289743  1.073962017  2.141717721 -1.162637354
##  [426]  1.428651266 -0.063304497  1.692137740 -1.296226954  1.502204350
##  [431] -0.278756047  1.543152219  0.708221722  0.416751170  1.828992545
##  [436]  2.520283626  0.867855971  2.080086253  2.161950643  0.513509686
##  [441]  0.450086611  0.351512569  2.255399439  1.818711708  0.357296201
##  [446]  1.522709542  1.426008570  2.468830833  1.554291179  2.257853718
##  [451]  2.540925952  0.688852908  1.054658788  0.686252148  2.408141725
##  [456]  2.016727861  1.760680607  1.145463737  1.122356883  2.118834349
##  [461]  0.187205738  2.207070799  1.622404467 -1.221521858  2.209498247
##  [466]  1.763628714  1.488527094  0.691064041  0.612147804  1.810438127
##  [471]  1.100873603  2.081930537  0.288233771  0.632212708  1.621867479
##  [476]  1.641323693  1.650880910  1.932854135  0.965588686  0.566506124
##  [481]  1.576700882  1.267477622 -0.509398155 -0.692900452  1.281940255
##  [486]  1.392027282 -0.627836208 -0.571223237  0.803006575  0.221379000
##  [491]  1.523228752  0.462351036  1.501312981  0.875059416  0.593930468
##  [496]  0.574258253  1.769799335 -1.123345704  2.102362058  0.593389705
##  [501]  2.132650672  0.254553265  1.104165720  1.846603659  0.007750072
##  [506]  1.755739515  0.198125588  2.105476284  0.356995853  1.110090667
##  [511]  1.294334119  2.014521247  2.161708725  1.857045289  1.472804239
##  [516]  2.466979633  1.397258827  1.039742600  1.015548848  1.012550712
##  [521] -0.910911889  0.939245057  2.238105593  1.773086603  0.915226728
##  [526]  1.786857936  0.602693935  2.210853879  1.462879099  0.149953706
##  [531]  1.011517876  0.197190080  0.172638087  1.810655416  0.501499249
##  [536]  1.177290493  1.077070672 -0.930602220  2.319317382  1.102449625
##  [541]  1.254801558  1.487250412  1.296265279  1.439087318  1.823393743
##  [546]  0.868476695  2.504673590  1.229927910  1.420454564  1.829993225
##  [551] -1.027474112  1.212671893  0.782877829  1.080258461  1.404228089
##  [556]  1.602947076 -0.701389583  1.598985895  1.264763266  1.511389385
##  [561]  1.350151308  1.900828667  1.780217138  0.087601090  0.479496018
##  [566]  0.435002894  1.815415819  1.288479992  0.905716467  2.157488358
##  [571]  1.454485345  1.984454587  1.122614021  1.822128115 -0.350729799
##  [576]  1.995649057 -0.355081325  1.808106221  1.249470660  0.713293899
##  [581]  0.268480534  0.724581534 -0.364623456  1.191129628  0.882203574
##  [586] -0.546816065  0.741720993  2.191597301  1.704317808  0.314410942
##  [591]  1.265518252  1.535551998  1.995875733  1.788153788  0.485266506
##  [596]  1.668496243 -0.084878424  1.824010134  1.918261685  1.014294115
##  [601]  1.308691423  2.519271938  0.143798000  1.777795202  0.609297448
##  [606] -0.112607842  0.938473650  1.348062433  1.219169363  2.127515923
##  [611]  0.912932390  0.435608856  1.362890261  0.856562241  1.074162292
##  [616]  1.142919080 -0.680573004 -0.615279569 -1.251756244  0.166761262
##  [621]  0.405587250  0.709419295  1.096068381  1.486491938  2.171248083
##  [626]  2.081568771  0.558112966 -1.197294347 -0.467103640  2.187252125
##  [631]  1.110346321  1.759408709  0.018061867  0.736069756 -0.420705231
##  [636]  1.266267442  0.560265381  1.280835981  1.879132823  1.801908002
##  [641]  1.296505568  1.806832813  1.081153243  1.225135766  1.705948252
##  [646]  0.058556314  2.158379765  2.143753724 -1.587227024  1.521734603
##  [651]  0.486918826  1.740075952  1.265135745  2.154459039  1.031266321
##  [656]  2.114783671  1.300606960  0.969114340  0.930442999  2.105149182
##  [661] -0.208541230 -0.767181962  2.101199829  0.181890939  1.328814374
##  [666]  0.930955362  1.464826783  2.017505132  0.015902296 -0.102016441
##  [671]  0.060116552  2.164008718  1.066622604  2.156428436  1.776131556
##  [676]  2.121783027  0.380749802  2.489908766  1.098603588  0.011004105
##  [681]  1.880248536  1.787412533  0.922708991  1.125772753  2.536723299
##  [686]  1.885540732  0.635942755  0.418263115 -0.999628766  0.432592981
##  [691]  1.754170509  0.253261162  0.994633874 -0.566389911  2.407967068
##  [696]  2.016449879  2.498750918  1.091242100  0.683029984  1.261931410
##  [701]  2.136403733  0.748367270  0.591767961  1.136818830  1.759089114
##  [706] -0.683450495  0.924593977 -1.148282147  1.808848398  1.754170509
##  [711]  2.310955956  1.264143401  1.274668825  1.092686702 -1.388418037
##  [716]  1.354187896  0.067829134  2.152968638  1.278616358  0.505708061
##  [721]  0.654003832  1.513941135  2.072291120  0.799395145  0.032004710
##  [726]  1.710507396  2.507399812  1.362874929  1.194403103  0.944394873
##  [731]  1.164566704  1.060272271 -1.242944336  1.786378648  0.514181734
##  [736]  1.747545871  1.249407970 -0.063239719  1.415343072  2.155737525
##  [741]  2.550863187 -1.304112853  0.626504101  1.415564274  0.624256569
##  [746]  0.198148488  1.234470538  2.096706556  0.463727732  2.213712008
##  [751]  0.495478913  1.269973205 -0.192693440  0.295501597  1.918245431
##  [756]  1.777795202  1.061608570 -0.610364041  0.066378835  0.705984794
##  [761]  0.126397237  1.664824255  1.520611538  0.322653485  0.838919855
##  [766]  0.281198180  2.412260562  1.064415568  0.738182299  0.240950231
##  [771]  1.883818350  1.806716762 -0.527774165  0.439009518  1.399171861
##  [776]  1.132792135  1.118720742  1.783778373  1.782112946  1.716489382
##  [781]  1.838157388  1.979085139  1.878434661  0.239914745  1.505470644
##  [786] -1.013983688  1.398590944  2.246136838  0.143366764  1.532573119
##  [791] -0.147122524  0.556192169  0.679457932  0.935280267  1.794669695
##  [796]  0.916880416  1.422756935  0.771575753  1.340375474  0.320928678
##  [801] -0.422442594  0.760124993  0.889159808  0.252050427 -0.901017360
##  [806]  1.226610324  0.610067170  2.085720293  1.202545524  2.128665124
##  [811]  0.735458841  1.257598653  2.166176326  0.812890814  1.353602233
##  [816]  0.160164165  0.499308412  0.738149730  0.628105140  2.129037754
##  [821] -0.893908463  2.088432731  0.234125109  2.396512369  1.241390675
##  [826]  1.973074895  0.494103155  2.502643271 -0.275897379  2.195075960
##  [831]  2.038509537  0.404808130  1.406135587  2.514806115  0.530713552
##  [836]  1.160641168  1.097004121  1.101020981  1.263311542 -1.478892065
##  [841]  1.858007210  2.151445266 -0.608306047  1.055982540  1.508849844
##  [846]  1.131105044  1.311468883  0.768557620  0.624457453  1.305534544
##  [851] -0.988041635  2.332291457  1.725626987  1.926194030  0.325327525
##  [856]  1.787554843  1.272006099  0.215560417  1.510455969 -0.722365648
##  [861]  1.261542862 -1.241190612  1.027339517  2.478836567  1.909249976
##  [866]  2.128972979  0.752896490  1.466371921  1.532548880  0.825706552
##  [871]  2.511617587  1.716334130  0.154461567 -0.706474400  0.685936463
##  [876] -1.257261507  1.835736049  0.008093252  0.092852481  1.459132715
##  [881]  1.174212581  1.382763613  0.531257515  1.761685977  1.846603659
##  [886]  2.176627960  2.495024237  2.217331011 -0.439566111  1.745510826
##  [891]  1.872416839  0.829694240  1.685000200  2.170531414  0.024866214
##  [896]  0.495810391  1.888996260  2.248786385  0.353913065  0.262119678
##  [901]  0.751341019  0.867226387  1.428935037 -0.056999002  1.510974555
##  [906]  0.527565862  0.669835138  1.370040231  1.255553808  1.105392856
##  [911]  1.090682686  1.278109310  0.402964905  0.071701300  2.455982029
##  [916]  2.177864708  1.284412699  0.808671808  1.284053012  1.498221368
##  [921]  1.627499744  0.963001230  2.107458187  2.111680842  1.271481923
##  [926]  2.128976799  1.089656093  1.251738221  2.121435382  2.116208332
##  [931]  1.610263527  1.826918651  1.604869432 -1.900321697 -1.672702725
##  [936]  1.241950553  1.252566664  1.258028978  2.041362052  1.342491572
##  [941] -0.734560883  1.483853920 -0.990879468  1.533112220  0.756461830
##  [946]  1.780585299 -0.227504132  1.236929388  1.602340630  1.077352151
##  [951]  1.482318077  0.366580200  0.124751681 -0.013348783  2.169927248
##  [956]  0.747972686  2.049244087  1.660929442  0.608825420  1.598230218
##  [961] -0.859396936  1.256345997  1.805457156  0.367954434  1.616008882
##  [966] -0.788732942  2.064843444  1.422630976  1.428612988 -1.726231297
##  [971]  1.241769001  2.098600177  2.204122897  2.179466084  1.746373255
##  [976]  1.418967970  0.828075259  0.031376518  1.421105954 -0.605982620
##  [981]  2.139585160 -1.205562867  1.254120811  1.750121504  1.636355409
##  [986]  1.562395173 -0.644463014  0.782342903  2.104429131  0.761871575
##  [991]  1.284042605  1.419596573  0.395905203  1.107832871  0.400042862
##  [996]  0.722761073  2.207695213  2.043336999  1.238398354  1.344813162
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
```

```r
fit2
```

```
##      mean         sd    
##   0.8971327   0.6068371 
##  (0.1918987) (0.1356900)
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
## [1]  0.20830892 -1.56446793 -0.10760717  0.08272764 -0.45352623 -0.43754721
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
## [1] -0.0145
```

```r
sd(results2$thetastar)  #the standard deviation of the theta stars is the SE of the statistic (in this case, the mean)
```

```
## [1] 0.8992933
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
## t1*      4.5 -0.01101101   0.9285706
```

One of the main advantages to the 'boot' package over the 'bootstrap' package is the nicer formatting of the output.

Going back to our original code, lets see how we could reproduce all of these numbers:


```r
table(sample(x,size=length(x),replace=T))
```

```
## 
## 2 3 4 5 6 9 
## 2 2 1 3 1 1
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
## [1] -0.0245
```

```r
se.boot
```

```
## [1] 0.9086228
```

Why do our numbers not agree exactly with those of the boot package? This is because our estimates of bias and standard error are just estimates, and they carry with them their own uncertainties. That is one of the reasons we might bother doing jackknife-after-bootstrap.

The 'boot' package has a LOT of functionality. If we have time, we will come back to some of these more complex functions later in the semester as we cover topics like regression and glm.

