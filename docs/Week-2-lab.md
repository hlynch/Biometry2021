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
## 0 1 2 4 6 7 9 
## 1 1 1 3 2 1 1
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
## [1] 0.0228
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
## [1] 2.796691
```

```r
UL.boot
```

```
## [1] 6.248909
```

Method #2: Simply take the quantiles of the bootstrap statistics


```r
quantile(xmeans,c(0.025,0.975))
```

```
##  2.5% 97.5% 
##   2.9   6.2
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
##    [1] 5.1 6.0 3.5 4.6 3.7 3.0 5.5 4.5 4.0 6.7 5.1 4.2 2.6 5.4 4.3 6.3 3.5 3.9
##   [19] 7.3 4.6 4.3 4.9 4.0 2.7 4.0 4.2 4.0 3.7 4.5 4.6 5.7 4.2 6.2 4.5 3.9 5.0
##   [37] 4.4 4.6 3.8 4.0 4.7 5.1 4.2 4.2 3.8 5.1 3.4 3.9 6.2 4.9 4.5 5.1 4.4 3.8
##   [55] 4.9 4.2 5.1 4.4 5.6 4.2 5.2 2.5 4.8 4.0 3.3 4.0 3.7 2.7 4.4 2.6 4.4 2.9
##   [73] 4.9 4.5 5.4 4.3 3.4 3.7 5.0 4.1 5.5 3.5 4.9 4.1 4.9 5.2 3.4 6.7 2.7 3.8
##   [91] 4.7 5.5 5.8 4.4 4.2 4.3 4.7 4.2 4.9 4.8 5.0 4.5 5.3 2.8 3.3 5.0 5.8 3.8
##  [109] 4.8 5.3 5.2 5.4 4.6 4.4 3.1 4.2 4.7 3.8 3.7 5.3 4.8 4.3 2.1 4.5 2.9 3.9
##  [127] 4.0 4.7 4.1 3.9 3.2 4.7 4.8 4.2 5.4 4.9 3.9 3.4 4.2 4.2 3.7 2.9 4.8 4.3
##  [145] 4.7 4.1 3.0 4.0 5.6 3.6 3.8 5.3 4.7 3.9 4.7 3.9 3.6 5.8 5.6 4.5 5.6 3.6
##  [163] 4.3 6.4 3.7 3.2 4.8 4.2 3.3 4.7 2.8 4.4 5.4 4.6 3.5 6.9 5.2 4.3 4.5 5.8
##  [181] 5.9 4.9 5.2 4.7 4.4 4.4 2.9 5.2 5.3 5.0 3.7 6.6 4.4 4.3 3.0 4.2 6.3 3.7
##  [199] 5.0 5.2 4.8 5.6 4.9 4.3 4.6 6.4 4.0 2.6 4.3 3.4 3.3 4.3 5.3 5.0 5.1 4.6
##  [217] 5.2 4.7 6.5 2.5 4.7 4.8 4.1 5.9 4.1 4.4 6.5 4.2 4.4 4.2 4.0 5.9 5.6 4.0
##  [235] 5.7 4.0 5.0 3.6 4.0 4.5 5.9 3.3 4.1 3.6 2.8 5.7 4.2 4.9 5.3 5.2 4.2 5.7
##  [253] 4.3 5.9 5.3 3.8 2.9 3.1 4.4 4.7 5.5 3.5 4.7 2.1 5.4 4.1 6.6 4.8 2.9 3.8
##  [271] 4.6 5.0 4.5 5.2 4.5 6.3 4.9 5.5 4.5 4.6 5.3 4.6 4.4 3.0 2.7 5.0 5.6 5.8
##  [289] 6.0 3.5 4.9 5.0 5.9 5.0 4.7 3.6 4.0 4.7 4.9 4.3 4.3 4.3 3.6 3.6 5.3 4.1
##  [307] 6.8 4.9 2.8 4.1 3.5 5.3 3.4 4.1 4.4 3.7 5.4 5.5 4.2 4.4 3.3 4.4 5.4 3.7
##  [325] 4.2 3.7 3.9 4.5 2.5 4.4 5.4 3.1 5.0 4.7 4.4 4.9 4.6 5.4 5.6 6.1 5.1 5.7
##  [343] 5.1 5.5 5.2 5.1 5.5 5.4 3.7 6.1 3.7 3.7 5.2 4.0 5.4 4.5 4.8 5.0 5.7 4.0
##  [361] 4.9 5.3 4.7 4.7 4.7 4.1 4.8 4.8 4.9 2.7 4.2 2.7 5.1 4.6 4.6 4.3 3.9 4.8
##  [379] 4.6 5.3 2.8 5.5 5.6 5.3 4.3 4.0 4.5 6.3 3.8 6.1 4.3 4.0 3.2 3.8 3.2 5.0
##  [397] 5.0 2.4 3.0 4.2 4.4 5.9 3.4 4.3 2.0 5.3 4.6 3.9 5.8 4.2 2.0 5.3 5.5 3.7
##  [415] 3.5 3.9 5.4 4.8 3.6 4.7 6.9 4.1 4.3 4.6 5.4 4.0 3.4 5.6 4.4 3.6 4.9 4.3
##  [433] 4.6 4.0 4.3 5.5 4.7 3.4 4.2 6.5 3.9 4.5 4.6 5.7 5.2 3.3 6.5 4.6 3.2 4.6
##  [451] 5.8 4.2 4.3 4.4 4.6 4.3 4.3 4.3 3.0 6.0 4.0 3.8 4.5 4.1 4.3 3.9 3.0 5.5
##  [469] 6.0 4.1 5.0 5.0 5.1 4.4 4.7 4.3 4.6 5.4 6.0 4.6 5.5 4.7 3.2 2.1 4.1 5.5
##  [487] 6.5 4.4 3.0 4.2 3.9 4.8 3.3 4.2 5.4 2.5 5.6 4.5 2.5 5.3 5.9 3.8 3.0 3.8
##  [505] 3.4 4.0 5.3 4.3 5.7 4.2 2.8 5.5 3.3 4.0 4.2 3.8 6.5 4.6 5.9 6.1 2.4 4.0
##  [523] 3.9 3.9 4.6 3.9 4.1 4.6 5.2 4.9 4.5 4.0 4.0 5.1 5.8 4.5 2.7 4.9 4.9 6.2
##  [541] 4.0 3.7 4.8 4.5 4.3 2.5 2.7 4.3 4.2 4.1 4.2 3.6 4.0 2.6 4.4 2.6 3.4 4.0
##  [559] 4.2 3.1 4.5 4.1 5.3 4.8 5.1 4.7 3.8 5.8 6.6 4.4 3.1 3.0 4.5 2.8 3.6 4.7
##  [577] 4.5 5.5 3.9 5.3 5.5 4.7 4.4 6.4 5.1 3.6 4.7 3.7 5.0 5.6 4.2 4.9 5.9 3.8
##  [595] 4.9 6.1 4.4 4.7 4.0 4.5 2.7 4.7 5.3 5.7 5.0 2.3 4.0 5.1 5.7 5.5 4.1 2.9
##  [613] 4.5 5.4 3.6 5.1 3.9 5.4 4.8 4.6 3.6 3.8 3.9 4.5 3.9 3.6 4.6 5.3 4.1 4.7
##  [631] 6.4 4.2 6.3 4.9 4.4 2.4 4.0 4.7 4.5 4.5 7.3 4.6 2.8 3.6 3.4 5.3 5.7 4.7
##  [649] 6.0 5.2 5.5 5.4 4.1 3.8 3.6 4.4 3.9 3.2 4.0 4.5 3.2 3.0 6.5 3.5 6.7 5.0
##  [667] 3.6 4.1 3.7 3.8 4.7 4.9 5.1 4.4 3.9 5.7 4.3 5.4 4.5 5.0 3.9 3.8 5.2 3.4
##  [685] 4.4 5.2 4.5 4.5 4.2 5.3 4.9 4.6 3.4 3.4 4.4 4.5 5.9 4.6 5.1 3.2 3.5 3.7
##  [703] 4.8 5.0 4.5 3.8 5.4 6.8 3.8 4.6 4.9 5.2 6.4 5.5 4.7 4.3 4.7 3.6 4.1 4.5
##  [721] 5.5 4.1 5.5 5.0 2.9 5.9 3.3 5.0 4.8 4.3 5.9 5.9 5.8 3.2 3.4 3.2 4.0 3.5
##  [739] 4.4 4.4 5.7 3.9 5.4 5.7 4.1 4.4 4.1 4.6 5.1 3.2 5.5 3.3 4.4 3.7 5.3 5.3
##  [757] 3.1 6.2 3.1 4.2 4.3 3.9 4.2 6.0 3.0 5.4 5.3 5.6 5.6 3.0 3.0 5.5 5.4 3.6
##  [775] 2.1 4.7 6.3 5.7 4.2 4.3 5.3 5.5 4.3 2.9 6.0 3.6 5.9 5.9 3.1 4.1 3.7 4.8
##  [793] 5.4 4.3 4.7 3.6 5.5 4.6 4.2 4.3 4.8 4.7 5.2 2.6 3.5 5.2 4.8 4.9 4.6 2.7
##  [811] 3.8 4.3 3.0 7.1 4.1 3.9 5.2 4.6 3.9 3.3 4.9 4.5 4.0 4.8 6.7 5.2 2.7 3.8
##  [829] 3.9 3.4 3.1 3.9 3.8 3.5 5.4 4.7 4.7 5.2 4.0 4.6 4.9 3.2 6.2 3.2 5.5 5.0
##  [847] 4.3 4.6 5.1 3.3 6.3 6.1 5.8 4.7 4.9 3.9 4.0 4.3 4.2 4.1 4.2 4.2 4.1 4.7
##  [865] 4.1 5.1 4.6 4.0 2.6 4.9 3.7 5.9 5.3 3.7 4.0 4.9 6.5 5.5 4.1 4.8 4.9 4.1
##  [883] 4.8 6.0 5.7 4.5 4.2 6.2 4.1 4.5 3.2 4.8 5.0 4.8 1.9 3.9 4.0 3.6 3.0 2.2
##  [901] 5.3 6.0 3.7 4.9 3.5 4.5 6.2 3.1 5.5 4.1 4.3 5.7 3.6 4.7 5.3 3.8 5.1 3.6
##  [919] 5.2 4.5 3.8 4.1 5.2 4.0 5.2 4.5 5.9 3.8 4.8 4.1 5.3 3.4 4.7 4.0 3.3 6.6
##  [937] 3.1 5.2 5.9 4.8 5.1 5.0 4.7 4.2 3.6 4.9 3.3 5.0 5.9 3.8 3.1 5.5 4.5 4.5
##  [955] 4.2 5.2 6.0 5.9 3.6 5.1 2.0 4.6 3.6 4.8 4.9 4.9 3.8 4.2 5.0 3.3 4.3 4.9
##  [973] 3.9 5.7 4.6 5.9 4.1 3.8 4.9 5.5 5.2 4.9 4.9 4.1 3.6 3.0 2.2 3.4 3.5 2.8
##  [991] 4.1 4.0 3.3 5.2 3.7 4.7 4.7 4.3 4.6 4.9
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
##    [1] 6.0 4.3 4.7 4.0 4.8 4.1 4.9 5.6 4.5 3.5 4.9 4.7 4.8 4.1 4.2 3.9 5.5 4.3
##   [19] 4.2 5.0 5.7 3.7 4.6 4.7 4.3 3.3 3.9 6.0 5.2 4.7 4.6 4.3 5.2 3.3 4.7 5.1
##   [37] 4.4 5.9 3.9 4.8 4.7 3.9 4.0 3.6 3.7 4.7 2.5 4.2 5.7 4.2 3.5 5.4 5.1 4.3
##   [55] 3.8 4.9 4.6 3.7 4.0 5.4 4.4 2.3 4.5 3.7 3.1 4.3 3.6 5.4 4.0 4.1 4.5 6.0
##   [73] 3.6 3.8 5.1 3.2 5.4 4.2 4.7 4.4 3.3 5.5 5.3 6.2 5.0 4.7 5.3 6.2 4.4 2.8
##   [91] 2.9 4.8 4.2 4.7 4.4 3.7 4.8 3.8 4.7 5.8 4.2 5.6 5.1 3.1 3.4 5.6 4.3 3.7
##  [109] 4.1 6.1 4.0 4.4 6.1 3.5 4.7 4.4 3.5 3.3 3.7 5.9 4.4 3.2 3.5 3.8 4.3 4.2
##  [127] 5.8 6.2 4.4 3.7 4.9 3.8 3.5 3.4 5.3 4.7 4.4 3.4 4.0 3.9 2.8 3.9 4.0 3.5
##  [145] 3.6 4.0 3.1 2.3 4.3 3.7 5.4 3.1 3.8 3.6 5.6 4.9 5.5 4.4 5.4 5.5 5.3 4.5
##  [163] 4.2 3.6 4.2 3.7 5.3 4.4 4.3 3.2 5.8 4.3 3.7 4.2 5.0 5.4 5.2 3.4 6.0 4.6
##  [181] 4.9 5.0 4.7 4.5 5.3 3.5 5.9 4.2 5.7 5.6 5.9 3.6 4.4 3.8 4.6 5.2 4.9 5.1
##  [199] 4.3 3.8 3.3 5.0 3.0 3.2 3.3 6.8 4.2 4.0 4.7 5.5 4.5 3.7 4.6 4.4 4.1 5.0
##  [217] 3.9 5.6 3.0 3.9 4.4 2.9 3.8 5.4 4.7 4.2 3.8 5.7 3.8 2.9 4.8 4.0 4.1 4.1
##  [235] 5.3 4.7 3.6 4.2 4.5 5.0 3.0 3.2 6.4 3.3 3.8 4.6 5.3 2.8 4.6 4.7 3.0 4.4
##  [253] 2.8 4.8 5.6 4.0 4.0 5.6 3.9 4.5 3.9 3.2 5.7 3.8 4.7 4.6 4.5 2.6 5.1 2.8
##  [271] 5.0 4.0 5.5 3.5 3.8 4.8 6.4 3.1 5.7 4.9 6.4 4.3 5.7 4.7 5.0 3.2 3.4 5.9
##  [289] 4.6 4.4 4.2 4.4 3.9 6.0 6.3 4.2 5.0 4.0 3.7 3.9 4.8 5.1 4.1 3.6 4.1 4.0
##  [307] 3.5 3.9 4.2 3.1 5.0 4.5 4.4 5.0 4.4 4.5 4.0 3.8 4.3 5.8 5.2 2.7 5.2 3.5
##  [325] 5.6 3.7 6.0 4.6 5.2 3.6 4.5 4.0 4.9 4.1 4.5 4.5 4.2 4.6 6.4 3.6 3.4 4.2
##  [343] 5.8 5.9 5.3 3.9 4.8 6.5 2.2 3.2 4.0 4.5 4.4 5.0 5.7 4.0 2.6 3.3 4.2 4.8
##  [361] 4.8 4.9 4.1 4.2 5.7 4.8 4.5 3.8 4.8 5.9 3.4 4.3 5.5 3.4 2.8 3.4 4.2 3.2
##  [379] 2.9 3.0 5.0 5.3 4.3 4.7 3.7 4.1 3.8 3.8 4.6 2.9 5.8 4.7 6.3 6.3 4.8 4.8
##  [397] 5.5 3.6 4.9 2.9 4.1 3.8 3.8 4.7 4.3 4.1 6.3 5.0 3.6 4.7 4.9 3.2 4.5 5.0
##  [415] 5.2 5.3 3.0 6.2 3.3 3.8 5.5 3.3 4.1 5.7 5.3 4.4 4.1 3.1 4.3 4.5 4.5 5.0
##  [433] 3.0 5.4 4.6 5.7 4.9 3.6 3.8 4.9 5.2 4.2 4.3 5.0 4.7 4.3 4.1 3.8 5.3 4.5
##  [451] 6.3 5.3 4.6 4.7 5.7 4.0 3.2 4.2 4.3 3.6 4.7 6.0 4.1 2.7 4.8 3.8 5.2 4.0
##  [469] 5.3 4.5 3.3 5.0 4.1 4.5 3.7 4.5 3.6 3.9 4.5 4.9 4.2 4.4 4.4 4.6 5.2 4.6
##  [487] 5.3 4.6 4.6 3.7 5.3 4.9 4.8 4.7 4.4 3.7 5.0 4.2 3.9 2.5 5.1 4.1 5.1 5.8
##  [505] 4.4 4.0 5.1 4.1 3.1 5.3 4.3 7.0 5.9 4.4 5.0 4.9 5.7 3.3 6.3 6.3 3.5 4.5
##  [523] 5.9 5.1 3.7 4.5 2.8 4.0 3.7 3.1 4.6 5.5 4.2 4.2 3.6 4.7 3.9 4.8 4.6 4.3
##  [541] 3.5 3.3 5.6 5.2 3.7 5.3 5.3 6.1 4.6 5.1 3.8 5.1 3.9 4.9 5.1 4.2 6.7 5.1
##  [559] 5.5 3.5 4.8 4.1 4.5 3.7 3.5 4.6 5.7 4.1 3.8 4.7 4.8 4.7 2.5 3.5 5.1 4.5
##  [577] 4.7 4.8 5.7 4.6 4.8 2.9 3.3 4.6 4.5 3.0 5.0 4.9 3.4 4.6 3.5 4.2 4.9 3.7
##  [595] 5.1 4.3 4.0 3.5 4.8 5.2 6.4 3.6 4.3 4.9 3.7 4.4 4.0 6.1 5.3 5.1 4.6 5.0
##  [613] 6.2 4.8 5.5 4.0 3.5 4.5 3.9 5.2 4.1 4.3 4.0 4.9 4.1 4.0 2.7 6.2 5.2 3.8
##  [631] 4.9 5.3 3.3 4.9 4.3 3.4 5.1 3.2 4.2 5.6 4.1 5.7 2.5 3.3 4.4 4.5 4.2 3.3
##  [649] 5.4 3.0 5.9 5.0 4.5 4.8 4.6 5.9 3.9 2.9 4.3 3.4 3.6 4.8 5.2 5.0 3.6 6.3
##  [667] 6.1 5.2 3.3 5.6 4.4 4.3 5.2 5.5 3.8 3.0 4.2 5.7 3.7 4.1 4.8 3.3 5.6 3.4
##  [685] 4.3 5.9 5.7 4.7 4.8 4.9 5.2 4.2 4.4 5.6 3.6 5.6 4.2 5.0 5.2 3.0 5.6 3.9
##  [703] 4.9 3.4 2.9 4.4 6.1 5.8 4.4 3.4 4.4 3.8 4.3 3.3 4.9 5.4 4.1 3.9 3.9 5.1
##  [721] 3.6 4.6 3.8 5.7 4.7 4.1 4.2 4.8 3.6 4.8 4.3 3.5 5.1 5.3 5.1 5.6 6.3 4.8
##  [739] 5.4 5.1 3.6 4.2 3.7 3.2 4.2 5.4 3.6 4.5 4.5 2.9 4.3 5.4 5.6 3.8 5.3 5.3
##  [757] 4.9 4.5 4.3 5.8 3.7 6.5 4.5 4.3 5.3 6.1 3.5 4.8 4.7 4.9 3.6 6.1 3.8 5.6
##  [775] 2.9 4.4 4.5 5.2 5.1 5.1 4.3 5.4 5.4 2.9 5.6 5.2 6.0 2.8 4.6 4.4 4.5 3.2
##  [793] 5.4 5.6 4.9 3.5 4.8 5.5 5.2 5.1 5.2 5.8 5.1 4.5 4.4 4.2 5.0 4.2 4.2 4.5
##  [811] 4.6 4.1 4.0 5.3 5.1 6.3 5.9 4.6 5.3 5.5 4.7 4.3 5.8 4.2 4.8 5.3 4.8 4.5
##  [829] 2.4 3.2 4.7 5.8 4.5 4.2 2.5 3.5 4.5 4.5 4.5 4.4 4.7 5.1 3.6 6.8 5.2 3.7
##  [847] 4.0 4.2 4.3 5.0 3.5 5.7 4.3 5.8 2.5 4.1 3.3 2.0 4.0 5.2 5.7 5.0 5.0 4.5
##  [865] 4.4 3.7 4.5 5.0 5.3 5.6 4.7 4.9 4.5 3.6 4.6 4.9 4.3 5.8 4.6 3.7 4.7 6.5
##  [883] 3.3 5.4 4.6 4.8 5.3 5.3 5.8 5.4 4.2 1.0 5.5 3.6 4.1 4.2 5.6 3.0 5.0 4.7
##  [901] 3.9 6.0 4.5 4.8 3.8 3.0 3.4 5.0 5.3 3.1 6.3 5.6 3.4 3.2 4.2 4.5 4.2 4.9
##  [919] 4.5 4.2 5.2 4.5 4.8 4.6 5.8 4.5 3.7 5.6 4.2 4.9 3.6 4.5 4.1 4.4 4.4 3.6
##  [937] 4.2 4.6 5.0 3.9 5.1 3.0 6.6 4.7 5.2 3.5 3.3 5.1 5.1 5.7 4.1 5.7 2.8 5.0
##  [955] 4.3 4.4 4.3 4.0 2.5 3.5 5.7 4.8 5.4 4.7 3.4 4.6 6.3 2.6 4.1 5.1 4.6 5.3
##  [973] 5.1 6.3 6.0 5.1 5.2 4.2 6.0 4.8 3.0 4.4 4.3 3.3 5.5 3.9 4.1 5.3 4.0 5.4
##  [991] 2.9 4.3 6.0 4.2 3.9 3.8 4.0 5.0 3.6 5.7
## 
## $func.thetastar
## [1] -0.0128
## 
## $jack.boot.val
##  [1]  0.49708455  0.38526012  0.21761364  0.16439169  0.02945946 -0.10395137
##  [7] -0.18011204 -0.27771261 -0.41049563 -0.49043210
## 
## $jack.boot.se
## [1] 0.9474088
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
##    [1] 4.3 4.6 4.5 3.9 3.0 4.5 5.1 5.5 3.1 4.1 3.6 5.5 4.7 4.5 4.9 5.5 3.1 4.1
##   [19] 6.4 3.7 4.7 3.7 4.8 2.6 3.5 5.0 5.1 3.5 6.6 5.0 5.0 5.1 3.0 4.7 4.9 5.3
##   [37] 5.4 5.3 6.2 4.6 3.6 4.7 6.5 4.0 5.1 5.1 4.8 5.2 4.2 3.9 6.5 4.9 6.2 5.4
##   [55] 3.7 4.8 5.4 5.0 4.7 4.2 5.7 3.1 4.3 3.9 4.0 4.3 5.4 4.1 5.6 4.1 4.8 4.0
##   [73] 4.7 6.0 4.4 4.5 3.3 4.7 3.6 4.6 3.9 5.2 4.9 3.9 6.3 4.6 5.0 4.1 4.6 4.7
##   [91] 3.5 5.9 4.6 3.5 3.8 5.2 3.7 6.2 5.5 3.5 4.9 4.6 3.5 4.1 4.7 3.9 4.4 5.0
##  [109] 3.1 4.9 4.7 4.5 5.8 4.4 2.7 5.4 4.9 4.3 2.9 4.5 5.0 4.2 3.5 5.1 6.1 5.0
##  [127] 6.4 3.8 5.1 5.4 4.1 5.1 5.5 4.4 6.2 3.8 4.2 5.3 5.6 4.8 5.5 3.5 4.7 3.7
##  [145] 4.2 4.3 2.9 5.6 4.7 3.3 3.3 3.7 5.9 4.2 5.8 4.0 4.6 3.5 5.1 4.3 5.7 4.2
##  [163] 6.2 4.5 3.7 3.7 3.3 5.5 4.3 3.6 4.8 2.9 3.9 5.3 4.4 5.4 3.1 3.1 4.7 4.3
##  [181] 4.1 4.4 6.3 4.4 4.0 6.3 5.1 3.3 4.5 4.7 4.2 6.1 5.3 4.8 3.4 4.9 4.6 2.5
##  [199] 4.2 5.6 6.0 3.6 4.4 5.2 3.7 5.6 3.3 6.8 4.2 5.3 4.2 4.9 5.7 5.2 3.4 6.1
##  [217] 5.9 3.9 4.6 6.0 3.9 5.0 4.8 5.8 5.4 3.5 2.8 4.2 5.5 4.4 2.7 4.4 4.0 5.0
##  [235] 5.2 6.2 3.8 4.9 3.7 3.6 4.3 4.3 3.8 5.9 3.8 4.4 4.8 3.4 5.6 4.6 4.3 5.2
##  [253] 5.5 4.2 3.0 4.1 5.0 3.5 3.8 5.3 3.0 4.1 4.4 6.3 4.6 5.3 3.8 5.2 5.8 4.6
##  [271] 6.1 5.7 2.9 4.7 4.1 5.5 4.6 4.1 4.1 4.0 2.6 3.6 4.9 6.9 2.4 4.6 4.9 5.1
##  [289] 5.8 3.7 5.5 4.9 5.0 4.6 3.6 4.7 3.4 3.5 4.9 4.2 2.6 3.5 4.6 2.8 4.1 4.2
##  [307] 5.9 3.5 4.4 4.3 5.1 3.1 3.7 3.1 3.7 4.9 4.9 4.5 5.7 4.6 4.6 4.4 5.9 4.4
##  [325] 2.5 4.3 3.9 5.6 5.6 5.0 5.3 5.7 3.6 3.9 5.3 5.0 5.4 4.2 5.7 4.9 2.6 4.5
##  [343] 5.5 3.1 5.0 4.4 4.4 5.7 4.4 4.0 4.1 4.2 4.7 3.9 6.6 4.0 4.8 3.2 4.6 4.5
##  [361] 6.1 3.9 4.5 4.4 4.7 4.8 4.1 2.1 4.5 3.0 3.7 3.8 5.0 5.2 5.1 4.5 3.3 4.9
##  [379] 4.3 5.4 5.3 4.3 5.5 4.9 4.8 3.5 5.0 3.5 4.8 3.7 4.6 3.3 5.1 3.6 4.4 3.3
##  [397] 4.5 4.7 4.2 3.6 5.2 4.6 4.7 4.8 3.3 3.4 4.7 4.0 6.7 3.1 4.4 5.6 3.9 3.0
##  [415] 4.5 3.9 4.0 4.8 4.3 4.8 5.0 5.1 4.7 3.8 5.2 3.0 5.9 4.3 4.2 4.9 2.9 4.5
##  [433] 4.6 4.4 4.6 3.8 3.2 3.6 4.2 5.3 4.6 5.8 5.3 4.1 7.3 4.4 6.1 5.0 3.7 5.2
##  [451] 4.2 4.0 4.1 5.7 3.6 3.8 3.8 5.4 4.6 3.8 5.6 3.7 4.0 2.9 3.2 5.5 5.6 5.1
##  [469] 4.6 4.5 2.2 4.4 4.9 5.0 5.3 4.1 4.3 3.0 5.1 3.5 5.5 3.4 4.2 3.6 4.5 4.9
##  [487] 5.4 3.6 5.5 4.6 4.5 6.1 5.4 5.6 3.9 5.7 4.5 3.8 5.4 6.0 5.1 5.3 4.4 5.3
##  [505] 4.0 5.8 5.0 4.4 4.5 5.0 3.5 3.7 4.3 4.8 4.2 5.4 4.4 4.4 5.3 3.7 4.9 5.1
##  [523] 2.9 4.9 5.1 4.1 4.8 4.2 5.0 4.4 4.8 3.6 1.7 4.9 5.7 4.7 3.6 5.0 5.0 4.1
##  [541] 5.7 3.2 4.2 4.3 4.2 3.4 3.8 5.1 5.2 3.1 3.6 4.4 2.5 4.0 3.4 3.4 4.5 4.1
##  [559] 4.9 4.0 4.3 6.2 2.5 3.5 2.0 5.4 4.0 4.7 5.7 4.8 4.3 3.3 5.3 3.9 2.8 3.6
##  [577] 4.7 4.0 4.5 4.6 5.7 4.9 4.3 4.0 4.2 3.8 5.2 5.2 2.6 4.0 2.6 4.8 5.1 2.4
##  [595] 4.8 2.8 5.9 4.2 6.0 4.2 4.8 4.7 5.9 4.2 5.7 3.2 5.3 5.8 2.8 5.1 4.5 3.4
##  [613] 4.2 4.5 4.4 4.7 4.1 5.2 4.5 3.7 4.1 3.9 4.5 5.0 5.3 3.5 5.6 5.5 5.9 3.9
##  [631] 3.9 5.1 5.3 3.4 3.6 5.7 3.3 5.3 3.3 4.4 5.8 4.8 5.7 5.6 5.5 4.8 3.6 4.3
##  [649] 4.7 3.7 4.6 4.0 3.4 3.5 5.1 5.0 3.9 3.4 3.8 3.6 4.2 4.2 4.8 5.2 4.5 5.0
##  [667] 4.3 4.7 5.0 3.1 4.2 4.6 3.7 2.9 6.4 3.8 4.3 4.3 6.3 5.4 4.2 4.9 5.3 4.0
##  [685] 3.0 5.0 3.5 3.9 4.8 5.0 2.7 4.9 5.7 4.1 4.9 5.3 4.3 5.1 4.9 3.3 4.5 4.5
##  [703] 3.9 3.9 4.6 5.6 4.2 2.7 3.9 5.8 4.1 5.0 4.9 5.0 5.8 4.7 4.1 3.1 5.0 4.8
##  [721] 4.0 4.8 4.3 4.1 4.4 6.1 4.5 4.1 3.4 3.1 4.3 4.5 3.9 4.3 6.4 3.7 4.7 2.7
##  [739] 4.7 5.6 5.8 6.1 3.8 4.8 3.7 4.4 5.3 5.5 5.0 3.9 3.9 6.0 4.7 6.2 3.5 3.7
##  [757] 3.9 4.6 5.6 4.0 3.5 4.7 4.7 5.8 4.5 4.8 4.9 4.0 4.1 4.6 4.5 5.3 6.1 6.5
##  [775] 6.1 4.6 4.3 4.3 3.5 4.6 4.7 5.3 4.0 5.4 4.6 5.0 4.8 5.5 3.8 4.3 4.8 4.1
##  [793] 4.5 2.7 3.6 5.3 5.6 4.0 4.5 4.2 3.3 3.4 4.4 4.3 5.8 3.7 5.3 4.9 5.8 4.9
##  [811] 5.7 4.0 5.5 6.4 4.8 3.5 5.3 5.0 4.2 5.2 4.1 6.0 4.8 5.1 5.1 3.9 4.6 5.2
##  [829] 5.3 6.0 4.1 4.7 4.7 4.2 4.6 4.5 4.3 4.2 3.9 3.7 5.0 6.1 2.9 5.0 5.3 2.6
##  [847] 4.1 4.7 3.1 4.9 5.5 4.8 4.8 3.7 4.7 6.2 4.2 4.6 5.6 3.9 3.5 4.6 5.2 4.1
##  [865] 4.0 4.5 4.9 6.0 5.1 4.5 4.6 5.6 4.8 4.6 4.5 3.6 4.7 3.0 3.5 4.9 5.1 5.7
##  [883] 4.5 4.1 4.1 3.5 3.9 3.6 2.9 4.2 4.7 6.0 5.6 3.2 3.4 4.2 3.0 3.3 4.2 1.8
##  [901] 4.6 4.7 4.5 2.9 5.5 3.7 4.2 3.8 5.6 4.5 3.6 5.2 4.4 4.4 5.1 4.5 5.7 3.9
##  [919] 5.0 2.3 4.5 5.1 4.1 4.1 4.5 3.1 5.9 4.5 3.4 2.7 3.7 3.8 4.2 5.6 6.1 4.7
##  [937] 3.6 2.7 4.2 3.8 6.7 5.7 4.5 3.3 4.4 4.3 5.6 5.1 5.4 4.7 5.7 2.1 3.3 3.2
##  [955] 4.4 2.8 3.3 5.1 4.2 6.9 3.0 4.5 4.8 4.5 3.7 4.2 4.4 3.9 5.3 3.3 4.6 4.6
##  [973] 3.9 3.6 4.7 3.7 5.4 5.1 6.3 4.7 5.0 4.6 3.8 4.2 3.4 5.1 4.8 3.4 4.7 5.1
##  [991] 4.6 3.7 5.6 2.7 4.1 5.1 3.0 6.2 4.7 4.8
## 
## $func.thetastar
## 72% 
##   5 
## 
## $jack.boot.val
##  [1] 5.5 5.4 5.2 5.1 5.2 5.0 4.7 4.8 4.6 4.6
## 
## $jack.boot.se
## [1] 0.9241753
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
## [1] 0.2532979
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
##   5.053832   9.596814 
##  (2.189618) (4.371560)
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
## [1]  0.60044066  0.01878966  1.49506937 -0.11784357  0.74036381  1.52959959
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
##    [1] -0.207065976  0.381483873  1.204655196 -0.321478317  0.319148279
##    [6] -0.213638047  0.303193721  0.370697037 -0.683371344  0.279199360
##   [11]  0.374368462 -0.372062033  0.583117940  0.248494640  0.848404531
##   [16]  0.151848403  0.499474758  0.272607730 -0.531431353  0.076194480
##   [21]  0.115405867  0.083104506  0.006824156 -0.306949702  0.316205457
##   [26] -0.260844298  0.410949745  0.707094510  0.507862112  0.684162814
##   [31] -0.573962589  0.307829208  0.930290245 -0.351392305  0.065781862
##   [36]  0.537858012  0.698396021  1.368134945 -0.681285332  0.497314612
##   [41]  1.254492831  0.531727657  0.880107221 -0.441073449  1.533163470
##   [46]  0.035289628  0.032878868  0.613690597  0.596190622 -0.040172762
##   [51]  0.503033112  0.031609665 -0.011433873 -0.022934602 -1.169107498
##   [56]  0.863876135  0.203239608  0.671630359 -0.127229090  0.134744241
##   [61] -0.430740486 -0.548579724  0.140421818 -0.344328262  0.172209056
##   [66] -0.408121557 -0.124248227  1.910250431  0.184889250 -0.275321019
##   [71]  0.459090242 -0.593736620  0.581899773  0.283154713  0.605058642
##   [76]  0.862110147 -0.536453821  0.909032001  0.219360181  0.266026719
##   [81] -0.279917696  0.547526142  0.730996559 -0.227728785 -0.069912521
##   [86] -0.159276013  0.527793714  0.993326936  0.227820098  0.619725258
##   [91] -0.161309340  0.725331084  0.059009945  0.039186757  0.529163956
##   [96]  0.607961691  0.158181377  0.176304082  0.778175423  0.426175644
##  [101]  0.192872175  0.686515948  0.456837432  0.830100268  0.259333393
##  [106]  0.597257571  0.105128557  0.458748142 -0.303905173 -0.784570878
##  [111]  0.387992325 -0.686568099  0.183693612  0.774031972  0.502305860
##  [116]  0.233107362 -0.363742848  0.599532127 -0.235201956 -0.423593525
##  [121]  0.921921691  0.059956643  0.097697423  0.136810840  0.934854720
##  [126]  0.197469639  1.603228334 -0.484354077  0.124379090 -0.111485784
##  [131]  0.846957888 -0.217531263  0.464979177  0.082531652  0.385691562
##  [136]  0.743413940 -0.569033580  0.242709874  0.895502367  0.736087598
##  [141]  0.376533413  1.026667913  0.185969494  0.381960762  0.624916903
##  [146]  0.169866999 -0.433257501 -0.171553782  0.275615206 -0.129085527
##  [151]  0.155779993 -0.300526232  0.875943107  1.148588605 -0.298090779
##  [156] -0.156592122 -0.309274140  0.168937235  0.105414743  1.024153548
##  [161] -0.351538638  0.173384280  1.149446569  0.620708300 -0.604409775
##  [166]  0.977223023 -0.188312961 -0.419265985  0.503259099  0.523826029
##  [171] -0.102808495  0.730402435  0.347642996 -0.461375457  0.412651326
##  [176]  0.958114027  0.649050003  0.531102894  0.775292494  0.026230049
##  [181]  0.177455316  0.345879716  0.156630944  0.114448958  0.255092042
##  [186]  0.444974817  0.385351878  0.929435093 -0.427016056  0.265188755
##  [191]  0.074250770  0.359719393  0.762544077  0.901869857  0.320970303
##  [196] -0.683548226 -0.104127151  0.589008515 -0.205836638  0.064658844
##  [201] -0.157466809  0.023580026 -0.109903346  0.015731756  0.448365494
##  [206]  0.341656499 -0.217400811 -0.067597237 -0.854895238  0.475880627
##  [211]  0.865824792 -0.070103933  0.082352060  0.666836749  0.883690840
##  [216]  0.377538632  0.894659149  0.752553633  1.118978255  0.535223836
##  [221]  0.421864957  1.106478884  0.613604151  0.601382738  0.804139881
##  [226]  0.767440551  0.409137726 -0.103663358 -0.935014771  0.394354912
##  [231] -0.134064407  0.573379694  0.871744422 -0.005794226 -0.078116207
##  [236] -0.352777089 -0.259006860 -0.094684273  0.393752092 -0.025969507
##  [241]  0.294106076  0.103340082 -0.770791607  0.505862122  0.071757310
##  [246]  0.318585631  0.148740204  0.273476566  0.746876540 -0.289204633
##  [251]  0.124387859  0.698922650  0.387707954  0.124038661 -0.343690697
##  [256] -0.446461123  0.915283156  0.149063458  0.293314037  0.126008918
##  [261]  0.393857961 -0.448309298  0.015021019  0.358041442 -0.533686964
##  [266]  0.282751943  0.347173300  0.446227845 -0.155568347  0.253127316
##  [271]  0.481010470  1.489437031  0.944080967 -0.055311026  0.258573591
##  [276] -0.159576232 -0.088195197 -0.037940189  0.278713982  0.511669302
##  [281]  1.204426814  0.488333002 -1.095740637  0.481025021  0.695464196
##  [286]  0.845679651  0.792158039  0.486905403 -0.116514217 -0.251330402
##  [291]  0.247544643  0.330242156  0.026690478  0.252017646  0.320729824
##  [296]  1.347473617 -0.105933105  0.105905978  1.365562537 -0.298391405
##  [301]  0.791910551 -0.029136113  0.442211957 -0.338772787  0.232726736
##  [306] -0.484640165  1.203316463  0.124955391  0.985133175  0.141668492
##  [311]  1.470606089  0.251186305 -0.375164265  0.362828094  0.088576420
##  [316]  0.538671368  0.475635389  1.363805501  0.403890428  0.494425337
##  [321] -0.080394886  0.823765567  0.467780409 -0.083197111 -0.392018542
##  [326]  0.649595577  0.169646476  0.619092389 -0.193810660  0.550991577
##  [331] -0.769656705  0.866483042 -0.408795745  0.602112390  0.605430911
##  [336] -0.320522369 -0.290115900  0.517989590 -0.035759175  0.056952416
##  [341]  0.767638829  0.228631006  0.461197250  0.028569180  0.184004053
##  [346]  0.133535759  0.642901901  0.484743694 -0.432361246  0.016115882
##  [351] -0.444457209  0.727219890  0.377713438 -0.304731558  1.132984854
##  [356]  0.823489694  0.384052002 -0.010940569 -0.539234310  0.266071857
##  [361]  0.295528883  0.284569545  0.567309125 -0.479914478  0.251770625
##  [366]  0.186736346 -0.273243380  0.683481956 -0.847220626  0.291544193
##  [371]  1.145944130  0.090426640  0.447400237 -0.115618584 -0.208986230
##  [376]  0.269148563  0.444360402  0.650576714  0.118036433  0.180306925
##  [381] -1.338822663  0.008190054  0.233609698 -0.050680392  0.553391663
##  [386]  0.057578511  0.670089439  0.613390625 -0.115284752  0.957595398
##  [391]  0.167279041 -0.256289731  0.588331196  0.187633190  0.943738628
##  [396]  0.025146916  0.832111426 -0.363212950  0.127739585  0.153865630
##  [401]  0.832091685  0.221347080 -0.085849453  0.299554306 -0.128261100
##  [406]  0.566612723  0.376533413  0.904728768  0.587200045  0.977095764
##  [411]  0.924873585 -0.100149800 -0.033389239 -0.452210260  0.219416280
##  [416]  0.565198189  0.462881919  1.125839738 -0.112129317  0.031411363
##  [421]  0.360817650  0.497191589 -0.133971559  0.278283158 -0.358377704
##  [426] -0.373867097  0.235220730  0.605430911  0.675088710  0.148418460
##  [431]  1.118352080 -0.079214237 -0.029102837  0.265120809  0.669117822
##  [436]  0.509564274  0.526129923  0.014808289  0.820011486 -0.024463257
##  [441]  0.670263801  0.078798355 -0.083310396 -0.120196523  1.703071825
##  [446]  0.534346315 -0.958511050  0.603029212 -0.638575513  0.088929744
##  [451]  1.286630060  0.598681053 -0.226191276  0.976228560  0.384610329
##  [456]  0.333747585 -0.183746571  0.808816954 -0.311285447  0.589478543
##  [461] -0.328047788  0.601084439  0.806786800 -0.013177327  1.612208264
##  [466] -0.193810660 -0.039176622 -0.661238974  0.716578810 -0.137387257
##  [471] -0.341224782 -0.407829781  1.324003537  0.691613617  0.112454540
##  [476]  0.113544983 -0.173798397 -0.213257144  0.978403206  0.988143740
##  [481]  0.787338597 -0.434491410  1.162229719  0.289182377  0.885194622
##  [486]  0.155549787  0.577587925  1.147528123  0.689716407  0.367401820
##  [491] -0.106116707  0.935206705  0.555404394  0.536606237 -0.152556073
##  [496] -0.494293457 -0.120823566  0.531173614  0.227738283  0.366537649
##  [501] -0.094448936  0.684481123 -0.019419813  0.494955990  0.632904788
##  [506] -0.022365124 -0.076490693  0.798007738 -0.301727299 -0.096392453
##  [511] -0.528779175 -0.569787188  0.285871251 -0.018908500  0.674753047
##  [516]  0.695201300  0.834407114  0.522470071  0.693028082 -0.007720622
##  [521]  0.422812643 -0.195301040  0.322134185  0.815948622  0.508458425
##  [526] -0.083215344  0.930527685 -0.317804431  0.169066144 -0.369266275
##  [531] -0.548791117  1.819395551 -0.263748417  0.785917885  0.639519634
##  [536] -0.249171145  0.549727049 -0.149734867  0.617194568  0.231578132
##  [541]  0.806786800  1.278203049  0.440529338 -1.095767357  0.186954073
##  [546]  1.433168694 -0.031231556  0.344346322 -0.202462792  0.644364514
##  [551]  0.016526342 -0.750736041  0.101684364  0.436165083  0.982641301
##  [556]  1.170529427 -0.197542027  0.325278873  0.425803262  1.351908910
##  [561]  0.666233934  0.200061584  0.783118127  0.225398569  1.148341474
##  [566]  0.444365469 -0.409438426 -0.005252835 -0.046568816  0.364904032
##  [571]  0.005031354 -0.381366262 -0.176840520  0.745169183 -1.148492585
##  [576]  0.523917713  0.644297327  0.132119878  0.293805832  1.376322785
##  [581]  0.110625372  0.831438848  0.119429516 -0.083596771 -0.014827703
##  [586]  0.351165356  0.233486090  0.260065978  0.151309372  0.715958455
##  [591]  0.164650407  0.113163229 -0.088370249 -0.003475054  0.252811370
##  [596]  1.186947233  0.109757964 -0.070942404  0.798805151 -0.187564373
##  [601] -0.635680292  0.833500660  0.632498918  0.240952787  0.224009543
##  [606]  0.235754739  0.327155569 -0.204837110 -0.543615029  0.607831914
##  [611]  0.559503022  0.267723214  0.320569303  1.038706351 -0.284404832
##  [616]  0.748709337  0.272818455  0.695885237 -0.342026450  0.654761359
##  [621]  0.366664071  1.475804656 -0.359063574 -0.110677869  0.193612114
##  [626]  0.755315699 -0.083477574  0.425350299  0.258749367  0.534106434
##  [631]  0.184456040  0.155208930  0.231449564  0.056318500 -0.097413501
##  [636]  0.948912439  0.068831245  0.024701277  0.658058863  0.907348309
##  [641]  0.131722297  0.371337705  0.298367699  0.385383449  0.428332539
##  [646]  0.187089437  0.729047775  1.119908237 -0.338601852  0.557234126
##  [651]  0.120458366  0.592093887  0.348175644  0.382948348  0.995216187
##  [656]  0.403576661  0.577607430  0.651888029  1.037765272  0.733276575
##  [661] -0.267726714  0.334256788 -0.829846298  0.645232863 -0.017468734
##  [666]  1.480838288  0.124722536 -0.177027006  0.572503679 -0.212603842
##  [671]  0.637358512 -0.932825158  0.017529809  0.180819794  0.325935693
##  [676] -0.442480505  0.720855675  0.357339913 -0.173926862  0.129259223
##  [681] -1.001746727  0.163540409  0.179365182 -0.157895506  0.257193045
##  [686]  0.577207862  0.731822634 -0.710535755  0.272718518  1.195273658
##  [691]  0.099710248 -0.335787993 -0.315951498  1.560839142  0.278641234
##  [696]  0.548437731  0.157556103  0.210758440  0.189337303  0.445337226
##  [701] -0.744715373  0.033068557  1.456024880  0.121768088 -0.462071189
##  [706]  0.898460926 -0.104008522 -0.095457119 -0.131372310  0.545020018
##  [711]  0.682140307  0.376497924  0.497994904  0.119391812  0.132208482
##  [716]  0.601084439  1.074658811  0.696256950  0.241790463  0.568622651
##  [721]  0.903981963  0.858311170  1.249846137  0.564202623  0.570648628
##  [726]  0.256591999  0.655269310  1.224183445  0.379242401  0.104601433
##  [731]  0.146602398  0.669738734  0.836303658  0.542774952 -1.180072682
##  [736]  0.193831763 -0.016811056 -0.064253102 -0.651444515 -0.211052922
##  [741]  0.285290828  0.703704690  0.150351860  0.023102518 -0.395645739
##  [746] -0.596802052 -0.530586299  0.446628716  1.017359249  1.275564636
##  [751] -0.113613866  0.336104116  2.195412971  0.392640764 -0.088886168
##  [756]  0.676933278 -0.167043149  0.645759111  0.047460630  0.267936725
##  [761] -0.320250536  0.479846407 -0.215605412  0.522757110 -1.141157247
##  [766] -0.193810660  0.582045941  0.416123030  0.473482257 -0.079099899
##  [771]  0.061742159  0.607961691  0.004368738  0.705752541  0.949730598
##  [776]  1.129007654 -0.045532535 -0.100172294 -0.151950673  0.193076371
##  [781]  0.812218496  0.483463874 -0.330892614  0.418616334  0.837886226
##  [786] -0.227977561  1.260791675  0.136466377  0.445829288  0.006013577
##  [791]  0.783588804  0.621729956  0.183209597 -0.071131177 -0.107493199
##  [796]  1.750197259  0.405787661  1.386971253  0.395250202  0.047253808
##  [801]  0.710249764 -0.816606619 -0.259523119  0.160876001  0.389117607
##  [806] -0.192667616 -0.224045903  0.686428683 -0.246013997  0.519043867
##  [811] -0.753610616  0.367163246  0.210670656  0.219407531 -0.219248836
##  [816] -0.359063574  0.197677377 -0.068138373  0.493498529  0.241917661
##  [821]  0.705710953  0.716576553 -0.043949187 -0.923470752  0.866146387
##  [826]  0.416899205  0.396932028 -0.324946374 -0.082871067  0.609326596
##  [831]  0.235220131  0.635797001 -0.213449289  0.720582429 -0.318202521
##  [836]  0.491027162  0.388044295  0.170075014  0.083019044  0.539489990
##  [841]  0.397924048  0.630164462  0.476747908  0.403890428  0.542494307
##  [846]  0.660017168 -0.023117002  0.409734641 -0.218266817  0.354596000
##  [851] -0.405758032 -0.235781150  0.866379273  0.378644728  0.037406176
##  [856]  0.528785030 -0.288590098  0.921609539  1.350312196  0.358429023
##  [861]  0.095300042  0.115356641 -0.217621202  0.161653588  0.168937235
##  [866]  0.540687383 -0.036338017 -0.338549395  0.669432478 -0.449447375
##  [871]  0.366185451  0.192079024  0.216982558 -0.060142486  1.211542247
##  [876]  0.106804371  0.177637981 -0.261570612  0.519578534 -0.021056650
##  [881] -0.259445142  0.056016495  0.510455624 -0.429799651  0.511515480
##  [886]  0.715682508  0.262262222  0.441241314  0.272818455  0.708358477
##  [891]  0.218052065  0.151649143  0.342128409 -0.369814324 -0.021375636
##  [896]  0.135649522  0.440166364  0.340273104  0.592093887  0.932158921
##  [901]  0.219666986 -0.054147202 -0.133641515  0.146597317 -0.252965030
##  [906]  0.421898735 -0.125872151  0.435125914 -0.339143677  0.235907040
##  [911]  0.457036345 -0.076855827 -0.133761147  0.096357949  0.031760928
##  [916]  0.316405539  0.680249442  0.460658627  0.513549179  0.228207498
##  [921]  0.549892200  0.239552524  1.122331300  1.577628675  0.145139949
##  [926]  0.179061550  0.190420460 -0.266775115  0.497220152  0.340003589
##  [931] -0.224291172 -0.605538431 -0.198409790  0.522217010  0.041075412
##  [936]  0.608705181  0.023675733 -0.009337507  0.408133075  0.958881348
##  [941] -0.607136103 -0.872754699  0.234244142  0.217682811  0.708431891
##  [946]  0.818035644  0.398391340 -0.723312638 -0.069912521 -0.519902788
##  [951] -0.170809309  0.036158134  0.008002082  0.808346240  1.486758435
##  [956]  0.784545080 -0.186910120  0.621729956  0.997425076  0.352567358
##  [961]  0.981000214  1.015690038  0.685315807  0.451974360  0.382086919
##  [966]  0.525459916  0.671647195 -0.314002152 -0.538213563  0.390475439
##  [971] -0.042500652  0.320245064  0.850451334  0.175474674 -0.005793090
##  [976]  1.019603940  1.105446842  0.455233229 -0.684618361  0.024556174
##  [981]  0.475376159  1.075444020  0.440098763  0.033736113  0.239656724
##  [986]  0.042439543  0.266189865  0.490382845 -0.449086682 -0.001610243
##  [991]  0.355862331  0.573593321 -0.222221093  1.103697979 -0.305560067
##  [996]  0.634146686  0.140455876  0.209672529 -0.191269619  0.522413467
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
##   0.52661249   0.22147274 
##  (0.07003583) (0.04951918)
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
## [1] -0.478762985 -0.026550594  0.009221888  0.032572747  0.051291803
## [6] -0.241285712
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
## [1] -0.0083
```

```r
sd(results2$thetastar)  #the standard deviation of the theta stars is the SE of the statistic (in this case, the mean)
```

```
## [1] 0.9021058
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
## t1*      4.5 0.05935936   0.9168922
```

One of the main advantages to the 'boot' package over the 'bootstrap' package is the nicer formatting of the output.

Going back to our original code, lets see how we could reproduce all of these numbers:


```r
table(sample(x,size=length(x),replace=T))
```

```
## 
## 0 2 3 4 5 6 8 
## 2 1 2 1 2 1 1
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
## [1] -0.0031
```

```r
se.boot
```

```
## [1] 0.9333872
```

Why do our numbers not agree exactly with those of the boot package? This is because our estimates of bias and standard error are just estimates, and they carry with them their own uncertainties. That is one of the reasons we might bother doing jackknife-after-bootstrap.

The 'boot' package has a LOT of functionality. If we have time, we will come back to some of these more complex functions later in the semester as we cover topics like regression and glm.

