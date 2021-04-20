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
## 0 1 2 3 4 5 6 8 9 
## 1 1 1 1 1 1 1 1 2
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
## [1] -0.0212
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
## [1] 2.631445
```

```r
UL.boot
```

```
## [1] 6.326155
```

Method #2: Simply take the quantiles of the bootstrap statistics


```r
quantile(xmeans,c(0.025,0.975))
```

```
##  2.5% 97.5% 
##   2.8   6.4
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
##    [1] 3.8 5.9 6.0 4.9 5.3 5.6 4.0 3.8 4.1 4.1 3.8 2.9 3.7 5.0 4.4 5.4 3.2 4.0
##   [19] 5.6 5.7 3.2 4.5 4.5 4.0 4.6 2.6 5.1 3.6 5.5 5.8 4.2 4.8 4.9 5.2 5.5 5.3
##   [37] 3.1 3.8 4.0 5.1 4.4 3.0 4.3 4.7 3.4 5.8 3.6 5.1 4.7 4.6 4.4 4.8 4.3 3.3
##   [55] 5.8 4.1 6.5 3.3 4.2 4.5 4.9 3.9 5.9 5.2 3.5 3.7 4.1 4.9 4.0 4.7 5.0 2.9
##   [73] 4.1 4.5 4.2 4.2 4.3 5.3 5.6 4.1 3.6 3.9 5.3 4.2 3.7 3.4 3.8 3.8 5.2 4.2
##   [91] 3.7 4.6 3.6 4.9 4.6 3.4 3.6 5.1 4.9 4.3 4.9 5.2 4.2 4.5 4.2 3.7 4.9 5.4
##  [109] 4.9 5.4 5.3 5.1 3.9 5.1 5.6 4.0 4.9 4.2 3.8 5.3 3.2 3.5 5.4 5.7 5.3 2.6
##  [127] 3.4 6.4 4.2 3.9 3.7 4.2 5.7 4.6 4.9 4.6 4.0 3.9 4.5 4.2 4.7 4.8 5.2 4.0
##  [145] 4.6 6.6 3.4 4.3 4.3 4.5 3.5 4.6 3.7 4.0 4.6 4.0 4.1 5.3 4.9 3.6 4.5 4.7
##  [163] 4.6 4.4 4.8 3.7 4.1 3.5 3.4 3.5 4.5 2.9 3.5 5.5 6.1 4.5 5.3 3.3 5.8 4.6
##  [181] 6.2 4.9 2.9 5.2 3.9 4.7 3.2 3.9 4.8 4.2 5.0 4.5 4.2 4.3 3.8 4.4 3.5 3.3
##  [199] 4.6 4.4 5.3 4.6 3.8 4.6 4.9 3.7 5.4 4.3 4.6 5.3 3.8 5.8 4.3 4.7 4.5 4.2
##  [217] 5.0 4.9 3.5 5.6 3.9 4.7 5.2 4.1 3.0 5.5 4.2 4.6 3.9 4.1 4.6 4.5 5.6 3.5
##  [235] 3.2 4.8 4.6 4.4 4.4 4.6 3.9 5.0 2.8 3.9 3.7 4.0 3.4 2.2 4.1 4.5 4.5 4.7
##  [253] 5.9 3.5 4.6 5.2 3.1 4.3 4.1 4.6 3.8 3.5 3.1 4.2 4.8 5.0 5.7 5.1 4.4 3.1
##  [271] 4.4 4.9 5.5 5.3 4.8 4.2 4.6 4.0 3.4 5.2 4.8 5.1 4.6 3.0 4.6 4.7 5.5 4.0
##  [289] 3.8 5.2 3.1 3.5 4.5 5.3 5.6 3.3 3.8 4.4 4.3 3.7 1.9 5.1 5.4 5.6 3.5 2.7
##  [307] 4.7 4.2 6.0 4.1 5.3 4.1 4.5 4.2 5.3 4.5 3.7 4.6 5.4 4.7 3.4 4.0 5.5 3.9
##  [325] 5.8 4.2 4.9 3.3 5.1 3.5 5.1 3.5 4.0 4.1 5.1 4.6 4.5 2.2 3.9 4.8 5.8 4.7
##  [343] 4.7 3.8 4.2 4.6 5.1 5.1 3.2 5.5 5.5 4.7 3.4 5.6 5.1 4.7 3.6 3.7 3.1 3.9
##  [361] 5.3 2.4 6.6 4.9 4.0 5.4 3.9 4.9 6.3 2.8 4.0 4.2 4.9 3.7 2.5 4.9 4.4 4.7
##  [379] 4.0 4.3 4.5 5.8 4.4 3.4 2.9 5.3 5.1 4.5 4.7 6.9 3.0 3.6 5.2 5.1 5.4 5.0
##  [397] 5.6 5.2 4.8 4.5 5.4 5.0 5.5 4.8 2.3 3.1 3.9 5.1 4.5 5.7 4.3 4.6 3.3 3.8
##  [415] 5.4 4.5 4.3 2.8 4.5 4.6 3.2 2.6 3.7 4.5 6.2 3.6 3.3 3.8 3.6 4.4 4.2 3.4
##  [433] 5.0 5.5 6.5 3.8 4.9 5.2 3.5 5.8 4.6 5.0 5.5 4.6 4.5 3.9 4.1 6.3 3.5 4.5
##  [451] 4.9 5.0 4.0 2.9 4.2 5.5 5.0 3.9 5.5 3.8 4.3 3.3 6.6 4.3 5.0 4.2 3.1 5.0
##  [469] 3.6 3.9 4.1 3.9 2.9 2.6 4.6 4.6 3.7 5.5 5.0 5.0 5.7 4.7 5.2 5.7 5.1 6.9
##  [487] 3.6 3.7 3.2 4.7 3.5 3.9 5.0 5.8 6.0 5.5 3.8 5.2 4.3 5.7 4.1 5.8 3.9 4.2
##  [505] 5.4 3.5 4.7 3.5 5.0 3.4 3.1 5.4 4.3 4.0 4.3 6.6 5.2 2.4 5.3 4.4 5.2 4.3
##  [523] 3.6 2.8 5.5 6.1 2.9 4.5 3.7 5.4 3.7 4.8 6.6 3.9 4.8 3.6 4.1 5.0 5.0 3.4
##  [541] 3.6 4.2 5.1 4.8 4.0 3.2 5.6 6.5 3.6 3.9 4.7 6.3 6.2 5.1 4.0 3.5 4.3 5.8
##  [559] 6.0 5.0 4.4 4.9 3.6 4.4 4.6 4.7 5.9 4.2 5.4 5.4 4.3 3.5 2.9 4.6 5.7 5.1
##  [577] 3.5 4.2 2.6 3.2 4.4 3.8 4.4 5.3 4.1 2.8 3.5 2.6 4.1 4.2 6.1 4.6 3.9 4.2
##  [595] 5.5 4.5 3.5 5.4 6.5 5.5 5.4 4.2 4.3 4.4 4.4 6.2 5.0 3.8 5.6 3.5 6.1 3.7
##  [613] 5.2 4.7 4.7 5.6 3.6 6.6 4.9 4.1 4.4 5.0 5.2 5.6 4.3 5.0 3.4 6.2 3.9 2.8
##  [631] 4.1 4.9 4.0 5.9 3.7 4.7 3.8 5.0 3.8 3.6 3.9 3.9 4.7 5.3 4.2 3.3 4.5 3.7
##  [649] 5.9 4.9 4.6 3.5 4.9 3.6 4.4 3.6 5.1 3.6 4.8 6.3 2.9 3.8 4.1 5.2 5.3 4.5
##  [667] 5.8 5.4 5.4 3.3 4.3 5.1 4.2 3.0 4.9 4.9 4.8 6.7 4.8 3.9 4.7 3.8 4.7 3.3
##  [685] 5.1 4.9 3.6 4.8 3.0 5.6 4.5 3.6 6.2 4.2 6.0 3.4 3.0 5.0 4.8 3.7 3.9 5.6
##  [703] 4.5 3.2 5.5 3.5 5.4 4.7 5.2 4.7 3.3 2.8 3.7 4.2 5.7 2.7 4.3 4.5 5.9 4.3
##  [721] 3.8 5.0 4.9 3.8 5.2 4.5 3.4 5.0 3.9 4.2 3.5 4.6 3.3 4.9 4.8 5.5 3.7 4.5
##  [739] 4.2 3.1 5.6 3.3 4.5 5.5 3.7 4.3 5.8 4.0 5.2 3.5 3.7 4.9 4.3 5.4 4.3 2.9
##  [757] 3.6 4.5 4.3 2.8 3.6 4.1 5.3 4.3 3.7 4.1 4.3 4.4 2.5 5.0 3.9 5.1 4.1 4.1
##  [775] 6.2 4.4 4.7 4.9 4.8 3.5 4.8 5.2 3.9 3.2 3.7 4.0 4.2 5.6 4.5 4.5 5.6 4.0
##  [793] 5.3 5.3 4.6 2.7 5.0 4.5 4.9 5.6 5.1 4.1 4.4 3.3 5.2 6.1 4.1 3.6 3.5 4.6
##  [811] 4.6 2.6 4.8 3.7 5.6 5.1 4.5 4.8 4.6 4.9 4.4 3.5 3.6 2.9 5.4 4.1 3.2 3.9
##  [829] 4.6 4.1 3.3 2.9 4.2 3.6 4.2 4.1 3.1 5.6 4.3 5.9 4.2 4.5 4.6 6.0 3.6 5.0
##  [847] 4.2 6.6 5.3 2.7 3.9 4.6 3.7 3.9 6.1 5.4 4.6 3.5 3.3 4.5 5.7 4.4 3.1 4.7
##  [865] 5.5 4.8 4.0 4.9 4.8 3.4 3.9 3.6 5.1 4.4 3.4 3.4 5.9 4.7 4.0 4.3 3.6 3.6
##  [883] 3.4 3.4 4.7 3.5 6.3 4.3 4.0 4.0 3.6 5.2 4.6 3.6 3.6 5.2 4.3 3.2 3.3 1.3
##  [901] 4.4 3.1 5.8 4.8 3.1 5.1 3.4 6.1 3.9 3.5 4.1 5.9 4.0 4.1 4.5 3.7 5.3 4.6
##  [919] 6.5 5.2 3.2 4.4 3.9 4.7 3.4 4.3 4.7 3.0 4.3 5.2 4.6 4.3 4.5 5.2 5.0 4.0
##  [937] 3.3 4.2 4.4 5.2 4.8 4.7 5.2 4.8 5.2 4.3 5.2 4.1 4.5 4.5 4.0 4.0 3.6 4.0
##  [955] 4.7 4.9 7.1 5.5 3.9 3.2 3.8 4.2 4.4 4.6 3.6 3.1 3.3 5.5 4.2 3.6 4.7 6.7
##  [973] 3.4 3.8 5.1 2.3 3.3 3.5 5.0 5.2 5.1 4.0 3.8 5.5 4.4 5.2 2.9 5.0 3.2 3.8
##  [991] 5.3 5.5 4.7 5.8 4.7 5.5 4.4 4.9 4.2 5.1
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
##    [1] 4.1 4.3 3.4 6.0 3.0 5.1 5.2 4.4 3.4 4.1 5.1 3.6 4.6 5.6 5.0 5.0 4.5 6.1
##   [19] 3.7 5.0 5.3 5.7 4.5 4.1 5.3 3.8 5.5 3.5 3.8 5.2 3.5 4.2 4.9 4.9 4.0 4.3
##   [37] 4.4 3.8 4.5 4.1 2.4 5.1 4.5 5.3 5.3 4.6 2.5 5.3 4.8 4.0 3.7 3.9 4.7 3.3
##   [55] 4.7 6.0 5.1 4.0 4.8 4.1 3.1 3.8 2.5 4.6 4.6 3.7 4.8 6.1 3.6 2.8 2.8 4.2
##   [73] 5.4 4.5 5.2 5.1 3.0 4.1 4.0 4.0 4.4 6.5 4.2 4.1 5.3 4.2 4.3 4.3 5.7 3.2
##   [91] 3.1 5.7 3.9 3.9 6.0 5.1 5.3 5.8 5.6 5.3 5.5 4.6 5.6 4.6 3.6 3.8 5.2 5.0
##  [109] 3.0 5.5 5.1 3.5 3.8 5.1 4.8 5.7 4.8 4.8 3.8 4.7 5.3 5.4 4.9 4.2 3.1 4.6
##  [127] 4.3 2.4 4.6 5.5 6.0 6.6 5.6 4.3 5.0 4.1 4.5 5.6 4.1 4.6 4.0 4.3 3.7 5.1
##  [145] 3.3 5.1 4.8 3.0 5.3 2.9 4.9 4.8 5.1 4.3 5.6 5.6 4.0 3.8 6.4 3.6 4.2 3.0
##  [163] 3.6 3.1 4.2 4.5 5.0 5.8 3.5 4.6 3.4 4.0 5.1 4.6 5.2 2.9 4.9 5.6 5.3 4.5
##  [181] 3.5 3.6 5.4 3.6 5.3 4.6 4.6 4.9 3.9 4.1 4.8 3.1 5.4 4.4 5.6 4.2 4.0 4.0
##  [199] 3.9 2.8 4.2 3.2 4.6 2.5 5.9 3.7 4.2 5.4 5.7 3.1 4.3 5.0 4.5 3.5 4.0 4.5
##  [217] 3.1 5.0 2.6 5.9 4.5 3.7 3.4 3.4 4.4 4.0 5.0 4.5 6.0 4.8 5.1 5.3 4.0 5.8
##  [235] 3.6 4.3 5.1 3.7 3.9 5.5 4.7 4.9 4.4 4.2 3.3 4.6 4.2 4.4 4.1 5.1 4.4 2.9
##  [253] 5.4 3.6 4.2 5.4 5.2 4.6 4.0 3.5 4.1 6.0 3.5 4.7 4.3 5.3 5.4 4.0 3.2 5.0
##  [271] 4.7 5.8 3.6 4.8 2.6 4.0 3.1 3.5 5.6 3.8 4.0 4.3 6.6 5.8 5.2 2.9 4.6 3.4
##  [289] 4.6 3.8 4.4 4.5 3.4 5.6 3.9 4.3 2.9 5.1 4.4 3.6 4.8 3.7 5.1 2.6 5.5 5.2
##  [307] 3.4 2.5 5.0 4.6 5.5 5.3 4.5 4.1 2.5 4.3 3.4 3.6 5.5 3.9 5.2 4.9 6.0 5.0
##  [325] 5.8 4.2 4.5 4.7 5.2 5.0 3.3 4.4 3.6 3.6 4.7 5.8 4.4 3.4 5.0 3.8 4.7 3.8
##  [343] 4.5 5.8 6.2 5.3 5.1 2.8 4.5 2.8 4.0 4.4 4.8 4.8 5.6 2.7 3.9 4.0 3.8 4.4
##  [361] 4.4 4.3 4.8 3.9 4.2 4.8 5.8 4.3 4.7 3.4 5.7 2.2 6.2 4.3 3.7 5.1 4.1 3.9
##  [379] 3.4 3.8 3.3 5.1 4.7 6.0 5.6 4.8 2.3 4.4 3.8 6.0 4.6 4.5 4.3 3.6 3.7 4.4
##  [397] 4.5 3.6 3.9 4.3 3.1 5.4 4.0 3.9 5.0 3.2 2.3 3.9 5.3 4.8 5.4 3.7 5.1 3.9
##  [415] 4.6 4.7 4.5 5.9 5.2 3.3 4.7 4.3 5.0 4.6 4.8 4.4 4.1 3.6 5.1 4.9 3.5 3.3
##  [433] 3.6 3.5 3.7 5.4 5.4 5.4 5.2 3.9 4.0 4.0 4.7 4.5 4.9 3.9 3.8 5.0 2.1 3.7
##  [451] 2.8 4.1 4.2 5.2 3.1 3.9 4.5 4.1 3.6 4.6 1.9 5.0 4.7 4.0 5.0 4.0 3.7 5.8
##  [469] 6.8 5.0 3.6 5.9 4.9 3.3 4.7 4.1 5.5 4.7 4.0 3.9 5.7 5.4 5.5 4.8 4.6 5.2
##  [487] 5.1 4.8 4.8 4.9 5.4 5.1 4.1 4.6 3.8 3.9 4.7 3.8 4.1 5.2 4.9 5.9 5.3 3.0
##  [505] 4.6 5.0 6.0 2.3 5.7 5.2 4.3 4.2 3.0 4.1 5.2 5.0 3.7 4.6 5.0 3.5 4.3 3.4
##  [523] 3.8 3.6 5.3 3.9 3.5 4.5 3.8 5.4 4.6 5.1 4.5 5.2 4.3 3.9 3.9 5.7 4.4 4.8
##  [541] 4.0 3.5 4.9 5.4 5.2 3.0 4.5 4.6 4.3 3.5 5.5 5.2 3.8 4.1 4.2 4.6 2.8 6.2
##  [559] 2.7 4.4 4.4 5.6 5.5 5.5 4.2 3.3 4.7 4.4 3.1 3.4 4.9 5.4 4.5 5.0 5.0 5.8
##  [577] 4.5 4.8 4.3 4.9 4.3 4.4 3.4 5.0 5.9 4.5 4.3 3.3 3.7 5.6 4.0 5.3 5.4 5.5
##  [595] 3.4 6.4 4.7 3.6 4.8 5.6 4.3 3.4 4.4 5.8 3.2 5.0 4.9 4.1 4.2 4.7 4.0 3.9
##  [613] 4.5 3.4 3.7 4.2 6.3 5.2 4.1 4.0 3.7 5.0 3.8 4.3 3.2 6.7 5.8 5.3 6.4 4.6
##  [631] 4.7 4.3 4.7 5.5 4.8 4.9 3.8 4.0 6.1 6.7 2.4 5.0 5.3 4.8 3.9 4.2 4.6 3.6
##  [649] 4.6 5.2 5.1 4.6 3.6 6.7 2.8 5.1 4.5 4.4 6.8 5.2 4.3 5.5 6.0 4.7 5.3 6.0
##  [667] 4.1 3.4 6.7 4.6 3.6 5.1 5.1 5.6 3.9 4.6 3.5 3.6 4.6 4.3 4.0 2.5 3.3 4.2
##  [685] 3.4 4.3 5.4 5.8 4.5 5.3 3.7 2.9 4.6 3.1 4.8 4.4 5.3 4.9 4.9 5.3 5.7 3.9
##  [703] 3.4 3.7 3.3 4.9 6.9 5.8 4.1 3.0 4.0 3.8 3.5 6.7 3.5 4.9 5.7 6.0 4.2 3.7
##  [721] 3.8 3.5 3.0 5.6 3.3 3.9 3.5 6.5 5.2 3.4 4.9 2.9 3.6 5.0 5.1 2.9 3.4 4.3
##  [739] 5.0 3.2 3.2 5.3 3.9 5.2 6.0 4.8 5.3 4.4 6.0 2.8 5.2 2.7 2.5 4.2 2.9 4.0
##  [757] 6.1 5.3 4.6 4.2 4.9 4.1 3.7 5.0 3.1 3.8 4.0 5.4 6.1 4.4 5.7 3.3 4.0 3.8
##  [775] 4.2 5.3 4.7 4.1 3.1 3.8 3.5 4.9 3.7 4.1 5.6 5.4 4.4 5.6 3.1 5.1 4.2 4.7
##  [793] 4.4 4.5 4.0 4.7 4.2 3.2 4.3 5.8 4.5 2.8 3.1 6.2 4.5 6.4 4.1 3.3 4.7 7.0
##  [811] 5.0 4.7 4.5 2.4 4.8 4.7 4.5 4.8 3.7 4.2 3.4 3.0 3.8 3.8 3.4 4.1 4.0 4.3
##  [829] 4.4 5.4 5.4 3.9 4.1 5.5 5.2 4.7 5.0 5.9 4.6 4.4 5.2 4.3 3.7 6.5 4.6 4.9
##  [847] 4.5 4.1 4.9 4.7 3.2 5.2 5.0 5.2 4.9 6.2 5.4 4.9 4.8 4.7 4.5 2.9 4.9 4.4
##  [865] 3.2 3.9 4.8 5.7 4.0 3.2 4.4 3.4 4.7 4.0 5.1 5.3 4.0 1.9 3.2 4.0 4.1 4.5
##  [883] 5.0 5.2 6.0 4.1 5.6 4.5 3.8 5.4 4.3 4.8 4.2 3.8 4.5 4.5 4.4 4.6 3.3 5.1
##  [901] 2.8 4.1 3.8 3.5 3.3 5.8 3.7 3.8 4.7 4.5 3.0 4.4 4.8 3.1 4.3 5.6 5.7 6.5
##  [919] 6.6 3.2 5.2 5.3 3.8 3.3 5.2 5.8 5.4 4.3 4.9 4.4 6.0 4.0 3.6 4.2 4.5 5.6
##  [937] 3.6 5.8 4.2 4.8 4.4 4.3 5.1 4.4 3.4 3.8 4.7 3.3 5.3 3.9 3.9 4.6 5.1 4.9
##  [955] 5.3 3.4 5.7 4.8 4.5 2.7 3.3 4.8 6.0 4.8 5.0 4.6 5.5 3.8 2.9 6.3 3.4 2.6
##  [973] 3.8 5.5 3.6 6.0 5.8 4.1 5.2 3.9 5.3 2.2 6.2 5.7 4.2 6.2 4.2 4.3 2.5 5.7
##  [991] 4.3 6.3 4.9 4.6 2.9 4.6 3.8 5.2 5.1 3.9
## 
## $func.thetastar
## [1] -0.0384
## 
## $jack.boot.val
##  [1]  0.488495575  0.328484848  0.271987952  0.053314121 -0.009638554
##  [6] -0.073262032 -0.221114370 -0.313424658 -0.480229226 -0.495988539
## 
## $jack.boot.se
## [1] 0.9638303
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
##    [1] 3.6 4.7 5.2 4.0 5.3 4.7 3.7 3.9 4.1 3.5 5.4 3.7 4.9 4.9 2.9 4.9 4.7 4.5
##   [19] 2.9 5.4 4.9 5.0 2.6 3.4 5.5 4.3 4.6 4.7 4.7 4.7 4.1 4.3 5.9 1.9 4.5 5.1
##   [37] 4.5 3.2 4.5 5.1 3.7 5.0 4.0 5.5 4.1 3.6 3.0 4.8 4.0 4.8 3.5 5.1 4.9 4.6
##   [55] 5.1 4.5 4.2 4.9 5.6 4.7 3.3 5.8 4.2 3.9 4.8 5.5 5.9 4.0 3.0 4.9 3.8 3.9
##   [73] 3.1 4.3 4.8 3.0 3.0 3.7 4.1 5.2 3.8 4.6 2.5 4.3 5.5 4.6 3.1 5.6 5.5 5.8
##   [91] 3.1 4.2 4.3 3.4 4.8 3.5 5.3 4.6 3.7 3.3 4.7 4.3 2.7 3.3 3.5 4.0 5.2 4.5
##  [109] 3.3 5.6 3.6 5.4 5.4 3.7 4.8 2.7 5.8 4.7 4.4 3.8 5.2 4.5 5.2 4.7 4.7 4.7
##  [127] 5.3 4.1 5.2 4.3 5.3 4.9 5.9 4.9 3.0 4.6 4.3 2.0 4.9 4.9 4.1 3.8 2.8 3.8
##  [145] 3.7 3.8 2.5 5.2 7.0 3.5 5.9 4.8 4.0 4.9 3.6 4.9 4.7 4.3 3.9 4.3 5.6 5.1
##  [163] 6.5 5.5 4.1 5.0 4.0 3.5 4.9 4.1 4.6 5.6 5.9 5.2 4.7 4.7 3.2 3.5 6.1 5.4
##  [181] 3.3 5.0 5.6 4.0 3.1 5.9 4.7 6.0 4.2 5.2 4.6 4.7 2.4 5.8 5.0 4.8 4.0 4.9
##  [199] 3.9 5.0 5.6 3.6 4.6 2.2 2.1 6.0 3.6 3.9 3.8 2.9 4.2 4.3 5.1 5.1 5.4 4.2
##  [217] 6.1 4.4 4.2 3.9 6.4 3.1 4.3 3.2 3.9 4.1 4.8 6.1 3.6 4.3 5.1 3.0 4.3 3.4
##  [235] 4.3 6.2 4.4 2.4 3.9 5.0 5.1 5.6 5.8 4.3 4.7 4.8 3.4 5.4 3.6 4.9 5.5 3.5
##  [253] 3.6 4.8 4.3 4.4 4.9 4.3 5.2 5.2 4.2 3.5 5.7 3.6 4.0 3.4 6.5 4.9 5.2 3.8
##  [271] 3.7 3.4 4.8 6.9 5.0 4.3 3.9 4.6 3.3 6.5 4.9 4.4 6.1 5.5 5.6 4.7 4.7 5.2
##  [289] 6.2 4.7 4.3 4.2 3.6 5.8 4.5 2.3 4.0 5.9 4.7 4.0 4.6 5.2 3.3 3.8 4.9 3.4
##  [307] 6.8 6.9 4.1 4.3 2.2 3.7 4.3 4.2 4.7 6.0 3.8 5.0 4.1 4.9 4.7 5.0 4.2 4.5
##  [325] 5.1 3.8 3.5 4.2 4.9 4.3 3.1 5.3 5.9 5.2 4.7 4.2 5.3 4.6 3.3 3.5 5.6 5.2
##  [343] 3.0 1.8 6.0 3.3 4.8 5.9 5.7 3.9 5.3 6.1 4.0 3.4 4.8 6.3 5.1 3.3 4.0 4.1
##  [361] 3.9 3.6 4.0 3.2 4.8 4.8 4.4 4.2 4.5 3.1 5.4 5.3 4.7 5.0 4.4 4.9 4.7 4.3
##  [379] 3.6 5.0 5.8 5.1 4.7 5.3 4.4 4.4 4.4 3.7 4.2 4.6 3.6 4.6 5.6 5.2 4.2 4.1
##  [397] 5.3 4.6 6.1 4.5 3.9 2.9 3.4 3.5 5.9 4.0 4.4 3.3 4.8 4.0 4.3 4.5 3.6 4.8
##  [415] 4.2 3.5 4.2 4.7 4.4 3.7 4.1 5.9 4.7 4.9 4.3 4.7 6.0 4.4 4.5 3.8 4.2 5.1
##  [433] 5.5 6.2 4.7 5.5 2.9 4.9 5.5 6.5 4.3 3.4 3.8 4.9 3.7 4.8 4.9 3.8 4.5 4.8
##  [451] 2.8 5.5 4.9 3.2 2.9 4.6 4.5 4.6 6.2 3.1 5.4 4.9 5.4 5.7 3.4 4.4 5.4 5.9
##  [469] 2.0 4.2 4.8 3.5 4.2 4.4 4.6 4.8 5.1 4.7 4.7 5.6 4.2 4.0 4.5 4.9 5.5 4.6
##  [487] 5.7 4.6 3.9 5.7 4.1 4.1 3.4 3.8 5.2 5.1 5.9 4.5 4.0 4.3 3.5 4.3 5.6 3.8
##  [505] 4.5 4.3 5.3 6.9 6.1 4.3 4.7 5.6 3.9 5.9 4.2 3.6 3.5 4.5 3.7 5.5 4.2 3.6
##  [523] 4.1 4.1 5.5 3.7 5.2 5.2 4.3 5.5 4.5 5.7 5.2 5.6 5.2 3.8 4.6 3.6 4.6 4.3
##  [541] 4.5 4.0 3.5 4.1 5.2 4.2 4.5 3.8 5.5 3.8 5.7 5.2 4.5 4.6 3.1 4.4 4.5 5.1
##  [559] 4.1 3.4 3.6 2.9 4.4 4.9 5.8 5.1 2.4 3.8 4.7 4.6 5.6 4.0 4.2 5.2 5.0 4.9
##  [577] 4.9 6.8 3.9 4.3 3.9 5.9 4.3 5.3 3.8 4.0 5.1 5.0 6.6 4.5 3.1 4.7 6.0 4.5
##  [595] 4.6 5.2 4.4 5.7 3.7 4.9 6.2 4.8 4.5 5.4 5.9 5.2 4.1 5.4 4.0 4.8 5.1 5.4
##  [613] 4.6 5.0 4.0 4.3 5.5 4.9 4.4 4.3 4.7 4.7 5.6 4.7 4.1 5.6 4.9 4.4 4.4 5.5
##  [631] 5.4 4.9 3.0 3.1 3.8 3.6 3.8 4.4 4.2 4.7 6.0 4.2 6.0 4.9 3.2 4.9 4.8 3.9
##  [649] 4.3 4.1 3.7 3.4 3.4 2.8 5.5 5.6 3.2 4.4 2.4 5.6 3.8 4.5 3.8 3.9 3.3 4.8
##  [667] 3.5 5.3 6.2 5.3 2.8 4.8 4.2 4.1 4.6 4.6 4.0 5.3 5.7 5.2 4.1 5.8 2.3 3.6
##  [685] 6.5 4.9 4.6 3.2 5.8 3.2 5.4 4.4 5.6 5.4 4.3 4.5 3.5 5.0 4.1 5.7 5.4 5.0
##  [703] 5.5 3.8 4.7 4.0 4.6 4.5 3.8 4.3 6.3 3.8 4.2 3.3 5.0 3.3 2.6 2.5 3.7 4.5
##  [721] 5.2 4.3 3.5 3.1 4.8 4.4 4.3 3.9 2.2 3.6 4.4 6.4 4.9 3.3 5.2 3.3 4.4 3.9
##  [739] 3.7 5.9 4.9 3.8 4.3 4.8 3.8 4.3 5.5 4.3 3.3 5.7 5.1 4.2 4.9 4.5 3.9 3.7
##  [757] 3.5 5.0 3.7 3.6 3.5 4.4 4.6 3.8 5.3 3.2 3.5 4.8 4.0 3.9 5.7 4.1 5.1 4.3
##  [775] 5.0 5.3 3.7 3.1 5.4 3.5 6.2 3.4 3.6 6.3 6.0 5.5 4.5 3.3 4.3 3.4 4.7 5.8
##  [793] 3.6 4.9 5.5 4.6 4.9 4.6 5.3 4.2 4.7 5.6 4.2 4.7 3.7 4.8 4.2 5.2 4.1 5.6
##  [811] 3.0 5.0 4.3 5.2 3.6 4.0 5.0 4.7 4.2 3.5 4.6 4.2 5.3 4.3 3.0 2.7 5.1 4.8
##  [829] 5.4 5.8 5.5 3.8 3.7 3.9 5.1 3.4 3.7 4.1 4.4 5.0 4.4 6.0 3.7 4.6 4.0 4.0
##  [847] 3.9 4.2 4.8 4.1 5.3 5.6 5.6 5.7 3.8 6.1 3.2 5.2 5.3 5.5 3.2 5.5 5.9 3.9
##  [865] 6.4 4.3 6.9 4.7 4.9 3.8 5.0 2.6 4.1 5.1 3.8 3.4 5.7 4.8 2.8 5.5 3.8 5.4
##  [883] 3.3 6.1 4.9 6.6 3.7 4.2 3.8 3.7 2.6 4.7 5.1 5.2 3.5 5.5 4.2 4.4 2.1 3.8
##  [901] 3.2 5.7 4.2 5.0 3.6 5.0 5.8 4.7 6.4 3.0 5.5 3.2 3.7 4.7 2.7 3.7 4.4 2.7
##  [919] 4.9 5.3 2.6 3.7 4.4 4.4 5.4 3.6 4.9 3.5 3.9 3.7 4.5 4.1 3.8 4.8 4.9 2.9
##  [937] 4.6 4.2 4.0 4.0 4.4 6.2 5.0 6.4 5.5 3.4 5.0 4.1 4.8 4.8 4.3 4.5 4.4 5.6
##  [955] 4.2 4.7 4.5 3.9 5.5 5.0 4.2 4.8 5.2 4.4 3.0 4.7 4.7 4.2 6.9 4.7 3.3 5.7
##  [973] 5.0 3.0 6.4 5.8 3.8 4.0 4.6 4.4 6.4 5.3 4.0 4.9 3.7 4.3 3.3 3.8 4.1 4.7
##  [991] 5.5 2.5 6.2 4.7 3.5 4.3 5.9 5.0 4.5 4.5
## 
## $func.thetastar
## 72% 
##   5 
## 
## $jack.boot.val
##  [1] 5.5 5.4 5.4 5.3 5.1 5.0 4.9 4.8 4.6 4.4
## 
## $jack.boot.se
## [1] 1.049571
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
## [1] 1.241525
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
##   4.214008   4.744052 
##  (1.814896) (2.169941)
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
## [1]  0.5227427 -0.9672462  0.6237564  0.1440538 -0.2006455  0.5917345
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
##    [1]  1.2518131499  0.2766875492  0.7567784910  1.8173282219  0.3828842597
##    [6]  1.0748147814  0.3663517493  1.9521076207  0.9888418641  1.6951003494
##   [11] -0.4139603199 -0.1508236545 -0.1678628212 -0.0009440321  0.4468890381
##   [16]  1.2231366992  1.8346119510  0.6819219746  0.0486260650  1.4134725323
##   [21] -0.4822776893  1.5721125031  1.2396371808  1.5481530848  0.5987134489
##   [26]  1.6629729380  0.0315546845  1.7527689794  1.1163851602 -0.6395886456
##   [31]  0.0710033894  1.0008457177  1.2558232703  0.8146826960 -0.1403571661
##   [36]  1.4649299048 -0.8307876017  0.9948441835  1.5506003805  1.6763000749
##   [41]  0.4132231206  0.7397492950  1.0360167191  0.2204433483  0.5809865158
##   [46]  0.1933341471  0.8088412062  0.8477970297  0.5668978317 -0.4068291321
##   [51]  0.5897869859 -0.7967174680  0.5919800591  1.0256560587  2.0266343822
##   [56]  1.5247768293  1.2679490491  0.0371714058  0.1883512871  1.4793355918
##   [61] -1.1151281658  1.2186833783  1.6443962847  0.5655274287  0.6588106120
##   [66]  1.3875295495  1.0292207459 -0.6192674404  1.4274956325  1.1536652106
##   [71] -0.1371586193  0.5858508814  1.1422806863  1.6116388073  0.5731754367
##   [76]  0.7193374623  2.0934701190  1.2187647967  1.1694984050  1.5499123434
##   [81]  0.1836273899  0.8698483507  0.7105825503  1.0576157001  1.3680707295
##   [86]  1.3803339200  1.1134330340 -0.1013523206  1.3506822982  1.0951046544
##   [91]  1.2204068115  1.2603580362  0.9264647244  1.0866440616  0.5394395776
##   [96]  0.2941485076  1.9918422379  0.4148856423  1.0059503564  2.1398354265
##  [101]  0.8816927553 -0.1334748216 -0.2366045585 -0.4205077423  0.3397375786
##  [106]  1.5154009487  1.0433323608  2.3221041661 -0.3431425220  1.4873840784
##  [111]  2.4677883816  1.5276200981 -0.3010941105  1.9415071824  0.7520248078
##  [116] -0.5792061010  0.5972435359  0.5768187237  1.2417369797  0.6665353383
##  [121]  1.2406442758  0.4050934840  1.2066782323  0.5775404498  0.5749207303
##  [126]  1.7459139987  0.0400915693  0.6203076278  2.1299857625 -0.6930615800
##  [131]  1.3460786852  1.3620297381  1.3742256176  1.6930157129  0.3892195872
##  [136]  0.7576850125 -0.7842935853  0.7182324890  1.6460448787  0.2928204878
##  [141]  0.7181989317  0.4675443380  1.6930157129  0.8591551775  0.3982330564
##  [146]  0.9715632547  0.8339683783  0.7577835668  1.0786074448  0.0710999521
##  [151]  0.4230736128  1.3988381821  0.5688416799 -0.7419324340  0.6750359074
##  [156]  0.6259629584  0.5419866280  2.0429069952  0.7850870737  0.5932916921
##  [161]  1.9862279700  0.9183441847  0.8807693736  1.1117693903  1.0560057501
##  [166]  0.3353747355  0.0926620024  1.7439198410  0.8822784246  1.2463145220
##  [171]  0.1466317480  1.1107731493  0.5649715875  0.8477970297  0.7027895888
##  [176]  1.4873264244  1.9394597893  2.5330043528  1.1292113371  1.3270252422
##  [181]  1.0580585936  1.2452672008  0.0966527517  0.6711656267  0.3286178844
##  [186] -0.5721951955  0.5435904202  0.6458734599  0.4967458114  1.0576935277
##  [191]  1.2241455276 -0.5914976102  0.4427061727  0.7191395470  1.4861162142
##  [196]  0.9340677565  1.3059781133  0.0434251626  0.5092665019 -0.1263600680
##  [201]  0.6005215601  0.0687061510 -0.5304933837 -0.3846065936  1.7996765341
##  [206]  1.6840674948  0.2849055143  1.9520068549 -0.6628969530 -0.1283724158
##  [211]  0.1628513989 -0.0742656905 -0.0841125568 -0.2406955871 -0.1315990084
##  [216]  0.0789228655  0.7788110232  1.0563939110  0.8464617765  1.9937385642
##  [221]  0.8844380696  0.1880558534  1.2302202698  0.6864932866  1.3747450544
##  [226] -0.0650859101  0.7970016318  1.2393915896 -0.8536134012  0.4382732285
##  [231] -0.1310037423 -0.3561372981  1.3456713766  1.4502830114 -0.3908253977
##  [236]  0.9273951135  0.1474373130  0.1323548467  0.5337761621  0.2363215130
##  [241]  1.8456064567  1.1026990601 -0.7392462577 -0.3941672547  1.2738252748
##  [246] -0.4558430572  0.7914197162  0.2115958740 -0.8703059354  0.6073677691
##  [251]  0.0094755870  0.2981417559  0.1523178066  1.5083769614  1.4311885973
##  [256]  0.7300437247  1.6325489282  1.3970717884  0.2253220866 -1.3040154435
##  [261]  1.2035876520  0.5651444230 -0.6705909347  0.1747322610  0.3521040778
##  [266]  1.1518268818  1.6710182260  0.7395570933  1.1092294445  1.4867954463
##  [271]  1.2514745974  1.7912126908  0.2216909791  1.6187901550  0.9563604829
##  [276]  1.6087034952  1.2978895907  0.4845507029  1.0738409929  1.6687676985
##  [281]  0.4447322286  1.0524084509  1.2097702491 -0.7214699101  0.3957398567
##  [286]  0.3413034415  0.6171913517  0.9590405457  1.5093908154  0.9465748981
##  [291]  1.0863403785  1.5650277845  1.2425234214  1.4165122090  1.2122930151
##  [296]  1.7303756551  1.2063043905  0.9479764152  0.3481793187  0.2592300358
##  [301] -0.2341269016  1.2702727594 -0.3992524521 -0.4742005270  1.1470599322
##  [306]  0.4996392998  0.4673535757  0.8307315562  1.8345791768  0.7133070217
##  [311] -0.1000232814 -0.4778710601  1.2227861441  1.4576179591  0.1137687918
##  [316]  1.9853791887  0.0420642752  0.6224456165  1.9512676968  1.0809366350
##  [321]  1.3725056710  1.6325501840  1.7613896579  1.3828396973  1.7113039034
##  [326]  0.5366732412  0.3710343091  2.1315032195  1.5212280032  1.0621638192
##  [331] -0.2141503507  0.4990921390 -0.7265194606 -0.2933270333  1.2554624844
##  [336]  1.1244594018  0.7261423902 -0.8495683557  1.0196699199  1.0799001819
##  [341]  1.6565616027  1.5665177526  1.0750997269  0.8975323591  1.1670505252
##  [346]  1.0112658652 -0.1377948511  0.7009919996  1.4143183252 -0.6628088504
##  [351]  1.5890896191  0.8782504894  0.0382785385  0.6499846743  0.7974275148
##  [356] -0.1383804233  0.3370790625  0.7348967419  0.7181495888  1.0877968258
##  [361]  0.9528612716  0.1592722899 -0.5459468173  1.4052285726  1.5071327718
##  [366]  1.0450963907  0.8737480224 -0.1770150740  0.5949666353 -0.5641687927
##  [371]  0.5482496723  0.4593127917  1.3631294472  0.2052591137  0.1284765667
##  [376]  0.5561307039  1.3671968410 -0.4018738155  0.8017068524 -0.2339969929
##  [381] -0.2127189707  0.1906989976 -0.1887614621  1.7472935684  1.5019991532
##  [386]  0.7692934189  1.4031326255  0.1975929055  1.3019992129  0.0132026477
##  [391]  1.0473255334  1.0360734645  0.4411493861  0.9472427220  1.4961286671
##  [396]  2.0844252227  0.3347039665  0.9536515438  0.6048234044  1.5964034589
##  [401]  1.1197276678  0.0334390003  0.9157623574  1.2736220658  0.7346251452
##  [406]  0.3755658034 -0.5057510546  1.7622522623  1.1145441766 -0.0860303004
##  [411]  1.6522779907 -0.0177065544  0.1771784830  0.7768872267 -0.1578924851
##  [416]  1.0047117202 -0.1896069051  0.5883339446 -0.1854876673 -0.5099310079
##  [421]  1.8472172995  0.8455933541 -0.4158920733  0.7146328071  1.0123628304
##  [426]  1.2345459219  0.6927369478  1.6972351879  0.1228816489 -1.1411954755
##  [431]  1.6311523160  1.1473368624 -0.4997417972  2.0177837379  0.1488668292
##  [436]  0.5258025272  0.0540181469  0.5429842711 -0.3229396758  1.2039680085
##  [441]  0.1506606351  0.5183465006  1.0653598025 -0.6409614396  1.9221717672
##  [446]  1.5890896191 -0.2485817330  0.2552569004 -0.6802046139 -0.3519749910
##  [451]  1.2531990589  1.6902887896 -0.1378747433 -0.0808527682  1.2029267720
##  [456]  0.1991430605 -0.3908253977  1.1693968065  0.6417818853  1.1043449961
##  [461]  1.3174657423  0.3090738551 -0.9332698092  0.9685585829  0.2526408605
##  [466]  0.7850760181  1.9876703053  1.9146611901  0.7713876899  1.2848484527
##  [471]  0.1264564327 -0.5242443100  0.9284062164  0.0753418833 -0.0545918291
##  [476]  0.1457808253  1.4395804557  0.3427727651  1.5721486986  0.1981937211
##  [481]  1.1497466412  1.4980271866  0.8877713278  0.4544006381  0.4925116561
##  [486]  1.8966290246 -0.7589533619  0.4634712420  0.9317229873  0.1497134787
##  [491] -0.5551406680  0.5126023238  0.7133259683  1.8324361989  1.5686052919
##  [496]  1.5270407926  0.2505369522  0.6999765839  1.3941338309  0.8809601294
##  [501]  1.5363957978  1.2054147461  1.2489774406  0.6699165950  0.1687918926
##  [506]  0.1724192925  1.5528227285  1.0396889246  1.1296200153  0.8773782487
##  [511]  0.6830683130 -0.1054574186  0.8155119246 -0.2107248618  0.1489639448
##  [516]  1.9157683991  0.6829579148 -0.2683626590  0.5764464157  1.7202192443
##  [521]  0.4283393067  0.9036638068  0.2460207287  0.3030305666  0.9581660921
##  [526]  1.4855362872  1.1581743274  1.9391250044  0.6905336357  1.0327025170
##  [531]  1.6539237748  0.8104350746  0.3922929182 -0.2072978791  0.2866618010
##  [536]  1.2330673260  0.7002513004 -0.1603539890  0.8664890013  0.6689876205
##  [541]  0.7576097273  1.4876096405  0.6569474241  1.7562239240  1.4448819327
##  [546]  0.9004016553 -0.7841095976 -0.0706827364  0.6978962800  1.6378779062
##  [551]  0.3927084050 -0.2385665532  0.3941712022  0.4698078817  1.0930385263
##  [556] -1.1218765587  1.5600065506  2.0746187707 -0.0297548986 -0.9380837683
##  [561]  0.8565608172  1.6972351879 -0.3956124648  2.1326594168  1.0861403928
##  [566]  1.8658829191  1.3186689160 -0.0646628935  1.1604344311  1.4516921123
##  [571]  0.3553211346  0.6050604531  0.7831219110  0.8014154415 -0.1355237139
##  [576]  0.6197136995  0.1227100194 -0.2995836163  0.7611027990  1.1518268818
##  [581]  1.2129676040  1.5873537271  0.5646458289  0.0257656656  1.4130477085
##  [586]  1.4181261878  1.2798261446  1.2063853081  1.4117404886  1.1934909593
##  [591] -0.3988388119  1.3764625255 -0.2535112097 -0.8340205729  0.8024483586
##  [596]  0.7105227975  1.8514334195  0.5423875454  0.1209430233  0.9645837846
##  [601]  1.0496375091 -1.1365029947  2.0837629516  1.2467627464 -0.0282442843
##  [606]  1.4019346285  0.8012259281  1.6515388048  1.7102786791  0.8309442015
##  [611]  1.6955974863  0.7291705099  0.4100879394  1.0361697348 -0.0760602518
##  [616]  0.1236649152  0.3736210913  1.2983733876  1.4507984676  1.7934293751
##  [621]  0.8180486889  1.0570383853  1.6544779248  1.6891749058  1.8198045769
##  [626]  0.3391825737  0.5563023440  0.7974462379  1.8235243800 -0.8596073186
##  [631]  2.0665915478  0.8005959857  1.1175240695  0.9766244722 -0.1055329147
##  [636]  0.5228660792  1.7616296957  0.5460968885  1.2840420181  0.0502503104
##  [641]  1.4031866703  0.2468629862  0.6348167928 -0.5738184972 -0.5394682428
##  [646] -0.3449367329  2.2018784193  1.3782880278  0.9778153980 -0.5031264579
##  [651]  0.8977985689  1.8045658940 -0.8632104949  1.1984241911  1.1374164951
##  [656] -0.1881278091  0.4437745877  1.5771097587  1.9310864314 -0.4960066898
##  [661]  1.0779863573  1.8738122902  0.0813560479  1.2904785389  0.6681826404
##  [666]  0.9815508111  0.6489307019  1.2936477825  0.1507884904  0.5957974573
##  [671]  0.8303075874  2.1743089710  0.7797270657 -0.5623129468  0.6751095163
##  [676]  1.3448264952  0.7099252583 -1.0484976991  1.8498814170  1.2262236053
##  [681]  1.2459542817  1.3104448793 -0.3475303159  1.6720471846  0.9692210588
##  [686]  1.1452533567  0.7485329518  0.6127565594  0.4629902383  2.0428704110
##  [691]  0.9553407906  0.8992757818  0.9760589030  1.0791426852  0.7538815902
##  [696]  0.1841612126  1.6402925882  0.5780322325  1.1035592332  0.6652441154
##  [701]  1.9676726677  1.4334891540  0.6917223883  0.9855481570 -0.1290570916
##  [706]  0.7903435962  0.9585491122  0.7415174803  0.1511314017  0.8139401740
##  [711]  0.3435563321  0.8050932012  0.4034857459  0.9769124135  0.7866481470
##  [716] -0.1737732537  1.7058167530  1.1368429171  0.9643150773  0.0449752892
##  [721]  0.3610354363  0.4133712422  1.2535006723 -0.3853565405  1.0733913083
##  [726]  0.2411158625  0.9086994981  0.0077200283 -0.7598071398  1.1202120285
##  [731]  0.9983309824  0.4815277146  1.0074837171  1.1899611027  1.4552291916
##  [736]  1.5099135965  1.4536232967  1.2071854541  1.1168910532  1.3420515220
##  [741]  0.0789785233  1.0083522400  1.3683264606  0.2007609251  0.8612727542
##  [746]  0.6464031587 -0.6740020695  1.6985771071  0.7642494601  1.5230113683
##  [751]  0.8422981829  0.6744074641  0.8414983473  0.1631653764  1.3193081991
##  [756]  1.2333425483  0.8250080443  0.3413034415  0.9352776640  0.5935993172
##  [761] -0.3640372517  1.0119928510 -0.0405123380  0.8127523768  1.4980088217
##  [766] -0.3907636221  1.2050728561  1.2228762458  1.4498286979  1.2334832413
##  [771]  0.5265355312  1.0319950221  0.9715907996  0.2688285047  0.1082019902
##  [776]  0.7322605743  0.8494184096 -0.1263727008  1.0412598773  1.9146797556
##  [781] -0.4907611049  1.1286860111 -0.7890845368  2.2696232799  1.6722868516
##  [786]  1.9669701450  0.7574025438  0.8356346756  0.5588188780  1.0781988534
##  [791]  1.1907719853  0.8228310380  0.7853586624  0.1938472969  0.4595785147
##  [796]  1.6331600051  0.1076439909 -0.1908950845  1.8320216140 -0.2106073677
##  [801]  0.0225195125 -0.7821674667 -0.0629674053  0.9600719163  1.2283776169
##  [806]  1.2757671455  0.2596619114  0.6828805851 -0.3538708878  0.5891846569
##  [811]  2.2924647633  0.8046022320  1.3583939016  0.9707917000 -0.6308844836
##  [816]  1.0863262184  1.1542345975  0.0674875699  1.1912458132  1.5902802098
##  [821]  1.0997153502 -0.9450394918  0.6246702285  0.5925542672  1.1815794862
##  [826]  1.2088858262  0.5475012902  0.5579048967  1.1107404421  1.3846576630
##  [831]  1.5223151388  0.7938780800  0.3140041169  1.2726400247  0.2234703465
##  [836]  1.1155440735  0.5639209243 -0.2470248116  0.9306432842  1.2196417773
##  [841]  0.8749343345 -0.5852674209  0.9672190422  2.2908545774  1.8856586694
##  [846]  0.3139090705  0.9703674763  1.2684086406  1.3671676932 -0.6762473302
##  [851]  1.9938952573  1.0415598043  1.2505119279  2.2563546411 -0.6525909540
##  [856]  0.4385222743  0.5367203737  1.3445141745  1.7074160861  1.1710213904
##  [861]  0.8499617882  1.6064622311  0.0876733473  0.4753602992  0.0915465451
##  [866]  0.5610905648 -0.1637812298  1.2724352482  0.3996727521  1.9305429170
##  [871]  1.1784765722  0.3396252371 -0.6561150703  1.1050822435  1.3524071105
##  [876]  0.6927377714  1.1877945955  0.2200109074  0.9460168081  0.6635402949
##  [881]  1.0368098194  1.3025786707  1.0217052290  1.8877788917  1.7138303441
##  [886]  0.0247954363  1.0328830132  1.7112895945 -0.2099857320  1.0397274563
##  [891]  1.1205978651  0.7864621512  0.9557186421  1.5941827781  0.3673394027
##  [896]  0.3838937676  0.9504831104  0.1562267766 -0.1648671263  0.9750670247
##  [901]  0.2487641018  1.1479929968 -0.1564552033  0.8081315039  1.0951190214
##  [906]  1.7422195957  0.6334627919 -0.0089337231 -0.1837530271  2.0747849145
##  [911]  0.5506012361  0.6311393786  0.4080965207 -1.0195910733  0.3274578528
##  [916]  0.4403429808 -0.4551047638 -0.0338420418  1.0415598043  0.7731654015
##  [921]  1.2251925285  0.2929812954 -0.4785909877  1.2771215314  1.7146644329
##  [926] -0.2597354142  1.2065605164  1.2440870243  0.8482826215  1.2929808622
##  [931]  0.5960757883  0.3762151510  0.6489587738 -0.3927632590  1.7587924293
##  [936]  0.8238875283  1.3694977663  0.9887290834  1.1405598498  1.1504105781
##  [941]  1.2619919133  0.2027445437 -1.1287956486 -0.1955418988  0.3058572265
##  [946]  0.9259601838  0.5359717999 -0.0646628935 -0.1377948511  0.7239571393
##  [951]  0.9420243963 -0.6092616378  0.0267608793  1.3240866328  0.4660349144
##  [956]  1.9876703053  0.3431095476  0.8340368225  0.7475431185  0.9603494755
##  [961]  1.4727378711 -0.9049636641 -0.0036696605  2.3954330135  1.4566385372
##  [966] -1.0111224953  0.5696184892  0.3423662054  1.1243643488  1.0295450680
##  [971]  0.3140117255  1.9654406872  0.8527492077  0.8274327723  2.1415351974
##  [976]  1.0624516375  1.4004454600  1.2053078684  0.0060794450  0.7117458072
##  [981]  0.7772629270  0.6672618839  1.3163790965 -0.2127189707  2.2174101086
##  [986]  1.3429623577  1.2780802342 -0.1432638769  1.0928094138  1.6586628464
##  [991]  0.7332225088  0.3240682376 -0.1623964552  1.8705352256  1.2814556409
##  [996]  0.2067782127  0.5290499068  0.1413136966  1.2099721963  0.5864836290
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
##   0.8882736   0.4561433 
##  (0.1442452) (0.1019947)
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
## [1] -0.52125767  0.63834524  0.86792658  1.26255150  0.20954244 -0.07081153
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
## [1] -0.0381
```

```r
sd(results2$thetastar)  #the standard deviation of the theta stars is the SE of the statistic (in this case, the mean)
```

```
## [1] 0.9162302
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
## t1*      4.5 0.006606607   0.8666248
```

One of the main advantages to the 'boot' package over the 'bootstrap' package is the nicer formatting of the output.

Going back to our original code, lets see how we could reproduce all of these numbers:


```r
table(sample(x,size=length(x),replace=T))
```

```
## 
## 0 1 2 3 4 6 7 8 
## 2 1 1 1 1 1 2 1
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
## [1] 0.0173
```

```r
se.boot
```

```
## [1] 0.833658
```

Why do our numbers not agree exactly with those of the boot package? This is because our estimates of bias and standard error are just estimates, and they carry with them their own uncertainties. That is one of the reasons we might bother doing jackknife-after-bootstrap.

The 'boot' package has a LOT of functionality. If we have time, we will come back to some of these more complex functions later in the semester as we cover topics like regression and glm.

