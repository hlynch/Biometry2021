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
## 0 3 4 5 8 9 
## 1 2 2 1 2 2
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
## [1] -0.0329
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
## [1] 2.6565
```

```r
UL.boot
```

```
## [1] 6.2777
```

Method #2: Simply take the quantiles of the bootstrap statistics


```r
quantile(xmeans,c(0.025,0.975))
```

```
##  2.5% 97.5% 
##   2.8   6.3
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
##    [1] 4.3 4.5 3.5 5.0 3.0 2.7 3.8 5.6 5.7 4.7 5.3 4.8 5.5 5.1 5.3 6.5 4.1 5.1
##   [19] 3.0 4.0 4.7 5.5 4.4 5.5 4.1 4.7 4.6 5.6 3.9 5.8 5.3 5.4 4.5 5.2 5.2 4.3
##   [37] 2.3 5.2 4.6 4.1 3.2 3.7 5.4 3.9 4.3 2.6 4.6 5.1 3.9 4.0 6.0 3.6 4.7 5.4
##   [55] 3.7 5.0 4.7 4.7 4.1 5.7 5.2 6.0 4.9 4.4 4.1 5.0 5.5 4.3 5.3 5.7 5.1 5.6
##   [73] 2.0 5.1 4.1 3.6 5.6 5.1 4.9 3.7 3.9 4.1 4.0 5.3 5.4 3.4 4.7 5.0 3.6 4.5
##   [91] 3.2 4.6 5.0 4.4 4.6 2.6 5.7 4.4 4.9 4.9 5.5 4.5 4.3 4.2 4.1 4.0 3.1 3.9
##  [109] 4.9 4.7 4.0 3.5 5.9 5.0 5.4 4.3 4.5 4.3 3.2 4.3 4.1 4.6 4.6 4.4 4.5 4.6
##  [127] 4.4 3.6 3.9 4.2 3.7 5.4 3.1 4.7 4.9 6.3 3.9 3.7 5.5 5.0 4.4 3.9 3.4 4.0
##  [145] 4.0 4.0 5.9 3.3 3.7 5.9 4.9 5.0 5.4 5.2 4.8 3.7 4.7 3.6 3.9 5.1 5.5 4.4
##  [163] 3.7 5.6 5.4 4.9 4.3 5.3 4.5 2.8 5.4 4.1 4.2 3.7 3.8 4.2 5.4 2.8 4.3 5.0
##  [181] 5.9 6.6 2.8 6.5 4.5 4.0 3.6 4.4 3.7 3.0 4.1 4.5 5.3 4.6 3.9 5.4 4.1 3.5
##  [199] 6.2 3.9 5.8 4.7 6.0 3.8 4.1 4.0 5.5 5.0 3.5 4.1 5.1 4.2 5.5 4.1 3.4 4.4
##  [217] 5.0 3.8 1.9 5.5 4.2 6.3 5.6 2.4 5.0 4.4 3.7 2.6 2.5 4.8 4.7 5.7 2.8 4.0
##  [235] 3.9 3.1 3.8 4.4 2.9 4.7 5.0 4.9 4.8 5.8 5.5 4.4 3.3 4.4 4.8 4.3 4.5 3.8
##  [253] 4.1 5.1 5.0 3.7 3.9 3.9 4.2 6.1 5.2 5.3 4.5 5.1 4.4 6.2 4.5 4.8 4.9 3.2
##  [271] 3.5 5.4 4.1 4.5 3.4 5.0 5.6 4.3 4.7 3.7 4.8 3.2 4.3 3.1 4.5 5.2 3.7 4.2
##  [289] 4.5 5.9 3.6 5.1 4.2 4.1 4.4 5.6 5.0 4.1 2.9 5.0 4.3 5.5 4.8 4.6 5.5 3.4
##  [307] 3.8 5.5 4.2 4.2 4.6 4.5 4.1 3.2 6.0 4.2 4.9 4.4 4.0 2.8 5.3 5.3 3.9 4.3
##  [325] 4.2 5.0 5.6 4.3 5.5 4.1 3.5 4.5 3.8 3.8 3.1 3.1 3.5 4.9 3.4 4.2 5.4 4.0
##  [343] 4.8 4.3 3.8 3.2 3.8 3.9 6.5 5.5 5.2 4.0 3.1 4.4 4.7 4.8 4.8 4.1 5.0 3.9
##  [361] 4.6 5.1 3.9 3.9 4.4 5.0 4.4 4.4 5.2 4.3 5.1 5.8 4.7 4.0 3.2 3.0 4.2 3.5
##  [379] 5.5 4.3 5.5 4.9 3.5 3.9 4.2 5.4 3.5 3.6 3.8 4.2 3.8 3.9 4.2 5.6 3.8 5.0
##  [397] 4.0 4.7 4.0 5.7 5.7 5.8 3.5 4.8 3.3 4.8 4.2 4.2 4.5 5.1 3.3 6.0 4.6 4.1
##  [415] 5.7 3.8 4.3 5.1 4.9 4.3 5.1 4.0 5.0 3.1 5.2 6.7 3.6 4.5 5.6 3.8 5.7 2.5
##  [433] 4.7 5.0 5.1 4.8 5.0 2.8 5.8 4.2 3.8 3.2 3.5 5.0 3.8 3.6 4.9 4.8 5.7 5.5
##  [451] 3.9 4.7 5.0 3.2 5.6 5.7 3.4 4.8 5.6 3.4 4.3 2.7 5.9 4.6 4.2 3.8 4.9 4.2
##  [469] 5.8 3.7 5.1 4.7 4.6 3.4 5.4 4.1 5.0 3.5 6.1 4.5 4.9 4.7 4.3 4.8 5.0 3.1
##  [487] 3.1 3.9 4.7 4.8 3.5 3.5 5.9 4.2 3.8 6.1 4.2 4.6 4.4 7.1 4.6 4.0 3.4 4.5
##  [505] 5.2 4.6 5.7 4.3 3.8 5.1 4.4 4.1 4.9 3.7 5.6 3.5 4.5 5.4 4.7 3.9 3.6 4.3
##  [523] 4.1 6.1 4.4 3.2 3.7 5.9 3.7 4.6 4.5 2.3 4.4 3.8 2.8 4.2 4.1 4.7 2.9 2.5
##  [541] 5.1 5.5 4.1 4.3 4.2 3.5 2.7 3.8 4.6 5.5 4.8 5.4 5.2 3.2 4.9 5.7 5.0 3.4
##  [559] 5.8 4.1 5.0 5.6 3.5 4.3 4.5 2.0 3.7 3.4 5.1 4.3 5.4 3.9 5.6 5.0 4.6 4.2
##  [577] 4.4 4.7 3.6 5.0 4.8 5.5 4.1 4.5 3.4 2.7 4.9 4.5 6.0 5.2 3.2 4.3 5.5 3.5
##  [595] 3.5 3.4 4.7 5.5 3.5 5.5 4.9 3.0 4.7 4.7 5.5 5.0 4.3 4.4 4.4 4.6 4.1 5.1
##  [613] 5.5 4.7 4.9 3.8 3.1 3.9 5.2 5.2 5.6 5.0 4.6 5.1 5.1 5.8 5.2 4.5 5.6 4.1
##  [631] 4.0 3.9 4.4 4.7 2.9 3.5 5.8 3.2 3.4 3.3 3.6 6.5 4.3 3.3 4.7 4.0 3.8 5.4
##  [649] 3.8 3.8 3.9 4.1 5.7 4.9 4.4 4.7 5.5 4.7 7.0 2.9 4.8 4.4 3.9 4.2 4.6 4.8
##  [667] 3.4 4.4 4.8 5.1 2.8 3.6 4.3 4.2 5.7 4.9 6.0 3.8 4.5 4.7 6.0 4.7 4.9 5.1
##  [685] 4.7 2.8 4.4 4.3 4.0 3.2 4.8 3.7 4.3 4.5 5.2 5.7 2.3 4.3 4.1 3.8 4.6 7.2
##  [703] 3.3 3.6 4.6 6.3 5.1 3.5 5.6 5.8 4.0 3.7 3.0 5.1 4.5 4.4 4.1 4.9 5.0 4.1
##  [721] 3.6 4.5 6.4 3.8 3.9 3.7 4.8 3.7 4.8 2.6 5.0 5.0 4.8 4.5 4.7 4.3 3.6 4.1
##  [739] 4.9 3.7 4.7 4.7 4.6 4.1 3.9 4.8 4.3 5.7 5.4 5.1 4.2 4.2 5.4 4.9 4.0 4.2
##  [757] 4.5 3.3 3.6 6.0 4.3 4.9 4.4 4.6 5.4 5.1 4.1 4.7 4.5 3.9 6.2 3.5 4.6 4.3
##  [775] 4.6 4.7 6.9 6.3 4.4 5.6 2.8 3.5 4.0 4.2 4.4 4.7 4.5 5.5 6.2 3.5 4.6 4.5
##  [793] 4.9 2.8 4.1 5.1 4.9 4.6 5.9 4.2 4.0 6.9 6.0 4.7 3.5 4.2 6.4 5.8 4.5 5.2
##  [811] 5.7 3.8 3.7 4.6 4.2 3.4 5.3 5.0 3.0 5.2 5.6 3.7 5.2 6.3 5.1 3.5 5.4 6.7
##  [829] 5.3 4.2 4.8 3.6 3.8 3.9 3.3 5.2 3.9 5.1 3.8 5.2 3.5 4.1 3.8 4.4 4.1 4.9
##  [847] 3.4 4.1 4.0 3.2 4.7 4.8 4.2 4.8 4.7 3.6 6.1 3.8 4.2 6.2 4.7 5.3 3.8 4.4
##  [865] 4.7 4.5 4.2 2.7 4.6 3.7 3.9 5.1 5.6 5.3 5.3 3.9 5.6 5.0 4.9 3.9 3.2 4.1
##  [883] 4.1 5.1 4.0 5.6 4.7 5.1 4.1 5.9 4.5 4.6 4.6 2.7 4.4 5.7 5.3 3.8 4.4 5.1
##  [901] 4.9 5.3 3.9 4.3 3.6 6.7 4.7 3.9 4.1 4.7 5.1 5.1 4.7 3.4 5.6 4.9 4.5 3.4
##  [919] 2.8 6.2 4.5 2.5 3.5 3.9 4.9 4.2 4.8 2.1 3.1 4.5 5.5 5.4 3.9 3.7 3.7 3.3
##  [937] 3.4 3.8 4.3 3.2 4.7 5.0 4.1 5.8 4.4 4.5 5.5 3.0 5.5 6.2 6.1 3.4 4.3 3.3
##  [955] 6.4 5.8 4.7 5.3 3.9 3.0 5.1 3.9 3.5 4.8 4.0 5.9 3.4 5.2 4.7 4.5 3.3 2.1
##  [973] 3.5 5.8 4.9 4.2 5.1 4.7 5.8 4.3 6.0 4.7 5.5 2.3 5.8 5.4 3.4 4.5 3.4 6.3
##  [991] 4.0 3.9 4.7 3.6 3.4 4.3 5.0 4.7 3.5 3.3
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
##    [1] 3.7 2.6 4.9 4.3 4.1 3.7 2.4 4.7 3.6 5.6 3.9 5.7 4.8 4.5 3.0 4.8 3.5 5.0
##   [19] 4.6 4.6 5.6 3.1 4.7 5.5 4.9 2.2 4.2 5.1 4.4 2.3 1.8 3.2 4.9 4.2 3.5 2.8
##   [37] 3.4 3.4 5.1 3.0 3.9 3.3 4.8 5.9 4.5 1.7 6.3 4.2 5.0 4.6 4.4 5.4 3.1 6.6
##   [55] 3.8 3.4 3.7 5.5 3.9 3.6 4.4 2.4 5.2 3.7 5.6 5.4 3.9 3.7 5.2 5.1 3.4 3.2
##   [73] 6.0 4.9 2.9 6.1 4.6 4.7 4.0 4.6 4.7 4.2 3.5 4.8 4.4 4.7 2.9 4.8 5.4 4.2
##   [91] 5.8 5.1 4.7 4.6 4.9 5.2 4.3 5.1 3.4 3.3 5.2 3.8 5.1 3.5 5.1 5.4 5.0 5.0
##  [109] 5.4 6.1 4.1 4.3 5.4 4.5 5.3 5.0 4.1 3.4 3.7 4.8 5.0 3.6 6.1 5.4 5.7 4.4
##  [127] 3.7 4.2 4.5 3.5 5.5 3.9 4.7 4.6 3.3 4.2 4.4 4.2 5.1 4.9 5.2 5.8 5.6 5.0
##  [145] 4.7 4.8 6.5 3.3 4.9 4.2 5.2 2.7 5.6 5.1 4.7 4.2 4.2 5.6 5.4 4.0 3.6 4.3
##  [163] 5.1 4.6 4.7 3.5 4.8 5.3 5.6 5.9 4.4 5.8 5.4 5.1 3.3 4.8 2.7 4.6 3.9 5.5
##  [181] 4.3 4.2 3.4 5.3 4.5 5.4 4.1 4.3 6.2 4.5 3.5 5.4 6.0 5.3 4.2 3.5 6.1 5.2
##  [199] 4.9 4.2 4.1 3.2 3.9 4.2 4.8 3.6 3.7 5.9 2.8 4.5 3.1 3.6 3.9 2.8 5.4 3.9
##  [217] 4.7 3.5 4.2 4.8 4.7 3.8 5.4 4.1 4.0 4.9 4.5 3.7 3.9 4.5 4.3 4.7 2.0 4.5
##  [235] 3.8 4.9 4.5 5.3 3.9 3.4 4.4 5.2 6.0 4.6 5.4 3.4 5.7 3.9 3.0 4.6 4.5 3.1
##  [253] 4.9 5.1 4.1 3.3 5.6 5.1 6.6 4.5 5.9 4.8 3.4 5.3 4.1 4.1 6.4 3.7 4.4 4.8
##  [271] 3.8 5.1 5.1 5.1 4.2 4.9 5.2 4.8 4.4 3.9 4.3 1.9 4.1 5.3 3.8 4.0 3.8 6.2
##  [289] 2.9 4.2 6.0 4.9 5.1 3.3 4.5 2.7 3.1 3.8 4.9 3.0 4.8 6.3 6.5 3.4 4.9 5.8
##  [307] 3.8 5.1 7.2 5.5 4.3 4.1 3.5 2.7 4.4 4.8 6.0 4.9 4.9 5.3 4.2 5.0 4.4 3.7
##  [325] 3.7 4.5 4.1 4.6 3.8 3.9 5.5 5.5 6.0 5.6 4.5 4.7 4.6 5.5 3.6 3.5 3.5 5.7
##  [343] 4.2 4.1 5.4 4.9 4.7 4.5 4.3 5.3 3.0 3.4 4.2 4.0 4.5 4.7 2.8 4.2 3.9 5.7
##  [361] 4.1 6.0 5.0 3.8 5.1 4.9 6.3 3.9 5.3 3.0 6.0 3.3 4.9 4.1 4.8 4.5 6.5 4.7
##  [379] 4.0 3.6 3.3 3.9 4.4 4.2 4.9 3.8 4.7 3.1 4.8 4.7 5.9 4.3 3.9 4.6 6.3 5.0
##  [397] 3.2 4.2 5.3 4.7 4.3 2.6 3.6 5.5 4.9 5.5 3.7 2.7 4.6 4.7 3.5 4.6 3.7 5.7
##  [415] 4.3 4.3 4.2 3.7 5.0 3.8 4.1 4.5 4.5 4.7 5.6 5.5 4.1 5.1 4.3 4.4 2.9 4.9
##  [433] 4.5 5.4 3.8 4.1 4.5 4.5 4.7 5.1 5.5 4.4 4.0 4.2 4.7 4.5 5.4 4.6 6.3 3.1
##  [451] 5.2 3.8 3.7 4.3 7.5 3.1 4.4 5.5 4.7 5.3 5.0 4.9 3.5 3.9 5.3 5.0 4.6 3.1
##  [469] 4.9 5.6 3.5 4.0 5.3 4.1 4.7 4.3 2.1 3.7 5.4 4.9 3.6 3.3 4.4 4.2 5.8 4.2
##  [487] 5.6 5.0 5.2 4.0 4.3 6.6 4.1 3.8 5.9 3.3 3.7 3.8 4.2 3.9 4.1 6.3 4.0 4.7
##  [505] 2.8 3.8 5.6 5.6 5.0 5.6 5.5 2.8 4.7 3.5 3.7 5.6 3.6 4.5 5.1 5.2 4.6 4.4
##  [523] 5.2 2.6 3.1 3.4 4.4 3.8 3.8 6.1 5.5 4.3 6.6 5.3 3.9 4.0 4.2 3.8 4.6 3.4
##  [541] 3.7 5.6 6.4 5.2 5.2 3.4 4.1 2.4 4.0 4.5 3.1 3.9 3.5 3.6 4.1 6.1 2.6 4.5
##  [559] 6.6 3.4 5.0 3.3 5.4 5.3 3.8 5.0 5.1 5.2 3.4 3.6 3.3 3.3 3.9 3.0 5.3 5.4
##  [577] 4.4 4.1 5.4 4.6 4.4 5.3 5.6 4.2 6.1 5.0 3.0 3.9 4.1 4.7 4.2 3.6 2.9 4.0
##  [595] 5.5 4.4 6.5 3.1 3.2 4.7 4.3 4.1 5.3 5.0 3.5 6.3 4.5 5.4 5.2 4.1 5.0 5.6
##  [613] 5.4 5.1 4.8 4.2 5.8 4.2 4.2 4.1 4.0 5.0 5.0 4.0 3.0 3.4 4.3 2.7 5.9 4.9
##  [631] 4.2 4.7 4.2 3.0 6.7 5.5 4.2 5.6 4.8 3.8 6.6 3.9 3.9 3.3 4.6 4.1 5.5 4.9
##  [649] 4.9 4.2 5.0 2.5 3.6 4.0 4.3 5.2 4.3 4.0 4.5 4.3 4.6 6.0 4.5 4.4 6.2 3.0
##  [667] 4.0 4.3 4.4 5.6 6.1 4.7 4.9 4.0 3.8 3.8 3.7 4.4 4.0 4.5 5.5 3.2 5.7 5.1
##  [685] 4.3 4.3 3.3 4.5 4.7 5.5 5.3 4.3 5.3 3.9 5.3 5.1 4.4 4.7 4.4 3.2 3.9 5.5
##  [703] 3.0 3.4 5.8 4.5 5.1 5.0 4.0 4.2 4.9 4.9 4.9 2.8 5.9 5.4 4.4 2.8 4.1 5.9
##  [721] 5.3 6.0 3.0 2.9 4.1 5.1 4.1 5.8 3.0 3.3 4.1 5.1 4.9 3.7 4.7 4.4 5.2 5.0
##  [739] 4.1 4.3 4.0 3.3 4.2 3.5 5.1 5.8 4.9 4.6 5.7 4.0 4.1 3.6 3.5 5.2 4.9 2.9
##  [757] 4.9 4.1 4.0 3.4 3.4 3.7 5.5 4.5 4.6 3.8 5.5 4.9 5.1 3.3 3.4 4.0 4.5 3.1
##  [775] 4.9 4.2 3.9 2.8 5.1 5.0 6.3 4.6 3.6 3.9 4.9 4.8 5.1 4.6 3.9 5.2 5.0 4.3
##  [793] 3.3 5.4 4.4 4.8 3.7 5.3 5.1 3.6 4.3 5.0 3.3 4.4 5.0 2.6 4.1 3.1 5.6 3.1
##  [811] 4.8 4.8 3.9 6.0 4.2 4.1 5.1 3.7 5.8 3.9 4.4 3.9 5.0 3.7 4.2 5.4 4.4 3.3
##  [829] 4.9 5.7 3.2 5.6 3.0 5.1 4.4 4.1 5.0 5.4 4.2 4.4 2.7 4.5 3.9 3.2 4.3 4.7
##  [847] 3.5 5.7 6.4 5.8 3.9 4.4 4.6 4.9 6.8 5.3 4.0 6.2 4.7 4.8 3.9 5.1 6.2 5.3
##  [865] 5.0 4.3 3.9 5.0 4.8 4.7 5.9 5.1 5.7 4.1 4.5 4.7 3.9 4.8 4.8 5.5 4.2 4.1
##  [883] 3.5 3.7 3.6 3.2 5.3 3.6 3.9 2.7 5.2 4.2 4.4 6.7 4.4 5.7 2.9 4.8 4.6 4.5
##  [901] 4.1 4.5 3.6 4.6 4.1 5.2 5.6 2.2 5.4 4.2 3.5 5.0 5.1 6.4 3.5 5.1 5.3 4.2
##  [919] 6.1 5.2 4.2 5.4 4.3 5.0 3.6 4.7 4.6 3.7 4.1 4.6 6.0 4.6 5.1 6.7 6.7 3.2
##  [937] 4.2 4.4 4.3 5.5 3.0 5.4 6.2 5.0 5.8 5.8 4.1 5.1 5.4 4.2 3.7 3.3 5.5 3.8
##  [955] 3.6 5.4 5.0 6.4 4.0 4.5 5.0 4.8 3.9 4.2 5.5 4.3 2.9 5.3 4.7 4.5 4.4 4.7
##  [973] 5.1 4.4 4.0 5.1 6.0 2.4 5.0 5.0 5.9 5.0 5.5 3.1 4.5 3.7 4.5 3.7 3.8 5.1
##  [991] 4.0 5.0 5.4 3.3 3.7 4.0 6.2 4.2 5.9 3.8
## 
## $func.thetastar
## [1] -0.0084
## 
## $jack.boot.val
##  [1]  0.522356495  0.381962865  0.321489971  0.103899721 -0.004359673
##  [6] -0.060795455 -0.196226415 -0.228977273 -0.452840909 -0.527616279
## 
## $jack.boot.se
## [1] 0.9990523
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
##    [1] 4.9 5.1 5.3 5.2 3.6 5.0 4.1 4.8 3.8 3.6 5.8 5.4 4.6 4.0 4.0 3.8 5.5 5.6
##   [19] 6.2 4.0 6.0 5.4 3.7 5.5 4.5 3.3 5.0 3.1 3.9 2.8 6.0 5.4 5.8 4.1 3.7 5.4
##   [37] 3.0 6.0 2.9 3.4 5.0 6.2 5.1 2.7 5.1 3.8 4.6 5.1 4.5 4.7 3.3 4.4 4.7 3.7
##   [55] 4.7 6.8 5.9 6.1 4.3 4.7 5.4 4.4 2.7 4.3 5.1 5.8 5.3 4.3 3.6 5.1 4.1 4.5
##   [73] 3.6 4.4 3.1 4.8 4.6 4.8 5.5 3.9 3.8 4.6 3.9 5.1 4.0 5.0 4.1 4.4 3.6 4.3
##   [91] 5.1 3.9 3.7 4.3 5.3 5.2 3.8 4.3 2.9 6.2 4.7 3.6 3.0 4.8 5.1 4.8 4.0 3.7
##  [109] 2.4 3.7 6.0 4.9 4.3 4.2 4.4 3.9 3.5 4.6 5.7 4.5 5.2 5.3 5.0 4.1 4.5 2.8
##  [127] 4.1 3.7 3.5 4.2 6.5 4.3 5.4 3.6 5.7 4.3 5.5 4.0 5.5 5.3 4.5 4.7 3.3 5.0
##  [145] 3.8 5.6 2.6 4.2 5.9 4.9 5.7 4.2 3.8 5.2 4.0 4.3 4.1 4.4 3.9 4.4 2.7 5.4
##  [163] 3.7 5.3 5.1 5.2 5.1 4.1 3.6 5.2 3.7 3.1 4.1 5.8 4.2 4.6 4.1 6.8 3.8 5.9
##  [181] 4.3 5.2 4.3 4.1 5.2 2.9 4.5 4.8 4.0 4.5 3.2 4.1 4.1 4.9 2.6 4.8 3.3 4.9
##  [199] 4.8 4.6 3.7 3.3 3.8 4.6 4.1 4.5 6.0 4.7 4.4 3.5 3.0 5.2 3.6 4.6 3.8 4.1
##  [217] 4.9 3.9 4.3 4.9 4.5 1.5 6.2 3.4 5.1 3.2 5.7 5.9 4.3 4.6 4.7 3.9 4.2 3.4
##  [235] 6.1 3.9 4.7 3.8 5.6 4.8 4.0 5.2 4.2 3.9 4.7 5.1 3.6 5.1 2.6 3.0 4.4 4.8
##  [253] 3.9 3.2 4.6 3.6 3.4 6.0 4.7 5.0 5.5 4.2 4.4 2.2 5.5 5.9 4.9 3.7 2.7 5.2
##  [271] 5.0 3.8 5.2 4.6 5.7 4.1 6.3 2.3 4.6 3.1 5.1 3.3 4.2 4.0 4.0 5.5 3.6 4.1
##  [289] 4.9 4.1 3.1 5.4 5.1 5.2 4.3 3.7 3.8 2.1 3.7 6.2 4.7 4.2 4.9 4.5 5.3 4.2
##  [307] 3.2 5.2 3.5 4.7 5.0 5.7 4.5 3.5 4.6 4.1 4.8 4.7 2.6 3.1 4.1 3.7 4.8 6.7
##  [325] 4.2 5.2 3.3 3.5 3.1 4.2 4.8 4.9 4.0 4.1 4.5 4.5 5.3 4.5 4.0 4.2 3.7 4.5
##  [343] 4.7 5.1 4.7 4.3 4.1 4.7 4.9 5.8 4.2 5.9 4.5 4.4 5.0 5.3 4.1 2.5 4.9 3.9
##  [361] 4.0 5.1 4.8 4.4 3.8 6.2 4.5 3.7 3.3 5.1 5.0 5.0 4.8 5.4 4.2 4.4 6.2 4.4
##  [379] 5.3 3.5 4.5 5.8 3.5 4.3 5.4 5.8 5.6 4.9 6.9 4.1 4.1 4.1 3.8 4.6 5.3 4.6
##  [397] 5.4 4.9 4.2 3.9 5.7 3.0 5.0 4.6 4.9 4.9 4.9 3.5 2.9 5.5 4.0 4.4 4.7 5.1
##  [415] 4.2 4.4 5.4 4.2 4.6 3.9 4.3 4.7 3.5 4.9 3.9 6.3 5.3 4.4 4.6 3.1 3.4 4.7
##  [433] 6.1 5.0 3.5 3.5 6.0 6.0 4.7 4.1 6.1 3.9 2.4 5.2 4.3 4.0 3.2 4.5 3.7 4.4
##  [451] 6.3 2.1 4.3 2.7 5.9 5.9 4.9 4.1 4.4 4.9 4.8 6.0 5.4 5.8 6.5 4.3 5.6 5.0
##  [469] 5.5 5.2 5.1 2.6 4.0 3.8 6.2 5.8 4.1 2.9 3.5 6.4 4.8 3.4 5.2 4.9 3.8 3.6
##  [487] 5.2 3.9 3.9 4.8 4.2 5.6 4.1 4.9 5.3 3.5 4.8 4.9 4.5 4.5 4.6 3.4 3.8 5.6
##  [505] 4.9 5.5 5.8 3.5 4.5 5.0 4.7 3.2 5.0 3.7 5.5 6.2 4.3 4.6 5.3 3.6 4.7 3.7
##  [523] 4.7 5.4 5.9 3.9 3.7 6.2 3.9 3.8 6.2 3.3 4.3 5.7 3.9 5.7 4.7 4.6 5.5 4.9
##  [541] 3.4 4.7 5.6 3.8 5.0 4.5 5.2 4.7 4.2 3.8 5.6 4.7 3.3 5.5 4.1 4.2 5.2 5.3
##  [559] 3.4 4.9 4.0 4.9 3.7 4.3 5.2 4.1 4.0 3.6 3.7 5.3 4.8 5.0 5.6 2.7 4.2 4.2
##  [577] 6.4 6.0 4.2 3.2 4.7 4.6 5.7 3.9 5.5 5.1 4.3 5.3 3.9 3.6 5.1 4.9 3.5 6.1
##  [595] 4.6 5.0 3.4 6.1 3.6 3.8 4.0 5.3 5.7 5.2 5.0 3.7 5.1 4.0 4.0 4.4 4.8 5.5
##  [613] 4.5 4.4 5.0 5.2 4.3 5.5 5.3 3.5 5.3 5.2 4.7 5.3 3.5 4.8 6.0 4.7 5.1 3.2
##  [631] 5.1 5.8 3.3 4.9 6.1 4.2 4.8 3.6 6.4 3.2 4.1 3.6 4.4 3.7 3.8 4.4 2.9 3.9
##  [649] 5.4 4.3 5.6 3.1 5.7 5.1 5.0 6.5 4.4 3.6 4.0 4.7 3.4 4.2 3.9 4.0 4.5 4.1
##  [667] 5.1 3.7 2.2 4.5 5.6 4.7 5.2 2.7 5.2 5.4 6.7 4.4 4.2 3.4 6.8 4.7 3.1 5.0
##  [685] 3.6 4.7 3.7 4.2 4.3 2.7 2.7 4.7 5.0 5.7 3.9 3.8 2.9 5.2 5.0 3.6 4.7 5.8
##  [703] 3.1 5.3 6.6 4.5 4.5 4.3 5.2 3.7 3.9 5.1 3.1 5.1 4.2 5.4 3.5 2.2 4.0 6.1
##  [721] 5.9 3.5 4.3 4.5 4.0 4.8 4.4 5.1 4.2 3.8 4.0 4.5 5.1 5.0 6.2 6.7 4.4 3.7
##  [739] 2.8 4.2 5.2 3.3 4.1 3.5 3.9 4.6 3.5 4.6 5.1 3.2 2.0 4.3 6.2 4.0 5.7 5.2
##  [757] 5.7 6.6 3.7 2.2 2.5 3.8 5.0 4.2 3.5 5.1 5.2 3.8 5.0 4.6 5.0 4.5 5.4 5.0
##  [775] 3.3 6.3 2.8 4.8 5.0 4.1 3.8 4.4 4.0 5.5 4.3 4.7 4.9 6.0 4.3 4.2 3.1 4.0
##  [793] 5.3 5.9 4.7 5.7 6.4 5.9 3.8 4.5 4.6 4.4 7.0 4.6 5.1 5.2 4.8 4.4 4.1 4.4
##  [811] 5.4 4.7 4.8 4.4 3.7 4.2 4.7 5.1 4.0 4.5 4.9 4.6 3.6 3.6 3.3 3.9 3.0 3.6
##  [829] 4.8 4.5 5.1 4.8 5.6 3.6 5.1 5.4 3.3 5.9 6.8 4.2 4.6 4.2 4.0 4.9 6.5 5.5
##  [847] 6.0 3.8 5.3 4.3 5.4 4.5 4.1 5.1 4.8 2.9 5.6 3.5 4.5 4.0 3.8 5.3 3.3 4.8
##  [865] 3.7 5.1 2.6 4.1 4.3 2.7 4.6 5.4 5.4 2.9 3.2 4.0 4.1 3.9 4.0 4.1 6.2 3.8
##  [883] 4.3 4.4 2.3 5.9 5.1 4.4 4.7 3.5 3.7 2.8 5.9 4.4 6.2 5.1 5.4 3.5 3.8 5.2
##  [901] 3.3 3.4 3.3 4.9 4.0 5.0 3.8 5.5 7.1 2.9 2.9 2.4 5.7 3.3 3.5 5.8 4.3 4.6
##  [919] 4.6 5.3 5.0 4.0 3.2 4.3 3.4 4.6 4.4 5.8 4.2 5.6 4.5 3.8 5.2 3.5 5.5 4.0
##  [937] 5.4 3.5 3.3 4.8 3.2 5.3 4.1 3.5 5.3 5.3 3.8 4.2 6.0 4.3 5.5 4.5 3.7 3.9
##  [955] 3.2 2.9 3.9 3.7 5.3 4.9 3.7 6.5 4.3 5.4 4.0 4.8 3.5 3.7 4.3 4.6 2.9 5.6
##  [973] 5.1 5.8 3.6 5.1 5.5 3.3 5.3 3.5 5.0 5.3 5.9 5.1 4.5 4.5 3.7 4.1 4.2 3.8
##  [991] 5.8 4.7 4.3 3.9 4.5 3.9 3.8 5.8 5.0 3.9
## 
## $func.thetastar
## 72% 
## 5.1 
## 
## $jack.boot.val
##  [1] 5.400 5.500 5.300 5.200 5.100 5.100 5.000 4.800 4.600 4.504
## 
## $jack.boot.se
## [1] 0.9394429
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
## [1] 0.8906098
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
##   11.192328   14.905943 
##  ( 4.932563) ( 6.718586)
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
## [1]  0.08714117  0.25708799 -0.76685139  0.44085688 -0.07230625 -0.52636009
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
##    [1]  0.163713276  0.499666102  1.554590047  0.559996887  0.897132931
##    [6]  0.200945427  1.536501628  0.814794659  0.412182920  0.498659234
##   [11]  0.497873383  0.736727402  0.236324500  0.141093844  0.147940263
##   [16] -0.145110320  0.452407232  1.901409657  0.819395368  0.960856224
##   [21]  1.401221931 -0.164602967  0.379938794  1.179919923  0.819149309
##   [26]  0.508402008  0.638336121  0.210452219  1.873097175 -0.301583188
##   [31]  0.447960146  2.095577050  0.727805758  0.498647772  0.962845659
##   [36]  0.110962973  1.046004369  0.937516387  1.019823652  0.806767346
##   [41]  0.609092323  0.845422675  1.554185679  0.780318770  1.542542921
##   [46]  1.361383153  1.604046227  1.827227802  1.089631327  1.526637068
##   [51] -0.117302494  0.443844927  1.113626854  0.558566846  1.050691114
##   [56]  1.219567231  0.652888568  1.598637947  0.123684931  0.929915701
##   [61]  0.099763587  0.585011715  0.537692746  0.618056907  1.377950777
##   [66]  1.113849663  0.615877454  0.629683190  0.861814685  0.223420961
##   [71]  0.417623081  0.617600447  0.516290686 -0.082649016  0.847559941
##   [76]  0.990688634  1.298619880  0.274383346  0.009063742  1.329963038
##   [81]  0.282619709  0.701447179  0.055839967  1.799934828  0.645379699
##   [86]  0.554803569  0.100543688  0.907406381  0.628829631  0.684549505
##   [91]  0.435233795  0.702563144  1.779485035  0.202309749  0.635816093
##   [96]  0.914842380  0.346927736  0.724669322 -0.160977519  1.278353887
##  [101]  0.731020281  0.205435002 -0.230829019  0.638923035 -0.120796783
##  [106]  0.887302429  0.899015018  1.212080105  0.689351851  0.393399076
##  [111]  0.436070198  0.606322479  2.027587980  0.716887429 -0.290394572
##  [116]  0.896435362 -0.279501828  0.681623386  0.335855645  0.362299571
##  [121]  0.662283340  0.622654433  1.815094777  0.952982085  0.345106779
##  [126] -0.008118390  1.224452104  1.097136030  0.146683537  0.592528991
##  [131]  0.478201767  1.938624888  1.009617480  0.699601739  0.085532909
##  [136]  0.413975616  0.401499696  0.816785965  1.616210863  0.249680969
##  [141]  1.108354440  0.842981282  0.776762856  0.428997857  0.502004534
##  [146] -0.060079769  0.857957968  2.170402727  1.215150051  0.971553943
##  [151]  1.475895548  1.217770689  1.345815253  0.137122440  0.653734198
##  [156]  1.062986277  0.528518120 -0.265857700  0.005712403 -0.538234193
##  [161]  1.000692434  0.670756622  0.554765786  0.598682291  0.275976272
##  [166]  0.655396264 -0.036981431  0.057959407  0.918831735  1.136715407
##  [171]  0.727772096  0.846643950  0.750114759  1.239977873  0.474889020
##  [176]  1.155441127  0.592980837  1.834994760  0.988450777  0.670154325
##  [181]  0.404638128  0.747682079  0.258661828  0.509459275 -0.199153064
##  [186]  0.240671024  1.394765590  0.688407462  0.650633108  1.056530817
##  [191]  1.277802332  0.004009465  1.057957185  1.753015588  0.720407087
##  [196]  1.852539729  0.466215369  0.517041849  0.988808656  0.062785462
##  [201]  0.503212059  0.578020530  1.282372327  0.807246397  0.848707425
##  [206]  0.365511807  0.627253855  0.184712685  0.915755432 -0.155933867
##  [211]  1.189900767  1.099222591  0.706106721  0.916977505  0.740093122
##  [216]  0.181509044  0.397841789  1.317175233  0.975944417  0.734184224
##  [221] -0.407681694  1.103559728  0.778238288  0.854418240  1.155417274
##  [226]  0.679491682  0.975956253  0.441068980  0.935489715  0.534186116
##  [231]  0.931162787  0.137775919  1.068648966  0.391348643 -0.807915507
##  [236]  0.846001925  0.271709477  0.156648948  0.251343603  1.070697180
##  [241]  1.062860107  0.636486125  0.684498615 -0.792995112  0.279935229
##  [246]  0.776000487  0.928836835  0.838158040  0.202307256  0.501124703
##  [251]  0.128192832  1.184787256  1.184229896  0.241045916  0.944020066
##  [256] -0.099291825  1.305295361  0.819984363  0.279449153  0.722081218
##  [261]  1.345732771  0.829275296  0.038081622  0.806034595 -0.176222167
##  [266]  0.401057809  1.879996153  0.323205825 -0.452508609  0.931677965
##  [271]  0.650061960  0.446985862  0.373060651  0.925922525  0.914872228
##  [276]  0.408444081  0.145238699 -0.419751518  0.508110722  1.002666963
##  [281]  0.302858501  0.727396754  0.762239969  0.932192031  0.413600497
##  [286]  1.475186593  0.951038310  0.132268844  0.979479534  1.430059766
##  [291]  0.776393321 -0.247113875 -0.042360481  0.391461890  1.449043524
##  [296]  0.326369242  0.061184029 -0.067621460  1.176450297  0.088181606
##  [301]  1.111491396  0.845273429  1.796915766  1.039887448  0.438148176
##  [306]  0.441817312 -0.281352062  0.688323644  0.663932847  1.032826082
##  [311]  0.355136317  1.023854865  0.937145736  0.129770743 -0.077282619
##  [316]  1.193826627  0.156648948  0.798188042  1.687855821  0.557321322
##  [321]  0.493298376  0.051178981  0.737907559 -0.054583410  1.185407110
##  [326]  0.398769734  0.772114979  0.667882213  1.135185380  1.634371761
##  [331]  0.364613363  1.576935779  0.521882073 -0.232644529  0.801317758
##  [336]  0.347082839  1.183670155  1.026077396  0.796130571  0.981649029
##  [341]  1.188685596  0.637385291  0.946846178  0.770708580  0.845828170
##  [346]  0.978667915  0.554737144  1.044519791  0.557240098 -0.283083793
##  [351]  0.994747715  1.202679622  0.501124703  0.924047140  1.048908129
##  [356]  0.532991509 -0.247113875  1.637038233  0.138654434  0.830571840
##  [361]  1.097662968 -0.380192639  1.243005135  0.555544759  1.611258454
##  [366]  0.302797072  0.896685081  0.983433228  0.732818724  0.325836582
##  [371]  0.647512046  0.121448728  1.197515934  0.837550261  0.746741385
##  [376]  0.908162180  0.745643286  0.557790036  0.824166427  0.600554557
##  [381] -0.154829525  0.653229000  1.344251766 -0.354572137  0.239793557
##  [386]  0.336073247  1.068994956  0.282976078  0.856459225  2.135848259
##  [391]  0.795340352  0.252065527  0.714598096  1.113183853  0.781873683
##  [396]  1.141915561  0.159991873  1.273550699  1.394403304  0.376738613
##  [401]  0.440006649  0.107518406  0.459930677  0.855081432  0.704614823
##  [406]  0.849996862  0.435165201  0.963561633 -0.109893752  1.417562397
##  [411]  0.814342154  0.675173380  0.433424042  0.700382281  0.552488370
##  [416]  1.274220497  0.884870738  0.890623803  0.830051676  0.699720387
##  [421]  1.664047634  0.998935233  0.960272189  0.936760416  0.960728730
##  [426]  1.926893779  0.463341232  0.860065291 -0.021225674  0.569763668
##  [431]  1.752501946  2.082013208  0.710441660  0.303725917  0.620458864
##  [436]  0.921548557  0.826872942  0.884298368  1.154199536  0.576276647
##  [441]  1.002714293  0.003497678  1.731005291  0.484351612  0.394753569
##  [446]  2.131714476  1.398021806  0.439698824  0.481211730  0.409531179
##  [451]  0.852120009  1.077290706 -0.072220737  0.267899089  1.019769549
##  [456]  1.045888717  1.293328335  1.018041619  0.213863902  0.738778171
##  [461]  1.793878349  0.586514214  0.590676971  0.712997700  1.111919237
##  [466] -0.125310127 -0.245257806  0.355385308  0.439103425  0.723384031
##  [471] -0.007050643  0.751128922 -0.064961536  1.091686643  0.413642765
##  [476]  0.032331039  1.452210453  0.393319249  1.101808224  0.891186032
##  [481]  1.512242651  0.592366331  1.931435915  0.131552079  1.449055109
##  [486]  0.464746018  0.671817858  0.297206541  0.463343496  0.274578136
##  [491]  0.192249202  0.210464667  1.087195512  0.860253171  0.660828426
##  [496]  0.717670005  0.410989466  0.811943298  1.008075752  0.326625684
##  [501]  0.482317851  0.124092743  0.694055704  1.024591675  0.106102621
##  [506]  1.572514506  0.509278224 -0.494808059  0.920697085  0.563839155
##  [511]  0.440348524  0.073275955  0.505126293  0.183294981  1.570641258
##  [516]  1.566448254  0.881426816  1.128841355  0.399827201  0.487629046
##  [521]  0.241242179  0.722014384  0.650307362  1.171422808  0.705744557
##  [526] -0.120007915  1.023976302  0.765857136  0.145860933  1.357988489
##  [531]  0.876812735  0.671612803  0.843241779  0.020964570  1.535274909
##  [536]  0.895323309  0.283493084  0.570190863  0.381672159  0.445095929
##  [541]  0.611276775  0.898122891  1.318816166  1.126580805  0.506700017
##  [546]  0.744434303  0.669997197  0.374289352  1.379184139  0.457329190
##  [551]  0.881398112  0.520871944  0.613869002  0.598902117  0.840185906
##  [556]  0.756571640  0.721591434  0.916873321  0.741200762  0.314851613
##  [561]  0.475767515  0.427692469  1.118469739  0.336581992  0.633729691
##  [566]  0.642126128  0.637805682  0.819501962  0.944795258 -0.285930142
##  [571]  1.086172500  0.560104517  0.509660073  1.352586081  0.141502027
##  [576]  0.285188575  0.941532875  0.173183635  0.602517905  0.587251413
##  [581]  0.501964459  1.061243152 -0.298141788  0.181509044 -0.058472911
##  [586]  0.882217639  1.257793770  0.733774105  0.381014657  1.586660879
##  [591]  0.717225816  0.921069806 -0.028228356  1.414796683  0.204106491
##  [596]  2.327069662  0.722557977  0.644286886  0.157286314  2.174382367
##  [601]  0.665687606  0.332175196  0.494160530  0.002492425  0.200118914
##  [606]  0.629683190  0.128020560  0.439675545  1.116597315  0.712759227
##  [611]  0.873329347 -0.123835255  1.289240405  0.564904612  0.626650364
##  [616]  0.840387635  0.575509769  0.744463930  0.307467474  0.860947462
##  [621]  0.208574532  1.492991763  0.194256455  0.194256455  0.635713175
##  [626]  1.220322469  0.781649736  0.396372235  1.668819950  1.025069001
##  [631]  1.061150041  0.757081359  0.779123683  0.038081622  1.157953986
##  [636] -0.066872435  0.335442188  0.222272962  0.877630293  0.650421165
##  [641]  0.850547727  1.141018049  0.389848755  1.421885730  1.143861789
##  [646]  0.181138901 -0.002800687  0.825372818  0.983669363  0.905061753
##  [651]  0.393488871  0.382742178  0.893589663  0.110054546  0.431168394
##  [656] -0.322579912  0.459930677  0.221264618  1.383116601  0.629258705
##  [661]  1.032818726  1.201933464  1.680300155  1.007243669  1.060726818
##  [666]  0.499556990  0.803111204  0.119440237  1.112879448  0.849544665
##  [671] -0.032832470  0.745445534  0.290364447  0.869025037  1.275608577
##  [676]  0.620325876  0.478678270  0.651778499  0.858490937 -0.012138416
##  [681]  0.722014384  0.219432926  0.997745216  0.926986543  0.397244268
##  [686]  1.173328389 -0.078460961  0.755426561  0.563693313  1.095951546
##  [691]  0.685511472  0.104549154  0.483918635  0.673823041 -0.064888683
##  [696]  0.861263461  0.073739265  1.029042002  0.373829544  0.365644046
##  [701]  1.108010671  1.123294802  0.200208803  0.720124447  0.251775558
##  [706]  0.569161952  1.271087835  0.119269239 -0.029567316  1.269636221
##  [711]  0.848636937 -0.237606859  0.619189699  0.794225256  0.699888667
##  [716]  0.351009422  0.528078720  0.296451928  1.697657499  0.415131904
##  [721]  1.318521911  0.716473318  0.582353369  0.630709300  1.370644622
##  [726]  0.358015314 -0.270722190  0.679758701  1.412396135  0.419157466
##  [731]  1.253500596  0.557298761  0.160560763  1.966284847  0.476031116
##  [736]  0.406063052  0.874538955  1.626563148  1.373201414 -0.011884340
##  [741]  0.801563143  0.341156869  1.163099653  0.733370272  1.071016664
##  [746]  0.886811045  1.257152309  0.447658728  0.990935310 -0.350324151
##  [751]  0.094022196  1.087857478  0.889570558  0.829953547  1.039655003
##  [756]  1.598161111  1.771274735  0.988768325  1.037374341  0.450247209
##  [761]  0.843156811  0.766769214  0.364699396  2.284364302  0.935489715
##  [766]  1.838803604  1.383135685  0.603945664  1.210773661  1.350277153
##  [771] -0.295613185  0.810523987  1.199845038  0.125648191  1.266306980
##  [776]  0.744810713  1.071889828  0.436304051  0.472512463  0.652681342
##  [781]  0.361070151  0.387788841  0.405791735  0.627316115  0.349862804
##  [786]  0.792280179  0.946222527  0.142022487  1.057983643  0.406994958
##  [791]  0.952588275  1.139483680  0.513465652  0.177807387  0.151557642
##  [796]  0.397668431  1.498860979  0.460529752  0.644933059  0.191441610
##  [801]  0.526428451  0.975027322  0.094173266  1.429131594  1.191880292
##  [806]  0.889569474  1.554675289  1.604046227  0.597674411  0.710272275
##  [811]  0.790305280  1.203029670 -0.097960846 -0.111617711  1.328609183
##  [816]  0.746588926  0.947013692  0.301296384  0.576307712  1.313276507
##  [821]  0.615134689  0.645243998  1.133020562  0.641637108  0.424253545
##  [826]  0.400855902  0.071191550  0.292606694  1.745588262  1.389675748
##  [831]  0.640658579 -0.011374023  0.617162530  0.676671486  1.509846051
##  [836]  1.058160014  0.909191548  1.439508804  0.796722718  0.692775626
##  [841] -0.094487170 -0.084965160  0.748779572  0.668165783  1.667736178
##  [846]  2.028647819  1.461782965  0.074393933  1.273394589 -0.574328656
##  [851] -0.043957386  1.112145715  0.158912418  0.807494120  1.375171249
##  [856]  1.220939712  0.340746726  1.560327256  1.015422720  0.805458611
##  [861] -0.041995909  0.545536584  1.108010671  0.147707857  1.410075288
##  [866]  0.648964858  0.835564726  0.670400067  1.005795570  0.062373098
##  [871]  0.756503402  1.542542921  0.988375611  0.028108444  1.593435432
##  [876]  0.258233081  0.735267791  1.595673513  0.468943409  1.985897395
##  [881]  0.651156043  1.145663841  0.904930866  0.327764193  1.033861768
##  [886] -0.049582760  0.636774677  0.881398112  0.246025170  0.937426648
##  [891]  0.904407266  0.238444412  0.710559402  1.139161319  0.304386953
##  [896]  0.677478297  0.983373596  0.232450261  0.537326572  0.840029933
##  [901]  1.455427637  0.248999284  1.188682387  0.430584613  0.306147657
##  [906] -0.487938034  0.205748376  0.510054994  1.459060493  1.227281525
##  [911]  0.353095765  0.726715726  0.883531274  0.125985764  1.857553447
##  [916]  0.928060513  0.662567258  1.037065422  0.379588231  0.421443416
##  [921]  0.049319820  1.101501383  0.990688634  0.294793995  0.412162852
##  [926]  0.570085387  0.474152424  0.478953038  0.746017518  0.078743421
##  [931]  0.037670954  0.505071853  0.516602478  0.356433127 -0.190762654
##  [936]  0.862745761  1.191781622  0.773739868  0.794856696  0.477109135
##  [941]  0.634926369  0.418010352  1.016906765  0.893589663  0.577965994
##  [946]  0.388528105  0.649399922  0.983018648  1.158012849  1.560310789
##  [951]  1.940698108  0.496818967 -0.057911173  0.836520768  0.533652001
##  [956]  0.767188063  0.598606809  1.017440712  1.352351978  0.624016927
##  [961]  0.575405385  0.187895366  1.162492333  0.269151572  1.059444694
##  [966]  0.636858148 -0.321076750  1.876572539  0.908906085  1.751330256
##  [971] -0.537382892  0.554528824  0.810256156  0.277791367  0.563450961
##  [976]  0.151238677  1.545112808  0.215785583  1.621429402  1.456987526
##  [981]  0.740994703  0.215366915  0.817694846  0.712629309  0.445595600
##  [986]  0.009256451  0.504148631  0.839020652 -0.094756554  1.059614648
##  [991]  0.952892139  0.537824690  0.742829719  2.152267850  0.158289939
##  [996]  1.058786130  0.528451954  0.923252076  0.385405711  0.016949928
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
##   0.75086241   0.23611388 
##  (0.07466577) (0.05279359)
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
## [1] -0.26442577 -0.02492872 -0.24124073 -0.30040588 -0.40713280 -0.21812159
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
## [1] 0.0494
```

```r
sd(results2$thetastar)  #the standard deviation of the theta stars is the SE of the statistic (in this case, the mean)
```

```
## [1] 0.925806
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
## t1*      4.5 0.01381381   0.9164109
```

One of the main advantages to the 'boot' package over the 'bootstrap' package is the nicer formatting of the output.

Going back to our original code, lets see how we could reproduce all of these numbers:


```r
table(sample(x,size=length(x),replace=T))
```

```
## 
## 1 2 3 4 8 9 
## 2 1 3 1 1 2
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
## [1] -0.0154
```

```r
se.boot
```

```
## [1] 0.901907
```

Why do our numbers not agree exactly with those of the boot package? This is because our estimates of bias and standard error are just estimates, and they carry with them their own uncertainties. That is one of the reasons we might bother doing jackknife-after-bootstrap.

The 'boot' package has a LOT of functionality. If we have time, we will come back to some of these more complex functions later in the semester as we cover topics like regression and glm.

