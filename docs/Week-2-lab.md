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
## 2 3 4 5 7 8 9 
## 1 1 1 1 2 2 2
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
## [1] 0.0195
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
## [1] 2.734588
```

```r
UL.boot
```

```
## [1] 6.304412
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
##    [1] 5.7 3.9 4.4 2.7 3.9 3.2 5.1 4.7 5.5 5.2 4.6 5.1 4.5 3.3 4.3 4.4 4.8 4.6
##   [19] 5.3 6.3 6.4 4.3 3.8 7.0 5.6 6.6 5.7 3.8 2.6 4.7 5.9 3.6 3.1 3.5 3.7 5.8
##   [37] 4.1 4.0 4.4 3.6 4.9 5.8 4.8 4.9 4.9 3.7 4.5 4.9 4.4 4.9 2.5 4.0 4.5 6.7
##   [55] 4.0 4.3 3.4 4.4 5.3 5.1 4.3 5.4 5.7 4.3 3.6 5.2 6.0 4.9 5.3 4.5 4.2 4.3
##   [73] 5.1 5.3 3.8 4.6 4.4 4.1 4.5 4.6 3.7 4.1 3.7 5.7 5.0 4.1 5.8 4.7 4.4 5.5
##   [91] 4.2 2.7 4.7 3.1 5.5 2.4 4.7 5.2 4.4 4.2 6.1 4.7 5.6 4.6 5.6 4.5 3.8 3.5
##  [109] 6.2 4.5 3.8 5.5 2.7 4.0 3.6 4.6 5.2 3.9 4.3 3.9 3.0 4.7 3.3 5.4 5.6 3.4
##  [127] 3.8 2.4 3.5 5.0 4.6 5.3 4.0 5.3 2.8 5.8 4.6 5.9 4.5 4.7 4.5 3.7 4.8 5.7
##  [145] 3.1 3.8 4.0 4.6 4.1 5.8 4.1 4.1 4.5 5.1 4.4 3.8 4.9 4.0 4.4 3.7 3.1 4.7
##  [163] 5.4 4.6 5.9 3.2 5.1 5.0 4.9 4.1 4.2 5.4 4.9 4.0 3.3 4.7 4.2 3.2 3.8 4.8
##  [181] 5.3 4.1 6.3 5.2 3.8 3.9 3.9 5.5 4.3 4.7 4.4 4.3 3.8 4.9 3.6 4.5 5.0 4.2
##  [199] 6.5 3.2 3.6 2.6 3.5 4.8 3.0 3.7 4.6 4.6 4.9 4.6 4.4 4.1 3.9 4.6 6.1 5.5
##  [217] 4.0 4.3 4.4 3.0 4.8 4.6 4.1 5.5 5.1 3.8 5.4 4.3 2.6 4.2 4.5 4.9 4.0 3.6
##  [235] 3.9 3.9 4.2 4.8 4.0 4.4 5.2 3.3 4.8 3.9 5.0 2.5 4.8 3.6 5.2 3.5 3.1 3.3
##  [253] 5.2 3.8 5.0 4.1 5.2 4.5 3.5 4.5 3.4 6.2 3.7 3.7 4.9 5.5 5.0 5.6 3.7 4.5
##  [271] 4.7 6.8 4.5 5.3 4.1 5.7 2.7 4.6 5.6 3.6 3.1 4.5 4.0 4.8 5.2 4.5 2.6 3.6
##  [289] 4.0 4.8 4.9 3.4 4.7 4.3 2.7 3.8 5.5 5.1 2.1 2.8 4.9 4.9 4.4 3.8 3.7 5.9
##  [307] 4.6 4.5 5.9 5.7 5.2 3.7 4.7 3.4 5.1 4.7 4.0 5.6 2.9 4.1 5.2 4.7 4.9 3.5
##  [325] 5.3 4.1 5.1 5.0 4.6 4.7 5.7 3.8 3.2 6.1 4.2 5.8 3.4 3.7 5.0 4.5 3.7 4.5
##  [343] 3.6 3.9 5.0 5.0 3.1 4.7 4.3 2.7 4.7 4.7 5.7 5.3 5.0 3.8 4.6 5.3 3.5 4.1
##  [361] 3.8 4.2 5.5 3.2 4.6 5.4 3.2 5.7 5.4 4.3 5.0 3.2 5.2 3.9 5.5 4.6 5.7 2.9
##  [379] 3.8 3.7 4.0 3.0 4.7 4.4 5.7 4.9 5.2 3.8 5.8 5.7 5.2 4.6 5.1 5.3 4.4 4.8
##  [397] 4.0 5.5 4.5 4.9 4.5 4.4 3.9 4.2 3.3 5.5 4.0 4.7 3.4 4.5 4.1 5.1 4.7 1.5
##  [415] 3.8 5.3 4.5 3.6 4.4 6.8 3.8 5.0 3.7 3.1 2.8 3.2 4.6 3.7 5.2 4.1 4.5 4.7
##  [433] 4.4 4.7 4.2 4.6 5.5 3.8 4.2 5.2 4.5 5.4 4.7 3.2 6.2 5.3 4.1 4.4 3.4 5.8
##  [451] 5.3 5.0 3.5 3.4 6.1 5.5 4.8 4.4 3.5 5.1 5.0 5.9 5.1 5.2 4.6 5.4 4.0 4.8
##  [469] 6.7 4.1 4.9 4.6 4.7 4.6 3.7 4.9 4.3 4.0 6.1 4.5 4.8 4.9 6.0 4.4 3.9 3.3
##  [487] 5.4 4.9 4.7 3.4 3.5 5.5 5.6 4.5 6.2 4.5 6.2 4.4 5.1 5.2 4.2 4.1 4.6 3.4
##  [505] 6.0 5.7 3.3 4.9 4.6 5.0 4.8 5.5 3.7 4.3 4.9 3.9 4.2 4.1 3.6 3.7 3.3 4.4
##  [523] 4.4 6.5 3.0 5.4 3.5 4.9 4.8 4.2 5.4 5.7 3.5 4.7 4.1 5.7 6.2 2.4 3.3 4.7
##  [541] 4.4 2.9 4.7 3.0 3.6 4.1 5.2 5.3 3.8 2.9 5.0 3.1 4.9 5.0 5.7 2.9 3.7 4.8
##  [559] 5.7 6.4 4.3 4.1 4.7 5.8 3.3 4.5 4.4 3.8 2.7 5.3 5.0 5.1 4.7 5.1 3.9 3.6
##  [577] 3.6 5.5 3.1 3.6 6.0 4.0 5.8 5.9 4.3 5.7 4.9 2.7 4.0 3.7 4.9 3.8 5.3 5.2
##  [595] 5.6 4.1 4.0 5.0 4.3 2.7 3.4 5.8 4.8 4.4 5.8 4.7 5.1 6.3 6.3 5.5 5.1 5.2
##  [613] 3.8 5.9 5.4 3.1 4.9 5.7 4.5 4.8 3.8 3.2 5.2 4.7 2.7 4.7 4.5 6.3 4.0 4.8
##  [631] 5.0 2.8 3.8 4.0 5.5 5.4 5.9 4.5 5.3 5.3 5.3 3.5 6.5 6.2 2.6 5.3 4.5 3.8
##  [649] 4.5 4.6 4.3 4.3 3.6 5.3 3.7 4.6 4.7 5.8 3.8 3.5 2.9 5.8 5.0 4.7 4.7 2.9
##  [667] 4.1 4.2 3.2 3.4 2.4 5.0 3.0 5.1 4.9 4.1 3.9 4.1 3.9 4.3 4.4 6.3 4.5 5.6
##  [685] 2.8 4.9 4.7 4.3 5.4 3.2 5.1 4.1 3.8 2.2 5.9 5.5 3.4 6.0 4.9 4.7 5.8 4.7
##  [703] 3.3 3.4 2.8 5.7 6.1 4.2 4.4 5.1 5.5 4.5 5.5 3.2 5.0 6.4 5.7 4.8 4.0 5.1
##  [721] 4.9 5.0 3.7 4.1 5.0 4.3 5.2 6.1 4.6 3.2 5.6 6.5 4.1 3.4 4.3 5.0 3.9 2.9
##  [739] 5.7 5.3 3.4 4.3 5.4 4.2 3.9 6.0 6.3 5.8 3.6 3.5 4.7 5.0 3.0 6.5 4.9 4.8
##  [757] 3.4 4.9 4.2 5.0 5.4 2.5 3.7 4.8 6.3 6.4 4.9 4.4 4.2 3.3 5.5 4.2 4.2 5.7
##  [775] 4.9 5.4 5.1 3.1 2.9 4.7 4.1 6.2 3.8 4.3 2.8 5.3 3.3 3.0 5.9 4.5 5.2 4.5
##  [793] 6.0 3.5 4.8 4.2 2.7 5.9 4.1 3.5 5.1 4.5 3.3 4.9 5.0 3.4 3.7 3.4 4.2 3.1
##  [811] 3.9 3.6 3.1 4.7 6.1 5.2 4.8 3.1 4.1 4.2 5.4 5.4 3.3 4.2 5.2 3.7 4.3 5.3
##  [829] 3.5 4.1 4.9 4.7 4.6 3.0 4.1 5.5 5.1 3.3 5.7 4.9 3.8 2.7 4.5 4.3 4.2 5.0
##  [847] 5.0 3.6 4.2 3.0 4.0 5.7 5.2 3.7 5.4 3.5 5.2 2.1 3.7 3.6 4.1 3.7 5.5 3.2
##  [865] 3.7 5.1 4.7 5.1 3.4 4.7 3.4 4.8 5.1 5.3 4.7 4.1 3.8 6.8 5.0 4.8 5.4 3.7
##  [883] 3.8 4.1 3.4 5.9 4.3 5.0 3.4 5.1 4.0 4.5 3.3 4.7 5.3 4.5 3.9 3.2 5.7 3.9
##  [901] 5.2 4.7 6.2 3.2 6.7 5.8 5.6 3.2 4.8 5.0 5.5 3.1 4.0 5.8 4.6 2.9 3.4 3.1
##  [919] 5.5 3.6 4.1 4.3 4.1 2.9 5.7 4.6 2.8 4.7 5.1 5.2 3.8 4.4 5.4 5.2 3.2 3.2
##  [937] 4.1 6.9 4.7 3.2 5.3 5.4 4.3 3.6 4.4 4.7 3.6 4.8 5.9 5.4 4.0 6.2 3.4 4.2
##  [955] 4.1 5.0 4.2 4.4 4.5 5.8 4.2 4.4 3.6 6.1 5.2 4.3 4.8 5.0 5.3 6.9 4.7 4.3
##  [973] 4.1 4.9 5.1 4.1 2.7 5.3 5.1 5.5 3.3 4.0 5.0 3.6 4.2 3.3 3.9 2.5 3.6 4.8
##  [991] 3.3 4.3 6.1 4.1 4.5 4.0 3.1 5.2 4.2 2.0
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
##   2.7   6.3
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
##    [1] 5.0 3.9 4.2 5.3 3.4 3.7 6.9 4.9 3.6 3.8 4.0 3.6 3.3 3.3 3.6 3.8 6.0 4.6
##   [19] 4.2 4.4 4.0 3.9 2.0 4.2 5.3 4.7 3.2 4.1 5.0 3.3 5.2 5.7 4.1 3.7 4.8 4.3
##   [37] 2.5 3.6 5.1 4.7 4.4 3.6 3.6 5.5 5.0 5.2 4.3 4.0 3.9 3.8 4.6 5.0 4.7 4.1
##   [55] 4.4 4.3 4.8 4.4 3.1 4.6 6.3 4.9 4.2 3.6 5.2 3.8 3.9 5.0 4.1 3.7 4.7 5.2
##   [73] 4.2 3.6 2.9 4.9 4.5 7.0 5.2 2.6 5.2 3.8 4.9 5.2 4.1 5.5 4.2 3.9 3.8 2.8
##   [91] 3.4 3.0 3.3 3.2 3.4 4.6 6.4 5.2 2.7 2.9 3.8 5.6 4.0 5.1 3.5 3.4 4.7 3.2
##  [109] 5.5 5.1 4.2 4.3 4.7 5.6 3.6 4.8 5.2 5.1 4.7 5.2 4.5 3.0 5.2 4.4 4.2 5.2
##  [127] 3.2 2.5 4.3 3.5 3.5 6.4 3.5 3.6 3.5 5.6 5.4 4.7 4.6 4.9 3.3 5.8 4.8 5.9
##  [145] 3.5 2.0 2.5 3.1 3.3 4.2 6.5 4.0 5.0 5.8 4.3 4.5 4.9 3.3 3.3 2.3 4.0 4.4
##  [163] 4.3 4.5 4.7 6.8 4.9 5.4 4.3 4.7 2.6 3.7 5.6 5.4 5.3 3.1 6.0 4.0 5.4 4.6
##  [181] 4.9 4.4 4.7 5.1 3.8 5.0 4.6 5.4 5.0 3.1 3.2 2.8 5.7 5.5 4.2 2.1 3.9 4.7
##  [199] 4.9 4.8 4.0 6.9 3.2 6.0 4.6 4.0 6.5 3.3 5.5 5.1 5.3 3.5 4.6 5.6 4.9 3.6
##  [217] 5.3 4.5 3.4 4.5 2.7 4.6 4.9 3.3 4.7 4.2 4.1 5.8 4.0 3.3 4.1 3.9 5.6 4.0
##  [235] 4.4 4.9 3.8 3.9 4.8 3.9 4.9 5.5 4.9 4.8 4.7 4.1 5.3 3.9 4.2 4.8 5.5 4.0
##  [253] 3.6 5.1 6.6 4.5 4.8 4.8 5.6 3.6 4.4 5.5 3.7 5.5 5.6 4.6 3.9 4.2 4.6 3.7
##  [271] 4.3 3.6 4.0 3.6 2.6 5.8 4.8 3.3 4.9 3.0 6.1 3.8 4.8 5.2 5.3 5.7 3.0 4.9
##  [289] 5.3 3.5 3.4 5.4 3.5 2.9 4.7 5.0 4.8 4.5 3.6 6.0 5.0 5.3 5.7 5.2 4.7 3.6
##  [307] 5.3 1.9 4.9 3.7 4.3 4.0 4.9 5.7 4.7 4.8 3.6 4.3 4.6 4.0 3.8 4.3 4.8 5.9
##  [325] 3.4 3.6 5.3 5.2 4.8 4.7 4.3 4.3 5.1 4.7 4.2 3.3 4.8 2.9 5.1 3.4 5.0 5.7
##  [343] 4.7 4.4 5.4 5.4 5.0 4.3 4.2 7.0 5.1 4.4 5.7 5.4 6.0 4.0 5.0 5.1 5.7 3.0
##  [361] 6.3 4.5 4.1 1.9 5.5 4.4 5.2 5.2 3.6 4.1 4.3 2.5 4.5 3.8 4.2 3.9 5.3 4.2
##  [379] 4.4 3.3 5.4 4.3 4.6 5.1 3.6 5.4 4.4 3.6 3.3 4.7 4.0 5.9 5.6 4.4 4.2 3.1
##  [397] 4.9 3.8 5.3 4.2 4.0 3.8 5.0 5.0 4.7 3.2 5.6 4.5 5.0 4.1 4.2 5.2 5.4 6.1
##  [415] 4.2 5.0 5.4 4.7 5.1 3.8 4.3 4.8 3.9 5.0 4.7 3.5 4.5 4.2 3.6 4.4 5.4 4.8
##  [433] 4.3 4.1 4.9 5.1 4.7 4.6 4.7 5.1 4.4 5.9 5.4 4.9 2.7 3.8 4.5 3.7 3.3 3.9
##  [451] 4.2 5.1 3.2 4.9 6.7 4.1 4.5 2.6 3.7 3.7 5.3 4.9 5.5 5.1 3.4 5.7 5.2 4.5
##  [469] 5.5 3.6 4.9 5.8 5.9 4.4 4.5 5.2 3.6 4.0 5.1 3.9 3.9 7.3 3.7 5.1 3.7 5.4
##  [487] 4.7 4.2 4.7 3.8 4.5 4.2 4.3 5.2 4.5 4.0 4.8 4.0 3.1 3.8 4.0 2.9 3.2 4.0
##  [505] 5.1 4.4 4.5 4.8 4.6 4.6 5.2 3.7 5.7 3.6 5.0 4.1 4.5 5.6 5.8 5.1 4.1 6.0
##  [523] 4.2 3.8 4.1 4.9 5.3 4.6 5.0 4.4 5.7 4.9 4.6 5.5 3.9 5.4 4.4 4.4 3.5 7.0
##  [541] 5.2 4.1 5.0 2.5 3.2 4.0 5.3 3.5 4.0 4.8 5.9 5.3 2.5 4.6 4.8 3.5 3.4 3.8
##  [559] 5.1 6.3 5.3 3.9 5.6 3.5 3.3 5.0 4.8 3.2 3.3 4.5 4.8 5.8 3.5 5.1 3.4 3.5
##  [577] 3.4 5.9 4.9 3.0 3.2 3.8 4.5 5.0 5.3 5.8 5.7 3.7 5.7 4.9 4.5 5.9 3.9 4.0
##  [595] 5.5 5.1 4.6 3.9 7.4 5.4 3.5 3.6 4.4 3.7 5.4 3.8 2.9 5.3 4.8 5.3 3.1 5.2
##  [613] 5.3 5.5 5.2 4.7 6.3 5.6 4.4 3.4 4.4 4.4 5.7 4.5 3.6 4.8 3.4 3.6 4.7 4.2
##  [631] 5.4 4.4 4.6 3.1 3.3 5.9 5.5 3.6 5.2 4.9 5.6 4.3 3.1 5.6 3.4 5.8 4.7 4.1
##  [649] 4.6 4.2 4.1 5.5 5.5 3.9 4.6 3.7 3.6 5.0 4.1 3.9 4.8 3.2 3.9 3.8 4.8 2.8
##  [667] 6.3 5.4 6.2 4.4 4.9 3.8 3.6 2.3 4.3 3.6 4.5 4.1 4.5 4.6 5.0 5.0 3.3 5.2
##  [685] 4.2 4.9 3.0 5.1 4.7 4.2 3.7 2.7 4.4 5.1 3.1 4.9 3.3 4.5 5.6 4.0 3.7 3.6
##  [703] 4.5 5.5 3.7 6.0 4.6 4.5 6.4 6.4 3.8 4.0 4.3 3.4 5.1 3.4 2.6 3.6 5.1 4.3
##  [721] 4.8 3.3 4.0 4.8 4.4 3.5 4.6 4.7 5.2 4.5 3.6 6.1 5.1 5.2 2.9 5.3 3.3 4.4
##  [739] 3.6 4.9 5.2 2.9 5.6 5.0 3.5 3.0 4.4 5.2 4.7 5.2 4.2 5.0 5.2 5.5 5.1 3.8
##  [757] 4.4 3.3 5.3 5.4 4.8 4.1 4.2 4.1 5.3 4.4 5.5 3.8 4.7 3.8 3.8 4.1 4.7 3.0
##  [775] 5.2 4.5 4.9 4.5 4.3 5.0 3.9 3.3 4.7 7.0 4.4 4.3 3.0 2.8 5.1 5.2 3.8 4.0
##  [793] 6.9 4.6 5.1 4.1 5.2 5.3 4.7 4.8 4.2 4.6 5.1 3.9 4.1 4.6 4.2 3.8 4.9 4.4
##  [811] 5.5 3.3 2.8 4.0 3.8 5.4 5.8 5.5 4.7 4.2 2.9 3.6 5.5 4.8 5.5 3.6 4.0 5.0
##  [829] 4.2 3.7 3.8 3.9 4.4 4.9 5.3 2.8 5.1 4.4 6.3 4.8 6.0 2.9 4.3 3.8 5.0 4.2
##  [847] 5.1 5.2 5.5 5.0 5.2 3.4 3.9 5.6 3.7 4.8 5.4 5.1 3.4 3.3 4.0 4.5 5.6 4.1
##  [865] 3.9 4.0 4.3 3.5 4.4 4.1 3.4 4.2 3.6 4.5 3.7 3.4 4.5 5.9 3.6 5.5 4.5 6.1
##  [883] 4.8 4.5 6.8 3.3 4.8 4.9 4.0 5.8 4.3 3.6 3.7 6.3 3.8 4.0 5.1 5.1 4.9 3.1
##  [901] 3.5 5.4 5.6 5.8 6.3 4.1 4.5 3.9 4.6 4.6 4.5 4.1 6.1 3.0 2.2 7.0 3.1 4.6
##  [919] 3.0 5.2 4.4 5.9 3.7 4.4 5.1 4.5 6.1 4.5 5.0 4.8 3.2 3.7 4.6 4.6 4.2 4.9
##  [937] 3.4 4.1 3.9 5.3 4.6 3.2 3.2 2.7 3.6 4.3 4.7 6.5 4.9 2.5 4.6 5.4 3.6 4.6
##  [955] 5.5 5.0 5.0 3.5 6.3 4.8 5.2 3.9 4.7 3.2 4.6 3.3 4.1 5.1 5.3 4.3 4.4 4.7
##  [973] 6.0 3.4 3.9 4.2 4.4 3.8 3.8 3.1 3.3 3.2 3.1 5.5 4.2 5.1 4.0 4.8 5.9 4.5
##  [991] 5.9 3.4 4.0 5.6 3.9 5.8 3.0 5.4 4.1 3.5
## 
## $func.thetastar
## [1] -0.0351
## 
## $jack.boot.val
##  [1]  0.48664773  0.38543417  0.28807339  0.11636905 -0.02149254 -0.12588556
##  [7] -0.19631902 -0.28296089 -0.45864865 -0.54735294
## 
## $jack.boot.se
## [1] 1.001345
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
##    [1] 5.3 5.3 5.4 4.5 3.9 4.6 5.0 5.6 6.1 4.2 5.3 3.5 4.1 2.7 4.2 3.5 3.5 5.4
##   [19] 4.5 4.6 4.9 4.8 3.2 3.4 3.3 4.8 4.5 4.2 4.5 4.2 3.2 3.7 6.4 3.9 3.3 5.2
##   [37] 6.1 5.1 4.8 4.2 4.4 4.3 3.2 4.7 5.5 4.7 3.3 4.4 4.3 4.5 4.5 3.9 5.3 6.8
##   [55] 5.2 3.9 3.3 4.7 4.9 5.1 5.8 4.2 6.9 4.5 5.8 3.9 5.6 5.0 5.2 4.7 4.2 5.8
##   [73] 3.2 6.5 5.1 4.9 2.7 4.8 3.6 5.0 5.1 3.1 2.6 5.5 4.6 3.7 7.0 5.0 4.0 4.2
##   [91] 5.6 3.3 6.2 4.7 4.3 3.8 3.8 4.8 3.6 5.3 4.3 5.5 5.5 2.8 4.2 4.7 4.1 4.1
##  [109] 3.9 4.6 3.9 4.4 4.2 6.0 4.2 2.3 5.3 4.9 3.4 5.6 3.6 4.9 2.6 3.6 4.2 4.9
##  [127] 3.6 4.1 4.8 3.6 5.9 4.8 4.7 4.3 5.6 2.8 3.7 4.1 4.0 4.0 4.8 4.3 5.9 4.1
##  [145] 4.4 4.3 3.1 4.1 5.2 2.8 3.0 3.5 4.1 3.7 3.0 4.0 4.1 3.9 4.5 4.2 4.2 4.9
##  [163] 4.7 4.4 4.1 5.3 4.9 4.6 4.1 3.5 3.2 3.4 5.7 3.5 2.6 4.5 4.3 3.9 3.5 6.0
##  [181] 5.1 5.4 3.1 5.3 4.5 5.7 4.6 4.7 5.6 4.5 2.7 4.7 3.3 4.4 4.5 6.9 4.6 3.6
##  [199] 4.5 5.3 3.9 3.3 4.0 4.9 4.4 3.6 4.5 4.9 4.9 5.5 3.5 5.4 4.1 3.0 4.0 5.1
##  [217] 5.6 3.1 5.1 5.5 3.4 5.6 5.7 4.1 5.2 6.2 4.0 3.3 5.3 5.4 5.3 2.7 5.3 2.9
##  [235] 4.4 5.4 4.5 4.6 4.3 4.2 4.4 3.4 3.2 5.3 4.6 4.2 4.0 5.4 3.1 4.8 4.8 4.4
##  [253] 3.6 4.4 4.5 2.9 4.8 4.1 5.0 4.7 4.4 4.5 3.6 4.3 4.4 4.8 3.3 6.5 4.9 3.6
##  [271] 3.5 4.9 3.1 5.8 5.2 6.2 4.2 4.2 2.8 6.4 3.9 5.7 5.0 4.2 4.1 4.6 5.5 2.4
##  [289] 4.5 3.8 4.4 4.0 5.2 4.8 5.7 5.2 4.3 3.5 3.7 5.5 5.6 5.3 6.2 2.9 4.3 4.3
##  [307] 4.3 6.2 4.1 5.1 6.4 4.2 4.9 6.2 3.1 6.0 3.7 4.9 4.9 4.9 4.8 3.5 4.1 5.7
##  [325] 3.7 3.3 3.2 4.8 4.6 3.4 4.0 4.7 4.6 3.3 4.6 5.1 4.1 5.3 3.8 3.8 3.0 5.6
##  [343] 3.7 3.8 4.6 4.7 5.4 5.1 3.1 4.7 4.8 5.6 4.8 3.5 5.4 4.3 6.5 4.0 5.8 4.0
##  [361] 4.4 5.5 3.0 3.8 4.4 6.4 4.5 3.6 4.5 4.9 7.1 4.5 5.1 3.8 5.5 5.0 5.4 4.7
##  [379] 6.3 3.8 3.9 3.8 4.9 5.3 4.9 4.6 6.2 3.7 5.1 4.0 4.6 4.3 4.5 1.2 5.1 4.9
##  [397] 5.6 4.1 3.8 4.8 5.2 4.9 4.3 4.9 5.0 4.7 5.4 4.8 3.1 3.9 4.7 4.9 4.5 4.8
##  [415] 4.9 5.5 6.3 5.4 2.4 3.6 5.2 3.9 4.3 4.0 4.1 4.2 2.2 5.2 3.9 5.2 5.0 6.4
##  [433] 3.9 5.7 4.3 4.1 5.4 4.5 3.8 3.4 3.7 4.6 6.4 5.3 4.5 5.3 4.4 4.4 3.5 4.6
##  [451] 4.0 5.2 3.4 4.0 5.0 3.5 3.4 3.1 4.4 4.5 4.4 3.2 4.5 4.5 5.2 3.8 4.5 5.0
##  [469] 4.5 4.4 3.7 5.1 3.1 6.0 3.0 4.7 3.4 3.7 6.6 4.3 4.1 4.9 4.9 6.0 4.3 4.7
##  [487] 5.3 5.5 5.6 5.4 3.9 5.1 4.9 3.8 4.0 4.6 4.4 4.1 4.7 2.3 4.7 5.2 5.4 2.5
##  [505] 3.7 4.0 4.7 2.6 4.1 7.4 4.3 4.6 3.1 5.1 3.9 2.6 5.6 5.1 3.2 4.0 3.2 5.7
##  [523] 3.6 4.5 5.2 3.7 4.9 5.4 5.2 4.2 4.7 3.1 4.6 4.5 4.0 4.4 3.7 4.2 3.8 4.1
##  [541] 4.2 4.1 3.8 4.8 5.9 4.0 5.5 4.4 3.6 3.5 4.0 3.8 4.4 2.9 5.0 3.7 3.7 4.1
##  [559] 3.7 3.8 4.3 4.2 5.2 3.7 4.9 6.0 5.4 3.2 4.9 4.6 3.5 4.5 4.4 3.6 5.5 4.5
##  [577] 6.1 2.7 3.9 5.1 4.3 4.3 4.8 3.4 4.1 3.7 4.3 3.0 5.7 3.5 5.2 3.5 5.5 5.0
##  [595] 4.5 4.9 3.3 6.8 4.0 5.2 3.7 4.8 5.1 6.7 4.7 4.2 4.2 3.3 3.4 5.6 4.4 4.1
##  [613] 3.7 5.1 5.3 5.1 4.3 5.2 5.8 4.5 2.5 4.1 4.2 6.2 5.0 4.9 4.0 2.8 4.7 5.3
##  [631] 6.1 6.0 3.8 3.2 3.0 4.6 5.9 4.1 5.5 5.0 6.3 4.3 4.7 2.1 5.7 4.7 2.8 3.1
##  [649] 3.6 6.3 4.4 5.0 5.0 4.9 3.6 2.8 4.3 4.9 2.2 3.8 4.6 5.4 4.8 4.6 4.3 4.4
##  [667] 6.0 3.2 3.7 6.0 3.3 5.6 5.2 4.4 5.1 6.6 4.7 2.3 4.9 5.0 3.8 2.8 3.7 3.7
##  [685] 3.4 4.6 6.5 2.0 5.5 4.9 5.4 4.0 3.3 5.2 5.7 4.1 5.4 5.5 4.8 3.2 7.2 4.7
##  [703] 4.4 4.0 3.7 5.5 3.4 4.4 4.8 4.2 4.1 5.2 5.0 5.1 4.9 3.6 4.7 5.1 4.5 4.3
##  [721] 3.6 3.1 6.0 2.5 4.2 2.4 4.5 4.9 3.5 5.5 5.2 4.3 3.8 3.9 3.1 5.7 4.8 6.4
##  [739] 5.0 3.4 3.9 5.8 3.1 4.9 3.9 4.0 5.4 5.5 4.6 3.5 5.6 4.0 5.2 4.0 6.2 4.5
##  [757] 4.3 2.9 3.2 5.3 3.3 4.4 4.3 7.1 3.1 4.2 3.9 4.3 5.9 3.8 4.3 4.8 5.4 5.6
##  [775] 4.3 4.5 4.3 3.9 3.9 5.7 5.0 4.4 5.6 4.8 4.1 4.4 4.7 5.2 3.1 4.8 2.9 5.8
##  [793] 4.8 4.1 4.5 4.2 6.0 3.6 4.8 4.2 3.6 5.5 4.8 4.9 4.1 3.8 4.8 3.6 2.7 4.0
##  [811] 4.2 5.5 5.5 4.4 4.5 5.6 3.2 6.0 3.6 2.1 4.0 3.9 4.4 4.6 3.9 3.8 3.9 5.3
##  [829] 4.2 4.7 4.0 6.0 5.5 3.6 4.9 3.9 5.7 3.8 4.6 6.4 4.0 5.7 5.1 3.2 3.9 6.2
##  [847] 5.4 6.0 4.7 3.8 5.1 2.9 4.2 5.2 4.3 4.8 3.3 3.6 4.4 4.4 5.4 4.6 3.4 4.4
##  [865] 4.3 3.5 3.8 3.2 2.7 4.9 3.9 3.9 4.6 5.3 5.0 4.2 5.1 4.4 4.7 3.1 6.1 5.6
##  [883] 5.5 4.8 5.3 6.1 4.7 4.8 5.3 4.8 3.6 5.4 3.5 4.9 4.5 5.2 4.3 4.0 5.4 4.1
##  [901] 3.5 4.3 4.9 5.2 3.9 4.7 4.9 4.5 4.5 4.6 2.8 5.4 3.2 3.0 2.5 3.5 3.4 5.1
##  [919] 5.5 3.7 3.6 5.1 3.5 4.7 4.3 4.2 6.0 3.8 3.9 3.4 5.2 5.1 4.4 4.1 6.0 3.8
##  [937] 5.5 5.3 5.1 5.3 5.1 2.8 5.1 4.0 4.6 3.3 4.9 4.9 4.6 4.1 3.6 1.7 4.2 3.8
##  [955] 4.5 3.8 4.5 5.2 4.5 4.7 3.7 5.3 4.5 5.3 5.8 5.4 5.5 4.7 4.7 4.9 4.4 5.4
##  [973] 5.3 5.3 6.7 4.2 5.1 6.1 6.0 4.4 4.7 4.6 4.4 4.5 5.5 6.1 5.7 3.8 3.6 4.5
##  [991] 4.0 2.6 4.9 4.2 4.6 4.3 3.7 4.3 3.8 4.7
## 
## $func.thetastar
## 72% 
##   5 
## 
## $jack.boot.val
##  [1] 5.500 5.400 5.300 5.208 5.200 4.900 4.900 4.644 4.700 4.400
## 
## $jack.boot.se
## [1] 1.029524
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
## [1] 0.7464885
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
##   1.835559   3.907881 
##  (0.758100) (1.853895)
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
## [1] -0.4561257  2.0850922  1.5857577  1.4703885  1.7973359  1.3926469
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
##    [1] -0.309222158  0.920113992  2.177881616  0.789130425  0.515020770
##    [6]  0.680808282  0.115107203  0.553323187  0.260083641  1.137215506
##   [11]  0.514722573  0.966955944  0.898916665  0.080560802  1.038618135
##   [16]  1.398403604  1.309314388  0.908162998  1.056601599 -0.120484884
##   [21]  0.776350091  1.650052880  0.027425848  0.844682468  0.265709147
##   [26]  0.288987620  1.216438528  1.221716048  0.814951307  1.995707571
##   [31]  0.282095266 -0.136479911  0.353432061  0.520595996  0.719201999
##   [36]  0.953784524  1.158891400  0.303894000  0.745749917  0.888041065
##   [41]  1.030801542  0.739152950  1.253735524  0.101442027  1.219578320
##   [46]  0.522580157  1.310422166  0.053197544  0.412719053  0.778864592
##   [51]  1.085659941  0.484648921 -0.039455466  0.913389307  0.040076528
##   [56]  0.736174961  1.012304295  0.419973283  0.676043171  1.630435046
##   [61]  1.666375690  0.903379591  0.431267462  0.468488576  1.104190670
##   [66]  0.633788295  0.467465151  0.202542406  0.476648150  0.380682257
##   [71]  0.226295306  0.518363578  0.560611284  0.432265701  0.398077903
##   [76]  1.137837979  1.269149792  1.865161857  0.228188923 -0.134947225
##   [81]  1.169284715  0.742615341  0.590453183  0.170812117  0.139312738
##   [86]  1.040917068  0.363344852  0.086059258  0.345545080  0.743145380
##   [91]  0.278484787  1.531916485  0.744672904  0.785793473  1.124739892
##   [96]  1.368963807  1.021288956 -0.017860527  0.927401757 -0.813493685
##  [101]  0.295461467  1.236139008  0.355518999  0.631649696  0.126325605
##  [106]  0.287883386  1.124351311  1.119592326  0.212268545  0.331918141
##  [111]  1.375268875  0.212302119  0.767566835  0.300692761 -2.073737117
##  [116]  0.326027740  1.523806131  0.714798947  1.047385659  0.038914638
##  [121]  1.688434528  0.885663447  1.006118534  0.747972968  1.047893688
##  [126]  1.781490441  0.695087441  0.374659371  0.851830758  0.741365467
##  [131]  0.368920839  0.165423740  0.682309275  1.643904528  0.487749508
##  [136]  0.712402128  0.954331288  0.760487270  0.199907889  0.493595117
##  [141]  0.091201938 -0.304563935  1.393441472  0.744909848  0.722334934
##  [146]  0.231790874  1.173464086  0.581704157  0.924067890  0.660975480
##  [151]  0.637420125  2.258064600  1.064621377  0.400985566  0.274725562
##  [156]  0.176759232  1.125262108  0.317612408  0.758866040  0.453098586
##  [161]  1.060497934  0.781001978  0.476998856 -0.426334611 -0.054907016
##  [166]  0.270803929 -0.368180874  0.409283337  0.616537653  1.754812404
##  [171]  0.913393262  0.260806259  0.396154086  0.992225181  0.654758230
##  [176]  1.102864167  0.164573370  0.339830272  0.594959155  1.287267588
##  [181]  0.393423649  0.015580315  0.635681194  1.411061891  0.361358960
##  [186]  0.442697801  0.727383568  0.289252725  1.265547309  0.948492010
##  [191] -0.144483241  0.489441039 -0.030455551  1.247870523 -0.287617750
##  [196]  0.033154265  0.476717586  1.231695766  0.759718681  0.794319454
##  [201]  0.646767016  1.047765386  0.024115808  1.765195271  0.445743362
##  [206]  0.402921689  0.658684791  1.944086472  1.452130519  0.259029539
##  [211]  0.509609721 -0.051280804  0.428049819  0.437830728  0.960204535
##  [216]  0.540530075 -0.331108369  0.032017649  1.757375746  1.669783178
##  [221]  0.361692877  1.247308207  1.013570173  1.117162696  0.639629243
##  [226]  0.307946654  0.806783223  0.480915053  0.634559314  0.633539254
##  [231]  0.614531835  0.116016762  0.328820508  0.751658283  0.802995979
##  [236]  1.119397076  1.320032084  0.888629760  0.410720430  0.186717659
##  [241]  0.728011463  0.557762606  0.491034186 -0.002800456 -0.306371203
##  [246]  0.588811454  0.251977266  0.239996966  1.044284462  1.299168603
##  [251]  0.201606889  0.395066826  1.117004670  0.560348454  1.432285835
##  [256] -0.035802194 -0.293867305 -0.243954354  0.161741055  1.687283360
##  [261]  1.411061891  1.184983867  0.614316856  0.449854628  0.703491687
##  [266]  0.853861183  1.028798127  1.653987012  0.295070672  0.876333743
##  [271] -0.316943099  0.943724671  0.634559314 -0.358778655  2.397417332
##  [276]  0.437835051  0.842357860  0.923411539  1.413988464  1.024072380
##  [281]  0.552952713  0.474666431  0.722318817  0.914155251  2.248683636
##  [286]  0.960357210  0.958815252  0.888779248  1.080015251  0.006797408
##  [291]  1.514990967  0.688559166  0.291686421  0.636117099  0.264779505
##  [296]  0.798557834  0.652294591  2.223792884  0.525329324  0.393061025
##  [301]  0.368576379  1.661378294  0.661226209  1.472033176  1.057462590
##  [306]  1.251049564  1.098649076  0.614777506  0.431578722  0.366577104
##  [311]  1.342949844  1.013389410  0.244466782  1.175542180 -0.271726355
##  [316]  0.563708151  1.364024732  0.605842661  0.329296615  1.021114919
##  [321]  1.881440109  0.129696284  0.307307030  0.766662783  0.467942993
##  [326] -0.344400322  0.510345859 -0.069606315  0.762905993  0.435605179
##  [331]  1.411161856  0.358468310  1.072424157  0.307692697  0.516640740
##  [336]  0.245649227  0.967282674  0.955320051  0.480915053  0.610485609
##  [341]  0.429572644  0.176787573 -0.367787622  0.642181377  0.725803845
##  [346]  0.386656259 -0.365309210  0.274023423  1.110575077 -0.881079747
##  [351]  0.277017273  0.857419122  1.294421297  1.801647161  0.070105764
##  [356]  0.283163010  1.254031709  0.632121473  0.187398027  1.105648517
##  [361]  0.916194523  2.342758632  0.501423219  0.792178451  0.614112164
##  [366]  1.377927808 -0.037817084  0.473884687 -0.429371872  0.237622513
##  [371]  1.232748396  0.223841791 -0.016584733  1.549197035  0.053915597
##  [376]  0.412379878  1.049252985 -0.138090471  0.364725854  0.639660070
##  [381]  0.642933935  0.379559582  0.638003043  0.738174502  2.345620905
##  [386] -0.002703213  0.734555813  0.919352741  0.830656278  0.503764092
##  [391]  0.031555756  0.053251784  0.575435083  0.608568979  0.677953057
##  [396]  0.650463561  0.762676446  0.552306535  0.743242030  0.557512755
##  [401]  1.131579032  1.051685694  1.067415124 -0.556318290  0.525585916
##  [406]  0.757167024  0.819295802  0.594396333  0.233506094  1.640997608
##  [411]  0.692251789  0.566428101  0.914329489  0.410642484  0.290903726
##  [416]  1.081650755  0.561963723  0.839264703 -0.022033491  0.178555724
##  [421]  0.133269228  0.238590772  0.908146805  0.420560039  0.471835885
##  [426]  0.554643076 -0.099842841  1.117126754  0.920141766  0.661987362
##  [431]  0.469835052  0.384329192  0.822250519  0.344088598  0.996486800
##  [436]  0.348938639  0.588811454  0.271617726  0.762404917 -0.404465528
##  [441]  0.913377325  0.962931948  0.781760164 -0.089144911  0.750571676
##  [446] -0.146491237  0.782383878  0.608504107  1.106127167  0.520980981
##  [451]  0.409595488  0.880245486  0.234681475  1.054162491  0.748453275
##  [456]  0.574900275  1.038267662  0.261719569  1.137366445  0.869147727
##  [461]  1.353235309  1.738809231  1.038445684  0.763881361  0.209076000
##  [466]  1.117448911 -0.148122254  2.267253267  0.630865202  0.231137544
##  [471]  0.474777378 -0.160595031  0.361772877 -0.058563179  1.739397902
##  [476]  0.771471661  1.175529424  1.148916304  0.777436240 -0.101634875
##  [481]  0.614777506  0.678868682  0.314721445  1.343023018  0.334722017
##  [486]  0.650688558  0.490690053  0.880837593  0.599550240  0.016831325
##  [491]  1.186863006  0.366372096 -0.086194342  1.646489273  0.104155436
##  [496]  1.234699401  0.047077535  1.471747215  0.568981478  0.293901305
##  [501]  0.137962046  1.457399900  0.348883871  0.856494164  0.385197579
##  [506]  0.279947608  0.983681314  0.793837846 -0.091982534  0.654764621
##  [511]  0.111960403  0.374131124  0.757282625  0.614102308  0.707156206
##  [516]  0.306579334 -0.014115497  0.810152613  1.293629921  0.745749917
##  [521] -0.100281054  1.221298069  0.317660615  0.790127121  1.202153928
##  [526]  0.437037768  1.239265632  0.736809545 -0.011602927  0.058500984
##  [531]  0.868857942  0.377024105  0.366084112  0.680479617  0.714485445
##  [536]  1.075928357  0.446148171  0.641687377  1.029630991  0.918932570
##  [541]  0.867534142  0.730108243  0.276082941  0.104951190  1.162566261
##  [546]  0.858928423  0.979631007  0.249710003  0.816072973  0.406207112
##  [551] -0.011602927 -0.086000316  1.077497893  0.422017146  0.870003476
##  [556]  0.328479629  0.168083184  0.888328085 -0.032679695  0.636238108
##  [561]  0.757521358  1.020623314  1.654386209  0.223024991  0.595844946
##  [566]  0.298542856  1.433468757  0.307720167 -0.252726781  1.258283636
##  [571]  0.641622968  0.740540466  0.323650349  1.085406044  0.614965178
##  [576]  0.224749540  0.664944392  2.310076017  1.097856105  0.815697424
##  [581]  0.407148072  0.650252897  1.301448766  0.113861526  0.474809480
##  [586]  1.352376746  1.281204358  0.728011463  1.086989674  0.848925856
##  [591]  0.454186097  0.041877585  0.774668476  0.265181047 -0.015149186
##  [596]  1.454313821  0.435389985  0.824053497  0.470067023  0.075768347
##  [601]  1.029173418  0.346394103 -0.004604438  0.316519453  0.514972168
##  [606]  1.701401660  0.750571676  0.455759474  0.515987755  0.559258832
##  [611]  0.191204152  1.132104019  1.024796209  1.665395447  0.177106526
##  [616]  0.199985139  0.186353429  1.409852495  0.911497731  0.634035517
##  [621]  0.689063641  0.839062717  0.158959359  1.461354901  1.346401085
##  [626]  0.306507226  0.731295486  0.134837185 -0.009525585  0.284382087
##  [631] -0.179793373  0.871511469  0.771896172 -0.397390386  0.353175821
##  [636]  2.553683505 -0.137687209  0.317982746 -0.570631729  0.241437619
##  [641]  0.931378030  1.073987430  0.063964726  0.031555756  0.295228583
##  [646]  2.277151767  1.207750249  0.889386836  0.084077165  0.132905338
##  [651]  0.728852312  1.193357008  1.711221420  0.888846859  0.888224190
##  [656]  0.908217770  1.122319102  0.908765669  0.570374330  0.755217734
##  [661]  0.029639853  0.053831177  1.112426581  1.244988640 -0.380290538
##  [666]  1.273720881  1.654746836 -0.464414685  0.301836915  1.295822698
##  [671]  2.271868069  0.158045981  1.132346403 -0.356349497 -0.135283438
##  [676]  0.435500536  1.133069124  0.933963303  0.769394330  0.148127670
##  [681]  1.154231164  1.098591761  0.324408260  0.670317209  1.091802272
##  [686]  0.491280089 -0.044024958 -0.005057424  0.960667903  1.090784944
##  [691]  1.194388132  1.704168447  0.135472962  1.257634444 -0.035705964
##  [696] -0.050756610  0.627746974  0.545033472  1.506380187  1.029680829
##  [701]  0.285091365 -0.673684175  0.517935379  0.293922793  1.103001180
##  [706]  0.089279845  1.073370994  1.471747215  0.774077202  0.644960277
##  [711]  0.367433674  0.808937012  0.732513342  1.092142989  0.225400538
##  [716]  0.495565056  1.112546260  0.366049167  0.788319964  1.273584545
##  [721]  0.724313930 -0.469393318  0.622117586  1.693025710  0.239996966
##  [726]  0.745147175  0.639807817  0.363336791  0.451350602  1.424785925
##  [731]  1.469580500  0.765698846  0.317612408  0.463337710  0.252373593
##  [736]  0.172096531  0.863723374  0.756222287  1.308083295  0.730009269
##  [741]  0.457549961 -0.156693427 -0.186368170  1.223104175  0.725844939
##  [746]  0.059841047  1.099769353  0.912107835  1.694782316 -0.057338240
##  [751] -0.035950195  0.264627951  0.690912228  0.566398301  0.637988216
##  [756]  0.830454654  0.422322566  0.642735261  0.503563384  1.274934581
##  [761]  0.238852825 -0.056187172  0.494962922  1.160853385  0.334674991
##  [766]  0.724997233 -0.062171980  0.386343264  0.435382525  0.619782680
##  [771]  0.111501436  0.558802737  0.304691139  1.202153928  0.024115808
##  [776]  0.147864066  2.323586746  0.428177553 -0.340674711  0.069726626
##  [781]  1.685955100  1.709225870  1.154017557  1.179242602  1.048610783
##  [786]  0.457704044  0.738689688  1.657846453  0.459340422  0.881241031
##  [791]  1.470666313 -0.330285278  0.762431798 -0.209159320  0.373424073
##  [796]  0.912497844  0.706451603  0.782231214  0.177699162  1.040026119
##  [801]  1.237624666  0.431906265  0.321787216  1.160607952  0.920141766
##  [806]  1.102380931  1.171376774  1.489606162  0.447625328  1.748401116
##  [811]  1.139194103  0.468482062  0.215229632  0.797813706  0.885592159
##  [816]  0.630479602  0.720620431  0.635787494  0.383404227  0.501831423
##  [821]  0.506379441  1.713303367  0.349661157 -0.387260399  0.122680428
##  [826]  0.530747858  0.051841407  1.072611008  0.895734876  0.686001277
##  [831]  0.527403269  1.204940817 -0.377223138  0.627385177  0.496849856
##  [836]  0.123761997  0.304517057  0.319468733  1.769800085  0.769752464
##  [841]  0.288278250 -0.047160128  1.759831496  0.514255313  0.223841791
##  [846]  0.150996341  1.599660849  1.459090629  1.086090814  1.142583975
##  [851]  0.473049605  0.982071329  1.104979449  0.424238260  0.319155959
##  [856]  0.865882641  0.363871564  1.106118721  1.228084844  0.871016771
##  [861]  1.062987715  0.675552162 -0.036667766  0.306729230  0.715180081
##  [866]  0.396895215  0.920586872 -0.172721633  0.674031021  0.683894927
##  [871]  0.461534342  1.163524078  0.489220888 -0.033750874  0.059672049
##  [876]  1.023831681  1.210015747 -0.017797957  0.913804842  0.091380360
##  [881]  1.525810820  0.002122518  0.745053780  0.398804609 -0.053206480
##  [886]  1.272547477  0.615002688  1.111602288 -0.434564958  0.448974955
##  [891] -0.102966199  0.226065361  0.715456309  0.279326862  0.156214605
##  [896]  0.291081825  0.714655955  1.449907205  0.597988067  1.024104944
##  [901]  0.628844508  0.588958365  1.185017156  0.151299638  0.541865654
##  [906]  0.730478978  0.645402710  0.939067173  0.920167357  0.594396333
##  [911]  1.029138239  0.757324189  0.441794314  0.360532420  0.459081458
##  [916]  0.656780901  0.616849488 -0.006477907  0.049613804  0.699212490
##  [921]  0.731686132  0.513575542  0.370692403 -0.050110992  0.323387982
##  [926]  0.157846834  1.177386072  1.070562782  0.358724410  0.351826715
##  [931]  0.746443756  1.326687014  0.456505723  0.582258943  0.805583633
##  [936]  2.385937926  1.245202416  0.289711695 -0.042235151  0.304469317
##  [941]  0.530414171  0.425398396  0.777315612  0.292932270  0.162280706
##  [946]  1.519071611  0.749334272  0.514512352  1.581965943  0.829421992
##  [951]  0.255446638  0.028708861  0.733683664  0.060983538  0.312366942
##  [956]  0.283963815  0.590890567  0.300531867 -0.337831395  0.354946438
##  [961]  0.328438361  1.056401488  0.661449049  0.058420839  0.831179118
##  [966]  0.555903858  0.725412564  0.090675421  0.567835444  0.786596455
##  [971] -0.520403444  1.032429941  0.783230024  0.597633328  0.417074555
##  [976]  1.078934072  1.781399958 -0.365775261  0.571669481  0.919010088
##  [981]  0.278743408  0.175666906  0.312774871  1.173740521  2.376717686
##  [986] -0.112905422  1.274297597  0.822788749  0.581214136  1.224633831
##  [991]  1.943895382  0.642477458 -0.017969660 -0.080964245  0.632121473
##  [996] -0.044032239  1.300658713  0.663008357  0.758761347  1.475193846
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
```

```r
fit2
```

```
##       mean          sd    
##   0.46970828   0.34188584 
##  (0.10811380) (0.07644467)
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
## [1]  0.263810602  0.469689270  0.131156183  0.086125364 -0.223486204
## [6] -0.000260287
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
## [1] 0.0169
```

```r
sd(results2$thetastar)  #the standard deviation of the theta stars is the SE of the statistic (in this case, the mean)
```

```
## [1] 0.9164957
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
##     original       bias    std. error
## t1*      4.5 0.0004004004   0.9037662
```

One of the main advantages to the 'boot' package over the 'bootstrap' package is the nicer formatting of the output.

Going back to our original code, lets see how we could reproduce all of these numbers:


```r
table(sample(x,size=length(x),replace=T))
```

```
## 
## 0 2 4 5 6 8 9 
## 1 2 1 1 2 2 1
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
## [1] -0.0191
```

```r
se.boot
```

```
## [1] 0.9277855
```

Why do our numbers not agree exactly with those of the boot package? This is because our estimates of bias and standard error are just estimates, and they carry with them their own uncertainties. That is one of the reasons we might bother doing jackknife-after-bootstrap.

The 'boot' package has a LOT of functionality. If we have time, we will come back to some of these more complex functions later in the semester as we cover topics like regression and glm.

