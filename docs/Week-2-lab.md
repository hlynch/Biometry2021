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
## 0 4 7 9 
## 1 1 3 5
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
## [1] -0.0079
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
## [1] 2.765938
```

```r
UL.boot
```

```
## [1] 6.218262
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
##    [1] 6.6 4.4 3.8 2.4 4.1 5.2 4.6 4.8 5.4 2.8 5.1 5.1 3.9 4.5 5.0 3.6 4.5 3.6
##   [19] 4.5 4.0 4.1 5.2 3.6 5.3 6.5 4.9 6.2 5.2 2.9 5.0 4.3 7.0 4.7 5.5 4.1 5.1
##   [37] 4.8 6.3 3.8 5.7 6.0 5.0 3.6 4.2 4.5 5.3 4.5 5.6 4.9 3.8 4.3 4.5 6.6 5.7
##   [55] 4.6 5.0 4.0 4.5 5.1 7.4 4.7 5.0 4.1 2.6 4.2 4.0 4.1 4.1 5.1 3.4 3.4 3.0
##   [73] 3.5 4.4 3.4 3.4 5.8 2.5 4.4 3.6 4.1 5.6 4.4 5.2 5.0 4.1 5.1 3.7 4.5 3.6
##   [91] 5.0 4.4 5.0 5.6 3.7 4.6 4.8 5.1 4.2 3.2 5.2 3.9 4.8 5.6 3.7 3.9 4.0 3.1
##  [109] 5.9 4.5 5.3 4.8 4.1 2.8 5.6 4.8 3.1 2.6 4.6 5.8 4.4 5.4 3.9 4.5 4.5 4.8
##  [127] 5.2 4.5 4.2 3.3 3.5 5.6 5.8 3.4 4.6 4.5 4.7 4.9 4.6 4.5 3.6 5.2 5.4 3.9
##  [145] 6.0 4.6 4.8 6.2 4.3 3.6 3.6 4.5 5.4 3.7 4.4 3.5 4.2 3.3 4.7 4.4 3.6 4.2
##  [163] 5.3 5.6 5.8 4.7 4.0 4.7 3.5 6.3 6.2 5.8 4.2 4.8 4.2 2.8 5.0 5.5 4.7 4.1
##  [181] 4.1 4.8 5.2 5.8 5.3 4.0 4.4 4.6 4.3 5.6 5.4 3.1 3.2 4.0 4.2 4.0 5.3 3.9
##  [199] 3.7 4.7 4.8 5.0 5.0 2.3 6.1 4.5 3.4 4.8 4.5 3.2 4.6 3.9 5.5 4.3 5.9 4.7
##  [217] 6.0 3.2 3.0 4.8 3.7 4.8 4.0 4.2 3.7 3.5 3.6 5.6 4.9 5.5 5.7 3.5 5.2 4.8
##  [235] 2.6 4.5 3.7 6.2 4.2 4.4 5.4 5.1 5.1 3.9 4.8 4.7 3.8 6.3 5.2 3.7 4.5 5.3
##  [253] 6.0 5.3 5.4 3.3 6.6 4.7 4.9 3.5 5.7 3.9 4.9 4.1 2.5 3.3 5.2 2.9 5.2 3.0
##  [271] 5.6 4.4 3.5 5.2 3.1 4.2 3.4 5.0 3.3 4.6 4.6 5.0 5.2 3.8 5.0 5.0 5.7 3.8
##  [289] 5.0 5.4 3.6 5.6 3.7 4.4 4.4 6.0 4.9 5.2 4.7 4.8 4.5 5.0 3.8 4.0 5.6 3.6
##  [307] 5.7 5.0 4.2 5.2 5.3 4.5 4.3 5.8 5.0 4.4 5.2 5.0 4.7 5.6 6.2 4.4 5.8 4.6
##  [325] 4.1 4.5 4.4 4.5 5.6 4.8 5.2 4.1 3.3 5.5 4.2 4.6 4.2 5.3 3.6 5.1 3.9 4.0
##  [343] 4.6 4.8 3.2 5.0 3.8 4.2 5.1 3.5 3.7 3.4 3.7 4.6 5.0 4.2 3.5 3.9 4.3 3.5
##  [361] 4.7 3.9 3.2 5.0 2.6 4.3 5.6 4.0 3.7 4.7 7.4 3.0 2.8 5.1 5.0 4.3 5.4 5.6
##  [379] 4.3 4.1 5.5 4.0 4.4 5.1 2.7 4.1 5.0 5.0 5.2 5.2 3.6 5.1 3.5 6.2 3.1 5.2
##  [397] 4.6 4.9 4.4 2.2 2.9 4.7 4.8 4.0 4.0 6.0 4.2 4.2 4.8 4.9 3.7 4.9 4.4 4.6
##  [415] 4.0 5.0 3.6 4.3 4.6 5.1 5.9 4.1 4.4 5.8 3.0 6.4 4.5 3.4 5.7 4.0 4.3 3.7
##  [433] 4.5 4.4 4.8 6.2 5.1 4.4 2.1 4.4 5.2 6.8 3.3 4.1 3.6 3.7 4.7 6.3 5.4 4.4
##  [451] 3.1 4.4 4.7 6.3 3.7 4.4 5.5 4.7 3.8 4.9 6.0 3.8 3.0 4.6 5.8 3.4 3.2 5.9
##  [469] 3.5 3.5 4.8 4.7 4.6 4.8 4.3 2.8 3.5 5.5 3.7 4.9 3.6 3.7 4.8 5.2 5.0 3.8
##  [487] 4.8 5.3 4.6 5.1 4.9 5.9 4.7 4.2 5.7 4.6 3.4 4.1 4.3 4.6 4.7 3.7 5.7 2.7
##  [505] 4.9 4.4 6.6 5.0 5.0 4.1 4.5 5.1 3.2 3.6 2.7 6.4 5.3 5.6 3.6 4.5 2.6 5.6
##  [523] 5.1 4.0 4.2 6.0 4.1 4.4 3.2 3.7 4.3 3.9 5.5 5.8 5.8 5.7 4.6 4.2 4.4 3.5
##  [541] 2.9 4.9 3.4 4.8 4.5 4.8 4.7 4.8 4.7 4.5 4.2 4.8 6.1 6.9 5.2 4.2 5.3 6.3
##  [559] 4.9 5.4 4.6 4.9 4.2 3.1 3.5 4.1 4.8 2.6 5.3 5.0 3.4 5.5 5.5 5.9 5.9 4.6
##  [577] 5.2 3.6 3.6 5.1 3.5 3.9 4.2 5.0 4.9 6.1 2.9 4.1 3.9 4.9 5.5 4.6 3.6 4.3
##  [595] 2.8 4.3 5.6 4.5 4.0 3.6 4.6 4.3 5.4 5.8 3.4 3.7 2.5 2.5 6.4 4.3 4.5 4.1
##  [613] 4.7 3.6 4.2 5.0 4.0 4.5 4.9 5.7 5.3 4.7 2.5 6.2 5.0 3.3 5.4 4.1 4.9 3.8
##  [631] 5.6 4.4 4.5 5.1 5.0 5.1 5.0 4.5 4.8 3.8 4.2 3.6 4.3 6.2 3.6 3.9 2.9 5.2
##  [649] 3.5 4.8 3.7 6.1 4.7 5.0 3.1 3.9 4.5 4.3 5.4 3.7 4.7 2.8 4.1 5.9 4.8 4.2
##  [667] 4.0 3.9 5.2 4.0 4.8 4.4 4.6 2.8 3.7 4.0 4.9 5.1 3.3 4.8 4.8 4.3 4.5 4.8
##  [685] 4.8 3.4 4.2 3.3 3.6 3.6 4.1 5.1 6.3 5.4 4.5 5.7 3.2 4.1 5.2 5.2 5.3 4.7
##  [703] 4.6 4.5 5.3 2.8 5.3 5.1 5.1 7.4 5.7 5.7 3.4 3.8 5.2 6.2 5.4 5.8 5.4 5.5
##  [721] 4.5 5.5 3.3 4.5 3.3 4.8 5.3 3.7 4.3 5.4 2.8 6.2 4.4 4.7 4.2 4.4 3.9 4.5
##  [739] 4.1 4.1 3.9 5.8 4.1 4.2 4.2 4.3 5.0 3.9 3.8 5.0 5.8 3.8 4.1 4.9 5.3 5.9
##  [757] 5.1 5.2 3.5 4.0 5.6 5.6 4.7 4.1 4.8 3.4 6.1 4.8 3.4 4.9 3.2 3.9 5.4 3.5
##  [775] 6.4 5.0 4.6 7.1 5.7 5.2 4.3 6.0 3.6 4.3 4.8 3.4 3.2 4.8 4.3 4.4 4.3 5.2
##  [793] 4.6 3.5 5.0 4.9 4.3 5.6 5.2 4.6 4.4 3.7 5.1 3.8 5.3 3.8 4.9 3.5 4.7 5.1
##  [811] 5.4 4.0 4.7 4.1 4.8 2.7 3.7 3.4 2.8 4.1 4.2 4.6 5.1 3.3 4.6 4.0 4.6 5.5
##  [829] 2.6 3.9 3.9 5.0 4.3 5.4 5.0 5.0 3.3 4.3 4.8 4.5 5.2 5.1 5.8 4.1 5.2 5.4
##  [847] 5.2 3.7 3.6 3.7 4.3 3.8 3.4 3.7 4.6 4.8 4.3 3.9 4.8 4.0 4.0 4.4 3.1 5.4
##  [865] 5.3 5.3 4.8 3.9 3.4 4.4 3.6 4.8 3.7 5.0 4.4 4.9 6.1 5.3 4.9 6.9 5.4 5.5
##  [883] 4.8 4.9 4.1 6.1 3.6 3.4 5.0 4.4 3.0 4.2 4.7 4.1 4.9 3.5 3.6 4.4 4.8 3.4
##  [901] 5.6 5.5 4.9 3.7 4.8 4.6 4.1 3.7 4.1 5.5 2.9 5.9 5.8 4.3 6.1 4.9 4.1 3.2
##  [919] 4.1 5.2 3.9 2.8 4.5 3.9 4.7 4.8 4.1 4.0 5.5 6.2 3.8 4.1 3.7 5.2 5.1 3.5
##  [937] 4.4 3.8 4.0 4.3 4.0 4.9 2.6 4.7 3.7 3.4 4.9 4.4 4.2 4.4 4.3 2.5 3.8 2.8
##  [955] 6.3 5.4 4.3 3.4 4.9 4.8 6.0 3.5 4.0 3.8 4.8 5.2 4.6 3.2 4.0 6.8 5.4 4.5
##  [973] 5.0 4.6 1.7 6.4 5.5 5.0 5.0 2.8 4.1 3.0 4.8 4.0 4.8 4.4 6.4 3.4 4.1 5.1
##  [991] 6.2 5.0 5.7 4.8 4.7 5.4 4.4 4.6 4.6 2.1
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
##   2.8   6.3
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
##    [1] 4.2 4.8 4.4 5.5 5.2 5.3 5.6 3.6 4.7 2.8 4.6 4.4 6.2 6.3 6.6 3.6 3.9 4.8
##   [19] 3.7 5.7 5.6 5.0 4.3 5.6 5.1 4.0 3.7 4.3 5.0 4.3 5.1 4.7 3.0 3.5 3.5 5.2
##   [37] 3.7 4.1 5.2 4.7 5.7 3.8 5.0 5.1 3.8 3.4 4.3 6.1 4.7 4.8 3.9 5.4 6.9 4.4
##   [55] 4.5 4.0 4.5 4.6 3.9 4.6 5.1 3.4 5.3 6.8 4.1 3.6 4.7 4.3 5.5 5.2 5.5 6.0
##   [73] 5.6 4.7 5.0 4.5 5.7 5.5 4.2 6.1 6.0 3.9 5.3 4.9 4.4 6.2 4.0 4.6 5.1 6.8
##   [91] 3.7 4.8 4.5 4.4 5.4 4.5 4.6 3.5 5.1 3.3 3.5 3.8 6.6 4.3 3.9 3.3 4.3 3.4
##  [109] 5.7 4.6 2.9 3.8 3.9 4.2 4.0 4.2 3.6 5.8 4.4 4.2 5.7 6.0 4.8 5.2 4.9 5.1
##  [127] 2.8 5.1 3.8 3.7 3.3 4.4 4.8 4.2 4.6 5.3 3.0 4.1 5.7 5.6 2.7 5.6 3.3 3.7
##  [145] 3.6 4.4 5.2 3.3 4.5 5.7 5.5 5.7 4.3 5.1 5.2 3.7 3.4 3.9 5.0 4.4 5.2 5.9
##  [163] 4.6 4.9 6.9 4.9 4.7 3.6 5.0 5.6 6.5 4.7 6.1 4.0 4.0 4.7 4.4 5.8 3.3 5.3
##  [181] 3.5 4.4 2.5 3.4 5.3 3.0 3.2 3.9 4.8 4.6 4.0 4.9 5.2 3.6 5.6 5.7 4.5 4.1
##  [199] 5.1 4.0 3.7 4.3 4.9 4.8 5.2 5.1 6.1 4.0 5.8 4.4 4.8 4.5 5.7 3.0 3.8 5.4
##  [217] 4.2 4.6 4.3 3.8 5.1 4.8 4.8 3.8 3.1 4.3 5.6 3.2 5.9 3.2 4.0 4.8 3.7 5.7
##  [235] 4.3 4.1 4.4 4.9 4.2 5.5 4.7 5.4 4.7 4.3 5.8 2.7 5.4 5.7 3.9 5.6 4.7 4.3
##  [253] 4.6 6.4 6.4 5.7 3.8 3.9 3.2 3.9 4.6 5.9 3.9 5.4 4.5 3.4 4.9 5.1 5.8 4.1
##  [271] 4.8 5.1 3.8 4.1 4.7 5.3 3.1 5.5 4.5 5.7 5.4 4.1 5.1 5.0 4.5 4.8 5.9 5.1
##  [289] 3.9 6.2 3.8 5.1 4.7 4.4 3.4 4.3 4.3 3.6 4.5 3.7 5.9 4.1 2.5 3.8 6.3 4.0
##  [307] 5.2 5.5 3.0 4.5 3.5 5.6 4.0 4.8 4.2 4.2 4.3 4.6 3.7 5.3 5.4 3.5 4.8 4.5
##  [325] 4.7 3.9 3.8 5.6 6.0 5.0 6.6 4.7 3.6 3.5 3.9 4.6 3.3 4.1 3.7 4.0 4.9 4.4
##  [343] 3.8 4.2 3.3 4.3 3.1 4.1 4.0 3.9 4.3 4.1 5.2 4.9 4.3 4.0 4.5 5.8 5.1 5.1
##  [361] 4.4 4.0 4.7 4.4 6.4 3.8 4.6 4.4 5.6 4.1 3.0 2.9 4.5 5.4 5.3 5.0 4.7 4.3
##  [379] 4.5 4.7 4.9 4.9 4.1 5.7 5.7 4.7 4.2 4.4 6.9 4.6 4.5 4.2 4.6 3.1 5.8 2.9
##  [397] 2.3 4.8 3.4 5.3 5.8 5.4 4.2 4.1 6.8 5.0 5.1 4.7 5.8 5.2 4.7 5.6 3.8 3.9
##  [415] 4.1 3.2 4.5 5.0 4.9 5.2 3.0 4.4 5.3 4.9 4.9 3.4 5.1 5.4 3.8 5.8 4.4 3.9
##  [433] 4.6 4.5 5.0 5.4 2.9 5.6 3.5 5.2 6.3 3.1 6.9 3.7 6.6 4.1 4.3 4.9 5.6 5.0
##  [451] 4.9 5.5 3.7 3.9 5.0 4.9 5.1 5.5 3.2 3.5 5.5 5.3 3.2 2.9 3.9 5.5 2.5 3.7
##  [469] 5.0 4.4 4.1 4.0 4.0 4.9 4.2 4.5 3.8 5.6 4.6 3.6 4.7 5.5 3.8 3.5 4.7 2.5
##  [487] 4.3 4.9 4.1 4.7 5.2 5.4 5.1 4.6 4.6 4.9 4.6 5.0 4.6 4.7 4.7 4.0 2.5 5.1
##  [505] 4.8 2.9 5.5 5.1 4.1 5.3 5.0 4.1 3.9 3.4 4.7 4.0 3.3 4.4 4.0 5.1 4.2 5.4
##  [523] 4.1 6.0 2.7 5.3 5.1 4.6 5.1 6.9 2.8 4.1 4.7 5.3 4.7 4.9 5.6 4.4 3.5 4.9
##  [541] 5.3 2.5 3.9 5.0 6.0 5.9 4.1 3.9 6.1 2.4 4.0 3.2 4.8 4.2 4.2 4.7 4.6 4.6
##  [559] 5.2 6.1 3.2 4.8 5.8 4.2 3.6 5.0 3.4 5.4 3.1 3.9 5.6 4.1 4.4 5.1 6.1 4.7
##  [577] 4.5 3.9 4.3 4.3 4.2 3.8 3.9 3.4 4.8 4.8 3.9 4.0 6.0 2.4 5.1 4.4 4.7 5.8
##  [595] 5.6 5.1 6.2 5.5 4.7 5.7 4.0 6.2 5.6 5.3 3.9 5.0 2.2 5.5 5.4 4.1 4.4 4.9
##  [613] 2.4 4.6 5.4 5.6 5.1 4.6 5.6 3.8 2.6 4.6 5.2 4.0 5.9 5.6 5.3 3.2 3.9 5.4
##  [631] 5.3 4.9 3.2 4.0 4.7 3.3 3.8 4.1 2.9 3.3 4.7 4.8 4.9 4.6 4.5 4.9 5.0 6.5
##  [649] 3.0 3.8 4.3 4.6 4.8 5.0 5.3 4.9 5.2 5.9 6.1 5.3 6.1 3.6 3.4 5.0 4.4 3.3
##  [667] 3.1 5.6 4.2 4.7 5.0 3.8 4.4 4.4 3.9 4.5 5.1 6.1 4.9 3.4 5.1 5.4 4.7 3.5
##  [685] 2.8 3.6 4.0 4.2 4.2 5.2 5.1 5.3 3.3 5.9 6.2 5.2 4.8 4.4 4.5 2.6 4.7 3.0
##  [703] 4.5 4.3 5.8 4.7 4.1 5.5 4.7 5.5 4.0 3.7 4.4 4.5 4.6 4.2 4.5 4.2 4.7 2.3
##  [721] 5.5 3.2 3.4 4.0 3.6 4.4 4.8 4.2 3.4 3.9 5.9 4.4 3.6 3.8 5.2 5.6 5.7 4.0
##  [739] 5.6 5.4 3.3 3.1 5.6 4.1 4.3 6.0 5.2 5.0 5.6 5.6 4.4 5.1 4.5 5.3 4.2 5.3
##  [757] 5.5 6.0 5.3 5.3 4.3 4.4 4.6 4.1 4.6 4.3 5.7 3.5 3.7 4.5 3.2 3.9 4.4 4.3
##  [775] 5.9 5.2 4.5 4.9 4.7 3.5 4.6 4.8 3.7 4.3 4.0 4.4 5.1 4.5 4.9 5.3 4.3 3.9
##  [793] 3.0 3.8 4.2 2.9 5.4 4.5 6.4 4.7 6.3 6.9 3.3 4.7 5.6 4.9 4.7 3.6 4.5 3.8
##  [811] 3.3 5.8 4.5 3.0 5.6 4.0 4.2 3.5 2.5 5.5 4.3 2.5 4.4 3.3 4.4 5.8 5.3 4.8
##  [829] 4.1 6.1 4.1 4.3 3.3 6.2 3.7 3.2 3.2 4.7 3.8 6.2 3.4 4.3 4.5 5.4 4.0 4.5
##  [847] 3.6 4.7 3.3 3.7 3.6 5.3 3.6 4.7 2.5 4.6 5.6 5.2 5.5 2.7 4.3 4.6 3.5 5.6
##  [865] 4.8 4.0 5.6 5.5 5.7 3.7 4.1 4.4 3.1 5.4 4.9 4.4 4.6 4.0 4.3 2.6 5.1 4.8
##  [883] 4.0 4.6 3.7 5.3 3.9 3.6 4.2 4.9 5.8 3.2 2.9 4.8 4.6 4.7 2.5 4.8 4.1 5.9
##  [901] 3.5 4.4 4.9 3.9 3.8 5.3 5.3 3.5 4.3 4.1 6.1 5.7 4.7 5.9 4.4 4.7 3.6 3.2
##  [919] 3.6 5.3 5.4 4.6 3.8 3.5 3.6 4.7 5.6 4.8 3.4 4.2 4.3 3.7 5.7 4.1 6.8 4.1
##  [937] 4.2 3.5 6.2 5.2 5.5 4.2 3.7 4.5 5.7 4.4 6.1 4.6 3.7 3.3 5.3 1.6 5.5 3.8
##  [955] 4.4 2.6 5.4 7.1 4.8 3.3 3.6 4.7 3.3 5.0 3.6 4.7 5.6 5.1 6.0 4.3 3.3 3.5
##  [973] 3.9 5.8 4.9 4.7 4.4 5.1 5.0 4.3 4.5 2.8 3.8 5.9 3.4 4.3 3.5 5.2 5.9 4.1
##  [991] 6.8 4.9 4.9 4.8 4.8 3.9 4.9 5.1 4.8 3.8
## 
## $func.thetastar
## [1] 0.0515
## 
## $jack.boot.val
##  [1]  0.56899441  0.40524862  0.32426667  0.25260274  0.05014663  0.02140845
##  [7] -0.02939481 -0.23524355 -0.39148936 -0.45575221
## 
## $jack.boot.se
## [1] 0.9725731
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
##    [1] 4.4 3.4 5.2 3.5 3.3 3.1 4.1 4.8 5.4 4.7 4.9 3.5 4.0 3.7 4.7 5.9 3.5 4.3
##   [19] 3.2 4.7 3.2 3.1 5.2 4.7 4.4 5.3 4.3 5.2 3.8 4.6 5.2 4.4 4.8 3.9 5.7 3.7
##   [37] 5.4 2.3 4.3 5.4 4.7 5.0 4.5 5.0 4.0 5.2 5.3 5.1 4.0 4.6 2.8 4.6 6.6 6.2
##   [55] 5.0 5.6 5.3 6.2 4.3 5.9 5.9 4.7 5.6 7.3 3.2 3.4 4.1 4.8 5.1 3.3 3.4 5.2
##   [73] 4.1 3.3 7.2 2.6 3.8 4.1 4.9 4.1 3.1 4.8 6.0 5.1 4.9 5.9 4.0 5.7 3.3 6.1
##   [91] 4.6 5.3 3.5 4.4 5.9 6.0 5.4 3.2 5.0 5.2 6.1 2.8 3.8 4.2 3.4 4.6 4.1 5.1
##  [109] 5.3 4.6 4.7 5.8 5.2 4.2 3.8 4.5 4.8 3.9 4.9 5.2 4.4 3.6 4.3 5.3 5.0 5.4
##  [127] 3.2 5.1 4.2 4.8 2.9 3.5 2.9 3.4 3.1 5.1 4.6 4.5 5.1 3.4 4.0 3.7 4.5 4.0
##  [145] 5.1 4.6 4.7 5.6 4.7 2.8 3.8 5.7 3.7 5.5 3.2 3.6 5.1 4.9 6.0 4.6 5.6 4.8
##  [163] 5.4 6.1 3.4 2.5 3.9 6.3 3.4 5.6 3.7 4.4 4.8 4.5 5.8 5.6 4.4 4.4 4.5 4.5
##  [181] 4.0 4.1 5.1 5.8 4.1 3.4 4.5 3.5 3.3 5.5 3.8 4.8 5.3 4.8 2.9 5.8 4.3 4.6
##  [199] 4.7 4.2 5.2 3.7 5.1 5.7 5.3 4.6 4.1 5.4 4.5 3.7 5.9 4.2 5.3 4.5 4.7 4.9
##  [217] 5.4 4.2 4.0 5.2 4.9 3.3 4.6 5.2 3.9 6.7 4.8 3.9 6.0 5.4 4.2 4.2 5.2 4.9
##  [235] 2.3 4.4 4.1 3.8 3.8 5.7 4.5 4.5 4.2 3.5 5.7 4.3 4.0 3.6 5.0 3.1 4.3 4.7
##  [253] 4.5 4.0 4.5 5.2 5.2 3.7 4.0 3.4 5.0 4.5 3.4 6.2 5.1 6.7 4.7 6.5 2.6 6.6
##  [271] 3.3 6.0 5.0 5.6 2.9 4.8 4.4 4.8 4.0 3.1 3.7 4.0 4.4 5.5 5.2 6.1 5.8 4.6
##  [289] 3.5 3.8 4.5 6.5 5.1 4.7 4.3 2.8 3.6 4.4 4.5 3.7 5.0 5.1 3.3 5.3 6.0 6.0
##  [307] 3.7 5.2 4.1 2.7 3.7 4.8 6.5 3.9 2.5 5.4 5.0 3.8 4.3 5.4 5.2 4.5 4.3 4.5
##  [325] 4.1 5.8 4.2 5.2 4.0 4.0 4.5 5.4 4.8 2.9 2.8 6.0 4.4 5.8 4.3 4.0 5.3 4.1
##  [343] 3.7 5.7 5.0 3.3 4.4 4.4 5.3 4.0 5.3 2.9 6.3 3.7 3.9 3.9 3.5 3.1 3.4 4.5
##  [361] 4.2 4.5 3.0 4.5 3.2 2.8 6.9 2.0 5.6 5.0 4.8 2.7 5.1 4.2 4.1 4.4 3.1 4.5
##  [379] 3.9 2.5 4.5 4.5 4.3 3.4 4.5 4.0 5.8 3.9 5.5 5.6 4.5 6.1 4.3 3.9 5.5 4.4
##  [397] 4.1 5.6 5.6 3.9 6.0 4.6 3.8 4.1 4.6 3.5 6.4 3.1 4.5 5.6 4.4 5.1 5.1 3.5
##  [415] 6.4 4.7 4.5 3.6 4.2 4.0 3.0 5.2 4.0 5.8 3.9 4.1 4.6 5.2 4.8 4.8 4.0 5.0
##  [433] 4.3 5.0 3.9 4.8 5.5 4.2 4.1 5.5 5.6 5.4 5.0 5.6 5.2 5.5 4.9 3.3 1.9 4.8
##  [451] 4.3 6.4 4.0 4.8 4.2 5.4 5.2 4.8 5.4 5.7 5.4 3.1 6.8 3.5 5.5 6.2 4.9 4.4
##  [469] 4.9 5.1 6.1 4.0 3.2 5.9 6.4 4.4 5.1 4.0 5.3 4.6 3.5 4.9 6.1 3.0 4.7 5.0
##  [487] 5.5 4.7 6.5 4.0 2.7 2.4 3.1 5.9 4.2 4.2 2.5 4.7 3.9 3.3 4.9 4.3 2.9 5.4
##  [505] 5.4 4.4 3.5 4.4 3.4 4.9 3.0 2.4 5.1 5.2 4.8 4.4 4.3 4.2 5.9 6.1 5.0 3.9
##  [523] 3.8 4.6 4.5 6.8 5.5 5.0 5.8 5.3 4.8 3.5 6.6 5.6 5.2 3.1 2.4 3.8 4.3 4.8
##  [541] 3.9 4.3 5.1 4.3 4.5 4.9 3.1 3.6 6.0 4.4 4.0 3.7 5.0 3.8 5.5 2.8 3.9 5.5
##  [559] 5.3 5.2 6.1 4.8 2.6 4.9 3.9 3.5 4.8 4.7 4.9 3.6 4.1 4.4 3.9 4.8 4.1 4.1
##  [577] 4.6 4.9 3.7 2.8 5.8 4.1 3.7 3.8 3.9 5.0 3.6 3.0 4.9 5.6 4.8 4.8 4.0 5.1
##  [595] 6.0 4.7 5.2 4.3 4.3 4.2 3.4 3.2 3.9 3.6 6.3 4.2 3.5 4.4 4.4 3.9 4.5 4.7
##  [613] 3.1 3.2 4.4 3.8 4.3 5.0 5.5 3.8 4.9 4.1 4.8 4.2 4.5 6.1 4.3 3.1 5.6 5.4
##  [631] 4.4 4.4 4.4 2.5 4.9 3.9 6.2 5.8 2.9 3.5 3.7 4.6 4.5 4.5 4.1 5.1 3.6 4.2
##  [649] 3.2 5.2 4.7 4.8 4.6 5.6 3.9 2.2 4.4 3.1 5.3 4.3 5.0 4.4 5.3 3.3 3.5 4.6
##  [667] 4.7 3.9 5.7 4.7 3.6 4.5 3.3 5.6 3.7 5.5 3.9 5.4 4.4 4.9 3.1 5.3 4.6 3.9
##  [685] 2.6 4.2 4.1 4.4 4.5 2.8 5.2 5.4 5.2 5.2 5.4 4.7 5.5 5.9 5.2 4.5 3.4 5.9
##  [703] 6.1 5.8 6.2 5.2 5.5 4.6 3.9 4.9 4.7 4.8 2.9 5.9 3.3 3.7 4.2 3.4 4.6 5.7
##  [721] 5.4 4.2 6.2 2.7 4.9 4.2 4.3 4.7 4.7 4.2 5.2 3.1 4.5 2.9 4.4 4.4 5.2 5.2
##  [739] 3.6 4.8 4.1 4.7 3.7 4.3 4.0 3.3 5.5 5.2 4.6 3.7 4.1 5.2 3.9 3.2 4.1 3.8
##  [757] 3.7 3.9 4.9 3.8 4.8 3.3 4.5 4.3 4.7 3.4 5.8 3.7 4.6 4.1 5.9 3.7 3.9 6.0
##  [775] 3.6 3.7 4.5 4.1 3.5 4.4 4.8 3.7 5.1 5.1 4.1 2.7 4.2 4.8 5.2 3.4 4.7 3.9
##  [793] 4.2 5.2 5.0 5.8 4.8 3.9 4.1 4.7 4.1 3.9 4.1 2.2 5.3 5.0 4.6 4.5 4.8 5.4
##  [811] 3.6 4.9 5.1 4.7 4.7 4.1 3.1 4.6 4.8 3.9 3.5 4.6 5.6 4.1 4.3 3.4 3.9 3.8
##  [829] 4.2 6.7 5.9 4.9 4.5 4.6 5.4 3.6 4.2 5.7 5.7 4.5 4.8 5.3 4.1 4.9 4.4 3.0
##  [847] 4.9 3.3 4.8 3.6 5.2 3.0 4.1 3.7 4.3 3.7 5.1 5.3 5.1 5.2 4.9 4.2 4.0 5.2
##  [865] 5.5 3.8 5.6 4.4 4.3 3.2 2.5 3.9 3.6 4.8 4.5 3.8 5.5 5.1 2.8 4.9 3.7 3.3
##  [883] 3.4 4.9 4.2 2.2 3.9 4.5 5.3 3.4 5.9 4.2 5.6 3.9 4.6 4.9 3.3 2.0 5.6 5.0
##  [901] 4.9 5.6 5.7 4.9 4.2 4.2 4.7 4.4 3.8 5.5 4.0 3.6 4.7 4.0 3.3 6.4 5.4 5.4
##  [919] 2.6 3.0 5.2 5.1 4.4 3.7 3.5 5.3 4.8 5.0 5.1 4.1 4.4 6.3 4.8 5.3 5.5 4.4
##  [937] 4.8 5.5 3.6 3.4 5.4 4.5 4.8 5.4 4.0 3.2 3.2 3.2 5.4 4.3 3.1 5.7 4.2 2.8
##  [955] 5.5 1.3 4.7 4.4 3.4 3.3 4.3 3.6 3.1 4.4 3.5 5.0 4.4 5.4 6.0 3.5 3.2 4.5
##  [973] 4.5 4.4 3.7 4.0 3.8 6.5 4.8 4.2 3.3 4.1 6.5 3.8 4.5 4.9 5.3 4.5 4.1 3.7
##  [991] 4.2 4.5 4.3 5.0 4.3 5.6 4.2 3.8 3.1 3.6
## 
## $func.thetastar
## 72% 
## 5.1 
## 
## $jack.boot.val
##  [1] 5.5 5.5 5.3 5.2 5.2 4.9 4.8 4.8 4.5 4.4
## 
## $jack.boot.se
## [1] 1.11
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
## [1] 0.6755591
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
##   4.813467   8.374548 
##  (2.082338) (3.818679)
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
## [1]  0.4390205  1.2348854  0.7241190  1.0870628 -0.3188306  0.4499389
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
##    [1]  3.741683e-01 -2.792834e-02 -7.291974e-02  4.561697e-01  8.450677e-01
##    [6]  7.504578e-02  6.230779e-01  9.772638e-01 -5.068049e-02 -9.183743e-03
##   [11] -6.341813e-01  1.737107e-01  1.187919e+00  4.200301e-01  3.038738e-01
##   [16]  4.396914e-01  1.384605e+00 -1.263293e+00  7.045854e-01  2.109139e-01
##   [21] -1.830407e-01  4.365923e-01 -6.539971e-02  9.320809e-01  3.739354e-01
##   [26]  8.527961e-01  1.074502e+00 -7.187732e-01  5.914176e-01  1.098156e+00
##   [31] -4.214522e-02 -1.029652e-01  2.936817e-02  5.488508e-01  7.863232e-01
##   [36]  1.050098e+00  4.383329e-01 -7.030954e-01  6.755591e-01  2.615758e-01
##   [41] -4.715067e-01  6.701140e-01  4.023917e-01  1.194713e+00  6.244787e-01
##   [46]  1.953497e-01 -1.847017e-01  7.833043e-01 -3.900629e-01  4.291596e-01
##   [51] -5.150186e-01  1.483044e-02  3.293727e-01  1.087837e+00  4.634807e-01
##   [56]  3.213548e-01  1.009409e+00 -4.495710e-01 -4.251071e-01  4.972704e-01
##   [61]  9.154516e-01  2.968299e-01 -1.206562e-01  1.887780e-01  6.336097e-01
##   [66]  2.218807e-02  8.491086e-01  6.593859e-01  3.854478e-01  8.318814e-01
##   [71]  8.682754e-01  8.988968e-01  4.272210e-01  8.542370e-01  3.564630e-01
##   [76]  7.237179e-01  9.570033e-01  8.961339e-01  4.529882e-01  1.339671e+00
##   [81]  4.915757e-01  1.175694e+00  9.287484e-01  4.077390e-01  7.467876e-02
##   [86]  2.612523e-01  7.796008e-01 -5.448209e-01 -1.974683e-01  1.069846e+00
##   [91]  6.857325e-01 -3.550117e-01  6.099600e-01 -9.678940e-01 -2.145018e-01
##   [96] -1.614815e-01 -1.580849e-01  9.130663e-01  5.227295e-01  3.190809e-01
##  [101]  9.825136e-01  1.562256e-01  3.859226e-01  1.013463e+00 -2.665407e-01
##  [106]  8.128943e-01  6.220514e-01  2.379973e-01  6.989144e-01 -6.630999e-02
##  [111]  7.525102e-01  8.687915e-01  1.391397e+00  1.003091e+00  9.788481e-01
##  [116]  2.853386e-01  3.217244e-01  4.927252e-01  1.887794e-01  1.374362e+00
##  [121]  1.140569e+00  8.444506e-01 -2.676133e-01  6.256706e-01 -4.901970e-02
##  [126] -5.016998e-01  3.670127e-02  4.344813e-01  5.853262e-01  5.664001e-01
##  [131]  1.621885e+00  9.175249e-02  2.588654e-01 -2.968635e-01  5.580767e-01
##  [136]  1.418332e+00  3.113887e-01  5.774364e-01  6.339336e-01  1.287988e+00
##  [141]  7.992898e-01  6.214430e-01 -1.099447e-03 -1.046829e+00 -4.432357e-02
##  [146]  1.656261e+00 -1.828722e-01  1.282086e+00  8.471729e-01  5.909342e-01
##  [151]  2.859519e-01  8.161750e-01  6.946724e-01  1.024870e+00  7.346933e-01
##  [156]  1.440317e+00  3.260195e-01  2.619075e-01  2.225446e-01  6.974901e-01
##  [161]  9.930175e-01  9.015214e-01  5.408911e-01  4.354086e-02 -1.847465e-05
##  [166]  6.350325e-03  1.455820e+00  6.857873e-01  1.629895e-01  1.250508e+00
##  [171]  1.370872e+00  8.361334e-01  1.809881e-01 -3.125336e-01 -4.086334e-01
##  [176]  6.588894e-01  4.878989e-01  7.313753e-01  2.590098e-01  1.255734e+00
##  [181]  4.900075e-01 -2.676133e-01  2.078061e+00  6.348916e-01  9.084717e-02
##  [186]  5.332316e-01  1.571755e-01  2.170659e-01  2.397974e-01  3.438989e-01
##  [191]  1.510880e-01 -6.810816e-01  3.653463e-01 -3.243891e-02  2.079824e-01
##  [196]  9.718185e-02 -1.827022e-01  3.881251e-01  9.831144e-01  6.170244e-02
##  [201]  1.009445e+00  6.543372e-01 -7.999348e-02  5.417520e-01  1.012077e-01
##  [206] -2.682367e-02  2.449424e-01  1.467426e+00  6.265696e-01  1.688134e-01
##  [211] -4.178540e-02  1.616559e-01  1.214023e+00  5.097786e-01  5.249725e-01
##  [216]  6.828369e-01  2.744823e-01  1.698541e-01 -1.697251e-01  2.822871e-01
##  [221]  1.202734e-01  8.889443e-01  2.712450e-02  1.394556e-01  1.133307e+00
##  [226]  1.059247e+00  7.169014e-01  3.839812e-01 -3.419391e-01 -1.330115e-02
##  [231]  4.959115e-01 -8.281945e-01 -6.872829e-01  7.440640e-02 -6.427881e-02
##  [236]  8.556170e-01  2.281132e-01  1.086949e-01  7.332911e-01 -3.532781e-01
##  [241] -6.985054e-01  5.375990e-01  1.785430e-01  5.315409e-01  7.893197e-01
##  [246]  6.964159e-01  1.184984e+00  5.271227e-01 -5.515739e-01 -3.180335e-01
##  [251]  3.740087e-01 -8.208971e-01  2.657290e-01 -1.974639e-01  9.520511e-02
##  [256]  9.036069e-01  2.647263e-01  3.840497e-01  1.837208e+00  6.210680e-01
##  [261] -4.879261e-01  6.356191e-01  6.011325e-01  4.913427e-01 -9.939929e-01
##  [266]  5.362020e-01  2.865551e-01  7.168970e-01 -4.032831e-01  9.281529e-01
##  [271]  2.501733e-01  8.550426e-01  6.688785e-01  4.724488e-01 -9.235940e-03
##  [276] -5.041068e-01  5.827545e-01 -2.650879e-01 -2.245542e-01  8.975711e-01
##  [281]  4.977473e-01  1.425413e+00 -4.348507e-01  7.943047e-01  3.172596e-01
##  [286]  4.005266e-01  9.960169e-01  8.610439e-01  1.296150e+00  1.397180e+00
##  [291] -4.270515e-01 -3.107771e-01  3.355641e-01  5.288519e-01  2.159059e-01
##  [296]  1.688291e-01  6.071179e-01  3.609345e-01  3.714329e-01  7.027975e-01
##  [301]  4.025930e-01 -3.587552e-01  9.782512e-02 -3.858461e-01 -1.129517e-01
##  [306] -4.132761e-02  8.038594e-01  7.201131e-01 -1.170622e-01  3.280413e-01
##  [311] -5.041068e-01  4.172509e-01  5.218335e-01  4.441863e-01  1.467622e-01
##  [316]  1.871696e-02 -2.755695e-02  7.743003e-01  8.105864e-01  9.659963e-01
##  [321]  8.355927e-02  8.103048e-01  5.493321e-01  1.282329e+00  4.660072e-01
##  [326]  2.414118e-01 -5.771602e-01  4.068255e-01  1.342179e-01  2.244822e-01
##  [331] -1.277318e+00  1.230082e+00  9.684198e-01  9.213563e-01  1.378430e+00
##  [336]  6.003023e-02 -3.202687e-01  7.871904e-01  1.236479e+00  6.525003e-01
##  [341]  2.013663e+00  6.723840e-01  1.156541e-01  1.215406e+00  2.012333e-01
##  [346]  7.463519e-01  6.674596e-01  2.402926e+00  6.419946e-01  1.712880e-01
##  [351]  6.677978e-01  1.221121e+00  1.412754e+00 -4.852360e-01  6.780999e-01
##  [356]  4.163373e-01 -2.377049e-01  1.265091e+00  5.099371e-01  3.518698e-01
##  [361]  1.186233e+00  8.226584e-01  1.017863e+00  9.613937e-01  8.824538e-02
##  [366] -1.170622e-01  9.933654e-01 -1.818715e-01  2.943763e-01  3.477599e-01
##  [371]  2.131674e-01  5.838732e-01  1.284851e+00  1.075349e+00  7.546864e-01
##  [376]  4.564530e-01 -5.370330e-01  1.630390e+00  1.022322e-01  4.074519e-01
##  [381] -4.035474e-01  2.720485e-01  4.416086e-01  1.362103e+00  9.605478e-01
##  [386]  7.045854e-01  9.376681e-01  1.280149e-01 -5.914928e-01  1.579698e+00
##  [391]  1.152612e+00  2.335903e-01  2.967486e-01  2.602376e-01  8.139291e-02
##  [396]  2.346912e-01  3.290458e-01  1.480554e+00  6.207436e-01  5.574261e-01
##  [401]  2.094610e-01  9.369166e-01  8.436029e-02  3.961853e-01  6.132372e-01
##  [406]  6.391865e-01  5.613777e-01  1.910632e-01 -5.051010e-01  5.458328e-01
##  [411]  7.471875e-01  1.777690e+00  9.991455e-02 -2.115162e-02  8.462555e-01
##  [416]  9.242622e-01  5.969293e-01  1.321414e+00  5.408082e-01  6.712167e-01
##  [421]  1.920168e+00  2.891305e-01  6.191797e-01  5.279655e-01 -8.248445e-02
##  [426]  1.071271e+00  4.319906e-01  6.528639e-01  6.316809e-01  4.168275e-01
##  [431]  5.034531e-01 -3.834578e-02  3.397873e-01  8.549806e-01  3.677237e-01
##  [436] -8.297894e-02  7.244551e-01  7.010202e-01 -7.550391e-01  7.126880e-01
##  [441]  4.163140e-01  7.043045e-01 -1.686710e-01 -5.243296e-01  8.654604e-01
##  [446]  6.404339e-01  8.449655e-01  1.717320e+00  7.646708e-01  9.782512e-02
##  [451]  3.638927e-01  8.960370e-01  9.570033e-01  5.871136e-01  1.190590e+00
##  [456]  9.957430e-02  1.495374e+00  3.237366e-01  3.714481e-01  5.992202e-01
##  [461]  2.385837e-01  1.789088e-01  3.136452e-01  4.958826e-01  8.536198e-01
##  [466]  3.257243e-01  6.327609e-02  4.690362e-01 -7.565739e-01 -5.555270e-01
##  [471]  5.928510e-02 -8.697704e-02 -3.437223e-01  4.982282e-01 -3.575415e-01
##  [476]  9.130487e-01 -8.445385e-01  9.830914e-02  1.981915e+00  1.585446e-01
##  [481] -2.423604e-01  3.384330e-01  1.307711e+00  9.335546e-01  5.288092e-01
##  [486]  7.224663e-01 -3.748279e-01  5.340710e-01 -9.515251e-02  5.919405e-01
##  [491]  1.048568e+00  1.395237e-01 -5.156791e-01  1.069846e+00  3.489136e-01
##  [496]  1.306626e+00  1.015116e+00  6.787463e-01  1.029063e+00  5.047722e-01
##  [501]  2.689427e-01  1.202734e-01  1.069979e+00  1.158319e-01 -2.433477e-01
##  [506]  3.618083e-01  1.344023e-01  3.400055e-01  8.311968e-01  3.621040e-01
##  [511]  1.424693e-01  7.994589e-01  5.218306e-01  8.230472e-01  6.637480e-01
##  [516]  7.937688e-01  1.062804e+00  6.033583e-01  4.296197e-01  4.453594e-01
##  [521] -5.282511e-01 -9.042203e-01  8.630502e-01  9.344724e-02  9.965878e-03
##  [526] -1.481636e-01  8.027203e-02  2.165063e-02  8.242488e-01 -1.782804e-01
##  [531] -1.710868e-02 -1.876603e-01  6.832408e-01  8.627611e-01  2.293759e-01
##  [536]  5.505641e-01  6.168342e-02  7.289925e-01  2.966816e-01  7.618903e-01
##  [541] -2.864992e-01  6.800724e-01  1.213533e+00  9.933654e-01 -6.078122e-02
##  [546]  1.517644e-01  6.042356e-01  6.150812e-01  4.296534e-01  6.503385e-01
##  [551]  3.719472e-01  3.659950e-01 -9.964940e-01  6.666519e-01  3.513390e-02
##  [556]  4.397296e-01  6.834294e-01  3.957954e-01  7.115380e-01  7.904562e-01
##  [561]  5.157779e-01  5.157876e-02  1.085147e+00  8.379339e-01  6.691360e-01
##  [566]  5.155059e-01  3.473542e-01  4.616129e-01  5.527396e-01  7.400910e-02
##  [571]  5.449577e-01 -4.219867e-02  4.729063e-01 -2.021249e-01 -5.197378e-02
##  [576] -6.378522e-02  1.398077e-01  1.095839e+00  9.048957e-01  1.727533e+00
##  [581] -4.972510e-01  1.444987e+00  1.200852e+00  8.212023e-01  1.756487e+00
##  [586]  1.774557e-01 -1.822906e-01  1.668747e-01 -7.176177e-01 -5.050289e-02
##  [591]  6.944558e-01  2.411169e-01  8.311968e-01  6.736078e-01  9.895242e-01
##  [596]  4.313637e-01  6.650695e-01  2.261361e-01  2.228487e+00  1.570547e+00
##  [601]  1.115959e+00  4.478922e-01  3.126973e-01  4.010500e-01  5.640008e-01
##  [606]  4.863252e-01  5.219103e-01 -4.939952e-02  3.169904e-01  1.334184e+00
##  [611] -3.538479e-01  5.400648e-01 -7.635182e-01  2.651878e-01 -4.845255e-02
##  [616]  7.569883e-01 -4.703624e-01  1.007199e+00  8.139034e-01  4.351873e-01
##  [621]  4.340733e-01  1.166559e+00  7.791383e-01  7.997347e-01 -3.553513e-01
##  [626]  7.852764e-02  2.958668e-01  1.447412e+00 -1.005726e-01  1.144936e+00
##  [631]  4.760084e-01  3.699153e-01 -4.886903e-01  4.040421e-01  5.292191e-01
##  [636]  9.509261e-01  9.411820e-02  7.849213e-01  4.823857e-01  7.009359e-01
##  [641]  8.492364e-01  9.943604e-02  5.433945e-01  1.529806e-01 -5.596682e-02
##  [646]  4.989684e-01 -6.516003e-03  1.351012e+00  8.847330e-01  6.800724e-01
##  [651]  3.253556e-01  4.153835e-01  1.370810e-01 -2.830306e-01  2.485980e-01
##  [656]  3.445961e-01 -5.151718e-01  4.335802e-01  4.672459e-01 -8.801097e-02
##  [661]  5.148876e-01  3.854131e-01  4.429949e-01  4.049397e-01  3.182815e-01
##  [666]  1.072301e+00 -1.897782e-01  2.335533e-01  1.069875e-01  1.074891e+00
##  [671] -5.228737e-02  1.038748e+00  1.746951e+00  4.846285e-01  5.891443e-01
##  [676]  5.034531e-01 -2.081967e-01  4.766630e-01 -2.826230e-01  4.647343e-01
##  [681]  1.086146e+00  9.480822e-01  3.041450e-02 -3.324055e-01  5.938099e-02
##  [686] -3.911010e-01  6.888625e-01  1.462650e+00  5.855756e-01  1.443146e-01
##  [691]  4.635876e-01 -2.559662e-01 -3.543818e-01  5.724768e-01  1.922045e-01
##  [696]  3.280216e-02  5.339513e-01  7.883767e-01  7.013036e-01  4.917470e-01
##  [701]  1.008740e+00  3.531833e-01 -2.815770e-02  1.118450e-01  2.059254e-01
##  [706]  5.201442e-01  6.941505e-01 -2.466065e-01 -2.948529e-01 -9.707169e-01
##  [711]  1.031648e+00  8.900590e-02 -1.317966e-01  5.175360e-01  1.271497e+00
##  [716]  1.280234e+00  7.046599e-01  1.179509e+00 -1.013407e-01  3.148404e-01
##  [721]  5.551867e-01  6.405676e-01  3.099405e-01  9.876402e-01  2.868119e-01
##  [726]  6.885037e-01  3.108677e-01  4.276572e-01  7.148519e-01  1.269852e+00
##  [731] -5.129628e-01  9.140288e-01 -1.146880e-01  4.177429e-01  1.154503e+00
##  [736]  1.495732e-02  9.270910e-01 -7.268646e-01 -6.756292e-01  7.548146e-01
##  [741]  7.538704e-01 -5.097542e-01  6.409904e-01  5.986041e-01  6.630188e-01
##  [746]  7.341030e-01  4.960651e-01  5.381882e-01 -3.720345e-01  1.591098e+00
##  [751]  9.658743e-01 -5.748027e-01  1.049721e-01  3.299350e-01  5.367069e-01
##  [756]  1.676888e+00 -2.047275e-01  8.606354e-01  5.959730e-01  5.052405e-01
##  [761]  8.414556e-01  1.483287e+00 -2.627788e-01  7.920520e-01 -5.995266e-02
##  [766]  3.406444e-01  9.567503e-01  1.629982e+00  2.467096e-02  6.705957e-01
##  [771]  8.523717e-01  7.277015e-01  6.832408e-01  5.935233e-01  5.928873e-01
##  [776]  3.393680e-01  5.962032e-01  1.063398e+00  3.259141e-01  1.124095e+00
##  [781]  4.650802e-01  7.824714e-01  8.167055e-01  6.926285e-01  8.442404e-01
##  [786]  4.420193e-02  1.826364e-01  1.742788e-03  7.727990e-01  1.360689e+00
##  [791]  1.067816e+00  1.857162e-01  5.057329e-01  8.979848e-02 -2.641898e-02
##  [796]  4.288450e-01  6.357239e-01  1.178538e+00  4.739456e-01  1.287605e-01
##  [801]  7.582091e-01  1.063247e+00  2.304596e-01  1.025270e+00  8.172880e-01
##  [806]  6.707517e-01  2.021552e-01  9.695663e-02  9.357620e-01  3.854352e-01
##  [811]  5.320578e-01  1.147138e+00  6.671448e-01  1.252557e-01  2.083959e+00
##  [816]  1.327219e+00  5.399147e-01  2.265780e-01  5.584461e-01  1.157545e+00
##  [821]  6.573588e-01  5.161506e-01  2.332534e-02 -4.722667e-01  8.434539e-02
##  [826]  7.553145e-01  1.510912e-01  9.717610e-01 -2.354268e-01  7.642813e-01
##  [831]  7.240874e-01  1.000189e+00 -3.886011e-01 -1.528816e-01  4.812594e-01
##  [836]  4.883982e-01  8.091993e-01  9.551676e-01  1.264272e+00  8.923306e-01
##  [841] -1.039848e-02  4.893771e-01  7.638688e-01  9.063101e-01 -6.953961e-02
##  [846] -2.457130e-01 -3.490125e-01  8.202821e-01  6.560479e-01  6.713645e-01
##  [851]  2.135653e+00 -9.145382e-02 -1.142907e-01  9.636903e-01 -4.858524e-01
##  [856] -6.516003e-03  7.104585e-01  4.473782e-01 -1.084643e-01  6.292802e-01
##  [861]  5.878700e-01  1.155849e-01 -3.031884e-01  6.874340e-01  7.636474e-01
##  [866] -4.953233e-01  2.190390e-01  4.504308e-01  5.717967e-01 -6.647181e-02
##  [871]  1.129259e-01  1.701632e-01  3.304503e-01  1.048713e+00 -3.207170e-01
##  [876]  1.813905e-01  4.406729e-01  2.220775e-01  5.944264e-01  5.853492e-01
##  [881]  9.255736e-01  1.479288e-01  8.800532e-01  4.897135e-01  1.187084e+00
##  [886]  5.797227e-01  5.039744e-02  8.226095e-01  6.589236e-01  1.051370e+00
##  [891]  1.033897e+00 -2.521396e-01 -2.536846e-01  8.358111e-01  5.555255e-01
##  [896] -2.757208e-01  1.515464e-01  1.058105e+00  8.848569e-01  1.654468e+00
##  [901]  6.524464e-01  1.112289e+00  3.745158e-01  5.724768e-01  6.011984e-01
##  [906]  2.357891e-01  8.961062e-01 -6.868595e-01  1.912486e-01  2.083584e-01
##  [911]  2.109648e-01  4.932838e-01  4.730130e-01  3.950414e-01  6.832662e-01
##  [916] -2.986808e-01  7.235933e-01  6.257274e-01  1.011359e+00  8.287908e-01
##  [921]  4.283729e-01 -2.730249e-01  3.245405e-01  6.800724e-01  1.962497e-01
##  [926]  1.201050e+00  8.053811e-01 -4.001643e-01  7.682166e-01  4.861383e-01
##  [931]  1.221003e+00  1.750821e+00  3.284586e-02  5.172580e-01  6.904698e-02
##  [936]  4.708921e-01 -1.909131e-01  2.659446e-01  1.089766e+00  7.763618e-02
##  [941]  7.815533e-01  3.454069e-01 -2.392251e-02 -3.777818e-01  2.813150e-01
##  [946] -3.175336e-01 -3.130396e-01  6.657847e-01  7.957082e-01  3.251379e-01
##  [951] -1.712111e-01 -8.606290e-01 -9.792115e-02  9.759621e-01 -1.005357e+00
##  [956]  5.303648e-01  3.798296e-01 -8.281945e-01  4.241137e-01  3.220551e-01
##  [961]  3.188969e-01  8.866765e-01  6.701297e-01  1.141073e-01  7.562483e-01
##  [966]  1.276157e-01 -4.603856e-01 -2.386552e-01  4.611433e-01 -3.001863e-01
##  [971]  1.324778e+00  1.475614e-01  6.828369e-01  6.755591e-01  1.936230e-01
##  [976]  8.659402e-01  5.961611e-01  6.551778e-01  8.048500e-01  4.940115e-01
##  [981]  9.193157e-01  1.000958e+00  7.016487e-02  7.503392e-01  1.263649e+00
##  [986] -1.512117e-01  4.870910e-01  7.889484e-01  9.410711e-01  8.202913e-01
##  [991]  3.991357e-01  6.667742e-01  7.506094e-01  1.700360e+00  3.590446e-01
##  [996]  8.993849e-01 -5.544123e-01  5.084850e-01  1.009625e+00  4.168536e-01
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
##       mean          sd    
##   0.57477277   0.26416560 
##  (0.08353650) (0.05906648)
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
## [1] -0.1023660  0.4629630 -0.1500307 -0.5043800 -0.2183728  0.6341872
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
## [1] 0.0188
```

```r
sd(results2$thetastar)  #the standard deviation of the theta stars is the SE of the statistic (in this case, the mean)
```

```
## [1] 0.9182209
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
## t1*      4.5 -0.007907908   0.9031605
```

One of the main advantages to the 'boot' package over the 'bootstrap' package is the nicer formatting of the output.

Going back to our original code, lets see how we could reproduce all of these numbers:


```r
table(sample(x,size=length(x),replace=T))
```

```
## 
## 1 3 4 5 6 7 9 
## 1 2 1 2 1 2 1
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
## [1] 0.0141
```

```r
se.boot
```

```
## [1] 0.9172309
```

Why do our numbers not agree exactly with those of the boot package? This is because our estimates of bias and standard error are just estimates, and they carry with them their own uncertainties. That is one of the reasons we might bother doing jackknife-after-bootstrap.

The 'boot' package has a LOT of functionality. If we have time, we will come back to some of these more complex functions later in the semester as we cover topics like regression and glm.

