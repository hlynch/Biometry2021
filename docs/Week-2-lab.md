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
## 0 1 2 3 7 8 9 
## 1 1 2 1 1 2 2
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
## [1] 0.0014
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
## [1] 2.707174
```

```r
UL.boot
```

```
## [1] 6.295626
```

Method #2: Simply take the quantiles of the bootstrap statistics


```r
quantile(xmeans,c(0.025,0.975))
```

```
##  2.5% 97.5% 
##   2.5   6.3
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
##    [1] 4.5 3.7 3.9 3.8 6.1 2.6 3.7 6.1 4.5 4.1 6.2 4.1 4.5 4.1 5.1 4.5 3.1 5.0
##   [19] 4.3 5.1 5.6 4.5 5.0 4.6 3.6 5.4 3.7 5.4 5.8 5.6 4.2 4.5 5.1 6.1 5.4 3.6
##   [37] 5.0 4.4 5.4 4.0 4.5 5.0 5.2 4.1 4.9 5.4 2.9 3.6 2.4 3.8 4.5 5.0 4.3 6.4
##   [55] 5.0 5.1 4.4 5.3 4.1 3.0 4.4 2.6 5.0 3.1 4.6 4.9 4.3 2.6 4.6 6.4 4.9 4.9
##   [73] 4.0 6.7 3.6 5.6 4.0 4.4 3.3 4.2 4.6 4.6 5.1 3.8 4.4 3.0 3.6 4.3 4.7 3.6
##   [91] 5.9 4.3 5.7 5.0 3.3 3.6 5.1 6.5 5.1 4.5 3.5 4.6 2.7 5.5 5.3 6.6 4.6 4.2
##  [109] 3.4 5.8 3.4 5.0 5.2 5.0 4.7 4.1 2.6 4.1 3.8 3.9 5.8 5.4 3.5 5.1 4.6 4.3
##  [127] 5.0 5.4 5.5 4.4 3.7 3.8 5.3 4.1 3.6 4.4 5.8 3.7 4.1 2.9 5.3 4.7 4.1 4.7
##  [145] 5.3 5.4 3.8 5.7 4.4 4.4 4.0 3.3 5.4 3.4 5.9 5.4 5.4 5.5 4.6 3.0 4.2 3.6
##  [163] 3.8 4.9 3.9 6.9 4.9 4.2 4.9 5.2 2.9 5.4 4.8 5.9 4.3 5.7 3.6 5.0 5.1 3.3
##  [181] 5.4 4.1 5.3 6.0 4.9 4.9 4.4 4.5 4.8 4.0 4.9 3.3 5.1 3.0 3.8 4.2 5.0 3.4
##  [199] 4.8 3.0 4.6 5.8 6.2 2.1 5.5 2.4 4.7 5.0 3.9 2.5 4.2 3.5 4.3 4.4 2.7 6.0
##  [217] 5.4 4.6 4.9 3.5 3.0 5.2 4.3 5.1 4.0 4.1 4.7 3.1 5.7 4.9 5.1 3.0 5.0 5.1
##  [235] 5.6 6.1 4.4 3.9 5.5 7.5 6.1 3.9 4.4 4.7 4.9 4.2 2.8 3.5 3.6 3.1 4.8 3.7
##  [253] 4.3 6.1 5.3 3.9 4.5 4.1 3.8 3.4 2.4 4.9 4.0 5.2 4.4 2.6 3.8 4.3 4.3 4.2
##  [271] 4.2 3.7 3.4 3.8 4.5 4.4 3.3 5.3 3.8 5.8 4.2 4.5 5.8 3.7 4.3 4.9 4.0 3.3
##  [289] 3.4 4.4 4.8 4.3 4.5 4.9 3.6 5.3 4.2 5.3 4.5 4.3 5.8 5.4 4.0 4.6 7.2 5.6
##  [307] 4.5 4.5 4.0 4.2 3.6 4.6 5.6 5.3 2.8 3.9 4.9 5.7 4.7 4.9 4.3 4.1 5.9 2.8
##  [325] 4.6 5.0 5.0 4.6 5.8 5.3 3.4 4.1 4.6 5.3 5.1 4.6 3.2 5.7 5.7 2.9 3.2 4.3
##  [343] 2.7 3.7 5.2 3.5 5.6 4.5 3.1 4.9 3.6 5.6 5.2 4.2 5.0 4.2 2.8 4.4 3.6 4.0
##  [361] 3.3 4.4 4.3 4.4 4.2 3.9 4.8 4.5 5.1 4.0 4.5 5.0 4.0 6.2 4.0 4.6 4.7 4.2
##  [379] 4.1 4.1 4.8 5.0 4.1 4.9 4.5 4.5 5.2 4.5 4.2 4.2 4.5 4.2 5.7 2.6 4.3 4.9
##  [397] 4.4 4.3 3.6 4.9 5.3 4.2 4.4 3.0 4.1 2.5 4.2 5.3 4.6 3.7 4.9 3.9 4.7 4.5
##  [415] 4.9 3.6 6.3 4.6 3.5 4.2 4.9 4.8 3.9 3.8 5.1 3.9 4.7 4.8 5.3 4.1 3.3 4.7
##  [433] 3.5 4.7 3.8 4.1 3.3 5.5 3.9 4.4 5.6 3.6 5.8 4.0 2.9 5.3 5.2 5.1 4.2 4.8
##  [451] 4.7 5.4 5.5 3.7 5.2 4.3 3.0 2.4 3.8 4.9 4.2 5.6 5.5 4.6 4.4 4.8 4.9 4.8
##  [469] 5.1 4.0 7.2 4.3 5.2 4.4 5.5 4.9 5.4 4.2 6.1 4.0 5.8 5.6 4.5 5.1 4.9 5.5
##  [487] 5.2 4.9 4.4 3.4 6.3 5.4 4.6 3.6 4.6 5.4 4.0 4.7 4.6 5.9 5.5 5.5 4.2 3.9
##  [505] 3.8 5.8 4.4 3.1 4.4 5.7 4.9 5.6 4.6 5.4 2.6 6.1 6.4 4.4 4.0 5.0 4.7 4.3
##  [523] 4.7 3.5 5.1 3.6 4.9 4.7 3.4 4.3 5.2 4.0 5.8 6.0 5.5 4.0 5.8 5.1 4.9 4.0
##  [541] 3.9 4.1 4.7 3.5 5.4 5.1 3.2 3.6 4.0 4.4 4.4 3.2 5.2 4.7 5.0 4.0 5.2 4.6
##  [559] 4.5 4.3 5.5 5.1 4.6 3.3 6.6 3.9 4.1 4.5 3.2 4.7 4.7 3.8 5.3 6.8 4.1 4.4
##  [577] 4.1 3.3 4.7 4.5 4.4 4.5 4.5 6.2 4.6 3.3 4.9 4.8 5.0 4.4 2.8 5.2 4.9 6.0
##  [595] 4.4 6.0 4.1 4.8 5.3 4.8 3.0 4.6 6.5 4.1 2.5 3.1 4.2 4.8 4.7 4.6 4.6 3.4
##  [613] 4.7 5.7 4.6 4.3 4.1 6.0 5.4 5.3 4.1 3.8 4.9 4.2 2.8 3.3 6.1 5.5 4.5 3.9
##  [631] 5.1 4.1 5.4 5.5 4.7 4.5 4.0 3.7 6.8 4.1 3.7 3.1 5.0 4.8 5.8 5.3 4.4 3.2
##  [649] 5.4 3.1 4.3 4.1 4.6 4.0 4.6 3.4 3.6 4.6 3.1 5.5 4.9 5.7 5.5 3.4 4.0 5.2
##  [667] 4.3 2.6 4.4 5.1 4.5 4.2 5.0 3.0 2.9 5.0 4.8 4.5 3.4 5.1 3.9 4.6 3.2 5.0
##  [685] 4.8 5.2 4.8 4.5 5.3 4.4 4.3 5.0 4.5 4.6 4.2 3.7 4.8 4.5 3.0 3.1 5.7 5.3
##  [703] 4.9 6.5 4.0 5.8 5.4 4.4 3.9 5.1 3.0 4.3 5.0 4.4 4.6 6.4 4.4 6.0 5.0 5.1
##  [721] 3.9 4.4 4.3 5.4 4.6 4.0 4.5 3.9 3.3 3.9 5.4 4.5 6.4 5.2 6.3 5.1 4.4 5.3
##  [739] 6.0 4.2 4.3 5.4 3.9 3.2 3.9 4.8 5.4 3.3 2.8 3.3 4.9 4.1 3.1 4.5 4.4 6.3
##  [757] 4.0 5.0 4.2 5.1 4.0 5.1 2.8 5.0 5.6 4.1 2.8 2.6 5.1 3.5 3.1 4.6 4.8 3.2
##  [775] 4.8 6.4 4.5 4.0 5.0 5.2 3.9 3.5 3.5 4.4 3.8 2.9 2.9 5.5 4.9 6.3 5.7 4.8
##  [793] 4.0 5.4 4.9 5.6 3.9 4.0 6.0 3.5 3.2 3.3 4.4 5.0 3.5 4.3 4.8 4.0 6.3 5.3
##  [811] 3.4 4.5 4.2 4.7 5.2 6.5 3.2 3.8 5.2 5.1 4.5 2.0 3.8 4.4 6.4 4.2 4.8 3.7
##  [829] 5.0 4.2 3.8 2.9 5.2 4.8 4.4 5.4 6.4 4.2 5.2 4.4 3.9 4.3 3.5 4.5 4.2 5.5
##  [847] 3.2 4.1 4.0 3.8 3.6 4.5 5.2 5.8 4.4 3.2 4.7 3.5 4.7 4.2 4.2 4.2 2.1 3.6
##  [865] 5.8 4.7 3.7 5.1 4.7 6.3 5.6 4.3 5.2 5.2 3.7 5.9 4.5 4.1 4.3 4.7 4.4 4.1
##  [883] 4.2 4.8 4.8 6.5 3.8 3.2 5.1 5.4 4.2 6.6 5.0 4.6 4.4 6.6 4.7 7.0 4.6 4.6
##  [901] 4.7 4.6 5.1 3.1 3.8 3.7 2.6 5.4 3.3 4.7 4.2 4.0 3.9 2.6 5.1 4.4 3.1 4.4
##  [919] 3.4 5.3 5.1 5.0 4.5 3.2 5.8 5.6 1.3 4.0 5.5 3.7 5.7 3.9 3.8 4.4 3.4 5.0
##  [937] 3.0 4.5 5.3 5.0 4.8 4.5 3.0 2.7 4.3 5.8 4.6 6.3 4.5 4.7 3.4 4.6 4.3 3.0
##  [955] 3.8 5.1 3.8 4.1 3.4 4.5 4.7 3.6 4.6 6.1 5.6 3.5 3.6 4.4 2.7 3.2 6.2 2.4
##  [973] 4.4 6.0 4.7 3.9 3.9 2.9 5.2 4.0 5.6 5.4 5.2 5.6 5.5 3.9 4.7 4.0 3.5 5.2
##  [991] 5.3 4.7 5.1 3.4 4.4 3.6 4.4 5.1 3.9 5.7
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
## 2.7000 6.3025
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
##    [1] 5.1 4.3 3.8 4.4 4.1 3.7 3.2 5.1 6.0 5.5 2.9 3.3 4.4 4.0 4.6 4.2 4.2 4.0
##   [19] 5.5 3.6 5.6 5.4 4.6 5.6 4.0 4.4 5.8 3.1 5.7 4.6 6.0 5.1 3.6 5.3 4.4 3.9
##   [37] 3.7 4.6 5.0 2.5 4.9 5.4 5.1 3.5 5.0 4.9 4.0 3.3 4.4 5.5 5.4 4.1 3.4 4.0
##   [55] 5.0 4.5 3.5 4.8 2.9 4.8 4.1 4.1 5.3 4.5 4.9 4.8 4.4 4.7 4.3 5.2 4.9 4.5
##   [73] 4.5 5.4 5.0 4.8 5.0 3.4 2.8 3.2 3.5 4.6 3.3 3.3 2.9 5.7 2.4 4.9 5.2 4.1
##   [91] 5.6 4.7 2.4 4.7 4.8 3.8 4.1 3.7 4.2 4.5 4.3 4.6 4.9 3.8 5.4 5.1 5.5 3.5
##  [109] 6.0 5.8 5.7 3.1 4.0 4.3 4.0 5.6 3.7 4.3 2.8 4.7 4.9 4.3 3.3 3.1 5.8 4.7
##  [127] 3.2 5.7 5.6 5.1 3.6 7.1 3.7 4.1 4.5 6.2 6.8 5.3 3.4 4.5 3.3 3.8 3.5 4.3
##  [145] 6.2 4.6 4.1 4.8 2.5 4.4 2.7 5.4 3.4 4.4 5.2 5.4 5.2 4.9 4.1 6.1 4.7 5.5
##  [163] 3.6 3.9 5.3 4.2 2.9 3.0 4.7 7.0 5.3 5.1 3.9 4.1 4.7 4.3 4.0 4.0 4.8 4.4
##  [181] 4.5 5.6 4.6 5.8 3.3 4.3 4.3 5.1 5.2 4.7 5.1 6.6 4.5 4.1 5.9 3.4 4.8 4.8
##  [199] 6.0 4.6 3.8 5.6 3.0 4.3 4.4 4.3 3.9 4.3 4.3 3.5 3.4 5.6 3.3 4.6 3.4 4.3
##  [217] 4.7 4.9 5.1 5.1 4.0 3.5 3.6 4.5 6.5 5.7 4.3 5.0 3.5 4.8 4.0 4.9 4.6 4.4
##  [235] 4.3 5.0 4.8 2.9 4.2 6.7 2.1 4.9 5.0 4.7 4.3 4.9 4.3 5.4 5.4 3.9 4.4 4.3
##  [253] 3.7 2.9 5.2 3.1 3.9 3.2 3.8 5.2 4.8 5.6 5.3 3.4 5.1 4.4 4.3 4.1 5.4 2.7
##  [271] 6.1 5.5 5.5 3.5 5.4 4.7 3.2 5.1 3.7 1.6 3.7 4.8 5.7 3.3 4.8 5.3 4.8 5.5
##  [289] 2.5 2.1 6.0 3.6 4.7 5.0 5.3 4.6 4.9 4.7 4.1 4.0 4.1 2.8 5.7 4.1 4.4 4.3
##  [307] 6.7 5.8 5.9 3.1 4.4 3.9 3.1 3.8 5.4 5.7 6.6 4.6 4.2 5.8 4.3 2.3 5.1 3.8
##  [325] 3.1 4.9 5.8 5.1 4.5 7.2 4.4 3.3 4.0 3.6 4.5 4.9 5.2 4.0 6.2 5.5 5.3 4.7
##  [343] 3.1 4.4 5.4 4.3 4.3 4.4 4.5 4.1 4.1 3.7 4.1 4.8 4.4 4.4 4.9 5.2 5.3 3.1
##  [361] 5.0 4.3 3.2 2.6 5.0 4.8 5.6 5.5 4.2 4.5 5.6 5.7 3.4 4.2 5.4 3.7 4.1 7.1
##  [379] 4.7 4.2 3.9 3.4 4.4 4.7 5.0 5.6 5.1 4.3 3.3 4.9 4.2 5.3 4.1 4.5 5.2 4.3
##  [397] 4.4 4.1 4.2 4.6 3.4 4.5 5.7 5.7 4.8 4.9 4.7 3.7 5.1 4.9 3.8 6.7 4.1 3.9
##  [415] 3.1 3.5 4.1 4.0 5.6 5.0 3.2 5.8 3.8 5.0 4.3 5.7 3.3 5.6 3.9 3.5 5.8 5.9
##  [433] 3.5 5.4 4.3 2.7 4.2 3.3 5.0 5.6 4.7 3.9 3.3 5.2 3.9 4.2 5.3 6.1 4.2 4.3
##  [451] 5.7 3.9 3.8 3.4 4.0 4.9 5.1 4.5 4.4 2.8 5.0 4.8 5.6 3.5 4.9 4.4 4.7 4.3
##  [469] 5.6 4.3 5.4 4.0 3.5 4.6 4.8 5.4 4.6 3.8 3.6 3.6 4.3 3.3 7.1 3.9 4.2 4.2
##  [487] 6.4 3.5 4.7 3.4 5.9 4.3 5.1 4.8 4.7 4.3 3.6 4.7 5.1 5.5 2.7 4.1 4.8 4.8
##  [505] 3.5 3.0 4.9 3.5 5.2 4.4 5.9 3.5 4.1 4.1 3.9 4.7 5.2 3.4 3.4 6.3 4.0 3.3
##  [523] 3.3 5.6 4.7 5.4 6.6 5.1 6.5 3.3 2.5 4.6 3.9 4.8 5.1 4.3 3.7 5.1 5.6 6.9
##  [541] 5.0 4.0 4.8 3.5 4.8 2.8 3.7 6.8 4.7 4.7 4.3 3.8 4.3 4.0 4.5 5.5 4.8 5.3
##  [559] 5.4 4.4 5.5 4.1 4.1 3.9 5.3 4.1 3.8 5.1 5.6 5.0 4.4 4.0 5.6 3.7 4.0 5.0
##  [577] 3.9 4.4 4.4 5.9 3.1 4.3 4.6 5.9 4.2 5.1 5.6 5.5 3.8 3.3 4.7 3.7 3.6 4.8
##  [595] 4.3 3.5 5.3 3.2 3.9 4.9 5.1 6.1 2.8 2.9 4.0 5.1 4.8 5.0 4.9 3.0 5.0 3.5
##  [613] 4.0 3.7 3.0 3.7 3.9 4.9 4.0 5.0 5.3 4.3 5.1 3.4 4.5 3.9 4.5 4.5 3.3 5.0
##  [631] 3.9 4.4 5.0 5.6 6.2 4.2 3.8 3.1 5.6 4.5 5.1 3.8 3.9 4.3 5.6 3.0 5.4 4.4
##  [649] 4.7 5.0 5.3 3.5 6.9 3.4 4.5 5.9 5.6 6.6 5.4 5.4 4.8 3.8 3.4 4.6 3.9 5.7
##  [667] 3.7 3.8 4.6 3.4 4.0 4.6 5.3 4.7 4.9 5.2 3.4 3.6 4.9 3.7 3.8 4.6 4.2 4.7
##  [685] 6.4 4.5 4.6 5.1 4.7 4.3 4.1 4.4 4.1 3.6 2.6 4.4 5.6 4.8 4.9 2.5 4.3 4.9
##  [703] 5.8 4.3 4.0 6.0 5.2 4.9 6.0 5.4 5.5 6.1 5.1 4.9 5.3 4.3 3.3 4.3 2.5 5.4
##  [721] 4.2 3.8 3.5 5.1 4.0 4.2 5.2 4.1 6.2 3.8 3.6 5.2 6.1 3.2 3.8 3.7 4.1 5.2
##  [739] 2.9 3.0 3.6 3.5 4.3 5.5 5.3 5.7 3.9 5.7 5.0 5.5 4.5 4.7 5.4 5.2 3.4 3.6
##  [757] 3.1 2.7 5.0 6.2 2.6 6.6 4.6 5.9 4.7 2.9 4.8 3.5 6.0 3.7 5.2 3.9 4.5 5.5
##  [775] 6.1 4.8 3.9 3.6 5.6 5.1 5.5 3.4 4.8 5.8 4.2 4.3 4.6 4.1 4.6 3.9 4.0 5.2
##  [793] 3.9 3.0 4.4 5.1 4.2 4.5 2.8 4.0 3.7 3.8 6.0 4.7 4.2 5.3 3.8 6.4 5.5 4.7
##  [811] 3.7 3.8 4.7 3.2 3.7 5.7 5.4 4.7 3.7 5.7 6.3 4.5 4.4 4.9 4.9 5.7 4.8 5.7
##  [829] 3.8 3.7 4.0 2.9 4.4 3.1 3.2 3.2 4.2 3.3 4.7 5.5 5.2 3.4 3.6 5.5 5.5 3.6
##  [847] 3.5 4.9 5.1 4.3 4.2 6.0 3.5 3.7 4.0 2.6 2.8 6.0 4.4 5.0 3.3 3.7 2.2 2.9
##  [865] 4.7 4.5 3.6 4.3 5.2 5.2 4.4 5.2 4.3 5.1 5.1 3.5 4.2 6.3 4.8 3.6 4.3 5.0
##  [883] 5.4 2.9 3.8 3.3 3.7 4.5 4.1 4.8 4.2 3.9 4.2 3.3 4.7 3.0 4.6 5.0 4.2 3.9
##  [901] 5.9 5.0 4.7 3.9 5.5 3.3 4.0 6.1 3.5 5.2 4.7 4.6 5.0 3.8 5.9 4.8 3.2 4.9
##  [919] 4.0 3.6 4.2 4.3 5.8 3.8 5.4 4.2 3.8 5.1 4.2 5.4 4.0 5.4 5.4 4.0 3.4 5.6
##  [937] 5.5 4.0 4.2 4.0 3.4 4.1 4.4 4.2 5.2 3.1 4.7 4.8 5.7 5.1 3.6 4.4 4.5 4.4
##  [955] 4.6 4.5 4.8 3.0 5.1 3.3 5.3 3.4 5.7 3.4 2.0 4.1 2.7 3.5 3.4 4.3 3.9 4.3
##  [973] 4.9 4.3 5.2 4.8 4.3 3.6 5.7 2.4 5.4 4.2 4.1 4.8 3.1 3.9 4.8 5.1 5.3 4.7
##  [991] 4.8 4.8 4.0 4.9 5.9 4.1 5.5 3.9 6.5 5.0
## 
## $func.thetastar
## [1] -0.0189
## 
## $jack.boot.val
##  [1]  0.51415929  0.37207447  0.21899441  0.13241758  0.10685358 -0.04763314
##  [7] -0.26318841 -0.30849315 -0.44171271 -0.51944444
## 
## $jack.boot.se
## [1] 0.9975978
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
##    [1] 4.4 3.8 5.3 4.0 5.8 4.4 3.6 4.8 4.3 4.0 3.2 3.8 4.2 4.0 4.5 4.9 4.6 3.3
##   [19] 3.3 4.1 4.9 3.7 5.3 5.4 2.6 5.4 4.5 5.5 5.1 3.6 6.0 4.9 5.6 4.1 4.8 4.4
##   [37] 4.6 3.9 4.8 5.4 4.3 5.2 5.0 4.4 2.7 6.0 5.0 4.5 4.6 3.4 3.6 2.6 5.6 2.8
##   [55] 4.5 5.8 2.9 5.3 4.9 4.2 5.9 4.3 3.4 4.7 2.6 4.5 4.4 4.0 4.1 4.9 5.0 5.1
##   [73] 3.5 4.9 3.1 6.0 4.6 4.3 7.3 3.5 4.0 3.8 5.0 2.5 4.1 4.2 5.3 5.2 2.5 3.7
##   [91] 5.8 3.0 3.8 3.6 3.4 5.9 4.3 3.1 3.2 4.2 4.2 4.0 3.2 3.0 3.7 5.5 5.3 5.0
##  [109] 6.2 4.1 5.3 6.5 4.5 4.5 5.6 5.5 5.2 4.4 5.2 5.8 3.9 6.5 4.7 6.4 4.1 4.7
##  [127] 4.5 4.8 4.8 5.6 3.8 3.2 4.9 4.5 4.9 3.2 4.3 5.0 3.5 4.7 3.7 5.9 4.7 4.3
##  [145] 2.9 5.3 2.8 3.7 4.0 4.5 5.2 4.9 5.4 4.3 4.7 3.3 5.4 7.3 4.7 3.7 4.1 5.8
##  [163] 3.1 4.1 3.5 5.3 4.3 5.0 5.0 4.8 3.0 5.1 7.4 3.5 3.7 3.0 4.4 4.1 4.0 4.9
##  [181] 6.3 3.4 5.1 4.4 5.2 3.9 3.6 5.2 5.0 4.5 4.1 4.1 5.4 4.7 4.3 4.4 4.3 4.5
##  [199] 4.1 4.8 4.8 4.4 5.9 5.2 4.2 6.5 4.7 2.8 5.9 3.3 5.1 6.1 4.1 5.6 3.5 5.4
##  [217] 5.3 4.5 5.8 5.1 3.5 4.8 2.8 5.1 3.8 4.4 4.8 4.7 4.4 3.9 3.9 5.6 3.8 4.2
##  [235] 4.8 6.1 3.4 4.6 4.2 5.8 3.1 4.2 4.4 5.0 5.2 5.3 3.7 5.6 4.6 4.8 4.3 3.8
##  [253] 5.3 4.0 6.2 4.1 4.0 6.4 5.6 4.8 4.6 3.7 5.0 4.3 3.1 5.3 3.9 3.4 5.0 5.2
##  [271] 3.3 3.9 3.8 4.4 6.1 4.9 4.9 4.5 4.4 3.8 4.0 4.8 4.6 4.3 6.1 3.1 5.5 4.6
##  [289] 4.3 4.4 4.7 2.9 4.6 4.3 2.0 4.8 5.6 4.7 6.1 4.8 5.3 4.9 4.3 4.3 4.4 5.0
##  [307] 5.8 4.1 5.6 4.9 5.4 6.5 4.7 3.3 4.1 3.7 5.8 3.8 6.6 3.9 4.3 4.2 4.2 4.5
##  [325] 4.6 4.8 6.6 5.5 4.9 3.5 5.0 5.1 4.9 4.2 4.7 4.1 3.9 4.8 4.2 6.0 3.4 5.3
##  [343] 6.8 4.9 4.3 4.5 5.7 4.5 5.1 4.9 4.1 5.7 4.7 3.7 4.0 4.8 5.7 3.5 7.8 4.2
##  [361] 3.6 3.4 3.4 5.3 2.5 1.6 5.2 5.2 3.4 4.3 3.6 6.1 5.6 3.8 5.1 3.7 2.6 4.8
##  [379] 4.4 5.8 4.6 3.7 5.0 5.4 5.2 4.4 2.3 4.8 4.8 7.0 4.2 2.7 4.6 4.8 5.0 3.4
##  [397] 5.0 5.3 5.4 5.4 4.0 4.9 5.7 5.3 4.8 4.2 5.1 3.6 4.5 2.9 3.7 3.8 4.1 5.3
##  [415] 4.5 4.3 3.8 5.8 3.6 5.4 3.6 3.9 4.2 5.5 4.2 4.5 3.6 4.4 4.4 4.3 3.9 5.6
##  [433] 4.3 3.5 4.3 5.5 4.4 3.8 4.3 4.3 3.9 3.0 4.2 4.3 3.9 4.7 4.5 5.9 3.8 3.4
##  [451] 4.5 3.7 3.9 4.1 3.2 3.7 4.8 4.7 4.8 5.6 4.6 4.6 4.3 5.4 5.9 5.3 3.1 6.7
##  [469] 2.7 6.1 5.7 5.9 4.6 5.0 2.4 3.3 3.2 3.6 4.4 3.5 5.2 4.2 5.2 5.4 5.4 3.1
##  [487] 4.9 3.9 4.6 5.2 4.7 4.1 4.0 5.3 4.3 3.5 3.9 3.5 5.0 5.5 4.5 3.1 6.3 3.9
##  [505] 5.7 5.4 5.2 3.6 4.6 4.8 4.9 5.6 5.5 3.4 3.4 4.9 2.0 4.1 5.4 5.9 3.3 4.7
##  [523] 2.0 4.7 3.5 3.5 4.4 4.2 3.9 4.0 2.9 4.3 3.5 5.3 5.5 4.1 5.3 6.1 5.1 4.9
##  [541] 3.4 3.7 4.1 4.4 5.5 4.4 4.8 5.4 5.5 5.9 4.5 3.6 5.3 5.2 4.0 3.9 6.3 1.8
##  [559] 4.7 5.2 4.7 4.5 3.7 3.6 3.0 5.3 4.9 5.9 5.5 4.5 4.3 3.1 4.9 6.2 4.7 3.4
##  [577] 3.6 4.0 3.1 6.1 4.4 5.5 3.6 3.6 3.5 5.4 5.0 4.4 2.3 4.1 5.1 4.1 4.6 6.9
##  [595] 4.5 2.8 3.7 4.6 4.2 4.9 5.3 4.8 5.9 4.5 4.0 3.8 4.4 4.0 4.1 3.7 3.2 5.1
##  [613] 3.6 5.2 4.1 4.1 5.3 3.8 5.0 5.0 5.6 5.4 3.3 2.4 4.6 4.8 2.8 3.9 5.6 3.4
##  [631] 4.9 3.0 4.2 3.9 4.9 3.8 4.3 5.2 6.2 2.7 5.5 5.2 5.1 4.4 5.0 4.7 3.5 6.1
##  [649] 3.5 4.8 5.0 5.2 4.0 6.2 6.0 4.5 4.5 3.6 4.6 5.0 4.7 6.4 4.8 4.9 5.6 4.6
##  [667] 3.6 5.2 5.5 4.9 5.0 3.8 4.7 6.0 4.0 5.0 2.8 6.5 3.2 2.8 4.2 4.1 4.0 3.5
##  [685] 2.7 5.5 5.1 5.4 3.5 3.9 4.6 5.5 5.1 5.6 4.6 4.1 2.7 3.9 3.8 5.9 4.2 2.7
##  [703] 4.7 4.1 3.3 6.2 4.0 5.7 4.3 3.0 4.8 3.9 5.4 5.4 3.9 4.8 4.9 4.7 5.4 4.6
##  [721] 4.1 3.8 4.3 5.2 4.4 4.5 2.9 4.8 3.9 5.2 4.3 5.2 4.9 4.4 4.2 3.9 2.7 1.9
##  [739] 4.0 4.2 4.7 4.1 3.7 3.4 5.7 3.9 2.7 4.6 5.0 3.9 4.9 3.9 5.8 5.4 4.6 3.3
##  [757] 5.9 4.3 4.3 4.8 3.3 5.2 5.2 2.4 5.1 4.1 5.0 4.0 4.0 4.8 5.7 4.4 4.7 5.6
##  [775] 4.2 3.2 5.0 4.9 4.3 4.3 4.9 5.0 5.0 5.4 5.4 4.4 4.7 5.3 5.2 4.4 4.0 1.8
##  [793] 4.2 4.5 4.7 4.1 5.6 2.6 3.9 3.7 4.1 3.9 5.1 2.9 4.4 4.8 3.7 4.5 5.8 4.6
##  [811] 3.6 3.4 6.1 4.9 5.0 2.8 6.0 4.1 4.7 4.4 7.4 5.8 3.9 4.0 4.2 4.3 5.6 4.5
##  [829] 5.1 4.0 6.0 5.6 3.3 3.6 6.0 3.4 4.5 3.7 3.6 3.9 4.4 3.1 3.7 4.6 3.4 2.4
##  [847] 3.9 6.3 4.1 3.8 2.2 4.6 4.0 4.4 5.6 4.2 3.4 5.0 6.7 3.5 7.1 3.5 3.1 4.5
##  [865] 5.3 6.5 3.9 2.1 4.5 4.0 4.8 6.1 6.0 3.2 4.8 3.6 3.2 4.5 4.8 3.6 5.4 3.4
##  [883] 4.5 4.7 5.4 4.8 5.0 4.1 4.6 5.0 3.9 5.2 3.4 3.0 4.2 3.9 3.5 3.2 6.4 4.9
##  [901] 4.5 6.1 4.8 4.1 4.3 4.2 3.5 5.6 3.5 5.3 3.6 5.3 5.0 4.1 3.7 3.9 3.5 4.5
##  [919] 4.1 3.0 3.9 4.5 3.9 5.0 4.4 6.0 3.7 3.6 3.5 5.7 3.9 5.6 5.0 5.0 4.1 4.4
##  [937] 5.4 5.0 6.1 4.9 4.6 3.6 6.8 6.2 4.0 3.7 4.4 5.2 4.4 5.5 4.4 4.1 5.9 5.3
##  [955] 3.6 5.4 4.1 5.4 4.2 4.7 3.7 5.1 5.7 4.7 4.7 5.2 5.1 4.6 4.9 4.5 5.9 4.9
##  [973] 4.7 2.9 5.8 5.3 5.0 3.1 4.3 4.1 4.3 4.8 6.0 5.3 4.8 5.5 4.2 3.5 5.5 4.4
##  [991] 3.5 3.2 5.3 3.1 3.8 4.2 4.0 3.1 3.9 6.3
## 
## $func.thetastar
## 72% 
##   5 
## 
## $jack.boot.val
##  [1] 5.500 5.500 5.400 5.300 5.100 5.200 4.800 4.800 4.656 4.400
## 
## $jack.boot.se
## [1] 1.088013
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
## [1] 0.2806456
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
##   2.617820   4.301649 
##  (1.104288) (1.999895)
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
## [1]  0.9285461  1.8130585  0.6209366 -1.1370996  0.5878266  0.9579254
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
##    [1]  0.2039859745 -0.1347469494  0.9593224599 -0.2870415809  0.9219015517
##    [6] -1.0832526426  0.0179367284  1.0151822362 -0.6056172157  0.1255310867
##   [11]  0.3188389589  0.7959531676  1.0368156102  0.0340455789  0.1654093299
##   [16]  0.3332924409  0.9260861345  0.0622523933  1.1049266323 -0.1296614200
##   [21]  0.0508499192  1.4112391020  0.0736134474  0.6066436683  0.7163497694
##   [26]  0.8575764547  0.0566984377 -0.2747385601  0.3467033446  0.0813557065
##   [31]  1.7708346879  0.6169771877 -0.4067451397 -0.2139072725  0.8628357293
##   [36]  0.3537585551 -0.0275122859  0.3456046591  0.5865569084  0.5812619188
##   [41]  1.1765773644  0.4347711820  0.7420776431  0.7922198630  0.5374858833
##   [46]  0.3567635309  0.9471968111  0.3160488436  0.2644707882  0.4371707021
##   [51]  0.7572259263  0.8435653227  0.0286624909 -0.2619766041  0.5450471846
##   [56]  0.1412290918  0.7127361530  0.4414644671  1.2026867628 -0.8996698801
##   [61]  0.2903122415  0.2344301010  0.4633320400  1.3295492628  0.8344482367
##   [66]  0.4119646684  0.2588165874  0.2899478822  0.0065781799  0.3486443189
##   [71]  2.2186633565  0.6047912911  0.4003764903  0.4361536793  0.0780648679
##   [76] -0.4745276322 -0.3987022537  0.5504782725  0.0758577315  0.7585884811
##   [81]  0.2855252417  0.2047222110 -0.0294190625  0.5720161552  0.0019898750
##   [86] -0.4080393857  0.7094023313  1.6755806951 -0.8646366188  0.1496709180
##   [91]  1.2881091375  0.7027377930  0.3109729299  0.1002577329  0.1738117514
##   [96] -0.0026656008  0.2055759793  0.1765601813  0.3737436174  0.6939030719
##  [101] -0.2555435780  0.5933434869  1.0269601105  0.0942596843 -0.1466905108
##  [106]  1.3934855734  0.2968479646 -0.8355824273  0.0019763488  0.1458270752
##  [111] -0.2609130034  0.4901165943  0.1928164329 -0.0844849616  0.6119454418
##  [116]  1.6275898678 -0.0848397414  0.8978255664  0.7995545890  0.0710124199
##  [121]  0.3314439566  0.1174214917  0.3830522200  0.1208925724 -0.8806522929
##  [126]  0.6571876996 -0.2040591630  0.5250344284  0.4524745520 -0.6357925257
##  [131]  1.0236496384  0.3939750848  1.1211539588 -0.2665669948 -0.2350275721
##  [136] -0.3243309139 -0.5211782310  0.0286364754  0.0197777433 -0.0701467434
##  [141] -0.0160478489  0.4678963728  0.1308888201  0.1303515750  0.5724469374
##  [146] -0.6765635338 -0.2575924061  0.3544812063  0.2496720571 -0.1938360767
##  [151]  0.9027099960  0.6152109188  0.0357209599  0.1588045837  0.1097579949
##  [156]  0.1462480267  0.5733105488 -0.0355317516  0.3619678784  1.3093155413
##  [161]  0.5056375044 -0.0043131335  1.0384841008  0.3404746237  1.1345121669
##  [166]  0.1263177725  1.4679847403  0.2667314807 -0.0225976316 -0.1616350495
##  [171] -0.0810889940  0.3140308695  0.1484949490  0.8368698823  0.4397626627
##  [176] -0.3581426483  0.7147109434  0.3413978456 -0.2967730009 -0.1180448999
##  [181]  0.1048216438  0.3937407744  0.2709739556  0.0359877276  0.1761377395
##  [186]  0.4211628900 -0.4471308911  0.7583874295  0.2199817228  0.9617616183
##  [191]  0.0466898495 -0.1951739291  1.0194526946 -0.6992510922  0.4042454388
##  [196]  0.0592652662 -0.3980659249  0.2974855140 -0.4322087228  0.5530821872
##  [201]  0.2753183802 -0.6219780078  0.1701943657  1.0227399064  1.2014332722
##  [206]  1.0096315810  0.6235658851  0.1245724049  0.4530886438  0.2663073294
##  [211]  0.7950773813  0.5304164210  0.3430581577 -0.0568115909  0.2940437626
##  [216] -0.5571031687 -0.0907547225  0.2693699677  0.0339028532  0.0766709702
##  [221]  0.0816120877  0.0942596843 -0.6315949850  0.1150155505  0.4948968326
##  [226]  1.6188809175 -0.1205002606  0.2359005139  0.2780761683  0.6332716534
##  [231] -0.0907742966 -0.8780383871  0.1054429490 -0.0341617386  0.2080974638
##  [236]  0.2996513556  0.2692047869 -0.3958447405  0.2231878296 -0.2259095171
##  [241] -0.5822223221 -0.6096076378  0.3600817608  0.3468445693  0.3199353491
##  [246] -0.6256069026  1.0810666704  0.2575042740 -1.0066077943  0.2211638316
##  [251] -0.3614683163  0.0820908099  0.0408016671  0.1367466298 -0.1071804457
##  [256]  0.4093394237  0.2056293280  1.3078774624  1.0186565551 -0.2296994004
##  [261]  0.5060602974  1.1828874816  0.3803526139  0.5133772766  0.0793645468
##  [266]  1.0419365396 -0.1890066337  0.6423710464  1.1941426961 -0.1103780335
##  [271]  0.3942469749  0.5679266270  0.0297806903  1.3271940118  0.8880525507
##  [276]  0.0126830117  0.6109056421  0.8743126002  0.0557622853  0.0419447236
##  [281]  0.1297127929  0.3717537644  0.2993403665 -0.6499627794  0.8449715269
##  [286]  0.4042997303  0.4620183330  0.2771991133  0.7591079821  0.4319641851
##  [291]  0.2273297735 -0.6203057039  0.4912452875  0.6335521653  0.3604819446
##  [296]  0.6269169075 -0.1611393344  0.6070858968 -0.1658401349 -0.5110382787
##  [301]  0.2424075914  0.4674922461  1.2546506662  0.2510545011 -0.0870354097
##  [306]  0.2973452281 -0.5401189967  0.2167307292  0.9680070012 -0.4554699601
##  [311]  0.0797541680  1.4329609116 -0.1050404358 -0.3508756999  0.5158411493
##  [316]  0.5815500064 -0.3180006712  1.2153143467 -0.0624665956  0.9188645053
##  [321]  0.3010855656 -0.0724581644  0.4802284203  0.0827557771 -0.2211658343
##  [326] -0.4404902322  0.5116124728  0.1911293650  0.0702005094  0.0006273787
##  [331]  0.6086957074  0.2523017657 -0.4637740559  1.5554518884 -0.4955609766
##  [336]  0.4699042559  0.3298022351  0.1165354524 -0.3032219532 -0.4521461124
##  [341] -0.6713072347 -0.2909430605 -0.0956582744 -1.0676585320  0.1851141726
##  [346] -0.5801362802 -0.1103780335 -0.1716405382  1.0008832578  1.7797601336
##  [351] -0.1179686190  0.6147143981 -0.3762690629  0.3932594264  0.0073704627
##  [356]  0.0112008951 -0.8133226505  0.4218868058  0.6939030719  0.3231634153
##  [361]  0.9908798305 -0.4857124463  1.1504472556  0.0590042482  0.2086114541
##  [366]  0.6296847066 -0.2242830283  0.6962196529 -0.1286918356  0.1762799125
##  [371]  0.8787608783  0.8045372625  0.0840415738  0.1079650998 -0.4641460875
##  [376] -0.3123658838 -0.0936246197 -0.1569209473 -0.0909340483 -0.6938642111
##  [381] -0.0308744911  0.0487385386 -0.0700294508  0.3301399109  0.6849589724
##  [386]  0.5223165438 -0.4947110846  0.3793732790  0.3783912843  0.9391390734
##  [391] -0.7431571723 -0.0986422131 -0.2234562075  0.5323303497 -0.0906269329
##  [396]  0.8149005870 -0.0351688154  0.8501220622  0.7765023738 -0.0659383372
##  [401]  0.3095666599  1.2012805185  1.2318279432 -0.0617841831  0.4831722685
##  [406] -0.3887556828  0.8753738213  0.2256652086  0.7131537093  0.4022854842
##  [411]  0.4666623953  0.6841494214  0.3631961660 -0.0923114708  0.3538388035
##  [416]  0.6193400703  0.2357936308  0.1609080773  0.0526333705  0.0937723839
##  [421]  0.2277594427  0.5675034680  0.0704999318  0.6237099390  0.7136048802
##  [426] -0.1005547379  0.3709217729  0.1404544504  1.1068368798 -0.0836570435
##  [431] -0.1344147857  0.7509393317  0.7529744746 -0.4117614202  1.7062140136
##  [436]  1.1069419368  0.4015491383  0.3109136519  1.1431273876  0.9359574504
##  [441]  0.0449406598 -0.5687525483  0.9090585245  0.7024251783  0.7851205246
##  [446]  1.1493252162  0.5864117191  0.5494180595 -0.8844271055  0.7166403709
##  [451]  0.2080974638  0.5116124728  1.0006934267  0.6956092591  1.8044667048
##  [456]  0.4155727488  0.4567643824  0.4876197098  0.7055297365  0.1286087442
##  [461]  0.5683835502  0.4126534902 -0.9280531600  0.6840246802  0.0803188026
##  [466]  0.8120126364  0.1765487223  0.0272213539  0.0852719199  0.4672511025
##  [471] -0.0049379385 -0.3074769602 -0.0844175069  0.2465679717  0.4150342558
##  [476]  0.3164433635 -0.3089671545  0.3495612817  0.5159686787 -0.1952571152
##  [481]  0.6699409607  0.5623706881  0.7459199928  0.2391861086  0.4144491674
##  [486]  0.7925202969 -1.2569848722  0.0295482912  1.0152560848 -0.0422786964
##  [491]  0.7648877770  1.0839411609 -0.2867906090  0.2239801441  1.0516823178
##  [496] -0.4276606322 -0.0532567510 -0.1712545309  0.7413497224  0.2742559032
##  [501] -0.0636243839  0.6367629208 -0.0350834694  0.7872465290  1.6743770239
##  [506]  0.4783685044 -0.1008479158  0.1717916244  0.5091167278  0.0577897947
##  [511]  0.4064768441  0.5717346076  0.2670864721 -0.5025599077  0.6118312138
##  [516]  0.8204275824  0.0678785190  0.8176032588  0.8149005870  0.0020879046
##  [521]  0.1835799695  0.0677932568  1.1561334419 -0.1214558988  1.1933905291
##  [526] -1.1371621764  0.1872697882  0.1895402899  0.1547222746  0.4017339188
##  [531]  1.1182370199 -0.0392075587  0.2770599556  0.6489272889  0.4285534078
##  [536]  0.2060687426  0.5302612443  0.0273723675  0.3020388576  0.0215614843
##  [541]  0.1239143491  0.9485309675  0.7881812475  0.6193893955  0.5352347746
##  [546]  1.5706289369 -1.0613339600 -0.4911470537  0.6030425309  0.7932450592
##  [551]  0.3725956274  0.5046909656  0.5714803674  0.0428602599  0.2420467125
##  [556]  1.4111521732  0.4568200963  0.0136590860 -0.3275654854  1.9098470049
##  [561]  0.2070473975  0.2047357802  0.3401776833  0.3010381433  0.2902351596
##  [566]  0.0770497329  0.2293170324 -0.4278632654 -0.3668969553  0.1240522152
##  [571]  0.1318356288 -0.0766809175  0.4462524666  1.0063343980  0.0275651995
##  [576] -0.0143614231  0.1585854759 -0.8734528681  0.6627363135 -0.1137509919
##  [581]  0.7158645099  0.2909603018  0.8590431553 -0.4646787736  0.1832975947
##  [586]  0.0136282565 -0.0869416527  0.4358315076  0.3091769894  0.9691105896
##  [591]  0.9356064835  1.0912473255  0.5481692631  0.3963163146  0.1477914050
##  [596] -0.0155054184  1.0709383336 -0.3773419864  0.7165430532  0.0354102281
##  [601] -0.1334321032  0.1878542347  0.9790009231 -0.0661387493  0.3023859363
##  [606]  0.8595275377 -0.0505402165  0.2996707180  0.1794980331  0.1417017625
##  [611]  0.4122830615  0.3842403519  0.1328848094  0.9253091841  0.1327980616
##  [616]  0.0105209038  1.1157786510  0.6860869926  0.5410226626  0.8990142822
##  [621]  0.9968175202  0.0815467880  0.4466567490  0.4346784163 -0.3562011534
##  [626]  0.4378869651  0.3509952167  0.1991651435  0.6964459374 -0.0423233283
##  [631] -0.1327583349  0.1649614651  0.6701496864 -0.2532064530  0.3414425743
##  [636]  0.4658460008 -0.4067451397 -0.3607306234  0.1637782350 -0.3381113384
##  [641]  0.1619228603  0.0680544076  0.8065722234  0.2711580946 -0.1866633535
##  [646] -0.0731274365  1.3567744922  0.0166656465 -0.4601880559 -0.0467682143
##  [651]  0.9405263129  0.6727741690  0.0242044447 -0.3870479341  0.3464671617
##  [656]  0.6575610678 -0.3277137793 -0.3157620649  0.6543012164  0.0428984083
##  [661]  1.0444379109  0.3203834864 -1.0141485686  0.1862570309 -0.0807414832
##  [666]  0.6975290198 -0.0294904919  0.4674581410  0.6309930903  0.4341049600
##  [671]  0.6749707789  1.2146845257  0.3331122685  0.0108937983 -0.1994392431
##  [676] -0.1099755062  0.6161450104  1.5888130859 -0.0860097425  0.2801495258
##  [681] -0.4271342881  0.2946616004  0.8441590935  0.2574991486 -0.0724117349
##  [686]  0.2429211329 -0.1755689280 -0.6481118885  0.3086820214  0.3551717064
##  [691]  0.2755413972 -0.6441560452  0.6485768306  1.6886045338 -0.0137586317
##  [696] -0.0592152378  0.9809288231  0.4623094121  0.3078949172  0.8032776514
##  [701]  0.4379619111  0.7479397776  0.2754406362  0.9270082971 -0.3866467750
##  [706]  0.8066323197  0.3181654684 -0.0752936734  0.4673840700  0.8343750694
##  [711]  0.4502980126 -0.1216643747  0.3330102759  0.2101215669  0.6048611349
##  [716] -0.0745450654  0.3990204884  0.2556059865  0.5698367731  0.9033389221
##  [721]  0.4527444905 -0.6083056535  0.1162853613 -1.0629140064 -0.0068891710
##  [726]  0.3785464853 -0.5346046915 -0.2637104154  1.4479347787  0.4051883478
##  [731]  0.4061896594  1.1287072391  0.0459115993  0.7579986014  0.8057332306
##  [736] -0.2833060259  0.5806781508  0.2406547647  0.4216304384  1.0649877196
##  [741]  0.1476479378  0.1518518954  0.2218522192  0.1620424846  0.2950625683
##  [746]  0.8251472608  0.5110674006  0.2309706851  0.8767252879  0.2968896694
##  [751]  0.5105834291  0.9941997635 -0.1124213839  0.5256844764 -0.2391647044
##  [756] -0.5346046915  0.6722782673 -0.2165813832  1.2216018086  0.9122594234
##  [761]  1.1786238386  0.1631186478  0.6853202142  0.0896577520  0.2792611628
##  [766]  0.6842072846 -0.1015704457  0.0966942139  0.2854946429  0.0346584741
##  [771] -0.7031158046  0.3837351810 -0.8876833623  0.5046909656  0.4869256218
##  [776] -0.1163584332 -0.6993877071  1.6755806951  0.4962675518  0.3441078855
##  [781]  0.6015800823  0.3266175726  0.0779494975  0.0991157880  0.6146984174
##  [786]  0.5446418042  0.0366605948  0.1709041046  0.1115755378  0.9552806251
##  [791]  0.4348483445  0.5336420070  0.8269136878  0.7647303183 -0.0271656972
##  [796]  0.5537983934  0.3747696794  0.6131175303  0.1734766111  0.4270050285
##  [801]  0.8797676575  0.1472085899  0.2733475087  0.5288921612  0.3057800800
##  [806]  0.6764173418  0.2664134755  1.1256727910  1.5644582110 -0.2037893078
##  [811]  1.7356024598  0.6352818313  0.5646480802  1.1461874788  0.5866746371
##  [816] -0.1952224511  0.2286043390 -0.7240802846  0.5727074098  0.1930632671
##  [821]  0.3657006208  0.7104021790  0.5747553264 -0.4292229653  0.3413722085
##  [826]  0.3050817848  0.3415275174  0.1421201572  0.0066245713  1.0973235779
##  [831] -0.4631336715  2.0886930131 -0.4545992689  0.7559511404  0.6465384628
##  [836]  0.2846269044 -0.1802504807  0.3148155287  0.1436267311  0.5733105488
##  [841] -0.4260791822  0.7476210697 -0.2725278760 -0.1135652727  0.5127267736
##  [846] -0.6752653814 -0.1182377357  0.6184843500  0.6201951761 -0.6489658925
##  [851] -0.1428293678 -0.4018289156  0.5269810287  0.4945021876 -0.8099033685
##  [856]  1.5609923612  0.2673007726  0.4202355648 -0.5722817228  0.5696508658
##  [861] -0.0844175069 -0.1936628308  1.1369904560  0.1054847816  0.0601108998
##  [866]  0.5401319343  0.1334179109 -0.2431979728  0.7629540681  0.4661027237
##  [871] -0.0673340160  0.2967013212  0.0068107692  1.7843042772  0.0959908076
##  [876]  0.5827114782  0.5137543441  1.3059367068  0.8968488135  0.2122101105
##  [881]  0.6053735124  0.4254966877  0.2070473975  0.1725899541  0.3425917110
##  [886]  0.4537936972 -0.1824922085  0.3688465964  0.5348657141 -0.8751288355
##  [891]  0.0984269759  0.2543071089  1.1437822104  0.1073977608  0.1651030724
##  [896]  1.2505462769 -0.0683872915  0.5951975329 -0.2343266447  0.5527831862
##  [901]  0.6178368154 -0.0310797727  0.3203119308  0.0415257571  0.0112943896
##  [906]  0.6823468815 -0.4040174091 -0.4057079414  0.1097579949  1.0700256359
##  [911]  0.8192703867  0.3032110309  0.1503447183 -0.0666404061  0.5651237287
##  [916]  0.7948876754  0.4748799870  0.7972074317  0.7319643416  0.0552366133
##  [921]  0.4128434768  1.0185247648 -0.4421545567  0.6198101315 -0.0241853773
##  [926] -0.0850202156  0.2394295784  0.5287250610  0.8825753895  0.1671983584
##  [931]  0.5027303545  0.3686591844 -0.1137865029  0.8089870986 -0.1505826331
##  [936]  0.4108481711  1.3852679017 -0.2237446797  0.0788606791  0.6011819116
##  [941]  0.1188869908  0.4042454388  0.2269883036  0.5061050583 -0.1698352891
##  [946] -0.1362576772  0.0392183194  0.3437710102 -0.4311451969 -0.5189454130
##  [951] -0.2058290522 -0.1281522278 -0.3359833008  0.2169028065 -0.1458752475
##  [956]  0.4512196631 -0.3381113384 -0.0871652019  0.8027635041  0.4298617682
##  [961]  1.2881372739  0.0537417677  0.1207411489  0.2955195410  0.1648403885
##  [966]  0.0732535552 -0.0300675434 -0.1006324575  0.1960769769  1.4123369528
##  [971]  0.2906506911  0.9106132597  0.6728575958  0.2700660312  0.4751625273
##  [976] -0.0341060123  0.7932206991 -0.5969572536  0.0257737869  0.1843529187
##  [981]  0.3596012529 -0.3389125960  0.7849640196  0.1281348121  0.5713777182
##  [986]  1.1273629304 -0.2556466267  0.4888557332 -0.0593967470  0.0431677086
##  [991] -0.2882623560 -0.0664503951  0.8306986959  1.2019813441  0.3010815833
##  [996]  0.2070275552  0.4748586231  0.5299459368 -0.1513095339  1.2223249578
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
##   0.6085631   0.3374187 
##  (0.1067012) (0.0754459)
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
## [1]  0.9657651  0.4097591  0.1306197  0.6391737  0.1028379 -0.3358236
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
## [1] 0.0107
```

```r
sd(results2$thetastar)  #the standard deviation of the theta stars is the SE of the statistic (in this case, the mean)
```

```
## [1] 0.9190322
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
## t1*      4.5 0.03873874   0.8948253
```

One of the main advantages to the 'boot' package over the 'bootstrap' package is the nicer formatting of the output.

Going back to our original code, lets see how we could reproduce all of these numbers:


```r
table(sample(x,size=length(x),replace=T))
```

```
## 
## 0 2 3 5 7 8 
## 1 1 1 1 5 1
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
## [1] 0.0351
```

```r
se.boot
```

```
## [1] 0.9340398
```

Why do our numbers not agree exactly with those of the boot package? This is because our estimates of bias and standard error are just estimates, and they carry with them their own uncertainties. That is one of the reasons we might bother doing jackknife-after-bootstrap.

The 'boot' package has a LOT of functionality. If we have time, we will come back to some of these more complex functions later in the semester as we cover topics like regression and glm.

