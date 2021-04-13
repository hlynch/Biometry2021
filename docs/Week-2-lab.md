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
## 1 2 5 6 9 
## 2 1 1 5 1
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
## [1] 0.0244
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
## [1] 2.75516
```

```r
UL.boot
```

```
## [1] 6.29364
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
##    [1] 3.3 3.8 4.4 4.5 4.5 4.6 3.7 4.1 5.7 5.5 5.0 3.9 5.5 4.2 6.7 4.7 4.5 4.0
##   [19] 2.1 6.8 4.6 3.8 3.4 4.7 5.3 6.1 4.0 4.1 2.1 4.3 5.0 6.8 3.9 4.3 6.5 3.6
##   [37] 4.4 5.1 5.5 4.7 4.5 4.8 6.0 3.6 5.9 3.5 3.7 5.0 3.7 5.2 5.5 3.2 4.6 5.0
##   [55] 5.0 3.8 4.8 4.3 4.3 5.7 4.4 5.6 4.2 4.9 3.7 3.2 5.8 3.6 3.8 6.0 4.1 4.8
##   [73] 5.0 1.4 5.2 4.0 4.3 4.0 6.3 4.7 3.3 5.5 4.2 4.8 2.7 5.5 5.2 3.7 4.9 5.8
##   [91] 5.3 5.3 6.1 4.7 4.2 6.0 5.4 3.9 4.1 4.3 2.9 4.6 3.9 4.9 3.4 6.0 4.6 3.5
##  [109] 3.7 5.5 4.5 5.5 4.4 4.6 4.5 3.8 4.7 6.5 5.8 6.0 5.1 5.1 4.0 5.4 6.5 4.8
##  [127] 3.7 3.6 4.6 4.1 4.1 3.6 6.4 5.2 4.3 5.7 3.0 4.4 5.5 5.1 4.7 4.0 3.5 4.9
##  [145] 3.5 3.5 5.3 4.4 5.4 4.8 4.3 3.4 4.7 4.7 4.6 4.3 3.2 6.0 5.3 4.7 5.1 2.4
##  [163] 3.3 3.7 3.3 4.9 6.2 6.6 4.2 4.5 3.4 4.8 6.1 5.4 4.6 4.3 5.7 4.7 6.0 4.4
##  [181] 3.1 4.8 4.8 4.2 4.9 5.0 6.1 4.3 4.5 3.9 3.6 4.7 5.6 3.5 3.2 5.2 4.4 3.5
##  [199] 5.4 5.9 6.2 3.8 4.0 5.3 5.1 3.0 4.9 5.4 3.9 4.7 5.1 4.7 4.9 5.5 4.0 4.4
##  [217] 3.7 4.4 5.1 4.1 2.8 3.3 5.4 2.9 4.9 4.2 5.4 5.5 3.9 4.1 3.1 4.0 5.1 4.4
##  [235] 3.0 5.4 4.9 4.6 3.9 4.2 5.2 5.9 5.3 6.2 5.9 4.3 4.8 5.5 4.2 4.1 4.8 5.0
##  [253] 3.7 4.7 3.0 3.4 4.3 3.4 3.7 3.7 4.1 5.2 5.1 3.8 4.7 5.9 4.0 4.3 5.6 6.0
##  [271] 2.9 3.3 3.2 4.7 4.2 3.9 4.6 4.6 5.1 3.7 3.6 3.9 2.9 3.8 5.5 4.5 5.6 5.9
##  [289] 3.9 4.0 4.7 3.1 3.4 4.7 3.6 4.3 4.7 6.5 4.2 6.1 4.9 3.9 5.0 3.9 5.6 5.6
##  [307] 5.5 3.4 4.7 3.6 4.3 3.5 4.4 4.8 3.5 2.8 4.4 4.1 2.8 5.0 3.5 4.3 2.6 5.4
##  [325] 4.8 3.5 4.3 3.7 2.6 3.5 4.1 5.7 4.2 4.5 6.1 5.1 4.6 3.6 3.9 4.7 2.5 4.1
##  [343] 3.5 3.3 3.5 3.0 4.6 3.8 4.5 4.0 6.2 5.7 5.0 5.1 5.6 5.7 5.6 5.2 1.9 3.2
##  [361] 5.7 3.8 4.4 6.8 4.8 4.7 4.5 5.7 5.1 4.9 5.5 3.3 5.7 5.2 5.2 5.4 4.5 3.9
##  [379] 3.3 4.4 3.6 4.4 2.8 4.2 2.9 4.9 5.5 4.4 4.0 2.8 5.0 4.7 4.5 5.0 3.2 3.3
##  [397] 4.2 3.3 6.2 4.1 3.4 6.1 4.1 5.5 3.4 5.7 3.4 4.7 5.2 5.4 4.2 3.3 4.5 3.5
##  [415] 4.2 3.7 4.7 4.6 3.9 4.7 4.9 3.9 5.5 5.1 5.3 3.9 4.4 5.6 4.6 4.6 5.1 4.2
##  [433] 4.4 2.2 4.5 2.8 5.4 4.5 4.3 4.4 4.9 3.7 2.6 4.6 5.3 6.3 5.4 6.0 5.4 3.9
##  [451] 4.9 4.1 3.4 6.0 4.9 5.1 6.4 3.2 5.7 4.1 5.9 5.6 4.2 5.4 6.5 3.1 5.0 4.2
##  [469] 5.9 5.3 3.8 4.5 4.6 3.5 4.2 5.2 4.1 5.5 5.5 3.0 5.0 4.4 4.7 3.3 5.2 5.2
##  [487] 6.4 3.1 4.9 3.4 5.6 5.6 5.3 4.9 4.9 6.2 5.3 3.0 4.6 4.8 3.4 5.8 4.1 4.5
##  [505] 4.5 3.9 4.4 4.0 4.3 4.1 4.8 5.2 4.2 6.1 5.1 4.9 2.3 5.0 5.2 4.0 3.5 6.5
##  [523] 4.8 4.4 4.1 4.1 4.2 3.1 7.1 2.8 4.0 2.8 4.8 5.0 5.1 3.2 4.5 2.2 5.6 5.4
##  [541] 3.3 4.5 5.3 4.9 3.9 4.3 5.3 4.7 4.2 3.1 4.2 3.8 6.4 4.7 5.1 4.6 6.1 4.6
##  [559] 4.9 4.2 3.1 2.7 3.7 4.6 4.9 5.1 5.2 4.3 5.6 4.6 3.3 5.4 6.0 4.4 5.1 5.2
##  [577] 4.6 4.1 4.3 5.5 4.4 4.4 5.8 4.0 5.1 4.1 4.1 4.3 5.3 3.6 6.1 1.5 5.2 4.3
##  [595] 5.1 2.0 4.0 3.6 2.7 3.9 3.3 5.6 4.3 5.1 6.4 5.6 3.8 4.6 4.0 5.9 3.4 4.1
##  [613] 3.1 3.1 3.7 2.7 2.7 3.8 5.3 4.1 4.9 3.9 5.1 4.1 5.4 4.8 5.8 3.0 4.7 3.4
##  [631] 4.1 3.9 4.2 4.8 3.9 3.7 6.4 4.7 4.5 4.6 3.9 3.4 4.8 2.6 5.5 4.0 4.5 4.0
##  [649] 4.7 6.8 4.2 4.4 4.1 4.4 3.1 5.1 3.2 6.2 4.0 4.8 3.2 5.1 3.7 4.2 4.7 4.8
##  [667] 5.3 5.7 5.2 4.4 3.7 4.4 3.8 5.3 3.5 6.1 5.8 3.7 5.0 4.0 3.3 4.0 6.3 3.8
##  [685] 4.1 4.2 5.8 2.8 5.3 5.9 4.3 3.7 4.1 3.6 3.7 4.3 2.9 5.9 4.2 3.6 4.4 4.2
##  [703] 3.8 3.6 5.4 5.1 4.2 5.3 4.6 4.2 3.7 4.8 4.2 4.7 4.0 3.8 4.5 4.8 3.6 6.1
##  [721] 4.2 5.2 5.9 4.8 4.5 4.2 4.9 4.0 5.2 5.6 4.9 3.5 3.7 4.6 6.5 4.8 5.0 5.9
##  [739] 3.7 3.9 5.4 3.8 5.0 4.1 5.8 4.3 4.4 5.0 2.8 5.4 5.0 3.1 4.3 5.8 2.8 3.1
##  [757] 2.9 5.1 4.0 4.6 4.3 3.9 5.1 4.1 5.8 4.6 4.4 2.9 5.3 4.8 4.3 4.5 5.8 5.5
##  [775] 3.1 5.9 5.0 5.1 3.5 4.6 3.3 5.0 4.2 4.5 4.3 5.3 3.8 5.0 6.1 2.4 3.8 3.8
##  [793] 5.9 4.2 3.8 3.6 3.8 4.9 5.4 3.8 5.2 5.2 4.6 2.9 6.0 6.0 2.5 3.4 6.5 3.9
##  [811] 3.9 4.6 5.0 3.7 5.2 3.4 3.7 5.7 5.1 3.0 4.4 5.8 5.8 6.1 3.3 4.1 4.0 5.3
##  [829] 5.1 4.5 4.4 4.2 6.2 6.4 3.9 4.7 4.7 5.4 5.7 5.0 4.3 5.9 4.6 6.0 3.3 5.8
##  [847] 5.5 5.4 6.1 2.6 5.4 5.5 4.2 5.6 4.8 5.6 3.9 3.8 5.6 3.7 5.7 4.0 2.9 4.9
##  [865] 4.9 4.9 5.0 4.6 3.0 5.7 3.4 4.4 5.2 2.7 5.4 5.4 3.5 4.0 4.3 3.9 3.6 5.1
##  [883] 4.6 2.7 4.7 5.4 4.2 4.2 4.6 6.0 3.2 4.9 5.0 5.1 5.1 4.4 3.4 3.2 3.4 4.3
##  [901] 3.9 4.9 5.4 3.9 2.4 4.0 2.6 4.5 4.6 4.0 2.2 3.8 5.9 4.6 4.4 4.5 4.2 3.4
##  [919] 4.8 6.5 5.3 3.0 4.2 5.3 5.4 5.0 5.4 5.3 4.6 3.8 5.9 5.7 3.5 3.0 5.1 5.2
##  [937] 3.0 4.5 5.6 3.9 3.5 5.0 4.5 4.3 4.7 5.2 4.5 2.7 4.5 4.8 4.2 4.6 4.7 4.7
##  [955] 3.0 4.5 4.8 4.4 4.3 5.6 5.0 5.3 4.9 3.9 4.3 3.6 3.2 4.3 5.9 3.7 3.3 3.4
##  [973] 4.3 5.8 5.8 4.8 4.1 3.5 2.7 4.5 4.8 4.4 5.2 5.5 4.8 4.6 4.8 5.0 4.2 4.5
##  [991] 4.8 3.9 3.4 3.2 6.1 4.8 3.2 4.3 3.0 5.2
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
##    [1] 4.9 4.4 4.0 5.6 2.8 3.6 4.5 3.7 3.6 2.2 4.2 4.2 4.1 3.5 3.4 3.6 4.6 4.5
##   [19] 3.5 2.9 3.3 4.2 5.8 5.0 4.9 3.6 4.2 3.0 5.8 4.5 4.7 4.3 4.5 3.6 5.7 5.9
##   [37] 6.0 5.9 4.7 4.2 4.2 5.9 4.1 4.7 4.3 3.5 3.5 4.2 4.2 4.8 6.2 4.1 3.3 5.3
##   [55] 4.0 5.1 3.5 4.9 4.6 3.9 5.4 4.2 4.7 3.8 4.2 4.6 4.0 5.0 4.2 3.0 5.6 3.5
##   [73] 6.6 3.6 3.3 3.6 4.9 2.5 2.8 4.9 5.3 4.8 6.5 3.7 4.6 5.5 4.0 2.7 3.7 3.2
##   [91] 4.6 2.6 3.1 3.7 3.6 3.8 4.3 5.4 5.2 4.5 6.5 5.8 4.3 4.6 3.9 3.8 4.7 5.4
##  [109] 2.5 4.8 5.0 4.2 4.5 4.8 4.1 5.6 5.8 4.6 5.2 4.1 4.7 3.3 3.4 2.6 4.1 3.5
##  [127] 4.2 3.9 6.8 3.3 4.1 4.5 4.5 5.8 3.5 5.1 5.2 3.2 3.8 3.2 4.6 5.1 4.7 3.9
##  [145] 5.3 4.3 4.8 5.0 5.6 5.2 5.2 4.4 3.7 3.4 4.5 5.8 4.4 4.1 4.7 5.5 4.5 5.4
##  [163] 5.5 6.0 4.0 4.7 3.0 5.0 4.7 4.6 4.2 4.7 5.3 4.2 4.4 5.1 3.9 5.1 4.3 5.1
##  [181] 3.0 6.4 5.2 4.7 4.8 5.6 4.2 4.0 4.9 4.2 4.6 4.6 4.2 3.0 6.0 4.3 4.8 5.6
##  [199] 2.3 4.7 4.2 4.8 5.9 6.1 5.7 4.7 7.3 4.6 6.7 3.5 4.6 4.1 2.6 3.9 4.7 4.9
##  [217] 4.7 5.2 4.1 5.0 4.5 3.6 5.5 5.7 4.4 5.1 3.4 4.4 3.7 3.0 4.5 3.8 4.5 4.1
##  [235] 5.0 4.1 3.1 5.1 4.8 4.7 5.3 4.4 3.4 2.9 4.4 5.3 5.0 3.6 4.3 3.8 5.8 4.4
##  [253] 6.2 3.2 6.0 6.1 5.5 6.2 3.6 4.9 4.5 2.4 3.4 4.5 4.7 5.4 5.4 5.1 2.9 6.2
##  [271] 5.2 3.6 4.8 5.2 4.9 5.8 5.0 4.9 5.1 5.9 4.7 4.5 5.5 4.8 3.8 5.3 4.0 4.2
##  [289] 5.2 4.8 4.0 3.6 5.3 4.4 3.0 4.4 5.0 3.5 3.3 5.0 3.7 3.2 4.5 5.2 4.3 2.8
##  [307] 3.6 5.6 4.5 4.0 4.1 3.4 5.7 5.3 4.8 4.6 4.3 2.2 3.9 4.6 3.6 3.6 4.6 4.7
##  [325] 4.0 4.3 2.2 5.0 5.6 4.6 4.2 3.9 4.2 4.4 4.1 5.6 3.5 6.6 5.3 5.9 4.1 5.0
##  [343] 3.2 4.4 6.7 6.0 4.4 4.5 5.5 5.5 4.0 3.8 3.9 5.1 6.4 4.4 4.6 3.7 4.2 4.3
##  [361] 3.6 4.6 4.9 4.0 4.7 3.5 4.4 4.7 4.0 4.9 3.1 4.3 6.4 3.3 4.2 3.9 4.9 5.4
##  [379] 4.9 4.9 4.0 4.3 4.8 5.3 4.0 3.0 5.7 5.2 4.2 4.1 3.8 4.6 5.6 3.8 5.3 5.2
##  [397] 5.1 4.6 5.8 3.0 5.2 3.5 3.7 2.4 4.8 5.4 6.5 5.0 4.6 3.5 4.4 4.8 2.6 5.1
##  [415] 6.6 6.5 6.2 2.6 4.6 5.6 3.6 4.7 4.6 3.8 4.2 4.4 4.2 3.4 3.9 5.0 4.8 4.1
##  [433] 3.7 4.2 5.1 4.0 4.9 3.5 5.1 5.1 5.2 4.4 4.7 4.3 5.0 5.0 4.3 4.5 4.0 3.9
##  [451] 4.2 3.9 3.1 3.0 5.1 5.4 3.3 6.1 3.9 4.5 4.1 4.7 5.0 4.6 4.0 4.9 4.5 4.3
##  [469] 5.6 3.6 5.2 3.7 3.0 3.9 4.8 4.8 4.6 5.0 2.3 5.4 4.0 4.8 4.0 3.5 4.1 6.0
##  [487] 5.2 5.5 3.3 4.6 4.7 6.5 3.9 4.0 3.4 4.6 5.1 4.7 4.2 3.5 5.5 3.2 5.1 2.8
##  [505] 3.3 3.8 3.1 4.9 4.6 4.4 3.1 4.2 5.2 6.3 5.3 4.6 3.0 4.0 5.2 3.8 4.9 5.9
##  [523] 4.9 3.8 4.3 4.4 4.0 6.3 4.2 4.6 4.1 2.5 2.9 3.3 3.9 3.8 5.2 4.5 4.0 4.2
##  [541] 3.5 5.9 4.4 5.9 3.5 4.3 4.6 4.1 3.5 4.7 4.7 3.4 3.4 6.0 6.5 5.3 4.0 4.0
##  [559] 3.8 4.1 5.3 5.1 4.6 4.1 4.1 4.6 4.3 3.2 4.5 5.5 5.2 5.2 4.3 3.9 3.7 3.2
##  [577] 3.7 5.9 5.5 4.3 4.8 4.0 4.2 5.1 5.1 3.3 4.1 4.6 3.1 5.8 5.9 6.6 5.6 3.8
##  [595] 4.2 3.9 3.7 4.2 4.8 3.3 4.9 2.9 3.4 4.7 6.1 5.1 4.3 3.7 3.9 4.4 5.2 5.3
##  [613] 5.3 5.1 3.7 4.8 4.6 5.1 4.7 5.6 5.6 4.0 5.0 3.6 5.0 4.1 4.9 5.1 5.0 3.7
##  [631] 4.2 5.0 3.7 4.2 4.9 5.4 5.3 5.2 4.9 4.2 3.7 5.3 4.1 4.4 6.0 5.2 5.9 5.3
##  [649] 4.2 4.8 4.4 3.5 5.4 6.0 2.5 5.4 5.2 5.4 3.6 3.6 4.0 5.4 4.4 3.1 4.5 5.7
##  [667] 4.9 3.5 4.1 4.4 7.5 5.1 4.7 2.1 3.6 5.6 4.7 2.8 4.8 5.5 3.8 4.3 4.7 5.4
##  [685] 4.0 5.0 4.5 3.7 4.1 3.6 4.4 4.2 3.3 2.8 6.0 4.0 4.1 4.4 4.3 5.7 5.9 5.0
##  [703] 4.8 5.3 4.9 5.1 4.8 5.8 4.5 3.9 5.6 5.9 4.1 4.6 3.7 4.4 3.0 5.4 4.5 3.2
##  [721] 3.3 4.4 4.6 6.2 5.2 4.6 5.9 5.4 6.4 3.1 5.1 2.7 5.1 4.8 4.3 5.8 5.5 5.2
##  [739] 4.1 5.3 4.1 4.4 4.0 5.4 5.2 5.6 4.9 3.5 4.7 2.9 6.0 5.2 3.7 5.1 3.4 4.0
##  [757] 5.4 2.8 3.5 3.8 5.7 3.7 2.5 5.3 4.7 4.1 4.1 3.7 5.3 4.6 4.3 5.7 5.1 3.0
##  [775] 4.5 4.6 5.5 4.8 5.0 5.4 4.3 5.1 6.3 4.0 5.1 4.9 4.3 4.6 3.9 5.7 4.1 4.0
##  [793] 3.9 3.5 6.4 4.0 3.4 4.9 5.0 3.8 3.3 4.9 5.2 5.7 3.0 3.9 4.6 4.3 5.6 4.2
##  [811] 3.6 4.6 5.4 3.7 3.1 5.1 4.7 5.8 4.9 4.1 5.1 3.5 3.8 4.7 5.4 4.4 4.4 4.4
##  [829] 3.8 2.8 4.0 4.6 3.8 3.7 4.2 4.8 5.1 3.7 3.2 4.8 5.5 3.4 2.9 4.5 6.3 3.5
##  [847] 4.7 2.9 4.4 5.6 3.7 5.5 5.7 4.0 2.9 4.1 5.4 4.4 4.1 4.9 5.4 5.0 4.3 4.7
##  [865] 4.3 4.1 4.8 3.8 5.2 5.1 4.5 5.1 2.8 4.3 5.8 3.7 3.5 6.1 5.1 4.3 5.4 5.6
##  [883] 4.8 6.4 3.3 4.7 5.8 5.3 3.5 5.7 5.6 2.6 4.0 5.3 7.0 3.0 5.1 4.0 3.1 3.7
##  [901] 4.0 4.6 4.3 3.3 5.0 4.6 4.0 4.1 2.8 4.8 4.5 4.9 3.9 6.3 5.1 5.8 3.9 4.2
##  [919] 5.0 4.1 6.2 4.1 3.4 2.5 4.0 4.2 4.8 4.5 3.8 4.1 5.0 2.8 4.6 4.6 5.0 2.5
##  [937] 4.8 5.6 4.3 3.9 5.7 3.5 4.5 4.6 4.4 5.2 4.2 6.4 5.4 4.8 3.6 3.8 3.1 3.0
##  [955] 4.1 5.0 3.5 5.1 5.7 4.9 5.9 4.4 5.4 4.5 4.5 3.6 6.6 4.8 5.0 5.5 5.8 4.1
##  [973] 3.8 4.1 4.7 5.0 4.9 3.8 3.3 5.4 6.0 4.0 4.1 5.7 3.9 6.5 5.3 5.5 3.8 3.6
##  [991] 5.5 4.9 4.8 4.5 4.1 4.6 3.4 4.6 4.4 5.1
## 
## $func.thetastar
## [1] 0.0011
## 
## $jack.boot.val
##  [1]  0.43600000  0.39604520  0.28618619  0.16586826  0.08542857 -0.03746479
##  [7] -0.13945205 -0.17750000 -0.41573034 -0.54023669
## 
## $jack.boot.se
## [1] 0.9393186
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
##    [1] 5.4 4.2 4.5 5.0 3.2 4.0 2.0 4.4 4.9 4.7 4.3 4.2 4.6 3.9 5.3 4.1 4.5 4.4
##   [19] 3.7 4.5 4.0 1.9 5.4 3.9 3.7 4.7 5.0 4.2 6.3 4.3 4.9 4.8 3.8 4.4 4.5 3.6
##   [37] 4.9 4.4 5.1 2.7 4.0 4.0 3.4 4.9 2.9 4.0 4.7 4.8 4.2 5.7 3.4 4.5 6.2 4.5
##   [55] 3.5 2.9 4.4 4.4 4.1 6.1 4.1 5.3 3.3 4.3 5.1 5.5 4.8 6.1 5.9 4.4 3.7 5.2
##   [73] 6.0 5.1 3.5 5.2 6.3 5.0 4.0 5.3 6.2 3.4 4.7 3.6 4.3 5.1 6.0 3.5 4.8 4.6
##   [91] 5.1 6.1 4.8 3.8 5.3 4.4 4.4 4.4 4.8 6.0 5.1 5.0 4.8 4.5 4.1 4.0 4.6 3.3
##  [109] 5.1 5.6 4.0 4.6 4.2 3.8 5.3 4.0 3.7 3.6 4.5 4.6 4.2 3.8 4.0 5.2 4.8 4.0
##  [127] 3.6 3.8 3.7 5.5 4.7 3.6 5.1 5.1 3.4 5.8 4.2 5.3 5.0 5.5 4.9 4.8 3.9 3.7
##  [145] 5.8 5.2 4.9 4.6 5.2 4.0 5.5 3.7 3.4 5.2 4.1 4.6 3.6 3.6 4.3 4.7 3.6 4.0
##  [163] 5.5 4.2 4.6 6.7 5.4 2.3 3.9 4.0 4.6 4.0 3.9 4.6 4.6 4.1 4.4 5.0 4.4 4.2
##  [181] 3.6 4.6 7.2 5.5 4.1 4.1 4.8 5.0 3.3 2.9 4.7 5.8 3.1 4.2 5.2 5.4 5.4 6.1
##  [199] 4.1 4.6 4.3 5.1 4.1 4.2 4.5 5.1 4.5 6.2 5.3 4.3 3.6 4.8 4.6 4.4 3.7 5.2
##  [217] 3.6 4.8 5.2 3.4 3.8 3.9 5.2 4.1 4.7 6.9 4.3 4.2 4.6 4.6 6.7 4.3 4.0 3.1
##  [235] 3.8 4.0 4.4 3.8 4.1 5.2 5.1 3.8 3.8 3.3 4.2 3.7 4.7 3.9 3.7 6.4 3.5 4.5
##  [253] 4.4 3.5 4.9 5.0 4.7 4.0 6.0 2.7 3.8 6.1 5.9 5.4 4.9 2.8 5.6 5.6 4.2 4.5
##  [271] 6.2 3.7 5.2 5.1 5.5 3.8 5.0 5.9 3.5 4.4 2.9 4.2 5.6 4.1 4.0 2.5 2.1 5.1
##  [289] 4.4 5.6 3.7 4.9 3.2 3.8 5.1 4.5 6.6 3.4 4.7 3.7 5.7 2.4 3.2 4.5 4.5 4.4
##  [307] 2.7 5.7 5.3 5.1 4.0 3.4 2.6 3.7 5.0 4.3 3.3 4.2 4.0 3.9 4.4 4.4 4.5 5.7
##  [325] 5.1 5.7 4.3 2.7 3.1 5.1 4.1 4.2 3.6 3.7 5.2 5.3 4.1 6.1 3.5 5.5 5.6 4.2
##  [343] 4.0 4.6 4.2 4.8 3.5 4.6 4.2 4.4 5.8 6.4 4.6 4.9 4.7 4.1 6.0 3.9 4.6 5.9
##  [361] 4.1 3.7 3.0 4.5 3.8 5.8 4.1 5.0 5.0 2.1 3.2 4.1 3.9 4.1 4.0 5.8 5.4 4.9
##  [379] 3.6 5.2 4.7 3.9 6.1 4.7 4.9 3.8 3.8 4.6 5.0 6.7 4.3 4.2 5.7 5.0 6.0 5.1
##  [397] 5.0 5.4 5.1 5.6 5.1 4.1 4.1 5.5 3.8 3.7 3.6 3.9 5.8 3.3 3.8 6.0 5.9 4.5
##  [415] 6.1 5.2 5.0 5.2 4.1 3.1 4.2 3.8 3.7 4.3 5.9 5.0 3.3 4.8 4.0 4.7 5.2 3.8
##  [433] 4.4 4.0 3.6 4.0 4.3 4.6 4.3 4.1 3.5 4.1 4.9 4.7 4.7 4.8 4.0 5.8 4.3 3.5
##  [451] 5.1 3.2 6.0 4.2 4.9 6.6 4.4 4.7 5.7 4.8 4.2 3.5 4.8 3.0 5.2 4.2 5.6 5.8
##  [469] 4.6 5.5 5.3 4.9 3.4 4.4 3.8 5.2 4.1 3.7 4.8 4.4 4.2 4.9 3.2 4.2 4.1 3.3
##  [487] 3.1 4.3 3.3 5.0 4.7 6.4 5.6 3.2 5.3 5.2 3.9 5.2 4.7 4.6 2.3 4.7 4.4 6.0
##  [505] 5.3 4.8 4.0 3.8 4.9 5.3 5.6 3.1 5.6 5.4 4.8 3.7 4.4 5.5 5.3 5.9 5.4 5.7
##  [523] 4.4 3.6 4.7 5.4 5.6 4.2 6.0 6.4 6.3 4.6 3.5 5.1 5.7 5.3 4.9 3.4 4.8 4.7
##  [541] 4.3 4.5 4.3 5.9 5.0 4.2 4.7 4.4 4.5 4.6 4.1 4.1 3.3 4.5 4.6 4.3 5.3 6.3
##  [559] 4.1 4.9 3.1 5.3 3.5 4.6 5.4 5.6 6.1 3.8 5.9 4.1 5.1 2.3 4.5 3.3 3.2 3.7
##  [577] 4.4 3.3 3.7 4.6 4.1 4.1 3.8 5.1 4.7 3.9 3.9 4.2 5.5 4.8 2.6 4.6 6.7 2.7
##  [595] 5.9 4.9 3.6 5.0 3.5 4.4 3.9 5.8 5.3 4.8 4.9 6.1 5.5 5.2 3.9 4.9 5.1 5.3
##  [613] 4.6 5.5 4.4 4.0 4.3 4.2 4.2 4.2 4.6 4.4 4.8 3.6 5.1 3.8 4.3 4.6 5.8 4.4
##  [631] 3.3 4.9 5.4 3.7 5.6 5.4 4.8 5.9 3.8 4.4 3.0 4.3 5.8 5.2 5.1 3.4 6.0 3.7
##  [649] 5.6 4.6 4.3 4.8 5.7 5.5 3.8 5.4 4.8 4.9 4.3 5.3 4.1 5.0 4.5 6.0 4.1 6.3
##  [667] 4.6 4.5 5.8 3.9 4.8 4.4 4.5 4.1 3.3 4.3 4.4 4.0 3.6 2.8 5.8 5.8 5.5 4.7
##  [685] 4.3 4.5 5.1 4.5 4.9 4.4 6.0 4.2 4.1 4.4 5.7 4.5 4.6 4.6 4.3 4.9 4.6 4.6
##  [703] 4.7 3.4 4.6 4.7 5.4 4.7 5.8 4.1 3.1 4.0 5.8 4.1 3.8 5.1 4.7 4.5 6.3 4.5
##  [721] 5.7 3.6 3.9 5.9 5.5 4.0 4.0 5.2 5.0 2.8 4.5 6.5 4.5 4.2 4.5 5.5 4.2 3.0
##  [739] 3.8 5.1 4.2 3.9 6.4 2.5 5.9 3.9 5.2 5.2 3.5 4.3 5.5 3.7 3.4 5.0 4.6 4.4
##  [757] 4.7 3.5 6.0 4.8 2.8 5.0 4.0 6.1 4.6 6.3 5.7 6.3 4.8 4.1 4.9 4.4 4.5 6.0
##  [775] 4.7 4.4 3.9 4.0 3.4 2.8 4.7 4.8 5.4 4.3 4.6 4.2 3.3 3.5 6.6 4.6 4.9 5.4
##  [793] 4.8 4.8 3.1 5.4 3.6 5.9 6.0 4.9 5.4 3.6 3.8 5.3 5.0 5.2 3.3 4.1 4.3 4.8
##  [811] 5.1 5.3 5.8 5.8 5.3 4.3 3.5 3.7 6.0 3.8 4.8 4.4 6.5 3.8 3.9 4.1 3.6 5.3
##  [829] 5.2 4.6 3.9 5.4 5.4 5.2 4.5 3.8 3.3 4.1 4.8 5.3 3.6 5.1 4.6 4.2 4.4 4.4
##  [847] 4.1 5.8 6.0 3.7 5.9 5.8 3.3 5.2 5.0 4.4 5.2 3.7 3.7 5.2 2.8 4.5 4.0 5.2
##  [865] 5.6 4.1 4.6 4.9 3.9 6.1 2.6 5.1 4.0 4.2 4.3 4.7 5.4 5.9 4.2 6.1 4.6 3.5
##  [883] 4.5 4.6 4.9 3.9 5.3 3.5 4.7 4.6 4.4 4.2 3.4 4.8 2.5 5.2 3.2 4.1 5.0 5.5
##  [901] 4.0 4.6 3.0 4.3 4.1 4.0 3.5 3.8 3.9 4.3 3.5 5.0 3.5 4.5 4.0 1.6 3.6 4.8
##  [919] 3.6 4.9 5.1 4.3 4.2 4.8 4.6 5.4 4.4 3.4 4.0 4.0 4.3 4.5 3.9 5.4 5.6 4.4
##  [937] 5.2 5.1 3.8 3.1 2.9 4.5 6.3 4.7 3.7 4.5 4.6 6.9 5.3 2.8 4.0 3.2 3.5 5.2
##  [955] 3.3 4.2 3.1 4.4 5.7 4.3 3.5 3.9 4.1 3.9 4.6 3.0 4.1 5.6 4.9 2.9 4.3 3.5
##  [973] 3.8 2.9 4.7 4.1 4.0 3.4 5.1 3.0 4.7 3.8 5.1 4.7 4.1 4.7 5.5 3.7 3.7 4.5
##  [991] 2.8 3.4 3.4 4.0 3.6 6.0 5.1 4.3 3.8 5.1
## 
## $func.thetastar
##   72% 
## 5.028 
## 
## $jack.boot.val
##  [1] 5.5 5.3 5.3 5.3 5.2 5.1 4.9 4.6 4.5 4.6
## 
## $jack.boot.se
## [1] 1.013361
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
## [1] 0.5558785
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
##   6.214924   8.281290 
##  (2.708083) (3.758269)
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
## [1]  0.1485309  0.7175859 -0.4277769  0.4864290 -0.3427941 -0.2221421
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
##    [1]  1.3963375426  0.0862701152 -0.0420689105  0.5240216073  1.5817441923
##    [6] -0.2439537190  1.6276433361  0.8661028837  1.1587262516  0.5924325069
##   [11]  0.5696298324 -0.6481470618  1.5108686384  0.7127244843  0.2196441207
##   [16]  0.2444347203  0.1579198402  0.5319620473  0.7575715828  1.6079687466
##   [21] -0.5662806486  0.9710167119 -0.7146398122  0.7532046440  0.9707397098
##   [26]  0.6009478592  0.2399500299  0.4197164241  1.5735149578 -0.0046525813
##   [31]  1.2174169361  1.1534728352  0.2217741963  1.0598497272  0.5314060236
##   [36]  0.3007496349  0.5437127014  0.1132732070 -0.3159517977  0.5213808687
##   [41]  0.0016063229  0.4162659035  0.3378253150 -1.0547917407 -0.1617573833
##   [46]  0.8116915266  0.2007714229  0.3548771614  0.1031563456 -0.3587580132
##   [51]  0.3987286222  0.8081275994 -0.3134352209  0.2480956163  0.3462414804
##   [56]  0.6336772351 -0.3657335395  0.3885541340  1.0840876416  0.3507891415
##   [61]  0.1098057996 -0.1724339584  0.3967277990  0.4596526817  0.3251742369
##   [66]  0.3463597163  1.2156861461  0.0885541059 -0.3341936590 -0.1638976759
##   [71] -0.4503129175  0.6660517807  1.2548679509  0.1682791225  1.0402279436
##   [76]  0.0836260715  1.5250726974  0.8910083446 -0.5148763256  0.9821038246
##   [81]  1.1992397648  0.1020168787  0.8520405572  0.4839952075  1.8084779568
##   [86] -0.2544239976  0.5301833709  0.6634665406  0.1229221687  0.6301421159
##   [91] -0.5746248422  0.1626501627  0.8839101006  0.8299428672  0.4169285889
##   [96] -0.1092179397 -0.1554321638  0.3734309758  0.3323156063  0.8937619223
##  [101]  0.3618266250  0.8671890411  0.9764018714  0.5244381299  0.4606232433
##  [106] -0.2308407010  0.9574784674  1.1523898999 -0.1642450713  0.7438779666
##  [111]  0.2972025366  0.9737379871  0.4519231594  1.0298223568  0.3174458450
##  [116]  0.3967038679 -0.0757873749  0.4343279056  0.6399697858  0.1634258550
##  [121]  0.6452027435  1.5135289027  0.7608490465  1.4451360660 -0.2116893910
##  [126] -0.2283430277  0.2974319961  0.5668723290  1.2097611104  0.8813202112
##  [131]  0.3548771614  0.8354200713 -0.6090324931 -0.7131470103  0.3825896439
##  [136]  0.5628025956 -0.0024148999 -0.2865822638  0.3534548700  0.7415306622
##  [141] -0.7579015060 -0.2858354043 -0.0504323068  0.1478997273 -0.5277408962
##  [146]  0.6481833238  0.3492776170  0.0117994347  0.1172707705  0.8415324307
##  [151]  0.2709004857  0.4499686830  0.0484969862  0.6967769019  0.8568722526
##  [156]  0.7486123213  0.7797317676 -0.8550496289 -0.0380646784  0.9390577428
##  [161]  0.1512641188 -0.3572376073 -0.6195995747 -0.2621251959  0.5834399878
##  [166]  0.0206012506  0.4728930418  0.7635111794  0.4538312606  0.2147120065
##  [171]  0.3440619535  0.7576275517 -0.8242606625  0.3592432134  0.3326827874
##  [176]  0.9662859494  1.9240800660  0.7705577135  0.9183760280  0.9743961950
##  [181]  0.2749683448 -0.4300375355  0.3713768623  0.0643149312  0.2309975756
##  [186] -0.2186824254  0.5399145609  0.9846376100  1.1785533811  1.1691125962
##  [191]  0.6078329969  1.0169414770  0.6478037170  0.0011211948  0.2024519998
##  [196] -0.0836490231 -0.4963018507  0.4755120493 -0.6548765526  0.8479318718
##  [201] -0.0910545821  0.1803327481  0.8684198528  0.5418815494 -0.0227972988
##  [206]  1.0112931053  1.4451360660  0.5121185316  0.3888617679  0.7024717888
##  [211]  0.7895598723 -0.6654421293  0.0489102579 -0.1228131050  0.5867898618
##  [216]  0.3777325624 -0.2251821146  0.0214832448  0.8590255719  0.6776926249
##  [221]  0.2772797346  0.5433138394  0.4368990766  0.4507566329  0.1731814399
##  [226]  0.2776132007  0.5962747340  0.1875268305  0.3985949096  0.5063932502
##  [231]  0.1151506337  0.3394000747  0.8451841245  0.8242770076 -0.1479063038
##  [236]  0.8026801282  2.2357770850  0.7381707631  1.0095694698  0.1825822977
##  [241]  1.1937607327  0.3654064731  0.9436259037 -0.2400262179  0.1507474177
##  [246]  0.5407873432  0.6426149087  0.7677812167  0.3554094515 -0.0903177320
##  [251]  0.4921842381  0.3907081239  0.2251476674  0.8104619545  1.4557692787
##  [256]  0.0933398186  0.2670341224  0.5837992736  0.7259944131  0.2335811259
##  [261]  0.0614520089 -0.5785522938  0.4644415411  1.1491738779  0.4155403854
##  [266]  0.2331811180  0.0593147587  0.7479603425  0.4193520204  1.1659809777
##  [271]  0.1458195355  0.0934379368  0.6498322007  0.3462414804 -0.1347520954
##  [276]  0.2208474799  0.9686603381  1.3981692109  0.4059098160  0.0115012757
##  [281]  0.3741128930  0.3193416913  0.5320703758  0.1228241411  0.0149991164
##  [286]  0.0637753313 -0.0210208318  0.3671136382  0.7001008134  2.3852465212
##  [291]  0.5834534509  0.4197164241  0.2898489326  1.0156002084  0.0289033930
##  [296] -0.1506589989 -0.0491656733  0.9885367265 -0.1295641749 -0.4604218126
##  [301]  0.2424740195  0.2101505545  1.3591081951  0.6189821642 -1.0001116319
##  [306]  0.0747509702  0.3991926593  0.6639446873  0.6032836090  0.9753651016
##  [311] -0.0457052750  0.5478313915  1.3855647202  0.7774950287  0.8174302411
##  [316]  0.0407288845  0.9878050536  0.8908943571  0.3026114549  0.7300810552
##  [321]  0.1466628151  1.7937757417  0.1088304049  1.0095694698  0.1048079606
##  [326]  0.4074043210  0.5988169356  0.4972346160  1.0065207501  0.7668073101
##  [331]  0.1326562735  0.2037520811  0.2984374069  0.4957553660  0.0464415713
##  [336]  0.2430323464  0.8498722034  0.8250397447  0.1712648326  0.7841135292
##  [341]  0.9809313166  0.0020316180 -0.1255741657  1.2000959758  0.5108113633
##  [346]  0.7408582528  0.0350304745 -0.1377085893  0.6362507355  0.4826747961
##  [351]  0.7278461462  0.7197331809 -0.3741105601  0.7689416559  0.1823563177
##  [356]  1.0658791710  0.3864123292  0.1492648817  1.7920406488  0.0682689806
##  [361]  0.5184513276  1.2045128622  1.1659857056  0.2716928416 -0.1442908703
##  [366]  0.1220771664  0.9824132452  0.2863957933  0.5911895967  0.4456890784
##  [371]  0.3646897372  0.2276784040 -0.0816304636  0.4708777464  0.1540181376
##  [376]  0.0417286929  0.3543843969 -0.3201652881 -0.1230172971  1.7692668237
##  [381]  0.5333157349  0.3942183071  0.1750487726  0.9251783396  0.3318189125
##  [386] -0.9322176515  0.7914959539  0.5693100550  0.3909335977 -0.3168549477
##  [391]  0.9999842384  0.2243828531  1.1140634850 -0.2324600836  0.0739851743
##  [396]  0.7614362669  0.6457740379  0.8296055576 -0.1202467144  0.3552913285
##  [401]  0.3713321502  1.5570324704  1.4756947639  0.1106760821  0.1654287261
##  [406]  0.2620209880  0.5056907675  0.4019286411  0.4176316808  0.0544064729
##  [411]  0.9800291112 -0.0387352882  0.8162803813  0.8110815757  1.1639102451
##  [416] -0.2961035324  0.0652604646 -0.1717873724  0.2867619682  0.9708137250
##  [421]  1.0398061607  1.3459345224  0.1988348085  1.2241349949  0.2647312921
##  [426]  0.1151506337 -0.5852415037  0.1883127582  0.0206012506  0.1787461143
##  [431]  0.7451964016  0.3679205715  1.3868788918  1.5621888735 -0.3709745123
##  [436] -0.2810158391  0.0644900578  0.2713156952  0.5183610740  0.0827174914
##  [441]  1.2245445696  0.1655091404  1.1223003902 -0.0338367758  0.4064953835
##  [446]  0.4525074667  0.7306182630  1.2650262785  0.3913816625 -0.0430035960
##  [451] -0.4454580781  0.3392272041  0.5582810557  0.4315392205  0.7636528207
##  [456]  0.2697384782  0.2584661242  1.0645259986  0.7642574915 -0.7092534263
##  [461]  1.2026062566 -0.3388268584 -0.2826038672  1.4413794312 -0.2882002193
##  [466]  0.3867116869  0.9464101137  0.2578766573  1.3702144226  0.7563375312
##  [471]  0.3327544308  0.8197854146 -0.4959333062  1.1787931819 -0.2068630547
##  [476] -0.2616479735  1.1734430340 -0.1425700953  0.5820402530  0.3805816745
##  [481]  1.0093139076 -0.2771499236  0.3515599944  0.8188452027 -0.1761101129
##  [486]  0.4369870728  0.4324078911 -0.0419857699 -0.4936536755 -0.2025557610
##  [491] -0.6113632690  2.0470879181  1.3254908569 -0.7714371747 -1.2234538308
##  [496]  0.0825009500  0.9471922835  0.2704034751  0.0548873430  0.3889004964
##  [501] -0.2849113021  1.1554292196  0.7129326794 -0.1916765052  0.0628966757
##  [506] -0.4726189547  0.6121442998  0.5867835833  0.3628222314 -0.0958037451
##  [511] -0.3689858355 -0.2163876405 -0.2710459854 -0.0433883813  0.3928348534
##  [516] -0.9029546233  0.0343839209 -0.1750121879  1.0759009456  1.5408572804
##  [521] -0.9399141027  0.3516660122 -0.0803698289 -0.0060816629  0.9431263725
##  [526]  0.0734529995  0.4167586589  0.3682374844  0.6634069899  0.9705231526
##  [531]  2.2443792690  0.5035778580  1.1859128107  0.5771426327 -0.4119522465
##  [536]  0.6766841264  0.8200226744  0.2945174811  0.5826582815  0.9066889266
##  [541]  0.5668227783  0.1402832381  0.0689086744 -0.3444212997  0.1610608799
##  [546]  0.0604561080  0.4029533941  0.4017673772  0.6771502560  0.6545302974
##  [551] -2.3000751209  0.2769671902  0.1798933469 -0.0622782195  0.4102262134
##  [556] -0.2170764047  0.1024869997  0.1680118324  0.7393114670 -0.0543630320
##  [561] -0.4875416187 -0.0631892745  1.2333294000  0.6750118599  0.5691732774
##  [566]  1.0430056294  0.1406367374  0.5634751916 -0.1074708522  0.6747193697
##  [571]  0.0545118166  0.7534337921  0.3233514503  0.3365366937  0.7452545425
##  [576]  0.9338404785  0.7300810552 -0.3821086427 -0.1244078592  0.7153442500
##  [581]  0.9846205474  0.0115158330 -0.3304345585  0.6062312274  1.0140559798
##  [586] -0.0947921672  0.3763983680  1.6499766207  0.8464599181  0.3165200563
##  [591] -0.2864643550  0.7932020236  0.3911086622  0.7052676099 -0.1017596575
##  [596]  0.1626876534  0.0484969862  1.0137681154  0.5944886573  0.5546127691
##  [601]  1.2971636887  0.7636384373 -0.2635855593  0.1323898736 -0.0552592304
##  [606] -0.4256753263  0.1769874185  0.7234724302 -0.2639336810 -0.1805681387
##  [611] -0.1023142281  0.2194985375  1.0158248871  0.6382081085  0.1581191374
##  [616] -0.7676060022 -0.0185061605  0.9641797369  0.9040173148  0.7009871129
##  [621]  0.1674248177  0.0409028604 -0.0048468646  0.8857047818  1.2644207085
##  [626] -0.2762297488 -0.2950800261  0.4386607959 -0.2771956689  1.1582938240
##  [631]  0.1983048860  0.8913776248  0.2450351916  0.5044115776  0.5281679632
##  [636]  0.0893812699  0.2061954154  0.8319808146 -0.6263198283  0.7972376215
##  [641]  0.3449281229 -0.1120995716  0.4868002266  0.2355487867 -0.1007644650
##  [646] -0.0130987817  0.6075950170 -0.3765457476  0.3216759746 -0.2258476498
##  [651] -0.0206375938  1.9327325017  0.2801155807  0.2048467134  0.2064185178
##  [656]  0.1169610153  1.1251232386  0.9950001609  0.0502000924  0.3186577014
##  [661]  0.6995978989 -0.1113147176 -0.3727869020 -0.1263496667  0.2889107561
##  [666]  0.5704838670  0.0512175738  0.1442440741  0.2146310561  0.4262199203
##  [671]  0.0592297204  0.2587390335 -0.2158042358  0.2463208363 -0.3721105274
##  [676]  0.7606323321  0.8552268576  1.2234221947  0.5603568974  1.7945277385
##  [681]  0.5964287367  1.4232023731  0.5175877176  0.8563021500  0.1989028892
##  [686] -0.2052903835 -0.3708114238  0.1533340257  0.6536720860 -0.4862063034
##  [691]  0.9968163203  1.2145505993  0.6179045983 -0.3009000054  0.7911097440
##  [696]  0.6652118163  1.1865046006  0.8384255678 -0.4267189524  0.4249891669
##  [701]  1.4049414145  0.1303369278  0.9176188968  0.7291836568  0.2065369398
##  [706]  0.3121223882  0.2826794170  1.7351903276 -0.1370589141 -0.0819602459
##  [711] -0.5826915552  0.1375435927  0.0966312635  0.6323989403  1.0810561228
##  [716]  0.5040935202  0.5308324844  0.0204955592  0.7510580428  0.6400562744
##  [721]  0.9679781676  0.4017759385 -0.5746248422  0.5544015504  0.6532529225
##  [726]  0.1641160350 -0.1556828865  0.0457788031  1.2938998841  0.1135663673
##  [731]  0.3205085813 -0.2208336788  0.7446227774  0.2449622153  0.6234431137
##  [736]  0.4588745899  1.2420065630  0.9479428198  0.5244381299  0.7593848800
##  [741]  0.2393961636  0.5990460834 -0.0158569241  0.9750451319 -0.2651665333
##  [746]  1.2761543249  0.3147324107  0.4078529301  0.3742959112  1.1491738779
##  [751]  0.2322208239 -0.2913958173  0.7175825662  0.5841593337 -1.1756410816
##  [756]  0.3596096444  0.4301986387  0.5771426327  0.0472801979  1.7705922146
##  [761]  0.1716950449 -0.4822786871  1.1499130436  0.4586288156  0.3482512648
##  [766]  0.5366545411  0.5441932484 -0.0666798648  1.2411174532  0.5210914175
##  [771]  0.2300706279 -0.5780420683  0.6616236012  1.1248752226  0.8949168106
##  [776]  0.6025351926  0.3401837184  1.4080982769  0.5292568396 -0.5473083399
##  [781]  0.3044076516  0.6022728378  0.1015812581 -0.0959760276  1.1895726976
##  [786] -0.9264898418  0.5819103559 -0.1931293399  0.6580790090  0.7297500954
##  [791]  0.0163394318  0.0157809499 -0.6376387004  0.5498594320 -0.1328214426
##  [796]  0.3767337283  0.3682374844  1.0894271946  0.3993016474  0.1636457958
##  [801]  0.2776132007  1.3981692109  0.5521632865 -0.0511076898 -0.1140011981
##  [806]  0.5941384265  0.7499052566  1.1333267021  0.7928648052  0.7415306622
##  [811] -0.2756763444  0.4376594996  0.5908474192  0.1978941145  0.0086137996
##  [816]  0.3911448445  0.5619436399  0.6047049088  0.8012244847  0.1695497627
##  [821]  0.3609795769  0.2191780188  0.5548371375 -0.0327134789  1.1162610169
##  [826] -0.6914014448  0.5110019450  0.9294747945 -0.3557207428 -0.1533543041
##  [831]  0.3396160126  0.4978088397 -0.0530254268  1.4236938154 -0.2513946084
##  [836]  0.0004688441  0.5545027338  0.1537620218 -0.4020772277 -0.6743822939
##  [841]  0.7770722653  0.9460649460  0.1131163924 -0.2265494032  0.7127249784
##  [846]  1.0468210109  0.2207473980 -0.1875916233  0.8586701561 -0.5103984172
##  [851] -0.0445392391  1.3003135691  0.7778984279  0.3155665915  0.2951478590
##  [856]  0.6462549467 -0.2737245476 -0.5931507370 -0.6530748328  0.5167377008
##  [861]  0.5121185316  0.1064829759  0.3230624572  0.6737717447  0.2999120047
##  [866]  0.8486748194  0.9938704510  0.4175430298  0.0003736385 -0.1988833670
##  [871] -0.5394044162  0.9546684051  0.4069901758  0.0792899034  0.2502697508
##  [876]  0.3672482289  0.7258725442  0.0499187624  1.2219310925 -0.1174694011
##  [881]  0.5377104493  0.3745861074  1.1358944291  0.2546801451  0.0345786316
##  [886]  0.0754029432  0.7404181518  0.3977236477 -0.2215647069 -0.5546781351
##  [891] -0.3287605626  0.0918915232  0.4414657280  0.0222587744  0.6892072578
##  [896] -0.2735160216  0.3392265194  0.1763927239  0.6651658663  0.4701365185
##  [901]  0.3853490524 -0.7863869922  1.5799655196  0.8203315625  0.7874729284
##  [906]  0.6634861461  0.7908050147 -0.3979224909  1.5555613202  0.4046192493
##  [911] -1.3386841361  0.4806400933  0.3127555384 -0.3329676297  0.2641351299
##  [916]  0.2584661242  0.5449774659  0.7745522102  0.2709894565  0.0522713515
##  [921]  0.3477858139  0.2978954482  1.1986262066  0.4553610063  0.9612957774
##  [926] -0.0212211182 -0.0428492486  0.5253165282 -0.5376613668  1.2472793761
##  [931]  0.0327363962  0.4473548894 -0.3646471936  0.7290508864  0.1274054026
##  [936]  0.7540883868  0.6105503630  0.7578823611  0.3708278956  0.3543843969
##  [941]  1.1261276882 -0.0835222939  0.2360945383  0.7335631203  0.3468223603
##  [946]  0.5448350734  0.0863386622  0.2595241269  0.0955394027  0.7682564078
##  [951]  1.0824670449 -0.3902621265  1.0993781693  0.8855868636  0.2168816093
##  [956]  0.1450656591  0.3127023058  0.4246701630  1.0358618157  1.0161306792
##  [961]  0.5799991132  0.5219895076 -0.8386761916  0.5141461560  0.7592021014
##  [966]  0.0056803372  1.3223243200  0.4383220057  1.1719166901  1.0205035213
##  [971]  0.5801722329  2.4223295401  0.4079870683 -0.3126361786  0.7398907358
##  [976]  0.0707439241  1.0095694698  0.1417226462  0.7082775811  0.2520287500
##  [981]  0.7724623750  0.3953167879 -0.1413437150  0.5367178793  0.0616546678
##  [986]  0.1396932467  0.7708244403  0.6852402208  0.3797968028  0.0373503015
##  [991]  1.6908434882  0.9263715803  1.1897386480  0.2820165407 -0.2492739691
##  [996]  0.3671985966 -0.2635451589  0.6352418156  0.4116758795  0.5468424444
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
##   0.75047692   0.29696714 
##  (0.09390925) (0.06640008)
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
## [1] -0.07250958  0.33807002 -0.46156986 -0.41597936  0.55902247  0.04628967
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
## [1] 0.0691
```

```r
sd(results2$thetastar)  #the standard deviation of the theta stars is the SE of the statistic (in this case, the mean)
```

```
## [1] 0.8936519
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
## t1*      4.5 0.03013013   0.9178715
```

One of the main advantages to the 'boot' package over the 'bootstrap' package is the nicer formatting of the output.

Going back to our original code, lets see how we could reproduce all of these numbers:


```r
table(sample(x,size=length(x),replace=T))
```

```
## 
## 1 2 3 4 5 7 9 
## 1 1 1 2 1 1 3
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
## [1] -0.0046
```

```r
se.boot
```

```
## [1] 0.9296468
```

Why do our numbers not agree exactly with those of the boot package? This is because our estimates of bias and standard error are just estimates, and they carry with them their own uncertainties. That is one of the reasons we might bother doing jackknife-after-bootstrap.

The 'boot' package has a LOT of functionality. If we have time, we will come back to some of these more complex functions later in the semester as we cover topics like regression and glm.

