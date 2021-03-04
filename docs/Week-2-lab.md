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
## 1 3 4 5 7 8 9 
## 1 1 2 1 2 2 1
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
## [1] -0.0106
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
## [1] 2.712241
```

```r
UL.boot
```

```
## [1] 6.266559
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
##    [1] 6.0 2.6 5.7 4.2 5.0 4.9 4.8 4.9 5.3 5.4 3.1 4.4 4.4 4.2 5.3 4.9 4.0 4.1
##   [19] 3.0 5.4 4.8 4.6 4.6 4.8 4.9 4.7 4.9 5.3 3.5 4.7 6.2 4.6 4.9 5.1 4.7 5.9
##   [37] 2.4 4.8 5.0 3.7 3.9 3.3 4.4 4.8 4.9 3.9 4.5 5.7 5.4 4.2 5.7 6.0 4.4 5.7
##   [55] 4.3 3.2 4.7 4.2 4.8 6.3 5.4 2.6 3.1 5.1 4.7 5.9 3.6 4.8 4.8 4.2 3.1 3.8
##   [73] 3.3 4.2 4.5 6.3 4.6 4.9 3.9 5.5 4.9 3.7 3.6 4.3 4.0 4.4 3.6 4.1 4.1 4.4
##   [91] 2.1 3.3 5.2 5.1 4.7 4.9 4.0 4.2 2.8 4.6 5.0 5.8 4.8 3.6 3.2 4.0 3.4 4.0
##  [109] 5.9 2.8 3.7 5.1 6.3 3.7 5.4 4.9 3.3 6.5 2.8 3.0 4.4 3.9 3.4 5.1 5.3 4.6
##  [127] 4.2 4.4 5.4 5.6 4.5 4.7 3.9 3.7 4.8 3.7 4.9 3.4 3.3 5.0 3.3 4.8 4.0 3.7
##  [145] 3.8 3.7 4.6 5.6 3.6 4.7 4.3 4.2 6.0 4.7 5.2 3.3 2.5 4.3 2.6 6.3 4.4 5.8
##  [163] 5.6 5.3 6.1 3.8 4.4 3.9 4.3 2.5 4.0 3.3 4.2 4.8 4.4 5.0 5.0 3.6 5.7 4.1
##  [181] 4.6 4.4 4.3 6.7 2.3 5.1 6.7 5.4 5.6 4.4 3.8 5.0 4.4 3.2 4.0 3.4 6.1 4.0
##  [199] 4.0 4.5 4.3 3.0 4.6 4.1 5.8 4.7 3.9 5.0 5.1 3.9 6.8 3.3 4.3 6.1 5.4 5.0
##  [217] 3.8 5.3 4.4 4.1 5.5 4.4 2.5 3.6 4.8 2.3 6.1 4.4 5.9 5.5 4.6 4.0 3.9 4.2
##  [235] 4.8 3.2 3.4 5.7 4.8 4.2 4.1 5.0 5.3 6.5 6.4 4.3 3.4 5.7 5.1 4.8 4.9 3.4
##  [253] 3.8 4.7 4.8 4.1 6.1 4.6 3.1 2.8 5.8 4.3 3.3 5.7 4.7 4.0 2.4 5.2 3.4 5.2
##  [271] 3.8 5.0 2.1 4.5 5.3 2.1 3.2 3.9 5.3 3.5 5.7 5.4 4.7 3.7 3.3 5.1 4.5 4.9
##  [289] 3.9 3.7 4.0 4.8 4.7 3.8 5.2 3.4 5.2 4.7 4.9 4.6 3.5 5.0 5.4 5.3 6.5 5.6
##  [307] 5.0 2.9 3.6 4.5 4.7 6.1 4.3 3.9 2.9 3.9 3.4 2.7 4.9 6.1 4.9 6.5 4.8 3.3
##  [325] 5.1 3.4 4.1 3.8 4.0 2.8 4.3 4.0 5.6 4.1 4.8 4.1 4.4 3.6 2.5 3.3 5.8 3.8
##  [343] 4.3 5.7 4.5 2.6 4.2 6.2 4.3 4.8 4.8 4.7 5.3 4.4 5.1 5.4 3.3 3.5 5.9 4.3
##  [361] 5.4 3.0 5.7 4.8 5.9 4.0 5.2 5.5 3.9 3.5 3.9 3.8 3.4 6.2 5.2 2.0 5.5 2.9
##  [379] 4.9 4.3 4.1 4.9 5.3 4.7 3.9 3.4 5.2 4.7 4.4 4.0 5.1 5.0 3.8 4.1 3.3 5.1
##  [397] 6.0 3.0 4.0 4.3 5.5 4.9 4.4 4.8 4.6 3.9 3.9 4.5 4.3 4.2 4.3 6.0 4.2 3.9
##  [415] 4.1 4.3 2.7 4.0 5.0 5.5 5.6 3.5 4.7 4.8 5.3 6.0 4.6 3.4 4.7 5.6 5.7 3.5
##  [433] 3.0 6.1 4.2 5.1 4.3 5.7 4.6 3.9 4.9 3.9 4.2 4.0 4.6 4.5 4.1 3.1 2.8 3.8
##  [451] 4.1 3.8 4.2 4.7 4.4 5.5 5.3 5.3 5.5 5.1 6.0 4.4 6.4 3.6 5.8 4.4 3.6 2.8
##  [469] 4.0 4.2 4.1 5.2 5.8 4.4 2.7 3.7 6.5 2.9 4.8 5.9 3.4 3.9 5.2 5.1 3.9 5.1
##  [487] 2.6 4.1 5.5 5.2 5.1 4.9 5.0 3.3 5.0 5.3 4.7 4.7 4.8 4.2 3.8 6.4 4.9 2.9
##  [505] 3.4 5.6 5.2 4.4 6.1 4.6 3.7 4.6 4.0 4.2 4.2 3.9 3.9 4.1 3.5 3.9 4.0 3.4
##  [523] 3.7 4.7 5.2 5.3 3.7 6.1 4.2 6.4 4.8 4.4 5.6 5.8 4.4 5.9 4.7 4.6 2.7 4.8
##  [541] 3.6 4.1 4.9 3.1 5.9 3.7 4.2 4.2 4.2 3.7 4.2 4.2 4.8 4.8 4.5 3.4 5.7 5.3
##  [559] 2.7 4.0 1.8 4.3 3.6 5.6 4.6 4.4 5.4 4.9 4.7 5.0 4.7 6.5 3.3 4.7 5.2 4.3
##  [577] 5.9 4.1 3.8 3.6 4.5 6.3 5.1 4.5 4.1 2.1 2.3 2.9 3.5 5.3 3.4 2.9 4.9 4.5
##  [595] 2.9 5.9 6.3 4.7 4.9 5.1 5.2 5.0 3.4 5.9 4.1 3.2 5.7 4.1 4.0 4.9 4.2 5.4
##  [613] 4.3 5.2 4.3 4.4 4.0 6.0 4.2 4.8 3.4 2.8 4.9 4.3 5.2 2.8 3.8 3.3 3.2 6.1
##  [631] 4.1 4.9 5.6 4.1 4.3 4.7 3.1 4.3 3.3 4.9 4.1 5.0 3.9 5.2 6.3 4.4 5.5 3.4
##  [649] 4.3 4.9 4.9 5.5 4.1 4.6 3.8 6.4 4.7 4.0 3.5 4.3 3.4 6.4 2.6 3.8 3.7 5.1
##  [667] 4.3 4.3 5.4 5.0 5.5 6.1 3.3 4.3 4.0 5.8 4.2 3.4 4.0 3.8 6.1 5.0 4.1 4.3
##  [685] 4.3 4.2 3.7 2.6 6.4 5.6 5.9 4.4 5.6 4.2 4.7 3.8 4.9 5.3 3.0 4.5 5.2 3.6
##  [703] 3.6 5.6 3.3 4.9 1.7 4.4 4.1 5.4 5.3 6.0 3.9 4.7 5.7 5.0 5.2 5.4 3.1 3.6
##  [721] 4.5 4.8 4.1 5.1 4.7 3.2 3.6 3.8 3.9 4.2 3.8 4.6 4.6 4.8 3.7 4.1 3.8 3.3
##  [739] 6.1 4.1 4.2 6.5 5.7 3.7 4.2 3.1 4.3 3.4 4.8 3.8 4.4 4.5 3.8 4.3 6.0 3.5
##  [757] 5.4 3.8 4.7 5.4 3.7 5.0 5.8 4.5 5.3 5.0 3.8 3.7 4.3 4.3 3.2 4.5 5.0 4.4
##  [775] 3.4 4.4 4.8 3.9 3.9 5.6 4.8 4.5 5.7 4.5 5.1 4.5 3.9 4.6 4.7 4.5 7.0 3.9
##  [793] 4.4 6.1 3.5 5.2 5.2 4.4 4.5 4.4 3.9 5.6 4.9 3.1 4.5 4.3 3.4 5.5 5.2 6.0
##  [811] 5.1 3.0 4.7 3.2 5.5 5.9 4.4 4.5 4.0 3.9 5.4 4.1 4.5 3.7 5.2 5.1 4.7 4.8
##  [829] 5.4 4.4 5.6 4.0 4.9 4.7 5.6 4.4 5.4 4.2 3.4 4.6 4.6 4.2 5.2 3.9 4.7 5.1
##  [847] 5.8 3.1 4.4 4.0 4.7 3.9 4.3 5.2 3.2 4.7 4.1 5.5 4.3 4.4 6.2 5.4 3.1 3.1
##  [865] 4.2 3.9 3.5 4.6 4.6 3.4 6.4 6.2 5.1 5.7 6.0 2.8 4.8 5.2 3.4 5.1 6.0 5.1
##  [883] 4.2 5.3 4.5 4.9 5.4 5.0 4.0 4.4 5.3 5.4 4.0 6.0 4.1 5.9 3.5 4.1 4.9 3.6
##  [901] 4.9 2.7 2.9 3.8 4.9 4.2 6.4 4.4 3.5 5.1 3.1 4.2 4.9 5.3 5.1 4.5 6.1 5.7
##  [919] 5.0 5.8 4.2 3.6 5.7 2.3 5.3 4.5 4.2 4.8 5.0 4.2 4.7 5.3 5.1 6.0 4.1 5.7
##  [937] 4.0 4.4 4.2 4.2 1.8 4.9 3.1 2.4 3.8 4.7 5.5 4.8 5.9 3.3 4.6 5.5 4.6 3.9
##  [955] 4.9 6.3 3.1 3.3 4.1 3.8 3.7 4.7 5.3 5.8 4.5 4.2 4.5 4.7 3.1 6.2 4.0 5.7
##  [973] 5.1 5.4 3.7 3.8 4.9 4.3 5.1 3.9 4.4 4.0 3.2 3.7 3.5 4.5 4.8 4.3 5.2 3.5
##  [991] 4.9 4.6 5.5 5.0 5.4 4.4 4.5 3.5 5.0 2.8
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
##   2.6   6.3
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
##    [1] 4.0 4.4 6.5 6.1 6.0 4.6 4.7 4.5 2.9 4.3 5.0 4.4 6.7 5.1 3.8 4.7 3.3 3.9
##   [19] 5.7 5.1 4.5 5.0 5.7 5.2 5.5 4.8 5.9 5.7 7.9 4.5 6.2 4.6 4.9 3.8 3.4 4.2
##   [37] 5.3 2.9 5.0 6.2 5.1 4.4 2.7 5.2 4.2 4.3 4.2 3.1 4.7 3.9 5.8 4.3 3.7 4.1
##   [55] 4.2 5.6 4.4 4.4 3.7 4.9 4.1 5.0 5.7 2.1 4.5 5.3 6.0 6.5 4.2 4.2 4.1 5.4
##   [73] 4.6 3.2 3.9 4.6 5.6 4.1 3.6 5.7 4.7 3.7 4.9 4.4 3.5 3.8 4.0 3.9 4.3 4.1
##   [91] 4.2 3.3 3.9 4.3 5.3 3.9 6.3 3.5 3.9 4.5 3.6 4.4 5.1 3.6 5.7 5.6 3.5 5.4
##  [109] 6.1 4.5 3.3 5.6 3.4 4.7 4.6 3.8 4.0 4.4 4.4 3.1 4.1 4.2 4.9 4.1 4.7 3.7
##  [127] 5.1 3.9 4.6 2.7 5.4 4.8 5.0 4.2 5.4 4.9 3.3 3.2 5.3 4.0 3.2 5.7 4.6 2.6
##  [145] 4.9 5.3 5.2 5.5 4.8 5.2 5.5 4.0 6.3 4.0 5.4 4.7 4.3 5.0 5.3 5.4 3.9 4.2
##  [163] 5.5 4.4 4.7 4.1 3.1 5.5 3.8 4.9 6.1 4.2 5.0 4.7 4.0 3.7 6.3 3.8 5.1 3.4
##  [181] 3.6 6.5 5.8 4.8 3.2 4.4 3.8 5.5 3.7 4.3 5.0 4.4 3.4 4.5 2.7 5.5 3.0 4.1
##  [199] 4.2 3.7 4.7 4.3 5.9 5.6 3.4 3.5 4.8 5.8 2.9 4.0 5.9 5.6 3.7 4.3 4.1 4.4
##  [217] 3.5 5.4 5.2 3.4 3.4 4.9 3.0 3.7 5.1 4.7 5.9 5.6 4.3 3.6 3.3 3.5 3.6 3.7
##  [235] 4.5 4.1 3.4 4.7 5.5 4.8 5.5 3.6 5.1 4.4 4.7 4.0 3.2 6.3 4.0 6.5 5.9 5.5
##  [253] 3.7 4.7 3.1 2.6 5.9 4.0 3.3 4.6 6.3 3.2 3.3 3.4 2.7 4.3 5.1 4.7 4.4 4.3
##  [271] 4.0 4.4 4.5 3.7 4.8 5.2 4.5 3.5 3.5 3.8 3.5 5.5 4.2 4.3 4.0 2.6 5.7 3.7
##  [289] 4.5 3.6 3.6 2.7 4.7 5.9 4.2 4.2 5.3 5.7 6.2 5.5 3.8 5.4 4.3 5.2 3.8 4.3
##  [307] 4.1 3.1 5.3 4.1 4.4 4.3 5.0 3.8 4.0 3.4 4.7 4.0 4.2 4.2 4.8 4.9 5.4 5.3
##  [325] 4.9 4.1 3.0 3.7 5.5 4.7 3.9 4.6 4.4 5.1 4.6 5.1 4.3 4.5 5.0 2.8 3.8 4.9
##  [343] 4.3 3.5 4.3 4.5 6.2 4.5 5.7 3.1 5.9 3.9 4.1 3.3 3.9 3.8 5.1 5.5 3.9 4.4
##  [361] 4.1 3.2 3.7 4.8 4.1 4.9 5.9 4.5 4.5 2.9 5.4 3.1 3.1 2.0 4.7 5.1 4.5 4.2
##  [379] 5.3 5.7 3.3 3.1 2.9 3.9 3.2 4.2 4.8 4.1 6.1 4.8 4.9 5.0 4.1 4.0 4.5 3.5
##  [397] 3.4 4.4 4.1 4.8 5.5 4.0 5.1 6.0 4.2 4.4 5.3 4.6 4.5 2.5 3.2 5.0 3.7 6.1
##  [415] 2.8 4.2 4.4 5.6 4.9 2.3 3.5 5.0 2.9 5.8 3.9 3.6 5.4 6.0 3.3 3.3 4.7 3.8
##  [433] 4.3 4.6 5.3 4.8 3.4 4.8 4.0 3.8 3.8 3.2 4.6 3.2 5.4 4.2 5.2 4.7 4.9 4.9
##  [451] 3.1 3.6 5.0 3.2 5.4 3.6 5.6 4.6 3.4 5.3 3.3 4.1 5.0 7.5 4.6 3.7 5.0 4.8
##  [469] 4.0 6.7 5.1 3.6 4.6 5.9 3.8 5.3 4.3 3.9 3.8 5.3 3.3 2.4 4.8 5.7 5.1 5.2
##  [487] 5.9 6.5 4.8 4.9 4.6 5.6 4.1 5.1 5.2 4.8 4.7 5.7 5.0 4.2 3.8 3.8 3.5 4.8
##  [505] 4.8 4.7 4.4 4.5 4.3 6.2 4.1 4.3 4.2 5.4 4.2 3.3 4.0 4.2 4.2 4.7 2.6 4.7
##  [523] 3.8 5.6 4.2 3.0 3.9 5.1 6.7 4.3 5.0 4.0 6.1 4.9 4.7 3.3 3.8 4.7 5.4 5.2
##  [541] 3.8 3.9 2.2 5.6 3.4 4.2 4.4 3.8 4.2 5.7 6.1 4.3 5.5 3.1 4.6 3.2 4.6 3.8
##  [559] 4.3 5.0 4.3 5.0 4.7 4.3 5.8 3.6 5.9 3.8 5.3 5.5 5.3 4.0 4.1 5.5 5.3 5.6
##  [577] 7.0 4.9 5.4 4.7 4.3 4.7 3.8 6.2 5.8 3.0 4.0 3.3 4.0 3.5 4.7 4.1 5.3 6.7
##  [595] 3.7 5.8 3.6 4.0 4.6 4.0 6.1 5.1 5.5 2.8 4.1 2.5 4.8 5.8 4.2 3.0 5.2 3.3
##  [613] 4.8 5.4 5.6 3.4 6.1 5.2 5.4 3.9 4.0 4.4 5.2 4.5 6.9 5.2 6.2 6.1 2.9 3.9
##  [631] 4.2 3.0 5.1 5.1 4.2 4.7 3.3 3.8 4.3 4.1 5.0 4.9 3.5 3.1 7.1 3.9 5.0 4.9
##  [649] 4.6 4.6 5.0 4.4 5.0 5.0 3.2 3.7 5.4 4.1 4.8 4.9 5.3 6.5 4.0 3.6 5.7 4.5
##  [667] 2.9 4.1 4.6 6.0 3.9 5.7 6.0 4.7 4.7 5.2 5.0 4.5 6.7 4.6 4.4 3.9 5.9 3.5
##  [685] 3.3 4.5 3.3 4.5 4.6 4.5 5.5 3.2 5.0 4.0 4.1 5.1 4.4 3.2 4.2 4.6 3.5 3.4
##  [703] 4.5 5.9 4.7 4.8 3.9 4.4 5.0 2.5 3.6 5.8 5.5 5.1 4.3 4.8 4.9 4.3 5.2 5.3
##  [721] 6.2 6.1 3.9 3.3 5.2 6.6 4.3 4.0 4.6 4.5 3.6 3.9 3.9 5.3 4.0 7.0 3.5 5.5
##  [739] 4.2 3.8 3.8 4.1 5.0 4.8 2.7 5.8 4.6 4.3 3.6 3.7 5.9 3.9 3.7 4.4 2.4 6.6
##  [757] 4.3 4.9 4.3 5.0 4.8 5.0 6.0 5.6 4.6 5.0 3.4 5.7 3.6 4.3 3.9 4.8 4.3 5.8
##  [775] 4.9 4.4 5.1 3.8 3.3 3.7 5.2 4.2 3.9 5.1 5.2 5.3 3.4 3.4 4.6 4.1 4.4 3.8
##  [793] 3.2 3.7 4.9 4.8 5.0 5.0 3.9 4.3 3.6 5.3 5.5 4.4 3.8 4.3 5.1 3.7 5.0 4.7
##  [811] 5.1 4.8 3.9 3.7 5.8 6.4 4.0 4.7 4.4 3.7 4.8 3.6 4.5 3.1 4.7 3.5 4.8 2.8
##  [829] 4.3 4.3 3.0 5.2 4.3 3.2 4.3 4.5 4.5 4.7 6.0 4.5 6.5 4.2 5.5 4.0 5.2 4.9
##  [847] 3.7 5.6 5.5 5.5 5.2 3.6 4.1 6.3 3.7 3.6 2.9 5.1 4.6 3.3 4.6 3.2 5.3 4.7
##  [865] 5.0 4.4 4.7 2.7 5.1 2.8 4.0 4.2 6.0 4.3 4.5 2.9 4.7 5.5 4.0 3.1 4.6 4.5
##  [883] 3.7 4.2 4.9 3.8 3.9 6.9 3.8 5.0 3.4 4.4 4.5 4.1 4.6 6.0 6.5 5.9 4.6 4.6
##  [901] 5.5 6.7 3.7 3.8 4.5 5.7 4.5 4.0 4.2 4.5 4.9 3.8 3.2 4.6 4.1 4.7 3.4 3.7
##  [919] 4.9 3.6 5.0 4.2 3.2 5.2 4.8 4.5 5.7 4.9 5.4 5.1 2.7 4.3 4.6 3.4 4.4 3.9
##  [937] 4.7 3.9 5.7 4.9 4.4 4.4 5.0 4.6 6.7 2.9 3.6 4.1 3.8 3.5 2.8 4.4 5.9 5.8
##  [955] 5.0 3.2 4.0 5.1 4.2 5.0 5.8 3.7 6.2 3.7 4.0 4.1 4.3 5.3 4.3 3.8 5.7 3.9
##  [973] 4.0 4.3 4.8 4.9 5.5 5.4 4.0 5.9 3.1 4.8 4.8 4.6 5.1 5.3 5.4 3.4 3.5 5.9
##  [991] 5.1 5.4 4.6 2.8 4.4 3.6 4.1 4.8 4.1 3.6
## 
## $func.thetastar
## [1] -0.001
## 
## $jack.boot.val
##  [1]  0.51485714  0.39029412  0.33293769  0.24689266  0.01578947 -0.03111702
##  [7] -0.18932927 -0.21095890 -0.47994100 -0.50758427
## 
## $jack.boot.se
## [1] 1.020923
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
##    [1] 4.5 3.2 3.9 2.9 4.9 5.2 4.9 4.0 6.4 5.0 3.6 6.8 5.1 3.5 4.3 4.8 4.7 3.7
##   [19] 5.3 3.8 4.5 4.5 4.9 4.5 3.7 4.0 5.3 5.0 4.3 3.3 4.7 5.0 4.3 3.4 5.7 4.8
##   [37] 2.5 3.8 4.7 5.9 5.2 5.5 4.4 4.9 3.2 4.5 3.1 6.8 4.8 4.3 6.0 5.1 4.2 3.9
##   [55] 5.7 5.5 5.0 4.6 5.7 5.4 5.9 3.8 5.4 5.0 4.5 4.1 4.3 3.8 4.6 4.3 4.3 5.7
##   [73] 5.1 4.7 4.3 3.8 4.3 5.0 5.6 4.2 3.6 5.6 6.1 3.9 5.1 3.9 6.3 4.4 5.3 6.0
##   [91] 5.7 3.6 3.2 4.7 5.0 3.4 5.1 3.8 4.9 4.3 5.3 4.6 3.2 5.6 4.0 5.7 5.3 4.0
##  [109] 4.8 4.5 3.1 5.2 3.5 4.7 5.7 3.4 3.9 4.7 3.2 4.0 3.8 4.3 5.2 4.1 5.5 4.9
##  [127] 4.0 4.5 5.0 4.7 4.0 5.4 3.4 5.4 4.8 5.1 5.1 3.7 5.2 3.6 4.3 3.5 4.0 4.4
##  [145] 6.2 5.6 6.5 7.2 4.8 5.3 5.6 4.5 4.9 3.9 5.0 3.0 5.3 6.5 2.4 3.0 4.5 4.6
##  [163] 4.7 3.5 4.0 4.9 5.3 3.9 5.5 4.1 3.4 4.0 4.5 4.8 3.9 3.6 3.4 4.4 6.3 3.6
##  [181] 3.8 3.6 4.1 6.8 5.5 4.3 4.1 4.3 3.5 4.7 3.8 3.7 5.6 3.3 3.0 3.3 5.9 3.2
##  [199] 4.2 4.4 3.9 4.3 4.7 4.3 6.0 5.1 2.2 5.0 3.2 4.1 3.0 4.1 3.3 2.4 5.4 4.1
##  [217] 4.1 4.7 4.9 6.0 3.5 4.5 4.1 5.6 4.5 4.6 4.0 5.5 6.5 5.0 6.4 5.8 4.2 3.0
##  [235] 5.9 3.4 3.5 4.0 2.9 2.7 6.3 6.8 2.9 4.1 2.6 3.7 5.3 4.3 4.0 3.9 4.7 3.8
##  [253] 2.4 4.7 3.7 5.0 6.4 2.7 6.4 5.4 6.9 5.9 4.9 4.5 4.8 4.4 4.4 4.7 3.2 3.7
##  [271] 5.6 5.5 5.9 5.5 3.5 4.9 5.7 3.4 4.1 4.6 5.5 4.5 6.1 4.2 3.0 3.6 3.1 3.3
##  [289] 4.6 3.0 3.2 4.9 3.4 4.9 5.7 4.2 4.1 4.8 4.8 4.8 3.0 5.8 3.3 5.8 4.0 3.7
##  [307] 5.6 5.7 5.5 6.0 2.8 5.4 5.8 3.3 3.2 5.1 3.5 5.4 4.6 4.4 4.1 4.1 3.3 3.4
##  [325] 4.3 6.5 2.8 4.3 5.3 6.0 4.6 5.6 4.0 4.8 3.4 4.1 5.6 3.6 4.7 4.7 4.5 3.3
##  [343] 4.0 2.8 5.8 3.1 5.2 4.7 5.0 4.5 4.1 5.4 4.8 3.7 5.0 3.9 4.6 4.8 4.6 4.5
##  [361] 5.7 3.7 2.9 2.8 4.9 2.9 5.4 4.2 6.8 3.5 4.2 5.4 2.6 5.3 6.4 3.8 4.6 6.8
##  [379] 5.4 5.4 4.3 5.8 4.1 4.3 7.2 4.7 3.3 3.8 5.3 3.5 3.0 5.6 4.9 4.4 4.3 5.4
##  [397] 5.4 4.6 4.9 4.9 4.9 5.5 4.5 4.3 5.3 4.0 5.0 4.3 4.3 4.1 3.8 4.0 4.8 4.2
##  [415] 6.5 2.7 2.9 4.4 4.1 4.8 4.1 4.2 5.0 4.2 4.0 5.2 4.7 5.9 3.4 2.3 5.8 5.2
##  [433] 5.5 5.6 4.5 2.8 4.0 4.4 4.2 5.7 4.4 4.9 3.6 6.2 4.9 4.0 3.6 4.8 4.2 5.2
##  [451] 5.2 5.5 4.4 3.5 4.0 2.8 5.3 4.8 4.0 3.7 5.5 5.8 3.1 5.1 4.6 3.4 4.5 5.2
##  [469] 5.3 5.4 4.7 5.0 4.5 3.5 4.6 4.6 3.3 5.0 3.5 4.4 4.5 4.1 3.5 4.3 3.4 4.9
##  [487] 5.6 4.9 2.7 6.5 5.5 3.9 4.5 2.1 5.8 4.6 5.4 4.2 4.6 5.3 4.6 4.6 2.7 4.4
##  [505] 3.8 2.9 3.4 6.2 5.5 3.1 5.3 4.2 4.4 3.1 4.5 4.4 4.2 6.1 4.6 5.1 4.4 4.0
##  [523] 3.4 3.6 3.3 5.7 6.7 2.7 4.4 3.8 3.9 5.9 3.7 4.8 5.0 2.9 2.1 3.8 5.2 5.0
##  [541] 4.5 4.2 3.3 4.4 4.0 3.5 5.7 3.7 5.2 4.1 4.2 4.7 4.1 6.7 3.6 4.2 4.7 3.7
##  [559] 3.7 4.4 2.6 1.9 4.8 4.8 5.7 2.7 3.3 4.2 4.2 5.0 6.2 3.9 5.2 5.1 4.1 3.6
##  [577] 3.4 4.6 3.8 4.4 4.0 2.3 5.3 5.4 3.1 4.1 4.2 4.1 4.7 4.0 4.5 4.7 5.1 4.1
##  [595] 4.8 3.7 3.5 5.1 4.9 4.8 3.4 3.9 4.6 3.6 4.4 5.4 5.0 5.3 3.6 4.0 4.1 2.9
##  [613] 5.8 5.7 4.4 4.0 5.4 4.4 6.0 5.1 5.3 4.2 4.1 5.2 3.3 5.0 3.3 3.6 3.5 5.8
##  [631] 4.7 3.6 4.5 3.5 4.8 4.8 4.3 2.9 5.6 3.8 4.5 4.9 4.7 5.5 5.2 5.4 4.4 5.5
##  [649] 3.5 3.9 5.3 5.3 5.3 5.5 5.6 5.8 4.8 5.5 4.3 5.7 4.8 5.2 5.2 5.4 5.0 4.2
##  [667] 4.3 2.8 3.8 3.7 5.6 5.9 2.6 4.0 2.2 3.6 4.4 4.4 4.6 5.0 4.3 3.9 3.5 4.6
##  [685] 3.5 4.5 4.6 3.1 4.2 4.1 3.2 5.0 3.6 4.1 4.1 4.7 4.9 2.6 3.8 4.8 5.8 3.9
##  [703] 2.9 4.6 4.4 3.8 3.5 4.5 3.4 4.2 5.2 5.0 5.4 5.2 3.5 3.6 3.5 5.4 3.9 6.3
##  [721] 3.4 5.2 5.8 4.1 3.4 3.7 5.5 4.3 3.3 3.9 4.3 6.3 4.2 4.7 2.9 5.7 5.3 6.1
##  [739] 3.8 6.1 5.6 5.7 5.0 3.0 2.7 5.5 5.4 5.0 4.5 5.2 5.9 4.4 3.2 6.2 4.3 4.3
##  [757] 5.1 4.0 6.1 4.3 4.6 3.8 4.1 4.8 3.4 3.9 6.4 4.4 5.3 3.7 4.9 3.7 4.4 5.1
##  [775] 3.8 5.3 2.5 2.9 3.8 4.1 4.8 4.2 5.1 4.8 5.6 4.7 4.9 4.9 4.3 4.3 3.3 4.3
##  [793] 6.1 4.3 5.2 4.9 5.5 4.5 5.8 6.4 3.9 3.7 5.8 4.1 6.1 4.1 5.1 4.8 5.1 5.4
##  [811] 5.7 5.3 4.9 5.7 4.3 5.1 4.0 5.6 3.4 4.4 5.5 2.5 5.1 3.8 2.6 4.1 4.1 2.9
##  [829] 5.4 5.1 3.8 4.9 4.9 4.3 4.0 4.6 5.9 5.2 4.1 5.9 3.1 4.6 4.1 3.1 4.3 5.0
##  [847] 5.1 2.7 4.6 4.5 3.8 4.8 4.3 4.6 3.6 5.1 3.2 5.1 4.7 4.7 5.2 5.5 3.8 4.5
##  [865] 4.1 4.3 4.9 4.9 3.8 3.9 4.8 5.6 4.4 4.7 5.1 4.3 3.4 4.1 5.9 3.8 3.8 3.7
##  [883] 4.2 4.4 4.3 4.7 4.2 4.9 4.7 5.3 5.8 4.6 5.3 5.5 5.8 5.4 5.3 4.3 6.1 4.8
##  [901] 3.4 5.6 6.0 4.1 5.1 4.5 6.1 3.3 4.6 5.5 5.1 4.0 3.2 4.6 5.0 5.8 2.7 5.1
##  [919] 5.5 4.4 5.3 2.9 3.6 4.9 4.0 5.1 4.8 5.7 5.1 3.5 5.1 3.8 4.6 4.6 4.2 4.2
##  [937] 3.1 4.6 4.9 4.0 2.8 5.0 2.5 5.1 4.5 5.4 4.1 5.0 4.9 4.5 3.6 4.7 5.3 4.3
##  [955] 3.6 5.7 4.8 5.9 6.0 5.0 4.4 5.1 6.0 3.0 3.8 3.4 4.3 4.4 5.8 6.0 4.4 5.4
##  [973] 3.2 4.8 3.7 4.4 2.4 4.5 4.1 3.5 4.1 2.4 2.6 3.8 6.0 3.9 4.2 3.7 4.8 4.3
##  [991] 4.7 4.1 3.1 3.8 4.3 4.7 5.1 3.9 4.8 3.5
## 
## $func.thetastar
## 72% 
## 5.1 
## 
## $jack.boot.val
##  [1] 5.5 5.5 5.4 5.2 5.1 5.0 4.9 4.8 4.7 4.5
## 
## $jack.boot.se
## [1] 0.9785704
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
## [1] 0.4438572
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
##      shape       rate   
##   2.1663976   4.4306117 
##  (0.9041979) (2.0798433)
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
## [1] 0.439434 1.388359 1.187570 1.627672 0.102301 1.249693
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
##    [1]  0.3297009614 -0.2770179120  0.5368274785  0.5911754929  1.0631512750
##    [6] -0.0521354130  0.1312919335  0.2880409816  0.3124277737 -0.1931118701
##   [11]  1.0447822727  0.4459597107  1.1092806157  0.6343296589  0.2295108138
##   [16]  0.1086998274  0.1594328068  0.2897658634  0.4336579828  0.3427625708
##   [21]  0.1415227035  0.4399970013  0.2470925966 -0.3135793645  0.2230482724
##   [26]  0.0820372625 -0.4373596960  0.0316086067  0.3535888039  0.2602769220
##   [31]  0.1555560751 -0.0290201677  0.1843835782 -0.8527589308 -0.0397130860
##   [36] -0.1064692027 -0.1601500076  0.6674038362 -0.1856038610  0.2491268799
##   [41]  0.7900001493  0.7239009829  0.6111509868  0.2126339518  0.2974726440
##   [46]  0.5427120592  0.2404998613  1.6894653461  0.0982098679  0.4887472283
##   [51]  0.9362019152  0.8824673682  0.5285056914  1.2587591582 -0.0890300247
##   [56] -0.4026773050  0.8337428367  1.5623809213  0.5292880025 -0.3235503156
##   [61]  1.6004463426  0.8667317644  0.5882750507 -1.0008743950  0.7159359162
##   [66] -0.4248718151  1.0833171528 -0.0699534935 -0.0939353380  0.9061084095
##   [71] -0.1345573335  0.7430054596  0.2987547059 -0.4371434810  0.2056226199
##   [76]  1.0291339955  0.4835987464  0.6187076981  0.4556819784  1.1256223929
##   [81]  0.3463513839  0.0921814685 -0.4952386152  0.4821049099  0.4399785975
##   [86]  0.4379220111  0.5926035963  0.4893662165  0.3133089371 -0.1171597104
##   [91]  0.2900958123  0.0786692806  0.8925041149  0.1362368510  0.4963985584
##   [96]  0.2951729484  0.0674985183  0.0023007682  0.8018385587 -0.0623077494
##  [101]  0.1045846804  0.4446429564  0.7604364711  0.4523103445  0.7722164732
##  [106]  0.2575951175  0.2199978805  0.0715718614  0.4202049102  0.4219635064
##  [111]  0.2245606501  0.1528753925  0.4624345905  1.0137922361  0.5210797706
##  [116] -0.6454240180 -0.3029200373 -0.2823202748 -0.0774653560 -0.2360996457
##  [121]  0.5325716526 -0.1270310884  0.5984018040  0.8793449242  0.7152442204
##  [126]  0.3099998873 -0.6412373590  0.2501049113  0.6662722763  1.0637225552
##  [131]  0.2737695062  0.0512465488  0.8675333121  0.6041134061 -0.0224258525
##  [136] -0.4938165069  1.1604622156  1.0433756280  0.2145566051  0.1128172879
##  [141]  0.3166699783  0.0709910476  0.8084810713  0.1480285259  0.2266059512
##  [146]  0.5054381465  0.6142526923  0.9669383644 -0.0336595143  0.5655941803
##  [151] -0.2822601934  0.0146375394  1.2509400152  0.1522398648  1.0131267846
##  [156]  0.5888998475 -0.0454162849  0.4335238327  0.4086148196  0.2306028931
##  [161]  0.5618532130  0.1752126423 -0.0248746848  0.3571707936  0.2751994052
##  [166]  0.3904305508 -0.0460899703  0.6314703539  1.0784110153  0.8129042855
##  [171]  1.0050574757  0.4401765975 -0.0135048295  0.3344581198 -0.0971035044
##  [176]  0.4309347127  0.5821396094  0.7399655535  0.6431532451  0.2283727100
##  [181]  0.5040516405  1.0067555849  1.0638512496  0.4718410929  0.1000483627
##  [186]  0.6539839425  0.6545050786  0.6301731715  0.0787115356 -0.1242380108
##  [191] -0.3781083477  0.1191835720 -0.7053876893  0.4825639540  0.2598968392
##  [196] -0.4525074841  0.9384641090  0.5304986971  0.1004394061  0.0824046565
##  [201] -0.5503353937  0.4549150988  0.3115042111  0.0817726189  0.7640083677
##  [206]  0.6320632034  0.9956363589  0.7852713164  0.5096531369  0.6493035751
##  [211]  0.2694496383  1.5476740257  0.2918600142  0.0594866268  0.1795873582
##  [216]  0.1336369630  0.1242702530  1.0247789846  1.0849162535  0.0930296130
##  [221]  0.3618652207 -0.0880349098  0.1956394888  0.3408814360  0.4705035974
##  [226]  0.4601367536  0.5288982416  0.6140230665 -0.2270077379 -0.2549400865
##  [231]  0.4543335503  0.5310620444 -0.0362736435  0.3020064231 -0.3418901683
##  [236]  0.7401972081  0.3690193293  0.1762634867  0.2354713761  0.4386941158
##  [241] -0.0633277693 -0.0431424235  0.7630926082  0.4062969067  0.7390908842
##  [246]  0.7190127030  0.5946684504  1.2246292213  0.0372323553  0.1097181000
##  [251]  0.7585786126  0.9377021604  0.6030811781  0.4894326757  0.3889617944
##  [256] -0.1784021355  0.7418879375 -0.2853032936 -0.3783388055 -0.0292166584
##  [261]  0.2928654819  0.3660918342  0.5548691367  1.3566468491  0.7515553519
##  [266]  0.0109557457  0.8199982687  0.9453898927  0.0451555224  1.4973049075
##  [271] -0.5367192153  0.2061075435 -0.1333940385  1.0073622066  0.0084226563
##  [276]  0.0580987385  0.1106404390  0.6690190920  1.0728002982  0.5999132005
##  [281]  0.0666842880  0.1181058014 -0.0943470220  0.2893688155  0.1036834141
##  [286]  0.5481740155  1.2412789477 -0.1217673922  1.0133303767  0.3944746746
##  [291]  1.1481638929  0.9700155024  0.4282612497  0.1946185093 -0.1279540592
##  [296]  0.7667025428 -0.0838076081  1.1635081458  0.6535329948 -0.5182760290
##  [301]  0.3971664830 -0.0737766774  0.3389510739  0.5999766440  0.5784641396
##  [306]  0.4365343787 -0.0920721304  0.5284274588  0.8792559831  0.3603590108
##  [311]  1.1382943652  1.4387600856 -0.0404076779 -0.2572341881  0.2130756082
##  [316]  1.5204470245 -0.6294857976  0.4440531436  0.4354546440  0.2784365326
##  [321]  1.1855246656  0.2001420484 -0.1539535075  0.0693208213  1.2701752336
##  [326]  0.9986858618  0.6766879068  0.2194854775  0.0625324902  0.7033350090
##  [331] -0.2408016641  0.8594279214 -0.4310492087  0.3803832310  0.7951315540
##  [336]  0.2966822451 -0.1181476093  0.3407396221  0.3506039522 -0.0368013530
##  [341]  0.4514740041  0.4840164130  0.1255408341  0.0248626039  0.4361760707
##  [346]  0.7086356066  0.5226523685  0.7304208510  0.5108497554  0.4717880187
##  [351]  0.3517821290  1.1013810777  0.8368766301  0.4198836317  0.5580655639
##  [356]  0.3556205682 -0.4098456167  0.0071049417 -0.0374040845  0.7791269341
##  [361]  0.1621819845  0.5101651098  0.0926777072  0.7767854444 -0.0488243685
##  [366]  0.7900620414  0.9334452885  0.7844437815  0.5153239446  0.2715291815
##  [371] -0.2013607633  1.4417393051 -0.3453381537  0.8332966210  0.5429589096
##  [376]  0.3506125284  0.9972395017  0.7390898385  0.8300631007  0.0624465965
##  [381]  0.3811447256  0.2486002066  1.1845702165  0.7075342880 -0.0948598970
##  [386] -0.1177674702  0.1773560201  1.0561763835 -0.1738902482  0.0720649455
##  [391]  0.2923769069 -0.2772588481  0.8026617929  0.2458594618  0.5244913449
##  [396]  0.3487904292 -0.1970037689  0.0540398727  0.7144698385  0.2116014324
##  [401]  0.6607665189 -0.0074642595  0.2940400149  0.4076533306 -0.1056786627
##  [406]  0.3002097315 -0.1088573009  0.4110214758  0.6610227851 -0.3296938755
##  [411]  0.2270358466  0.2864941814 -0.2929764426 -0.3073037286  0.8899578356
##  [416]  0.4178278473 -0.4493997495  0.8667308187  0.5281544815  0.2506577675
##  [421] -0.2189737963  0.3187512738 -0.4653734438  0.3990763155  0.6729032661
##  [426]  1.0115003429  1.5916062343 -0.1073969996  0.0194351784  0.3853699456
##  [431] -0.0625986894  0.3195821056  0.4053780990  1.1037958705 -0.1031787771
##  [436]  0.3639497175  0.4942501492  0.3029981593  0.5267806994  0.6321440781
##  [441]  0.2230482724  0.2292490048  0.3816827050 -0.7683779581  0.6549759100
##  [446]  0.8067189120 -0.1609439030  0.7490126693  0.1583202739  0.5887587628
##  [451]  0.7374942136  0.3870258636  0.6324073584  0.8062881657 -0.1915033370
##  [456] -0.4251035711  0.7114080047  0.3971664830  0.5511985816  0.8076696393
##  [461]  0.4793678261  0.5710310176 -0.0070644012 -0.3069686815  1.5377517237
##  [466]  0.0884262720 -0.1755457258  0.3593400496 -0.3633956331  0.5768050378
##  [471]  0.4234517713  0.4199419993  0.9132645953  0.8663454060  0.2018455925
##  [476]  0.3767577006  0.0745917854 -0.0876238461  0.3433224834  0.6411731196
##  [481] -0.1837263270  0.8005963465 -0.1570013287  0.5525786059  0.8771417499
##  [486]  0.6352707656  0.2686666615  0.3547815393  0.5686477237  0.4216606284
##  [491]  0.2992397528  0.3710208292  0.6841402615  1.2900933674  0.9556478344
##  [496]  0.0839962982  0.7612396477  0.7827376650  0.4016545935  0.9412153623
##  [501]  0.9395745286  0.1919910410 -0.3025405319  1.0524354373  1.7691668608
##  [506]  1.0865662057  0.2855077045  0.2466959572  0.3338303073  0.8777374935
##  [511]  0.7614904612  0.7468819405  0.4573662236  0.4323834009  1.3531736223
##  [516]  0.7618459714  0.0351261872  0.3776200753  0.1155475232  0.3735599335
##  [521] -0.0009278322  0.5367723533  0.2856119957 -0.2959080701  1.2536758451
##  [526]  0.0114591511  0.1381101524  0.8831695644  1.0358167890  0.5861090319
##  [531]  0.2977661751  0.4879019035  1.1003027528  0.1745562593  0.2330699031
##  [536]  0.6165698970 -0.2238939045  0.1597097206  0.6013722136  0.1620765516
##  [541]  0.0435605135 -0.2623142124 -0.0319053692  0.3851890596  0.5205950900
##  [546]  0.3881266809  0.6729530299 -0.1384618956 -0.3148906235  0.2891323905
##  [551]  0.2451959628  0.5675753463  0.0849184667 -0.2851628661 -0.2692968699
##  [556]  0.6725055303  0.0438272838 -0.1511078926 -0.1308228098  0.1364840381
##  [561] -0.3258165160  0.9003689272  0.1762542812  0.7657614264  0.8504832548
##  [566]  0.3729048359  0.5427908181  0.7471399011  0.3290063651  1.2167074098
##  [571] -0.7179599631  1.2076622124 -0.3188265098  0.5240444374 -0.1984675772
##  [576]  1.4077343049  0.1971871849 -0.5734974815  0.3990836718  0.1245060670
##  [581]  0.1876071868  0.6673641142 -0.3554176789  0.0929054834  0.3315859719
##  [586] -0.4497067497  0.8191189942  0.4233125956 -0.3877238650 -0.3192397804
##  [591]  0.5411716277  1.3976251688  0.9260192215  0.3162137777  0.4729076178
##  [596]  0.3531184537  0.0003767408  0.0444527635  0.0825486703  0.5092115493
##  [601]  0.0366843127 -0.3181838941  1.0132867788  0.4144992055 -0.2197894849
##  [606]  0.7922468546  0.5155235052  0.6736943031  0.8566886949 -0.3611997767
##  [611]  0.7519069671 -1.0619118623  0.9075852379 -0.2478115857  0.5838687775
##  [616]  0.3712827477 -0.0978442968  0.2638155664  0.4375051045 -0.7441775919
##  [621]  0.4875282429 -0.0072996023  0.0073279829  0.2350747527 -0.3909920198
##  [626] -0.1573952513  0.9803294896  0.1157960108  0.5382696806  1.3754514401
##  [631] -0.0486450961  0.2594032084 -0.0837567168 -0.4223678152  1.1143948348
##  [636]  0.9389562629  0.0186096126  0.0838361374 -0.0472108099 -0.1732034366
##  [641]  0.2060703030  0.5860638591  0.6068011146  0.2779069232  1.4085267548
##  [646]  0.0554335299  0.3759423971 -0.1202111009 -0.2021200307  0.0902965176
##  [651]  0.1895490616  0.0048500756  0.5756478085  0.3881263254  0.3361411785
##  [656]  0.5103094173  0.3063242935  0.2760581630 -0.2025622700  0.6148565568
##  [661]  0.2498934558 -0.0678426987 -0.2608213528  0.1706786317 -0.0939353380
##  [666]  0.0168329357 -0.0863550070  0.8976926687  0.8872091189  1.6807909600
##  [671]  0.0548980291  0.1514114744 -0.0139545961  0.6136350092  1.0358839349
##  [676] -0.2875633943  0.4680952017 -0.2851542758 -0.0806629802  0.5789851384
##  [681]  0.5103094173  0.6316741984  0.8213547489  0.3330979163  0.6030811781
##  [686] -0.2501697389 -0.0884187236 -0.0701486107 -0.5238950230 -0.3520660806
##  [691] -0.0580381810  0.0244249961  0.2347479151  0.0873056268  0.5130477687
##  [696] -0.0021369231  0.6649332688  0.5553838121 -0.2062683598  1.3786641161
##  [701]  0.1715532057 -0.1291269859  0.1043773073  0.2432200568  0.6727726712
##  [706]  0.4930749389  0.8875757912  0.9257829233  1.2787154845  0.6354943459
##  [711] -0.9347961205  0.1605056371  0.1124844469  1.1391446760  0.2477213775
##  [716]  0.3155336585  1.1880262276  0.8172585467 -0.1454056605  0.3365133472
##  [721]  0.3852810493  0.6511019574  0.3884863785 -0.0689733294  0.2295840169
##  [726] -0.2241862380 -0.2616365461  1.1560720882 -0.0165241297  0.4745327685
##  [731]  0.2165686948  0.9298757700  0.5370064804 -0.2401680516  0.5561971947
##  [736] -0.9896995099  0.2713767284  0.8922636314  0.8912544575  0.0653718345
##  [741]  0.8243417021  0.6400739879  1.0653444263  0.1376172563  0.0837780470
##  [746] -0.3396330026 -1.0340961315 -0.0538806640  0.4191477980  0.4654958087
##  [751] -0.0073457734  0.1189691015 -0.2264458256  0.7293903213 -0.0724302389
##  [756] -0.1330669869 -0.2365904329  1.5053730823  1.0141072595  0.4814925884
##  [761] -0.1046422647  0.0142554696  0.5094008747  0.5392282592  0.4015243543
##  [766]  0.4438086343  0.4489159218  0.2377068697  0.4258178338  0.1799187955
##  [771]  0.2818201481  0.7954592005  0.6966469778 -0.0601142106  0.5698760030
##  [776]  0.1031794468  1.3523847840  0.2157061823 -0.1931859567  0.6939963575
##  [781]  0.0119466213 -0.0785264577  0.3541399624  0.2176818481  0.2249209821
##  [786]  1.4541168724  0.0327966759  0.6169093028  0.3689078636 -0.1855670703
##  [791]  0.1400107339  0.3998117202 -0.1575953856  0.4203885683  0.2234979948
##  [796]  0.7380609989  0.5334849954 -0.1350175348 -0.0267795338  0.3872068432
##  [801]  0.6267727027  0.4694041549  0.4000690604  0.0202227892  0.2776364809
##  [806]  0.2754118970  0.3300306286  0.7409115569  1.0406734519  1.4390069713
##  [811]  1.2265289616  0.1260591440  0.1396353833  1.0755648192  0.3935572780
##  [816] -0.1692676029  0.0770196276  0.4540464111  0.7966660822  0.4486251896
##  [821] -0.6238905939 -0.0180236382  0.6318775261  0.5714786562  0.8409548625
##  [826]  0.7515330448  0.6521473574  0.1985672842 -0.2880968126 -0.1971516254
##  [831]  1.1951379331  0.1276772260  0.4844872471  0.2671503950 -0.4231899115
##  [836]  0.6501881774  0.0529068881  0.4181244143  0.5088618798 -0.5911300478
##  [841] -0.1827857122  0.1453584643  0.0069074027  0.5737318882  0.0627460746
##  [846]  0.5124559732  0.2826480181  1.2512924012  0.1821354979  0.3958682681
##  [851]  0.8104279638  0.3471330185 -0.0562605634 -0.2324638110 -0.0024666436
##  [856] -0.3313066435  0.9154839290 -0.7522841843  0.1908531619  1.2721967099
##  [861]  0.1717597700  0.2369898627  0.0107717748  0.0264545803  1.1497063753
##  [866] -0.0562320643  0.2181379696  1.2375775925  0.7508712304  0.6131856508
##  [871]  0.6521976075  0.0140148654  0.4351660715  0.2798115489  0.2370400840
##  [876]  0.0703010175  0.6760865576  0.3349645898  0.5275450568 -0.2126794110
##  [881]  1.0132269330 -0.0431202778  0.0246510696  0.2405834324  0.5644606507
##  [886]  0.5959508948  0.5052133748  0.2714722258  0.1704812883 -0.1079796738
##  [891]  0.7485285505  0.3540166794 -0.1377412294  0.6579648163  0.6410706165
##  [896]  0.1449893109  0.7717800884  0.4821049099  0.5607611619  0.8575133154
##  [901]  0.1968538198  0.4236527296 -0.1302786687  0.3544692434  1.0021407164
##  [906]  0.4425436836  0.2180180400  0.3476280777  0.1661745603  1.4528526684
##  [911] -0.3929807014  0.1956482535  0.2821360506  0.4011741286  1.0527487087
##  [916]  0.7404016967 -0.7386298994  0.4708820932  0.3717692946 -0.0745723361
##  [921]  0.5573498238  0.0506545143  0.6839339903  0.1123610905 -0.2496895682
##  [926]  0.9805090464  1.4381425049  0.9717962900  1.0017810764 -0.0459254643
##  [931] -0.5351473549  0.2692542714  0.0463389463 -0.0819545284 -0.0475636284
##  [936]  0.4665461240 -0.0326527068  0.8779393794  0.5870374945  1.2167074098
##  [941] -0.2367236852  0.6883697650  0.7205036751 -0.0731809926 -0.0494489863
##  [946]  0.3701073945  0.0685099997  0.8222576749  1.2633773467  0.1981514785
##  [951] -0.0767238876 -0.1744958624  0.1411255929  0.2431812870 -0.3920867436
##  [956]  0.1812242756 -0.1127827313  0.2428225376 -0.2421821835 -0.0593970257
##  [961]  0.4941358757 -0.1578936489  0.9782058910  0.4464484572  1.0811361117
##  [966]  0.4546375962  0.0880920752  0.2656335331  0.8159946696  0.8619294051
##  [971] -0.0555876575  0.7844097342  0.7363359320  0.1143652729  0.9915318507
##  [976]  0.8793449242  0.3656031708  0.7270822628  1.0676989871 -0.0658003831
##  [981] -0.3051838390 -0.2781516149  0.2922720149  1.1518633600  0.6805993404
##  [986]  0.0810679495  0.9084524956  0.2563575599  0.5110656487  0.3529621828
##  [991] -0.1226776873  0.4680886969  0.5285209788 -0.2996891596  1.1079982660
##  [996] -0.4387927726  1.4333397588  0.4699832074 -0.2801583835  0.5555834826
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
##   0.48896149   0.30423764 
##  (0.09620839) (0.06802725)
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
## [1]  0.5159975 -0.9457983 -0.8982600  0.4452285 -0.3312530  0.5204171
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
## [1] 0.0251
```

```r
sd(results2$thetastar)  #the standard deviation of the theta stars is the SE of the statistic (in this case, the mean)
```

```
## [1] 0.9059806
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
## t1*      4.5 0.007707708   0.9387935
```

One of the main advantages to the 'boot' package over the 'bootstrap' package is the nicer formatting of the output.

Going back to our original code, lets see how we could reproduce all of these numbers:


```r
table(sample(x,size=length(x),replace=T))
```

```
## 
## 0 1 3 4 5 7 8 
## 1 2 2 1 1 2 1
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
## [1] 0.0368
```

```r
se.boot
```

```
## [1] 0.9310062
```

Why do our numbers not agree exactly with those of the boot package? This is because our estimates of bias and standard error are just estimates, and they carry with them their own uncertainties. That is one of the reasons we might bother doing jackknife-after-bootstrap.

The 'boot' package has a LOT of functionality. If we have time, we will come back to some of these more complex functions later in the semester as we cover topics like regression and glm.

