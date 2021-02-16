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
## 0 1 3 6 8 
## 1 3 1 4 1
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
## [1] 0.0173
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
## [1] 2.73547
```

```r
UL.boot
```

```
## [1] 6.29913
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
##    [1] 4.6 3.8 4.6 4.5 4.7 3.6 5.2 5.4 6.3 4.2 3.5 5.8 5.3 3.4 5.6 5.2 5.2 5.0
##   [19] 2.7 4.4 4.4 4.4 3.5 3.7 5.3 3.6 4.6 5.1 5.1 3.2 5.3 3.8 3.9 4.6 4.8 2.2
##   [37] 3.7 4.7 4.9 3.8 4.7 3.9 4.5 5.2 3.5 5.9 5.2 4.4 4.2 3.4 4.7 4.8 3.7 5.5
##   [55] 4.8 4.8 4.5 3.3 4.6 3.9 4.3 6.1 4.6 4.4 3.9 4.1 5.7 5.9 4.7 4.8 4.8 3.6
##   [73] 4.5 6.2 6.0 6.0 3.3 5.6 4.7 4.1 5.7 4.9 5.4 4.3 3.6 4.8 5.0 5.7 5.8 5.2
##   [91] 3.1 2.0 4.7 4.0 4.9 5.0 5.3 3.1 4.3 4.3 5.6 2.9 3.8 5.5 4.8 4.7 4.6 4.5
##  [109] 4.7 3.7 6.9 4.7 4.5 5.4 5.3 6.0 4.5 3.0 5.8 4.5 4.7 6.4 5.0 2.7 4.3 3.9
##  [127] 4.4 5.8 4.9 3.5 5.3 5.2 5.4 4.4 4.2 3.5 5.4 5.1 5.9 4.0 5.7 4.5 5.8 6.1
##  [145] 5.4 2.9 5.3 4.8 4.7 4.8 4.1 4.4 5.4 5.8 4.6 4.6 5.9 4.5 4.0 3.4 5.8 4.6
##  [163] 5.5 7.2 4.8 4.1 5.2 3.9 6.1 4.5 4.3 3.8 4.9 4.6 5.7 4.2 5.0 5.4 4.1 3.1
##  [181] 3.5 4.8 4.7 5.6 4.4 2.7 4.8 4.7 4.6 4.6 4.2 5.7 4.1 4.4 3.4 4.1 5.2 6.2
##  [199] 4.1 5.0 4.0 4.3 4.9 5.5 5.7 4.0 4.0 4.5 5.8 4.6 5.2 3.5 4.1 6.3 4.8 6.0
##  [217] 4.6 4.6 5.3 4.9 5.5 4.8 4.5 3.9 5.0 3.7 4.8 3.4 4.2 4.1 2.7 4.1 5.5 5.0
##  [235] 4.7 4.3 4.5 3.6 5.9 3.7 5.4 5.6 4.4 6.1 4.5 4.9 4.9 4.1 5.3 3.5 5.6 2.8
##  [253] 4.8 3.2 3.9 3.2 5.4 3.7 3.4 3.9 5.0 3.9 5.2 3.1 5.1 4.6 4.5 3.7 5.3 3.3
##  [271] 4.2 3.5 3.6 5.5 4.7 6.0 5.0 3.8 1.6 6.1 3.5 5.3 5.4 4.7 3.9 4.6 4.0 4.2
##  [289] 3.7 2.8 4.3 3.8 5.3 4.2 4.5 4.7 3.8 4.1 5.5 5.0 4.4 4.0 4.8 4.7 5.0 4.6
##  [307] 4.3 4.4 5.9 3.8 5.7 3.5 5.4 4.0 5.5 3.9 4.5 4.8 4.3 4.5 3.7 4.9 4.4 5.0
##  [325] 3.5 5.2 4.4 5.8 5.8 4.7 6.3 4.6 3.5 3.9 3.9 4.4 3.9 5.2 4.9 4.5 2.2 4.3
##  [343] 6.2 5.7 4.6 4.8 3.8 2.9 3.6 5.6 5.6 4.2 6.6 4.8 4.0 6.6 5.2 5.4 4.0 5.0
##  [361] 3.9 5.1 4.2 3.6 4.3 5.3 4.7 4.5 3.9 5.4 3.6 5.6 5.6 5.6 3.6 4.8 4.0 3.8
##  [379] 4.4 1.0 4.0 3.7 5.4 5.4 5.1 5.7 5.0 4.8 5.0 6.0 5.2 5.4 5.2 4.3 5.3 3.0
##  [397] 3.9 3.9 3.9 4.0 5.1 4.9 4.8 5.2 4.4 3.3 4.6 3.4 4.5 4.8 4.7 4.2 3.3 3.8
##  [415] 4.7 3.8 4.1 4.1 4.9 5.1 4.4 5.0 2.4 3.9 5.9 5.2 5.3 2.8 5.1 4.9 3.5 2.8
##  [433] 6.6 5.8 5.1 2.4 5.0 2.9 4.4 4.5 4.8 5.5 5.1 4.8 4.0 5.6 4.3 4.1 3.6 4.0
##  [451] 3.2 4.2 4.7 4.9 4.9 3.1 6.2 4.7 4.7 5.0 6.6 5.3 5.3 4.8 5.0 3.4 4.6 4.6
##  [469] 5.3 4.3 4.8 6.4 4.4 4.1 4.6 3.9 5.5 4.0 4.7 3.9 5.3 4.8 4.3 4.3 5.2 4.6
##  [487] 4.3 6.4 5.3 3.7 6.5 4.2 5.8 4.1 5.1 5.6 4.9 5.1 3.8 5.7 4.5 5.6 4.5 5.2
##  [505] 3.6 4.4 5.5 4.4 4.8 4.6 4.5 4.0 3.5 4.1 4.6 4.2 4.3 4.1 5.2 3.8 5.4 6.1
##  [523] 3.2 3.0 3.7 6.0 6.6 5.4 4.6 4.5 4.6 5.3 5.6 4.8 5.8 4.5 3.2 5.2 4.3 6.1
##  [541] 4.5 3.2 4.2 5.6 4.7 4.2 3.8 5.2 5.2 3.5 4.7 4.9 4.3 5.1 4.2 4.2 4.4 4.0
##  [559] 4.5 4.7 3.8 4.2 4.6 3.4 4.6 5.0 5.0 6.6 6.6 3.4 6.3 3.4 4.4 4.5 5.3 4.0
##  [577] 6.1 3.9 4.7 6.3 4.5 4.2 5.1 2.8 6.4 3.9 3.4 4.3 5.1 5.1 3.6 3.3 5.5 4.7
##  [595] 4.6 3.1 3.3 3.9 5.2 5.8 5.0 3.7 4.3 4.8 4.0 5.4 3.9 4.0 4.4 4.2 3.6 4.2
##  [613] 3.9 4.8 4.9 4.0 5.5 4.1 4.4 4.7 4.5 5.2 6.2 4.4 4.2 5.2 5.8 4.8 4.6 3.2
##  [631] 4.8 4.3 4.7 3.2 3.9 5.2 3.5 4.8 4.9 5.3 4.4 4.5 4.9 3.0 3.3 4.2 5.6 3.8
##  [649] 4.9 4.7 5.6 3.6 3.5 5.2 3.8 3.2 4.4 3.6 4.0 4.6 4.9 4.5 4.2 4.6 3.8 5.7
##  [667] 5.4 4.9 5.6 4.4 4.1 5.1 5.0 2.9 5.1 4.4 3.2 6.0 4.1 4.1 5.2 4.5 5.4 4.3
##  [685] 6.1 3.9 3.2 5.6 3.1 3.2 3.8 2.9 4.1 4.6 4.8 4.6 5.1 3.3 3.1 5.3 3.1 5.1
##  [703] 5.1 4.7 4.8 3.4 4.5 5.0 3.8 4.7 4.2 5.2 3.9 4.2 4.6 5.3 6.5 4.7 5.9 4.2
##  [721] 4.2 6.0 4.8 6.0 5.4 4.6 4.5 4.9 3.7 6.7 3.7 3.8 4.7 5.1 3.6 5.2 5.4 5.1
##  [739] 4.0 4.6 4.6 4.7 4.5 4.9 5.4 4.4 4.1 4.4 4.7 5.8 3.8 3.2 3.5 4.1 5.5 6.1
##  [757] 5.1 5.7 5.6 3.8 5.0 5.3 3.8 4.6 4.8 6.5 4.1 4.3 4.2 5.2 5.5 5.4 4.4 4.0
##  [775] 4.2 3.7 4.6 3.6 3.8 4.8 6.0 4.1 5.6 4.1 6.1 5.1 3.2 4.3 3.4 5.6 5.4 5.1
##  [793] 4.6 4.8 5.4 3.2 4.5 4.7 4.9 2.4 4.6 6.4 4.7 5.3 4.6 3.8 3.2 4.5 4.1 3.3
##  [811] 4.8 4.7 3.3 5.0 4.6 4.0 5.5 3.5 4.6 3.6 4.2 6.4 5.2 3.5 3.9 4.8 4.4 3.7
##  [829] 4.2 5.2 3.7 3.5 5.8 3.8 4.0 6.6 4.9 4.4 3.0 4.0 4.1 2.7 4.9 5.5 4.8 5.4
##  [847] 5.8 3.6 3.5 3.4 3.9 4.1 4.9 5.7 4.7 4.6 3.1 5.1 2.7 6.9 3.8 4.2 4.1 3.8
##  [865] 4.9 4.3 5.0 4.1 5.6 5.3 6.2 4.2 5.7 6.1 5.9 5.0 6.5 4.8 5.2 5.1 5.3 5.8
##  [883] 3.8 4.4 5.1 5.0 6.4 4.9 4.3 4.2 5.5 3.6 7.7 4.4 4.2 3.3 3.8 4.4 4.3 4.8
##  [901] 5.1 5.8 2.8 3.4 4.5 5.0 4.3 4.9 3.0 4.7 4.3 3.5 5.7 4.1 4.4 3.3 4.4 4.7
##  [919] 3.5 3.4 4.2 4.4 4.4 4.8 6.1 5.5 3.9 5.3 4.5 6.4 4.7 4.4 5.2 4.7 4.9 3.7
##  [937] 4.2 4.7 2.3 5.3 3.8 5.5 6.1 5.0 4.0 3.2 4.8 4.0 4.6 5.2 5.2 4.8 3.6 3.6
##  [955] 4.0 4.0 4.9 4.9 6.0 3.7 4.7 3.8 6.3 4.1 6.5 5.1 2.9 5.5 5.9 5.5 5.6 3.7
##  [973] 2.6 4.7 2.4 5.5 5.7 2.7 5.4 5.0 3.5 3.4 6.7 3.5 4.9 5.3 4.2 4.4 4.8 4.4
##  [991] 5.3 3.4 4.0 4.2 5.8 5.1 4.3 5.0 5.9 4.1
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
##   2.9   6.4
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
##    [1] 3.6 4.5 3.4 5.3 3.9 4.8 5.1 5.1 4.1 4.0 4.7 3.9 2.2 3.9 5.5 5.0 3.1 5.0
##   [19] 4.5 3.3 5.8 3.8 5.9 4.5 4.8 3.9 5.0 5.3 4.7 2.9 2.5 4.9 3.8 5.1 4.7 4.2
##   [37] 4.6 6.3 4.2 3.5 4.9 2.1 2.9 4.0 3.7 4.8 5.4 4.6 4.8 5.2 5.7 4.7 4.4 4.4
##   [55] 5.0 3.9 3.6 6.5 3.1 3.8 5.4 4.4 5.2 3.5 4.7 4.4 4.9 4.3 6.0 4.7 3.8 3.9
##   [73] 3.7 5.6 3.4 5.8 5.2 4.1 3.2 4.7 4.4 3.8 4.8 4.3 4.8 2.2 3.9 4.6 5.1 4.2
##   [91] 4.6 4.2 5.3 3.7 4.2 4.1 4.9 3.0 4.9 4.6 2.1 4.8 5.4 6.3 6.1 3.2 5.4 5.3
##  [109] 2.7 2.0 3.0 5.7 4.8 3.9 4.0 4.4 3.8 3.4 5.1 6.6 5.5 5.1 3.8 4.4 3.8 4.0
##  [127] 4.1 3.3 4.4 5.3 5.3 4.0 4.9 5.9 3.9 4.9 3.8 4.7 5.5 4.6 3.3 4.9 3.6 3.8
##  [145] 3.6 3.2 3.4 3.6 4.7 5.6 3.8 5.1 3.8 4.9 3.1 3.8 4.6 5.2 4.0 5.7 3.2 4.8
##  [163] 4.4 2.4 4.3 5.1 5.2 1.8 3.5 5.2 4.2 4.8 4.9 4.4 5.7 3.5 5.3 5.4 3.9 4.4
##  [181] 5.4 4.8 4.5 5.5 5.2 5.5 5.5 3.5 4.9 5.6 4.9 4.2 4.2 3.8 5.5 4.1 4.2 5.9
##  [199] 5.6 3.3 2.5 5.8 3.5 3.1 4.2 4.3 3.9 4.6 3.3 2.4 4.7 3.6 4.5 3.7 5.9 5.0
##  [217] 4.0 4.4 4.7 2.6 5.5 4.1 4.2 5.1 3.1 3.0 5.0 2.7 5.4 4.4 5.1 4.7 2.2 5.2
##  [235] 5.4 3.3 3.7 4.6 4.7 4.6 5.2 4.1 5.8 4.4 4.5 6.2 2.9 4.9 3.3 2.4 5.6 2.9
##  [253] 3.9 5.1 3.4 4.7 4.6 4.0 4.8 5.3 4.5 5.3 6.0 6.2 3.3 5.2 2.9 5.1 4.5 6.3
##  [271] 6.0 3.4 5.6 4.9 4.2 4.0 3.7 5.1 5.0 4.4 4.7 5.0 4.6 4.7 5.0 3.1 5.5 4.9
##  [289] 3.3 5.0 4.7 5.2 4.6 5.3 5.7 4.9 5.4 5.0 4.2 4.7 4.1 2.4 4.0 5.4 5.8 5.2
##  [307] 6.7 4.3 5.0 5.8 3.6 5.1 3.3 3.5 5.3 4.2 2.7 4.4 5.3 2.3 4.6 3.8 5.1 4.4
##  [325] 6.2 4.4 5.5 2.8 5.1 4.9 2.3 5.3 5.6 3.2 3.8 5.1 3.9 4.0 4.0 4.3 3.9 5.6
##  [343] 6.1 4.4 4.8 3.3 4.7 4.5 4.4 4.2 4.2 4.3 4.6 3.6 3.8 5.6 4.6 4.6 3.7 4.3
##  [361] 5.4 6.3 4.5 4.7 4.3 4.4 4.0 5.7 3.5 4.2 4.6 5.7 4.5 6.3 4.8 5.4 4.9 5.0
##  [379] 3.7 7.0 5.0 4.9 6.0 4.6 6.2 4.6 3.8 2.5 4.6 3.4 5.9 4.1 3.7 5.2 6.3 3.1
##  [397] 4.5 5.8 4.2 4.9 5.1 5.5 5.2 3.3 3.6 4.2 5.0 4.1 5.0 3.8 6.1 5.5 4.3 3.8
##  [415] 4.1 3.6 5.1 4.0 5.2 4.2 2.2 6.1 3.7 5.8 5.4 3.6 3.8 6.0 3.9 5.3 4.3 4.0
##  [433] 4.1 5.5 6.3 4.5 4.6 3.8 5.8 5.7 5.4 3.7 3.5 4.0 3.4 4.2 4.5 4.5 5.0 3.5
##  [451] 4.4 5.1 4.3 4.7 5.1 3.3 3.4 3.6 5.1 2.9 4.6 6.0 5.5 5.7 4.3 5.4 4.3 5.9
##  [469] 5.9 4.0 5.4 3.4 4.6 4.3 4.7 5.3 5.2 5.6 4.0 3.2 4.6 4.0 2.9 3.8 3.2 4.4
##  [487] 3.8 4.2 5.3 3.7 5.4 5.4 5.5 5.2 5.5 3.5 5.4 5.6 3.5 4.6 4.0 4.6 4.1 4.9
##  [505] 5.2 4.1 3.3 4.9 5.3 3.1 1.8 4.5 4.9 3.8 5.3 4.8 4.9 5.1 4.3 3.3 5.8 4.2
##  [523] 4.9 4.2 5.4 4.6 5.0 4.5 4.0 5.4 4.6 4.5 4.7 5.1 3.5 4.6 4.4 4.4 5.0 4.2
##  [541] 3.2 3.7 3.9 5.2 3.5 5.7 4.7 5.3 3.3 4.8 3.2 4.9 5.4 6.3 3.3 5.2 5.0 5.1
##  [559] 5.4 5.6 4.3 4.5 4.5 5.3 4.3 3.8 4.8 4.8 5.6 4.7 5.3 4.5 4.5 3.8 3.0 5.7
##  [577] 4.8 4.5 4.9 4.1 3.6 5.6 5.2 3.9 4.3 4.8 3.9 3.4 5.3 4.3 3.3 4.0 6.2 4.8
##  [595] 6.8 4.4 4.3 5.0 5.8 4.6 4.3 4.4 4.2 2.5 5.5 5.1 6.5 3.0 4.4 4.2 4.2 3.5
##  [613] 4.7 4.3 3.3 4.2 5.7 3.6 4.2 3.9 5.0 5.7 3.6 5.5 4.8 3.2 5.1 3.9 5.0 3.9
##  [631] 4.6 4.2 4.1 3.0 3.6 3.8 3.4 5.3 4.4 3.5 5.6 3.9 4.8 3.8 4.4 6.4 4.8 5.1
##  [649] 4.8 4.9 4.2 3.8 3.7 3.1 4.3 4.2 4.0 6.1 4.7 5.2 3.1 4.8 5.2 4.8 3.8 2.9
##  [667] 5.6 3.1 3.1 4.0 4.0 6.2 4.9 3.1 5.7 3.2 3.9 4.4 3.9 4.0 4.4 4.3 4.5 3.1
##  [685] 5.7 4.7 4.9 4.6 4.7 3.9 4.2 4.5 4.3 4.8 5.2 2.6 4.5 4.7 3.4 5.7 4.0 3.8
##  [703] 4.8 5.6 5.5 2.9 4.3 4.9 2.7 3.3 2.5 5.4 6.2 4.6 4.4 3.7 4.0 4.0 4.2 5.9
##  [721] 3.2 5.1 5.7 4.1 5.3 3.3 3.8 4.8 4.6 3.7 5.0 4.7 5.5 3.0 4.5 5.1 3.6 5.0
##  [739] 4.1 4.6 4.5 4.7 4.3 3.7 5.2 4.4 4.5 4.0 4.8 4.2 5.1 4.2 4.4 5.9 4.1 3.9
##  [757] 3.6 4.6 5.8 5.1 5.1 4.3 5.2 5.3 5.5 4.1 3.9 5.0 3.1 4.0 4.8 4.0 3.2 3.7
##  [775] 3.7 3.8 5.2 3.1 4.9 3.6 4.0 2.7 4.4 5.9 5.7 6.1 2.3 4.0 5.0 5.8 5.1 5.2
##  [793] 4.0 4.7 5.1 3.9 3.7 4.8 4.6 3.9 3.7 3.5 5.7 4.9 4.9 3.0 6.4 4.3 5.8 2.8
##  [811] 3.8 4.8 3.9 3.9 4.5 4.3 5.2 4.7 5.9 3.0 4.3 4.9 4.6 5.1 6.2 4.3 6.1 4.3
##  [829] 5.7 4.7 5.0 3.4 6.3 5.5 4.3 4.0 5.9 4.1 4.4 4.6 3.5 4.3 3.5 5.5 5.2 5.1
##  [847] 4.7 2.6 5.0 5.9 5.2 3.3 5.0 4.6 5.4 5.7 5.4 2.6 5.4 4.4 4.6 4.3 3.5 5.4
##  [865] 4.0 5.2 3.2 4.3 4.1 2.6 4.8 6.1 3.9 5.6 3.4 3.3 6.5 2.7 5.1 4.5 3.4 3.1
##  [883] 3.5 5.8 3.6 2.9 4.1 6.5 2.1 6.1 3.7 4.6 4.6 3.6 5.0 5.6 3.2 5.1 4.5 6.0
##  [901] 4.2 5.7 5.9 4.9 3.3 4.5 3.8 2.1 4.5 4.8 5.6 4.4 6.0 3.3 3.3 3.8 3.5 3.4
##  [919] 4.9 5.4 3.6 6.2 4.0 4.8 5.5 3.0 6.5 4.6 6.3 4.6 5.5 4.3 4.7 4.5 4.2 5.1
##  [937] 4.8 2.1 4.7 4.0 5.2 1.7 6.0 4.1 4.9 4.4 4.1 3.6 4.9 5.0 4.5 4.9 5.6 3.9
##  [955] 5.4 2.9 3.9 3.7 4.6 4.9 4.0 4.4 6.3 5.5 4.4 3.5 4.9 3.3 4.7 4.4 4.9 3.6
##  [973] 5.1 6.1 5.0 4.1 5.8 2.8 3.4 5.1 5.7 3.8 2.4 4.3 4.2 3.4 4.5 5.0 4.6 4.8
##  [991] 3.0 5.6 4.4 4.5 6.0 6.2 5.7 5.5 4.2 2.8
## 
## $func.thetastar
## [1] -0.017
## 
## $jack.boot.val
##  [1]  0.49181287  0.41250000  0.28066465  0.12339833 -0.01855956 -0.06880223
##  [7] -0.15890805 -0.33707865 -0.44081633 -0.54894260
## 
## $jack.boot.se
## [1] 1.012311
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
##    [1] 4.5 2.4 2.8 2.8 4.6 4.7 4.3 3.1 4.4 4.9 4.3 4.3 4.2 3.5 5.4 3.0 4.0 5.4
##   [19] 4.6 5.7 4.3 3.5 4.8 5.0 3.1 4.2 4.0 4.6 3.3 5.5 4.6 5.4 2.3 6.2 2.6 5.4
##   [37] 4.2 3.8 4.0 4.2 3.3 5.0 5.4 5.8 4.0 5.1 4.2 3.7 4.7 3.9 5.4 5.3 4.8 5.0
##   [55] 6.7 4.7 5.0 5.3 3.2 4.9 4.8 3.8 5.0 5.7 4.1 3.0 4.4 4.4 4.6 5.0 5.7 4.5
##   [73] 5.4 5.0 3.7 5.1 4.7 3.8 6.0 4.4 5.0 4.3 5.2 2.9 5.5 5.0 4.7 4.1 4.1 3.8
##   [91] 3.2 4.2 3.9 4.9 3.2 3.8 3.3 4.4 3.0 5.1 3.3 4.6 3.4 6.3 5.2 4.5 3.6 4.3
##  [109] 4.3 3.9 4.1 4.8 3.6 4.6 4.3 4.4 4.9 3.0 4.1 4.5 4.0 3.4 6.2 4.3 4.2 6.4
##  [127] 4.9 4.6 5.7 4.3 4.0 5.7 3.7 5.3 4.1 3.7 4.3 4.0 4.2 4.8 4.5 6.1 5.5 5.1
##  [145] 4.8 5.0 5.1 3.3 4.3 4.2 5.6 5.3 6.3 4.1 5.6 4.5 6.0 5.1 4.4 6.1 3.1 4.5
##  [163] 5.6 5.8 5.1 4.1 4.2 5.9 4.5 3.5 4.6 5.4 5.8 4.4 3.4 5.5 4.9 4.6 4.0 5.5
##  [181] 3.1 5.6 5.0 6.6 4.6 4.3 5.1 4.6 5.9 5.3 5.4 5.7 4.4 5.0 4.3 3.4 3.7 3.0
##  [199] 5.3 6.2 3.0 4.3 3.9 5.1 3.7 4.0 3.4 6.0 4.7 4.8 5.1 3.5 2.7 5.3 5.0 4.0
##  [217] 5.1 5.0 5.5 4.9 3.9 5.5 4.6 5.4 4.3 4.6 4.5 3.8 3.9 4.4 5.2 4.2 5.4 4.1
##  [235] 5.2 4.0 5.2 3.5 4.0 5.7 5.8 5.2 6.0 3.1 6.9 4.3 3.7 3.3 6.2 3.4 2.5 4.2
##  [253] 5.5 4.7 5.6 5.5 3.3 2.8 4.9 4.8 6.7 4.5 4.7 4.4 3.6 5.8 4.0 4.6 3.5 5.0
##  [271] 5.5 4.0 4.8 5.1 4.9 5.4 5.1 4.4 4.5 2.3 6.4 4.7 3.0 4.3 3.6 4.0 4.0 5.1
##  [289] 5.6 4.0 3.5 3.2 6.6 3.7 6.2 5.2 3.7 4.1 4.1 4.0 3.1 3.2 5.4 2.8 4.7 4.2
##  [307] 5.0 5.6 3.4 5.0 4.1 6.6 2.2 4.2 4.4 2.1 3.2 4.5 3.6 3.6 6.7 4.3 4.3 5.9
##  [325] 4.7 3.5 4.9 5.2 6.4 4.3 4.1 5.3 2.8 3.8 4.2 4.5 5.1 3.2 3.8 3.2 2.4 4.0
##  [343] 3.7 5.4 3.5 4.3 4.0 3.6 6.2 4.4 3.5 4.1 3.7 3.5 4.6 5.0 5.1 5.3 4.2 3.6
##  [361] 6.8 3.6 3.8 2.3 4.7 4.5 4.0 5.1 4.7 5.4 5.0 7.1 4.5 4.3 3.2 4.9 4.3 4.1
##  [379] 4.9 3.9 4.9 4.0 4.1 3.8 3.1 5.1 4.0 4.2 3.1 4.0 5.2 3.0 3.1 5.3 4.9 3.5
##  [397] 3.1 4.0 4.3 5.7 5.4 4.9 5.3 4.6 3.1 3.9 3.7 4.1 4.0 6.1 4.5 3.7 5.7 5.8
##  [415] 4.1 4.6 3.6 3.9 6.1 5.7 3.3 5.8 5.0 2.7 3.7 4.4 5.3 3.8 4.1 5.4 5.7 5.7
##  [433] 5.6 5.7 4.9 3.8 4.1 5.4 4.3 3.9 3.3 4.5 4.3 4.4 5.3 3.9 3.3 4.4 4.9 4.4
##  [451] 4.6 4.8 5.5 3.8 3.4 4.2 4.0 4.5 3.8 4.6 5.1 6.2 4.9 4.0 4.5 3.7 5.2 2.9
##  [469] 4.8 4.2 5.5 4.2 5.7 4.3 4.7 5.2 5.1 6.2 4.0 4.7 6.5 5.1 5.7 5.0 4.7 4.3
##  [487] 4.5 5.5 4.3 4.0 5.2 4.5 5.7 3.6 4.0 4.8 4.2 4.4 3.0 5.4 5.8 4.7 2.8 4.7
##  [505] 4.8 5.0 5.0 4.5 5.6 5.0 4.9 4.6 3.6 4.9 6.0 4.3 4.5 2.9 3.9 3.7 6.0 3.1
##  [523] 6.3 5.7 3.8 4.0 4.8 4.5 4.6 4.2 4.8 6.2 3.7 4.3 3.5 4.9 6.5 4.5 4.3 4.3
##  [541] 5.2 4.5 4.6 2.9 4.7 4.2 3.1 4.2 4.0 4.1 5.6 5.1 4.5 2.8 6.1 2.9 4.1 5.3
##  [559] 2.6 5.2 3.8 3.8 6.5 3.5 5.1 5.1 2.7 2.3 5.1 5.0 4.3 5.0 5.7 4.7 5.0 4.5
##  [577] 4.3 5.4 4.9 1.8 5.5 4.3 4.7 4.7 4.4 4.4 5.4 4.7 4.4 5.3 2.8 4.7 4.4 5.7
##  [595] 3.9 5.1 4.5 5.7 4.6 4.2 6.0 3.8 5.5 5.6 4.6 4.2 4.2 4.7 3.7 4.2 3.4 5.4
##  [613] 5.2 6.1 3.2 6.1 5.7 6.3 4.9 2.2 4.8 3.5 3.6 3.7 4.6 4.9 5.8 5.1 5.7 3.7
##  [631] 4.4 4.9 2.7 3.8 4.2 3.2 5.3 3.4 5.4 4.8 3.0 5.0 3.9 3.9 5.5 5.0 3.8 5.4
##  [649] 5.0 4.7 3.7 3.6 3.4 2.6 2.9 3.5 4.3 6.2 6.5 4.7 4.1 4.3 5.4 4.5 5.6 3.4
##  [667] 4.1 5.6 6.3 3.5 3.7 3.2 4.5 4.0 4.2 5.2 5.6 5.9 4.0 5.2 3.4 5.2 5.2 5.3
##  [685] 4.7 5.0 5.8 4.9 3.6 3.6 4.4 5.3 5.2 4.4 4.6 4.4 4.8 5.7 4.2 3.8 5.4 5.2
##  [703] 3.3 5.8 3.6 4.3 4.4 4.1 4.8 3.6 5.3 4.7 4.2 6.0 3.5 5.3 4.4 6.3 3.8 4.7
##  [721] 5.7 3.0 4.5 4.1 5.5 4.9 5.0 5.7 5.0 3.5 3.7 4.2 2.0 4.5 4.2 3.7 5.0 5.0
##  [739] 5.3 5.0 4.8 4.2 3.9 5.0 5.1 6.5 4.3 4.3 5.6 4.7 3.7 3.9 4.5 4.4 4.2 5.0
##  [757] 4.2 3.9 4.3 4.4 4.9 4.9 5.8 4.8 4.9 4.1 4.2 4.3 6.1 5.8 3.7 4.4 4.6 3.1
##  [775] 3.3 5.5 3.4 4.7 5.5 5.0 3.6 5.4 4.0 5.8 5.3 6.3 3.7 3.5 4.6 3.8 5.0 4.4
##  [793] 4.4 4.9 4.5 5.0 4.1 5.4 3.1 5.0 3.5 4.3 6.1 3.5 4.2 6.2 5.2 4.4 3.6 3.7
##  [811] 4.4 5.1 7.0 4.5 3.7 4.9 4.8 4.5 4.4 3.4 2.7 6.0 3.7 4.6 3.2 4.0 5.4 2.8
##  [829] 5.0 4.4 4.5 5.8 3.3 4.0 5.7 3.8 4.3 2.9 4.0 5.1 6.0 3.0 3.0 4.5 5.7 4.7
##  [847] 4.2 6.0 5.1 4.9 5.3 4.4 5.4 6.0 5.4 4.3 4.6 3.6 3.8 5.8 4.7 4.4 4.8 5.3
##  [865] 5.5 5.3 4.3 5.9 3.4 5.1 3.9 5.3 5.3 4.5 3.7 2.8 3.7 4.2 5.1 6.1 4.9 3.3
##  [883] 5.6 3.5 4.2 5.3 5.8 4.3 4.3 4.8 4.6 4.7 5.4 5.8 4.4 5.3 6.4 3.9 2.9 4.4
##  [901] 4.6 2.4 4.4 3.5 5.4 4.6 3.5 4.1 6.0 4.5 4.7 3.6 4.1 4.8 4.4 3.8 4.8 5.2
##  [919] 3.4 3.0 5.9 4.9 3.3 4.2 5.2 3.9 5.2 3.3 4.7 3.0 6.1 5.8 3.7 4.8 2.5 4.2
##  [937] 4.5 3.6 5.6 6.4 4.6 5.1 5.1 5.3 4.1 3.2 4.4 4.7 4.0 4.5 4.5 1.7 3.3 5.7
##  [955] 6.3 4.2 5.3 6.3 3.5 4.7 4.7 5.8 4.0 4.6 4.7 4.4 3.7 3.6 4.0 6.3 4.3 6.2
##  [973] 4.3 4.4 4.5 4.9 4.6 5.4 5.4 3.2 3.6 5.6 3.1 4.7 3.5 5.4 4.7 5.2 4.8 4.5
##  [991] 3.8 4.9 2.7 3.6 5.0 3.4 4.8 4.9 5.3 4.9
## 
## $func.thetastar
## 72% 
## 5.1 
## 
## $jack.boot.val
##  [1] 5.500 5.500 5.400 5.300 5.200 5.000 5.000 4.900 4.600 4.464
## 
## $jack.boot.se
## [1] 1.027318
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
## [1] 2.201447
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
##    8.038008   13.022080 
##  ( 3.522626) ( 5.888811)
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
## [1] -0.09221914  1.03657485  0.04843033  0.61887907 -0.43275047  1.03367385
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
##    [1]  2.208491556  2.180584016  2.520456255  2.533432593  2.535480901
##    [6]  1.292846428 -0.560395583  0.969623913  0.837808722  1.927433853
##   [11]  1.092198577  2.539083419  2.475416464  2.550564366  2.628985061
##   [16]  1.287708398  2.219990273  1.457703456  1.174637931  2.225439839
##   [21]  1.768513214  1.251828867  0.925966471  1.638903630  1.953954033
##   [26]  1.653090714  1.335987994  1.296329464  2.137471264  1.324832638
##   [31]  0.571944771  1.427209950  2.224638198 -0.582292061  0.843132616
##   [36]  1.349150922  2.523971578  2.609505296  1.453846666  1.970783272
##   [41]  2.581158831  2.178146452  1.165265868  2.572342245 -0.727610942
##   [46] -0.631144488  0.661497520  2.557504262  1.296705706  2.272375318
##   [51]  0.705016739  2.576073825  0.061136512 -0.155659721  1.227220618
##   [56]  1.316478911  0.251665256  2.113518437  2.568157884  2.541996769
##   [61]  1.339986457  1.264732107  2.234614215  1.166641092  0.513762694
##   [66]  1.321723665  1.302642789  2.192478398  1.212829872  2.006190755
##   [71]  0.902077062  1.272088275  0.842388520  2.198481663  1.422699207
##   [76]  0.897627900  1.575357996  1.928186398  2.182117634  1.434903662
##   [81]  0.838707565 -0.315358494  1.887609805  2.176012291  2.015439270
##   [86]  2.177423609  2.613520757  1.142278364  2.517902429  0.141710544
##   [91]  0.852930830 -0.070464502  1.145531096  2.541565997  0.758046772
##   [96] -0.240048610  1.946603237  1.706610816  0.856805281  1.164325215
##  [101]  2.590962038  1.447473343  2.559783524  0.171987667  1.468464134
##  [106]  0.272400544  0.164814180 -0.713254028  0.629080432  2.632856431
##  [111]  1.305204998  0.844030906  0.001111365  0.959261130  1.210494665
##  [116]  1.213117275  1.136944848  2.605247128  1.042320872  1.448341842
##  [121]  1.899890687  1.160538665 -1.081183016  1.327115113  1.167948054
##  [126]  1.329489662  1.809826404  2.284321183  2.254767421  1.215722661
##  [131]  1.832507595 -0.569055493  1.342862757  0.548497095 -0.023013667
##  [136]  2.032530815  0.346114339  1.916811555  1.342539546  1.843959833
##  [141]  1.386115995  1.304717369  2.582281344  2.532562209  1.154726100
##  [146]  0.740931412  1.435266866  0.685362583  2.532562209  2.211017549
##  [151]  1.305816799  2.439413475  1.177397817  2.146429181  1.209043285
##  [156]  2.213882552  2.227906384  1.034566661  2.529645737  0.766021413
##  [161]  2.217529387  1.048716608  2.616163609  2.266227129  1.325657697
##  [166]  0.415876981  1.428906963  2.524610786 -0.582292061  1.752506694
##  [171]  0.842092418  2.576484015  1.310884079  0.949463787  2.184822037
##  [176]  0.738581427  1.469877647  2.551247507  1.936822082  0.246547150
##  [181]  1.437447498  1.173658296  1.944921523  2.256179916  1.449147264
##  [186]  2.556135806  1.061649649 -0.041420181 -1.346862848  1.465539842
##  [191]  1.466703615  0.641876187 -1.520771786  1.459275271  2.233765369
##  [196]  1.291490501  1.315500793  1.830282086  1.445381028  1.926904327
##  [201]  2.518592185  2.588987903  1.312872802  1.475871840  2.257033807
##  [206] -0.229809265  2.440463048  2.600098780  0.733258735  2.288565203
##  [211]  0.689157455  1.338985481  1.447623511  1.273405541  0.761803690
##  [216]  0.831212576 -0.359915295  2.539365005  2.012419260  2.149593183
##  [221]  1.930547977  2.234583253  1.602061988  1.681564854  1.917997539
##  [226] -1.571699874  1.653816084  1.322479657 -0.062579184  1.306240758
##  [231]  1.207875570  1.200535852  2.570898928  1.285250487  1.290159168
##  [236]  1.184348051  1.342645988  0.666204839  1.900804679  1.290109562
##  [241]  1.437451283  2.566514146  2.588684693  1.437237003 -0.564763674
##  [246]  1.297134987  0.748570892  2.219624357  1.263451490  2.487236047
##  [251]  2.179076358  1.349760034 -0.088065887  1.251954957  1.953342865
##  [256]  1.946784417  2.212945068  0.839689924  1.137289742 -0.014643278
##  [261]  1.342642894  2.206800291  1.231890831  2.171591591  2.213882552
##  [266]  2.260590088  1.329014303  2.540867292  2.033751060  0.862130570
##  [271]  1.093742702  2.149107444  1.662402505  1.445042569  1.308865612
##  [276]  1.926893142  1.450102207  2.576981663  0.846247369  1.328231128
##  [281]  1.449776503  1.760675422  1.324049483  1.692731963  0.754199110
##  [286]  2.181780976  1.791906721  1.590644250  1.082253111  2.216273969
##  [291]  1.320653808  1.955048873  2.178317464  1.292079432 -0.056574249
##  [296]  1.764206075  1.333593522  1.175427306  2.316572295 -0.376575029
##  [301]  1.523162215  0.836152864  1.428965716  2.134483979  2.591057370
##  [306]  2.218739328  1.975076649  0.728709907  1.308952197  1.998395038
##  [311]  2.208143073  1.298946059  1.072591235  1.485748668  1.506012081
##  [316]  2.139974245  2.506467660  0.213679537  1.296891488  2.292523036
##  [321]  0.625821188  2.257085727  2.518041375  2.218024685  0.312690500
##  [326]  2.206014787  2.136806286  2.189100672  0.740488368  1.450076704
##  [331]  2.230979787  1.688935999 -0.300599134  2.476702903  0.352377818
##  [336]  1.985369832  1.114162725  2.203083602  1.978726448  1.759939389
##  [341]  1.159395064  2.566152974  2.283713695  1.196038219  1.439823171
##  [346] -0.567752242  2.546015636  0.381612914  1.781429961  1.274008579
##  [351]  1.452023375  0.878647560  2.193325791  1.148328878  2.263177895
##  [356]  1.447615228  1.922948719  2.281603363  2.458943546  1.292305576
##  [361]  2.284321183  0.840512525  1.188509449  1.849636961  2.528735480
##  [366]  1.603884566  0.717230515  1.035472951 -0.466623404  2.533699811
##  [371]  1.030626964  0.591532176  0.846465338  1.337987242  1.778155834
##  [376]  1.446718243  1.285290283  2.496406012  1.781429961  2.493270275
##  [381]  1.480566365  1.946250532  1.946288123  0.840560882  2.224668719
##  [386]  2.019922470  2.501428968  1.934385021  0.675937638  0.652362341
##  [391]  0.737384024  1.196150120  0.983787138  1.421687035  1.940985758
##  [396]  1.917772944 -0.330212476  1.780545518  2.294809529  0.745251782
##  [401]  2.172844342  2.597003658  1.364857358 -0.766563644  1.169253297
##  [406]  2.002996579  1.293739071 -0.418295375  1.286904948  0.851832237
##  [411]  2.525400147  2.544146919  1.760634977  2.080542635  0.857619633
##  [416]  1.468195630  0.749435739  1.018507234  1.325951578  2.012089177
##  [421]  1.592609319  2.186365725  1.126566350  2.160646637  0.735096838
##  [426]  1.443534226  2.265815717  1.075306699  0.047681485  2.580975862
##  [431]  1.297164596  1.337810693  0.619985136  1.840503266  1.711374380
##  [436]  0.941718053  2.102577372  2.216587719  1.462893997  2.597017605
##  [441]  0.302652210  2.181345648  0.767215627  1.882203079 -0.381372564
##  [446]  1.291752346  1.947929414  2.175462903  2.529946578  2.156920445
##  [451]  1.801954509  2.259448803  0.545910727  1.766774383  2.172492939
##  [456] -0.139891084  1.522654402  0.719567265  2.129914994  2.503668087
##  [461]  2.244545809  2.229805606  0.763071558  1.630550326  0.576507803
##  [466]  1.847233519  2.185624418  1.470499530  1.350391754  1.610791481
##  [471]  1.327436723  1.324213851  1.930454907  2.218292812  1.488835136
##  [476]  1.302731569  2.502511846  1.158189962 -0.563676957  2.176183047
##  [481]  2.202528454  2.584808041  1.671096779  0.205776498  1.171315747
##  [486]  1.768457238  0.508943481  2.291647421  2.333450253  1.478692551
##  [491]  2.578308549 -0.636334345  2.311655092  2.285893775  1.345114983
##  [496]  2.314215814  1.515238014  1.438638797 -0.004116973 -0.817501933
##  [501]  2.219764486  0.720795335  1.469686378  0.760200632 -0.518768102
##  [506]  1.478906077  1.129317233  2.231268524  2.053147732  2.538098787
##  [511]  1.336011271  2.184365400  0.623062912  1.608192348  1.984003522
##  [516]  1.670095463  1.983104446  1.112958300 -0.439887778 -0.111024586
##  [521]  2.546684184  2.272847554 -1.506995319  2.564275551  2.640708164
##  [526]  2.448187446  1.026159617  2.213754420  1.868229688  2.146753584
##  [531]  1.442217421  1.322652390  1.875306077 -1.013790699 -0.464068234
##  [536]  2.577592790  2.608660329  2.124695784  1.443161200  1.273808821
##  [541]  1.465366484  1.106920186  1.189416869 -0.220012767  0.078406584
##  [546]  1.804982267  1.238339331  1.328765418  1.500769517  2.343012485
##  [551]  1.351420112 -0.443142049  2.611399724  1.121244995  1.434990221
##  [556]  1.938234923  1.316227481  1.239601979  2.205122850  1.332849830
##  [561]  2.541917020  2.103945795  1.788148747  1.905082814  2.211768036
##  [566]  2.204933055  1.244773840  1.340586603  1.456181661 -0.453735802
##  [571]  1.128138265  1.379169355  1.743612788  1.325859263 -0.977845372
##  [576]  2.525646907  2.520683790  1.484414965  1.942387187  1.191490509
##  [581] -0.638175918  1.126245455  2.541159116 -0.403582355  2.198705110
##  [586]  1.497323985  1.836562591  1.755150760  1.095297463  2.243783890
##  [591]  2.139904222  1.020496194  2.081200939  1.960549056  1.068657932
##  [596]  1.711708826  1.907996331  1.405998502  2.563841051  2.281226252
##  [601]  1.198214067  1.872770331  2.232412737  1.466885028  1.307582330
##  [606]  0.733336934  0.067706477  2.553608228  1.862434906  2.483638884
##  [611]  2.213958698  1.324888223  1.762683031 -0.523129846  0.180915052
##  [616]  2.302913504  1.286763807  1.754504795  1.651565003  0.071868427
##  [621]  0.336906180  2.310731929  1.873841588  1.789603118  2.536748970
##  [626]  0.171524327  1.323095903  2.230440281  1.458515594  1.676432325
##  [631]  0.272111744  2.327737359  2.051376012  2.559394790  1.738358789
##  [636] -0.064278593  0.305300345  0.005328101 -0.321548879  0.025394684
##  [641]  1.474902595  0.760681492  1.212829872  1.007151555  1.282734391
##  [646]  2.575252369  1.307286774  2.192438319  1.890742566  2.245550629
##  [651]  2.406505862  0.998010511  1.243411211  0.922429558 -0.575021116
##  [656] -0.821846966  0.048912673  1.752867069  2.469429518  0.714424123
##  [661]  1.572877692  1.265256670  0.406235099  2.488619527  1.223770597
##  [666]  1.476390808  0.581347124  0.843215109  0.561541045 -1.009790908
##  [671]  0.204771138  2.532095376  1.402751752  0.399572782  2.205274264
##  [676]  1.230748119 -0.128431229 -0.387456679  2.008062383 -1.503196798
##  [681]  0.580734794  2.205010886  2.086843434  1.494772443  2.460938734
##  [686]  1.444615100  1.254140373  1.406355329  2.507048808  2.574871862
##  [691]  0.892720274  1.948110757  1.320897403  0.833727614  1.430881926
##  [696]  1.805283230  2.217630177  2.575099449 -0.017089376  2.192922729
##  [701]  2.146522512  2.488584689  1.874724822  1.443715870  2.264072234
##  [706] -0.083865104  1.143503040  2.195064039  2.195765677  1.467731586
##  [711] -0.259173934  2.278605172  1.891063303  1.477982815  1.357168230
##  [716]  2.360868780  2.595564643  1.123465256  1.779292949  0.878324450
##  [721]  2.561523891  1.305403217  0.829068947  0.799919401  1.423779681
##  [726]  0.811201696  1.472963662  2.546043022  1.736820181  2.291037509
##  [731]  1.664677581  1.942573955  2.294702616  1.984888357 -0.534352005
##  [736]  2.236978150  2.315936223  2.584116155  1.305611199  2.621232634
##  [741]  1.296351826  1.310328993  1.597382883  1.071787644  1.941189495
##  [746]  1.280820476 -0.185558504 -0.032302042  2.137984166  2.001116445
##  [751]  0.840074952  0.166833971  1.466147967  0.194309438  1.445751202
##  [756]  1.138645024  1.320510409 -0.789323088 -0.291819410  1.566741820
##  [761]  1.921315718  2.214625874  1.226550841  2.519954493  2.061103289
##  [766]  2.585791288  1.474902595 -0.791094616  2.471061825  2.611391994
##  [771]  1.382311048  2.038613600  1.465103875  1.892398858  2.530785922
##  [776]  2.147574179  1.455949889  2.329665644  1.021398306  2.594383490
##  [781]  2.546894925  2.211512077  2.147456736  2.530750117  2.545640252
##  [786]  0.783597416  2.605708950  1.210138985  1.674534874  0.642761106
##  [791]  0.761742877  1.581749638  0.560973936  1.257961469  1.292487569
##  [796] -0.165763994  2.119754541  1.450292870 -0.037621222  2.558180465
##  [801]  2.588715439  0.761299192  1.443901968  1.602787425  0.846771609
##  [806]  1.465539842  0.687640729  1.169253297  2.626959660  1.976348482
##  [811]  1.042074684  1.166913843  2.182322746  0.937193988  0.370184569
##  [816] -0.463191961  1.549144375  1.459340418  0.124004127  1.238645994
##  [821]  1.767330462  0.824963500  2.541160146  2.020586678  2.626511503
##  [826] -0.683666739  0.145338432  2.135805092  0.833543187  1.084676558
##  [831]  0.296115496  2.220669722  1.447911760  1.443634143  1.985126683
##  [836]  1.942040184  1.979008288  2.529022459  1.099864948  1.522246166
##  [841]  1.472211845  2.041643416  1.583478307  1.654541164  2.039670879
##  [846]  2.226411394  1.804691028  1.680840786  1.170208800 -0.019683016
##  [851]  2.194279222  1.553484482  2.216789741  2.287829685  1.493003800
##  [856]  1.959601253  1.477869417  2.608750854  2.473211009  1.917772944
##  [861]  1.211445035  1.307366801  2.148798120  1.983856541  2.130142459
##  [866]  0.859175830  2.207869990 -0.336932510  2.553282221  1.301960488
##  [871]  1.687002408  1.925671953  2.175951040 -0.572605962  1.165960978
##  [876]  0.918910556  2.584971015  0.342662741  1.464609807  0.110744261
##  [881]  2.557314698  0.997769128  1.922925961  2.271713277  2.268586549
##  [886]  2.547804894  2.612681432  2.215551769  2.545485889  1.930959244
##  [891]  1.834629731  2.517866652  2.559343756  1.205104234  0.597688965
##  [896]  1.310675636  1.339704685  1.760501088  2.312904360 -0.770273983
##  [901]  1.662770097  1.336682654  0.297128770  0.758398836  1.936163217
##  [906]  2.300659694  2.555361949  0.638388305  1.688173520  0.521817698
##  [911]  0.843506828  1.180803359  1.463925199  1.946784417  1.805001257
##  [916]  2.558191427  2.233463646  2.272201860  1.983550798  2.142759509
##  [921]  1.861410762  2.185624418  1.292348802  1.816954971  1.365484002
##  [926]  2.218149986  2.269709880  2.559825105  1.648054373  2.381165417
##  [931]  0.681101844  1.950296635 -0.300599134  1.487088457  2.262759924
##  [936]  2.095121801  1.062902722  1.210107442  0.466558796  1.494329977
##  [941]  1.312714273  0.849597974  0.863204153  1.592390920  2.203392356
##  [946]  1.300165071  2.559246625  2.245765864  2.440552117  0.996124217
##  [951]  1.800365356  2.547793239  1.572127886  1.959791128  1.342334695
##  [956]  1.158379308 -0.175140687  2.349801103  1.448468045  1.325839358
##  [961]  1.477376627  1.153239391  1.066867248  1.633605553  1.104518819
##  [966]  2.449166583  1.851884079  1.422539018  2.028440993  1.433833370
##  [971]  1.521287683  1.424842308  1.173818691  2.064979696  0.518883362
##  [976]  0.985612486  0.741722684  1.322558477  2.490093657  1.920819805
##  [981]  1.917443057  0.569614787  0.743221381  2.501243132  2.274185298
##  [986]  2.486089958  0.647016653  1.181184816  1.170763609  0.049382184
##  [991]  2.500533507  1.183934355 -0.446890377  2.549105964  2.305096922
##  [996]  0.994100241  1.533416679  1.143468395  1.355189616  1.957468002
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
##   0.61725959   0.26652318 
##  (0.08428203) (0.05959363)
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
## [1]  0.27495840 -0.56382340 -0.33750352  0.38054255 -0.03441034  1.07925833
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
## [1] 0.0078
```

```r
sd(results2$thetastar)  #the standard deviation of the theta stars is the SE of the statistic (in this case, the mean)
```

```
## [1] 0.9161651
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
## t1*      4.5 0.01981982   0.9333541
```

One of the main advantages to the 'boot' package over the 'bootstrap' package is the nicer formatting of the output.

Going back to our original code, lets see how we could reproduce all of these numbers:


```r
table(sample(x,size=length(x),replace=T))
```

```
## 
## 0 2 3 5 6 7 
## 1 2 1 2 2 2
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
## [1] 0.0158
```

```r
se.boot
```

```
## [1] 0.903098
```

Why do our numbers not agree exactly with those of the boot package? This is because our estimates of bias and standard error are just estimates, and they carry with them their own uncertainties. That is one of the reasons we might bother doing jackknife-after-bootstrap.

The 'boot' package has a LOT of functionality. If we have time, we will come back to some of these more complex functions later in the semester as we cover topics like regression and glm.

