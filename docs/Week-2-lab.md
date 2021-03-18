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
## 0 1 2 3 5 7 8 9 
## 1 1 1 2 2 1 1 1
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
## [1] -0.0207
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
## [1] 2.726462
```

```r
UL.boot
```

```
## [1] 6.232138
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
##    [1] 5.3 3.2 4.9 5.3 3.9 3.7 3.4 2.8 4.4 4.8 4.9 4.2 3.9 5.0 2.2 4.0 4.7 3.7
##   [19] 4.0 5.2 6.3 3.4 2.7 4.6 4.9 4.6 5.4 6.1 5.2 5.4 4.4 4.4 3.7 3.7 4.5 2.9
##   [37] 3.5 5.3 4.7 5.3 5.0 3.5 4.1 3.7 5.0 5.4 5.1 3.6 4.6 3.7 6.2 5.1 5.6 4.0
##   [55] 5.0 4.6 3.4 4.6 3.8 5.7 3.4 4.1 4.6 5.4 2.2 4.5 5.0 3.0 4.2 2.9 5.1 4.5
##   [73] 4.3 5.0 6.1 5.2 4.2 4.0 5.8 3.8 4.9 3.6 4.5 2.2 4.6 4.5 4.9 4.7 5.2 4.2
##   [91] 3.5 3.4 5.3 3.4 5.3 3.3 3.5 4.3 4.5 4.5 3.8 4.8 5.6 5.7 4.9 5.0 6.5 4.1
##  [109] 4.2 4.5 3.3 4.3 5.5 4.6 3.2 4.2 3.2 4.4 4.8 4.4 4.4 4.9 5.0 4.2 5.7 4.6
##  [127] 4.4 3.4 4.8 3.4 3.8 4.4 4.7 3.7 5.7 4.1 5.0 5.1 4.2 5.5 4.7 6.3 3.8 3.0
##  [145] 4.4 4.1 4.4 1.7 3.5 5.0 4.0 4.5 4.3 4.2 4.3 4.6 6.3 3.7 5.0 4.7 3.9 3.9
##  [163] 3.2 4.8 4.5 3.6 4.0 4.5 5.3 5.2 6.6 4.8 5.1 3.8 2.6 3.1 4.7 3.3 3.9 4.0
##  [181] 6.2 4.0 4.3 3.9 4.9 4.0 5.6 3.7 4.0 3.5 5.3 3.7 4.5 4.1 4.5 5.4 4.1 5.0
##  [199] 3.4 4.7 4.8 3.1 5.1 5.4 6.4 4.1 4.4 3.5 4.2 5.8 4.1 6.1 2.9 4.0 3.8 4.5
##  [217] 2.8 4.2 4.4 5.6 4.9 5.4 4.3 4.1 4.8 5.4 5.4 5.4 3.0 5.4 3.3 5.4 4.0 4.0
##  [235] 5.1 3.7 5.0 5.0 4.9 2.8 3.1 3.1 5.5 5.8 6.2 3.7 5.1 7.1 4.2 5.8 4.8 4.6
##  [253] 3.4 4.7 3.8 5.1 4.5 3.8 5.9 4.4 4.8 4.3 5.6 4.4 4.2 5.2 3.7 4.5 5.9 5.2
##  [271] 4.7 6.6 6.4 1.7 3.5 5.1 3.3 2.7 5.5 4.5 5.5 5.2 4.5 3.1 3.8 3.2 2.9 4.4
##  [289] 4.1 4.6 5.1 4.2 3.5 3.5 3.8 5.2 3.8 6.0 4.2 5.4 4.2 4.7 5.9 3.5 3.7 5.5
##  [307] 3.3 4.4 3.9 4.8 4.4 3.9 4.9 3.8 3.9 4.6 3.2 4.3 4.0 6.4 4.6 3.8 4.5 3.9
##  [325] 6.1 3.7 4.5 4.9 4.1 4.4 4.4 3.5 5.5 4.0 5.7 6.1 4.5 5.6 5.1 2.9 5.0 4.2
##  [343] 3.3 4.4 4.8 4.9 4.5 5.5 4.4 5.5 5.9 6.1 5.7 3.6 3.8 4.4 5.3 4.9 4.7 5.3
##  [361] 4.8 2.5 4.0 3.6 5.1 4.8 4.0 5.3 3.7 4.1 4.2 3.9 4.6 5.3 2.6 3.5 4.9 3.6
##  [379] 4.1 4.0 5.7 4.3 3.9 5.6 3.2 4.1 4.2 4.3 4.2 4.7 3.5 4.1 5.0 4.5 4.5 3.6
##  [397] 5.3 3.8 4.8 4.5 4.0 4.8 5.4 4.5 5.5 3.5 4.5 4.7 3.9 4.2 5.1 5.4 3.9 4.9
##  [415] 5.0 4.4 3.4 5.2 5.0 4.0 3.7 2.8 2.2 5.4 5.3 3.6 4.6 5.5 5.1 5.4 6.2 4.4
##  [433] 4.6 4.0 5.4 4.8 4.7 3.4 5.4 4.5 2.8 5.5 3.7 3.0 3.0 3.7 2.4 4.9 3.7 4.4
##  [451] 3.1 4.8 4.6 3.9 3.5 4.1 5.0 6.5 5.2 2.4 4.9 3.8 4.8 2.4 5.6 3.2 6.2 5.1
##  [469] 3.3 4.8 5.6 5.0 5.4 2.5 4.2 5.0 5.1 5.3 5.2 4.1 5.0 4.2 3.8 4.1 4.3 3.0
##  [487] 6.0 4.8 6.1 5.3 5.2 3.9 5.5 4.9 4.0 5.4 4.4 3.1 4.7 4.0 5.9 4.7 3.1 4.4
##  [505] 4.2 5.4 4.0 3.7 5.9 3.4 4.3 4.0 5.0 3.5 5.6 3.7 3.7 4.0 4.1 4.2 4.3 4.8
##  [523] 4.2 5.6 5.3 4.9 2.8 4.6 4.0 5.2 4.0 4.2 3.5 4.0 4.1 4.9 3.3 3.4 4.6 4.4
##  [541] 4.8 4.6 3.0 6.7 5.3 2.8 3.8 6.2 4.5 3.6 5.7 6.8 6.4 4.3 5.4 4.7 4.4 5.4
##  [559] 4.2 4.0 5.7 5.4 5.4 3.5 4.0 4.7 4.1 5.0 4.8 4.4 2.7 3.8 4.4 3.1 3.2 5.0
##  [577] 5.2 4.5 4.5 4.1 5.3 4.9 3.2 4.4 3.9 4.9 5.3 4.4 6.3 3.7 3.2 3.2 4.5 3.9
##  [595] 4.9 5.5 2.9 5.3 4.2 4.5 3.7 3.6 4.7 4.3 3.6 5.9 4.8 6.0 3.1 4.8 3.8 4.9
##  [613] 4.3 5.5 5.8 4.8 3.9 4.7 4.5 6.6 5.1 4.5 3.5 4.4 5.0 5.0 3.8 5.2 6.1 5.0
##  [631] 2.6 2.4 5.1 6.0 4.5 3.1 3.9 5.7 5.0 4.5 4.3 2.9 4.0 2.7 6.2 5.5 3.7 4.1
##  [649] 4.6 5.0 4.9 5.0 6.4 3.0 4.1 5.2 5.1 5.0 4.0 6.2 5.0 3.1 5.7 5.2 4.4 4.9
##  [667] 4.9 5.3 5.6 3.7 3.0 3.9 6.1 4.8 3.3 3.7 4.6 4.2 3.3 3.7 2.9 4.8 4.4 5.5
##  [685] 3.8 4.0 4.8 5.4 4.0 4.2 5.3 3.9 3.1 4.5 5.0 3.2 4.7 4.4 4.0 3.8 6.3 5.4
##  [703] 3.9 4.8 5.8 3.5 3.8 3.4 3.5 4.1 3.6 4.1 5.2 6.2 4.7 6.1 5.1 4.6 4.1 3.9
##  [721] 3.5 4.2 4.0 5.5 5.7 2.7 4.2 5.0 4.1 5.1 4.0 3.7 3.9 4.4 5.7 4.0 5.0 3.3
##  [739] 5.5 3.5 4.2 3.3 3.5 5.4 4.3 4.9 3.8 3.5 4.7 5.3 4.5 5.6 3.9 3.3 4.5 5.2
##  [757] 4.4 5.1 4.1 4.6 3.1 4.8 4.5 3.1 4.1 5.2 4.1 4.1 3.3 5.5 5.0 5.1 4.4 4.6
##  [775] 3.6 3.1 5.4 3.8 3.6 4.4 5.0 3.5 5.3 3.5 3.9 4.9 4.3 4.0 4.1 5.2 3.7 4.9
##  [793] 3.9 5.2 3.5 3.9 5.4 4.2 5.1 5.6 3.7 5.3 3.5 2.9 6.1 3.4 3.6 4.8 4.2 5.7
##  [811] 2.6 6.6 4.5 4.4 6.5 3.3 5.5 4.3 4.8 3.6 5.3 5.2 4.5 3.6 6.5 5.6 3.3 4.2
##  [829] 3.6 3.7 5.7 4.6 4.9 3.7 5.0 2.9 5.0 4.4 2.9 3.1 4.2 4.6 6.0 4.0 3.4 5.1
##  [847] 4.2 6.3 4.1 4.1 6.4 3.7 4.7 4.8 2.8 4.6 3.5 3.3 4.3 4.8 3.9 3.7 3.9 4.1
##  [865] 6.6 5.4 6.3 2.7 2.8 5.1 5.5 7.0 5.3 2.9 5.4 4.9 3.8 5.1 4.5 5.2 3.5 4.7
##  [883] 4.8 5.0 4.8 4.8 4.0 5.6 5.4 5.6 3.9 3.9 5.3 5.5 5.3 5.4 5.0 5.3 5.2 7.5
##  [901] 4.9 4.4 3.4 4.6 4.7 4.5 4.6 4.1 4.9 5.9 4.9 4.8 5.0 4.6 2.5 4.7 5.6 2.8
##  [919] 3.7 3.5 5.4 4.1 5.4 4.3 4.6 4.3 3.7 4.3 5.0 4.1 6.2 4.6 3.3 4.6 4.2 4.0
##  [937] 2.8 4.1 4.5 4.4 5.0 4.0 3.9 4.1 5.0 3.2 3.2 2.5 6.4 4.1 5.0 5.0 5.8 6.4
##  [955] 3.6 3.9 4.8 4.3 4.3 5.5 3.7 4.8 5.2 4.4 4.1 4.3 5.8 3.2 4.8 3.4 2.6 5.1
##  [973] 4.5 4.6 4.3 5.3 3.2 3.2 3.0 5.0 3.8 5.5 5.9 4.4 5.3 5.5 5.0 4.4 4.4 4.9
##  [991] 5.8 4.1 4.1 5.2 5.1 5.1 3.3 4.2 4.8 4.0
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
## 2.7975 6.3000
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
##    [1] 6.1 5.6 4.0 4.7 4.6 4.7 2.8 4.2 3.8 4.0 5.3 3.8 4.6 4.1 3.5 6.3 5.3 3.9
##   [19] 5.0 6.0 5.1 4.2 4.0 4.4 4.6 3.6 4.7 1.9 4.8 2.9 4.1 5.0 2.7 5.8 4.7 6.1
##   [37] 3.8 3.9 3.8 4.0 3.8 6.0 4.2 3.4 4.4 4.1 4.7 3.9 4.6 4.4 4.2 4.6 2.7 5.1
##   [55] 4.5 2.8 3.7 4.0 4.4 3.5 5.0 4.6 4.8 5.9 4.6 3.4 5.1 3.7 4.6 5.1 5.3 4.5
##   [73] 5.1 5.1 3.9 4.8 4.7 4.3 3.8 5.9 5.0 4.5 5.7 4.4 3.6 2.9 3.9 4.8 4.6 4.7
##   [91] 4.6 4.6 3.8 5.6 4.7 4.9 5.0 6.2 4.9 5.7 4.2 5.1 5.6 4.0 2.8 5.2 5.5 6.8
##  [109] 5.2 3.9 4.5 4.4 5.1 4.7 3.9 4.7 4.0 4.3 6.1 4.7 5.3 3.2 3.6 5.1 3.6 3.8
##  [127] 4.7 3.9 2.8 4.5 5.0 5.4 4.6 2.7 5.6 4.3 4.5 3.2 4.0 3.4 5.9 4.2 4.9 5.1
##  [145] 5.1 3.1 3.3 4.3 3.4 4.5 4.4 3.9 4.6 5.2 6.8 5.6 3.9 4.9 3.4 5.4 4.1 5.0
##  [163] 3.3 6.2 3.3 6.4 4.3 2.9 5.0 4.4 3.8 5.0 2.8 4.6 3.4 3.2 5.5 5.0 3.5 4.7
##  [181] 3.1 4.5 3.2 5.7 4.9 6.5 4.5 3.4 3.7 4.7 4.1 6.0 5.6 6.5 2.7 4.4 3.5 3.9
##  [199] 3.4 2.8 5.4 4.3 5.5 3.7 4.5 3.7 4.7 3.1 3.8 4.5 4.3 4.4 4.3 5.2 4.4 5.6
##  [217] 4.8 4.8 3.0 3.7 4.1 3.7 4.7 6.4 5.4 4.9 5.0 5.3 5.3 4.8 4.8 4.7 4.4 5.3
##  [235] 5.1 3.8 4.1 4.6 4.3 4.1 5.6 4.8 3.7 4.7 5.1 3.8 3.9 5.8 4.0 3.6 5.1 5.5
##  [253] 4.9 4.2 6.3 4.8 4.5 3.1 4.9 5.2 4.1 4.1 5.7 5.6 6.2 4.5 4.9 5.7 5.2 3.6
##  [271] 4.2 4.1 5.2 5.0 4.1 5.3 4.0 4.3 5.4 4.6 4.6 4.3 5.5 5.5 4.5 4.7 2.7 3.6
##  [289] 5.2 4.9 3.3 2.6 3.8 5.1 6.1 6.0 5.8 4.1 6.0 4.2 7.2 3.6 3.2 7.8 4.9 4.6
##  [307] 4.4 5.9 4.3 5.1 6.2 3.9 2.3 6.0 5.9 5.7 4.3 4.0 6.9 4.9 4.7 3.8 4.5 4.4
##  [325] 5.3 4.1 3.5 5.7 5.6 5.0 4.7 6.3 3.3 4.1 4.7 5.0 5.4 5.9 3.8 3.8 4.2 4.6
##  [343] 3.5 5.4 4.1 3.8 4.8 3.4 5.1 3.0 5.0 4.1 5.1 4.7 5.8 5.4 3.6 5.7 4.8 3.1
##  [361] 4.2 4.3 4.8 4.1 4.0 4.7 4.6 4.1 3.7 5.0 3.8 4.2 5.0 5.6 5.0 4.9 5.4 4.5
##  [379] 4.7 4.3 3.4 3.4 3.5 3.8 2.9 3.0 5.3 5.9 5.5 5.9 3.9 3.1 3.8 4.4 3.3 5.6
##  [397] 3.8 5.5 5.4 5.0 5.6 2.6 5.5 4.9 4.2 3.5 4.9 4.3 4.4 5.4 4.6 3.9 3.7 4.3
##  [415] 5.0 5.4 3.6 3.4 4.9 3.4 3.9 4.7 5.2 4.3 4.4 5.5 3.3 4.5 6.4 4.5 5.4 3.7
##  [433] 5.6 4.2 4.3 4.6 5.6 4.0 4.5 5.3 3.2 5.0 3.6 4.5 4.1 4.7 4.4 4.4 6.4 4.2
##  [451] 5.4 4.1 4.8 4.8 4.5 3.5 6.3 3.6 4.5 4.3 5.2 5.5 4.2 5.4 5.4 3.4 6.5 4.1
##  [469] 4.5 3.0 4.6 4.2 3.5 4.8 4.5 1.8 3.5 4.9 5.3 2.9 3.7 5.2 4.7 4.7 5.7 3.2
##  [487] 4.9 5.1 5.2 5.6 5.9 7.3 4.3 3.7 3.8 4.6 3.1 4.9 4.3 5.9 3.3 4.9 3.8 7.0
##  [505] 3.4 4.7 2.9 5.3 7.1 4.7 3.6 4.0 4.8 4.1 3.3 3.8 3.9 4.4 3.9 3.3 2.8 4.0
##  [523] 4.7 4.0 5.8 6.0 4.6 4.9 6.0 5.2 4.0 4.4 4.1 6.4 3.1 3.3 5.5 5.5 3.6 6.0
##  [541] 3.0 4.8 4.5 3.3 2.3 5.3 5.0 3.7 5.6 4.5 5.4 3.5 3.5 3.7 2.9 4.6 3.5 5.5
##  [559] 2.9 4.3 5.2 5.1 4.3 5.0 5.4 4.8 4.6 3.1 4.0 5.6 5.0 3.7 3.8 5.7 6.0 4.6
##  [577] 5.1 2.6 4.5 5.7 4.4 4.0 4.7 6.8 6.3 5.2 5.6 3.9 4.2 6.0 3.5 5.4 3.5 3.9
##  [595] 3.6 3.3 5.3 4.4 5.7 5.1 3.0 3.9 3.0 4.6 5.0 6.0 3.8 4.4 6.1 3.5 3.1 4.4
##  [613] 4.7 3.2 3.7 4.7 4.5 3.1 4.4 4.1 4.9 5.6 3.9 4.1 5.0 4.6 5.1 4.1 4.7 4.5
##  [631] 4.0 5.6 5.1 4.1 4.7 4.3 4.4 5.3 2.6 4.6 3.1 4.9 5.2 5.6 3.8 5.5 4.2 4.9
##  [649] 4.7 4.5 4.7 5.8 4.3 3.8 5.6 6.2 3.9 3.7 4.2 4.3 4.5 5.4 4.9 3.5 5.2 4.5
##  [667] 3.1 5.6 5.2 5.4 3.7 3.0 4.9 5.8 5.6 4.2 5.6 5.5 3.1 4.1 4.1 4.4 4.0 3.4
##  [685] 4.7 5.1 3.7 5.4 5.2 4.7 3.1 4.2 4.3 4.6 5.5 3.4 3.9 5.0 4.5 6.4 4.1 4.8
##  [703] 5.2 3.3 2.0 5.4 4.6 5.0 6.0 3.8 4.0 4.1 2.9 5.3 4.6 2.7 4.1 4.4 5.7 5.2
##  [721] 5.6 3.2 4.1 3.7 4.7 5.9 4.6 4.6 5.6 2.7 4.7 4.4 2.7 5.6 5.3 5.5 4.8 4.5
##  [739] 2.3 3.4 5.3 3.5 4.4 5.0 3.7 3.7 3.3 4.7 4.1 4.4 3.0 3.7 3.7 4.7 4.3 4.1
##  [757] 3.9 2.5 3.9 5.6 4.0 5.5 4.0 4.8 4.3 5.3 2.9 5.9 4.5 4.8 5.8 5.7 5.2 5.3
##  [775] 5.1 3.9 3.7 3.8 4.7 5.4 5.0 4.1 4.0 4.1 4.6 4.2 5.1 4.4 3.4 3.4 5.0 3.4
##  [793] 6.6 4.1 5.7 2.8 3.4 6.2 5.7 4.6 4.4 3.2 4.3 3.5 4.6 2.7 4.2 3.3 2.7 3.7
##  [811] 5.1 6.0 4.2 4.0 6.3 4.2 4.6 4.6 5.0 4.0 5.0 4.8 4.7 3.6 5.6 3.9 3.0 4.3
##  [829] 3.9 4.0 3.8 4.8 3.6 6.2 3.4 4.1 4.9 5.3 2.7 5.0 4.6 4.5 4.4 3.7 5.6 3.8
##  [847] 4.4 5.9 4.9 4.8 2.9 3.8 4.7 5.1 4.2 3.0 3.8 5.1 4.4 3.8 4.3 5.9 4.7 4.8
##  [865] 5.2 5.7 5.5 6.0 3.4 4.0 5.9 3.9 3.9 3.5 3.5 5.4 5.4 4.4 3.7 4.7 4.3 3.6
##  [883] 2.6 5.3 4.6 3.8 4.5 6.5 4.8 4.2 4.4 4.5 4.9 4.6 5.9 4.3 3.3 4.3 4.5 5.3
##  [901] 3.7 3.9 5.3 4.1 2.9 4.4 5.3 4.3 3.6 3.9 4.2 4.4 5.7 3.2 5.0 6.1 4.5 4.7
##  [919] 4.6 4.9 4.2 4.8 5.5 4.8 4.1 2.9 3.6 5.8 4.8 4.6 4.0 5.3 5.2 5.2 3.9 5.2
##  [937] 4.3 3.7 5.0 3.5 4.2 4.0 4.5 3.0 4.3 4.3 5.4 5.6 4.9 4.3 4.9 5.6 3.4 4.9
##  [955] 5.2 3.8 6.0 3.6 3.2 5.3 6.6 5.4 4.8 4.3 4.1 3.8 3.0 4.1 4.4 3.7 3.0 4.9
##  [973] 4.0 5.3 5.8 4.7 5.9 3.9 4.8 6.8 4.6 4.0 4.7 3.6 4.2 4.1 3.7 4.0 6.2 4.8
##  [991] 5.3 5.2 5.7 3.7 4.5 3.8 4.6 4.1 4.1 3.5
## 
## $func.thetastar
## [1] 0.0062
## 
## $jack.boot.val
##  [1]  0.50201729  0.38476454  0.32057143  0.15694051  0.07529412 -0.06283988
##  [7] -0.17235294 -0.24808743 -0.41400560 -0.50088235
## 
## $jack.boot.se
## [1] 0.9723169
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
##    [1] 5.3 4.5 4.3 5.3 5.6 4.6 4.0 4.4 4.9 4.3 3.9 2.8 3.9 5.8 5.4 6.0 4.4 5.2
##   [19] 5.2 3.4 5.0 5.2 5.3 4.8 3.5 5.4 3.4 6.1 3.6 3.4 5.2 4.9 3.6 5.8 4.1 6.7
##   [37] 4.7 5.2 2.9 2.9 3.9 3.9 4.7 3.6 3.5 5.0 4.9 5.0 4.1 2.6 4.7 3.3 4.2 4.5
##   [55] 4.1 4.8 6.4 5.3 5.6 4.4 5.7 5.1 4.0 3.8 4.7 4.0 3.5 3.1 4.2 4.2 5.3 4.3
##   [73] 5.4 5.1 4.3 4.2 5.3 5.3 3.6 4.5 3.6 4.6 5.8 3.2 6.0 4.6 5.3 4.2 4.0 5.4
##   [91] 3.8 5.2 2.3 2.9 4.6 2.8 4.3 5.4 4.8 3.0 3.7 3.6 4.5 4.8 3.6 3.4 4.0 5.5
##  [109] 5.7 6.3 3.3 4.2 3.8 5.8 4.4 4.4 4.1 4.0 2.9 5.0 4.0 4.3 5.3 4.4 3.3 4.6
##  [127] 5.6 4.3 4.2 2.4 3.0 2.0 4.5 3.4 5.6 5.6 5.8 4.6 4.2 3.9 4.8 4.6 3.1 4.6
##  [145] 6.0 4.6 4.9 4.8 4.6 6.0 2.8 5.3 4.0 4.1 4.0 5.3 3.3 4.4 4.1 3.8 3.8 5.1
##  [163] 3.4 3.4 5.0 3.8 4.5 3.4 3.7 4.1 4.4 4.3 5.2 5.6 5.1 4.1 4.5 4.9 4.4 3.5
##  [181] 4.8 3.9 5.0 4.7 3.4 4.4 5.6 6.0 5.1 5.1 4.2 4.0 3.5 3.9 5.1 3.8 3.9 3.8
##  [199] 6.4 5.4 4.8 3.2 5.5 5.0 4.4 4.8 4.4 5.9 6.1 3.8 4.8 4.4 5.1 3.9 4.5 3.8
##  [217] 5.9 4.1 3.5 4.4 3.5 3.9 5.5 3.3 4.0 5.9 3.5 4.3 5.2 3.9 3.6 6.3 3.7 4.3
##  [235] 3.0 6.3 3.8 4.7 2.8 5.1 4.8 2.7 3.6 5.9 5.3 4.9 3.7 3.0 5.1 4.8 5.2 4.1
##  [253] 6.3 3.6 4.7 3.1 3.8 2.2 5.0 4.2 4.7 3.8 3.0 3.7 4.6 4.6 3.9 3.9 4.1 4.3
##  [271] 4.3 3.4 3.8 3.6 6.4 3.0 5.2 4.5 3.8 5.0 4.0 4.8 5.4 5.7 4.9 4.6 6.0 5.7
##  [289] 5.7 5.3 4.5 5.0 3.9 3.7 3.5 4.4 4.9 3.7 3.3 3.2 4.4 4.6 3.2 3.3 5.0 5.2
##  [307] 6.4 3.9 5.9 4.9 3.8 5.1 5.8 4.2 3.7 5.4 4.4 5.8 5.4 4.6 3.0 3.6 3.1 5.6
##  [325] 4.6 5.7 5.0 4.9 5.4 4.3 4.0 4.1 4.5 5.5 4.0 4.1 5.3 4.0 4.6 5.6 4.9 5.3
##  [343] 2.6 4.1 3.8 5.0 5.8 4.8 5.0 5.2 3.1 4.3 4.8 5.5 4.9 3.9 4.7 3.4 2.9 3.8
##  [361] 3.6 3.6 2.8 4.1 6.5 4.9 2.9 5.3 4.7 6.7 3.7 5.2 6.3 3.9 3.5 4.1 5.8 4.3
##  [379] 4.4 5.6 6.0 4.6 3.2 3.7 3.1 4.8 5.5 4.9 4.0 6.5 4.1 3.1 3.4 4.7 5.5 4.7
##  [397] 4.2 3.2 6.0 4.5 4.1 3.5 4.1 3.5 5.1 4.5 4.8 5.9 4.1 3.7 3.8 4.1 3.1 4.4
##  [415] 3.7 4.6 4.1 4.4 3.4 5.7 4.0 5.0 4.4 5.5 3.4 5.4 3.9 5.6 5.4 3.6 5.3 4.8
##  [433] 5.1 5.0 5.8 3.2 4.4 5.2 3.2 4.2 3.7 4.5 4.4 5.5 4.0 5.1 4.4 4.0 3.7 5.0
##  [451] 4.2 5.6 3.6 2.9 4.8 5.1 4.4 3.7 3.2 4.8 3.8 3.9 4.5 4.1 3.2 4.1 5.4 4.3
##  [469] 3.9 4.4 5.2 4.9 2.7 4.5 2.6 4.2 3.5 4.0 5.1 4.7 4.2 4.3 2.8 4.6 4.0 4.6
##  [487] 6.0 3.7 4.5 4.9 4.0 4.5 4.3 5.8 3.6 4.5 5.5 2.9 4.5 5.5 4.6 4.8 3.5 3.8
##  [505] 3.6 4.5 3.6 3.8 4.3 5.2 4.5 5.0 5.6 2.9 3.6 3.7 4.7 5.5 4.3 3.2 3.9 3.1
##  [523] 2.6 4.3 4.9 5.5 4.2 5.1 5.4 2.8 3.5 5.4 5.3 6.2 3.2 3.5 3.7 4.4 3.2 4.3
##  [541] 5.1 4.6 5.1 3.5 3.9 3.3 3.0 5.0 6.1 4.0 3.7 4.1 2.9 3.4 4.9 2.9 5.2 3.8
##  [559] 5.3 4.8 4.5 3.5 5.0 4.8 5.1 5.9 5.0 4.9 4.2 3.7 4.8 4.7 4.3 6.3 6.2 5.8
##  [577] 4.4 5.6 5.3 3.3 3.4 5.1 4.6 2.7 5.5 4.3 4.7 4.9 3.6 4.9 4.2 5.3 4.4 2.9
##  [595] 4.5 4.2 5.5 4.8 4.1 5.1 4.8 6.3 3.2 3.1 5.0 4.3 4.7 4.6 3.6 2.1 4.1 3.9
##  [613] 4.5 2.8 4.2 5.8 4.9 5.2 5.1 2.6 3.4 4.5 5.1 3.8 3.6 4.0 4.7 3.9 5.9 3.6
##  [631] 6.6 3.7 3.5 4.0 4.5 6.3 3.7 3.9 4.3 5.5 4.5 4.3 2.4 5.2 4.5 4.7 3.9 3.8
##  [649] 3.7 4.4 3.9 3.5 4.6 3.4 6.2 5.5 4.7 4.0 4.0 5.6 3.1 4.3 4.0 4.3 3.4 3.7
##  [667] 5.2 5.3 3.5 4.8 6.0 3.6 3.5 5.4 3.5 4.5 4.7 4.3 3.4 5.3 5.0 3.6 3.8 4.3
##  [685] 5.6 4.6 4.0 4.7 3.9 4.2 4.7 3.3 4.8 4.4 3.1 4.9 4.3 4.5 4.2 4.6 5.6 4.1
##  [703] 4.7 5.3 3.3 4.5 4.1 5.2 3.4 4.2 4.1 4.7 4.4 2.9 4.8 4.9 6.5 4.5 6.5 3.1
##  [721] 5.5 5.6 4.8 5.4 4.9 5.3 4.4 4.7 3.6 5.3 5.5 4.3 4.2 3.2 4.7 5.4 3.8 5.0
##  [739] 4.9 3.8 3.7 2.7 4.9 5.3 4.3 4.6 5.4 4.0 3.8 4.7 5.0 5.4 4.4 4.2 4.9 4.8
##  [757] 4.4 5.4 3.5 3.2 5.3 4.1 5.6 4.1 2.8 5.0 5.8 5.5 4.6 4.7 5.7 5.0 4.8 5.7
##  [775] 4.1 3.8 5.4 3.5 5.8 4.0 5.1 4.1 4.6 4.7 5.3 3.8 6.6 3.6 4.9 4.1 3.5 5.4
##  [793] 4.3 4.7 4.6 3.8 2.9 5.6 5.4 4.5 4.4 4.1 3.5 3.8 4.5 4.8 4.2 2.5 4.4 5.0
##  [811] 2.9 3.8 5.0 6.5 4.4 4.8 4.6 4.7 3.2 3.7 4.4 4.0 4.1 5.0 4.3 4.2 4.8 5.2
##  [829] 2.6 4.9 6.6 5.0 4.6 6.3 4.2 4.0 3.6 3.4 5.0 4.9 3.9 4.6 2.5 4.0 4.1 4.1
##  [847] 5.4 5.2 4.8 3.8 6.2 3.4 4.9 4.7 2.9 4.5 2.3 5.9 5.5 4.0 6.1 4.8 3.9 5.6
##  [865] 3.6 5.1 3.5 5.1 4.3 4.4 5.5 4.8 4.9 3.6 5.0 3.0 4.2 5.6 6.2 4.5 4.7 5.8
##  [883] 4.9 5.1 6.9 4.5 5.9 6.0 5.1 4.0 4.5 3.4 3.5 5.3 4.9 4.1 4.0 6.0 3.2 5.8
##  [901] 4.7 4.4 4.6 4.1 4.2 4.4 5.4 3.7 4.0 6.5 3.6 4.6 3.6 3.6 4.2 5.6 4.9 4.0
##  [919] 4.9 6.2 2.7 5.4 3.4 4.8 4.8 5.2 4.0 3.8 6.5 4.6 4.9 2.1 3.9 4.2 4.2 5.2
##  [937] 3.1 4.6 3.5 4.5 5.3 5.8 5.4 4.1 4.6 3.7 4.7 4.9 5.1 5.9 4.4 3.3 3.4 4.5
##  [955] 4.0 4.1 4.3 3.8 4.7 5.6 3.7 3.6 6.1 4.9 4.6 3.5 4.3 4.8 5.2 5.3 5.0 4.8
##  [973] 5.0 4.1 5.0 3.3 3.1 4.5 3.8 4.8 2.8 5.2 5.6 5.3 2.7 2.9 3.6 6.2 4.3 5.6
##  [991] 6.2 4.6 5.3 4.5 5.1 4.6 4.9 4.3 3.5 2.8
## 
## $func.thetastar
## 72% 
##   5 
## 
## $jack.boot.val
##  [1] 5.500 5.300 5.300 5.212 5.100 4.900 4.900 4.700 4.500 4.400
## 
## $jack.boot.se
## [1] 1.039841
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
## [1] 0.9532624
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
##   3.493847   6.157406 
##  (1.493912) (2.831399)
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
## [1] 2.2684127 0.4555944 0.7727185 1.7698869 0.3032580 1.1627536
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
##    [1] -0.3649279288  0.9487727306 -0.4982573593  0.9581659818  0.9652086803
##    [6]  0.9080561198  0.7046178191  0.8564485329  0.1943736916  1.0816483173
##   [11]  0.8493330368  0.0072876324  0.4119503276  1.1438553106  0.6924260627
##   [16]  1.5220156965  0.8605164609  0.1241343087  0.3834477214  1.0128311564
##   [21]  1.0324791491  0.9339708450  1.0425294641  1.1912850185  1.1392605624
##   [26]  0.1451143209  1.3325573060  1.1326468734  1.0823856968 -0.4658282587
##   [31]  0.1092408529  0.0940982670  0.5359117328  0.7618106935 -0.0124781352
##   [36]  0.4974620156  0.7154663698  0.3863281185  0.4131288896  1.3210911914
##   [41]  0.7903610679  0.9283458527  1.2061614700 -0.3878534699  0.6087136596
##   [46]  1.9954903065  0.8775808051 -0.1586853171  0.4090686736  1.4931553122
##   [51]  0.9293851254  0.4139620407  0.9540624287  1.0586478031  0.7635430309
##   [56]  0.8895175509  0.9658668289  0.8437349343  0.6304461374  0.5786454250
##   [61] -0.0334081039  1.8680212490  1.5372112235  0.2795932944  0.7237483613
##   [66]  0.5627045108  0.5562349459  0.4071024375  0.2482922553  0.6542616558
##   [71]  0.8906824496  0.3494767069  1.3221958751  0.5111648027  1.4903540999
##   [76]  1.4079783832 -0.3621360994  0.7168232894  0.7392634673 -0.0827626433
##   [81]  0.5098243008  0.6211698084  1.3523369124  0.8496178824  0.5421680215
##   [86] -1.1759007045  1.9208836761  1.2207428095  0.1448829185  1.2324651047
##   [91]  0.1695198880  0.3375821623  0.9290064511  0.9277305633  1.2844546934
##   [96]  0.8798969801  1.3094777642  0.2961693984 -0.6668299338 -0.0098601810
##  [101]  1.4785098621  0.4646546678  1.0446831685  0.5071170287  1.1958717780
##  [106]  0.7083231969  0.6888073187 -0.1998658729  0.1289537848  0.1677947036
##  [111]  0.7262925480  0.7253020719  1.5324852979  0.9410225286  0.8932649122
##  [116]  0.4814796580  2.0753037047  1.7842013464  1.3342518727  0.0596164644
##  [121]  0.9028313290  0.2267100039  0.0490520727  0.6137742250  2.0435225779
##  [126]  0.5304775958  1.8109516832  1.3910366267  1.4694953616  0.4539586382
##  [131]  1.1918769132  0.4474056694  0.8092993319  1.1455437461  0.4878931604
##  [136]  1.3545458439  0.3590503973  1.0667214764  0.2764679443  1.0981770114
##  [141]  1.0914664315  0.7679737375  1.8652320873  0.6227968159  0.3207660445
##  [146]  0.5546886449  1.1087351962  1.3508320293  0.4178330273  0.7521856090
##  [151]  0.9200026866  0.9942592178  1.7417696431  1.3565513991  0.7667008759
##  [156] -0.1701741375  1.4196204268  0.6786628293 -0.0155791535  1.1234076945
##  [161] -0.0069787281  0.6161429756 -0.1835039957  0.9553209721  0.9670890278
##  [166]  0.9042786954  0.0140472513  0.1264113955 -0.2432393304  0.7655470479
##  [171]  0.9127585743  0.4017976558  0.8011018546  1.1900275932  1.2343330700
##  [176]  0.2188663990  0.7587679390  1.2370801244  1.0298696985  0.4569386808
##  [181]  1.6977076669  0.3888595241  0.2193146878  0.7512843320  0.1302284116
##  [186]  0.7329136914  1.4109540156  1.2461453036  0.7625437437  0.8766052169
##  [191]  1.7965475605  1.1070647623 -1.1160439577 -0.6810674795  0.6382767326
##  [196]  1.0565671524  0.3067958712  0.8169325048  0.9842133565  0.2645269291
##  [201]  1.0351806615  0.4592539287  0.9403336651  1.0861254995  1.2876698820
##  [206]  0.4657652855  0.2410877089  0.7548885743 -0.0049984382  1.6999405062
##  [211] -0.0173141966  0.1164325140 -0.0415212207  0.0267810447  0.3842978243
##  [216]  0.3413371411  0.2964338864  0.4445131866  0.5034102061  1.7973484879
##  [221]  1.1593884127  0.6177313533  0.4620873255  0.3919861439  0.4277272495
##  [226]  1.1979845452  0.2100060241  1.0083609867  1.0310299440  0.8546437044
##  [231] -0.1511201902  1.1234663162  0.4471699018  0.8691311251  0.9648099386
##  [236]  1.6234065552  0.9960588325 -0.2320393886 -0.0435902990  0.9886537489
##  [241]  0.6268915952  0.6467845611  0.9888285344  0.4151559545  0.0109073905
##  [246]  0.1371558281 -0.0784542803  0.8373813352  0.4993129674  0.4641099630
##  [251]  1.6515919262  0.8800616723  0.4046190236  0.1740508066  0.6991507405
##  [256]  1.5093266783  0.6412075084  0.2594571623  0.9307954234  1.0772716061
##  [261]  0.6228369628  0.4383869448  0.7260040922  0.1468262463  0.9940926610
##  [266]  0.4984389381  1.2381217208  0.1728948828  0.6256322819  0.9801922383
##  [271]  1.0695465622  1.4917882416  0.2411135495  0.4083354376 -0.1111987206
##  [276]  1.0656248097  1.3059142599  0.1585129835 -0.1876441708  0.6818524492
##  [281]  1.1011027172  0.7326466669  1.0546660178 -0.1652613245  1.5204149058
##  [286] -0.1573018025  0.4526378759  1.2412872808  1.0615328822  0.2367804473
##  [291]  0.3243517793  0.4027079358  0.0380603701  0.6583877197  0.5416974532
##  [296]  0.1132245995  1.2877569060  0.0955001297  0.5505736949 -0.1829568746
##  [301]  0.8423674111  0.1341951266  1.2333946427  0.2698965772  0.2711389807
##  [306]  0.7253065227  0.1962239595  0.1448829185  0.3061721160 -0.2026821383
##  [311]  0.8292608492  0.7786505498 -0.2371285386  0.7155525316  0.6416823116
##  [316]  1.5693142465  0.8275134412  1.5703627437  2.4239124849  0.5066903730
##  [321] -0.2194036760  1.2080311836  0.3867765215  0.1971803796  0.6957414777
##  [326]  1.0121405770  0.0762972654  0.5694036415  1.4901595867 -0.1042477922
##  [331] -0.1716136469  0.8275134412  0.6286450242  1.0117035211  1.6789601472
##  [336] -0.1063362436  0.2534190207  0.9479804471  1.0920081425  1.3678276117
##  [341]  0.4294626356  1.1616881002  1.1667835062  0.8933502264  1.7386326085
##  [346]  2.1360761339  0.1470840149  1.3678371608  0.3523893089  1.5676624264
##  [351]  1.1100692614  2.3803694800  0.8624175064  0.7676973642  0.5346308843
##  [356]  0.3666187739  0.3162769918  0.9582550318  1.1307222912  0.5531797031
##  [361]  0.0419548783  0.5503514345  1.1260504912 -0.0415463685  1.1996328177
##  [366]  0.0261777481 -0.0387099797  1.0209411869  1.0520391791  0.0872135011
##  [371]  0.3224609965 -0.4155903671  0.5634575145  0.1842874403  0.9713187846
##  [376]  1.0986725811  1.5533424068  0.7525314712  0.0610326086  1.0340054677
##  [381]  0.6396613336  0.3362865344  1.4410540441  0.9012835874  0.7487760804
##  [386] -0.2876876237 -0.6349288640  0.9664844089  0.2782505268 -0.2401201040
##  [391]  0.0083004258  0.9916147243 -0.0997089292  0.5839888736  1.5425196452
##  [396]  0.8180405631  1.1902784684  0.3554310565  0.8364634924  1.0427313535
##  [401]  0.9396139671  0.1054231941  0.2440476500  0.8825949416  0.8515747327
##  [406]  1.0422456535  0.7223348221  1.1796387806  0.5888034039  1.2008847860
##  [411]  1.3556497198  1.6553733710  0.5597057385  0.5247465582  0.7784616347
##  [416]  0.7337232442  1.0249726345  0.9518618324  0.1803038540 -0.1434211310
##  [421]  0.3394585196  0.6597224182  0.0217085348  0.6606091952  1.2463741258
##  [426]  0.0881593041  1.9312732662  0.7693541512  2.2421226994  0.6832075608
##  [431]  0.0453410153  1.1611581379  1.1261004639 -0.0445840852  0.9284013123
##  [436] -0.0660961029  0.5285475986  0.4442212063  1.8083031064  0.4996884591
##  [441]  1.7064824248  0.8712019592  0.6477292009  0.3387861041  1.6083606864
##  [446]  0.5138515965  0.7780768822  1.7263663169 -0.9114767460 -0.5235188858
##  [451]  0.9411222689  0.5464491182  1.1679065561  1.8969610353  0.4342816148
##  [456]  0.7692495549  0.0442107942  1.5828116404  0.5268453574  1.1940914214
##  [461] -0.0187128569  0.4535456450  0.6630116599  0.8835992226 -0.4042576352
##  [466] -0.1216270508  0.1748479871  1.2987234715  1.4424857238  1.1631075758
##  [471]  1.6683587553  0.4918476363  1.5703839609  0.2767447378  1.4239605441
##  [476]  0.1029803994 -0.1076474718  0.2892524646  1.5610995519  0.9422123906
##  [481]  0.4754676279  0.1271841846  0.5012257420  1.4222131390  1.6523731954
##  [486] -0.0143578606  1.5895900433  1.7573510154  0.3028844384 -0.8240948638
##  [491] -0.2601193897  1.4171685104  0.8773241552  0.8343715159  0.9155626638
##  [496]  0.9984063673 -0.1535733790  0.1566203988  1.1052274263  0.3939962323
##  [501]  0.8879164867  0.5627015082  0.8472735593  0.6497569997  0.3833406370
##  [506]  0.3928679562  0.0632258823  0.7189966419 -0.4554983067 -0.1161131793
##  [511] -0.7873854375  0.8872545577  0.5428477403  0.8710298453  0.8153735312
##  [516]  0.6941673315  1.0389296985  0.6524883370 -0.1398897221 -0.7400943255
##  [521]  0.1418543116  1.1051360958  1.7114537779  0.8586691573  0.7830665475
##  [526]  0.8406438490  1.0095794947  0.7470799970  0.3862297873  1.6356979081
##  [531]  0.4842652261  0.6249419881  0.6640477209  0.8683095670  0.5084662185
##  [536]  0.7216417095  0.4627900879  1.7583872458  0.6184364924  0.2937164127
##  [541]  1.0585384161  0.9257014320  0.7311733268  1.8597678722  1.2631195098
##  [546]  0.5743334348  0.8706547853  1.2151026724 -0.6833037289  0.4140031643
##  [551]  0.0361709802  1.0509305682  1.0491620556  0.8869996797 -0.0352169238
##  [556]  0.6097406770  0.6453901717  0.3224609965  0.5198133601  2.1637039830
##  [561] -0.2706720659  0.2714084414  0.8103352052  0.5709422457  0.8215024476
##  [566]  1.4951867571 -0.2662994868  1.0844576828 -0.4234834835  0.9566987130
##  [571]  0.6528780302  0.5524104202  0.2576504066  1.5449207963  1.0970635236
##  [576]  0.7447216373 -0.3675146377  0.9150916051  0.1667912667  0.0816124800
##  [581]  0.4718210988 -0.3736503589  1.6129214269  0.4510057238  0.1894646572
##  [586]  1.3259795685  0.7049525037  1.2334855082  1.4955248362  0.1948206530
##  [591]  0.7134234676  0.7272429754  0.3141649043  0.5211370032  0.1587757449
##  [596]  1.1157662833  0.6322999237  0.9073508231  1.0525565944  1.1501073340
##  [601]  1.6753520812  0.1981683359  0.4628033814  1.4978067512  0.1261854681
##  [606]  0.6569368946  1.2720021341  1.3239941181  0.8955787625  0.6689054123
##  [611]  0.2895544448  0.1204305063  0.4103530359  1.8182050137  0.6305806326
##  [616]  0.9055711191  1.3880891542  1.0713068605 -0.3576661420  0.3607484572
##  [621]  0.6696788785  1.7953672772  0.3874866608  0.7356655530  0.5029225190
##  [626]  1.2337996237  0.4194282616  0.9016299906  0.3291580452  0.5355309780
##  [631]  0.0562757629  0.3248582185  1.0066457749 -0.2402035581  0.4067351606
##  [636]  0.4678792006  0.8477384001  0.0874501201 -0.5928230169  0.7385222485
##  [641]  1.0722009996  0.6414170280  0.0460680379  0.9226272878  1.1534411061
##  [646]  0.5065380676  1.3724918158  1.0480095041  1.3544700134  1.0093807051
##  [651]  0.7183551152  1.3210911914  1.3224667358  1.2286852881  1.0004580301
##  [656]  1.0425294641  0.4686407540  1.8617729840  1.1430737835  1.2012207943
##  [661]  0.6728979717  0.2973082808  0.9248143666  0.4460601520  1.1674855166
##  [666]  0.9586009396  0.9074276638  0.2958954393  0.7782778142  0.0236857675
##  [671]  0.6884769575  0.5689669215 -0.0274811962 -0.1092530416  1.0835759691
##  [676]  0.6763369599  0.5450214677  0.2712292835  0.2411358714  0.5867238528
##  [681]  1.2284409389  0.9185657203  0.9836139622 -0.0274836574  0.9331976394
##  [686]  0.6738390958  0.8592254792  1.1083830487  0.7358056373  1.0943789494
##  [691]  1.1467276998  0.1985227616  0.9910830486  0.0023747717  1.2470462704
##  [696]  1.0972441459 -0.9019864353 -0.0115663836  0.8243431681  1.1475179027
##  [701]  0.1419773315  0.3691194481  0.4441866442  0.4929578036  0.3707993569
##  [706]  0.8094469208  1.1083830487  0.6529510315  1.2984828564  0.0400468412
##  [711]  1.6238193533  1.6028524927  1.3834042969  1.0551236644  0.0518562779
##  [716]  1.9902238631  0.7493009020  0.5743334348  0.8158320429  0.7523010315
##  [721]  0.7707078518  0.1631972906  0.3395445323  1.5071304485  0.3111019333
##  [726]  0.9364545955  0.3440249449  1.5734635226  0.0305131186 -0.5395485395
##  [731]  1.6003027945  0.5792125743  0.3804695547 -0.4254037488 -0.4132439864
##  [736]  1.1070160176  0.8697811156  2.4572700830  0.8978297356  0.7302200359
##  [741] -0.0814831324  0.8181099923  0.1109049464  0.2829284929  0.9832106735
##  [746]  0.4581093350  1.0833965190  1.2305677573  1.0230909752  0.6249419881
##  [751]  0.3663712647  0.0981318858  0.7645780450  1.1511265508  0.7371651542
##  [756]  0.3742822963 -0.0049674449  0.5001910569 -0.2104917695 -0.2180488327
##  [761]  0.8032617768  0.5263275756  0.4338260325  1.0261084352  0.8249074022
##  [766]  0.9507500141 -0.2921638698  0.6156931201  0.1360498995 -0.0468423678
##  [771]  0.8671174702  0.9254145148  1.0352124334  0.6887409895  0.5296407819
##  [776]  0.6178070049  1.5023455268  1.0359121597  0.2902803677  1.5552020020
##  [781]  0.1569174286  0.8720059289  0.6967396656  0.3152655392  1.5450087657
##  [786]  0.6059937608  0.3864840739  0.8280439878  0.5440819735  0.0764057058
##  [791]  1.3351881421  1.9029065213  1.0371462066  0.7461294502  0.7163058376
##  [796]  0.6681607650  0.5691171141  0.5729584271 -0.3612074723  0.8809399968
##  [801]  0.3836545299  0.5692823058  0.9113577076  0.6029889913  0.6033047993
##  [806]  1.5870664260  1.2984890400  0.4580713133  0.5064240392  0.7859997796
##  [811]  0.1592722607  0.4567793599  1.6769010996 -0.2488387871  0.3945730503
##  [816]  0.6208396193  0.6035514725  0.8499904585  0.8873995594  0.9690229287
##  [821]  1.5637211510  1.1044635751  0.6076128712 -0.0365519104  1.4953426228
##  [826]  0.4225467889  0.1009098032  0.1645753368  0.7951105535  0.5325781618
##  [831] -0.1067177411  0.6069750968  1.2170695424  0.1740508066  1.4323937424
##  [836]  0.2762306569  1.5068577479  0.2461569072  0.5476438054  1.5329904524
##  [841] -0.0744041805 -0.1288738541  0.6340755207  0.7453883317  0.2267566064
##  [846]  0.8462518861  1.6776930406  0.5300110535  1.2153073241  0.5380865881
##  [851]  0.5826916686  0.8806636593  0.5498166613  0.1556033498  1.1051243350
##  [856]  1.0878342963  1.1011027172  0.2918855083  1.0011869756  0.6051114232
##  [861]  0.7483100574  0.1575319801  0.8585579051  0.8946099483 -0.1890961680
##  [866]  0.5627015082  1.1046835333  0.9878061593  0.9422013136  1.1145012889
##  [871]  0.9488849063  0.8267011834 -0.0130790597 -0.1165856902  1.2781260611
##  [876]  1.7649495154  0.3312975225  0.7602220088  0.5176595850  0.3977122211
##  [881]  0.4669046316  1.0610667104  0.6465636859  0.3948635347  0.2728419230
##  [886]  0.7944939202 -0.1486085808 -0.2229120783  1.1543524851  0.9892697232
##  [891]  0.6053863398 -0.6052028314  0.3417623582  1.0197977262  0.6127016190
##  [896]  1.6064154707  0.4568600719  0.0439982019  1.3259795685  1.3949032950
##  [901]  1.0261084352  0.1353585482  0.3705044885  0.0684670602  1.0002379258
##  [906]  1.9300980938  0.4398703534  0.1790104898  0.9807578196  0.0090469362
##  [911]  0.6316871472  1.0641554948  0.5852621422  0.2142802996  1.6791104302
##  [916]  0.2709178855 -0.0816704424  0.3703860979  0.0316181182  0.9582115527
##  [921] -0.1308652329  0.7816753754 -0.4942286361  0.5350896365 -0.2989366350
##  [926]  0.1202955163  1.0607905869  0.5930674331 -0.2091948101  0.5594038632
##  [931]  0.3041807977  0.3717594881  0.7368811662  0.5771609445  0.1793511463
##  [936]  0.7652728819 -0.1728258631  0.0130882611  0.6053562215  0.4862392636
##  [941]  0.2464418009 -0.3581429407  0.4686622009  0.9884225530 -0.2695970127
##  [946] -0.1331758464  0.6778034223 -0.7202611096  1.0914867991  0.5538028099
##  [951]  0.7643163793  1.3681402740  1.0005560281  0.3376127629 -0.1960661256
##  [956]  0.1545325911 -0.1111987206  1.8111447236  0.9312471161 -0.4341369448
##  [961]  0.7128809877  0.7067350232 -0.5522053904  0.7785779189  1.0102390788
##  [966]  1.2563132617  1.4583447268  0.3844344249  2.2162360552  0.3941183326
##  [971]  0.8762774737  0.0907967699 -0.5693814333 -0.0001899499  1.4869225008
##  [976]  0.2284030254  0.1327006351  1.0005800715 -0.3042758647  1.1200748503
##  [981]  0.0834930910  0.5041786956  0.0967565868 -0.1851835700  0.5365268270
##  [986]  1.2084073752  0.5115392074  0.5009532056  0.0141333368  0.6816569395
##  [991]  0.4146757104  1.1329974282  1.4112955775  0.4344931558  0.4855708348
##  [996]  0.4141051366  0.9497032885  0.2231135869  0.6759750987  0.5238574827
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
##   0.56742183   0.30990597 
##  (0.09800087) (0.06929473)
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
## [1]  0.40560076  0.52562250 -0.95059730  0.14320770 -0.20386863 -0.07649019
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
## [1] 0.0328
```

```r
sd(results2$thetastar)  #the standard deviation of the theta stars is the SE of the statistic (in this case, the mean)
```

```
## [1] 0.9267702
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
## t1*      4.5 -0.009409409   0.9171226
```

One of the main advantages to the 'boot' package over the 'bootstrap' package is the nicer formatting of the output.

Going back to our original code, lets see how we could reproduce all of these numbers:


```r
table(sample(x,size=length(x),replace=T))
```

```
## 
## 0 3 5 7 8 9 
## 4 1 2 1 1 1
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
## [1] -0.0142
```

```r
se.boot
```

```
## [1] 0.9080875
```

Why do our numbers not agree exactly with those of the boot package? This is because our estimates of bias and standard error are just estimates, and they carry with them their own uncertainties. That is one of the reasons we might bother doing jackknife-after-bootstrap.

The 'boot' package has a LOT of functionality. If we have time, we will come back to some of these more complex functions later in the semester as we cover topics like regression and glm.

