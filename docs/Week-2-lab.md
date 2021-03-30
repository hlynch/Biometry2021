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
## 1 2 3 4 5 6 7 9 
## 1 1 1 2 1 1 2 1
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
## [1] 0.0095
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
## [1] 2.739744
```

```r
UL.boot
```

```
## [1] 6.279256
```

Method #2: Simply take the quantiles of the bootstrap statistics


```r
quantile(xmeans,c(0.025,0.975))
```

```
##  2.5% 97.5% 
##   2.8   6.2
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
##    [1] 3.4 4.9 3.1 5.8 5.2 6.2 4.2 2.3 4.1 3.7 4.2 2.7 5.3 4.2 4.1 5.2 4.4 5.0
##   [19] 3.9 5.3 4.1 3.4 5.3 4.0 3.8 4.2 5.3 3.6 2.8 4.3 5.1 5.1 5.1 3.9 4.2 4.1
##   [37] 6.2 4.3 3.7 3.4 4.2 5.8 4.6 3.9 3.9 3.5 5.1 4.2 3.1 4.6 5.2 5.6 4.2 5.4
##   [55] 4.1 4.2 3.5 3.7 4.6 3.3 4.1 5.1 5.5 4.3 4.5 4.7 4.8 5.5 6.0 6.2 4.8 1.4
##   [73] 4.2 5.5 4.2 5.1 6.0 4.0 3.5 3.5 4.5 4.8 4.4 4.2 3.0 3.8 3.9 3.4 4.9 4.1
##   [91] 5.0 3.2 2.5 5.0 3.6 3.5 5.2 4.0 5.1 5.5 4.9 5.1 5.3 5.1 4.4 6.2 5.4 4.9
##  [109] 3.9 5.4 5.1 3.4 5.2 3.7 4.4 5.3 5.3 3.6 4.5 5.4 2.9 3.4 5.0 3.5 2.8 4.8
##  [127] 6.3 5.5 5.2 4.3 4.2 4.1 3.7 4.2 3.4 3.7 3.3 4.6 5.8 3.3 3.9 5.1 3.1 5.0
##  [145] 4.3 4.7 5.0 4.9 3.5 4.2 4.5 5.8 5.2 5.1 5.9 5.0 4.9 4.0 3.3 4.4 3.8 5.0
##  [163] 4.2 4.4 2.4 4.5 2.9 4.8 3.3 4.7 5.7 4.7 2.9 4.8 5.4 4.7 5.0 6.2 5.5 3.5
##  [181] 4.5 5.5 5.2 5.2 4.0 5.7 4.7 3.1 4.9 5.9 3.9 4.7 4.1 5.4 3.7 4.0 5.7 4.9
##  [199] 4.7 4.7 4.1 4.6 5.3 3.9 5.1 6.2 6.4 4.9 5.5 3.3 4.8 4.8 3.8 3.9 3.9 5.0
##  [217] 5.7 3.2 3.1 5.1 3.3 3.9 4.1 5.1 4.8 3.1 5.9 4.5 2.6 2.9 5.4 4.9 3.7 3.6
##  [235] 4.7 2.6 4.1 5.0 4.5 3.6 5.3 5.3 5.4 4.1 4.8 3.5 3.6 6.2 5.4 3.7 4.6 5.4
##  [253] 2.3 5.5 5.6 4.7 6.6 5.0 3.2 5.9 4.8 4.7 5.1 5.0 3.8 2.7 3.9 4.0 5.1 4.2
##  [271] 5.2 4.3 2.9 5.1 5.0 5.4 4.8 4.2 6.1 4.5 4.4 3.4 4.3 3.6 4.5 5.6 4.8 5.1
##  [289] 4.3 5.0 6.1 4.6 4.8 4.0 5.7 4.3 4.1 6.9 5.5 4.8 4.1 6.1 4.8 4.4 4.7 3.1
##  [307] 4.7 4.6 5.5 3.9 5.4 4.7 4.9 3.9 5.3 5.1 3.6 2.7 3.9 2.4 3.4 5.1 5.1 5.1
##  [325] 6.1 5.6 2.3 3.3 3.4 2.9 6.1 4.7 3.2 4.0 5.6 4.9 3.4 4.5 5.3 3.2 3.9 3.9
##  [343] 4.1 4.7 3.7 5.4 3.7 4.9 4.5 4.1 6.0 4.7 4.5 5.9 3.9 5.8 4.7 3.8 4.3 4.1
##  [361] 6.5 3.6 4.2 5.7 4.3 5.9 3.0 5.3 3.7 4.6 3.5 5.3 4.6 3.9 4.3 2.6 4.1 4.4
##  [379] 5.3 4.0 4.7 2.5 4.5 5.0 5.1 5.2 2.7 4.7 4.5 5.0 4.5 4.6 3.7 4.5 5.2 4.9
##  [397] 3.9 2.9 4.6 3.7 4.7 4.7 4.9 4.3 4.4 3.6 5.6 5.8 4.1 4.6 4.6 6.4 6.3 4.9
##  [415] 5.1 3.5 6.1 3.8 3.8 6.2 3.8 5.2 4.9 4.4 5.2 4.6 3.7 4.7 4.8 6.0 4.8 5.8
##  [433] 2.6 4.8 5.0 3.8 4.7 4.1 4.1 3.8 4.4 3.8 3.1 5.0 4.8 4.1 4.3 5.9 3.9 4.4
##  [451] 3.2 5.9 3.9 4.2 4.6 5.4 4.8 4.2 3.0 4.2 5.3 5.7 4.2 4.6 5.3 3.6 5.1 3.6
##  [469] 5.4 5.3 5.5 6.0 4.1 5.3 5.1 4.6 4.8 4.5 4.1 6.0 3.9 3.8 3.0 5.6 4.6 4.7
##  [487] 4.1 4.3 5.8 6.2 4.7 4.9 4.1 4.6 5.8 3.9 5.6 3.8 4.9 6.0 5.8 5.4 3.9 6.5
##  [505] 2.8 5.6 3.4 5.7 4.7 4.0 3.9 3.6 4.1 4.8 4.9 3.8 5.6 4.4 3.4 5.5 2.9 4.3
##  [523] 5.7 4.0 4.5 4.1 5.0 4.3 5.3 5.4 4.0 4.5 4.9 4.9 5.1 4.6 4.2 4.8 3.6 3.6
##  [541] 3.0 5.5 3.6 4.3 5.5 4.6 3.0 5.6 4.0 4.2 4.6 6.0 5.6 5.1 4.6 2.8 3.8 4.5
##  [559] 3.9 4.8 6.2 4.2 4.8 3.7 3.7 3.9 4.2 3.0 5.1 4.2 4.5 5.1 4.6 5.1 3.0 4.7
##  [577] 4.4 4.3 4.8 5.6 3.7 4.6 3.5 5.0 5.6 5.2 3.6 3.8 4.9 2.5 2.8 4.6 5.0 2.9
##  [595] 4.4 3.3 4.6 4.3 4.2 3.9 4.5 3.9 6.0 4.2 3.6 3.8 4.2 4.3 4.1 6.0 4.8 4.0
##  [613] 5.0 3.8 5.0 6.0 4.8 5.3 2.7 4.5 3.9 3.1 2.7 5.1 4.5 5.0 4.0 4.0 4.7 5.6
##  [631] 5.0 4.3 4.9 4.3 4.0 2.5 3.6 4.4 4.3 4.7 6.7 5.1 4.7 5.0 4.3 4.1 4.3 2.7
##  [649] 4.4 4.0 3.3 6.1 4.2 2.9 5.2 2.5 4.8 5.1 2.9 3.2 5.0 4.2 6.2 4.4 3.3 5.0
##  [667] 5.0 4.1 5.5 3.7 6.0 4.8 3.7 4.3 4.8 5.7 4.4 4.2 4.7 4.0 3.7 2.9 3.7 3.9
##  [685] 5.4 5.9 4.9 3.3 5.7 2.9 4.0 5.6 2.6 6.2 3.6 4.0 3.3 4.1 6.4 4.5 5.2 4.5
##  [703] 5.2 5.3 4.8 3.8 4.2 5.2 4.1 4.1 6.1 4.3 2.9 3.2 5.1 4.3 3.6 1.8 4.0 3.3
##  [721] 4.4 3.7 4.2 2.9 3.0 4.5 3.2 3.2 2.4 3.2 6.3 4.7 4.3 2.5 5.3 4.6 4.4 2.8
##  [739] 3.8 4.0 6.0 4.4 5.3 4.4 6.1 4.6 6.0 3.2 4.5 3.9 4.2 3.1 4.2 5.9 4.2 5.3
##  [757] 3.7 5.6 4.2 2.9 5.0 3.6 5.0 5.2 3.6 3.8 3.4 4.1 5.0 4.2 4.5 5.3 4.3 3.8
##  [775] 4.4 5.4 4.1 3.5 3.4 4.9 5.1 3.9 5.3 4.3 4.8 4.7 3.4 5.8 3.9 4.0 5.5 5.1
##  [793] 4.3 4.0 3.1 4.1 5.7 4.6 4.4 4.6 4.4 5.3 4.9 4.1 3.5 6.4 5.6 4.9 4.9 4.0
##  [811] 4.0 3.8 4.6 4.4 5.2 5.5 5.1 5.8 4.4 6.2 3.6 4.6 3.4 4.1 2.8 5.5 3.6 4.4
##  [829] 3.5 4.4 4.1 5.4 3.3 5.4 5.6 5.9 4.2 4.5 5.1 3.7 5.2 3.7 2.2 5.0 5.4 4.8
##  [847] 4.2 4.1 3.5 6.7 4.9 4.0 6.0 3.9 4.6 6.3 6.1 3.9 4.6 5.3 5.3 6.4 5.0 5.0
##  [865] 5.0 5.0 4.4 5.5 5.3 3.9 5.1 4.1 4.9 5.3 4.3 3.4 4.6 4.2 5.9 4.1 3.7 6.3
##  [883] 4.0 6.3 3.8 4.4 3.3 4.9 4.2 5.2 4.4 5.3 5.3 5.6 4.3 4.0 5.1 4.9 4.2 5.0
##  [901] 3.4 3.3 5.1 4.0 3.8 5.6 3.8 5.2 3.3 2.8 4.6 4.3 2.9 4.6 3.6 5.8 3.9 3.3
##  [919] 3.9 4.7 4.3 4.6 5.6 4.2 5.4 4.2 5.2 3.8 4.2 5.2 5.4 4.5 3.4 2.6 3.5 3.6
##  [937] 5.0 3.9 5.0 3.7 5.1 4.6 4.2 3.9 4.7 4.9 4.8 4.7 2.9 4.4 5.0 4.1 5.4 3.5
##  [955] 4.9 5.0 5.6 3.7 3.7 2.9 5.0 3.5 3.4 4.1 2.9 4.4 3.9 4.9 3.3 4.5 3.3 4.1
##  [973] 4.1 4.9 6.2 5.3 4.1 5.4 4.3 5.5 4.5 4.0 6.0 3.8 2.7 4.2 4.0 3.4 5.9 4.1
##  [991] 3.3 3.3 6.0 4.9 4.0 3.9 2.9 3.8 4.6 6.2
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
##   2.7   6.2
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
##    [1] 4.7 4.4 3.6 4.6 4.0 3.5 4.8 5.1 3.9 6.2 3.3 4.5 3.3 5.4 3.9 6.7 4.8 4.1
##   [19] 5.2 5.8 5.3 5.9 4.9 4.8 5.2 2.8 4.8 3.2 4.2 4.2 3.2 2.6 3.4 4.4 4.3 3.7
##   [37] 4.2 4.9 5.4 4.9 3.9 4.7 3.7 4.9 4.3 3.5 3.6 4.1 4.9 3.0 3.2 5.0 3.0 4.8
##   [55] 5.4 5.1 3.7 5.4 3.7 4.4 6.0 5.0 4.6 3.3 3.9 3.4 5.4 3.9 4.6 4.4 5.1 4.9
##   [73] 4.3 6.2 3.8 4.5 3.6 5.0 4.6 4.3 4.4 5.3 4.7 5.2 4.1 4.3 2.8 4.8 6.3 5.5
##   [91] 4.5 5.5 6.6 6.0 3.8 4.4 4.6 3.9 4.3 3.9 3.0 5.0 4.5 3.8 3.7 4.4 3.8 4.5
##  [109] 5.9 3.1 5.2 5.7 4.8 4.8 5.4 4.0 5.7 5.1 6.1 5.5 3.0 5.3 3.0 4.5 5.0 3.8
##  [127] 4.2 6.1 3.4 4.4 6.0 4.5 5.0 4.6 3.5 6.0 6.0 3.1 5.5 4.2 5.0 2.6 5.6 4.0
##  [145] 6.4 3.8 5.7 4.1 3.9 6.0 4.0 4.7 4.5 4.3 3.2 3.4 4.5 3.5 4.8 5.8 4.6 3.9
##  [163] 3.4 5.0 5.0 5.3 4.8 5.1 4.5 2.7 4.4 5.1 2.3 4.6 4.3 4.4 4.2 4.5 4.8 3.9
##  [181] 4.1 4.0 4.6 5.8 3.8 6.0 3.3 5.1 4.3 4.8 4.2 4.9 3.4 3.7 3.5 3.5 5.6 4.2
##  [199] 3.2 4.1 4.3 4.6 3.0 5.2 6.2 5.4 4.4 3.1 4.9 5.7 5.1 5.5 5.5 3.9 5.1 4.2
##  [217] 3.2 5.0 3.2 5.9 5.7 5.5 4.8 4.2 6.1 4.6 4.1 6.1 5.9 4.6 4.8 3.6 6.1 3.7
##  [235] 4.0 3.7 5.5 2.4 3.5 4.2 3.2 3.7 4.7 6.0 2.5 5.0 4.2 4.4 5.4 3.5 4.5 4.8
##  [253] 4.0 4.5 4.0 2.9 3.6 3.7 4.3 2.9 4.2 3.2 6.6 5.0 3.2 4.4 4.1 4.6 4.3 4.9
##  [271] 5.7 5.0 4.1 4.1 4.6 5.5 5.0 4.6 2.8 4.9 2.9 4.7 5.2 3.8 4.3 4.7 4.4 3.3
##  [289] 3.9 5.4 5.8 3.8 3.9 5.2 4.8 4.9 3.6 3.4 5.0 3.7 5.0 2.5 3.7 3.3 5.7 4.2
##  [307] 5.3 3.4 2.4 4.9 3.9 6.2 3.8 3.7 5.3 5.3 1.9 4.9 4.7 4.0 4.6 5.8 3.9 3.9
##  [325] 4.5 5.1 3.0 3.2 4.3 5.4 4.4 3.2 3.4 4.9 4.7 3.6 5.3 4.6 5.1 4.2 4.2 5.5
##  [343] 3.4 3.7 2.4 3.2 3.7 5.2 5.8 4.6 4.2 4.7 5.3 5.3 5.4 3.6 5.1 5.3 3.5 4.6
##  [361] 3.3 5.5 5.7 5.1 4.9 4.4 4.1 4.3 3.3 4.1 4.3 4.2 3.2 4.6 5.7 4.8 4.0 4.4
##  [379] 5.2 5.7 4.0 5.8 4.5 5.7 4.7 5.9 2.5 4.9 5.4 4.9 5.3 4.3 5.1 4.1 5.2 3.9
##  [397] 5.2 4.3 2.8 5.6 3.3 3.7 6.8 4.2 4.3 3.4 4.9 4.6 4.5 4.7 3.9 4.2 4.3 4.6
##  [415] 3.5 4.5 3.4 3.7 4.9 4.4 4.7 5.0 5.1 5.1 3.4 3.8 5.3 3.8 4.0 4.7 3.7 4.8
##  [433] 5.0 4.0 3.0 3.6 5.0 3.0 3.1 3.9 4.3 5.9 4.8 6.8 2.5 5.7 3.4 6.0 5.2 5.2
##  [451] 2.8 6.7 4.9 5.5 3.6 4.4 5.3 3.4 4.8 4.9 4.6 3.9 4.8 4.1 4.8 3.4 4.7 2.9
##  [469] 5.2 5.3 4.6 5.2 4.7 4.2 3.0 4.8 5.5 4.5 5.9 2.9 4.6 4.2 4.4 3.7 3.9 4.0
##  [487] 3.7 5.7 4.5 4.9 5.4 4.0 3.5 5.0 5.3 5.3 4.3 5.6 4.6 4.4 5.8 3.0 3.1 4.5
##  [505] 5.7 5.4 4.3 2.4 3.7 3.9 5.2 3.3 5.0 4.4 2.9 4.8 4.6 3.2 4.0 4.9 4.6 4.3
##  [523] 3.6 5.6 3.4 4.0 5.1 4.1 3.9 5.7 4.6 3.8 4.2 3.6 5.1 5.5 4.2 3.4 5.5 5.5
##  [541] 4.7 4.8 5.7 4.8 6.4 4.1 4.6 2.7 4.6 5.4 4.9 3.4 6.4 5.1 4.8 4.0 4.5 4.2
##  [559] 3.9 3.7 3.5 4.8 4.2 4.7 2.3 6.1 5.0 5.9 5.0 4.8 3.1 6.0 4.1 6.8 4.9 4.4
##  [577] 4.5 5.1 5.0 4.4 5.2 2.3 5.3 4.1 5.5 5.4 5.4 1.5 4.2 5.9 4.3 6.4 6.2 6.0
##  [595] 5.0 3.5 4.0 5.3 4.7 6.6 6.8 4.8 5.1 3.0 5.6 4.7 4.2 3.6 5.1 4.4 5.1 3.7
##  [613] 6.2 4.1 5.5 3.4 4.7 6.0 3.8 5.7 4.9 2.3 4.2 4.7 2.5 4.1 4.4 3.2 3.7 5.1
##  [631] 3.6 4.1 4.3 4.6 4.0 4.0 3.5 5.7 5.6 4.1 5.9 4.4 5.8 5.8 5.7 4.0 6.3 5.5
##  [649] 3.6 4.3 4.5 5.2 4.4 4.0 5.2 4.1 3.7 5.2 4.2 4.9 6.4 5.8 5.0 5.7 2.9 5.1
##  [667] 4.1 3.8 5.7 3.7 5.3 4.2 4.4 4.8 4.5 4.7 5.6 4.7 2.4 6.4 4.8 4.5 3.9 5.3
##  [685] 4.7 5.3 4.5 4.6 4.1 4.6 3.1 4.2 3.9 5.5 3.4 2.7 3.2 4.6 4.6 6.3 4.4 4.8
##  [703] 4.8 5.2 4.0 4.6 6.0 4.4 5.0 5.6 4.6 3.9 4.7 4.0 4.2 3.9 4.6 5.1 2.9 4.6
##  [721] 5.5 5.2 5.2 5.5 3.0 6.2 3.3 4.8 5.9 5.1 3.8 6.5 5.8 4.6 6.7 3.0 3.8 4.8
##  [739] 3.0 2.2 3.5 3.4 3.8 4.3 3.0 4.5 5.0 5.2 4.4 4.7 4.7 4.7 4.1 6.6 4.3 3.7
##  [757] 4.3 4.9 3.1 4.6 4.0 5.6 4.8 4.8 4.7 4.8 4.7 4.5 5.5 6.9 3.0 2.7 4.1 3.7
##  [775] 3.8 4.2 4.5 6.3 4.5 4.8 5.9 3.2 4.4 4.6 5.6 4.8 4.3 5.3 3.5 3.2 5.5 5.1
##  [793] 5.9 4.2 3.9 3.9 4.8 5.8 4.2 4.4 5.4 4.2 4.6 4.5 5.0 3.7 5.2 6.3 3.9 4.3
##  [811] 5.9 5.1 2.9 3.7 3.6 4.6 5.2 4.5 5.2 4.2 3.7 4.0 4.4 4.5 4.6 3.7 4.4 3.6
##  [829] 3.7 3.7 6.0 3.9 4.5 3.5 2.7 5.2 3.5 4.2 4.5 3.8 6.0 4.5 5.0 3.3 3.9 4.0
##  [847] 4.0 4.8 5.6 4.7 4.4 4.2 4.8 4.2 3.3 3.6 4.4 5.0 5.3 2.6 4.8 4.6 5.4 5.5
##  [865] 6.2 5.0 4.9 3.7 3.5 3.4 3.3 5.4 3.9 4.2 5.5 4.4 4.9 5.2 4.0 5.0 4.5 4.6
##  [883] 3.9 3.6 4.7 5.9 3.4 4.3 4.3 3.6 4.5 4.9 3.4 4.7 2.7 3.6 4.3 4.0 5.6 3.8
##  [901] 3.9 4.1 3.9 4.1 3.7 5.1 3.2 6.4 4.1 4.3 5.0 5.6 5.1 4.7 3.2 5.4 5.3 4.1
##  [919] 4.0 5.7 3.9 5.1 4.6 2.9 5.1 3.3 4.9 5.4 3.7 4.0 3.5 5.5 3.2 3.2 4.8 4.8
##  [937] 5.2 4.5 4.3 4.8 5.3 5.0 3.8 3.4 5.4 5.2 5.6 4.0 3.6 6.1 5.4 5.9 5.7 4.3
##  [955] 4.8 4.5 5.2 3.8 5.5 4.0 5.5 4.8 3.2 4.3 4.4 4.9 4.1 4.9 6.4 4.0 3.6 5.3
##  [973] 6.9 6.0 6.1 6.2 4.1 4.3 4.1 5.0 5.5 5.0 4.6 4.7 4.0 6.0 3.9 3.5 4.0 6.4
##  [991] 5.2 4.8 3.4 4.2 5.3 5.8 4.6 4.0 4.0 4.0
## 
## $func.thetastar
## [1] 0.0084
## 
## $jack.boot.val
##  [1]  0.54940120  0.33176796  0.36657143  0.22962963  0.11893491 -0.09129213
##  [7] -0.13982036 -0.27382353 -0.46253521 -0.49553571
## 
## $jack.boot.se
## [1] 1.027798
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
##    [1] 4.1 3.8 3.2 3.3 3.9 4.5 6.6 3.8 6.2 5.3 4.3 5.0 3.6 5.5 4.3 4.7 3.6 5.2
##   [19] 5.2 4.8 5.5 4.5 4.9 4.0 3.6 3.5 5.0 4.5 4.3 5.5 4.5 3.0 4.1 4.4 3.8 5.8
##   [37] 3.8 4.8 3.5 4.2 3.9 3.7 3.2 2.9 4.3 5.2 4.2 4.0 3.9 5.1 4.0 4.2 4.8 3.2
##   [55] 6.1 3.6 6.8 3.7 5.7 5.0 3.3 6.2 4.1 5.2 2.9 3.5 4.8 4.9 3.9 4.1 4.0 2.8
##   [73] 4.3 4.9 4.4 4.7 3.4 4.6 5.0 5.3 4.0 4.5 3.7 3.1 5.0 5.1 4.9 6.0 4.0 3.7
##   [91] 3.7 4.9 4.5 4.3 2.7 3.4 3.4 4.9 5.1 5.0 4.4 4.1 3.8 4.8 4.2 2.6 4.9 5.5
##  [109] 3.6 4.2 5.6 3.2 5.6 5.0 5.0 4.9 4.6 4.6 4.6 3.9 4.0 4.8 5.3 4.4 3.9 3.9
##  [127] 3.3 4.6 6.3 4.5 3.6 3.9 3.2 5.2 5.1 5.0 4.4 5.6 4.6 4.4 6.6 3.9 5.1 4.6
##  [145] 5.2 5.8 5.8 3.6 4.9 4.4 4.1 3.4 4.5 4.4 4.0 4.3 5.8 3.0 4.3 4.7 3.5 4.8
##  [163] 3.0 5.6 5.1 3.6 5.1 3.8 6.4 4.7 2.9 4.0 5.4 3.6 4.2 4.2 3.6 3.6 4.8 4.4
##  [181] 4.8 4.1 4.9 4.5 4.3 5.5 4.5 4.7 4.9 4.0 5.4 2.7 3.1 4.2 4.7 4.5 5.2 4.7
##  [199] 5.8 5.0 5.3 3.9 3.6 5.0 3.9 4.3 5.4 3.1 4.8 5.5 5.1 5.7 5.1 3.1 2.7 3.7
##  [217] 4.4 5.2 4.2 4.8 4.3 4.9 5.1 4.4 5.1 6.3 5.0 4.8 3.8 4.3 3.7 4.2 3.9 5.2
##  [235] 3.9 4.3 2.7 4.1 4.8 6.3 4.8 4.6 4.0 3.8 4.9 4.6 5.3 3.7 4.7 5.3 6.5 5.0
##  [253] 4.5 5.5 4.6 4.5 4.6 5.7 2.6 4.3 3.3 3.4 4.6 4.1 5.5 5.7 5.6 4.3 4.4 5.5
##  [271] 5.8 3.3 4.5 5.4 4.2 3.5 3.6 4.2 2.7 4.7 4.6 5.3 4.6 3.7 5.4 6.1 5.0 5.2
##  [289] 3.1 4.9 5.9 3.9 4.6 4.9 2.8 4.0 3.0 4.2 2.8 5.2 4.3 4.9 4.9 6.2 4.4 4.9
##  [307] 5.3 5.0 5.8 3.4 5.3 4.6 4.9 3.5 5.3 4.9 4.8 4.0 3.7 4.5 4.3 5.2 4.2 3.4
##  [325] 3.1 5.2 7.2 4.8 3.9 5.1 3.4 6.0 4.7 5.9 4.8 4.7 6.0 3.8 3.5 2.0 3.5 5.5
##  [343] 4.5 3.5 4.3 4.6 5.0 3.2 4.7 3.4 3.2 4.2 3.5 4.4 3.8 4.5 4.7 4.4 3.5 3.6
##  [361] 3.8 3.7 3.4 3.7 4.1 4.4 4.3 4.9 4.7 4.4 3.9 5.3 2.6 4.9 4.8 5.8 5.2 5.2
##  [379] 4.7 4.2 3.7 5.5 4.2 7.0 5.4 4.4 4.4 4.1 5.3 5.2 3.8 5.4 2.6 4.0 5.0 3.8
##  [397] 4.6 4.6 4.7 4.5 5.5 4.6 3.9 4.3 3.5 3.7 4.4 5.3 3.1 5.2 4.7 5.1 5.0 5.1
##  [415] 3.8 4.1 4.4 5.5 4.3 5.0 5.3 2.5 2.9 3.8 4.2 4.0 5.1 4.5 3.7 3.6 4.6 3.3
##  [433] 3.9 4.8 4.5 3.9 4.6 3.7 5.9 4.1 5.3 2.3 3.2 4.6 3.5 4.5 4.3 4.2 3.9 5.8
##  [451] 3.9 4.8 4.1 4.6 3.2 5.0 3.3 4.6 5.0 4.3 4.8 4.3 6.3 5.3 5.1 4.0 2.8 3.4
##  [469] 4.2 6.4 2.0 4.2 4.6 3.9 5.2 4.9 4.9 3.3 4.6 6.5 4.6 3.6 4.3 4.0 4.5 5.7
##  [487] 3.2 5.0 5.0 4.5 5.7 4.7 5.2 4.5 4.3 2.8 5.5 4.5 4.6 5.4 4.1 4.9 5.2 4.5
##  [505] 5.2 4.5 5.5 2.8 5.0 5.4 5.5 3.7 3.3 4.3 5.2 4.6 6.0 4.3 5.4 4.6 5.5 4.9
##  [523] 6.9 3.4 3.8 5.6 4.4 3.8 5.8 3.1 3.7 4.5 3.1 3.6 5.5 5.0 4.8 4.3 4.7 4.3
##  [541] 4.9 5.3 3.0 4.1 5.3 3.2 5.6 3.9 3.4 2.9 4.0 4.1 4.1 5.6 4.3 6.6 3.8 5.1
##  [559] 4.7 3.7 4.1 5.3 4.3 5.1 4.2 4.7 4.6 3.7 4.6 5.0 5.3 3.9 4.3 4.6 6.2 3.2
##  [577] 4.5 4.9 4.7 4.4 4.7 5.9 5.0 4.9 4.2 4.9 3.3 2.7 5.0 5.3 3.4 3.9 4.6 3.3
##  [595] 5.2 3.8 5.4 3.0 6.3 4.8 4.4 4.3 5.0 4.9 4.5 5.0 3.9 4.8 3.2 4.4 3.4 3.0
##  [613] 4.2 4.0 3.8 4.9 4.8 5.3 4.7 5.0 5.3 3.1 3.6 3.9 4.4 4.8 5.7 4.3 5.6 4.4
##  [631] 5.5 3.7 4.8 5.1 4.4 4.6 5.1 2.9 4.6 4.4 5.2 3.7 5.8 4.7 5.6 4.4 4.0 5.1
##  [649] 6.0 3.9 5.6 4.8 4.6 3.7 4.7 5.3 3.9 4.5 4.1 2.9 3.3 4.0 2.6 5.4 3.9 2.0
##  [667] 4.2 3.6 3.9 3.7 3.2 2.9 4.1 4.2 5.2 4.5 3.4 5.3 5.6 3.6 4.3 4.5 3.8 4.6
##  [685] 5.4 5.5 3.1 3.7 5.5 3.6 4.2 4.9 4.9 5.0 5.6 3.8 3.7 3.8 6.3 4.6 5.1 5.5
##  [703] 5.4 5.2 4.9 4.2 4.4 3.4 3.8 4.2 3.6 4.6 5.3 4.0 4.9 5.3 4.6 4.4 4.7 4.7
##  [721] 4.5 4.4 4.0 5.3 3.9 4.6 2.9 5.4 5.7 5.5 5.3 4.1 4.1 3.4 4.2 3.4 4.0 4.7
##  [739] 2.4 4.0 4.9 4.6 4.4 4.4 4.5 3.9 5.3 2.4 4.9 5.0 3.3 3.9 4.2 2.7 3.1 3.0
##  [757] 5.8 4.9 4.0 3.9 4.7 5.0 5.3 5.8 4.8 4.7 5.2 3.8 4.1 3.5 6.1 3.8 3.5 5.0
##  [775] 4.2 4.6 3.5 5.0 4.6 5.8 4.9 4.2 4.4 5.2 5.1 2.7 4.6 4.5 3.6 4.1 4.3 4.4
##  [793] 3.0 4.4 4.5 3.8 4.7 4.0 4.0 4.6 5.9 5.3 5.0 3.1 6.1 4.8 6.7 4.7 2.8 4.9
##  [811] 5.3 5.2 5.3 2.9 3.9 6.7 5.0 3.4 5.8 4.7 5.4 5.3 4.8 5.8 6.5 3.4 4.8 2.1
##  [829] 3.9 5.4 5.1 4.6 4.1 4.1 3.6 3.7 5.2 4.0 4.2 4.9 4.1 5.2 6.7 4.2 3.3 4.3
##  [847] 3.9 4.2 3.2 3.7 4.0 3.0 3.2 5.5 4.2 4.7 3.9 5.9 6.2 5.9 5.0 3.3 4.9 4.1
##  [865] 5.2 4.8 4.4 6.0 3.5 3.9 5.2 5.1 3.4 5.5 3.7 4.8 5.1 4.0 4.7 3.4 3.7 5.4
##  [883] 5.2 3.7 3.8 5.0 5.6 5.5 5.8 5.5 5.3 4.9 5.9 3.1 5.0 3.4 3.2 4.4 3.3 5.0
##  [901] 3.3 4.3 4.4 3.9 5.1 4.1 3.7 3.4 4.9 4.6 5.5 2.7 6.0 1.5 4.1 2.9 3.8 3.9
##  [919] 5.6 5.4 4.5 5.5 3.4 4.6 5.6 4.6 5.6 3.8 5.1 5.6 5.1 5.5 3.9 4.0 3.6 5.1
##  [937] 5.1 3.9 3.7 4.8 3.9 3.2 3.1 3.8 5.7 3.4 5.9 4.8 4.8 3.6 3.8 2.8 4.2 4.1
##  [955] 5.4 4.2 4.8 5.9 4.2 5.1 4.2 4.0 4.6 4.8 4.5 5.8 2.8 4.7 5.1 4.5 5.9 4.3
##  [973] 4.6 5.5 4.7 4.2 3.9 3.0 4.4 3.9 3.1 5.8 4.3 4.1 5.6 4.6 3.9 4.9 3.0 5.3
##  [991] 4.8 5.1 5.8 4.7 4.8 3.8 4.9 3.7 2.9 4.2
## 
## $func.thetastar
## 72% 
##   5 
## 
## $jack.boot.val
##  [1] 5.400 5.300 5.200 5.100 5.068 4.900 4.800 4.800 4.600 4.500
## 
## $jack.boot.se
## [1] 0.8451281
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
## [1] 0.7294112
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
##   4.086516   6.308636 
##  (1.758043) (2.887865)
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
## [1]  0.2967908 -0.3820790 -0.1524023  0.1830765  0.7357794  0.6335143
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
##    [1]  0.0882405873  0.1684601399 -0.0863682082  1.4690751259  1.1752623965
##    [6]  2.2300357581 -0.3747574273 -0.3842470169  1.1430760497  2.2107474462
##   [11]  0.1694435258  0.2896953152  1.1467186607  0.3318879972  0.8191247435
##   [16]  0.9342695121  1.8736405201  0.2533248937 -0.0317440412  0.0602776146
##   [21] -0.4604933159  1.4027095616  0.6287625327  0.3527934913  0.3502424320
##   [26]  0.3186261967  0.6995477301  1.7097298974  0.8320598509  1.2347575631
##   [31]  0.8545451108 -0.0078477445 -0.3404990935  0.4103004754  1.1719915395
##   [36]  0.8128535482  0.0217588579  0.3746630150  0.9067946441  0.7910862491
##   [41]  0.5079715692  0.3184931010  0.7861871205  0.0219093208  1.1860153161
##   [46]  0.4300903831  1.1188374100  1.3771304478  2.1968729173  0.0284074034
##   [51]  0.3513607820  0.6920724044  0.8028535621  1.7918431592  0.8479803431
##   [56] -0.3620658130 -0.6586326820  1.3274216453  0.8075776125  0.3259216292
##   [61]  0.3587924547  0.7299971014  0.3795709484  0.4721655266  0.4275769901
##   [66]  0.0736783281  0.4530913333  0.7774520408  0.4897962561  1.1366548511
##   [71] -0.0498011530  0.7112141760 -0.4697114583  0.7687260786  0.7670363413
##   [76] -0.0100855601  0.0380250116  1.2353739775  0.6317143685  0.8221615979
##   [81]  1.1647476843  0.6202785426  1.0699994804  0.6408592180  0.6421639660
##   [86]  0.6597184761  0.2673913358  1.0712496053  0.3928048983  1.2466626944
##   [91]  1.6295763917  0.1993349497 -0.3552550776  0.3674859702  1.0690550852
##   [96]  0.9319933432  0.4450908367  0.2845021193  1.3371125448  1.1213241390
##  [101]  1.1020396643  2.1846218306  0.0246925702 -0.4098186585  0.4275769901
##  [106]  0.0820611045  0.6460892865  0.8563905333  1.3769137659  0.0583884739
##  [111]  0.2502988629 -0.8266433938  1.3847121846 -0.0529772602 -0.3377708814
##  [116]  0.8289188132  0.4330091591 -0.0415487177  1.2165613959  0.2774539793
##  [121]  1.3234091130  0.2993292446  0.7331708823  0.7395219999  1.2636300964
##  [126]  0.4045613189 -0.0212611751  1.2173410604 -0.1047603201  0.4909076978
##  [131]  1.5914521006  0.4063952296  0.7902110715  1.1023122330  0.7755724294
##  [136] -0.0345503822  1.1548658901  1.3835485268  0.3860581266  0.9154922471
##  [141]  0.2158261361  0.2615315387 -0.2440052515  1.6022384745  0.6428599631
##  [146]  0.3894105796 -0.0883059081 -0.8353558177  0.3817997365  1.0752350158
##  [151]  1.1144288324  0.3860581266  1.0882472124  0.4893539179  0.0404988074
##  [156]  1.1908067485  0.4319211508  0.0378065245  0.7688557746  1.6646185303
##  [161]  0.2909276739  0.5080759203 -0.7501288964  1.0764013406  0.8593190991
##  [166]  1.3072038806  0.2403390507  1.1690060061  0.3665730273  0.3668969804
##  [171]  0.3807205698  1.2919794391  0.9087398695  0.6119535054  0.2878025541
##  [176]  0.3881548865  1.8335524153  0.4193006775  1.0658690108  1.3497113375
##  [181]  2.0347557055  0.5876055163  0.7362209306  0.3174723992  1.1319475189
##  [186]  2.5183628310 -0.1915994659  0.1610913879  0.0302667311 -0.0510345951
##  [191]  0.4676315141  0.3857606782  1.3251944196  0.6261324224  0.8189591908
##  [196]  0.7336573300  0.8283834173  0.9942727934  0.2791370307  0.0231561189
##  [201] -0.0180944562 -0.0286389932  0.4016618627  0.6847824088 -0.2567454711
##  [206]  1.3278349382  0.1644833580  0.3349799754 -0.8164758272  0.2427257348
##  [211]  1.6137260501  0.4944475875  0.6936290114  0.7834906870 -0.4444328850
##  [216]  0.3086166227  0.4727392755  1.2240003728  0.1769737426  0.6813391832
##  [221] -0.2424794182  0.3597457980  0.4629450808  1.0671893542  0.3732283700
##  [226]  0.3246846771  0.6839659048  1.1319475189  0.8195835971  1.1801413823
##  [231]  0.0973905247  0.6067264055  0.7416647585  0.3242734320  0.8194779803
##  [236]  0.6438831032  0.9784631024  1.2076957070  0.8109602364  1.3841732046
##  [241]  0.3517103156 -0.0187525596  0.3914700333  1.2034253910 -0.1993206310
##  [246]  0.3228446987  0.6866096387  0.3384285241  0.2757109446  0.6702397249
##  [251]  0.8878722184  0.4954809163  2.0847555741  1.8247974385  1.1028704218
##  [256]  1.4541782687 -0.0314357727  1.0046769879  1.1448271466  1.1742369420
##  [261]  0.7647865962 -0.3544375147  0.7647865962  1.3841922937  1.6670504079
##  [266]  1.5875470694  2.1525940634 -0.0108910731  0.8868473775  0.8239469452
##  [271]  0.7653961188 -1.3333824922  1.4408889012  0.3759020996  1.1157918547
##  [276]  0.4863620589  0.7731014954  0.8253886859  0.0570967387  0.6456594747
##  [281]  0.8176118790  0.4611886847  0.5834223036 -0.0740690214  0.8391925094
##  [286]  1.3340812949  0.4267699537  0.1311150588  1.0616595009  1.2689411957
##  [291]  1.6574793806  1.1427448831  2.4914103593  2.1791819275 -0.4439742151
##  [296]  0.4415026708  0.8275291870  1.6650782553 -0.4787540280  0.8152637141
##  [301]  1.3794810270  0.3925085989  0.4258172604  0.3467632963  0.4395586822
##  [306]  0.3983423480 -0.0833174125  0.4370996809  1.2313088460  0.5520143004
##  [311] -0.0989940179  0.7258251034  0.8227031352  1.1078887459  0.1291894478
##  [316] -0.8808553226  0.5627931138  1.6666467849  0.3599226968  0.9949843946
##  [321]  0.8539881998  0.8184653420  1.2789432463  2.1505706742  2.4311405292
##  [326]  1.3342356711  0.3259915331 -0.2711497759  1.3371125448  0.0438082095
##  [331]  1.9965770043  0.5992655442  0.5658024498  0.9946579853  0.7826095591
##  [336]  0.3571555237  0.4096170881  1.2152253731  0.3005496817  1.1801413823
##  [341] -0.0883059081  1.2951956459  0.0280709221 -0.0206305534  1.4652021151
##  [346]  0.3773966393 -2.2089122999  0.3809917160  0.3275450241  0.3771525872
##  [351]  0.6653551290  0.8698424524  0.2311696007  0.3831443030  1.2609181165
##  [356]  0.7906748922  0.6821406110  0.7289129827  0.5434910549  0.7639044707
##  [361]  0.4044587221  0.4066532583  1.8448118146  2.0798178774  1.3512497234
##  [366]  1.0577423166  0.6093698398  1.2946174552  0.4280671212 -0.1149810420
##  [371] -0.4136513924  0.2838112638  0.7073512981  1.2002045271  0.4256571886
##  [376]  0.8569458369  0.3549591896  0.8086398130  0.6438831032  0.7463881779
##  [381]  0.5031829159  0.3732506719  0.8337262362 -0.0061623833  1.4048634706
##  [386] -0.0124826999  1.1826897725 -0.0972040783  0.7286586769  0.3977586117
##  [391]  1.1385764007  1.2632156659  0.3652399093  0.6652442215 -0.0777315427
##  [396]  1.2355534761  1.5765439385  0.7663119949  1.2188241650  1.3460521187
##  [401]  0.3934879217  0.4058754285  1.2866799228  0.7648801188  0.3938823373
##  [406]  0.0758097101  1.5870498823  0.9986119006  1.1842144518  0.1778277341
##  [411]  1.3779431063 -0.2594654926  0.5855842890 -0.7728298722  0.3970234023
##  [416]  0.7933988248 -0.2931806585  0.9583871013  0.8595792688  0.4516007472
##  [421] -0.3741852235  0.6192027340  0.6434502011  0.9853576815  0.0359539929
##  [426]  0.3433778107  0.4305412222  0.7870562779  0.7263488466  1.1434961353
##  [431]  1.8891403860  0.8276065739  1.1392701947  1.2726029785  0.3234441734
##  [436]  1.1172166608  0.9041743440  0.0245933702  0.7008772614  0.7283021959
##  [441]  0.3564708796  1.5742255169  1.0191082980 -0.0531243213  0.7957551444
##  [446]  0.0503969912  1.0983926170  0.0236448058  0.8071507202  0.1120573371
##  [451]  1.8914580298  0.3873673591  0.4043804283  1.6382415348  0.4666514952
##  [456]  0.7763853368  0.4896793850  1.3155081091 -0.0228608286  0.7516186233
##  [461]  0.2352476854  1.2539079808  1.5381399726  0.6983442759  0.4399874635
##  [466]  0.8309429551  1.3651182432  2.0056229063  0.9508100994 -0.0354963800
##  [471] -0.0108291672 -0.3856900928  0.2948838993  0.4096340063  0.3759224659
##  [476]  0.4909076978  0.3096209397  0.3303016430  0.4377349840  1.0020545181
##  [481]  1.4270452480  0.8271098382 -0.0429795154  0.2524289373  1.3496371751
##  [486]  1.2065559511  0.8310238156  0.3635347799  0.8398224699  0.6276800856
##  [491] -0.0873134179  1.0672362332  0.9998093880  0.3845685817  0.3651923054
##  [496] -0.0072761345  0.3264549378  0.0473131414  0.0832164771 -0.0504605622
##  [501]  0.7950186355  0.6562468633 -0.2084308916  0.6617587318  0.0584258474
##  [506]  1.5351302632  0.3949323773 -0.0206200653  0.3270420450  1.4499053677
##  [511]  0.1784799190  0.7488964918  0.6862456344  0.4444629674  0.0416382294
##  [516]  1.2545908450  1.4185154357  1.2803013100  0.6910001435  0.5982469866
##  [521]  0.3675681535  0.4448227229  2.2275733206  0.0214377120  0.6053596690
##  [526]  0.7657097693  0.0242863736  0.0185195960  0.7382894881  0.4097381426
##  [531]  0.7696087848 -0.1227151661  0.0191578268  1.5675220501  1.0265722055
##  [536]  0.3484844236  0.3869507897  0.2989318034  1.0626420355  0.6423885481
##  [541]  1.3375280831 -0.2087249489  0.2201445674  1.1321722872  0.0054715614
##  [546] -0.0008361423  0.9826372394  0.8481634998 -0.3736592345  1.0045008820
##  [551]  1.4073681573  1.1680219011  0.9561486893 -0.0892530120  0.2339494205
##  [556]  0.3469690353  0.5472098450  0.7289890519  0.4593208357  0.5678084761
##  [561]  0.0050318211  0.9811950594  0.7645210569 -0.0362703023  1.6779913087
##  [566]  0.8742687175 -0.1251008859  0.7646626772  1.0237576926  0.6101858715
##  [571]  0.7074926864  0.5082178622  0.4594322885  0.9539208136  2.4794097146
##  [576]  0.7684125778  0.2189909899 -0.0097673884  0.4025148544  1.3408487532
##  [581]  0.4394695217  1.1096045996  1.0005618067  0.7530355954  0.4115892098
##  [586]  0.7538867033  0.7856970738  0.6628343083  1.2159496721  1.3781176583
##  [591]  0.5213227739 -0.0815410877  0.2222566675  0.5566478632  0.0392274618
##  [596]  0.7639225168  0.6971258658  1.0518228513  0.6878467919  0.4108511047
##  [601]  0.8558831217  0.7693995809  1.5893393986  0.8647702572  0.4233659677
##  [606]  0.8322590466  2.5818741228  0.2961797603  0.1035095845  0.4927640332
##  [611]  0.3608484655  1.8244852616  0.2028259357  0.7100197833  0.4235164873
##  [616]  2.2175546642  1.3480427181  1.0695616602  0.9324758014  0.3895521503
##  [621]  0.2632903162  0.4151866062  0.8719274113  1.2813938711  1.0841107904
##  [626]  0.8117366738 -0.0124646899  0.3433778107  1.7046724175  0.5293566057
##  [631]  0.3186884695  0.6040592083 -0.0120680414  0.0226992687 -0.0632629577
##  [636]  1.3264123475  0.2171115171  1.2818909899  0.0199973388  0.2314809424
##  [641]  0.1617093532  0.7762995657  0.8220160102 -0.0494907006 -0.7069331956
##  [646]  1.3404554648  1.4010136833  1.7128432334  1.1806549264  0.4000318870
##  [651] -1.2011126584  1.0160135162  0.7893449896  0.0641566881  0.8649759423
##  [656] -0.3743236613  1.4691395806  1.3258800656  1.5635541529 -0.2227738351
##  [661]  0.0205007931  0.0225134983  0.9324439442  2.1685305056  0.2046248711
##  [666]  0.9304057691  0.7902549718  1.5754561822 -0.0208274028  1.9253529288
##  [671]  0.8186044588  0.7114279622  1.4652021151  1.0169702261  0.1694020240
##  [676]  0.3963857783  0.3534435031  0.9903947635 -0.0361685172 -0.4020945877
##  [681]  1.9500355115  1.0582162600  0.4348107081 -0.3191534089  1.0638758723
##  [686]  1.3571208018  0.8369235511  1.1995917690 -0.0200948722  0.3953043250
##  [691]  0.4372006384 -0.3694136359  1.1540133283  0.4110163211  0.7652621107
##  [696]  0.7644917205  0.4242220690  0.0244909734  0.9258144310  0.6604649315
##  [701]  0.8522859644  0.8446001331  1.3631054548 -0.0223918375 -0.3394522100
##  [706]  0.0048054428  0.4123364027  0.7096728430  0.5878666995 -0.0106263807
##  [711]  0.7297539123  0.3759224659  1.6000576201  0.9324439442 -0.0867473707
##  [716]  0.7113086389  0.4115767113  1.0169702261  1.1326975958  0.1085110603
##  [721]  0.5562893557  1.3961241146  0.4368890316  0.7289027834 -0.3255114193
##  [726]  1.0964776603  1.1075347334  1.3554879228  0.8258809547  2.2298396855
##  [731]  0.0036515802  0.5030441586  0.3267905105  2.6113161292  0.7098972934
##  [736]  0.4909076978  2.1517698181  0.0045051971  1.6720558896 -0.7980373902
##  [741]  1.1244814114  1.4051219467 -0.0223345056  1.1763302545  1.2366672083
##  [746]  0.8037203335  0.6895734466  0.7461903386  1.1704274315  0.5038232398
##  [751]  0.4296242032  1.1149645813  0.0612753430  0.8320040182  0.8641024379
##  [756]  0.0571647186 -0.1423301780  0.7255373785  0.6032026680  1.2959298383
##  [761]  2.5830016940  0.9484464854  0.2246854639  0.3384986591  2.2689146818
##  [766]  0.0577218066  0.0628681791  0.4265097489  0.0482974442 -0.2249042120
##  [771]  0.0258491032 -0.0439795514  0.3608484655  0.7330136239  0.6698039124
##  [776]  0.0194115216  0.6495586120  1.0735210955  0.7841437828  0.6976078819
##  [781]  0.3076348550  0.2520531119  0.7645484199 -0.0182324463  0.1011954658
##  [786]  0.2478926080  1.1611544797  0.2492044387 -0.0339088813 -0.0122276911
##  [791]  0.0917029649  0.0793940668  0.0056161447  0.3613301940  0.4471760959
##  [796]  1.4264522867  0.7036292037  1.0642251807  0.5676116211  0.7329609339
##  [801]  0.0125889289  0.6579553600 -0.5521169895  0.3569449114  0.4059625986
##  [806]  0.5797446503  0.1323652059  1.2341511886  0.7338635565  0.3495836820
##  [811]  1.0918343808  1.9608082573  1.6778697102  0.3905447453  1.1801968527
##  [816]  0.3768417563  0.3926270909 -0.1122621994  1.5917247665  0.0340233994
##  [821]  0.0223387448  0.0456604519  0.7411545208  1.3956464986  1.2156909349
##  [826]  0.3839603090  1.1147404171  0.8115699811  0.5212661624  0.7018276119
##  [831]  0.5049578828  0.7612885687  0.7959399558  1.0740312442  0.7256070155
##  [836]  0.0065509662 -0.3210118212  0.3814800537  0.3607958378 -0.0209036866
##  [841]  0.3906884145  2.2712200992  2.0689613166  0.2371103466  0.6111266067
##  [846]  0.2667827551  0.8168375002  1.4015008027  0.4058064547  0.0322551096
##  [851]  0.4866998858  0.1760942139  1.2367983657  0.8180070212  1.4071446950
##  [856]  0.7117547193  0.8147435533  0.3196985417  1.6616855958  0.7657097693
##  [861]  0.7287290707  1.0086387856  0.7731021560  1.5214888242  1.0729009432
##  [866]  2.3764986215  0.6893542853  0.2486194893  2.3452142480  0.6063872707
##  [871]  0.9726459319  0.4211803549  1.3774345079  0.3545506617  0.3700095170
##  [876]  0.0527147979  0.0297067853  0.7931542645  0.4231946035  1.7177747724
##  [881]  0.0273276154  1.0331497205  1.3158929905  0.3888916522  1.6425855553
##  [886]  0.7870587765 -0.0218110495  0.8755693319 -0.3917348952  0.3192288331
##  [891]  0.2995060599  0.8814688112  0.9242272957  0.0962172717  0.0212919903
##  [896]  1.7040748104  0.2817571519  0.3752053733 -0.3644333198  0.4323071507
##  [901]  1.6681131413  0.5616043974  0.0896102146  2.3403816254  2.3594527889
##  [906]  1.5108192138  1.3658993917  0.6958196627  0.3087877272  0.4020739746
##  [911] -0.2445991143  0.5954121764  0.6989565291  0.8365034412  1.3171099986
##  [916]  0.7040292951  0.7730680798  0.3410370829  0.0042136795 -0.0256681408
##  [921]  1.4915081120  0.1005202097  0.3816584099  0.7013798784  0.3820694602
##  [926]  0.0337003035  0.8117204247  1.2189041500  1.0269700743  0.2177269548
##  [931]  0.7173974490  1.4562062281  2.1077468761  0.4287987973  1.2327873914
##  [936]  1.1101864472  1.0046008272 -0.3441129893  0.7532216043  0.8210872022
##  [941]  0.4311724394  0.8310342871  0.7289890519  1.1705819911  0.0981157272
##  [946]  1.4114916123  0.3817855432  0.3799726251  1.2522219940 -0.1177732035
##  [951]  0.7838871747 -0.1669904096  2.2925635489  1.2551598272  0.8294411747
##  [956]  0.1254515171 -0.0521127532  0.3190192596  1.3304792315  0.1825665306
##  [961]  1.4196314802  0.6701745710  1.5845156967 -0.0961042478  0.3794798788
##  [966]  0.0371664660  0.4270704291  0.7022107563  0.7255373785  1.4554807918
##  [971]  1.6815893430  0.3659557569  0.3723704655 -0.4080816296  0.8171070284
##  [976] -0.0870852466  0.4407237986  0.3969893537  0.7309905684  1.3633796216
##  [981]  1.7775807552  0.8664006379  0.3353752687  0.1954277262  1.0092347709
##  [986]  0.6983442759  0.7248525493  1.5952829278  0.3554361795  0.2965297740
##  [991] -0.0189971826  0.8945290239  1.3471354866  1.6659991797  1.0168262259
##  [996]  0.3795894036  0.7595838392  0.5043626536  0.0802346097  1.8792371855
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
##   0.64776237   0.33525418 
##  (0.10601668) (0.07496469)
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
## [1]  0.1208904 -0.1886559 -0.1810988  0.4774924  0.2964451 -0.5605078
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
## [1] -0.0036
```

```r
sd(results2$thetastar)  #the standard deviation of the theta stars is the SE of the statistic (in this case, the mean)
```

```
## [1] 0.9017429
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
## t1*      4.5 -0.06626627   0.8879471
```

One of the main advantages to the 'boot' package over the 'bootstrap' package is the nicer formatting of the output.

Going back to our original code, lets see how we could reproduce all of these numbers:


```r
table(sample(x,size=length(x),replace=T))
```

```
## 
## 1 2 3 6 8 9 
## 1 4 1 1 1 2
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
## [1] 0.0314
```

```r
se.boot
```

```
## [1] 0.8901159
```

Why do our numbers not agree exactly with those of the boot package? This is because our estimates of bias and standard error are just estimates, and they carry with them their own uncertainties. That is one of the reasons we might bother doing jackknife-after-bootstrap.

The 'boot' package has a LOT of functionality. If we have time, we will come back to some of these more complex functions later in the semester as we cover topics like regression and glm.

