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
## 0 1 3 6 7 9 
## 1 2 1 2 3 1
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
## [1] 0.0065
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
## [1] 2.847033
```

```r
UL.boot
```

```
## [1] 6.165967
```

Method #2: Simply take the quantiles of the bootstrap statistics


```r
quantile(xmeans,c(0.025,0.975))
```

```
##  2.5% 97.5% 
##   2.9   6.2
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
##    [1] 4.1 3.7 5.5 3.9 4.7 5.8 7.4 4.7 2.8 4.9 1.6 5.4 5.3 4.7 3.1 4.7 4.5 5.4
##   [19] 5.0 3.6 5.4 3.4 5.1 4.5 6.0 4.9 5.9 5.6 5.9 4.4 4.8 3.5 4.5 3.3 3.5 4.2
##   [37] 5.6 3.9 3.6 3.4 4.4 5.1 4.2 4.1 4.2 3.5 4.3 4.0 4.0 4.9 5.3 5.1 3.0 4.4
##   [55] 3.8 3.6 3.0 4.3 6.1 3.8 5.6 3.9 5.7 4.2 3.6 4.1 3.6 5.1 5.5 5.1 4.1 4.5
##   [73] 3.3 3.0 4.4 4.4 4.4 5.4 4.6 5.3 5.1 5.7 3.7 3.4 4.2 6.0 4.2 4.2 4.8 4.6
##   [91] 4.1 7.3 5.0 4.4 5.7 3.5 4.2 3.3 5.2 4.2 3.5 5.5 4.2 4.1 4.1 4.0 5.0 2.5
##  [109] 5.3 3.1 5.5 4.7 4.0 4.3 3.4 3.8 4.0 5.1 5.0 4.8 3.2 4.6 4.0 4.5 4.9 4.0
##  [127] 5.4 4.5 3.9 5.6 4.0 4.4 2.9 4.2 5.4 4.1 5.2 5.3 6.6 3.9 5.4 3.8 4.1 5.7
##  [145] 3.3 4.4 5.7 2.3 4.8 4.1 4.0 4.6 4.9 4.3 3.9 2.8 4.7 4.4 2.6 4.7 4.7 5.3
##  [163] 3.2 4.8 4.6 6.6 3.7 5.2 4.2 5.7 3.4 4.3 5.0 2.3 4.4 6.4 4.3 5.4 3.3 3.0
##  [181] 5.1 5.6 5.3 5.5 5.3 4.6 4.8 4.7 4.0 4.3 3.4 5.1 3.0 3.5 5.1 4.9 2.7 2.1
##  [199] 4.1 5.3 5.0 2.8 5.7 4.1 5.6 3.9 2.5 4.4 4.2 5.9 3.8 3.5 4.0 3.6 2.9 4.8
##  [217] 5.9 4.0 4.2 4.7 5.2 5.7 6.0 4.9 3.2 5.4 4.3 6.6 5.2 3.6 4.9 4.2 5.0 5.3
##  [235] 6.0 3.9 4.7 5.3 3.5 4.7 3.3 6.3 5.3 6.2 3.8 4.0 3.4 2.9 5.4 6.0 4.6 2.5
##  [253] 4.2 5.3 3.9 3.5 3.8 4.1 5.5 3.8 5.5 5.1 5.4 3.8 4.1 4.4 4.1 3.3 3.4 3.5
##  [271] 4.7 6.0 4.5 4.2 2.5 5.8 4.8 3.4 5.8 6.1 4.1 3.0 3.7 3.7 4.3 3.7 4.3 4.6
##  [289] 4.7 3.6 5.8 3.0 5.9 3.5 5.2 3.7 3.8 4.3 4.8 4.2 4.1 5.1 3.6 5.4 2.0 6.2
##  [307] 5.3 3.5 3.9 4.3 3.7 4.4 4.7 4.7 3.9 5.0 3.4 4.5 3.3 3.5 5.2 4.5 5.6 4.9
##  [325] 4.3 5.0 3.0 4.1 5.7 3.6 4.4 3.9 3.4 5.8 4.3 4.2 4.2 4.0 5.6 4.3 3.9 4.5
##  [343] 5.4 3.9 3.5 2.6 3.3 4.8 3.0 5.1 3.3 4.4 3.2 4.0 3.9 5.0 6.5 5.0 5.1 3.2
##  [361] 3.5 3.0 4.1 4.2 4.5 3.5 5.2 4.4 4.7 4.9 5.5 4.6 3.6 4.6 5.9 5.1 4.3 3.5
##  [379] 2.2 4.6 4.3 4.9 5.2 3.6 5.2 3.5 5.2 3.8 4.1 4.6 3.6 4.4 4.6 4.5 3.6 4.4
##  [397] 4.1 3.1 4.9 5.9 5.3 5.7 4.4 5.8 5.3 5.6 5.2 3.6 5.8 4.2 4.4 3.8 4.1 4.6
##  [415] 4.4 4.5 3.4 4.8 5.5 5.0 4.3 3.2 5.1 4.7 5.4 3.4 4.3 4.4 4.6 4.0 3.6 3.5
##  [433] 3.7 5.0 3.9 3.3 4.5 4.5 5.5 4.3 4.5 5.4 4.9 5.6 3.2 3.8 4.9 4.1 4.3 4.3
##  [451] 5.1 6.2 4.5 3.6 4.1 5.0 4.5 4.4 3.7 4.4 4.1 5.4 5.3 5.2 3.4 5.4 5.0 4.6
##  [469] 5.1 3.6 4.8 5.0 5.7 3.6 4.1 4.5 4.2 4.7 4.3 3.2 3.4 4.8 4.0 4.8 5.3 4.8
##  [487] 3.9 2.6 5.3 3.8 4.2 3.7 5.7 4.3 4.7 5.5 3.7 3.7 5.2 3.9 4.1 4.0 3.9 5.4
##  [505] 3.5 3.1 4.7 5.0 2.8 4.9 7.5 4.2 5.1 3.7 4.1 4.6 5.1 4.9 5.4 3.3 6.0 4.4
##  [523] 4.1 4.3 5.6 5.1 4.7 4.6 4.4 4.4 2.9 4.0 3.3 3.9 5.2 4.9 4.7 4.8 3.8 6.0
##  [541] 4.6 4.1 4.5 4.9 5.7 4.7 4.4 3.5 4.4 2.9 4.8 2.8 3.7 5.3 3.8 5.9 4.2 2.3
##  [559] 3.8 3.6 3.1 4.6 3.3 4.5 4.0 4.2 3.2 4.6 4.8 4.8 4.9 5.7 4.7 5.5 4.8 4.2
##  [577] 5.2 3.6 4.1 4.5 6.3 5.9 5.6 5.0 5.0 2.8 4.8 4.2 4.4 3.7 4.7 4.5 3.9 7.3
##  [595] 3.9 4.3 4.2 3.9 3.9 4.2 4.9 5.0 4.4 3.0 7.1 4.1 2.5 6.1 5.1 5.3 4.3 5.8
##  [613] 3.2 3.2 5.7 5.7 3.6 4.8 6.2 3.1 4.0 5.9 5.2 3.2 4.8 5.0 4.5 6.1 3.5 2.9
##  [631] 3.1 5.0 6.0 5.0 5.1 3.8 3.5 4.9 5.7 4.3 3.9 5.1 3.8 4.7 4.7 2.9 4.8 5.0
##  [649] 5.2 5.8 4.8 6.2 3.8 2.7 5.1 4.0 5.4 5.8 3.1 5.3 4.1 3.4 5.0 5.2 3.5 6.8
##  [667] 4.3 5.8 5.2 5.3 4.1 4.4 5.4 3.5 4.6 4.2 4.6 5.9 5.0 4.4 5.3 4.1 5.8 2.9
##  [685] 4.7 4.7 4.7 5.5 2.6 4.8 4.5 5.2 5.8 3.9 4.4 2.1 3.2 5.7 3.9 4.4 4.0 3.8
##  [703] 4.6 5.7 5.1 5.3 5.0 3.6 3.3 4.3 5.2 3.7 4.8 3.9 4.0 4.1 5.2 4.8 5.2 3.6
##  [721] 4.7 4.2 5.7 4.4 4.0 4.3 4.0 5.1 4.6 5.1 3.4 5.0 5.4 4.2 5.5 2.6 4.6 4.0
##  [739] 5.3 5.0 4.6 5.0 4.8 4.6 6.5 4.2 4.5 5.7 3.4 5.3 3.4 4.8 4.2 3.5 5.6 3.4
##  [757] 4.6 2.2 4.2 4.7 3.4 4.1 4.7 4.2 4.7 5.6 6.5 4.2 6.6 4.5 3.1 4.2 2.4 4.7
##  [775] 3.9 3.9 3.6 2.4 4.8 4.3 4.4 4.3 3.7 3.9 5.1 4.4 5.5 6.4 5.3 3.1 4.5 5.6
##  [793] 5.4 4.6 4.4 2.4 4.5 4.0 5.7 2.7 2.9 5.3 5.4 5.5 4.5 5.5 4.4 3.9 4.0 2.9
##  [811] 5.2 5.1 5.5 2.9 5.0 4.8 4.9 2.9 2.6 5.3 4.9 5.3 5.2 4.5 3.7 5.4 3.4 5.6
##  [829] 3.6 5.0 4.3 5.1 5.5 4.6 4.0 4.8 4.9 3.2 4.6 4.7 4.9 4.7 5.5 3.0 2.7 6.2
##  [847] 5.1 4.0 4.7 4.0 6.1 4.3 3.7 2.8 3.0 5.9 4.3 3.4 5.1 5.0 7.1 4.4 4.8 3.9
##  [865] 3.5 5.2 4.2 3.2 4.1 5.7 5.4 4.8 5.4 3.7 5.7 3.9 4.4 4.3 5.9 5.2 5.3 5.9
##  [883] 3.2 4.9 2.5 4.8 2.8 3.1 4.6 2.6 4.4 3.9 4.6 3.3 4.6 4.3 4.7 5.8 3.5 3.6
##  [901] 5.3 4.4 4.7 3.4 4.7 2.9 5.2 5.4 4.2 4.1 5.5 5.7 5.7 5.1 3.4 4.6 5.1 4.2
##  [919] 3.8 6.4 3.3 4.4 3.8 5.0 4.0 5.1 5.2 5.2 4.2 5.1 4.8 5.1 3.4 4.6 5.5 4.0
##  [937] 5.5 6.1 3.6 4.2 3.3 4.1 6.2 4.0 4.4 4.7 3.4 4.4 5.5 4.2 6.1 5.6 4.4 4.7
##  [955] 3.7 4.6 4.5 3.9 5.3 3.7 3.3 3.0 3.2 3.8 4.8 4.3 4.5 4.8 5.1 2.4 4.1 4.3
##  [973] 3.9 3.9 3.6 2.4 6.0 3.9 4.7 5.3 3.4 6.2 4.9 3.8 4.7 4.7 5.4 3.2 3.6 5.2
##  [991] 3.7 4.8 2.4 3.9 5.0 4.6 4.5 4.4 4.6 3.3
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
##   2.6   6.2
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
##    [1] 5.1 4.2 4.0 4.2 4.9 4.7 6.0 3.8 3.3 5.7 5.8 4.7 5.8 4.2 5.8 3.9 4.2 3.8
##   [19] 4.5 4.9 3.5 4.9 5.3 4.9 3.5 2.6 5.4 3.9 3.8 5.6 4.8 3.1 4.5 4.9 4.5 4.3
##   [37] 4.2 4.5 4.7 3.3 4.5 2.9 5.2 4.1 3.3 2.9 5.7 5.4 5.0 4.3 4.1 4.5 4.1 3.6
##   [55] 3.9 6.2 4.7 6.0 4.6 4.5 4.0 3.8 5.5 5.0 4.6 4.0 2.7 5.2 3.2 3.1 5.3 6.3
##   [73] 5.5 4.1 3.2 4.7 3.6 5.0 5.4 5.4 4.2 3.1 4.7 5.0 4.7 4.8 4.5 3.8 3.7 4.7
##   [91] 4.2 4.3 4.8 5.9 4.3 4.7 5.8 4.9 3.2 4.0 5.4 3.6 3.0 3.9 5.4 3.6 4.8 4.8
##  [109] 4.8 4.2 4.1 6.2 4.3 3.8 3.7 5.3 3.3 2.8 4.7 4.7 3.0 5.2 4.2 6.1 3.5 4.7
##  [127] 4.3 4.3 4.7 4.9 4.0 5.2 5.3 4.6 5.0 5.3 4.4 3.6 3.5 4.7 3.0 4.6 7.0 4.7
##  [145] 2.8 3.7 4.5 3.1 4.7 4.4 5.0 5.0 4.4 4.9 3.9 4.0 5.6 3.8 4.7 4.9 5.8 2.8
##  [163] 4.2 4.5 4.2 3.3 3.1 3.9 4.4 4.5 5.7 5.1 4.4 5.2 3.3 3.1 3.4 4.7 5.3 4.4
##  [181] 3.2 4.4 3.2 4.7 5.6 4.9 5.5 3.1 4.0 4.7 1.8 4.9 4.5 4.2 4.6 5.1 4.2 4.5
##  [199] 3.6 5.9 6.0 5.7 4.4 4.4 3.3 4.8 3.5 4.4 5.6 4.6 4.0 4.8 4.8 3.6 5.1 4.6
##  [217] 4.4 3.7 5.1 5.2 5.7 3.3 4.8 4.9 4.2 5.3 3.3 5.0 4.8 5.0 3.0 5.8 5.5 5.1
##  [235] 6.6 3.9 5.6 4.5 5.5 2.7 2.7 3.9 4.8 3.8 4.2 4.9 5.1 5.0 4.7 5.2 4.3 4.6
##  [253] 6.0 3.7 5.2 3.1 3.9 4.4 3.5 2.0 5.2 3.5 3.9 2.9 3.7 5.9 5.2 5.1 6.9 4.4
##  [271] 4.0 4.1 5.7 3.7 3.5 4.9 4.7 5.6 5.9 3.9 5.0 3.3 7.0 5.4 3.8 4.1 4.8 3.6
##  [289] 4.3 3.1 4.5 4.9 6.2 5.9 5.5 4.0 4.2 3.8 4.6 3.1 3.3 4.5 5.4 3.4 4.8 4.9
##  [307] 4.0 3.1 4.1 5.3 3.6 3.7 3.5 4.1 4.2 4.4 2.5 4.0 4.8 3.6 5.1 3.8 3.5 4.5
##  [325] 4.8 6.1 5.6 5.3 4.9 4.8 3.7 4.9 4.8 5.2 4.4 3.6 5.7 4.8 5.0 3.2 3.9 3.8
##  [343] 5.6 5.0 4.4 5.0 3.8 5.6 5.5 5.8 4.9 6.6 4.3 4.2 4.4 3.9 3.8 6.1 6.1 3.9
##  [361] 5.3 4.7 4.8 3.5 5.2 4.7 4.1 4.6 3.0 3.4 4.1 5.6 4.3 5.0 4.2 4.4 3.5 2.9
##  [379] 5.0 5.7 4.3 5.5 4.5 6.5 4.3 4.3 3.7 4.4 5.9 3.2 4.6 5.6 4.8 5.5 3.4 4.3
##  [397] 5.1 5.0 4.8 4.4 5.1 5.5 5.0 5.7 5.2 4.2 4.3 4.3 5.1 3.6 1.9 4.5 4.6 4.3
##  [415] 5.0 5.6 4.9 4.0 4.6 4.1 3.4 5.6 5.0 5.4 3.3 5.3 3.0 4.9 3.7 4.8 4.7 4.6
##  [433] 4.4 3.4 4.5 6.0 3.8 6.7 4.2 4.8 5.1 4.4 4.0 4.3 5.3 3.4 4.2 4.0 4.3 4.4
##  [451] 3.9 4.6 3.6 4.7 4.6 4.2 3.6 4.0 4.3 4.4 5.7 4.9 6.0 5.8 4.9 5.6 4.6 3.8
##  [469] 4.1 5.5 4.9 4.8 4.4 4.5 4.3 4.9 5.2 4.6 4.3 2.9 4.2 4.2 4.3 5.2 5.5 4.3
##  [487] 4.8 6.0 4.9 4.8 3.8 3.0 4.7 5.4 5.5 5.2 5.1 5.9 3.4 5.7 5.0 5.1 4.0 5.0
##  [505] 3.4 3.2 2.8 5.4 4.1 5.4 2.9 4.5 4.0 4.7 3.6 4.3 4.9 3.8 5.8 6.5 3.6 4.5
##  [523] 4.3 3.7 4.4 5.0 5.4 4.3 5.1 4.6 4.2 2.1 3.7 3.1 4.2 5.0 4.7 5.3 4.8 3.4
##  [541] 3.4 4.5 5.4 4.1 5.1 5.7 4.9 4.9 4.0 4.1 3.8 5.1 3.8 4.3 3.2 4.3 3.1 4.4
##  [559] 5.6 3.6 5.1 4.4 5.8 4.8 5.4 6.3 5.0 3.4 4.7 3.8 5.6 3.4 4.2 4.5 3.4 4.7
##  [577] 3.5 4.6 4.9 4.8 4.1 4.8 6.0 4.1 5.0 5.4 7.5 5.2 4.7 5.0 3.7 5.3 4.8 2.5
##  [595] 5.0 6.4 5.5 4.4 4.5 4.0 5.0 5.9 5.3 4.5 5.2 6.1 2.9 3.6 4.2 4.9 4.6 3.3
##  [613] 4.9 4.5 5.2 4.0 5.8 4.6 4.8 4.5 4.9 4.0 3.7 4.1 4.8 4.7 5.3 3.5 6.4 5.0
##  [631] 4.3 6.3 4.6 4.0 6.2 2.4 5.6 4.2 4.8 4.1 4.1 3.7 5.1 5.5 3.9 4.7 5.1 5.5
##  [649] 3.5 4.9 3.5 6.4 3.2 4.1 5.1 5.3 5.7 2.3 4.1 5.1 5.1 5.8 5.3 4.9 5.1 2.0
##  [667] 5.0 3.8 6.4 4.7 4.8 2.5 4.9 3.3 4.9 3.7 4.1 5.4 4.9 5.3 5.2 4.5 2.4 4.0
##  [685] 4.4 2.5 4.2 4.9 3.6 4.2 2.7 3.7 3.2 5.0 3.7 4.9 3.2 5.3 5.3 4.4 5.9 2.4
##  [703] 6.0 4.3 5.2 4.7 3.2 3.6 5.4 4.0 5.9 4.8 5.6 4.6 4.4 4.1 4.9 4.1 4.2 2.8
##  [721] 5.3 3.3 4.8 3.0 3.7 4.4 4.2 5.4 4.5 4.6 3.4 3.8 6.4 3.9 5.0 2.4 4.8 3.4
##  [739] 5.1 6.1 4.9 2.4 4.9 3.8 5.3 5.5 3.8 4.7 5.5 4.7 5.4 4.2 4.4 4.1 5.5 4.4
##  [757] 5.8 4.1 3.8 3.8 5.5 3.9 4.0 4.6 4.0 3.5 4.6 3.8 6.8 3.8 5.5 4.5 3.7 5.1
##  [775] 2.8 4.9 5.6 5.5 4.2 4.4 5.9 4.5 4.1 4.2 4.4 4.4 5.5 4.6 4.9 3.8 3.4 4.0
##  [793] 5.2 6.4 4.0 4.3 4.1 5.3 3.8 5.2 4.7 3.2 4.6 4.6 5.4 4.4 3.8 5.6 3.4 6.0
##  [811] 5.0 3.8 5.5 5.1 4.6 5.8 4.8 3.7 5.4 4.8 4.4 4.3 3.9 4.5 5.1 5.2 3.9 5.3
##  [829] 5.6 4.7 4.9 5.3 3.8 4.9 4.3 4.2 5.3 5.6 6.3 4.9 3.3 4.9 4.3 5.6 4.5 5.7
##  [847] 4.3 4.3 5.3 4.8 5.9 5.8 4.7 4.0 3.0 3.8 4.1 5.6 3.4 5.2 3.3 5.1 4.4 3.6
##  [865] 4.7 4.5 4.6 5.7 4.9 4.6 4.5 4.3 5.6 3.6 5.0 6.0 5.4 5.5 3.9 5.1 5.7 5.3
##  [883] 5.3 5.2 4.6 4.3 5.5 3.1 4.2 3.0 3.8 3.7 3.6 5.0 3.6 3.2 4.0 5.3 3.8 3.4
##  [901] 4.8 4.5 4.4 3.5 3.5 5.0 3.9 4.8 2.6 5.5 4.6 5.9 4.8 5.1 3.6 3.6 4.1 4.2
##  [919] 5.3 5.2 6.0 4.6 5.3 4.9 5.1 4.0 4.1 3.9 5.4 4.0 4.0 3.8 4.8 4.7 4.8 4.2
##  [937] 5.7 5.8 4.3 4.9 4.9 5.3 4.9 5.2 3.6 3.6 2.8 5.6 4.9 3.7 4.7 4.9 5.8 4.3
##  [955] 5.1 5.7 4.4 3.9 6.0 5.4 4.4 5.2 5.3 3.4 5.7 5.9 6.0 4.4 3.5 4.5 4.1 4.2
##  [973] 6.2 4.4 2.1 5.3 2.8 3.7 3.7 3.7 5.7 3.2 5.6 3.5 3.6 6.0 3.8 3.7 5.5 4.5
##  [991] 4.6 4.5 6.0 5.8 4.1 4.2 4.3 3.6 5.1 2.6
## 
## $func.thetastar
## [1] 0.027
## 
## $jack.boot.val
##  [1]  0.51117479  0.42215385  0.29002695  0.17598870  0.06950147 -0.07804154
##  [7] -0.10802292 -0.19732938 -0.33537604 -0.47753165
## 
## $jack.boot.se
## [1] 0.9242802
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
##    [1] 3.6 4.3 5.2 3.6 4.2 4.2 4.1 4.8 3.2 3.9 5.6 4.1 4.4 3.5 4.7 4.3 4.5 4.3
##   [19] 3.8 4.5 5.7 4.6 4.6 4.3 4.5 5.2 4.6 5.1 4.1 4.6 4.3 4.6 6.1 6.0 5.6 4.4
##   [37] 5.2 5.3 5.1 4.9 4.3 3.8 5.3 4.4 4.6 5.2 5.1 6.8 5.9 4.2 4.7 3.2 6.4 5.2
##   [55] 3.3 4.1 5.1 3.9 5.9 5.0 4.4 6.2 6.2 3.5 3.1 4.7 6.3 4.9 3.2 5.4 3.8 4.5
##   [73] 2.6 3.4 4.4 4.7 5.5 4.9 4.6 3.9 3.6 5.2 4.1 2.8 3.6 3.9 2.8 3.9 3.4 4.6
##   [91] 4.6 4.6 3.5 3.8 2.9 3.7 5.4 5.0 5.5 3.6 3.7 5.5 3.6 3.2 4.7 5.0 4.4 4.5
##  [109] 4.0 5.4 4.8 4.5 4.1 4.5 4.9 3.9 4.3 2.0 3.4 3.8 3.7 3.1 4.8 3.9 3.4 5.7
##  [127] 4.0 5.2 3.5 2.1 4.5 4.9 5.8 4.3 5.6 3.8 3.1 6.1 4.1 6.0 4.2 3.9 6.0 5.0
##  [145] 4.6 5.0 4.7 6.0 4.9 4.9 5.9 4.3 5.2 4.3 4.2 5.6 4.8 3.2 5.4 3.8 5.7 4.9
##  [163] 3.6 4.0 5.1 3.4 5.5 4.3 4.4 4.3 4.9 3.4 4.1 3.9 4.5 5.5 4.8 3.8 5.4 5.0
##  [181] 6.4 4.6 4.9 3.8 4.0 4.4 4.9 2.3 4.7 3.8 4.6 5.8 5.4 5.5 4.6 6.0 4.8 6.7
##  [199] 4.8 5.1 4.3 5.8 6.2 6.4 4.6 4.1 4.3 3.5 5.2 6.0 4.5 4.9 3.1 3.7 3.1 3.3
##  [217] 3.9 5.4 4.8 4.4 6.4 4.5 3.2 3.5 4.2 5.5 4.5 2.0 4.0 3.7 3.5 3.5 4.2 5.1
##  [235] 4.7 4.5 3.7 4.9 2.6 4.6 6.9 2.9 4.6 4.9 4.3 6.6 5.6 3.5 6.5 5.4 4.9 4.5
##  [253] 3.4 4.7 5.7 3.2 2.7 4.4 4.3 4.6 4.3 3.6 3.6 5.4 4.6 4.7 4.2 3.8 5.8 4.1
##  [271] 3.3 3.4 3.6 5.1 3.4 5.2 5.2 5.5 5.6 4.9 4.5 3.3 7.0 2.4 4.7 4.3 3.1 4.4
##  [289] 4.2 3.5 5.1 6.1 4.6 6.5 3.0 4.7 4.1 4.6 4.7 5.7 6.3 4.7 6.1 3.6 4.7 5.0
##  [307] 3.5 6.3 4.4 3.3 5.1 6.5 4.0 4.4 2.7 5.4 4.3 4.4 4.2 4.4 6.0 4.9 4.6 4.3
##  [325] 4.3 2.7 3.4 5.3 3.8 4.6 4.6 4.0 4.6 4.7 4.8 4.6 4.7 5.1 3.2 3.3 3.9 4.0
##  [343] 6.1 3.6 4.1 5.3 4.4 4.3 3.9 4.5 3.6 5.4 3.9 4.7 3.7 5.2 4.5 4.5 3.6 4.8
##  [361] 4.9 4.7 3.5 4.1 4.5 4.5 5.3 4.2 5.4 4.2 4.2 4.4 4.9 5.7 2.9 6.9 5.0 4.2
##  [379] 4.5 4.5 4.0 5.5 3.3 3.6 5.6 4.5 4.5 3.2 4.4 4.0 4.7 4.7 4.4 5.8 4.8 5.3
##  [397] 3.3 5.0 3.4 4.1 3.7 3.9 5.9 5.3 2.8 3.4 3.8 4.5 3.1 3.5 3.3 5.4 4.5 3.1
##  [415] 5.7 6.0 6.8 5.0 4.8 4.5 5.1 3.6 3.3 3.9 4.3 4.3 3.4 5.8 5.2 5.2 3.9 6.0
##  [433] 3.2 3.7 5.9 3.4 3.9 6.1 4.3 2.8 4.8 4.7 4.1 4.9 6.2 5.7 5.4 4.4 5.0 4.1
##  [451] 4.6 5.4 4.0 4.4 4.7 3.3 5.5 4.4 4.5 3.4 5.2 4.4 5.7 3.9 3.5 4.8 4.1 4.2
##  [469] 3.3 4.0 5.1 5.1 5.6 5.2 5.4 2.9 6.2 5.1 3.8 4.3 4.9 2.7 4.6 4.5 4.2 4.3
##  [487] 4.5 6.5 3.9 5.1 4.6 3.6 5.5 5.1 1.7 4.6 4.9 5.1 5.0 4.3 4.8 4.8 5.2 4.2
##  [505] 5.6 5.0 2.2 5.0 5.5 3.8 3.8 4.5 3.5 4.7 4.1 4.8 3.3 2.4 4.4 3.7 4.4 5.5
##  [523] 4.8 3.3 4.7 4.0 3.8 4.1 5.0 5.1 5.8 4.8 4.9 4.4 5.1 4.4 5.0 4.0 4.1 6.1
##  [541] 3.6 3.9 5.4 4.7 5.7 4.3 5.4 3.9 5.7 3.5 4.5 4.4 4.7 3.7 5.4 5.3 6.0 4.1
##  [559] 6.4 3.6 4.4 4.7 5.3 4.7 6.2 4.5 4.5 6.5 4.5 4.5 4.1 4.0 4.9 3.3 5.2 3.9
##  [577] 3.5 5.7 6.4 5.7 3.8 3.8 2.7 3.4 4.5 4.3 5.1 4.8 4.4 3.4 4.3 6.0 6.2 4.1
##  [595] 3.5 3.8 4.6 3.1 4.6 3.7 3.0 5.4 4.2 2.3 3.6 2.2 4.4 5.1 5.5 3.5 5.6 4.3
##  [613] 5.2 3.7 4.4 5.9 3.6 3.9 4.0 4.3 4.1 4.1 3.0 4.9 5.0 6.0 4.9 3.9 2.6 4.4
##  [631] 5.2 6.4 4.3 5.3 3.8 4.0 4.2 3.3 4.4 5.0 5.9 5.3 6.1 4.0 4.6 5.5 3.2 4.0
##  [649] 3.5 3.9 4.8 5.1 4.2 4.0 4.4 5.6 3.1 3.4 3.5 4.6 3.6 4.5 4.8 3.6 5.0 3.0
##  [667] 3.1 2.8 3.8 4.8 4.4 4.1 4.2 6.2 4.2 4.4 4.0 4.4 5.7 4.9 4.9 4.4 6.4 5.4
##  [685] 4.2 4.8 3.2 4.9 3.3 5.5 4.3 4.7 5.2 5.2 6.0 5.0 5.4 5.8 4.5 5.1 4.4 4.3
##  [703] 4.0 4.0 5.1 4.1 4.9 5.5 5.0 4.9 5.8 4.7 4.6 4.9 4.4 4.7 4.2 3.6 4.7 6.1
##  [721] 5.5 4.9 3.3 4.4 6.0 3.4 3.7 3.3 4.9 5.8 6.7 4.5 5.4 4.1 4.4 5.4 4.3 3.8
##  [739] 6.1 5.7 4.5 5.1 5.9 5.8 4.0 3.6 4.5 4.7 6.5 4.7 4.9 3.9 4.1 4.8 3.6 3.4
##  [757] 4.2 4.3 4.5 4.1 5.4 4.5 2.4 3.5 5.9 4.4 5.1 5.3 2.9 4.5 4.8 4.6 3.4 6.7
##  [775] 3.8 3.3 7.0 4.5 4.1 4.8 6.1 2.3 5.3 3.4 5.3 4.2 5.0 4.9 6.1 4.4 3.7 4.5
##  [793] 4.4 4.6 6.6 3.2 5.4 3.8 3.9 4.2 2.7 5.5 5.3 3.6 5.6 3.9 4.3 4.5 5.2 3.8
##  [811] 3.7 3.8 3.8 4.1 5.8 7.1 5.4 4.7 4.2 4.9 4.0 6.3 5.7 4.8 4.2 5.3 4.1 4.8
##  [829] 4.1 4.3 6.1 4.0 4.2 5.6 4.4 6.5 3.7 4.9 5.3 4.0 3.6 3.2 2.3 4.9 4.4 4.4
##  [847] 2.8 3.4 4.4 5.0 3.2 5.2 4.5 4.0 5.4 5.1 4.4 3.7 3.8 4.7 4.8 4.1 3.9 3.7
##  [865] 4.6 4.6 4.9 4.4 4.4 5.5 3.1 4.8 5.0 6.0 3.7 2.3 3.2 5.8 4.5 4.5 4.7 4.0
##  [883] 6.3 4.2 5.1 3.5 5.1 6.1 5.5 3.1 6.5 3.9 4.2 5.5 4.5 4.5 4.3 5.1 2.8 4.8
##  [901] 3.4 4.5 4.9 6.0 3.3 4.7 3.3 3.6 4.2 3.9 4.9 5.8 2.4 5.3 5.6 4.3 5.0 3.7
##  [919] 4.1 6.0 4.4 5.8 5.1 5.5 3.9 1.5 5.5 4.9 4.9 4.1 4.6 3.9 3.2 3.7 4.3 3.2
##  [937] 4.6 5.2 4.5 6.0 4.2 3.0 3.8 4.7 5.9 5.2 4.4 5.3 3.3 4.4 4.5 3.5 3.5 3.7
##  [955] 5.2 3.5 5.0 6.1 4.5 4.7 4.5 5.2 5.5 3.7 4.9 4.5 4.3 5.7 5.2 4.3 7.1 2.4
##  [973] 3.8 6.7 3.9 6.0 4.5 3.2 5.3 5.7 4.9 4.0 2.5 5.3 5.2 4.3 3.6 2.5 4.2 4.8
##  [991] 5.4 5.1 4.7 4.7 5.9 4.7 4.3 4.1 6.3 2.9
## 
## $func.thetastar
## 72% 
##   5 
## 
## $jack.boot.val
##  [1] 5.50 5.60 5.22 5.40 5.20 5.10 4.80 4.70 4.60 4.50
## 
## $jack.boot.se
## [1] 1.110389
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
## [1] 0.8686341
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
##   6.207697   9.166483 
##  (2.704855) (4.160071)
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
## [1] 0.2852967 0.4702354 1.2131310 0.1419910 0.7088664 0.4370629
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
##    [1] -0.2281600535  0.7526121156  0.3939417583  0.6775455697  0.8352716907
##    [6]  0.6825368006  0.9705013694  0.7067122693  0.6261517734  0.9716859874
##   [11]  0.5113976099  0.6466320606  0.6620855929  1.7194246161  0.5894667928
##   [16]  0.2935179992  0.5568293571 -0.0780601793  0.4906979423  0.5406743598
##   [21]  0.7855487773  0.3394269216  0.5043831460  0.6765422539  0.5022681294
##   [26]  1.6669431194  1.1535408900  0.7590905105  0.7146204699  2.0118967447
##   [31] -0.3515523289  0.6618121379  0.9139955202  0.8811346861  1.1578297611
##   [36]  0.0658074540  0.1305522222  0.6207732756 -0.1892755306  0.5877573151
##   [41]  0.8674189654  1.4177188837  2.1101978839  2.1962180579  0.1758235964
##   [46]  1.1165863133  0.7501507248 -0.4417990950  0.0685091602  0.6488003943
##   [51]  0.5109836142  0.1180237300  0.9643172434 -0.3895127054 -0.3837864794
##   [56]  0.7187723650  0.8210636169  1.2647437070  0.7501555697  1.6121740177
##   [61]  0.3147384170  0.0343700400  0.6225433441  1.2753086922  1.3307561757
##   [66]  1.2346032717  1.2108755264  0.3656786321  0.3032837254  1.7788490189
##   [71]  0.2562407344  1.3393311516  0.2603298902  0.5164243453  1.3223373028
##   [76]  0.6396587411  0.8906398835  0.6114575786  0.6065598963  0.1943373342
##   [81]  0.4600363688  1.2169451997  0.0401479891  0.3631686853  0.9184259908
##   [86]  0.8519566519  0.7621593756  0.1501744736  2.0895941943  1.0753808294
##   [91]  1.0317605077  0.5723826114  0.6837510363  0.6546927318  1.2789988062
##   [96]  0.2205653538  0.7943600198  0.6866155824  0.6673815576  1.1557114148
##  [101] -0.1798946282  2.4767800488 -0.0307045457  1.2625486035  0.4300399528
##  [106] -0.4709958099  0.3680773519  0.4431202130  0.2313663397  0.9419527991
##  [111]  1.8410652959  0.7044753003  0.5245498291  0.9397258587  1.1446117143
##  [116]  0.3291187760  1.2708146499  1.7679879787  1.2669700461  0.9487993428
##  [121]  0.5064083849  0.4160260023  0.6278982319  0.5162647927  0.9495338686
##  [126]  0.7511571215  1.2700344914  1.5617686574  1.0713732592  1.2031768066
##  [131]  2.0356121075  0.6554784955  1.2068498761  0.4564850393  1.2021791213
##  [136]  0.3298950983  0.4645144643  0.2012148496  0.1624441546  0.8921887438
##  [141]  1.5865902840  1.3006948006  1.0270519491  0.3735661791  0.7232387545
##  [146]  0.4527921699  1.3052142636  0.4626205927  0.8473883349  1.4627154774
##  [151]  1.4522446629  0.8042909139  0.1947456820  0.8599680847  1.3941586686
##  [156]  1.5660768245  1.2659879248  0.1488339655  1.0875672445  0.5291562054
##  [161]  0.2253021065  0.4132678123  1.2752311003  0.9180365400  0.8788009037
##  [166]  0.4288635764  0.7073377753  0.3559836907  1.6282200538  1.2547890644
##  [171]  0.6416947490  0.5962471119  0.8542091734  0.3932584422  0.3548475014
##  [176]  0.8086504066  0.5811837121  1.4478089350  0.8159025527  1.1075588171
##  [181]  1.4221880929  0.3671452861  1.2203629269  0.8306048396  1.4307774315
##  [186]  2.3005123950  0.2500395580  1.1575802592 -0.0112827746  0.4725101520
##  [191]  0.6433667122  1.1760967231  0.2026278669  0.7139409905  1.0911742845
##  [196]  0.5618755203  0.0522646767  1.2492537452  0.7163320107  0.2096831523
##  [201]  0.3410603900 -0.1105850684  0.4921911745  0.1220311209 -0.8399459291
##  [206]  1.2660715749  1.3291380649 -0.1326619650  0.7403328565  0.8755446107
##  [211]  0.5171450295  0.4378889809  1.8006387756  1.3332050502  0.7435933242
##  [216]  0.9743526599  0.4090463067  1.1030995349  0.4265007759 -0.2245169158
##  [221]  0.8563118296  0.1771469152  0.9657809488  2.5619605307  0.7973170926
##  [226]  1.2605648879  0.9555823377  0.4426873425  0.7988519014  0.8321070794
##  [231]  1.2923877580  1.7612031630  1.1399629075  0.1976510628  0.1514469943
##  [236]  0.4321503250  0.1188862360  1.4612040892  0.6820619341  2.2314632917
##  [241]  0.0892243021  0.0988298739  0.7538580205 -0.1096745278  2.4389616600
##  [246]  0.4927440844  0.4204806061  0.0868288259  0.6369764987  2.0784977131
##  [251]  0.4232817103  1.6581514128  0.8051164854  0.7216551206  0.8935541243
##  [256]  1.3734777719  0.4626031909  0.8501866645  0.7484748611  1.0956887740
##  [261]  0.3417034884  1.2622517055  0.8685041838  1.3501093723  0.6265856540
##  [266]  0.8541897969  1.2597511415 -0.6721621468  0.2627412113  0.7386246733
##  [271] -0.0636839037  0.5116578369  1.9922623221  0.4704047312  0.0445884196
##  [276]  0.3186697093  0.4052958283  1.4213519534  0.2159540354  1.6868721463
##  [281]  0.0640512307  0.5246386364  1.0848396491  0.8724021220  0.5120400829
##  [286]  1.8745376392  0.6285615122  1.4481205550  0.4614807275  1.2306407410
##  [291]  0.3198288526  2.0425843622  0.7441030164  1.4434105523 -0.4133267155
##  [296]  0.2665412345  0.8997159929  0.1237091697  0.5344489086  0.6784965856
##  [301]  0.1160353608 -0.2682096687  1.1726761603  0.1895355268  0.9609370207
##  [306]  1.4043610359  1.7813010866  0.0039500388  2.4127801290  1.8027859012
##  [311]  0.4260053212  1.0413731928  0.8110422182  1.2405467361  0.9014890414
##  [316]  0.8687451754  0.4506628038  1.5710254697  2.5540162352  1.1218924114
##  [321]  0.5341536457  1.0756867317  1.1395026301  1.8124397007  0.8606812293
##  [326]  0.7628595313  1.0860187653  0.7234905992  0.3432623892  0.9826864789
##  [331]  0.4859812510  0.1705200420  1.5796575435  0.8495410331  0.2972481281
##  [336]  0.4877956940  1.0657567939  0.0252699800  0.8023052577  1.1739367714
##  [341]  0.7783335492  0.5591213389  0.6155380253  0.7235541472  0.5068120229
##  [346]  1.2311060599  0.5243013937  0.3852953382  0.8707080161  1.2531947155
##  [351]  0.2614999484  0.0551959766  0.9865897460  1.3240080575  1.0596231372
##  [356]  1.9973536589  0.5292301489  0.4455868077 -0.6437573801  0.8209455853
##  [361]  0.4532887189  2.1756692508  0.7411018403  1.9866448448  1.2896906248
##  [366]  1.3555117351  1.2844094801  0.4837612986  0.0393665124  0.5262623640
##  [371]  0.9543746486  0.7655584793  1.6500013556  0.7626620866  0.5356563276
##  [376]  0.3710330246  0.4257788407  0.2233147059  0.3707360233 -0.0892963222
##  [381]  0.8630944284  1.1975093190 -0.4079508673  0.7433158810  0.4304759648
##  [386] -0.4506450709  0.9924469367  0.2339678038  0.9075452174  0.5254592028
##  [391]  0.9761169417  0.7998522893  1.1501993632  1.0014431782  0.2783255313
##  [396]  0.4686558832  1.6647452086  1.3642809848  0.8885383942  0.9758172708
##  [401]  0.5199288318  0.3626334486  0.8100537547  1.4538472319  0.7663200618
##  [406]  0.3563497376  1.2425330028  0.9649249109  0.2895919301  0.9682739063
##  [411]  0.7388995134  2.2256947562  0.9020848616  0.9562498311  0.6811621607
##  [416]  0.5604165528  0.4965444826  0.5459194207  0.6710264032  0.7658081096
##  [421]  0.3355929353  1.2019420835  1.0941290416  0.4537992769  0.7409324870
##  [426]  0.6984556309 -0.0649348611  0.2953031059  0.9770791809  0.2102922759
##  [431]  1.4361460137  0.3892185466  1.1014054205  0.9129107363  0.5485327741
##  [436]  0.2952204724  2.2558795729  0.7220669874  0.8100820787  2.5094831009
##  [441]  1.4030103106 -0.0439279369  1.3928531471  1.3142173520  0.8837340797
##  [446]  0.9667016647  1.4585662289  0.5292301489 -0.0852013070  0.9733067742
##  [451]  0.5889449955  0.4587368756  0.9827451502  0.7117609803  1.6149453127
##  [456]  0.5123368912  0.2314442207  0.3167440401 -0.4302258609  0.3382579558
##  [461]  1.0754931902 -0.3703191285 -0.2600029442  0.8874981310  0.7457041016
##  [466]  1.4400362445  0.8866969517  0.5274954734  0.3892128320  0.3550259970
##  [471]  0.4145226647  0.2925958069 -0.4352965882  1.1790790587  0.9801451691
##  [476]  0.0525611016  1.3119937121  0.9360077014  0.4273512816  1.6665660822
##  [481]  0.9166258407  1.3893880180  0.6877998564  0.9643172434  0.4129930087
##  [486]  1.1102696282  0.6344952500  1.5824246199  2.5630157696  0.9437609038
##  [491]  0.2715194852  0.9576031975  0.4956661421  1.1443359237  0.6774771825
##  [496]  1.3720630821  0.2614083508  0.8478599919  0.2956589162  0.9119597788
##  [501]  0.4945732775  1.2555271383  0.4610981937  0.4947660546  1.5199020757
##  [506]  0.9897940386  0.8331324223  0.5536877552 -0.1758251141  1.1636918544
##  [511]  1.0428996289  1.0223324062  0.3534675589  0.7907204371  0.0547707314
##  [516]  0.5680829991  0.7587434189 -0.2107643899  0.8892223335  0.9775843384
##  [521]  1.3148144526  0.8295133415  0.5580044399  0.8056611963  0.7075364155
##  [526]  0.7911820301  1.6294896499 -0.1802865693  0.9296094714  1.9507424372
##  [531]  1.1315446511  0.4213691485  1.1411004740  0.7681541125  0.2939836585
##  [536]  0.4457260321 -0.4966086030  0.4990838751  1.6752490367  0.7921341556
##  [541]  0.3222093300  0.1570459461  0.6213029544  0.3152927677  0.5863785494
##  [546]  0.6202715974  0.9449104720  0.4853607471  2.2372398705  0.7632778639
##  [551]  0.7656161009  0.1499529766  1.2145775914  0.5389041452  0.6824854723
##  [556]  1.2397616916  0.7957399223  0.2243179878  0.6826689011  0.8738091767
##  [561]  0.6078924822  1.5453747763  1.1834262701  0.6998677261  0.5124715088
##  [566]  1.0053431724  0.3912537630  0.4062885262  1.3923487625 -0.0956754284
##  [571]  0.0220948297  0.2772316929  1.4020183120  0.4286946607  1.7476828588
##  [576]  0.8573062287  0.3609495030 -0.1226116133  0.9311884358  0.5224991344
##  [581]  0.9758172708  0.8448149884 -0.3196619576  2.0244252081 -0.0467024779
##  [586]  1.0222680721  1.3501093723  1.5736312359  0.6782153324  0.0461082115
##  [591]  1.5773191891  0.4770213665  0.3518270403  0.5543113977  0.8280789306
##  [596]  0.5624199781  0.2413880104  1.0377558769  0.0654489521  0.9662185902
##  [601]  1.0028822325  0.0630099240  0.0867201125  1.1058840497  1.5765982426
##  [606]  1.8218359951  1.2292439012  0.7995219690  0.2764127793  0.8429476798
##  [611]  0.4903961400  0.3718703447 -0.6120984787  0.7072114686  0.9458476066
##  [616] -0.2156816450 -0.3216919391  0.5036739008  0.1238713368 -0.4716262720
##  [621]  0.6931174202  0.0852371600 -0.3699907821  0.4871693545  1.0840358429
##  [626]  1.9726144271  0.8401473466  1.2409511489  0.3906677472  0.4852980629
##  [631]  0.3841108311  0.2233147059  1.5992094491  0.9144587351  0.5174684690
##  [636]  0.4389985677  1.3959631934  0.4067095835  0.8606812293 -0.2473901249
##  [641]  1.1296478646  0.4141057691  1.0837024783  1.0756867317  1.0055621467
##  [646]  1.0046670802  0.8922091973  0.6184027639  0.2927061141  1.1480290611
##  [651]  0.4144514079  0.3669209558  1.1606523399  1.3956959643  0.7535069242
##  [656]  2.1207926568  0.5399267761  0.7610945940  0.7761235336  0.4103438910
##  [661]  1.7699300739  0.3994503199  1.4076396874  0.7825049305  0.7703968064
##  [666]  1.1177791306  0.9657701414  1.1865713525  1.2253066958  0.3995878636
##  [671]  0.0828293571  1.1146738508  0.2571634400  0.1989800519  0.3452616536
##  [676] -0.2895047196  0.6933610851  0.8608963841  0.8357846478  0.6445546136
##  [681]  0.7018089517  0.4887669184  2.0629250859  1.1350761938  1.3707166089
##  [686] -0.4019586296  0.3049858051  0.7594153099  0.7568363268  0.2982497863
##  [691]  0.3543155853  1.4213496803  0.9647593680  1.3246385559  0.5140766829
##  [696]  0.9152225582  0.0899665393  0.9623301963  1.2829077074  0.9611292930
##  [701]  0.6114575786  0.8330498722  0.7501555697  0.7337699418  0.4205654107
##  [706]  0.2839940219  0.3719713234  0.3884208260  1.5706321430  0.6969310255
##  [711] -0.1453085749  1.6765867287  0.0236267519 -0.0941409952  1.4510650042
##  [716]  0.2044183654  0.2875180333  0.3159828897  0.6283609615  1.5052409239
##  [721]  1.4760470739  1.0442527967  1.9105312062  0.5810917068 -0.1867603420
##  [726]  1.1268786424  1.3223633941  1.4095805360 -0.1963633668  0.9711963437
##  [731]  1.2844418516  0.5241059984  0.3581003672  0.3561026144  0.3268524478
##  [736]  0.1102416415  1.1241414666  1.7156163684  0.3674207972  0.1734901187
##  [741]  0.0536784697 -0.0003414102  0.6937101343  0.3962362020  0.5113316369
##  [746]  0.6241837446  0.6182524187  0.5250425117  0.4853500690  1.1272465340
##  [751]  0.9271584442  0.8975339032  0.5377043922  0.0975950112  0.7575595393
##  [756]  1.1390152035  0.7246136594  0.2257003239  1.5944424738  0.4568254561
##  [761]  0.7623249562  0.6251374569  0.1702716219  2.0936857150  0.8868047361
##  [766]  0.3134503208 -0.0036663956  0.4053645786  0.5141434491  0.5069443175
##  [771] -0.1242465259  0.5015323788  0.4424911918  0.1405016927  1.6123263528
##  [776]  0.9633019076  1.2378699742  0.8947930085  1.4303715510  0.8866969517
##  [781]  2.0667463813  0.4319684741  0.8546344841  0.4929329548  0.8394949631
##  [786]  1.1825481891  0.2971487428  0.8495410331  1.1140344127  0.8821965990
##  [791]  0.4515379378  1.8695348136  0.4587855909  0.7515588889  0.2362311696
##  [796]  1.3792758131  0.5404896561  0.6774809779  0.9865897460  1.3820500126
##  [801]  0.7348450337  1.5003208831  1.2605648879  0.9614984829  0.6603678605
##  [806]  0.5358898617  0.4910529178  0.9549456861  0.7066911871  0.8543926249
##  [811]  2.2813379383  0.7491843853  1.1921960147  0.7027765407 -0.0188724590
##  [816]  0.4955165490  0.7013604061 -0.0568365558  0.9858876978  0.9356509475
##  [821]  1.7906166415  0.9681149750  0.3300620126  0.8503479438  1.3045002286
##  [826]  0.9066272750  0.6787316072  0.7005588624  0.9690264816  0.9728767413
##  [831]  1.0320023194  0.5024853501  0.5293370261  1.2057591646  0.4531504658
##  [836]  0.0228074000 -0.1235198545  1.0999962861  0.6339824462  0.9227555534
##  [841]  0.5168714652  1.1651790033  0.8210102259  0.6972244509  0.8353002483
##  [846]  1.3052711628  1.4854627767  0.1401701308  0.6805251677  0.8311082074
##  [851]  0.5267046789  1.3941902178  0.5514552516  1.5975284257  1.3114184999
##  [856]  0.2970307321  0.2196380919  0.1806374713  0.6446902496  0.5165718960
##  [861]  0.3698648569  0.7821634980  0.5514552516  0.7112420514 -0.0787565803
##  [866]  1.2009427594  1.0570364077  0.1360849983  1.0848396491  0.6890322275
##  [871]  1.2752311003  0.7774004714  0.4976865722  1.8066657178  0.3957873669
##  [876] -0.0112676987  0.0692007456  0.2173996826  0.9018151160  0.8130993136
##  [881]  0.6815370902  0.9659944769  1.3370835443  1.2353256406  0.3550259970
##  [886]  0.9539098981  0.2266697638  0.0324667738  0.7420186024  0.9932799716
##  [891]  0.3314542230  0.2106863615  0.0821379040  1.8082746348  0.7122216548
##  [896]  0.5310800696  0.5889449955  0.9528400410  0.7981023893  0.4440415152
##  [901]  1.0362451076  0.8352716907  0.0576148837 -0.4990951035  0.1515414079
##  [906]  0.0237414302 -0.3698023002  1.0882605647  1.1312638436  0.8764791199
##  [911]  0.2543826609  0.5850618095  1.1648751223  1.3083822388  0.3260962552
##  [916]  0.6070663627  1.1247798378  0.9367447666  0.8943235434  0.3478527915
##  [921]  1.1581831685  1.1032073699  0.2940443321  0.1998272639 -0.4084185812
##  [926] -0.1413870359  1.3006948006  0.9564638452  0.7164270201  0.4079851270
##  [931]  0.0620821249  0.4827095532  0.5051030780  0.2817762523  0.7070033106
##  [936]  0.5045124099  0.9666702435  0.5849454876  1.7783924127  0.5181166270
##  [941]  1.4204962077  0.0592309737  1.9462865808  0.5257915951  1.4033621356
##  [946]  0.9351646788  0.6730735639  0.6077514424  0.2084404905  0.3177866674
##  [951]  0.3321987462  1.5700150559  0.2960550030  2.1924228712  0.4673215024
##  [956]  0.8346851678  0.6700508577  2.3890054364  0.5200476301  1.2178417618
##  [961]  0.0525611016  0.6045524782  0.5389555945  1.5448552482  0.9717411329
##  [966]  1.4626687154  2.1017735691  0.8213469765  0.6366670891  0.6120494589
##  [971]  0.5478985774  1.1568285755  1.0636504801  0.1553088078  0.7377706079
##  [976] -0.4837220349  0.6352836519  0.3440085332  1.2880145097  1.3837601599
##  [981]  1.2138769562  1.6682297941  1.3612087240  0.8570906413  0.4424911918
##  [986]  1.1234582788  0.6821864733  0.6282233879  1.2110998276  0.6891976904
##  [991]  1.1216523121  0.4920145188  1.2341317795  1.2420875758  0.8539247894
##  [996]  0.9700225320  0.6750151786  0.9577899526  0.7218871113  1.2462680894
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
##   0.67722034   0.28862343 
##  (0.09127074) (0.06453565)
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
## [1]  0.52098331  0.73301437 -0.07399979 -0.12712176 -0.08375173 -0.29291740
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
## [1] 0.0014
```

```r
sd(results2$thetastar)  #the standard deviation of the theta stars is the SE of the statistic (in this case, the mean)
```

```
## [1] 0.8936983
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
## t1*      4.5 -0.02762763   0.9072506
```

One of the main advantages to the 'boot' package over the 'bootstrap' package is the nicer formatting of the output.

Going back to our original code, lets see how we could reproduce all of these numbers:


```r
table(sample(x,size=length(x),replace=T))
```

```
## 
## 1 2 3 4 6 9 
## 1 2 3 2 1 1
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
## [1] 0.0083
```

```r
se.boot
```

```
## [1] 0.9242919
```

Why do our numbers not agree exactly with those of the boot package? This is because our estimates of bias and standard error are just estimates, and they carry with them their own uncertainties. That is one of the reasons we might bother doing jackknife-after-bootstrap.

The 'boot' package has a LOT of functionality. If we have time, we will come back to some of these more complex functions later in the semester as we cover topics like regression and glm.

