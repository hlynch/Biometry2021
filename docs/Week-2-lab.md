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
## 0 3 6 8 9 
## 2 3 1 3 1
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
## [1] 0.0402
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
## [1] 2.710505
```

```r
UL.boot
```

```
## [1] 6.369895
```

Method #2: Simply take the quantiles of the bootstrap statistics


```r
quantile(xmeans,c(0.025,0.975))
```

```
##  2.5% 97.5% 
##   2.8   6.4
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
##    [1] 3.4 6.2 3.4 4.1 5.7 4.3 3.6 3.8 3.3 3.6 5.0 5.0 5.3 5.4 5.0 5.2 4.6 3.4
##   [19] 3.4 3.8 6.7 5.5 4.7 5.5 4.7 4.8 6.0 5.0 5.2 5.4 4.7 5.0 4.4 2.6 5.8 4.3
##   [37] 4.8 5.2 4.4 4.0 6.0 4.5 3.2 6.4 3.9 4.3 5.5 4.2 3.8 3.4 4.3 3.4 3.3 4.6
##   [55] 5.2 3.3 3.4 5.3 5.4 4.4 5.9 3.8 5.5 4.0 4.3 4.1 3.6 5.7 3.5 5.4 4.4 4.8
##   [73] 4.7 4.6 4.2 5.0 4.6 3.7 3.6 5.5 5.4 3.8 4.8 5.6 4.8 2.9 4.5 3.3 3.8 3.1
##   [91] 3.4 5.7 5.3 4.6 5.9 3.9 2.8 4.0 4.1 5.1 4.9 5.0 4.5 6.1 4.6 4.8 4.1 4.4
##  [109] 5.3 4.1 4.0 3.8 3.8 4.1 4.1 5.2 5.3 4.7 5.4 3.2 4.0 5.7 3.9 3.5 4.9 4.0
##  [127] 6.1 3.3 5.6 4.0 3.2 5.1 5.2 4.8 4.9 6.5 4.1 5.8 6.0 5.5 3.4 3.9 5.8 4.8
##  [145] 4.5 5.1 5.0 5.1 5.6 3.9 4.0 4.4 4.1 3.5 3.2 5.3 6.8 5.5 4.2 4.5 5.0 5.8
##  [163] 4.0 4.8 4.9 4.7 3.5 3.7 4.6 4.5 4.0 4.8 3.5 4.2 4.4 2.3 2.9 4.4 5.0 4.7
##  [181] 3.9 5.8 5.4 5.6 3.3 3.0 3.9 5.2 5.8 5.2 3.5 5.0 4.6 4.9 4.6 4.8 5.5 4.7
##  [199] 5.3 6.2 4.5 3.4 4.6 4.0 4.7 2.2 2.9 4.5 3.9 4.1 4.2 4.2 3.4 5.2 5.2 4.5
##  [217] 4.8 4.3 4.3 5.4 5.5 4.7 4.3 3.9 6.1 3.6 3.7 3.6 3.3 3.2 4.3 5.7 3.9 4.4
##  [235] 4.4 3.5 3.1 4.0 4.6 4.9 4.7 5.4 4.6 4.6 4.7 3.2 2.5 4.0 4.2 5.0 5.6 4.8
##  [253] 4.8 4.0 4.5 4.7 2.8 3.0 5.1 4.1 4.3 5.3 3.7 4.5 3.2 5.2 4.6 4.7 6.2 4.6
##  [271] 5.1 5.6 4.5 6.1 5.5 3.8 4.0 4.2 6.0 4.1 4.6 4.5 4.2 6.9 2.6 4.6 3.8 0.5
##  [289] 3.5 4.3 5.0 3.9 3.6 4.2 3.6 3.2 4.6 4.1 3.5 3.6 5.4 4.8 3.7 4.6 4.0 5.7
##  [307] 5.0 3.9 4.1 5.0 3.7 3.4 5.2 5.2 4.8 3.4 5.7 6.1 4.7 5.3 5.0 5.7 3.6 4.8
##  [325] 4.9 4.8 4.1 5.3 5.3 4.7 3.9 4.6 5.4 4.2 4.5 5.5 3.5 4.9 4.9 3.7 5.5 4.6
##  [343] 4.7 4.9 4.9 4.4 4.6 3.3 4.6 3.5 5.3 5.2 4.1 3.7 4.8 4.5 5.0 3.0 3.7 4.9
##  [361] 4.9 5.5 5.5 3.6 2.9 2.8 4.6 5.3 5.8 2.7 4.8 4.1 4.7 5.2 3.5 4.7 3.5 6.1
##  [379] 3.5 5.0 4.6 4.3 4.5 5.1 5.1 2.8 4.8 4.6 4.8 4.8 4.7 5.5 4.8 4.0 4.7 4.9
##  [397] 3.6 2.0 3.3 4.2 5.1 4.5 4.1 6.8 4.7 5.6 5.3 2.8 5.5 4.7 6.4 4.3 2.3 4.7
##  [415] 5.7 4.7 4.9 4.7 2.9 5.8 4.6 4.3 5.4 4.8 3.3 5.1 3.8 3.8 5.5 3.8 4.0 4.7
##  [433] 4.2 4.8 5.9 5.2 4.1 6.2 5.6 4.5 5.3 3.2 4.6 4.7 4.5 5.0 3.8 3.7 3.9 5.5
##  [451] 4.4 3.3 4.3 3.5 4.8 4.1 4.6 3.6 4.3 5.3 4.8 4.9 5.0 3.7 5.1 6.0 5.0 5.6
##  [469] 3.3 4.9 3.6 6.0 3.2 4.6 5.4 4.7 3.8 5.5 4.2 3.5 5.2 6.7 5.4 3.3 5.1 4.9
##  [487] 4.7 5.4 4.5 4.2 3.7 4.0 5.1 3.1 5.5 5.0 3.0 3.8 4.7 4.1 3.5 4.3 3.6 5.0
##  [505] 3.7 5.2 5.0 4.2 4.3 3.8 5.8 3.6 6.2 3.3 6.1 3.3 3.2 3.0 4.5 4.9 5.5 4.2
##  [523] 4.1 3.6 3.7 5.4 4.1 5.4 5.3 4.5 4.3 4.3 5.8 5.2 4.4 5.0 4.4 3.0 4.5 5.1
##  [541] 4.4 5.4 3.5 4.7 4.4 3.7 3.7 3.7 4.0 5.1 4.5 6.0 6.0 3.4 5.3 4.2 5.7 3.2
##  [559] 5.0 4.4 5.3 5.2 4.5 3.7 4.4 4.0 4.6 2.5 5.4 3.9 3.6 3.7 4.3 4.0 5.1 5.0
##  [577] 3.2 5.4 5.5 3.1 3.6 4.4 2.9 4.1 4.8 4.6 4.1 4.8 5.1 4.6 5.0 4.9 3.5 3.3
##  [595] 3.5 5.3 5.9 4.9 4.3 4.4 2.5 4.3 4.9 4.7 4.9 4.3 3.5 5.1 4.7 2.8 4.9 4.1
##  [613] 4.3 5.4 5.3 4.4 2.9 3.6 3.5 5.4 5.3 3.3 4.6 3.6 4.0 4.8 5.3 5.7 5.7 6.0
##  [631] 5.7 4.4 4.8 3.7 5.4 4.6 3.2 2.7 3.3 4.6 6.0 5.0 4.7 5.6 3.3 3.8 3.5 5.4
##  [649] 4.1 4.4 4.5 4.9 5.4 3.4 3.9 3.5 5.6 4.7 3.5 5.8 3.6 5.7 5.3 4.7 6.3 3.5
##  [667] 3.2 5.3 3.3 4.3 3.8 5.7 5.4 6.0 5.1 5.4 3.9 4.9 3.6 3.8 4.0 2.7 6.6 3.5
##  [685] 6.1 4.8 3.9 4.8 4.7 4.4 4.9 4.4 4.1 5.8 4.4 4.9 4.0 3.7 5.0 4.0 5.8 5.3
##  [703] 4.9 4.1 6.1 5.6 3.8 4.0 5.0 4.6 3.5 6.1 4.5 3.6 4.2 4.3 3.3 5.7 3.3 4.6
##  [721] 4.0 4.8 5.6 2.7 3.4 5.3 3.6 4.3 5.9 4.9 5.5 3.5 4.8 3.7 4.2 4.0 5.5 5.0
##  [739] 4.2 3.5 4.2 5.7 4.7 5.2 4.1 4.5 3.6 4.3 3.1 4.1 3.5 4.0 4.7 5.6 6.0 3.7
##  [757] 5.4 4.9 3.2 4.0 4.9 4.3 3.9 6.2 4.4 5.5 5.5 4.4 4.3 5.3 4.2 3.8 4.7 3.9
##  [775] 4.0 3.7 4.5 3.3 3.4 3.9 6.1 4.6 3.5 2.2 3.4 4.6 2.5 4.2 5.4 3.7 5.6 5.2
##  [793] 5.7 5.9 5.2 5.2 4.6 4.1 5.3 3.0 4.6 3.7 5.8 4.9 2.9 4.5 5.1 5.4 4.8 4.1
##  [811] 6.4 3.4 5.2 4.5 4.9 3.7 5.3 5.0 5.5 4.3 4.4 2.9 3.1 5.7 5.2 5.4 4.5 1.8
##  [829] 4.9 3.9 3.2 3.1 2.1 6.2 3.5 3.5 4.0 6.2 4.0 4.4 4.5 3.9 4.7 3.9 5.2 3.9
##  [847] 5.0 3.5 3.7 5.1 4.5 6.1 4.8 4.3 5.2 5.3 3.8 4.8 5.4 5.1 4.1 3.7 4.8 5.3
##  [865] 3.3 4.8 5.8 4.3 2.7 4.9 5.3 6.2 2.7 4.5 5.8 3.9 5.6 5.5 5.4 3.5 4.9 4.0
##  [883] 3.9 4.0 5.2 4.8 3.6 3.8 4.5 3.8 6.2 4.7 6.5 4.3 6.0 3.6 4.9 5.7 3.3 6.2
##  [901] 4.7 4.7 4.2 4.3 4.2 3.7 2.9 3.4 6.5 5.3 4.6 3.8 5.6 4.7 4.4 5.0 4.8 3.6
##  [919] 4.4 5.0 5.1 4.9 4.6 4.2 4.6 4.7 4.9 4.4 4.2 3.4 4.3 4.3 3.8 5.2 4.4 5.8
##  [937] 4.1 4.6 4.0 5.0 5.2 6.6 4.6 3.8 4.5 3.4 5.5 5.3 4.3 5.9 4.5 2.5 6.2 5.4
##  [955] 4.2 5.1 4.5 6.5 3.6 5.9 4.9 5.4 4.0 5.2 4.4 5.7 5.2 4.1 3.1 3.9 4.4 5.3
##  [973] 5.5 3.4 4.8 6.9 3.9 3.3 4.2 3.5 5.4 2.9 2.9 4.2 3.7 4.7 4.3 4.3 2.9 3.7
##  [991] 4.6 3.8 3.4 4.0 4.2 3.8 5.1 6.0 5.1 3.9
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
##    [1] 4.2 4.2 5.3 5.6 4.3 4.1 3.0 4.4 4.2 4.2 3.7 4.8 6.4 4.8 4.2 4.7 4.1 4.1
##   [19] 5.2 4.3 3.9 4.1 3.9 3.5 5.9 4.2 3.6 4.9 4.6 4.3 4.5 5.6 3.5 4.6 4.3 5.3
##   [37] 4.2 4.2 4.4 4.9 4.9 4.7 4.4 3.5 3.0 4.5 2.7 6.0 5.7 4.2 5.2 5.7 3.0 3.9
##   [55] 5.6 5.3 5.0 4.8 4.0 4.7 3.2 5.3 5.1 5.3 4.3 3.8 3.2 4.4 3.6 3.7 4.6 5.7
##   [73] 3.8 3.1 2.8 5.1 3.4 4.6 3.8 3.4 3.1 4.3 4.7 4.7 3.7 3.5 5.0 4.4 4.6 5.0
##   [91] 5.5 4.3 5.6 4.5 4.0 4.5 5.5 3.4 4.5 5.6 5.5 5.1 5.1 5.4 3.2 3.4 4.7 3.6
##  [109] 4.4 5.9 5.2 4.5 5.2 3.3 4.3 4.9 6.1 4.8 4.0 3.4 5.2 4.1 3.9 4.0 5.3 6.2
##  [127] 5.1 5.5 4.2 3.7 6.3 5.2 3.8 5.6 4.8 3.8 4.5 5.1 3.8 3.3 4.1 4.7 3.8 3.9
##  [145] 4.1 5.6 4.5 4.3 3.9 3.8 5.5 4.3 3.2 4.6 3.8 4.7 5.1 5.2 4.7 5.8 3.9 3.2
##  [163] 5.0 5.7 4.1 4.2 3.7 4.3 2.8 3.8 5.3 4.1 5.0 4.7 4.1 5.1 4.2 3.7 3.0 5.9
##  [181] 3.6 5.5 5.3 4.5 4.6 4.9 2.6 5.5 6.0 4.9 4.6 4.8 3.9 4.7 4.7 3.0 5.7 3.7
##  [199] 3.8 4.6 4.6 3.8 4.8 3.6 5.1 4.8 4.7 3.8 4.8 4.3 4.5 4.3 4.2 5.5 5.8 4.2
##  [217] 6.0 4.0 4.9 4.3 2.9 4.4 4.8 4.0 5.1 4.3 5.1 4.8 3.9 2.6 5.1 3.3 4.6 5.0
##  [235] 4.2 5.7 5.0 4.9 4.5 3.9 4.9 4.0 4.6 5.1 4.5 5.0 4.8 4.2 4.1 2.6 3.3 5.0
##  [253] 5.2 5.7 5.1 3.8 5.7 5.4 3.3 2.7 5.5 3.4 4.7 5.0 3.4 3.4 3.8 3.1 4.7 3.6
##  [271] 4.5 5.6 3.7 4.8 6.1 4.8 5.2 3.9 3.3 4.1 4.5 5.1 4.4 4.3 5.3 5.7 3.8 4.8
##  [289] 6.3 6.7 5.0 4.4 3.4 2.6 3.6 3.5 3.5 4.6 4.5 4.7 3.8 4.6 5.3 3.6 5.0 5.5
##  [307] 5.3 5.0 4.5 5.8 4.3 5.2 1.7 3.9 5.6 3.6 6.6 4.5 5.4 3.3 3.3 4.1 4.4 5.2
##  [325] 3.7 4.1 5.0 4.7 5.0 5.1 3.7 5.7 4.5 5.8 4.0 7.5 4.1 4.2 5.3 7.1 4.1 3.9
##  [343] 4.7 5.4 4.8 3.6 4.1 3.9 5.6 4.4 4.2 5.2 6.1 3.7 4.8 4.9 4.4 4.4 5.2 4.8
##  [361] 4.9 3.8 4.1 5.3 3.9 4.1 3.4 4.1 5.6 3.5 4.1 2.5 5.3 3.6 3.8 6.0 4.0 3.9
##  [379] 5.0 4.6 2.8 4.5 4.4 5.4 4.5 4.5 3.9 6.6 4.2 4.4 5.1 5.1 4.7 5.2 5.2 5.7
##  [397] 3.1 5.2 4.3 3.8 4.5 4.2 3.6 4.8 5.4 4.3 5.2 4.2 4.5 4.9 3.0 4.4 5.1 5.2
##  [415] 5.2 4.0 5.6 5.2 2.6 5.1 3.5 5.7 3.7 4.5 5.2 4.2 6.4 5.1 4.3 4.6 3.5 5.3
##  [433] 3.7 5.2 5.7 2.2 4.0 3.8 5.8 4.8 3.8 4.2 4.8 4.5 3.0 4.8 3.0 4.9 3.7 3.6
##  [451] 5.1 4.3 3.9 6.5 5.0 4.1 4.6 5.2 3.4 3.2 4.3 5.1 4.3 3.6 5.7 6.2 3.7 2.9
##  [469] 3.4 6.2 5.2 3.6 4.0 5.6 5.2 4.0 4.2 4.2 5.2 2.1 4.0 5.1 4.3 5.7 4.7 4.3
##  [487] 5.0 3.3 3.0 4.7 5.5 5.6 5.4 3.8 3.4 5.5 3.5 6.3 5.8 6.7 4.6 4.5 2.7 4.0
##  [505] 2.3 4.8 4.9 4.3 4.8 4.7 5.4 3.9 4.2 6.0 3.5 4.9 4.9 4.4 4.2 5.2 4.8 4.8
##  [523] 4.7 6.6 5.0 5.5 3.9 4.0 4.2 5.0 4.7 4.7 4.3 4.2 4.3 5.5 5.2 4.2 4.1 4.4
##  [541] 5.5 5.5 3.6 3.8 6.2 4.8 4.1 5.7 5.8 4.2 4.7 4.7 4.2 5.1 3.1 3.1 4.7 4.0
##  [559] 4.8 6.6 4.9 4.0 5.3 4.7 5.2 1.7 5.5 5.7 2.9 6.3 4.2 6.2 3.1 3.8 4.4 3.5
##  [577] 5.8 4.1 4.3 5.2 4.1 4.0 2.3 5.1 3.1 4.6 5.2 5.5 3.5 4.8 5.8 4.5 3.1 4.7
##  [595] 4.7 3.3 3.8 5.3 5.7 5.9 3.4 4.7 3.0 3.0 5.5 4.9 5.7 4.2 4.1 5.4 4.6 4.9
##  [613] 4.6 4.0 5.0 3.6 6.4 4.1 5.7 5.0 3.9 3.7 4.6 3.0 5.3 3.9 6.0 4.4 4.8 4.0
##  [631] 4.5 5.2 5.1 4.6 4.4 4.6 3.7 4.3 4.8 5.2 4.4 5.4 5.0 5.6 5.4 6.1 3.9 3.5
##  [649] 5.4 4.4 4.9 5.4 5.1 4.8 4.4 6.7 5.1 5.2 5.2 4.6 5.3 3.6 5.2 5.6 3.2 5.2
##  [667] 4.9 3.6 4.5 3.9 2.5 7.1 5.9 4.9 4.5 5.6 4.1 5.5 5.8 4.6 4.6 4.9 4.7 4.4
##  [685] 4.5 5.1 3.9 4.4 4.0 4.4 4.2 5.5 5.3 4.7 4.2 4.5 4.4 5.4 5.4 4.7 2.8 4.8
##  [703] 3.2 4.5 4.9 4.4 4.8 3.9 3.7 4.2 5.1 5.6 5.0 5.7 4.9 2.7 5.1 6.1 3.9 4.8
##  [721] 4.7 4.3 2.7 2.6 2.9 5.6 4.7 2.6 4.1 6.6 4.4 4.3 2.8 3.2 3.9 3.8 6.0 4.4
##  [739] 5.2 4.7 5.1 6.0 3.9 5.4 4.5 4.5 4.4 4.7 6.0 4.5 4.2 4.6 4.3 4.9 2.7 5.9
##  [757] 4.8 4.1 5.1 6.4 5.0 5.6 4.8 5.9 5.9 5.1 5.1 4.5 5.0 3.0 5.1 6.7 5.2 4.4
##  [775] 5.3 4.2 5.3 4.5 4.0 3.6 3.4 4.3 4.9 3.7 4.1 6.1 3.5 2.8 4.1 4.5 3.9 5.1
##  [793] 5.4 4.8 4.5 3.8 4.5 3.6 4.0 5.4 6.0 4.4 3.5 3.8 5.3 4.3 4.8 5.8 3.9 3.5
##  [811] 4.9 5.2 4.5 4.1 4.5 5.8 4.0 3.6 6.4 5.1 4.0 4.6 3.3 3.6 6.2 5.9 5.3 5.6
##  [829] 4.7 4.5 3.9 5.3 4.0 4.4 2.8 5.7 4.3 4.2 4.8 4.0 4.1 4.3 5.3 5.6 4.9 5.0
##  [847] 3.7 5.4 5.3 3.9 5.4 4.7 4.4 5.3 5.2 3.8 2.5 5.1 5.1 5.6 4.2 3.8 4.0 6.1
##  [865] 3.9 3.8 4.5 4.6 5.6 3.7 5.1 4.4 4.0 4.4 5.1 5.6 5.2 5.0 5.2 5.3 5.4 3.8
##  [883] 5.2 5.8 3.8 4.1 5.3 6.1 3.8 4.4 4.7 4.9 4.9 4.1 3.2 2.0 4.9 3.8 4.5 5.4
##  [901] 4.0 4.0 4.9 3.5 4.6 4.0 3.1 4.9 4.3 6.3 3.5 6.0 4.4 5.0 2.4 4.2 3.5 5.2
##  [919] 5.0 4.2 4.1 4.2 4.7 4.5 4.9 3.0 4.3 4.3 4.4 4.9 4.4 4.5 4.6 4.3 4.6 4.9
##  [937] 5.2 4.5 3.4 4.0 5.8 3.6 5.4 5.0 4.0 5.0 4.7 4.8 4.7 4.1 4.3 3.8 2.9 4.1
##  [955] 4.9 6.2 5.5 5.3 3.5 3.2 4.6 5.5 4.4 5.3 4.3 4.3 4.3 6.2 4.3 4.2 3.0 4.4
##  [973] 5.4 3.3 3.4 5.6 3.8 4.8 4.9 3.3 4.3 5.0 3.6 5.4 4.5 3.8 5.5 5.9 2.5 4.7
##  [991] 5.3 4.3 5.1 2.8 4.6 3.8 5.6 4.1 3.4 3.6
## 
## $func.thetastar
## [1] 0.0379
## 
## $jack.boot.val
##  [1]  0.49892473  0.39942029  0.28206522  0.17159420  0.04752747  0.01927711
##  [7] -0.09298780 -0.26810345 -0.27863501 -0.46338462
## 
## $jack.boot.se
## [1] 0.8914699
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
##    [1] 4.6 4.1 5.2 3.9 5.3 2.9 3.9 5.2 4.6 5.2 5.4 4.7 3.5 4.2 5.1 5.6 4.9 3.3
##   [19] 4.5 4.5 4.6 4.6 5.6 4.5 5.6 3.2 4.2 3.5 5.7 4.1 3.3 5.2 5.6 3.5 4.5 4.4
##   [37] 5.1 5.1 5.0 2.0 3.4 4.5 4.8 4.5 5.1 5.4 4.5 4.2 2.2 4.1 5.0 4.8 5.7 4.4
##   [55] 6.2 3.7 4.6 4.5 4.3 4.2 3.7 5.1 5.4 4.1 2.5 5.1 4.8 4.8 4.0 5.6 4.2 3.2
##   [73] 5.1 4.0 6.0 3.7 5.2 4.8 4.4 4.7 3.3 3.6 4.7 6.3 3.7 5.4 5.6 4.2 4.4 4.4
##   [91] 3.9 2.1 5.1 3.9 3.4 4.7 4.7 3.6 4.5 4.0 2.4 3.2 5.6 4.2 5.8 5.5 4.6 5.4
##  [109] 5.2 4.2 4.0 4.6 4.1 4.2 2.9 4.7 5.5 4.9 4.9 5.3 4.0 5.5 2.8 3.4 3.9 3.8
##  [127] 5.1 3.8 4.2 4.0 4.8 4.5 4.9 6.4 4.4 5.2 4.0 4.4 5.3 4.5 2.9 5.8 5.3 4.7
##  [145] 4.2 4.0 4.5 2.1 2.8 6.0 5.1 4.4 4.1 3.6 5.0 4.1 3.8 4.1 4.6 5.2 2.7 4.4
##  [163] 5.3 3.3 5.7 5.4 4.0 4.4 3.1 3.8 3.8 5.5 5.4 4.6 5.6 3.7 4.1 4.2 6.0 5.3
##  [181] 5.3 3.0 2.6 1.1 4.3 4.4 3.8 4.4 4.8 5.5 5.0 5.3 4.8 4.1 6.5 5.8 4.7 4.7
##  [199] 6.0 5.3 4.5 3.7 4.9 4.7 3.0 4.7 3.6 6.3 4.3 5.9 3.3 5.4 3.7 6.1 4.5 4.1
##  [217] 4.6 6.5 4.8 5.4 5.2 4.7 3.2 4.1 5.8 3.9 4.6 5.2 3.2 4.4 4.7 4.7 4.3 5.6
##  [235] 3.3 4.0 5.0 5.5 3.9 4.2 3.0 5.8 3.2 3.8 4.1 4.3 4.0 4.1 5.5 3.8 3.8 4.2
##  [253] 4.2 3.3 4.5 4.6 4.3 2.4 5.6 4.9 5.9 2.9 5.2 5.0 3.4 5.1 3.2 4.6 6.1 4.5
##  [271] 3.0 5.0 3.9 3.2 3.2 4.6 5.9 2.7 3.4 4.9 5.4 5.3 4.8 5.7 3.3 4.8 4.4 4.4
##  [289] 4.6 4.3 5.8 5.6 4.0 4.4 5.4 5.7 5.5 4.7 2.4 5.9 3.4 5.7 5.4 5.9 5.1 4.7
##  [307] 4.4 3.6 4.4 3.6 3.9 7.0 4.5 3.4 5.2 4.0 4.6 3.5 6.0 5.3 5.1 4.3 3.9 4.4
##  [325] 4.4 4.2 4.1 4.7 6.0 3.8 3.6 4.0 6.9 4.7 5.1 4.5 4.7 5.4 3.6 4.8 4.2 5.0
##  [343] 5.6 5.6 5.0 4.1 5.8 4.9 3.6 3.1 5.0 3.7 6.9 4.5 4.7 2.9 4.8 2.8 4.3 3.5
##  [361] 3.7 5.2 3.6 4.5 4.5 2.7 3.2 3.1 3.8 3.2 6.6 4.0 5.4 3.3 4.3 3.0 4.5 4.0
##  [379] 6.4 3.0 4.2 4.2 4.0 4.6 5.4 2.3 5.0 6.3 4.7 3.7 3.8 4.3 5.2 4.4 3.4 6.1
##  [397] 2.7 3.8 3.9 4.5 3.9 5.3 6.1 5.8 3.6 4.2 4.3 5.0 4.0 4.4 6.1 4.7 4.9 6.2
##  [415] 4.2 5.1 5.4 6.1 4.9 3.1 4.1 5.0 4.2 4.7 4.3 5.0 3.3 3.9 3.7 4.3 5.0 6.2
##  [433] 4.2 5.5 3.1 4.5 3.8 4.0 6.4 5.8 5.2 4.0 3.1 4.5 3.4 4.9 5.4 5.6 4.9 4.5
##  [451] 4.5 4.9 5.1 4.5 3.9 6.4 3.3 4.7 5.1 3.2 5.4 5.0 5.6 3.7 4.9 6.1 5.4 5.0
##  [469] 5.4 4.9 6.3 5.9 4.1 4.8 4.9 4.6 5.1 4.8 5.3 5.2 4.3 4.7 3.9 5.7 4.3 4.2
##  [487] 5.5 4.4 3.8 5.1 4.1 4.1 3.3 4.2 6.6 3.8 5.2 4.2 5.9 6.1 4.8 3.9 3.7 5.3
##  [505] 4.5 3.9 3.6 5.4 4.2 6.3 3.9 3.9 6.1 5.6 4.9 4.1 5.1 4.7 6.3 3.1 4.2 3.3
##  [523] 4.4 6.2 3.3 4.6 4.8 3.5 4.6 4.4 4.5 4.8 3.9 5.3 3.9 4.4 4.7 4.5 3.8 4.3
##  [541] 5.1 3.8 3.8 3.2 4.0 4.4 5.5 5.0 4.1 4.0 5.3 4.8 3.2 4.5 5.4 4.5 2.4 3.4
##  [559] 4.3 5.4 4.8 5.4 4.7 4.5 4.0 5.3 3.9 4.4 4.9 4.8 5.1 5.6 4.4 5.0 3.1 6.1
##  [577] 5.5 4.3 5.5 4.9 4.2 5.0 5.3 3.9 4.5 4.2 4.6 3.7 5.2 4.3 5.2 4.0 2.9 6.6
##  [595] 3.9 4.4 4.6 4.0 4.7 3.7 3.9 5.3 5.9 4.2 4.5 3.9 3.3 3.7 5.4 5.7 4.0 4.6
##  [613] 5.1 5.0 3.9 3.8 5.1 4.2 5.3 5.2 5.0 5.7 3.9 5.0 5.7 5.1 3.5 5.7 5.5 5.1
##  [631] 6.0 4.9 4.8 4.0 3.7 5.4 3.0 4.3 6.8 3.7 4.0 4.1 5.3 5.0 5.4 3.8 4.6 6.5
##  [649] 3.1 3.3 5.1 3.7 4.7 4.0 5.6 4.0 4.6 6.2 4.4 4.2 4.8 4.2 4.2 5.5 4.6 4.2
##  [667] 5.2 3.6 4.3 3.9 2.6 5.1 5.4 4.5 5.9 3.8 3.1 4.8 4.5 4.5 5.2 3.7 4.7 6.1
##  [685] 5.2 4.0 3.5 4.1 3.7 4.7 5.1 5.9 3.7 5.0 4.5 5.2 5.3 4.0 5.1 3.5 5.0 5.8
##  [703] 6.9 2.9 5.3 5.3 5.5 6.1 2.7 3.6 4.7 5.3 5.6 4.2 4.7 4.5 4.4 5.1 3.6 4.6
##  [721] 4.8 5.8 5.7 4.8 3.3 4.8 4.6 4.7 4.4 3.6 4.5 3.0 4.2 6.0 3.3 4.3 5.6 3.2
##  [739] 4.1 4.1 4.9 3.7 5.1 4.2 6.0 3.2 3.3 5.1 4.9 7.3 5.0 3.1 5.3 5.8 4.8 5.3
##  [757] 5.5 4.6 5.8 4.4 5.5 5.9 3.9 2.5 4.0 5.1 5.4 5.5 7.5 2.6 4.2 4.0 5.1 3.4
##  [775] 3.9 3.9 3.9 4.0 5.2 4.4 3.8 5.1 5.3 4.0 4.9 2.6 4.1 3.4 4.2 4.8 5.0 3.5
##  [793] 5.7 4.0 4.9 4.8 4.7 3.9 3.8 4.9 3.1 3.5 4.6 6.8 4.8 3.6 5.4 3.9 6.9 3.9
##  [811] 5.1 3.8 4.5 5.7 4.7 5.6 3.0 2.9 2.3 4.5 3.9 5.3 5.3 2.8 2.9 2.9 3.5 4.7
##  [829] 3.7 3.6 5.4 5.0 5.1 4.8 5.0 5.7 5.7 3.0 5.2 5.2 5.9 3.7 3.7 4.2 3.6 5.2
##  [847] 3.6 3.9 4.4 3.2 5.0 3.8 2.5 4.0 3.4 4.1 3.3 4.6 4.6 5.0 5.4 6.1 5.6 5.0
##  [865] 3.7 4.9 7.5 4.5 3.7 4.3 5.4 5.8 3.5 4.8 3.2 4.5 3.6 5.9 4.4 5.3 5.9 3.1
##  [883] 5.0 2.5 4.9 4.6 4.0 4.4 4.7 4.8 4.2 4.7 5.7 4.2 4.0 4.8 5.0 5.0 4.9 5.0
##  [901] 3.8 4.7 5.4 4.4 3.0 5.6 6.6 5.4 5.9 4.9 3.9 5.1 4.2 3.0 4.5 5.6 4.8 3.9
##  [919] 3.4 5.4 5.8 5.6 4.0 4.2 5.9 4.4 4.2 5.1 2.9 4.8 4.8 3.8 4.0 5.9 6.3 4.0
##  [937] 5.4 2.5 4.9 5.0 3.8 6.6 4.5 3.1 4.2 4.5 5.0 5.1 4.3 4.7 4.0 4.0 3.9 5.5
##  [955] 3.0 4.9 3.6 4.1 6.8 4.6 3.5 4.6 4.8 4.0 3.4 4.8 5.1 4.0 4.9 3.8 5.0 5.3
##  [973] 4.0 4.7 5.2 4.8 4.0 4.5 4.5 2.9 4.6 3.6 5.1 3.8 4.3 3.5 4.9 4.2 6.4 4.9
##  [991] 4.9 4.2 4.4 5.5 4.5 3.5 4.2 4.0 4.2 4.1
## 
## $func.thetastar
## 72% 
## 5.1 
## 
## $jack.boot.val
##  [1] 5.5 5.4 5.4 5.2 5.2 5.0 4.9 4.8 4.6 4.5
## 
## $jack.boot.se
## [1] 0.9881801
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
## [1] 0.1065881
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
##    6.639282   17.208156 
##  ( 2.897644) ( 7.801603)
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
## [1] 0.14665634 0.26187412 0.59498699 0.29052991 0.01739739 0.52149440
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
##    [1]  1.1982178026 -0.4045941693 -0.2290926311 -0.7921052717  0.6686758697
##    [6]  0.3389657454  0.3074699574  0.0293934466 -0.1389224745 -0.1067255177
##   [11]  1.0962731560 -0.3668544386 -0.1908662448  0.0971565749  0.0123838477
##   [16] -0.5817026779  0.9222879894  0.4029561549 -0.1578659066 -0.9632529993
##   [21] -0.2890495419  0.3018828417 -0.0897408551  0.1480648135  0.7236980917
##   [26] -0.1172818212  0.3580270979  0.1065880550 -0.0733627219  0.2029771070
##   [31]  0.1821779431  0.4033845148 -0.6869801304 -0.3184853183  0.5702578777
##   [36]  0.7823729860  0.4602603769 -0.1786768967  0.7790033483 -0.2582830521
##   [41]  0.3292258296  0.3857757859 -0.0558483441 -0.0604086964  0.0256435543
##   [46]  0.6066811681  0.1140820148 -0.3022663072  0.1513144579  0.4911726343
##   [51] -0.2476641835  0.2076112748  1.3111762742  0.5087058940  0.1670733911
##   [56]  0.1065880550 -0.2012474809 -0.2337394714  1.3899335753  0.0066053035
##   [61] -0.5039535125  0.3532374648 -0.6421591366  0.1919652373 -0.0767956289
##   [66]  0.3739324493  0.3879444559 -0.3327346603  0.4161284045  0.7494254767
##   [71] -0.0218237139  0.0129419502  0.8463724221  0.5565914431 -0.3223853054
##   [76]  0.6701576291  0.1104253428  0.0288773771  0.1702058727  0.4351220194
##   [81]  0.6327665917 -0.2830133116 -0.1670184387  0.2609421486  0.0865080624
##   [86]  0.2819333460 -0.0468710648  0.2756368112 -0.5912349198  0.2561750651
##   [91]  0.1950375257  0.4030189334  0.1499532215  0.1591598423  0.1712432918
##   [96]  0.3433949249 -0.0635089738  0.2061358817  0.0957675034  0.6760954375
##  [101]  0.0532134693  1.6894979693  0.1081174772 -0.1362624382  0.1198219181
##  [106] -0.2712120862  0.4008705953 -0.3730906688 -0.4824400430  1.0948245885
##  [111]  1.1439450361  0.5668328142  0.2574118562  0.6574133890 -1.4517179303
##  [116]  0.0082531382 -0.6914471062  0.0517708921  0.0405424245  0.0783390162
##  [121] -0.3377442650  0.7585536052  0.2225567191  0.2347914813  0.3450474702
##  [126]  1.3269717285 -0.2243469934 -0.7031279777  0.5022488994 -1.0502252519
##  [131]  0.4477347803  0.4089635843  0.7228961141 -0.0658959404  0.0395391719
##  [136]  0.4986000044  0.3496173944  0.5910138726 -0.6836855902 -0.9377238151
##  [141]  0.6351913203 -1.4520810932  0.6929748847 -0.1899428463  0.1475459626
##  [146]  0.2000800338 -0.2907226929  0.3995330959  0.1919929531 -0.1791337570
##  [151]  0.0174191999  0.8077410496 -0.1043879015 -0.0462654214  0.4670321120
##  [156]  0.1081129073  0.1744627581  0.0751194262  0.5030064962  0.2311582693
##  [161]  0.3316005264 -0.2730497826  0.0117999959 -0.7219939802 -0.0552893608
##  [166]  0.3582861855 -0.9304662122 -0.1234950650  0.2514787248  0.7480686075
##  [171] -0.5595884788  0.3684168176  0.2398673371  0.5553658388  0.4082449714
##  [176]  1.2254450370 -0.4926163444  0.0241354662  0.9026327967  0.5073072149
##  [181]  0.1414451915  0.2362549738 -0.2170281085  0.3487228520  0.1699271558
##  [186] -0.2494049300 -0.1541575219  0.5883486707  0.6070616488 -0.5788218762
##  [191]  0.1887532879  0.3656335368  0.0590681173 -0.2029600232  0.6407452972
##  [196]  0.0626870418 -0.2594871531  0.2000490545  0.2179859311  0.2009622394
##  [201]  0.8816612038  0.0716166022  0.3177999853 -0.7490352865  0.0617229977
##  [206]  0.2312787646  0.2619123854 -0.2876967165 -0.7752693513 -0.1903922689
##  [211]  0.9045438511  0.2509267985 -0.2274732892 -0.7155012734  0.3156051659
##  [216]  0.3203052418  0.2896732364  0.3311946164 -0.1786768967  0.3976411015
##  [221] -0.1493407506  0.8013399852 -0.4009228855 -0.9372685345 -0.1451065018
##  [226]  1.0573291047 -0.4830371821  2.0679214743 -0.2398229859  1.6143934963
##  [231]  0.1126671738 -0.4168933454 -0.0604130222  0.3182305052  0.6584419171
##  [236] -0.6382578181  0.4895104686 -0.0572954961 -0.2008469821  0.3947075219
##  [241]  0.0445351997  0.9152736659 -0.1587673291 -0.4494842495 -0.4025600906
##  [246]  0.5484520774  0.2353645543 -0.4282096666  0.1471328573  0.3672231907
##  [251]  0.3037417278 -0.1726472017  0.2897347783  0.1743404348 -0.0332432610
##  [256]  1.1144745709  0.3720565493  0.3830561993 -0.1264599527 -0.2542125861
##  [261] -0.4478015099 -1.0182112305 -0.0267756963  0.5006195817  0.0565284505
##  [266]  0.1827718611  0.4277672162 -0.2583337275 -0.7084480242  0.4970278696
##  [271]  0.2295735409  0.3839492629  1.1843518658  1.3415393396 -0.0540972468
##  [276] -0.1677783268 -0.6136098526  0.0305020976 -0.2468612291  0.6239842964
##  [281] -0.4600300220  0.5689423147  0.6814566943 -0.1605599080  0.3414102403
##  [286]  0.2599814735 -0.0178701512  0.6079635452 -0.1506464624  0.2392084714
##  [291] -0.7696983452 -0.3176808623 -0.0446193550  0.4634528172 -0.0060358344
##  [296] -0.0994058675  0.0413577793  0.1129656461 -0.7965258825  0.9582340356
##  [301] -0.1660092440  0.3200177773 -0.3365896846 -0.2174153118  0.0844796698
##  [306] -0.6612253404  0.1112705161 -0.0638229083  0.1312505540 -0.4126481644
##  [311] -0.2090224760 -0.7664725887  1.3312470217  0.7203932784 -0.5718432293
##  [316]  0.4360649313 -0.3267119726 -0.9829351210 -0.1963149680  0.3164916815
##  [321] -0.3404896133 -0.0890581014 -0.2481510218 -0.0507153210  0.6543762794
##  [326]  0.0194418419  0.6428195718 -0.1589007716  0.1995642655 -0.5297809316
##  [331] -0.1400069992 -0.2845960366  0.4194939946  0.5480810703  0.1101789571
##  [336] -0.4478015099  0.4354335112 -0.9485597527  0.0504148319 -0.1232798544
##  [341] -0.2289840362  0.3884104998  0.1063036512 -0.0765013695  0.2611262554
##  [346] -0.2588719502  0.6893048935  0.1543551224 -0.2634490831  0.5475821977
##  [351]  0.9555916181  0.2492573953 -0.1819914389  0.1748214893 -0.0898738948
##  [356]  0.8479915254  0.3039377749 -0.4153818818  0.9936814646 -0.2855201082
##  [361]  0.2464764519  0.3529871786 -1.1505107250 -0.5221689034  0.3307685485
##  [366]  0.1030032738 -0.4711025064  0.6546041463  0.1463300159  0.3452484650
##  [371]  0.0342031348  0.7306435530 -0.1027816089 -0.0204502905  0.0825048211
##  [376] -0.4952515749  0.1800406885 -0.6279925365 -0.3342169582  0.4339900162
##  [381] -0.0939775002 -0.1039843079 -0.0668421479 -0.8800985030 -0.7447558760
##  [386]  0.3975995831 -0.5090695500  0.5183574793  0.3062522477  0.4300160423
##  [391] -0.3446893947  0.0789260748  0.9935099614 -0.0453163113  0.2281012994
##  [396]  0.4139412166  0.7403062819  0.7413421219 -0.3841349008 -0.1684944932
##  [401] -0.2914422661  1.2141763400 -0.3923059620  0.4508115677 -0.4217204197
##  [406]  0.0673548299  0.3360661099  0.1073929116  0.3489626817 -0.5908260022
##  [411] -0.2911625344 -0.5832940571  1.0958614379  0.1411609162 -0.2701770796
##  [416] -0.6619956198  0.6486781989  0.3079075180  0.0064971525  0.1738593013
##  [421] -0.3990017524 -0.4108502696 -0.1898097173  0.3879086182 -0.3573875899
##  [426]  0.6852908830  0.0423550051  0.1754014147 -0.1078978616  0.4347024867
##  [431]  0.2634596472  0.1375077535  0.7708794082  0.1509288980 -0.1195758506
##  [436]  0.4618167346  0.1535084186  1.2229486496  0.1710026419  0.2917778449
##  [441] -0.5874190023 -0.0996857536  0.0053555980  0.2897397423  0.2677083762
##  [446]  0.3088094899  0.4056106448  0.4341538442 -0.1539537386  0.6382244906
##  [451] -0.4791540089 -0.0102503316  0.2283533343 -0.1121213078  0.1368708810
##  [456] -0.4643442566  0.4435224145  0.7927664248  0.0826253097  0.1545395051
##  [461]  0.0671270776 -0.6667085572 -0.1664218257  0.4051608018 -0.1328515805
##  [466]  1.4338739486  0.1059050533  0.1400543300  0.4959630459 -0.0078914537
##  [471] -0.1866608928  0.0193553326 -0.0708345674  0.4156041174 -0.8134889887
##  [476]  0.3185825453  0.0146919414  0.1531491779  0.5155719474 -0.2198656652
##  [481]  0.8312951119 -0.4531779762  0.6829664341 -0.5792384708  0.3559848713
##  [486]  1.0138586994  0.2059535688 -0.1726434772 -0.3671848043 -0.9798684485
##  [491]  0.3388307821  0.5036251983 -0.7596795949  0.8321681079 -0.1812119806
##  [496]  0.5944011416 -0.0585736205  0.3513839193  0.4593888749  0.4349051934
##  [501] -0.1018152533 -0.1954130031  0.7432814777  0.1588684783  0.3575655581
##  [506]  0.1287785186 -0.7683129930 -0.0131516045  0.1874801560  0.1915359144
##  [511] -0.0138442311  0.5460373833 -0.2954972336  0.6785467278  0.2276548446
##  [516]  0.2472383355 -0.6821064614  0.2113071843  1.4838974188 -0.9549170736
##  [521]  0.7010942343 -0.1124154756  0.0524189090  0.6869641773 -0.8358524216
##  [526]  0.5961295120 -0.4965311663  0.2669981688  0.1089636721  0.4819947314
##  [531]  0.2503528358 -0.4899107212  0.0236685719  0.7330810283 -0.0959099390
##  [536] -0.3249602290 -0.9609291657 -0.8364772947  0.0998590748 -0.2060876846
##  [541] -0.5005322129  0.2012915819  0.1287086877  0.0445132601 -0.4736371995
##  [546]  0.0710056310  0.5874696896 -0.1262626625  0.1829124590  0.2235419660
##  [551]  1.7108570689  0.0977696402  0.7589859631 -0.0006550867  0.2802631825
##  [556] -0.0550878909  1.2274383286 -0.1827426651  0.6600931640  0.3594671942
##  [561] -0.2018609666 -0.3951554398 -0.4349708495 -0.0709180891  0.4384080498
##  [566] -0.1681435088 -0.2968218716  0.1168979790  0.1129656461 -0.2175873039
##  [571] -0.2632436329  0.4414736478  0.4347213077  0.0756953367  0.1822261727
##  [576]  0.3077743353  0.3179474215 -0.1344688032  0.0306788000 -0.2489876919
##  [581]  0.2659004047 -0.3489804324  0.4478303719  0.4997339670  1.1448290564
##  [586]  0.5789238865 -0.0550716422  0.1782764488 -0.2596020189  1.0468834827
##  [591]  0.3154884213 -0.3557994770  0.3881129919 -0.5268783497 -0.1324922486
##  [596] -0.6080675394 -0.1927047808  0.8073304308 -0.1232798544  0.0598935987
##  [601]  0.4881505687  0.2642930388  0.0007578665  0.0661432660  0.3523880012
##  [606]  0.4589576534 -0.1204065406  0.3320272590  0.6996416047 -0.8993660447
##  [611]  0.6016197740 -0.3069119830  0.6034408842  0.7909471657  0.6683330968
##  [616] -0.2005659503 -0.0009495990 -0.2367567895 -0.2599566638  0.2897347783
##  [621]  0.4448584390 -0.7859090020 -0.2853339160  1.0074112831  0.4051036358
##  [626] -0.0824487176  0.5528391173 -0.2544347956 -0.1046085446  0.1085165193
##  [631] -0.1528790596 -0.3569880440  0.5044676227 -0.0950565543 -0.3769622249
##  [636]  0.0601904209  0.0882209930  0.0289679561  0.5165726722 -0.2207037792
##  [641]  0.4208226813  0.2853349822  0.3592250486 -0.0229653727  0.5097533399
##  [646]  1.2691696040  0.2603607531  0.5959882210  0.3875997512 -0.0612696048
##  [651]  0.6388703229 -1.5371721038  0.2017408433 -0.2328350846  0.6087656391
##  [656]  0.7421817875 -0.3361244218  0.0830042960 -0.0740640315  0.2147720325
##  [661]  0.5996955959  0.6693439122  0.0455374075 -0.1730090661  0.0941129519
##  [666] -0.1619472425 -0.0682156148  0.2280351697  0.6683330968 -0.0293714396
##  [671]  0.2943093459 -0.6476459695 -0.4506655238 -0.2854976021 -0.9203111597
##  [676] -0.5019734552  0.0516954001 -0.0497488901  1.3598263046 -0.4379995781
##  [681]  0.2641480039  0.2688169659  0.0416676904 -0.8254114688  0.5294638732
##  [686] -0.0202725634 -0.0042697057  0.7103920268  0.2309212106  0.1215558250
##  [691]  1.0716070022 -0.2578092213  0.0306958561  1.0847249259  0.4625364003
##  [696]  0.9187784471  0.3945052000  0.3215149904 -0.0572117249 -0.0183604103
##  [701]  0.1778311517  0.0262732987 -0.2897946245  0.9499508127 -0.5228003908
##  [706]  0.1277773626 -0.1728929743  0.3859973402  0.0323549712  0.0017596373
##  [711]  0.4057306284 -0.2075349209  0.2693948702  0.4055498000  0.2085855654
##  [716] -0.5582976634  1.4402914463 -0.5210066551 -0.1459777609  0.1519706768
##  [721]  0.4965265712  0.2212853562  0.1601939728  0.0760213461 -0.6201361371
##  [726]  0.6258804293  0.5625858411  0.1532597086  0.2332285450 -0.3440577377
##  [731] -0.0179203189  0.8403901953  0.2593469968  0.1513144579  0.2584786565
##  [736] -0.0009201759 -0.4493599074 -0.4739877295 -0.7253831156  0.5333182091
##  [741] -0.8436882774 -0.1899151107  0.1076578915 -0.0903014588 -0.3540766918
##  [746]  0.1305001988 -0.1568966685  0.2863030475 -0.1675484951  0.0598935987
##  [751]  0.7891861629  0.3103312768  0.1553815256 -0.1188195481 -1.1937866541
##  [756] -0.0362930223  0.0004788235  0.0266878405 -0.4162652569 -0.1698810476
##  [761] -0.1899154157 -0.4527451724  1.0847249259 -0.6268193427  0.4716896587
##  [766] -0.1487750337  0.2046483000  0.0975652752 -0.3597347190 -0.2436819957
##  [771]  0.1751346094 -0.6634796269 -0.0218237139  0.8046577081  0.3769561960
##  [776] -0.9763763463 -0.3140562523 -0.3767288047  0.0383841515  0.6461779848
##  [781] -0.0982186753  0.4939662499  1.0873330811 -1.1130886578  0.6198498526
##  [786] -0.6331032723  0.3452993181  0.0137280162  0.9788209488 -0.8034109915
##  [791]  0.2912051627 -0.8239969087 -0.2885531458  0.6864689382 -0.9130809180
##  [796] -0.7752693513  0.1670732637  0.2978303295  0.7612724829  0.4892411946
##  [801]  0.0655190238 -0.9311377054 -0.1401047005  0.7256764406  0.9951742390
##  [806] -0.0740730559  0.2104114439  0.2034844027  0.4291791636 -0.1199644039
##  [811] -0.1606493239  0.0133662397 -0.2010187210 -0.1212945312 -0.0244538926
##  [816]  0.2680323523  0.3381649923 -0.2999535157  0.1751035454  0.4478681815
##  [821] -0.1688798533  0.0007555803  1.1512088380 -0.2392101587  0.9318449827
##  [826] -0.4410729068  0.1584386289  0.4964344569 -0.0367395350 -0.1039292920
##  [831] -1.1603337858  0.3910343530 -0.2706434742 -0.1283341483  0.3693833675
##  [836]  0.2670941153  1.2213935424 -0.3319877979  0.4881505687 -0.1567855207
##  [841]  0.2174767533  0.8819099452  0.8312459165 -0.4498274958  0.6079050827
##  [846] -0.1010270001 -0.6102044673  0.9188500076  0.5212925137  0.5789853845
##  [851]  0.5136640100 -0.1350263404 -0.3941868862  0.3117102788  0.3075746426
##  [856]  0.8804366113  0.3373816823 -0.5390970592  0.0475851564  0.6564739179
##  [861]  0.8399418578  0.2104225421  0.0666019003 -0.3026788094 -0.4957714734
##  [866] -0.0204502905  0.0537527317  0.0014058508  0.5049664082 -0.7808505664
##  [871] -0.3421374608 -0.3732176927  0.1321849416  0.4801700822 -0.5489254558
##  [876]  0.9021391320  0.1795477531  0.5162174714 -0.3895201400  0.8052910231
##  [881]  0.1060385909  0.9844212367  0.2182169287  0.0819253576  0.2531489740
##  [886]  0.9094524087  0.5849382152 -0.1602326823 -0.1476158285  0.1292759909
##  [891]  0.0103441369  0.8737400508  0.2441346890 -0.6761382263 -0.9297627613
##  [896]  0.0286265133 -0.8537858942  0.1589620282 -0.8254114688 -0.0255385542
##  [901]  0.7134747859  0.7134289364 -1.2009283743  0.4258571210  1.3859600635
##  [906]  0.2372940786  0.7134289364  0.4965265712  0.2758055819 -0.4516850547
##  [911] -0.2218273199  0.3548958330  0.1485593314 -0.6265886380  1.1139113448
##  [916] -0.3345961267  0.0531920344  0.5474437849  0.2802707458  0.4331901610
##  [921]  0.4534083924 -0.9441901618 -0.3251510890 -0.0571755204  0.7865972394
##  [926] -1.1921001412 -0.4500019829  0.4123130749  0.2322822698 -0.6235082074
##  [931]  0.2672214312  0.5728916737  0.4122556746  0.3379641337  0.1629677197
##  [936] -0.0769147744  0.2243033499  0.3541160977 -0.1239096691 -1.0737219978
##  [941] -0.4251723077  0.0929807233  0.6070616488  0.4842261997  1.2708349176
##  [946] -0.7373249721  0.1150518115 -0.1352612684  0.1673864351  0.4412457262
##  [951]  0.4594343956  0.3294026327  0.0950376471 -0.2822119157 -0.5442432019
##  [956] -0.2079518511  0.0051534484 -0.3660151140  0.3331215214  0.3383367801
##  [961] -0.1371423142  0.7446346221  1.0600890933 -0.0002831273  0.5065937258
##  [966]  0.5067960977  0.5837301382  1.3710353478  0.5336542292  0.7387129417
##  [971]  1.2158304149 -0.3285098157  0.5402035905  0.3534627189  0.1947315546
##  [976]  0.5995441752  0.2151879684  0.4659157044  0.2040318898  0.7799243364
##  [981]  0.1867945149  0.4521206549 -0.0440543321  0.2836071154  0.1065880550
##  [986] -0.8851940239 -0.1225284335  0.5086535327  1.7134073484  0.3714047327
##  [991]  0.7486696394  0.4349625305  0.8118418313 -0.0028095241 -0.3058630891
##  [996]  0.4471617095 -0.0872333536  0.0950376471 -0.1743298489 -0.3832752221
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
##   0.38581063   0.14529033 
##  (0.04594483) (0.03248207)
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
## [1] -0.95435442 -0.39153472 -0.23016692  0.07357029  0.79932298  0.78465806
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
## [1] 0.0221
```

```r
sd(results2$thetastar)  #the standard deviation of the theta stars is the SE of the statistic (in this case, the mean)
```

```
## [1] 0.9029047
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
## t1*      4.5 0.02862863   0.8829044
```

One of the main advantages to the 'boot' package over the 'bootstrap' package is the nicer formatting of the output.

Going back to our original code, lets see how we could reproduce all of these numbers:


```r
table(sample(x,size=length(x),replace=T))
```

```
## 
## 0 2 5 6 9 
## 3 1 1 4 1
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
## [1] 0.0217
```

```r
se.boot
```

```
## [1] 0.9205469
```

Why do our numbers not agree exactly with those of the boot package? This is because our estimates of bias and standard error are just estimates, and they carry with them their own uncertainties. That is one of the reasons we might bother doing jackknife-after-bootstrap.

The 'boot' package has a LOT of functionality. If we have time, we will come back to some of these more complex functions later in the semester as we cover topics like regression and glm.

