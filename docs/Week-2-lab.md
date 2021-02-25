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
## 1 2 3 6 7 8 9 
## 2 1 3 1 1 1 1
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
## [1] 0.0073
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
## [1] 2.700943
```

```r
UL.boot
```

```
## [1] 6.313657
```

Method #2: Simply take the quantiles of the bootstrap statistics


```r
quantile(xmeans,c(0.025,0.975))
```

```
##  2.5% 97.5% 
##   2.6   6.3
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
##    [1] 4.7 4.0 4.0 4.2 5.2 5.5 3.9 5.0 5.0 5.1 5.1 4.5 6.1 3.6 5.4 5.2 5.0 3.6
##   [19] 2.9 5.2 5.8 3.8 5.1 2.8 3.6 5.7 2.1 4.6 5.2 3.9 3.9 2.4 4.3 3.5 5.2 4.1
##   [37] 5.0 6.2 5.0 4.2 4.7 4.7 6.2 3.0 6.3 4.0 3.5 3.6 6.1 3.9 4.3 4.6 4.2 5.0
##   [55] 4.1 3.4 2.0 4.4 4.0 4.1 5.1 3.7 5.8 5.0 5.9 3.5 4.9 5.5 4.7 4.8 3.6 4.2
##   [73] 4.0 5.1 4.6 5.0 5.9 5.5 4.7 3.1 4.7 5.8 6.3 5.8 5.5 5.2 4.5 4.2 4.1 5.2
##   [91] 5.8 5.0 3.5 5.9 3.9 5.4 4.0 4.0 4.8 4.2 4.1 4.2 4.5 2.9 6.2 6.0 5.2 5.2
##  [109] 5.4 5.0 4.3 5.8 5.3 4.7 3.3 4.7 4.0 5.7 2.4 4.9 4.4 5.9 5.1 4.3 3.5 4.6
##  [127] 4.3 5.9 4.4 3.4 3.4 3.9 3.1 4.7 4.0 5.0 6.0 5.0 6.5 5.7 4.3 4.9 4.9 3.9
##  [145] 4.4 5.2 4.7 5.2 3.4 3.4 5.7 4.0 4.8 4.0 4.4 4.1 4.2 3.4 3.6 4.0 5.9 4.9
##  [163] 5.1 4.7 7.7 3.9 4.6 6.7 5.0 4.2 4.0 4.7 5.4 3.9 3.3 4.5 2.9 4.9 3.7 3.2
##  [181] 4.2 2.1 4.9 3.8 5.5 3.6 4.5 3.9 5.2 3.6 5.5 5.6 4.6 3.9 4.4 4.3 5.3 4.7
##  [199] 4.4 5.1 3.9 5.8 4.6 3.7 3.6 4.3 4.4 4.3 4.5 3.7 2.9 4.4 5.3 3.2 4.4 4.4
##  [217] 4.2 4.8 4.5 4.7 4.6 4.7 4.8 3.8 2.9 4.2 4.5 4.2 5.6 4.9 4.6 4.9 4.1 4.4
##  [235] 3.7 2.2 3.6 4.3 4.6 3.4 3.5 4.5 5.5 4.6 3.5 4.0 2.9 3.8 5.5 4.7 4.4 5.2
##  [253] 4.7 5.4 4.5 2.6 5.3 5.7 4.8 4.9 5.0 6.8 5.2 5.0 4.4 5.6 3.4 3.8 2.8 4.9
##  [271] 5.3 5.0 4.8 4.0 3.9 5.3 5.3 4.3 3.8 4.6 4.9 4.4 4.7 5.0 5.0 4.5 3.4 5.8
##  [289] 4.9 5.7 4.1 5.8 7.2 4.4 3.8 3.9 4.6 4.6 4.3 3.9 3.8 4.4 4.3 4.1 3.8 4.0
##  [307] 3.0 3.8 5.6 3.6 5.7 4.6 4.1 3.8 3.8 6.2 5.3 4.2 2.1 5.3 3.4 3.0 4.4 5.5
##  [325] 3.8 4.8 3.6 3.2 4.7 4.5 6.0 2.6 4.7 6.8 5.8 5.1 5.9 4.1 2.5 4.9 4.7 5.0
##  [343] 4.4 5.7 5.3 5.0 3.8 4.2 5.3 3.8 5.2 4.8 3.4 4.6 3.6 3.3 4.1 4.9 4.2 4.3
##  [361] 4.7 3.8 4.3 5.4 4.4 4.8 3.4 3.7 4.0 3.8 4.1 4.5 5.0 3.4 5.2 4.7 4.5 3.9
##  [379] 4.3 4.2 4.3 4.9 4.4 5.3 3.6 6.0 2.8 4.0 6.7 4.1 3.6 2.8 4.5 5.7 4.5 5.1
##  [397] 4.6 5.5 6.9 5.3 5.0 5.9 4.4 5.9 4.2 3.8 3.5 3.6 4.2 3.8 6.0 6.2 4.0 7.0
##  [415] 2.9 3.6 4.3 4.7 2.9 4.7 5.0 3.7 5.0 4.2 4.3 4.4 5.8 5.7 5.8 3.8 4.3 3.9
##  [433] 5.1 3.4 3.6 4.9 5.2 4.7 5.4 5.0 4.9 5.3 6.4 3.6 6.2 6.4 4.5 6.4 3.6 4.6
##  [451] 4.6 2.8 3.1 3.4 5.6 4.9 5.7 4.3 5.4 6.0 3.9 5.5 3.8 5.0 3.3 4.3 4.7 4.9
##  [469] 3.8 3.3 3.2 4.3 3.1 4.8 5.0 5.2 3.1 5.8 2.7 4.2 6.1 3.9 4.1 3.8 4.5 5.1
##  [487] 3.5 3.8 5.1 5.1 4.2 4.4 4.5 2.3 3.7 6.3 4.1 4.6 4.0 4.2 7.3 3.9 4.0 3.0
##  [505] 5.3 6.4 4.5 5.3 4.0 4.5 6.0 2.8 4.0 2.7 2.9 3.8 4.9 5.2 5.1 4.8 3.4 6.1
##  [523] 3.6 4.7 6.9 6.4 5.2 3.1 5.1 5.0 5.8 3.6 4.2 5.2 5.1 5.1 4.8 5.5 4.6 4.2
##  [541] 2.9 4.4 6.4 5.1 5.1 4.4 3.8 3.3 3.9 4.7 4.4 3.9 4.3 4.8 5.6 5.2 5.4 3.5
##  [559] 5.4 5.7 4.8 5.5 3.4 5.1 6.1 4.6 4.7 4.7 4.6 2.9 5.8 3.9 3.4 5.1 3.8 3.4
##  [577] 3.7 3.0 4.8 5.3 4.1 3.8 4.3 6.8 3.6 4.7 5.5 5.6 3.9 5.3 3.5 2.7 4.6 4.9
##  [595] 3.9 3.5 2.2 4.7 3.4 4.2 3.4 6.2 5.1 5.6 4.0 3.0 5.5 5.2 3.9 5.6 3.1 5.2
##  [613] 4.3 5.9 5.1 5.0 3.9 5.2 5.4 5.2 5.4 4.6 4.0 5.0 4.1 3.5 5.6 3.7 3.8 4.4
##  [631] 3.4 4.8 6.0 3.4 4.5 5.1 4.0 4.8 4.8 4.6 3.2 4.5 3.1 4.8 4.4 3.7 4.1 6.6
##  [649] 3.4 3.6 4.8 5.4 5.5 4.5 4.1 3.0 4.0 4.6 4.5 4.8 4.7 4.4 4.8 5.0 4.2 3.1
##  [667] 3.2 4.6 5.9 5.0 4.5 3.6 3.3 4.4 3.1 5.3 3.7 6.1 3.6 4.4 4.3 5.4 6.1 5.8
##  [685] 4.5 5.1 3.5 5.0 4.2 3.8 5.1 4.6 4.3 5.6 5.9 5.7 4.6 3.8 4.4 3.7 5.4 5.1
##  [703] 4.1 4.0 5.8 4.8 5.2 5.9 4.0 5.3 4.8 2.8 4.2 4.0 4.7 5.5 3.2 4.2 5.7 2.9
##  [721] 3.7 3.8 3.8 3.2 5.5 4.2 4.9 3.4 5.0 4.5 3.3 4.5 4.3 5.4 5.8 4.6 5.1 3.8
##  [739] 4.9 5.2 4.9 5.6 4.0 3.6 4.8 3.9 7.0 4.6 4.7 3.1 4.6 4.8 4.1 5.0 4.3 3.9
##  [757] 6.3 3.2 6.2 5.6 3.8 5.6 4.4 4.2 3.0 5.6 5.0 5.0 2.4 6.5 4.3 3.9 4.6 4.9
##  [775] 3.3 4.4 4.2 5.8 2.6 6.1 5.3 4.3 3.0 4.8 3.5 4.6 5.1 4.1 3.7 4.8 5.3 6.8
##  [793] 4.8 5.6 4.5 5.3 5.4 4.3 4.1 4.9 5.5 4.0 4.2 4.8 4.9 6.5 4.4 4.6 3.6 2.4
##  [811] 4.3 3.2 4.6 3.3 4.7 3.5 4.0 5.6 3.4 4.9 6.2 3.9 4.3 4.5 5.8 3.7 4.3 4.6
##  [829] 4.4 4.5 2.9 4.9 5.0 2.6 4.1 4.3 5.8 4.3 5.9 4.1 5.4 5.1 3.8 4.1 2.9 4.9
##  [847] 5.2 4.6 5.3 4.8 4.1 3.7 4.7 3.0 3.4 4.5 5.5 4.5 3.0 6.3 4.2 5.3 4.1 6.8
##  [865] 4.7 4.6 3.7 4.1 4.1 3.9 5.6 5.3 4.2 5.8 4.1 4.3 3.9 5.7 5.2 5.0 4.0 4.7
##  [883] 3.5 3.0 4.4 5.8 5.0 6.0 5.6 5.3 4.6 3.5 6.6 2.8 4.3 5.7 4.2 5.1 4.7 5.7
##  [901] 5.8 5.2 4.5 3.7 4.8 3.2 3.1 3.9 3.4 4.8 5.7 3.5 5.5 5.5 5.6 5.4 3.9 5.3
##  [919] 3.4 5.6 3.2 4.9 3.5 4.9 4.5 5.2 4.6 5.4 4.4 5.0 5.0 4.4 4.3 6.8 4.7 5.1
##  [937] 3.5 5.7 5.7 5.2 4.8 4.9 4.7 3.5 4.8 3.9 4.2 3.4 4.1 3.9 5.9 4.5 5.4 3.2
##  [955] 4.2 3.7 6.1 4.5 4.6 5.0 3.4 2.2 4.6 3.7 3.0 5.1 3.0 4.9 3.7 5.2 5.4 6.0
##  [973] 4.9 4.8 3.3 4.9 3.7 4.5 6.3 3.7 5.2 5.1 4.5 5.1 4.2 5.9 3.8 2.0 4.5 3.1
##  [991] 5.5 4.8 4.9 5.6 3.6 3.1 2.8 3.8 5.3 4.6
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
##   2.8   6.4
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
##    [1] 4.6 4.0 4.0 3.9 4.5 5.1 4.4 4.7 2.4 4.9 4.4 5.0 4.6 5.2 5.3 3.2 4.9 4.3
##   [19] 3.2 4.1 4.0 5.4 5.3 4.7 5.1 5.3 5.1 4.9 4.7 4.3 5.0 6.0 5.2 4.7 4.2 5.1
##   [37] 4.0 5.4 2.9 4.2 3.2 3.9 5.0 4.4 4.6 3.6 4.7 4.7 5.8 4.1 4.3 5.3 4.7 3.4
##   [55] 5.7 3.6 4.1 2.5 5.1 4.1 6.3 3.2 4.6 3.6 3.5 4.8 5.0 4.5 4.7 5.3 5.2 5.3
##   [73] 5.1 5.5 3.0 4.6 2.2 5.2 4.8 5.3 3.9 5.9 4.8 4.2 4.3 4.6 4.4 5.9 4.2 5.5
##   [91] 3.6 2.9 4.1 5.7 5.2 4.8 4.6 3.8 6.5 3.5 4.8 3.6 3.5 5.6 4.2 5.6 4.3 3.9
##  [109] 3.8 3.5 3.7 2.4 4.5 3.9 3.6 4.5 4.6 4.1 2.6 4.4 2.9 4.6 4.6 3.7 3.1 4.0
##  [127] 4.8 4.2 4.9 5.6 3.8 4.1 5.0 3.8 2.7 3.0 4.3 5.1 6.0 4.8 4.8 5.0 5.4 4.1
##  [145] 4.7 4.3 6.0 3.6 5.0 7.2 4.5 3.7 4.6 2.9 3.6 3.5 5.1 3.8 3.7 6.2 3.8 3.1
##  [163] 4.8 4.6 4.3 3.6 2.8 4.7 3.6 3.9 4.9 5.2 5.9 4.9 3.8 4.8 5.1 5.3 4.6 4.7
##  [181] 3.5 3.3 5.2 4.3 4.0 7.1 5.2 4.0 4.7 4.9 5.8 4.8 5.5 4.2 4.7 4.5 4.1 5.4
##  [199] 4.1 5.3 3.7 6.6 6.2 4.3 3.3 4.6 3.4 5.1 5.4 4.8 5.2 4.3 3.9 5.7 5.9 3.0
##  [217] 4.7 4.6 4.3 4.4 4.4 4.4 4.4 4.1 6.0 3.0 3.3 6.1 5.5 4.5 4.0 5.2 4.4 7.1
##  [235] 4.3 4.8 4.3 4.3 4.4 3.6 4.8 3.9 5.9 4.9 5.0 2.9 3.8 4.2 5.4 4.5 5.3 3.6
##  [253] 4.2 3.4 6.2 4.3 4.5 5.0 4.4 4.5 3.4 3.4 5.5 2.9 4.0 4.7 3.1 4.2 5.6 4.1
##  [271] 4.4 4.8 4.0 4.8 6.1 4.8 4.8 4.1 4.2 5.6 3.8 4.6 4.9 3.6 4.4 5.3 4.5 4.9
##  [289] 4.3 4.8 3.7 4.3 4.7 4.4 5.6 5.0 5.1 5.5 3.7 4.2 4.6 3.6 7.0 3.4 5.2 3.9
##  [307] 4.5 2.3 5.0 5.4 3.1 4.8 5.8 4.2 5.1 3.9 4.1 3.7 4.9 3.1 5.0 3.9 4.0 4.5
##  [325] 3.7 3.8 4.8 5.7 5.8 3.4 3.3 5.7 3.9 4.4 6.0 3.0 2.8 4.9 5.0 5.1 4.3 4.4
##  [343] 2.6 4.6 4.4 4.0 4.5 2.2 3.0 3.5 3.6 5.4 4.8 5.0 5.1 3.7 4.8 4.9 4.1 5.5
##  [361] 5.2 5.5 4.7 6.8 4.5 6.0 3.2 4.6 5.1 4.7 5.4 5.1 5.0 5.6 4.0 4.7 3.5 5.1
##  [379] 5.6 2.3 1.9 5.4 3.6 5.3 5.9 3.9 4.5 4.2 3.9 4.2 4.1 4.1 4.9 5.7 3.1 5.3
##  [397] 5.3 3.9 4.9 3.7 4.3 3.6 5.6 4.1 6.1 5.0 6.3 4.4 4.5 4.1 3.6 5.8 5.4 6.1
##  [415] 4.4 5.3 6.0 4.2 5.1 5.0 4.7 5.6 3.5 6.1 4.3 4.4 5.5 2.4 6.0 4.5 5.5 6.2
##  [433] 4.3 3.9 5.1 3.4 3.6 5.7 3.5 5.4 4.6 5.1 5.2 4.7 4.4 4.2 4.4 3.8 3.8 3.6
##  [451] 5.0 5.3 4.4 4.5 2.9 3.6 3.2 3.8 5.4 4.3 4.1 4.9 3.7 5.5 4.9 4.5 2.9 4.9
##  [469] 6.2 5.4 3.9 3.9 3.8 4.8 3.5 4.3 5.7 5.8 3.0 5.1 3.7 4.9 3.6 6.3 5.9 4.2
##  [487] 5.4 5.0 5.1 4.6 3.8 4.1 3.2 4.5 4.2 3.5 4.1 4.0 6.0 4.1 3.3 2.5 3.9 4.0
##  [505] 3.3 6.1 4.4 5.7 5.5 4.1 3.9 4.3 6.2 3.9 3.3 5.2 3.9 2.5 3.7 4.9 5.9 4.4
##  [523] 5.4 4.8 5.1 4.6 4.9 6.4 3.4 4.9 4.6 3.2 5.7 4.6 5.5 3.7 4.9 4.4 6.2 5.9
##  [541] 6.5 4.5 4.4 4.6 3.8 3.9 4.4 3.4 5.3 5.5 5.9 5.3 5.2 3.1 5.8 4.5 4.2 5.7
##  [559] 4.4 3.3 4.4 3.4 2.9 3.7 4.1 4.9 4.7 5.1 4.4 4.0 3.6 4.4 4.3 3.4 3.7 5.6
##  [577] 3.8 3.0 3.9 4.0 5.9 3.6 3.9 5.1 4.5 4.6 4.2 6.2 3.9 5.4 5.0 4.4 4.1 4.5
##  [595] 5.5 3.8 5.4 4.4 5.5 4.7 3.9 4.3 5.6 3.2 4.8 3.5 4.4 6.0 3.0 4.6 3.1 3.6
##  [613] 5.8 3.3 3.5 4.3 4.1 3.7 4.6 3.4 4.4 5.4 4.8 4.2 6.6 5.1 3.9 3.7 3.4 4.5
##  [631] 3.7 5.1 4.1 5.9 6.1 5.4 4.3 5.5 4.4 6.7 4.6 3.7 4.0 4.9 3.9 3.4 5.6 3.2
##  [649] 5.1 4.3 5.1 5.4 4.4 3.8 3.5 5.2 3.3 4.8 4.8 5.0 5.5 3.9 4.6 4.4 2.9 5.3
##  [667] 4.1 4.8 4.2 3.3 3.3 4.7 4.4 4.2 2.9 4.5 5.5 4.9 4.6 5.2 4.2 5.1 5.5 4.3
##  [685] 6.1 5.9 4.3 6.0 3.5 4.6 5.6 5.0 4.7 4.3 5.7 5.0 4.7 4.5 6.1 5.1 5.4 5.0
##  [703] 4.2 4.9 4.7 4.1 4.0 3.7 3.8 5.7 3.3 4.7 3.4 4.5 4.1 4.2 3.3 4.4 4.5 4.2
##  [721] 5.4 5.7 3.5 4.6 4.6 5.4 5.4 3.5 4.6 4.3 3.7 4.5 4.8 3.4 4.3 5.0 3.2 5.1
##  [739] 5.2 3.6 3.6 4.8 3.3 5.4 5.3 4.9 4.8 6.9 4.4 4.2 4.4 4.1 4.9 5.1 2.7 3.5
##  [757] 5.5 4.6 4.1 3.2 5.7 4.9 4.0 4.2 3.9 4.6 4.0 3.9 5.0 6.7 5.5 4.8 3.2 5.8
##  [775] 3.3 3.5 3.9 3.6 3.8 3.2 3.2 4.4 5.0 4.4 3.8 4.7 3.4 4.0 4.6 3.9 5.1 4.9
##  [793] 2.4 4.8 3.6 2.5 5.6 4.4 4.6 5.2 5.4 4.0 4.0 4.6 4.0 3.7 2.8 4.3 3.4 3.7
##  [811] 4.9 4.3 4.5 6.9 6.4 5.0 4.2 3.7 3.7 5.3 4.4 4.1 3.8 4.8 4.3 3.6 3.9 5.2
##  [829] 5.4 4.9 5.2 3.8 3.7 4.2 4.7 3.9 3.5 5.3 5.3 4.7 3.7 4.2 5.8 3.7 4.2 5.1
##  [847] 2.8 5.0 3.9 3.9 5.0 5.2 6.0 5.1 5.9 3.9 5.7 7.6 3.9 2.1 3.9 4.2 5.2 4.5
##  [865] 5.7 5.0 2.7 5.4 5.1 5.5 5.1 6.3 3.2 4.7 4.8 4.3 5.5 3.9 4.8 3.9 5.7 3.6
##  [883] 5.5 5.1 4.5 4.0 5.3 3.9 4.4 5.9 4.2 3.8 4.3 3.9 3.5 4.9 4.2 4.6 5.9 6.3
##  [901] 4.4 3.2 3.7 4.7 3.8 5.0 3.1 3.7 3.9 4.0 4.1 4.9 5.5 5.2 4.5 4.6 4.5 3.9
##  [919] 5.0 3.8 4.3 6.0 5.7 4.8 3.3 4.8 4.8 4.1 6.4 4.5 3.7 4.6 4.9 3.9 3.5 4.4
##  [937] 5.0 3.8 6.4 5.6 5.5 3.9 4.1 5.2 5.2 4.3 3.9 4.6 4.8 6.1 4.9 3.6 2.8 4.2
##  [955] 5.8 4.2 4.7 5.4 6.5 6.2 3.0 5.5 4.3 5.5 5.6 2.9 4.0 4.1 7.0 4.5 2.0 5.1
##  [973] 3.9 4.4 4.9 4.1 3.9 4.4 4.9 4.9 3.7 3.5 5.8 5.4 4.4 3.8 5.2 4.3 5.3 5.0
##  [991] 5.0 5.3 4.9 4.7 5.4 3.2 4.4 3.2 4.5 6.5
## 
## $func.thetastar
## [1] 0.016
## 
## $jack.boot.val
##  [1]  0.49531680  0.40000000  0.29799427  0.17323529  0.09147727  0.03684211
##  [7] -0.20306407 -0.20862534 -0.31929825 -0.50426136
## 
## $jack.boot.se
## [1] 0.9334415
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
##    [1] 5.2 3.2 4.2 5.0 3.9 3.8 3.3 4.6 5.4 4.4 3.4 5.4 3.6 3.1 3.8 5.5 4.5 3.8
##   [19] 6.4 5.3 4.6 4.4 3.4 5.8 4.8 4.1 6.2 3.2 3.5 5.1 3.3 4.3 4.3 2.9 4.2 2.6
##   [37] 5.8 3.4 5.2 5.4 4.5 4.3 5.5 3.0 3.2 3.9 4.7 4.3 6.7 4.1 4.8 4.4 4.9 2.8
##   [55] 4.9 4.8 4.9 4.0 5.0 4.5 4.4 6.1 3.5 2.8 4.6 4.0 5.0 4.9 5.8 5.1 4.2 5.4
##   [73] 5.1 4.8 3.7 6.1 2.9 5.6 3.2 5.8 2.9 5.0 3.8 4.8 2.7 3.7 4.3 4.1 3.6 5.4
##   [91] 6.1 4.1 4.5 5.4 4.8 5.4 6.0 3.6 3.9 5.4 4.6 3.0 4.8 4.6 4.1 3.7 5.8 4.5
##  [109] 3.7 4.5 4.2 5.2 3.9 3.3 5.8 4.8 5.2 3.8 4.3 1.8 6.3 3.8 4.0 3.9 4.3 4.9
##  [127] 5.2 3.0 7.1 3.1 4.3 4.8 6.8 4.3 4.5 4.8 5.2 3.4 4.6 5.5 2.9 3.8 4.3 3.7
##  [145] 4.9 5.5 4.8 4.8 5.1 4.0 3.4 5.2 4.6 5.7 4.2 4.3 3.8 5.5 6.0 4.4 4.4 5.7
##  [163] 5.7 4.5 4.4 5.4 5.0 4.7 5.4 3.3 4.0 6.0 4.8 3.7 5.7 5.2 4.3 2.7 4.6 5.1
##  [181] 4.3 4.6 3.5 5.3 5.5 4.1 2.9 2.7 3.6 3.8 3.8 4.6 4.6 3.9 5.2 3.9 4.9 5.0
##  [199] 5.4 4.3 3.9 4.1 5.5 5.4 6.1 5.7 6.4 4.6 4.5 4.6 6.3 2.9 3.8 5.2 3.8 4.1
##  [217] 5.0 4.1 6.0 4.3 4.3 4.4 2.9 4.8 4.3 4.7 4.4 3.3 5.7 5.4 2.2 4.5 4.4 4.8
##  [235] 5.1 5.3 2.2 5.2 4.2 4.9 2.9 4.0 5.2 3.8 5.2 2.4 4.9 4.9 5.3 5.0 6.3 5.6
##  [253] 5.5 4.1 5.4 4.5 3.1 4.7 3.9 4.3 4.9 5.9 4.9 4.5 4.5 2.9 3.6 4.7 4.1 6.0
##  [271] 4.5 4.2 3.8 3.3 4.1 3.9 4.4 5.3 3.8 5.0 5.0 3.1 5.8 5.6 4.2 3.8 2.6 3.4
##  [289] 5.5 3.2 3.8 4.0 4.5 3.8 5.2 3.5 3.1 4.0 4.1 4.8 6.2 4.1 3.7 3.7 4.8 3.8
##  [307] 3.2 4.7 3.5 3.1 4.7 5.1 4.8 4.7 4.3 5.6 4.9 3.1 3.2 4.9 4.6 4.7 5.7 4.6
##  [325] 4.2 4.5 4.3 5.5 4.3 5.1 3.4 5.0 4.0 4.0 4.3 2.7 4.6 5.8 4.6 5.0 4.8 5.7
##  [343] 3.6 4.8 4.6 5.4 5.0 6.0 5.7 4.2 4.5 3.0 5.1 4.3 4.9 3.9 5.2 4.4 3.4 4.9
##  [361] 4.3 4.0 5.6 5.3 5.4 1.9 3.1 4.3 5.6 3.8 6.5 5.7 6.8 4.8 3.7 3.4 4.7 3.4
##  [379] 4.0 3.3 2.6 4.5 3.7 2.4 3.8 4.4 3.3 4.2 4.2 3.4 4.2 4.6 3.3 5.1 4.5 5.4
##  [397] 4.8 4.1 4.8 4.9 5.7 6.4 4.2 4.7 4.0 5.8 5.9 4.6 3.5 5.2 5.4 4.0 3.4 4.8
##  [415] 3.4 3.7 3.2 3.4 5.5 4.1 3.7 4.9 3.2 3.4 2.9 4.5 5.9 3.8 4.2 4.9 6.7 4.2
##  [433] 5.0 4.2 4.8 5.0 3.7 4.5 4.9 6.1 3.8 5.3 3.2 4.7 5.4 4.1 4.0 3.9 5.2 3.5
##  [451] 3.9 6.0 4.6 5.0 3.4 5.8 3.4 3.1 3.8 4.2 4.0 3.5 4.3 5.1 4.5 4.9 3.1 4.9
##  [469] 7.0 4.8 4.5 5.2 5.3 4.2 5.2 4.0 3.9 5.0 5.0 3.7 4.9 4.4 5.4 6.1 4.4 3.9
##  [487] 5.2 4.2 2.6 6.3 5.5 5.0 4.8 3.6 4.7 4.7 4.9 4.8 4.3 4.8 4.9 4.2 5.0 4.9
##  [505] 4.0 3.9 4.2 4.5 4.5 4.2 4.3 4.8 3.0 4.1 3.3 3.7 4.6 3.0 6.2 5.8 4.1 3.8
##  [523] 7.0 4.9 4.2 5.4 2.4 4.8 4.6 4.6 5.6 3.8 4.6 3.8 6.2 5.5 6.5 4.5 5.7 5.1
##  [541] 6.3 3.9 3.8 2.4 4.3 3.7 4.6 3.6 3.7 5.2 3.1 6.5 5.3 5.5 4.9 5.1 3.9 3.2
##  [559] 5.6 2.8 3.6 2.9 3.4 5.6 5.4 4.2 6.7 4.7 3.9 4.4 5.1 5.0 5.2 5.0 4.0 4.8
##  [577] 6.5 4.7 4.1 4.8 4.9 4.1 4.8 3.8 5.6 4.5 3.9 3.9 4.4 3.7 4.2 5.4 4.6 3.9
##  [595] 3.7 3.6 5.1 4.9 3.7 6.1 4.7 4.9 4.2 3.4 4.2 4.8 5.9 3.6 5.3 3.0 4.3 4.6
##  [613] 4.2 5.1 4.3 4.4 5.4 5.6 4.7 3.2 3.9 4.8 4.2 6.3 5.0 3.2 4.5 3.6 4.2 5.9
##  [631] 6.3 3.3 3.9 5.0 4.5 2.2 3.8 5.5 5.5 4.0 5.5 3.4 4.0 5.0 4.8 4.5 4.9 4.0
##  [649] 4.6 5.0 3.9 2.8 4.5 6.1 4.8 4.5 3.6 4.7 5.4 3.7 5.2 4.2 5.4 5.7 3.3 2.2
##  [667] 2.9 4.1 4.5 6.2 5.0 4.2 4.9 3.6 5.6 4.8 5.7 5.9 4.0 2.6 4.5 5.6 5.1 4.8
##  [685] 4.9 5.3 4.0 4.1 5.4 4.6 5.7 4.7 3.3 5.5 3.9 4.8 4.5 4.4 5.3 4.8 6.1 4.6
##  [703] 3.2 5.8 6.7 4.3 4.7 4.2 5.8 5.0 6.2 3.1 6.2 4.4 3.4 3.9 4.6 4.1 6.2 4.5
##  [721] 4.3 4.0 4.4 2.8 5.0 4.4 4.2 6.3 3.4 5.7 3.5 3.9 3.0 4.3 3.8 4.6 4.0 4.8
##  [739] 4.7 4.9 4.2 4.5 3.5 5.2 3.9 5.3 4.6 5.3 2.5 4.9 4.0 3.9 3.1 4.2 3.8 6.1
##  [757] 4.5 5.6 5.0 4.7 3.0 4.8 3.7 4.3 5.5 5.8 4.5 3.6 5.3 4.1 3.7 5.3 4.1 4.9
##  [775] 5.8 4.6 3.0 3.0 4.7 3.9 4.2 5.5 3.8 4.9 4.1 5.1 5.0 5.8 5.0 5.1 4.8 4.5
##  [793] 3.6 5.3 5.7 3.7 4.2 4.1 2.9 5.2 3.3 5.3 4.4 5.0 6.0 5.0 5.2 5.7 5.4 5.3
##  [811] 4.9 3.9 3.9 5.3 3.7 3.8 6.2 3.7 3.4 4.8 4.3 5.4 5.9 6.3 4.5 2.3 4.1 3.6
##  [829] 6.7 5.6 4.4 5.2 3.7 4.6 2.9 3.6 3.9 3.2 4.6 5.8 4.0 4.7 3.5 2.3 4.1 5.1
##  [847] 3.6 5.3 4.6 3.5 2.7 4.9 5.0 5.3 4.2 4.0 4.6 5.9 3.7 4.4 4.6 6.4 4.7 4.4
##  [865] 3.5 3.5 4.0 3.1 5.3 5.6 3.3 5.7 3.4 4.3 5.3 4.1 4.0 5.4 5.5 5.9 3.7 3.3
##  [883] 4.5 5.4 5.0 6.1 5.1 5.3 3.6 5.4 3.6 4.0 6.1 4.8 3.9 3.9 6.2 5.9 5.7 5.6
##  [901] 5.4 5.7 5.0 5.0 4.6 5.5 4.6 4.5 5.1 4.3 4.9 4.4 3.2 4.9 5.6 4.7 4.3 4.3
##  [919] 4.5 3.9 6.1 5.4 5.9 4.7 2.9 4.7 4.7 4.8 5.2 3.4 6.1 3.1 3.7 3.2 5.0 3.0
##  [937] 5.4 3.9 5.1 4.5 4.8 5.1 5.5 4.3 4.6 4.8 5.5 5.5 4.1 4.4 4.8 2.8 4.4 4.7
##  [955] 4.4 4.1 6.4 5.3 5.3 4.4 5.6 4.4 4.6 5.1 4.7 5.9 5.3 4.7 4.8 4.6 3.8 6.1
##  [973] 4.0 3.7 5.1 5.3 4.6 5.2 4.1 3.5 4.4 5.0 3.8 4.1 4.0 4.9 5.3 3.3 4.6 4.5
##  [991] 3.7 5.3 3.3 3.7 5.1 4.5 5.0 4.7 5.0 3.6
## 
## $func.thetastar
## 72% 
##   5 
## 
## $jack.boot.val
##  [1] 5.50 5.44 5.20 5.30 5.10 5.10 5.00 4.80 4.60 4.50
## 
## $jack.boot.se
## [1] 0.9555082
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
## [1] 1.295119
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
##   4.556370   9.572778 
##  (1.967615) (4.370386)
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
## [1] 1.8792439 0.3487101 0.7505952 1.1844423 0.5622747 1.1736204
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
##    [1]  0.063544101  0.850958255  1.000295264  2.200893989  0.629688519
##    [6]  1.824446003  0.818022809  1.271600457  1.092251082  0.626078015
##   [11] -0.247760537  1.634521424  0.074658534  0.353570103 -0.461250427
##   [16]  0.738996599  1.462486016 -0.562731135  0.895493651  1.129055641
##   [21]  0.361400328  1.170166910  1.924267053  0.502083069 -0.255102987
##   [26]  1.432479750  1.829007795  0.358483839  1.232548854  1.132928970
##   [31]  1.312295570  2.413857591  2.082589419  0.285487869  1.184305683
##   [36]  0.589234153  0.789845499  1.835701120  1.084335123  0.548106840
##   [41]  0.945417508 -0.236621134  0.970912075  1.202577706  2.522538939
##   [46]  1.412000794  1.110524120  1.004444000 -0.408943248  1.472816012
##   [51]  1.001158901  1.215285309  1.277489402  1.704823276  1.781202196
##   [56]  1.334576053  1.141250388  0.977384821 -0.077316162 -0.169181341
##   [61]  0.984359139  0.909543182  1.583547305  0.700667426  0.319191852
##   [66]  0.817860311  0.659642602  0.177235842  1.213188459  2.198800538
##   [71] -0.037470801  1.243276559  1.116026797  1.247419902  0.801768968
##   [76]  0.222669145  2.003489579  1.012393676  1.262994281 -0.121094036
##   [81]  1.715352415  1.534769983  0.342450838  2.418610940  0.706989111
##   [86]  0.355316354  0.433982748  0.729899338  0.974606724  0.963217521
##   [91]  0.858758498  1.771957091 -0.354201848  1.228011565 -0.192977662
##   [96]  0.430748226  1.416470683  1.866262117  1.324197090  1.699299510
##  [101]  1.719367469  0.852870035  1.642258146  0.220330032  1.553802162
##  [106]  0.101847545  0.378824492  1.129591633  0.634478197  0.825652503
##  [111] -0.390328313  1.434934905  0.904458671  1.202577706  1.242484774
##  [116]  1.480210240  0.776144406  0.582859484  1.215565326  1.282351384
##  [121]  0.193455291 -0.101245376  0.855910975  1.499274160  0.931482197
##  [126] -0.222959049  1.780848986  1.357353813  0.905138439  1.606001657
##  [131]  0.548989066 -0.146588901  0.499640785  1.545081615 -0.059478339
##  [136]  1.711118871  1.636373230  1.280925895  0.792176402  1.548237468
##  [141]  1.546878205  1.272804554  1.460096125  2.079339512  1.524002773
##  [146]  0.516197056 -0.683988155  0.001182501  0.394148523  1.062235029
##  [151]  0.562658653  1.136652686  0.512049322  0.236671833  1.801179292
##  [156]  1.186181723  1.329015513  1.725740058  0.933827880  0.857191668
##  [161] -0.328184436 -0.143071898  0.190591980  0.073515419  0.473498042
##  [166]  2.148703635  0.163524215  1.005599781  1.566369895 -0.254889279
##  [171]  0.205436192  1.275500302  1.118408848  1.039653474 -0.058071756
##  [176]  0.434257015  0.190749984 -0.187974422  0.753121386  1.210299357
##  [181] -0.164223118  1.276530291  0.808001301  0.099723875  1.207847507
##  [186]  1.169963829  0.959995711 -0.626237602 -0.058907958 -0.212795692
##  [191] -0.343766694  1.953176489  1.490124982  0.746934036 -0.039525376
##  [196]  0.796931504  1.747713473  1.130881300  0.410292878  0.471107141
##  [201] -0.249035892  1.228434162  0.218318919 -0.221637344  1.181824772
##  [206]  1.384383849  0.299942254  1.169346060  1.438319524  0.940725563
##  [211]  0.925947570  1.408107926  1.148903704  1.464786644  1.478893441
##  [216]  1.567431239  0.076626387  0.045029101  1.202296572  1.481259480
##  [221]  1.604476173  1.651045901  0.723725436  0.752207552  1.643336937
##  [226]  1.228467022  0.694001980  0.628517424  1.328953724  1.417455944
##  [231]  1.771524625  1.090016156  0.819985352  1.308457606  2.276323568
##  [236]  0.613097161 -0.126591283  1.390882966  0.426025771  1.273886135
##  [241]  1.268845992  1.314444160  0.797104766  1.968420860  0.798170466
##  [246]  1.466428942  0.791955709  1.204888678  0.087558134  1.277592666
##  [251]  1.306796410  0.390130876  1.817013917  0.135712257  2.035857029
##  [256]  1.012495617 -0.207976666 -0.646799971 -0.446140461  1.123851905
##  [261]  0.464724734  0.759895114  1.042529546 -0.322525780  1.087163484
##  [266]  1.067093286  0.243301878  1.256780790  0.713387599  0.462268676
##  [271]  1.858680789  0.576852048  0.063870869  1.160427481  0.723013988
##  [276] -0.242851296  1.092848691  1.382191623  1.093269003  0.855100523
##  [281]  1.078525956  0.179212611  1.083884822  0.025427653  1.408954805
##  [286]  1.717552861  0.804969094  0.763092085  0.769307164  0.888535025
##  [291]  0.502760427  2.080133256  0.110923124  1.070227771  1.041329015
##  [296]  1.524710594  1.480822306  0.829952538  0.086470392  0.838818349
##  [301]  0.421341630  0.594752549  1.051106941  1.227781846 -1.899134473
##  [306]  0.338544992  1.170086574  1.228441901  0.492876473  0.001277120
##  [311]  1.149619398  0.231085694  1.692163068  0.170339109  2.367686261
##  [316]  0.316052110 -0.147518435  1.558656957  0.241514286  0.774729601
##  [321]  1.721434041  1.356872226  1.054967997  1.042529546  0.725605168
##  [326]  1.893219594  1.318060220  1.343053717  1.289918414  0.068141274
##  [331] -0.508840796  0.517005889  0.562343935 -0.045016222  0.776790044
##  [336]  0.994332549 -0.486980547 -0.230581220  0.854122117  2.252842060
##  [341]  0.790318291  1.007650198  1.145284705  0.306990569  0.428596375
##  [346]  0.968046378  1.025119071  0.925517682  0.023078649  1.229063638
##  [351]  0.850500050 -0.527120571 -0.130079250  1.792096610  1.192057927
##  [356]  0.385179763 -0.243187130  0.697597135  1.448116428  1.348972568
##  [361]  0.723269343  0.105279661  0.501485547  1.887443645  0.736719042
##  [366]  1.457977184  1.448160281  0.827485214  0.138700283  0.677076463
##  [371]  1.312640755  0.214262800  1.303510272  1.285306146 -0.038244709
##  [376] -0.073393437  0.745198299  0.086707975  1.315731315 -0.144343532
##  [381]  0.802133583  1.291381691  0.885663895  1.237879490  0.668747476
##  [386]  1.003875990 -0.013018945  1.229789626 -0.171249026  0.418861285
##  [391] -0.387248825  0.256666538  1.637459936  0.387601331  0.561504479
##  [396]  1.089175126  1.541093270  1.220827079  1.233523160  1.883695801
##  [401]  1.219433644  0.618388518  1.532342831  0.971595335  1.064582753
##  [406]  0.093685742  1.115660830  0.149419618  1.454138711  0.745654000
##  [411]  1.755441228  1.535791943  0.574266877  0.567866140  1.678416108
##  [416]  0.898175767 -0.072869651  0.216501035  0.252164713  0.211751592
##  [421]  0.365472886 -0.399845673  1.690263130  0.167613871  0.267898159
##  [426]  0.590557357  1.376296933  0.324935807  1.134779637  0.519334085
##  [431]  0.423774585  0.707911025  1.512799520  0.966484330  1.518151265
##  [436]  1.062611363  0.245153849  0.371056470  0.375173201  1.616189819
##  [441]  1.605437369 -0.152783358  1.080315342  2.304271248  0.008214810
##  [446]  0.651737956  1.495900601  0.601919982  2.018684682 -0.174669061
##  [451]  0.855757654  0.825352753  1.170784134  1.536272556  0.927471345
##  [456]  1.222770200  1.016834212  1.323790721  1.408363350  0.208959662
##  [461]  1.735298077  0.913921888  0.362786766  0.708977260 -0.347998472
##  [466]  1.151552143  0.423003210  1.357401355  1.270476480  0.152834178
##  [471] -0.279991663  0.134152381  1.394224172  1.143028216  0.021311372
##  [476]  1.677103176  1.384951597  0.536311973  0.540006497  1.431465075
##  [481]  0.217183844  1.159978441  1.256302616  0.917726621  0.011825120
##  [486]  1.869244304  1.693302026  0.174364969 -0.068924846  0.707796032
##  [491]  1.594421801  0.218969350  1.126088386  0.039824648  0.003244240
##  [496]  1.329713082  1.243476241  1.877506703  1.037088655  1.078843384
##  [501]  1.410747804  1.302780883  0.920628959  1.055517786  0.674237816
##  [506]  1.186143412  1.067005167  0.550560263 -0.164498968  1.133683084
##  [511]  1.929235209 -0.760709161  1.145765833  0.209275651  0.316197549
##  [516]  1.026817613  0.482832903  1.225034754  1.347183281  0.317721429
##  [521]  0.413089144  1.134079630  1.633911814  1.027166387  1.664966840
##  [526]  0.835306522  1.247396220  0.692159954  0.551582549  0.125697421
##  [531] -0.041942686  1.337038638  2.089702648 -0.025576830  0.659441696
##  [536]  1.610102624  0.979120474  0.798917167  0.122890005  1.044128590
##  [541]  0.946718700  0.778739588  0.528660250  1.277407499  0.222588025
##  [546] -0.034741892  1.698697732  1.338157547  0.597701878  1.519339933
##  [551]  1.004293237  0.631042478  1.234843935  0.955835548 -0.197723352
##  [556] -0.077602915  1.412000794  1.843727783  0.188732944  1.748187416
##  [561]  1.748689530  0.314462559  0.167039622 -0.046490676  1.514989359
##  [566]  1.567595550  0.599457143  0.773209580  1.332238544  1.269884104
##  [571]  0.922743992  0.527744566  0.631482296  0.559189937  0.553630840
##  [576]  0.640373085  0.089911629  0.975922105  0.923166329  0.574793488
##  [581]  0.460076496 -0.356532865  0.353538666  1.347392302  1.293418201
##  [586]  0.897985746  1.384165594  0.658365213  1.555517649  0.609071687
##  [591]  1.045116473  2.099455614  1.825791984  0.006938723  0.396929774
##  [596]  1.401143521  0.238737681  0.409881239  0.320476594  1.623994256
##  [601]  1.368506610  0.971988249  0.025854856 -0.234077719  2.096664881
##  [606]  1.263144353  1.748750774  0.650459181 -0.421809166  0.679808170
##  [611]  1.003217193  1.266663577  0.790317943  1.018111656  1.299388171
##  [616] -0.459482768  1.785870267  0.569323596  1.401524196  1.063214980
##  [621]  1.616126038  1.065170292  2.212924712  1.661654244 -0.128768991
##  [626]  1.069402908  0.920712090  0.526161024  1.178970888  1.014105301
##  [631]  0.689932559  1.557663244  1.184507897  1.465292402 -0.170725194
##  [636]  0.789635338  0.338009757  0.218100449  2.319377269  0.686194722
##  [641]  1.586256396 -0.251264912 -0.185078890  1.792984283 -0.189178370
##  [646]  1.123466998  1.282258371  1.684403890  0.212123587  1.741449448
##  [651] -0.391203538  1.103875074  0.104164890  0.806398130  1.185849561
##  [656]  1.612751984  1.046761790  2.093349036  1.535556382  0.930600953
##  [661]  0.874697354  0.419694964  0.254153264  1.171228401  0.917073247
##  [666]  0.705491102  0.987481098 -0.469269798  0.425458085  1.574903624
##  [671]  0.294910643  1.619462702  1.602855814  0.778963834  1.167546038
##  [676]  1.325160733  1.259579307  0.978773292  0.075486121  1.567861252
##  [681]  0.931732687  1.381157999  1.530813236  0.033695873  0.845136866
##  [686]  1.153817392 -0.237942969  1.518253490  0.483412761  1.083905026
##  [691]  1.117840859  0.501181659  0.104164890  1.347213239 -0.379691462
##  [696]  1.469272707  1.832093511  1.109549462  0.588855082  1.954827989
##  [701]  0.935476103  0.568209815  0.070421113  1.234356680  0.856819019
##  [706]  0.588940977  1.614655797 -0.050612521  1.626412772  0.437810462
##  [711]  0.194901854  0.856825081  0.452194087  0.161924303  1.316538001
##  [716] -0.028154151 -0.092239480  1.259579307  1.010788742  1.676969296
##  [721]  0.411334596  1.711253988  1.362903704 -0.097552623  0.545694811
##  [726]  1.012495617  1.466682187  0.037869287 -0.010026777  1.161932112
##  [731]  0.543924679 -0.032658000  0.536448471  1.741760970  0.468529098
##  [736]  0.673999845  0.206407245  1.822135874  0.185781053  0.821038081
##  [741]  0.517043432  0.903666174  1.198247501  1.719415344  0.644810767
##  [746]  1.523437231 -0.009693135  0.920328574  0.601937430  2.142772725
##  [751]  1.225821446  1.303841638  0.436316852  0.678447886  1.265608246
##  [756]  1.029892569  1.293595300  0.894657258  0.944156454  1.511533068
##  [761] -0.505893119  1.216492410  0.812143801  1.634946133  0.495125004
##  [766] -0.339073821  0.925774461  0.712442769  1.578335310  0.398971191
##  [771]  1.655692541  1.685202708  0.536767121  1.680852241  0.832013742
##  [776]  0.971404510  0.073777302  1.674110175  1.667663718  1.318681489
##  [781] -0.533108316  0.486695644  0.408403439  0.215657533  0.483340053
##  [786] -0.347005505  0.559923614  1.588729401  0.496597424  1.094102423
##  [791]  0.101734994  0.440587781  0.729166217  0.987363445  1.034346758
##  [796]  0.425849667  2.119442647  0.324037633  2.127280389 -0.052041439
##  [801] -0.675101438  1.579547090  0.267438275  0.774595063  1.200839648
##  [806]  0.529674676  1.175045804  0.539692431  1.172432945  1.030176342
##  [811]  1.161883221  0.925517682  0.893573916  1.277407499 -0.056914559
##  [816]  0.184453847  0.558401967  0.488478414  1.190881769  1.948010833
##  [821]  1.047825547  0.205566651  1.235238127  0.890328280  0.239853232
##  [826]  1.380698916  0.977664920  0.935515706  1.071109016  1.098357563
##  [831]  0.382802114  0.759643195  0.744526797  0.589821308  0.840925226
##  [836]  0.748287683  1.210050283  1.942270607  1.926863879  1.827048026
##  [841]  2.453627667 -0.185952067  1.447098103  0.890028462  0.424443240
##  [846]  0.095050822  1.417353142  0.326905568 -0.273022472 -0.964646384
##  [851]  0.222823245  1.123056004  1.561060764  0.928868787  1.219892319
##  [856]  1.588729401  0.023654101  0.430908704  0.663446111  1.710444452
##  [861]  0.407734076  1.866471438  0.029801841  0.115209988  1.272804554
##  [866]  0.860005033 -0.179584827  0.807232145  0.573136702  1.226931466
##  [871]  1.270069157  1.080252948 -0.008005020  0.945402677  0.277078643
##  [876]  1.688015665  0.838916591  0.168484727 -0.034559455 -0.232669378
##  [881]  1.073683490  0.420117534  0.431576662  1.406221753  0.018238706
##  [886]  1.458959314  1.271615514  1.020418199  1.812674372  1.190958734
##  [891]  0.735559661  0.993462668  1.266069502  1.532183624  1.333090834
##  [896]  1.172010690  1.577132428  0.842762339  0.044211940  1.208214872
##  [901]  0.783189037  0.678650793  1.367062803  1.216767144  0.562586660
##  [906]  1.650688977  0.027287658  1.331450574  0.626673775  0.964359244
##  [911]  1.552302200 -0.410929339  0.948022721  1.233632623  0.665702834
##  [916]  0.717190430  0.902466373  0.346895415  0.042776406  0.036668988
##  [921]  1.548457781  0.350034977  1.303510272  1.339413694  1.195467759
##  [926]  1.305939387  0.367885349 -0.102854649  0.232254708  0.504591488
##  [931]  0.262124678  1.086080640  1.108542281  1.115542471  1.610449250
##  [936] -1.033450748  0.108671211  1.060697824  1.628648249  1.028835339
##  [941]  0.335517261 -0.135351222  0.267539964  1.574292913  0.678447886
##  [946]  0.903147359  1.234328063 -0.008270625  1.155181251  0.932714250
##  [951]  1.458208399  1.537546920  0.396281704  0.792852784 -0.554853531
##  [956]  1.373578056  0.107654168  0.859101325  1.475551769  1.562561385
##  [961]  1.480319526  0.453044371  1.289983673  0.173073384  0.630745568
##  [966]  1.031647165  1.128723757  0.576674667  0.510440230  1.777291897
##  [971] -0.187795637  1.033129697  1.139856065  0.680425900  1.225852976
##  [976]  1.359848622  0.832887238  0.719098256  0.483554642  0.270346956
##  [981]  1.314754025  1.368067886 -0.278179564 -0.006412753  1.523119177
##  [986]  0.660463742  0.426450170  1.335768050 -0.365124523 -0.630649343
##  [991]  0.639285236  0.367468463  0.388870090  0.171629144  1.010788742
##  [996]  2.347430100  1.219433644  0.046649712  0.260636549  0.906776310
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
##   0.47597439   0.24069872 
##  (0.07611562) (0.05381783)
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
## [1]  0.3939440  0.4098091 -0.7690112  0.8811942  0.1542702  0.2283423
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
## [1] 0.009
```

```r
sd(results2$thetastar)  #the standard deviation of the theta stars is the SE of the statistic (in this case, the mean)
```

```
## [1] 0.9083634
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
## t1*      4.5 0.005205205   0.9046711
```

One of the main advantages to the 'boot' package over the 'bootstrap' package is the nicer formatting of the output.

Going back to our original code, lets see how we could reproduce all of these numbers:


```r
table(sample(x,size=length(x),replace=T))
```

```
## 
## 0 1 2 3 4 6 7 8 
## 1 1 1 2 1 2 1 1
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
## [1] -0.0312
```

```r
se.boot
```

```
## [1] 0.9220286
```

Why do our numbers not agree exactly with those of the boot package? This is because our estimates of bias and standard error are just estimates, and they carry with them their own uncertainties. That is one of the reasons we might bother doing jackknife-after-bootstrap.

The 'boot' package has a LOT of functionality. If we have time, we will come back to some of these more complex functions later in the semester as we cover topics like regression and glm.

