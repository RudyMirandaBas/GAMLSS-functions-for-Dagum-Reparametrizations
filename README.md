# Alternative Parameterizations of the Dagum Distribution and Their Implementation in GAMLSS
 
Scripts accompanying the paper *Alternative Parameterizations of the Dagum Distribution: Inference and Regression Modeling*. This repository includes the simulation code and results, the code to reproduce the three real-data applications, and, of course, the code for the three proposed models.
 
## Using the RDQ1 and RDQ2 Models
 
The RDQ1 and RDQ2 model scripts were written as **4-parameter** functions rather than 3-parameter ones, treating the fourth parameter ($\tau$) as the quantile of interest. This was the simplest way to implement them.
 
So, for example, to model the 0.3 quantile, you need to specify it as both the **starting value** and a **fixed value**, as follows:
 
```r
gamlss(Y ~ X1 + X2 + ..., sigma.fo = ~ ...,
       ...,
       tau.fix = TRUE, tau.start = 0.3)
```
 
## Disclaimers
 
**Disclaimer 1:** Unlike what's described above, the simulation scripts for RDQ1 and RDQ2 use the **3-parameter** version, since the quantile in that case was fixed at 0.5, making the general 4-parameter formulation unnecessary. The underlying logic, however, is essentially the same.
 
**Disclaimer 2:** Although the derivatives for the models are included, the **CG** method (and consequently **mixed**) does **not** work — only **RS** does. Since RS is also the default method, be sure to use it explicitly, as shown below:
 
```r
gamlss(Y ~ X1 + X2 + ..., sigma.fo = ~ ...,
       ...,
       tau.fix = TRUE, tau.start = 0.3,
       method = RS())
```
 
**Disclaimer 3:** The applications use a model called **BS2**. Since this model was not developed by us, it is **not included** in this repository. You can find it in the [RelDists](https://github.com/fhernanb/RelDists/) package.
