# Logistic regression model II
##### Last updates: 07/15/2020

Summary
======

Finalize the model
======
## Step1: Explore the feastures that have predictive power on phasing error

##### four features distribution by class
![alt text](figures/4features.png "Logo Title Text 1")
##### four features distribution
![alt text](figures/4features_all.png "Logo Title Text 1")
##### dicatomize features
![alt text](figures/feature_dicatomize.png "Logo Title Text 1")

## Step2: Check mean-variance relationship & mean-square-error relationship
binned data: data was binned by certain characteristics and then splitted randomly into 50 subgroups


\overline{R}

![alt text](figures/logreg_step2_plan.png "process of calculating the statistics")
**sample_varaince**
![alt text](figures/sample_variance.png "sample_variance")
**standard error for binomial**
![alt text](figures/se_binomial.png "standard error for binomial")

## Step3: Check sample variance and standard error for binomial on mean phasing error rate for binned data


## Step3: logistic regression
**1. Intercept model**

**2. full model**
##### logistic regression model shows that all four features (IVs) have significant affect on error (DV)
![alt text](figures/logistic_full.png "logistic_full")
**3. overfitted model (with all intxn terms)**
