# Logistic regression for features on phasing error
##### Weekly updates: 06/17/2020

Background
======
#### What we learned from BEASTIE model performance already
![alt text](figures/model_lowrate.jpg "BEASTIE model phasing error rate 𝜋")
When we assume data with 10% phasing error rate by setting 𝜋 ~ beta(1,10), we observe:
1. BEASTIE outperforms ADAM when true error rate belows 10%.
2. BEASTIE outperforms ADAM when true error rate equal/above 10% with prior distribution centers at true error rate
