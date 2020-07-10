# GSD patients summary
##### Weekly updates: 06/24/2020

**P.S.**

*mp: mpileup*

*mc: manual correction*


**Data description:**
Parsing/filtering threshold: total allele counts >= `0`


**Table: ID matching**
|aliqot ID| patient ID| new identifier from DNA VCF|
|--|--|--|
|122687|KG-24|390449|
|122698|CG-46|383579|
|125249|JP-06|383582|
|125260|CP-0416|383583|
|123667|NW-10|383581|
|123375|MA-1|390450|

125249 and 125260 siblings
======
**Target gene - chr16: PHKG2**

Notes from Bill: Found 13 hets near TSS of PHKG2 in the 2 siblings

@Patient: 125249

|type| source VCF|Chr| varianst pos |rs ID |gnomad AF | genotype|location|varaint type        |ALT(mp)|REF(mp)|total(mp)|ALT (mc)|REF(mc)|total(mc)|
|--| --------- |-- |:------------:|:----:|:-:|:-------:|:------:|:------------------:|:----: |:----: |:------: |:-----: |:----: |:------: |
|Non-coding| DNA only | 16 |  30761285   | [rs11309088](http://grch37.ensembl.org/Homo_sapiens/Variation/Explore?db=core;r=16:30760786-30761786;source=dbSNP;tl=0575iecWJsdVkwHU-6360597;v=rs11309088;vdb=variation;vf=312211312)| [0.1613](https://gnomad.broadinstitute.org/variant/16-30761285-CA-C?dataset=gnomad_r2_1) or [0.00003196](https://gnomad.broadinstitute.org/variant/16-30761285-CAG-C?dataset=gnomad_r2_1)|0/1  |  intron varaint |deletion (**CA>C**) or deletion (CAG>C) |||||||
|Non-coding| DNA only | 16 |  30761534   |VCF: rs11406717 NO INFO but VEP: [rs535247818](http://grch37.ensembl.org/Homo_sapiens/Variation/Explore?db=core;r=16:30761035-30762044;source=dbSNP;tl=0575iecWJsdVkwHU-6360597;v=rs535247818;vdb=variation;vf=339107695) |[0.0001710](https://gnomad.broadinstitute.org/variant/16-30761534-TAAA-T?dataset=gnomad_r2_1) or [**0.01481**](https://gnomad.broadinstitute.org/variant/16-30761534-T-TA?dataset=gnomad_r2_1) or [0.00003421]()| 0/1   |  intron variant| Indel TAAA>T or **T>TA** or TAAAA>T|||||||
|Non-coding| DNA only | 16 |  `30762416`    |[rs751600886](http://grch37.ensembl.org/Homo_sapiens/Variation/Explore?db=core;r=16:30761916-30762916;tl=v018KV05GmSzgbuL-6360299;v=rs751600886;vdb=variation;vf=318505695) | [0.00002482](https://gnomad.broadinstitute.org/variant/16-30762416-G-A?dataset=gnomad_r2_1) |0/1   |intron variant |**`pathogenic (splice site) G>A,T`** ||||`23`|`24`|`47`|
|Non-coding|DNA/RNA | 16 |  30765947    |[rs1433965292](http://grch37.ensembl.org/Homo_sapiens/Variation/Explore?db=core;r=16:30765447-30766447;source=dbSNP;tl=q83OcqvHRlC0yN0N-6360591;v=rs1433965292;vdb=variation;vf=344664191) |no information found| 0/1   |intron variant  | (T>G) |13|2|15|||| 
|Non-coding| DNA only| 16|  ~~30768492~~  | [rs71820772]()|no information found|[**1/1**]()   |exon    |**3' UTR varaint , Indel (GCTCTGGC>G)**|||||||0|41|41|
|Non-coding|VEP| 16|  ~~30768493~~  | [rs56107910](http://grch37.ensembl.org/Homo_sapiens/Variation/Explore?db=core;r=16:30767993-30769009;tl=6K6eHdpe5RrkChj1-6360435;v=rs56107910;vdb=variation;vf=313798550)|0.16|  |exon  |**3' UTR varaint , Indel (CTCTGGCC>C)**||||||||||
|Non-coding|DNA/RNA| 16 |  30768499    | [rs77175815](http://grch37.ensembl.org/Homo_sapiens/Variation/Explore?db=core;r=16:30767999-30768999;tl=h0XsP9wDDEB4sf59-6360504;v=rs77175815;vdb=variation;vf=315567850)|no information found|0/1   |exon    |3' UTR variant (C>G)|115|60|175||||
|Non-coding|DNA/RNA | 16 |  30768503    |[rs867637932](http://grch37.ensembl.org/Homo_sapiens/Variation/Explore?db=core;r=16:30768003-30769003;source=dbSNP;tl=hTBYYcWlacR0SAuk-6360507;v=rs867637932;vdb=variation;vf=318873162)| no information found|0/1   |exon    |3' UTR variant (substitution T>A)|44|294|338||||
|Non-coding|DNA only | 16 |  30770074    |[rs11645367](http://grch37.ensembl.org/Homo_sapiens/Variation/Explore?db=core;r=16:30769574-30770574;tl=EuTHLgOybcYGgppD-6360513;v=rs11645367;vdb=variation;vf=312336262)|[0.06613](https://gnomad.broadinstitute.org/variant/16-30770074-C-T?dataset=gnomad_r2_1)|[0/1]()   |   exon  |3' UTR variant (C>T)|||||||
|Non-coding|DNA/RNA | 16 |  30770950    |[rs62622830](http://grch37.ensembl.org/Homo_sapiens/Variation/Explore?db=core;r=16:30770450-30771450;tl=6RRSTBGSeJWtDc0s-6360579;v=rs62622830;vdb=variation;vf=314470340)|[0.1498](https://gnomad.broadinstitute.org/variant/16-30770950-T-C?dataset=gnomad_r2_1)|0/1   |exon    |downstream_gene_variant (missence T>C)|4|6|10||||

conclusion: we should not consider deletion as ASE, or variants close to deletion. SM's filter may not be applicable to our context. We have different goals. 

@Patient: 125260
| source VCF|Chr| varianst pos |rs ID|phasing|location|varaint type        |ALT(mp)|REF(mp)|total(mp)|ALT (mc)|REF(mc)|total(mc)|ALT(ab)|REF(ab)|total(ab)|
| -- | -- |:------------:|:----:|:----: |:------:|:------------------:|:----: |:----: |:----: |:----: |:----: |:----: |:----: |:----: |:----: |
| [DNA only]() | [16]() |  ~~[30761285]()~~   | [rs11309088]()|[0/1]()   |   |[CA>C]() ||||||||||
| [DNA only]() | [16]() |  ~~[30761534]()~~    |[rs11406717]() |[0/1]()   |   |[T>TA]() ||||||||||
| DNA only | 16 |  `30762416`    |. |0/1   |intron  |**`pathogenic (splice site) G>A`** ||||`13`|`28`|`41`||||
| DNA/RNA | 16 |  30765947    |rs1433965292 | 0/1   |intron  |intron variant (T>G) |8|4|12|||| |||
| [DNA only]()| [16]() |  ~~[30768492]()~~  | [rs71820772]()|[**1/1**]()   |[exon]()    |[**deletion (GCTCTGGC>G)**]()||||||||||
| DNA/RNA| 16 |  30768499    | rs77175815|0/1   |exon    |3' UTR variant (C>G)|50|16|66||||
| DNA/RNA | 16 |  30768503    | .| 0/1   |exon    |3' UTR variant (substitution T>A)|10|107|117|||||||
| [DNA only]() | [16]() |  [30770074]()    |[rs11645367]()|[0/1]()   |     |[(C>T)]()||||||||||
| DNA/RNA | 16 |  30770950    |rs62622830|0/1   |exon    |3' UTR variant (missence T>C)|2|3|5|||||||



**chr16: 30768492 deletion (1/1); 30768499 30768503 (1/0)**

Note: those three variants are not in the coding region: 30768499, 30768503, 30770950
![alt text](figures/125249_125260_3vari.png "3 close variants")


122687 and 122698 siblings
======
**Target gene: chr11: SLC37A4**

Patient: 122687

| Chr| varianst pos |rs ID|phasing|location|varaint type        |ALT|REF|total|ALT (mc)|REF(mc)|total(mc)|
| -- |:------------:|:----:|:----: |:------:|:------------------:|:----: |:----: |:----: |:----: |:----: |:----: |
| 11 |  30765947    |rs1433965292 | 0/1   |intron  |intron variant (T>G) |13|2|15|13|2|15|
| 11 |  30768499    | rs77175815|0/1   |exon    |3' UTR variant (C>G)|115|60|175|115|60|175|
| 11 |  30768503    | .| 0/1   |exon    |3' UTR variant (substitution T>A)|44|294|338|44|294|338|
| 11 |  30770950    |rs62622830|0/1   |exon    |3' UTR variant (missence T>C)|4|6|10|4|6|10|
| 11 |  `30762416`    | |0/1   |intron  |`pathogenic (splice site)` ||||`23`|`24`|`47`|

Patient: 122698



123667
======
**Target gene: chr1: GBA**

Patient: 123667

| Chr| varianst pos |rs ID|phasing|location|varaint type        |ALT|REF|total|ALT (mc)|REF(mc)|total(mc)|
| -- |:------------:|:----:|:----: |:------:|:------------------:|:----: |:----: |:----: |:----: |:----: |:----: |
| 1 |  155210451    |rs387906315 | 0/1   |exon  |`pathogenic (frameshift variant, inserG)` |0|8|8|`4`|`1`|`5`|



123375
======
**Target gene: chr1: AGL**

Patient: 123667

| Chr| varianst pos |rs ID|phasing|location|varaint type        |ALT|REF|total|ALT (mc)|REF(mc)|total(mc)|
| -- |:------------:|:----:|:----: |:------:|:------------------:|:----: |:----: |:----: |:----: |:----: |:----: |
| 1 |  100341003    |. | 0/1   |exon  |`pathogenic (deletion of G)` |0|10|10|`5`|`5`|`10`|

### Future
* Bill will do GSNAP for 125249, 125260, not allow mismatch, but allow mismatch in deletion. 
* look at all coding variants, VEP, AF information
* homozygous alternatives
* Bill will look at deletion with VCF shared with siblings, manually remove the deletion.