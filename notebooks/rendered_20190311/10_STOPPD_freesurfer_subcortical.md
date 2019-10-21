# Freesurfer Derived Subcortical Volumes


```r
library(tidyverse)
```

```
## -- Attaching packages --------------------------------------------------------------------------------------------- tidyverse 1.2.1 --
```

```
## v ggplot2 3.1.0     v purrr   0.2.5
## v tibble  1.4.2     v dplyr   0.7.8
## v tidyr   0.8.2     v stringr 1.3.1
## v readr   1.1.1     v forcats 0.2.0
```

```
## -- Conflicts ------------------------------------------------------------------------------------------------ tidyverse_conflicts() --
## x dplyr::filter() masks stats::filter()
## x dplyr::lag()    masks stats::lag()
```

```r
library(lme4)
```

```
## Loading required package: Matrix
```

```
## 
## Attaching package: 'Matrix'
```

```
## The following object is masked from 'package:tidyr':
## 
##     expand
```

```
## Loading required package: methods
```

```r
library(lmerTest)
```

```
## 
## Attaching package: 'lmerTest'
```

```
## The following object is masked from 'package:lme4':
## 
##     lmer
```

```
## The following object is masked from 'package:stats':
## 
##     step
```

```r
library(broom)

df <- read_csv("../generated_csvs/STOPPD_masterDF_2018-11-05.csv",na = "empty") #spreadsheet created by 03_STOPPD_masterDF.rmd
```

```
## Parsed with column specification:
## cols(
##   .default = col_character(),
##   STUDYID = col_integer()
## )
```

```
## See spec(...) for full column specifications.
```

```r
FS <- read_csv('../data/fs-enigma-long_201811/LandRvolumes.csv') #bring in subcortical data, from pipelines
```

```
## Parsed with column specification:
## cols(
##   SubjID = col_character(),
##   LLatVent = col_double(),
##   RLatVent = col_double(),
##   Lthal = col_double(),
##   Rthal = col_double(),
##   Lcaud = col_double(),
##   Rcaud = col_double(),
##   Lput = col_double(),
##   Rput = col_double(),
##   Lpal = col_double(),
##   Rpal = col_double(),
##   Lhippo = col_double(),
##   Rhippo = col_double(),
##   Lamyg = col_double(),
##   Ramyg = col_double(),
##   Laccumb = col_double(),
##   Raccumb = col_double(),
##   ICV = col_double()
## )
```

```r
RandomArmColors = c( "#FFC200", "#007aa3")
```


```r
# remove participants that did not complete first and second scan (n=74)
# then add offlabel and dateDiff (in days columns)
# + a scan is by definition offlabel if it is the third scan
# then select the cols for analysis
df <- df %>%
  filter(first_complete == "Yes", 
         second_complete == "Yes",
         MR_exclusion == "No") %>%
  mutate(offLabel  = if_else(third_complete == "Yes", "Yes", ''),
         dateDiff = round(difftime(second_date, first_date, units = "days"), 0),
         STUDYID = parse_character(STUDYID),
         age = parse_number(age),
         category = factor(second_timepoint, levels = c("RCT","Relapse", "Off protocol"))) %>%
  select(STUDYID, randomization, sex, age, category, offLabel, dateDiff)
```

## cleaning the CT data



```r
# separating the subject id and anything afterwards to identify the longtudinal pipeline participants
# separating the subject id into site, "STUDYID" and timepoint columns
# filtering (two steps) to only include the longitudinal pipeline data
FS_long <- FS %>%
  separate(SubjID, into = c("subid", "longitudinal_pipe"), sep = '\\.', extra = "drop", fill = "right") %>%
  separate(subid, into = c("study", "site", "STUDYID", "timepoint"), fill = "right") %>%
  filter(longitudinal_pipe == "long") %>%
  filter(timepoint != "00", timepoint != "03", timepoint != "")

# adding columns that combine L and R
FS_long_plus <- FS_long %>%
  mutate(Thalamus = Lthal + Rthal,
         Hippocampus = Lhippo + Rhippo,
         Striatum = Lcaud + Rcaud + Lput + Rput)


# move CT from long to wide format
FS_wide <- FS_long_plus %>%
  gather(region, volume, -study, -site, -timepoint, -STUDYID, -longitudinal_pipe) %>%
  spread(timepoint, volume) %>%
  mutate(change = `02` - `01`,
         percchange = (`02`-`01`)/`01`) %>%
  gather(timepoint, volume, `01`, `02`, change, percchange) %>%
  unite(newcolnames, region, timepoint) %>%
  spread(newcolnames, volume)
```



```r
# merge CT values with df
ana_df <- inner_join(df, FS_wide, by='STUDYID') %>%
    mutate(STUDYID = as.character(STUDYID),
         dateDiff = as.numeric(dateDiff),
         RandomArm = factor(randomization, 
                       levels = c("O", "P"),
                       labels = c("Olanzapine", "Placebo"))) 

# write.csv
write_csv(ana_df, '../generated_csvs/STOPPD_participants_LandRVolumes_20181116.csv')
```

## report any mising values from clinical trial sample


```r
anti_join(df, FS_wide, by='STUDYID') %>%
  summarise(`Number of participants missing` = n()) %>%
  knitr::kable()
```


\begin{tabular}{r}
\hline
Number of participants missing\\
\hline
0\\
\hline
\end{tabular}

```r
ana_df %>%
  filter(is.na(ICV_01)) %>%
  summarise(`Number of participants missing timepoint 01` = n()) %>%
  knitr::kable()
```


\begin{tabular}{r}
\hline
Number of participants missing timepoint 01\\
\hline
0\\
\hline
\end{tabular}


```r
ana_df %>%
  filter(is.na(ICV_02)) %>%
  summarise(`Number of participants missing timepoint 02` = n()) %>%
  knitr::kable()
```


\begin{tabular}{r}
\hline
Number of participants missing timepoint 02\\
\hline
0\\
\hline
\end{tabular}


```r
library(tableone)

print(CreateTableOne(data = ana_df,
               strata = c("RandomArm"),
               vars = c("category", "Hippocampus_01", "Striatum_01", 'Thalamus_01')))
```

```
##                             Stratified by RandomArm
##                              Olanzapine         Placebo            p     
##   n                                38                 34                 
##   category (%)                                                      0.008
##      RCT                           26 (68.4)          14 (41.2)          
##      Relapse                        8 (21.1)          19 (55.9)          
##      Off protocol                   4 (10.5)           1 ( 2.9)          
##   Hippocampus_01 (mean (sd))  7538.32 (871.37)   7390.00 (1099.58)  0.526
##   Striatum_01 (mean (sd))    16931.60 (1825.84) 16610.18 (2077.11)  0.487
##   Thalamus_01 (mean (sd))    13326.38 (1834.30) 12989.71 (1916.03)  0.449
##                             Stratified by RandomArm
##                              test
##   n                              
##   category (%)                   
##      RCT                         
##      Relapse                     
##      Off protocol                
##   Hippocampus_01 (mean (sd))     
##   Striatum_01 (mean (sd))        
##   Thalamus_01 (mean (sd))
```

```r
print(CreateTableOne(data = ana_df,
               vars = c("category", "Hippocampus_01", "Striatum_01", 'Thalamus_01')))
```

```
##                             
##                              Overall           
##   n                                72          
##   category (%)                                 
##      RCT                           40 (55.6)   
##      Relapse                       27 (37.5)   
##      Off protocol                   5 ( 6.9)   
##   Hippocampus_01 (mean (sd))  7468.28 (981.44) 
##   Striatum_01 (mean (sd))    16779.82 (1941.31)
##   Thalamus_01 (mean (sd))    13167.39 (1867.72)
```

```r
ana_df %>%
  select(RandomArm, Hippocampus_01, Striatum_01, Thalamus_01) %>%
  gather(brain, mm, -RandomArm) %>%
  group_by(brain) %>%
  do(tidy(t.test(mm~RandomArm, data = .))) %>%
  select(brain, statistic, parameter, p.value, method) %>%
  rename(t = statistic, df = parameter) %>%
  knitr::kable(caption = "t.test for baseline group differences", digits = 2)
```

\begin{table}[t]

\caption{(\#tab:unnamed-chunk-5)t.test for baseline group differences}
\centering
\begin{tabular}{l|r|r|r|l}
\hline
brain & t & df & p.value & method\\
\hline
Hippocampus\_01 & 0.63 & 62.82 & 0.53 & Welch Two Sample t-test\\
\hline
Striatum\_01 & 0.69 & 66.19 & 0.49 & Welch Two Sample t-test\\
\hline
Thalamus\_01 & 0.76 & 68.33 & 0.45 & Welch Two Sample t-test\\
\hline
\end{tabular}
\end{table}

## creating an control error term calculating data frame


```r
## identify the repeat control in a column and mangle the STUDYID to match in a new column
FS_long1 <- FS_long_plus %>%
  mutate(repeat_run = if_else(str_sub(STUDYID,1,1)=="R", "02", "01"),
         STUDYID = str_replace(STUDYID, 'R',"")) 

## extra the repeat study ids as a character vector
repeat_ids <- filter(FS_long1, repeat_run == "02")$STUDYID

## filter for only the subjects who are in the repeats list then switch to wide format
FS_wide_controls <- FS_long1 %>%
  filter(STUDYID %in% repeat_ids) %>% 
  gather(region, volume, -study, -site, -timepoint, -STUDYID, -longitudinal_pipe, -repeat_run) %>%
  unite(newcolnames, region, repeat_run) %>%
  spread(newcolnames, volume)

#write.csv
  write.csv(FS_wide_controls, '../generated_csvs/STOPPD_errorControls_LandRVolumes_2018-11-05.csv', row.names = FALSE)
  
rm(FS_long1, repeat_ids)
```

## run RCT analysis (because it's simpler across volumes)


```r
# make sure that STUDYID is an character not a number
# make sure that dateDiff is a number, not an interger
# label the RandomArm variable  
RCT_SubCort <- ana_df %>%

  filter(category == "RCT")
```



```r
#boxplot of difference in thickness (y axis) by RandomArm group (x axis)
RCT_SubCort  %>%
  gather(region, volume_change, Thalamus_change, Hippocampus_change, Striatum_change) %>%
  mutate(Region = str_replace(region, '_change','')) %>%
ggplot(aes(x= RandomArm, y = volume_change, fill = RandomArm)) + 
     geom_boxplot(outlier.shape = NA, alpha = 0.0001) + 
     geom_dotplot(binaxis = 'y', stackdir = 'center') +
     geom_hline(yintercept = 0) +
     ggtitle("Freesurfer Subcortical Volume Changes") +
     xlab(NULL) +
     ylab("Change in Volume") +
     scale_fill_manual(values = RandomArmColors) +
     scale_shape_manual(values = c(21)) +
     facet_wrap(~Region) +
     theme_bw()
```

```
## `stat_bindot()` using `bins = 30`. Pick better value with `binwidth`.
```

![](10_STOPPD_freesurfer_subcortical_files/figure-latex/10-boxplot-ROIs_supplfig-1.pdf)<!-- --> 


```r
#boxplot of difference in thickness (y axis) by RandomArm group (x axis)
RCT_SubCort  %>%
  gather(region, volume_percchange, Thalamus_percchange, Hippocampus_percchange, Striatum_percchange) %>%
  mutate(Region = str_replace(region, '_percchange','')) %>%
ggplot(aes(x= RandomArm, y = volume_percchange, fill = RandomArm)) + 
     geom_boxplot(outlier.shape = NA, alpha = 0.0001) + 
     geom_dotplot(binaxis = 'y', stackdir = 'center') +
     geom_hline(yintercept = 0) +
     ggtitle("Freesurfer Subcortical Percent Volume Changes") +
     xlab(NULL) +
     ylab("Percent Change in Volume") +
     scale_fill_brewer(palette = "Dark2", direction = -1) +
     scale_shape_manual(values = c(21)) +
     facet_wrap(~Region) +
     theme_bw()
```

```
## `stat_bindot()` using `bins = 30`. Pick better value with `binwidth`.
```

![](10_STOPPD_freesurfer_subcortical_files/figure-latex/10-boxplot-ROIs-pchange-1.pdf)<!-- --> 

### Running RCT Linear Models

#### Thalamus


```r
#run linear model with covariates of sex and age
  fit_rct <- lmer(Thalamus_change ~ RandomArm + sex + age + (1|site), data= RCT_SubCort)
  summary(fit_rct)
```

```
## Linear mixed model fit by REML. t-tests use Satterthwaite's method [
## lmerModLmerTest]
## Formula: Thalamus_change ~ RandomArm + sex + age + (1 | site)
##    Data: RCT_SubCort
## 
## REML criterion at convergence: 504.2
## 
## Scaled residuals: 
##     Min      1Q  Median      3Q     Max 
## -2.1283 -0.6307  0.1721  0.5984  1.4686 
## 
## Random effects:
##  Groups   Name        Variance Std.Dev.
##  site     (Intercept)  4782     69.15  
##  Residual             41571    203.89  
## Number of obs: 40, groups:  site, 4
## 
## Fixed effects:
##                  Estimate Std. Error       df t value Pr(>|t|)  
## (Intercept)      -182.712    139.263   29.080  -1.312   0.1998  
## RandomArmPlacebo  156.222     69.401   34.664   2.251   0.0308 *
## sexM               16.157     66.777   34.519   0.242   0.8102  
## age                -1.784      2.470   35.941  -0.722   0.4748  
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Correlation of Fixed Effects:
##             (Intr) RndmAP sexM  
## RndmArmPlcb -0.011              
## sexM        -0.063  0.075       
## age         -0.892 -0.192 -0.183
```

```r
#run linear model with covariates of sex and age and site intercept
  fit_rct <- lmer(Thalamus_percchange ~ RandomArm + sex + age + (1|site), data= RCT_SubCort)
  summary(fit_rct)
```

```
## Linear mixed model fit by REML. t-tests use Satterthwaite's method [
## lmerModLmerTest]
## Formula: Thalamus_percchange ~ RandomArm + sex + age + (1 | site)
##    Data: RCT_SubCort
## 
## REML criterion at convergence: -182.1
## 
## Scaled residuals: 
##     Min      1Q  Median      3Q     Max 
## -1.9666 -0.6349  0.1623  0.6530  1.5828 
## 
## Random effects:
##  Groups   Name        Variance  Std.Dev.
##  site     (Intercept) 2.737e-05 0.005232
##  Residual             2.181e-04 0.014767
## Number of obs: 40, groups:  site, 4
## 
## Fixed effects:
##                    Estimate Std. Error         df t value Pr(>|t|)  
## (Intercept)      -0.0071511  0.0101398 29.1678166  -0.705   0.4862  
## RandomArmPlacebo  0.0116258  0.0050286 34.6583454   2.312   0.0268 *
## sexM              0.0035373  0.0048381 34.5125265   0.731   0.4696  
## age              -0.0002695  0.0001793 35.9629119  -1.503   0.1415  
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Correlation of Fixed Effects:
##             (Intr) RndmAP sexM  
## RndmArmPlcb -0.010              
## sexM        -0.064  0.075       
## age         -0.890 -0.192 -0.182
```

#### Striatum


```r
#run linear model with covariates of sex and age
  fit_rct <- lmer(Striatum_change ~ RandomArm + sex + age + (1|site), data= RCT_SubCort)
  summary(fit_rct)
```

```
## Linear mixed model fit by REML. t-tests use Satterthwaite's method [
## lmerModLmerTest]
## Formula: Striatum_change ~ RandomArm + sex + age + (1 | site)
##    Data: RCT_SubCort
## 
## REML criterion at convergence: 526.8
## 
## Scaled residuals: 
##      Min       1Q   Median       3Q      Max 
## -2.68942 -0.49396  0.03391  0.53963  1.92558 
## 
## Random effects:
##  Groups   Name        Variance Std.Dev.
##  site     (Intercept)     0      0.0   
##  Residual             82431    287.1   
## Number of obs: 40, groups:  site, 4
## 
## Fixed effects:
##                  Estimate Std. Error       df t value Pr(>|t|)
## (Intercept)      -169.875    180.214   36.000  -0.943    0.352
## RandomArmPlacebo  -80.518     96.821   36.000  -0.832    0.411
## sexM              -12.233     93.203   36.000  -0.131    0.896
## age                -2.318      3.343   36.000  -0.693    0.492
## 
## Correlation of Fixed Effects:
##             (Intr) RndmAP sexM  
## RndmArmPlcb -0.022              
## sexM        -0.044  0.067       
## age         -0.921 -0.181 -0.201
```

```r
#run linear model with covariates of sex and age
  fit_rct <- lmer(Striatum_percchange ~ RandomArm + sex + age + (1|site), data= RCT_SubCort)
  summary(fit_rct)
```

```
## Linear mixed model fit by REML. t-tests use Satterthwaite's method [
## lmerModLmerTest]
## Formula: Striatum_percchange ~ RandomArm + sex + age + (1 | site)
##    Data: RCT_SubCort
## 
## REML criterion at convergence: -172.9
## 
## Scaled residuals: 
##      Min       1Q   Median       3Q      Max 
## -2.55786 -0.37884  0.03767  0.58473  2.04779 
## 
## Random effects:
##  Groups   Name        Variance  Std.Dev.
##  site     (Intercept) 0.0000000 0.00000 
##  Residual             0.0002991 0.01729 
## Number of obs: 40, groups:  site, 4
## 
## Fixed effects:
##                    Estimate Std. Error         df t value Pr(>|t|)
## (Intercept)      -0.0066369  0.0108551 36.0000000  -0.611    0.545
## RandomArmPlacebo -0.0048423  0.0058319 36.0000000  -0.830    0.412
## sexM              0.0018895  0.0056140 36.0000000   0.337    0.738
## age              -0.0002271  0.0002013 36.0000000  -1.128    0.267
## 
## Correlation of Fixed Effects:
##             (Intr) RndmAP sexM  
## RndmArmPlcb -0.022              
## sexM        -0.044  0.067       
## age         -0.921 -0.181 -0.201
```

#### Hippocampus


```r
#run linear model with covariates of sex and age
  fit_rct <- lmer(Hippocampus_change ~ RandomArm + sex + age + (1|site), data= RCT_SubCort)
  summary(fit_rct)
```

```
## Linear mixed model fit by REML. t-tests use Satterthwaite's method [
## lmerModLmerTest]
## Formula: Hippocampus_change ~ RandomArm + sex + age + (1 | site)
##    Data: RCT_SubCort
## 
## REML criterion at convergence: 476.5
## 
## Scaled residuals: 
##     Min      1Q  Median      3Q     Max 
## -3.4916 -0.6315 -0.0032  0.5931  1.8058 
## 
## Random effects:
##  Groups   Name        Variance Std.Dev.
##  site     (Intercept)     0      0.0   
##  Residual             20379    142.8   
## Number of obs: 40, groups:  site, 4
## 
## Fixed effects:
##                  Estimate Std. Error     df t value Pr(>|t|)  
## (Intercept)         5.198     89.606 36.000   0.058   0.9541  
## RandomArmPlacebo   97.320     48.141 36.000   2.022   0.0507 .
## sexM               -6.910     46.342 36.000  -0.149   0.8823  
## age                -2.171      1.662 36.000  -1.306   0.1998  
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Correlation of Fixed Effects:
##             (Intr) RndmAP sexM  
## RndmArmPlcb -0.022              
## sexM        -0.044  0.067       
## age         -0.921 -0.181 -0.201
```


```r
#run linear model with covariates of sex and age
  fit_rct <- lm(Hippocampus_percchange ~ RandomArm + sex + age, data= RCT_SubCort)
  summary(fit_rct)
```

```
## 
## Call:
## lm(formula = Hippocampus_percchange ~ RandomArm + sex + age, 
##     data = RCT_SubCort)
## 
## Residuals:
##      Min       1Q   Median       3Q      Max 
## -0.06398 -0.01338  0.00013  0.01082  0.03488 
## 
## Coefficients:
##                    Estimate Std. Error t value Pr(>|t|)  
## (Intercept)       0.0039836  0.0119300   0.334   0.7404  
## RandomArmPlacebo  0.0124023  0.0064094   1.935   0.0609 .
## sexM             -0.0001382  0.0061699  -0.022   0.9823  
## age              -0.0003643  0.0002213  -1.646   0.1084  
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Residual standard error: 0.01901 on 36 degrees of freedom
## Multiple R-squared:  0.1357,	Adjusted R-squared:  0.0637 
## F-statistic: 1.884 on 3 and 36 DF,  p-value: 0.1497
```

```r
#run linear model with covariates of sex and age
  fit_rct <- lm(Hippocampus_percchange ~ RandomArm + sex + age + site, data= RCT_SubCort)
  summary(fit_rct)
```

```
## 
## Call:
## lm(formula = Hippocampus_percchange ~ RandomArm + sex + age + 
##     site, data = RCT_SubCort)
## 
## Residuals:
##       Min        1Q    Median        3Q       Max 
## -0.062098 -0.012941  0.001371  0.014436  0.035558 
## 
## Coefficients:
##                    Estimate Std. Error t value Pr(>|t|)  
## (Intercept)       0.0073077  0.0137456   0.532   0.5985  
## RandomArmPlacebo  0.0123147  0.0067210   1.832   0.0760 .
## sexM             -0.0009176  0.0064477  -0.142   0.8877  
## age              -0.0004477  0.0002510  -1.784   0.0836 .
## siteMAS           0.0040824  0.0085714   0.476   0.6370  
## siteNKI          -0.0020478  0.0081914  -0.250   0.8041  
## sitePMC           0.0080461  0.0100633   0.800   0.4297  
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Residual standard error: 0.01955 on 33 degrees of freedom
## Multiple R-squared:  0.1618,	Adjusted R-squared:  0.009442 
## F-statistic: 1.062 on 6 and 33 DF,  p-value: 0.4047
```
----

## RCT & Relapse (with time as factor)

### Thalamus



```r
#restructure data for RCT & Relapse participants (N=72)
  RCTRelapse_Thalamus <- ana_df %>%
    gather(oldcolname, volume, Thalamus_01, Thalamus_02) %>%
    mutate(model_days = if_else(oldcolname == "Thalamus_01", 1, dateDiff))

RCTRelapse_Thalamus %>% filter(model_days == 1) %>% count(RandomArm) %>% knitr::kable() 
```


\begin{tabular}{l|r}
\hline
RandomArm & n\\
\hline
Olanzapine & 38\\
\hline
Placebo & 34\\
\hline
\end{tabular}

```r
RCTRelapse_Thalamus_sense <- RCTRelapse_Thalamus %>% filter(category != "Off protocol")
```


```r
 RCTRelapse_Thalamus %>%
  mutate(roi = "Thalamus") %>%
  ggplot(aes(x=model_days, y=volume, fill = RandomArm)) + 
  geom_point(aes(shape = category)) + 
  geom_line(aes(group=STUDYID, color = RandomArm), alpha = 0.5) + 
  geom_smooth(aes(color = RandomArm), method="lm") +
  labs(x = "Days between MRIs", y = "Volume", colour = NULL) +
  scale_colour_manual(values = RandomArmColors) +
  scale_fill_manual(values = RandomArmColors) +
  scale_shape_manual(values = c(21:23)) +
  theme_bw()  +
  facet_wrap(~roi)
```

![](10_STOPPD_freesurfer_subcortical_files/figure-latex/RCTRelapse_Thalamus_plot_suppF-1.pdf)<!-- --> 


```r
#run mixed linear model, with covariates
  fit_all <- lmer(volume ~ RandomArm*model_days + sex + age + (1|site) + (1|STUDYID), data= RCTRelapse_Thalamus)
  summary(fit_all)
```

```
## Linear mixed model fit by REML. t-tests use Satterthwaite's method [
## lmerModLmerTest]
## Formula: volume ~ RandomArm * model_days + sex + age + (1 | site) + (1 |  
##     STUDYID)
##    Data: RCTRelapse_Thalamus
## 
## REML criterion at convergence: 2189.4
## 
## Scaled residuals: 
##     Min      1Q  Median      3Q     Max 
## -3.5900 -0.4277  0.0075  0.3954  3.5698 
## 
## Random effects:
##  Groups   Name        Variance Std.Dev.
##  STUDYID  (Intercept) 1725112  1313.4  
##  site     (Intercept)  329065   573.6  
##  Residual               29304   171.2  
## Number of obs: 144, groups:  STUDYID, 72; site, 4
## 
## Fixed effects:
##                               Estimate Std. Error         df t value
## (Intercept)                 15820.7399   678.5347    34.2363  23.316
## RandomArmPlacebo             -154.1418   316.5539    66.4583  -0.487
## model_days                     -1.1867     0.1804    70.1254  -6.580
## sexM                         1860.4791   316.4113    65.6780   5.880
## age                           -61.0926    10.2816    65.3437  -5.942
## RandomArmPlacebo:model_days     0.8101     0.3046    70.3093   2.659
##                             Pr(>|t|)    
## (Intercept)                  < 2e-16 ***
## RandomArmPlacebo             0.62791    
## model_days                  7.18e-09 ***
## sexM                        1.51e-07 ***
## age                         1.20e-07 ***
## RandomArmPlacebo:model_days  0.00969 ** 
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Correlation of Fixed Effects:
##             (Intr) RndmAP mdl_dy sexM   age   
## RndmArmPlcb -0.187                            
## model_days  -0.031  0.059                     
## sexM        -0.171  0.052  0.001              
## age         -0.810 -0.064  0.004 -0.078       
## RndmArmPl:_  0.017 -0.078 -0.592  0.001 -0.001
```


```r
#run mixed linear model, with covariates
  fit_all <- lmer(volume ~ RandomArm*model_days + sex + age + (1|site) + (1|STUDYID), data= RCTRelapse_Thalamus_sense)
  summary(fit_all)  
```

```
## Linear mixed model fit by REML. t-tests use Satterthwaite's method [
## lmerModLmerTest]
## Formula: volume ~ RandomArm * model_days + sex + age + (1 | site) + (1 |  
##     STUDYID)
##    Data: RCTRelapse_Thalamus_sense
## 
## REML criterion at convergence: 2039.4
## 
## Scaled residuals: 
##     Min      1Q  Median      3Q     Max 
## -3.4937 -0.4158  0.0095  0.3908  3.4791 
## 
## Random effects:
##  Groups   Name        Variance Std.Dev.
##  STUDYID  (Intercept) 1786865  1336.7  
##  site     (Intercept)  334546   578.4  
##  Residual               30876   175.7  
## Number of obs: 134, groups:  STUDYID, 67; site, 4
## 
## Fixed effects:
##                               Estimate Std. Error         df t value
## (Intercept)                 15947.4510   700.6484    35.0741  22.761
## RandomArmPlacebo             -196.6842   332.4836    61.3571  -0.592
## model_days                     -1.1802     0.1920    65.1012  -6.146
## sexM                         1979.3127   335.5814    60.6584   5.898
## age                           -63.4241    10.8760    60.3263  -5.832
## RandomArmPlacebo:model_days     0.8061     0.3168    65.2702   2.544
##                             Pr(>|t|)    
## (Intercept)                  < 2e-16 ***
## RandomArmPlacebo              0.5563    
## model_days                  5.40e-08 ***
## sexM                        1.76e-07 ***
## age                         2.31e-07 ***
## RandomArmPlacebo:model_days   0.0133 *  
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Correlation of Fixed Effects:
##             (Intr) RndmAP mdl_dy sexM   age   
## RndmArmPlcb -0.191                            
## model_days  -0.033  0.061                     
## sexM        -0.119  0.008 -0.001              
## age         -0.815 -0.059  0.005 -0.125       
## RndmArmPl:_  0.020 -0.080 -0.606  0.004 -0.004
```

### Striatum



```r
#restructure data for RCT & Relapse participants (N=72)
  RCTRelapse_Striatum <- ana_df %>%
    gather(oldcolname, volume, Striatum_01, Striatum_02) %>%
    mutate(model_days = if_else(oldcolname == "Striatum_01", 1, dateDiff)) %>%
  mutate(age_centered = age - mean(age),
         model_days_centered = model_days - mean(model_days))

RCTRelapse_Striatum %>% filter(model_days == 1) %>% count(RandomArm) %>% knitr::kable() 
```


\begin{tabular}{l|r}
\hline
RandomArm & n\\
\hline
Olanzapine & 38\\
\hline
Placebo & 34\\
\hline
\end{tabular}

```r
RCTRelapse_Striatum_sense <- RCTRelapse_Striatum %>% filter(category != "Off protocol")
```


```r
 RCTRelapse_Striatum %>%
  mutate(roi = "Striatum") %>%
  ggplot(aes(x=model_days, y=volume, fill = RandomArm)) + 
  geom_point(aes(shape = category)) + 
  geom_line(aes(group=STUDYID, color = RandomArm), alpha = 0.5) + 
  geom_smooth(aes(color = RandomArm), method="lm") +
  labs(x = "Days between MRIs", y = "Volume", colour = NULL) +
  scale_colour_manual(values = RandomArmColors) +
  scale_fill_manual(values = RandomArmColors) +
  scale_shape_manual(values = c(21:23)) +
  theme_bw()  +
  facet_wrap(~roi)
```

![](10_STOPPD_freesurfer_subcortical_files/figure-latex/RCTRelapse_Striatum_plot_supplE-1.pdf)<!-- --> 


```r
  fit_all <- lmer(volume ~ RandomArm*model_days + age + sex + (1|site) + (1|STUDYID), data= RCTRelapse_Striatum)
  print(fit_all)
```

```
## Linear mixed model fit by REML ['lmerModLmerTest']
## Formula: volume ~ RandomArm * model_days + age + sex + (1 | site) + (1 |  
##     STUDYID)
##    Data: RCTRelapse_Striatum
## REML criterion at convergence: 2250.494
## Random effects:
##  Groups   Name        Std.Dev.
##  STUDYID  (Intercept) 1673.98 
##  site     (Intercept)   50.84 
##  Residual              215.57 
## Number of obs: 144, groups:  STUDYID, 72; site, 4
## Fixed Effects:
##                 (Intercept)             RandomArmPlacebo  
##                  18785.0451                    -197.5670  
##                  model_days                          age  
##                     -1.1427                     -47.6647  
##                        sexM  RandomArmPlacebo:model_days  
##                   1604.6240                      -0.1886
```

```r
  summary(fit_all)
```

```
## Linear mixed model fit by REML. t-tests use Satterthwaite's method [
## lmerModLmerTest]
## Formula: volume ~ RandomArm * model_days + age + sex + (1 | site) + (1 |  
##     STUDYID)
##    Data: RCTRelapse_Striatum
## 
## REML criterion at convergence: 2250.5
## 
## Scaled residuals: 
##      Min       1Q   Median       3Q      Max 
## -2.37598 -0.36232  0.03756  0.33443  2.82794 
## 
## Random effects:
##  Groups   Name        Variance Std.Dev.
##  STUDYID  (Intercept) 2802195  1673.98 
##  site     (Intercept)    2585    50.84 
##  Residual               46470   215.57 
## Number of obs: 144, groups:  STUDYID, 72; site, 4
## 
## Fixed effects:
##                               Estimate Std. Error         df t value
## (Intercept)                 18785.0451   774.3552    64.4702  24.259
## RandomArmPlacebo             -197.5670   398.9866    68.8067  -0.495
## model_days                     -1.1427     0.2271    70.1337  -5.031
## age                           -47.6647    13.0310    67.1314  -3.658
## sexM                         1604.6240   399.1139    67.7920   4.020
## RandomArmPlacebo:model_days    -0.1886     0.3836    70.3100  -0.492
##                             Pr(>|t|)    
## (Intercept)                  < 2e-16 ***
## RandomArmPlacebo            0.622056    
## model_days                  3.61e-06 ***
## age                         0.000501 ***
## sexM                        0.000148 ***
## RandomArmPlacebo:model_days 0.624398    
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Correlation of Fixed Effects:
##             (Intr) RndmAP mdl_dy age    sexM  
## RndmArmPlcb -0.202                            
## model_days  -0.034  0.059                     
## age         -0.902 -0.055  0.003              
## sexM        -0.172  0.037  0.000 -0.079       
## RndmArmPl:_  0.019 -0.078 -0.592 -0.001  0.002
```


```r
    fit_all <- lmer(volume ~ RandomArm*model_days + age + sex + (1|site) + (1|STUDYID), data= RCTRelapse_Striatum_sense)
  print(fit_all)
```

```
## Linear mixed model fit by REML ['lmerModLmerTest']
## Formula: volume ~ RandomArm * model_days + age + sex + (1 | site) + (1 |  
##     STUDYID)
##    Data: RCTRelapse_Striatum_sense
## REML criterion at convergence: 2069.416
## Random effects:
##  Groups   Name        Std.Dev. 
##  STUDYID  (Intercept) 1.534e+03
##  site     (Intercept) 5.355e-04
##  Residual             1.999e+02
## Number of obs: 134, groups:  STUDYID, 67; site, 4
## Fixed Effects:
##                 (Intercept)             RandomArmPlacebo  
##                   1.921e+04                   -1.834e+01  
##                  model_days                          age  
##                  -1.245e+00                   -6.041e+01  
##                        sexM  RandomArmPlacebo:model_days  
##                   1.809e+03                   -8.351e-02
```

```r
  summary(fit_all)
```

```
## Linear mixed model fit by REML. t-tests use Satterthwaite's method [
## lmerModLmerTest]
## Formula: volume ~ RandomArm * model_days + age + sex + (1 | site) + (1 |  
##     STUDYID)
##    Data: RCTRelapse_Striatum_sense
## 
## REML criterion at convergence: 2069.4
## 
## Scaled residuals: 
##      Min       1Q   Median       3Q      Max 
## -1.87831 -0.38523  0.02011  0.35407  1.90856 
## 
## Random effects:
##  Groups   Name        Variance  Std.Dev. 
##  STUDYID  (Intercept) 2.353e+06 1.534e+03
##  site     (Intercept) 2.868e-07 5.355e-04
##  Residual             3.996e+04 1.999e+02
## Number of obs: 134, groups:  STUDYID, 67; site, 4
## 
## Fixed effects:
##                               Estimate Std. Error         df t value
## (Intercept)                  1.921e+04  7.235e+02  6.316e+01  26.558
## RandomArmPlacebo            -1.834e+01  3.783e+02  6.384e+01  -0.048
## model_days                  -1.245e+00  2.184e-01  6.511e+01  -5.700
## age                         -6.041e+01  1.243e+01  6.300e+01  -4.862
## sexM                         1.809e+03  3.815e+02  6.300e+01   4.742
## RandomArmPlacebo:model_days -8.351e-02  3.604e-01  6.527e+01  -0.232
##                             Pr(>|t|)    
## (Intercept)                  < 2e-16 ***
## RandomArmPlacebo               0.961    
## model_days                  3.14e-07 ***
## age                         8.09e-06 ***
## sexM                        1.25e-05 ***
## RandomArmPlacebo:model_days    0.817    
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Correlation of Fixed Effects:
##             (Intr) RndmAP mdl_dy age    sexM  
## RndmArmPlcb -0.205                            
## model_days  -0.036  0.061                     
## age         -0.901 -0.054  0.005              
## sexM        -0.115 -0.007 -0.001 -0.126       
## RndmArmPl:_  0.022 -0.080 -0.606 -0.003  0.004
```

### Hippocampus



```r
#restructure data for RCT & Relapse participants (N=72)
  RCTRelapse_Hippocampus <- ana_df %>%
    gather(oldcolname, volume, Hippocampus_01, Hippocampus_02) %>%
    mutate(model_days = if_else(oldcolname == "Hippocampus_01", 1, dateDiff)) %>%
  mutate(age_centered = age - mean(age),
         model_days_centered = model_days - mean(model_days))

RCTRelapse_Hippocampus %>% filter(model_days == 1) %>% count(RandomArm) %>% knitr::kable() 
```


\begin{tabular}{l|r}
\hline
RandomArm & n\\
\hline
Olanzapine & 38\\
\hline
Placebo & 34\\
\hline
\end{tabular}

```r
RCTRelapse_Hippocampus_sense <- RCTRelapse_Hippocampus %>% filter(category != "Off protocol")
```


```r
#plot all data, including outlier (participant 210030)
 RCTRelapse_Hippocampus %>%
  mutate(roi = "Hippocampus") %>%
  ggplot(aes(x=model_days, y=volume, fill = RandomArm)) + 
  geom_point(aes(shape = category)) + 
  geom_line(aes(group=STUDYID, color = RandomArm), alpha = 0.5) + 
  geom_smooth(aes(color = RandomArm), method="lm") +
  labs(x = "Days between MRIs", y = "Volume", colour = NULL) +
  scale_colour_manual(values = RandomArmColors) +
  scale_fill_manual(values = RandomArmColors) +
  scale_shape_manual(values = c(21:23)) +
  theme_bw()  +
  facet_wrap(~roi)
```

![](10_STOPPD_freesurfer_subcortical_files/figure-latex/RCTRelapse_Hippocampus_plot_supplD-1.pdf)<!-- --> 


```r
#run mixed linear model, with covariates
  fit_all <- lmer(volume ~ RandomArm*model_days_centered + sex + age_centered + (1|site) + (1|STUDYID), data= RCTRelapse_Hippocampus)
  print(fit_all)
```

```
## Linear mixed model fit by REML ['lmerModLmerTest']
## Formula: volume ~ RandomArm * model_days_centered + sex + age_centered +  
##     (1 | site) + (1 | STUDYID)
##    Data: RCTRelapse_Hippocampus
## REML criterion at convergence: 2049.592
## Random effects:
##  Groups   Name        Std.Dev.
##  STUDYID  (Intercept) 829.3   
##  site     (Intercept) 164.5   
##  Residual             100.6   
## Number of obs: 144, groups:  STUDYID, 72; site, 4
## Fixed Effects:
##                          (Intercept)  
##                            7238.7852  
##                     RandomArmPlacebo  
##                             -72.6709  
##                  model_days_centered  
##                              -0.4047  
##                                 sexM  
##                             584.5259  
##                         age_centered  
##                             -31.7712  
## RandomArmPlacebo:model_days_centered  
##                               0.2636
```

```r
  summary(fit_all)
```

```
## Linear mixed model fit by REML. t-tests use Satterthwaite's method [
## lmerModLmerTest]
## Formula: volume ~ RandomArm * model_days_centered + sex + age_centered +  
##     (1 | site) + (1 | STUDYID)
##    Data: RCTRelapse_Hippocampus
## 
## REML criterion at convergence: 2049.6
## 
## Scaled residuals: 
##      Min       1Q   Median       3Q      Max 
## -2.40824 -0.41230 -0.00417  0.40063  2.36527 
## 
## Random effects:
##  Groups   Name        Variance Std.Dev.
##  STUDYID  (Intercept) 687763   829.3   
##  site     (Intercept)  27053   164.5   
##  Residual              10111   100.6   
## Number of obs: 144, groups:  STUDYID, 72; site, 4
## 
## Fixed effects:
##                                       Estimate Std. Error        df
## (Intercept)                          7238.7852   188.6898    7.7620
## RandomArmPlacebo                      -72.6709   198.2076   66.6777
## model_days_centered                    -0.4047     0.1060   70.1144
## sexM                                  584.5259   198.7201   66.5371
## age_centered                          -31.7712     6.4714   65.6248
## RandomArmPlacebo:model_days_centered    0.2636     0.1790   70.2723
##                                      t value Pr(>|t|)    
## (Intercept)                           38.363 3.93e-10 ***
## RandomArmPlacebo                      -0.367 0.715048    
## model_days_centered                   -3.820 0.000286 ***
## sexM                                   2.941 0.004491 ** 
## age_centered                          -4.909 6.36e-06 ***
## RandomArmPlacebo:model_days_centered   1.473 0.145272    
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Correlation of Fixed Effects:
##             (Intr) RndmAP mdl_d_ sexM   ag_cnt
## RndmArmPlcb -0.529                            
## mdl_dys_cnt -0.009  0.008                     
## sexM        -0.518  0.046  0.000              
## age_centerd  0.069 -0.060  0.003 -0.079       
## RndmArmP:__  0.005  0.005 -0.592  0.001 -0.001
```

```r
#run mixed linear model, with covariates
  fit_all <- lmer(volume ~ RandomArm*model_days_centered*age_centered + sex + (1|site) + (1|STUDYID), data= RCTRelapse_Hippocampus)
  print(fit_all)
```

```
## Linear mixed model fit by REML ['lmerModLmerTest']
## Formula: volume ~ RandomArm * model_days_centered * age_centered + sex +  
##     (1 | site) + (1 | STUDYID)
##    Data: RCTRelapse_Hippocampus
## REML criterion at convergence: 2055.319
## Random effects:
##  Groups   Name        Std.Dev.
##  STUDYID  (Intercept) 831.9   
##  site     (Intercept) 189.1   
##  Residual             100.5   
## Number of obs: 144, groups:  STUDYID, 72; site, 4
## Fixed Effects:
##                                       (Intercept)  
##                                         7.239e+03  
##                                  RandomArmPlacebo  
##                                        -7.478e+01  
##                               model_days_centered  
##                                        -4.255e-01  
##                                      age_centered  
##                                        -3.460e+01  
##                                              sexM  
##                                         5.897e+02  
##              RandomArmPlacebo:model_days_centered  
##                                         2.785e-01  
##                     RandomArmPlacebo:age_centered  
##                                         6.427e+00  
##                  model_days_centered:age_centered  
##                                        -8.378e-03  
## RandomArmPlacebo:model_days_centered:age_centered  
##                                         1.802e-02  
## fit warnings:
## Some predictor variables are on very different scales: consider rescaling
```

```r
  summary(fit_all)  
```

```
## Linear mixed model fit by REML. t-tests use Satterthwaite's method [
## lmerModLmerTest]
## Formula: volume ~ RandomArm * model_days_centered * age_centered + sex +  
##     (1 | site) + (1 | STUDYID)
##    Data: RCTRelapse_Hippocampus
## 
## REML criterion at convergence: 2055.3
## 
## Scaled residuals: 
##      Min       1Q   Median       3Q      Max 
## -2.35631 -0.44217  0.01185  0.39228  2.31761 
## 
## Random effects:
##  Groups   Name        Variance Std.Dev.
##  STUDYID  (Intercept) 692108   831.9   
##  site     (Intercept)  35772   189.1   
##  Residual              10110   100.5   
## Number of obs: 144, groups:  STUDYID, 72; site, 4
## 
## Fixed effects:
##                                                     Estimate Std. Error
## (Intercept)                                        7.239e+03  1.958e+02
## RandomArmPlacebo                                  -7.478e+01  1.990e+02
## model_days_centered                               -4.255e-01  1.073e-01
## age_centered                                      -3.460e+01  9.044e+00
## sexM                                               5.897e+02  2.000e+02
## RandomArmPlacebo:model_days_centered               2.785e-01  1.799e-01
## RandomArmPlacebo:age_centered                      6.427e+00  1.328e+01
## model_days_centered:age_centered                  -8.378e-03  6.862e-03
## RandomArmPlacebo:model_days_centered:age_centered  1.802e-02  1.470e-02
##                                                           df t value
## (Intercept)                                        7.246e+00  36.977
## RandomArmPlacebo                                   6.544e+01  -0.376
## model_days_centered                                6.812e+01  -3.966
## age_centered                                       6.636e+01  -3.825
## sexM                                               6.530e+01   2.948
## RandomArmPlacebo:model_days_centered               6.827e+01   1.548
## RandomArmPlacebo:age_centered                      6.657e+01   0.484
## model_days_centered:age_centered                   6.810e+01  -1.221
## RandomArmPlacebo:model_days_centered:age_centered  6.836e+01   1.226
##                                                   Pr(>|t|)    
## (Intercept)                                        1.6e-09 ***
## RandomArmPlacebo                                  0.708339    
## model_days_centered                               0.000178 ***
## age_centered                                      0.000291 ***
## sexM                                              0.004430 ** 
## RandomArmPlacebo:model_days_centered              0.126208    
## RandomArmPlacebo:age_centered                     0.629926    
## model_days_centered:age_centered                  0.226334    
## RandomArmPlacebo:model_days_centered:age_centered 0.224457    
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Correlation of Fixed Effects:
##             (Intr) RndmAP mdl_d_ ag_cnt sexM   RnAP:__ RnAP:_ md__:_
## RndmArmPlcb -0.514                                                  
## mdl_dys_cnt -0.008  0.008                                           
## age_centerd  0.082 -0.036  0.002                                    
## sexM        -0.507  0.046 -0.001 -0.104                             
## RndmArmP:__  0.004  0.006 -0.596 -0.002  0.002                      
## RndmArmPl:_ -0.049 -0.011 -0.001 -0.695  0.069  0.000               
## mdl_dys_c:_  0.005 -0.002  0.157 -0.012 -0.007 -0.094   0.008       
## RndmAP:__:_ -0.004  0.001 -0.074  0.006  0.008  0.005   0.031 -0.467
## fit warnings:
## Some predictor variables are on very different scales: consider rescaling
```

```r
#run mixed linear model, with covariates
  fit_all <- lmer(volume ~ RandomArm*model_days_centered + sex + age_centered + (1|site) + (1|STUDYID), data= RCTRelapse_Hippocampus_sense)
  print(fit_all)
```

```
## Linear mixed model fit by REML ['lmerModLmerTest']
## Formula: volume ~ RandomArm * model_days_centered + sex + age_centered +  
##     (1 | site) + (1 | STUDYID)
##    Data: RCTRelapse_Hippocampus_sense
## REML criterion at convergence: 1903.747
## Random effects:
##  Groups   Name        Std.Dev.
##  STUDYID  (Intercept) 846.51  
##  site     (Intercept)  92.00  
##  Residual              99.13  
## Number of obs: 134, groups:  STUDYID, 67; site, 4
## Fixed Effects:
##                          (Intercept)  
##                            7254.3186  
##                     RandomArmPlacebo  
##                            -110.9660  
##                  model_days_centered  
##                              -0.4211  
##                                 sexM  
##                             583.4918  
##                         age_centered  
##                             -31.2064  
## RandomArmPlacebo:model_days_centered  
##                               0.2795
```

```r
  summary(fit_all)
```

```
## Linear mixed model fit by REML. t-tests use Satterthwaite's method [
## lmerModLmerTest]
## Formula: volume ~ RandomArm * model_days_centered + sex + age_centered +  
##     (1 | site) + (1 | STUDYID)
##    Data: RCTRelapse_Hippocampus_sense
## 
## REML criterion at convergence: 1903.7
## 
## Scaled residuals: 
##     Min      1Q  Median      3Q     Max 
## -2.4428 -0.4317 -0.0118  0.4159  2.3980 
## 
## Random effects:
##  Groups   Name        Variance Std.Dev.
##  STUDYID  (Intercept) 716581   846.51  
##  site     (Intercept)   8465    92.00  
##  Residual               9827    99.13  
## Number of obs: 134, groups:  STUDYID, 67; site, 4
## 
## Fixed effects:
##                                       Estimate Std. Error        df
## (Intercept)                          7254.3186   181.3540    7.4248
## RandomArmPlacebo                     -110.9660   208.3516   61.6101
## model_days_centered                    -0.4211     0.1083   65.0858
## sexM                                  583.4918   210.8283   61.6355
## age_centered                          -31.2064     6.8588   60.0991
## RandomArmPlacebo:model_days_centered    0.2795     0.1788   65.2177
##                                      t value Pr(>|t|)    
## (Intercept)                           40.001 6.06e-10 ***
## RandomArmPlacebo                      -0.533  0.59623    
## model_days_centered                   -3.887  0.00024 ***
## sexM                                   2.768  0.00745 ** 
## age_centered                          -4.550 2.66e-05 ***
## RandomArmPlacebo:model_days_centered   1.563  0.12278    
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Correlation of Fixed Effects:
##             (Intr) RndmAP mdl_d_ sexM   ag_cnt
## RndmArmPlcb -0.569                            
## mdl_dys_cnt -0.011  0.010                     
## sexM        -0.526 -0.003 -0.001              
## age_centerd  0.116 -0.056  0.004 -0.126       
## RndmArmP:__  0.005  0.003 -0.606  0.003 -0.003
```

