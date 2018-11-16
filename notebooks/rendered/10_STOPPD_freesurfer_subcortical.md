# Freesurfer Derived Subcortical Volumes


```r
library(tidyverse)
```

```
## -- Attaching packages --------------------------------------------------------------------------------------------- tidyverse 1.2.1 --
```

```
## √ ggplot2 3.1.0     √ purrr   0.2.4
## √ tibble  1.4.1     √ dplyr   0.7.7
## √ tidyr   0.7.2     √ stringr 1.2.0
## √ readr   1.1.1     √ forcats 0.2.0
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
         age = parse_number(age)) %>%
  rename(category = "second_timepoint") %>%
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
  mutate(change = `02` - `01`) %>%
  gather(timepoint, volume, `01`, `02`, change) %>%
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
# label the randomization variable  
RCT_SubCort <- ana_df %>%

  filter(category == "RCT")
```



```r
#boxplot of difference in thickness (y axis) by randomization group (x axis)
RCT_SubCort  %>%
  gather(region, volume_change, Thalamus_change, Hippocampus_change, Striatum_change) %>%
  mutate(Region = str_replace(region, '_change','')) %>%
ggplot(aes(x= RandomArm, y = volume_change)) + 
     geom_boxplot(outlier.shape = NA) + 
     geom_dotplot(binaxis = 'y', stackdir = 'center') +
     geom_hline(yintercept = 0) +
     ggtitle("Freesurfer Subcortical Volume Changes") +
     xlab(NULL) +
     ylab("Change in Volume") +
     facet_wrap(~Region) +
     theme_bw()
```

```
## `stat_bindot()` using `bins = 30`. Pick better value with `binwidth`.
```

![](10_STOPPD_freesurfer_subcortical_files/figure-latex/10-boxplot-ROIs-1.pdf)<!-- --> 
### Running RCT Linear Models

#### Thalamus


```r
#run linear model without covariates
  fit_rct <- lm(Thalamus_change ~ RandomArm, data= RCT_SubCort)
  summary(fit_rct)
```

```
## 
## Call:
## lm(formula = Thalamus_change ~ RandomArm, data = RCT_SubCort)
## 
## Residuals:
##     Min      1Q  Median      3Q     Max 
## -454.67 -105.09    3.09  155.71  366.04 
## 
## Coefficients:
##                  Estimate Std. Error t value Pr(>|t|)    
## (Intercept)       -286.13      41.25  -6.937    3e-08 ***
## RandomArmPlacebo   157.29      69.73   2.256   0.0299 *  
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Residual standard error: 210.3 on 38 degrees of freedom
## Multiple R-squared:  0.1181,	Adjusted R-squared:  0.09489 
## F-statistic: 5.089 on 1 and 38 DF,  p-value: 0.02992
```

```r
#run linear model with covariates of sex and age
  fit_rct <- lm(Thalamus_change ~ RandomArm + sex + age, data= RCT_SubCort)
  summary(fit_rct)
```

```
## 
## Call:
## lm(formula = Thalamus_change ~ RandomArm + sex + age, data = RCT_SubCort)
## 
## Residuals:
##     Min      1Q  Median      3Q     Max 
## -484.06 -102.67   17.14  139.25  326.59 
## 
## Coefficients:
##                  Estimate Std. Error t value Pr(>|t|)  
## (Intercept)      -188.053    134.129  -1.402   0.1695  
## RandomArmPlacebo  169.144     72.061   2.347   0.0245 *
## sexM               30.398     69.369   0.438   0.6639  
## age                -2.146      2.488  -0.863   0.3940  
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Residual standard error: 213.7 on 36 degrees of freedom
## Multiple R-squared:  0.1377,	Adjusted R-squared:  0.06582 
## F-statistic: 1.916 on 3 and 36 DF,  p-value: 0.1444
```

```r
#run linear model with covariates of sex and age
  fit_rct <- lm(Thalamus_change ~ RandomArm + sex + age + site, data= RCT_SubCort)
  summary(fit_rct)
```

```
## 
## Call:
## lm(formula = Thalamus_change ~ RandomArm + sex + age + site, 
##     data = RCT_SubCort)
## 
## Residuals:
##     Min      1Q  Median      3Q     Max 
## -396.26 -147.24   65.82  121.84  294.44 
## 
## Coefficients:
##                  Estimate Std. Error t value Pr(>|t|)  
## (Intercept)      -302.140    144.384  -2.093   0.0442 *
## RandomArmPlacebo  145.518     70.598   2.061   0.0472 *
## sexM                7.308     67.727   0.108   0.9147  
## age                -1.325      2.636  -0.503   0.6186  
## siteMAS           171.904     90.034   1.909   0.0649 .
## siteNKI           179.824     86.043   2.090   0.0444 *
## sitePMC            90.282    105.706   0.854   0.3992  
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Residual standard error: 205.3 on 33 degrees of freedom
## Multiple R-squared:   0.27,	Adjusted R-squared:  0.1373 
## F-statistic: 2.035 on 6 and 33 DF,  p-value: 0.08872
```
#### Striatum


```r
#run linear model without covariates
  fit_rct <- lm(Striatum_change ~ RandomArm, data= RCT_SubCort)
  print(fit_rct)
```

```
## 
## Call:
## lm(formula = Striatum_change ~ RandomArm, data = RCT_SubCort)
## 
## Coefficients:
##      (Intercept)  RandomArmPlacebo  
##          -296.60            -91.83
```

```r
  summary(fit_rct)
```

```
## 
## Call:
## lm(formula = Striatum_change ~ RandomArm, data = RCT_SubCort)
## 
## Residuals:
##     Min      1Q  Median      3Q     Max 
## -756.70 -141.60    8.11  153.80  542.53 
## 
## Coefficients:
##                  Estimate Std. Error t value Pr(>|t|)    
## (Intercept)       -296.60      55.23  -5.370 4.15e-06 ***
## RandomArmPlacebo   -91.83      93.35  -0.984    0.331    
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Residual standard error: 281.6 on 38 degrees of freedom
## Multiple R-squared:  0.02483,	Adjusted R-squared:  -0.0008284 
## F-statistic: 0.9677 on 1 and 38 DF,  p-value: 0.3315
```

```r
#run linear model with covariates of sex and age
  fit_rct <- lm(Striatum_change ~ RandomArm + sex + age, data= RCT_SubCort)
  print(fit_rct)
```

```
## 
## Call:
## lm(formula = Striatum_change ~ RandomArm + sex + age, data = RCT_SubCort)
## 
## Coefficients:
##      (Intercept)  RandomArmPlacebo              sexM               age  
##         -169.875           -80.518           -12.233            -2.318
```

```r
  summary(fit_rct)
```

```
## 
## Call:
## lm(formula = Striatum_change ~ RandomArm + sex + age, data = RCT_SubCort)
## 
## Residuals:
##     Min      1Q  Median      3Q     Max 
## -772.16 -141.82    9.74  154.93  552.85 
## 
## Coefficients:
##                  Estimate Std. Error t value Pr(>|t|)
## (Intercept)      -169.875    180.214  -0.943    0.352
## RandomArmPlacebo  -80.518     96.821  -0.832    0.411
## sexM              -12.233     93.203  -0.131    0.896
## age                -2.318      3.343  -0.693    0.492
## 
## Residual standard error: 287.1 on 36 degrees of freedom
## Multiple R-squared:  0.0397,	Adjusted R-squared:  -0.04032 
## F-statistic: 0.4961 on 3 and 36 DF,  p-value: 0.6873
```

```r
#run linear model with covariates of sex and age
  fit_rct <- lm(Striatum_change ~ RandomArm + sex + age + site, data= RCT_SubCort)
  print(fit_rct)
```

```
## 
## Call:
## lm(formula = Striatum_change ~ RandomArm + sex + age + site, 
##     data = RCT_SubCort)
## 
## Coefficients:
##      (Intercept)  RandomArmPlacebo              sexM               age  
##         -147.226           -69.518           -25.309            -3.412  
##          siteMAS           siteNKI           sitePMC  
##          -21.034            84.805           157.298
```

```r
  summary(fit_rct)
```

```
## 
## Call:
## lm(formula = Striatum_change ~ RandomArm + sex + age + site, 
##     data = RCT_SubCort)
## 
## Residuals:
##     Min      1Q  Median      3Q     Max 
## -742.29 -116.85   26.69  199.56  431.93 
## 
## Coefficients:
##                  Estimate Std. Error t value Pr(>|t|)
## (Intercept)      -147.226    205.644  -0.716    0.479
## RandomArmPlacebo  -69.518    100.552  -0.691    0.494
## sexM              -25.309     96.463  -0.262    0.795
## age                -3.412      3.755  -0.909    0.370
## siteMAS           -21.034    128.235  -0.164    0.871
## siteNKI            84.805    122.550   0.692    0.494
## sitePMC           157.298    150.555   1.045    0.304
## 
## Residual standard error: 292.5 on 33 degrees of freedom
## Multiple R-squared:  0.08654,	Adjusted R-squared:  -0.07955 
## F-statistic: 0.521 on 6 and 33 DF,  p-value: 0.7881
```
#### Hippocampus


```r
#run linear model without covariates
  fit_rct <- lm(Hippocampus_change ~ RandomArm, data= RCT_SubCort)
  summary(fit_rct)
```

```
## 
## Call:
## lm(formula = Hippocampus_change ~ RandomArm, data = RCT_SubCort)
## 
## Residuals:
##     Min      1Q  Median      3Q     Max 
## -490.91  -89.61    7.16   94.86  270.16 
## 
## Coefficients:
##                  Estimate Std. Error t value Pr(>|t|)    
## (Intercept)       -111.36      27.95  -3.984 0.000296 ***
## RandomArmPlacebo    86.58      47.25   1.832 0.074752 .  
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Residual standard error: 142.5 on 38 degrees of freedom
## Multiple R-squared:  0.08118,	Adjusted R-squared:  0.057 
## F-statistic: 3.357 on 1 and 38 DF,  p-value: 0.07475
```

```r
#run linear model with covariates of sex and age
  fit_rct <- lm(Hippocampus_change ~ RandomArm + sex + age, data= RCT_SubCort)
  summary(fit_rct)
```

```
## 
## Call:
## lm(formula = Hippocampus_change ~ RandomArm + sex + age, data = RCT_SubCort)
## 
## Residuals:
##     Min      1Q  Median      3Q     Max 
## -498.44  -90.15   -0.46   84.67  257.79 
## 
## Coefficients:
##                  Estimate Std. Error t value Pr(>|t|)  
## (Intercept)         5.198     89.606   0.058   0.9541  
## RandomArmPlacebo   97.320     48.141   2.022   0.0507 .
## sexM               -6.910     46.342  -0.149   0.8823  
## age                -2.171      1.662  -1.306   0.1998  
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Residual standard error: 142.8 on 36 degrees of freedom
## Multiple R-squared:  0.1268,	Adjusted R-squared:  0.05408 
## F-statistic: 1.743 on 3 and 36 DF,  p-value: 0.1756
```

```r
#run linear model with covariates of sex and age
  fit_rct <- lm(Hippocampus_change ~ RandomArm + sex + age + site, data= RCT_SubCort)
  summary(fit_rct)
```

```
## 
## Call:
## lm(formula = Hippocampus_change ~ RandomArm + sex + age + site, 
##     data = RCT_SubCort)
## 
## Residuals:
##     Min      1Q  Median      3Q     Max 
## -491.31  -98.13    3.17   97.16  258.53 
## 
## Coefficients:
##                  Estimate Std. Error t value Pr(>|t|)  
## (Intercept)        32.179    103.540   0.311   0.7579  
## RandomArmPlacebo   97.302     50.627   1.922   0.0633 .
## sexM              -10.971     48.568  -0.226   0.8227  
## age                -2.748      1.891  -1.454   0.1555  
## siteMAS            22.156     64.565   0.343   0.7337  
## siteNKI           -24.291     61.703  -0.394   0.6964  
## sitePMC            47.155     75.803   0.622   0.5382  
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Residual standard error: 147.3 on 33 degrees of freedom
## Multiple R-squared:  0.1483,	Adjusted R-squared:  -0.006514 
## F-statistic: 0.9579 on 6 and 33 DF,  p-value: 0.4684
```

----

## RCT & Relapse (with time as factor)

### Thalamus



```r
#restructure data for RCT & Relapse participants (N=72)
  RCTRelapse_Thalamus <- ana_df %>%
    gather(oldcolname, volume, Thalamus_01, Thalamus_02) %>%
    mutate(model_days = if_else(oldcolname == "Thalamus_01", 1, dateDiff))

RCTRelapse_Thalamus %>% filter(model_days == 1) %>% count(randomization) %>% knitr::kable() 
```


\begin{tabular}{l|r}
\hline
randomization & n\\
\hline
O & 38\\
\hline
P & 34\\
\hline
\end{tabular}


```r
#plot all data, including outlier (participant 210030)
  RCTRelapse_Thalamus %>%
   ggplot(aes(x=model_days, y=volume, colour=RandomArm)) + 
   geom_point() + 
   geom_line(aes(group=STUDYID), alpha = 0.5) + 
   geom_smooth(method="lm", formula=y~poly(x,1)) +
   ggtitle("Volume of Thalamus over time") +
   labs(x = "Days between MRIs", y = "Volume", colour = NULL) +
   theme_minimal()
```

![](10_STOPPD_freesurfer_subcortical_files/figure-latex/RCTRelapse_Thalamus_plot-1.pdf)<!-- --> 


```r
#run mixed linear model, with covariates
  fit_all <- lmer(volume ~ RandomArm*model_days + sex + age + (1|STUDYID), data= RCTRelapse_Thalamus)
  print(fit_all)
```

```
## Linear mixed model fit by REML ['lmerModLmerTest']
## Formula: volume ~ RandomArm * model_days + sex + age + (1 | STUDYID)
##    Data: RCTRelapse_Thalamus
## REML criterion at convergence: 2193.44
## Random effects:
##  Groups   Name        Std.Dev.
##  STUDYID  (Intercept) 1396.1  
##  Residual              171.2  
## Number of obs: 144, groups:  STUDYID, 72
## Fixed Effects:
##                 (Intercept)             RandomArmPlacebo  
##                  15637.4008                    -181.4066  
##                  model_days                         sexM  
##                     -1.1851                    1898.4711  
##                         age  RandomArmPlacebo:model_days  
##                    -58.6352                       0.8093
```

```r
  summary(fit_all)
```

```
## Linear mixed model fit by REML. t-tests use Satterthwaite's method [
## lmerModLmerTest]
## Formula: volume ~ RandomArm * model_days + sex + age + (1 | STUDYID)
##    Data: RCTRelapse_Thalamus
## 
## REML criterion at convergence: 2193.4
## 
## Scaled residuals: 
##     Min      1Q  Median      3Q     Max 
## -3.6119 -0.4282  0.0145  0.3924  3.5470 
## 
## Random effects:
##  Groups   Name        Variance Std.Dev.
##  STUDYID  (Intercept) 1949229  1396.1  
##  Residual               29305   171.2  
## Number of obs: 144, groups:  STUDYID, 72
## 
## Fixed effects:
##                               Estimate Std. Error         df t value
## (Intercept)                 15637.4008   645.0482    68.1290  24.242
## RandomArmPlacebo             -181.4066   332.4480    68.7735  -0.546
## model_days                     -1.1851     0.1804    70.1127  -6.570
## sexM                         1898.4711   332.6682    67.9899   5.707
## age                           -58.6352    10.8628    67.9907  -5.398
## RandomArmPlacebo:model_days     0.8093     0.3046    70.2726   2.657
##                             Pr(>|t|)    
## (Intercept)                  < 2e-16 ***
## RandomArmPlacebo             0.58706    
## model_days                  7.46e-09 ***
## sexM                        2.75e-07 ***
## age                         9.26e-07 ***
## RandomArmPlacebo:model_days  0.00976 ** 
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Correlation of Fixed Effects:
##             (Intr) RndmAP mdl_dy sexM   age   
## RndmArmPlcb -0.201                            
## model_days  -0.032  0.056                     
## sexM        -0.172  0.037  0.000              
## age         -0.903 -0.054  0.003 -0.079       
## RndmArmPl:_  0.018 -0.074 -0.592  0.002 -0.001
```

```r
#run mixed linear model, with covariates
  fit_all <- lmer(volume ~ RandomArm*model_days + sex + age + site + (1|STUDYID), data= RCTRelapse_Thalamus)
  print(fit_all)
```

```
## Linear mixed model fit by REML ['lmerModLmerTest']
## Formula: 
## volume ~ RandomArm * model_days + sex + age + site + (1 | STUDYID)
##    Data: RCTRelapse_Thalamus
## REML criterion at convergence: 2140.5
## Random effects:
##  Groups   Name        Std.Dev.
##  STUDYID  (Intercept) 1313.0  
##  Residual              171.2  
## Number of obs: 144, groups:  STUDYID, 72
## Fixed Effects:
##                 (Intercept)             RandomArmPlacebo  
##                  15877.4970                    -149.6319  
##                  model_days                         sexM  
##                     -1.1878                    1844.5154  
##                         age                      siteMAS  
##                    -61.8649                    -822.5404  
##                     siteNKI                      sitePMC  
##                    821.9333                      37.9774  
## RandomArmPlacebo:model_days  
##                      0.8111
```

```r
  summary(fit_all)
```

```
## Linear mixed model fit by REML. t-tests use Satterthwaite's method [
## lmerModLmerTest]
## Formula: 
## volume ~ RandomArm * model_days + sex + age + site + (1 | STUDYID)
##    Data: RCTRelapse_Thalamus
## 
## REML criterion at convergence: 2140.5
## 
## Scaled residuals: 
##     Min      1Q  Median      3Q     Max 
## -3.5816 -0.4404  0.0083  0.3873  3.5784 
## 
## Random effects:
##  Groups   Name        Variance Std.Dev.
##  STUDYID  (Intercept) 1724006  1313.0  
##  Residual               29304   171.2  
## Number of obs: 144, groups:  STUDYID, 72
## 
## Fixed effects:
##                               Estimate Std. Error         df t value
## (Intercept)                 15877.4970   637.1250    65.1359  24.921
## RandomArmPlacebo             -149.6319   317.2828    65.8206  -0.472
## model_days                     -1.1878     0.1804    70.1251  -6.586
## sexM                         1844.5154   317.2625    64.9936   5.814
## age                           -61.8649    10.2955    64.9952  -6.009
## siteMAS                      -822.5404   402.5790    64.9957  -2.043
## siteNKI                       821.9333   446.1341    64.9957   1.842
## sitePMC                        37.9774   459.0477    64.9931   0.083
## RandomArmPlacebo:model_days     0.8111     0.3046    70.3103   2.663
##                             Pr(>|t|)    
## (Intercept)                  < 2e-16 ***
## RandomArmPlacebo              0.6388    
## model_days                  7.00e-09 ***
## sexM                        2.02e-07 ***
## age                         9.35e-08 ***
## siteMAS                       0.0451 *  
## siteNKI                       0.0700 .  
## sitePMC                       0.9343    
## RandomArmPlacebo:model_days   0.0096 ** 
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Correlation of Fixed Effects:
##             (Intr) RndmAP mdl_dy sexM   age    sitMAS sitNKI sitPMC
## RndmArmPlcb -0.140                                                 
## model_days  -0.033  0.059                                          
## sexM        -0.130  0.055  0.001                                   
## age         -0.882 -0.066  0.004 -0.076                            
## siteMAS     -0.292 -0.147  0.004 -0.066  0.108                     
## siteNKI     -0.175 -0.119 -0.004 -0.119  0.010  0.357              
## sitePMC     -0.153 -0.089  0.000 -0.144 -0.009  0.343  0.319       
## RndmArmPl:_  0.018 -0.078 -0.592  0.001 -0.001  0.000  0.004  0.001
```

### Striatum



```r
#restructure data for RCT & Relapse participants (N=72)
  RCTRelapse_Striatum <- ana_df %>%
    gather(oldcolname, volume, Striatum_01, Striatum_02) %>%
    mutate(model_days = if_else(oldcolname == "Striatum_01", 1, dateDiff))

RCTRelapse_Striatum %>% filter(model_days == 1) %>% count(randomization) %>% knitr::kable() 
```


\begin{tabular}{l|r}
\hline
randomization & n\\
\hline
O & 38\\
\hline
P & 34\\
\hline
\end{tabular}


```r
#plot all data, including outlier (participant 210030)
  RCTRelapse_Striatum %>%
   ggplot(aes(x=model_days, y=volume, colour=RandomArm)) + 
   geom_point() + 
   geom_line(aes(group=STUDYID), alpha = 0.5) + 
   geom_smooth(method="lm", formula=y~poly(x,1)) +
   ggtitle("Volume of Striatum over time") +
   labs(x = "Days between MRIs", y = "Volume", colour = NULL) +
   theme_minimal()
```

![](10_STOPPD_freesurfer_subcortical_files/figure-latex/RCTRelapse_Striatum_plot-1.pdf)<!-- --> 


```r
#run mixed linear model, with covariates
  fit_all <- lmer(volume ~ RandomArm*model_days + sex + age + (1|STUDYID), data= RCTRelapse_Striatum)
  print(fit_all)
```

```
## Linear mixed model fit by REML ['lmerModLmerTest']
## Formula: volume ~ RandomArm * model_days + sex + age + (1 | STUDYID)
##    Data: RCTRelapse_Striatum
## REML criterion at convergence: 2250.495
## Random effects:
##  Groups   Name        Std.Dev.
##  STUDYID  (Intercept) 1674.5  
##  Residual              215.6  
## Number of obs: 144, groups:  STUDYID, 72
## Fixed Effects:
##                 (Intercept)             RandomArmPlacebo  
##                  18786.7886                    -199.5106  
##                  model_days                         sexM  
##                     -1.1427                    1603.1890  
##                         age  RandomArmPlacebo:model_days  
##                    -47.6384                      -0.1886
```

```r
  summary(fit_all)
```

```
## Linear mixed model fit by REML. t-tests use Satterthwaite's method [
## lmerModLmerTest]
## Formula: volume ~ RandomArm * model_days + sex + age + (1 | STUDYID)
##    Data: RCTRelapse_Striatum
## 
## REML criterion at convergence: 2250.5
## 
## Scaled residuals: 
##      Min       1Q   Median       3Q      Max 
## -2.37595 -0.36241  0.03762  0.33458  2.82797 
## 
## Random effects:
##  Groups   Name        Variance Std.Dev.
##  STUDYID  (Intercept) 2804106  1674.5  
##  Residual               46470   215.6  
## Number of obs: 144, groups:  STUDYID, 72
## 
## Fixed effects:
##                               Estimate Std. Error         df t value
## (Intercept)                 18786.7886   774.0095    68.1517  24.272
## RandomArmPlacebo             -199.5106   399.0078    68.8613  -0.500
## model_days                     -1.1427     0.2271    70.1337  -5.031
## sexM                         1603.1890   399.1561    67.9984   4.016
## age                           -47.6384    13.0339    67.9994  -3.655
## RandomArmPlacebo:model_days    -0.1886     0.3836    70.3099  -0.492
##                             Pr(>|t|)    
## (Intercept)                  < 2e-16 ***
## RandomArmPlacebo            0.618657    
## model_days                  3.61e-06 ***
## sexM                        0.000150 ***
## age                         0.000502 ***
## RandomArmPlacebo:model_days 0.624440    
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Correlation of Fixed Effects:
##             (Intr) RndmAP mdl_dy sexM   age   
## RndmArmPlcb -0.202                            
## model_days  -0.034  0.059                     
## sexM        -0.172  0.036  0.000              
## age         -0.903 -0.054  0.003 -0.079       
## RndmArmPl:_  0.019 -0.078 -0.592  0.002 -0.001
```

```r
#run mixed linear model, with covariates
  fit_all <- lmer(volume ~ RandomArm*model_days + sex + age + site + (1|STUDYID), data= RCTRelapse_Striatum)
  print(fit_all)
```

```
## Linear mixed model fit by REML ['lmerModLmerTest']
## Formula: 
## volume ~ RandomArm * model_days + sex + age + site + (1 | STUDYID)
##    Data: RCTRelapse_Striatum
## REML criterion at convergence: 2204.864
## Random effects:
##  Groups   Name        Std.Dev.
##  STUDYID  (Intercept) 1680.8  
##  Residual              215.6  
## Number of obs: 144, groups:  STUDYID, 72
## Fixed Effects:
##                 (Intercept)             RandomArmPlacebo  
##                  19124.4744                    -100.8250  
##                  model_days                         sexM  
##                     -1.1423                    1673.9180  
##                         age                      siteMAS  
##                    -48.9143                    -595.3613  
##                     siteNKI                      sitePMC  
##                   -807.5124                    -307.2624  
## RandomArmPlacebo:model_days  
##                     -0.1903
```

```r
  summary(fit_all)
```

```
## Linear mixed model fit by REML. t-tests use Satterthwaite's method [
## lmerModLmerTest]
## Formula: 
## volume ~ RandomArm * model_days + sex + age + site + (1 | STUDYID)
##    Data: RCTRelapse_Striatum
## 
## REML criterion at convergence: 2204.9
## 
## Scaled residuals: 
##      Min       1Q   Median       3Q      Max 
## -2.38638 -0.35547  0.03791  0.32934  2.81721 
## 
## Random effects:
##  Groups   Name        Variance Std.Dev.
##  STUDYID  (Intercept) 2825222  1680.8  
##  Residual               46470   215.6  
## Number of obs: 144, groups:  STUDYID, 72
## 
## Fixed effects:
##                               Estimate Std. Error         df t value
## (Intercept)                 19124.4744   815.4829    65.1363  23.452
## RandomArmPlacebo             -100.8250   406.0694    65.7992  -0.248
## model_days                     -1.1423     0.2271    70.1259  -5.029
## sexM                         1673.9180   406.0846    64.9986   4.122
## age                           -48.9143    13.1779    65.0002  -3.712
## siteMAS                      -595.3613   515.2865    65.0006  -1.155
## siteNKI                      -807.5124   571.0354    65.0006  -1.414
## sitePMC                      -307.2624   587.5646    64.9981  -0.523
## RandomArmPlacebo:model_days    -0.1903     0.3836    70.3052  -0.496
##                             Pr(>|t|)    
## (Intercept)                  < 2e-16 ***
## RandomArmPlacebo            0.804680    
## model_days                  3.63e-06 ***
## sexM                        0.000109 ***
## age                         0.000429 ***
## siteMAS                     0.252157    
## siteNKI                     0.162099    
## sitePMC                     0.602793    
## RandomArmPlacebo:model_days 0.621329    
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Correlation of Fixed Effects:
##             (Intr) RndmAP mdl_dy sexM   age    sitMAS sitNKI sitPMC
## RndmArmPlcb -0.140                                                 
## model_days  -0.033  0.058                                          
## sexM        -0.130  0.055  0.001                                   
## age         -0.882 -0.066  0.004 -0.076                            
## siteMAS     -0.292 -0.147  0.004 -0.066  0.108                     
## siteNKI     -0.175 -0.119 -0.004 -0.119  0.010  0.357              
## sitePMC     -0.153 -0.089  0.000 -0.144 -0.009  0.343  0.319       
## RndmArmPl:_  0.017 -0.077 -0.592  0.001 -0.001  0.000  0.004  0.001
```

### Hippocampus



```r
#restructure data for RCT & Relapse participants (N=72)
  RCTRelapse_Hippocampus <- ana_df %>%
    gather(oldcolname, volume, Hippocampus_01, Hippocampus_02) %>%
    mutate(model_days = if_else(oldcolname == "Hippocampus_01", 1, dateDiff))

RCTRelapse_Hippocampus %>% filter(model_days == 1) %>% count(randomization) %>% knitr::kable() 
```


\begin{tabular}{l|r}
\hline
randomization & n\\
\hline
O & 38\\
\hline
P & 34\\
\hline
\end{tabular}


```r
#plot all data, including outlier (participant 210030)
  RCTRelapse_Hippocampus %>%
   ggplot(aes(x=model_days, y=volume, colour=RandomArm)) + 
   geom_point() + 
   geom_line(aes(group=STUDYID), alpha = 0.5) + 
   geom_smooth(method="lm", formula=y~poly(x,1)) +
   ggtitle("Volume of Hippocampus over time") +
   labs(x = "Days between MRIs", y = "Volume", colour = NULL) +
   theme_minimal()
```

![](10_STOPPD_freesurfer_subcortical_files/figure-latex/RCTRelapse_Hippocampus_plot-1.pdf)<!-- --> 


```r
#run mixed linear model, with covariates
  fit_all <- lmer(volume ~ RandomArm*model_days + sex + age + (1|STUDYID), data= RCTRelapse_Hippocampus)
  print(fit_all)
```

```
## Linear mixed model fit by REML ['lmerModLmerTest']
## Formula: volume ~ RandomArm * model_days + sex + age + (1 | STUDYID)
##    Data: RCTRelapse_Hippocampus
## REML criterion at convergence: 2049.905
## Random effects:
##  Groups   Name        Std.Dev.
##  STUDYID  (Intercept) 840.5   
##  Residual             100.6   
## Number of obs: 144, groups:  STUDYID, 72
## Fixed Effects:
##                 (Intercept)             RandomArmPlacebo  
##                   8989.3972                     -84.5737  
##                  model_days                         sexM  
##                     -0.4047                     607.4115  
##                         age  RandomArmPlacebo:model_days  
##                    -31.6143                       0.2634
```

```r
  summary(fit_all)
```

```
## Linear mixed model fit by REML. t-tests use Satterthwaite's method [
## lmerModLmerTest]
## Formula: volume ~ RandomArm * model_days + sex + age + (1 | STUDYID)
##    Data: RCTRelapse_Hippocampus
## 
## REML criterion at convergence: 2049.9
## 
## Scaled residuals: 
##      Min       1Q   Median       3Q      Max 
## -2.41294 -0.41410 -0.00566  0.40252  2.36012 
## 
## Random effects:
##  Groups   Name        Variance Std.Dev.
##  STUDYID  (Intercept) 706389   840.5   
##  Residual              10111   100.6   
## Number of obs: 144, groups:  STUDYID, 72
## 
## Fixed effects:
##                              Estimate Std. Error        df t value
## (Intercept)                 8989.3972   388.2347   68.1293  23.155
## RandomArmPlacebo             -84.5737   200.0678   68.7432  -0.423
## model_days                    -0.4047     0.1060   70.1138  -3.820
## sexM                         607.4115   200.2276   67.9968   3.034
## age                          -31.6143     6.5381   67.9975  -4.835
## RandomArmPlacebo:model_days    0.2634     0.1790   70.2661   1.472
##                             Pr(>|t|)    
## (Intercept)                  < 2e-16 ***
## RandomArmPlacebo            0.673815    
## model_days                  0.000285 ***
## sexM                        0.003420 ** 
## age                         7.94e-06 ***
## RandomArmPlacebo:model_days 0.145501    
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Correlation of Fixed Effects:
##             (Intr) RndmAP mdl_dy sexM   age   
## RndmArmPlcb -0.201                            
## model_days  -0.031  0.055                     
## sexM        -0.172  0.037  0.000              
## age         -0.903 -0.054  0.003 -0.079       
## RndmArmPl:_  0.017 -0.073 -0.592  0.002 -0.001
```

```r
#run mixed linear model, with covariates
  fit_all <- lmer(volume ~ RandomArm*model_days + sex + age + site + (1|STUDYID), data= RCTRelapse_Hippocampus)
  print(fit_all)
```

```
## Linear mixed model fit by REML ['lmerModLmerTest']
## Formula: 
## volume ~ RandomArm * model_days + sex + age + site + (1 | STUDYID)
##    Data: RCTRelapse_Hippocampus
## REML criterion at convergence: 2005.898
## Random effects:
##  Groups   Name        Std.Dev.
##  STUDYID  (Intercept) 827.3   
##  Residual             100.6   
## Number of obs: 144, groups:  STUDYID, 72
## Fixed Effects:
##                 (Intercept)             RandomArmPlacebo  
##                   8920.7798                    -106.7365  
##                  model_days                         sexM  
##                     -0.4048                     548.9133  
##                         age                      siteMAS  
##                    -32.1879                      46.0261  
##                     siteNKI                      sitePMC  
##                    116.0375                     631.6705  
## RandomArmPlacebo:model_days  
##                      0.2637
```

```r
  summary(fit_all)
```

```
## Linear mixed model fit by REML. t-tests use Satterthwaite's method [
## lmerModLmerTest]
## Formula: 
## volume ~ RandomArm * model_days + sex + age + site + (1 | STUDYID)
##    Data: RCTRelapse_Hippocampus
## 
## REML criterion at convergence: 2005.9
## 
## Scaled residuals: 
##      Min       1Q   Median       3Q      Max 
## -2.40104 -0.40874 -0.00025  0.39616  2.37266 
## 
## Random effects:
##  Groups   Name        Variance Std.Dev.
##  STUDYID  (Intercept) 684488   827.3   
##  Residual              10111   100.6   
## Number of obs: 144, groups:  STUDYID, 72
## 
## Fixed effects:
##                              Estimate Std. Error        df t value
## (Intercept)                 8920.7798   401.2060   65.1204  22.235
## RandomArmPlacebo            -106.7365   199.7292   65.7164  -0.534
## model_days                    -0.4048     0.1060   70.1110  -3.820
## sexM                         548.9133   199.7986   64.9966   2.747
## age                          -32.1879     6.4837   64.9980  -4.964
## siteMAS                       46.0261   253.5270   64.9984   0.182
## siteNKI                      116.0375   280.9562   64.9984   0.413
## sitePMC                      631.6705   289.0890   64.9961   2.185
## RandomArmPlacebo:model_days    0.2637     0.1790   70.2722   1.474
##                             Pr(>|t|)    
## (Intercept)                  < 2e-16 ***
## RandomArmPlacebo            0.594864    
## model_days                  0.000285 ***
## sexM                        0.007766 ** 
## age                         5.26e-06 ***
## siteMAS                     0.856506    
## siteNKI                     0.680959    
## sitePMC                     0.032497 *  
## RandomArmPlacebo:model_days 0.145021    
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Correlation of Fixed Effects:
##             (Intr) RndmAP mdl_dy sexM   age    sitMAS sitNKI sitPMC
## RndmArmPlcb -0.140                                                 
## model_days  -0.031  0.055                                          
## sexM        -0.130  0.055  0.001                                   
## age         -0.882 -0.066  0.004 -0.076                            
## siteMAS     -0.292 -0.147  0.004 -0.066  0.108                     
## siteNKI     -0.175 -0.119 -0.004 -0.119  0.010  0.357              
## sitePMC     -0.153 -0.089  0.000 -0.144 -0.009  0.343  0.319       
## RndmArmPl:_  0.017 -0.073 -0.592  0.001 -0.001  0.000  0.004  0.001
```
