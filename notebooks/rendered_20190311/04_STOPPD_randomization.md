# Report Randomization numbers

This script identifies the number of participants in olanzapine vs. placebo by scan timepoint, using the logic of group inclusion that Judy and Dielle provided, and Nick and Aristotle have agreed to.

__Note:__ this script includes data from all participants with data in Judy's master log and our file system. It has not excluded participants on any other basis (e.g., QC fail, processing fail, post-hoc clinical trial ineligibility, etc.)





## Identify baseline scans

First - identify the number of baseline scans (i.e., scans completed at week 20).


```r
#count the number of participants that have a 'yes' for 'completed' in "Scan.completed.1"
n_first_complete = sum(na.omit(df$first_complete == "Yes")) #88 participants completed week 20 scan

#for clarity, print the IDs of the N=88 participants that completed week 20 scans
(df %>% filter(first_complete == "Yes"))$STUDYID 
```

```
##  [1] 110008 110009 110013 110016 110022 110025 110028 110030 110031 110034
## [11] 120011 120012 120015 120016 120017 120021 120026 210012 210013 210014
## [21] 210017 210020 210022 210024 210026 210030 210033 210036 210038 210042
## [31] 210048 210049 210051 220002 220003 220004 220006 220008 220009 310010
## [41] 310015 310025 310037 310051 310070 320006 320013 320021 320022 320032
## [51] 320041 320042 320043 320045 410004 410008 410009 410010 410011 410012
## [61] 410013 410015 410019 410022 410023 410029 410030 410031 410037 410039
## [71] 410040 410043 410045 410047 420005 420007 420013 420016 420018 420019
## [81] 420020 420023 420029 420032 420039 420042 420043 420044
```

**The number of participants who completed their first scan is 88**

**RANDOMIZATION** - as expected, there's no difference in first scan completion between those randomized to O vs. P group


```r
#RANDOMIZATION - as expected, there's no difference in first scan completion between those randomized to O vs. P group
(R <- addmargins(table(df$first_complete == 'Yes', df$randomization))) #O = 45; P = 43 (total = 88)
```

```
##       
##         O  P Sum
##   TRUE 45 43  88
##   Sum  45 43  88
```

## Identify week 56 scans

Second - identify the number of week 56 scans (i.e., 36 weeks after week 20).



```r
#make sure that all the participants that completed week 56 scan also completed week 20
all_second_complete <- all((df$second_complete == "Yes") %in% (df$first_complete== "Yes")) #all TRUE

#count the number of participants that have a 'yes' for 'completed' in "Scan.completed" - but this includes 'relapse' and 'off protocol', as well as RCT or "true completers"
(n_second_complete <- sum(na.omit(df$second_complete == "Yes"))) #74 completed week 56 scan
```

```
## [1] 74
```

Subject ids of the n = 74 who completed their second scan. Note: it is TRUE that all participants who completed their second scan have baseline data.


```r
#for clarity, print the IDs of the N=74 participants that completed week 56 scans
(df %>% filter(second_complete == "Yes"))$STUDYID
```

```
##  [1] 110008 110009 110013 110022 110031 110034 120011 120012 120015 120016
## [11] 120017 120021 120026 210012 210013 210014 210017 210020 210022 210026
## [21] 210030 210033 210038 210042 210049 210051 220002 220003 220004 220006
## [31] 220009 310010 310015 310025 310037 310051 320006 320013 320021 320022
## [41] 320032 320042 320043 320045 410004 410008 410009 410010 410011 410012
## [51] 410013 410015 410019 410022 410023 410029 410030 410031 410037 410039
## [61] 410040 410043 420007 420013 420016 420018 420019 420020 420023 420029
## [71] 420032 420039 420042 420043
```



```r
#count how many participants that completed week 56 scan are classified as RCT 
sum(na.omit(df$second_complete == 'Yes' & df$second_timepoint == 'RCT')) #RCT = 41
```

```
## [1] 41
```

```r
  (as.vector(na.omit(df$STUDYID[df$second_complete == "Yes" & df$second_timepoint == 'RCT']))) #for clarity, print N=41 RCT participant IDs
```

```
##  [1] 110008 110009 110013 110022 110031 110034 120011 120012 120015 210012
## [11] 210013 210014 210017 210020 210030 210051 220004 310051 320006 320021
## [21] 320032 320042 320043 320045 410004 410008 410010 410013 410022 410023
## [31] 410029 410030 410037 410039 410043 420013 420020 420029 420039 420042
## [41] 420043
```

```r
sum(na.omit(df$second_complete == 'Yes' & df$second_timepoint == 'Relapse')) #Relapse = 28
```

```
## [1] 28
```

```r
  (as.vector(na.omit(df$STUDYID[df$second_complete == "Yes" & df$second_timepoint == 'Relapse']))) #for clarity, print N=28 Relapse participant IDs
```

```
##  [1] 120016 120017 120021 120026 210022 210026 210033 210038 210042 210049
## [11] 220002 220003 220006 220009 310010 310015 310025 310037 320013 410009
## [21] 410011 410012 410031 410040 420007 420016 420023 420032
```

```r
sum(na.omit(df$second_complete == 'Yes' & df$second_timepoint == 'Off protocol')) #Off protocol = 5
```

```
## [1] 5
```

```r
  (as.vector(na.omit(df$STUDYID[df$second_complete == "Yes" & df$second_timepoint == 'Off protocol']))) #for clarity, print N=5 Off protocol participant IDs
```

```
## [1] 320022 410015 410019 420018 420019
```

```r
df %>% 
  filter(second_complete == "Yes") %>%
  count(second_timepoint) %>%
  kable(caption = "breakdown of those who where scanned at two timepoints")
```

\begin{table}[t]

\caption{(\#tab:SecondScan-subcounts)breakdown of those who where scanned at two timepoints}
\centering
\begin{tabular}{l|r}
\hline
second\_timepoint & n\\
\hline
Off protocol & 5\\
\hline
RCT & 41\\
\hline
Relapse & 28\\
\hline
\end{tabular}
\end{table}

```r
#RANDOMIZATION- look at randomization info for those who completed a second timepoint RCT scan
(R <- addmargins(table(df$second_complete == 'Yes' & df$second_timepoint == 'RCT', df$randomization))) #O = 27; P = 14 (total = 41)
```

```
##        
##          O  P Sum
##   FALSE 14 24  38
##   TRUE  27 14  41
##   Sum   41 38  79
```

```r
df %>% 
  filter(second_complete == "Yes") %>%
  count(second_timepoint, randomization) %>%
  kable(caption = "breakdown of those who where scanned at two timepoints, by arm")
```

\begin{table}[t]

\caption{(\#tab:SecondScan-subcounts)breakdown of those who where scanned at two timepoints, by arm}
\centering
\begin{tabular}{l|l|r}
\hline
second\_timepoint & randomization & n\\
\hline
Off protocol & O & 4\\
\hline
Off protocol & P & 1\\
\hline
RCT & O & 27\\
\hline
RCT & P & 14\\
\hline
Relapse & O & 8\\
\hline
Relapse & P & 20\\
\hline
\end{tabular}
\end{table}

## Identify off label scans

Third - identify the number of "off label" scans also at week 56.


```r
#make sure timepoint is a character
df$second_timepoint <- as.character(df$second_timepoint)

#count the number of scans completed at *third* timepoint, which are by definition "off label"
n_offlable <- sum(na.omit(df$third_complete == 'Yes')) #8 off-label scans

#for clarity, print the IDs of the N=8 participants that completed off-label scans
(as.vector(na.omit(df$STUDYID[df$third_complete == "Yes"])))
```

```
## [1] 110016 210033 210049 220006 310037 320022 410019 420032
```

```r
#of these, determine how many "off protocol" vs. "relapse", based on second timepoint scan
sum(na.omit(df$third_complete == 'Yes' & df$second_timepoint  == 'Off protocol')) #2 "off protocol" scans
```

```
## [1] 2
```

```r
  (as.vector(na.omit(df$STUDYID[df$third_complete == "Yes" & df$second_timepoint  == 'Off protocol'])))
```

```
## [1] 320022 410019
```

```r
sum(na.omit(df$third_complete == 'Yes' & df$second_timepoint == 'Relapse')) #6 relapse scans
```

```
## [1] 6
```

```r
  (as.vector(na.omit(df$STUDYID[df$third_complete == "Yes" & df$second_timepoint  == 'Relapse'])))
```

```
## [1] 110016 210033 210049 220006 310037 420032
```

```r
#RANDOMIZATION 
df %>%
  filter(third_complete == "Yes") %>%
  count(randomization) %>%
  kable(caption = str_c("Breakdown of thrid timepoint off-label scans ", n_offlable, " total"))
```

\begin{table}[t]

\caption{(\#tab:ThirdScanOffLabel)Breakdown of thrid timepoint off-label scans 8 total}
\centering
\begin{tabular}{l|r}
\hline
randomization & n\\
\hline
O & 3\\
\hline
P & 5\\
\hline
\end{tabular}
\end{table}

```r
df %>%
  filter(df$second_timepoint == 'Off protocol') %>%
  count(randomization, third_complete) %>%
  kable(caption = str_c("Breakdown of off-protocol scans by presence of third timepoint"))
```

\begin{table}[t]

\caption{(\#tab:ThirdScanOffLabel)Breakdown of off-protocol scans by presence of third timepoint}
\centering
\begin{tabular}{l|l|r}
\hline
randomization & third\_complete & n\\
\hline
O & Yes & 2\\
\hline
O & NA & 2\\
\hline
P & NA & 1\\
\hline
\end{tabular}
\end{table}

```r
df %>%
  filter(df$second_timepoint == 'Relapse') %>%
  count(randomization, third_complete) %>%
  kable(caption = str_c("Breakdown of thrid timepoint 'Relapse' scans by presence of third timepoint"))
```

\begin{table}[t]

\caption{(\#tab:ThirdScanOffLabel)Breakdown of thrid timepoint 'Relapse' scans by presence of third timepoint}
\centering
\begin{tabular}{l|l|r}
\hline
randomization & third\_complete & n\\
\hline
O & Yes & 1\\
\hline
O & NA & 9\\
\hline
P & Yes & 5\\
\hline
P & NA & 18\\
\hline
\end{tabular}
\end{table}

## Identify "Relapse" Scans

Identify the scans completed between week 20 and week 56 which are the relapse scans (and in a small minority of cases may be a scan when somebody is moving or wants out of the study despite being well). 


```r
#count relapse - note: both 'relapse' and 'off protocol' is included here (everything other than 'RCT')
sum(na.omit((df$second_timepoint == 'Relapse' | df$second_timepoint == 'Off protocol') & df$second_complete == 'Yes')) #33 participants relapsed/off protocol
```

```
## [1] 33
```

```r
#of these, count how many were "relapse" and how many were "off protocol"
sum(na.omit(df$second_timepoint == 'Relapse' & df$second_complete == 'Yes'))# 28 relapse
```

```
## [1] 28
```

```r
sum(na.omit(df$second_timepoint == 'Off protocol' & df$second_complete == 'Yes'))#5 off protocol
```

```
## [1] 5
```

```r
#RANDOMIZATION 
(R <- addmargins(table((df$second_timepoint == 'Relapse' | df$second_timepoint == 'Off protocol') & df$second_complete == 'Yes', df$randomization))) #relapse & off-protocol : O = 12; P = 21 (total = 33)
```

```
##        
##          O  P Sum
##   FALSE 28 15  43
##   TRUE  12 21  33
##   Sum   40 36  76
```

```r
(R <- addmargins(table(df$second_timepoint == 'Relapse' & df$second_complete == 'Yes', df$randomization))) #relapse: O = 8; P = 20 (total = 28)
```

```
##        
##          O  P Sum
##   FALSE 32 16  48
##   TRUE   8 20  28
##   Sum   40 36  76
```

```r
(R <- addmargins(table(df$second_timepoint == 'Off protocol' & df$second_complete == 'Yes', df$randomization))) #relapse: O = 4; P = 1 (total = 5)
```

```
##        
##          O  P Sum
##   FALSE 38 38  76
##   TRUE   4  1   5
##   Sum   42 39  81
```

```r
rm(df, R)
```




