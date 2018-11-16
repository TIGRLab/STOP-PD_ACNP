# Mangle Freesurfer Outputs

This script pulls together completion information alongside cortical thickness (CT) values and demographic information, for statistical purposes (error calculations). It is required for subsequent CT analyses. It was made in preparation for, and discussed at, the meeting with Jason Lerch. 


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
CT <- read_csv('../data/fs-enigma-long_201811/CorticalMeasuresENIGMA_ThickAvg.csv') #bring in CT data, from pipelines
```

```
## Parsed with column specification:
## cols(
##   .default = col_double(),
##   SubjID = col_character()
## )
## See spec(...) for full column specifications.
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
         STUDYID = parse_character(STUDYID)) %>%
  rename(category = "second_timepoint") %>%
  select(STUDYID, randomization, sex, age, category, offLabel, dateDiff)
```

## cleaning the CT data



```r
# separating the subject id and anything afterwards to identify the longtudinal pipeline participants
# separating the subject id into site, "STUDYID" and timepoint columns
# filtering (two steps) to only include the longitudinal pipeline data
CT_long <- CT %>%
  separate(SubjID, into = c("subid", "longitudinal_pipe"), sep = '\\.', extra = "drop", fill = "right") %>%
  separate(subid, into = c("study", "site", "STUDYID", "timepoint"), fill = "right") %>%
  filter(longitudinal_pipe == "long") %>%
  filter(timepoint != "00", timepoint != "03", timepoint != "")


# move CT from long to wide format
CT_wide <- CT_long %>%
  gather(region, thickness, ends_with('thickavg'), LThickness, RThickness, LSurfArea, RSurfArea, ICV) %>%
  spread(timepoint, thickness) %>%
  mutate(change = `02` - `01`) %>%
  gather(timepoint, thickness, `01`, `02`, change) %>%
  unite(newcolnames, region, timepoint) %>%
  spread(newcolnames, thickness)
```



```r
# merge CT values with df
ana_df <- inner_join(df, CT_wide, by='STUDYID')

# write.csv
write_csv(ana_df, '../generated_csvs/STOPPD_participantsCT_20181111.csv')
```

## report any mising values from clinical trial sample


```r
anti_join(df, CT_wide, by='STUDYID') %>%
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
  filter(is.na(LThickness_01)) %>%
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
  filter(is.na(LThickness_02)) %>%
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

## creating an control error term calculating spreadsheet


```r
## identify the repeat control in a column and mangle the STUDYID to match in a new column
CT_long1 <- CT_long %>%
  mutate(repeat_run = if_else(str_sub(STUDYID,1,1)=="R", "02", "01"),
         STUDYID = str_replace(STUDYID, 'R',"")) 

## extra the repeat study ids as a character vector
repeat_ids <- filter(CT_long1, repeat_run == "02")$STUDYID

## filter for only the subjects who are in the repeats list then switch to wide format
CT_wide_controls <- CT_long1 %>%
  filter(STUDYID %in% repeat_ids) %>% 
  gather(region, thickness, ends_with('thickavg'), LThickness, RThickness, LSurfArea, RSurfArea, ICV) %>%
  unite(newcolnames, region, repeat_run) %>%
  spread(newcolnames, thickness)

#write.csv
  write.csv(CT_wide_controls, '../generated_csvs/STOPPD_errorControls_2018-11-05.csv', row.names = FALSE)
```

