# Verifying number of scans

## Checking the TIGRLab "/archive/data"

This script pulls in and cleans up the naming of STOPPD scans as they exist in the Kimel lab file system. At earlier stages, this script helped us identify naming errors in the file system (all have since been fixed).

**Purpose:** The contents of the file system will, in other scripts, be checked against (1) the scans we have in XNAT, to ensure that there are no discrepancies between these databases, and also against (2) our subject inclusion list.


```r
library('stringi')
library('stringr')
library('plyr')
library('tidyr')
```


```r
#import spreadsheet ('ls' of file system)
terminal <- read.csv('../data/stoppd_NiiFolderContents_2018-01-25.csv', header = TRUE, stringsAsFactors = FALSE)
```


```r
#make a new column for site component of ID
terminal$site <- str_sub(terminal$scan_id, 8, 10)

#cut out study and site component from ID (first 11 characters)
terminal$scan_id <- substring(terminal$scan_id, 12)

#make a new column for session component of ID
terminal$session <- str_sub(terminal$scan_id, -2)

#cut out session information from ID (last 3 characters)
terminal$scan_id <- stri_sub(terminal$scan_id, 1, -4)

#make a new column that captures alphabetic component of ID ('R')
terminal$contains_R <- grepl('R', terminal$scan_id, fixed=TRUE) #36 participants

#cut out the 'R' in some participant IDs (indicates repeat for controls)
terminal$scan_id <- gsub("[R]", "", terminal$scan_id)

#make a 'group' column to capture case vs. control information
terminal$group <- stri_sub(terminal$scan_id, 2, 2) #note: 1 or 2 is patient, 6 is control

#for clarity, change values in 'group' column to labels for clarity
terminal$group[terminal$group == 1] <- "patient"
terminal$group[terminal$group == 2] <- "patient"
terminal$group[terminal$group == 6] <- "control"

#make a variable that combines unique ID and session number
terminal$id_session <- paste(terminal$scan_id, '_', terminal$session, sep='')

#write csv
write.csv(terminal, '../generated_csvs/terminal_clean_2018-01-25.csv', row.names=FALSE)

#cleanup
rm(terminal)
```

## Checking XNAT

This script pulls in and cleans up the naming of STOPPD scans as they exist in XNAT. At earlier stages, this script helped us identify naming errors in XNAT (all have since been fixed).

**Purpose:** The contents of XNAT will, in other scripts, be checked against (1) the scans we have in our file system, to ensure that there are no discrepancies between these databases, and also against (2) our subject inclusion list.



```r
#import spreadsheets (exported from XNAT)
xnat_camh <- read.csv('../data/xnat_records/xnat_cmh_2018-01-25.csv')
xnat_nki <- read.csv('../data/xnat_records/xnat_nki_2018-01-25.csv')
xnat_pitt <- read.csv('../data/xnat_records/xnat_pmc_2018-01-25.csv')
xnat_umass <- read.csv('../data/xnat_records/xnat_umas_2018-01-25.csv')

#combine XNAT spreadsheets, take only ID and date columns
xnat <- Reduce(function(x, y) merge(x, y, all=TRUE), list(xnat_camh, xnat_nki, xnat_pitt, xnat_umass))
xnat <- xnat[c('MR.ID', 'Date') ]

#cleanup
rm (xnat_camh, xnat_nki, xnat_pitt, xnat_umass)

#import spreadsheet of data in file system (made in script 01_STOPPD_terminal)
terminal <- read.csv('../generated_csvs/terminal_clean_2018-01-25.csv')
```


```r
#remove all CAMH scans with '00' as timepoint (NOTE: '00' this is a consequence of creative naming to account for MRS scans)
xnat$timepoint <- str_sub(xnat$MR.ID, start= -2) #make column with timepoint data
xnat <- xnat[-grep('00', xnat$timepoint),] #remove those with 00

#cut out timepoint info from subject ID string - now meaningless - and remove timepoint column
xnat$MR.ID  <- str_sub(xnat$MR.ID, 1, -4)
xnat <- xnat[, -grep('timepoint', colnames(xnat))]

#cut out study and site info from subject ID string - not needed
xnat$MR.ID <- substring(xnat$MR.ID, 12)

#make a new column for session component of ID
xnat$session <- str_sub(xnat$MR.ID, -2)
table(xnat$session)
```

```
## 
##  00  01  02  03 
##  17 222  77   7
```

```r
#cut out session from subject ID string - not needed
xnat$MR.ID <- str_sub(xnat$MR.ID, 1, -4)

#make a new column that captures alphabetic component of ID ('R')
xnat$contains_R <- grepl('R', xnat$MR.ID, fixed=TRUE) 

#cut out the 'R' in some participant IDs (indicates repeat for controls)
xnat$MR.ID <- gsub("[R]", "", xnat$MR.ID)

#make a variable that combines unique ID and session number
xnat$id_session <- paste(xnat$MR.ID, '_', xnat$session, sep='')

#check for consistency between file system and XNAT 
X <- terminal$id_session %in% xnat$id_session 
  which(X == FALSE) #identical
```

```
## integer(0)
```

```r
Y <- xnat$id_session %in% terminal$id_session 
  which(Y == FALSE) #identical
```

```
## integer(0)
```

```r
#write csv
write.csv(xnat, '../generated_csvs/xnat_clean_2018-01-25.csv', row.names=FALSE)

#cleanup
rm(terminal, xnat)
```


