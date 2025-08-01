# Process_NHANES_2003-2006_7164_semi-raw_accelerometer_data
 My workflow on how to derive meaningful parameters (Physical Activity, part of Circadian Rhythm etc...) from NHANES 2003-2006 SEMI-RAW Accelerometer data

The NHANES 2003 and 2005 cycles released semi-raw accelerometer data nearly two decades ago:
- [2003 Cycle](https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2003/DataFiles/PAXRAW_C.htm)
- [2005 Cycle](https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2005/DataFiles/PAXRAW_D.htm)

NHANES is an excellent dataset because it is nationally representative (as long as survey weights are handled correctly 😄), free, and open access. Using these two cycles, you can generate some impactful papers ~~(Association between A and B among C / Prevalence and correlates of A among B / Trends in A among B)~~!

However, the semi-raw accelerometer data must be processed to yield meaningful results. While NCI has already released an [algorithm](https://epi.grants.cancer.gov/nhanes-pam/create.html) (why there is some hardcoding in the SAS code??), there's also a [great R package](https://github.com/vandomed/nhanesaccel) for processing the 2003–2006 data (I really like the C part, fast, simple, and straightforward). But my goal is a bit different: I already processed the 2011–2014 cycles, and I want the results to be comparable aka processed using the same pipeline. I prefer GGIR for that. Below is a step-by-step guide to processing these files using the same pipeline as the 2011–2014 cycles.

---

## Steps for Processing NHANES Semi-Raw Accelerometer Data

### Step 0: Download All Files

Just download [PAXRAW_C](https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2003/DataFiles/PAXRAW_C.zip) and [PAXRAW_D](https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2005/DataFiles/PAXRAW_D.zip).

### Step 1: Go Through NCI Calibration

I suggest going through the NCI calibration process using [create.pam_perminute.sas](https://epi.grants.cancer.gov/nhanes-pam/create.pam_perminute.sas). This file is quite interesting since NCI appears to have manually reviewed all 2003–2004 data and performed some QC (see lines 122–156). I also recommend excluding records with PAXSTAT = 2.

However, this file is only available for 2003–2004. You can easily modify it to work for 2005–2006, just skip the manual QC part. After this step, you’ll have two giant SAS files:
pam_perminute_c.sas7bdat and pam_perminute_d.sas7bdat, totaling about 8.43 GB.


### Step 2: Import DEMO and Merge with PAM
Next, merge age information from the DEMO files into the PAM dataset. NCI used age-specific cutoffs (see [create.pam_perday.sas](https://epi.grants.cancer.gov/nhanes-pam/create.pam_perday.sas), lines 206–250). I’m lazy, so I’m using the tidyverse.

```bash
pam_perminute_C <- read_sas("pam_perminute_c.sas7bdat")
pam_perminute_D <- read_sas("pam_perminute_d.sas7bdat")
demo_C <- read_xpt("DEMO_C.XPT")
demo_D <- read_xpt("DEMO_D.XPT")
demo <- bind_rows(demo_C, demo_D) %>% dplyr::select(SEQN, RIDAGEYR)
pam_perminute <- bind_rows(pam_perminute_C, pam_perminute_D %>% dplyr::select(-PAXSTEP)) 
pam_perminute <- inner_join(pam_perminute, demo, by = "SEQN")
```

### Step 3: Generate Age-Specific CSV Files

As mentioned, there are age-specific cutoffs, so you'll need to process the data by age group. Let's focus on age ≥ 18 first, since that's the majority.

Then, format the result into ActiGraph GT3X-style CSV files, using VM = count, and match the data format used in the 2011–2014 cycles. This saves time and ensures comparability.

```bash
# Working dataset
pam_perminuteWorking <- filter(pam_perminute, RIDAGEYR>=18) %>% dplyr::select(-RIDAGEYR)
allSEQN <- unique(pam_perminuteWorking$SEQN)

# Map PAXDAY (1=Sunday…7=Saturday) to a date , silimar to 2011-2014 format (1st week of 2000-01-01)
first_week <- as.Date("2000-01-01") + 0:6
# POSIXlt$wday: 0=Sunday…6=Saturday, so add 1 to align with PAXDAY
day_codes  <- as.POSIXlt(first_week)$wday + 1
date_map   <- setNames(first_week, day_codes)

# 3) loop
for (seq in allSEQN) {
  
  # subset
  temp <- filter(pam_perminuteWorking, SEQN == seq)
  
  # transform body exactly as before
  out <- temp %>%
    mutate(
      Date = date_map[as.character(PAXDAY)],
      TimeStamp = as.POSIXct(
        paste0(Date, " ", sprintf("%02d:%02d:00", PAXHOUR, PAXMINUT)),
        tz = "UTC"
      ),
      TimeStamp = format(TimeStamp, "%Y-%m-%dT%H:%M:%SZ"),
      vm    = PAXINTEN,
      axis1 = vm, axis2 = 0, axis3 = 0 # Those are dummy values , GGIR won't use XYZ axis acc to do any of the calculation
    ) %>%
    select(TimeStamp, axis1, axis2, axis3, vm) # VM is derived count data
  
  # rebuild acc header using this participant’s first row
  hdr <- c(
    "------------ Data Table File Created By Actigraph Link ActiLife v6.11.9 date format dd/MM/yyyy Filter Normal -----------,,,,,",
    "Serial Number: TAS1D48140206,,,,,", # Just a random SN
    paste0(
      "Start Time ",
      format(
        as.POSIXct(
          paste0(
            date_map[ as.character(temp$PAXDAY[1]) ], " ",
            sprintf("%02d:%02d:00", temp$PAXHOUR[1], temp$PAXMINUT[1])
          ),
          tz = "UTC"
        ),
        "%H:%M:%S"
      ),
      ",,,,,"
    ),
    paste0(
      "Start Date ",
      format(
        as.Date(date_map[ as.character(temp$PAXDAY[1]) ]),
        "%m-%d-%Y"
      ),
      ",,,,,"
    ),
    "Epoch Period (hh:mm:ss) 00:00:60,,,,,",
    "Download Time 15:29:44,,,,,",
    "Download Date 09-19-2017,,,,,",
    "Current Memory Address: 0,,,,,",
    "Current Battery Voltage: 3.85     Mode = 13,,,,,",
    "--------------------------------------------------,,,,,"
  )
  
  # file path
  file_path <- file.path("/NHANESGGIR/NHANES03040506Adult/", paste0(seq, ".csv"))

  # write header + body
  writeLines(hdr, file_path)
  write.table(
    out, file_path,
    sep        = ",",
    row.names  = FALSE,
    col.names  = TRUE,
    quote      = FALSE,
    append     = TRUE
  )
}
```
### Step 4: Process the CSV Files

The core idea is to use GGIR to process the 2003–2006 accelerometer data without calibration and imputation. The 2003–2006 data were collected using the ActiGraph 7164, a single-axis accelerometer with 60-second epoch count data. Participants were instructed to remove the device while showering, swimming, or sleeping ([see protocol](https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2003/DataFiles/PAXRAW_C.htm)). 

As a result:
1. Most established sleep detection and step count algorithms (e.g., Verisense step count algorithm) do not work on this data.
2. You won’t get reliable circadian rhythm estimates because the wear time is generally insufficient. Yes I know there are some papers report RA, M10, and L5. However, L5 is most likely to coincide with the non-wear period (during sleep), resulting in extremely high RA values, often near 1 in 90% of cases which doesn’t make sense at all.
3. In ~10% of participants who did not follow protocol, it may be possible to estimate RA and other circadian rhythm parameters.

I wrote code to match the NCI algorithm as closely as possible. Key features:
1. Calibration and imputation are turned off; these steps are not feasible with single-axis count-based data.
2. Windowsizes are set to c(60, 3600, 3600). This is matches the NCI algorithm. Here's an explanation:
```
window1 = 60: Matches the native 60-second epoch of the ActiGraph 7164
window2 = 3600: 1-hour windows for non-wear detection
window3 = 3600: 1-hour assessment windows with overlapping detection
```
4. Uses GGIR’s non-wear detection (NotWorn algorithm) to identify sleep. I tried L5+/-12 but it gave me some pretty wild results.
5. A 10-hour minimum wear time defines a valid day (same as NCI).
6. Uses a 100-count non-wear threshold (though I’m not sure it works well).
7. Applies NCI-defined thresholds for sedentary, light and MVPA.
8. Uses already-derived ActiGraph 7164 counts (VM) as the acc metric.
9. You have to use NotWorn for both HASPT.algo and HASIB.algo.
10. Any columns containing NeishabouriCount_x|y|z are not usable because they contain dummy values.

My GGIR parameter:

```bash
GGIR(
  idloc                    = 6, # SEQN before dot
  datadir                  = "/NHANESGGIR/NHANES03040506Adult",
  outputdir                = "/NHANESGGIR/ResultNHANES03040506Adult",
  do.cal                   = FALSE,
  dataFormat               = "actigraph_csv",
  extEpochData_timeformat  = "%m-%d-%Y %H:%M:%S",
  windowsizes              = c(60, 3600, 3600),
  do.neishabouricounts     = FALSE,
  acc.metric               = "NeishabouriCount_vm",
  studyname                = "NHANES03040506",
  HASPT.algo               = "NotWorn",
  HASIB.algo               = "NotWorn",
  do.imp                   = FALSE,
  nonwear_range_threshold  = 100, 
  HASPT.ignore.invalid     = NA,
  ignorenonwear            = FALSE,
  threshold.lig            = 100, # Below 100 cpm is sedentary
  threshold.mod            = 2020, # Change this to age specific cut off
  threshold.vig            = 5999, # Change this to age specific cut off
  includedaycrit           = 10, # NCI consider >=10 hr per day is valid
  includenightcrit         = 10)
```
I compared MVPA (total duration, no need to consider bouts at the 60-second epoch level, as the 60s resolution already smooths out random acceleration) and sedentary behavior time against the NCI SAS algorithm. Among participants aged ≥18, the correlations were 97.12% for MVPA and 85.43% for sedentary time (raw correlation without survey weights). I would consider this safe to use.

# If you have any questions or would like access to the processed data (it’s way too large to upload to GitHub), feel free to contact me on LinkedIn!
