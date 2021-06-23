###############################################################################

# Prepare the clinical data

###############################################################################

library(stringr)

###############################################################################

# Load the data
data <- read.csv("~/School/WUR/SSB-80336 - Thesis/Provided Data/csv_INFECT_export_180604.csv",
                 sep = ';', dec = ',', check.names = FALSE)
rownames(data) <- data$`Patient ID`

###############################################################################

# Remove columns with too many levels
for (word in c('date','Date','Description','Hospital','immunodef_specified',
               'replaced','surgical_number')) { #'Other','Microbiological.species.found.in.',
  to_remove <- grepl(word, colnames(data))
  data <- data[,-which(to_remove==TRUE)]
}

# Replace unknown by NA
for (column in which(grepl('unknown', colnames(data))==TRUE)) {
  splitted<-strsplit(colnames(data)[column], split = '[=]|[(]|[)]|[,]|[ ]')[[1]]
  # Find which number represents 'unknown' (usually 2, but not always)
  for (i in 1:length(splitted)){
    if (splitted[i] == 'unknown') {
      number <- splitted[i-1]
    }
  }
  # If the column contains this number, change it to NA
  for (sample in 1:nrow(data)){
    if (data[sample,column] == as.integer(number)) {
      data[sample,column] <- NA
    }
  }
}
# Not all colnames have the word 'unknown' in it
for (column in which(grepl('0=no', colnames(data))==TRUE)) {
  # if (any(is.na(data[,column]))==FALSE) {
  if (length(unique(data[,column]))>2) {
    for (sample in 1:nrow(data)){
      # If it is 2, change it to NA
      if (is.na(data[sample,column]) == FALSE) {
        if (data[sample,column] == '2') {
        data[sample,column] <- NA
        }
      }
    }
  }
}
# Not all colnames have the word 'unknown' in it
for (column in which(grepl('_amputated', colnames(data))==TRUE)) {
  # print(colnames(data)[column])
  if (any(is.na(data[,column]))==FALSE) {
    if (length(unique(data[,column]))>2) {
      for (sample in 1:nrow(data)){
        # If it is 2, change it to NA
        if (data[sample,column] == '2') {
          data[sample,column] <- NA
        }
      }
    }
  }
}
# When there is no entry, change it to NA
for (column in 1:ncol(data)) {
  levels(data[,column])[levels(data[,column]) == ""]<-NA
  # New R version:
  empty <- which(data[, column]=="")
  data[empty, column] <- NA
}

# Remove variables with only 1 level
remove_level <- NULL
for (column in 1:ncol(data)) {
  level <- length(unique(data[,column]))
  if (level==1) {
    remove_level <- c(remove_level, column)
  }
}
data <- data[,-remove_level]

###############################################################################

# Change some names to improve interpretation
colnames(data)[which(colnames(data)
                     =='2=Copenhagen, 3=Stockholm, 4=Karlskrona, 5=Gothenburg, 6=Bergen')]<-
  'City (2=Copenhagen, 3=Stockholm, 4=Karlskrona, 5=Gothenburg, 6=Bergen)'
colnames(data)[which(colnames(data)
                     =='Comorbidity, 1=yes, 0=no')]<-'Comorbidity (1=yes, 0=no)'
colnames(data)[which(colnames(data)
                     =='Septic shock Baseline, 1=yes, 0=no')]<-'Septic shock Baseline (1=yes, 0=no)'
colnames(data)[which(colnames(data)
                     == "BMI WHO classification: 1=underweight, 2=normal, 3=pre-obesity, 4=obesity class I, 5=obesity class II, 6=obesity class III")] <- 
  "BMI WHO classification (1=underweight, 2=normal, 3=pre-obesity, 4=obesity class I, 5=obesity class II, 6=obesity class III)"
colnames(data)[which(colnames(data) =='Mechanical ventilation Baseline')]<-
  'Mechanical ventilation Baseline (1=yes, 0=no)'
colnames(data)[which(colnames(data) =='Need For reconstructive Surgery ((if not in hospital day 90)')]<-
  'Need For reconstructive Surgery (1=yes, 0=no) (if not in hospital day 90)'
colnames(data)[which(colnames(data) =='Vancomycin(0=no, 1=yes) Day 1')]<-
  'Vancomycin (0=no, 1=yes) Day 1'

# colnames(data)[which(colnames(data) =='toe_amputated')]<-
#   'toe_amputated (0=no, 1=yes)'
# colnames(data)[which(colnames(data) =='foot_amputated')]<-
#   'foot_amputated (0=no, 1=yes)'
# colnames(data)[which(colnames(data) =='hand_amputated')]<-
#   'hand_amputated (0=no, 1=yes)'
# colnames(data)[which(colnames(data) =='penis_amputated')]<-
#   'penis_amputated (0=no, 1=yes)'
# colnames(data)[which(colnames(data) =='finger_amputated')]<-
#   'finger_amputated (0=no, 1=yes)'
# colnames(data)[which(colnames(data) =='lower_arm_amputated')]<-
#   'lower_arm_amputated (0=no, 1=yes)'
# colnames(data)[which(colnames(data) =='upper_arm_amputated')]<-
#   'upper_arm_amputated (0=no, 1=yes)'
# colnames(data)[which(colnames(data) =='lower_leg_amputated')]<-
#   'lower_leg_amputated (0=no, 1=yes)'
# colnames(data)[which(colnames(data) =='upper_leg_amputated')]<-
#   'upper_leg_amputated (0=no, 1=yes)'


loc <- which(colnames(data)=='Amputation of Limb (Time: from diagnosis to ICU day 7)')
data[which(data[,loc] == 'yes'),loc] <- 1
data[which(data[,loc] == 'no'),loc] <- 0
class(data[,loc]) <- 'integer'

# colnames(data)[which(colnames(data) =='Amputation of Limb (Time: from diagnosis to ICU day 7)')]<-
  # 'Amputation of Limb (Time: from diagnosis to ICU day 7) (0=no, 1=yes)'

###############################################################################

# Make sure the spelling of the different levels is the same

simpleCap <- function(x) {
  s <- strsplit(x, " ")[[1]]
  paste(toupper(substring(s, 1, 1)), substring(s, 2),
        sep = "", collapse = " ")
}
for (column in c(46,274:280)) {
  for (sample in 1:nrow(data)) {
    data[sample, column] <- simpleCap(data[sample, column])
  }
  
}

for (antibiotic in c(46,275:280)) {
  data[which(data[,antibiotic]=='Abboticin (Erythromycin)'),antibiotic] <- 'Erythromycin'
  data[which(data[,antibiotic]=='Erytromycin'),antibiotic] <- 'Erythromycin'
  data[which(data[,antibiotic]=='Abboticin'),antibiotic] <- 'Erythromycin'
  data[which(data[,antibiotic]=='Erythromycin (Abboticin)'),antibiotic] <- 'Erythromycin'
  data[which(data[,antibiotic]=='Abboticin, Fungizone'),antibiotic] <- 'Erythromycin, Fungizone'
  data[which(data[,antibiotic]=='Diflucan, Abboticin'),antibiotic] <- 'Erythromycin, Abboticin'
  data[which(data[,antibiotic]=='Colimycin, Abboticin, Fungizone'),antibiotic] <- 'Colimycin, Erythromycin, Fungizone'
  data[which(data[,antibiotic]=='Erythromycin, Abboticin'),antibiotic] <- 'Erythromycin'
  
  data[which(data[,antibiotic]=='Bioclavid (amoxicillin+clavulansyre)'),antibiotic] <- 'Bioclavid'
  
  data[which(data[,antibiotic]=='Flukloxacillin'),antibiotic] <- 'Flucoxacillin'
  data[which(data[,antibiotic]=='Flukloxacillin'),antibiotic] <- 'Flucoxacillin'
  data[which(data[,antibiotic]=='Flukloxacillin'),antibiotic] <- 'Flucoxacillin'
  
  data[which(data[,antibiotic]=='Flukonazol'),antibiotic] <- 'Fluconazol'
  data[which(data[,antibiotic]=='Fluconazole'),antibiotic] <- 'Fluconazol'
  data[which(data[,antibiotic]=='Flucanazol'),antibiotic] <- 'Fluconazol'
    
  data[which(data[,antibiotic]=='Tienam(IMIPENEM)'),antibiotic] <- 'Tienam'
  data[which(data[,antibiotic]=='Tienam (imipenem)'),antibiotic] <- 'Tienam'
  data[which(data[,antibiotic]=='Tienam = Imipenem Cilastin'),antibiotic] <- 'Tienam'
  data[which(data[,antibiotic]=='Cilastatin/Imipenem (Tienam)'),antibiotic] <- 'Tienam'
  data[which(data[,antibiotic]=='Tienam (Imipenem + Cilastattin)'),antibiotic] <- 'Tienam'

  data[which(data[,antibiotic]=='Clarithromycin.'),antibiotic] <- 'Clarithromycin'
  
  data[which(data[,antibiotic]=='Colimycin Lundbeck'),antibiotic] <- 'Colimycin'

  data[which(data[,antibiotic]=='Bensylpenicillin'),antibiotic] <- 'Benzylpenicillin'
  data[which(data[,antibiotic]=='Besylpenicillin'),antibiotic] <- 'Benzylpenicillin'
  
  data[which(data[,antibiotic]=='Epivir, Ziagen Kaletra'),antibiotic] <- 'Epivir, Ziagen, Kaletra'
  data[which(data[,antibiotic]=='Epivir, Ziagen,kaletra'),antibiotic] <- 'Epivir, Ziagen, Kaletra'
  data[which(data[,antibiotic]=='Epivir,ziagen, Kaletra'),antibiotic] <- 'Epivir, Ziagen, Kaletra'
}

# sort(unique(c(unique(data[,46]),unique(data[,275]),unique(data[,276]),unique(data[,277]),unique(data[,278]),unique(data[,279]),unique(data[,280]))))

###############################################################################

# Remove some variables that are not of interest

for (word in c('HOME','Given Before ICU','korrektion','Day 2','Day 3','Day 4',
               'Day 5','Day 6','Day 7','Nr. 2','Nr. 3','Nr. 4','Nr. 5','Nr. 6',
               'Nr. 7','Nr. of Surgical Procedures','Data Collected Through',
               'Tissue Biopsy For Histology','Information Collected',
               'Discharged To','Certainty of tissue sample',
               'Bodypart affected at arrival - conclusion')) { #'Other','Microbiological.species.found.in.',
  to_remove <- grepl(word, colnames(data))
  data <- data[,-which(to_remove==TRUE)]
}

###############################################################################

# Misc.

# Use the rownames
data$Row.names <- rownames(data)

# Get rid of quotation marks
for (name in colnames(data)) {
  colnames(data)[which(colnames(data)==name)]<-str_remove_all(name, "'")
}

###############################################################################

# Combine the amputation into one column

amputation_columns <- which(grepl('putat', colnames(data)))

data$amputation <- rep(NA, nrow(data))

for (column in amputation_columns) {
  for (sample in 1:nrow(data)) {
    if (is.na(data[sample,column])==FALSE) {
      if (data[sample,column]==1) {
        data$amputation[sample] <- 1
        }
      if (data[sample,column]==0) {
        if (is.na(data$amputation[sample])==TRUE) {
          data$amputation[sample] <- 0
        }
      }
    }
  }
}

data <- data[,-amputation_columns]

# Combine the temperature columns into one column

temperature_columns <- which(grepl('emperature', colnames(data)))

data$`Temperature >38 Degrees Celcius (0=no, 1=yes) Day 1` <- rep(NA,nrow(data))

for (column in temperature_columns) {
  for (sample in 1:nrow(data)) {
    if (is.na(data[sample,column])==FALSE) {
      if (data[sample,column]==1) {
        data$`Temperature >38 Degrees Celcius (0=no, 1=yes) Day 1`[sample] <- 1
      }
      if (data[sample,column]==0) {
        if (is.na(data$`Temperature >38 Degrees Celcius (0=no, 1=yes) Day 1`[sample])==TRUE) {
          data$`Temperature >38 Degrees Celcius (0=no, 1=yes) Day 1`[sample] <- 0
        }
      }
    }
  }
}

data <- data[,-temperature_columns]

data$`Temperature >38 Degrees Celcius (0=no, 1=yes) Day 1`

###############################################################################

# Save the data
setwd("~/School/WUR/SSB-80336 - Thesis/Provided data/Merged/")
write.csv(data, 'Prepped_INFECT.csv')
