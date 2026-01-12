#example EDA, PCA, etc. code using several postgraduate datasets
#W.J. Matthaeus 7/10/2026
#this includes some data quality control on the way in
#and interpretation

#----packages-------------------------------------------------------------------
#this is to load in required packages
pkgs <- c("readr", "ggplot2", "dplyr", "readxl",
          "tibble", "tidyr", "stringr",
          "PerformanceAnalytics", "psych",
          "FactoMineR", "factoextra","rstatix",
          "nsprcomp","compositions","FactoInvestigate")
vapply(pkgs, library, logical(1), character.only = TRUE, logical.return = TRUE)

#----wd + input swtiches--------------------------------------------------------
#change the working directory. update the string "/Users.." to contain the 
#path to the directory containing this script and your data file
#setwd("/Users/willmatthaeus/Dropbox/didactics/didactic_xrf/")
'%!in%' <- function(x,y)!('%in%'(x,y))
#only change this here, then the rest of the code should work for *your* data
# input_switch <-"Ant"
# input_switch <- "Sid"
input_switch <- "Cat_Penn"
#input_switch <- "Cat_Modern"

#change the working directory. update the string "/Users.." to contain the 
#path to the directory containing this script and your data file
if(input_switch=="Cat_Penn"){
  setwd("C:/Users/catar/OneDrive/Ambiente de Trabalho/Pennsylvanina/Pennsylvanian XRF data")
}

if(input_switch=="Cat_Modern"){
  setwd("C:/Users/catar/OneDrive/Ambiente de Trabalho/modern")}

#----read in data---------------------------------------------------------------
#read in data
if(input_switch=="Sid"){
#read in the data
xrf <- read_excel(path = "data/Mastersheet.xlsx")
#find the group ids from the sample name and make a new column
xrf$group<-NA
xrf[grep("CSD",xrf$SAMPLE),]$group <- "new"
xrf[is.na(xrf$group),]$group <- "old"
xrf$group <- as.factor(xrf$group)
#
elements<-xrf[seq(2,76,2)]
element_names <- colnames(elements)
errors<-xrf[seq(3,77,2)]
error_names <- colnames(errors)

#grab some more info from the sample names
xrf<-xrf %>% select(all_of(c("SAMPLE", element_names))) %>%
              mutate(splitted = str_split_fixed(SAMPLE,"_",n = 2),
                     site_time = splitted[,1],
                     time = str_match(splitted[,1],"[A-Z]+"),
                     site = str_match(splitted[,1],"[0-9]+"),
                     species = splitted[,2],
                     )

# %>%
#               rename(time = "time[, 1]")
# 
# 
# #prototypes
# xrf$t
# str_split_fixed(xrf$SAMPLE[2],"_",n = 2)[,2]
#separate out the data columns and the error columns

}

if(input_switch=="Ant"){
  
  #read in the data
  xrf <- read_excel(path = "data/STANAST_XRF_GINK.xlsx")
  #find the group ids from the sample name and make a new column
  which(xrf$SAMPLE=="b4-51258-F")
  xrf[57,]$SAMPLE<-"b4-51258-F_2"
  which(xrf$SAMPLE=="b4-51408-F")
  xrf[77,]$SAMPLE<-"b4-51408-F_2"
  
  #grouping variables
  xrf$locality<-NA
  xrf[grep("AST",xrf$SAMPLE),]$locality <- "Astartekløft"
  xrf[is.na(xrf$locality),]$locality <- "South Tankrediakløft"
  xrf$locality <- as.factor(xrf$locality)
  which(is.na(xrf$locality))
  
  xrf$taxon <- factor(xrf$TAXON)
  
  #separate out the data columns and the error columns
  elements<-xrf[3:32]
  element_names <- colnames(elements)
  
}

##----Cat Penn------------------------------------------------------------------          
if(input_switch=="Cat_Penn"){
  #read in the data
  xrf <- read_excel(path = "C:/Users/catar/OneDrive/Ambiente de Trabalho/Pennsylvanina/Pennsylvanian XRF data/Carboniferous XRF.xlsx",sheet = 2)
 
  #find the group ids from the sample name and make a new column
  sdarm2 <- xrf %>% filter(SAMPLE%in%c("sdarm2","sdarm"))

###-----QC------------------------
  #compare measured+-error v standard+-20%
  
  #clean up data before testing
  #separate out the data columns and the error columns
  element.colnums<-seq(17,85,2)
  #this next line turns the column into numeric data
  #if any values can’t be converted, they become NA
  sdarm2[element.colnums] <- sapply(sdarm2[element.colnums],as.numeric)
  elements<-sdarm2[element.colnums]
  element_names <- colnames(elements)
  
  sdarm2 <-  sdarm2 %>% select(
    -`Reading No`,
    -Time,
    -Type,
    -Duration,
    -Units,
    -Sequence,
    -Flags,
    -LOCATION,
    -`Leaf/Axis`,
    -`Structure #`,
    -taxon,
    -INSPECTOR,
    -MISC,
    -NOTE,
    -`User Login`
      )
  
  #Standardized values
  standards <- c(
    Si = 343411, Al = 65995, K = 41507, Fe = 18395, Ca = 6003,
    Mg = 2955, Ti = 1798, Ba = 990, Pb = 808, Mn = 1038,
    Zn = 772, p = 345, Zr = 259, Cu = 239, sr = 144,
    cr = 51, sb = 111, v = 25.2, rb = 149, as = 80, mo = 13.1
  )
  
  # Identify element + error columns
  element_cols <- names(sdarm2)[
    names(sdarm2) %in% names(standards) |
      names(sdarm2) %in% paste0(names(standards), " Error")
  ]
  
  # Force element + error columns to numeric
  # (character artifacts become NA)
  sdarm2_clean <- sdarm2 %>%
    mutate(
      across(
        all_of(element_cols),
        ~ readr::parse_number(as.character(.))
      )
    )
  
  
  # Reshape and evaluate QC
  qc_results <- sdarm2_clean %>%
    mutate(row_id = row_number()) %>%
    pivot_longer(
      cols = all_of(element_cols),
      names_to = "variable",
      values_to = "value"
    ) %>%
    mutate(
      Element = str_remove(variable, " Error$"),
      type = if_else(str_detect(variable, " Error$"), "error", "value")
    ) %>%
    select(-variable) %>%
    pivot_wider(
      names_from = type,
      values_from = value
    ) %>%
    mutate(
      standard = standards[Element],
      measured_min = value - error,
      measured_max = value + error,
      acceptable_min = standard * 0.8,
      acceptable_max = standard * 1.2,
      QC_pass = measured_max >= acceptable_min &
        measured_min <= acceptable_max
    )
  
  #if QC_pass=NA then elemeent fails

##-----data finagleing-----------------------------

  #removes the standards from the xrf dataframe
  xrf<- xrf %>% filter(SAMPLE%!in%c("sdarm2","sdarm"))

  #remove non vegetative taxa from df
  xrf <- filter(xrf, taxon!=c("Lepidostrobophyllum"))
  xrf <- filter(xrf, taxon!=c("?"))
          
  #prints every unique name in SAMPLE
  xrf %>% select(SAMPLE)%>%unique%>%print(n=126)
          
   #creates two dataframes
    #one with the names that are already formatted correctly
  goodnames_xrf <- xrf[grep(x = xrf$SAMPLE,pattern = "_"),]
    #and one with Bad Formatting
  badnames_xrf <- xrf[-grep(x = xrf$SAMPLE,pattern = "_"),]

          
  goodnames_xrf <- goodnames_xrf %>% 
    #creates new column "name_split" and splits SAMPLE in 3 strings divided by "-"
      #name_split is a 3 column matrix inside my df      
    mutate(name_split = str_split_fixed(SAMPLE,'-',n=3)) %>%
      #create a new column "name" from the third element of name split
    mutate(name = name_split[,3]) %>%
      #changes "name_split" to be divided by "_"
    mutate(name_split = str_split_fixed(name,'_',n=3))%>%
    #create column "core" from the first element of "name_split"
    #              "side"          second
    #              "spot"          third
    mutate(core= name_split[,1], side = name_split[,2], spot = name_split[,3])%>%
    #removes name_split and name
    select(-name_split, -name)

  #in the column "side" - change all the entries matching "front" to be upper case
  goodnames_xrf$side[which(goodnames_xrf$side == 'front')]<-"FRONT"
  #fill the column substrate w NAs
  goodnames_xrf$substrate <- NA
  #fill the column substrate w Matrix if the spot column has an M or an m
  goodnames_xrf$substrate[grep(x = goodnames_xrf$spot, pattern = '[mM]')] <- 'Matrix'
  #fill the column substrate w Fossil if the spot column has an F or an f
  goodnames_xrf$substrate[grep(x = goodnames_xrf$spot, pattern = '[fF]')] <- 'Fossil'
  
#check for all the different spellings in the LOCATION column
  unique(goodnames_xrf$LOCATION)
  # [1] "womac"              "murphysboro"        "new haven"          "seline,springfield"
  # [5] "chapel"             "Murphysboro"        "Chapel"
          
  #find replace names in the LOCATION column´
    #== means that they're exactly the same
  goodnames_xrf<-goodnames_xrf%>%
  mutate(
    location = case_when(
    LOCATION == "seline,springfield" ~ "Seline",
    LOCATION == "womac"  ~ "Womac",
    LOCATION == "murphysboro"  ~ "Murphysboro",
    LOCATION == "new haven"  ~ "New haven",
    LOCATION == "chapel"  ~ "Chapel",
    .default = as.character(LOCATION)
    )
  )
            
  #
  # unique(goodnames_xrf$core)
  # unique(goodnames_xrf$side)
  # unique(goodnames_xrf$spot)
  # unique(goodnames_xrf$substrate)
  # unique(goodnames_xrf$location)

  #
  badnames_xrf <- badnames_xrf %>% 
    mutate(name_split = str_split_fixed(SAMPLE,'-',n=2)) %>% 
    mutate(core = name_split[,1], spot=name_split[,2])%>% 
    select(-name_split)
  #
  badnames_xrf$substrate <- NA
  badnames_xrf$substrate[grep(x = badnames_xrf$spot, pattern = '[mM]')] <- 'Matrix'
  badnames_xrf$substrate[grep(x = badnames_xrf$spot, pattern = '[fF]')] <- 'Fossil'
  badnames_xrf$substrate[which(is.na(badnames_xrf$substrate))] <-'Fossil'
  #
  badnames_xrf$side<-"FRONT"
  #
  # unique(badnames_xrf$LOCATION)
  # [1] "paum mine" "womac"     "Paum mine"
  badnames_xrf<-badnames_xrf%>%
    mutate(
      location = case_when(
        LOCATION == "paum mine" ~ "Paum mine",
        LOCATION == "womac"  ~ "Womac",
        .default = as.character(LOCATION)
      )
    )
  #
  # unique(badnames_xrf$spot) 
  # unique(badnames_xrf$substrate) 
  # unique(badnames_xrf$core)
  # unique(badnames_xrf$location)
  
  #marge
  xrf<-rbind(goodnames_xrf, badnames_xrf)
  
  # colnames(xrf)
  # [1] "Reading No" "Time"       "Type"       "Duration"   "Units"      "Sequence"   "Flags"     
  # [8] "SAMPLE"     "LOCATION"   "INSPECTOR"  "MISC"       "NOTE"       "User Login" "Ba"        
  # [15] "Ba Error"   "Sb"         "Sb Error"   "Sn"         "Sn Error"   "Cd"         "Cd Error"  
  # [22] "Pd"         "Pd Error"   "Ag"         "Ag Error"   "Bal"        "Bal Error"  "Mo"        
  # [29] "Mo Error"   "Nb"         "Nb Error"   "Zr"         "Zr Error"   "Sr"         "Sr Error"  
  # [36] "Rb"         "Rb Error"   "Bi"         "Bi Error"   "As"         "As Error"   "Se"        
  # [43] "Se Error"   "Au"         "Au Error"   "Pb"         "Pb Error"   "W"          "W Error"   
  # [50] "Zn"         "Zn Error"   "Cu"         "Cu Error"   "Ni"         "Ni Error"   "Co"        
  # [57] "Co Error"   "Fe"         "Fe Error"   "Mn"         "Mn Error"   "Cr"         "Cr Error"  
  # [64] "V"          "V Error"    "Ti"         "Ti Error"   "Ca"         "Ca Error"   "K"         
  # [71] "K Error"    "Al"         "Al Error"   "P"          "P Error"    "Si"         "Si Error"  
  # [78] "Cl"         "Cl Error"   "S"          "S Error"    "Mg"         "Mg Error"   "core"      
  # [85] "side"       "spot"       "substrate"  "location" 
  
  #separate out the data columns and the error columns
  element.colnums<-seq(17,85,2)
  #this next line turns the column into numeric data
  #if any values can’t be converted, they become NA
  xrf[element.colnums] <- sapply(xrf[element.colnums],as.numeric)
  elements<-xrf[element.colnums]
  element_names <- colnames(elements)

          
xrf <- xrf %>%
    mutate(LOCATION =
             factor(LOCATION, levels = 
                      c("Murphysboro",
                        "Paum Mine",
                        "Springfield",
                        "Chapel",
                        "Womac",
                        "New Haven"
                        )))

          
}

##----Cat Modern----------------------------------------------------------------
if(input_switch=="Cat_Modern"){
  #read in the data
  xrf <- read_xlsx(path = "C:/Users/catar/OneDrive/Ambiente de Trabalho/modern/modernXRF.xlsx",sheet = 2)
  #creates a separate table dataframe w the standards 
  standards <- xrf %>% filter(SAMPLE%in%c("sdarm2","sdarm"))
  #removes the standards from the xrf dataframe
  xrf<- xrf %>% filter(SAMPLE%!in%c("sdarm2","sdarm"))
  #prints all the unique values in SAMPLE
  xrf %>% select(SAMPLE)%>%unique%>%print(n=126)

  
  #changes location from ttec to the trinity botanic gardens
  xrf<-xrf%>%
    mutate(
      LOCATION = case_when(
                    #if      #then
        LOCATION == "ttec" ~ "tbg",
        #else: remain the same
        .default = as.character(LOCATION)
      )
    )
  #creates dataframes
    #one that's just the Ti tests
 # Ti_names <- xrf[grep(x = xrf$SAMPLE,pattern = "-Ti"),]
    #one with the names that are formatted good
      #excludes samples that include -t-
  goodnames_xrf <- xrf[-grep(x = xrf$SAMPLE,pattern = "-t-"),]
    #and one for the names that are formatted bad
  badnames_xrf <- xrf[grep(x = xrf$SAMPLE,pattern = "-t-"),]

###----goodnames----------------------------------------------------------------
  goodnames_xrf <- goodnames_xrf %>%
    #creates new column "name_split" and splits SAMPLE in 3 strings divided by "-"
    #name_split is a 3 column matrix inside my df
        #if Mgr-5 fill method column with CF2
    mutate(method = case_when(
      #CF for coffee grinder  
      str_detect(SAMPLE, pattern = "Mgr-5") ~ "CF2",
      #else fill with MP - for mortar and pestle
      .default = "MP")) %>%

    mutate(name_split = str_split_fixed(SAMPLE,'-',n=3)) %>%
    #create column "genus" from the first element of "name_split"
    #              "replicate"      second
    #              "plate"           third
    mutate(genus = name_split[,1], replicate = name_split[,2], plate = name_split[,3])%>%
    #removes name_split and name
    select(-name_split)

#check that your column names are right  
colnames(goodnames_xrf)  

  #fills the plate rows without Ti to have NAs instead? i think
  goodnames_xrf$plate[-grep(x = goodnames_xrf$plate, pattern = "Ti")] <- NA

  #changes "genus" from the codes to the full genus names
  goodnames_xrf<-goodnames_xrf%>%
    mutate(
      genus = case_when(
        #ferns
        str_detect(genus, pattern = "Ore") ~ "Osmunda	regalis",
        str_detect(genus, pattern = "Dic") ~ "Dicksonia? fibrosa?",
        str_detect(genus, pattern = "Dcy") ~ "Dryopteris cycadena",
        str_detect(genus, pattern = "Dsq") ~ "Dicksonia squarrosa",
        str_detect(genus, pattern = "Wfi") ~ "Woodwardia fimbriata",
        str_detect(genus, pattern = "Gdr") ~ "Gymnocarpium dryopteris",
        str_detect(genus, pattern = "Ose") ~ "Onoclea sensibilis",
        str_detect(genus, pattern = "Cfa") ~ "Cyrtomium	falcatum",
        str_detect(genus, pattern = "Efl") ~ "Equisetum	fluviata",
        str_detect(genus, pattern = "Cco") ~ "Cyathea	cooperi",
        #gymnos
        str_detect(genus, pattern = "Gbi")   ~ "Ginkgo biloba",
        str_detect(genus, pattern = "Sse")   ~ "Sequoia	sempervirens",
       # str_detect(genus, pattern = "Sse-2") ~ "Sequoia sempervirens-2",
        str_detect(genus, pattern = "Tta") ~ "Torreya	taxifolia",
      #  str_detect(genus, pattern = "Tta-2") ~ "Torreya taxifolia-2",
       # str_detect(genus, pattern = "Tta-3") ~ "Torreya taxifolia-3",
        str_detect(genus, pattern = "LPE")   ~ "Lepidozamia	peroffskyana",
        str_detect(genus, pattern = "Efr") ~ "Ephedra fragilis",
        #str_detect(genus, pattern = "Efr-2") ~ "Ephedra	fragilis-2",
        str_detect(genus, pattern = "Wno")   ~ "Wollemia nobilis",
        str_detect(genus, pattern = "Tcr") ~ "Taiwania crypotomeroides",
       # str_detect(genus, pattern = "Tcr-2") ~ "Taiwania crypotomeroides-2",
        str_detect(genus, pattern = "Aau")   ~ "Agathis	australis",
        str_detect(genus, pattern = "Ptr")   ~ "Phyllocarpus tricomanoides",
        str_detect(genus, pattern = "Ppa")   ~ "Pinus	patula",       
        #angio
        str_detect(genus, pattern = "Isi")   ~ "Illicium simonsi",
        str_detect(genus, pattern = "Ihe")   ~ "Illicium henryi",
        str_detect(genus, pattern = "Dwi")   ~ "Drimys winteri",
        str_detect(genus, pattern = "Dta")   ~ "Drimys tasmanica?",
        str_detect(genus, pattern = "Pco")   ~ "Pseudowintera	colorata",
        str_detect(genus, pattern = "Pda")   ~ "Phoenix	dactilifera",
        str_detect(genus, pattern = "Shi")   ~ "Sorbus hibernica",
        str_detect(genus, pattern = "Mgr") ~ "Magnolia grandifolia",
        #str_detect(genus, pattern = "Mgr-2") ~ "Magnolia grandifolia-2",
        #str_detect(genus, pattern = "Mgr-3") ~ "Magnolia grandifolia-3",
        #str_detect(genus, pattern = "Mgr-4") ~ "Magnolia grandifolia-4",
        #str_detect(genus, pattern = "Mgr5") ~ "Magnolia grandifolia-5",
        str_detect(genus, pattern = "Mdo")   ~ "Malus	domestica",
        str_detect(genus, pattern = "Aun")   ~ "Arbutus unido",
        str_detect(genus, pattern = "Jni") ~ "Juglans nigra",
      #  str_detect(genus, pattern = "Jni-2") ~ "Juglans nigra-2"
      )
    )
  
##-----------------------------badnames-----------------------------------------
    badnames_xrf<-badnames_xrf %>%
    #create method column, assigns them values
    mutate(method = case_when(
      #CF for coffee grinder
      str_detect(SAMPLE, pattern = "-cf") ~ "CF1",
      #or MP for mortar and pestle     
      str_detect(SAMPLE, pattern = "-mp") ~ "MP")) %>%
    #fixes the names
    mutate(
      SAMPLE = case_when(
        str_detect(SAMPLE, pattern = "t-mp1") ~ "Mgr-1.1",
        str_detect(SAMPLE, pattern = "t-mp2") ~ "Mgr-1.2",
        str_detect(SAMPLE, pattern = "t-mp3") ~ "Mgr-1.3",
        str_detect(SAMPLE, pattern = "t-cf1") ~ "Mgr-1.4",
        str_detect(SAMPLE, pattern = "t-cf2") ~ "Mgr-1.5",
        str_detect(SAMPLE, pattern = "t-cf3") ~ "Mgr-1.6",
        .default = as.character(SAMPLE)
      )
    ) %>%
    #creates column "name_split" and splits SAMPLE in 3 strings divided by "-"
    #name_split is a 3 column matrix inside my df
    mutate(name_split = str_split_fixed(SAMPLE,'-',n=3)) %>%
    #create column "genus" from the first element of "name_split"
    #              "replicate"      second
    #              "plate"           third
    mutate(genus = name_split[,1],
           replicate = name_split[,2],
           plate = name_split[,3])%>%
    #changes genus from the code to the genus name
    mutate(
      genus = case_when(
        str_detect(genus, pattern = "Mgr") ~ "Magnolia	grandifolia",
        .default = genus
      )) %>%
   # mutate(method = case_when(
      #CF for coffee grinder
    #  str_detect(SAMPLE, pattern = "-1.4") ~ "CF1",
     # str_detect(SAMPLE, pattern = "-1.5") ~ "CF1",
      #str_detect(SAMPLE, pattern = "-1.6") ~ "CF1",
      #or MP for mortar and pestle     
  #    str_detect(SAMPLE, pattern = "-1.1") ~ "MP",
   #   str_detect(SAMPLE, pattern = "-1.2") ~ "MP",
    #  str_detect(SAMPLE, pattern = "-1.3") ~ "MP",
     # )) %>%
    #removes name_split
    select(-name_split)
  
  badnames_xrf$plate[-grep(x = badnames_xrf$plate, pattern = "Ti")] <- NA
  
  #check that your column names are right
  colnames(badnames_xrf)
  
  #merge
  xrf<-rbind(goodnames_xrf, badnames_xrf)
  
colnames(xrf)  
  
  #add fern/gymno/angio column
  xrf<-xrf%>%
    mutate(
      group = case_when(
        #ferns
        str_detect(SAMPLE, pattern = "Ore") ~ "Fern",
        str_detect(SAMPLE, pattern = "Dic") ~ "Fern",
        str_detect(SAMPLE, pattern = "Dcy") ~ "Fern",
        str_detect(SAMPLE, pattern = "Dsq") ~ "Fern",
        str_detect(SAMPLE, pattern = "Wfi") ~ "Fern",
        str_detect(SAMPLE, pattern = "Gdr") ~ "Fern",
        str_detect(SAMPLE, pattern = "Ose") ~ "Fern",
        str_detect(SAMPLE, pattern = "Cfa") ~ "Fern",
        str_detect(SAMPLE, pattern = "Efl") ~ "Fern",
        str_detect(SAMPLE, pattern = "Cco") ~ "Fern",
        #gymnos
        str_detect(SAMPLE, pattern = "Gbi") ~ "Gymno",
        str_detect(SAMPLE, pattern = "Sse") ~ "Gymno",
        str_detect(SAMPLE, pattern = "Tta") ~ "Gymno",
        #str_detect(SAMPLE, pattern = "Tta-2") ~ "Gymno",
        #str_detect(SAMPLE, pattern = "Tta-3") ~ "Gymno",
        str_detect(SAMPLE, pattern = "LPE") ~ "Gymno",
        str_detect(SAMPLE, pattern = "Efr") ~ "Gymno",
        #str_detect(SAMPLE, pattern = "Efr-2") ~ "Gymno",
        str_detect(SAMPLE, pattern = "Wno") ~ "Gymno",
        str_detect(SAMPLE, pattern = "Tcr") ~ "Gymno",
        #str_detect(SAMPLE, pattern = "Tcr-2") ~ "Gymno",
        str_detect(SAMPLE, pattern = "Aau") ~ "Gymno",
        str_detect(SAMPLE, pattern = "Ptr") ~ "Gymnos",
        str_detect(SAMPLE, pattern = "Ppa") ~ "Gymno",       
        #angio
        str_detect(SAMPLE, pattern = "Isi") ~ "Angio",
        str_detect(SAMPLE, pattern = "Ihe") ~ "Angio",
        str_detect(SAMPLE, pattern = "Dwi") ~ "Angio",
        str_detect(SAMPLE, pattern = "Dta") ~ "Angio",
        str_detect(SAMPLE, pattern = "Pco") ~ "Angio",
        str_detect(SAMPLE, pattern = "Pda") ~ "Angio",
        str_detect(SAMPLE, pattern = "Shi") ~ "Angio",
        str_detect(SAMPLE, pattern = "Mgr") ~ "Angio",
        str_detect(SAMPLE, pattern = "Mdo") ~ "Angio",
        str_detect(SAMPLE, pattern = "Aun") ~ "Angio",
        str_detect(SAMPLE, pattern = "Jni") ~ "Angio",
        #str_detect(SAMPLE, pattern = "Jni-2") ~ "Angio",    
        )
    )
  
  #xrf<-xrf%>%
  #mutate(replicate_split = str_split_fixed(replicate,'.',n=2)) %>%
    #create column "samplenumber" from the first element of "replicate"
    #              "rep"      second
  #  mutate(samplenumber = replicate_split[,1],
     #      rep = replicate_split[,2])%>%
   # select(-replicate_split)
colnames(xrf)
  
#   colnames(xrf)
#  [1] "Reading No" "Time"       "Type"       "Duration"   "Units"      "Sequence"   "Flags"     
#  [8] "SAMPLE"     "LOCATION"   "INSPECTOR"  "MISC"       "NOTE"       "User Login" "Ba"        
#  [15] "Ba Error"   "Sb"         "Sb Error"   "Sn"         "Sn Error"   "Cd"         "Cd Error"  
#  [22] "Pd"         "Pd Error"   "Ag"         "Ag Error"   "Bal"        "Bal Error"  "Mo"        
#  [29] "Mo Error"   "Nb"         "Nb Error"   "Zr"         "Zr Error"   "Sr"         "Sr Error"  
#  [36] "Rb"         "Rb Error"   "Bi"         "Bi Error"   "As"         "As Error"   "Se"        
#  [43] "Se Error"   "Au"         "Au Error"   "Pb"         "Pb Error"   "W"          "W Error"   
#  [50] "Zn"         "Zn Error"   "Cu"         "Cu Error"   "Ni"         "Ni Error"   "Co"        
#  [57] "Co Error"   "Fe"         "Fe Error"   "Mn"         "Mn Error"   "Cr"         "Cr Error"  
#  [64] "V"          "V Error"    "Ti"         "Ti Error"   "Ca"         "Ca Error"   "K"         
#  [71] "K Error"    "Al"         "Al Error"   "P"          "P Error"    "Si"         "Si Error"  
#  [78] "Cl"         "Cl Error"   "S"          "S Error"    "Mg"         "Mg Error"   "genus"     
#  [85] "replicate"  "plate"      "method"


  
  #separate out the data columns and the error columns
  #check the column names using the colnames
  #the first element is the 14th column, the last is the 83rd
  element.colnums<-seq(14,83,2)
  xrf[element.colnums] <- sapply(xrf[element.colnums],as.numeric)
  elements<-xrf[element.colnums]
  element_names <- colnames(elements)
}

#-----boxplots------------------------------------------------------------------
##Boxplots
if(input_switch=="Sid"){
  
  
  #adjust the 'shape' of the data from a matrix of elements to a table where 
  #each row is a single element concentration with a other categorical variables
  #like taxon and locality
  xrf_long <- xrf%>%
    pivot_longer(cols = all_of(element_names),names_to = c("element"))%>%
    rename(ppm=value)
  
  
  #if you want to filter rows, keep certain ones or drop others
  # xrf_long <- filter(locality=="Astartekløft", 
  #                    taxon!=c("GINKGOITES MINUTA", "GINKGOITES")) 
  
  a_big_boxplot <- xrf_long %>%  
    ggplot()+#just to get things started
    geom_boxplot(aes(x=species, y=ppm, color=time, fill = time))+#make boxplot shapes, separate in space using taxon, and color using locality
    facet_grid(element~site, scales="free")+#separate element plots out into separate panels
    theme(axis.text.x=element_text(angle=45, hjust = 1))
  
  #some plot saving code
  #useful for making nice plots with high res at a particular size  
  #i just played with the size til it looked good
  #it saves into whatever directory you set at the the top
  plotName<-"Sid_element_boxplots.png"
  ggsave(plotName, plot = a_big_boxplot, device = "png", path = ".",
         scale = 1, height = 16, width = 14, units = c("in"),
         dpi = 300, limitsize = TRUE)
  
  #simplified a bit, focused on locality
  locality_boxplot<-xrf_long %>%  
    ggplot()+#just to get things started
    geom_boxplot(aes(x=locality, y=ppm, color=locality))+#make boxplot shapes, separate in space using taxon, and color using locality
    facet_wrap(element~., scales="free")+#separate element plots out into separate panels
    theme(axis.text.x=element_text(angle=45, hjust = 1))
  
  #some plot saving code
  #useful for making nice plots with high res at a particular size  
  #i just played with the size til it looked good
  #it saves into whatever directory you set at the the top
  plotName<-"Anto_elements_locality_boxplots.png"
  ggsave(plotName, plot = locality_boxplot, device = "png", path = ".",
         scale = 1, height = 16, width = 14, units = c("in"),
         dpi = 300, limitsize = TRUE)
  
  
  #simplified a bit, focused on taxon
  taxon_boxplot<-xrf_long %>%  
    ggplot()+#just to get things started
    geom_boxplot(aes(x=taxon, y=ppm, color=taxon))+#make boxplot shapes, separate in space using taxon, and color using locality
    facet_wrap(element~., scales="free")+#separate element plots out into separate panels
    theme(axis.text.x=element_text(angle=45, hjust = 1))
  
  #some plot saving code
  #useful for making nice plots with high res at a particular size  
  #i just played with the size til it looked good
  #it saves into whatever directory you set at the the top
  plotName<-"Anto_elements_taxon_boxplots.png"
  ggsave(plotName, plot = taxon_boxplot, device = "png", path = ".",
         scale = 1, height = 16, width = 14, units = c("in"),
         dpi = 300, limitsize = TRUE)
  
  
  
  
}

if(input_switch=="Ant"){
  #adjust the 'shape' of the data from a matrix of elements to a table where 
  #each row is a single element concentration with a other categorical variables
  #like taxon and locality
  xrf_long <- xrf%>%
    pivot_longer(cols = all_of(element_names),names_to = c("element"))%>%
    rename(ppm=value)
  
  #if you want to filter rows, keep certain ones or drop others
  # xrf_long <- filter(locality=="Astartekløft", 
  #                    taxon!=c("GINKGOITES MINUTA", "GINKGOITES")) 
  
  a_big_boxplot <- xrf_long %>%  
    ggplot()+#just to get things started
    geom_boxplot(aes(x=taxon, y=ppm, color=element, fill = locality))+#make boxplot shapes, separate in space using taxon, and color using locality
    facet_wrap(element~., scales="free")+#separate element plots out into separate panels
    theme(axis.text.x=element_text(angle=45, hjust = 1))
  
  #some plot saving code
  #useful for making nice plots with high res at a particular size  
  #i just played with the size til it looked good
  #it saves into whatever directory you set at the the top
  plotName<-"Anto_element_boxplots.png"
  ggsave(plotName, plot = a_big_boxplot, device = "png", path = ".",
         scale = 1, height = 16, width = 14, units = c("in"),
         dpi = 300, limitsize = TRUE)
  
  #simplified a bit, focused on locality
  locality_boxplot<-xrf_long %>%  
    ggplot()+#just to get things started
    geom_boxplot(aes(x=locality, y=ppm, color=locality))+#make boxplot shapes, separate in space using taxon, and color using locality
    facet_wrap(element~., scales="free")+#separate element plots out into separate panels
    theme(axis.text.x=element_text(angle=45, hjust = 1))
  
  #some plot saving code
  #useful for making nice plots with high res at a particular size  
  #i just played with the size til it looked good
  #it saves into whatever directory you set at the the top
  plotName<-"Anto_elements_locality_boxplots.png"
  ggsave(plotName, plot = locality_boxplot, device = "png", path = ".",
         scale = 1, height = 16, width = 14, units = c("in"),
         dpi = 300, limitsize = TRUE)
  
  
  #simplified a bit, focused on taxon
  taxon_boxplot<-xrf_long %>%  
    ggplot()+#just to get things started
    geom_boxplot(aes(x=taxon, y=ppm, color=taxon))+#make boxplot shapes, separate in space using taxon, and color using locality
    facet_wrap(element~., scales="free")+#separate element plots out into separate panels
    theme(axis.text.x=element_text(angle=45, hjust = 1))
  
  #some plot saving code
  #useful for making nice plots with high res at a particular size  
  #i just played with the size til it looked good
  #it saves into whatever directory you set at the the top
  plotName<-"Anto_elements_taxon_boxplots.png"
  ggsave(plotName, plot = taxon_boxplot, device = "png", path = ".",
         scale = 1, height = 16, width = 14, units = c("in"),
         dpi = 300, limitsize = TRUE)
  
  
  
  
}

##----Cat Penn----
if(input_switch=="Cat_Penn"){
  #adjust the 'shape' of the data from a matrix of elements to a table where 
  #each row is a single element concentration with a other categorical variables
  #like taxon and locality
  categories<-c( "spot", "substrate", "core", "LOCATION", "taxon")
  
  xrf_long <- xrf%>%select(all_of(c(categories,element_names)))%>%
    pivot_longer(cols = all_of(element_names),names_to = c("element"))%>%
    rename(ppm=value) %>% filter(!is.na(ppm))
  
  #if you want to filter rows, keep certain ones or drop others
  #xrf_long <- filter(xrf_long, element==c("As","Pb","Cu"))
  #                    taxon!=c("GINKGOITES MINUTA", "GINKGOITES")) 
  
  a_big_boxplot <- xrf_long %>%  
    ggplot()+ #just to get things started
    geom_boxplot(aes(x=LOCATION, y=ppm, color=substrate))+ #make boxplot shapes, separate in space using taxon, and color using locality
    facet_wrap(element~., scales="free")+ #separate element plots out into separate panels
    theme(axis.text.x=element_text(angle=45, hjust = 1))
  
  #some plot saving code
  #useful for making nice plots with high res at a particular size  
  #i just played with the size til it looked good
  #it saves into whatever directory you set at the the top
  plotName<-"Cat_element_boxplots.png"
  ggsave(plotName, plot = a_big_boxplot, device = "png", path = ".",
         scale = 1, height = 16, width = 25, units = c("in"),
         dpi = 300, limitsize = TRUE)
  
  #simplified a bit, focused on locality
  locality_boxplot <- xrf_long %>%
    ggplot()+#just to get things started
    geom_boxplot(aes(x=locality, y=ppm, color=LOCATION))+#make boxplot shapes, separate in space using taxon, and color using locality
    facet_wrap(element~., scales="free")+#separate element plots out into separate panels
    theme(axis.text.x=element_text(angle=45, hjust = 1))
  
  #some plot saving code
  #useful for making nice plots with high res at a particular size  
  #i just played with the size til it looked good
  #it saves into whatever directory you set at the the top
  plotName<-"Cat_fire_elements_locality_boxplots.png"
  ggsave(plotName, plot = locality_boxplot, device = "png", path = ".",
         scale = 1, height = 16, width = 14, units = c("in"),
         dpi = 300, limitsize = TRUE)
  
  
  #simplified a bit, focused on taxon
  taxon_boxplot<-xrf_long %>%  
    ggplot()+#just to get things started
    geom_boxplot(aes(x=taxon, y=ppm, color=taxon))+#make boxplot shapes, separate in space using taxon, and color using locality
    facet_wrap(element~., scales="free")+#separate element plots out into separate panels
    theme(axis.text.x=element_text(angle=45, hjust = 1))
  
  #some plot saving code
  #useful for making nice plots with high res at a particular size  
  #i just played with the size til it looked good
  #it saves into whatever directory you set at the the top
  plotName<-"Cat_elements_taxon_boxplots.png"
  ggsave(plotName, plot = taxon_boxplot, device = "png", path = ".",
         scale = 1, height = 16, width = 14, units = c("in"),
         dpi = 300, limitsize = TRUE)
}

##----Cat Modern----
  if(input_switch=="Cat_Modern"){
    #adjust the 'shape' of the data from a matrix of elements to a table where 
    #each row is a single element concentration with a other categorical variables
    #like taxon and locality
    categories<-c( "genus", "plate", "method", "group")
    
    xrf_long <- xrf%>%
      select(all_of(c(categories,element_names)))%>%
      pivot_longer(cols = all_of(element_names),names_to = c("element"))%>%
      rename(ppm=value) %>%
      filter(!is.na(ppm))
    
    #if you want to filter rows, keep certain ones or drop others
    #xrf_long <- filter(xrf_long, is.na(plate))
    xrf_long <- filter (xrf_long, !genus %in% c("Arbutus unido",
                                            "Dicksonia? fibrosa?",
                                            "Dryopteris cycadena",
                                            "Equisetum fluviata",
                                            "Magnolia grandifolia",
                                            "Woodwardia fimbriata"
    ))
    
    
    xrf_long <- filter(xrf_long, element==c("Ti"))
                 

###----genus--------------------------------------------------------------------
    a_big_boxplot <- xrf_long %>%  
      ggplot()+#just to get things started
      #make boxplot shapes
      geom_boxplot(aes(x=genus,
                       y=ppm,
                       color=genus
                         ))+
      facet_wrap(element~., scales="free")+#separate element plots out into separate panels
      theme(axis.text.x=element_text(angle=45, hjust = 1))
    
    #some plot saving code
    #useful for making nice plots with high res at a particular size  
    #i just played with the size til it looked good
    #it saves into whatever directory you set at the the top
    plotName<-"Cat_element_boxplots.png"
    ggsave(plotName, plot = a_big_boxplot, device = "png", path = ".",
           scale = 1, height = 16, width = 25, units = c("in"),
           dpi = 300, limitsize = TRUE)
    
###----group------------------------------------------------------------------
    group_boxplot <- xrf_long %>%
      ggplot()+#just to get things started
      #make boxplot shapes, separate in space using taxon, and color using locality
      geom_boxplot(aes(x=group, y=ppm, color=group))+ 
      facet_wrap(element~., scales="free")+ #separate element plots out into separate panels
      theme(axis.text.x=element_text(angle=45, hjust = 1))
    
    #some plot saving code
    #useful for making nice plots with high res at a particular size  
    #i just played with the size til it looked good
    #it saves into whatever directory you set at the the top
    plotName<-"Cat_group_elements_boxplots.png"
    ggsave(plotName, plot = group_boxplot, device = "png", path = ".",
           scale = 1, height = 15, width = 14, units = c("in"),
           dpi = 300, limitsize = TRUE)
    
###----method------------------------------------------------------------------
    #simplified a bit, focused on group
    Ti_boxplot<-xrf_long %>%  
      ggplot()+#just to get things started
      #make boxplot shapes, separate in space using method, and color using locality
      geom_boxplot(aes(x= genus, y=ppm, color=plate))+
      facet_wrap(element~., scales="free")+#separate element plots out into separate panels
      theme(axis.text.x=element_text(angle=45, hjust = 1))
    
    #some plot saving code
    #useful for making nice plots with high res at a particular size  
    #i just played with the size til it looked good
    #it saves into whatever directory you set at the the top
    plotName<-"Ti_depthtest_boxplots.png"
    ggsave(plotName, plot = Ti_boxplot, device = "png", path = ".",
           scale = 1, height = 10, width = 20, units = c("in"),
           dpi = 300, limitsize = TRUE)
}

#----PCA and correlation matrix----  


#if doing pre v post mlpb
#run untitled 1.2 here

##remove columns with zero variance, which causes errors for subsequent analyses
#find zero variance columns
var_nz <- function(x) !is.na(var(x[x != 0]))
varying_elements_map<-apply(elements,2,var_nz)
#update data coulmns and names to only those with varianc
elements <- elements[,varying_elements_map]
element_names <- colnames(elements)


#individual tests of normality
elements %>% ungroup %>%  
  shapiro_test(element_names) %>% 
  arrange(variable)
#guide for interprting this and the next one:
#https://www.datanovia.com/en/lessons/normality-test-in-r/#shapiro-wilks-normality-test
#if p > 0.05 then data is normal

#test of multivariate normality
elements %>% mshapiro_test
#the data are non-normal individually and as a multivariate group
#this is probably ok on it's own for PCA

#computes the correlations between the columns of the elements df
ele_cor <- cor(elements) %>%
  #computes their absolute value
  abs


#correlation matrix with significance and correlation coefficients in top 
#right triangle, histograms on diagonal, and biplots with best fit lines in red
#this is ugly, but save it in a large format, zoom in and look at the 
#relationships between variables, are they linear?
chart.Correlation(elements)

#couple of different options for removing some columns
#things that don't have correlations about 0.3 with anything
#this is the more important thing for PCA

if(input_switch=="Sid"){
low_cor <- c("Ba","Sb","Rb")
#elements with normal-ish distributions based on the histograms
normal <- c("Bal","Mo","Nb","Cr","V","Ti","Ca","K")
}

if(input_switch=="Ant"){
  low_cor <- c("Ag","Mo")
  #elements with normal-ish distributions based on the histograms
  normal <- c("S","Cl","Si","P","K","Ca","Cr")
}

if(input_switch=="Cat_Penn"){
  low_cor <- c("Zr")
  #elements with normal-ish distributions based on the histograms
  normal <- c("Bal","Nb","Rb","Cr","V","K","Si")
}

#normal elements
ele_norm <- elements%>%select(normal)#select(!low_cor) 
ele_norm_names <- colnames(ele_norm)
#individual tests of normality ... again
ele_norm %>% ungroup %>%  
  shapiro_test(ele_norm_names) %>% 
  arrange(variable)


#test of multivariate normality ... again
ele_norm %>% mshapiro_test
#still not normal

#to understand the linear relationship assumption let's break the correlation 
#matrix down into the normal and non-normal elements
chart.Correlation(ele_norm)
ele_notnorm <- elements%>%select(!normal)#select(!low_cor) 
chart.Correlation(ele_notnorm)

#filtering out the few elements with very low correlations
#this is just going to clarify things in the PCA
ele_cor <- elements%>%select(!low_cor) 
ele_cor_names <- colnames(ele_cor)
corr_mat_cor <- ele_cor %>% cor

#Conclusion: the outputs of this PCA are suspect because the underlying 
#relationships may not be linear, particularly for the non-normal elements
#you can still use PCA as an exploratory technique, but refer back to the
#correlation matrix when making element-level interpretations


ele_cor<-data.frame(ele_cor)

#if doing pre vs post mlpb 
#run Untitled 2 here

#rownames(ele_cor) <- xrf$SAMPLE
# https://www.sthda.com/english/wiki/wiki.php?id_contents=7851
res.pca <- PCA(ele_cor, graph = FALSE)
eigenvalues <- res.pca$eig
head(eigenvalues[, 1:2])
#the variance explained in pc1 is 39%
#the variance explained in pc2 is 20%

#contributions of each element to each dimension
res.pca$var$contrib
#pc1 is mostly K, si, Bal, Al
#pc2 is mostly Fe, Ca, Nb, Cr

#fviz_pca_var(res.pca)
#the first two dimensions of the PCA don't single out any elements
#as driving variation
#fviz_pca_var(res.pca, axes = c(1,3))
#the third two dimension (vertical there) 
#collapses a bit

#fviz_pca_ind(res.pca, label="none", habillage = xrf$locality)
#no separation of old and new

#fviz_pca_ind(res.pca, axes = c(1,3), label="none", habillage = xrf$locality)
#no separation of old and new
#but remember, the PCA (1) didn't work very well and (2) is not designed
#to test for differences between groups

###try nsprcomp for constrained PCA
#nn_pca <- nsprcomp(elements, ncomp = 4, scale.=T, center = F,  nneg = T)
#fviz_pca_var(nn_pca)
#fviz_pca_ind(nn_pca)
#summary(nn_pca)

#nn_comp_pca <- nscumcomp(elements, ncomp=4, k=150, scale.= T, nneg=TRUE, gamma=1)
#fviz_pca_var(nn_comp_pca)
#fviz_pca_ind(nn_comp_pca)

###try compositions package
#x <- rcomp(elements)
#compo_pca<-PCA(x) #removes missing values


# compo_pca <- princomp(x)
#plot(compo_pca,habillage = xrf$locality)

#Investigate(compo_pca)


#--------Cat_Penn x=time y=pc1 and 2-----------------------------------------------------

ggplot(pca_df)+
#  geom_point(aes(
 #   x = xrf$Age,
  #  y = Dim.1,
#    color = xrf$LOCATION,
 #   shape = xrf$taxon)
  #  ) +
  geom_point(aes(
    x = xrf$Age,
    y = Dim.2,
    color = xrf$LOCATION,
    shape = xrf$taxon)
  )
