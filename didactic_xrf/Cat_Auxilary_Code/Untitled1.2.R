xrf<-xrf%>%
  mutate(
    Age = case_when(
      str_detect(LOCATION, pattern = "Murphysboro") ~ "311.5",
      str_detect(LOCATION, pattern = "Paum Mine") ~ "311.5",
      str_detect(LOCATION, pattern = "Springfield") ~ "307.5",
      str_detect(LOCATION, pattern = "Chapel") ~ "305.79",
      str_detect(LOCATION, pattern = "Womac") ~ "305.49",
      str_detect(LOCATION, pattern = "New Haven") ~ "305.07",
      .default = as.character(LOCATION))) %>%
      
  mutate(
    LOCATION = case_when(
      #murphysboro and springfield to pre
      str_detect(LOCATION, pattern = "Murphysboro") ~ "Pre-MLPB",
      str_detect(LOCATION, pattern = "Paum Mine") ~ "Pre-MLPB",
      str_detect(LOCATION, pattern = "Springfield") ~ "Pre-MLPB",
      #chapel, womac, and new haven to post
      str_detect(LOCATION, pattern = "Chapel") ~ "Post-MLPB",
      str_detect(LOCATION, pattern = "Womac") ~ "Post-MLPB",
      str_detect(LOCATION, pattern = "New Haven") ~ "Post-MLPB",
      .default = as.character(LOCATION)
    ))
