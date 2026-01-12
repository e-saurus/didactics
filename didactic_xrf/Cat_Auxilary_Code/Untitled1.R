#code for splitting Pennsylvanian data into into Pre Post MLPB


  #adjust the 'shape' of the data from a matrix of elements to a table where 
  #each row is a single element concentration with a other categorical variables
  #like taxon and locality
  categories<-c( "spot", "substrate", "core", "LOCATION","taxon")
  
  xrf_long <- xrf%>%select(all_of(c(categories,element_names)))%>%
    pivot_longer(cols = all_of(element_names),names_to = c("element"))%>%
    rename(ppm=value) %>% filter(!is.na(ppm))

  xrf_long<-xrf_long%>%
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

  a_big_boxplot <- xrf_long %>%  
    ggplot()+#just to get things started
    geom_boxplot(aes(x=LOCATION, y=ppm, color=substrate))+#make box plot shapes, separate in space using taxon, and color using locality
    facet_wrap(element~., scales="free")+#separate element plots out into separate panels
    theme(axis.text.x=element_text(angle=45, hjust = 1))  
  
  plotName<-"Cat_prevpost_element_boxplots.png"
  ggsave(plotName, plot = a_big_boxplot, device = "png", path = ".",
         scale = 1, height = 16, width = 25, units = c("in"),
         dpi = 300, limitsize = TRUE)  
  