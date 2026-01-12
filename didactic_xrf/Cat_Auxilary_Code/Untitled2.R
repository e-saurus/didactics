#Untitled 2
xrf <- xrf %>%
  mutate(SAMPLE_ID = str_c(SAMPLE, `Reading No`, sep = "_")) %>%
  select(-SAMPLE) %>%
  select(-`Reading No`)

rownames(ele_cor) <- xrf$SAMPLE_ID

res.pca <- PCA(ele_cor, graph = FALSE)
eigenvalues <- res.pca$eig
head(eigenvalues[, 1:2])
#the variance explained in pc1 is 33%
#the variance explained in pc2 is 24%

#contributions of each element to each dimension
res.pca$var$contrib
#pc1 is mostly K, si, Bal, Al
#pc2 is mostly Rb, Fe, Ti


#pca plot
pca_df <- as.data.frame(res.pca$ind$coord) %>%
  mutate(
    LOCATION = xrf$LOCATION,
    taxon = xrf$taxon
  )
pca_plt <-  
 ggplot(pca_df, aes(Dim.1, Dim.2)) +
  geom_point(aes(color = LOCATION, fill = LOCATION, shape = taxon), size = 3) +
  labs(
    x = paste0("PC1 (", round(res.pca$eig[1,2], 1), "%)"),
    y = paste0("PC2 (", round(res.pca$eig[2,2], 1), "%)")
  ) +
  scale_shape_manual(values=c(1,2,3,4,22,6,7,8)) +
  theme_classic()

plotName<-"Cat_prevpost_pca_plot.svg"
ggsave(plotName, plot = pca_plt, device = "svg", path = ".",
       scale = 1, height = 16, width = 25, units = c("in"),
       dpi = 300, limitsize = TRUE)  
