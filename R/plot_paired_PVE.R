paired_PVE_plot <- function(PVE, mean_name="Mean", var_name = "Variance"){
  meanPVE <- PVE[grepl(paste0(mean_name, "$"), PVE$PaperName), ]
  varPVE <- PVE[grepl(paste0(var_name, "$"), PVE$PaperName), ]
  # Add non-mean/variance to mean table
  meanPVE <- rbind(meanPVE, PVE[!(PVE$PaperName %in% c(meanPVE$PaperName, varPVE$PaperName)),])
  varPVE <- rbind(varPVE, PVE[!(PVE$PaperName %in% c(meanPVE$PaperName, varPVE$PaperName)),])
  meanPVE$basename <- gsub(paste0(mean_name, "$"), "", meanPVE$PaperName)
  varPVE$basename <- gsub(paste0(var_name, "$"), "", varPVE$PaperName)
  # Adding false variance rows for non-mean/variance phenotypes
  novar <- grepl(paste0(var_name, "$"), varPVE$PaperName)
  varPVE[novar, "PVE"] <- 0
  varPVE[novar, "PaperName"] <- ""
  varPVE[novar, "PVESE"] <- 0

  pvorder <- sort(meanPVE$PVE, decreasing = T)
  # Plot the plots
  mean_pvep <-
    ggplot(meanPVE, aes(reorder(PaperName, pvorder), PVE, fill = Group)) + geom_bar(color =
                                                                               "black", stat = "identity") +
    geom_errorbar(aes(ymin = PVE - PVESE, ymax = PVE + PVESE), width = .2) +
    xlab("Phenotype") + coord_flip()
  var_pvep <-
    ggplot(varPVE, aes(reorder(PaperName, pvorder), PVE, fill = Group)) + geom_bar(color =
                                                                                      "black", stat = "identity") +
    geom_errorbar(aes(ymin = PVE - PVESE, ymax = PVE + PVESE), width = .2) +
    xlab("Phenotype") + coord_flip()

  return(list(mean_plot = mean_pvep, var_plot = var_pvep))

  }
