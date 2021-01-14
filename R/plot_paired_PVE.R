paired_PVE_plot <- function(PVE, var_name = "Variance", minherit=0){
  meanPVE <- as.data.frame(PVE[!grepl(paste0(var_name), PVE$PaperName), ])
  varPVE <- as.data.frame(PVE[grepl(paste0(var_name), PVE$PaperName), ])
  # Add non-mean/variance to mean table
  varPVE <- rbind(varPVE, PVE[!(PVE$PaperName %in% c(meanPVE$PaperName, varPVE$PaperName)),])
  varPVE$basename <- gsub(paste0(var_name), "", varPVE$PaperName)
  # Adding false variance rows for non-mean/variance phenotypes
  novar <- setdiff(meanPVE$PaperName, varPVE$basename)
  addv <- meanPVE[meanPVE$PaperName %in% novar,]
  addv$basename <- addv$PaperName
  addv[,"PVE"] <- NA
  addv[, "PVESE"] <- NA
  varPVE <- rbind(varPVE, addv)
  rownames(varPVE) <- varPVE$basename
  varPVE <- varPVE[meanPVE$PaperName, ]
  # Plot the plots
  mean_pvep <-
    ggplot(meanPVE, aes(reorder(PaperName, -meanPVE$PVE), PVE, fill = Group, alpha = PVE>minherit)) + geom_bar(color =
                                                                               "black", stat = "identity") +
    geom_errorbar(aes(ymin = PVE - PVESE, ymax = PVE + PVESE), width = .2) + ylim(0,1) +
    scale_alpha_manual(values=c(0.6, 1)) +
    xlab("Phenotype") + coord_flip() + ggtitle("Mean")
  var_pvep <-
    ggplot(varPVE, aes(reorder(basename, -meanPVE$PVE), PVE, fill = Group, alpha = PVE>minherit)) + geom_bar(color =
                                                                                      "black", stat = "identity") +
    geom_errorbar(aes(ymin = PVE - PVESE, ymax = PVE + PVESE), width = .2) + ylim(0,1) +
    scale_alpha_manual(values=c(0.6, 1))
    xlab("") + coord_flip() + ggtitle("Variance")

  return(list(mean_plot = mean_pvep, var_plot = var_pvep))

  }
