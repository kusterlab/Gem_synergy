# To run this code, please download the files "03_Kinobeads_curves_forplotting.txt" and "03_Kinobeads_proteingroups_forplotting.txt" from Zenodo (see README). 
# These files should be imported as dataframes called "curves" and "proteingroups", respectively.
# Before staring the script, run the following lines:

LL.4.function <- function(B, x){
  (B[2]-B[3])/(1+(x/B[1])^B[4])+B[3] 
}
xlim <- c(0.4, 40000)
ylim <- c(0, 1.5)
ylabels <- seq(0,1.6, by = 0.2)
many.x.values <- c(seq(0.5, 1, by = 0.05), seq(1, 10, by = 0.1),seq(10, 100, by = 1), seq(100, 1000, by = 2),seq(1000, 35000, by = 100))


doses.withoutDMSO = c( 1, 3, 10, 30, 100, 300, 1000, 3000, 30000)
doses.toplot = c(0, 1, 3, 10, 30, 100, 300, 1000, 3000, 10000, 30000)




############################################################
# Fig. 3b (heatmap)
############################################################

# Select relevant columns, filter for Targets, arrange genes by alphabet
curves_targets <- curves %>% dplyr::select(Target, Genes, pKdapp, Experiment) %>% dplyr::filter(Target ==T & Kinase == T)
curves_targets$Gene.names <- factor(curves_targets$Genes, levels = sort(unique(curves_targets$Genes)))

# Arrange data
curves_targets$Experiment <- factor(curves_targets$Experiment, levels = rev(c('Gartisertib','Berzosertib','Ceralasertib','Elimusertib')))

# vertical heatmap
ggplot(curves_targets, aes(Experiment, factor(Genes, sort(unique(Genes), decreasing = T)))) +
  geom_tile(aes(fill = pKdapp), colour = "white") +
  theme_classic()+
  scale_fill_viridis_c(option = "magma", na.value = 'white',limits=c(5.5,10.5),direction=-1,oob = scales::squish)+ #oob=squish
  theme(axis.text.x = element_text(angle = 35, vjust = 1, hjust=1, size =10), panel.border = element_rect(colour = "black", fill=NA, size=1),
        axis.line = element_line(colour = "black"), axis.title = element_blank(),
        axis.text = element_text(size =10, colour = "black"), legend.text= element_text(size = 12), legend.title = element_text(size = 12, face = "bold"))




############################################################
# Fig. 3c (curve plots)
############################################################

colors = c("#54309c", "#E08766")

dev.off()
par(tck = -.02, mgp = c(1.5, 0.5, 0), mar = c(6,5,3.5,3.5))
plot(x = 'n', y = 'n',
     log = 'x',
     xlab = "", ylab = "", 
     xlim = xlim, ylim = ylim, 
     axes = F)

param <- curves %>% dplyr::filter(Experiment =="Elimusertib" & Target == T & Kinase == T)
for (i in 1:nrow(param)){
  parameters <- param[i,] %>% dplyr::select(c('Curve.R2','EC50','Curve.Front', 'Curve.Back', 'Curve.Slope', 'Kdapp', 'Curve.Fold.Change'))  # EC50 (in nM), top, bottom, slope
  
  if(!is.na(parameters$Curve.R2)){
    Gene_ID <- as.character(param$Genes[i])
    cols <- colors[i]
    
    y.response <- param %>% dplyr::filter(Genes == Gene_ID) %>% dplyr::select(starts_with("Ratio")) # subset row by row
    
    colnames(y.response)[-1] <- doses.withoutDMSO
    colnames(y.response)[1] <- 0.3
    
    points(x = c(0.3, doses.withoutDMSO), y = c(1, t(y.response)[2:ncol(y.response)]),
           log = 'x',
           legend = F, 
           cex = 1.2,
           cex.main = 1.6,
           bg = cols,
           lty =1,
           pch = 16, 
           xlab = "", ylab = "", 
           xlim = xlim, ylim = ylim,
           col = cols,
           axes = F)
    
    four.parameters <-t(parameters)[2:5] 
    y.fitted.curve <- LL.4.function(four.parameters, many.x.values)
    lines(many.x.values, y.fitted.curve, type = "l", log = "x", lwd = 3, col= cols)
  }
  add=T
}
title(ylab = "Residual Binding", cex.lab = 1.2, line = 2)
title(xlab = "Concentration [nM]", cex.lab = 1.2, line = 3.6)

axis(1, at =c(0.3, doses.toplot), labels=c(0, doses.toplot), lwd = 0, lwd.ticks = 1, las =2)
axis(2, at =seq(0, 1.4, 0.2), labels =format(seq(0, 1.4, 0.2),digits=2), lwd = 0, lwd.ticks = 1, las = 2)
abline(h = 1, lty = 'dotted', lwd=1)

legendtext <- param$Genes
lgd <- legend('topright', legendtext, plot = F)
legend('topright',box.lty=0,  bty="n",
       legend = legendtext, cex =0.8, pch = 19,  lwd=3, col = colors[1:i],
       plot = T)
box(bty="L")




############################################################
# Supplementary Fig. 2a (curve plots)
############################################################

colors = list(c("#005293", "#999999"), c("grey30", "#64a0c8"))
pchs = c(16,16)
Experiment_sets <- list(c("Elimusertib", "Ceralasertib"), c("Gartisertib", "Berzosertib"))

for(n in 1:length(Experiment_sets)){
  par(tck = -.02, mgp = c(1.5, 0.5, 0), mar = c(6,5,3.5,3.5))
  plot(x = 'n', y = 'n',
       log = 'x',
       xlab = "", ylab = "", # Labels
       xlim = xlim, ylim = ylim, # Limits
       axes = F)
  
  param <- curves %>% dplyr::filter(Experiment %in% Experiment_sets[[n]] & Genes == "ATR") #%>% arrange(EC50, desc(r2))
  for (i in 1:2){
    parameters <- param[i,] %>% dplyr::select(c('Curve.R2','EC50','Curve.Front', 'Curve.Back', 'Curve.Slope', 'Kdapp', 'Curve.Fold.Change'))  # EC50 (in nM), top, bottom, slope
    cols <- colors[[n]][i]

    y.response <- param %>% dplyr::filter(Experiment == param[i,]$Experiment & Genes == "ATR") %>% dplyr::select(starts_with("Ratio")) # subset row by row
    
    colnames(y.response)[-1] <- doses.withoutDMSO
    colnames(y.response)[1] <- 0.3
    
    points(x = c(0.3, doses.withoutDMSO), y = c(1, t(y.response)[2:ncol(y.response)]),
           log = 'x',
           legend = F, 
           cex = 1.5,
           cex.main = 1.6,
           bg = cols,
           lty =1,
           pch = pchs[i], 
           xlab = "", ylab = "", 
           xlim = xlim, ylim = ylim,
           col = cols,
           axes = F)
    
    box(bty="L")
    
    four.parameters <-t(parameters)[2:5] 
    y.fitted.curve <- LL.4.function(four.parameters, many.x.values)
    lines(many.x.values, y.fitted.curve, type = "l", log = "x", lwd = 2, col= cols)
  }
  title(ylab = "Residual Binding", cex.lab = 1.2, line = 2)
  title(xlab = "Concentration (nM)", cex.lab = 1.2, line = 3.6)
  title("ATR", line =0.2, cex.main = 1.5)
  axis(1, at =c(0.3, doses.toplot), labels=c(0, doses.toplot), lwd = 0, lwd.ticks = 1, las =2)
  axis(2, at =seq(0, 1.4, 0.2), labels =format(seq(0, 1.4, 0.2),digits=2), lwd = 0, lwd.ticks = 1, las = 2)
  abline(h = 1, lty = 'dotted', lwd=1)
  
  legendtext <- as.character(param$Experiment)
  lgd <- legend('topright', legendtext, plot = F)
  legend('topright',box.lty=0,  bty="n",
         legend = legendtext, pt.cex =1.5, cex =1, pch = pchs, lwd=2, col = colors[[n]],
         plot = T)
}
      
    
  

############################################################
# Supplementary Fig. 2b (barplots)
############################################################

dev.off()
y.response <- proteingroups %>% dplyr::filter(Genes == "ATR") %>% dplyr::select(starts_with("MS.MS.count."), Experiment) %>% 
  pivot_longer(starts_with("MS.MS"), names_to = "dose1", values_to = "count" ) %>% mutate(readout = "MSMScounts") 
y.response$dose1 <- gsub("MS.MS.count.", "",  y.response$dose1)
y.response2 <- proteingroups %>% dplyr::filter(Genes == "ATR") %>% dplyr::select(starts_with("Unique.peptides."), Experiment) %>%  
  pivot_longer(starts_with("Unique"), names_to = "dose1", values_to = "count" ) %>% mutate(readout = "uniquepeptides")

y.response2$dose1 <- gsub("Unique.peptides.", "",  y.response2$dose1)
y.response <- rbind(y.response, y.response2) %>% dplyr::filter(dose1 != "PDPD")
y.response$dose1 <- gsub("DMSO", "0",  y.response$dose1)

y.response$dose1 <- factor(y.response$dose1, levels = sort(unique(as.numeric(y.response$dose1))))
y.response$Experiment <- factor(y.response$Experiment, levels = c("Elimusertib", "Ceralasertib", "Berzosertib", "Gartisertib"))

ggplot(y.response, aes(fill=readout, y=count, x=dose1)) + 
  geom_bar(position="dodge", stat="identity")+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), axis.text = element_text(colour = "black"),
        axis.ticks = element_line(color = "black"))+ 
  scale_fill_manual(values = c("grey40","grey70")) + 
  scale_y_continuous(breaks= pretty_breaks(),expand = c(0, NA))+
  facet_wrap(~Experiment, scales = "free")+
  theme(strip.background = element_blank(),
        strip.placement = "outside")




############################################################
# Supplementary Plots (Here: Curve plots of all targets per experiment)
############################################################

experiments <- as.character(unique(curves$Experiment))

dev.off()
for(n in 1:length(experiments)){
  df_plots <- curves %>% dplyr::filter(Experiment == experiments[n] & Target ==T & Kinase == T)
  nf <- layout(matrix(c(1:15), ncol = 3, nrow = 5, byrow=TRUE), 
               widths=c(1.5,0.5,1.5), heights =c(1,1,1,1,1))
  
  for (i in 1:nrow(df_plots)){
    parameters <- df_plots[i,] %>% dplyr::select(c('Curve.R2','EC50','Curve.Front', 'Curve.Back', 'Curve.Slope', 'Kdapp', 'Curve.Fold.Change'))  # R2, EC50 (in nM), top, bottom, slope, Kdapp, log2FC
    
    if(!is.na(parameters$Curve.R2)){
      Gene_ID <- as.character(df_plots$Genes[i])
      cols <- "black"
      par(tck = -.03, mgp = c(1.5, 0.5, 0), mar = c(5,4,3.5,3.5))

      y.response <- curves %>% dplyr::filter(Experiment == experiments[n] & Genes == df_plots$Genes[i]) %>% dplyr::select(starts_with("Ratio")) # subset row by row
      
      colnames(y.response)[-1] <- doses.withoutDMSO
      colnames(y.response)[1] <- 0.3
      
      plot(x = c(doses.withoutDMSO), y = c(t(y.response)[2:ncol(y.response)]),
           log = 'x',
           legend = F, 
           cex = 1,
           cex.main = 1.6,
           bg = cols,
           lty =1,
           pch = 19, 
           xlab = "", ylab = "", 
           xlim = xlim, ylim = ylim, 
           col = cols,
           axes = F)
      box(bty="L")

      title(ylab = "Residual Binding", cex.lab = 1.2, line = 2)
      title(xlab = "Concentration (nM)", cex.lab = 1.2, line = 3.6)
      title(Gene_ID, line =0.2, cex.main = 1.5)
      
      axis(1, at =c(doses.toplot)[-1], labels=c(doses.toplot)[-1], lwd = 0, lwd.ticks = 1, las =2)
      axis(2, at =seq(0, 1.4, 0.2), labels =format(seq(0, 1.4, 0.2),digits=2), lwd = 0, lwd.ticks = 1, las = 2)

      four.parameters <-t(parameters)[2:5] # Make numeric vector out of it. Leave out R2
      y.fitted.curve <- LL.4.function(four.parameters, many.x.values)
      lines(many.x.values, y.fitted.curve, type = "l", log = "x", lwd = 1.5)
      
      # Add legend
      parameters <- round(parameters,2)
      legendtext <- c(paste("R2: ", parameters[1], sep=''),
                      paste("EC50: ", parameters[2], ' nM', sep=''),
                      paste("KDapp: ", parameters[6], ' nM', sep=''),
                      paste("Log2 FC: ", parameters[7], sep=''))

      lgd <- legend('bottomleft', legendtext, plot = F)
      legend('bottomleft',box.lty=0, text.width = strwidth(legendtext)[1]*1.5, bty="n", 
             legend = legendtext, cex =0.7, pch = c(20),  col = "transparent",
             plot = T)
      
      # 2. density plots
      par(tck = -.03, mgp = c(2, 0.5, 0),mar = c(4.5,0.6,3.5,2.5))
      df <- proteingroups %>% dplyr::filter(Experiment == experiments[n])
      colnames(df) <- gsub("Raw 0", "LFQ.intensity.0", colnames(df))
      df <- df %>% filter(LFQ.intensity.0 != 0) 
      d <- log2(df$LFQ.intensity.0) 
      
      density.d <- density(d)
      plot(density.d$y,density.d$x, main = "", yaxt = 'n', xaxt='n', xlim=rev(range(density.d$y)), type='l',
           ylab= '', xlab = "", cex.main=1.2, )
      axis(4)
      mtext("Log2 Intensity (Ctrl)", side=2, line=0.7, cex = 0.8)
      polygon(density.d$y,density.d$x, col="grey", border="grey")
      abline(h = log2(df[which(df$Genes == Gene_ID),'LFQ.intensity.0']), 
             col= 'grey30', lwd=3)
      
      # 3. unique peptides and ms ms scans
      proteingroups_temp <- proteingroups %>% dplyr::filter(Experiment == experiments[n])
      unique.peptides.df <- dplyr::select(proteingroups_temp, Genes | starts_with("Unique.peptides.") & -contains("PD"))
      colnames(unique.peptides.df) <- gsub(".0", ".0.3", colnames(unique.peptides.df), fixed = T)
      colnames(unique.peptides.df) <- gsub("Unique.peptides.", "", colnames(unique.peptides.df))
      colnames(unique.peptides.df) <- gsub("nM", "", colnames(unique.peptides.df))
      
      
      unique.peptides.df.melt <- melt(unique.peptides.df, id.vars = "Genes")
      b <- unique.peptides.df.melt[which(unique.peptides.df.melt$Genes==Gene_ID),]
      
      par(tck = -.03, mgp = c(2, 0.5, 0),mar = c(4.5,4,3.5,4.5))
      plot(as.numeric(levels(b$variable))[b$variable], b$value, 
           log="x",
           cex = 1,
           cex.main = 1.6,
           bg = "black",
           xlab = "", ylab = "",
           lty =1,
           pch = 19, 
           xlim = xlim, 
           ylim = c(0, max(b$value)),
           col = "black",
           axes= F)
      
      title(ylab = "Unique peptides", cex.lab = 1.2, line = 1.5)
      axis(2, at = 0:max(b$value), labels= 0:max(b$value), las =1, col.ticks = "black")
      
      # Get MS MS counts for protein from filtered proteingroups.txt
      MSMScounts.df <- dplyr::select(proteingroups_temp, Genes | starts_with("MS.MS.Count.") & -contains("PD"))
      colnames(unique.peptides.df) <- gsub(".0", ".0.3", colnames(unique.peptides.df), fixed = T)
      colnames(MSMScounts.df) <- gsub("Unique.peptides.", "", colnames(unique.peptides.df))
      colnames(MSMScounts.df) <- gsub("nM", "", colnames(unique.peptides.df))
      
      MSMScounts.df.melt <- melt(MSMScounts.df, id.vars = "Genes")
      b <- MSMScounts.df.melt[which(MSMScounts.df.melt$Genes==Gene_ID),]
      
      par(new=TRUE)
      plot(as.numeric(levels(b$variable))[b$variable], b$value, main ="", log="x",
           cex = 1,
           cex.main = 1.4,
           bg = 'red',
           xlab = "", ylab = "",
           lty =1,
           pch = 19, 
           xlim = xlim, 
           ylim = c(0, max(b$value)),
           col = "red",
           axes= F)
      box(bty = 'L')
      
      title(ylab = "MS/MS counts", cex.lab = 1.2, line = 2.5, col.lab = 'red')
      axis(4, at = 0:max(b$value), labels= 0:max(b$value),las =1, col.ticks = 'red' , col.axis = 'red', col = 'red')
      title(xlab = "Concentration [nM]", cex.lab = 1.2, line = 3.6)
      title(Gene_ID, line =0.2, cex.main = 1.5)
      axis(1, at =c(0.3, doses.toplot[-1]), labels=c(doses.toplot), lwd = 0, lwd.ticks = 1, las =2)
      
    }
    if(i %in% seq(1,200,5)){
      mtext(paste0("Target of ", experiments[n]), outer=TRUE, line=-2)
    }
  }
}