# To run this code, please download the files "06_curves_forplotting.txt" from Zenodo (see README). 
# This file should be imported as a dataframe called "curves".


############################################################
# Supplementary Fig. 3b + Supplementary Fig. 4a (barplots)
############################################################

# Prepare the data: summarise counts for ATRi
curves$reg_exp <- paste0(curves$Experiment,"_", curves$`Curve Regulation`)
barplot_df <- curves %>% group_by(Experiment, `Curve Regulation`, reg_exp) %>% 
  dplyr::summarise(n=n()) %>% arrange(desc(n))
barplot_df$`Curve Regulation` <- as.character(barplot_df$`Curve Regulation`)
barplot_df$`Curve Regulation`[1:4] <- "all"
barplot_df$reg_exp <- factor(barplot_df$reg_exp, levels = c('Elimusertib_','Elimusertib_up','Elimusertib_down',
                                                            'Ceralasertib_','Ceralasertib_up','Ceralasertib_down',
                                                            'Berzosertib_','Berzosertib_up','Berzosertib_down',
                                                            'Gartisertib_','Gartisertib_up','Gartisertib_down'))

df_temp <-  curves %>% select(Gene_SitePos, GEM_updown, GEM_sign) %>% na.omit %>% unique

# Add GEM to datatable, needed for Supplementary Fig. 4a later
vector1 <- c("GEM", "all", "GEM_", length((df_temp$Gene_SitePos))) %>% as.data.frame %>% mutate(across(everything(as.character))) %>% t
vector2 <- c("GEM", "up", "GEM_up", sum(df_temp$GEM_updown == 'up')) %>% as.data.frame %>% mutate(across(everything(as.character))) %>% t
vector3 <- c("GEM", "down", "GEM_down", sum(df_temp$GEM_updown == 'down')) %>% as.data.frame %>% mutate(across(everything(as.character))) %>% t

colnames(vector1) <- colnames(barplot_df)
colnames(vector2) <- colnames(barplot_df)
colnames(vector3) <- colnames(barplot_df)

barplot_df <- barplot_df %>% as.data.frame %>% mutate(across(everything(as.character)))
barplot_df <- rbind(barplot_df, (vector1), (vector2), (vector3))

barplot_df$reg_exp <- factor(barplot_df$reg_exp, levels = c('GEM_', 'GEM_up', 'GEM_down',
                                                            'Elimusertib_','Elimusertib_up','Elimusertib_down',
                                                            'Ceralasertib_','Ceralasertib_up', 'Ceralasertib_down',
                                                            'Berzosertib_','Berzosertib_up','Berzosertib_down',
                                                            'Gartisertib_','Gartisertib_up','Gartisertib_down'))
barplot_df$n <- as.numeric(barplot_df$n)

# prepare table with SQ/TQ counts, first for ATRi and then for GEM
barplot_df_SQTQ <- curves %>% dplyr::filter(SQTQ==T) %>%  group_by(Experiment, `Curve Regulation`, reg_exp) %>% 
  dplyr::summarise(n=n())%>% arrange(desc(n))
barplot_df_SQTQ$`Curve Regulation` <- as.character(barplot_df_SQTQ$`Curve Regulation`)
barplot_df_SQTQ$`Curve Regulation`[1:4] <- "all"

barplot_df_SQTQ$reg_exp <- factor(barplot_df_SQTQ$reg_exp, levels = c('Elimusertib_','Elimusertib_up','Elimusertib_down',
                                                                      'Ceralasertib_','Ceralasertib_up','Ceralasertib_down',
                                                                      'Berzosertib_','Berzosertib_up','Berzosertib_down',
                                                                      'Gartisertib_','Gartisertib_up','Gartisertib_down'))

df_temp <-  curves %>% dplyr::filter(SQTQ==T) %>% select(Gene_SitePos, GEM_updown, GEM_sign) %>% na.omit %>% unique
vector1 <- c("GEM", "all", "GEM_", length((df_temp$Gene_SitePos))) %>% t
vector2 <- c("GEM", "up", "GEM_up", sum(df_temp$GEM_updown == 'up')) %>% as.data.frame %>% mutate(across(everything(as.character))) %>% t
vector3 <- c("GEM", "down", "GEM_down", sum(df_temp$GEM_updown == 'down')) %>% as.data.frame %>% mutate(across(everything(as.character))) %>% t

colnames(vector1) <- colnames(barplot_df_SQTQ)
colnames(vector2) <- colnames(barplot_df_SQTQ)
colnames(vector3) <- colnames(barplot_df_SQTQ)

barplot_df_SQTQ <- barplot_df_SQTQ %>% as.data.frame %>% mutate(across(everything(as.character)))
barplot_df_SQTQ <- rbind(barplot_df_SQTQ, (vector1), (vector2), (vector3))

barplot_df_SQTQ$reg_exp <- factor(barplot_df_SQTQ$reg_exp, levels = levels(barplot_df$reg_exp))
barplot_df_SQTQ$n <- as.numeric(barplot_df_SQTQ$n)

# Plot 4 ATRi barplots
ggplot((barplot_df %>% dplyr::filter(Experiment != "GEM")), aes(fill=`Curve Regulation`, y=n, x=reg_exp)) + 
  geom_bar(position="stack", stat="identity", colour="grey30")+ scale_y_break(c(800,14000), scales =0.4 ) + ylim(0,18000)+
  theme(legend.title=element_blank())+
  scale_fill_manual(values = c("grey85",'#2171b5','#6baed6' ))+#c("grey85", "#0065bd", "grey40"))+
  geom_bar(data=barplot_df_SQTQ%>% dplyr::filter(Experiment != "GEM"), aes(fill=`Curve Regulation`, y=n, x=reg_exp), position="stack", stat="identity", fill = 'red', colour="grey30")

# Plot GEM only barplots
ggplot((barplot_df %>% dplyr::filter(Experiment == "GEM")), aes(fill=`Curve Regulation`, y=n, x=reg_exp)) + 
  geom_bar(position="stack", stat="identity", colour="grey30")+ scale_y_break(c(800,14000), scales =0.4) + ylim(0,25000)+
  theme(legend.title=element_blank())+
  scale_fill_manual(values = c("grey85",'#2171b5','#6baed6' ))+#c("grey85", "#0065bd", "grey40"))+
  geom_bar(data=barplot_df_SQTQ%>% dplyr::filter(Experiment == "GEM"), aes(fill=`Curve Regulation`, y=n, x=reg_exp), position="stack", stat="identity", fill = 'red', colour="grey30")




############################################################
# Supplementary Fig. 3b + Supplementary Fig. 4a (pie charts)
############################################################

pie_temp <- barplot_df %>% arrange(reg_exp)
pie_temp2 <- barplot_df_SQTQ %>% arrange(reg_exp)
pie_temp$n_SQTQ <- pie_temp2$n
pie_temp$fraction_SQTQ <- pie_temp$n_SQTQ/pie_temp$n
pie_temp$fraction_nonSQTQ <- 1-pie_temp$fraction_SQTQ

ratios_for_pieplot <- data.frame(as.list(pie_temp$fraction_SQTQ)) %>% t %>% as.data.frame()
ratios_for_pieplot$V2 <- 1-ratios_for_pieplot$V1
ratios_for_pieplot$experiment <- rownames(ratios_for_pieplot)

# plot for four ATRi and GEM
ggplot(pie_temp %>% pivot_longer(c(fraction_nonSQTQ, fraction_SQTQ)), aes(x=1, y=value, fill = name)) +
  geom_col() +
  coord_polar("y", start=0)+
  scale_fill_manual(values=c("grey70", "red"))+
  theme_void()+
  theme(text=element_text(size=12, color = 'black'),
        legend.key.size = unit(0.2,"line"),)+
  facet_wrap('reg_exp', ncol = 3) 




############################################################
# Supplementary Fig. 3c (barplots)
############################################################

# plot all p-sites 
barplot_df <- curves %>% dplyr::filter(ATRi_count_regulated >0) %>% dplyr::select(ATRi_count_regulated, Gene_SitePos) %>% distinct %>% group_by(ATRi_count_regulated) %>% dplyr::summarise(n=n())%>% arrange(desc(n))
barplot_df <- rbind(barplot_df, setNames(data.frame("3_4", barplot_df$n[3]+barplot_df$n[4]), names(barplot_df))) %>% dplyr::filter(ATRi_count_regulated %in% c(1,2, "3_4"))

ggplot(barplot_df, aes(y=n, x=ATRi_count_regulated)) + 
  geom_bar(position="stack", stat="identity", fill="grey80", color = 'black')+ 
  theme(legend.title=element_blank())+ scale_y_continuous(expand = c(0, 0), limits = c(0,1500)) +
  geom_text(aes(label=n), vjust=0, angle = 90, hjust = 0)

# plot SQTQ p-sites only
barplot_df <- curves %>% dplyr::filter(ATRi_count_regulated >0 & SQTQ==T) %>% dplyr::select(ATRi_count_regulated, Gene_SitePos) %>% distinct %>% group_by(ATRi_count_regulated) %>% dplyr::summarise(n=n())%>% arrange(desc(n))
barplot_df <- rbind(barplot_df, setNames(data.frame("3_4", barplot_df$n[3]+barplot_df$n[4]), names(barplot_df))) %>% dplyr::filter(ATRi_count_regulated %in% c(1,2, "3_4"))

ggplot(barplot_df, aes(y=n, x=ATRi_count_regulated)) + 
  geom_bar(position="stack", stat="identity", fill="grey80", color = 'black')+ 
  theme(legend.title=element_blank())+ scale_y_continuous(expand = c(0, 0), limits = c(0,100)) +
  geom_text(aes(label=n), vjust=0, angle = 90, hjust = 0)




############################################################
# Supplementary Fig. 4c (Venn diagrams)
############################################################

# plot all p-sites 
venn_temp <- curves %>% dplyr::select(Gene_SitePos, ATRi_updown_atleast3, GEM_sign)%>% unique #dplyr::filter(ATRi_count_regulated >0) %>% unique
venn_temp$GEM_sign[is.na(venn_temp$GEM_sign)] <- FALSE
venn_temp %>% group_by(ATRi_updown_atleast3, GEM_sign) %>% count

set1 <- venn_temp %>% dplyr::filter(ATRi_updown_atleast3 != "") %>% dplyr::pull(`Gene_SitePos`)
set2 <- venn_temp %>% dplyr::filter(GEM_sign ==T) %>% dplyr::pull(`Gene_SitePos`)

dev.off()
venn_plot <- venn.diagram(
  x = list(set1, set2),
  category.names = c("atleast3ATRi", "GEM"),
  filename = NULL,
  output=T,
  lwd = 1.5,
  lty = 'blank',
  fill = c( '#3070B3','grey'),
  cex = 2,
  fontface = "bold",
  fontfamily = "sans",
  cat.cex = 0.6,
  cat.default.pos = "outer",
  cat.fontfamily = "sans"
)
grid.draw(venn_plot)


# plot SQTQ p-sites only
venn_temp <- curves %>% dplyr::filter(SQTQ ==T) %>% dplyr::select(Gene_SitePos, ATRi_updown_atleast3, GEM_sign)%>% unique #dplyr::filter(ATRi_count_regulated >0) %>% unique
venn_temp$GEM_sign[is.na(venn_temp$GEM_sign)] <- FALSE
venn_temp %>% group_by(ATRi_updown_atleast3, GEM_sign) %>% count

set1 <- venn_temp %>% dplyr::filter(ATRi_updown_atleast3 != "") %>% dplyr::pull(`Gene_SitePos`)
set2 <- venn_temp %>% dplyr::filter(GEM_sign ==T) %>% dplyr::pull(`Gene_SitePos`)

dev.off()
venn_plot <- venn.diagram(
  x = list(set1, set2),
  category.names = c("atleast3ATRi", "GEM"),
  filename = NULL,
  output=T,
  lwd = 1.5,
  lty = 'blank',
  fill =  c( '#3070B3','grey'),
  cex = 2,
  fontface = "bold",
  fontfamily = "sans",
  cat.cex = 0.6,
  cat.default.pos = "outer",
  cat.fontfamily = "sans"
)
grid.draw(venn_plot)
dev.off()




############################################################
# Fig. 4a (scatterplot, pEC50 vs FC four ATRi)
############################################################

curves$Experiment <- factor(curves$Experiment, levels = c("Elimusertib",  "Gartisertib" ,  "Berzosertib", "Ceralasertib"))

ggplot(curves %>% filter(ATRi_regulated == T & SQTQ ==F), aes(y=pEC50, x=`Curve Fold Change`)) +
  geom_point(size = 0.75, color='grey70', pch = 16)+
  geom_point(size = 0.75, data = curves %>% filter(ATRi_regulated == T & SQTQ ==T), aes(y=pEC50, x=`Curve Fold Change`), color='red', pch = 16)+
  scale_x_continuous(limits=c(-6.6,6.6), breaks = seq(-6,6,2))+
  scale_y_continuous(limits=c(5,9))+
  theme_classic()+
  theme(text=element_text(size=6, color = 'black'),
    axis.text.x=element_text(angle = 0, hjust = 0.5),
    axis.line = element_line(colour = 'black', size = 0.1), 
    axis.ticks = element_line(size = 0.1),
    axis.ticks.length=unit(.05, "cm"),
    panel.border = element_rect(colour = "black", fill=NA, size=0,1),
    strip.background = element_blank())+
  geom_text_repel(data = curves %>% dplyr::filter(Gene_SitePos %in% c("CHEK1_S317", "BRCA1_S1239","FANCD2_S319",  "TOPBP1_S1504")),
                  aes(label=Gene_SitePos),
                  nudge_x = -4.5,
                  size          = 1,
                  box.padding   = 0.5,
                  force         = 100,
                  segment.size  = 0.2,
                  direction     = "y",
                  hjust = "left")+
  facet_wrap('Experiment', nrow=2, scales='free')




############################################################
# Fig. 4c (pEC50 distribution four ATRi)
############################################################

col_assignment <- c("#005293", "#999999", "grey30", "#64a0c8")
names(col_assignment) <- c( "Elimusertib" , "Ceralasertib" ,"Berzosertib", "Gartisertib" )

curves %>% filter(ATRi_regulated ==T & SQTQ == T) %>% group_by(Experiment) %>% tally
ggplot(curves %>% filter(ATRi_regulated ==T & SQTQ == T), aes(x=pEC50, colour = Experiment)) +
  geom_density(alpha = 0.5,size = 1)+
  theme_classic()+
  theme(legend.position=c(0.8,0.85), legend.title=element_blank(),
        legend.key.size = unit(0.5,"line"),
        text=element_text(size=10, color = 'black'),
        axis.text.x=element_text(angle = 0, hjust = 0.5, size=10),
        axis.text.y=element_text( size=10),
        axis.line = element_line(colour = 'black', size = 0.1), 
        axis.ticks = element_line(size = 0.1),
        axis.ticks.length=unit(.05, "cm"),
        panel.border = element_rect(colour = "black", fill=NA, size=0,1))+
  xlab("pEC50 [M]")+ ylab("Density")+ 
  scale_color_manual(values=col_assignment) +
  scale_x_reverse(expand = c(0,0),limits = c(10.2,3.8),  breaks = seq(4,10, by = 1))+
  scale_y_continuous(expand = c(0, 0), breaks = seq(0, 1, by = 0.5), limits = c(0,1.2))




############################################################
# Fig. 5a (scatterplot, pEC50 vs FC, GEM and >= 3 ATRi)
############################################################

ggplot(curves %>% filter(GEM_sign ==T & SQTQ ==F & ATRi_count_regulated <= 2) %>% distinct(GEM_qvalue, GEM_log2FC, Gene_SitePos),
       aes(y=-log10(GEM_qvalue), x=GEM_log2FC)) +
  theme_classic()+
  geom_hline(yintercept=-log10(0.01), col = 'black', size = 0.1)+
  geom_vline(xintercept=c(-1,1), col = 'black', size = 0.1)+
  geom_point(size = 1.5, color='grey80', pch = 16)+
  geom_point(size = 1.5, data = curves %>% filter(GEM_sign ==T & SQTQ ==T & ATRi_count_regulated <= 2) %>% distinct(GEM_qvalue, GEM_log2FC, Gene_SitePos), 
             aes(y=-log10(GEM_qvalue), x=GEM_log2FC), color='#f9a19a', pch = 16)+
  scale_x_continuous(limits=c(-5,5), breaks = seq(-5,5,1))+
  scale_y_continuous(limits=c(1.5,8), breaks = seq(2,8,1))+
  theme(text=element_text(size=8, color = 'black'),
        axis.text.x=element_text(angle = 0, hjust = 0.5, size=10),
        axis.text.y=element_text( size=10),
        axis.line = element_line(colour = 'black', size = 0.1), 
        axis.ticks = element_line(size = 0.1),
        axis.ticks.length=unit(.05, "cm"),
        panel.border = element_rect(colour = "black", fill=NA, size=0,1))+
  geom_point(curves %>% filter(GEM_sign ==T & ATRi_count_regulated >2 & SQTQ ==F) %>% distinct(GEM_qvalue, GEM_log2FC, Gene_SitePos), 
             mapping = aes(y=-log10(GEM_qvalue), x=GEM_log2FC),  
             size = 1.5,    color = 'grey40', pch = 16)+
  geom_point(curves %>% filter(GEM_sign ==T & ATRi_count_regulated >2 & SQTQ ==T) %>% distinct(GEM_qvalue, GEM_log2FC, Gene_SitePos), 
             mapping = aes(y=-log10(GEM_qvalue), x=GEM_log2FC),  
             size = 1.5,   color = 'red', pch = 16) +
  geom_text_repel(data = curves %>% 
                    dplyr::filter(Gene_SitePos %in% c("CHEK1_S317", "CHEK1_S468", "FANCD2_S319", "TOPBP1_S1504", "NBN_S397", "NBN_S615","UIMC1_S171", "BRCA1_1239", "EXO1_S714", "H2AFX_S139")) %>%
                    distinct(GEM_qvalue, GEM_log2FC, Gene_SitePos),
                  aes(label=Gene_SitePos),
                  nudge_x = 4.5,
                  size          = 2,
                  box.padding   = 0.5,
                  force         = 100,
                  segment.size  = 0.2,
                  direction     = "y",
                  hjust = "right")+
  geom_text_repel(data = curves %>% 
                    dplyr::filter(Gene_SitePos %in% c("HIST1H1E_T146", "HIST1H1D_T147", "HIST1H1C_T146", "HIST1H1D_T18", "HIST1H1E_T18", "MYBL2_T494","MKI67_T761","NIFK_T223", "PRC1_T481", "TPX_T72")) %>%
                    distinct(GEM_qvalue, GEM_log2FC, Gene_SitePos),
                  aes(label=Gene_SitePos),
                  nudge_x = -4.5,
                  size          = 2,
                  box.padding   = 0.5,
                  force         = 100,
                  segment.size  = 0.2,
                  direction     = "y",
                  hjust = "left")




############################################################
# Supplementary Fig. 4b (Heatmap)
############################################################

heatmap_df <- curves %>% dplyr::filter(combi_regulation == 'GEMup_ATRidown' & SQTQ ==T) 
hotlist <- unique(heatmap_df$Gene_SitePos)
heatmap_df <-  heatmap_df %>% dplyr::filter(Gene_SitePos %in% hotlist)
heatmap_df$GEM <- heatmap_df$GEM_ratio_mean
heatmap_df$GEM <- 2^heatmap_df$GEM_log2FC
heatmap_df$vehicle <- 1
heatmap_df$ATRi <- 2^(heatmap_df$`Curve Fold Change`)*heatmap_df$GEM

regulated_annotation <-  heatmap_df %>% dplyr::select(Experiment, Gene_SitePos, ATRi_regulated)
regulated_annotation$regulated <- ifelse(regulated_annotation$ATRi_regulated ==F, "*", "")

heatmap_df <- heatmap_df %>% dplyr::select(vehicle, GEM, ATRi, Experiment, Gene_SitePos) %>% 
  pivot_wider(names_from = Experiment, values_from = ATRi) %>%  
  pivot_longer(-c(Gene_SitePos), names_to = "Experiment", values_to = "values")

heatmap_df <- merge(heatmap_df, regulated_annotation, by=c('Experiment', 'Gene_SitePos'), all.x = T)
heatmap_df$Gene_SitePos <- gsub("_", "-",heatmap_df$Gene_SitePos)
Gene_sorted <- heatmap_df %>% dplyr::filter(Experiment == "GEM") %>% arrange((values)) %>% pull(Gene_SitePos)
Experiments_sorted <- c("vehicle", "GEM", "Elimusertib", "Gartisertib", "Berzosertib", "Ceralasertib")

ggplot(heatmap_df, aes(factor(Experiment, levels = Experiments_sorted), factor(Gene_SitePos, levels= Gene_sorted))) +
  geom_tile(aes(fill = log2(values)), colour = "black") +
  theme_minimal()+
  scale_fill_gradient2(low = "blue",
                       mid = "white",
                       high = "#FF0000",midpoint = 0) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), #, panel.border = element_rect(colour = "grey", fill=NA, size=1), axis.line = element_line(colour = "black"),
        axis.title = element_blank(),axis.text = element_text(size =10, color="black"), 
        legend.text= element_text(size = 12, color="black"), 
        legend.title = element_text(size = 12, face = "bold"))+ 
  coord_fixed()+
  geom_text(aes(label=regulated)) 





############################################################
# Fig. 4b, Fig. 5b, and Supplementary Fig. 4d (decryptM PTM response curves)
############################################################


col_assignment <- c("#005293", "#999999", "grey30", "#64a0c8")
names(col_assignment) <- c( "Elimusertib" , "Ceralasertib" ,"Berzosertib", "Gartisertib" )

symbols_assignment <- c(16,16,16,16)
names(symbols_assignment) <- c( "Elimusertib" , "Ceralasertib" ,"Berzosertib", "Gartisertib" )

# Function for curve fit: B is vector with four parameters. x is many x values
LL.4.function <- function(B, x) {
  (B[2]-B[3])/(1 + 10^(B[4] * (x - log10(B[1])))) + B[3]
}

xlim.max <- 10000*1.05
xlim.min <- 0.5
doses = c(0,  1, 3, 10, 30, 100, 300, 1000, 3000, 10000)
doses_rev <- rev(doses)

curves_targets <- curves %>% dplyr::filter(combi_regulation == 'GEMup_ATRidown' & SQTQ==T) 
peptides <- unique(curves_targets$Gene_SitePos)
proteins <- unique(curves_targets$`Gene names`)

nf <- layout(matrix(c(1:8), ncol = 2, nrow = 4, byrow=TRUE), 
             widths=c(1), heights =c(1,1,1,1))

for(k in 1:length(peptides)){
  df_plots <- curves %>% dplyr::filter(Gene_SitePos %in% peptides[k]) 
  df_plots_ratios <- df_plots %>% dplyr::select(contains("Ratio "), Experiment) #, starts_with("Curve")
  colnames(df_plots_ratios) <- c(rev(doses), "Experiment")
  df_plots_ratios <- df_plots_ratios %>% pivot_longer(cols = colnames(df_plots_ratios)[1:length(doses)])
  colnames(df_plots_ratios) <- c("Experiment", "dose", "response")
  
  legendtext <- ""
  
  if(sum(is.na(df_plots_ratios$response))<nrow(df_plots_ratios)){
    
    ylim_max <- 1
    ylim <- c(0,1.3)
    ylabels <- seq(0, 1, by = 0.5)
    drugs <- as.character(unique(df_plots_ratios$Experiment)) %>% sort
    
    par(tck = -.025, mgp = c(5, 0.8, 0), mar =c(5.5,5.5,6,8.5))
    plot(1, type="n", xlab="", ylab="", xlim = c(xlim.min, xlim.max), ylim = ylim,  log = 'x',
         axes = F, legend = T, main = "", cex.main = 1.2)
    abline(h = 1/(2^df_plots$GEM_log2FC[1]), col= "red")
    
    for(n in 1:length(drugs)){
      legendtext[n] <- drugs[n]
      cols <- col_assignment[drugs[n]]
      symbols <- symbols_assignment[drugs[n]]
      
      add = TRUE
      df_plots_ratios_temp <- df_plots_ratios %>% dplyr::filter(Experiment %in% drugs[n])
      df_plots_ratios_temp$dose[df_plots_ratios_temp$dose == 0] <- "0.03"
      
      points(x = df_plots_ratios_temp$dose %>% as.numeric, y = (df_plots_ratios_temp$response), col = cols, 
             cex = 2, lwd = 2, 
             bg = cols,
             lty =1,
             pch = symbols,
             xlim = c(xlim.min, xlim.max), ylim = ylim)
      
      df_plots_param <- df_plots %>% dplyr::filter(Experiment %in% drugs[n])
      
      if(!is.na(df_plots_param$`Curve R2`)){
        four.parameters <- df_plots_param %>% dplyr::select(c('pEC50','Curve Front', 'Curve Back', 'Curve Slope'))
        four.parameters[1,1] <- 10^(-four.parameters[1,1])*1000000000
        four.parameters <- t(four.parameters) # Make numeric vector out of it. Leave out R2
        many.x.values <- c(seq(0.5, 1, by = 0.05), seq(1, 10, by = 0.1),seq(10, 100, by = 1), seq(100, 1000, by = 2),seq(1000, 13000, by = 100))
        y.fitted.curve <- LL.4.function(four.parameters, log10(many.x.values))
        lines(many.x.values, (y.fitted.curve), type = "l", log = "x", col = cols, lwd = 3)
      }
      
    }
    mtext(side=3, line=2, cex=1, gsub("_", "-p",unique(as.character(df_plots_param$`Gene_SitePos`))), font = 2)
    mtext(side=3, line=0.6, cex=0.8, unique(as.character(df_plots_param$`ModSeqShort`)))
    
    box(bty="l")
    unique_doses <- df_plots_ratios_temp$dose
    axis(1, at =unique_doses[c(1,3,5,7,9,10)], labels=c(10000,1000,100,10,1, 0), cex.axis =1.2)
    axis(2, at =ylabels, labels=ylabels, cex.axis =1.2, las = 2)
    title(xlab="1 ??M GEM + ATRi [nM]", ylab="PTM Response\n normalized to 1 ??M GEM", mgp = c(2.6, 1, 0), cex.lab = 1.2)
    
    legend('topleft',inset =c(1.02,0), legend = c(legendtext, "Vehicle"),
           plot = T, xpd = T, title = as.expression(bquote(bold("ATR inhibitor"))),
           lty =1, cex=1, pch = c(symbols_assignment[drugs],NA), bty="n", horiz = F, col = c(col_assignment[drugs], "red"))
    }
}
