
# This function removes samples that are not part of the original Roux et al paper.
select_only_roux_samples<-function(data){
  tmp<- data %>%
    mutate(station_number=gsub('Station([0-9]+)_.*', '\\1', sample, perl=TRUE ))
  
  tmp$station_number[tmp$station_number=='L4']<-0
  tmp$station_number<-as.numeric(tmp$station_number)
  tmp<- tmp %>%
    filter(station_number<=138) %>%
    dplyr::select(-station_number)
  return(tmp)
}

#set the colors for the contig type
get_contig_type_colors<-function(){
  return(c('L4 hybrid'="#ffd700", 
           'L4 short'="#ffb14e", 
           'L4 pilon'="#fa8775",
           'GOV'="#5fadea",
           'fosmid'="#003697",
           'RefSeq'="#9d02d7",
           'vSAG'="#774be5",
           'Luo 2017'='#000000'))
}


add_contig_types<-function(data){
  fosmid_ids<-read.table('../../data/L4.virome/abundance/fosmid.ids')
  return(data %>% mutate(contig_type=ifelse(grepl('pilon', contig), 'L4 pilon', 
                            ifelse(grepl('H_NODE', contig), 'L4 hybrid',
                                   ifelse(grepl('S_NODE', contig), 'L4 short',
                                          ifelse(grepl('vSAG', contig), 'vSAG', 
                                                 ifelse(grepl('GOV', contig), 'GOV', 
                                                        ifelse(grepl('Tp1', contig), 'GOV', 
                                                               ifelse(grepl('Mp1', contig), 'GOV',
                                                                      ifelse(grepl('luo_2017', contig), 'Luo 2017',
                                                                      ifelse(contig %in% fosmid_ids$V1, 'fosmid', 'RefSeq')))))))))))
}



#This returns a matrix of abundances
get_abundance_matrix<-function(sample_type='ALL', rank_per_sample=50, min_pct_samples=10){
  abundance_mtx<-abundance_df %>%
    select_only_roux_samples() %>%
    mutate(contig=gsub('VIRSorter_', '', contig, perl=TRUE)) %>%
    mutate(contig=gsub('_cov.*', '', contig, perl=TRUE)) %>%
    filter(rank.per.sample<=rank_per_sample) %>%
    mutate(sample=gsub('Station', '', sample, perl=TRUE)) %>%
    filter(grepl(ifelse(sample_type =='ALL', '.*',
                        paste(sample_type, '|L4', sep='')), sample)) %>%
    dplyr::select(contig, sample, estimated_abundance) %>%
    spread(contig, estimated_abundance,fill=0)
  
  abundance_mtx<-data.frame(abundance_mtx)
  rownames(abundance_mtx)<-abundance_mtx$sample
  
  abundance_mtx<-abundance_mtx %>% 
    dplyr::select(-sample) %>% 
    as.matrix()
  
  contigs_to_choose<-names(which(colSums((abundance_mtx>0))>nrow(abundance_mtx)/min_pct_samples))
  
  
  
  abundance_mtx_filtered<- as.data.frame(abundance_mtx) %>% 
    dplyr::select(contigs_to_choose) %>%
    as.matrix()
  
  return(abundance_mtx_filtered)
}



## This function plots a heatmap of samples
plot_heatmap<-function(sample_type='ALL', km=1, rank_per_sample=50, min_pct_samples=10){
  
  abundance_mtx_filtered<- get_abundance_matrix(sample_type = sample_type, rank_per_sample = rank_per_sample, min_pct_samples = min_pct_samples)
  
  color_df<-data.frame(get_contig_type_colors()) %>% set_colnames(c('contig_type_color'))
  color_df$contig_type<-rownames(color_df)
  
  barplot_df<-data.frame(contig=colnames(abundance_mtx_filtered)) %>%
    add_contig_types() %>%
    mutate(totals=colSums(abundance_mtx_filtered)) %>%
    left_join(color_df, by='contig_type')
  
  barplot_df$contig_type<-as.factor(barplot_df$contig_type)
  barplot_df$contig_type_color<-as.character(barplot_df$contig_type_color)
  
  
  ha1 = HeatmapAnnotation(dist1 = anno_barplot(barplot_df$totals, bar_width = 1, 
                                               gp = gpar(col = NA, 
                                                        fill = barplot_df$contig_type_color), 
                                               border = FALSE, axis = TRUE, axis_side='right'), show_legend = FALSE)
  
  p<-Heatmap(abundance_mtx_filtered, name='Abundance',
             col = colorRamp2(c(0, 10000, 20000, 30000), c("white", "cornflowerblue", "yellow", "red")),
             row_title=paste(sample_type, 'Samples'),
             row_names_gp = gpar(fontsize = 7),
             column_names_gp = gpar(fontsize =7),
             show_column_dend = FALSE,
             row_dend_side = "right",
             row_dend_width = unit(2, "cm"),
             km=km,
             top_annotation = ha1, 
             top_annotation_height = unit(1.5, "cm"),
             cluster_columns = FALSE,
             column_order=rev(order(colSums(abundance_mtx_filtered))))
  
  
  return(p)
}

slideFunct<- function(data, name='test', window=500, step=1, cutoff=0.2){
  
  total <- length(data)
  spots <- seq(from=1, to=(total-window), by=step)
  result <- vector(length = length(spots))
  start <- vector(length = length(spots))
  end <- vector(length = length(spots))
  for(i in 1:length(spots)){
    start[i] <- i * step - step +1
    result[i] <- median(data[spots[i]:(spots[i]+window)])
    end[i] <- i * step - step + window
  }
  median_value=median(data)
  df<-data.frame('contig'=name, 
                 'start'=start, 'end'=end,
                 'window_median_coverage'=result,
                 'contig_median_coverage'=median_value,
                 'potential_GI'= result < median_value * cutoff)
  return(df)
}


analyse_GIs<-function(x){
  tmp<-suppressMessages(read_csv(x))
  
  return(list('count'=sum(rle(tmp$potential_GI)$values),
              'total_length'=sum(
                rle(tmp$potential_GI)$lengths[which(
                  rle(tmp$potential_GI)$values)]),
              'max_length'=max(
                rle(tmp$potential_GI)$lengths[which(
                  rle(tmp$potential_GI)$values)])))
}

calculate_GI_lengths<-function(x){
  tmp<-suppressMessages(read_csv(x))
  return(sum(
    rle(tmp$potential_GI)$lengths[which(
      rle(tmp$potential_GI)$values)]))
}


