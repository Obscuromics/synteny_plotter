suppressPackageStartupMessages(library("argparse"))
suppressPackageStartupMessages(library("dplyr"))

parser <- ArgumentParser()
parser$add_argument("--busco_list", nargs = '+',
                    help = "A list of Busco \"full_table.tsv\" of the two or more query species ")
parser$add_argument("--chrom_list", nargs = '+',
                    help = "A list of .tsv tables with chromosomes of the two or more query species (5 columns expected chromosome, length, order, direction, annotation)")
parser$add_argument("-o", "--output_prefix", default = "synteny_plot",
                    help = "Name pattern for the output")
parser$add_argument("-g", "--gap", type = "integer", default = 6,
                    help = "Gap between two chromosomal sets")
parser$add_argument("-f", "--filter", type = "integer", default = 5,
                    help = "The minimal number of BUSCOs on a chromosome to include")
parser$add_argument("--alpha", type = "integer", default = 0,
                    help = "Set transparency to colours [%]")
parser$add_argument("--colour_by", default = "chromosomes", help = "chromosomes, algs")
parser$add_argument("--alg_file", default = NULL, help = "List of ALGs")
parser$add_argument("--alg_colours", default = NULL, help = "A list of ALGs and assigned colours")

args <- parser$parse_args()

#source('scripts/interactive_args.R')

### specify arguments and parameters ###
busco_list <- readLines(file(args$busco_list)) 
#busco_list <- readLines(file("/Users/ab66/Documents/sanger_work/Tools/synteny_plotter/busco_list.tsv"))
chrom_list <- readLines(file(args$chrom_list))
#chrom_list <- readLines(file("/Users/ab66/Documents/sanger_work/Tools/synteny_plotter/chrom_list.tsv"))

minimum_buscos <- args$f
#minimum_buscos <- 5

print('[+] Processing list of file(s):')
print(busco_list)

### import functions ###
source('scripts/helper_functions.R')
#source('/Users/ab66/Documents/sanger_work/Tools/synteny_plotter/scripts/helper_functions.R')

### read ALGs ###
if (args$colour_by == "algs"){
  if(is.null(args$alg_file)){
    print("Warnin: provide ALG file with the option --alg_file")
  }else{
    algs <- read.csv(args$alg_file, sep='\t', header=FALSE)[,c(1,2,3)]
    colnames(algs) <- c('busco', 'status', 'alg')
    algs <- algs %>% filter(alg != "unassigned")
  }
}

# Read ref files
ref_df <- read_buscos(busco_list[1], 'R')
ref_chroms <- read.table(chrom_list[1], sep = '\t', header = TRUE)
ref_chroms <- ref_chroms %>% arrange(order)

chr_offset <- max(ref_chroms$length) / 2 # make chr offset 50% of the largest chr size in the ref chr set
# NB: using ref genome as a proxy for chr lengths in other genomes too

### generate alignments ###

make_alignment_table <- function(R_df, R_chroms, Q_df, Q_chroms, chr_offset, algs = NULL){
  alignments <- merge(Q_df, R_df, by = 'busco')
  alignments <- alignments %>% group_by(chrR) %>% filter(n() > minimum_buscos) %>% ungroup()
  
  # apply any filters
  # ref chromosomes
  R_chroms <- R_chroms %>% filter(chr %in% alignments$chrR)
  R_chroms <- R_chroms %>% arrange(order)
  chr_order_R <- R_chroms[,c("chr", "length")] # extract order and length of chr
  
  #query chromosomes
  Q_chroms <- Q_chroms %>% filter(chr %in% alignments$chrQ)
  Q_chroms <- Q_chroms %>% arrange(order)
  chr_order_Q <- Q_chroms[,c("chr", "length")] #extract order and length of chr
  
  alignments <- perform_inverts(alignments, Q_chroms)
  offset_alignments_Q <- offset_chr(alignments, 'Q', chr_offset, chr_order_Q)
  offset_alignments_RQ <- offset_chr(offset_alignments_Q$df, 'R', chr_offset, chr_order_R)
  alignments <- offset_alignments_RQ$df
  offset_list_R <- offset_alignments_RQ$offset_list
  offset_list_Q <- offset_alignments_Q$offset_list
  
  if(!is.null(args$alg_file)){
    alignments <- merge(alignments, algs, by='busco')
  }else{
    alignments$alg <- NA
  }
  
  output_list <- list('alignments' = alignments, 
                      'chr_order_R' = chr_order_R, 'chr_order_Q' = chr_order_Q,
                      'offset_list_R' = offset_list_R, 'offset_list_Q' = offset_list_Q)
  return(output_list)
}

processed_Q_list <- list()
max_ends <- list()

temp_ref_chroms <- ref_chroms
temp_ref_df <- ref_df

for (file in busco_list[-1]){
  #print(file)
  i <- match(file, busco_list)
  query_df <- read_buscos(file, 'Q')
  query_chroms <- read.table(chrom_list[i], sep = '\t', header = TRUE)
  processed_Q <- make_alignment_table(temp_ref_df, temp_ref_chroms, 
                                      query_df, query_chroms, chr_offset, algs = algs)
  alignments <- processed_Q$alignments
  processed_Q_list <- append(processed_Q_list, processed_Q)
  max_ends <- append(max_ends, max(alignments$Rend))
  max_ends <- append(max_ends, max(alignments$Qend))
  temp_ref_df <- query_df
  colnames(temp_ref_df) <- c('busco', 'chrR', 'Rstart', 'Rend', 'Rstrand')
  temp_ref_chroms <- query_chroms
}

### plotting ###

# define colours

colour_by <- args$colour_by
#colour_by = "algs"

if(colour_by == "chromosomes"){
  if(nrow(ref_chroms) <= 9){
    col_list <- c("#FFC759","#FF7B9C", "#607196", "#BABFD1", '#BACDB0', '#C6E2E9', '#F3D8C7', '#47A8BD','#FFAD69')
  }else{
    col_list <- c('#577590', '#617A8B', '#6B8086', '#748581', '#7E8A7C', '#889077', '#929572', '#9B9A6D',
                  '#A5A068', '#AFA563', '#B9AA5E', '#C2AF59', '#CCB554', '#D6BA4F', '#E0BF4A', '#E9C545', 
                  '#F3CA40', '#F1C544', '#EFC148', '#EDBC4C', '#EBB84F', '#E9B353', '#E7AF57', '#E5AA5B',
                  '#E3A55F', '#E1A163', '#DF9C67', '#DD986A', '#DB936E', '#D98F72', '#D78A76')
    num_remove <- length(col_list) - nrow(ref_chroms)
    if(num_remove > 0){
      set.seed(123) # set a random seed for reproducibility
      indices_to_remove <- sample(length(col_list), num_remove) # get random indices of elements to remove
      my_list_cleaned <- col_list[-indices_to_remove] # remove elements from the list using negative indexing
      col_list <- my_list_cleaned
    }
  }
  
  # for more than 31 chromosomes
  if(nrow(ref_chroms) > 31){
    col_list <- c('#577590', '#617A8B', '#6B8086', '#748581', '#7E8A7C', '#889077', '#929572', '#9B9A6D',
                  '#A5A068', '#AFA563', '#B9AA5E', '#C2AF59', '#CCB554', '#D6BA4F', '#E0BF4A', '#E9C545', 
                  '#F3CA40', '#F1C544', '#EFC148', '#EDBC4C', '#EBB84F', '#E9B353', '#E7AF57', '#E5AA5B',
                  '#E3A55F', '#E1A163', '#DF9C67', '#DD986A', '#DB936E', '#D98F72', '#D78A76',
                  '#577590', '#617A8B', '#6B8086', '#748581', '#7E8A7C', '#889077', '#929572', '#9B9A6D',
                  '#A5A068', '#AFA563', '#B9AA5E', '#C2AF59', '#CCB554', '#D6BA4F', '#E0BF4A', '#E9C545', 
                  '#F3CA40', '#F1C544', '#EFC148', '#EDBC4C', '#EBB84F', '#E9B353', '#E7AF57', '#E5AA5B',
                  '#E3A55F', '#E1A163', '#DF9C67', '#DD986A', '#DB936E', '#D98F72', '#D78A76',
                  '#577590', '#617A8B', '#6B8086', '#748581', '#7E8A7C', '#889077', '#929572', '#9B9A6D',
                  '#A5A068', '#AFA563', '#B9AA5E', '#C2AF59', '#CCB554', '#D6BA4F', '#E0BF4A', '#E9C545', 
                  '#F3CA40', '#F1C544', '#EFC148', '#EDBC4C', '#EBB84F', '#E9B353', '#E7AF57', '#E5AA5B',
                  '#E3A55F', '#E1A163', '#DF9C67', '#DD986A', '#DB936E', '#D98F72', '#D78A76',
                  '#577590', '#617A8B', '#6B8086', '#748581', '#7E8A7C', '#889077', '#929572', '#9B9A6D',
                  '#A5A068', '#AFA563', '#B9AA5E', '#C2AF59', '#CCB554', '#D6BA4F', '#E0BF4A', '#E9C545', 
                  '#F3CA40', '#F1C544', '#EFC148', '#EDBC4C', '#EBB84F', '#E9B353', '#E7AF57', '#E5AA5B',
                  '#E3A55F', '#E1A163', '#DF9C67', '#DD986A', '#DB936E', '#D98F72', '#D78A76')
  }
  
  col_list_final <- col_list[1:nrow(ref_chroms)] # subset col_list to number needed based on number of ref chromosomes
  ref_chroms$colour <- col_list_final
  busco2colour <- ref_chroms[,c(1,6)]
  temp <- ref_df[,c(1,2)]
  colnames(temp) <- c('busco', 'chr')
  busco2colour <- merge(busco2colour, temp, by = 'chr')
}

if(colour_by == "algs"){
  if(is.null(args$alg_colours)){
    algs2colour <- data.frame(alg = unique(algs$alg))
    col_list <- c("#882255", "#AA4499", "#CC6677", "#DDCC77", "#88CCEE", "#44AA99", "#117733", "#332288")
    col_list_final <- col_list[1:nrow(algs2colour)]
    algs2colour$colour <- col_list_final
  }else{
    algs2colour <- read.table(args$alg_colours, sep = '\t', header = FALSE)
    colnames(algs2colour) <- c("alg", "colour")
  }
  busco2colour <- merge(algs, algs2colour, by = 'alg')
}

max_end <- max(unlist(max_ends))
plot_length <- max_end # make plot_length the max of the longest chr set
gap <- args$gap
#gap <- 6
alpha = 0.6
show_outline = TRUE

ref_df <- read_buscos(busco_list[1], 'R')
ref_chroms <- read.table(chrom_list[1], sep = '\t', header = TRUE)
ref_chroms <- ref_chroms %>% arrange(order)

pdf(paste0(args$output_prefix, '.pdf'))
#pdf(paste0('/Users/ab66/Documents/sanger_work/Tools/test', '.pdf'))
print('[+] Generating plot')
plot(0,cex = 0, xlim = c(1, plot_length), 
     #ylim = c(((gap+1)*-1*length(busco_list)*2),((gap+1)*length(busco_list)*2)),
     ylim = c(((gap+1)*-1*2*2),((gap+1)*2*2)),
     xlab = "", ylab = "", bty = "n", yaxt="n", xaxt="n")

main_counter <- 1
y_offset <- 0
y_increment <- 11.65

for (file in busco_list[-1]){
  j <- match(file, busco_list)
  query_chroms <- read.table(chrom_list[j], sep = '\t', header = TRUE)
  
  alignments <- processed_Q_list[[main_counter]]
  chr_order_R <- processed_Q_list[[main_counter+1]]
  chr_order_Q <- processed_Q_list[[main_counter+2]]
  offset_list_R <- processed_Q_list[[main_counter+3]]
  offset_list_Q <- processed_Q_list[[main_counter+4]]
  
  if (max(alignments$Qend) != max_end){ # i.e. if this is the longest chr_set
    if (max(alignments$Rend) != max_end){
      adjustment_length_R <- (max_end - max(alignments$Rend)) / 2 
      adjustment_length_Q <- (max_end - max(alignments$Qend)) / 2 
    }
    else{
      adjustment_length_R <- 0
      adjustment_length_Q <- (max_end - max(alignments$Qend)) / 2 
    }
  }
  else{
    adjustment_length_Q <- 0
    adjustment_length_R <- (max_end - max(alignments$Rend)) / 2 
  }
  
  # align everything to the left
  #adjustment_length_Q <- 0
  #adjustment_length_R <- 0
  
  # plot alignments
  counter <- 1
  for (i in chr_order_R$chr){
    temp <- alignments[alignments$chrR == i,]
    plot_one_ref_chr(temp, adjustment_length_R, adjustment_length_Q, y_offset, busco2colour, alpha)
    counter <- counter + 1
  }
  
  # plotting query chromosomes
  if (main_counter == (length(processed_Q_list) - 4)){
    counter <- 0
    offset <- 0
    for (i in chr_order_Q$chr){
      chr_length <- chr_order_Q[chr_order_Q$chr == i,]$length
      Qfirst <- offset
      Qlast <- chr_order_Q[chr_order_Q$chr == i,]$length + offset
    
      if (counter != 0){ # only need to offset start/end if this is not the first chr
        Qfirst <- offset  # allows for accumulative chr positions
        Qlast <- chr_length + offset # allows for accumulative chr positions
      }
    
      offset <- offset + chr_length + chr_offset # accumulative offset
      counter <- counter + 1
    
      segments(Qfirst+adjustment_length_Q, 1-gap-y_offset, 
               Qlast+adjustment_length_Q, 1-gap-y_offset, lwd = 10)
    
      text(x = ((Qlast+Qfirst+1)/2)+adjustment_length_Q, y = 1-gap-y_offset, 
           label = query_chroms[query_chroms$chr == i,]$annot,
           srt = 0, cex = 0.5, col = "grey")
    
    }
  }
  
  # plotting reference chromosomes
  counter <- 0
  offset <- 0
  for (i in chr_order_R$chr){
    chr_length <- chr_order_R[chr_order_R$chr == i,]$length
    Rfirst <- offset
    Rlast <- chr_order_R[chr_order_R$chr == i,]$length + offset
    
    if (counter != 0){ # only need to offset start/end if this is not the first chr
      Rfirst <- offset  # allows for accumulative chr positions
      Rlast <- chr_length + offset # allows for accumulative chr positions
    }
    
    offset <- offset + chr_length + chr_offset # accumulative offset
    counter <- counter + 1
    
    segments(Rfirst+adjustment_length_R, gap-y_offset, 
             Rlast+adjustment_length_R, gap-y_offset, lwd = 10)
    
    text(x = ((Rlast+Rfirst+1)/2)+adjustment_length_R, y = gap-y_offset, 
         label = ref_chroms[ref_chroms$chr == i,]$annot,
         srt = 0, cex = 0.5, col = "grey")
    }
  
  main_counter <- main_counter + 5
  y_offset <- y_offset + y_increment
  ref_chroms <- query_chroms
}

dev.off()

q()