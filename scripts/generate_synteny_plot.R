suppressPackageStartupMessages(library("argparse"))
suppressPackageStartupMessages(library("dplyr"))

parser <- ArgumentParser()
parser$add_argument("-busco_list", nargs = '+',
                    help = "A list of Busco \"full_table.tsv\" of the two or more query species ")
parser$add_argument("-chrom_list", nargs = '+',
                    help = "A list of .tsv tables with chromosomes of the two or more query species (5 columns expected chromosome, length, order, direction, annotation)")
parser$add_argument("-o", "--output_prefix", default = "synteny_plot",
                    help = "Name pattern for the output")
parser$add_argument("-g", "-gap", type = "integer", default = 6,
                    help = "Gap between two chromosomal sets")
parser$add_argument("-f", "-filter", type = "integer", default = 5,
                    help = "The minimal number of BUSCOs on a chromosome to include")
parser$add_argument("-alpha", type = "integer", default = 0,
                    help = "Set transparency to colours [%]")
args <- parser$parse_args()

#source('scripts/interactive_args.R')

### specify arguments and parameters ###
#busco_list = args$busco_list
busco_list <- readLines(file("/Users/ab66/Documents/sanger_work/Tools/synteny_plotter/busco_list.tsv"))
#chrom_list = args$chrom_list
chrom_list <- readLines(file("/Users/ab66/Documents/sanger_work/Tools/synteny_plotter/chrom_list.tsv"))

#minimum_buscos <- args$f
minimum_buscos <- 5

print('[+] Processing list of file(s):')
print(args$busco_list)

### import functions ###
#source('scripts/helper_functions.R')
source('/Users/ab66/Documents/sanger_work/Tools/synteny_plotter/scripts/helper_functions.R')

### read in data ###

# TODO: allow for algs in the arguments
# algs <- read.csv(args$alg_file, sep='\t', header=FALSE)[,c(1,3)]
# colnames(algs) <- c('busco', 'alg')
# alignments <- merge(alignments, algs, by='busco')

# Read ref files
ref_df <- read_buscos(busco_list[1], 'R')
ref_chroms <- read.table(chrom_list[1], sep = '\t', header = TRUE)
ref_chroms <- ref_chroms %>% arrange(order)

chr_offset <- max(ref_chroms$length) / 2 # make chr offset 50% of the largest chr size in the ref chr set
# NB: using ref genome as a proxy for chr lengths in other genomes too

### generate alignments ###

make_alignment_table <- function(R_df, R_chroms, Q_df, Q_chroms){
  alignments <- merge(Q_df, R_df, by = 'busco')
  alignments <- alignments %>% group_by(chrR) %>% filter(n() > minimum_buscos) %>% ungroup()
  # apply any filters
  # ref chromosomes
  R_chroms <- R_chroms %>% filter(chr %in% alignments$chrR)
  R_chroms <- R_chroms %>% arrange(order)
  chr_order_R <- R_chroms$chr # extract order of chr
  
  #query chromosomes
  Q_chroms <- Q_chroms %>% filter(chr %in% alignments$chrQ)
  Q_chroms <- Q_chroms %>% arrange(order)
  chr_order_Q <- Q_chroms$chr #extract order of chr
  
  alignments$alg <- NA
  alignments <- perform_inverts(alignments, Q_chroms)
  offset_alignments_Q <- offset_chr(alignments, 'Q', chr_offset, chr_order_Q)
  offset_alignments_RQ <- offset_chr(offset_alignments_Q$df, 'R', chr_offset, chr_order_R)
  alignments <- offset_alignments_RQ$df
  offset_list_R <- offset_alignments_RQ$offset_list
  offset_list_Q <- offset_alignments_Q$offset_list
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
  processed_Q <- make_alignment_table(temp_ref_df, temp_ref_chroms, query_df, query_chroms)
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
if (nrow(ref_chroms) <= 9){
  col_list <- c("#FFC759","#FF7B9C", "#607196", "#BABFD1", '#BACDB0', '#C6E2E9', '#F3D8C7', '#47A8BD','#FFAD69')
} else { 
  col_list <- c('#577590', '#617A8B', '#6B8086', '#748581', '#7E8A7C', '#889077', '#929572', '#9B9A6D',
                '#A5A068', '#AFA563', '#B9AA5E', '#C2AF59', '#CCB554', '#D6BA4F', '#E0BF4A', '#E9C545', 
                '#F3CA40', '#F1C544', '#EFC148', '#EDBC4C', '#EBB84F', '#E9B353', '#E7AF57', '#E5AA5B',
                '#E3A55F', '#E1A163', '#DF9C67', '#DD986A', '#DB936E', '#D98F72', '#D78A76')
  num_remove <- length(col_list) - nrow(ref_chroms)
  if (num_remove > 0){
    set.seed(123) # set a random seed for reproducibility
    indices_to_remove <- sample(length(col_list), num_remove) # get random indices of elements to remove
    my_list_cleaned <- col_list[-indices_to_remove] # remove elements from the list using negative indexing
    col_list <- my_list_cleaned
  }
}

# for more than 31 chromosomes
if (nrow(ref_chroms) > 31){
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
ref_chroms$chr_colour <- col_list_final
busco2colour <- ref_chroms[,c(1,6)]
temp <- ref_df[,c(1,2)]
colnames(temp) <- c('busco', 'chr')
busco2colour <- merge(busco2colour, temp, by='chr') # assigns a colour to each busco based on which chr its in in ref

max_end <- max(unlist(max_ends))
plot_length <- max_end # make plot_length the max of the longest chr set
#gap <- args$gap
gap <- 6
alpha = 0.6
show_outline = TRUE

#pdf(paste0(args$output_prefix, '.pdf'))
pdf(paste0('/Users/ab66/Documents/sanger_work/Tools/synteny_plotter/test', '.pdf'))
print('[+] Generating plot')
plot(0,cex = 0, xlim = c(1, plot_length), 
     ylim = c(((gap+1)*-1*length(busco_list)*2),((gap+1)*length(busco_list)*2)),
     #ylim = c(((gap+1)*-1*2*2),((gap+1)*2*2)),
     xlab = "", ylab = "", bty = "n", yaxt="n", xaxt="n")

main_counter <- 1
y_offset <- 0
y_increment <- 17

for (i in busco_list[-1]){
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
  
  # plot alignments
  counter <- 1
  for (i in chr_order_R){
    temp <- alignments[alignments$chrR == i,]
    plot_one_ref_chr(temp, adjustment_length_R, adjustment_length_Q, y_offset, busco2colour, alpha)
    counter <- counter + 1
  }
  
  counter <- 1  # chr outlines for query
  for (i in chr_order_Q){
    temp <- alignments[alignments$chrQ == i,]
    Qfirst <- min(temp$Qstart)
    Qlast <- max(temp$Qend)
    segments(Qfirst+adjustment_length_Q, 1-gap-y_offset, 
             Qlast+adjustment_length_Q, 1-gap-y_offset, lwd = 5)
  }
  
  counter <- 1   # chr outlines for ref
  for (i in chr_order_R){
    temp <- alignments[alignments$chrR == i,]
    Rfirst <- min(temp$Rstart)
    Rlast <- max(temp$Rend)
    segments(Rfirst+adjustment_length_R, gap-y_offset, 
             Rlast+adjustment_length_R, gap-y_offset, lwd = 5)
  }
  
  # add text labels
  if (main_counter == 1){ # only need to plot text labels ref if this is the first chr set being plotted, else get duplicates
    counter <- 1    # text labels for ref
    for (i in chr_order_R){
      temp <- alignments[alignments$chrR == i,]
      Rfirst <- min(temp$Rstart)
      Rlast <- max(temp$Rend)
      offset <- offset_list_R[[counter]]
      text(x = ((Rlast+Rfirst+1)/2)+adjustment_length_R, y = gap+3.5, 
           label = i,
           srt = 90, cex=0.5) # Rotation
      counter <- counter + 1
    }
  }
  
  counter <- 1  # text labels for query:
  for (i in chr_order_Q){
    temp <- alignments[alignments$chrQ == i,]
    Qfirst <- min(temp$Qstart)
    Qlast <- max(temp$Qend)
    text(x = ((Qlast+Qfirst+1)/2)+adjustment_length_Q, y = -gap-2-y_offset, label = i,
         srt = 90, cex=0.5) # Rotation
    counter <- counter + 1
  }
  
  main_counter <- main_counter + 5
  y_offset <- y_offset + y_increment
}

dev.off()

q()