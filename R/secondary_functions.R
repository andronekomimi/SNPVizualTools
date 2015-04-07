
create_gr_from_df = function(DF,start_idx,end_idx,label) {
  #label = ""
  left = c()
  right = c()
  base_width = c()
  color = c()
  alpha = c()  
  
  for (i in (1:nrow(DF)))
  {
    start_value = as.numeric(DF[i,start_idx])
    end_value = as.numeric(DF[i,end_idx])
    
    if((start_value > (start(current_range) - 1) 
        && start_value < (end(current_range) - 1)) 
       || (end_value > (start(current_range) - 1) 
           && end_value < (end(current_range) - 1))) {
      
      left = c(left, start_value)
      right = c(right, end_value)
      base_width = c(base_width ,as.numeric(DF[i,start_idx+1])-start_value)
      color = c(color,label)
    }
  }
  
  ranges = IRanges(left,right)
  
  if(length(left) > 0) {
    g_ranges <- GRanges(seqnames = current_chr, ranges = ranges, size = base_width, color = color, alpha = rep(0.5,length(left)))
    print(paste0("-> Find ",length(g_ranges)," interaction for ", label))
    invisible(g_ranges)
  }
  else
  {
    print(paste0("No relation found for ", label))
    invisible(GRanges())
  }
}

create_gr_from_microTfile = function(DF,start_idx,end_idx,label) {
#   label = ""
  left = c()
  right = c()
  base_width = c()
  color = c()
  alpha = c()
  transcript_id = c()
  annotated = c()
  
  for (i in (1:nrow(DF)))
  {
    start_value = as.numeric(DF[i,start_idx])
    end_value = as.numeric(DF[i,end_idx])
    
    if((start_value > (start(current_range) - 1) 
        && start_value < (end(current_range) - 1)) 
       || (end_value > (start(current_range) - 1) 
           && end_value < (end(current_range) - 1))) {
      
      left = c(left, start_value)
      right = c(right, end_value)
      base_width = c(base_width ,as.numeric(DF[i,start_idx+1])-start_value)
      color = c(color,label)
      
      # transcript_id
      str = strsplit(DF[i,9], split = "; ")[[1]]
      transcript_id_idx = grep(str, pattern = "transcript_id")
      tId = gsub(pattern = "transcript_id ",x = str[transcript_id_idx], replacement = "")
      transcript_id = c(transcript_id, tId)
      
      # transcript annotation
      annotation_idx = grep(str, pattern = "tstatus")
      annotation = ! grepl(x = str[annotation_idx], pattern = "unannotated")
      annotated = c(annotated, annotation)
    }
  }
  
  ranges = IRanges(left,right)
  
  if(length(left) > 0) {
    g_ranges <- GRanges(seqnames = current_chr, ranges = ranges, 
                        base_width = base_width, color = color, alpha = rep(0.5,length(left)),
                        annotated = annotated, tsID = transcript_id)
    print(paste0("-> Find ",length(g_ranges)," interaction for ", label))
    invisible(g_ranges)
  }
  else
  {
    print(paste0("No relation found for ", label))
    invisible(GRanges())
  }
}

create_gr_from_df_4_circos_range = function(DF,label) {
  
  chrom_start_idx = 1
  left_start_idx = 2
  chrom_end_idx = 4
  right_start_idx = 5
  left_end_idx = 3
  right_end_idx = 6
  
  chrom_start = c()
  chrom_end =  c()
  left_start = c()
  left_end = c()
  right_start = c()
  right_end = c()
  base_width = c()
  color = c()
  
  for (i in (1:nrow(DF)))
  {
    chrom_start_value =  sub("chr","",DF[i,chrom_start_idx])
    chrom_end_value =  sub("chr","",DF[i,chrom_end_idx])
    
    left_start_value = as.numeric(DF[i,left_start_idx])
    left_end_value = as.numeric(DF[i,left_end_idx])
    
    right_start_value = as.numeric(DF[i,right_start_idx])
    right_end_value = as.numeric(DF[i,right_end_idx])
    
    if(chrom_start_value == as.numeric(sub("chr","",seqlevels(current_range)))) {
      if((left_start_value >= start(current_range) && left_start_value <= end(current_range)) || (left_end_value >= start(current_range) && left_end_value <= end(current_range)))
      {
        chrom_start = c(chrom_start, chrom_start_value)
        chrom_end = c(chrom_end, chrom_end_value)
        left_start = c(left_start, left_start_value)
        left_end = c(left_end, left_end_value)
        right_start = c(right_start, right_start_value)
        right_end = c(right_end, right_end_value)
        base_width = c(base_width , left_end_value - left_start_value)
        color = c(color,label)
      }
    }
    else
    {
      if(chrom_end_value == as.numeric(sub("chr","",seqlevels(current_range)))) {
        if((right_start_value >= start(current_range) && right_start_value <= end(current_range)) || (right_end_value >= start(current_range) && right_end_value <= end(current_range)))
        {
          chrom_start = c(chrom_start, chrom_start_value)
          chrom_end = c(chrom_end, chrom_end_value)
          left_start = c(left_start, left_start_value)
          left_end = c(left_end, left_end_value)
          right_start = c(right_start, right_start_value)
          right_end = c(right_end, right_end_value)
          base_width = c(base_width , left_end_value - left_start_value)
          color = c(color,label)
        }
      }
    }
    
  }
  
  out = list()
  left_ranges = IRanges(left_start,left_end)
  right_ranges = IRanges(right_start,right_end)
  
  if(length(left_ranges) != 0) {
    
    left_range <- GRanges(seqnames = chrom_start, ranges = left_ranges, base_width = base_width, color = color)
    right_range <- GRanges(seqnames = chrom_end, ranges = right_ranges, base_width = base_width, color = color)
    print(paste0("-> Find ",length(left_ranges)," interaction for ", label))
    out = list(left_range, right_range)
  }
  else
  {
    print(paste0("No relation found for ", label))
  }
  
  invisible(out)
}

get_chiapet_arch <- function(rep) {
  
  g = ggplot(rep) +
    geom_arch() +
    theme_bw() +
    #aes(size=base_width, color=color) +
    #aes(size = alpha) +
    aes(color=color, alpha = alpha) +
#     aes(color=color) +
#     scale_colour_manual(values = c("gray","promoter"="black")) +
    xlim(current_range) +
    theme(axis.ticks = element_blank(), axis.text.y = element_blank()) +
    ylab("") + guides(alpha=FALSE, color=FALSE)
  
#   if (!is.null(highlight_region)) {
#     
#     highlight_region = strsplit(x = highlight_region, fixed = T, split = ":")[[1]]
#     d = data.frame(x1=as.numeric(highlight_region[1]), x2=as.numeric(highlight_region[2]), y1=0, y2=10)
#     g = g + geom_rect(d, mapping=aes(xmin=x1, xmax=x2,ymin=y1, ymax=y2), color="black", alpha = 0.4) + xlim(current_range)
#   }
#   
  g
  
}


make_emphasis = function(range, range_degree) {
  
  for (i in(1 : length(range))) {
    
    range_row = range[i]
    start_value = start(range_row)
    end_value = end(range_row)
    # in the spec zone?
    if((start_value > (start(range_degree) - 1) 
        && start_value < (end(range_degree) - 1)) 
       || (end_value > (start(range_degree) - 1) 
           && end_value < (end(range_degree) - 1))) {
      range[i]$alpha = range_degree$alpha
      range[i]$color = range_degree$color
    }
  }
  
  range
}

