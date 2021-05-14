
library(readr)
disease <- "Urate_EA"
annota <- read_tsv("/Users/jnz/Document/JHU/Research/PWAS/Analysis/4_PWAS_results/prot.anno_autosomal.txt")
dat.pwas <- read_tsv(paste0("/Users/jnz/Document/JHU/Research/PWAS/Analysis/4_PWAS_results/", disease, "/allchr.pwas"))
p.pwas <- 0.05/nrow(dat.pwas)
dat.pwas <- dat.pwas[!(is.na(dat.pwas$TWAS.P)),]
pos_lookup <- read_tsv("/Users/jnz/Document/JHU/Research/PWAS/Analysis/4_PWAS_results/Plasma_Protein_hg19.pos")
pos_lookup <- unique(pos_lookup[,c(3,5,6)])
m <- match(dat.pwas$ID,pos_lookup$ID)
dat.pwas$P0 <- pos_lookup$P0[m]
dat.pwas$P1 <- pos_lookup$P1[m]
dat_Urate <- dat.pwas

disease <- "Gout_EA"
annota <- read_tsv("/Users/jnz/Document/JHU/Research/PWAS/Analysis/4_PWAS_results/prot.anno_autosomal.txt")
dat.pwas <- read_tsv(paste0("/Users/jnz/Document/JHU/Research/PWAS/Analysis/4_PWAS_results/", disease, "/allchr.pwas"))
p.pwas <- 0.05/nrow(dat.pwas)
dat.pwas <- dat.pwas[!(is.na(dat.pwas$TWAS.P)),]
pos_lookup <- read_tsv("/Users/jnz/Document/JHU/Research/PWAS/Analysis/4_PWAS_results/Plasma_Protein_hg19.pos")
pos_lookup <- unique(pos_lookup[,c(3,5,6)])
m <- match(dat.pwas$ID,pos_lookup$ID)
dat.pwas$P0 <- pos_lookup$P0[m]
dat.pwas$P1 <- pos_lookup$P1[m]
dat_Gout <- dat.pwas

dat_Urate$TWAS.P[dat_Urate$ID=="INHBC"] = 10^(-30)
dat_Urate$ID[dat_Urate$ID=="INHBC"] = "INHBC(7.95e-63)"

dat <- rbind(cbind(dat_Gout, Gout=rep(-1, nrow(dat_Gout))), cbind(dat_Urate, Gout=rep(1, nrow(dat_Urate))))


ggmiami <- function(
  data,
  split_by,
  split_at,
  chr = "chr",
  pos = "pos",
  p  =  "p",
  chr_colors = c("black", "grey"),
  upper_ylab = "-log10(p)",
  lower_ylab = "-log10(p)",
  genome_line = 5e-8,
  genome_line_color = "red",
  suggestive_line = 1e-5,
  suggestive_line_color = "blue",
  hits_label_col = NULL,
  hits_label = NULL,
  top_n_hits1 = 5,
  top_n_hits2 = 5,
  upper_labels_df = NULL,
  lower_labels_df = NULL,
  upper_highlight = NULL,
  upper_highlight_col = NULL,
  upper_highlight_color = "green",
  lower_highlight = NULL,
  lower_highlight_col = NULL,
  lower_highlight_color = "green") {
  
  # Prepare the data
  plot_data <- prep_miami_data(data = data, split_by = split_by,
                               split_at = split_at, chr = chr, pos = pos, p = p)
  
  # Prepare the colors
  if (length(chr_colors) == 2) {
    chr_colors <- rep(chr_colors, length.out = nrow(plot_data$axis))
  } else if (length(chr_colors) == 1) {
    chr_colors <- rep(chr_colors, length.out = nrow(plot_data$axis))
  } else if (length(chr_colors) == nrow(plot_data$axis)) {
    chr_colors <- chr_colors
  } else {
    stop("The number of colors specified in {chr_colors} does not match the
         number of chromosomes to be displayed.")
  }
  
  # Prepare the axis titles
  if (upper_ylab == "-log10(p)") {
    upper_ylab <- expression("-log" [10]* "(p)")
  } else {
    upper_ylab <- bquote(atop(.(upper_ylab), "-log" [10]* "(p)"))
    # upper_ylab <- bquote(atop(NA, atop(textstyle(.(upper_ylab)),
    # textstyle("-log" [10]* "(p)"))))
    
  }
  
  # Prepare the axis titles
  if (lower_ylab == "-log10(p)") {
    lower_ylab <- expression("-log" [10]* "(p)")
  } else {
    lower_ylab <- bquote(atop(.(lower_ylab), "-log" [10]* "(p)"))
    # lower_ylab <- bquote(atop(atop(textstyle(.(lower_ylab)),
    # textstyle("-log" [10]* "(p)")), NA))
  }
  
  # When ggtext is published on CRAN, this is a better solution for axis text.
  # lower_ylab <- paste0(lower_ylab, "<br>-log<sub>10</sub>(p)")
  # Then in ggplot2 call: axis.title.y = ggtext::element_markdown() and un-set
  # the l margin.
  
  # Create base upper plot.
  upper_plot <- ggplot2::ggplot() +
    ggplot2::geom_point(data = plot_data$upper,
                        aes(x = .data$rel_pos, y = .data$logged_p,
                            color = as.factor(.data$chr)), size = 0.25) +
    ggplot2::scale_color_manual(values = chr_colors) +
    ggplot2::scale_x_continuous(labels = plot_data$axis$chr,
                                breaks = plot_data$axis$chr_center,
                                expand = ggplot2::expansion(mult = 0.01),
                                guide = ggplot2::guide_axis(check.overlap =
                                                              TRUE)) +
    ggplot2::scale_y_continuous(limits = c(0, plot_data$maxp),
                                expand =
                                  ggplot2::expansion(mult = c(0.02, 0))) +
    ggplot2::labs(x = "", y = upper_ylab) +
    ggplot2::theme_classic() +
    ggplot2::theme(legend.position = "none",
                   axis.title.x = ggplot2::element_blank(),
                   plot.margin = ggplot2::margin(t = 5, r = 3, b = 5, l = 10))
  
  # Create base lower plot
  lower_plot <- ggplot2::ggplot() +
    ggplot2::geom_point(data = plot_data$lower,
                        aes(x = .data$rel_pos, y = .data$logged_p,
                            color = as.factor(.data$chr)), size = 0.25) +
    ggplot2::scale_color_manual(values = chr_colors) +
    ggplot2::scale_x_continuous(breaks = plot_data$axis$chr_center,
                                position = "top",
                                expand = ggplot2::expansion(mult = 0.01)) +
    ggplot2::scale_y_reverse(limits = c(plot_data$maxp, 0),
                             expand = ggplot2::expansion(mult = c(0, 0.02))) +
    ggplot2::labs(x = "", y = lower_ylab) +
    ggplot2::theme_classic() +
    ggplot2::theme(legend.position = "none",
                   axis.text.x = ggplot2::element_blank(),
                   axis.title.x = ggplot2::element_blank(),
                   plot.margin = ggplot2::margin(t = 5, r = 3, b = 5, l = 10))
  
  # If the user has requested a suggetive line, add:
  if (!is.null(suggestive_line)) {
    upper_plot <- upper_plot +
      ggplot2::geom_hline(yintercept = -log10(suggestive_line),
                          color = suggestive_line_color,
                          linetype = "solid", size = 0.4)
    lower_plot <- lower_plot +
      ggplot2::geom_hline(yintercept = -log10(suggestive_line),
                          color = suggestive_line_color,
                          linetype = "solid", size = 0.4)
  }
  
  # If the user has requested a genome-wide line, add:
  if (!is.null(genome_line)) {
    upper_plot <- upper_plot +
      ggplot2::geom_hline(yintercept = -log10(genome_line),
                          color = genome_line_color,
                          linetype = "dashed", size = 0.4)
    lower_plot <- lower_plot +
      ggplot2::geom_hline(yintercept = -log10(genome_line),
                          color = genome_line_color,
                          linetype = "dashed", size = 0.4)
  }
  
  # If they try to specify both hits_label_col and dataframe(s), return an
  # error message.
  if (all(!is.null(hits_label_col), any(!is.null(upper_labels_df),
                                        !is.null(lower_labels_df)))) {
    stop("You have specified both hits_label_col and a *_labels_df. This
         package does not know how to use both information simultaneously.
         Please only use one method for labelling: either hits_label_col (with
         or without hits_label), or *_labels_df.")
  }
  
  # If the user requests labels based on the hits_label_col
  if (all(!is.null(hits_label_col), is.null(upper_labels_df),
          is.null(lower_labels_df))) {
    # Create the labels for the upper and lower plot
    upper_labels_df <- make_miami_labels(data = plot_data$upper,
                                         hits_label_col = hits_label_col,
                                         hits_label = hits_label,
                                         top_n_hits = top_n_hits1)
    
    lower_labels_df <- make_miami_labels(data = plot_data$lower,
                                         hits_label_col = hits_label_col,
                                         hits_label = hits_label,
                                         top_n_hits = top_n_hits2)
    
    # Add to the plots
    upper_plot <- upper_plot +
      ggrepel::geom_label_repel(data = upper_labels_df,
                                aes(x = .data$rel_pos,
                                    y = .data$logged_p,
                                    label = .data$label),
                                size = 2, segment.size = 0.2,
                                point.padding = 0.3,
                                ylim = c(plot_data$maxp / 2, NA),
                                min.segment.length = 0, force = 2,
                                box.padding = 0.5)
    
    lower_plot <- lower_plot +
      ggrepel::geom_label_repel(data = lower_labels_df,
                                aes(x = .data$rel_pos,
                                    y = .data$logged_p,
                                    label = .data$label),
                                size = 2, segment.size = 0.2,
                                point.padding = 0.3,
                                ylim = c(NA, -(plot_data$maxp / 2)),
                                min.segment.length = 0, force = 2,
                                box.padding = 0.5)
  }
  
  # If the user requests labels based on a dataframe for the upper plot
  if (all(is.null(hits_label_col), !is.null(upper_labels_df))) {
    # Make sure that the column names in upper_labels_df are correct
    checkmate::assertNames(colnames(upper_labels_df),
                           identical.to = c("rel_pos", "logged_p", "label"))
    
    # Add to plot
    upper_plot <- upper_plot +
      ggrepel::geom_label_repel(data = upper_labels_df,
                                aes(x = .data$rel_pos,
                                    y = .data$logged_p,
                                    label = .data$label),
                                size = 2, segment.size = 0.2,
                                point.padding = 0.3,
                                ylim = c(plot_data$maxp / 2, NA),
                                min.segment.length = 0, force = 2,
                                box.padding = 0.5)
  }
  
  # If the user requests labels based on a dataframe for the lower plot
  if (all(is.null(hits_label_col), !is.null(lower_labels_df))) {
    # Make sure that the column names in lower_labels_df are correct
    checkmate::assertNames(colnames(lower_labels_df),
                           identical.to = c("rel_pos", "logged_p", "label"))
    
    # Add to plot
    lower_plot <- lower_plot +
      ggrepel::geom_label_repel(data = lower_labels_df,
                                aes(x = .data$rel_pos,
                                    y = .data$logged_p,
                                    label = .data$label),
                                size = 2, segment.size = 0.2,
                                point.padding = 0.3,
                                ylim = c(NA, -(plot_data$maxp / 2)),
                                min.segment.length = 0, force = 2,
                                box.padding = 0.5)
  }
  
  # Highlight upper snps
  if (all(!is.null(upper_highlight), !is.null(upper_highlight_col))) {
    # Make highlighting df
    upper_highlight_df <- highlight_miami(data = plot_data$upper,
                                          highlight = upper_highlight,
                                          highlight_col = upper_highlight_col,
                                          highlight_color =
                                            upper_highlight_color)
    
    # Add to plot
    upper_plot <- upper_plot +
      ggplot2::geom_point(data = upper_highlight_df,
                          aes(x = .data$rel_pos, y = .data$logged_p),
                          color = upper_highlight_df$color,
                          size = 0.25)
  }
  
  # Highlight lower snps
  if (all(!is.null(lower_highlight), !is.null(lower_highlight_col))) {
    # Make highlighting df
    lower_highlight_df <- highlight_miami(data = plot_data$lower,
                                          highlight = lower_highlight,
                                          highlight_col = lower_highlight_col,
                                          highlight_color =
                                            lower_highlight_color)
    
    # Add to plot
    lower_plot <- lower_plot +
      ggplot2::geom_point(data = lower_highlight_df,
                          aes(x = .data$rel_pos, y = .data$logged_p),
                          color = lower_highlight_df$color,
                          size = 0.25)
  }
  
  # Put the two together
  gridExtra::grid.arrange(upper_plot, lower_plot, nrow = 2)
  
}

# colorRampPalette(brewer.pal(12, "Paired"))(22)
mycolors <- c("#A6CEE3", "#5FA0CA", "#257CB2", "#72B29C", "#A5D981", "#63B84F", 
              "#4F9F3B", "#B89B74", "#F68181", "#E93E3F", "#E9412F", "#F6975B", 
              "#FDAC4F", "#FE8B15", "#ED8F47", "#D1AAB7", "#A585BF", "#73489F", 
              "#A99099", "#F7F599", "#D9AF63", "#B15928")
p <- ggmiami(data = dat, 
             chr_colors = mycolors, 
             split_by = "Gout", split_at = 0, p = "TWAS.P",
             chr = "CHR", top_n_hits1=8, top_n_hits2=4,
             pos = "P0", suggestive_line=NULL, genome_line = p.pwas, genome_line_color = "black", 
             hits_label_col="ID", hits_label=c("INHBB","ITIH1","SPP1","BTN3A3","C11orf68","B3GAT3","INHBC","INHBC(7.95e-63)","SNUPN","IL1RN"),
             upper_ylab = "Urate",
             lower_ylab = "Gout")
ggsave(filename="PWAS.png", 
       plot=p, device="png",
       path="/Users/jnz/Document/JHU/Research/PWAS/Analysis/*Figures/", 
       width=8, height=4, units="in", dpi=500)

