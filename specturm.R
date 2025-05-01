args=commandArgs(trailingOnly=TRUE)

library(MutationalPatterns)
library(ggplot2)
library(BSgenome)
library(reshape2)
library(dplyr)

plot_indel_profile <- function(counts, colors = NA) {
   colors <- c(
    "#FDBE6F", "#FF8001", "#B0DD8B", "#36A12E", "#FDCAB5", "#FC8A6A",
    "#F14432", "#BC141A", "#D0E1F2", "#94C4DF", "#4A98C9", "#1764AB",
    "#E2E2EF", "#B6B6D8", "#8683BD", "#61409B"
  )

  count <- muttype <- muttype_sub <- muttype_total <- sample <- NULL

  # Separate muttype and muttype_sub. Then make data long
  counts <- counts %>%
    as.data.frame() %>%
    tibble::rownames_to_column("muttype_total") %>%
    tidyr::separate(muttype_total, c("muttype", "muttype_sub"), sep = "_(?=[0-9])") %>%
    dplyr::mutate(muttype = factor(muttype, levels = unique(muttype))) %>%
    tidyr::gather(key = "sample", value = "count", -muttype, -muttype_sub) %>% 
    dplyr::mutate(sample = factor(sample, levels = unique(sample)))

  # Count nr mutations. (This is used for the facets)
  nr_muts <- counts %>%
    dplyr::group_by(sample) %>%
    dplyr::summarise(nr_muts = round(sum(count)))

  # Create facet texts
  facet_labs_y <- stringr::str_c(nr_muts$sample)
  names(facet_labs_y) <- nr_muts$sample
  facet_labs_x <- c("1: C", "1: T", "1: C", "1: T", 2, 3, 4, "5+", 2, 3, 4, "5+", 2, 3, 4, "5+")
  names(facet_labs_x) <- levels(counts$muttype)

  width <- 0.8
  spacing <- 0.8

  # Create figure
  fig <- ggplot(counts, aes(x = muttype_sub, y = count, fill = muttype, width = width)) + 
    geom_bar(stat = "identity") +
    facet_grid(sample ~ muttype,
      scales = 'free_x', space = "free_x",
      labeller = labeller(muttype = facet_labs_x, sample = facet_labs_y)
    ) +
    scale_fill_manual(values = colors) +
    xlab('') +
    ylab("Percentage") +
    guides(fill = FALSE) +
    theme_minimal() +
    theme(
      axis.title.y = element_text(size = 10, vjust = 1),
      axis.text.y = element_text(size = 10),
      axis.title.x = element_text(size = 10),
      axis.text.x = element_text(size = 6),
      strip.text.x = element_text(size = 10),
      strip.text.y = element_text(size = 10),
      strip.background = element_blank(),
      panel.spacing.x = unit(spacing, "lines") 
   )

  return(fig)
}

plot_snv_profile <- function(mut_matrix, colors = NA) {
  library(dplyr)
  colors <- c(
  "#2EBAED", "#000000", "#DE1C14",
  "#D4D2D2", "#ADCC54", "#F0D0CE"
  )

  # Get substitution and context from rownames and make long.
  tb <- mut_matrix %>%
    as.data.frame() %>%
    tibble::rownames_to_column("full_context") %>%
    dplyr::mutate(
      substitution = stringr::str_replace(full_context, "\\w\\[(.*)\\]\\w", "\\1"),
      context = stringr::str_replace(full_context, "\\[.*\\]", "\\.")
    ) %>%
    dplyr::select(-full_context) %>%
    tidyr::pivot_longer(c(-substitution, -context), names_to = "sample", values_to = "count") %>%
    dplyr::mutate(sample = factor(sample, levels = unique(sample)))

    width <- 1
    spacing <- 0

  # Create figure
  plot <- ggplot(data = tb, aes(
    x = context,
    y = count,
    fill = substitution,
    width = width
  )) +
    geom_bar(stat = "identity", colour = "black", size = .2) +
    scale_fill_manual(values = colors) +
    facet_grid(sample ~ substitution) +
    ylab("Percentage") +
    ylim(0,0.1) +
#    ylab("Absolute burden") +
    guides(fill = FALSE) +
    theme_bw() +
    theme(
      axis.title.y = element_text(size = 12, vjust = 1),
      axis.text.y = element_text(size = 8),
      axis.title.x = element_text(size = 12),
      axis.text.x = element_text(size = 5, angle = 90, vjust = 0.5),
      strip.text.x = element_text(size = 9),
      strip.text.y = element_text(size = 9),
      strip.background = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.spacing.x = unit(spacing, "lines")
    )

  return(plot)
}


