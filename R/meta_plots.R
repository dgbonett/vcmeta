# replicate.plot
#' Plot to compare estimates from original and follow-up studies
#'
#'
#' @description 
#' Generates a basic plot using ggplot2 to visualize the estimates from
#' and original and follow-up studies.
#'  
#'  
#' @param result - a result matrix from any of the replicate functions in vcmeta
#' @param focus - Optional specification of the focus of the plot; 
#'   defaults to 'Both'
#' * Both - plots each estimate, differencence, and average
#' * Difference - plot each estimate and difference between them
#' * Average - plot each estimate and the average effect size
#' @param reference_line - Optional x-value for a reference line. Only applies
#'   if focus is 'Difference' or 'Both'. Defaults to NULL, in which case a 
#'   reference line is not drawn.
#' @param diamond_height - Optional height of the diamond representing average
#'  effect size. Only applies if focus is 'Average' or 'Both'.
#'  Defaults to 0.2
#' @param difference_axis_ticks - Optional requested number of ticks on the
#'   difference axis. Only applies if focus is 'Difference' or 'Both'. 
#'   Defaults to 5.
#' @param ggtheme - optional ggplot2 theme object; defaults to theme_classic()
#' 
#' 
#' @return 
#' Returns a ggplot object.  If stored, can be further customized via
#'   the ggplot API
#'  
#'    
#' @examples    
#' # Compare Damisch et al., 2010 to Calin-Jageman & Caldwell 2014
#' # Damisch et al., 2010, Exp 1, German participants made 10 mini-golf putts.
#' # Half were told they had a 'lucky' golf ball; half were not.
#' # Found a large but uncertain improvement in shots made in the luck condition
#' # Calin-Jageman & Caldwell, 2014, Exp 1, was a pre-registered replication with
#' # input from Damisch, though with English-speaking participants.
#' #
#' # Here we compare the effect sizes, in original units, for the two studies.
#' # Use the replicate.mean2 function because the design is a 2-group design.
#' 
#' library(ggplot2)
#' damisch_v_calinjageman_raw <- replicate.mean2(
#'   alpha = 0.05,
#'   m11 = 6.42,
#'   m12 = 4.75,
#'   sd11 = 1.88,
#'   sd12 = 2.15,
#'   n11 = 14,
#'   n12 = 14,
#'   m21 = 4.73,
#'   m22 = 4.62,
#'   sd21 = 1.958,
#'   sd22 = 2.12,
#'   n21 = 66,
#'   n22 = 58
#' )
#' 
#' # View the comparison:
#' damisch_v_calinjageman_raw
#' 
#' 
#' # Now plot the comparison, focusing on the difference
#' replicate.plot(damisch_v_calinjageman_raw, focus = "Difference")
#' 
#' # Plot the comparison, focusing on the average
#' replicate.plot(damisch_v_calinjageman_raw, 
#'   focus = "Average", 
#'   reference_line = 0,
#'   diamond_height = 0.1
#' )
#' 
#' 
#' # Plot the comparison with both difference and average.
#' # In this case, store the plot for manipulation
#' myplot <- replicate.plot(
#'   damisch_v_calinjageman_raw,
#'   focus = "Both",
#'   reference_line = 0
#' )
#' 
#' # View the stored plot
#' myplot
#' 
#' # Change x-labels and study labels
#' myplot <- myplot + xlab("Difference in Putts Made, Lucky - Control")
#' myplot <- myplot + scale_y_discrete(
#'     labels = c(
#'       "Average",
#'       "Difference",
#'       "Calin-Jageman & Caldwell, 2014",
#'       "Damisch et al., 2010"
#'       )
#'   )
#'   
#' # View the updated plot
#' myplot
#' 
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 aes_string
#' @importFrom ggplot2 geom_linerange
#' @importFrom ggplot2 geom_segment
#' @importFrom ggplot2 geom_vline
#' @importFrom ggplot2 geom_point
#' @importFrom ggplot2 geom_polygon
#' @importFrom ggplot2 sec_axis
#' @importFrom ggplot2 theme
#' @importFrom ggplot2 ylab
#' @importFrom ggplot2 element_blank
#' @importFrom ggplot2 scale_x_continuous
#' @export
replicate.plot <- function(
  result, 
  focus = c("Both", "Difference", "Average"),
  reference_line = NULL,
  diamond_height = 0.2,
  difference_axis_ticks = 5,
  ggtheme = ggplot2::theme_classic()
) {
  
  # Options ---------------------------------------------
  focus <- match.arg(focus)
  plot_average <- focus != "Difference"
  plot_difference <- focus != "Average"
  is_log <- "exp(Estimate)" %in% colnames(result)
  diff_axis_y <- 0
  
  # Definitions ------------------------------------------
  # Row names
  comparison_name <- "Original:"
  ref_name <- "Follow-up:"
  avg_name <- "Average:"
  diff_name <- "Original - Follow-up:"
  
  # Column names
  se_name <- "SE"
  if (is_log) {
    es_name <- "exp(Estimate)"
    ll_name <- "exp(LL)"
    ul_name <- "exp(UL)"
  } else {
    es_name <- "Estimate"
    ll_name <- "LL"
    ul_name <- "UL"
  }
  col_adj <- c(es_name, ll_name, ul_name)
  
  
  # Data prep --------------------------------------------
  # Convert matrix to data frame, drop unused columns, save names to name
  as_df <- as.data.frame(result)
  as_df <- as_df[ , c(col_adj, se_name)]
  as_df$name <- row.names(result)
  
  # Filter out average or difference based on focus
  if(!plot_average) {
    as_df <- as_df[as_df$name != avg_name, ]
  }
  if(!plot_difference) {
    as_df <- as_df[as_df$name != diff_name, ]
  }
  
  # Convert names to factor with levels in rev order for plotting
  as_df$name <- factor(
    as_df$name,
    levels = rev(as_df$name)
  )
  
  # Set y values and size
  as_df$y_value <- as.integer(as_df$name)
  as_df$size <- max(as_df[[se_name]]) - as_df[[se_name]] + min(as_df[[se_name]])
  as_df$diff_axis_y <- diff_axis_y
  
  if (plot_difference) {
    # Shift difference by reference value
    ref_es <- as_df[ref_name, es_name]
    as_df[diff_name, col_adj] <- as_df[diff_name, col_adj] + ref_es
    
    # Calculate floating axis breaks
    diff_start <- min(0, as_df[[ll_name]]-ref_es)
    diff_end <-  max(0, as_df[[ul_name]]-ref_es)
    diff_breaks <- pretty(diff_start:diff_end, n = difference_axis_ticks)
  }
   
  if (plot_average) {
    # Generate polygon data
    row_avg <- as.list(as_df[avg_name, ])
    diamond_xs <- c(
      row_avg[[ll_name]], 
      row_avg[[es_name]], 
      row_avg[[ul_name]], 
      row_avg[[es_name]]
    )
    d_y <- row_avg$y_value
    diamond_ys <- c(d_y, d_y - diamond_height , d_y, d_y + diamond_height)
    poly_data <- data.frame(x = diamond_xs, y = diamond_ys)
  }
  
  
  # Make graph ----------------------------------------------
  # Basic graph
  myplot <- ggplot(
    data = as_df, 
    aes_string(x = es_name, y = "name")
  )
  myplot <- myplot + ggtheme

  # CIs
  myplot <- myplot + geom_linerange(aes_string(xmin = ll_name, xmax = ul_name))
  
  # For differences, plot reference lines
  if (plot_difference) {
    # Comparison line
    myplot <- myplot + geom_segment(
      data = as_df[diff_name, ], 
      linetype = "dotted",
      color = "black",
      aes_string(
        x = es_name, 
        xend = es_name, 
        y = "y_value", 
        yend = "diff_axis_y"
      )
    )
    
    # Reference line
    myplot <- myplot + geom_segment(
      data = as_df[ref_name, ], 
      linetype = "dashed",
      color = "black",
      aes_string(
        x = es_name, 
        xend = es_name, 
        y = "y_value", 
        yend = "diff_axis_y"
      )
    )

  } 
  
  if (plot_average & !is.null(reference_line)) {
    myplot <- myplot + geom_vline(
      xintercept = reference_line, 
      linetype = "dotted"
    )
  }

  # Effect sizes
  myplot <- myplot + geom_point(
    aes_string(
      colour = "name", 
      fill = "name",
      size = "size"
    ),
    shape = "square filled"
  )
  
  # For averages, plot diamond for effect size
  if (plot_average) {
    myplot <- myplot + geom_polygon(
      data = poly_data, 
      aes_string(x = "x", y = "y")
    )    
  }

  # For differences, display floating difference axis
  if (plot_difference) {
    # Specify axis
    myplot <- myplot + scale_x_continuous(
      position = "top",
      sec.axis = sec_axis(
        name = "Difference",
        trans = ~.-ref_es,
        breaks = diff_breaks
      )
    )
    
    # Floating difference axis
    myplot <- myplot + geom_segment(
      linetype = "solid", 
      color = "black", 
      aes(
        x = min(diff_breaks)+ref_es, 
        xend = max(diff_breaks)+ref_es, 
        y = diff_axis_y, 
        yend = diff_axis_y
      )
    )
  }
  
  # Clean up axis lines and hide legends
  myplot <- myplot + theme(axis.line.y.left = element_blank())
  myplot <- myplot + theme(axis.ticks.y.left = element_blank())
  if(plot_difference) {
    myplot <- myplot + theme(axis.line.x.bottom = element_blank())
  }
  myplot <- myplot + ylab("")
  myplot <- myplot + theme(legend.position = "none")
 
  return(myplot) 
}


# meta.ave.plot 
#' Forest plot for average effect sizes
#'
#'
#' @description 
#' Generates a forest plot to visualize effect sizes estimates and overall
#' averages from the meta.ave functions in vcmeta. If the column
#' exp(Estimate) is present, this function plots the exponentiated
#' effect size and CI found in columns exp(Estimate), exp(LL), and exp(UL).
#' Otherwise, this function plots the effect size and CI found in
#' the columns Estimate, LL, and UL.
#'  
#'  
#' @param result - a result matrix from any of the replicate functions in vcmeta
#' @param reference_line - Optional x-value for a reference line.  Only applies
#'   if focus is 'Difference' or 'Both'.  Defaults to NULL, in which case a 
#'   reference line is not drawn.
#' @param diamond_height - Optional height of the diamond representing average
#'  effect size. Only applies if focus is 'Average' or 'Both'.
#'  Defaults to 0.2
#' @param ggtheme - optional ggplot2 theme object; defaults to theme_classic()
#' 
#' 
#' @return 
#' Returns a ggplot object.  If stored, can be further customized via
#'   the ggplot API
#'  
#' @examples
#' # Plot results from meta.ave.mean2
#' m1 <- c(7.4, 6.9)
#' m2 <- c(6.3, 5.7)
#' sd1 <- c(1.72, 1.53)
#' sd2 <- c(2.35, 2.04)
#' n1 <- c(40, 60)
#' n2 <- c(40, 60)
#' result <- meta.ave.mean2(.05, m1, m2, sd1, sd2, n1, n2, bystudy = TRUE)
#' meta.ave.plot(result, reference_line = 0)
#' 
#' 
#' # Plot results from meta.ave.meanratio2
#' # Note that this plots the exponentiated effect size and CI
#' m1 <- c(53, 60, 53, 57)
#' m2 <- c(55, 62, 58, 61)
#' sd1 <- c(4.1, 4.2, 4.5, 4.0)
#' sd2 <- c(4.2, 4.7, 4.9, 4.8)
#' cor <- c(.7, .7, .8, .85)
#' n <- c(30, 50, 30, 70)
#' result <- meta.ave.meanratio.ps(.05, m1, m2, sd1, sd2, cor, n, bystudy = TRUE)
#' myplot <- meta.ave.plot(result, reference_line = 1)
#' myplot
#' 
#' # Change x-scale to log2
#' library(ggplot2)
#' myplot <- myplot + scale_x_continuous(
#'   trans = 'log2', 
#'   limits = c(0.75, 1.25), 
#'   name = "Estimated Ratio of Means, Log2 Scale"
#' )
#' myplot
#'
#'
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 aes_string
#' @importFrom ggplot2 geom_linerange
#' @importFrom ggplot2 geom_segment
#' @importFrom ggplot2 geom_vline
#' @importFrom ggplot2 geom_point
#' @importFrom ggplot2 geom_polygon
#' @importFrom ggplot2 sec_axis
#' @importFrom ggplot2 theme
#' @importFrom ggplot2 ylab
#' @importFrom ggplot2 xlab
#' @importFrom ggplot2 element_blank
#' @importFrom ggplot2 scale_x_continuous
#' @importFrom ggplot2 scale_y_continuous
#' @importFrom utils head
#' @export    
meta.ave.plot <- function(
  result, 
  reference_line = NULL,
  diamond_height = 0.2,
  ggtheme = ggplot2::theme_classic()
) {

  
  # Options ----------------------------------------------
  is_log <- "exp(Estimate)" %in% colnames(result)
  

  # Definitions ------------------------------------------
  avg_name <- "Average"
  se_name <- "SE"
  if (is_log) {
    es_name <- "exp(Estimate)"
    ll_name <- "exp(LL)"
    ul_name <- "exp(UL)"
  } else {
    es_name <- "Estimate"
    ll_name <- "LL"
    ul_name <- "UL"
  }
  
  
  # Data prep --------------------------------------------
  # Convert matrix to data frame
  as_df <- as.data.frame(result)
  
  # Move average to bottom
  as_df <- rbind(
    as_df[-1, ],
    head(as_df, 1)
  )
  
  # Set name column and levels for proper order
  as_df$name <- factor(
    row.names(as_df), 
    levels = rev(row.names(as_df))
  )

  # Set y values and size
  as_df$y_value <- as.integer(as_df$name)
  as_df$size <- max(as_df[ , se_name]) - as_df[, se_name] + min(as_df[ , se_name])
  
  # Generate polygon data
  # Generate polygon data
  row_avg <- as.list(as_df[avg_name, ])
  diamond_xs <- c(
    row_avg[[ll_name]], 
    row_avg[[es_name]], 
    row_avg[[ul_name]], 
    row_avg[[es_name]]
  )
  d_y <- row_avg$y_value
  diamond_ys <- c(d_y, d_y - diamond_height , d_y, d_y + diamond_height)
  poly_data <- data.frame(x = diamond_xs, y = diamond_ys)
  
  
  # Make graph ----------------------------------------------
  # Basic graph
  myplot <- ggplot(
    data = as_df, 
    aes_string(x = es_name, y = "name")
  )
  myplot <- myplot + ggtheme

  # Optional reference line  
  if (!is.null(reference_line)) {
    myplot <- myplot + geom_vline(
      xintercept = reference_line, 
      linetype = "dotted"
    )
    
  }
    
  # CIs
  myplot <- myplot + geom_linerange(aes_string(xmin = ll_name, xmax = ul_name))
  
  # Effect sizes
  myplot <- myplot + geom_point(
    aes_string(
      colour = "name", 
      fill = "name",
      size = "size"
    ),
    shape = "square filled"
  )
  
  # Diamond for average effect size
  myplot <- myplot + geom_polygon(
    data = poly_data, 
    aes_string(x = "x", y = "y")
  )   
  
  # Clean up axis lines and hide legends
  myplot <- myplot + ylab("")
  myplot <- myplot + theme(legend.position = "none")
  
  return(myplot) 
}
