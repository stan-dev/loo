#' Compare models across domains
#' 
#' The LOO difference plot shows how the ELPD of two different models
#' changes when a predictor is varied. This can is useful for identifying
#' opportunities for model stacking or expansion.
#' 
#' @param y A vector of observations. See Details.
#' @param psis_object_1,psis_object_2 If using loo version 2.0.0 or greater, 
#' an object returned by the `[loo::psis()]` function (or by the 
#' `[loo::loo()]` function with argument `save_psis` set to `TRUE`).
#' @param ... Currently unused.
#' @param group A grouping variable (a vector or factor) the same length 
#' as `y`. Each value in group is interpreted as the group level pertaining 
#' to the corresponding value of `y`. If `FALSE`, ignored.
#' @param outlier_thresh Flag values when the difference in the ELPD exceeds
#' this threshold. Defaults to `NULL`, in which case no values are flagged.
#' @param size,alpha,jitter Passed to `[ggplot2::geom_point()]` to control 
#' aesthetics. `size` and `alpha` are passed to to the `size` and `alpha` 
#' arguments of `[ggplot2::geom_jitter()]` to control the appearance of 
#' points. `jitter` can be either a number or a vector of numbers.
#' Passing a single number will jitter variables along the x axis only, while 
#' passing a vector will jitter along both axes.
#' @param sort_by_group Sort observations by `group`, then plot against an 
#' arbitrary index. Plotting by index can be useful when categories have 
#' very different sample sizes. 
#' 
#' 
#' @template return-ggplot
#' 
#' @template reference-vis-paper
#' 
#' @examples 
#' 
#' library(loo)
#' 
#' cbPalette <- c("#636363", "#E69F00", "#56B4E9", "#009E73", 
#'                "#F0E442", "#0072B2","#CC79A7")
#' 
#' # Plot using groups from WHO
#' 
#' plot_loo_dif(factor(GM@data$super_region_name), loo3, loo2, 
#'              group = GM@data$super_region_name, alpha = .5, 
#'              jitter = c(.45, .2)
#'              ) + 
#'              xlab("Region") + scale_colour_manual(values=cbPalette)
#' 
#' # Plot using groups identified with clustering
#' 
#' plot_loo_dif(factor(GM@data$cluster_region), loo3, loo2, 
#'              group = GM@data$super_region_name, alpha = .5, 
#'              jitter = c(.45, .2)
#'              ) + 
#'              xlab("Cluster Group") + scale_colour_manual(values=cbPalette)
#'              
#' # Plot using an index variable to reduce crowding
#' 
#' plot_loo_dif(1:2980, loo3, loo2, group = GM@data$super_region_name, 
#'              alpha = .5, sort_by_group = TRUE, 
#'              ) + 
#'              xlab("Index") + scale_colour_manual(values=cbPalette)
#'
#' 
#' # Example using kid IQ Dataset with a continuous predictor
#'
#' data(kidiq)
#'
#' t_prior <- student_t(df = 10, location = 0, scale = .5)
#' coef_prior <- student_t(df = 10, location = .5, scale = .25)
#' kidiq$kid_std <- (kidiq$kid_score - 100) / 15
#' kidiq$mom_std <- (kidiq$mom_iq - 100) / 15
#' kidiq$age_std <- (kidiq$mom_age - mean(kidiq$mom_age)) / sd(kidiq$mom_age)
#' kidiq$hs_cent <- kidiq$mom_hs - mean(kidiq$mom_hs)
#'
#' coFit <- stan_glm(kid_std ~ hs_cent, data = kidiq,
#'                   family = gaussian(), prior = coef_prior,
#'                   prior_intercept = t_prior,
#'                   seed = 1776, chains = 2
#' )
#' iqFit <- stan_glm(kid_std ~ mom_std + hs_cent, data = kidiq,
#'                   family = gaussian(),
#'                   prior = coef_prior, prior_intercept = t_prior,
#'                   seed = 1776, chains = 2
#' )
#'
#'
#' coLoo <- loo(iqFit, save_psis = TRUE)
#' iqLoo <- loo(coFit, save_psis = TRUE)
#'
#'
#' plot_loo_dif(kidiq$mom_iq, coLoo, iqLoo, group = kidiq$mom_hs,
#'              alpha = .5, jitter = c(.1, .1)
#'              ) +
#'   ggplot2::geom_smooth() +
#'   ggplot2::xlab("IQ of Mother") +
#'   ggplot2::scale_colour_manual(values=cbPalette)
#' 

plot_loo_dif <- 
  function(y,
           psis_object_1,
           psis_object_2,
           ...,
           group = NULL,
           outlier_thresh = NULL,
           size = 1,
           alpha = 1,
           jitter = 0,
           sort_by_group = FALSE
  ){
    
    # Adding a 0 at the end lets users provide a single number as input.
    # In this case, only horizontal jitter is applied.
    jitter <- c(jitter, 0) 
    
    elpdDif <- psis_object_1$pointwise[, "elpd_loo"] - 
               psis_object_2$pointwise[, "elpd_loo"]

    
    if (sort_by_group){
      if (identical(group, NULL) || !identical(y, 1:length(y))){
        stop("ERROR: sort_by_group should only be used for grouping categorical 
             variables, then plotting them with an arbitrary index. You can
             create such an index using `1:length(data)`.
             ")
      }
      
      ordering <- order(group)
      elpdDif <- elpdDif[ordering]
      group <- group[ordering]
      
    }
    
    
    plot <- ggplot2::ggplot(mapping=aes(y, elpdDif)) +
            ggplot2::geom_hline(yintercept=0) + 
            ggplot2::xlab(ifelse(sort_by_group, "y", "Index")) +
            ggplot2::ylab(expression(ELPD[i][1] - ELPD[i][2])) + 
            ggplot2::labs(color = "Groups")
    

    
    if (identical(group, FALSE)){ 
      # Don't color by group if no groups are passed
      plot <- plot +  
              ggplot2::geom_jitter(width = jitter[1], height = jitter[2], 
                                   alpha = alpha, size = size
                                  ) 
    }
    else{
      # If group is passed, use color
      plot <- plot +  
          ggplot2::geom_jitter(aes(color = factor(group)),
                               width = jitter[1], height = jitter[2], 
                               alpha = alpha, size = size
                              ) 
    }
    
    if (!identical(outlier_thresh, NULL)){
      # Flag outliers
      is_outlier <- elpdDif > outlier_thresh
      index <- 1:length(y)
      outlier_labs <- index[is_outlier]
      
      plot <- plot + ggplot2::annotate("text",
                                       x = y[is_outlier],
                                       y = elpdDif[outlier_labs],
                                       label = outlier_labs,
                                       size = 4
                                      )
              
      
    }
    
    return(plot)
  }

