## Code to conduct FR-BCA 

## library(dplyr)
## library(magrittr)

## Inputs (list)
## - Model name(s)
## - Analysis parameters {delta, T, losses}
##    - Baseline
##    - Ranges for sensitivity analysis
## - Cost data for each model
## - Data from performance assessment for each model
##    - Format?

## Outputs
## - Table with baseline BCR + high/low values, for each parameter included in
## the sensitvity analysis
## - Columns:
##    - model/variant
##    - sensitivity parameter name
##    - bcr
##    - label (base, low, high)
## - Ready to print tables and plots?


## for each parameter in inputs$parameters$sensitivity:
## - update "base" list with {low, high} values for that parameter
## - compute the BCRs and add a column indicating the parameter that is varied

## So there will be a single function to compute pv_cost and pv_loss/pv_benefit
## The inputs will either be base "as is" or base with a sensitivity parameter varied
## There will be a wrapper function that calls each function iteratively for each sensitivity

## The other thing to note is this will be done one model at a time
## so the input tables need to be separated, so that there is only one status quo
## and all other rows are compared to that status quo
## There will be a wrapper that separates the models in the background,
## does all the computations for each model, and then assembles back together


preprocess_model <- function(eal, cost, p) {
    ## Purpose:
    ## join eal and cost tables
    ## compute height (floors) and total area
    ## separate models by height
    ## 
    models <- list()
    dat <- eal %>%
        dplyr::left_join(cost, by=c('model', 'intervention')) %>%
        dplyr::mutate(total_floors=as.numeric(gsub('(.*-)(\\d{1,2})$', '\\2', model))) %>%
        dplyr::mutate(total_area=p$floor_area*total_floors)
    floors <- dat %>% dplyr::distinct(total_floors) %>% pull()
    for (i in 1:length(floors)) {
        models[[i]] <- dat %>%
            dplyr::filter(total_floors == floors[i])
    }
    return(models)
}


preprocess_cost <- function(model) {
    ## Purpose:
    ## ignore rows without nonstructural costs
    return(model %>% dplyr::filter(!is.na(c_ns)))
}


## PV(total cost) = structural cost + PV(nonstructural cost)
pv_cost <- function(model, params) {
    ## Purpose:
    ## Calculate present value ns and total costs 
    p = params$parameters$base
    return(
        model %>%
        ## filter out missing NS cost
        ## preprocess_cost() %>%
        dplyr::mutate(
                   pv_s=c_s,
                   pv_ns=c_ns*(1-p$delta)^(2018-2011)) %>%
        dplyr::mutate(pv_total=pv_s+pv_ns+pv_ns/(1+p$delta)^25)
    )
}

## compute cost delta
## pv cost formula:
## pv_cost = pv_s + pv_ns + (pv_ns / (1+delta))
## pv_s = structural cost
## pv_ns = nonstructural cost * (1-delta)^(2018-2011)
pv_dcost <- function(model, params) {
    ## Purpose:
    ## Calculate present value cost deltas, relative to status quo
    ## NB: assumes cost table filtered to only include rows with ns cost
    return(
        model %>%
        pv_cost(params) %>%
        dplyr::mutate(
                   cost_diff=pv_total - pv_total[intervention == 0],
                   cost_delta=(pv_total/pv_total[intervention == 0]) - 1
                      ) 
    )
}


## loss formulas:
## total_area = (floor area) * number of stories
## displacement = displace_per_area * tenant_per_area * reocc_days * total_area
## bi = (1 - recapture) * bi_per_area * fr_days * total_area
## ri = ri_per_area * fr_days * total_area
pv_loss <- function(model, p) {
    ## Purpose:
    ## Calculate losses for each model
    return(
        model %>%
        ## NB: assumes total_area column has been created
        dplyr::mutate(
                   displacement=p$displacement*p$tenant*re_occupancy_time*total_area,
                   business_income=(1 - p$recapture)*p$bi*functional_recovery_time*total_area,
                   rental_income=p$ri*functional_recovery_time*total_area
                   )
        )
}

## pv benefits formula:
## (avoided losses) * ( (1 - (1+delta)^(-T)) / delta)
pv_benefit <- function(model, params, label='base') {
    ## Purpose:
    ## Calculate present value avoided losses, relative to status quo
    join_cols = c('model', 'intervention')
    loss_cols = c('repair_cost', 'displacement', 'business_income', 'rental_income')
    p = params$parameters$base
    m = model %>%
        ## NB: assumes total_area column has been created
        pv_loss(p) %>%
        dplyr::select(all_of(c(join_cols, loss_cols))) %>%
        dplyr::rowwise() %>%
        dplyr::mutate(loss_total=sum(across(all_of(loss_cols)))) %>%
        dplyr::ungroup() %>%
        dplyr::mutate(delta_loss=loss_total[intervention == 0] - loss_total) %>%
        dplyr::mutate(benefit=delta_loss * ((1 - (1+p$delta)^(-p$T))/p$delta))
    return(
        model %>%
        dplyr::select(!repair_cost) %>%
        dplyr::left_join(m, by=join_cols)
    )
}


bcr <- function(model, params, label='base') {
    ## Purpose:
    ## Compute BCR and NPV
    ## calls pv_benefit and pv_dcost
    ## TODO: add column with label {base, high, low}
    model <- pv_benefit(model, params)
    model <- pv_dcost(model, params)
    return(model %>%
           dplyr::mutate(
                      bcr=benefit/cost_diff,
                      npv=benefit-cost_diff,
                      label=label)
           )
}


set_params <- function(params, param, bound='low') {
    ## Purpose:
    ## Reset baseline parameter to one of {low, high}
    p <- params
    p[['parameters']][['base']][[param]] <- p[['parameters']][['sensitivity']][[param]][[bound]]
    return(p)
}


sensitivity <- function(model, params) {
    ## Purpose:
    ## Compute BCR and NPV, using low/high values of parameters
    ## --- ##
    ## get list of sensitivity parameters
    s <- params[['parameters']][['sensitivity']]
    ## store calculations
    m <- list()
    ## iterate low/high
    for (hi_low in c('low', 'high')) {
        ## iterate over parameters
        for (n in names(s)) {          
            p <- set_params(params, n, bound=hi_low)
            m[[paste(n, hi_low, sep='-')]] <- bcr(model, p, label=hi_low) %>%
                dplyr::mutate(parameter=n)
        }
    }
    return(dplyr::bind_rows(m))
}


bca <- function(model, params) {
    ## Purpose:
    ## Wrapper that computes (1) baseline BCA and (2) sensitivity analysis
    ## --- ##
    ## first pass baseline bca
    m_b <- bcr(model, params)
    ## second pass sensitivity
    m_s <- sensitivity(model, params)
    ## TODO: pivot_wider on bcr -> bcr_base, bcr_low, bcr_high
    return(dplyr::bind_rows(m_b, m_s))
}


frbca <- function(eal, cost, params) {
    ## Purpose:
    ## Wrapper that conducts BCA for each set of models in list
    models <- preprocess_model(eal, cost, params[['parameters']][['base']])
    for (i in 1:length(models)) {
        models[[i]] <- bca(models[[i]], params)
    }
    ## TODO: filter out NaN as base case?
    ## TODO: filter out NA for missing cost?
    models <- dplyr::bind_rows(models)
    return(models)
}

plot_frbca <- function(outout, n_floors=4, system='RCMF') {
  ## Purpose: post-process data and generate plot for sensitivity analysis
  plot_df <- output %>%
    dplyr::filter(!is.na(bcr)) %>%
    dplyr::filter(total_floors == n_floors) %>%
    dplyr::select(model, bcr, label, parameter)
  base <- plot_df %>%
    dplyr::filter(label == 'base') %>%
    dplyr::select(!c(label, parameter))
  sen <-plot_df %>%
    dplyr::filter(label != 'base') %>%
    tidyr::pivot_wider(names_from=label, values_from=bcr) %>%
    dplyr::left_join(base, by='model') %>%
    dplyr::rename(bcr_low=low, bcr_high=high)
  ## generate plot
  label_begin <- 'Sensitivity Analysis: Benefit-cost ratios for'
  label_end <- 'archetypes, relative to baseline ASCE 7-16 design.'
  plot.sen <- sen %>%
    ggplot() +
    geom_segment(aes(x=parameter, xend=parameter, y=bcr_low, yend=bcr_high),
                 linewidth = 5, colour = "red", alpha = 0.6) +
    geom_segment(aes(x=parameter, xend=parameter, y=bcr-0.001, yend=bcr+0.001),
                 linewidth = 5, colour = "black") +
    geom_hline(yintercept=1, colour='red') +
    coord_flip() +
    facet_wrap(~model, ncol=1) +
    ## geom_hline(data=rcmf, aes(yintercept=bcr)) +
    theme_light() +
    theme(legend.position='bottom') +
    labs(
      title=paste(label_begin, paste0(total_floors, '-story'), system, label_end),
      x='Parameter',
      y='Benefit-cost ratio')
return(plot.sen)
}
