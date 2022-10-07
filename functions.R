# script to load which contains all the new functions used in the project

# this function is based on plotBiomassObservedVsModel and adds a lower boundary to selected size
plotBiomassObservedVsModelCustom <- function (object, species = NULL, ratio = FALSE, log_scale = TRUE, 
                                              return_data = FALSE, labels = TRUE, show_unobserved = FALSE) 
{
  if (is(object, "MizerSim")) {
    params = object@params
    n <- finalN(object)
  }
  else if (is(object, "MizerParams")) {
    params = object
    n <- initialN(params)
  }
  else {
    stop("You have not provided a valid mizerSim or mizerParams object.")
  }
  sp_params <- params@species_params
  species = valid_species_arg(object, species)
  if (length(species) == 0) 
    stop("No species selected, please fix.")
  row_select = match(species, sp_params$species)
  if (!"biomass_observed" %in% names(sp_params)) {
    stop("You have not provided values for the column 'biomass_observed' ", 
         "in the mizerParams/mizerSim object.")
  }
  else if (!is.numeric(sp_params$biomass_observed)) {
    stop("The column 'biomass_observed' in the mizerParams/mizerSim object", 
         " is not numeric, please fix.")
  }
  else {
    biomass_observed = sp_params$biomass_observed
  }
  
  cutoffLow <- sp_params$biomass_cutoffLow[row_select]
  if (is.null(cutoffLow)) {
    cutoffLow = rep(0, length(species))
  }
  else if (!is.numeric(cutoffLow)) {
    stop("params@species_params$biomass_cutoffLow is not numeric, \",\n                 \"please fix.")
  }
  cutoffLow[is.na(cutoffLow)] <- 0
  
  cutoffHigh <- sp_params$biomass_cutoffHigh[row_select]
  if (is.null(cutoffHigh)) {
    cutoffHigh = rep(0, length(species))
  }
  else if (!is.numeric(cutoffHigh)) {
    stop("params@species_params$biomass_cutoffHigh is not numeric, \",\n                 \"please fix.")
  }
  cutoffHigh[is.na(cutoffHigh)] <- 0
  
  sim_biomass = rep(0, length(species))
  for (j in 1:length(species)) {
    sim_biomass[j] = sum((n[row_select[j], ] * params@w * 
                            params@dw)[params@w >= cutoffLow[j] & cutoffHigh[j] >= params@w])
  }
  dummy = data.frame(species = species, model = sim_biomass, 
                     observed = biomass_observed[row_select]) %>% mutate(species = factor(species, 
                                                                                          levels = species), is_observed = !is.na(observed) & observed > 
                                                                           0, observed = case_when(is_observed ~ observed, !is_observed ~ 
                                                                                                     model), ratio = model/observed)
  if (sum(dummy$is_observed) == 0) {
    cat(paste("There are no observed biomasses to compare to model,", 
              "only plotting model biomasses.", sep = "\n"))
  }
  if (!show_unobserved) {
    dummy <- filter(dummy, is_observed)
  }
  if (return_data == TRUE) 
    return(dummy)
  tre <- round(sum(abs(1 - dummy$ratio)), digits = 3)
  caption <- paste0("Total relative error = ", tre)
  if (any(!dummy$is_observed)) {
    caption <- paste(caption, "\n Open circles represent species without biomass observation.")
  }
  if (ratio == FALSE) {
    gg <- ggplot(data = dummy, aes(x = observed, y = model, 
                                   colour = species, shape = is_observed)) + geom_abline(aes(intercept = 0, 
                                                                                             slope = 1), colour = "purple", linetype = "dashed", 
                                                                                         size = 1.3) + geom_point(size = 3) + labs(y = "model biomass [g]") + 
      coord_cartesian(ylim = range(dummy$model, dummy$observed))
  }
  else {
    gg <- ggplot(data = dummy, aes(x = observed, y = ratio, 
                                   colour = species, shape = is_observed)) + geom_hline(aes(yintercept = 1), 
                                                                                        linetype = "dashed", colour = "purple", 
                                                                                        size = 1.3) + geom_point(size = 3) + labs(y = "model biomass / observed biomass") + 
      coord_cartesian(ylim = range(dummy$ratio))
  }
  gg <- gg + labs(x = "observed biomass [g]", caption = caption) + 
    scale_colour_manual(values = getColours(params)[dummy$species]) + 
    scale_shape_manual(values = c(`TRUE` = 19, `FALSE` = 1)) + 
    guides(shape = "none")
  if (log_scale == TRUE & ratio == FALSE) {
    gg = gg + scale_x_log10() + scale_y_log10()
  }
  if (log_scale == TRUE & ratio == TRUE) {
    gg = gg + scale_x_log10()
  }
  if (labels == TRUE) {
    gg = gg + ggrepel::geom_label_repel(aes(label = species), 
                                        box.padding = 0.35, point.padding = 0.5, segment.color = "grey50", 
                                        show.legend = FALSE, max.overlaps = Inf, seed = 42)
  }
  gg
}

# adapting cutoff here too

calibrateBiomassCustom <- function (params) 
{
  if ((!("biomass_observed" %in% names(params@species_params))) || 
      all(is.na(params@species_params$biomass_observed))) {
    return(params)
  }
  no_sp <- nrow(params@species_params)
  
  cutoffLow <- params@species_params$biomass_cutoffLow
  if (is.null(cutoffLow)) 
    cutoffLow <- rep(0, no_sp)
  cutoffLow[is.na(cutoffLow)] <- 0
  
  cutoffHigh <- params@species_params$biomass_cutoffHigh
  if (is.null(cutoffHigh)) 
    cutoffHigh <- rep(0, no_sp)
  cutoffHigh[is.na(cutoffHigh)] <- 0
  
  observed <- params@species_params$biomass_observed
  observed_total <- sum(observed, na.rm = TRUE)
  sp_observed <- which(!is.na(observed))
  model_total <- 0
  for (sp_idx in sp_observed) {
    model_total <- model_total + sum((params@initial_n[sp_idx, 
    ] * params@w * params@dw)[params@w >= cutoffLow[sp_idx] & cutoffHigh[sp_idx] >= params@w])
  }
  scaleModel(params, factor = observed_total/model_total)
}


# adding the cutoff when calculating error

getErrorCustom <- function(vary, params, dat, tol = 0.001, 
                           timetorun = 10)
{
  params@species_params$R_max[1:9]<-10^vary[1:9]
  params@species_params$erepro[1:9]<-vary[10:18]
  params@species_params$interaction_resource[1:9] <- vary[19:27]
  
  params <- setParams(params)
  
  interaction <- params@interaction
  interaction[] <- matrix(vary[28:108],nrow = 9) # stop at 54 if looking only at 3 biggest species
  
  params <- setInteraction(params,interaction)
  
  params <- projectToSteady(params, distance_func = distanceSSLogN, 
                            tol = tol, t_max = 200, return_sim = F)
  
  sim <- project(params, t_max = timetorun, progress_bar = F)
  
  sim_biomass = rep(0, length(params@species_params$species))
  
  cutoffLow <- params@species_params$biomass_cutoffLow
  if (is.null(cutoffLow)) 
    cutoffLow <- rep(0, no_sp)
  cutoffLow[is.na(cutoffLow)] <- 0
  
  cutoffHigh <- params@species_params$biomass_cutoffHigh
  if (is.null(cutoffHigh)) 
    cutoffHigh <- rep(0, no_sp)
  cutoffHigh[is.na(cutoffHigh)] <- 0
  
  for (j in 1:length(sim_biomass)) {
    sim_biomass[j] = sum((sim@n[dim(sim@n)[1],j,] * params@w * 
                            params@dw)[params@w >= cutoffLow[j] & cutoffHigh[j] >= params@w])
  }
  
  
  pred <- log(sim_biomass)
  dat <- log(dat)
  discrep <- pred - dat
  discrep <- (sum(discrep^2))
  return(discrep)
}
