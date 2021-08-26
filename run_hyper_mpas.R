library(marlin)

library(tidyverse)

library(here)

library(furrr)

options(dplyr.summarise.inform = FALSE)

theme_set(marlin::theme_marlin())



# global options ----------------------------------------------------------


resolution <-
  10 # resolution is in squared patches, so 20 implies a 20X20 system, i.e. 400 patches

years <- 200

seasons <- 1

time_step <- 1 / seasons

steps <- years * seasons

workers <- 8

future::plan(future::multisession, workers = workers)

on.exit(future::plan(future::sequential))
run_experiments <- FALSE


# do gadus case study -----------------------------------------------------



# setup experiments -------------------------------------------------------

n <- 250

experiments <- tibble(
  scientific_name = sample(c("Gadus morhua","Thunnus obesus"), n, replace = TRUE),
  adult_move = runif(n, 0, 2),
  rec_move = runif(n, 0, 2),
  hyper = runif(n, 1, 3),
  steepness = runif(n, .21, .8),
  fished_depletion = runif(n, .1, .75),
  rec_form = sample(c(0:3), n, replace = TRUE),
  p_sel = sample(c(0.1,1), n, replace = TRUE)
) %>%
  mutate(xid = 1:nrow(.)) %>%
  mutate(marshall = TRUE)

experiments <- experiments %>%
  bind_rows(experiments %>% mutate(marshall = FALSE))

find_refs <-
  function(log_emult = 0,
           fauna,
           fleet,
           mpa,
           years = 100,
           use = "opt") {
    emult <- exp(log_emult)
    
    tmp_fleet <- fleet
    
    tmp_fleet[[1]]$base_effort <-
      tmp_fleet[[1]]$base_effort  * emult
    
    
    sim <- simmar(
      fauna = fauna,
      fleets = tmp_fleet,
      years = years,
      mpas = list(locations = mpa,
                  mpa_year = floor(years * .5))
    )
    
    eq <- sim[[length(sim)]][[1]]
    
    yield <- sum(eq$c_p_a)
    
    f <-
      mean((eq$e_p_fl * tmp_fleet[[1]]$metiers[[1]]$catchability)[, 1])
    refs <-
      tibble(
        e_msy = tmp_fleet[[1]]$base_effort,
        fmsy = f,
        msy = yield,
        ssb_msy = sum(eq$ssb_p_a),
        b_msy = sum(eq$b_p_a, na.rm = TRUE),
        ssb_msy_to_ssb0 = sum(eq$ssb_p_a, na.rm = TRUE) / eq$ssb0,
        bmsy_to_b0 = sum(eq$b_p_a, na.rm = TRUE) / fauna[[1]]$b0
      )
    
    if (use == "opt") {
      out <- -yield
    } else {
      out <- refs
    }
    
    return(out)
    
  }



run_experiment <-
  function(adult_move,
           rec_move,
           hyper,
           steepness,
           fished_depletion,
           rec_form,
           marshall,
           p_sel,
           scientific_name,
           xid) {
    mpas <- expand_grid(x = 1:resolution, y = 1:resolution) %>%
      mutate(mpa = FALSE)
    
    iso_critter <-
      list(
        "bigeye" = create_critter(
          scientific_name = scientific_name,
          adult_movement = 0,
          adult_movement_sigma = adult_move * resolution / 2,
          recruit_movement = 0,
          recruit_movement_sigma = rec_move * resolution / 2,
          rec_form = rec_form,
          seasons = seasons,
          fished_depletion = fished_depletion,
          resolution = resolution,
          steepness = steepness,
          fec_expo = 1
        )
      )

    iso_fleet <- list("longline" = create_fleet(list(
      "bigeye" = Metier$new(
        critter = iso_critter$bigeye,
        price = 10,
        sel_form = "logistic",
        sel_start = p_sel,
        sel_delta = .01,
        catchability = 0,
        p_explt = 1
      )
    ),
    base_effort = 10000 * resolution ^ 2))
    
    iso_fleet <-
      tune_fleets(iso_critter, iso_fleet, tune_type = "depletion")
    
    ref_opt <-
      nlminb(
        log(mean(iso_critter[[1]]$m_at_age) * (sqrt((4 * iso_critter[[1]]$steepness) / (1 - iso_critter[[1]]$steepness)
        ) - 1)),
        find_refs,
        fauna = iso_critter,
        fleet = iso_fleet,
        mpa = mpas
      )
    
    base_iso_refpoints <-
      find_refs(
        ref_opt$par,
        fauna = iso_critter,
        fleet = iso_fleet,
        mpa = mpas,
        use = "refs"
      )
    
    
    # repeat but for hyper
    
    
    hyper_critter <-
      list(
        "bigeye" = create_critter(
          scientific_name = scientific_name,
          adult_movement = 0,
          adult_movement_sigma = adult_move * resolution / 2,
          recruit_movement = 0,
          recruit_movement_sigma = rec_move * resolution / 2,
          rec_form = rec_form,
          seasons = seasons,
          fished_depletion = fished_depletion,
          resolution = resolution,
          steepness = steepness,
          fec_expo = hyper
        )
      )

    hyper_fleet <- list("longline" = create_fleet(list(
      "bigeye" = Metier$new(
        critter = hyper_critter$bigeye,
        price = 10,
        sel_form = "logistic",
        sel_start = p_sel,
        sel_delta = .01,
        catchability = 0,
        p_explt = 1
      )
    ),
    base_effort = 10000 * resolution ^ 2))
    
    hyper_fleet <-
      tune_fleets(hyper_critter, hyper_fleet, tune_type = "depletion")
    
    ref_opt <-
      nlminb(
        log(1e-3),
        find_refs,
        fauna = hyper_critter,
        fleet = hyper_fleet,
        mpa = mpas
      )
    
    
    base_hyper_refpoints <-
      find_refs(
        ref_opt$par,
        fauna = hyper_critter,
        fleet = hyper_fleet,
        mpa = mpas,
        use = "refs"
      )
    
    # do MPAs
    #
    run_mpa <-
      function(prop_mpa,
               resolution,
               iso_critter,
               iso_fleet,
               hyper_critter,
               hyper_fleet,
               marshall = FALSE) {
        mpa <- expand_grid(x = 1:resolution, y = 1:resolution) %>%
          mutate(mpa = (1:nrow(.)) < (prop_mpa * resolution ^ 2))
        
        if (marshall) {
          # following Marshall redefine F to be Fmsy conditional on MPA
          ref_opt <-
            nlminb(
              ref_opt$par,
              find_refs,
              fauna = iso_critter,
              fleet = iso_fleet,
              mpa = mpa
            )
          
          iso_fleet[[1]]$base_effort <-
            iso_fleet[[1]]$base_effort  * exp(ref_opt$par)
          
          ref_opt <-
            nlminb(
              ref_opt$par,
              find_refs,
              fauna = hyper_critter,
              fleet = hyper_fleet,
              mpa = mpa
            )
          
          hyper_fleet[[1]]$base_effort <-
            hyper_fleet[[1]]$base_effort  * exp(ref_opt$par)
          
        }
        
        iso_sim <- simmar(
          fauna = iso_critter,
          fleets = iso_fleet,
          years = years,
          mpas = list(locations = mpa,
                      mpa_year = floor(years * .5))
        )
        
        
        hyper_sim <- simmar(
          fauna = hyper_critter,
          fleets = hyper_fleet,
          years = years,
          mpas = list(locations = mpa,
                      mpa_year = floor(years * .5))
        )
        
        iso_end <- iso_sim[[length(iso_sim)]]
        
        fs <-
          (iso_end[[1]]$e_p_fl[mpa$mpa == FALSE, 1] * iso_fleet[[1]]$metiers[[1]]$spatial_catchability[mpa$mpa == FALSE])
        
        fs <- mean(fs[, 1])
        
        iso_res <-
          tibble(
            biomass = sum(iso_end[[1]]$b_p_a, na.rm = TRUE),
            ssb = sum(iso_end[[1]]$ssb_p_a, na.rm = TRUE),
            catch = sum(iso_end[[1]]$c_p_a, na.rm = TRUE)
          ) %>%
          mutate(
            catch_msy = catch / base_iso_refpoints$msy,
            b_bmsy = biomass / base_iso_refpoints$b_msy,
            ssb_ssbmsy = ssb / base_iso_refpoints$ssb_msy,
            f = fs,
            f_fmsy = fs / base_iso_refpoints$fmsy,
            type = "Isometric"
          )
        
        # process hyper
        
        hyper_end <- hyper_sim[[length(hyper_sim)]]
        
        fs <-
          (hyper_end[[1]]$e_p_fl[mpa$mpa == FALSE, 1] * hyper_fleet[[1]]$metiers[[1]]$spatial_catchability[mpa$mpa == FALSE])
        
        fs <- mean(fs[, 1])
        
        hyper_res <-
          tibble(
            biomass = sum(hyper_end[[1]]$b_p_a, na.rm = TRUE),
            ssb = sum(hyper_end[[1]]$ssb_p_a, na.rm = TRUE),
            catch = sum(hyper_end[[1]]$c_p_a, na.rm = TRUE)
          ) %>%
          mutate(
            catch_msy = catch / base_hyper_refpoints$msy,
            b_bmsy = biomass / base_hyper_refpoints$b_msy,
            ssb_ssbmsy = ssb / base_hyper_refpoints$ssb_msy,
            f = fs,
            f_fmsy = fs / base_hyper_refpoints$fmsy,
            type = "Hyperallometric"
          )
        
        mpa_result <- iso_res %>%
          bind_rows(hyper_res)
        
        return(mpa_result)
      } # close run_mpa
    
    # a <- Sys.time()
    mpa_results <- tibble(prop_mpa = seq(0, .9, by = .1)) %>%
      mutate(
        tmp = future_map(
          prop_mpa,
          run_mpa,
          resolution = resolution,
          iso_critter = iso_critter,
          iso_fleet = iso_fleet,
          hyper_critter = hyper_critter,
          hyper_fleet = hyper_fleet,
          marshall = marshall,
          .options = furrr_options(seed = 42)
        )
      )
    
    # Sys.time() - a
    mpa_results <- mpa_results %>%
      unnest(cols = tmp)
    
    # mpa_results %>%
    #   ggplot(aes(prop_mpa, catch_msy, color = type)) +
    #   geom_line()
    message(glue::glue("{scales::percent(xid/(n * 2))} Done"))
    return(mpa_results)
    
  }

if (run_experiments) {
  a <- Sys.time()
  experiment_results <- experiments %>%
    mutate(tmp = pmap(
      list(
        adult_move = adult_move,
        rec_move = rec_move,
        hyper = hyper,
        steepness = steepness,
        fished_depletion = fished_depletion,
        rec_form = rec_form,
        marshall = marshall,
        p_sel = p_sel,
        scientific_name = scientific_name,
        xid = xid
      ),
      safely(run_experiment)
    ))
  
  b <- Sys.time() - a
  
  message(glue::glue("{n} Experiments took {round(b)} minutes"))
  
  write_rds(experiment_results,
            here("results", "hyper_experiments.rds"))
} # run experiments
experiment_results <-
  read_rds(here("results", "hyper_experiments.rds"))

# process results ---------------------------------------------------------


experiment_worked <-
  map(experiment_results$tmp, "error") %>% map_lgl(is.null)

experiment_results <- experiment_results %>%
  filter(experiment_worked) %>%
  mutate(tmp = map(tmp, "result")) %>%
  unnest(cols = tmp)

tmp <- experiment_results %>%
  select(marshall, xid, prop_mpa, catch_msy, type) %>%
  left_join(experiments, by = c("xid", "marshall")) %>%
  pivot_wider(values_from = catch_msy, names_from = "type")

experiment_results %>% 
  mutate(marshall = ifelse(marshall == TRUE, "F Optimized for MPA", "Constant Effort")) %>% 
  # filter(marshall == FALSE) %>% 
  ggplot(aes(prop_mpa, catch_msy, color = 1 - exp(-f))) + 
  geom_hline(aes(yintercept = 1)) +
  geom_jitter(alpha = 0.65) +
  scale_color_viridis_c(labels = scales::percent, name= "Mortality Rate Outside MPA") + 
  scale_x_continuous(name = "Percent MPA", labels = scales::percent) +
  scale_y_continuous("Catch relative to no-MPA MSY")+
  theme(legend.position = "top") + 
  facet_grid(type~marshall) + 
  labs(caption = "F set to MPA conditional Fmsy")


experiment_results %>% 
  mutate(marshall = ifelse(marshall == TRUE, "F Optimized for MPA", "Constant Effort")) %>% 
  # filter(marshall == FALSE) %>% 
  ggplot(aes(prop_mpa, catch_msy, color = adult_move)) + 
  geom_hline(aes(yintercept = 1)) +
  geom_jitter(alpha = 0.65) +
  scale_color_viridis_c(name = "Adult Movement") + 
  scale_x_continuous(name = "Percent MPA", labels = scales::percent) +
  scale_y_continuous("Catch relative to no-MPA MSY")+
  theme(legend.position = "top") + 
  facet_grid(type~marshall) + 
  labs(caption = "F set to MPA conditional Fmsy")

experiment_results %>% 
  mutate(marshall = ifelse(marshall == TRUE, "F Optimized for MPA", "Constant Effort")) %>% 
  # filter(marshall == FALSE) %>% 
  ggplot(aes(prop_mpa, catch_msy, color = rec_move)) + 
  geom_hline(aes(yintercept = 1)) +
  geom_jitter(alpha = 0.65) +
  scale_color_viridis_c(name = "Recruit Movement") + 
  scale_x_continuous(name = "Percent MPA", labels = scales::percent) +
  scale_y_continuous("Catch relative to no-MPA MSY")+
  theme(legend.position = "top") + 
  facet_grid(type~marshall) + 
  labs(caption = "F set to MPA conditional Fmsy")


experiment_results %>% 
  mutate(marshall = ifelse(marshall == TRUE, "F Optimized for MPA", "Constant Effort")) %>% 
  # filter(marshall == FALSE) %>% 
  ggplot(aes(prop_mpa, catch_msy, color = hyper)) + 
  geom_hline(aes(yintercept = 1)) +
  geom_jitter(alpha = 0.65) +
  scale_color_viridis_c(name = "Hyperallometric Exponent") + 
  scale_x_continuous(name = "Percent MPA", labels = scales::percent) +
  scale_y_continuous("Catch relative to no-MPA MSY")+
  theme(legend.position = "top") + 
  facet_grid(type~marshall) + 
  labs(caption = "F set to MPA conditional Fmsy")



experiment_results %>% 
  filter(marshall == TRUE) %>% 
  ggplot(aes(prop_mpa, catch_msy, color = hyper)) + 
  geom_hline(aes(yintercept = 1)) +
  geom_jitter(alpha = 0.65) +
  scale_color_viridis_c(name= "Recruit Movement Rate") + 
  scale_x_continuous(name = "Percent MPA", labels = scales::percent) +
  scale_y_continuous("Hyperallometric Catch relative to no-MPA MSY")+
  theme(legend.position = "top") + 
  facet_wrap(~type)

  
tmp %>%
  filter(marshall == FALSE) %>% 
  ggplot(aes(prop_mpa, Hyperallometric / Isometric, color = marshall)) +
  geom_hline(aes(yintercept = 1)) +
  geom_point() +
  scale_x_continuous(name = "Percent MPA", labels = scales::percent) +
  scale_y_continuous("Hyperallometric C/MSY / Isometric C/MSY")

tmp %>%
  ggplot(aes(prop_mpa, Hyperallometric, color = rec_move)) +
  geom_hline(aes(yintercept = 1)) +
  geom_point() +
  scale_x_continuous(name = "Percent MPA", labels = scales::percent) +
  scale_y_continuous("Hyperallometric C/MSY") +
  facet_wrap( ~ marshall) + 
  scale_color_viridis_c()

tmp %>%
  ggplot(aes(prop_mpa, adult_move, color = Hyperallometric)) +
  geom_hline(aes(yintercept = 1)) +
  geom_point() +
  scale_x_continuous(name = "Percent MPA", labels = scales::percent) +
  scale_y_continuous("Hyperallometric C/MSY") +
  facet_wrap( ~ marshall) + 
  scale_color_viridis_c()



tmp %>%
  ggplot(aes(prop_mpa, Isometric, color = marshall)) +
  geom_hline(aes(yintercept = 1)) +
  geom_point() +
  scale_x_continuous(name = "Percent MPA", labels = scales::percent) +
  scale_y_continuous("Isometric C/MSY")


