
<!-- README.md is generated from README.Rmd. Please edit that file -->

# hyper-mpas

<!-- badges: start -->
<!-- badges: end -->

The goal of hyper-mpas is to expand on the ideas of Marshall et
al. (2021) and explore how hyperallometry affects the outcomes of MPA in
a range of circumstances.

See run\_hyper\_mpas.R for current work.

This project is set up with
[`renv`](https://rstudio.github.io/renv/articles/renv.html) to manage
package dependencies. Inside R (and with your working directory set
correctly) run `renv::restore()`. Follow all prompts. This will install
the correct versions of all the packages needed to replicate our
results. Packages are installed in a stand-alone project library for
this paper, and will not affect your installed R packages anywhere else.

Setup

``` r
library(marlin)

library(tidyverse)
#> ── Attaching packages ─────────────────────────────────────── tidyverse 1.3.1 ──
#> ✓ ggplot2 3.3.5     ✓ purrr   0.3.4
#> ✓ tibble  3.1.3     ✓ dplyr   1.0.7
#> ✓ tidyr   1.1.3     ✓ stringr 1.4.0
#> ✓ readr   2.0.1     ✓ forcats 0.5.1
#> ── Conflicts ────────────────────────────────────────── tidyverse_conflicts() ──
#> x dplyr::filter() masks stats::filter()
#> x dplyr::lag()    masks stats::lag()

options(dplyr.summarise.inform = FALSE)

theme_set(marlin::theme_marlin())

resolution <- 20 # resolution is in squared patches, so 20 implies a 20X20 system, i.e. 400 patches 

years <- 200

seasons <- 1

time_step <- 1 / seasons

steps <- years * seasons

adult_movement_sigma <- resolution


recruit_movement_sigma <- resolution * 10

rec_form <- 1

hyper_expo <- 1.5

steepness <- 0.6

fished_depletion <- 0.15

mpa <- expand_grid(x = 1:resolution, y= 1:resolution) %>% 
  mutate(mpa = (1:nrow(.)) < (0.7 * resolution^2))

mpa %>% 
  ggplot(aes(x,y,fill = mpa)) + 
  geom_tile()
```

![](README_files/figure-gfm/unnamed-chunk-2-1.png)<!-- -->

``` r
find_refs <- function(log_emult = 0, fauna, fleet, years = 100, use = "opt"){

emult <- exp(log_emult)

tmp_fleet <- fleet

tmp_fleet[[1]]$base_effort <- tmp_fleet[[1]]$base_effort  * emult
  
sim <- simmar(fauna = fauna,
              fleets = tmp_fleet,
              years = years)

eq <- sim[[length(sim)]][[1]]

yield <- sum(eq$c_p_a)

f <- mean((eq$e_p_fl * tmp_fleet[[1]]$metiers[[1]]$catchability)[,1])
refs <- tibble(e_msy = tmp_fleet[[1]]$base_effort, fmsy = f, msy = yield, ssb_msy = sum(eq$ssb_p_a),
               ssb_msy_to_ssb0 = sum(eq$ssb_p_a) / eq$ssb0, bmsy_to_b0 = sum(eq$b_p_a) / fauna[[1]]$b0)

if (use == "opt"){
  out <- -yield
} else {
  
  out <- refs
}

return(out)

}
```

First a simulation with isometry

``` r

iso_critter <- 
  list(
    "bigeye" = create_critter(
      common_name = "bigeye tuna",
      adult_movement = 0,
      adult_movement_sigma = adult_movement_sigma,
      recruit_movement = 0,
      recruit_movement_sigma = recruit_movement_sigma,
      rec_form = rec_form,
      seasons = seasons,
      fished_depletion = fished_depletion,
      resolution = resolution,
      steepness = steepness,
      fec_expo = 1,
      ssb0 = 1000
    )
  )
#> ══  1 queries  ═══════════════
#> 
#> Retrieving data for taxon 'bigeye tuna'
#> ✔  Found:  bigeye+tuna[Common Name]
#> ══  Results  ═════════════════
#> 
#> ● Total: 1 
#> ● Found: 1 
#> ● Not Found: 0

iso_critter$bigeye$plot()
```

![](README_files/figure-gfm/unnamed-chunk-4-1.png)<!-- -->

``` r

iso_fleet <- list("longline" = create_fleet(
    list("bigeye" = Metier$new(
        critter = iso_critter$bigeye,
        price = 10,
        sel_form = "logistic",
        sel_start = .1,
        sel_delta = .01,
        catchability = 0,
        p_explt = 1
      )
    ),
    base_effort = 10000 * resolution ^ 2
  )
)

a <- Sys.time()

iso_fleet <- tune_fleets(iso_critter, iso_fleet, tune_type = "depletion") 

Sys.time() - a
#> Time difference of 3.54111 secs

ref_opt <- nlminb(log(1e-3), find_refs, fauna = iso_critter, fleet = iso_fleet)

iso_refpoints <- find_refs(ref_opt$par, fauna = iso_critter, fleet = iso_fleet, use = "refs")
 
iso_critter$bigeye$plot()
```

![](README_files/figure-gfm/unnamed-chunk-4-2.png)<!-- -->

``` r
iso_ref_check <- tibble(emult = seq(1e-3,2, length.out = 40)) %>% 
  mutate(tmp = map(log(emult), ~find_refs(.x, fauna = iso_critter, fleet = iso_fleet, use = "refs"))) %>% 
  mutate(type = "Isometric")

iso_ref_check %>% 
  unnest(cols = tmp) %>% 
  ggplot(aes(bmsy_to_b0, msy)) + 
  geom_point() + 
  geom_vline(aes(xintercept = exp(ref_opt$par)))
```

![](README_files/figure-gfm/unnamed-chunk-4-3.png)<!-- -->

``` r
a <- Sys.time()

iso_sim <- simmar(fauna = iso_critter,
                  fleets = iso_fleet,
                  years = years,
                  mpas = list(locations = mpa,
                              mpa_year = floor(years * .5)))

Sys.time() - a
#> Time difference of 0.6312251 secs


proc_iso_sim <- process_marlin(sim = iso_sim, time_step = time_step, keep_age = FALSE)

plot_marlin(proc_iso_sim)
```

![](README_files/figure-gfm/unnamed-chunk-4-4.png)<!-- -->

``` r
plot_marlin(proc_iso_sim, plot_var = "c")
```

![](README_files/figure-gfm/unnamed-chunk-4-5.png)<!-- -->

Now, a simulation with hyperallometry

``` r
hyper_critter <- 
  list(
    "bigeye" = create_critter(
      common_name = "bigeye tuna",
      adult_movement = 0,
      adult_movement_sigma = adult_movement_sigma,
      recruit_movement = 0,
      recruit_movement_sigma = recruit_movement_sigma,
      rec_form = rec_form,
      seasons = seasons,
      fished_depletion = fished_depletion,
      resolution = resolution,
      steepness = steepness,
      fec_expo = hyper_expo,
      ssb0 = 1000
    )
  )
#> ══  1 queries  ═══════════════
#> 
#> Retrieving data for taxon 'bigeye tuna'
#> ✔  Found:  bigeye+tuna[Common Name]
#> ══  Results  ═════════════════
#> 
#> ● Total: 1 
#> ● Found: 1 
#> ● Not Found: 0

hyper_critter$bigeye$plot()
```

![](README_files/figure-gfm/unnamed-chunk-5-1.png)<!-- -->

``` r
hyper_fleet <- list("longline" = create_fleet(
    list("bigeye" = Metier$new(
        critter = hyper_critter$bigeye,
        price = 10,
        sel_form = "logistic",
        sel_start = .1,
        sel_delta = .01,
        catchability = 0,
        p_explt = 1
      )
    ),
    base_effort = 10000*resolution ^ 2
  )
)

a <- Sys.time()

hyper_fleet <- tune_fleets(hyper_critter, hyper_fleet, tune_type = "depletion") 

Sys.time() - a
#> Time difference of 3.152861 secs

ref_opt <- nlminb(log(1e-3), find_refs, fauna = hyper_critter, fleet = hyper_fleet)

hyper_refpoints <- find_refs(ref_opt$par, fauna = hyper_critter, fleet = hyper_fleet, use = "refs")

hyper_ref_check <- tibble(emult = seq(1e-3,2, length.out = 40)) %>% 
  mutate(tmp = map(log(emult), ~find_refs(.x, fauna = hyper_critter, fleet = hyper_fleet, use = "refs"))) %>%
  mutate(type = "Hyperallometric")

hyper_ref_check %>% 
  unnest(cols = tmp) %>% 
  ggplot(aes(emult, msy)) + 
  geom_point() + 
  geom_vline(aes(xintercept = exp(ref_opt$par)))
```

![](README_files/figure-gfm/unnamed-chunk-5-2.png)<!-- -->

``` r

a <- Sys.time()

hyper_sim <- simmar(fauna = hyper_critter,
                  fleets = hyper_fleet,
                  years = years,
                  mpas = list(locations = mpa,
                              mpa_year = floor(years * .5)))

Sys.time() - a
#> Time difference of 0.451395 secs


proc_hyper_sim <- process_marlin(sim = hyper_sim, time_step = time_step, keep_age = FALSE)

plot_marlin(proc_hyper_sim)
```

![](README_files/figure-gfm/unnamed-chunk-5-3.png)<!-- -->

``` r
plot_marlin(proc_hyper_sim, plot_var = "c")
```

![](README_files/figure-gfm/unnamed-chunk-5-4.png)<!-- -->

And compare

``` r
fmsy_ratio <- round(iso_refpoints$fmsy / hyper_refpoints$fmsy,2)

yield_curve <- iso_ref_check %>% 
  bind_rows(hyper_ref_check) %>% 
  unnest(cols = tmp)


yield_curve %>% 
  ggplot(aes(fmsy, msy, color = type)) + 
  geom_line() + 
  scale_x_continuous(name = "F") + 
  scale_y_continuous(name  = "Yield") + 
  labs(caption = glue::glue("Isometric FMSY is {fmsy_ratio} times Hyperallometric FMSY"))
```

![](README_files/figure-gfm/unnamed-chunk-6-1.png)<!-- -->

``` r
iso_critter$bigeye$plot()
```

![](README_files/figure-gfm/unnamed-chunk-8-1.png)<!-- -->

``` r
hyper_critter$bigeye$plot()
```

![](README_files/figure-gfm/unnamed-chunk-8-2.png)<!-- -->

``` r
plot_marlin("Isometric" = proc_iso_sim, "Hyperallometric" = proc_hyper_sim)
```

![](README_files/figure-gfm/unnamed-chunk-8-3.png)<!-- -->

``` r
plot_marlin("Isometric" = proc_iso_sim, "Hyperallometric" = proc_hyper_sim, plot_var = "c")
```

![](README_files/figure-gfm/unnamed-chunk-9-1.png)<!-- -->

``` r
iso_cmsy <- proc_iso_sim$fleets %>% 
  group_by(step, critter) %>% 
  summarise(catch = sum(catch, na.rm = TRUE)) %>% 
  mutate(cmsy = catch / iso_refpoints$msy) %>% 
  mutate(type = "Isometric")

hyper_cmsy <- proc_hyper_sim$fleets %>% 
  group_by(step, critter) %>% 
  summarise(catch = sum(catch, na.rm = TRUE)) %>% 
  mutate(cmsy = catch / hyper_refpoints$msy) %>% 
    mutate(type = "Hyperallometric")


cmsys <- iso_cmsy %>% 
  bind_rows(hyper_cmsy)


cmsys %>% 
  ggplot(aes(step, cmsy, color = type)) + 
  geom_point()
```

![](README_files/figure-gfm/unnamed-chunk-10-1.png)<!-- -->
