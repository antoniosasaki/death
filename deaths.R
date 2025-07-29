# 1) Install and load packages
if (!requireNamespace(c("demography","readxl","writexl","dplyr","tidyr","purrr","forecast"), quietly=TRUE)) {
  install.packages(c("demography", "readxl", "writexl", "dplyr", "tidyr", "purrr", "forecast"))
}
library(demography)
library(readxl)
library(writexl)
library(dplyr)
library(tidyr)
library(purrr)
library(forecast)

# 2) Import data (adjust paths as needed)
pop <- read_excel("~/data/popdata.xlsx")
is  <- read_excel("~/data/isdata.xlsx")
cvd <- read_excel("~/data/cvddata.xlsx")

# 3) Compute totals (male + female)
pop <- pop  %>% mutate(total_alive = m_alive + w_alive)
is  <- is   %>% mutate(total_death = m_death + w_death)
cvd <- cvd  %>% mutate(total_death = m_death + w_death)

# 4) Build exposure matrices by sex and total (rows=ages, cols=years)
exp_m <- with(pop,    tapply(m_alive,    list(l_age, year), sum))
exp_w <- with(pop,    tapply(w_alive,    list(l_age, year), sum))
exp_t <- with(pop,    tapply(total_alive, list(l_age, year), sum))

# 5) Build death matrices similarly
dth_m_IS  <- with(is,  tapply(m_death,    list(l_age, year), sum))
dth_w_IS  <- with(is,  tapply(w_death,    list(l_age, year), sum))
dth_t_IS  <- with(is,  tapply(total_death, list(l_age, year), sum))

dth_m_CVD <- with(cvd, tapply(m_death,    list(l_age, year), sum))
dth_w_CVD <- with(cvd, tapply(w_death,    list(l_age, year), sum))
dth_t_CVD <- with(cvd, tapply(total_death, list(l_age, year), sum))

# 6) Helper function to fit Lee–Carter model
enable_lr_fit <- function(deaths_matrix, exposure_matrix, disease, sex) {
  # Mortality rates
  mx <- deaths_matrix / exposure_matrix
  # Create demogdata object
  demo <- demogdata(
    data  = mx,
    pop   = exposure_matrix,
    ages  = as.numeric(rownames(mx)),
    years = as.numeric(colnames(mx)),
    type  = "mortality",
    label = paste(disease, sex),
    name  = paste(disease, sex)
  )
  # Fit Lee–Carter
  lca(demo)
}

# 7) Fit models for each disease and sex
fit_list <- list(
  IS_male    = enable_lr_fit(dth_m_IS,  exp_m,  "IS",  "male"),
  IS_female  = enable_lr_fit(dth_w_IS,  exp_w,  "IS",  "female"),
  IS_total   = enable_lr_fit(dth_t_IS,  exp_t,  "IS",  "total"),
  CVD_male   = enable_lr_fit(dth_m_CVD, exp_m,  "CVD", "male"),
  CVD_female = enable_lr_fit(dth_w_CVD, exp_w,  "CVD", "female"),
  CVD_total  = enable_lr_fit(dth_t_CVD, exp_t,  "CVD", "total")
)

# 8) Projection settings
duration_horizon <- 21  # years ahead (e.g., 2019 + 21 = 2040)
last_year <- as.numeric(colnames(exp_t))[ncol(exp_t)]
proj_years <- (last_year + 1):(last_year + duration_horizon)

# 9) Forecast mortality rates for each fit
rate_fc <- map(fit_list, ~ forecast(.x, h = duration_horizon, level = 95)$rate)

# 10) Forecast exposures using ETS and transpose to match rates
exp_fc <- list(
  male   = t(apply(exp_m,  1, function(x) forecast(ets(x), h = duration_horizon)$mean)),
  female = t(apply(exp_w,  1, function(x) forecast(ets(x), h = duration_horizon)$mean)),
  total  = t(apply(exp_t,  1, function(x) forecast(ets(x), h = duration_horizon)$mean))
)
# Assign dimnames: rows=ages, cols=proj_years
dimnames(exp_fc$male)   <- list(rownames(exp_m),   proj_years)
dimnames(exp_fc$female) <- list(rownames(exp_w),   proj_years)
dimnames(exp_fc$total)  <- list(rownames(exp_t),   proj_years)

# 11) Calculate absolute deaths by age/sex/disease/year
proj_list <- imap(rate_fc, function(rate_list, label) {
  parts   <- strsplit(label, "_")[[1]]
  disease <- parts[1]
  sex     <- parts[2]
  
  # Extract rate matrices
  rate_pt <- rate_list[[1]]
  rate_lo <- rate_list[[2]]
  rate_hi <- rate_list[[3]]
  
  # Corresponding exposure matrix
  exp_mat <- exp_fc[[tolower(sex)]]
  
  # Compute absolute deaths
  d_pt <- round(rate_pt * exp_mat)
  d_lo <- round(rate_lo * exp_mat)
  d_hi <- round(rate_hi * exp_mat)
  
  # Reshape into tibble
  tibble(
    disease = disease,
    sex     = sex,
    age     = rep(rownames(d_pt), times = ncol(d_pt)),
    year    = rep(proj_years, each = nrow(d_pt)),
    death   = as.vector(d_pt),
    lower95 = as.vector(d_lo),
    upper95 = as.vector(d_hi)
  )
})

# 12) Combine and export
proj_age <- bind_rows(proj_list) %>% arrange(disease, sex, age, year)
write_xlsx(proj_age, "~/Desktop/deaths2040.xlsx")
