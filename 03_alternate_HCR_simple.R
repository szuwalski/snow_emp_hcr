
# ---- Libraries --------------------------------------------------------------
library(crabpack)                        # survey pull + bioabundance
library(dplyr); library(tidyr); library(reshape2)   # data wrangling
library(ggplot2); library(ggridges); library(patchwork)  # figures
library(png); library(grid)              # raster/image helpers

# ---- Helpers ----------------------------------------------------------------
# annotation_custom2(): place a grob (e.g. a background image) on a ggplot facet.
# NOTE: currently UNUSED in this script (kept from a prior version that overlaid
#       a snow-crab image; see the commented readPNG below). Safe to delete if
#       the background image is not coming back.
annotation_custom2 <- function(grob, xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, data) {
  layer(data = data,
        stat = StatIdentity,
        position = PositionIdentity,
        geom = ggplot2:::GeomCustomAnn,
        inherit.aes = TRUE,
        params = list(grob = grob, xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax))
}
# in_png_opie <- readPNG('continuity/snowcrab.png')   # (background image; unused)


# =============================================================================
# 1. PULL SURVEY SPECIMEN DATA  (crabpack API)
# =============================================================================
# years = 1982:2026  ->  terminal year = the 2026 summer survey (advance each
# cycle). Requires API access; everything below derives from `specimen_data`.
specimen_data <- crabpack::get_specimen_data(species = "SNOW",
                                             region  = "EBS",
                                             years   = c(1982:2026),
                                             channel = 'API')


# =============================================================================
# 2. MALE ABUNDANCE / BIOMASS INDICES
# =============================================================================
# Legal males (used only as a rough MMB-CV proxy, below) and large + preferred
# males (the >101 mm index written out for the report/model).
male_snow_ind <- crabpack::calc_bioabund(crab_data = specimen_data,
                                         species = "SNOW",
                                         region  = "EBS",
                                         crab_category = c("legal_male"))

big_male_snow_ind <- crabpack::calc_bioabund(crab_data = specimen_data,
                                             species = "SNOW",
                                             region  = "EBS",
                                             crab_category = c("large_male", "preferred_male"))

male_snow <- crabpack::calc_bioabund(crab_data = specimen_data,
                                     species = "SNOW",
                                     region  = "EBS",
                                     sex = 'male',
                                     shell_condition = "all_categories",
                                     size_min = 25,
                                     bin_1mm  = TRUE)

tmp<-filter(specimen_data$specimen,SIZE_1MM>101)
in_tmp<-merge(tmp,specimen_data$haul)
temp_occupied_lg<-in_tmp%>%
  group_by(YEAR)%>%
  summarize(temp_occ=weighted.mean(GEAR_TEMPERATURE ,SAMPLING_FACTOR ,na.rm=T))

ggplot(temp_occupied_lg)+
  geom_line(aes(x=YEAR,y=temp_occ))

new_dat <- seq(27.5, 132.5, 5)                       

new_male_mat_dat_in<-get_male_maturity(
  species = "SNOW",
  region = c("EBS", "NBS")[1])


new_male_mat_dat <- read.csv("data/snow_ogives.csv") %>%
  filter(SIZE_5MM %in% new_dat) %>%
  select(YEAR, SIZE_5MM, PROP_MATURE) %>%
  tidyr::pivot_wider(names_from = SIZE_5MM, values_from = PROP_MATURE) %>%
  arrange(YEAR) %>%
  as.data.frame()

tmp<-melt(new_male_mat_dat,id.vars=c("YEAR"))
colnames(tmp)<-c("Year","Size","p_mat")
ggplot(tmp)+
  geom_line(aes(x=Size,y=p_mat,col=as.factor(Year),group=Year),lwd=1.3,alpha=.8)+
  theme_bw()


# new_male_mat_dat <- new_male_mat_dat_in$male_mat_ratio%>%
#   filter(SIZE_BIN  %in% new_dat) %>%
#   select(YEAR, SIZE_BIN, PROP_MATURE) %>%
#   tidyr::pivot_wider(names_from = SIZE_BIN, values_from = PROP_MATURE) %>%
#   arrange(YEAR) %>%
#   as.data.frame()
# 
# new_male_mat_dat <- new_male_mat_dat[, c("YEAR", as.character(new_dat))]
# fill model bins outside each year's observed size range (0 below, 1 above)
for (i in seq_len(nrow(new_male_mat_dat))) {
  v   <- as.numeric(new_male_mat_dat[i, -1])
  obs <- new_dat[!is.na(v)]
  if (length(obs)) {
    v[new_dat < min(obs) & is.na(v)] <- 0
    v[new_dat > max(obs) & is.na(v)] <- 1
    new_male_mat_dat[i, -1] <- v
  }
}

# build a (year x size-bin) maturity matrix, filling missing years with the mean
all_yr <- seq(1982, max(new_male_mat_dat[, 1]))
allmat <- matrix(nrow = length(all_yr), ncol = length(new_dat))
rownames(allmat) <- all_yr
colnames(allmat) <- new_dat
yrs <- unique(new_male_mat_dat[, 1])
for (x in 1:length(yrs))
  allmat[match(yrs[x], all_yr), ] <- unlist(new_male_mat_dat[x, 2:ncol(new_male_mat_dat)])

mean_mat <- apply(allmat, 2, mean, na.rm = T)
for (x in which(is.na(allmat[, 1])))
  allmat[x, ] <- mean_mat

male_maturity_ogive <- data.frame(year = as.numeric(rownames(allmat)), allmat,
                                  check.names = FALSE)
names(male_maturity_ogive) <- c("year", paste0("m", new_dat))
# =============================================================================
# 4. MALE SIZE COMPOSITIONS BY SHELL CONDITION  (5-mm bins, 27.5-132.5 + plus grp)
# =============================================================================
# Re-bin the 1-mm male numbers to 5-mm bins separately for new-shell
# (new_hardshell + soft_molting) and old-shell (oldshell + very_oldshell), then
# reshape wide, constrain to the model's 27.5-132.5 bins with a 132.5+ plus
# group, and convert to millions of crab -> MaleNew, MaleOld.
# NB: right = FALSE here (the "previous survey approach"): a crab on a 5-mm cutoff
#     goes to the UPPER bin. Consistent with the fishery comps (2026-08, per Grant).

# ---- new-shell males --------------------------------------------------------
new_male_snow <- filter(male_snow, SHELL_TEXT %in% c('new_hardshell', 'soft_molting')) %>%
  group_by(YEAR, group = cut(SIZE_1MM, breaks = seq(0, max(SIZE_1MM), 5), right = FALSE)) %>%
  summarise(n = sum(ABUNDANCE),
            data = 'crabpack')

# derive the numeric bin midpoint from the "[lo,hi)" factor label
tmp <- strsplit(as.character(new_male_snow$group), ',')
new_male_snow$mid_pts <- NA
for (x in 1:length(tmp)) {
  in1 <- as.numeric(substr(unlist(tmp[x])[1], 2, nchar(unlist(tmp[x])[1])))
  in2 <- as.numeric(substr(unlist(tmp[x])[2], 1, nchar(unlist(tmp[x])[2]) - 1))
  new_male_snow$mid_pts[x] <- sum(in1, in2) / 2
}

# ---- old-shell males --------------------------------------------------------
old_male_snow <- filter(male_snow, SHELL_TEXT %in% c('oldshell', 'very_oldshell')) %>%
  group_by(YEAR, group = cut(SIZE_1MM, breaks = seq(0, max(SIZE_1MM), 5), right = FALSE)) %>%
  summarise(n = sum(ABUNDANCE),
            data = 'crabpack')

tmp <- strsplit(as.character(old_male_snow$group), ',')
old_male_snow$mid_pts <- NA
for (x in 1:length(tmp)) {
  in1 <- as.numeric(substr(unlist(tmp[x])[1], 2, nchar(unlist(tmp[x])[1])))
  in2 <- as.numeric(substr(unlist(tmp[x])[2], 1, nchar(unlist(tmp[x])[2]) - 1))
  old_male_snow$mid_pts[x] <- sum(in1, in2) / 2
}

# ---- reshape to year x size-bin matrices ------------------------------------
in_old <- old_male_snow[, c(1, 3, 5)]
old_male_wide <- dcast(in_old, YEAR ~ mid_pts, value.var = c("n"))

in_new <- new_male_snow[, c(1, 3, 5)]
new_male_wide <- dcast(in_new, YEAR ~ mid_pts, value.var = c("n"))

# ---- constrain to the model size bins + 132.5 plus group --------------------
# (Cody's note: this bin-constraining is awkward and should eventually move into
#  the data pull.) Everything above 132.5 mm is summed into the top bin.
sizes <- seq(27.5, 132.5, 5)

use_old_male <- old_male_wide[, which(!is.na(match(as.numeric(colnames(old_male_wide)), sizes)))]
rownames(use_old_male) <- (old_male_wide[, 1])
use_old_male[, ncol(use_old_male)] <- use_old_male[, ncol(use_old_male)] +
  apply(old_male_wide[, which(as.numeric(colnames(old_male_wide)) > 132.5)], 1, sum)

use_new_male <- new_male_wide[, which(!is.na(match(as.numeric(colnames(new_male_wide)), sizes)))]
rownames(use_new_male) <- (new_male_wide[, 1])
use_new_male[, ncol(use_new_male)] <- use_new_male[, ncol(use_new_male)] +
  apply(new_male_wide[, which(as.numeric(colnames(new_male_wide)) > 132.5)], 1, sum)

MaleOld <- use_old_male / 1000000   # -> millions of crab
MaleNew <- use_new_male / 1000000

# =============================================================================
# 6. MATURE / IMMATURE SIZE COMPOSITIONS + BIOMASS INDICES
# =============================================================================
# Split new-shell males into mature/immature via the ogive (dropping 2020, the
# cancelled-survey COVID year, so `in_mat` aligns with MaleNew's years), add all
# old-shell males to the mature pool, then normalize to size comps and weight by
# weight-at-size to get biomass indices.
in_mat        <- allmat[-which(rownames(allmat) == 2020), ]
male_new_mature <- MaleNew * in_mat
male_immature <- MaleNew * (1 - in_mat)
male_mature   <- male_new_mature + MaleOld

# normalized MALE size compositions (rows sum to 1)
sc_male_mat <- sweep(male_mature,   1, apply(male_mature,   1, sum, na.rm = T), FUN = "/")
sc_male_imm <- sweep(male_immature, 1, apply(male_immature, 1, sum, na.rm = T), FUN = "/")

# ---- biomass indices (numbers-at-size x weight-at-size) ---------------------
wt_at_size   <- read.csv("data/wt_at_size.csv")
male_mat_bio <- apply(sweep(male_mature,   2, wt_at_size[, 1], FUN = "*"), 1, sum)
male_imm_bio <- apply(sweep(male_immature, 2, wt_at_size[, 1], FUN = "*"), 1, sum)

# quick interactive sanity plot (mature vs immature biomass); not saved.
plot(male_mat_bio, type = 'b', ylim = c(0, 400))
lines(male_imm_bio, type = 'b', col = 2)



#========================================
# calculate the incoming immature biomass given assumed survival
size_trans<-as.matrix(read.csv("data/size_transition.csv",header=F))
proj_m<-0.3

#==projected under no fishing
male_immature_incoming<-as.matrix(male_immature)
for(x in 1:nrow(male_immature_incoming))
  male_immature_incoming[x,]<-exp(-proj_m)*male_immature_incoming[x,]%*%size_trans
male_mature_remaining<-as.matrix(male_mature)*exp(-proj_m)
all_male_proj_nofish<-male_immature_incoming + male_mature_remaining
all_male_obs<-male_immature + male_mature

proj_long<-melt(all_male_proj_nofish)
colnames(proj_long)<-c("Year","Size","N")
proj_long$type<-"proj"
proj_long$Year<-proj_long$Year+1

#==projected under observed fishing
ret_cat<-read.csv("data/retained_n_at_size.csv")[-1,]
disc_cat<-read.csv("data/discard_n_at_size.csv")[-1,]
div_n<-1000000
disc_mort<-0.3

#==just use from 1991 (remove first year)
ret_cat[,2:ncol(ret_cat)]<-ret_cat[,2:ncol(ret_cat)]/div_n
disc_cat[,2:ncol(disc_cat)]<-disc_mort*disc_cat[,2:ncol(disc_cat)]/div_n

tot_mort<-ret_cat
tot_mort[,2:ncol(tot_mort)]<-tot_mort[,2:ncol(tot_mort)] + disc_cat[,2:ncol(disc_cat)]

#==remove catches
male_mature_cat_removed<-male_mature
male_immature_cat_removed<-male_immature

for(x in 1:nrow(tot_mort))
{
 ind<-which(!is.na(match(rownames(male_mature_cat_removed),tot_mort[x,1])))
 prop_imm<-male_immature_cat_removed[ind,]/(male_immature_cat_removed[ind,] + male_mature_cat_removed[ind,])
 
 if(length(ind)>0)
 {
  male_mature_cat_removed[ind,]<- male_mature_cat_removed[ind,] - ((1-prop_imm)*tot_mort[x,-1])
  male_immature_cat_removed[ind,]<- male_immature_cat_removed[ind,] - (prop_imm*tot_mort[x,-1])
 }
}

#==then growth happens
male_immature_incoming_cat<-as.matrix(male_immature_cat_removed)
for(x in 1:nrow(male_immature_incoming_cat))
  male_immature_incoming_cat[x,]<-exp(-proj_m)*male_immature_incoming_cat[x,]%*%size_trans
male_mature_remaining_cat<-as.matrix(male_mature_cat_removed)*exp(-proj_m)
all_male_proj_fish<-male_immature_incoming_cat + male_mature_remaining_cat

cat_long<-melt(all_male_proj_fish)
colnames(cat_long)<-c("Year","Size","N")
cat_long$type<-"proj + catch"
cat_long$Year<-cat_long$Year+1

library(tibble)
df_with_rownames <- rownames_to_column(all_male_obs, var = "RowName")
obs_long <- melt(df_with_rownames, id.vars = "RowName")
colnames(obs_long)<-c("Year","Size","N")
obs_long$type<-"obs"
inny<-rbind(cat_long,proj_long,obs_long)

#========================
# estimate yearly survival
library(RTMB)

# 1. PRE-SLICE YOUR DATA LIST OUTSIDE THE FUNCTION
# Apply all row and column modifications here so the shapes match perfectly

library(RTMB)

# 1. PRE-SLICE YOUR DATA LIST OUTSIDE THE FUNCTION
# Apply all row and column modifications here so the shapes match perfectly

pars<-list(
  proj_m=rep(log(0.3), nrow(fit_obs)),
  log_m_mu= -1.2,
  mu_sd= log(0.6)
)


dat <- list(
  size_trans = size_trans,
  
  # Pre-slice the input matrices to drop the last row so they align with fit_obs
  male_immature_cat_removed = as.matrix(male_immature_cat_removed[which(rownames(male_immature_cat_removed) > 1990), ]),
  male_mature_cat_removed   = as.matrix(male_mature_cat_removed[which(rownames(male_immature_cat_removed) > 1990), ]),
  
  # Pre-slice fit_obs to only include columns 16 onwards
  fit_obs = fit_obs[-1, 16:ncol(fit_obs)] 
)

# 2. DEFINE THE CLEAN OBJECTIVE FUNCTION
small_pop_dy <- function(pars) {
  getAll(pars, dat)
  
  # Wrap the pre-sliced observation matrix cleanly
  fit_obs <- OBS(fit_obs) 
  
  # Initialize tracking matrices
  male_immature_incoming_cat <- male_immature_cat_removed
  male_mature_remaining_cat  <- male_mature_cat_removed
  
  # Population projection loop
  for(x in 1:nrow(male_immature_incoming_cat)) 
  {
    in_m<-exp(proj_m[x])
    male_immature_incoming_cat[x, ] <- exp(-in_m) * male_immature_cat_removed[x, ] %*% size_trans
    male_mature_remaining_cat[x, ]  <- male_mature_cat_removed[x, ] * exp(-in_m)
  }
  
  all_male_proj_fish <- male_immature_incoming_cat + male_mature_remaining_cat
  
  nll <- 0
  
  # A. Random Effect Component
  nll <- nll - sum(dnorm(proj_m, mean = log_m_mu, sd = exp(mu_sd), log = TRUE))
  
  # B. Numbers at Size Component
  # Slicing the prediction matrix: drop the last row, and take columns 15+ to match fit_obs columns (16:ncol)
  # (Note: since fit_obs was sliced from column 16 of the original matrix, 
  # pred_matrix must take the matching columns 16 to ncol)
  pred_matrix <- all_male_proj_fish[1:(nrow(all_male_proj_fish) - 1), 16:ncol(all_male_proj_fish)]
  
  # Force both sides explicitly back to a tracked matrix structure to be safe
  log_res <- (as.matrix(fit_obs) + 1e-6) - (as.matrix(pred_matrix) + 1e-6)
  
  nll <- nll - sum(dnorm(log_res, mean = 0, sd = 0.1, log = TRUE))
  REPORT(all_male_proj_fish)
  return(nll)
}

# 3. BUILD THE OBJECTIVE FUNCTION
obj <- MakeADFun(
  func = small_pop_dy, 
  parameters = pars, 
  random = "proj_m"
)

# 4. ESTIMATE PARAMETERS
opt <- nlminb(obj$par, obj$fn, obj$gr)
print(opt)

library(ggplot2)

# 1. Extract summaries for both random effects and fixed effects
rep <- sdreport(obj)
rep_rand  <- summary(rep, select = "random")
rep_fixed <- summary(rep, select = "fixed")

# 2. Extract the estimated log-mean and calculate its normal-scale value
# (We grab the value from the row labeled "log_m_mu")
log_M_mean_est <- rep_fixed["log_m_mu", "Estimate"]
M_mean_normal  <- exp(log_M_mean_est)

# 3. Build the plotting dataset for the yearly random effects
years <- as.numeric(rownames(dat$male_immature_cat_removed))

plot_data <- data.frame(
  Year  = years,
  M_log = rep_rand[, "Estimate"],
  SE    = rep_rand[, "Std. Error"]
)

# 4. Transform log random effects to strictly positive normal space
plot_data$M_est <- exp(plot_data$M_log)
plot_data$lower <- exp(plot_data$M_log - 1.96 * plot_data$SE)
plot_data$upper <- exp(plot_data$M_log + 1.96 * plot_data$SE)

# 5. Generate the ggplot with the average horizontal reference line
est_m_plot<-ggplot(plot_data, aes(x = Year)) +
  # Shaded uncertainty band for yearly variations
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = "skyblue", alpha = 0.4) +
  
  # --- ADDED: Horizontal line representing the overall average M ---
  geom_hline(aes(yintercept = M_mean_normal), color = "darkred", 
             linetype = "dashed", linewidth = 1, show.legend = TRUE) +
  
  # Yearly estimated point lines
  geom_line(aes(y = M_est), color = "blue", linewidth = 1) +
  geom_point(aes(y = M_est), color = "blue", size = 2) +
  
  # Add text label for the average line to make it easily readable
  annotate("text", x = min(plot_data$Year) + 2, y = M_mean_normal, 
           label = paste("Average M =", round(M_mean_normal, 3)), 
           vjust = -0.7, color = "darkred", fontface = "bold", size = 4) +
  
  # Labels and styling
  labs(
    x = "Year",
    y = "Natural Mortality"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    panel.grid.minor = element_blank()
  )+
  ylim(0,2.25)


est_numbers_at_size <- rep_values$all_male_proj_fish
rownames(est_numbers_at_size) <- rownames(dat$male_immature_cat_removed)
colnames(est_numbers_at_size) <- colnames(dat$male_immature_cat_removed)





#====================================

mod_long <- melt(est_numbers_at_size, id.vars = "RowName")
colnames(mod_long)<-c("Year","Size","N")
mod_long$type<-"proj + catch + est_surv"
mod_long$Year<-mod_long$Year+1
inny<-rbind(cat_long,proj_long,obs_long,mod_long)

ggplot(filter(inny,as.numeric(Size)>100 & Year>1991))+
  geom_line(aes(x=as.numeric(Size),y=N,lty=type,col=type))+
  facet_wrap(~Year)+theme_bw()
  
ggplot(filter(inny,as.numeric(Size)>100 & Year>2013))+
  geom_line(aes(x=as.numeric(Size),y=N,lty=type,col=type))+
  facet_wrap(~Year)+theme_bw()

big_bois<-filter(inny,as.numeric(Size)>100)

merr<-big_bois%>%
  group_by(type,Year)%>%
  summarize(tot_n=sum(N))

ggplot(merr)+
  geom_line(aes(x=as.numeric(Year),y=tot_n,col=type))+
  geom_point(aes(x=as.numeric(Year),y=tot_n,col=type))+
  theme_bw()

#======================================
# GAMS for temp and density on survival
#======================================
colnames(temp_occupied_lg)[1]<-"Year"
uhm<-merge(filter(merr,type=='obs')[,-1],temp_occupied_lg)
in_dat<-merge(uhm,plot_data)
in_dat$mortality<- 1-exp(-in_dat$M_est)

library(mgcv)
mod <- gam(mortality ~ s(tot_n,k=3) + s(temp_occ,k=3),
           family = betar(link = "logit"), 
           method = "REML", 
           data = in_dat)

# mod <- gam(mortality ~ s(tot_n,temp_occ,k=12),
#            family = betar(link = "logit"), 
#            method = "REML", 
#            data = in_dat)

predict(mod,newdata=data.frame(tot_n=100,temp_occ=3),type='response')

summary(mod)
plot(mod,pages=1,too.far=1,scheme=2)

#==estimate time-varying proj_m for immature + mature, but only include size range that survey selectivity
#===could be assumed to be 1 AND represents one molt away from commercial?


#===========================================================================
#==I need to keep shell condition in here so I can estimate mature mortality
#==without confounding with immature mortality
