########################
# Load data sets
########################

library(dplyr)
library(readr)
# Load environmental data; starting with species richness to ensure as many pixels are recovered as possible
env <- list.files(path="./environmental_data", pattern = "richness", full.names = TRUE) %>% lapply(read_csv) %>% bind_rows 
env <- data.frame(env)
env[env == -9999] <- NA
env$x <- round(env$x, digit = 1)
env$y <- round(env$y, digit = 1)
env %>% group_by(x, y) %>% summarize_if(is.numeric, mean, na.rm = TRUE) -> env
env <- as.data.frame(env)

# Add RPD randomizations
rand_RPD <- read.csv("./Mimosoid_CSVs_ToShare_WGS84/rand_RPD_20km.csv")
rand_RPD$x <- round(rand_RPD$x, digit = 1)
rand_RPD$y <- round(rand_RPD$y, digit = 1)
rand_RPD %>% group_by(x, y) %>% summarize_if(is.numeric, mean, na.rm = TRUE) -> rand_RPD
rand_RPD <- as.data.frame(rand_RPD)
# Add significance column
rand_RPD$RPD_significance <- as.factor(ifelse(rand_RPD$value < 0.05, "Low", ifelse(rand_RPD$value > 0.95, "High", "NS")))
rand_RPD$value <- NULL

combined <- merge(rand_RPD, env, by = c("x", "y"))

# Add RPD
RPD <- read.csv("./Mimosoid_CSVs_ToShare_WGS84/RPD_20km.csv")
names(RPD) <- c("x", "y", "RPD")
RPD$x <- round(RPD$x, digit = 1)
RPD$y <- round(RPD$y, digit = 1)
RPD %>% group_by(x, y) %>% summarize_if(is.numeric, mean, na.rm = TRUE) -> RPD
RPD <- as.data.frame(RPD)

combined <- merge(combined, RPD, by = c("x", "y"))

# Add CANAPE
CANAPE <- read.csv("./Mimosoid_CSVs_ToShare_WGS84/CANAPE_20km.csv")
names(CANAPE) <- c("x", "y", "CANAPE")
CANAPE$x <- round(CANAPE$x, digit = 1)
CANAPE$y <- round(CANAPE$y, digit = 1)
CANAPE %>% group_by(x, y) %>% summarize_if(is.character, max) -> CANAPE
CANAPE <- as.data.frame(CANAPE)
CANAPE$CANAPE <- as.factor(CANAPE$CANAPE)

combined <- merge(combined, CANAPE, by = c("x", "y"))

## Add SR
SR <- read.csv("./Mimosoid_CSVs_ToShare_WGS84/richness_20km.csv")
names(SR) <- c("x", "y", "SR")
SR$x <- round(SR$x, digit = 1)
SR$y <- round(SR$y, digit = 1)
SR %>% group_by(x, y) %>% summarize_if(is.numeric, mean, na.rm = TRUE) -> SR
SR <- as.data.frame(SR)

combined <- merge(combined, SR, by = c("x", "y"))

# Normalize entire data frame
combined.temp <- combined
combined.scaled <- rapply(combined.temp, scale, c("numeric","integer"), how="replace")
combined.scaled$SR <- combined$SR
combined.scaled <- as.data.frame(combined.scaled)
combined.scaled$y <- combined.temp$y
combined.scaled$x <- combined.temp$x


########################
# Model
########################


# RPD model
linear_model_complex <- lm(RPD ~ aridity_index_UNEP + BIOCLIM_1 + BIOCLIM_12 + BIOCLIM_7 + BIOCLIM_17 + ISRICSOILGRIDS_new_average_nitrogen_reduced + ISRICSOILGRIDS_new_average_phx10percent_reduced + ISRICSOILGRIDS_new_average_soilorganiccarboncontent_reduced, data = combined.scaled)
# Top 5 predictors by GLM normalized coefficient... note that nitrogen is a poorly performing predictor
linear_model_simple <- lm(RPD ~ aridity_index_UNEP + BIOCLIM_1 + BIOCLIM_12 + BIOCLIM_17 + ISRICSOILGRIDS_new_average_soilorganiccarboncontent_reduced, data = combined.scaled)
library(lme4)
mixed_model_complex <- lmer(RPD ~ aridity_index_UNEP + BIOCLIM_1 + BIOCLIM_12 + BIOCLIM_7 + BIOCLIM_17 + ISRICSOILGRIDS_new_average_nitrogen_reduced + ISRICSOILGRIDS_new_average_phx10percent_reduced + ISRICSOILGRIDS_new_average_soilorganiccarboncontent_reduced + (1 | y) + (1 | x), na.action = na.omit, data = combined.scaled)
# Top 5 predictors by LMM normalized coefficient... note that nitrogen is a poorly performing predictor
mixed_model_simple <- lmer(RPD ~ aridity_index_UNEP + BIOCLIM_12 + BIOCLIM_7 + ISRICSOILGRIDS_new_average_phx10percent_reduced + ISRICSOILGRIDS_new_average_soilorganiccarboncontent_reduced + (1 | y) + (1 | x), na.action = na.omit, data = combined.scaled)
mixed_model_noenvironment <- lmer(RPD ~ (1 | y) + (1 | x), na.action = na.omit, data = combined.scaled)

AIC(linear_model_complex)
AIC(linear_model_simple)
AIC(mixed_model_complex)
AIC(mixed_model_simple)
AIC(mixed_model_noenvironment)
# Complex mixed model favored

summary(mixed_model_complex)

library(MuMIn)
r.squaredGLMM(mixed_model_complex)


## RPD significance model -- if included need to redo the simple model variable set
#levels(combined.scaled$RPD_significance)[levels(combined.scaled$RPD_significance)=="NS"] <- NA # NS excluded by the metrics
## Top 5 variables in simple models chosen based on coefficients in complex logit model
#logit_complex <- glm(RPD_significance ~ aridity_index_UNEP + BIOCLIM_1 + BIOCLIM_12 + BIOCLIM_7 + BIOCLIM_17 + ISRICSOILGRIDS_new_average_nitrogen_reduced + ISRICSOILGRIDS_new_average_phx10percent_reduced + ISRICSOILGRIDS_new_average_soilorganiccarboncontent_reduced, family = binomial(link='logit'), data = combined.scaled)
#logit_simple <- glm(RPD_significance ~ aridity_index_UNEP + BIOCLIM_1 + BIOCLIM_17 + ISRICSOILGRIDS_new_average_nitrogen_reduced + ISRICSOILGRIDS_new_average_phx10percent_reduced, family = binomial(link='logit'), data = combined.scaled)
## Be warned, these logit models will take about 5 minutes and 2.8 GB RAM.
#mixed_model_complex <- glmer(RPD_significance ~ aridity_index_UNEP + BIOCLIM_1 + BIOCLIM_12 + BIOCLIM_7 + BIOCLIM_17 + ISRICSOILGRIDS_new_average_nitrogen_reduced + ISRICSOILGRIDS_new_average_phx10percent_reduced + ISRICSOILGRIDS_new_average_soilorganiccarboncontent_reduced + (1 | y) + (1 | x), family=binomial(link='logit'), na.action = na.omit, data = combined.scaled, control = glmerControl(optimizer = "bobyqa",optCtrl = list(maxfun = 2e5))) # Stronger likelihood search options per warnings + documentation
#mixed_model_simple <- glmer(RPD_significance ~ aridity_index_UNEP + BIOCLIM_1 + BIOCLIM_17 + ISRICSOILGRIDS_new_average_nitrogen_reduced + ISRICSOILGRIDS_new_average_phx10percent_reduced + (1 | y) + (1 | x), family=binomial(link='logit'), na.action = na.omit, data = combined.scaled, control = glmerControl(optimizer = "bobyqa",optCtrl = list(maxfun = 2e5))) # Stronger likelihood search options per warnings + documentation
#mixed_model_noenvironment <- glmer(RPD_significance ~ (1 | y) + (1 | x), family=binomial(link='logit'), na.action = na.omit, data = combined.scaled, control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun = 2e5)))
#
#AIC(logit_complex)
#AIC(logit_simple)
#AIC(mixed_model_complex)
#AIC(mixed_model_simple)
#AIC(mixed_model_noenvironment)
#
## Mixed model complex favored
#summary(mixed_model_complex)
#
#library(MuMIn)
#r.squaredGLMM(mixed_model_complex)



# CANAPE significance model
combined.scaled$CANAPE_significant <- combined.scaled$CANAPE
levels(combined.scaled$CANAPE_significant)[levels(combined.scaled$CANAPE_significant)=="Neo"] <-"Sig"
levels(combined.scaled$CANAPE_significant)[levels(combined.scaled$CANAPE_significant)=="Mixed"] <-"Sig"
levels(combined.scaled$CANAPE_significant)[levels(combined.scaled$CANAPE_significant)=="Paleo"] <-"Sig"
levels(combined.scaled$CANAPE_significant)
# Reduce factor to 2 levels; otherwise glm guesses and does this silently
# Top 5 variables in simple models chosen based on coefficients in complex logit model
logit_complex <- glm(CANAPE_significant ~ aridity_index_UNEP + BIOCLIM_1 + BIOCLIM_12 + BIOCLIM_7 + BIOCLIM_17 + ISRICSOILGRIDS_new_average_nitrogen_reduced + ISRICSOILGRIDS_new_average_phx10percent_reduced + ISRICSOILGRIDS_new_average_soilorganiccarboncontent_reduced, family = binomial(link='logit'), data = combined.scaled)
logit_simple <- glm(CANAPE_significant ~ aridity_index_UNEP + BIOCLIM_1 + BIOCLIM_17 + ISRICSOILGRIDS_new_average_nitrogen_reduced + ISRICSOILGRIDS_new_average_phx10percent_reduced, family = binomial(link='logit'), data = combined.scaled)
# Be warned, these logit models will take about 5 minutes and 2.8 GB RAM.
mixed_model_complex <- glmer(CANAPE_significant ~ aridity_index_UNEP + BIOCLIM_1 + BIOCLIM_12 + BIOCLIM_7 + BIOCLIM_17 + ISRICSOILGRIDS_new_average_nitrogen_reduced + ISRICSOILGRIDS_new_average_phx10percent_reduced + ISRICSOILGRIDS_new_average_soilorganiccarboncontent_reduced + (1 | y) + (1 | x), family=binomial(link='logit'), na.action = na.omit, data = combined.scaled, control = glmerControl(optimizer = "bobyqa",optCtrl = list(maxfun = 2e5))) # Stronger likelihood search options per warnings + documentation
mixed_model_simple <- glmer(CANAPE_significant ~ aridity_index_UNEP + BIOCLIM_1 + BIOCLIM_17 + ISRICSOILGRIDS_new_average_nitrogen_reduced + ISRICSOILGRIDS_new_average_phx10percent_reduced + (1 | y) + (1 | x), family=binomial(link='logit'), na.action = na.omit, data = combined.scaled, control = glmerControl(optimizer = "bobyqa",optCtrl = list(maxfun = 2e5))) # Stronger likelihood search options per warnings + documentation
mixed_model_noenvironment <- glmer(CANAPE_significant ~ (1 | y) + (1 | x), family=binomial(link='logit'), na.action = na.omit, data = combined.scaled, control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun = 2e5)))

AIC(logit_complex)
AIC(logit_simple)
AIC(mixed_model_complex)
AIC(mixed_model_simple)
AIC(mixed_model_noenvironment)

# Mixed model simple favored
summary(mixed_model_simple)

library(MuMIn)
r.squaredGLMM(mixed_model_simple)



########################
# Some plots
########################

library(ggplot2)
levels(combined.scaled$RPD_significance)[levels(combined.scaled$RPD_significance)=="NS"] <- NA # NS excluded by the metrics
levels(combined$RPD_significance)[levels(combined$RPD_significance)=="NS"] <- NA # NS excluded by the metrics

# Violin plots
ggplot(combined[!is.na(combined$RPD_significance), ], aes(x = RPD_significance, y = aridity_index_UNEP, fill = RPD_significance)) + geom_violin(trim=FALSE) + geom_boxplot(width=0.1, fill="white") + labs(title="Aridity vs. RPD significance", x="RPD significance", y = "UNEP aridity index") + ylim(quantile(combined$aridity_index_UNEP, 0.025, na.rm = TRUE), quantile(combined$aridity_index_UNEP, 0.975, na.rm = TRUE))
ggplot(combined[!is.na(combined$RPD_significance), ], aes(x = RPD_significance, y = BIOCLIM_1/10, fill = RPD_significance)) + geom_violin(trim=FALSE) + geom_boxplot(width=0.1, fill="white") + labs(title="Temperature vs. RPD significance", x="RPD significance", y = "Mean annual temperature (°C)") + ylim(quantile(combined$BIOCLIM_1/10, 0.025, na.rm = TRUE), quantile(combined$BIOCLIM_1/10, 0.975, na.rm = TRUE))
ggplot(combined[!is.na(combined$RPD_significance), ], aes(x = RPD_significance, y = ISRICSOILGRIDS_new_average_nitrogen_reduced, fill = RPD_significance)) + geom_violin(trim=FALSE) + geom_boxplot(width=0.1, fill="white") + labs(title="Nitrogen vs. RPD significance", x="RPD significance", y = "Nitrogen content") + ylim(quantile(combined$ISRICSOILGRIDS_new_average_nitrogen_reduced, 0.025, na.rm = TRUE), quantile(combined$ISRICSOILGRIDS_new_average_nitrogen_reduced, 0.975, na.rm = TRUE))
ggplot(combined[!is.na(combined$RPD_significance), ], aes(x = RPD_significance, y = ISRICSOILGRIDS_new_average_phx10percent_reduced/10, fill = RPD_significance)) + geom_violin(trim=FALSE) + geom_boxplot(width=0.1, fill="white") + labs(title="pH vs. RPD significance", x="RPD significance", y = "Soil pH") + ylim(quantile(combined$ISRICSOILGRIDS_new_average_phx10percent_reduced/10, 0.025, na.rm = TRUE), quantile(combined$ISRICSOILGRIDS_new_average_phx10percent_reduced/10, 0.975, na.rm = TRUE))
# Only BIO1 shows clear separation

ggplot(combined, aes(x = CANAPE, y = aridity_index_UNEP, fill = CANAPE)) + geom_violin(trim=FALSE) + geom_boxplot(width=0.1, fill="white") + labs(title="Aridity vs. CANAPE significance", x="CANAPE significance", y = "UNEP aridity index") + ylim(quantile(combined$aridity_index_UNEP, 0.025, na.rm = TRUE), quantile(combined$aridity_index_UNEP, 0.975, na.rm = TRUE))
ggplot(combined, aes(x = CANAPE, y = BIOCLIM_1/10, fill = CANAPE)) + geom_violin(trim=FALSE) + geom_boxplot(width=0.1, fill="white") + labs(title="Temperature vs. CANAPE significance", x="CANAPE significance", y = "Mean annual temperature (°C)") + ylim(quantile(combined$BIOCLIM_1/10, 0.025, na.rm = TRUE), quantile(combined$BIOCLIM_1/10, 0.975, na.rm = TRUE))
ggplot(combined, aes(x = CANAPE, y = ISRICSOILGRIDS_new_average_nitrogen_reduced, fill = CANAPE)) + geom_violin(trim=FALSE) + geom_boxplot(width=0.1, fill="white") + labs(title="Nitrogen vs. CANAPE significance", x="CANAPE significance", y = "Nitrogen content") + ylim(quantile(combined$ISRICSOILGRIDS_new_average_nitrogen_reduced, 0.025, na.rm = TRUE), quantile(combined$ISRICSOILGRIDS_new_average_nitrogen_reduced, 0.975, na.rm = TRUE))
ggplot(combined, aes(x = CANAPE, y = ISRICSOILGRIDS_new_average_phx10percent_reduced/10, fill = CANAPE)) + geom_violin(trim=FALSE) + geom_boxplot(width=0.1, fill="white") + labs(title="pH vs. CANAPE significance", x="CANAPE significance", y = "Soil pH") + ylim(quantile(combined$ISRICSOILGRIDS_new_average_phx10percent_reduced/10, 0.025, na.rm = TRUE), quantile(combined$ISRICSOILGRIDS_new_average_phx10percent_reduced/10, 0.975, na.rm = TRUE))


# Hex plots
# Species richness response
ggplot(combined, aes(x = aridity_index_UNEP, y = SR)) + geom_hex(bins = 35) + scale_fill_continuous(type = "viridis") + theme_bw() + geom_smooth(method='lm', formula = y ~ x) + labs(title="Aridity vs. species richness", x="Aridity index", y = "SR")
ggplot(combined, aes(x = BIOCLIM_1/10, y = SR)) + geom_hex(bins = 35) + scale_fill_continuous(type = "viridis") + theme_bw() + geom_smooth(method='lm', formula = y ~ x) + labs(title="Mean annual temperature (°C) vs. species richness", x="Aridity index", y = "SR")

# Latitude
ggplot(combined, aes(x = y, y = SR)) + geom_hex(bins = 30) + scale_fill_continuous(type = "viridis") + theme_bw() + geom_smooth(method='lm', formula = y ~ x) + labs(title="Species richness vs. latitude", y="Species richness", x = "Latitude")
ggplot(combined, aes(x = y, y = RPD)) + geom_hex(bins = 30) + scale_fill_continuous(type = "viridis") + theme_bw() + geom_smooth(method='lm', formula = y ~ x) + labs(title="RPD vs. latitude", y="RPD", x = "Latitude")


# North South Hemisphere RPD
# Note stronger mimosoid RPD separation in the northern hemisphere
ggplot(combined[combined$y > 0, ], aes(x = RPD_significance, y = aridity_index_UNEP, fill = RPD_significance)) + geom_violin(trim=FALSE) + geom_boxplot(width=0.1, fill="white") + labs(title="Aridity vs. RPD significance, Northern Hemisphere", x="RPD_significance", y = "UNEP aridity index") + ylim(quantile(combined$aridity_index_UNEP, 0.025, na.rm = TRUE), quantile(combined$aridity_index_UNEP, 0.975, na.rm = TRUE))
ggplot(combined[combined$y < 0, ], aes(x = RPD_significance, y = aridity_index_UNEP, fill = RPD_significance)) + geom_violin(trim=FALSE) + geom_boxplot(width=0.1, fill="white") + labs(title="Aridity vs. RPD significance, Southern Hemisphere", x="RPD_significance", y = "UNEP aridity index") + ylim(quantile(combined$aridity_index_UNEP, 0.025, na.rm = TRUE), quantile(combined$aridity_index_UNEP, 0.975, na.rm = TRUE))

ggplot(combined[combined$y > 0, ], aes(x = RPD_significance, y = BIOCLIM_1, fill = RPD_significance)) + geom_violin(trim=FALSE) + geom_boxplot(width=0.1, fill="white") + labs(title="Mean annual temperature vs. \nRPD significance", x="RPD_significance", y = "Mean annual temperature") + ylim(quantile(combined$BIOCLIM_1, 0.025, na.rm = TRUE), quantile(combined$BIOCLIM_1, 0.975, na.rm = TRUE))
ggplot(combined[combined$y < 0, ], aes(x = RPD_significance, y = BIOCLIM_1, fill = RPD_significance)) + geom_violin(trim=FALSE) + geom_boxplot(width=0.1, fill="white") + labs(title="Mean annual temperature vs. \nRPD significance", x="RPD_significance", y = "Mean annual temperature") + ylim(quantile(combined$BIOCLIM_1, 0.025, na.rm = TRUE), quantile(combined$BIOCLIM_1, 0.975, na.rm = TRUE))

# CANAPE is not so straightforward. Segregating by E-W or N-S hemisphere, or tropics vs. non-tropics, doesn't remove bimodality
ggplot(combined_nod[combined_nod$x > -30 | combined_nod$x < -168, ], aes(x = CANAPE, y = aridity_index_UNEP, fill = CANAPE)) + geom_violin(trim=FALSE) + geom_boxplot(width=0.1, fill="white") + labs(title="Aridity vs. CANAPE significance", x="CANAPE significance", y = "UNEP aridity index") + ylim(quantile(combined$aridity_index_UNEP, 0.025, na.rm = TRUE), quantile(combined$aridity_index_UNEP, 0.975, na.rm = TRUE))
ggplot(combined_nod[combined_nod$x < -30 & combined_nod$x > -168, ], aes(x = CANAPE, y = aridity_index_UNEP, fill = CANAPE)) + geom_violin(trim=FALSE) + geom_boxplot(width=0.1, fill="white") + labs(title="Aridity vs. CANAPE significance", x="CANAPE significance", y = "UNEP aridity index") + ylim(quantile(combined$aridity_index_UNEP, 0.025, na.rm = TRUE), quantile(combined$aridity_index_UNEP, 0.975, na.rm = TRUE))
ggplot(combined[combined$y < 0, ], aes(x = CANAPE, y = aridity_index_UNEP, fill = CANAPE)) + geom_violin(trim=FALSE) + geom_boxplot(width=0.1, fill="white") + labs(title="Aridity vs. CANAPE significance", x="CANAPE significance", y = "UNEP aridity index") + ylim(quantile(combined$aridity_index_UNEP, 0.025, na.rm = TRUE), quantile(combined$aridity_index_UNEP, 0.975, na.rm = TRUE))
ggplot(combined[combined$y < 0, ], aes(x = CANAPE, y = aridity_index_UNEP, fill = CANAPE)) + geom_violin(trim=FALSE) + geom_boxplot(width=0.1, fill="white") + labs(title="Aridity vs. CANAPE significance", x="CANAPE significance", y = "UNEP aridity index") + ylim(quantile(combined$aridity_index_UNEP, 0.025, na.rm = TRUE), quantile(combined$aridity_index_UNEP, 0.975, na.rm = TRUE))
ggplot(combined[combined$y > 23.43629 | combined$y < -23.43629, ], aes(x = CANAPE, y = aridity_index_UNEP, fill = CANAPE)) + geom_violin(trim=FALSE) + geom_boxplot(width=0.1, fill="white") + labs(title="Aridity vs. CANAPE significance", x="CANAPE significance", y = "UNEP aridity index") + ylim(quantile(combined$aridity_index_UNEP, 0.025, na.rm = TRUE), quantile(combined$aridity_index_UNEP, 0.975, na.rm = TRUE))
ggplot(combined[combined$y < 23.43629 && combined$y > -23.43629, ], aes(x = CANAPE, y = aridity_index_UNEP, fill = CANAPE)) + geom_violin(trim=FALSE) + geom_boxplot(width=0.1, fill="white") + labs(title="Aridity vs. CANAPE significance", x="CANAPE significance", y = "UNEP aridity index") + ylim(quantile(combined$aridity_index_UNEP, 0.025, na.rm = TRUE), quantile(combined$aridity_index_UNEP, 0.975, na.rm = TRUE))



########################
# Add nodulation response
########################

## Add nodulation -- doing this separately to minimize cell dropout for above
proportion_nodulating <- read.csv("./Mimosoid_CSVs_ToShare_WGS84/propNodulating_20km.csv")
proportion_nodulating$richness <- NULL
names(proportion_nodulating) <- c("x", "y", "prop_nod")
proportion_nodulating$x <- round(proportion_nodulating$x, digit = 1)
proportion_nodulating$y <- round(proportion_nodulating$y, digit = 1)
proportion_nodulating %>% group_by(x, y) %>% summarize_if(is.numeric, mean, na.rm = TRUE) -> proportion_nodulating
proportion_nodulating <- as.data.frame(proportion_nodulating)

combined_nod <- merge(combined, proportion_nodulating, by = c("x", "y"))

# Normalize entire data frame
combined_nod.temp <- combined_nod
combined_nod.scaled <- rapply(combined_nod.temp, scale, c("numeric","integer"), how="replace")
combined_nod.scaled$SR <- combined_nod$SR
combined_nod.scaled <- as.data.frame(combined_nod.scaled)
combined_nod.scaled$y <- combined_nod.temp$y
combined_nod.scaled$x <- combined_nod.temp$x

########################
# Nodule model
########################

linear_model_complex <- lm(prop_nod ~ aridity_index_UNEP + BIOCLIM_1 + BIOCLIM_12 + BIOCLIM_7 + BIOCLIM_17 + ISRICSOILGRIDS_new_average_nitrogen_reduced + ISRICSOILGRIDS_new_average_phx10percent_reduced + ISRICSOILGRIDS_new_average_soilorganiccarboncontent_reduced, data = combined_nod.scaled)
# Top 5 predictors by GLM normalized coefficient
linear_model_simple <- lm(prop_nod ~ aridity_index_UNEP + BIOCLIM_12 + BIOCLIM_17 + ISRICSOILGRIDS_new_average_nitrogen_reduced + ISRICSOILGRIDS_new_average_phx10percent_reduced, data = combined_nod.scaled)
library(lme4)
mixed_model_complex <- lmer(prop_nod ~ aridity_index_UNEP + BIOCLIM_1 + BIOCLIM_12 + BIOCLIM_7 + BIOCLIM_17 + ISRICSOILGRIDS_new_average_nitrogen_reduced + ISRICSOILGRIDS_new_average_phx10percent_reduced + ISRICSOILGRIDS_new_average_soilorganiccarboncontent_reduced + (1 | y) + (1 | x), na.action = na.omit, data = combined_nod.scaled)
# Top 5 predictors by LMM normalized coefficient
mixed_model_simple <- lmer(prop_nod ~ aridity_index_UNEP + BIOCLIM_12 + BIOCLIM_17 + ISRICSOILGRIDS_new_average_nitrogen_reduced + ISRICSOILGRIDS_new_average_soilorganiccarboncontent_reduced + (1 | y) + (1 | x), na.action = na.omit, data = combined_nod.scaled)

mixed_model_noenvironment <- lmer(prop_nod ~ (1 | y) + (1 | x), na.action = na.omit, data = combined_nod.scaled)

AIC(linear_model_complex)
AIC(linear_model_simple)
AIC(mixed_model_complex)
AIC(mixed_model_simple)
AIC(mixed_model_noenvironment)
# Simple mixed model favored

summary(mixed_model_simple)

library(MuMIn)
r.squaredGLMM(mixed_model_simple)
# First number is marginal (fixed effects only) and conditional (entire model)

########################
# Nodulation plots
########################

# BIO12 looks very similar to aridity, just show that in place of other precipitation variables
library(ggplot2)
ggplot(combined_nod, aes(x = aridity_index_UNEP, y = prop_nod) ) + geom_hex(bins = 28) + scale_fill_continuous(type = "viridis") + theme_bw() + ylim(0.5, 1) + geom_smooth(method='lm', formula = y ~ x) + labs(title="Aridity vs.\nproportion nodulating", x="Aridity index", y = "Proportion nodulating")
ggplot(combined_nod, aes(x = BIOCLIM_12, y = prop_nod) ) + geom_hex(bins = 28) + scale_fill_continuous(type = "viridis") + theme_bw() + ylim(0.5, 1) + geom_smooth(method='lm', formula = y ~ x) + labs(title="Precipitation vs.\nproportion nodulating", x="Annual precipitation", y = "Proportion nodulating")
# Exclude 1 here, see 0-inflation in text
ggplot(combined_nod, aes(x = BIOCLIM_7, y = prop_nod) ) + geom_hex(bins = 38) + scale_fill_continuous(type = "viridis") + theme_bw() + ylim(0.5, 0.99) + geom_smooth(method='lm', formula = y ~ x) + labs(title="Temperature annual range vs.\nproportion nodulating", x="BIO7", y = "Proportion nodulating")
ggplot(combined_nod, aes(x = ISRICSOILGRIDS_new_average_soilorganiccarboncontent_reduced, y = prop_nod) ) + geom_hex(bins = 38) + scale_fill_continuous(type = "viridis") + theme_bw() + ylim(0.5, 0.99) + geom_smooth(method='lm', formula = y ~ x) + labs(title="Carbon vs.\nproportion nodulating", x="Carbon content", y = "Proportion nodulating")



#ggplot(combined_nod, aes(x = ISRICSOILGRIDS_new_average_phx10percent_reduced, y = prop_nod) ) + geom_hex(bins = 70) + scale_fill_continuous(type = "viridis") + theme_bw() + ylim(0.01, 0.75) + geom_smooth(method='lm', formula = y ~ x) + labs(title="pH vs.\nproportion nodulating", x="pH", y = "Proportion nodulating")



combined_nod$non_nod <- as.factor(ifelse(combined_nod$prop_nod < 1, "Present", "Absent"))
# Aridity relationship not affected by hemisphere
ggplot(combined_nod, aes(x = non_nod, y = aridity_index_UNEP, fill = non_nod)) + geom_violin(trim=FALSE) + geom_boxplot(width=0.1, fill="white") + labs(title="Aridity vs. presence of nodulators", x="Nodulator presence", y = "UNEP aridity index") + ylim(quantile(combined$aridity_index_UNEP, 0.025, na.rm = TRUE), quantile(combined$aridity_index_UNEP, 0.975, na.rm = TRUE))
ggplot(combined_nod[combined_nod$y > 0, ], aes(x = non_nod, y = aridity_index_UNEP, fill = non_nod)) + geom_violin(trim=FALSE) + geom_boxplot(width=0.1, fill="white") + labs(title="Aridity vs. presence of nodulators,\nNorthern Hemisphere", x="Nodulator presence", y = "UNEP aridity index") + ylim(quantile(combined$aridity_index_UNEP, 0.025, na.rm = TRUE), quantile(combined$aridity_index_UNEP, 0.975, na.rm = TRUE))
ggplot(combined_nod[combined_nod$y < 0, ], aes(x = non_nod, y = aridity_index_UNEP, fill = non_nod)) + geom_violin(trim=FALSE) + geom_boxplot(width=0.1, fill="white") + labs(title="Aridity vs. presence of nodulators,\nSouthern Hemisphere", x="Nodulator presence", y = "UNEP aridity index") + ylim(quantile(combined$aridity_index_UNEP, 0.025, na.rm = TRUE), quantile(combined$aridity_index_UNEP, 0.975, na.rm = TRUE))
ggplot(combined_nod[combined_nod$x > -30 | combined_nod$x < -168, ], aes(x = non_nod, y = aridity_index_UNEP, fill = non_nod)) + geom_violin(trim=FALSE) + geom_boxplot(width=0.1, fill="white") + labs(title="Aridity vs. presence of nodulators,\nWestern Hemisphere", x="Nodulator presence", y = "UNEP aridity index") + ylim(quantile(combined$aridity_index_UNEP, 0.025, na.rm = TRUE), quantile(combined$aridity_index_UNEP, 0.975, na.rm = TRUE))
ggplot(combined_nod[combined_nod$x < -30 & combined_nod$x > -168, ], aes(x = non_nod, y = aridity_index_UNEP, fill = non_nod)) + geom_violin(trim=FALSE) + geom_boxplot(width=0.1, fill="white") + labs(title="Aridity vs. presence of nodulators,\nEastern Hemisphere", x="Nodulator presence", y = "UNEP aridity index") + ylim(quantile(combined$aridity_index_UNEP, 0.025, na.rm = TRUE), quantile(combined$aridity_index_UNEP, 0.975, na.rm = TRUE))


# Temperature is though
ggplot(combined_nod, aes(x = non_nod, y = BIOCLIM_1/10, fill = non_nod)) + geom_violin(trim=FALSE) + geom_boxplot(width=0.1, fill="white") + labs(title="Mean annual temperature vs. \npresence of nodulators", x="Nodulator presence", y = "Mean annual temperature")
ggplot(combined_nod[combined_nod$y > 0, ], aes(x = non_nod, y = BIOCLIM_1/10, fill = non_nod)) + geom_violin(trim=FALSE) + geom_boxplot(width=0.1, fill="white") + labs(title="Mean annual temperature vs. \npresence of nodulators,\nNorthern Hemisphere", x="Nodulator presence", y = "Mean annual temperature")
ggplot(combined_nod[combined_nod$y < 0, ], aes(x = non_nod, y = BIOCLIM_1/10, fill = non_nod)) + geom_violin(trim=FALSE) + geom_boxplot(width=0.1, fill="white") + labs(title="Mean annual temperature vs. \npresence of nodulators,\nSouthern Hemisphere", x="Nodulator presence", y = "Mean annual temperature")
ggplot(combined_nod[combined_nod$x > -30 | combined_nod$x < -168, ], aes(x = non_nod, y = BIOCLIM_1/10, fill = non_nod)) + geom_violin(trim=FALSE) + geom_boxplot(width=0.1, fill="white") + labs(title="Mean annual temperature vs. \npresence of nodulators,\nWestern Hemisphere", x="Nodulator presence", y = "Mean annual temperature")
ggplot(combined_nod[combined_nod$x < -30 & combined_nod$x > -168, ], aes(x = non_nod, y = BIOCLIM_1/10, fill = non_nod)) + geom_violin(trim=FALSE) + geom_boxplot(width=0.1, fill="white") + labs(title="Mean annual temperature vs. \npresence of nodulators,\nEastern Hemisphere", x="Nodulator presence", y = "Mean annual temperature")


