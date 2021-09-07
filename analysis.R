# Meta-Analysis of Gut Passage Review studies

# NOTE: We present two versions of the dataset. The raw version has data (and
# many less-relevant columns) that were not included in the analyses presented
# in the paper. For example, experiments involving non-fleshy fruits or
# instances where other less common experimental treatments were applied. We
# share these raw data for potential use in future analyses. The code in the
# section "Load raw data and manipulate" show how certain studies were removed
# from analysis and how we filtered the raw data for analysis and handled
# missing or alternately formatted raw data.
# 
# The clean version shows the data that fit the conditions for inclusion in the
# analysis and the most relevant columns.



### Functions ---------------------------------------------------------

ipak <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg)) 
    install.packages(new.pkg, dependencies = TRUE)
  sapply(pkg, require, character.only = TRUE)
}

factor.by.frequency <- function(x){
  x <- trimws(tolower(x))
  tb <- table(x)
  factor(x, levels = names(tb[order(tb, decreasing = T)]))
}



### Packages ---------------------------------------------------------


# Load libraries
packages <- c("metafor", "dplyr", "tidyr", "Taxonstand", "plotrix", "plyr",
              "ape", "httr", "brranching",
              "multcomp")
ipak(packages)



### Load raw data and manipulate -----------------------------------------------
# Data on plant-frugivore combinations

gut.dat <- read.csv("gut_passage_germination_data_raw.csv", stringsAsFactors = T)


# Want to produce a citation for a references table
gut.dat$citation <- paste0(gut.dat$firstauth, " (",
                           gut.dat$Year,") ",
                           gut.dat$Title, ". ",
                           gut.dat$Source.title, ". " ,
                           gut.dat$Volume, ". ",
                           gut.dat$Page.start, "-",
                           gut.dat$Page.end, ".")
gut.dat$citation <- gsub("NA. ", "", gut.dat$citation)
gut.dat$citation <- gsub("NA-NA", "", gut.dat$citation) %>% trimws()
gut.dat$citation <- gsub(" .", "", gut.dat$citation, fixed = T)


# Deal with data being read in as factors
gut.dat$n.gut <- as.numeric(as.character(gut.dat$n.gut)) # These NAs are okay - these are true NAs
gut.dat$n.control <- as.numeric(as.character(gut.dat$n.control)) # These NAs are okay - these are true NAs
gut.dat$gut.portion <- as.numeric(as.character(gut.dat$gut.portion)) # These NAs are okay - these are true NAs
gut.dat$control.portion <- as.numeric(as.character(gut.dat$control.portion)) # These NAs are okay - these are true NAs

# Remove records that don't fit our requirements #
levels(as.factor(gut.dat$compare.against))
gut.dat<-gut.dat %>% filter(compare.against ==  "mechanically cleaned" | compare.against == "whole fruit") 
gut.dat<-gut.dat %>% filter(plant.genus != "Vine" & plant.genus != "Plant")
# Keep only studies on fleshy fruits
gut.dat <- gut.dat[which(gut.dat$fleshy == "yes"), ]

# There are 21 recorded from on riverine islands - we have decided to treat 
# these as mainland sites for analysis.
gut.dat$mainisland <- revalue(gut.dat$mainisland, c("riverine island" = "mainland"))

# There's a study that involved seeds from both feeding trials and field 
# collection, and I will treat that as a "yes" for feedtrial
gut.dat$feedtrial <- revalue(gut.dat$feedtrial, c("both" = "yes"))

# Make sure paper ID is a factor
gut.dat$numid <- factor(gut.dat$numid)

# Make latitude region more consistent
gut.dat$lat_region <- gsub("\\s*\\([^\\)]+\\)","",as.character(gut.dat$lat_region)) %>% factor()

# Want a single column for plant species
gut.dat$plant.spp <- paste(gut.dat$plant.genus, gut.dat$plant.species)
gut.dat$frug.spp <- paste(gut.dat$frug.genus, gut.dat$frug.species)




### Get phylogenetic data ---------------------------------------------------------

# Plants
# https://cran.r-project.org/web/packages/rotl/vignettes/data_mashups.html
plant.vec <- unique(paste(gut.dat$plant.genus, gut.dat$plant.species, sep = " "))
plant.tree <- brranching::phylomatic(plant.vec, get = "POST")

# There's also something with class phylomatic in here, so writing
# and reading back in so this talks well with rbladj below
write.tree(plant.tree, "plant.tree")
plant.tree <- read.tree("plant.tree")

# Get node ages
wikstrom.ages <- read.table("https://raw.githubusercontent.com/phylocom/phylocom/master/example_data/bladj/wikstrom.ages")
plant.tree <- rbladj(tree = plant.tree, ages = wikstrom.ages)


### Prep effect size data for analysis -----------------------------------------

# We will assume median sample sizes when not reported
median.n.gut <- median(as.numeric(gut.dat$n.gut), na.rm=T)
median.n.control <- median(as.numeric(gut.dat$n.control), na.rm=T)

gut.dat$n.gut <- ifelse(is.na(gut.dat$n.gut), median.n.gut, gut.dat$n.gut)
gut.dat$n.control <- ifelse(is.na(gut.dat$n.control), median.n.control, gut.dat$n.control)


# Need to change to whole counts
# Note that we will record percentages because studies typically use
# proportions/percentages and dont always report sample sizes.

# ai is number of germ, gut-passed
# bi is number of dead, gut-passed
# ci is number of germ, control
# di is number of dead, control

# Get these in the escalc description, although we could use more easy to interpret colnames
gut.dat$ai <- round(gut.dat$gut.portion * gut.dat$n.gut)
gut.dat$bi <- gut.dat$n.gut - gut.dat$ai

gut.dat$ci <- round(gut.dat$control.portion * gut.dat$n.control)
gut.dat$di <- gut.dat$n.control - gut.dat$ci


### Prep sleeker dataset for analysis ----------------------------------------

colnames(gut.dat) <- tolower(colnames(gut.dat))

# Keep only studies where proportions were reported
gut.dat <- gut.dat[complete.cases(gut.dat[c("ai", "bi", "ci", "di")]), ] # There just weren't data for about 245 rows

combos <- gut.dat #rename, so we can run code below. 



# Remove columns relevant to data entry notes
combos <- combos[, c("numid", "plant.genus", "plant.species", "plant.spp",
                     "frug.taxon", "family", "frug.genus", "frug.species", "frug.spp",
                     "lat_region", "mainisland", "old_new_world", "year", "country",
                     "plant.invasive", "frug.invasive",
                     "plant.iucn.category", "frug.iucn.category",
                     "feedtrial", "germtrial",
                     "planting.location",
                     "germ.effect", "rate.effect",
                     "compare.against",
                     "fleshy",
                     "gut.portion", "n.gut",
                     "control.portion", "n.control",
                     "maybe.regurg",
                     "ai", "bi", "ci", "di",
                     "citation")]


# Get values for effect size (yi) and sampling variances (vi)

combos <- escalc(measure = "OR", 
                 ai = ai, bi = bi, ci = ci, di = di,
                 data = combos)



# Need complete cases for the variables used in the analysis

combos.set <- combos[complete.cases(combos[,c("yi", "vi", 
                                              "compare.against", 
                                              "planting.location", # There are NAs for planting location...
                                              "feedtrial", 
                                              "frug.taxon",
                                              "lat_region",
                                              "mainisland",
                                              "old_new_world")]),]

# Will save the analyzed version here, and a table of references
write.csv(combos.set, "gut_passage_germination_data_clean.csv", row.names = F)
tibble(Reference = combos.set$citation %>% unique() %>% sort()) %>% 
  write.csv("references.csv")



# Sort by factor level frequency ----------------------------

combos.set$compare.against <- combos.set$compare.against %>% factor.by.frequency()
combos.set$feedtrial <- combos.set$feedtrial %>% factor.by.frequency()
combos.set$planting.location <- combos.set$planting.location %>% factor.by.frequency()
combos.set$frug.taxon <- combos.set$frug.taxon %>% factor.by.frequency()
combos.set$frug.invasive <- combos.set$frug.invasive %>% factor.by.frequency()
combos.set$plant.invasive <- combos.set$plant.invasive %>% factor.by.frequency()
combos.set$lat_region <- combos.set$lat_region %>% factor.by.frequency()
combos.set$mainisland <- combos.set$mainisland %>% factor.by.frequency()



# Run (non-phylo) models ----------------------------

# Note that these take a very long time - you can simply pull in the model
# outputs just below

# best.mod.ran <- rma.mv(yi = yi, V = vi,
#                        mods = ~ compare.against + feedtrial + planting.location + frug.taxon + mainisland + lat_region,
#                        method = "ML",
#                        random = list(~ 1 | plant.spp, ~ 1 | frug.spp), 
#                        data = combos.set)
# 
# full.mod.ran <- rma.mv(yi = yi, V = vi,
#                        mods = ~ compare.against + feedtrial + planting.location + frug.taxon + frug.invasive + plant.invasive + lat_region + mainisland,
#                        method = "ML",
#                        random = list(~ 1 | plant.spp, ~ 1 | frug.spp), 
#                        data = combos.set)
# 
# dredged.full.mod.ran <- MuMIn::dredge(full.mod.ran)
# 
# save(best.mod.ran, full.mod.ran, file = "model_outputs.RData")

# Will pull saved versions of these in:
load("model_outputs.RData")

AIC(best.mod.ran, full.mod.ran)
summary(best.mod.ran)
summary(full.mod.ran)




# Run phylogenetic model ----------------------------

combos.set$plant.phylo <- tolower(paste(combos.set$plant.genus, combos.set$plant.species, sep = "_"))
combos.set$plant.tax <- tolower(paste(combos.set$plant.genus, combos.set$plant.species, sep = "_"))

combos.set.phy <- combos.set[combos.set$plant.phylo %in% plant.tree$tip.label, ]
dim(combos.set.phy)
dim(combos.set)

R.plant <- vcv(plant.tree, corr=T) # Gets var covar structure from tree

time0 <- Sys.time()
full.phy.mod <-rma.mv(yi, vi,
                      mods = ~ compare.against + feedtrial + planting.location + frug.taxon + frug.invasive + plant.invasive + lat_region + mainisland,
                      random = list(~ 1 | plant.phylo, ~ 1 | frug.spp),
                      R = list(plant.phylo = R.plant),
                      data = combos.set.phy)
time1 <- Sys.time()
full.phy.mod





# Linear hypothesis testing ----------------------

# Reference levels = tropical mainland bird using petri dish and feeding trials.
c1.petri <- c(1, # compare.againstwhole fruit
              1, # feedtrialyes
              0, # planting.locationgreenhouse/nursery soil
              0, # planting.locationfield soil
              0, # planting.locationother
              0, # frug.taxonnon-flying mammal other than primate
              0, # frug.taxonprimate
              0, # frug.taxonreptile
              0, # frug.taxonbat
              0, # frug.taxoninvertebrate
              0, # frug.taxonfish
              0, # mainislandisland
              0, # lat_regiontemperate
              0) # lat_regionsubtropical
c1.field <- c(1, # compare.againstwhole fruit
              1, # feedtrialyes
              0, # planting.locationgreenhouse/nursery soil
              1, # planting.locationfield soil
              0, # planting.locationother
              0, # frug.taxonnon-flying mammal other than primate
              0, # frug.taxonprimate
              0, # frug.taxonreptile
              0, # frug.taxonbat
              0, # frug.taxoninvertebrate
              0, # frug.taxonfish
              0, # mainislandisland
              0, # lat_regiontemperate
              0) # lat_regionsubtropical
c1.nursery <- c(1, # compare.againstwhole fruit
                1, # feedtrialyes
                1, # planting.locationgreenhouse/nursery soil
                0, # planting.locationfield soil
                0, # planting.locationother
                0, # frug.taxonnon-flying mammal other than primate
                0, # frug.taxonprimate
                0, # frug.taxonreptile
                0, # frug.taxonbat
                0, # frug.taxoninvertebrate
                0, # frug.taxonfish
                0, # mainislandisland
                0, # lat_regiontemperate
                0) # lat_regionsubtropical
c1.other <- c(1, # compare.againstwhole fruit
              1, # feedtrialyes
              0, # planting.locationgreenhouse/nursery soil
              0, # planting.locationfield soil
              1, # planting.locationother
              0, # frug.taxonnon-flying mammal other than primate
              0, # frug.taxonprimate
              0, # frug.taxonreptile
              0, # frug.taxonbat
              0, # frug.taxoninvertebrate
              0, # frug.taxonfish
              0, # mainislandisland
              0, # lat_regiontemperate
              0) # lat_regionsubtropical


c1.feedtrialyes <- c(1, # compare.againstwhole fruit
                     1, # feedtrialyes
                     0, # planting.locationgreenhouse/nursery soil
                     0, # planting.locationfield soil
                     0, # planting.locationother
                     0, # frug.taxonnon-flying mammal other than primate
                     0, # frug.taxonprimate
                     0, # frug.taxonreptile
                     0, # frug.taxonbat
                     0, # frug.taxoninvertebrate
                     0, # frug.taxonfish
                     0, # mainislandisland
                     0, # lat_regiontemperate
                     0) # lat_regionsubtropical
c1.feedtrialno <- c(1, # compare.againstwhole fruit
                    0, # feedtrialyes
                    0, # planting.locationgreenhouse/nursery soil
                    0, # planting.locationfield soil
                    0, # planting.locationother
                    0, # frug.taxonnon-flying mammal other than primate
                    0, # frug.taxonprimate
                    0, # frug.taxonreptile
                    0, # frug.taxonbat
                    0, # frug.taxoninvertebrate
                    0, # frug.taxonfish
                    0, # mainislandisland
                    0, # lat_regiontemperate
                    0) # lat_regionsubtropical



c2.scar <- c(0, # compare.againstwhole fruit
             1, # feedtrialyes
             0, # planting.locationgreenhouse/nursery soil
             0, # planting.locationfield soil
             0, # planting.locationother
             0, # frug.taxonnon-flying mammal other than primate
             0, # frug.taxonprimate
             0, # frug.taxonreptile
             0, # frug.taxonbat
             0, # frug.taxoninvertebrate
             0, # frug.taxonfish
             0, # mainislandisland
             0, # lat_regiontemperate
             0) # lat_regionsubtropical
c2.deinhib <- c(1, # compare.againstwhole fruit
                1, # feedtrialyes
                0, # planting.locationgreenhouse/nursery soil
                0, # planting.locationfield soil
                0, # planting.locationother
                0, # frug.taxonnon-flying mammal other than primate
                0, # frug.taxonprimate
                0, # frug.taxonreptile
                0, # frug.taxonbat
                0, # frug.taxoninvertebrate
                0, # frug.taxonfish
                0, # mainislandisland
                0, # lat_regiontemperate
                0) # lat_regionsubtropical

c3.bird <- c(1, # compare.againstwhole fruit
             1, # feedtrialyes
             0, # planting.locationgreenhouse/nursery soil
             0, # planting.locationfield soil
             0, # planting.locationother
             0, # frug.taxonnon-flying mammal other than primate
             0, # frug.taxonprimate
             0, # frug.taxonreptile
             0, # frug.taxonbat
             0, # frug.taxoninvertebrate
             0, # frug.taxonfish
             0, # mainislandisland
             0, # lat_regiontemperate
             0) # lat_regionsubtropical
c3.bat <- c(1, # compare.againstwhole fruit
            1, # feedtrialyes
            0, # planting.locationgreenhouse/nursery soil
            0, # planting.locationfield soil
            0, # planting.locationother
            0, # frug.taxonnon-flying mammal other than primate
            0, # frug.taxonprimate
            0, # frug.taxonreptile
            1, # frug.taxonbat
            0, # frug.taxoninvertebrate
            0, # frug.taxonfish
            0, # mainislandisland
            0, # lat_regiontemperate
            0) # lat_regionsubtropical
c3.primate <- c(1, # compare.againstwhole fruit
                1, # feedtrialyes
                0, # planting.locationgreenhouse/nursery soil
                0, # planting.locationfield soil
                0, # planting.locationother
                0, # frug.taxonnon-flying mammal other than primate
                1, # frug.taxonprimate
                0, # frug.taxonreptile
                0, # frug.taxonbat
                0, # frug.taxoninvertebrate
                0, # frug.taxonfish
                0, # mainislandisland
                0, # lat_regiontemperate
                0) # lat_regionsubtropical
c3.other <- c(1, # compare.againstwhole fruit
              1, # feedtrialyes
              0, # planting.locationgreenhouse/nursery soil
              0, # planting.locationfield soil
              0, # planting.locationother
              1, # frug.taxonnon-flying mammal other than primate
              0, # frug.taxonprimate
              0, # frug.taxonreptile
              0, # frug.taxonbat
              0, # frug.taxoninvertebrate
              0, # frug.taxonfish
              0, # mainislandisland
              0, # lat_regiontemperate
              0) # lat_regionsubtropical
c3.reptile <- c(1, # compare.againstwhole fruit
                1, # feedtrialyes
                0, # planting.locationgreenhouse/nursery soil
                0, # planting.locationfield soil
                0, # planting.locationother
                0, # frug.taxonnon-flying mammal other than primate
                0, # frug.taxonprimate
                1, # frug.taxonreptile
                0, # frug.taxonbat
                0, # frug.taxoninvertebrate
                0, # frug.taxonfish
                0, # mainislandisland
                0, # lat_regiontemperate
                0) # lat_regionsubtropical
c3.fish <- c(1, # compare.againstwhole fruit
             1, # feedtrialyes
             0, # planting.locationgreenhouse/nursery soil
             0, # planting.locationfield soil
             0, # planting.locationother
             0, # frug.taxonnon-flying mammal other than primate
             0, # frug.taxonprimate
             0, # frug.taxonreptile
             0, # frug.taxonbat
             0, # frug.taxoninvertebrate
             1, # frug.taxonfish
             0, # mainislandisland
             0, # lat_regiontemperate
             0) # lat_regionsubtropical
c3.invertebrate <- c(1, # compare.againstwhole fruit
                     1, # feedtrialyes
                     0, # planting.locationgreenhouse/nursery soil
                     0, # planting.locationfield soil
                     0, # planting.locationother
                     0, # frug.taxonnon-flying mammal other than primate
                     0, # frug.taxonprimate
                     0, # frug.taxonreptile
                     0, # frug.taxonbat
                     1, # frug.taxoninvertebrate
                     0, # frug.taxonfish
                     0, # mainislandisland
                     0, # lat_regiontemperate
                     0) # lat_regionsubtropical



c4.tropical <- c(1, # compare.againstwhole fruit
                 1, # feedtrialyes
                 0, # planting.locationgreenhouse/nursery soil
                 0, # planting.locationfield soil
                 0, # planting.locationother
                 0, # frug.taxonnon-flying mammal other than primate
                 0, # frug.taxonprimate
                 0, # frug.taxonreptile
                 0, # frug.taxonbat
                 0, # frug.taxoninvertebrate
                 0, # frug.taxonfish
                 0, # mainislandisland
                 0, # lat_regiontemperate
                 0) # lat_regionsubtropical
c4.subtropical <- c(1, # compare.againstwhole fruit
                    1, # feedtrialyes
                    0, # planting.locationgreenhouse/nursery soil
                    0, # planting.locationfield soil
                    0, # planting.locationother
                    0, # frug.taxonnon-flying mammal other than primate
                    0, # frug.taxonprimate
                    0, # frug.taxonreptile
                    0, # frug.taxonbat
                    0, # frug.taxoninvertebrate
                    0, # frug.taxonfish
                    0, # mainislandisland
                    0, # lat_regiontemperate
                    1) # lat_regionsubtropical
c4.temperate <- c(1, # compare.againstwhole fruit
                  1, # feedtrialyes
                  0, # planting.locationgreenhouse/nursery soil
                  0, # planting.locationfield soil
                  0, # planting.locationother
                  0, # frug.taxonnon-flying mammal other than primate
                  0, # frug.taxonprimate
                  0, # frug.taxonreptile
                  0, # frug.taxonbat
                  0, # frug.taxoninvertebrate
                  0, # frug.taxonfish
                  0, # mainislandisland
                  1, # lat_regiontemperate
                  0) # lat_regionsubtropical


c4.mainland <- c(1, # compare.againstwhole fruit
                 1, # feedtrialyes
                 0, # planting.locationgreenhouse/nursery soil
                 0, # planting.locationfield soil
                 0, # planting.locationother
                 0, # frug.taxonnon-flying mammal other than primate
                 0, # frug.taxonprimate
                 0, # frug.taxonreptile
                 0, # frug.taxonbat
                 0, # frug.taxoninvertebrate
                 0, # frug.taxonfish
                 0, # mainislandisland
                 0, # lat_regiontemperate
                 0) # lat_regionsubtropical
c4.island <- c(1, # compare.againstwhole fruit
               1, # feedtrialyes
               0, # planting.locationgreenhouse/nursery soil
               0, # planting.locationfield soil
               0, # planting.locationother
               0, # frug.taxonnon-flying mammal other than primate
               0, # frug.taxonprimate
               0, # frug.taxonreptile
               0, # frug.taxonbat
               0, # frug.taxoninvertebrate
               0, # frug.taxonfish
               1, # mainislandisland
               0, # lat_regiontemperate
               0) # lat_regionsubtropical



add0 <- function(x) c(0,x)

# Fig 1a
c1a <- rbind(c1.petri, c1.field, c1.nursery, c1.other)
c1a.preds <- predict(best.mod.ran, c1a) %>% as.data.frame()
c1a.glht <- glht(best.mod.ran, 
                 linfct=rbind(add0(c1.petri) - add0(c1.field),
                              add0(c1.petri) - add0(c1.nursery),
                              add0(c1.petri) - add0(c1.other),
                              add0(c1.field) - add0(c1.nursery),
                              add0(c1.field) - add0(c1.other),
                              add0(c1.nursery) - add0(c1.other))) %>% summary(test = adjusted("none"))
# Fig 1b
c1b <- rbind(c1.feedtrialyes, c1.feedtrialno)
c1b.preds <- predict(best.mod.ran, c1b) %>% as.data.frame()
c1b.glht <- glht(best.mod.ran, 
                 linfct=rbind(add0(c1.feedtrialyes) - add0(c1.feedtrialno))) %>% summary(test = adjusted("none"))


# Fig 2
c2.preds <- c1b.preds[1:2,]
c2.preds[1:2,] <- NA
c2.preds$pred <- summary(best.mod.ran)$b[c("intrcpt", "compare.againstwhole fruit"),]
c2.preds$ci.lb <- summary(best.mod.ran)$ci.lb[1:2]
c2.preds$ci.ub <- summary(best.mod.ran)$ci.ub[1:2]
c2.glht <- glht(best.mod.ran, 
                linfct=rbind(add0(c2.scar) - add0(c2.deinhib))) %>% summary(test = adjusted("none"))


# Fig 3
c3 <- rbind(c3.bird, c3.bat, c3.primate, c3.other, c3.reptile, c3.fish, c3.invertebrate)
c3.preds <- predict(best.mod.ran, c3) %>% as.data.frame()
c3.glht <- glht(best.mod.ran, 
                linfct=rbind(add0(c3.bird) - add0(c3.bat), #1
                             add0(c3.bird) - add0(c3.primate), #2
                             add0(c3.bird) - add0(c3.other), #3
                             add0(c3.bird) - add0(c3.reptile), #4
                             add0(c3.bird) - add0(c3.fish), #5
                             add0(c3.bird) - add0(c3.invertebrate), #6
                             add0(c3.bat) - add0(c3.primate), #7
                             add0(c3.bat) - add0(c3.other), #8
                             add0(c3.bat) - add0(c3.reptile), #9
                             add0(c3.bat) - add0(c3.fish), #10
                             add0(c3.bat) - add0(c3.invertebrate), #11
                             add0(c3.primate) - add0(c3.other), #12
                             add0(c3.primate) - add0(c3.reptile), #13
                             add0(c3.primate) - add0(c3.fish), # 14
                             add0(c3.primate) - add0(c3.invertebrate), #15
                             add0(c3.other) - add0(c3.reptile), #16
                             add0(c3.other) - add0(c3.fish), #17
                             add0(c3.other) - add0(c3.invertebrate), #18
                             add0(c3.reptile) - add0(c3.fish), #19
                             add0(c3.reptile) - add0(c3.invertebrate), #20
                             add0(c3.fish) - add0(c3.invertebrate))) %>% summary(test = adjusted("none"))

c3.cld <- matrix(nrow=nrow(c3), ncol=nrow(c3)) %>% as.data.frame()
c3.cld[t(combn(nrow(c3), 2))] <- ifelse(c3.glht$test$pvalues < 0.05,
                                        ifelse(c3.glht$test$coefficients > 0,1,-1),0) %>% as.vector()
colnames(c3.cld) <- rownames(c3.cld) <- c("bird", "bat", "primate", "other", "reptile", "fish", "invertebrate")


# Fig 4a
c4a <- rbind(c4.tropical, c4.subtropical, c4.temperate)
c4a.preds <- predict(best.mod.ran, c4a) %>% as.data.frame()
c4a.glht <- glht(best.mod.ran, 
                 linfct=rbind(add0(c4.tropical) - add0(c4.subtropical),
                              add0(c4.tropical) - add0(c4.temperate),
                              add0(c4.subtropical) - add0(c4.temperate))) %>% summary(test = adjusted("none"))

c4a.cld <- matrix(nrow=nrow(c4a), ncol=nrow(c4a)) %>% as.data.frame()
c4a.cld[t(combn(nrow(c4a), 2))] <- ifelse(c4a.glht$test$pvalues < 0.05,
                                          ifelse(c4a.glht$test$coefficients > 0,1,-1),0) %>% as.vector()
colnames(c4a.cld) <- rownames(c4a.cld) <- c("tropical", "subtropical", "temperate")


# Fig 4b
c4b <- rbind(c4.mainland, c4.island)
c4b.preds <- predict(best.mod.ran, c4b) %>% as.data.frame()
c4b.glht <- glht(best.mod.ran, 
                 linfct=rbind(add0(c4.mainland) - add0(c4.island))) %>% summary(test = adjusted("none"))



# Make table for the linear hypothesis testing

get_glht_df <- function(res_glht){
  df <- cbind(res_glht$test$coefficients, 
              res_glht$test$sigma,
              res_glht$test$tstat,
              res_glht$test$pvalues)
  df <- df %>% round(digits = 3)
  df <- rbind(rep(NA, 4), df)
  colnames(df) <- c("Estimate", "Std. Error", "z value", "P value")
  df
}

table_glht <- rbind(get_glht_df(c1a.glht),
                  get_glht_df(c1b.glht),
                  get_glht_df(c2.glht),
                  get_glht_df(c3.glht),
                  get_glht_df(c4a.glht),
                  get_glht_df(c4b.glht))

test_col <- c(NA,
              "Petri - Field",
              "Petri - Nursery",
              "Petri - Other",
              "Field - Nursery",
              "Field - Other",
              "Nursery - Other",
              NA,
              "Feeding trial - Field collection",
              NA,
              "Scarification - Deinhibition",
              NA,
              "Bird - Bat",
              "Bird - Primate",
              "Bird - Other mammal",
              "Bird - Reptile",
              "Bird - Fish",
              "Bird - Invertebrate",
              "Bat - Primate",
              "Bat - Other",
              "Bat - Reptile",
              "Bat - Fish",
              "Bat - Invertebrate",
              "Primate - Other",
              "Primate - Reptile",
              "Primate - Fish",
              "Primate - Invertebrate",
              "Other - Reptile",
              "Other - Fish",
              "Other - Invertebrate",
              "Reptile - Fish",
              "Reptile - Invertebrate",
              "Fish - Invertebrate",
              NA,
              "Tropical - Subtropical",
              "Tropical - Temperate",
              "Subtropical - Temperate",
              NA,
              "Mainland - Island")

table_glht <- table_glht %>% 
  as.data.frame() %>% 
  #tibble %>% 
  dplyr::mutate(Test = test_col, .before = 1) %>% 
  dplyr::mutate(Factor = NA, .before = 1)

table_glht$Factor[which(is.na(table_glht$Test))] <- c("Sowing medium", 
                                                  "Feeding trial",
                                                  "Scarification vs. deinhibition",
                                                  "Frugivore taxon",
                                                  "Latitude region",
                                                  "Mainland vs. island")
table_glht












# Make table for the model outputs

level.names <- c("",
                 "whole fruit", 
                 "field-collected", #"captive & field",
                 "greenhouse soil", "field soil", "other",
                 "other mammals", "primates", "reptiles", "bats", "invertebrates", "fish",
                 "yes", 
                 "yes",
                 "temperate", "subtropical",
                 "island")

coef.names <- c("Intercept","Control type", "Feeding trial", "Sowing medium",
                "Frugivore taxon", "Frugivore invasive?", "Plant invasive?",
                "Latitude region", "Mainland vs. island")

n.coef <- c(1,1,1,3,6,1,1,2,1)

full.df <- coef(full.mod.ran) %>% round(2) %>% as.data.frame()
colnames(full.df) <- "estimate"
full.df$ci.lb <- summary(full.mod.ran)$ci.lb
full.df$ci.ub <- summary(full.mod.ran)$ci.ub
full.df$level.names <- level.names
full.df[cumsum(n.coef),"coef.names"] <- coef.names

full.df <- full.df[,c("coef.names", "level.names", "estimate", "ci.lb", "ci.ub")]

full.df$coef.names <- ifelse(is.na(full.df$coef.names),
                             "", full.df$coef.names)
full.df$ci <- paste("(", 
                    round(full.df$ci.lb, 2), 
                    " - ",
                    round(full.df$ci.ub, 2), 
                    ")",
                    sep = "")

full.df <- full.df[,c("coef.names", "level.names", "estimate", "ci")]


phylo.df <- coef(full.phy.mod) %>% round(2) %>% as.data.frame()
colnames(phylo.df) <- "estimate"
phylo.df$ci.lb <- summary(full.phy.mod)$ci.lb
phylo.df$ci.ub <- summary(full.phy.mod)$ci.ub
phylo.df$level.names <- level.names
phylo.df[cumsum(n.coef),"coef.names"] <- coef.names

phylo.df <- phylo.df[,c("coef.names", "level.names", "estimate", "ci.lb", "ci.ub")]

phylo.df$coef.names <- ifelse(is.na(phylo.df$coef.names),
                              "", phylo.df$coef.names)
phylo.df$ci <- paste("(", 
                     round(phylo.df$ci.lb, 2), 
                     " - ",
                     round(phylo.df$ci.ub, 2), 
                     ")",
                     sep = "")

phylo.df <- phylo.df[,c("coef.names", "level.names", "estimate", "ci")]

