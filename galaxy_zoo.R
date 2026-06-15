library(treemap)
library(tidyverse)
library(RColorBrewer)

DATA    <- "../Research/Data/galaxies_training/"
FIGURES <- "../Research/Figures/MIA Project/"
my_palette <- c(brewer.pal(9, "Pastel1"), brewer.pal(8, "Pastel2"))

# --------------------------------------------------------------------------- #
# Load and reshape classification data
# --------------------------------------------------------------------------- #

classes <- read_csv(paste0(DATA, "training_solutions_rev1.csv"),
                    show_col_types = FALSE) %>%
  pivot_longer(cols = -GalaxyID,
               names_to = "Class",
               values_to = "Probability") %>%
  mutate(Class    = str_replace(Class, "Class", ""),
         Question = gsub("\\..*", "", Class))

classes_detail <- classes %>%
  filter(
    !Class %in% c("1.1", "1.2"),
    !Question %in% c("6", "8")
  ) %>%
  mutate(
    Galaxy = ifelse(Class %in% c("7.1","7.2","7.3"), "Smooth", ""),
    Galaxy = ifelse(!Class %in% c("7.1","7.2","7.3"), "Disk", Galaxy),
    Galaxy = ifelse(Class == "1.3", "Artifact", Galaxy)
  ) %>%
  mutate(
    Question = ifelse(Question == 7,  "7. How rounded?", Question),
    Question = ifelse(Question == 2,  "2. Could this disk be viewed edge on?", Question),
    Question = ifelse(Question == 9,  "9. Does the galaxy have a bulge?", Question),
    Question = ifelse(Question == 3,  "3. Is there a bar through the center?", Question),
    Question = ifelse(Question == 4,  "4. Is there any sign of a spiral arm?", Question),
    Question = ifelse(Question == 10, "10. How tightly wound are they?", Question),
    Question = ifelse(Question == 11, "11. How many spiral arms ar there?", Question),
    Question = ifelse(Question == 5,  "5. How prominent is the bulge?", Question),
    Question = ifelse(Question == 1,  "Artifact", Question)
  ) %>%
  mutate(
    Class = ifelse(Class == "7.1",  "Completely", Class),
    Class = ifelse(Class == "7.2",  "In-between", Class),
    Class = ifelse(Class == "7.3",  "Cigar",      Class),
    Class = ifelse(Class == "2.1",  "Yes",        Class),
    Class = ifelse(Class == "2.2",  "No",         Class),
    Class = ifelse(Class == "9.1",  "Yes, round", Class),
    Class = ifelse(Class == "9.2",  "Yes, boxy",  Class),
    Class = ifelse(Class == "9.3",  "No bulge",   Class),
    Class = ifelse(Class == "3.1",  "Yes",        Class),
    Class = ifelse(Class == "3.2",  "No",         Class),
    Class = ifelse(Class == "4.1",  "Yes",        Class),
    Class = ifelse(Class == "4.2",  "No",         Class),
    Class = ifelse(Class == "10.1", "Tight",      Class),
    Class = ifelse(Class == "10.2", "Medium",     Class),
    Class = ifelse(Class == "10.3", "Loose",      Class),
    Class = ifelse(Class == "11.1", "1",          Class),
    Class = ifelse(Class == "11.2", "2",          Class),
    Class = ifelse(Class == "11.3", "3",          Class),
    Class = ifelse(Class == "11.4", "4",          Class),
    Class = ifelse(Class == "11.5", ">4",         Class),
    Class = ifelse(Class == "11.6", "N/A",        Class),
    Class = ifelse(Class == "5.1",  "No bulge",   Class),
    Class = ifelse(Class == "5.2",  "Barely",     Class),
    Class = ifelse(Class == "5.3",  "Clearly",    Class),
    Class = ifelse(Class == "5.4",  "Dominant",   Class),
    Class = ifelse(Class == "1.3",  "",           Class)
  )

classes_oddities <- classes %>%
  filter(
    Class != "6.1",
    Question %in% c("6", "8")
  ) %>%
  mutate(
    Oddity = ifelse(Class == "6.2", "Not odd", ""),
    Oddity = ifelse(Class != "6.2", "Odd",     Oddity)
  ) %>%
  mutate(
    Class = ifelse(Class == "8.1", "Ring",      Class),
    Class = ifelse(Class == "8.2", "Lens",      Class),
    Class = ifelse(Class == "8.3", "Disturbed", Class),
    Class = ifelse(Class == "8.4", "Irregular", Class),
    Class = ifelse(Class == "8.5", "Other",     Class),
    Class = ifelse(Class == "8.6", "Merger",    Class),
    Class = ifelse(Class == "8.7", "Dust lane", Class),
    Class = ifelse(Class == "6.2", "",          Class)
  )

# --------------------------------------------------------------------------- #
# Save dataframes
# --------------------------------------------------------------------------- #

write.csv(filter(classes_detail, Question != "Artifact"),
          "../Research/Data/galaxy_classifications.csv", row.names = FALSE)

write.csv(select(filter(classes_oddities, Question == "8"), -Question),
          "../Research/Data/galaxy_oddities.csv", row.names = FALSE)

# --------------------------------------------------------------------------- #
# Plot: overall galaxy classification treemap
# --------------------------------------------------------------------------- #

png(paste0(FIGURES, "classification_treemap_all.png"), width = 1200, height = 800)
classes_detail %>%
  treemap(
    index        = c("Galaxy", "Question", "Class"),
    vSize        = "Probability",
    title        = "Galaxy Classifications",
    fontsize.title = 20,
    type         = "index",
    palette      = "Pastel2",
    fontcolor.labels = c("black", "black", "white"),
    overlap.labels   = 0.1,
    bg.labels        = c("transparent"),
    align.labels = list(
      c("center", "center"),
      c("center", "center"),
      c("right",  "bottom")
    )
  )
dev.off()

# --------------------------------------------------------------------------- #
# Plot: overall disk classification treemap
# --------------------------------------------------------------------------- #

png(paste0(FIGURES, "classification_treemap_disk.png"), width = 1200, height = 800)
classes_detail %>%
  filter(Galaxy == "Disk") %>%
  treemap(
    index        = c("Question", "Class"),
    vSize        = "Probability",
    title        = "Disk Classifications",
    fontsize.title = 20,
    type         = "index",
    palette      = "Pastel1",
    fontcolor.labels = c("black", "black", "white"),
    overlap.labels   = 0.1,
    bg.labels        = c("transparent"),
    align.labels = list(
      c("center", "center"),
      c("center", "center"),
      c("right",  "bottom")
    )
  )
dev.off()

# --------------------------------------------------------------------------- #
# Plot: overall oddities treemap
# --------------------------------------------------------------------------- #

png(paste0(FIGURES, "oddities_treemap_all.png"), width = 1200, height = 800)
classes_oddities %>%
  treemap(
    index        = c("Oddity", "Class"),
    vSize        = "Probability",
    title        = "Galaxy Oddities",
    fontsize.title = 20,
    type         = "index",
    palette      = "Pastel1",
    fontcolor.labels = c("black", "white"),
    overlap.labels   = 0.1,
    bg.labels        = c("transparent"),
    align.labels = list(
      c("center", "center"),
      c("right",  "bottom")
    )
  )
dev.off()

# --------------------------------------------------------------------------- #
# Plot: Galaxy 161279 – disk subclassifications
# --------------------------------------------------------------------------- #

png(paste0(FIGURES, "galaxy_161279_classifications.png"), width = 1200, height = 800)
classes_detail %>%
  filter(GalaxyID == 161279, Galaxy == "Disk", Probability > 0) %>%
  treemap(
    index        = c("Question", "Class"),
    vSize        = "Probability",
    title        = "Galaxy 161279 Classifications",
    fontsize.title = 20,
    type         = "index",
    palette      = "Pastel2",
    fontcolor.labels = c("black", "white"),
    overlap.labels   = 0.1,
    bg.labels        = c("transparent"),
    align.labels = list(
      c("center", "center"),
      c("right",  "bottom")
    )
  )
dev.off()

# --------------------------------------------------------------------------- #
# Plot: Galaxy 161279 – oddities
# --------------------------------------------------------------------------- #

png(paste0(FIGURES, "galaxy_161279_oddities.png"), width = 1200, height = 800)
classes_oddities %>%
  filter(GalaxyID == 161279, Probability > 0) %>%
  treemap(
    index        = c("Oddity", "Class"),
    vSize        = "Probability",
    title        = "Galaxy 161279 Oddities",
    fontsize.title = 20,
    type         = "index",
    palette      = "Pastel1",
    fontcolor.labels = c("black", "white"),
    overlap.labels   = 0.1,
    bg.labels        = c("transparent"),
    align.labels = list(
      c("center", "center"),
      c("right",  "bottom")
    )
  )
dev.off()

message("Done. PNGs written to ", FIGURES)
message("CSVs written to ../Research/Data/")
