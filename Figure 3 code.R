#Figure 3 A 1








# Figure 3 A 3
# PT
library(Seurat)
library(ggplot2)
library(dplyr)
library(tidyr)

df_mouse <- PT_subcl_mouse_res_0.5_dims_10@meta.data %>%
  group_by(ID_to_plot, celltype_updated) %>%
  summarise(count = n(), .groups = 'drop') %>%
  ungroup() %>%
  group_by(ID_to_plot) %>%
  mutate(total = sum(count), 
         percent = count / total) %>%
  ungroup()

ordered_ids <- c("Balbc to Balbc 1", "Balbc to Balbc 2", "Blck6 to Blck6 1",
                 "Balbc to Blck6 1", "Balbc to Blck6 2", 
                 "Blck6 to Balbc 1", "Blck6 to Balbc 2", "Blck6 to Balbc 3")

# Factor the ID with the levels in the order you want
df_mouse$ID_to_plot <- factor(df_mouse$ID_to_plot, levels = ordered_ids)

# Pivot the dataframe for mouse data
df_mouse_long <- df_mouse %>%
  select(ID_to_plot, celltype_updated, percent) %>%
  pivot_longer(
    cols = percent,
    names_to = "name",
    values_to = "value"
  )

# Define custom colors for ID samples
ID_colors <- c("Balbc to Balbc 1" = "grey70", 
               "Balbc to Balbc 2" = "grey50", 
               "Blck6 to Blck6 1" = "grey30",
               "Balbc to Blck6 1" = "darkolivegreen2",
               "Balbc to Blck6 2" = "darkolivegreen3",
               "Blck6 to Balbc 1" = "lightblue3",
               "Blck6 to Balbc 2" = "lightblue4",
               "Blck6 to Balbc 3" = "lightblue")


# Plot with ggplot2 for mouse data with custom colors and ordered ID
ggplot(df_mouse_long, aes(x = celltype_updated, y = value, fill = ID_to_plot)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  scale_fill_manual(values = ID_colors) +
  scale_y_continuous(labels = scales::percent_format()) +
  theme_minimal() +
  theme(
    
    axis.text.y = element_text(size = 13),  # Behalte die y-Achsenbeschriftungen
# Entferne die Legende
  )

# Plot with ggplot2 for mouse data with custom colors and ordered ID
ggplot(df_mouse_long, aes(x = celltype_updated, y = value, fill = ID_to_plot)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  scale_fill_manual(values = ID_colors) +
  scale_y_continuous(labels = scales::percent_format()) +
  theme_minimal() +
  theme(
    axis.text.x = element_blank(),    # Entferne die x-Achsenbeschriftungen
    axis.text.y = element_blank(),  # Behalte die y-Achsenbeschriftungen
    axis.title.x = element_blank(),   # Entferne den Titel der x-Achse
    axis.title.y = element_blank(),   # Entferne den Titel der y-Achse
    legend.position = "none",         # Entferne die Legende
  )


# TAL
library(Seurat)
library(ggplot2)
library(dplyr)
library(tidyr)

# Data preparation for TAL subclustering
df_tal <- Ktx_TAL_mouse_subclustering@meta.data %>%
  group_by(ID_to_plot, celltype_updated) %>%
  summarise(count = n(), .groups = 'drop') %>%
  ungroup() %>%
  group_by(ID_to_plot) %>%
  mutate(total = sum(count), 
         percent = count / total) %>%
  ungroup()

# Define the order for the IDs
ordered_ids_tal <- c("Balbc to Balbc 1", "Balbc to Balbc 2", "Blck6 to Blck6 1",
                     "Balbc to Blck6 1", "Balbc to Blck6 2", 
                     "Blck6 to Balbc 1", "Blck6 to Balbc 2", "Blck6 to Balbc 3")

# Factor the ID with the levels in the order you want
df_tal$ID_to_plot <- factor(df_tal$ID_to_plot, levels = ordered_ids_tal)

# Pivot the dataframe for TAL data
df_tal_long <- df_tal %>%
  select(ID_to_plot, celltype_updated, percent) %>%
  pivot_longer(
    cols = percent,
    names_to = "name",
    values_to = "value"
  )

# Define custom colors for ID samples
ID_colors_tal <- c("Balbc to Balbc 1" = "grey70", 
                   "Balbc to Balbc 2" = "grey50", 
                   "Blck6 to Blck6 1" = "grey30",
                   "Balbc to Blck6 1" = "darkolivegreen2",
                   "Balbc to Blck6 2" = "darkolivegreen3",
                   "Blck6 to Balbc 1" = "lightblue3",
                   "Blck6 to Balbc 2" = "lightblue4",
                   "Blck6 to Balbc 3" = "lightblue")

# Plot with ggplot2 for TAL data with custom colors and ordered ID
ggplot(df_tal_long, aes(x = celltype_updated, y = value, fill = ID_to_plot)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  scale_fill_manual(values = ID_colors_tal) +
  scale_y_continuous(labels = scales::percent_format()) +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 13),  # Behalte die y-Achsenbeschriftungen
  )

# Plot with ggplot2 for TAL data with custom colors and ordered ID without labels
ggplot(df_tal_long, aes(x = celltype_updated, y = value, fill = ID_to_plot)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  scale_fill_manual(values = ID_colors_tal) +
  scale_y_continuous(labels = scales::percent_format()) +
  theme_minimal() +
  theme(
    axis.text.x = element_blank(),    # Entferne die x-Achsenbeschriftungen
    axis.text.y = element_blank(),    # Entferne die y-Achsenbeschriftungen
    axis.title.x = element_blank(),   # Entferne den Titel der x-Achse
    axis.title.y = element_blank(),   # Entferne den Titel der y-Achse
    legend.position = "none"          # Entferne die Legende
  )
