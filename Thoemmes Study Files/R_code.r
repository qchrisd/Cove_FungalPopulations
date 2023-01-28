#load packages
library(ggplot2)
library(grid)
library(vegan)
library(mctoolsr)
library(PMCMRplus)
library(stats)
library(dplyr)
library(FSA)

#import filtered dataset

input_woodrat = load_taxa_table("otu-table-w-taxonomy-deblur_edit.txt", "WoodratMap_MST.txt")

sort <- sort(colSums(input_woodrat$data_loaded))

print(sort)


#rarefy

woodrat_rar = single_rarefy(input_woodrat, 4000)

sort(colSums(woodrat_rar$data_loaded))


#export rarefied dataset

export_taxa_table(woodrat_rar, "OTU_woodrat_rar.txt")


#import rarefied dataset

woodrat_rar = load_taxa_table("OTU_woodrat_rar.txt", "WoodratMap_MST.txt")


#filtered datasets

input_wr_bodies = filter_data(woodrat_rar, 'Location', keep_vals = 'Key Largo Woodrat')

input_wr_nests = filter_data(woodrat_rar, 'Type', keep_vals = 'InNest')

input_wr_envir = filter_data(woodrat_rar, 'Location', keep_vals = c('Supplemental In','Natural In','Forest Floor'))

input_wr_all = filter_data(woodrat_rar, 'Location', filter_vals = c('NA','Key Largo Cotton Mouse','Black Rat'))

input_wr_bod_nests = filter_data(woodrat_rar, 'Location', keep_vals = c('Key Largo Woodrat','Natural In','Supplemental In'))


#all locations and bodies

dm_wr_all = calc_dm(input_wr_all$data_loaded)

ord_wr_all = calc_ordination(dm_wr_all, 'nmds')

plot_ordination(input_wr_all, ord_wr_all, 'Location', hulls = FALSE)


#nests

dm_wr_nests = calc_dm(input_wr_nests$data_loaded)

ord_wr_nests = calc_ordination(dm_wr_nests, 'nmds')

plot_ordination(input_wr_nests, ord_wr_nests, 'Location', hulls = FALSE)

ado_nests <- adonis(dm_wr_nests ~ Location, data=input_wr_nests$map_loaded, permutations=999, method="euclidean")

print(ado_nests)


#nests and forest

dm_wr_envir = calc_dm(input_wr_envir$data_loaded)

ord_wr_envir = calc_ordination(dm_wr_envir, 'nmds')

plot_ordination(input_wr_envir, ord_wr_envir, 'Type', hulls = FALSE)

ado_envir <- adonis(dm_wr_envir ~ Type, data=input_wr_envir$map_loaded, permutations=999, method="euclidean")

print(ado_envir)


#bodies and nests

dm_wr_bod_nests = calc_dm(input_wr_bod_nests$data_loaded)

ord_wr_bod_nests = calc_ordination(dm_wr_bod_nests, 'nmds')

plot_ordination(input_wr_bod_nests, ord_wr_bod_nests, 'Type', hulls = FALSE)

ado_bod_nests <- adonis(dm_wr_bod_nests ~ Type, data=input_wr_bod_nests$map_loaded, permutations=999, method="euclidean")

print(ado_bod_nests)


#diversity in nests

richness_wr_nests <- calc_diversity(input_wr_nests$data_loaded, metric = "richness")

rich.map_wr_nests <- merge(richness_wr_nests, input_wr_nests$map_loaded, by=0, all=TRUE)

rich_nests <- aggregate(rich.map_wr_nests$`x`~Location, rich.map_wr_nests, mean)

print(rich_nests)

plot_diversity(input_wr_nests, variable = "Location", metric = "richness")

kruskal.test(`x`~Location, data=rich.map_wr_nests)


#diversity for all samples

richness_wr_all <- calc_diversity(input_wr_all$data_loaded, metric = "richness")

rich.map_wr_all <- merge(richness_wr_all, input_wr_all$map_loaded, by=0, all=TRUE)

rich_all <- aggregate(rich.map_wr_all$`x`~Location, rich.map_wr_all, mean)

print(rich_all)

kruskal.test(`x`~Location, data=rich.map_wr_all)


#diversity for natural nests and forests

input_wr_nat_for = filter_data(input_wr_envir, 'Location', filter_vals = 'Supplemental In')

richness_wr_nat_for <- calc_diversity(input_wr_nat_for$data_loaded, metric = "richness")

rich.map_wr_nat_for <- merge(richness_wr_nat_for, input_wr_nat_for$map_loaded, by=0, all=TRUE)

rich_nat_for <- aggregate(rich.map_wr_nat_for$`x`~Location, rich.map_wr_nat_for, mean)

print(rich_nat_for)

kruskal.test(`x`~Location, data=rich.map_wr_nat_for)


#shannon diversity in nests

shannon_wr_nests <- calc_diversity(input_wr_nests$data_loaded, metric = "shannon")

shan.map_wr_nests <- merge(shannon_wr_nests, input_wr_nests$map_loaded, by=0, all=TRUE)

plot_diversity(input_wr_nests, variable = "Location", metric = "shannon")

kruskal.test(`x`~Location, data=shan.map_wr_nests)


#diversity for nests and forests

richness_wr_envir <- calc_diversity(input_wr_envir$data_loaded, metric = "richness")

rich.map_wr_envir <- merge(richness_wr_envir, input_wr_envir$map_loaded, by=0, all=TRUE)

rich_envir <- aggregate(rich.map_wr_envir$`x`~Location, rich.map_wr_envir, mean)

print(rich_envir)

kruskal.test(`x`~Location, data=rich.map_wr_envir)

plot_diversity(input_wr_envir, variable = "Location", metric = "richness")


#shannon diversity in nests

shannon_wr_envir <- calc_diversity(input_wr_envir$data_loaded, metric = "shannon")

shan.map_wr_envir <- merge(shannon_wr_envir, input_wr_envir$map_loaded, by=0, all=TRUE)

plot_diversity(input_wr_envir, variable = "Location", metric = "shannon")

kruskal.test(`x`~Location, data=shan.map_wr_envir)


#shannon diversity for natural nests and forests

shannon_wr_nat_for <- calc_diversity(input_wr_nat_for$data_loaded, metric = "shannon")

shan.map_wr_nat_for <- merge(shannon_wr_nat_for, input_wr_nat_for$map_loaded, by=0, all=TRUE)

shan_nat_for <- aggregate(shan.map_wr_nat_for$`x`~Location, shan.map_wr_nat_for, mean)

print(shan_nat_for)

kruskal.test(`x`~Location, data=shan.map_wr_nat_for)


#taxa summary phyla

tax_sum_phyla_nests = summarize_taxonomy(input_wr_nests, level = 2, report_higher_tax = TRUE)

sum_phy_nests <- taxa_summary_by_sample_type(taxa_smry_df = tax_sum_phyla_nests, 'Location', metadata_map = input_wr_nests$map_loaded, smry_fun = sum)

print(sum_phy_nests)

tax_sum_phyla_bodies = summarize_taxonomy(input_wr_bodies, level = 2, report_higher_tax = TRUE)

sum_phy_bodies <- taxa_summary_by_sample_type(taxa_smry_df = tax_sum_phyla_bodies, 'Location', metadata_map = input_wr_bodies$map_loaded, smry_fun = sum)

print(sum_phy_bodies)

sum_phy_bod <- t(sum_phy_bodies)


#taxa summary families

tax_sum_family_nests = summarize_taxonomy(input_wr_nests, level = 5, report_higher_tax = TRUE)

sum_fam_nests <- taxa_summary_by_sample_type(taxa_smry_df = tax_sum_family_nests, 'Location', metadata_map = input_wr_nests$map_loaded, smry_fun = mean)

tax_fam_nests_plot <- plot_taxa_bars(tax_sum_family_nests, input_wr_nests$map_loaded, type_header = "Location", num_taxa = 10)

plot(tax_fam_nests_plot)

tax_sum_family_bodies = summarize_taxonomy(input_wr_bodies, level = 5, report_higher_tax = TRUE)

sum_fam_bodies <- taxa_summary_by_sample_type(taxa_smry_df = tax_sum_family_bodies, 'Location', metadata_map = input_wr_bodies$map_loaded, smry_fun = mean)


#taxa summary genera - all

tax_sum_genera_all = summarize_taxonomy(input_wr_all, level = 6, report_higher_tax = TRUE)

sum_gen_all <- taxa_summary_by_sample_type(taxa_smry_df = tax_sum_genera_all, 'Location', metadata_map = input_wr_all$map_loaded, smry_fun = mean)

tax_gen_all_plot <- plot_taxa_bars(tax_sum_genera_all, input_wr_all$map_loaded, type_header = "Location", num_taxa = 20)

genera.map <- merge(sum_gen_nests, input_wr_nests$map_loaded, by=0, all=TRUE)

aggregate(genera.map$`Natural In`~Location, genera.map, mean)

tax_sum_genera_bodies = summarize_taxonomy(input_wr_bodies, level = 6, report_higher_tax = TRUE)

sum_gen_bodies <- taxa_summary_by_sample_type(taxa_smry_df = tax_sum_genera_bodies, 'Location', metadata_map = input_wr_bodies$map_loaded, smry_fun = mean)


#taxa summary genera

tax_sum_genera_nests = summarize_taxonomy(input_wr_nests, level = 6, report_higher_tax = TRUE)

sum_gen_nests <- taxa_summary_by_sample_type(taxa_smry_df = tax_sum_genera_nests, 'Location', metadata_map = input_wr_nests$map_loaded, smry_fun = mean)

tax_gen_nests_plot <- plot_taxa_bars(tax_sum_genera_nests, input_wr_nests$map_loaded, type_header = "Location", num_taxa = 10)

plot(tax_gen_nests_plot)

genera.map <- merge(sum_gen_nests, input_wr_nests$map_loaded, by=0, all=TRUE)

aggregate(genera.map$`Natural In`~Location, genera.map, mean)

tax_sum_genera_bodies = summarize_taxonomy(input_wr_bodies, level = 6, report_higher_tax = TRUE)

sum_gen_bodies <- taxa_summary_by_sample_type(taxa_smry_df = tax_sum_genera_bodies, 'Location', metadata_map = input_wr_bodies$map_loaded, smry_fun = mean)


#taxa summary species

tax_sum_species_nests = summarize_taxonomy(input_wr_nests, level = 7, report_higher_tax = TRUE)

sum_spp_nests <- taxa_summary_by_sample_type(taxa_smry_df = tax_sum_species_nests, 'Location', metadata_map = input_wr_nests$map_loaded, smry_fun = mean)

tax_spp_nests_plot <- plot_taxa_bars(tax_sum_species_nests, input_wr_nests$map_loaded, type_header = "Location", num_taxa = 20)

plot(tax_spp_nests_plot)


#Pseudonocardiaceae in nests

pseudo.input <- filter_taxa_from_input(input_wr_nests, taxa_to_keep = 'f__Pseudonocardiaceae')

pseudo.tot<- as.data.frame(colSums(pseudo.input$data_loaded))

names(pseudo.tot)[names(pseudo.tot)=="colSums(pseudo.input$data_loaded)"]

pseudo.per <- ((pseudo.tot)/4000)*100

pseudo.per.map <- merge(pseudo.per, input_wr_nests$map_loaded, by=0, all=TRUE)

aggregate(pseudo.per.map$`colSums(pseudo.input$data_loaded)`~Location, pseudo.per.map, mean)


#Pseudo on bodies

pseudo_bodies.input <- filter_taxa_from_input(input_wr_bodies, taxa_to_keep = 'f__Pseudonocardiaceae')

pseudo_bod.tot<- as.data.frame(colSums(pseudo_bodies.input$data_loaded))

names(pseudo_bod.tot)[names(pseudo_bod.tot)=="colSums(pseudo_bodies.input$data_loaded)"]

pseudo_bod.per <- ((pseudo_bod.tot)/4000)*100

pseudo_bod.per.map <- merge(pseudo_bod.per, input_wr_bodies$map_loaded, by=0, all=TRUE)

aggregate(pseudo_bod.per.map$`colSums(pseudo_bodies.input$data_loaded)`~Location, pseudo_bod.per.map, mean)


#Strep on bodies

strep_bodies.input <- filter_taxa_from_input(input_wr_bodies, taxa_to_keep = 'f__Streptomycetaceae')

strep_bod.tot<- as.data.frame(colSums(strep_bodies.input$data_loaded))

names(strep_bod.tot)[names(strep_bod.tot)=="colSums(strep_bodies.input$data_loaded)"]

strep_bod.per <- ((strep_bod.tot)/4000)*100

strep_bod.per.map <- merge(strep_bod.per, input_wr_bodies$map_loaded, by=0, all=TRUE)

aggregate(strep_bod.per.map$`colSums(strep_bodies.input$data_loaded)`~Location, strep_bod.per.map, mean)


#Pseudonocardiaceae and Streptomycetaceae in nest and forest samples

pseudo_strep.input <- filter_taxa_from_input(input_wr_envir, taxa_to_keep = c('f__Pseudonocardiaceae','f__Streptomycetaceae'))

pseudo_strep.tot<- as.data.frame(colSums(pseudo_strep.input$data_loaded))

names(pseudo_strep.tot)[names(pseudo_strep.tot)=="colSums(pseudo_strep.input$data_loaded)"]

pseudo_strep.per <- ((pseudo_strep.tot)/4000)*100

pseudo_strep.per.map <- merge(pseudo_strep.per, input_wr_envir$map_loaded, by=0, all=TRUE)

aggregate(pseudo_strep.per.map$`colSums(pseudo_strep.input$data_loaded)`~Location, pseudo_strep.per.map, mean)

kruskal.test(`colSums(pseudo_strep.input$data_loaded)`~Location, data=pseudo_strep.per.map)


#Pseudonocardiaceae and Streptomycetaceae on bodies

pseudo_strep_bodies.input <- filter_taxa_from_input(input_wr_bodies, taxa_to_keep = c('f__Pseudonocardiaceae','f__Streptomycetaceae'))

pseudo_strep_bodies.tot<- as.data.frame(colSums(pseudo_strep_bodies.input$data_loaded))

names(pseudo_strep_bodies.tot)[names(pseudo_strep_bodies.tot)=="colSums(pseudo_strep_bodies.input$data_loaded)"]

pseudo_strep_bodies.per <- ((pseudo_strep_bodies.tot)/4000)*100

pseudo_strep_bodies.per.map <- merge(pseudo_strep_bodies.per, input_wr_bodies$map_loaded, by=0, all=TRUE)

aggregate(pseudo_strep_bodies.per.map$`colSums(pseudo_strep_bodies.input$data_loaded)`~Type, pseudo_strep_bodies.per.map, mean)

kruskal.test(`colSums(pseudo_strep.input$data_loaded)`~Type, data=pseudo_strep.per.map)


#pathogenic bacteria nests - genera

pathogen.input_nests_gen <- filter_taxa_from_input(input_wr_nests, taxa_to_keep = c('g__Bacillus|g__Bartonella|g__Bordetella|g__Brucella|g__Campylobacter|g__CAR|g__Chlamydia|g__Citrobacter|g__Clostridium|g__Corynebacterium|g__Coxiella|g__Ehrlichia|g__Elmeria|g__Escherichia|g__Flexispira|g__Francisella|g__Helicobacter|g__Kleibsiella|g__Lawsonia|g__Mycobacterium|g__Mycoplasma|g__Pasteurella|g__Pneumocystis|g__Proteus|g__Pseudomonas|g__Salmonella|g__Spirillum|g__Staphylococcus|g__Streptobacillus|g__Streptococcus|g__Treponema|g__Yersinia|g__Enterococcus|g__Burkholderia|g__Spiroplasma|g__Pasteurella|g__Anaplasma|g__Leptospira'))

pathogen.tot_nests_gen <- as.data.frame(colSums(pathogen.input_nests_gen$data_loaded))

names(pathogen.tot_nests_gen)[names(pathogen.tot_nests_gen)=="colSums(pathogen.input_nests_gen$data_loaded)"]

pathogen.per_nests_gen <- ((pathogen.tot_nests_gen)/4000)*100

pathogen.per.map_nests_gen <- merge(pathogen.per_nests_gen, input_wr_nests$map_loaded, by=0, all=TRUE)

aggregate(pathogen.per.map_nests_gen$`colSums(pathogen.input_nests_gen$data_loaded)`~Location, pathogen.per.map_nests_gen, mean)

kruskal.test(`colSums(pathogen.input_nests_gen$data_loaded)`~NestType, data=pathogen.per.map_nests_gen)


#pathogenic bacteria nests - species

pathogen.input_nests <- filter_taxa_from_input(pathogen.input_nests_gen, taxa_to_keep = c('s__anthracis|g__Bartonella|s__bronchiseptica|s__tularensis|s__hinzii|g__Brucella|g__Campylobacter|s__bacillus|s__trachomatis|s__muridarum|s__psittaci|s__rodentium|s__piliforme|s__difficile|s__perfringens|s__kutscheri|s__pseudodiptheriticum|s__burnetii|s__chaffeenis|s__muris|g__Elmeria|s__coli|s__rappini|s__tularensis|s__hepaticus|s__bilis|s__muridarum|s__trogontum|s__rodentium|s__typhlonius|s__pneumoniae|s__oxytoca|s__intracellularis|s__avium-intracellulare|s__lepraemurium|s__tuberculosis|s__microti|s__bovis|s__africanum|s__canetti|s__pulmonis|s__neurolyticum|s__hemomuris|s__coccoides|s__pneumotropica|s__murina|s__mirabilis|s__aeruginosa|s__enteritidis|s__enterica|s__minus|s__aureus|s__moniliformis|s__pneumoniae|s__dysgalactiae|g__Treponema|s__pestis|s__faecalis|s__faecium|s__pseudomallei|s__mallei|s__thailandensis|s__oklahomensis|s__humptydooensis|s__mirum|s__pneumotropica|s__phagocytophilum|g__Leptospira'))

pathogen.tot_nests <- as.data.frame(colSums(pathogen.input_nests$data_loaded))

names(pathogen.tot_nests)[names(pathogen.tot_nests)=="colSums(pathogen.input_nests$data_loaded)"]

pathogen.per_nests <- ((pathogen.tot_nests)/4000)*100

pathogen.per.map_nests <- merge(pathogen.per_nests, input_wr_nests$map_loaded, by=0, all=TRUE)

aggregate(pathogen.per.map_nests$`colSums(pathogen.input_nests$data_loaded)`~Location, pathogen.per.map_nests, mean)

kruskal.test(`colSums(pathogen.input_nests$data_loaded)`~NestType, data=pathogen.per.map_nests)

View(pathogen.input_nests$taxonomy_loaded)


#pathogenic bacteria in nest and forest samples - genera

pathogen.input_envir_gen <- filter_taxa_from_input(input_wr_envir, taxa_to_keep = c('g__Bacillus|g__Bartonella|g__Bordetella|g__Brucella|g__Campylobacter|g__CAR|g__Chlamydia|g__Citrobacter|g__Clostridium|g__Corynebacterium|g__Coxiella|g__Ehrlichia|g__Elmeria|g__Escherichia|g__Flexispira|g__Francisella|g__Helicobacter|g__Kleibsiella|g__Lawsonia|g__Mycobacterium|g__Mycoplasma|g__Pasteurella|g__Pneumocystis|g__Proteus|g__Pseudomonas|g__Salmonella|g__Spirillum|g__Staphylococcus|g__Streptobacillus|g__Streptococcus|g__Treponema|g__Yersinia|g__Enterococcus|g__Burkholderia|g__Spiroplasma|g__Pasteurella|g__Anaplasma|g__Leptospira'))

pathogen.tot_envir_gen <- as.data.frame(colSums(pathogen.input_envir_gen$data_loaded))

names(pathogen.tot_envir_gen)[names(pathogen.tot_envir_gen)=="colSums(pathogen.input_envir_gen$data_loaded)"]

pathogen.per_envir_gen <- ((pathogen.tot_envir_gen)/4000)*100

pathogen.per.map_envir_gen <- merge(pathogen.per_envir_gen, input_wr_envir$map_loaded, by=0, all=TRUE)

aggregate(pathogen.per.map_envir_gen$`colSums(pathogen.input_envir_gen$data_loaded)`~Location, pathogen.per.map_envir_gen, mean)

kruskal.test(`colSums(pathogen.input_envir_gen$data_loaded)`~Location, pathogen.per.map_envir_gen)


#pathogenic bacteria in nest and forest samples - species

pathogen.input_envir <- filter_taxa_from_input(pathogen.input_envir_gen, taxa_to_keep = c('s__anthracis|g__Bartonella|s__bronchiseptica|s__tularensis|s__hinzii|g__Brucella|g__Campylobacter|s__bacillus|s__trachomatis|s__muridarum|s__psittaci|s__rodentium|s__piliforme|s__difficile|s__perfringens|s__kutscheri|s__pseudodiptheriticum|s__burnetii|s__chaffeenis|s__muris|g__Elmeria|s__coli|s__rappini|s__tularensis|s__hepaticus|s__bilis|s__muridarum|s__trogontum|s__rodentium|s__typhlonius|s__pneumoniae|s__oxytoca|s__intracellularis|s__avium-intracellulare|s__lepraemurium|s__tuberculosis|s__microti|s__bovis|s__africanum|s__canetti|s__pulmonis|s__neurolyticum|s__hemomuris|s__coccoides|s__pneumotropica|s__murina|s__mirabilis|s__aeruginosa|s__enteritidis|s__enterica|s__minus|s__aureus|s__moniliformis|s__pneumoniae|s__dysgalactiae|g__Treponema|s__pestis|s__faecalis|s__faecium|s__pseudomallei|s__mallei|s__thailandensis|s__oklahomensis|s__humptydooensis|s__mirum|s__pneumotropica|s__phagocytophilum|g__Leptospira'))

pathogen.tot_envir <- as.data.frame(colSums(pathogen.input_envir$data_loaded))

names(pathogen.tot_envir)[names(pathogen.tot_envir)=="colSums(pathogen.input_envir$data_loaded)"]

pathogen.per_envir <- ((pathogen.tot_envir)/4000)*100

pathogen.per.map_envir <- merge(pathogen.per_envir, input_wr_envir$map_loaded, by=0, all=TRUE)

aggregate(pathogen.per.map_envir$`colSums(pathogen.input_envir$data_loaded)`~Location, pathogen.per.map_envir, mean)

kruskal.test(`colSums(pathogen.input_envir$data_loaded)`~Location, pathogen.per.map_envir)

View(pathogen.input_envir$taxonomy_loaded)


#pathogenic bacteria on woodrats - genera

pathogen.input_bod_gen <- filter_taxa_from_input(input_wr_bodies, taxa_to_keep = c('g__Bacillus|g__Bartonella|g__Bordetella|g__Brucella|g__Campylobacter|g__CAR|g__Chlamydia|g__Citrobacter|g__Clostridium|g__Corynebacterium|g__Coxiella|g__Ehrlichia|g__Elmeria|g__Escherichia|g__Flexispira|g__Francisella|g__Helicobacter|g__Kleibsiella|g__Lawsonia|g__Mycobacterium|g__Mycoplasma|g__Pasteurella|g__Pneumocystis|g__Proteus|g__Pseudomonas|g__Salmonella|g__Spirillum|g__Staphylococcus|g__Streptobacillus|g__Streptococcus|g__Treponema|g__Yersinia|g__Enterococcus|g__Burkholderia|g__Spiroplasma|g__Pasteurella|g__Anaplasma|g__Leptospira'))

pathogen.tot_bod_gen <- as.data.frame(colSums(pathogen.input_bod_gen$data_loaded))

names(pathogen.tot_bod_gen)[names(pathogen.tot_bod_gen)=="colSums(pathogen.input_bod_gen$data_loaded)"]

pathogen.per_bod_gen <- ((pathogen.tot_bod_gen)/4000)*100

pathogen.per.map_bod_gen <- merge(pathogen.per_bod_gen, input_wr_bodies$map_loaded, by=0, all=TRUE)

aggregate(pathogen.per.map_bod_gen$`colSums(pathogen.input_bod_gen$data_loaded)`~Location, pathogen.per.map_bod_gen, mean)


#pathogenic bacteria on woodrats - species

pathogen.input_bod <- filter_taxa_from_input(pathogen.input_bod_gen, taxa_to_keep = c('s__anthracis|g__Bartonella|s__bronchiseptica|s__tularensis|s__hinzii|g__Brucella|g__Campylobacter|s__bacillus|s__trachomatis|s__muridarum|s__psittaci|s__rodentium|s__piliforme|s__difficile|s__perfringens|s__kutscheri|s__pseudodiptheriticum|s__burnetii|s__chaffeenis|s__muris|g__Elmeria|s__coli|s__rappini|s__tularensis|s__hepaticus|s__bilis|s__muridarum|s__trogontum|s__rodentium|s__typhlonius|s__pneumoniae|s__oxytoca|s__intracellularis|s__avium-intracellulare|s__lepraemurium|s__tuberculosis|s__microti|s__bovis|s__africanum|s__canetti|s__pulmonis|s__neurolyticum|s__hemomuris|s__coccoides|s__pneumotropica|s__murina|s__mirabilis|s__aeruginosa|s__enteritidis|s__enterica|s__minus|s__aureus|s__moniliformis|s__pneumoniae|s__dysgalactiae|g__Treponema|s__pestis|s__faecalis|s__faecium|s__pseudomallei|s__mallei|s__thailandensis|s__oklahomensis|s__humptydooensis|s__mirum|s__pneumotropica|s__phagocytophilum|g__Leptospira'))

pathogen.tot_bod <- as.data.frame(colSums(pathogen.input_bod$data_loaded))

names(pathogen.tot_bod)[names(pathogen.tot_bod)=="colSums(pathogen.input_bod$data_loaded)"]

pathogen.per_bod <- ((pathogen.tot_bod)/4000)*100

pathogen.per.map_bod <- merge(pathogen.per_bod, input_wr_bodies$map_loaded, by=0, all=TRUE)

aggregate(pathogen.per.map_bod$`colSums(pathogen.input_bod$data_loaded)`~Location, pathogen.per.map_bod, mean)

View(pathogen.input_bod$taxonomy_loaded)


#pathogenic bacteria in nests and on woodrats - genera

pathogen.input_bod_nests_gen <- filter_taxa_from_input(input_wr_bod_nests, taxa_to_keep = c('g__Bacillus|g__Bartonella|g__Bordetella|g__Brucella|g__Campylobacter|g__CAR|g__Chlamydia|g__Citrobacter|g__Clostridium|g__Corynebacterium|g__Coxiella|g__Ehrlichia|g__Elmeria|g__Escherichia|g__Flexispira|g__Francisella|g__Helicobacter|g__Kleibsiella|g__Lawsonia|g__Mycobacterium|g__Mycoplasma|g__Pasteurella|g__Pneumocystis|g__Proteus|g__Pseudomonas|g__Salmonella|g__Spirillum|g__Staphylococcus|g__Streptobacillus|g__Streptococcus|g__Treponema|g__Yersinia|g__Enterococcus|g__Burkholderia|g__Spiroplasma|g__Pasteurella|g__Anaplasma|g__Leptospira'))

pathogen.tot_bod_nests_gen <- as.data.frame(colSums(pathogen.input_bod_nests_gen$data_loaded))

names(pathogen.tot_bod_nests_gen)[names(pathogen.tot_bod_nests_gen)=="colSums(pathogen.input_bod_nests_gen$data_loaded)"]

pathogen.per_bod_nests_gen <- ((pathogen.tot_bod_nests_gen)/4000)*100

pathogen.per.map_bod_nests_gen <- merge(pathogen.per_bod_nests_gen, input_wr_bod_nests$map_loaded, by=0, all=TRUE)

aggregate(pathogen.per.map_bod_nests_gen$`colSums(pathogen.input_bod_nests_gen$data_loaded)`~Location, pathogen.per.map_bod_nests_gen, mean)


#pathogenic bacteria in nests and on woodrats - species

pathogen.input_bod_nests <- filter_taxa_from_input(pathogen.input_bod_nests_gen, taxa_to_keep = c('s__anthracis|g__Bartonella|s__bronchiseptica|s__tularensis|s__hinzii|g__Brucella|g__Campylobacter|s__bacillus|s__trachomatis|s__muridarum|s__psittaci|s__rodentium|s__piliforme|s__difficile|s__perfringens|s__kutscheri|s__pseudodiptheriticum|s__burnetii|s__chaffeenis|s__muris|g__Elmeria|s__coli|s__rappini|s__tularensis|s__hepaticus|s__bilis|s__muridarum|s__trogontum|s__rodentium|s__typhlonius|s__pneumoniae|s__oxytoca|s__intracellularis|s__avium-intracellulare|s__lepraemurium|s__tuberculosis|s__microti|s__bovis|s__africanum|s__canetti|s__pulmonis|s__neurolyticum|s__hemomuris|s__coccoides|s__pneumotropica|s__murina|s__mirabilis|s__aeruginosa|s__enteritidis|s__enterica|s__minus|s__aureus|s__moniliformis|s__pneumoniae|s__dysgalactiae|g__Treponema|s__pestis|s__faecalis|s__faecium|s__pseudomallei|s__mallei|s__thailandensis|s__oklahomensis|s__humptydooensis|s__mirum|s__pneumotropica|s__phagocytophilum|g__Leptospira'))

pathogen.tot_bod_nests <- as.data.frame(colSums(pathogen.input_bod_nests$data_loaded))

names(pathogen.tot_bod_nests)[names(pathogen.tot_bod_nests)=="colSums(pathogen.input_bod_nests$data_loaded)"]

pathogen.per_bod_nests <- ((pathogen.tot_bod_nests)/4000)*100

pathogen.per.map_bod_nests <- merge(pathogen.per_bod_nests, input_wr_bod_nests$map_loaded, by=0, all=TRUE)

aggregate(pathogen.per.map_bod_nests$`colSums(pathogen.input_bod_nests$data_loaded)`~Location, pathogen.per.map_bod_nests, mean)

kruskal.test(`colSums(pathogen.input_bod_nests$data_loaded)`~Type, pathogen.per.map_bod_nests)

View(pathogen.input_bod_nests$taxonomy_loaded)
