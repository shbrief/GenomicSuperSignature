## The universe of MeSH terms
# githubURL <- "https://github.com/shbrief/GenomicSuperSignaturePaper/raw/master/inst/extdata/MeSH_terms_1399refinebio.rda"
# load(url(githubURL))
# x <- mesh_table$name
# length(x)   # 18,798 MeSH terms
# length(unique(x))   # 4,430 unique MeSH terms

## Top 9 most common MeSH terms in RAVmodel_536
## Because weighting process of drawWordcloud handles these terms, we didn't
## add them in the droplist.
general = c("Humans", "RNA-Seq", "Transcriptome", "RNA",
            "Whole Exome Sequencing", "Sequencing, RNA", "Animals", "mRNA",
            "Base Sequence")

## Non-specific, less-informative terms appeared more than once
ns = c("Goals", "Phenotype", "DNA", "Phenotypes", "Deep Sequencing",
       "Sequence Analysis, RNA", "Mutation", "Mutations", "Gene Expression",
       "Computational Biology", "Genome", "Records", "Data Analysis",
       "Sequencing, Whole Genome", "Research Subjects", "Publication",
       "Publishing", "Cues", "Filing", "Survey", "Workflow", "Prevalence",
       "Reagents", "Biases", "Data Files", "Gold", "Genotypes", "Nutrient",
       "Proteins", "Respect", "Seasons", "Solvents", "Survival Rate",
       "Data Accuracy", "High-Throughput RNA Sequencing", "Human Body",
       "ROC Curve", "Sequence Alignment", "Sequence Analysis", "Universities",
       "Intelligence", "Internet", "Islands", "Registry", "Running", "Silver",
       "Steel", "Biological Products", "Pilot Projects", "Research Design")

## Non-specific, less-informative terms appeared only once
ns_only1 = c("Access to Information", "Bar Codes", "Benchmark", "Big Data",
             "Biologics", "Biology", "Cataloging", "Checklist",
             "DNA Sequencing", "Drugs", "Dry Ice", "Dyes", "Electricity",
             "Electrode", "Electrodes", "Emigrant","Emigrants and Immigrants",
             "Emigrations", "Holidays", "Incubator", "Individuality",
             "Inpatients", "Intensive Care Units", "Intuition", "Logic",
             "Maintenance", "Motivation", "Personality", "Privacy",
             "Product Labeling", "PubMed", "Stainings", "Surveys", "Walking")

## Default drop list for drawWordcloud function
droplist = c(ns, ns_only1)
save(droplist, file = "data/droplist.RData")
