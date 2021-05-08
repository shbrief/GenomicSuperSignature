## Top 9 most common MeSH terms in RAVmodel_536
# Because weighting process of drawWordcloud handles these terms, we didn't add
# them in the droplist. These are most likely showing up in the word cloud representing
# only small scale.
general = c("Humans", "RNA-Seq", "Transcriptome", "RNA",
            "Whole Exome Sequencing", "Sequencing, RNA", "Animals", "mRNA",
            "Base Sequence")

## Non-specific, less-informative terms appeared more than once
nonspecific = c("Goals", "Phenotype", "DNA", "Phenotypes", "Deep Sequencing",
                "Sequence Analysis, RNA", "Mutation", "Mutations", "Gene Expression",
                "Computational Biology", "Genome", "Records", "Data Analysis",
                "Sequencing, Whole Genome", "Research Subjects", "Publication",
                "Publishing", "Cues", "Filing", "Survey", "Workflow", "Prevalence",
                "Reagents", "Biases", "Data Files", "Gold", "Genotypes", "Nutrient",
                "Proteins", "Respect", "Seasons", "Solvents", "Survival Rate",
                "Data Accuracy", "High-Throughput RNA Sequencing", "Human Body",
                "ROC Curve", "Sequence Alignment", "Sequence Analysis",
                "Universities", "Intelligence", "Internet", "Islands", "Registry",
                "Running", "Silver", "Steel")

## Non-specific, less-informative terms appeared only once
nonspecific_only1 = c("Access to Information", "Bar Codes", "Benchmark", "Big Data",
                      "Biologics", "Biology", "Cataloging", "Checklist", "DNA Sequencing",
                      "Drugs", "Dry Ice", "Dyes", "Electricity", "Electrode", "Electrodes",
                      "Emigrant","Emigrants and Immigrants", "Emigrations", "Holidays",
                      "Incubator", "Individuality", "Inpatients", "Intensive Care Units",
                      "Intuition", "Logic", "Maintenance", "Motivation", "Personality",
                      "Privacy", "Product Labeling", "PubMed", "Stainings", "Surveys",
                      "Walking")

# Default drop list for drawWordcloud function
droplist = c(nonspecific, nonspecific_only1)
save(droplist, file = "data/droplist.RData")
