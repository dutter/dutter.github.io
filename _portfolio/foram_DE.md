---
title: "Santa Barbara Basin foraminiferal transcriptomes"
excerpt: From Gomaa et al., 2021 *Science Advances**
permalink: projects/SBB_forams

toc: true
toc_sticky: true
highlight: pygments
---


Overview
--------

This document is meant to provide the code necessary to reproduce the
differential expression analyses performed in the accompanying paper
(Gomaa et al. 2021 Science Advances).

A full description of the sample collection is given in the manuscript.
Briefly, the study focused on two foraminiferal species from the Santa
Barbara Basin, *Nonionella stella* and *Bolivina argentea*. Sediment
samples were taken, from which some *N. stella* and *B. argentea* were
immediately isolated and preserved onboard the ship, while the rest of
the sample was transported to the laboratory for incubation under
varying atmospheric conditions (aerated, hypoxia, anoxia, and
[euxinia](https://en.wikipedia.org/wiki/Euxinia)) with or without
resupplying nitrate and peroxide, both of which are present in their
natural sediments and become depleted over time. From all samples, two
libraries were prepared for transcriptome sequencing -
metatranscriptomes prepared with [NuGEN Trio
kits](https://www.nugen.com/products/trio-rna-seq-library-preparation-kit)
to capture the bulk community transcription of the forams and any
associated bacteria or archaea, and host-enriched transcriptomes using
polyA capture to enrich for eukaryotic mRNAs. Host transcriptomes were
pooled and sequenced on an Illumina NovaSeq (two lanes), and
metatranscriptomes were sequenced on an Illumina HiSeq (two runs).
Transcripts were assembled de novo with Trinity according to the
developer’s recommendations.

Differential expression
-----------------------

The transcript abundances, functional annotations, and RapClust
transcript/cluster associations wereanalyzed in R with
[edgeR](https://bioconductor.org/packages/release/bioc/html/edgeR.html)
using
[limma](https://bioconductor.org/packages/release/bioc/html/limma.html).

The inputs can be found at [this FigShare
DOI.](10.6084/m9.figshare.14183567)

The following libraries were used at various points:

    library(ggplot2)
    library(scales)
    library(reshape2)
    library(vegan)
    library(gtools)
    library(tidyr)

    library(limma)

    ## Warning: package 'limma' was built under R version 3.6.2

    library(edgeR)

    ## Warning: package 'edgeR' was built under R version 3.6.2

First, the various data were read into R, and the separate Interpro and
EggNOG functions were formatted for combination:

    # get functional info
    target <- "NstellaConcat"

    genecalls <- read.csv(paste0(target,'TrinityAll-genecalls.txt'), sep = '\t')
    dups <- read.csv(paste0(target,'TrinityAll-DUPLICATES.txt'), sep = '\t')  ## ONLY FOR NSTELLA 

    rapclust <- read.csv(paste0(target,'_mag.flat.clust'), sep = "\t", header = F)

    eggnog <- read.csv(paste0(target,"_eggnog.emapper.annotations"), sep="\t", header=F, skip = 4,
                       col.names = c("query_name","seed_eggNOG_ortholog","seed_ortholog_evalue","seed_ortholog_score","best_tax_level","Preferred_name","GOs","EC",
                                     "KEGG_ko","KEGG_Pathway","KEGG_Module","KEGG_Reaction","KEGG_rclass", "BRITE","KEGG_TC","CAZy","BiGG_Reaction","tax_scope",
                                     "eggNOG_OGs","bestOG","COG_Functional_Category","eggNOG_free_text_description"))
    eggnog$query_name <- gsub("_.*$","", as.character(eggnog$query_name))
    eggnog <- eggnog[-grep("^#", eggnog$query_name),] # drop rows with comments, usually last few lines

    # drop duplicate annotations for duplicate/identical sequences ID'd by salmon and skipped
    eggnog <- eggnog[!(eggnog$query_name %in% genecalls$gene_callers_id[genecalls$contig %in% dups$DuplicateTxp]),]      ## ONLY FOR NSTELLA

    counts_salmon_raw <- read.csv(paste0(target,'TrinityAll-salmon-counts.txt'), sep = "\t", header = F)
    colnames(counts_salmon_raw) <- c("sample","ctg","length","effective_length","tpm", "num_reads_mapped")

    funcs <- read.csv(paste0('interpro-results-fmt-',target,'TrinityAll.tsv'), sep = "\t")

    funcs <- funcs[order(funcs$source, funcs$e_value),]  # sort for grabbing the top of multiples
    funcs$gene_callers_id <- gsub("_.*$","", as.character(funcs$gene_callers_id))
    # drop duplicate annotations for duplicate/identical sequences ID'd by salmon and skipped
    funcs <- funcs[!(funcs$gene_callers_id %in% genecalls$gene_callers_id[genecalls$contig %in% dups$DuplicateTxp]),]   # ONLY FOR NSTELLA

    funcs <- distinct(funcs, source, gene_callers_id, .keep_all=T)  # take only the best e-value hit per gene per program if multiple

    # get names of the transcripts on which each gene is found
    funcs$transcript <- as.character(genecalls$contig)[as.numeric(as.character(funcs$gene_callers_id))]
    eggnog$transcript <- as.character(genecalls$contig)[as.numeric(as.character(eggnog$query_name))]

Note that this script processes one species + transcriptome combination
at a time (e.g., *N. stella* polyA or *N. stella* Trio), as specified by
the `target` variable. Also, for *N. stella*, duplicate transcripts were
collapsed by running the lines with the comment `## ONLY FOR NSTELLA`;
for *B. argentea* there were no duplicates and so these lines were not
run.

At this point, the RapClust data are in the format
`transcript   cluster`. This information was then used to combine
transcripts, abundances, functions:

    # make a dict of transcript IDs and corresponding RapClust clusters
    rapDict <- setNames(rapclust$V1, rapclust$V2)

    # lump by rapClusters
    funcs$cluster <- as.character(funcs$transcript) # copy transcripts to new column to use in renaming
    clustersForFunctions <- as.character(rapDict[as.character(funcs$cluster)]) # get clusters, in the order in which they occur in funcs
    funcs$cluster[funcs$cluster %in% names(rapDict)] <- clustersForFunctions[!is.na(clustersForFunctions)]  # order is the same, so just drop the NAs and slot in
    eggnog$cluster <- as.character(eggnog$transcript) # copy transcripts to new column to use in renaming
    clustersForFunctions <- as.character(rapDict[as.character(eggnog$cluster)]) # get clusters, in the order in which they occur in funcs
    eggnog$cluster[eggnog$cluster %in% names(rapDict)] <- clustersForFunctions[!is.na(clustersForFunctions)]  # order is the same, so just drop the NAs and slot in

    rm(clustersForFunctions)

    # spread out by InterPro source (Pfam, etc) so can combine with eggnog
    funcs <- spread(funcs[,!(colnames(funcs) %in% c('e_value', 'accession'))], source, function.) # drop accessions (e.g. PF0000) and evalue to avoid needless dups
    funcs <- merge(funcs, eggnog, by.x = 'gene_callers_id', by.y = 'query_name', all = T)  # 1st 3 rows are from the eggnog header junk that was skipped CHECK
    funcs <- mutate(funcs, cluster=coalesce(cluster.x, cluster.y), transcript=coalesce(transcript.x, transcript.y))  # not all genes got annotated by both, so combine
    funcs <- funcs[,-which(colnames(funcs) %in% c('cluster.x','cluster.y','transcript.x','transcript.y'))]

    # convert factors to characters, then group by RapClust cluster and flatten all the character columns (transcript, Pfam, etc) for each cluster into |separated list
    clustHomogeneity <- mutate_if(funcs, is.factor, as.character) %>% group_by(cluster) %>% summarise_if(is.character, function(x) paste(unique(x), collapse="|"))
    clustHomogeneity <- as.data.frame(clustHomogeneity)
    clustHomogeneity <- clustHomogeneity[-nrow(clustHomogeneity),]  # drop last row as this has no clusterID


    ##### COUNTED BY SALMON

    # overwrite transcript IDs, in column 'ctg' with clusters
    counts_salmon_raw$ctg <- as.character(counts_salmon_raw$ctg)
    counts_salmon_raw$ctg[counts_salmon_raw$ctg %in% names(rapDict)] <-
      as.character(rapDict[as.character(counts_salmon_raw$ctg[counts_salmon_raw$ctg %in% names(rapDict)])])

    # print number of transcript clusters for the record
    length(unique(as.character(counts_salmon_raw$ctg)))

The metadata for each sample (e.g., oxygen condition, nitrate or
peroxide addition, etc) was then extracted from the sample names and
used to construct a treatment table to use for model creation.

    metadata <- unique(as.character(counts_salmon_raw$sample))

    # make a dataframe of all the samples, to have all the metadata from each sample
    metadata <- data.frame(oxygen = metadata,  condition = metadata, check = metadata, stringsAsFactors = F)

    # now just extract the treatment/condition info from each name with gsub
    metadata$oxygen[grepl("[Aa]noxia", metadata$oxygen) & grepl("[Ss]ulfide", metadata$oxygen)] <- "anoxiaSulfide"
    metadata$oxygen[grepl("[Aa]noxia", metadata$oxygen) & !grepl("[Ss]ulfide", metadata$oxygen)] <- "anoxiaOnly"
    metadata$oxygen[grep("[Hh]ypoxia", metadata$oxygen)] <- "hypoxia"
    metadata$oxygen[grep("[Ss]hip", metadata$oxygen)] <- "ship"
    metadata$oxygen[grep("[Aa][e]*r[e]*ated", metadata$oxygen)] <- "aerated"

    metadata$condition[grep("[Bb]oth", metadata$condition)] <- "both"
    metadata$condition[grep("NO3|[Nn]itrate", metadata$condition)] <- "nitrate"
    metadata$condition[grep("H2O2|[Pp]eroxide", metadata$condition)] <- "peroxide"
    metadata$condition[grep("[Cc]ontrol", metadata$condition)] <- "control"
    metadata$condition[grep("[Aa][e]*r[e]*ate|[Ss]hip",metadata$condition)] <- "NONE"

    # get the order of factor levels right for each treatment, and set the baseline treatment as the first level to be the control later on
    metadata <- metadata
    metadata$oxygen <- factor(metadata$oxygen)
    metadata$oxygen <- relevel(metadata$oxygen, "ship")   ## IMPORTANT this sets the baseline for OXYGEN
    metadata$condition <- factor(metadata$condition)
    metadata$condition <- relevel(metadata$condition, "both")   ## IMPORTANT this sets baseline for CONDITION
    rownames(metadata) <- metadata$check

    metadata$concat <- factor(paste0(metadata$oxygen, "-", metadata$condition))

A function was then written to take the counts for each transcript and
combine them by RapClust cluster. During this step, RapClust clusters
consisting entirely of unannotated transcripts were dropped as they
complicated the analysis but could not be directly interpreted.

    # get matrix for deseq
    processCounts <- function(counts, dropUnannotated=F){
      
      counts <- dcast(counts, ctg ~ sample, value.var = 'num_reads_mapped', fun.aggregate = sum)  # sum adds together rapclust-binned transcripts
      rownames(counts) <- counts$ctg
      counts <- as.matrix(counts[,-1])
      
      #status
      print(paste0(nrow(counts)," - this many transcripts/bins"))
      
      if (dropUnannotated){
        counts <- counts[unique(as.character(clustHomogeneity$cluster)),]
        
        print(paste0("Left with ", nrow(counts)," transcripts/bins after dropping unannotateds"))
      }
      
      # fix any NAs introduced from contigs not having any reads mapped and so being introduced as NAs when reshaping
      counts[is.na(counts)] <- 0
      print('NAs fixed')
      
      return(counts)
    }

Note that since the experimental design was not fully factorial (i.e.,
shipboard and aerated samples do not have nitrate/peroxide amendments),
the design can be thought of as two separate experiments: one of various
oxygen treatments and one of amendment.

    # make rectangular (clusters/tx by samples)
    counts_salmon <- processCounts(counts_salmon_raw, dropUnannotated = T)


    counts_salmon_oxygen <- counts_salmon[,metadata$check[!(metadata$condition %in% c('nitrate','peroxide'))]]
    counts_salmon_amendment <- counts_salmon[,metadata$check[metadata$condition %in% c('both','control')]]

Functions were then written to create a model matrix for samples
included in a counts matrix, turn the counts data into a DGE structure
used by edgeR and normalize the counts, and fit the model to the
normalized data with `lmFit()` using mean/variance estimation by
`voom()`:

    # function to set up the model
    makeModel <- function(countsMatrix, m = metadata, metadataCol="condition", legendPos='bottomright'){
      # limma design
      limma_design_bolovina <- model.matrix(as.formula(paste0("~",metadataCol)), data = m[colnames(countsMatrix),])  # treatments relative to "control"
      
      limma_design_bolovina <- limma_design_bolovina[,colSums(limma_design_bolovina) > 0]
      
      return(limma_design_bolovina)
    }

    makeModel(counts_salmon_nitrate, metadataCol = 'oxygen')

    makeDGE <- function(counts = counts_salmon_hypoxia, mCol="condition", design = makeModel(counts, metadataCol = mCol)) {
      counts <- counts[,rownames(design)]
      
      # reload the trimmed set
      d <- DGEList(counts = counts)
      # normalize with TMM
      d <- calcNormFactors(d)
      
      return(d)
    }

    fitModel <- function(counts = counts_salmon_hypoxia, mCol="condition"){ 
      design <- makeModel(counts, metadataCol = mCol)
      
      d <- makeDGE(counts = counts, design = design)
      
      v <- voom(d, design, plot = T)
      
      # linear fit the voom estimation
      f <- lmFit(v, design)
      f <- eBayes(f)
      
      return(f)
    }

A model was then fit for the oxygen treatments, and separately, for the
nitrate/peroxide amendments.

    fit_salmon_oxygen <- fitModel(counts_salmon_oxygen, mCol='oxygen')
    fit_salmon_amendment <- fitModel(counts_salmon_amendment, mCol='condition')

Two custom plotting functions were then written to take the normalized
counts data, subset to speficied functions, and aggregate the normalized
expression by taking the mean of samples per treatment level (e.g. the
mean of all ‘anoxia control’ samples).

The first function, `plotCpmByTreat()` plots the CPM values for the
specified data as a heatmap.

The second function, `plotPieByTreat()` plots the CPM values as dot
size, and apportions the dot by color as a pie chart based on taxonomy.

    plotCpmByTreat <- function(f = fit, dge=dge_salmon, functions = funcs, m = metadata, contrastCols = c(2), keepN = 100, func.source = 'Pfam', clust=F,
                                   customGenes=NULL, sumByFunction=T, collapseTo='concat', minLFC=2, abbrevFunctions=T, customScale=NULL) {
        # get log cpm differential expression
        de <- cpm(dge, log = F, prior.count = 0)
        
        if (is.null(customGenes)) {
          genesToKeep <- rownames(topTable(f, coef = contrastCols, lfc=minLFC, number = keepN))
        } else {
          genesToKeep <- unique(customGenes)
          genesToKeep <- genesToKeep[genesToKeep %in% rownames(de)]
          if (length(genesToKeep) < length(customGenes)) {
            print(paste0("WARNING!!! You asked for ", length(customGenes), " but ", length(customGenes) - length(genesToKeep), " didn't exist. Is there a reason?"))
          }
        }
        
        plotCPM <- melt(de[genesToKeep,])
        colnames(plotCPM) <- c('Transcript','Sample','CPM')
        
        plotCPM$Treatment <- m[as.character(plotCPM$Sample), collapseTo]
        
        deDist <- vegdist(de[genesToKeep,], method = 'euc')
        deClust <- hclust(deDist, method = 'average')
        
        sampDist <- vegdist(t(de[genesToKeep,]), method='euc')
        sampClust <- hclust(sampDist, method='average')
        
        # group samples by treatment
        if (clust) {
          plotCPM$Sample <- factor(plotCPM$Sample, levels=sampClust$labels[sampClust$order])
        } else {
          plotCPM$Treatment <- factor(plotCPM$Treatment, levels = sort(unique(as.character(plotCPM$Treatment))))
          
          plotCPM$Sample <- factor(plotCPM$Sample, levels = unique(as.character(plotCPM$Sample[order(plotCPM$Treatment)])))
        }

        # subset functions and make rownames
        rownames(functions) <- as.character(functions$cluster)
        
        # get the functions in the right order
        functions <- as.character(functions[as.character(plotCPM$Transcript),func.source])
        functions[is.na(functions)] <- ""  # make na's blank strings so doesn't skip over
        if (abbrevFunctions == T) {
          functions <- gsub("\\|.*$","| ...",functions)
        }
        

        if (sumByFunction) {
          print('SUMMING!!')

          plotCPM$functions <- as.factor(functions)  # already in right order so paste in
          plotCPM <- plotCPM %>% group_by(Treatment, functions) %>% summarise(CPM=log2(mean(CPM) + 1))
          
          # group transcripts
          summedCPM <- as.data.frame(plotCPM)
        }

        if (!is.null(names(customGenes)) | all(customGenes %in% plotCPM$functions)) {
          plotCPM$functions <- factor(plotCPM$functions, levels=as.character(customGenes))
        }
        
        p <- ggplot(plotCPM, aes(x = functions, y = Treatment)) +
          geom_tile(aes(fill = CPM)) + theme_minimal() +
          scale_x_discrete(expand = c(0,0)) + scale_y_discrete(expand = c(0,0)) +
          labs(size = "CPM (log2)", fill="CPM (log2)") +
          scale_fill_gradient(low = 'lightgoldenrodyellow', high = 'red', na.value = 'grey82') +
          scale_size_continuous(breaks= function(x) seq(round(x[1]), round(x[2]), length.out = 4)) +
          theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 6), panel.grid = element_blank(), panel.background = element_rect(fill='grey82'))
        
        if (!is.null(names(customGenes))) {
          p <- p + scale_x_discrete(expand = c(0,0), labels=names(customGenes))#[order(customGenes)])
        }
        if (!is.null(customScale)) {
          if (is.na(customScale[3])) {
            customScale <- c(customScale,4)
          }
          p <- p + scale_fill_gradient(low = 'lightgoldenrodyellow', high='red', na.value = 'grey82', limits = c(min(customScale), max(customScale)),
                                       breaks = seq(min(customScale), max(customScale), length.out = customScale[3]))
        }
        print(p)
    }



    plotPieByTreat <- function(f = fit, dge=dge_salmon, functions = funcs, m = metadata, contrastCols = c(2), keepN = 100, func.source = 'Pfam', clust=F,
                                   customGenes=NULL, sumByFunction=T, targetTaxLevel='best_tax_level', abbrevFunctions=T, abbrevTax=F, collapseTo='concat') {
      # get log cpm differential expression
      de <- cpm(dge, log = F, prior.count = 0)
      
      if (is.null(customGenes)) {
        genesToKeep <- rownames(topTable(f, coef = contrastCols, lfc=minLFC, number = keepN))
      } else {
        genesToKeep <- unique(customGenes)
        genesToKeep <- genesToKeep[genesToKeep %in% rownames(de)]
        if (length(genesToKeep) < length(customGenes)) {
          print(paste0("WARNING!!! You asked for ", length(customGenes), " but ", length(customGenes) - length(genesToKeep), " didn't exist. Is there a reason?"))
        }
      }
      
      plotCPM <- melt(de[genesToKeep,])
      colnames(plotCPM) <- c('Transcript','Sample','CPM')
      
      plotCPM$Treatment <- m[as.character(plotCPM$Sample), collapseTo]
      
      # group samples by treatment
      if (clust) {
        deDist <- vegdist(de[genesToKeep,], method = 'euc')
        deClust <- hclust(deDist, method = 'average')
        
        sampDist <- vegdist(t(de[genesToKeep,]), method='euc')
        sampClust <- hclust(sampDist, method='average')
      
        plotCPM$Sample <- factor(plotCPM$Sample, levels=sampClust$labels[sampClust$order])
      } else {
        plotCPM$Treatment <- factor(plotCPM$Treatment, levels = sort(unique(as.character(plotCPM$Treatment))))
        
        plotCPM$Sample <- factor(plotCPM$Sample, levels = unique(as.character(plotCPM$Sample[order(plotCPM$Treatment)])))
      }
      
      # subset functions and make rownames
      rownames(functions) <- as.character(functions$cluster)
      
      # get the functions in the right order
      if (length(unique(as.character(functions[unique(as.character(plotCPM$Transcript)), 'best_tax_level']))) < 99) {
        tax <- as.character(functions[unique(as.character(plotCPM$Transcript)), targetTaxLevel])
      } else if (targetTaxLevel == 'best_tax_level'){
        print("You have >99 different tax levels, but stubbornly want them (change targetTaxLevel to  be 'tax_scope' if you repent)")
        tax <- as.character(functions[unique(as.character(plotCPM$Transcript)), 'best_tax_level'])
      } else {
        print("You have >99 different tax levels, aggregating to tax_scope instead of best_tax_level")
        tax <- as.character(functions[unique(as.character(plotCPM$Transcript)), 'tax_scope'])
      }
      
      functions <- as.character(functions[unique(as.character(plotCPM$Transcript)),func.source])
      
      ##functions[is.na(functions)] <- ""  # make na's blank strings so doesn't skip over
      tax[is.na(tax)] <- "Unknown"  # make na's blank strings so doesn't skip over
      if (abbrevFunctions == T){
        functions <- gsub("\\|.*$","| ...",functions)
      }
      
      if (abbrevTax == T){
        tax <- gsub("\\|.*$","| ...",tax)
      }
      
      if (sumByFunction) {
        print('SUMMING!!')
        
        plotCPM$functions <- as.factor(functions)  # already in right order so paste in
        plotCPM$tax <- as.factor(tax)  # already in right order so paste in
        
        plotCPM <- plotCPM %>% #group_by(Treatment, functions, tax) %>% summarise(CPM=sum(CPM)) %>%  
          group_by(Treatment, functions, tax) %>% summarise(CPM=sum(CPM)) %>% mutate(normCPM = CPM/sum(CPM), cumCPM = log2(sum(CPM) + 1))

        #print(plotCPM)
      }
      
      
      # make a very inelegant legend
      plotCPM <- rbind(as.data.frame(plotCPM), 
                       data.frame(Treatment=paste0('LEGEND',seq(round(min(plotCPM$cumCPM)), round(max(plotCPM$cumCPM)), length.out = 5)), 
                                  functions=rep('LEGEND',5), tax=rep('LEGEND',5),
                                  CPM=seq(round(min(plotCPM$cumCPM)), round(max(plotCPM$cumCPM)), length.out = 5), 
                                  normCPM=seq(round(min(plotCPM$cumCPM)), round(max(plotCPM$cumCPM)), length.out = 5), 
                                  cumCPM=seq(round(min(plotCPM$cumCPM)), round(max(plotCPM$cumCPM)), length.out = 5)))
      
      plotCPM[is.na(plotCPM)] <- 0  # fix cases where it's 0/0 for normCPM

      ggplot(plotCPM, aes(x = cumCPM/2, y = normCPM, fill=tax, width=cumCPM)) + #x=functions
        coord_polar("y", start = 0) + theme_void() +
        geom_bar(position = 'fill', stat = 'identity') + 
        facet_grid(Treatment~functions) +
        scale_fill_hue(na.value = 'black') +
        scale_size_continuous(breaks= function(x) seq(round(x[1]), round(x[2]), length.out = 4)) +
        labs(fill= 'Best taxonomy')
    }

The cluster IDs for transcripts of interest were manually identified
based on high coverage and predicted function and are specified here:

    nstella_polyA <- setNames(c("cluster1117","cluster4836","cluster34930", "cluster13217", "cluster59869",
                                "cluster2582",
                                "cluster5632", "cluster57","cluster29558","cluster19494"), 
                              c("GOGAT","GS","GDH", 'NNP', "NRT", 
                                "NR",
                                "nirK", "nor", "Fe-hydrogenase","PFOR"))

    barg_polyA <- setNames(c("cluster3668","cluster26749","cluster206540", "cluster8555", "cluster56220", 
                             "cluster15120",
                             "cluster1205",
                             "cluster31395", "cluster249","cluster62445", "cluster155620"), 
                           c("GS","GOGAT","GDH", "NNP", "NRT",
                             "NR1",
                             "NR2",
                             "nirK", "nor", "Fe-hydrogenase","PFOR"))

For each host, the plot and statistics were then generated as follows:

    dset <- 'nstella_polyA'
    plotCpmByTreat(f = fit_salmon_oxygen, functions = clustHomogeneity, dge = makeDGE(counts_salmon_oxygen, mCol = 'oxygen'), m = metadata, sumByFunction = T,
                       contrastCols = c("oxygenanoxiaOnly","oxygenhypoxia","oxygenaerated","oxygenanoxiaSulfide"), keepN = Inf, func.source = 'cluster',
                       customGenes = get(dset), abbrevFunctions = F, customScale = c(0,14,5))
    ggsave(paste0("plots/Figure2B_", dset, "_heatmap.pdf"), width=5, height=3)
    write.table(x=topTable(fit_salmon_oxygen, number = Inf)[get(dset),], file = paste0("plots/Figure2B_", dset, "_DE_oxygen.tsv"), 
                sep = "\t", quote = F, row.names = paste0(get(dset), " (", names(get(dset)), ")"), col.names = NA)
    write.table(x=topTable(fit_salmon_amendment, number = Inf)[get(dset),], file = paste0("plots/Figure2B_", dset, "_DE_Amendment.tsv"), 
                sep = "\t", quote = F, row.names = paste0(get(dset), " (", names(get(dset)), ")"), col.names = NA)

We then compared the expression and predicted taxonomy of the de novo
assembled nuclear- and plastid-encoded transcripts from the Trio data as
below:

    # GAPDH, PGK, and rubisco via Pfam, FCP via EggNOG free text
    plotPieByTreat(f = fit_salmon_oxygen, functions = clustHomogeneity, dge = makeDGE(counts_salmon_oxygen, mCol = 'oxygen'), m = metadata, sumByFunction = T,
                       contrastCols = c("oxygenanoxiaOnly","oxygenhypoxia","oxygenaerated","oxygenanoxiaSulfide"), keepN = Inf, func.source = 'Pfam', collapseTo = 'oxygen',
                       customGenes = c(clustHomogeneity$cluster[grep("^Glyceraldehyde 3-phosphate dehydrogenase, C-terminal domain$|Phosphoglycerate kinase|^Ribulose bisphosphate carboxylase large chain, [a-z ]*$|^Ribulose bisphosphate carboxylase, small", clustHomogeneity$Pfam)],
                                       clustHomogeneity$cluster[grep("fucoxanthin chlorophyll", clustHomogeneity$eggNOG_free_text_description)]) ,
                   abbrevFunctions = F)

Which produced the raw plot for Figure 1B that was manually colored and
arragned.
