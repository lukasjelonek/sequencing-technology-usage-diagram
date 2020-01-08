
params.out = "results"
Channel.value([2007, 2019]).set{ch_years}


process lookup_data {
    publishDir params.out, mode: 'copy'
    
    container "quay.io/biocontainers/entrez-direct:10.2--pl526_0"

    input:
    set val(start), val(end) from ch_years

    output:
    file 'data.csv' into ch_data

    shell:
    '''
    #!/bin/bash
    PLATFORM=("abi solid" "bgiseq" "capillary" "complete genomics" "helicos" "illumina" "ion torrent" "ls454" "oxford nanopore" "pacbio smrt")
    TYPE=("genomic" "genomic single cell" "metagenomic" "metatranscriptomic" "other" "synthetic" "transcriptomic" "transcriptomic single cell" "viral rna")
    echo "technology,type,year,count" > data.csv
    for p in "${PLATFORM[@]}"; do
      for t in "${TYPE[@]}"; do
        for y in $(seq !{start} !{end}); do
          COUNT=$(esearch -db sra -query "(((\\"$y/01/01\\"[Modification Date] : \\"$y/12/31\\"[Modification Date])) AND \\"$p\\"[Platform]) AND \\"$t\\"[Source]" | xtract -pattern ENTREZ_DIRECT -element Count)
          echo ${p},${t},${y},${COUNT} >> data.csv
        done
      done
    done
    '''
}

process plot {
    publishDir params.out, mode: 'copy'
    
    container "r-ggplot2:latest"

    input:
    file 'data.csv' from ch_data

    output:
    file "*" into ch_end

    script:
    '''
    #!/usr/bin/env Rscript
    library(ggplot2)
    library(scales)
    # prepare data
    data = read.csv("data.csv")
    df = data.frame(data)
    df$year = as.character(df$year)


    # prepare common plot elements 
    p_palette = scale_fill_brewer(palette="Set1")
    p_scale_log = scale_y_continuous(labels=comma, trans="log10")
    p_scale_norm = scale_y_continuous(labels=comma)
    
    # plot 
    ggplot(data = df[which(df$count>0),], aes(x=year, y=count, fill=type)) +
        theme_bw() +
        geom_bar(stat='identity', position=position_dodge(preserve = "single")) +
        p_palette +
        facet_grid(facets =  technology ~ .) +
        p_scale_log + 
        theme(axis.text.x = element_text(angle=45, hjust=1), legend.position=c(0.15,0.87), legend.text = element_text(size=8), legend.box.background=element_rect(colour = "black")) + 
        labs(title="New public experiments in SRA per year by platform and source", x="Year", y="#Experiments", fill="Sample source")
    ggsave("experiment_count_log_by_platform_and_source.svg", width=16, height=24, units="cm")
    ggsave("experiment_count_log_by_platform_and_source.png", width=16, height=24, units="cm")

    aggregate_by_technology = aggregate.data.frame(df$count, by=list(technology=df$technology, year=df$year), FUN=sum)
    ggplot(data = aggregate_by_technology[which(aggregate_by_technology$x>0),], aes(x=year, y=x, fill=technology)) +
        theme_bw() +
        geom_bar(stat='identity', position=position_dodge()) +
        p_palette +
        p_scale_norm +
        theme(axis.text.x = element_text(angle=45, hjust=1), legend.position=c(0.15,0.6), legend.text = element_text(size=8), legend.background=element_blank()) + 
        labs(title="New public experiments in SRA per year by platform", x="Year", y="#Experiments", fill="Sequencing platform")
    ggsave("experiment_count_by_platform.svg", width=16, height=12, units="cm")
    ggsave("experiment_count_by_platform.png", width=16, height=12, units="cm")

    aggregate_by_type = aggregate.data.frame(df$count, by=list(type=df$type, year=df$year), FUN=sum)
    ggplot(data = aggregate_by_type[which(aggregate_by_type$x>0),], aes(x=year, y=x, fill=type)) +
        theme_bw() +
        geom_bar(stat='identity', position=position_dodge()) +
        p_palette +
        p_scale_norm +
        theme(axis.text.x = element_text(angle=45, hjust=1), legend.position=c(0.15,0.7), legend.text = element_text(size=8), legend.background=element_blank()) + 
        labs(title="New public experiments in SRA per year by source", x="Year", y="#Experiments", fill="Sample source")
  ggsave("experiment_count_by_source.svg", width=16, height=12, units="cm")
  ggsave("experiment_count_by_source.png", width=16, height=12, units="cm")
    '''
}
