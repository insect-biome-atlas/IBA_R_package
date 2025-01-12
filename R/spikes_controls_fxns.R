library(data.table)

#' Stats for OTU/ASV partitions
#' 
#' This function computes various statistics for the reads of a set of OTUs
#' or ASVs identified for potential removal (for instance, because they occur
#' in too many control samples). The set is specified as an index (a vector of
#' `TRUE` and `FALSE` values of the same length as the number of rows in the 
#' `counts`table, and should specify the rows identified for potential removal.
#' Statistics are provided for both the rows identified for removal and the rows
#' that would remain after removal.
#' 
#' @param counts The OTU/ASV table as a data.frame or data.table object.
#' @param index Vector of booleans indicating which rows are considered for
#' removal.
#' @returns A list of statistics on the clusters to be removed or kept.
#' @export
mean_max <- function(counts,index) {

    is_missing <- rowSums(counts[,-1])==0

    remove_stats <- index & !is_missing
    remove_mean <- rowSums(counts[remove_stats,-1]) / rowSums(counts[remove_stats,-1]>0)
    remove_max <- apply(as.matrix(counts[remove_stats,-1]),1,max)
    remove_prop <- rowMeans(counts[remove_stats,-1]>0)

    is_missing_keep <- rowSums(counts[!index,-1])==0
    keep_stats <- !index & !is_missing
    keep_mean <- rowSums(counts[keep_stats,-1]) / rowSums(counts[keep_stats,-1]>0)
    keep_max <- apply(as.matrix(counts[keep_stats,-1]),1,max)
    keep_prop <- rowMeans(counts[keep_stats,-1]>0)

    list(remove_mean=remove_mean,
         remove_max=remove_max,
         remove_prop=remove_prop,
         remove_missing=sum(is_missing & index),
         keep_mean=keep_mean,
         keep_max=keep_max,
         keep_prop=keep_prop,
         keep_missing=sum(is_missing & !index))
}


# Function to identify control clusters
#
#   counts:     cluster counts for all samples (including controls)
#   samples:    names of samples
#   controls:   names of controls
#   cutoff:     threshold for removing control clusters
#' Identify control_clusters
#'
#' Identifies control clusters according to specified criteria.
#' 
#' @param counts OTU table, with first column being cluster name.
#' @param taxonomy The taxonomic annotation (assumed to contain cluster and
#' taxonomic annotations for all ASVs)
#' @param samples Names (column headers) of the real samples.
#' @param controls Names (column headers) of the control samples.
#' @param cutoff Threshold for identifying control clusters, given in terms of the
#' proportion of control samples in which a cluster is allowed to occur before
#' it is removed. The default threshold is set to 0.05, meaning that clusters
#' occurring in more than 5% of control samples will be removed.
#' @returns A list containing a vector of identified control clusters and
#' a data frame containing detailed information on these clusters.
identify_control_clusters <- function(counts, taxonomy, samples, controls, cutoff=0.05) {

    # Identify clusters that appear in controls above cutoff
    orig_tax <- taxonomy
    taxonomy <- taxonomy[taxonomy$representative==1,]
    # Extract sample and control reads
    idx <- which(colnames(counts) %in% samples)
    idx <- c(1,idx) # keep cluster name
    if (class(counts)[1]=="data.table")
        sample_counts <- counts[,..idx]
     else
        sample_counts <- counts[,idx]

    idx <- which(colnames(counts) %in% controls)
    idx <- c(1,idx) # keep cluster name
    if (class(counts)[1]=="data.table")
        control_counts <- counts[,..idx]
     else
        control_counts <- counts[,idx]

    # identify clusters that have at least one read in a control
    include_rows <- rowSums(control_counts[, 2:ncol(control_counts)])>0
    # get counts of the control clusters in the controls
    control_counts <- control_counts[include_rows,]
    # get counts of the control clusters in the samples
    sample_counts <- sample_counts[include_rows,]

    # Output some diagnostics
    tot_clusters <- nrow(control_counts)
    # prop_controls is the proportion of control samples in which each cluster occurs
    prop_controls <- rowMeans(control_counts[,2:ncol(control_counts)]>0)

    cat("There are", tot_clusters, "clusters in the control samples\n")
    cat("Of these,", sum(prop_controls>cutoff), "occur in more than", cutoff, "of samples\n")
    # identify clusters that occur in more than cutoff proportion of control samples
    remove_clusters <- control_counts$cluster[prop_controls>cutoff]
    
    res <- mean_max(control_counts, prop_controls>cutoff)
    cat("Reads in controls of removed clusters:\n")
    remove_tax <- data.frame(list(cluster=remove_clusters,
                                  prop_controls=res$remove_prop,
                                  mean_reads=res$remove_mean,
                                  max_reads=res$remove_max,
                                  Genus=taxonomy$Genus[match(remove_clusters,taxonomy$cluster)],
                                  Species=taxonomy$Species[match(remove_clusters,taxonomy$cluster)],
                                  BOLD_bin=taxonomy$BOLD_bin[match(remove_clusters,taxonomy$cluster)]))
    #print(remove_tax)
    cat("Summary of reads in controls of kept clusters:\n")
    cat("prop_samples:\n")
    print(summary(res$keep_prop))
    cat("mean:\n")
    print(summary(res$keep_mean))
    cat("max:\n")
    print(summary(res$keep_max))

    prop_samples <- rowMeans(sample_counts[,2:ncol(sample_counts)]>0)
    res <- mean_max(sample_counts, prop_controls>cutoff)
    cat("Reads in samples of removed clusters:\n")
    remove_tax <- data.frame(list(cluster=remove_clusters,
                                  prop_samples=res$remove_prop,
                                  mean_reads=res$remove_mean,
                                  max_reads=res$remove_max,
                                  Genus=taxonomy$Genus[match(remove_clusters,taxonomy$cluster)],
                                  Species=taxonomy$Species[match(remove_clusters,taxonomy$cluster)],
                                  BOLD_bin=taxonomy$BOLD_bin[match(remove_clusters,taxonomy$cluster)]))
    print(remove_tax)
    cat("Summary of reads in samples of kept clusters:\n")
    cat("prop_samples:\n")
    print(summary(res$keep_prop))
    cat("mean:\n")
    print(summary(res$keep_mean))
    cat("max:\n")
    print(summary(res$keep_max))

    list(remove_clusters=remove_clusters, remove_tax=remove_tax)
}


#' Remove control_clusters
#'
#' The function removes the control clusters according to specified criteria
#' 
#' @param filtered_counts Filtered cluster counts; the first
#' column should be the cluster name.
#' @param all_cluster_counts Cluster counts for all samples (real samples and
#' control samples).
#' @param taxonomy The taxonomic annotation (assumed to contain cluster and
#' taxonomic annotations for all ASVs)
#' @param samples Names (column headers) of the real samples.
#' @param controls Names (column headers) of the control samples.
#' @param cutoff Threshold for removing control clusters, given in terms of the
#' proportion of control samples in which a cluster is allowed to occur before
#' it is removed. The default threshold is set to 0.05, meaning that clusters
#' occurring in more than 5% of control samples will be removed.
#' @returns A list containing the OTU table with filtered counts, with control
#' clusters removed. The list also includes the taxonomy with control clusters
#' removed, as well as a data frame containing information on the removed
#' control clusters.
remove_control_clusters <- function(filtered_counts,
                                    all_cluster_counts,
                                    taxonomy,
                                    samples,
                                    controls,
                                    cutoff=0.05) {

    res <- identify_control_clusters(all_cluster_counts, taxonomy, samples, controls, cutoff)
    remove_clusters <- res$remove_clusters
    remove_tax <- res$remove_tax

    res <- list(
        counts=filtered_counts[!filtered_counts$cluster %in% remove_clusters,], 
        filtered_tax=taxonomy[!taxonomy$cluster %in% remove_clusters,], 
        remove_tax=remove_tax)
    res
}


#' Identify spike-ins
#'
#' The function identifies biological and synthetic spike-ins using occurrence
#' pattern. Specifically, spike-ins are identified as clusters occurring in
#' more than the `cutoff` proportion of samples. For biological spike-ins, it
#' is also required that they are annotated as belonging to Insecta. Synthetic
#' spike-ins are identified as those candidate spike-in clusters missing a
#' taxonomic annotation.
#' 
#' @param counts OTU table with first column being OTU cluster name
#' @param spikein_samples The samples (column headers in the OTU table)
#' supposed to contain spike-in reads.
#' @param taxonomy The taxonomic annotations for the OTU clusters. If a complete
#' ASV taxonomy is provided, the annotation of the representative ASVs is used.
#' @param cutoff Cutoff for identifying spike-in clusters, given in terms of
#' the proportion of samples. The default cutoff is set to 0.8, meaning that
#' clusters occurring in more than 80% of `spikein_samples` are identified as
#' candidate spike-in clusters.
#' @returns A list giving all spikein clusters identified, as well as separate
#' vectors of biological and synthetic spikeins. The taxonomic annotation of
#' the identified biological spike-in clusters is also returned.
identify_spikes <- function(counts, spikein_samples, taxonomy, cutoff=0.8) {

    # Get a rep asv taxonomy in case a complete cluster taxonomy is provided
    if ("representative" %in% colnames(taxonomy))
        taxonomy <- taxonomy[taxonomy$representative==1,]

    # Get counts for the samples containing spikeins
    # Account for counts being either data.frame or data.table
    idx <- which(colnames(counts) %in% spikein_samples)
    idx <- c(1,idx)
    if (class(counts)[1]=="data.table")
        counts <- counts[,..idx]
     else
        counts <- counts[,idx]

    # Identify spikeins
    prop_samples <- rowMeans(counts[,2:ncol(counts)]>0)
    spikein_candidates <- counts$cluster[prop_samples>cutoff]
    synthetic_spikeins <- spikein_candidates[is.na(match(spikein_candidates,taxonomy$cluster))]
    if (length(synthetic_spikeins>0)) {
        cat("\nFound", length(synthetic_spikeins), "synthetic spikein clusters: ")
        cat(synthetic_spikeins,sep=",")
        cat("\n")
    }
    bio_spikein_candidates <- spikein_candidates[!(spikein_candidates %in% synthetic_spikeins)]
    bio_spikeins <- bio_spikein_candidates[taxonomy$Class[match(bio_spikein_candidates,taxonomy$cluster)]=="Insecta"]
    cat("\nFound", length(bio_spikeins), "biological spikein clusters:\n")
    spike_tax <- taxonomy[taxonomy$cluster %in% bio_spikeins,c("cluster","Genus","Species","BOLD_bin")]
    spike_tax$prop_samples <- rowMeans(counts[match(spike_tax$cluster,counts$cluster),2:ncol(counts)]>0)
    print(spike_tax)

    # Return clusters
    res <- list(spikein_clusters=c(synthetic_spikeins, bio_spikeins), synthetic_spikeins=synthetic_spikeins, bio_spikeins=bio_spikeins, spike_tax=spike_tax)
    res
}


#' Remove spike-ins
#'
#' Identifies the spikeins by calling the `identify_spikes` function, and then
#' removes them.
#' 
#' @param counts OTU table with first column being OTU cluster name
#' @param spikein_samples The samples (column headers in the OTU table)
#' supposed to contain spike-in reads.
#' @param taxonomy The taxonomic annotations for the OTU clusters. If a complete
#' ASV taxonomy is provided, the annotation of the representative ASVs is used.
#' @param cutoff Cutoff for identifying spike-in clusters, given in terms of
#' the proportion of samples. The default cutoff is set to 0.8, meaning that
#' clusters occurring in more than 80% of `spikein_samples` are identified as
#' candidate spike-in clusters.
#' @returns A list containing the counts table with spikeins removed, as well
#' as the taxonomic annotation of the biological spikeins identified.
remove_spikes <- function(counts, spikein_samples, taxonomy, cutoff=0.8) {

    res <- identify_spikes(counts, spikein_samples, taxonomy, cutoff)
    spikein_clusters <- res$spikein_clusters
    spike_tax <- res$spike_tax

    res <- list(counts=counts[!(counts$cluster %in% spikein_clusters),], spike_tax=spike_tax)
    res
}

