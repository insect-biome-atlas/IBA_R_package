# Code relies on data.table
library(data.table)

#' Getting IBA CO1 data
#'
#' Retrieves the requested CO1 data from the IBA project, given paths to the
#' raw data and processed data. The retrieved data can be calibrated using
#' biological spike-ins, if desired, and spike-in data can either be removed
#' or kept, as desired.
#' 
#' The `dataset` parameter allows the user to select the desired data subset.
#' The available options are `lysate`, `homogenate`, `soil`, `litter` and
#' `ethanol`. For "MG" (Madagascar), only `lysate` and `litter` data are
#' available. It is possible to request several subsets at once by concatenating
#' the dataset names. For instance, `lysate|homogenate` would retrieve both
#' lysate and homogenate data.
#' 
#' Note that only biological spike-in data are used to
#' calibrate reads. Spike-ins are identified by occurrence patterns (occurring
#' in more than 80% of samples) and taxonomic annotation (belonging to Insecta
#' if bioliogical spike-ins), calling the function `identify_spikes`.
#' 
#' The return value is an OTU read table that contains the full taxonomic
#' annotation of the OTU (using the recommended IBA annotation combining SINTAX
#' and EPA-NG annotations). The headers of the read columns are the
#' `sampleID_NGI` identifiers associated with the sequencing run in which the
#' data were generated. Thus, in a merged dataset it is not possible to
#' distinguish the data subsets without looking into the IBA metadata files.
#' 
#' Use the function `get_iba_sample_data` to get the associated field sample
#' and site metadata for those sequencing run identifiers. The function
#' `get_iba_taxon_data` can be used to retrieve such data together with read
#' numbers for specific taxa (individual OTUs or groups of OTUs with the same
#' taxonomic annotation at some rank).
#'
#' @param country The country for which data is requested (`"MG"` or `"SG"`).
#' @param dataset The dataset requested: `"lysate"`, `"homogenate"`, `"soil"`,
#' `"litter"` or `"ethanol"`. See *Details* for more information.
#' @param data_path The path to the processed data. Defaults to
#' `"~/dev/figshare-repos/iba/processed_data/v3/"`
#' @param metadata_path The path to the metadata (raw data). Defaults to
#' `"~/dev/figshare-repos/iba/raw_data/v6/"`
#' @param calibrate Whether to calibrate data using biological spike-ins.
#' Defaults to `FALSE`.
#' @param remove_spikes Whether to remove spike-ins. Defaults to `TRUE`.
#' @param format The desired data format ("data.table" or "data.frame")
#' @returns OTU read table in desired format with taxonomic annotation of OTUs
#' @examples
#' get_iba_co1_data(country="SE",dataset="lysate",calibrate=TRUE)
#' get_iba_co1_data(country="MG",dataset="lysate|homogenate")
#' @export
get_iba_co1_data <- function(country,
                             dataset="lysate|homogenate|soil|litter|ethanol",
                             data_path='~/dev/figshare-repos/iba/processed_data/v3/',
                             metadata_path='~/dev/figshare-repos/iba/raw_data/v6/',
                             calibrate=FALSE,
                             remove_spikes=TRUE,
                             format="data.frame"
                             ) {

    # Convert country to upper case
    country <- toupper(country)

    # Check that country is set correctly
    if (country!="MG" && country !="SE") {
        cat ("ERROR: 'country' must be one of 'MG' or 'SE'\n")
        return (NA)
    }

    # Initialize index into columns we want to keep
    index <- numeric()

    # Get data files
    if (format=="data.frame") {
       cat("\nImporting counts into a 'data.frame' object. If this takes more than a few minutes,\n")
       cat("or your system runs out of memory, try changing the format into 'data.table' using\n")
       cat ("'format=\"data.table\"'.\n\n")
    } else if (calibate) {
       cat("Importing counts into a 'data.table' object. Be prepared that the 'data.table' format\n")
       cat("will require a long time for read calibration, at least on MacOS systems. If this is\n")
       cat("the case, you can try using the default 'data.frame' format instead.\n\n")
    }

    if (country=="SE") {

      cat("Reading data files, this may take a while...\n")
      if (format=="data.table")
            counts <- fread(paste0(data_path,"cleaned_noise_filtered_cluster_counts_SE.tsv"),header=TRUE,sep="\t")
        else
            counts <- read.delim(paste0(data_path,"cleaned_noise_filtered_cluster_counts_SE.tsv"))
        if (colnames(counts)[1]!="cluster")
            colnames(counts)[1]<-"cluster"
        taxonomy <- read.delim(paste0(data_path,"cleaned_noise_filtered_cluster_taxonomy_SE.tsv"))
        taxonomy <- taxonomy[taxonomy$representative==1,]

        # Get metadata file and remove non-samples and sequencing failures
        meta <- read.delim(paste0(metadata_path,"CO1_sequencing_metadata_SE.tsv"))
        meta <- meta[meta$lab_sample_type=="sample" & meta$sequencing_successful==TRUE,]

        # Get and possibly adjust malaise trap data. We use the fact here
        # that the spike-ins are the same in both cases.
        if (grepl("lysate", dataset) || grepl("homogenate",dataset)) {
            malaise_samples <- character()
            if (grepl("lysate", dataset))
                malaise_samples <- with(meta, sampleID_NGI[dataset=="CO1_lysate_2019_SE"])
            if (grepl("homogenate", dataset)) {
                homogenates <- with(meta, sampleID_NGI[dataset=="CO1_homogenate_2019_SE"])
                malaise_samples <- c(malaise_samples, homogenates)
            }
            counts <- handle_spikes(counts, malaise_samples, taxonomy, calibrate, remove_spikes)
            index <- c(index, match(malaise_samples,colnames(counts)))
        }

        # Get ethanol data if requested
        if (grepl("ethanol", dataset)) {

            ethanol_samples <- with(meta, sampleID_NGI[dataset=="CO1_ethanol_2019_SE"])
            index <- c(index, match(ethanol_samples,colnames(counts)))
        }
 
        # Get soil and/or litter data if requested
        if (grepl("soil", dataset)) {

            soil_samples <- with(meta, sampleID_NGI[grepl("_S",sampleID_FIELD)])
            index <- c(index, match(soil_samples,colnames(counts)))
        }
        if (grepl("litter", dataset)) {

            litter_samples <- with(meta, sampleID_NGI[grepl("_L",sampleID_FIELD) & lab_sample_type=="sample"])
            index <- c(index, match(litter_samples,colnames(counts)))
        }
    }
    if (country=="MG") {

        # Get data files
        cat("Reading data files, this may take a while...\n")
        if (format=="data.table")
            counts <- fread(paste0(data_path,"cleaned_noise_filtered_cluster_counts_MG.tsv"),header=TRUE,sep="\t")
        else
            counts <- read.delim(paste0(data_path,"cleaned_noise_filtered_cluster_counts_MG.tsv"))
        if (colnames(counts)[1]!="cluster")
            colnames(counts)[1]<-"cluster"
        taxonomy <- read.delim(paste0(data_path,"cleaned_noise_filtered_cluster_taxonomy_MG.tsv"))
        taxonomy <- taxonomy[taxonomy$representative==1,]

        # Get metadata file and remove non-samples and sequencing failures
        meta <- read.delim(paste0(metadata_path,"CO1_sequencing_metadata_MG.tsv"))
        meta <- meta[meta$lab_sample_type=="sample" & meta$sequencing_successful==TRUE,]

        # Adjust lysate data if requested
        if (grepl("lysate", dataset)) {

            lysates <- with(meta, sampleID_NGI[dataset=="CO1_lysate_2019_MG"])
            counts <- handle_spikes(counts, lysates, taxonomy, calibrate, remove_spikes)
            index <- c(index, match(lysates,colnames(counts)))
        }

        # Get litter data if requested
        if (grepl("litter", dataset)) {

            litter_samples <- with(meta, sampleID_NGI[dataset=="CO1_litter_2019_MG"])
            index <- c(index, match(litter_samples,colnames(counts)))
        }
    }
 
    # Extract indices, do not forget to keep the cluster column. Also, safeguard against
    # samples that have no counts data despite the metadata claiming the contrary
    index <- c(1, index)
    index <- index[!is.na(index)]
    if (class(counts)[1]=="data.table")
        counts <- counts[,..index]
    else
        counts <- counts[,index]

    # Remove clusters that are not encountered in this dataset
    tot <- rowSums(counts[,2:ncol(counts)])
    counts <- counts[tot!=0,]

    # Remove samples with no data
    tot <- colSums(counts[,2:ncol(counts)])
    index <- c(TRUE,tot!=0)
    if (class(counts)[1]=="data.table")
        counts <- counts[,..index]
    else
        counts <- counts[,index]

    # Make sure artificial spikeins are included (they are not in the taxonomy file by default)
    if (sum(counts$cluster %in% taxonomy$cluster) < length(counts$cluster)) {
        art_spikes <- counts$cluster[!(counts$cluster %in% taxonomy$cluster)]
        if (remove_spikes==TRUE || !grepl("homogenate",dataset)) {
            cat ("ERROR: Unexpected clusters in counts file:\n")
            cat (art_spikes,"\n")
            return (NA)
        }
        for (i in 1:length(art_spikes)) {
            taxonomy <- rbind(taxonomy, rep("",times=ncol(taxonomy)))
            taxonomy[nrow(taxonomy),"cluster"] <- art_spikes[i]
        }
    }

    # Add in taxonomy data
    if (class(counts)[1]=="data.table")
        ret <- merge(data.table(taxonomy), counts, by="cluster")
    else
        ret <- merge(taxonomy, counts, by="cluster")

    # Make absolutely sure we remove Zoarces gillii (Species==NA for artificial spike-ins, if present)
    ret[is.na(ret$Species) | ret$Species!="Zoarces gillii",]
}


#' Handle spike-ins
#'
#' This function identifies biological and synthetic spike-ins. It then
#' calibrates the data using the biological spike-in data, if this is requested.
#' Finally, it removes the spike-ins, if this is requested.
#' 
#' In the OTU table passed in as `counts`, the first column is expected to
#' contain the OTU cluster names. The function handles both `data.frame` and
#' `data.table` formats for the OTU table. Note that data in `data.table` format
#' takes a long time to calibrate, at least on MacIntosh systems.
#' 
#' If some of the specified `samples` columns do not contain spike-in reads at
#' all, the read numbers are not changed.
#' 
#' The calibrated read numbers are returned as integers. The numbers are rounded
#' upwards to the nearest integer to make sure that presences are not turned
#' into absences for samples with few reads.
#' 
#' @param counts The OTU table as a data.frame or data.table object.
#' @param samples The names of the columns in the OTU table containing the
#' spike-in data, and for which the reads should (potentially) be calibrated.
#' @param taxonomy The taxonomic annotation of the OTUs.
#' @param calibrate Whether to calibrate data using biological spike-ins.
#' @param remove_spikes Whether to remove spike-ins.
#' @returns OTU read table in desired format with taxonomic annotation of OTUs
#' @examples
#' handle_spikes(counts, samples, taxonomy, calibrate=TRUE, remove_spikes=TRUE)
#' @export
handle_spikes <- function(counts, samples, taxonomy, calibrate, remove_spikes) {
  
    # Make sure first coumn header is set correctly
    colnames(counts)[1] <- "cluster"

    # get column indices in the counts data table
    # safe_guard against samples that are missing in the counts table
    idx <- match(samples,colnames(counts))
    idx <- idx[!is.na(idx)]

    # identify spikeins
    if (remove_spikes || calibrate) {
        res <- identify_spikes(counts, samples, taxonomy)
        spikein_clusters <- res$spikein_clusters
        if (length(spikein_clusters)==0) {
            cat("ERROR: Could not find any spikein clusters\n")
            return (NA)
        }
    }

    # calibrate
    # Note that there are occasional samples without spike-ins; we simply do
    # not correct these read counts (what else can we do?). Presumably, spike-ins
    # were never added to these samples.
    # We only calibrate based on biological spike-ins, and based on the mean of the
    # log of the sum of spike-in counts 
    if (calibrate && length(spikein_clusters) > 0) {
        if (class(counts)[1]=="data.table") {
            cat("Trying to calibrate counts in a 'data.table' format.If this takes exceptionally long,\n")
            cat("you can try to read in or convert the OTU table to the standard 'data.frame' instead.\n")
            cat("If you have enough memory for this operation, the calibration will run much faster,.\n")
            cat("at least on MacOS systems.\n")
        }
        cat("\nCalibrating...\n")
        if (class(counts)[1]=="data.table")
            spike_counts <- colSums(counts[counts$cluster %in% res$bio_spikeins,..idx])
        else
            spike_counts <- colSums(counts[counts$cluster %in% res$bio_spikeins,idx])
        correction <- log10(spike_counts) - mean(log10(spike_counts[spike_counts!=0]))
        correction[spike_counts==0] <- 0.0
        pb <- txtProgressBar(min = 0, max = length(idx), style = 3, width = 50)
        for (i in 1:length(idx)) {
            j <- idx[i]
            if (class(counts)[1]=="data.table")
                counts[,j] <- ceiling(counts[,..j] / 10^(correction[i]))
            else
                counts[,j] <- ceiling(counts[,j] / 10^(correction[i]))
            setTxtProgressBar(pb, i)
        }
        cat("\n")
    }

    # Remove the spike-ins
    if (calibrate || remove_spikes) {
        counts <- counts[!(counts$cluster %in% spikein_clusters),]
    }

    counts
} 

