# Code for getting IBA samples and sites for specific taxa


# Function for getting the IBA samples containing a specific taxon
#
# Parameters
#  -     iba_data: the IBA data (from the get_iba_co1_data fxn)
#  -        taxon: the taxon (or taxa) we are interested in (ASVs are also allowed)
#  -         rank: the rank of the taxon ("ASV" for ASVs, "cluster" for clusters)
#
# Return value: data.frame with read numbers and sample data for all individual
#               clusters matching the taxon
get_taxon_samples <- function(iba_data,
                              taxon,
                              rank,
                              metadata_path) {

    # Make sure the rank is in column names
    colnames(iba_data) -> names
    idx <- which(names==rank)
    if (is.na(idx)) {
        cat ("ERROR: A rank named ", rank, "was not found\n")
        return (NA)
    }

    # Extract the relevant rows in the data table 
    x <- D[,..idx]
    include_rows <- x[[1]]==taxon
    D <- D[include_rows,]

    # Extract the columns with non-zero reads
    sample_cols <- grepl("_",colnames(D)) & !grepl("BOLD",colnames(D))
    tot_reads <- colSums(D[,..sample_cols])
    include_cols <- c(rep(TRUE,times=sum(!sample_cols)),tot_reads>0)
    D <- D[,..include_cols]

    # Now melt the data table and only keep samples with reads
    D <- reshape2::melt(D,id.vars=1:sum(!sample_cols),variable.name="sampleID_NGI",value.name="reads")
    
    # Return as data frame
    data.frame(D[D$reads>0,])
}

