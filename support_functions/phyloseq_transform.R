## Functions to transform phyloseq OTU tables

## TO DO:
# - add warning about negative values
# - add function description
# - add documentation and links (e.g., to DESeq2 manual)
# - make single wrapper function for all methods ??


#' @title Replace zeros in OTU abundance
#'
#' @param physeq A phyloseq-class object
#' @param method Which method to apply for zero replacement ("pseudocount")
#' @param pseudocount Constant value
#'
#' @return Phyloseq object with transformed counts in OTU table.
#' @export
#'
#' @examples
#'
phyloseq_replace_zero <- function(physeq, method = "pseudocount", pseudocount = 0.65){

  ## Extract OTU table
  tmp <- phyloseq::otu_table(physeq)

  ## Replace zeros with pseudocount
  if(method == "pseudocount"){
    tmp[ tmp == 0 ] <- pseudocount
  }

  ## Replace zeros with minimum observed OTU abundance
  if(method == "min"){
    trows <- taxa_are_rows(physeq)

    if(trows == TRUE){
      ## Find minimum non-zero abundances for each sample
      mins <- apply(X = tmp, MARGIN = 2, FUN = function(z){ z <- z[z > 0]; min(z) })

      ## Replace zeros with sample-specific count
      for(i in 1:ncol(tmp)){
        tmp[ tmp[, i] == 0, i ] <- mins[i]
      }
    }

    if(trows == FALSE){
      mins <- apply(X = tmp, MARGIN = 1, FUN = function(z){ z <- z[z > 0]; min(z) })
      for(i in 1:nrow(tmp)){
        tmp[ i, tmp[i, ] == 0 ] <- mins[i]
      }
    }
  } # end of 'min' method

  ## Replace OTU table
  phyloseq::otu_table(physeq) <- tmp
  return(physeq)
}



#' @title ALDEx2-based centred log-ratio transformation of OTU table
#'
#' @param physeq A phyloseq-class object
#' @param variable Variable name (from \code{\link[phyloseq]{sample_data}}) containing a descriptor for the samples, allowing them to be grouped and compared
#' @param iter Number of Monte Carlo samples to use
#' @details
#' ALDEx2 generates Monte Carlo samples of the Dirichlet distribution for each sample and
#' converts each instance using the centred log-ratio transformation.
#' @return Phyloseq object with CLR-transformed OTU counts (or a list of rarefied phyloseq-objects if iter > 1).
#' @export
#' @seealso \code{\link[ALDEx2]{aldex.clr}}
#' @references
#' Fernandes DA, Reid J, Macklaim MJ, McMurrough TA, Edgell DR, Gloor BG (2014) Unifying the analysis of high-throughput sequencing datasets: characterizing RNA-seq, 16S rRNA gene sequencing and selective growth experiments by compositional data analysis // Microbiome 2:15. DOI: 10.1186/2049-2618-2-15
#' @examples
#'
phyloseq_transform_aldex_clr <- function(physeq, variable = NULL, iter = 1){

  ## Extract OTU abundance table
  OTUS <- as.data.frame(phyloseq::otu_table(physeq))

  ## Transpose data
  trow <- phyloseq::taxa_are_rows(physeq)
  if(trow == FALSE){
    OTUS <- t(OTUS)
  }

  ## Extract meta-data
  if(is.null(variable)){
    ## Create dummy sample descriptor
    conds <- rep("A", phyloseq::nsamples(physeq))
    warning("Dummy sample descriptor was used. Please consider to use one of the groupping variables from the sample metadata.\n")
  } else {
    ## Or use the provided groupping variable
    conds <- as.character( phyloseq::get_variable(physeq, varName = variable) )
  }

  ## Generate Monte Carlo samples of the Dirichlet distribution for each sample
  ## Convert each instance using the centred log-ratio transform
  ## Rows should contain OTUs and columns should contain sequencing read counts
  CLRs <- ALDEx2::aldex.clr(reads = OTUS, conds = conds, mc.samples = iter, denom = "iqlr", verbose = TRUE)

  ## Extract CLR-transformed abundances
  exract_aldex_clr <- function(ald){
    ## ald = result of ALDEx2::aldex.clr

    ## Extract CLR-transformed Monte Carlo samples of the Dirichlet distribution for each sample.
    dtt <- ALDEx2::getMonteCarloInstances(ald)
    # dtt <- ald@analysisData                      # the same

    ## Combine MC instances for each sample
    res <- plyr::mlply(
      .data = data.frame(i = 1:ald@mc.samples),
      .fun = function(i){
        mcs <- plyr::llply(.data = dtt, .fun = function(z){ z[, i, drop=FALSE] })
        rez <- do.call(cbind, mcs)
        colnames(rez) <- names(mcs)
        return(rez)
      })

    ## List with the same number of transformed abundandance tables as the number of Monte Carlo samples used
    return(res)
  }

  ## Extract CLR-transformed abundances
  CLRs_ab <- exract_aldex_clr(CLRs)

  ## Transpose OTU tables
  if(trow == FALSE){
    CLRs_ab <- plyr::llply(.data = CLRs_ab, .fun = function(z){ t(z) })
  }

  ## Function to replace OTU table in phyloseq object with new table
  replace_otu_table <- function(phys, new_tab, trow = TRUE){
    phyloseq::otu_table(phys) <- phyloseq::otu_table(new_tab, taxa_are_rows = trow)
    return(phys)
  }

  ## Replace phyloseq tables
  if(iter == 1){
    ## Take only the first MC sample
    physeq_CLR <- replace_otu_table(physeq, CLRs_ab[[ 1 ]], trow)
  } else {
    ## Use the specified number of MC samples
    physeq_CLR <- list()
    for(i in 1:iter){
      physeq_CLR[[ i ]] <- replace_otu_table(physeq, CLRs_ab[[ i ]], trow)
    }
  }

  ## Add aldex.clr results as attributes to the phyloseq object
  ## Could be useful to check which OTUs were used in CLR-transformation as denominator
  attr(physeq_CLR, which = "ALDEx2") <- CLRs

  return(physeq_CLR)
}



#' @title Cumulative sum scaling (CSS) normalization of OTU abundance table.
#'
#' @param physeq A phyloseq-class object
#' @param norm Logical, return normalized counts
#' @param log Logical, apply a logarithmic transform (log2) to the normalized count data
#' @param ... Additional arguments will be passed to \code{\link{MRcounts}}
#'
#' @details Median scaling factor across samples will be used as default.
#' @return Phyloseq object with transformed counts in OTU table.
#' @export
#' @references
#' Paulson et al. Nature Methods 10, 1200???1202 (2013) doi:10.1038/nmeth.2658. http://www.nature.com/nmeth/journal/v10/n12/full/nmeth.2658.html
#' @seealso \code{\link{phyloseq_transform_vst_blind}}, \code{\link{phyloseq_transform_rlog_blind}}, \code{\link{rlog}}, \code{\link{varianceStabilizingTransformation}}
#' @examples
#' data("GlobalPatterns")
#' gp_tr <- phyloseq_transform_css(GlobalPatterns)
#' head(otu_table(gp_tr))
#'
#' # Ordination on CSS-transformed data
#' plot_ordination(physeq = gp_tr, ordination = ordinate(gp_tr, method = "PCoA", distance = "bray"), color = "SampleType")
#' # Ordination on raw data
#' plot_ordination(physeq = GlobalPatterns, ordination = ordinate(GlobalPatterns, method = "PCoA", distance = "bray"), color = "SampleType")
#'
phyloseq_transform_css <- function(physeq, norm = TRUE, log = TRUE, ...){

  # require(metagenomeSeq)

  MGS <- phyloseq::phyloseq_to_metagenomeSeq(physeq)

  # get the normalized count matrix
  otu_norm <- metagenomeSeq::MRcounts(MGS, norm = norm, log = log, ...)
  # exportMat(datt.norm, file = "tmp.txt")    # norm = TRUE, log = TRUE
  ## save sample statistics (sample scaling factor, quantile value, number of identified features and library size):
  # exportStats(dattM.sf, file = "tmp_stats.txt")

  # Substitue raw abundance to the css-normalized data
  physeq.tr <- physeq
  phyloseq::otu_table(physeq.tr) <- phyloseq::otu_table(otu_norm, taxa_are_rows = T)
  return(physeq.tr)
}


#' @title Variance stabilizing transformation (VST) of OTU abundance table.
#'
#' @param physeq A phyloseq-class object
#' @param dropneg Logical, replace negative transformed values with 0
#' @param dropmissing Logical, remove missing data
#' @param ... Not yet implemented
#@param ... Additional arguments will be passed to \code{\link{varianceStabilizingTransformation}
#'
#' @details For downstream analysis it could be better to use sample covariate information (blind = FALSE in \code{\link{varianceStabilizingTransformation}}).
#' @return Phyloseq object with transformed counts in OTU table.
#' @export
#' @seealso \code{\link{phyloseq_transform_rlog_blind}}, \code{\link{phyloseq_transform_css}}, \code{\link{varianceStabilizingTransformation}}, \code{\link{rlog}}
#' @examples
#'
#'
phyloseq_transform_vst_blind <- function(physeq, dropneg = F, dropmissing = T, ...){

  # require(DESeq2)

  ## Add dummy sample data (phyloseq_to_deseq2 doesn't work without sample_data)
  if(is.null( phyloseq::sample_data(physeq, errorIfNULL = F) )){
    smpdat_nul <- TRUE
    smpdat <- data.frame(TMP = rep(1, times = phyloseq::nsamples(physeq)))
    rownames(smpdat) <- phyloseq::sample_names(physeq)
    phyloseq::sample_data(physeq) <- smpdat
  }

  ## Convert phyloseq data to DESeq2 object
  dsc <- phyloseq::phyloseq_to_deseq2(physeq, design = formula(~ 1))
  # otu_norm <- varianceStabilizingTransformation(dsc, blind = T, ...)
  # otu_norm <- assay(otu_norm)

  ## Estimate the size factors and dispersions
  dsc_tr <- try( DESeq2::estimateSizeFactors(dsc) )
  ## If it fails (probably because of excessive zeros) try a workaround
  ## see https://github.com/joey711/phyloseq/issues/445
  if("try-error" %in% class(dsc_tr)){

    ## Zero-tolerant version of geometric mean
    gm_mean <- function(x, na.rm=TRUE){ exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x)) }

    geoMeans <- apply(counts(dsc), 1, gm_mean)     # geometric means of the counts
    dsc <- DESeq2::estimateSizeFactors(dsc, geoMeans = geoMeans)
  } else {
    dsc <- dsc_tr    # if everithing is OK
    rm(dsc_tr)
  }

  ## Estimate dispersions
  # dsc <- DESeq2::estimateDispersions(dsc)
  dsc <- DESeq(dsc, fitType="local")

  ## Extract variance stabilized data
  otu_norm <- DESeq2::getVarianceStabilizedData(dsc)

  ## Negative values probably correspond to "less than one count"
  ## Set to zero all values less than zero
  if(dropneg == TRUE){
    otu_norm[otu_norm < 0.0] <- 0.0
  }

  ## Substitue raw abundance to the variance stabilized data
  physeq.tr <- physeq
  phyloseq::otu_table(physeq.tr) <- phyloseq::otu_table(otu_norm, taxa_are_rows = T)

  ## Remove missing OTUs
  if(dropneg == TRUE & dropmissing == TRUE){
    physeq.tr <- phyloseq::prune_taxa(phyloseq::taxa_sums(physeq.tr) > 0, physeq.tr)
  }

  ## Remove dummy sample data if present
  if(smpdat_nul == TRUE){
    pp@sam_data <- NULL
  }

  return(physeq.tr)
}


#' @title Regularized-log (rlog) transformation of OTU abundance table.
#'
#' @param physeq A phyloseq-class object
#' @param dropneg Logical, replace negative transformed values with 0
#' @param dropmissing Logical, remove missing data
#' @param ... Additional arguments will be passed to \code{\link{rlogTransformation}}
#'
#' @details rlog transformation (this function) is preferable to the vst (\code{\link{phyloseq_transform_vst_blind}}) if the size factors vary widely.
#' @return Phyloseq object with transformed counts in OTU table.
#' @export
#' @seealso \code{\link{phyloseq_transform_vst_blind}}, \code{\link{phyloseq_transform_css}}, \code{\link{rlog}}, \code{\link{varianceStabilizingTransformation}}
#' @examples
#'
#'
phyloseq_transform_rlog_blind <- function(physeq, dropneg = F, dropmissing = T, ...){
  require(DESeq2)

  # Add dummy sample data (phyloseq_to_deseq2 doesn't work without sample_data)
  if(is.null( sample_data(physeq, errorIfNULL = F) )){
    smpdat_nul <- TRUE
    smpdat <- data.frame(TMP = rep(1, times = nsamples(physeq)))
    rownames(smpdat) <- sample_names(physeq)
    sample_data(physeq) <- smpdat
  }

  dsc <- phyloseq_to_deseq2(physeq, design = formula(~ 1))
  otu_norm <- rlogTransformation(dsc, blind = T, ...)
  otu_norm <- assay(otu_norm)

  # Set to zero all values less than zero
  if(dropneg == TRUE){
    otu_norm[otu_norm < 0.0] <- 0.0
  }

  # Substitue raw abundance to the rlog transformed data
  physeq.tr <- physeq
  otu_table(physeq.tr) <- otu_table(otu_norm, taxa_are_rows = T)

  # Remove missing OTUs
  if(dropneg == TRUE & dropmissing == TRUE){
    physeq.tr <- prune_taxa(taxa_sums(physeq.tr) > 0, physeq.tr)
  }

  # Remove dummy sample data if present
  if(smpdat_nul == TRUE){
    pp@sam_data <- NULL
  }

  return(physeq.tr)
}



#' @title Log-transformation of OTU abundance table.
#'
#' @param physeq A phyloseq-class object
#' @param ... Additional arguments (e.g., "logbase" for the logarithm base) will be passed to \code{\link{decostand}}
#'
#' @details
#' Logarithmic transformation as suggested by Anderson et al. (2006) will be applied to the OTU table in phyloseq-object. First, non-integer data will be divided by smallest positive value. Second, log(x) + 1 for x > 0.
#' Default value of the logarithm base ("logbase" parameter) is 2.
#' This function is a wrapper to \code{\link{decostand}} in vegan.
#' @return Phyloseq object with log-transformed counts in OTU table.
#' @export
#' @seealso \code{\link{decostand}}
#' @references  Anderson MJ, Ellingsen KE and McArdle BH. (2006) Multivariate dispersion as a measure of beta diversity. Ecology Letters 9, P. 683???693.
#'
#' @examples
#' data(enterotype)
#' physeq_transform_anderson_log(enterotype)
#'
physeq_transform_anderson_log <- function(physeq, ...){
  require(vegan)
  otus <- as.data.frame( otu_table(physeq) )
  otu_log <- decostand(otus, method = "log", ...)
  rownames(otu_log) <- taxa_names(physeq)
  otu_table(physeq) <- otu_table(otu_log, taxa_are_rows = T)
  return(physeq)
}
