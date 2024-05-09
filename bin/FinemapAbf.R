#!/usr/bin/env Rscript

library(data.table)
library(coloc)
library(arrow)
library(stringr)
library(dplyr)
library(argparse)

parser <- ArgumentParser()
parser$add_argument("--pheno", 
                    help = "Phenotype name")

parser$add_argument("--maf", type = "double", default = 0.01, 
                    help = "MAF threshold for the SNPs to include [default %(default)s]")
parser$add_argument("--info", type = "double", default = 0.8, 
                    help = "INFO score threshold for the SNPs to include [default %(default)s]")
parser$add_argument("--pvalue", type = "double", default = 5e-8, 
                    help = "P-value threshold for the locus definition [default %(default)s]")
parser$add_argument("--window", type = "integer", default = 1000000, 
                    help = "Window size for defining the locus boundaries [default %(default)s]")
parser$add_argument("--prior", type = "double", default =  1e-04, 
                    help = "Prior used for ABF finemapping [default %(default)s]")

parser$add_argument("--input",
                    help = "Path to input parquet file")
parser$add_argument("--reference",
                    help = "Path to SNP reference parquet file")

parser$add_argument("--outputpip",
                    help = "Path to output file with per-snp PIPs")
parser$add_argument("--outputsummary",
                    help = "Path to output file with summary over loci")

args <- parser$parse_args()

# args <- c(
#   "/Users/urmovosa/Documents/projects/2023/sampo_pipelines/HeritabilitySampoNF/tests/input/sumstats/results_concat_H52.1.parquet.snappy",
#   "/Users/urmovosa/Documents/projects/2023/sampo_pipelines/temp_finemap_abf/temp/snp_reference.parquet",
#   0.01,
#   0.8,
#   5e-8,
#   1000000,
#   1e-04,
#   "/Users/urmovosa/Documents/projects/2023/sampo_pipelines/temp_finemap_abf/temp/finemapping_results.txt"
# )

# functions
ZtoP <- function(Z, largeZ = FALSE, log10P = TRUE) {
  if (!is.numeric(Z)) {
    message("Some of the Z-scores are not numbers! Please check why!")
    message("Converting the non-numeric vector to numeric vector.")
    Z <- as.numeric(Z)
  }

  if (largeZ == TRUE) {
    P <- log(2) + pnorm(abs(Z), lower.tail = FALSE, log.p = TRUE)

    if (largeZ == TRUE & log10P == TRUE) {
      P <- -(P * log10(exp(1)))
    }
  } else {
    P <- 2 * pnorm(abs(Z), lower.tail = FALSE)

    if (min(P) == 0) {
      P[P == 0] <- .Machine$double.xmin
      message("Some Z-score indicates very significant effect and P-value is truncated on 2.22e-308. If relevant, consider using largeZ = TRUE argument and logarithmed P-values instead.")
    }
  }

  return(P)
}

IdentifyLeadSNPs <- function(data,
                             window = 1000000,
                             Pthresh = 5e-8,
                             snp_id_col = "snp",
                             snp_chr_col = "chr",
                             snp_pos_col = "pos",
                             eff_all_col = "ea",
                             other_all_col = "nea",
                             beta_col = "beta",
                             se_col = "se",
                             p_col = NULL,
                             loci = FALSE,
                             verbose = TRUE) {
  if (is.null(p_col) & !is.null(beta_col) & !is.null(se_col)) {
    data <- data.table(
      SNP = data[[snp_id_col]],
      chr = data[[snp_chr_col]],
      pos = data[[snp_pos_col]],
      ea = data[[eff_all_col]],
      nea = data[[other_all_col]],
      beta = data[[beta_col]],
      se = data[[se_col]]
    )

    data$P <- ZtoP(data$beta / data$se)
    data$sig_ind <- data$beta / data$se
  } else if (!is.null(p_col) & is.null(beta_col) & is.null(se_col)) {
    data <- data.table(
      SNP = data[[snp_id_col]],
      chr = data[[snp_chr_col]],
      pos = data[[snp_pos_col]],
      ea = data[[eff_all_col]],
      nea = data[[other_all_col]],
      beta = data[[beta_col]],
      se = data[[se_col]],
      P = data[[p_col]]
    )

    data$sig_ind <- -log10(data$P)
  } else {
    message("There is no beta/se nor P-value in the data")
    stop()
  }

  data_f <- data[data$P < as.numeric(Pthresh), ]
  message(paste("Sig. rows remaining:", nrow(data_f)))
  if (nrow(data_f) == 0) {
    message(paste0("No locus meets significance threshold ", Pthresh, "!"))
    return(data[-c(1:nrow(data)), ])
  }

  # Iteratively identify most significant SNP, and remove all other SNPs in the window
  res <- data_f[-c(1:nrow(data_f)), ]

  while (min(data_f$P) <= Pthresh) {
    lead_snp <- data_f[abs(data_f$sig_ind) == max(abs(data_f$sig_ind)), ]
    if (nrow(lead_snp) > 1) {
      lead_snp <- lead_snp[1, ]
    }
    res <- rbind(res, lead_snp)
    data_f <- data_f[!(data_f$chr == lead_snp$chr & data_f$pos > lead_snp$pos - window & data_f$pos < lead_snp$pos + window), ]
    if (isTRUE(verbose)) {
      message(paste0("Added: ", lead_snp$SNP, " | ", lead_snp$chr, ":", lead_snp$pos))
    }
    if (nrow(data_f) == 0) {
      break
    }
  }

  if (isTRUE(verbose)) {
    
    message(paste(nrow(res), "loci found!"))
    
  }
  
  res <- res[order(chr, pos), -ncol(res), with = FALSE]

  if (isTRUE(loci)) {
    if (isTRUE(verbose)) {
    message("Extracting full loci...")
    }
    temp_res <- data.table(locus = 1, res, lead = TRUE)
    temp_res <- temp_res[-c(1:nrow(temp_res)), ]
    for (i in 1:nrow(res)) {
      locus <- data[(data$chr == res[i]$chr & data$pos > res[i]$pos - window & data$pos < res[i]$pos + window), -ncol(data), with = FALSE]
      locus <- data.table(locus = i, locus, lead = FALSE)
      locus[SNP == res[i]$SNP]$lead <- TRUE
      temp_res <- rbind(temp_res, locus)
    }
    res <- temp_res
    if (isTRUE(verbose)) {
    message("Extracting full loci...done")
    }
  }

  return(res)
}


# Read the input parquet file
message("Reading input...")
sumstats <- read_parquet(args$input)
# Convert to a data.table and calculate P from LOG10P
sumstats <- data.table(SNP = sumstats$ID, BETA = sumstats$BETA, SE = sumstats$SE, P = sumstats$LOG10P)
message("Reading input...done!")
sumstats$P <- 10^(-sumstats$P)


if (nrow(sumstats[P < as.numeric(args$pvalue)]) == 0){



  abi <- data.table(locus = NA, 
                    SNP = NA, 
                    chr = NA, 
                    pos = NA, 
                    ea = NA, 
                    nea = NA, 
                    beta = NA, 
                    se = NA, 
                    P = NA, 
                    lead = NA, 
                    PP = NA, 
                    CS_95percent = NA, 
                    CS_99percent = NA, 
                    cumPP = NA)[-1]
  
  fwrite(abi, args$outputpip, sep = "\t")

  abi2 <- data.table(
            locus = NA,
            LeadVariant = NA,
            chr = NA,
            pos = NA,
            ea = NA,
            nea = NA,
            beta = NA,
            se = NA,
            P = NA,
            NrVariants = NA,
            CS_95_size = NA,
            CS_95_CI_start = NA,
            CS_95_CI_end = NA,
            CS_99_size = NA,
            CS_99_CI_start = NA,
            CS_99_CI_end = NA,
            LargestPipVariant = NA,
            LargestPipVariantChr = NA,
            LargestPipVariantPos = NA,
            LargestPipVariantEa = NA,
            LargestPipVariantNea = NA,
            LargestPipVariantBeta = NA,
            LargestPipVariantSe = NA,
            LargestPipVariantP = NA,
            LargestPip = NA
  )[-1]
  
  fwrite(abi2, args$outputsummary, sep = "\t")
  message("No significant loci!")
  quit(save = "no", status = 0, runLast = FALSE)
}

message("Indexing...")
setkey(sumstats, "SNP")
message("Indexed...done!")

# Add needed info
message("Reading referece...")
snpref <- read_parquet(args$reference)
colnames(snpref)[1] <- "SNP"
snpref <- as.data.table(snpref)
message("Reading reference...done!")
message("Indexing...")
setkey(snpref, SNP)
message("Indexing...done!")
# Merge
sumstats <- merge(sumstats, snpref, by = "SNP")

# Filter on INFO and MAF
sumstats <- sumstats[A1FREQ > as.numeric(args$maf) & A1FREQ < 1 - as.numeric(args$maf) &
  INFO > as.numeric(args$info)]

if (nrow(sumstats[P < as.numeric(args$pvalue)]) == 0){
  
  abi <- data.table(locus = NA, 
                    SNP = NA, 
                    chr = NA, 
                    pos = NA, 
                    ea = NA, 
                    nea = NA, 
                    beta = NA, 
                    se = NA, 
                    P = NA, 
                    lead = NA, 
                    PP = NA, 
                    CS_95percent = NA, 
                    CS_99percent = NA, 
                    cumPP = NA)[-1]
  
  fwrite(abi, args$outputpip, sep = "\t")

  abi2 <- data.table(
          locus = NA,
          LeadVariant = NA,
          chr = NA,
          pos = NA,
          ea = NA,
          nea = NA,
          beta = NA,
          se = NA,
          P = NA,
          NrVariants = NA,
          CS_95_size = NA,
          CS_95_CI_start = NA,
          CS_95_CI_end = NA,
          CS_99_size = NA,
          CS_99_CI_start = NA,
          CS_99_CI_end = NA,
          LargestPipVariant = NA,
          LargestPipVariantChr = NA,
          LargestPipVariantPos = NA,
          LargestPipVariantEa = NA,
          LargestPipVariantNea = NA,
          LargestPipVariantBeta = NA,
          LargestPipVariantSe = NA,
          LargestPipVariantP = NA,
          LargestPip = NA
)[-1]
  
  fwrite(abi2, args$outputsummary, sep = "\t")
  message("No significant loci!")
  quit(save = "no", status = 0, runLast = FALSE)
}

# Find loci
message("Finding loci...")
loci <- IdentifyLeadSNPs(sumstats,
  window = as.numeric(args$window),
  Pthresh = as.numeric(args$pvalue),
  snp_id_col = "SNP",
  snp_chr_col = "CHR",
  snp_pos_col = "POS",
  eff_all_col = "ALLELE1",
  other_all_col = "ALLELE0",
  beta_col = "BETA",
  se_col = "SE",
  p_col = NULL,
  loci = TRUE
)
message("Finding loci...done")

message("Running fine-mapping iteratively over loci...")
nr_loci <- length(unique(loci$locus))
res <- data.table(loci[-c(1:nrow(loci)), ], PP = NA, CS_95percent = NA, CS_99percent = NA, cumPP = NA)[-1, ]
for (l in unique(loci$locus)) {
  abi <- loci[locus == l, ]
  abi <- abi[order(pos)]
  inp <- list(snp = abi$SNP, beta = abi$beta, varbeta = abi$se * abi$se, type = "cc")

  finemapped_res <- finemap.abf(inp, p1 = as.numeric(args$prior))
  finemapped_res <- finemapped_res[, c(5, 7)]
  colnames(finemapped_res) <- c("SNP", "PP")
  abi <- merge(abi, finemapped_res, by = "SNP")

  # calculate credible sets
  abi$CS_95percent <- "no"
  abi$CS_99percent <- "no"

  abi <- abi[order(PP, decreasing = TRUE)]
  abi$cumPP <- cumsum(abi$PP)

  abi$CS_95percent <- "no"
  abi$CS_99percent <- "no"
    abi$CS_95percent[1] <- "yes"
    for (i in 1:(nrow(abi) - 1)) {
      if (abi$cumPP[i] < 0.95) {
        abi$CS_95percent[i + 1] <- "yes"
      } else {
        abi$CS_95percent[i + 1] <- "no"
      }
    }

    abi$CS_99percent[1] <- "yes"
    for (i in 1:(nrow(abi) - 1)) {
      if (abi$cumPP[i] < 0.99) {
        abi$CS_99percent[i + 1] <- "yes"
      } else {
        abi$CS_99percent[i + 1] <- "no"
      }
    }

  res <- rbind(res, abi)

  message(paste0("Locus ", l, " fine-mapped (", round((as.numeric(l) / nr_loci) * 100, digits = 1)), "%)")
}

message("Running fine-mapping iteratively over loci...done")

res_summary <- res %>% 
  group_by(locus) %>%
  mutate(
    NrVariants = length(SNP),
    CS_95_size = length(SNP[CS_95percent == "yes"]),
    CS_95_CI_start = min(pos[CS_95percent == "yes"]),
    CS_95_CI_end = max(pos[CS_95percent == "yes"]),
    CS_99_size = length(SNP[CS_99percent == "yes"]),
    CS_99_CI_start = min(pos[CS_99percent == "yes"]),
    CS_99_CI_end = max(pos[CS_99percent == "yes"])) %>% 
    summarise(LeadVariant = SNP[lead == TRUE][1],
            chr = chr[lead == TRUE][1],
            pos = pos[lead == TRUE][1],
            ea = ea[lead == TRUE][1],
            nea = nea[lead == TRUE][1],
            beta = beta[lead == TRUE][1],
            se = se[lead == TRUE][1],
            P = P[lead == TRUE][1],
            NrVariants = unique(NrVariants),
            CS_95_size = unique(CS_95_size),
            CS_95_CI_start = unique(CS_95_CI_start),
            CS_95_CI_end = unique(CS_95_CI_end),
            CS_99_size = unique(CS_99_size),
            CS_99_CI_start = unique(CS_99_CI_start),
            CS_99_CI_end = unique(CS_99_CI_end),
            LargestPipVariant = SNP[PP == max(PP)][1],
            LargestPipVariantChr = chr[PP == max(PP)][1],
            LargestPipVariantPos = pos[PP == max(PP)][1],
            LargestPipVariantEa = ea[PP == max(PP)][1],
            LargestPipVariantNea = nea[PP == max(PP)][1],
            LargestPipVariantBeta = beta[PP == max(PP)][1],
            LargestPipVariantSe = se[PP == max(PP)][1],
            LargestPipVariantP = P[PP == max(PP)][1],
            LargestPip = PP[PP == max(PP)][1]
            )

res_summary <- data.table(Phenotype = args$pheno, res_summary)
res <- data.table(Phenotype = args$pheno, res)

fwrite(res, args$outputpip, sep = "\t")
fwrite(res_summary, args$outputsummary, sep = "\t")

