#' @export
Fold_Change <- function(data,
               asterisk = "t-test",
               reverse = F,
               FC_range = c(1,1),
               pval_intercept = 0.05,
               directory = getwd(),
               subdirectory = "Foldchange") {
  ifelse(!dir.exists(file.path(directory, subdirectory)), dir.create(file.path(directory, subdirectory)), FALSE)
  data[["Data_renamed"]] <-
    data[["Data_renamed"]] %>% plyr::mutate(ZZZZ = data[["Data_renamed"]][, 2])
  data[["Data_renamed"]] <- data[["Data_renamed"]][, c(-1, -2)]
  data[["Data_renamed"]][, 1:(ncol(data[["Data_renamed"]]) - 1)] <-
    sapply(
      data[["Data_renamed"]][, 1:(ncol(data[["Data_renamed"]]) - 1)],
      function(x) {
        as.numeric(x)
      }
    )
  NAMES <- colnames(data[["Data"]][3:ncol(data[["Data"]])])
  data[["Data_renamed"]]$ZZZZ <-
    as.factor(data[["Data_renamed"]]$ZZZZ)
  Comb <-
    gtools::combinations(length(unique(data[["Data"]]$Group)), 2, unique(data[["Data"]]$Group))

  {
    if (asterisk == "Scheffe") {
      p_val_data <- data[["Anova_PostHoc"]]
    } else if (asterisk == "t-test") {
      p_val_data <- data[["t_test"]]
    } else if (asterisk == "u-test") {
      p_val_data <- data[["u_test"]]
    } else if (asterisk == "Dunn") {
      p_val_data <- data[["Dunn"]]
    } else {
      print("Wrong asterisk input must be one of (Dunn,Scheffe,u_test,t_test)")
    }
    }

  {
    if (length(unique(data[["Data"]]$Group)) != 2) {
      for (a in 1:nrow(p_val_data)) {
        for (b in 1:ncol(p_val_data)) {
          if (is.nan(p_val_data[a, b]) == T) {
            p_val_data[a, b] <- 1
          }
        }
      }
    } else if (length(unique(data[["Data"]]$Group)) == 2) {
      for (b in 1:length(p_val_data)) {
        if (is.nan(p_val_data[b]) == T) {
          p_val_data[b] <- 1
        }
      }
    }
  }

  for (a in 1:nrow(Comb)) {{
    if (reverse == T) {
      G1 <- Comb[a, 2]
      G2 <- Comb[a, 1]
    } else if (reverse == F) {
      G1 <- Comb[a, 1]
      G2 <- Comb[a, 2]
    }
    data[["Data"]][, 3:ncol(data[["Data"]])] <-
      apply(data[["Data"]][, 3:ncol(data[["Data"]])], 2, function(x) {
        as.numeric(x)
      })
    Part_dat_1 <- data[["Data"]] %>% dplyr::filter(Group == G1)
    Part_dat_2 <- data[["Data"]] %>% dplyr::filter(Group == G2)
    Volc_Dat <- matrix(nrow = (ncol(data[["Data"]]) - 2), ncol = 3)
  }
    for (i in 3:ncol(Part_dat_1)) {
      Volc_Dat[(i - 2), 1] <-
        mean(Part_dat_1[, i]) / mean(Part_dat_2[, i])
    }
    if (length(unique(data[["Data"]]$Group)) > 2) {
      if (asterisk == "Scheffe") {{ if (paste0(Comb[a, 1], "-", Comb[a, 2], "___", "ANO_posthoc") %in% colnames(p_val_data) == T) {
        part_p_val <-
          p_val_data[, paste0(Comb[a, 1], "-", Comb[a, 2], "___", "ANO_posthoc")]
      } else {
        part_p_val <-
          p_val_data[, paste0(Comb[a, 2], "-", Comb[a, 1], "___", "ANO_posthoc")]
      } }} else if (asterisk == "Dunn") {{ if (paste0(Comb[a, 1], " - ", Comb[a, 2], "___", "Kru_posthoc(Dunn)") %in% colnames(p_val_data) == T) {
        part_p_val <-
          p_val_data[, paste0(Comb[a, 1], " - ", Comb[a, 2], "___", "Kru_posthoc(Dunn)")]
      } else {
        part_p_val <-
          p_val_data[, paste0(Comb[a, 2], " - ", Comb[a, 1], "___", "Kru_posthoc(Dunn)")]
      } }} else {{ if (paste0(Comb[a, 1], "-", Comb[a, 2], "___", asterisk) %in% colnames(p_val_data) == T) {
        part_p_val <-
          p_val_data[, paste0(Comb[a, 1], "-", Comb[a, 2], "___", asterisk)]
      } else {
        part_p_val <-
          p_val_data[, paste0(Comb[a, 2], "-", Comb[a, 1], "___", asterisk)]
      } }}
    } else if (length(unique(data[["Data"]]$Group)) == 2) {{ part_p_val <- p_val_data }}

    Volc_Dat[, 2] <- part_p_val
    rownames(Volc_Dat) <-
      colnames(data[["Data"]])[3:ncol(data[["Data"]])]
    Volc_Dat <- as.data.frame(Volc_Dat)
    colnames(Volc_Dat) <- c("FoldChange", "pvalue", "Direction")
    Volc_Dat$Direction <- "Not Significant"
    Volc_Dat$Direction[Volc_Dat$FoldChange > FC_range[2] &
                         Volc_Dat$pvalue < pval_intercept] <- "UP"
    Volc_Dat$Direction[Volc_Dat$FoldChange < FC_range[1] &
                         Volc_Dat$pvalue < pval_intercept] <- "DOWN"

    write.csv(Volc_Dat, paste0(directory,"/", subdirectory, "/", G1, " divided by ", G2, ".csv", collapse = ""), row.names= TRUE )
  }
  return("Fold change files were saved")
}
