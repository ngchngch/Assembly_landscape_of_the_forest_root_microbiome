library(dplyr)






#########################################################################
save.dir <- "02_06_02_02_basincomposition_freq_summary"
dir.create(save.dir)

fb <- "Fungi"

true_sa <- readRDS(sprintf("%s/%s/runSA_ST_full_%s_Family.rds",dir_02_06,fb,fb))
true_ssid <- readRDS(sprintf("%s/%s/graph_param_SStab_full_%s_Family.rds",dir_02_06,fb,fb))
ls <- list.files(dir_02_06_02,pattern = sprintf("basincomposition_freq_%s",fb),
           recursive = TRUE, full.names = TRUE)

df <- do.call(rbind,
        lapply(ls, readRDS))

#df <- rbind(as,readRDS("/Users/ngch/Desktop/Now_using/data/sugaFB_ELA/revise_251212/output_spacom/result_for_revise/basincomposition_freq_Fungi_boot52.rds"))
info <- cbind(new_ssid=true_ssid[df[,2],"rename_SS"],df[,1:3])
freq <- apply(df[,4:ncol(df)],2,as.numeric)

colnames(freq) <- rownames(true_sa[[1]])

smry <- cbind(info,as.data.frame(freq)) %>%
  group_by(new_ssid,plant, ssid) %>%
  summarise(across(where(is.numeric), mean))%>%
  as.data.frame()

saveRDS(smry,sprintf("%s/basincomposition_freq_%s_summary.rds",save.dir,fb))


fb <- "Prokaryota"

true_sa <- readRDS(sprintf("%s/%s/runSA_ST_full_%s_Family.rds",dir_02_06,fb,fb))
true_ssid <- readRDS(sprintf("%s/%s/graph_param_SStab_full_%s_Family.rds",dir_02_06,fb,fb))
ls <- list.files(dir_02_06_02,pattern = sprintf("basincomposition_freq_%s",fb),
                 recursive = TRUE, full.names = TRUE)

df <- do.call(rbind,
              lapply(ls, readRDS))

#df <- rbind(as,readRDS("/Users/ngch/Desktop/Now_using/data/sugaFB_ELA/revise_251212/output_spacom/result_for_revise/basincomposition_freq_Fungi_boot52.rds"))
info <- cbind(new_ssid=true_ssid[match(df[,2], rownames(true_ssid)),"rename_SS"],df[,1:3])
freq <- apply(df[,4:ncol(df)],2,as.numeric)

colnames(freq) <- rownames(true_sa[[1]])

smry <- cbind(info,as.data.frame(freq))[!is.na(info[,"new_ssid"]),] %>%
  group_by(new_ssid,plant, ssid) %>%
  summarise(across(where(is.numeric), mean))%>%
  as.data.frame()

saveRDS(smry,sprintf("%s/basincomposition_freq_%s_summary.rds",save.dir,fb))