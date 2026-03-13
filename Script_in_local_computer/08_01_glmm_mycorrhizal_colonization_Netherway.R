setwd("/Users/ngch/Desktop/Now_using/data/sugaFB_ELA/revise_251212/output_spacom/result_for_revise")
# %% loading library and functions
library(tidyverse)
library(dplyr)
library(compositions)
library(glmmTMB)

save.dir <- "analysis_in_local_computer/Output/08_01_glmm_mycorrhizal_colonization_Netherway"
dir.create(save.dir)

#######
df_f <- read.csv("analysis_in_local_computer/SourceData_08/root_FungiOTU_matrix.csv",row.names = 1)
tx_f <- read.csv("analysis_in_local_computer/SourceData_08/FungiOTU_taxonomy.csv",row.names = 1)
df_b <- read.csv("analysis_in_local_computer/SourceData_08/root_BacteriaOTU_matrix.csv",row.names = 1)
tx_b <- read.csv("analysis_in_local_computer/SourceData_08/BacteriaOTU_taxonomy.csv",row.names = 1)

info <- read.csv("analysis_in_local_computer/SourceData_08/env_metadata.csv",row.names = 1)
info$Row.names <- rownames(info)
#######
dir_03_02 <- "03_02_Zconv_randELA_withRA_fixP"
# read data
z_f <- readRDS(sprintf("%s/zvalue_Fungi.rds",dir_03_02))
z_p <- readRDS(sprintf("%s/zvalue_Prokaryote.rds",dir_03_02))

z_f$p_land_BH <- p.adjust(c(z_f$p_land,z_p$p_land), method = "BH")[1:nrow(z_f)]

z_f$p_even_BH <- p.adjust(c(z_f$p_even,z_p$p_even), method = "BH")[1:nrow(z_f)]

df2 <- z_f %>%
  filter(ra=="median")　%>%
  group_by(plant) %>%
  mutate(rank_desc = rank(-z_land, ties.method = "first")) %>%
  filter(rank_desc <= 10) %>%
  arrange(rank_desc, .by_group = TRUE) %>%
  select(Taxa, plant, z_land,p_land_BH,z_even,p_even_BH,rank_desc) %>%
  as.data.frame()

write.csv(df2, sprintf("%s/top10_taxa.csv",save.dir),row.names = FALSE)

target <- unique(df2$Taxa)

dff_g <- Taxa.mat(t(df_f),tx_f,"Genus")

clr_dff <- clr(dff_g+1)

dfb_g <- Taxa.mat(t(df_b),tx_b,"Genus")

clr_dfb <- clr(dfb_g+1)

merge_df <- merge(clr_dff[,which(colnames(dff_g) %in% target)],
                  clr_dfb[,which(colnames(dfb_g) %in% target)],
                  by="row.names")　%>%
  merge(info, by.x="Row.names") 

merge_df2 <- merge_df %>%
  filter(Tree_species %in% c("Alnus_glutinosa", "Alnus_incana"))

ntarget <- ncol(merge_df2)-ncol(info)

saveRDS(merge_df2, sprintf("%s/data_for_glmm.rds",save.dir))
total <- 100

res <- NULL

set.seed(1234)
for(i in 1:ntarget){
  merge_df2$abundance <- merge_df2[,i+1]
  
  model_am <- glmmTMB::glmmTMB(
    cbind(AM_colonisation, total - AM_colonisation) ~ Tree_species+scale(Soil_pH)+scale(abundance) + (1 | Site),
    family = binomial(link = "logit"),
    data = merge_df2
  )
  
  sm_am <- summary(model_am)
  
  model_ecm <- glmmTMB::glmmTMB(
    EcM_colonisation/total ~ Tree_species+scale(Soil_moisture)+scale(abundance) + (1 | Site),
    family = beta_family(link = "logit"),
    data = merge_df2
  )
  
  sm_ecm <- summary(model_ecm)
  
  res <- rbind(res,rbind(cbind(Var="AM_colonization",
                               target=colnames(merge_df)[i+1],
                               as.data.frame(matrix(sm_am$coefficients$cond["scale(abundance)",],nrow=1))
              ),
        cbind(Var="EcM_colonization",
              target=colnames(merge_df)[i+1],
              as.data.frame(matrix(sm_ecm$coefficients$cond["scale(abundance)",],nrow=1))
        )
        ))
  
}

colnames(res)[3:6] <- c("Estimate","se","zvalue","Pr")

res$p_BH <- NA
res[res$Var=="AM_colonization","p_BH"] <- p.adjust(res[res$Var=="AM_colonization","Pr"], method = "BH")
res[res$Var=="EcM_colonization","p_BH"] <- p.adjust(res[res$Var=="EcM_colonization","Pr"], method = "BH")

res$sig <- res$p_BH < 0.05

res[grep("Rhizobium",res$target),"target"] <- "Rhizobium"

write.csv(res, sprintf("%s/glmm_results.csv",save.dir),row.names = FALSE)

res$target2 <- factor(sprintf("*%s*",res$target),
                     levels=sprintf("*%s*",res[res$Var=="AM_colonization","target"])[order(res[res$Var=="AM_colonization","Estimate"],decreasing = TRUE)])


res$Var2 <- ifelse(res$Var=="AM_colonization",
                   "Root colonization rates\n(AM)",
                   "Root colonization rates\n(EcM)")

g <- ggplot(res,aes(x=target2,y=Estimate))+
  geom_hline(yintercept = 0)+
  geom_errorbar(aes(ymin=Estimate-se,ymax=Estimate+se))+
  geom_point(aes(fill=sig),shape=21)+
  labs(y="Standardized regression coefficients\n(Abundance of each highlighted genus [CLR-transformed])",
       x="Fungal/Proaryotic genera",
       fill="*P*(FDR) < 0.05")+
  facet_wrap(~Var2,scales = "free")+
  scale_fill_manual(values=c("gray","red"))+
  theme_bw()+
  theme(aspect.ratio = 2.5,
        strip.text = element_text(size=14),
        axis.text.y = element_markdown(size=12),
        axis.text.x = element_text(size=12),
        axis.title = element_text(size=14),
        legend.title = element_markdown(size=14),
        legend.text = element_text(size=12))+
  coord_flip()



g

ggsave(plot = g,
       sprintf("%s/clonizationrates_GLMM_results.pdf",save.dir),
       h=8,w=10)

ggsave(plot = g,
       sprintf("%s/clonizationrates_GLMM_results.png",save.dir),
       h=8,w=10,dpi=300)

