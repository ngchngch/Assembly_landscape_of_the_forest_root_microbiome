library(ggplot2)
library(ggtext)
library(patchwork)

dir_03_02_05_02 <- "Output_supercomputer/03_10_graphics_states_flow_flow_spl_250508"

####
save.dir <- "Output/03_02_05_03_graphics_boot_zdL_zdE_keystone"
dir.create(save.dir)

########
fb_mat <- cbind(taxa=c("Oidiodendron",
                  "Cladophialophora",
                  "Candidatus Udaeobacter",
                  "Meliniomyces"),
                fb=c("Fungal landscape<br>(*Acer*)","Fungal landscape<br>(*Populus*)",
                     "Prokaryotic landscape<br>(*Juglans*)","Prokaryotic landscape<br>(*Pinus*)"))

boot_res <- readRDS(sprintf("%s/Zconv_ELA_withRA_4step_boot_summary.rds",dir_03_02_05_02))
boot_res2 <- merge(boot_res,fb_mat,by="taxa")
boot_res2$taxa2 <- factor(sprintf("*%s*",boot_res2$taxa),levels=sprintf("*%s*",
                                                             c("Oidiodendron",
                                                               "Cladophialophora",
                                                               "Candidatus Udaeobacter",
                                                               "Meliniomyces")))

g_zdL <- ggplot(boot_res2[boot_res2$ra=="median",],aes(x=z_dland))+
  geom_vline(xintercept = 0,linetype="dashed",color="gray30")+
  geom_histogram()+
  facet_wrap(~taxa2,ncol=1)+
  labs(x="Z-standardized *\u0394landscape*",y="Count")+
  theme_bw()+
  theme(strip.background = element_rect(fill="white"),
        aspect.ratio = 0.5,
        strip.text = element_markdown(size=12),
        axis.text = element_text(size=10),
        axis.title.y = element_text(size=12),
        axis.title.x = element_markdown(size=12))

g_zdE <-ggplot(boot_res2[boot_res2$ra=="median",],aes(x=z_deven))+
  geom_vline(xintercept = 0,linetype="dashed",color="gray30")+
  geom_histogram()+
  facet_wrap(~taxa2,ncol=1)+
  labs(x="Z-standardized *\u0394evenness*",y="Count")+
  theme_bw()+
  theme(strip.background = element_rect(fill="white"),
        aspect.ratio = 0.5,
        strip.text = element_markdown(size=12),
        axis.text = element_text(size=10),
        axis.title.y = element_text(size=12),
        axis.title.x = element_markdown(size=12))


g <- g_zdL + g_zdE

ggsave(sprintf("%s/histogram_zdL_zdE.pdf",save.dir),plot=g,
       width=6,height=7,device = cairo_pdf)
