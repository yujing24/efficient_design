rm(list=ls())
library(tidyr)
library(tidyverse)
library(paletteer)
library(RColorBrewer)
library(ggplot2)

{
  # display.brewer.all()
  # display.brewer.pal(2, 'Dark2')
  # brewer.pal(2, 'Dark2')
}

# Consider a simple case 
# The treatment assignment for the RPCT data is completely at random 
# There is no unmeasured confounders

## Power
{
  
  # Power  -  02
  {
    setwd("/Users/jane.gao/Desktop/Hybrid experimental design/5. simulation")
    Res <- readRDS("res_02_08.RData")
    
    Res <- bind_rows(Res)
    Res$n <- seq(6, 800, by = 2) 
    
    
    Res_plot <-  Res %>% 
      pivot_longer(
        cols = "power_IF_full":"power_IF_simple", 
        names_to = "estimator",
        values_to = "power"
      )
    
    Res_plot$estimator <- factor(Res_plot$estimator, 
                                 levels = c("power_IF_simple",
                                            "power_IF_full",
                                            "power_IF_RCT",
                                            "power_naive"
                                 ))
    levels(Res_plot$estimator) <- c("Estimator 1", "Estimator 2", "Estimator 3", "Estimator 4")
    
    
    p1 <- ggplot(data = Res_plot, aes(x = n, y = power, colour = estimator)) +
      geom_line(linewidth = 0.65) +
      scale_colour_manual(values=c("#207CF5FF", "#038C33FF", "#FDC718FF", "#F27127FF"),
                          breaks=c("Estimator 1", "Estimator 2", "Estimator 3", "Estimator 4"))+
      geom_hline(yintercept = 0.8,
                 linetype="dashed",
                 color = "#084594",
                 size = 0.8) +
      # geom_vline(xintercept= 22,
      #            linetype="dashed",
      #            color = "#011C40FF",
      #            size = 0.8)+
      geom_vline(xintercept= 62,
                 linetype="dashed",
                 color =  "#005A32",
                 size = 0.8)+
      # geom_vline(xintercept= 30,
      #            linetype="dashed",
      #            color = "#05C7F2FF",
      #            size = 0.8)+
      geom_vline(xintercept= 384,
                 linetype="dashed",
                 color = "#99000D",
                 size = 0.8)+
      theme(legend.position = "none")+
      labs(x="n", y="Power")+
      annotate("text",x=600, y=0.7, label="beta = 0.8",size = 5)+
      labs(title= expression(paste("Power When ", tau, " = 0.8, ", pi, " = 0.2")))
    p1
    
    
    ggsave(plot = p1, 
           width = 5.5, 
           height = 5, 
           dpi = 500, 
           filename = "/Users/jane.gao/Desktop/Hybrid experimental design/5. simulation/plots/power_02.png")
    
    
  }
  
  
  # Power  -  04
  {
    setwd("/Users/jane.gao/Desktop/Hybrid experimental design/5. simulation")
    Res <- readRDS("res_04_08.RData")
    
    Res <- bind_rows(Res)
    Res$n <- seq(6, 800, by = 2) 
    
    
    Res_plot <-  Res %>% 
      pivot_longer(
        cols = "power_IF_full":"power_IF_simple", 
        names_to = "estimator",
        values_to = "power"
      )
    
    Res_plot$estimator <- factor(Res_plot$estimator, 
                                 levels = c("power_IF_simple",
                                            "power_IF_full",
                                            "power_IF_RCT",
                                            "power_naive"
                                 ))
    levels(Res_plot$estimator) <- c("Estimator 1", "Estimator 2", "Estimator 3", "Estimator 4")
    
    
    p1 <- ggplot(data = Res_plot, aes(x = n, y = power, colour = estimator)) +
      geom_line(linewidth = 0.65) +
      scale_colour_manual(values=c("#207CF5FF", "#038C33FF", "#FDC718FF", "#F27127FF"),
                          breaks=c("Estimator 1", "Estimator 2", "Estimator 3", "Estimator 4"))+
      geom_hline(yintercept = 0.8,
                 linetype="dashed",
                 color = "#084594",
                 size = 0.8) +
      # geom_vline(xintercept= 22,
      #            linetype="dashed",
      #            color = "#011C40FF",
      #            size = 0.8)+
      geom_vline(xintercept= 31,
                 linetype="dashed",
                 color =  "#005A32",
                 size = 0.8)+
      # geom_vline(xintercept= 30,
      #            linetype="dashed",
      #            color = "#05C7F2FF",
      #            size = 0.8)+
      geom_vline(xintercept= 257,
                 linetype="dashed",
                 color = "#99000D",
                 size = 0.8)+
      theme(legend.position = "none")+
      labs(x="n", y="Power")+
      annotate("text",x=600, y=0.7, label="beta = 0.8",size = 5)+
      labs(title= expression(paste("Power When ", tau, " = 0.8, ", pi, " = 0.4")))
    p1
    
    
    ggsave(plot = p1, 
           width = 5.5, 
           height = 5, 
           dpi = 500, 
           filename = "/Users/jane.gao/Desktop/Hybrid experimental design/5. simulation/plots/power_04.png")
    
    
  }
  
  
  # Power  -  05
  {
    setwd("/Users/jane.gao/Desktop/Hybrid experimental design/5. simulation")
    Res <- readRDS("res_05_08.RData")
    
    Res <- bind_rows(Res)
    Res$n <- seq(6, 800, by = 2) 
    
    
    Res_plot <-  Res %>% 
      pivot_longer(
        cols = "power_IF_full":"power_IF_simple", 
        names_to = "estimator",
        values_to = "power"
      )
    
    Res_plot$estimator <- factor(Res_plot$estimator, 
                                 levels = c("power_IF_simple",
                                            "power_IF_full",
                                            "power_IF_RCT",
                                            "power_naive"
                                 ))
    levels(Res_plot$estimator) <- c("Estimator 1", "Estimator 2", "Estimator 3", "Estimator 4")
    
    
    # check color
    {
      # display.brewer.all()
      
      # display.brewer.pal(n = 8, 'Reds')
      # brewer.pal(n = 8, 'Reds')
      # # "#FB6A4A"
      # 
      # display.brewer.pal(n = 8, 'Greens')
      # brewer.pal(n = 8, 'Greens')
      # # "#41AB5D"
      # 
      # display.brewer.pal(n = 8, 'Blues')
      # brewer.pal(n = 8, 'Blues')
      # # "#4292C6"
      
      # display.brewer.pal(n = 8, 'YlOrRd')
      # brewer.pal(n = 8, 'YlOrRd')
      # # "#FEB24C"
      
      
      # paletteer_d("fishualize::Acanthurus_leucosternon")
      # paletteer_d("fishualize::Acanthostracion_polygonius")
      # paletteer_d("fishualize::Acanthurus_olivaceus")
      # paletteer_d("fishualize::Acanthostracion_polygonius_y")
    }
    
    
    p1 <- ggplot(data = Res_plot, aes(x = n, y = power, colour = estimator)) +
      geom_line(linewidth = 0.65) +
      scale_colour_manual(values=c("#207CF5FF", "#038C33FF", "#FDC718FF", "#F27127FF"),
                          breaks=c("Estimator 1", "Estimator 2", "Estimator 3", "Estimator 4"))+
      geom_hline(yintercept = 0.8,
                 linetype="dashed",
                 color = "#084594",
                 size = 0.8) +
      # geom_vline(xintercept= 22,
      #            linetype="dashed",
      #            color = "#011C40FF",
      #            size = 0.8)+
      geom_vline(xintercept= 25,
                 linetype="dashed",
                 color =  "#005A32",
                 size = 0.8)+
      # geom_vline(xintercept= 30,
      #            linetype="dashed",
      #            color = "#05C7F2FF",
      #            size = 0.8)+
      geom_vline(xintercept= 246,
                 linetype="dashed",
                 color = "#99000D",
                 size = 0.8)+
      theme(legend.position = "none")+
      labs(x="n", y="Power")+
      annotate("text",x=600, y=0.7, label="beta = 0.8",size = 5)+
      labs(title= expression(paste("Power When ", tau, " = 0.8, ", pi, " = 0.5")))
    p1
    
    
    ggsave(plot = p1, 
           width = 5.5, 
           height = 5, 
           dpi = 500, 
           filename = "/Users/jane.gao/Desktop/Hybrid experimental design/5. simulation/plots/power_05.png")
    
    
  }
  
  
  # Power  -  06
  {
    setwd("/Users/jane.gao/Desktop/Hybrid experimental design/5. simulation")
    Res <- readRDS("res_06_08.RData")
    
    Res <- bind_rows(Res)
    Res$n <- seq(6, 800, by = 2) 
    
    
    Res_plot <-  Res %>% 
      pivot_longer(
        cols = "power_IF_full":"power_IF_simple", 
        names_to = "estimator",
        values_to = "power"
      )
    
    Res_plot$estimator <- factor(Res_plot$estimator, 
                                 levels = c("power_IF_simple",
                                            "power_IF_full",
                                            "power_IF_RCT",
                                            "power_naive"
                                 ))
    levels(Res_plot$estimator) <- c("Estimator 1", "Estimator 2", "Estimator 3", "Estimator 4")
    
    p1 <- ggplot(data = Res_plot, aes(x = n, y = power, colour = estimator)) +
      geom_line(linewidth = 0.65) +
      scale_colour_manual(values=c("#207CF5FF", "#038C33FF", "#FDC718FF", "#F27127FF"),
                          breaks=c("Estimator 1", "Estimator 2", "Estimator 3", "Estimator 4"))+
      geom_hline(yintercept = 0.8,
                 linetype="dashed",
                 color = "#084594",
                 size = 0.8) +
      # geom_vline(xintercept= 22,
      #            linetype="dashed",
      #            color = "#011C40FF",
      #            size = 0.8)+
      geom_vline(xintercept= 21,
                 linetype="dashed",
                 color =  "#005A32",
                 size = 0.8)+
      # geom_vline(xintercept= 30,
      #            linetype="dashed",
      #            color = "#05C7F2FF",
      #            size = 0.8)+
      geom_vline(xintercept= 257,
                 linetype="dashed",
                 color = "#99000D",
                 size = 0.8)+
      theme(legend.position = "none")+
      labs(x="n", y="Power")+
      annotate("text",x=600, y=0.7, label="beta = 0.8",size = 5)+
      labs(title= expression(paste("Power When ", tau, " = 0.8, ", pi, " = 0.6")))
    p1
    
    
    ggsave(plot = p1, 
           width = 5.5, 
           height = 5, 
           dpi = 500, 
           filename = "/Users/jane.gao/Desktop/Hybrid experimental design/5. simulation/plots/power_06.png")
    
    
  }
  
  
  # Power  -  08
  {
    setwd("/Users/jane.gao/Desktop/Hybrid experimental design/5. simulation")
    Res <- readRDS("res_08_08.RData")
    
    Res <- bind_rows(Res)
    Res$n <- seq(6, 800, by = 2) 
    
    
    Res_plot <-  Res %>% 
      pivot_longer(
        cols = "power_IF_full":"power_IF_simple", 
        names_to = "estimator",
        values_to = "power"
      )
    
    Res_plot$estimator <- factor(Res_plot$estimator, 
                                 levels = c("power_IF_simple",
                                            "power_IF_full",
                                            "power_IF_RCT",
                                            "power_naive"
                                 ))
    levels(Res_plot$estimator) <- c("Estimator 1", "Estimator 2", "Estimator 3", "Estimator 4")
    
    p1 <- ggplot(data = Res_plot, aes(x = n, y = power, colour = estimator)) +
      geom_line(linewidth = 0.65) +
      scale_colour_manual(values=c("#207CF5FF", "#038C33FF", "#FDC718FF", "#F27127FF"),
                          breaks=c("Estimator 1", "Estimator 2", "Estimator 3", "Estimator 4"))+
      geom_hline(yintercept = 0.8,
                 linetype="dashed",
                 color = "#084594",
                 size = 0.8) +
      # geom_vline(xintercept= 22,
      #            linetype="dashed",
      #            color = "#011C40FF",
      #            size = 0.8)+
      geom_vline(xintercept= 16,
                 linetype="dashed",
                 color =  "#005A32",
                 size = 0.8)+
      # geom_vline(xintercept= 30,
      #            linetype="dashed",
      #            color = "#05C7F2FF",
      #            size = 0.8)+
      geom_vline(xintercept= 384,
                 linetype="dashed",
                 color = "#99000D",
                 size = 0.8)+
      theme(legend.position = "none")+
      labs(x="n", y="Power")+
      annotate("text",x=600, y=0.7, label="beta = 0.8",size = 5)+
      labs(title= expression(paste("Power When ", tau, " = 0.8, ", pi, " = 0.8")))
    p1
    
    
    ggsave(plot = p1, 
           width = 5.5, 
           height = 5, 
           dpi = 500, 
           filename = "/Users/jane.gao/Desktop/Hybrid experimental design/5. simulation/plots/power_08.png")
    
    
  }
  
  
}


## Type I error
{
  
  # Type I error - 02
  {
    setwd("/Users/jane.gao/Desktop/Hybrid experimental design/5. simulation")
    Res <- readRDS("res_02_00.RData")
    
    Res <- bind_rows(Res)
    Res$n <- seq(6, 800, by = 2) 
    
    Res_plot <-  Res %>% 
      pivot_longer(
        cols = "power_IF_full":"power_IF_simple", 
        names_to = "estimator",
        values_to = "Type_I_Error"
      )
    
    Res_plot$estimator <- factor(Res_plot$estimator, 
                                 levels = c("power_IF_simple",
                                            "power_IF_full",
                                            "power_IF_RCT",
                                            "power_naive"
                                 ))
    levels(Res_plot$estimator) <- c("Estimator 1", "Estimator 2", "Estimator 3", "Estimator 4")
    
    
    p2 <- ggplot(data = Res_plot, aes(x = n, y = Type_I_Error, colour = estimator)) +
      geom_line(linewidth = 0.65) +
      scale_colour_manual(values=c("#207CF5FF", "#038C33FF", "#FDC718FF", "#F27127FF"),
                          breaks=c("Estimator 1", "Estimator 2", "Estimator 3", "Estimator 4"))+
      geom_hline(yintercept = 0.05,
                 linetype="dashed",
                 color = "#084594",
                 size = 0.8) +
      geom_vline(xintercept= 62,
                 linetype="dashed",
                 color =  "#005A32",
                 size = 0.8)+
      geom_vline(xintercept= 384,
                 linetype="dashed",
                 color = "#99000D",
                 size = 0.8)+
      labs(x="n", y="Type I Error")+
      annotate("text",x=600, y=0.1,label="alpha = 0.05",size=5)+
      labs(title= expression(paste("Type I Error When ", tau, " = 0, ", pi, " = 0.2")))+
      expand_limits(y=c(0, 0.546))+
      scale_y_continuous(breaks = seq(0, 0.546, by=0.1))
    p2
    
    
    ggsave(plot = p2, 
           width = 7, 
           height = 5, 
           dpi = 500, 
           filename = "/Users/jane.gao/Desktop/Hybrid experimental design/5. simulation/plots/TypeI_02.png")
    
    
  }
  
  # Type I error - 04
  {
    
    setwd("/Users/jane.gao/Desktop/Hybrid experimental design/5. simulation")
    Res <- readRDS("res_04_00.RData")
    
    Res <- bind_rows(Res)
    Res$n <- seq(6, 800, by = 2) 
    
    Res_plot <-  Res %>% 
      pivot_longer(
        cols = "power_IF_full":"power_IF_simple", 
        names_to = "estimator",
        values_to = "Type_I_Error"
      )
    
    Res_plot$estimator <- factor(Res_plot$estimator, 
                                 levels = c("power_IF_simple",
                                            "power_IF_full",
                                            "power_IF_RCT",
                                            "power_naive"
                                 ))
    levels(Res_plot$estimator) <- c("Estimator 1", "Estimator 2", "Estimator 3", "Estimator 4")
    
    
    p2 <- ggplot(data = Res_plot, aes(x = n, y = Type_I_Error, colour = estimator)) +
      geom_line(linewidth = 0.65) +
      scale_colour_manual(values=c("#207CF5FF", "#038C33FF", "#FDC718FF", "#F27127FF"),
                          breaks=c("Estimator 1", "Estimator 2", "Estimator 3", "Estimator 4"))+
      geom_hline(yintercept = 0.05,
                 linetype="dashed",
                 color = "#084594",
                 size = 0.8) +
      geom_vline(xintercept= 31,
                 linetype="dashed",
                 color =  "#005A32",
                 size = 0.8)+
      geom_vline(xintercept= 257,
                 linetype="dashed",
                 color = "#99000D",
                 size = 0.8)+
      labs(x="n", y="Type I Error")+
      annotate("text",x=600, y=0.1,label="alpha = 0.05",size=5)+
      labs(title= expression(paste("Type I Error When ", tau, " = 0, ", pi, " = 0.4")))+
      expand_limits(y=c(0, 0.546))+
      scale_y_continuous(breaks = seq(0, 0.546, by=0.1))
    p2
    
    
    ggsave(plot = p2, 
           width = 7, 
           height = 5, 
           dpi = 500, 
           filename = "/Users/jane.gao/Desktop/Hybrid experimental design/5. simulation/plots/TypeI_04.png")
    
    
  }
  
  # Type I error - 05
  {
    setwd("/Users/jane.gao/Desktop/Hybrid experimental design/5. simulation")
    Res <- readRDS("res_05_00.RData")
    
    Res <- bind_rows(Res)
    Res$n <- seq(6, 800, by = 2) 
    
    Res_plot <-  Res %>% 
    pivot_longer(
      cols = "power_IF_full":"power_IF_simple", 
      names_to = "estimator",
      values_to = "Type_I_Error"
    )
    
    Res_plot$estimator <- factor(Res_plot$estimator, 
                                 levels = c("power_IF_simple",
                                            "power_IF_full",
                                            "power_IF_RCT",
                                            "power_naive"
                                 ))
    levels(Res_plot$estimator) <- c("Estimator 1", "Estimator 2", "Estimator 3", "Estimator 4")
    
    
    p2 <- ggplot(data = Res_plot, aes(x = n, y = Type_I_Error, colour = estimator)) +
      geom_line(linewidth = 0.65) +
      scale_colour_manual(values=c("#207CF5FF", "#038C33FF", "#FDC718FF", "#F27127FF"),
                          breaks=c("Estimator 1", "Estimator 2", "Estimator 3", "Estimator 4"))+
      geom_hline(yintercept = 0.05,
                 linetype="dashed",
                 color = "#084594",
                 size = 0.8) +
      geom_vline(xintercept= 25,
                 linetype="dashed",
                 color =  "#005A32",
                 size = 0.8)+
      geom_vline(xintercept= 246,
                 linetype="dashed",
                 color = "#99000D",
                 size = 0.8)+
      labs(x="n", y="Type I Error")+
      annotate("text",x=600, y=0.1,label="alpha = 0.05",size=5)+
      labs(title= expression(paste("Type I Error When ", tau, " = 0, ", pi, " = 0.5")))+
      expand_limits(y=c(0, 0.546))+
      scale_y_continuous(breaks = seq(0, 0.546, by=0.1))
    p2
    
    
    ggsave(plot = p2, 
           width = 7, 
           height = 5, 
           dpi = 500, 
           filename = "/Users/jane.gao/Desktop/Hybrid experimental design/5. simulation/plots/TypeI_05.png")
    
    
  }
  
  # Type I error - 06
  {
    setwd("/Users/jane.gao/Desktop/Hybrid experimental design/5. simulation")
    Res <- readRDS("res_06_00.RData")
    
    Res <- bind_rows(Res)
    Res$n <- seq(6, 800, by = 2) 
    
    Res_plot <-  Res %>% 
      pivot_longer(
        cols = "power_IF_full":"power_IF_simple", 
        names_to = "estimator",
        values_to = "Type_I_Error"
      )
    
    Res_plot$estimator <- factor(Res_plot$estimator, 
                                 levels = c("power_IF_simple",
                                            "power_IF_full",
                                            "power_IF_RCT",
                                            "power_naive"
                                 ))
    levels(Res_plot$estimator) <- c("Estimator 1", "Estimator 2", "Estimator 3", "Estimator 4")
    
    
    p2 <- ggplot(data = Res_plot, aes(x = n, y = Type_I_Error, colour = estimator)) +
      geom_line(linewidth = 0.65) +
      scale_colour_manual(values=c("#207CF5FF", "#038C33FF", "#FDC718FF", "#F27127FF"),
                          breaks=c("Estimator 1", "Estimator 2", "Estimator 3", "Estimator 4"))+
      geom_hline(yintercept = 0.05,
                 linetype="dashed",
                 color = "#084594",
                 size = 0.8) +
      geom_vline(xintercept= 21,
                 linetype="dashed",
                 color =  "#005A32",
                 size = 0.8)+
      geom_vline(xintercept= 257,
                 linetype="dashed",
                 color = "#99000D",
                 size = 0.8)+
      labs(x="n", y="Type I Error")+
      annotate("text",x=600, y=0.1,label="alpha = 0.05",size=5)+
      labs(title= expression(paste("Type I Error When ", tau, " = 0, ", pi, " = 0.6")))+
      expand_limits(y=c(0, 0.546))+
      scale_y_continuous(breaks = seq(0, 0.546, by=0.1))
    p2
    
    
    ggsave(plot = p2, 
           width = 7, 
           height = 5, 
           dpi = 500, 
           filename = "/Users/jane.gao/Desktop/Hybrid experimental design/5. simulation/plots/TypeI_06.png")
    
    
  }
  
  # Type I error - 08
  {
    setwd("/Users/jane.gao/Desktop/Hybrid experimental design/5. simulation")
    Res <- readRDS("res_08_00.RData")
    
    Res <- bind_rows(Res)
    Res$n <- seq(6, 800, by = 2) 
    
    Res_plot <-  Res %>% 
      pivot_longer(
        cols = "power_IF_full":"power_IF_simple", 
        names_to = "estimator",
        values_to = "Type_I_Error"
      )
    
    Res_plot$estimator <- factor(Res_plot$estimator, 
                                 levels = c("power_IF_simple",
                                            "power_IF_full",
                                            "power_IF_RCT",
                                            "power_naive"
                                 ))
    levels(Res_plot$estimator) <- c("Estimator 1", "Estimator 2", "Estimator 3", "Estimator 4")
    
    
    p2 <- ggplot(data = Res_plot, aes(x = n, y = Type_I_Error, colour = estimator)) +
      geom_line(linewidth = 0.65) +
      scale_colour_manual(values=c("#207CF5FF", "#038C33FF", "#FDC718FF", "#F27127FF"),
                          breaks=c("Estimator 1", "Estimator 2", "Estimator 3", "Estimator 4"))+
      geom_hline(yintercept = 0.05,
                 linetype="dashed",
                 color = "#084594",
                 size = 0.8) +
      geom_vline(xintercept= 16,
                 linetype="dashed",
                 color =  "#005A32",
                 size = 0.8)+
      geom_vline(xintercept= 384,
                 linetype="dashed",
                 color = "#99000D",
                 size = 0.8)+
      labs(x="n", y="Type I Error")+
      annotate("text",x=600, y=0.1,label="alpha = 0.05",size=5)+
      labs(title= expression(paste("Type I Error When ", tau, " = 0, ", pi, " = 0.8")))+
      expand_limits(y=c(0, 0.546))+
      scale_y_continuous(breaks = seq(0, 0.546, by=0.1))
    p2
    
    
    ggsave(plot = p2, 
           width = 7, 
           height = 5, 
           dpi = 500, 
           filename = "/Users/jane.gao/Desktop/Hybrid experimental design/5. simulation/plots/TypeI_08.png")
    
    
  }
  
  
  
}


