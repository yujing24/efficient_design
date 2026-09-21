rm(list=ls())
library(tidyr)
library(tidyverse)
library(paletteer)
library(RColorBrewer)
library(ggplot2)
library(ggpubr)

setwd("/code/simulation/results")


### Main paper - Sufficient EC (N_E = 1000)
{
  ## Power sufficient
  {
    
    # Power  -  05
    {
      
      Res <- readRDS("res_suff/res_05_power_suff.RData")
      
      Res <- bind_rows(Res)
      Res$n <- seq(6, 500, by = 1) 
      
      Res_plot <-  Res %>% 
        pivot_longer(
          cols = c("power_EC", "power_AIPW", "power_naive"), 
          names_to = "Estimator",
          values_to = "power"
        )
      
      # add single-arm case
      Res_sa <- readRDS("res_suff/res_sa_power_suff.RData")
      
      Res_sa <- bind_rows(Res_sa)
      Res_sa$n <- seq(6, 500, by = 1) 
      
      Res_plot_sa <-  Res_sa %>% 
        pivot_longer(
          cols = c("power_SA"), 
          names_to = "Estimator",
          values_to = "power"
        )
      
      Res_plot <- rbind(Res_plot, Res_plot_sa)
      
      Res_plot$Estimator <- factor(Res_plot$Estimator, 
                                   levels = c(
                                     "power_EC",
                                     "power_AIPW",
                                     "power_naive",
                                     "power_SA"
                                   ))
      levels(Res_plot$Estimator) <- c(
        "Hybrid", 
        "AIPW",
        "Diff-in-mean",
        "Single-arm")
      
      p_suff_power_05 <- ggplot(data = Res_plot, aes(x = n, y = power, colour = Estimator)) +
        geom_line(linewidth = 1) +
        scale_colour_manual(values=c("#F27127FF","#FDC718FF",
                                     "#038C33FF", "#207CF5FF"),
                            breaks=c("Diff-in-mean", "AIPW",
                                     "Hybrid", "Single-arm"))+
        geom_hline(yintercept = 0.8,
                   linetype="dashed",
                   color = "#084594",
                   linewidth = 1) +
        geom_vline(xintercept= 42,
                   linetype="dashed",
                   color = "#011C40FF",
                   linewidth = 1)+
        geom_vline(xintercept= 83,
                   linetype="dashed",
                   color =  "#005A32",
                   linewidth = 1) +
        geom_vline(xintercept= 157,
                   linetype="dashed",
                   color = "#FEB24C",
                   linewidth = 1)+
        theme(legend.position = "none")+
        # labs(x=expression(RCT~Sample~Size~N[R]), y="Power")+
        labs(x=expression(Sample~Size~N[R]), y="Power")+
        scale_x_continuous(breaks=seq(0, 500, by=100))+
        scale_y_continuous(breaks=seq(0, 1, by = 0.2))+
        theme(title = element_text(size = 18),
              axis.title=element_text(size = 16),
              legend.title=element_text(size = 18),
              axis.text=element_text(size = 18),
              legend.text=element_text(size = 18))+
        labs(title = expression(Sufficient~EC~(N[E]==1000)~","~pi[A] ==0.5))
      
      p_suff_power_05
      
    }
    
    # Power  -  06
    {
      
      Res <- readRDS("res_suff/res_06_power_suff.RData")
      
      Res <- bind_rows(Res) 
      Res$n <- seq(6, 500, by = 1) 
      
      Res_plot <-  Res %>% 
        pivot_longer(
          cols = c("power_EC", "power_AIPW", "power_naive"), 
          names_to = "Estimator",
          values_to = "power"
        )
      
      # add single-arm case
      Res_sa <- readRDS("res_suff/res_sa_power_suff.RData")
      
      Res_sa <- bind_rows(Res_sa)
      Res_sa$n <- seq(6, 500, by = 1) 
      
      Res_plot_sa <-  Res_sa %>% 
        pivot_longer(
          cols = c("power_SA"), 
          names_to = "Estimator",
          values_to = "power"
        )
      
      Res_plot <- rbind(Res_plot, Res_plot_sa)
      
      Res_plot$Estimator <- factor(Res_plot$Estimator, 
                                   levels = c(
                                     "power_EC",
                                     "power_AIPW",
                                     "power_naive",
                                     "power_SA"
                                   ))
      levels(Res_plot$Estimator) <- c(
        "Hybrid", 
        "AIPW",
        "Diff-in-mean",
        "Single-arm")
      
      p_suff_power_06 <- ggplot(data = Res_plot, aes(x = n, y = power, colour = Estimator)) +
        geom_line(linewidth = 1) +
        scale_colour_manual(values=c("#F27127FF","#FDC718FF",
                                     "#038C33FF", "#207CF5FF"),
                            breaks=c("Diff-in-mean", "AIPW",
                                     "Hybrid", "Single-arm"))+
        geom_hline(yintercept = 0.8,
                   linetype="dashed",
                   color = "#084594",
                   linewidth = 1) +
        geom_vline(xintercept= 42,
                   linetype="dashed",
                   color = "#011C40FF",
                   linewidth = 1)+
        geom_vline(xintercept= 69,
                   linetype="dashed",
                   color =  "#005A32",
                   linewidth = 1)+
        geom_vline(xintercept= 164,
                   linetype="dashed",
                   color = "#FEB24C",
                   linewidth = 1)+
        theme(legend.position = "none")+
        labs(x=expression(Sample~Size~N[R]), y="Power")+
        scale_x_continuous(breaks=seq(0, 500, by=100))+
        scale_y_continuous(breaks=seq(0, 1, by = 0.2))+
        theme(title = element_text(size = 18),
              axis.title=element_text(size = 16),
              legend.title=element_text(size = 18),
              axis.text=element_text(size = 18),
              legend.text=element_text(size = 18))+
        labs(title = expression(Sufficient~EC~(N[E]==1000)~","~pi[A] ==0.6))
      
      p_suff_power_06
      
    }
    
    # Power  -  07
    {
      
      Res <- readRDS("res_suff/res_07_power_suff.RData")
      
      Res <- bind_rows(Res)
      Res$n <- seq(6, 500, by = 1) 
      
      Res_plot <-  Res %>% 
        pivot_longer(
          cols = c("power_EC", "power_AIPW", "power_naive"), 
          names_to = "Estimator",
          values_to = "power"
        )
      
      # add single-arm case
      Res_sa <- readRDS("res_suff/res_sa_power_suff.RData")
      
      Res_sa <- bind_rows(Res_sa)
      Res_sa$n <- seq(6, 500, by = 1) 
      
      Res_plot_sa <-  Res_sa %>% 
        pivot_longer(
          cols = c("power_SA"), 
          names_to = "Estimator",
          values_to = "power"
        )
      
      Res_plot <- rbind(Res_plot, Res_plot_sa)
      
      Res_plot$Estimator <- factor(Res_plot$Estimator, 
                                   levels = c(
                                     "power_EC",
                                     "power_AIPW",
                                     "power_naive",
                                     "power_SA"
                                   ))
      levels(Res_plot$Estimator) <- c(
        "Hybrid", 
        "AIPW",
        "Diff-in-mean",
        "Single-arm")
      
      p_suff_power_07 <- ggplot(data = Res_plot, aes(x = n, y = power, colour = Estimator)) +
        geom_line(linewidth = 1) +
        scale_colour_manual(values=c("#F27127FF","#FDC718FF",
                                     "#038C33FF", "#207CF5FF"),
                            breaks=c("Diff-in-mean", "AIPW",
                                     "Hybrid", "Single-arm"))+
        geom_hline(yintercept = 0.8,
                   linetype="dashed",
                   color = "#084594",
                   linewidth = 1) +
        geom_vline(xintercept= 42,
                   linetype="dashed",
                   color = "#011C40FF",
                   linewidth = 1)+
        geom_vline(xintercept= 59,
                   linetype="dashed",
                   color =  "#005A32",
                   linewidth = 1)+
        geom_vline(xintercept= 187,
                   linetype="dashed",
                   color = "#FEB24C",
                   linewidth = 1)+
        theme(legend.position = "none")+
        labs(x=expression(Sample~Size~N[R]), y="Power")+
        scale_x_continuous(breaks=seq(0, 500, by=100))+
        scale_y_continuous(breaks=seq(0, 1, by = 0.2))+
        theme(title = element_text(size = 18),
              axis.title=element_text(size = 16),
              legend.title=element_text(size = 18),
              axis.text=element_text(size = 18),
              legend.text=element_text(size = 18))+
        labs(title = expression(Sufficient~EC~(N[E]==1000)~","~pi[A] ==0.7))
      
      p_suff_power_07
      
    }
    
    # Power  -  08
    {
      
      Res <- readRDS("res_suff/res_08_power_suff.RData")
      
      Res <- bind_rows(Res) 
      Res$n <- seq(6, 500, by = 1) 
      
      Res_plot <-  Res %>% 
        pivot_longer(
          cols = c("power_EC", "power_AIPW", "power_naive"), 
          names_to = "Estimator",
          values_to = "power"
        )
      
      # add single-arm case
      Res_sa <- readRDS("res_suff/res_sa_power_suff.RData")
      
      Res_sa <- bind_rows(Res_sa)
      Res_sa$n <- seq(6, 500, by = 1) 
      
      Res_plot_sa <-  Res_sa %>% 
        pivot_longer(
          cols = c("power_SA"), 
          names_to = "Estimator",
          values_to = "power"
        )
      
      Res_plot <- rbind(Res_plot, Res_plot_sa)
      
      Res_plot$Estimator <- factor(Res_plot$Estimator, 
                                   levels = c(
                                     "power_EC",
                                     "power_AIPW",
                                     "power_naive",
                                     "power_SA"
                                   ))
      levels(Res_plot$Estimator) <- c(
        "Hybrid", 
        "AIPW",
        "Diff-in-mean",
        "Single-arm")
    
      p_suff_power_08 <- ggplot(data = Res_plot, aes(x = n, y = power, colour = Estimator)) +
        geom_line(linewidth = 1) +
        scale_colour_manual(values=c("#F27127FF","#FDC718FF",
                                     "#038C33FF", "#207CF5FF"),
                            breaks=c("Diff-in-mean", "AIPW",
                                     "Hybrid", "Single-arm"))+
        geom_hline(yintercept = 0.8,
                   linetype="dashed",
                   color = "#084594",
                   linewidth = 1) +
        geom_vline(xintercept= 42,
                   linetype="dashed",
                   color = "#011C40FF",
                   linewidth = 1)+
        geom_vline(xintercept= 52,
                   linetype="dashed",
                   color =  "#005A32",
                   linewidth = 1)+
        geom_vline(xintercept= 246,
                   linetype="dashed",
                   color = "#FEB24C",
                   linewidth = 1)+
        theme(legend.position = "none")+
        labs(x=expression(Sample~Size~N[R]), y="Power")+
        scale_x_continuous(breaks=seq(0, 500, by=100))+
        scale_y_continuous(breaks=seq(0, 1, by = 0.2))+
        theme(title = element_text(size = 18),
              axis.title=element_text(size = 16),
              legend.title=element_text(size = 18),
              axis.text=element_text(size = 18),
              legend.text=element_text(size = 18))+
        labs(title = expression(Sufficient~EC~(N[E]==1000)~","~pi[A] ==0.8))
      
      p_suff_power_08
      
    }
    
    # Power  -  09
    {
      
      Res <- readRDS("res_suff/res_09_power_suff.RData")
      
      Res <- bind_rows(Res)
      # Res$n <- seq(6, 800, by = 2) 
      Res$n <- seq(6, 500, by = 1) 
      
      Res_plot <-  Res %>% 
        pivot_longer(
          cols = c("power_EC", "power_AIPW", "power_naive"), 
          names_to = "Estimator",
          values_to = "power"
        )
      
      # add single-arm case
      Res_sa <- readRDS("res_suff/res_sa_power_suff.RData")
      
      Res_sa <- bind_rows(Res_sa)
      Res_sa$n <- seq(6, 500, by = 1) 
      
      Res_plot_sa <-  Res_sa %>% 
        pivot_longer(
          cols = c("power_SA"), 
          names_to = "Estimator",
          values_to = "power"
        )
      
      Res_plot <- rbind(Res_plot, Res_plot_sa)
      
      Res_plot$Estimator <- factor(Res_plot$Estimator, 
                                   levels = c(
                                     "power_EC",
                                     "power_AIPW",
                                     "power_naive",
                                     "power_SA"
                                   ))
      levels(Res_plot$Estimator) <- c(
        "Hybrid", 
        "AIPW",
        "Diff-in-mean",
        "Single-arm")
      
      p_suff_power_09 <- ggplot(data = Res_plot, aes(x = n, y = power, colour = Estimator)) +
        geom_line(linewidth = 1) +
        scale_colour_manual(values=c("#F27127FF","#FDC718FF",
                                     "#038C33FF", "#207CF5FF"),
                            breaks=c("Diff-in-mean", "AIPW",
                                     "Hybrid", "Single-arm"))+
        geom_hline(yintercept = 0.8,
                   linetype="dashed",
                   color = "#084594",
                   linewidth = 1) +
        geom_vline(xintercept= 42,
                   linetype="dashed",
                   color = "#011C40FF",
                   linewidth = 1)+
        geom_vline(xintercept= 46,
                   linetype="dashed",
                   color =  "#005A32",
                   linewidth = 1)+
        geom_vline(xintercept= 437,
                   linetype="dashed",
                   color = "#FEB24C",
                   linewidth = 1)+
        # theme(legend.position = "none")+
        labs(x=expression(Sample~Size~N[R]), y="Power")+
        scale_x_continuous(breaks=seq(0, 500, by=100))+
        scale_y_continuous(breaks=seq(0, 1, by = 0.2))+
        theme(title = element_text(size = 18),
              axis.title=element_text(size = 16),
              legend.title=element_text(size = 18),
              axis.text=element_text(size = 18),
              legend.text=element_text(size = 18))+
        labs(title = expression(Sufficient~EC~(N[E]==1000)~","~pi[A] ==0.9))
      
      p_suff_power_09
      
    }
    
  }
  
  ## Type I error sufficient
  {
    
    # Type I error - 05
    {
      
      Res <- readRDS("res_suff/res_05_TypeI_suff.RData")
      
      Res <- bind_rows(Res)
      Res$n <- seq(6, 500, by = 1) 
      
      Res_plot <-  Res %>% 
        pivot_longer(
          cols = c("power_EC", "power_AIPW", "power_naive"), 
          names_to = "Estimator",
          values_to = "Type_I_Error"
        )
      
      # add single-arm case
      Res_sa <- readRDS("res_suff/res_sa_TypeI_suff.RData")
      
      Res_sa <- bind_rows(Res_sa)
      Res_sa$n <- seq(6, 500, by = 1) 
      
      Res_plot_sa <-  Res_sa %>% 
        pivot_longer(
          cols = c("power_SA"), 
          names_to = "Estimator",
          values_to = "Type_I_Error"
        )
      
      Res_plot <- rbind(Res_plot, Res_plot_sa)
      
      Res_plot$Estimator <- factor(Res_plot$Estimator, 
                                   levels = c(
                                     "power_EC",
                                     "power_AIPW",
                                     "power_naive",
                                     "power_SA"
                                   ))
      levels(Res_plot$Estimator) <- c(
        "Hybrid", 
        "AIPW",
        "Diff-in-mean",
        "Single-arm")
      
      p_suff_TypeI_05 <- ggplot(data = Res_plot, aes(x = n, y = Type_I_Error, colour = Estimator)) +
        geom_line(linewidth = 1) +
        scale_colour_manual(values=c("#F27127FF","#FDC718FF",
                                     "#038C33FF", "#207CF5FF"), 
                            breaks=c("Diff-in-mean", "AIPW",
                                     "Hybrid", "Single-arm"))+
        geom_hline(yintercept = 0.05,
                   linetype="dashed",
                   color = "#084594",
                   linewidth = 1) +
        geom_vline(xintercept= 42,
                   linetype="dashed",
                   color = "#011C40FF",
                   linewidth = 1)+
        geom_vline(xintercept= 83,
                   linetype="dashed",
                   color =  "#005A32",
                   linewidth = 1)+
        geom_vline(xintercept= 157,
                   linetype="dashed",
                   color = "#FEB24C",
                   linewidth = 1)+
        theme(legend.position = "none")+
        labs(x=expression(Sample~Size~N[R]), y="Type I Error")+
        scale_x_continuous(breaks=seq(0, 500, by=100))+
        scale_y_continuous(breaks = seq(0, 0.546, by=0.1))+
        expand_limits(y=c(0, 0.546))+
        theme(title = element_text(size = 18),
              axis.title=element_text(size = 16),
              legend.title=element_text(size = 18),
              axis.text=element_text(size = 18),
              legend.text=element_text(size = 18))+
        labs(title = expression(Sufficient~EC~(N[E]==1000)~","~pi[A] ==0.5))
      
      p_suff_TypeI_05

    }
    
    # Type I error - 06
    {
      
      Res <- readRDS("res_suff/res_06_TypeI_suff.RData")
      
      Res <- bind_rows(Res)
      Res$n <- seq(6, 500, by = 1) 
      
      Res_plot <-  Res %>% 
        pivot_longer(
          cols = c("power_EC", "power_AIPW", "power_naive"), 
          names_to = "Estimator",
          values_to = "Type_I_Error"
        )
      
      # add single-arm case
      Res_sa <- readRDS("res_suff/res_sa_TypeI_suff.RData")
      
      Res_sa <- bind_rows(Res_sa)
      Res_sa$n <- seq(6, 500, by = 1) 
      
      Res_plot_sa <-  Res_sa %>% 
        pivot_longer(
          cols = c("power_SA"), 
          names_to = "Estimator",
          values_to = "Type_I_Error"
        )
      
      Res_plot <- rbind(Res_plot, Res_plot_sa)
      
      Res_plot$Estimator <- factor(Res_plot$Estimator, 
                                   levels = c(
                                     "power_EC",
                                     "power_AIPW",
                                     "power_naive",
                                     "power_SA"
                                   ))
      levels(Res_plot$Estimator) <- c(
        "Hybrid", 
        "AIPW",
        "Diff-in-mean",
        "Single-arm")
      
      
      p_suff_TypeI_06 <- ggplot(data = Res_plot, aes(x = n, y = Type_I_Error, colour = Estimator)) +
        geom_line(linewidth = 1) +
        scale_colour_manual(values=c("#F27127FF","#FDC718FF",
                                     "#038C33FF", "#207CF5FF"), 
                            breaks=c("Diff-in-mean", "AIPW",
                                     "Hybrid", "Single-arm"))+
        geom_hline(yintercept = 0.05,
                   linetype="dashed",
                   color = "#084594",
                   linewidth = 1) +
        geom_vline(xintercept= 42,
                   linetype="dashed",
                   color = "#011C40FF",
                   linewidth = 1)+
        geom_vline(xintercept= 69,
                   linetype="dashed",
                   color =  "#005A32",
                   linewidth = 1)+
        geom_vline(xintercept= 164,
                   linetype="dashed",
                   color = "#FEB24C",
                   linewidth = 1)+
        theme(legend.position = "none")+
        labs(x=expression(Sample~Size~N[R]), y="Type I Error")+
        scale_x_continuous(breaks=seq(0, 500, by=100))+
        scale_y_continuous(breaks = seq(0, 0.546, by=0.1))+
        expand_limits(y=c(0, 0.546))+
        theme(title = element_text(size = 18),
              axis.title=element_text(size = 16),
              legend.title=element_text(size = 18),
              axis.text=element_text(size = 18),
              legend.text=element_text(size = 18))+
        labs(title = expression(Sufficient~EC~(N[E]==1000)~","~pi[A] ==0.6))
      
      p_suff_TypeI_06
      
    }
    
    # Type I error - 07
    {
      
      Res <- readRDS("res_suff/res_07_TypeI_suff.RData")
      
      Res <- bind_rows(Res)
      Res$n <- seq(6, 500, by = 1) 
      
      Res_plot <-  Res %>% 
        pivot_longer(
          cols = c("power_EC", "power_AIPW", "power_naive"), 
          names_to = "Estimator",
          values_to = "Type_I_Error"
        )
      
      # add single-arm case
      Res_sa <- readRDS("res_suff/res_sa_TypeI_suff.RData")
      
      Res_sa <- bind_rows(Res_sa)
      Res_sa$n <- seq(6, 500, by = 1) 
      
      Res_plot_sa <-  Res_sa %>% 
        pivot_longer(
          cols = c("power_SA"), 
          names_to = "Estimator",
          values_to = "Type_I_Error"
        )
      
      Res_plot <- rbind(Res_plot, Res_plot_sa)
      
      Res_plot$Estimator <- factor(Res_plot$Estimator, 
                                   levels = c(
                                     "power_EC",
                                     "power_AIPW",
                                     "power_naive",
                                     "power_SA"
                                   ))
      levels(Res_plot$Estimator) <- c(
        "Hybrid", 
        "AIPW",
        "Diff-in-mean",
        "Single-arm")
      
      
      p_suff_TypeI_07 <- ggplot(data = Res_plot, aes(x = n, y = Type_I_Error, colour = Estimator)) +
        geom_line(linewidth = 1) +
        scale_colour_manual(values=c("#F27127FF","#FDC718FF",
                                     "#038C33FF", "#207CF5FF"), 
                            breaks=c("Diff-in-mean", "AIPW",
                                     "Hybrid", "Single-arm"))+
        geom_hline(yintercept = 0.05,
                   linetype="dashed",
                   color = "#084594",
                   linewidth = 1) +
        geom_vline(xintercept= 42,
                   linetype="dashed",
                   color = "#011C40FF",
                   linewidth = 1)+
        geom_vline(xintercept= 59,
                   linetype="dashed",
                   color =  "#005A32",
                   linewidth = 1)+
        geom_vline(xintercept= 187,
                   linetype="dashed",
                   color = "#FEB24C",
                   linewidth = 1)+
        theme(legend.position = "none")+
        labs(x=expression(Sample~Size~N[R]), y="Type I Error")+
        scale_x_continuous(breaks=seq(0, 500, by=100))+
        scale_y_continuous(breaks = seq(0, 0.546, by=0.1))+
        expand_limits(y=c(0, 0.546))+
        theme(title = element_text(size = 18),
              axis.title=element_text(size = 16),
              legend.title=element_text(size = 18),
              axis.text=element_text(size = 18),
              legend.text=element_text(size = 18))+
        labs(title = expression(Sufficient~EC~(N[E]==1000)~","~pi[A] ==0.7))
      
      p_suff_TypeI_07
      
    }
    
    # Type I error - 08
    {
      
      Res <- readRDS("res_suff/res_08_TypeI_suff.RData")
      
      Res <- bind_rows(Res)
      Res$n <- seq(6, 500, by = 1) 
      
      Res_plot <-  Res %>% 
        pivot_longer(
          cols = c("power_EC", "power_AIPW", "power_naive"), 
          names_to = "Estimator",
          values_to = "Type_I_Error"
        )
      
      # add single-arm case
      Res_sa <- readRDS("res_suff/res_sa_TypeI_suff.RData")
      
      Res_sa <- bind_rows(Res_sa)
      Res_sa$n <- seq(6, 500, by = 1) 
      
      Res_plot_sa <-  Res_sa %>% 
        pivot_longer(
          cols = c("power_SA"), 
          names_to = "Estimator",
          values_to = "Type_I_Error"
        )
      
      Res_plot <- rbind(Res_plot, Res_plot_sa)
      
      Res_plot$Estimator <- factor(Res_plot$Estimator, 
                                   levels = c(
                                     "power_EC",
                                     "power_AIPW",
                                     "power_naive",
                                     "power_SA"
                                   ))
      levels(Res_plot$Estimator) <- c(
        "Hybrid", 
        "AIPW",
        "Diff-in-mean",
        "Single-arm")
      
      
      p_suff_TypeI_08 <- ggplot(data = Res_plot, aes(x = n, y = Type_I_Error, colour = Estimator)) +
        geom_line(linewidth = 1) +
        scale_colour_manual(values=c("#F27127FF","#FDC718FF",
                                     "#038C33FF", "#207CF5FF"), 
                            breaks=c("Diff-in-mean", "AIPW",
                                     "Hybrid", "Single-arm"))+
        geom_hline(yintercept = 0.05,
                   linetype="dashed",
                   color = "#084594",
                   linewidth = 1) +
        geom_vline(xintercept= 42,
                   linetype="dashed",
                   color = "#011C40FF",
                   linewidth = 1)+
        geom_vline(xintercept= 52,
                   linetype="dashed",
                   color =  "#005A32",
                   linewidth = 1)+
        geom_vline(xintercept= 246,
                   linetype="dashed",
                   color = "#FEB24C",
                   linewidth = 1)+
        # theme(legend.position = "none")+
        labs(x=expression(Sample~Size~N[R]), y="Type I Error")+
        scale_x_continuous(breaks=seq(0, 500, by=100))+
        scale_y_continuous(breaks = seq(0, 0.546, by=0.1))+
        expand_limits(y=c(0, 0.546))+
        theme(title = element_text(size = 18),
              axis.title=element_text(size = 16),
              legend.title=element_text(size = 18),
              axis.text=element_text(size = 18),
              legend.text=element_text(size = 18))+
        labs(title = expression(Sufficient~EC~(N[E]==1000)~","~pi[A] ==0.8))
      
      p_suff_TypeI_08
      
    }
    
    # Type I error - 09
    {
      
      Res <- readRDS("res_suff/res_09_TypeI_suff.RData")
      
      Res <- bind_rows(Res)
      Res$n <- seq(6, 500, by = 1) 
      
      Res_plot <-  Res %>% 
        pivot_longer(
          cols = c("power_EC", "power_AIPW", "power_naive"), 
          names_to = "Estimator",
          values_to = "Type_I_Error"
        )
      
      # add single-arm case
      Res_sa <- readRDS("res_suff/res_sa_TypeI_suff.RData")
      
      Res_sa <- bind_rows(Res_sa)
      Res_sa$n <- seq(6, 500, by = 1) 
      
      Res_plot_sa <-  Res_sa %>% 
        pivot_longer(
          cols = c("power_SA"), 
          names_to = "Estimator",
          values_to = "Type_I_Error"
        )
      
      Res_plot <- rbind(Res_plot, Res_plot_sa)
      
      Res_plot$Estimator <- factor(Res_plot$Estimator, 
                                   levels = c(
                                     "power_EC",
                                     "power_AIPW",
                                     "power_naive",
                                     "power_SA"
                                   ))
      levels(Res_plot$Estimator) <- c(
        "Hybrid", 
        "AIPW",
        "Diff-in-mean",
        "Single-arm")
      
      
      p_suff_TypeI_09 <- ggplot(data = Res_plot, aes(x = n, y = Type_I_Error, colour = Estimator)) +
        geom_line(linewidth = 1) +
        scale_colour_manual(values=c("#F27127FF","#FDC718FF",
                                     "#038C33FF", "#207CF5FF"), 
                            breaks=c("Diff-in-mean", "AIPW",
                                     "Hybrid", "Single-arm"))+
        geom_hline(yintercept = 0.05,
                   linetype="dashed",
                   color = "#084594",
                   linewidth = 1) +
        geom_vline(xintercept= 42,
                   linetype="dashed",
                   color = "#011C40FF",
                   linewidth = 1)+
        geom_vline(xintercept= 46,
                   linetype="dashed",
                   color =  "#005A32",
                   linewidth = 1)+
        geom_vline(xintercept= 437,
                   linetype="dashed",
                   color = "#FEB24C",
                   linewidth = 1)+
        # theme(legend.position = "none")+
        labs(x=expression(Sample~Size~N[R]), y="Type I Error")+
        scale_x_continuous(breaks=seq(0, 500, by=100))+
        scale_y_continuous(breaks = seq(0, 0.546, by=0.1))+
        expand_limits(y=c(0, 0.546))+
        theme(title = element_text(size = 18),
              axis.title=element_text(size = 16),
              legend.title=element_text(size = 18),
              axis.text=element_text(size = 18),
              legend.text=element_text(size = 18))+
        labs(title = expression(Sufficient~EC~(N[E]==1000)~","~pi[A] ==0.9))
      
      p_suff_TypeI_09
      
    }
    
  }
  
}

### Main paper - Insufficient EC (N_E = 60)
{
  ## Power insufficient
  {
    
    # Power  -  05
    {
      
      Res <- readRDS("res_insuff_60/res_05_power_insuff.RData")
      
      Res <- bind_rows(Res)
      Res$n <- seq(6, 500, by = 1) 
      
      Res_plot <-  Res %>% 
        pivot_longer(
          cols = c("power_EC", "power_AIPW", "power_naive"), 
          names_to = "Estimator",
          values_to = "power"
        )
      
      Res_plot$Estimator <- factor(Res_plot$Estimator, 
                                   levels = c("power_EC",
                                              "power_AIPW",
                                              "power_naive"
                                   ))
      levels(Res_plot$Estimator) <- c("Hybrid", 
                                      "AIPW",
                                      "Diff-in-mean")
      
      p_insuff_60_power_05 <- ggplot(data = Res_plot, aes(x = n, y = power, colour = Estimator)) +
        geom_line(linewidth = 1) +
        scale_colour_manual(values=c("#F27127FF",
                                     "#FDC718FF",
                                     "#038C33FF"),
                            breaks=c("Diff-in-mean",
                                     "AIPW",
                                     "Hybrid"))+
        geom_hline(yintercept = 0.8,
                   linetype="dashed",
                   color = "#084594",
                   linewidth = 1) +
        # geom_vline(xintercept= 119,
        #            linetype="dashed",
        #            color = "#011C40FF",
        #            linewidth = 1)+
        geom_vline(xintercept= 126,
                   linetype="dashed",
                   color =  "#005A32",
                   linewidth = 1) +
        # geom_vline(xintercept= 145,
        #            linetype="dashed",
        #            color = "#99000D",
        #            linewidth = 1)+
        geom_vline(xintercept= 157,
                   linetype="dashed",
                   color = "#FEB24C",
                   linewidth = 1)+
        theme(legend.position = "none")+
        labs(x=expression(Sample~Size~N[R]), y="Power")+
        scale_x_continuous(breaks=seq(0, 500, by=100))+
        scale_y_continuous(breaks=seq(0, 1, by = 0.2))+
        theme(title = element_text(size = 18),
              axis.title=element_text(size = 16),
              legend.title=element_text(size = 18),
              axis.text=element_text(size = 18),
              legend.text=element_text(size = 18))+
        labs(title = expression(Insufficient~EC~(N[E]==60)~","~pi[A] ==0.5))
      
      p_insuff_60_power_05
      
    }
    
    # Power  -  06
    {
      
      Res <- readRDS("res_insuff_60/res_06_power_insuff.RData")
      
      Res <- bind_rows(Res)
      Res$n <- seq(6, 500, by = 1) 
      
      Res_plot <-  Res %>% 
        pivot_longer(
          cols = c("power_EC", "power_AIPW", "power_naive"), 
          names_to = "Estimator",
          values_to = "power"
        )
      
      Res_plot$Estimator <- factor(Res_plot$Estimator, 
                                   levels = c("power_EC",
                                              "power_AIPW",
                                              "power_naive"
                                   ))
      levels(Res_plot$Estimator) <- c("Hybrid", 
                                      "AIPW",
                                      "Diff-in-mean")
      
      
      p_insuff_60_power_06 <- ggplot(data = Res_plot, aes(x = n, y = power, colour = Estimator)) +
        geom_line(linewidth = 1) +
        scale_colour_manual(values=c("#F27127FF",
                                     "#FDC718FF",
                                     "#038C33FF"),
                            breaks=c("Diff-in-mean",
                                     "AIPW",
                                     "Hybrid"))+
        geom_hline(yintercept = 0.8,
                   linetype="dashed",
                   color = "#084594",
                   linewidth = 1) +
        geom_vline(xintercept= 118,
                   linetype="dashed",
                   color =  "#005A32",
                   linewidth = 1)+
        # geom_vline(xintercept= 134,
        #            linetype="dashed",
        #            color = "#99000D",
        #            linewidth = 1)+
        geom_vline(xintercept= 164,
                   linetype="dashed",
                   color = "#FEB24C",
                   linewidth = 1)+
        theme(legend.position = "none")+
        labs(x=expression(Sample~Size~N[R]), y="Power")+
        scale_x_continuous(breaks=seq(0, 500, by=100))+
        scale_y_continuous(breaks=seq(0, 1, by = 0.2))+
        theme(title = element_text(size = 18),
              axis.title=element_text(size = 16),
              legend.title=element_text(size = 18),
              axis.text=element_text(size = 18),
              legend.text=element_text(size = 18))+
        labs(title = expression(Insufficient~EC~(N[E]==60)~","~pi[A] ==0.6))
      
      p_insuff_60_power_06
      
    }
    
    # Power  -  07
    {
      
      Res <- readRDS("res_insuff_60/res_07_power_insuff.RData")
      
      Res <- bind_rows(Res)
      Res$n <- seq(6, 500, by = 1) 
      
      Res_plot <-  Res %>% 
        pivot_longer(
          cols = c("power_EC", "power_AIPW", "power_naive"), 
          names_to = "Estimator",
          values_to = "power"
        )
      
      Res_plot$Estimator <- factor(Res_plot$Estimator, 
                                   levels = c("power_EC",
                                              "power_AIPW",
                                              "power_naive"
                                   ))
      levels(Res_plot$Estimator) <- c("Hybrid", 
                                      "AIPW",
                                      "Diff-in-mean")
      
      
      p_insuff_60_power_07 <- ggplot(data = Res_plot, aes(x = n, y = power, colour = Estimator)) +
        geom_line(linewidth = 1) +
        scale_colour_manual(values=c("#F27127FF",
                                     "#FDC718FF",
                                     "#038C33FF"),
                            breaks=c("Diff-in-mean",
                                     "AIPW",
                                     "Hybrid"))+
        geom_hline(yintercept = 0.8,
                   linetype="dashed",
                   color = "#084594",
                   linewidth = 1) +
        # geom_vline(xintercept= 108,
        #            linetype="dashed",
        #            color = "#011C40FF",
        #            linewidth = 1)+
        geom_vline(xintercept= 116,
                   linetype="dashed",
                   color =  "#005A32",
                   linewidth = 1)+
        # geom_vline(xintercept= 128,
        #            linetype="dashed",
        #            color = "#99000D",
        #            linewidth = 1)+
        geom_vline(xintercept= 187,
                   linetype="dashed",
                   color = "#FEB24C",
                   linewidth = 1)+
        theme(legend.position = "none")+
        labs(x=expression(Sample~Size~N[R]), y="Power")+
        scale_x_continuous(breaks=seq(0, 500, by=100))+
        scale_y_continuous(breaks=seq(0, 1, by = 0.2))+
        theme(title = element_text(size = 18),
              axis.title=element_text(size = 16),
              legend.title=element_text(size = 18),
              axis.text=element_text(size = 18),
              legend.text=element_text(size = 18))+
        labs(title = expression(Insufficient~EC~(N[E]==60)~","~pi[A] ==0.7))
      
      p_insuff_60_power_07
      
      
    }
    
    # Power  -  08
    {
      
      Res <- readRDS("res_insuff_60/res_08_power_insuff.RData")
      
      Res <- bind_rows(Res)
      Res$n <- seq(6, 500, by = 1) 
      
      Res_plot <-  Res %>% 
        pivot_longer(
          cols = c("power_EC", "power_AIPW", "power_naive"), 
          names_to = "Estimator",
          values_to = "power"
        )
      
      Res_plot$Estimator <- factor(Res_plot$Estimator, 
                                   levels = c("power_EC",
                                              "power_AIPW",
                                              "power_naive"
                                   ))
      levels(Res_plot$Estimator) <- c("Hybrid", 
                                      "AIPW",
                                      "Diff-in-mean")
      
      
      p_insuff_60_power_08 <- ggplot(data = Res_plot, aes(x = n, y = power, colour = Estimator)) +
        geom_line(linewidth = 1) +
        scale_colour_manual(values=c("#F27127FF",
                                     "#FDC718FF",
                                     "#038C33FF"),
                            breaks=c("Diff-in-mean",
                                     "AIPW",
                                     "Hybrid"))+
        geom_hline(yintercept = 0.8,
                   linetype="dashed",
                   color = "#084594",
                   linewidth = 1) +
        # geom_vline(xintercept= 116,
        #            linetype="dashed",
        #            color = "#011C40FF",
        #            linewidth = 1)+
        geom_vline(xintercept= 124,
                   linetype="dashed",
                   color =  "#005A32",
                   linewidth = 1)+
        # geom_vline(xintercept= 131,
        #            linetype="dashed",
        #            color = "#99000D",
        #            linewidth = 1)+
        geom_vline(xintercept= 246,
                   linetype="dashed",
                   color = "#FEB24C",
                   linewidth = 1)+
        theme(legend.position = "none")+
        labs(x=expression(Sample~Size~N[R]), y="Power")+
        scale_x_continuous(breaks=seq(0, 500, by=100))+
        scale_y_continuous(breaks=seq(0, 1, by = 0.2))+
        theme(title = element_text(size = 18),
              axis.title=element_text(size = 16),
              legend.title=element_text(size = 18),
              axis.text=element_text(size = 18),
              legend.text=element_text(size = 18))+ 
        labs(title = expression(Insufficient~EC~(N[E]==60)~","~pi[A] ==0.8))
      
      p_insuff_60_power_08
      
      
    }
    
    # Power  -  09
    {
      
      Res <- readRDS("res_insuff_60/res_09_power_insuff.RData")
      
      Res <- bind_rows(Res)
      Res$n <- seq(6, 500, by = 1) 
      
      Res_plot <-  Res %>% 
        pivot_longer(
          cols = c("power_EC", "power_AIPW", "power_naive"), 
          names_to = "Estimator",
          values_to = "power"
        )
      
      Res_plot$Estimator <- factor(Res_plot$Estimator, 
                                   levels = c("power_EC",
                                              "power_AIPW",
                                              "power_naive"
                                   ))
      levels(Res_plot$Estimator) <- c("Hybrid", 
                                      "AIPW",
                                      "Diff-in-mean")
      
      
      p_insuff_60_power_09 <- ggplot(data = Res_plot, aes(x = n, y = power, colour = Estimator)) +
        geom_line(linewidth = 1) +
        scale_colour_manual(values=c("#F27127FF",
                                     "#FDC718FF",
                                     "#038C33FF"),
                            breaks=c("Diff-in-mean",
                                     "AIPW",
                                     "Hybrid"))+
        geom_hline(yintercept = 0.8,
                   linetype="dashed",
                   color = "#084594",
                   linewidth = 1) +
        # geom_vline(xintercept= 141,
        #            linetype="dashed",
        #            color = "#011C40FF",
        #            linewidth = 1)+
        geom_vline(xintercept= 153,
                   linetype="dashed",
                   color =  "#005A32",
                   linewidth = 1)+
        # geom_vline(xintercept= 145,
        #            linetype="dashed",
        #            color = "#99000D",
        #            linewidth = 1)+
        geom_vline(xintercept= 437,
                   linetype="dashed",
                   color = "#FEB24C",
                   linewidth = 1)+
        # theme(legend.position = "none")+
        labs(x=expression(Sample~Size~N[R]), y="Power")+
        scale_x_continuous(breaks=seq(0, 500, by=100))+
        scale_y_continuous(breaks=seq(0, 1, by = 0.2))+
        theme(title = element_text(size = 18),
              axis.title=element_text(size = 16),
              legend.title=element_text(size = 18),
              axis.text=element_text(size = 18),
              legend.text=element_text(size = 18))+
        labs(title = expression(Insufficient~EC~(N[E]==60)~","~pi[A] ==0.9))
      p_insuff_60_power_09
      
    }
    
    
  }
  
  ## Type I error insufficient
  {
    
    # Type I error - 05
    {
      
      Res <- readRDS("res_insuff_60/res_05_TypeI_insuff.RData")
      
      Res <- bind_rows(Res)
      Res$n <- seq(6, 500, by = 1) 
      
      Res_plot <-  Res %>% 
        pivot_longer(
          cols = c("power_EC", "power_AIPW", "power_naive"),
          names_to = "Estimator",
          values_to = "Type_I_Error"
        )
      
      Res_plot$Estimator <- factor(Res_plot$Estimator, 
                                   levels = c("power_EC",
                                              "power_AIPW",
                                              "power_naive"
                                   ))
      levels(Res_plot$Estimator) <- c("Hybrid", 
                                      "AIPW",
                                      "Diff-in-mean")
      
      
      p_insuff_60_TypeI_05 <- ggplot(data = Res_plot, aes(x = n, y = Type_I_Error, colour = Estimator)) +
        geom_line(linewidth = 1) +
        scale_colour_manual(values=c("#F27127FF",
                                      "#FDC718FF",
                                      "#038C33FF"),
                            breaks=c("Diff-in-mean",
                                     "AIPW",
                                     "Hybrid"))+
        geom_hline(yintercept = 0.05,
                   linetype="dashed",
                   color = "#084594",
                   linewidth = 1) +
        # geom_vline(xintercept= 119,
        #            linetype="dashed",
        #            color = "#011C40FF",
        #            linewidth = 1)+
        geom_vline(xintercept= 126,
                   linetype="dashed",
                   color =  "#005A32",
                   linewidth = 1)+
        # geom_vline(xintercept= 145,
        #            linetype="dashed",
        #            color = "#99000D",
        #            linewidth = 1)+
        geom_vline(xintercept= 157,
                   linetype="dashed",
                   color = "#FEB24C",
                   linewidth = 1)+
        theme(legend.position = "none")+
        labs(x=expression(Sample~Size~N[R]), y="Type I Error")+
        scale_x_continuous(breaks=seq(0, 500, by=100))+
        scale_y_continuous(breaks = seq(0, 0.546, by=0.1))+
        expand_limits(y=c(0, 0.546))+
        theme(title = element_text(size = 18),
              axis.title=element_text(size = 16),
              legend.title=element_text(size = 18),
              axis.text=element_text(size = 18),
              legend.text=element_text(size = 18))+
        labs(title = expression(Insufficient~EC~(N[E]==60)~","~pi[A] ==0.5))
      
      p_insuff_60_TypeI_05
      
    }
    
    # Type I error - 06
    {
      
      Res <- readRDS("res_insuff_60/res_06_TypeI_insuff.RData")
      
      Res <- bind_rows(Res)
      Res$n <- seq(6, 500, by = 1)  
      
      
      Res_plot <-  Res %>% 
        pivot_longer(
          cols = c("power_EC", "power_AIPW", "power_naive"),
          names_to = "Estimator",
          values_to = "Type_I_Error"
        )
      
      Res_plot$Estimator <- factor(Res_plot$Estimator, 
                                   levels = c("power_EC",
                                              "power_AIPW",
                                              "power_naive"
                                   ))
      levels(Res_plot$Estimator) <- c("Hybrid", 
                                      "AIPW",
                                      "Diff-in-mean")
      
      
      p_insuff_60_TypeI_06 <- ggplot(data = Res_plot, aes(x = n, y = Type_I_Error, colour = Estimator)) +
        geom_line(linewidth = 1) +
        scale_colour_manual(values=c("#F27127FF",
                                     "#FDC718FF",
                                     "#038C33FF"),
                            breaks=c("Diff-in-mean",
                                     "AIPW",
                                     "Hybrid"))+
        geom_hline(yintercept = 0.05,
                   linetype="dashed",
                   color = "#084594",
                   linewidth = 1) +
        # geom_vline(xintercept= 111,
        #            linetype="dashed",
        #            color = "#011C40FF",
        #            linewidth = 1)+
        geom_vline(xintercept= 118,
                   linetype="dashed",
                   color =  "#005A32",
                   linewidth = 1)+
        # geom_vline(xintercept= 134,
        #            linetype="dashed",
        #            color = "#99000D",
        #            linewidth = 1)+
        geom_vline(xintercept= 164,
                   linetype="dashed",
                   color = "#FEB24C",
                   linewidth = 1)+
        theme(legend.position = "none")+
        labs(x=expression(Sample~Size~N[R]), y="Type I Error")+
        scale_x_continuous(breaks=seq(0, 500, by=100))+
        scale_y_continuous(breaks = seq(0, 0.546, by=0.1))+
        expand_limits(y=c(0, 0.546))+
        theme(title = element_text(size = 18),
              axis.title=element_text(size = 16),
              legend.title=element_text(size = 18),
              axis.text=element_text(size = 18),
              legend.text=element_text(size = 18))+
        labs(title = expression(Insufficient~EC~(N[E]==60)~","~pi[A] ==0.6))
      
      p_insuff_60_TypeI_06
      
    }
    
    # Type I error - 07
    {
      
      Res <- readRDS("res_insuff_60/res_07_TypeI_insuff.RData")
      
      Res <- bind_rows(Res)
      Res$n <- seq(6, 500, by = 1) 
      
      
      Res_plot <-  Res %>% 
        pivot_longer(
          cols = c("power_EC", "power_AIPW", "power_naive"),
          names_to = "Estimator",
          values_to = "Type_I_Error"
        )
      
      Res_plot$Estimator <- factor(Res_plot$Estimator, 
                                   levels = c("power_EC",
                                              "power_AIPW",
                                              "power_naive"
                                   ))
      levels(Res_plot$Estimator) <- c("Hybrid", 
                                      "AIPW",
                                      "Diff-in-mean")
      
      p_insuff_60_TypeI_07 <- ggplot(data = Res_plot, aes(x = n, y = Type_I_Error, colour = Estimator)) +
        geom_line(linewidth = 1) +
        scale_colour_manual(values=c("#F27127FF",
                                     "#FDC718FF",
                                     "#038C33FF"),
                            breaks=c("Diff-in-mean",
                                     "AIPW",
                                     "Hybrid"))+
        geom_hline(yintercept = 0.05,
                   linetype="dashed",
                   color = "#084594",
                   linewidth = 1) +
        # geom_vline(xintercept= 108,
        #            linetype="dashed",
        #            color = "#011C40FF",
        #            linewidth = 1)+
        geom_vline(xintercept= 116,
                   linetype="dashed",
                   color =  "#005A32",
                   linewidth = 1)+
        # geom_vline(xintercept= 128,
        #            linetype="dashed",
        #            color = "#99000D",
        #            linewidth = 1)+
        geom_vline(xintercept= 187,
                   linetype="dashed",
                   color = "#FEB24C",
                   linewidth = 1)+
        theme(legend.position = "none")+
        labs(x=expression(Sample~Size~N[R]), y="Type I Error")+
        scale_x_continuous(breaks=seq(0, 500, by=100))+
        scale_y_continuous(breaks = seq(0, 0.546, by=0.1))+
        expand_limits(y=c(0, 0.546))+
        theme(title = element_text(size = 18),
              axis.title=element_text(size = 16),
              legend.title=element_text(size = 18),
              axis.text=element_text(size = 18),
              legend.text=element_text(size = 18))+
        labs(title = expression(Insufficient~EC~(N[E]==60)~","~pi[A] ==0.7))
      
      p_insuff_60_TypeI_07
      
    }
    
    # Type I error - 08
    {
      
      Res <- readRDS("res_insuff_60/res_08_TypeI_insuff.RData")
      
      Res <- bind_rows(Res)
      Res$n <- seq(6, 500, by = 1)  
      
      
      Res_plot <-  Res %>% 
        pivot_longer(
          cols = c("power_EC", "power_AIPW", "power_naive"),
          names_to = "Estimator",
          values_to = "Type_I_Error"
        )
      
      Res_plot$Estimator <- factor(Res_plot$Estimator, 
                                   levels = c("power_EC",
                                              "power_AIPW",
                                              "power_naive"
                                   ))
      levels(Res_plot$Estimator) <- c("Hybrid", 
                                      "AIPW",
                                      "Diff-in-mean")
      
      
      p_insuff_60_TypeI_08 <- ggplot(data = Res_plot, aes(x = n, y = Type_I_Error, colour = Estimator)) +
        geom_line(linewidth = 1) +
        scale_colour_manual(values=c("#F27127FF",
                                     "#FDC718FF",
                                     "#038C33FF"),
                            breaks=c("Diff-in-mean",
                                     "AIPW",
                                     "Hybrid"))+
        geom_hline(yintercept = 0.05,
                   linetype="dashed",
                   color = "#084594",
                   linewidth = 1) +
        # geom_vline(xintercept= 116,
        #            linetype="dashed",
        #            color = "#011C40FF",
        #            linewidth = 1)+
        geom_vline(xintercept= 124,
                   linetype="dashed",
                   color =  "#005A32",
                   linewidth = 1)+
        # geom_vline(xintercept= 131,
        #            linetype="dashed",
        #            color = "#99000D",
        #            linewidth = 1)+
        geom_vline(xintercept= 246,
                   linetype="dashed",
                   color = "#FEB24C",
                   linewidth = 1)+
        # theme(legend.position = "none")+
        labs(x=expression(Sample~Size~N[R]), y="Type I Error")+
        scale_x_continuous(breaks=seq(0, 500, by=100))+
        scale_y_continuous(breaks = seq(0, 0.546, by=0.1))+
        expand_limits(y=c(0, 0.546))+
        theme(title = element_text(size = 18),
              axis.title=element_text(size = 16),
              legend.title=element_text(size = 18),
              axis.text=element_text(size = 18),
              legend.text=element_text(size = 18))+
        labs(title = expression(Insufficient~EC~(N[E]==60)~","~pi[A] ==0.8))
      
      p_insuff_60_TypeI_08
      
    }
    
    # Type I error - 09
    {
      Res <- readRDS("res_insuff_60/res_09_TypeI_insuff.RData")
      
      Res <- bind_rows(Res)
      Res$n <- seq(6, 500, by = 1) 
      
      
      Res_plot <-  Res %>% 
        pivot_longer(
          cols = c("power_EC", "power_AIPW", "power_naive"),
          names_to = "Estimator",
          values_to = "Type_I_Error"
        )
      
      Res_plot$Estimator <- factor(Res_plot$Estimator, 
                                   levels = c("power_EC",
                                              "power_AIPW",
                                              "power_naive"
                                   ))
      levels(Res_plot$Estimator) <- c("Hybrid", 
                                      "AIPW",
                                      "Diff-in-mean")
      
      
      p_insuff_60_TypeI_09 <- ggplot(data = Res_plot, aes(x = n, y = Type_I_Error, colour = Estimator)) +
        geom_line(linewidth = 1) +
        scale_colour_manual(values=c("#F27127FF",
                                     "#FDC718FF",
                                     "#038C33FF"),
                            breaks=c("Diff-in-mean",
                                     "AIPW",
                                     "Hybrid"))+
        geom_hline(yintercept = 0.05,
                   linetype="dashed",
                   color = "#084594",
                   linewidth = 1) +
        # geom_vline(xintercept= 141,
        #            linetype="dashed",
        #            color = "#011C40FF",
        #            linewidth = 1)+
        geom_vline(xintercept= 153,
                   linetype="dashed",
                   color =  "#005A32",
                   linewidth = 1)+
        # geom_vline(xintercept= 145,
        #            linetype="dashed",
        #            color = "#99000D",
        #            linewidth = 1)+
        geom_vline(xintercept= 437,
                   linetype="dashed",
                   color = "#FEB24C",
                   linewidth = 1)+
        # theme(legend.position = "none")+
        labs(x=expression(Sample~Size~N[R]), y="Type I Error")+
        scale_x_continuous(breaks=seq(0, 500, by=100))+
        scale_y_continuous(breaks = seq(0, 0.546, by=0.1))+
        expand_limits(y=c(0, 0.546))+
        theme(title = element_text(size = 18),
              axis.title=element_text(size = 16),
              legend.title=element_text(size = 18),
              axis.text=element_text(size = 18),
              legend.text=element_text(size = 18))+
        labs(title = expression(Insufficient~EC~(N[E]==60)~","~pi[A] ==0.9))
      
      p_insuff_60_TypeI_09
      
    }
    
  }
}

### Main paper - Insufficient EC (N_E = 30)
{
  ## Power_not sufficient
  {
    
    # Power  -  05
    {
      
      Res <- readRDS("res_insuff_30/res_05_power_insuff.RData")
      
      Res <- bind_rows(Res)
      Res$n <- seq(6, 500, by = 1) 
      
      Res_plot <-  Res %>% 
        pivot_longer(
          cols = c("power_EC", "power_AIPW", "power_naive"), 
          names_to = "Estimator",
          values_to = "power"
        )
      
      Res_plot$Estimator <- factor(Res_plot$Estimator, 
                                   levels = c("power_EC",
                                              "power_AIPW",
                                              "power_naive"
                                   ))
      levels(Res_plot$Estimator) <- c("Hybrid", 
                                      "AIPW",
                                      "Diff-in-mean")
      
      
      p_insuff_30_power_05 <- ggplot(data = Res_plot, aes(x = n, y = power, colour = Estimator)) +
        geom_line(linewidth = 1) +
        scale_colour_manual(values=c("#F27127FF",
                                     "#FDC718FF",
                                     "#038C33FF"),
                            breaks=c("Diff-in-mean",
                                     "AIPW",
                                     "Hybrid"))+
        geom_hline(yintercept = 0.8,
                   linetype="dashed",
                   color = "#084594",
                   linewidth = 1) +
        # geom_vline(xintercept= 119,
        #            linetype="dashed",
        #            color = "#011C40FF",
        #            linewidth = 1)+
        geom_vline(xintercept= 138,
                   linetype="dashed",
                   color =  "#005A32",
                   linewidth = 1) +
        # geom_vline(xintercept= 145,
        #            linetype="dashed",
        #            color = "#99000D",
        #            linewidth = 1)+
        geom_vline(xintercept= 157,
                   linetype="dashed",
                   color = "#FEB24C",
                   linewidth = 1)+
        theme(legend.position = "none")+
        labs(x=expression(Sample~Size~N[R]), y="Power")+
        scale_x_continuous(breaks=seq(0, 500, by=100))+
        scale_y_continuous(breaks=seq(0, 1, by = 0.2))+
        theme(title = element_text(size = 20),
              axis.title=element_text(size = 20),
              legend.title=element_text(size = 20),
              axis.text=element_text(size = 20),
              legend.text=element_text(size = 20))+
        labs(title = expression(Insufficient~EC~(N[E]==30)~","~pi[A] ==0.5))
      p_insuff_30_power_05
      
    }
    
    
    # Power  -  06
    {
      
      Res <- readRDS("res_insuff_30/res_06_power_insuff.RData")
      
      Res <- bind_rows(Res)
      Res$n <- seq(6, 500, by = 1) 
      
      Res_plot <-  Res %>% 
        pivot_longer(
          cols = c("power_EC", "power_AIPW", "power_naive"), 
          names_to = "Estimator",
          values_to = "power"
        )
      
      Res_plot$Estimator <- factor(Res_plot$Estimator, 
                                   levels = c("power_EC",
                                              "power_AIPW",
                                              "power_naive"
                                   ))
      levels(Res_plot$Estimator) <- c("Hybrid", 
                                      "AIPW",
                                      "Diff-in-mean")
      
      
      p_insuff_30_power_06 <- ggplot(data = Res_plot, aes(x = n, y = power, colour = Estimator)) +
        geom_line(linewidth = 1) +
        scale_colour_manual(values=c("#F27127FF",
                                     "#FDC718FF",
                                     "#038C33FF"),
                            breaks=c("Diff-in-mean",
                                     "AIPW",
                                     "Hybrid"))+
        geom_hline(yintercept = 0.8,
                   linetype="dashed",
                   color = "#084594",
                   linewidth = 1) +
        geom_vline(xintercept= 136,
                   linetype="dashed",
                   color =  "#005A32",
                   linewidth = 1)+
        # geom_vline(xintercept= 134,
        #            linetype="dashed",
        #            color = "#99000D",
        #            linewidth = 1)+
        geom_vline(xintercept= 164,
                   linetype="dashed",
                   color = "#FEB24C",
                   linewidth = 1)+
        theme(legend.position = "none")+
        labs(x=expression(Sample~Size~N[R]), y="Power")+
        scale_x_continuous(breaks=seq(0, 500, by=100))+
        scale_y_continuous(breaks=seq(0, 1, by = 0.2))+
        theme(title = element_text(size = 20),
              axis.title=element_text(size = 20),
              legend.title=element_text(size = 20),
              axis.text=element_text(size = 20),
              legend.text=element_text(size = 20))+
        labs(title = expression(Insufficient~EC~(N[E]==30)~","~pi[A] ==0.6))
      p_insuff_30_power_06
      
    }
    
    
    # Power  -  07
    {
      
      Res <- readRDS("res_insuff_30/res_07_power_insuff.RData")
      
      Res <- bind_rows(Res)
      Res$n <- seq(6, 500, by = 1) 
      
      Res_plot <-  Res %>% 
        pivot_longer(
          cols = c("power_EC", "power_AIPW", "power_naive"), 
          names_to = "Estimator",
          values_to = "power"
        )
      
      Res_plot$Estimator <- factor(Res_plot$Estimator, 
                                   levels = c("power_EC",
                                              "power_AIPW",
                                              "power_naive"
                                   ))
      levels(Res_plot$Estimator) <- c("Hybrid", 
                                      "AIPW",
                                      "Diff-in-mean")
      
      
      p_insuff_30_power_07 <- ggplot(data = Res_plot, aes(x = n, y = power, colour = Estimator)) +
        geom_line(linewidth = 1) +
        scale_colour_manual(values=c("#F27127FF",
                                     "#FDC718FF",
                                     "#038C33FF"),
                            breaks=c("Diff-in-mean",
                                     "AIPW",
                                     "Hybrid"))+
        geom_hline(yintercept = 0.8,
                   linetype="dashed",
                   color = "#084594",
                   linewidth = 1) +
        # geom_vline(xintercept= 108,
        #            linetype="dashed",
        #            color = "#011C40FF",
        #            linewidth = 1)+
        geom_vline(xintercept= 143,
                   linetype="dashed",
                   color =  "#005A32",
                   linewidth = 1)+
        # geom_vline(xintercept= 128,
        #            linetype="dashed",
        #            color = "#99000D",
        #            linewidth = 1)+
        geom_vline(xintercept= 187,
                   linetype="dashed",
                   color = "#FEB24C",
                   linewidth = 1)+
        theme(legend.position = "none")+
        labs(x=expression(Sample~Size~N[R]), y="Power")+
        scale_x_continuous(breaks=seq(0, 500, by=100))+
        scale_y_continuous(breaks=seq(0, 1, by = 0.2))+
        theme(title = element_text(size = 20),
              axis.title=element_text(size = 20),
              legend.title=element_text(size = 20),
              axis.text=element_text(size = 20),
              legend.text=element_text(size = 20))+
        labs(title = expression(Insufficient~EC~(N[E]==30)~","~pi[A] ==0.7))
      p_insuff_30_power_07
      
      
    }
    
    
    # Power  -  08
    {
      
      Res <- readRDS("res_insuff_30/res_08_power_insuff.RData")
      
      Res <- bind_rows(Res)
      Res$n <- seq(6, 500, by = 1) 
      
      Res_plot <-  Res %>% 
        pivot_longer(
          cols = c("power_EC", "power_AIPW", "power_naive"), 
          names_to = "Estimator",
          values_to = "power"
        )
      
      Res_plot$Estimator <- factor(Res_plot$Estimator, 
                                   levels = c("power_EC",
                                              "power_AIPW",
                                              "power_naive"
                                   ))
      levels(Res_plot$Estimator) <- c("Hybrid", 
                                      "AIPW",
                                      "Diff-in-mean")
      
      
      p_insuff_30_power_08 <- ggplot(data = Res_plot, aes(x = n, y = power, colour = Estimator)) +
        geom_line(linewidth = 1) +
        scale_colour_manual(values=c("#F27127FF",
                                     "#FDC718FF",
                                     "#038C33FF"),
                            breaks=c("Diff-in-mean",
                                     "AIPW",
                                     "Hybrid"))+
        geom_hline(yintercept = 0.8,
                   linetype="dashed",
                   color = "#084594",
                   linewidth = 1) +
        # geom_vline(xintercept= 116,
        #            linetype="dashed",
        #            color = "#011C40FF",
        #            linewidth = 1)+
        geom_vline(xintercept= 170,
                   linetype="dashed",
                   color =  "#005A32",
                   linewidth = 1)+
        # geom_vline(xintercept= 131,
        #            linetype="dashed",
        #            color = "#99000D",
        #            linewidth = 1)+
        geom_vline(xintercept= 246,
                   linetype="dashed",
                   color = "#FEB24C",
                   linewidth = 1)+
        theme(legend.position = "none")+
        labs(x=expression(Sample~Size~N[R]), y="Power")+
        scale_x_continuous(breaks=seq(0, 500, by=100))+
        scale_y_continuous(breaks=seq(0, 1, by = 0.2))+
        theme(title = element_text(size = 20),
              axis.title=element_text(size = 20),
              legend.title=element_text(size = 20),
              axis.text=element_text(size = 20),
              legend.text=element_text(size = 20))+
        labs(title = expression(Insufficient~EC~(N[E]==30)~","~pi[A] ==0.8))
      p_insuff_30_power_08
      
      
    }
    
    
    # Power  -  09
    {
      
      Res <- readRDS("res_insuff_30/res_09_power_insuff.RData")
      
      Res <- bind_rows(Res)
      Res$n <- seq(6, 500, by = 1) 
      
      Res_plot <-  Res %>% 
        pivot_longer(
          cols = c("power_EC", "power_AIPW", "power_naive"), 
          names_to = "Estimator",
          values_to = "power"
        )
      
      Res_plot$Estimator <- factor(Res_plot$Estimator, 
                                   levels = c("power_EC",
                                              "power_AIPW",
                                              "power_naive"
                                   ))
      levels(Res_plot$Estimator) <- c("Hybrid", 
                                      "AIPW",
                                      "Diff-in-mean")
      
      
      p_insuff_30_power_09 <- ggplot(data = Res_plot, aes(x = n, y = power, colour = Estimator)) +
        geom_line(linewidth = 1) +
        scale_colour_manual(values=c("#F27127FF",
                                     "#FDC718FF",
                                     "#038C33FF"),
                            breaks=c("Diff-in-mean",
                                     "AIPW",
                                     "Hybrid"))+
        geom_hline(yintercept = 0.8,
                   linetype="dashed",
                   color = "#084594",
                   linewidth = 1) +
        # geom_vline(xintercept= 141,
        #            linetype="dashed",
        #            color = "#011C40FF",
        #            linewidth = 1)+
        geom_vline(xintercept= 261,
                   linetype="dashed",
                   color =  "#005A32",
                   linewidth = 1)+
        # geom_vline(xintercept= 145,
        #            linetype="dashed",
        #            color = "#99000D",
        #            linewidth = 1)+
        geom_vline(xintercept= 437,
                   linetype="dashed",
                   color = "#FEB24C",
                   linewidth = 1)+
        # theme(legend.position = "none")+
        labs(x=expression(Sample~Size~N[R]), y="Power")+
        scale_x_continuous(breaks=seq(0, 500, by=100))+
        scale_y_continuous(breaks=seq(0, 1, by = 0.2))+
        theme(title = element_text(size = 20),
              axis.title=element_text(size = 20),
              legend.title=element_text(size = 20),
              axis.text=element_text(size = 20),
              legend.text=element_text(size = 20))+
        labs(title = expression(Insufficient~EC~(N[E]==30)~","~pi[A] ==0.9))
      p_insuff_30_power_09
      
      
    }
    
    
  }
  
  ## Type I error_not sufficient
  {
    
    # Type I error - 05
    {
      
      Res <- readRDS("res_insuff_30/res_05_TypeI_insuff.RData")
      
      Res <- bind_rows(Res)
      Res$n <- seq(6, 500, by = 1) 
      
      Res_plot <-  Res %>% 
        pivot_longer(
          cols = c("power_EC", "power_AIPW", "power_naive"),
          names_to = "Estimator",
          values_to = "Type_I_Error"
        )
      
      Res_plot$Estimator <- factor(Res_plot$Estimator, 
                                   levels = c("power_EC",
                                              "power_AIPW",
                                              "power_naive"
                                   ))
      levels(Res_plot$Estimator) <- c("Hybrid", 
                                      "AIPW",
                                      "Diff-in-mean")
      
      
      p_insuff_30_TypeI_05 <- ggplot(data = Res_plot, aes(x = n, y = Type_I_Error, colour = Estimator)) +
        geom_line(linewidth = 1) +
        scale_colour_manual(values=c("#F27127FF",
                                     "#FDC718FF",
                                     "#038C33FF"),
                            breaks=c("Diff-in-mean",
                                     "AIPW",
                                     "Hybrid"))+
        geom_hline(yintercept = 0.05,
                   linetype="dashed",
                   color = "#084594",
                   linewidth = 1) +
        # geom_vline(xintercept= 119,
        #            linetype="dashed",
        #            color = "#011C40FF",
        #            linewidth = 1)+
        geom_vline(xintercept= 138,
                   linetype="dashed",
                   color =  "#005A32",
                   linewidth = 1)+
        # geom_vline(xintercept= 145,
        #            linetype="dashed",
        #            color = "#99000D",
        #            linewidth = 1)+
        geom_vline(xintercept= 157,
                   linetype="dashed",
                   color = "#FEB24C",
                   linewidth = 1)+
        theme(legend.position = "none")+
        labs(x=expression(Sample~Size~N[R]), y="Type I Error")+
        scale_x_continuous(breaks=seq(0, 500, by=100))+
        scale_y_continuous(breaks = seq(0, 0.546, by=0.1))+
        expand_limits(y=c(0, 0.546))+
        theme(title = element_text(size = 20),
              axis.title=element_text(size = 20),
              legend.title=element_text(size = 20),
              axis.text=element_text(size = 20),
              legend.text=element_text(size = 20))+
        labs(title = expression(Insufficient~EC~(N[E]==30)~","~pi[A] ==0.5))
      
      p_insuff_30_TypeI_05
      
    }
    
    # Type I error - 06
    {
      
      Res <- readRDS("res_insuff_30/res_06_TypeI_insuff.RData")
      
      Res <- bind_rows(Res)
      Res$n <- seq(6, 500, by = 1)  
      
      
      Res_plot <-  Res %>% 
        pivot_longer(
          cols = c("power_EC", "power_AIPW", "power_naive"),
          names_to = "Estimator",
          values_to = "Type_I_Error"
        )
      
      Res_plot$Estimator <- factor(Res_plot$Estimator, 
                                   levels = c("power_EC",
                                              "power_AIPW",
                                              "power_naive"
                                   ))
      levels(Res_plot$Estimator) <- c("Hybrid", 
                                      "AIPW",
                                      "Diff-in-mean")
      
      
      p_insuff_30_TypeI_06 <- ggplot(data = Res_plot, aes(x = n, y = Type_I_Error, colour = Estimator)) +
        geom_line(linewidth = 1) +
        scale_colour_manual(values=c("#F27127FF",
                                     "#FDC718FF",
                                     "#038C33FF"),
                            breaks=c("Diff-in-mean",
                                     "AIPW",
                                     "Hybrid"))+
        geom_hline(yintercept = 0.05,
                   linetype="dashed",
                   color = "#084594",
                   linewidth = 1) +
        # geom_vline(xintercept= 111,
        #            linetype="dashed",
        #            color = "#011C40FF",
        #            linewidth = 1)+
        geom_vline(xintercept= 136,
                   linetype="dashed",
                   color =  "#005A32",
                   linewidth = 1)+
        # geom_vline(xintercept= 134,
        #            linetype="dashed",
        #            color = "#99000D",
        #            linewidth = 1)+
        geom_vline(xintercept= 164,
                   linetype="dashed",
                   color = "#FEB24C",
                   linewidth = 1)+
        theme(legend.position = "none")+
        labs(x=expression(Sample~Size~N[R]), y="Type I Error")+
        scale_x_continuous(breaks=seq(0, 500, by=100))+
        scale_y_continuous(breaks = seq(0, 0.546, by=0.1))+
        expand_limits(y=c(0, 0.546))+
        theme(title = element_text(size = 20),
              axis.title=element_text(size = 20),
              legend.title=element_text(size = 20),
              axis.text=element_text(size = 20),
              legend.text=element_text(size = 20))+
        labs(title = expression(Insufficient~EC~(N[E]==30)~","~pi[A] ==0.6))
      
      p_insuff_30_TypeI_06
      
    }
    
    # Type I error - 07
    {
      
      Res <- readRDS("res_insuff_30/res_07_TypeI_insuff.RData")
      
      Res <- bind_rows(Res)
      Res$n <- seq(6, 500, by = 1) 
      
      
      Res_plot <-  Res %>% 
        pivot_longer(
          cols = c("power_EC", "power_AIPW", "power_naive"),
          names_to = "Estimator",
          values_to = "Type_I_Error"
        )
      
      Res_plot$Estimator <- factor(Res_plot$Estimator, 
                                   levels = c("power_EC",
                                              "power_AIPW",
                                              "power_naive"
                                   ))
      levels(Res_plot$Estimator) <- c("Hybrid", 
                                      "AIPW",
                                      "Diff-in-mean")
      
      p_insuff_30_TypeI_07 <- ggplot(data = Res_plot, aes(x = n, y = Type_I_Error, colour = Estimator)) +
        geom_line(linewidth = 1) +
        scale_colour_manual(values=c("#F27127FF",
                                     "#FDC718FF",
                                     "#038C33FF"),
                            breaks=c("Diff-in-mean",
                                     "AIPW",
                                     "Hybrid"))+
        geom_hline(yintercept = 0.05,
                   linetype="dashed",
                   color = "#084594",
                   linewidth = 1) +
        # geom_vline(xintercept= 108,
        #            linetype="dashed",
        #            color = "#011C40FF",
        #            linewidth = 1)+
        geom_vline(xintercept= 143,
                   linetype="dashed",
                   color =  "#005A32",
                   linewidth = 1)+
        # geom_vline(xintercept= 128,
        #            linetype="dashed",
        #            color = "#99000D",
        #            linewidth = 1)+
        geom_vline(xintercept= 187,
                   linetype="dashed",
                   color = "#FEB24C",
                   linewidth = 1)+
        theme(legend.position = "none")+
        labs(x=expression(Sample~Size~N[R]), y="Type I Error")+
        scale_x_continuous(breaks=seq(0, 500, by=100))+
        scale_y_continuous(breaks = seq(0, 0.546, by=0.1))+
        expand_limits(y=c(0, 0.546))+
        theme(title = element_text(size = 20),
              axis.title=element_text(size = 20),
              legend.title=element_text(size = 20),
              axis.text=element_text(size = 20),
              legend.text=element_text(size = 20))+
        labs(title = expression(Insufficient~EC~(N[E]==30)~","~pi[A] ==0.7))
      
      p_insuff_30_TypeI_07
      
    }
    
    # Type I error - 08
    {
      
      Res <- readRDS("res_insuff_30/res_08_TypeI_insuff.RData")
      
      Res <- bind_rows(Res)
      Res$n <- seq(6, 500, by = 1)  
      
      
      Res_plot <-  Res %>% 
        pivot_longer(
          cols = c("power_EC", "power_AIPW", "power_naive"),
          names_to = "Estimator",
          values_to = "Type_I_Error"
        )
      
      Res_plot$Estimator <- factor(Res_plot$Estimator, 
                                   levels = c("power_EC",
                                              "power_AIPW",
                                              "power_naive"
                                   ))
      levels(Res_plot$Estimator) <- c("Hybrid", 
                                      "AIPW",
                                      "Diff-in-mean")
      
      
      p_insuff_30_TypeI_08 <- ggplot(data = Res_plot, aes(x = n, y = Type_I_Error, colour = Estimator)) +
        geom_line(linewidth = 1) +
        scale_colour_manual(values=c("#F27127FF",
                                     "#FDC718FF",
                                     "#038C33FF"),
                            breaks=c("Diff-in-mean",
                                     "AIPW",
                                     "Hybrid"))+
        geom_hline(yintercept = 0.05,
                   linetype="dashed",
                   color = "#084594",
                   linewidth = 1) +
        # geom_vline(xintercept= 116,
        #            linetype="dashed",
        #            color = "#011C40FF",
        #            linewidth = 1)+
        geom_vline(xintercept= 170,
                   linetype="dashed",
                   color =  "#005A32",
                   linewidth = 1)+
        # geom_vline(xintercept= 131,
        #            linetype="dashed",
        #            color = "#99000D",
        #            linewidth = 1)+
        geom_vline(xintercept= 246,
                   linetype="dashed",
                   color = "#FEB24C",
                   linewidth = 1)+
        theme(legend.position = "none")+
        labs(x=expression(Sample~Size~N[R]), y="Type I Error")+
        scale_x_continuous(breaks=seq(0, 500, by=100))+
        scale_y_continuous(breaks = seq(0, 0.546, by=0.1))+
        expand_limits(y=c(0, 0.546))+
        theme(title = element_text(size = 20),
              axis.title=element_text(size = 20),
              legend.title=element_text(size = 20),
              axis.text=element_text(size = 20),
              legend.text=element_text(size = 20))+
        labs(title = expression(Insufficient~EC~(N[E]==30)~","~pi[A] ==0.8))
      
      p_insuff_30_TypeI_08
      
    }
    
    # Type I error - 09
    {
      
      Res <- readRDS("res_insuff_30/res_09_TypeI_insuff.RData")
      
      Res <- bind_rows(Res)
      Res$n <- seq(6, 500, by = 1) 
      
      
      Res_plot <-  Res %>% 
        pivot_longer(
          cols = c("power_EC", "power_AIPW", "power_naive"),
          names_to = "Estimator",
          values_to = "Type_I_Error"
        )
      
      Res_plot$Estimator <- factor(Res_plot$Estimator, 
                                   levels = c("power_EC",
                                              "power_AIPW",
                                              "power_naive"
                                   ))
      levels(Res_plot$Estimator) <- c("Hybrid", 
                                      "AIPW",
                                      "Diff-in-mean")
      
      
      p_insuff_30_TypeI_09 <- ggplot(data = Res_plot, aes(x = n, y = Type_I_Error, colour = Estimator)) +
        geom_line(linewidth = 1) +
        scale_colour_manual(values=c("#F27127FF",
                                     "#FDC718FF",
                                     "#038C33FF"),
                            breaks=c("Diff-in-mean",
                                     "AIPW",
                                     "Hybrid"))+
        geom_hline(yintercept = 0.05,
                   linetype="dashed",
                   color = "#084594",
                   linewidth = 1) +
        # geom_vline(xintercept= 141,
        #            linetype="dashed",
        #            color = "#011C40FF",
        #            linewidth = 1)+
        geom_vline(xintercept= 261,
                   linetype="dashed",
                   color =  "#005A32",
                   linewidth = 1)+
        # geom_vline(xintercept= 145,
        #            linetype="dashed",
        #            color = "#99000D",
        #            linewidth = 1)+
        geom_vline(xintercept= 437,
                   linetype="dashed",
                   color = "#FEB24C",
                   linewidth = 1)+
        # theme(legend.position = "none")+
        labs(x=expression(Sample~Size~N[R]), y="Type I Error")+
        scale_x_continuous(breaks=seq(0, 500, by=100))+
        scale_y_continuous(breaks = seq(0, 0.546, by=0.1))+
        expand_limits(y=c(0, 0.546))+
        theme(title = element_text(size = 20),
              axis.title=element_text(size = 20),
              legend.title=element_text(size = 20),
              axis.text=element_text(size = 20),
              legend.text=element_text(size = 20))+
        labs(title = expression(Insufficient~EC~(N[E]==30)~","~pi[A] ==0.9))
      
      p_insuff_30_TypeI_09
      
    }
    
  }
}

### Combine plots 
{
  
  ##########################################
  ## For main paper simulation section
  p1_power <- ggarrange(p_suff_power_06, p_suff_power_08,
                        p_insuff_60_power_06, p_insuff_60_power_08,
                        # widths = c(1, 1.4),
                        ncol = 2, nrow = 2)
  p1_power
  ggsave(plot = p1_power,
         width = 12,
         height = 11,
         dpi = 500,
         filename = "/code/simulation/plots/p1_power.pdf")

  p1_TypeI <- ggarrange(p_suff_TypeI_06, p_suff_TypeI_08,
                        p_insuff_60_TypeI_06, p_insuff_60_TypeI_08,
                        widths = c(1, 1.4),
                        ncol = 2, nrow = 2)
  p1_TypeI
  ggsave(plot = p1_TypeI,
         width = 14,
         height = 11,
         dpi = 500,
         filename = "/code/simulation/plots/p1_TypeI.pdf")
  
  
  p_simu <- ggarrange(p_suff_power_06, p_suff_power_08, p_suff_TypeI_06, p_suff_TypeI_08,
                      p_insuff_60_power_06, p_insuff_60_power_08, p_insuff_60_TypeI_06, p_insuff_60_TypeI_08,
                      widths = c(1, 1, 1, 1.4),
                      ncol = 4, nrow = 2)
  p_simu
  ggsave(plot = p_simu,
         width = 25,
         height = 11,
         dpi = 800,
         filename = "/code/simulation/plots/p_simu.pdf")
  
  
  ##########################################
  ## For supplementary material B.1
  p2_power <- ggarrange(p_suff_power_05, p_suff_power_07, p_suff_power_09,
                        p_insuff_60_power_05, p_insuff_60_power_07, p_insuff_60_power_09,
                        widths = c(1, 1, 1.38),
                        ncol = 3, nrow = 2)
  p2_power
  p2_power_new <- annotate_figure(p2_power,
                                  bottom = text_grob(expression((a)~Power~(tau==0.4)), face = "bold", size = 22)
                                  )
  p2_power_new
  
  ggsave(plot = p2_power,
         width = 21,
         height = 11,
         dpi = 500,
         filename = "/code/simulation/plots/p2_power.pdf")
  
  
  p2_TypeI <- ggarrange(p_suff_TypeI_05, p_suff_TypeI_07, p_suff_TypeI_09,
                        p_insuff_60_TypeI_05, p_insuff_60_TypeI_07, p_insuff_60_TypeI_09,
                        widths = c(1, 1, 1.38),
                        ncol = 3, nrow = 2)
  p2_TypeI
  p2_TypeI_new <- annotate_figure(p2_TypeI,
                                  bottom = text_grob(expression((b)~Type~I~Error~(tau==0)), face = "bold", size = 22)
  )
  p2_TypeI_new

  ggsave(plot = p2_TypeI,
         width = 21,
         height = 11,
         dpi = 500,
         filename = "/code/simulation/plots/p2_TypeI.pdf")
  
  
  p2_supp <- ggarrange(p2_power_new, 
                       NULL,   # this row will be pure whitespace
                       p2_TypeI_new, 
                       heights = c(1, 0.06, 1),
                       ncol = 1, nrow = 3)
  p2_supp
  ggsave(plot = p2_supp,
         width = 21,
         height = 24,
         dpi = 800,
         filename = "/code/simulation/plots/p2_supp.pdf")
  
  
  ##########################################
  ## For supplementary material B.2
  p_insuff_30 <- ggarrange(p_insuff_30_power_05, p_insuff_30_power_06, p_insuff_30_power_07, 
                           p_insuff_30_power_08, p_insuff_30_power_09,
                           p_insuff_30_TypeI_05, p_insuff_30_TypeI_06, p_insuff_30_TypeI_07,
                           p_insuff_30_TypeI_08, p_insuff_30_TypeI_09,
                         widths = c(1, 1, 1, 1, 1.38),
                         ncol = 5, nrow = 2)
  p_insuff_30
 
  ggsave(plot = p_insuff_30,
         width = 37.5,
         height = 13.2,
         dpi = 800,
         filename = "/code/simulation/plots/p_insuff_30.pdf")
  
  
}


### Sensitivity analysis 1 (N_E = 1000)
{
  
  ## Power_sufficient
  {
    
    # Power  -  05
    {
      
      Res <- readRDS("res_sensi1/res_05_power_sensi1.RData")
      
      Res <- bind_rows(Res)
      Res$n <- seq(6, 600, by = 1) 
      
      Res_plot <-  Res %>% 
        pivot_longer(
          cols = c("power_EC", "power_AIPW", "power_naive"), 
          names_to = "Estimator",
          values_to = "power"
        )
      
      # add single-arm case
      Res_sa <- readRDS("res_sensi1/res_sa_power_sensi1.RData")
      
      Res_sa <- bind_rows(Res_sa)
      Res_sa$n <- seq(6, 600, by = 1) 
      
      Res_plot_sa <-  Res_sa %>% 
        pivot_longer(
          cols = c("power_SA"), 
          names_to = "Estimator",
          values_to = "power"
        )
      
      Res_plot <- rbind(Res_plot, Res_plot_sa)
      
      Res_plot$Estimator <- factor(Res_plot$Estimator, 
                                   levels = c(
                                     "power_EC",
                                     "power_AIPW",
                                     "power_naive",
                                     "power_SA"
                                   ))
      levels(Res_plot$Estimator) <- c(
        "Hybrid", 
        "AIPW",
        "Diff-in-mean",
        "Single-arm")
      
      p_sensi1_power_05 <- ggplot(data = Res_plot, aes(x = n, y = power, colour = Estimator)) +
        geom_line(linewidth = 1) +
        scale_colour_manual(values=c("#F27127FF","#FDC718FF",
                                     "#038C33FF", "#207CF5FF"),
                            breaks=c("Diff-in-mean", "AIPW",
                                     "Hybrid", "Single-arm"))+
        geom_hline(yintercept = 0.8,
                   linetype="dashed",
                   color = "#084594",
                   linewidth = 1) +
        geom_vline(xintercept= 47,
                   linetype="dashed",
                   color = "#011C40FF",
                   linewidth = 1)+
        geom_vline(xintercept= 88,
                   linetype="dashed",
                   color =  "#005A32",
                   linewidth = 1) +
        geom_vline(xintercept= 162,
                   linetype="dashed",
                   color = "#FEB24C",
                   linewidth = 1)+
        theme(legend.position = "none")+
        labs(x=expression(Sample~Size~N[R]), y="Power")+
        scale_x_continuous(breaks=seq(0, 600, by=100))+
        scale_y_continuous(breaks=seq(0, 1, by = 0.2))+
        theme(title = element_text(size = 20),
              axis.title=element_text(size = 20),
              legend.title=element_text(size = 20),
              axis.text=element_text(size = 20),
              legend.text=element_text(size = 20))+
        labs(title = expression(Heterogeneous~treatment~effect~","~pi[A] ==0.5))
      
      p_sensi1_power_05
       
    }
    
    # Power  -  06
    {
      
      Res <- readRDS("res_sensi1/res_06_power_sensi1.RData")
      
      Res <- bind_rows(Res) 
      Res$n <- seq(6, 600, by = 1) 
      
      Res_plot <-  Res %>% 
        pivot_longer(
          cols = c("power_EC", "power_AIPW", "power_naive"), 
          names_to = "Estimator",
          values_to = "power"
        )
      
      # add single-arm case
      Res_sa <- readRDS("res_sensi1/res_sa_power_sensi1.RData")
      
      Res_sa <- bind_rows(Res_sa)
      Res_sa$n <- seq(6, 600, by = 1) 
      
      Res_plot_sa <-  Res_sa %>% 
        pivot_longer(
          cols = c("power_SA"), 
          names_to = "Estimator",
          values_to = "power"
        )
      
      Res_plot <- rbind(Res_plot, Res_plot_sa)
      
      Res_plot$Estimator <- factor(Res_plot$Estimator, 
                                   levels = c(
                                     "power_EC",
                                     "power_AIPW",
                                     "power_naive",
                                     "power_SA"
                                   ))
      levels(Res_plot$Estimator) <- c(
        "Hybrid", 
        "AIPW",
        "Diff-in-mean",
        "Single-arm")
      
      p_sensi1_power_06 <- ggplot(data = Res_plot, aes(x = n, y = power, colour = Estimator)) +
        geom_line(linewidth = 1) +
        scale_colour_manual(values=c("#F27127FF","#FDC718FF",
                                     "#038C33FF", "#207CF5FF"),
                            breaks=c("Diff-in-mean", "AIPW",
                                     "Hybrid", "Single-arm"))+
        geom_hline(yintercept = 0.8,
                   linetype="dashed",
                   color = "#084594",
                   linewidth = 1) +
        geom_vline(xintercept= 47,
                   linetype="dashed",
                   color = "#011C40FF",
                   linewidth = 1)+
        geom_vline(xintercept= 74,
                   linetype="dashed",
                   color =  "#005A32",
                   linewidth = 1)+
        geom_vline(xintercept= 169,
                   linetype="dashed",
                   color = "#FEB24C",
                   linewidth = 1)+
        theme(legend.position = "none")+
        labs(x=expression(Sample~Size~N[R]), y="Power")+
        scale_x_continuous(breaks=seq(0, 600, by=100))+
        scale_y_continuous(breaks=seq(0, 1, by = 0.2))+
        theme(title = element_text(size = 20),
              axis.title=element_text(size = 20),
              legend.title=element_text(size = 20),
              axis.text=element_text(size = 20),
              legend.text=element_text(size = 20))+
        labs(title = expression(Heterogeneous~treatment~effect~","~pi[A] ==0.6))
      
      p_sensi1_power_06
      
    }
    
    # Power  -  07
    {
      
      Res <- readRDS("res_sensi1/res_07_power_sensi1.RData")
      
      Res <- bind_rows(Res)
      Res$n <- seq(6, 600, by = 1) 
      
      Res_plot <-  Res %>% 
        pivot_longer(
          cols = c("power_EC", "power_AIPW", "power_naive"), 
          names_to = "Estimator",
          values_to = "power"
        )
       
      # add single-arm case
      Res_sa <- readRDS("res_sensi1/res_sa_power_sensi1.RData")
      
      Res_sa <- bind_rows(Res_sa)
      Res_sa$n <- seq(6, 600, by = 1) 
      
      Res_plot_sa <-  Res_sa %>% 
        pivot_longer(
          cols = c("power_SA"), 
          names_to = "Estimator",
          values_to = "power"
        )
      
      Res_plot <- rbind(Res_plot, Res_plot_sa)
      
      Res_plot$Estimator <- factor(Res_plot$Estimator, 
                                   levels = c(
                                     "power_EC",
                                     "power_AIPW",
                                     "power_naive",
                                     "power_SA"
                                   ))
      levels(Res_plot$Estimator) <- c(
        "Hybrid", 
        "AIPW",
        "Diff-in-mean",
        "Single-arm")
      
      
      p_sensi1_power_07 <- ggplot(data = Res_plot, aes(x = n, y = power, colour = Estimator)) +
        geom_line(linewidth = 1) +
        scale_colour_manual(values=c("#F27127FF","#FDC718FF",
                                     "#038C33FF", "#207CF5FF"),
                            breaks=c("Diff-in-mean", "AIPW",
                                     "Hybrid", "Single-arm"))+
        geom_hline(yintercept = 0.8,
                   linetype="dashed",
                   color = "#084594",
                   linewidth = 1) +
        geom_vline(xintercept= 47,
                   linetype="dashed",
                   color = "#011C40FF",
                   linewidth = 1)+
        geom_vline(xintercept= 65,
                   linetype="dashed",
                   color =  "#005A32",
                   linewidth = 1)+
        geom_vline(xintercept= 192,
                   linetype="dashed",
                   color = "#FEB24C",
                   linewidth = 1)+
        theme(legend.position = "none")+
        labs(x=expression(Sample~Size~N[R]), y="Power")+
        scale_x_continuous(breaks=seq(0, 600, by=100))+
        scale_y_continuous(breaks=seq(0, 1, by = 0.2))+
        theme(title = element_text(size = 20),
              axis.title=element_text(size = 20),
              legend.title=element_text(size = 20),
              axis.text=element_text(size = 20),
              legend.text=element_text(size = 20))+
        labs(title = expression(Heterogeneous~treatment~effect~","~pi[A] ==0.7))
      
      p_sensi1_power_07
      
    }
    
    # Power  -  08
    {
      
      Res <- readRDS("res_sensi1/res_08_power_sensi1.RData")
      
      Res <- bind_rows(Res) 
      Res$n <- seq(6, 600, by = 1) 
      
      Res_plot <-  Res %>% 
        pivot_longer(
          cols = c("power_EC", "power_AIPW", "power_naive"), 
          names_to = "Estimator",
          values_to = "power"
        )
      
      # add single-arm case
      Res_sa <- readRDS("res_sensi1/res_sa_power_sensi1.RData")
      
      Res_sa <- bind_rows(Res_sa)
      Res_sa$n <- seq(6, 600, by = 1) 
      
      Res_plot_sa <-  Res_sa %>% 
        pivot_longer(
          cols = c("power_SA"), 
          names_to = "Estimator",
          values_to = "power"
        )
      
      Res_plot <- rbind(Res_plot, Res_plot_sa)
      
      Res_plot$Estimator <- factor(Res_plot$Estimator, 
                                   levels = c(
                                     "power_EC",
                                     "power_AIPW",
                                     "power_naive",
                                     "power_SA"
                                   ))
      levels(Res_plot$Estimator) <- c(
        "Hybrid", 
        "AIPW",
        "Diff-in-mean",
        "Single-arm")
      
      
      p_sensi1_power_08 <- ggplot(data = Res_plot, aes(x = n, y = power, colour = Estimator)) +
        geom_line(linewidth = 1) +
        scale_colour_manual(values=c("#F27127FF","#FDC718FF",
                                     "#038C33FF", "#207CF5FF"),
                            breaks=c("Diff-in-mean", "AIPW",
                                     "Hybrid", "Single-arm"))+
        geom_hline(yintercept = 0.8,
                   linetype="dashed",
                   color = "#084594",
                   linewidth = 1) +
        geom_vline(xintercept= 47,
                   linetype="dashed",
                   color = "#011C40FF",
                   linewidth = 1)+
        geom_vline(xintercept= 57,
                   linetype="dashed",
                   color =  "#005A32",
                   linewidth = 1)+
        geom_vline(xintercept= 251,
                   linetype="dashed",
                   color = "#FEB24C",
                   linewidth = 1)+
        theme(legend.position = "none")+
        labs(x=expression(Sample~Size~N[R]), y="Power")+
        scale_x_continuous(breaks=seq(0, 600, by=100))+
        scale_y_continuous(breaks=seq(0, 1, by = 0.2))+
        theme(title = element_text(size = 20),
              axis.title=element_text(size = 20),
              legend.title=element_text(size = 20),
              axis.text=element_text(size = 20),
              legend.text=element_text(size = 20))+
        labs(title = expression(Heterogeneous~treatment~effect~","~pi[A] ==0.8))
      
      p_sensi1_power_08
    
      
    }
    
    # Power  -  09
    {
      
      Res <- readRDS("res_sensi1/res_09_power_sensi1.RData")
      
      Res <- bind_rows(Res)
      Res$n <- seq(6, 600, by = 1) 
      
      Res_plot <-  Res %>% 
        pivot_longer(
          cols = c("power_EC", "power_AIPW", "power_naive"), 
          names_to = "Estimator",
          values_to = "power"
        )
      
      # add single-arm case
      Res_sa <- readRDS("res_sensi1/res_sa_power_sensi1.RData")
      
      Res_sa <- bind_rows(Res_sa)
      Res_sa$n <- seq(6, 600, by = 1) 
      
      Res_plot_sa <-  Res_sa %>% 
        pivot_longer(
          cols = c("power_SA"), 
          names_to = "Estimator",
          values_to = "power"
        )
      
      Res_plot <- rbind(Res_plot, Res_plot_sa)
      
      Res_plot$Estimator <- factor(Res_plot$Estimator, 
                                   levels = c(
                                     "power_EC",
                                     "power_AIPW",
                                     "power_naive",
                                     "power_SA"
                                   ))
      levels(Res_plot$Estimator) <- c(
        "Hybrid", 
        "AIPW",
        "Diff-in-mean",
        "Single-arm")
      
      
      p_sensi1_power_09 <- ggplot(data = Res_plot, aes(x = n, y = power, colour = Estimator)) +
        geom_line(linewidth = 1) +
        scale_colour_manual(values=c("#F27127FF","#FDC718FF",
                                     "#038C33FF", "#207CF5FF"),
                            breaks=c("Diff-in-mean", "AIPW",
                                     "Hybrid", "Single-arm"))+
        geom_hline(yintercept = 0.8,
                   linetype="dashed",
                   color = "#084594",
                   linewidth = 1) +
        geom_vline(xintercept= 47,
                   linetype="dashed",
                   color = "#011C40FF",
                   linewidth = 1)+
        geom_vline(xintercept= 51,
                   linetype="dashed",
                   color =  "#005A32",
                   linewidth = 1)+
        geom_vline(xintercept= 441,
                   linetype="dashed",
                   color = "#FEB24C",
                   linewidth = 1)+
        # theme(legend.position = "none")+
        labs(x=expression(Sample~Size~N[R]), y="Power")+
        scale_x_continuous(breaks=seq(0, 600, by=100))+
        scale_y_continuous(breaks=seq(0, 1, by = 0.2))+
        theme(title = element_text(size = 20),
              axis.title=element_text(size = 20),
              legend.title=element_text(size = 20),
              axis.text=element_text(size = 20),
              legend.text=element_text(size = 20))+
        labs(title = expression(Heterogeneous~treatment~effect~","~pi[A] ==0.9))
      
      p_sensi1_power_09
      
    }
    
  }
  
  ## Type I error_sufficient
  {
    
    # Type I error - 05
    {
      
      Res <- readRDS("res_sensi1/res_05_TypeI_sensi1.RData")
      
      Res <- bind_rows(Res)
      Res$n <- seq(6, 600, by = 1) 
      
      Res_plot <-  Res %>% 
        pivot_longer(
          cols = c("power_EC", "power_AIPW", "power_naive"), 
          names_to = "Estimator",
          values_to = "Type_I_Error"
        )
      
      # add single-arm case
      Res_sa <- readRDS("res_sensi1/res_sa_TypeI_sensi1.RData")
      
      Res_sa <- bind_rows(Res_sa)
      Res_sa$n <- seq(6, 600, by = 1) 
      
      Res_plot_sa <-  Res_sa %>% 
        pivot_longer(
          cols = c("power_SA"), 
          names_to = "Estimator",
          values_to = "Type_I_Error"
        )
      
      Res_plot <- rbind(Res_plot, Res_plot_sa)
      
      Res_plot$Estimator <- factor(Res_plot$Estimator, 
                                   levels = c(
                                     "power_EC",
                                     "power_AIPW",
                                     "power_naive",
                                     "power_SA"
                                   ))
      levels(Res_plot$Estimator) <- c(
        "Hybrid", 
        "AIPW",
        "Diff-in-mean",
        "Single-arm")
      
      
      p_sensi1_TypeI_05 <- ggplot(data = Res_plot, aes(x = n, y = Type_I_Error, colour = Estimator)) +
        geom_line(linewidth = 1) +
        scale_colour_manual(values=c("#F27127FF","#FDC718FF",
                                     "#038C33FF", "#207CF5FF"), 
                            breaks=c("Diff-in-mean", "AIPW",
                                     "Hybrid", "Single-arm"))+
        geom_hline(yintercept = 0.05,
                   linetype="dashed",
                   color = "#084594",
                   linewidth = 1) +
        geom_vline(xintercept= 47,
                   linetype="dashed",
                   color = "#011C40FF",
                   linewidth = 1)+
        geom_vline(xintercept= 88,
                   linetype="dashed",
                   color =  "#005A32",
                   linewidth = 1)+
        geom_vline(xintercept= 162,
                   linetype="dashed",
                   color = "#FEB24C",
                   linewidth = 1)+
        theme(legend.position = "none")+
        labs(x=expression(Sample~Size~N[R]), y="Type I Error")+
        scale_x_continuous(breaks=seq(0, 600, by=100))+
        scale_y_continuous(breaks = seq(0, 0.546, by=0.1))+
        expand_limits(y=c(0, 0.546))+
        theme(title = element_text(size = 20),
              axis.title=element_text(size = 20),
              legend.title=element_text(size = 20),
              axis.text=element_text(size = 20),
              legend.text=element_text(size = 20))+
        labs(title = expression(Heterogeneous~treatment~effect~","~pi[A] ==0.5))
      
      p_sensi1_TypeI_05
        
    }
    
    # Type I error - 06
    {
      
      Res <- readRDS("res_sensi1/res_06_TypeI_sensi1.RData")
      
      Res <- bind_rows(Res)
      Res$n <- seq(6, 600, by = 1) 
      
      Res_plot <-  Res %>% 
        pivot_longer(
          cols = c("power_EC", "power_AIPW", "power_naive"), 
          names_to = "Estimator",
          values_to = "Type_I_Error"
        )
      
      # add single-arm case
      Res_sa <- readRDS("res_sensi1/res_sa_TypeI_sensi1.RData")
      
      Res_sa <- bind_rows(Res_sa)
      Res_sa$n <- seq(6, 600, by = 1) 
      
      Res_plot_sa <-  Res_sa %>% 
        pivot_longer(
          cols = c("power_SA"), 
          names_to = "Estimator",
          values_to = "Type_I_Error"
        )
      
      Res_plot <- rbind(Res_plot, Res_plot_sa)
      
      Res_plot$Estimator <- factor(Res_plot$Estimator, 
                                   levels = c(
                                     "power_EC",
                                     "power_AIPW",
                                     "power_naive",
                                     "power_SA"
                                   ))
      levels(Res_plot$Estimator) <- c(
        "Hybrid", 
        "AIPW",
        "Diff-in-mean",
        "Single-arm")
      
      
      p_sensi1_TypeI_06 <- ggplot(data = Res_plot, aes(x = n, y = Type_I_Error, colour = Estimator)) +
        geom_line(linewidth = 1) +
        scale_colour_manual(values=c("#F27127FF","#FDC718FF",
                                     "#038C33FF", "#207CF5FF"), 
                            breaks=c("Diff-in-mean", "AIPW",
                                     "Hybrid", "Single-arm"))+
        geom_hline(yintercept = 0.05,
                   linetype="dashed",
                   color = "#084594",
                   linewidth = 1) +
        geom_vline(xintercept= 47,
                   linetype="dashed",
                   color = "#011C40FF",
                   linewidth = 1)+
        geom_vline(xintercept= 74,
                   linetype="dashed",
                   color =  "#005A32",
                   linewidth = 1)+
        geom_vline(xintercept= 169,
                   linetype="dashed",
                   color = "#FEB24C",
                   linewidth = 1)+
        theme(legend.position = "none")+
        labs(x=expression(Sample~Size~N[R]), y="Type I Error")+
        scale_x_continuous(breaks=seq(0, 600, by=100))+
        scale_y_continuous(breaks = seq(0, 0.546, by=0.1))+
        expand_limits(y=c(0, 0.546))+
        theme(title = element_text(size = 20),
              axis.title=element_text(size = 20),
              legend.title=element_text(size = 20),
              axis.text=element_text(size = 20),
              legend.text=element_text(size = 20))+
        labs(title = expression(Heterogeneous~treatment~effect~","~pi[A] ==0.6))
      
      p_sensi1_TypeI_06
       
      
    }
    
    # Type I error - 07
    {
     
       Res <- readRDS("res_sensi1/res_07_TypeI_sensi1.RData")
      
      Res <- bind_rows(Res)
      Res$n <- seq(6, 600, by = 1) 
      
      Res_plot <-  Res %>% 
        pivot_longer(
          cols = c("power_EC", "power_AIPW", "power_naive"), 
          names_to = "Estimator",
          values_to = "Type_I_Error"
        )
      
      # Res_plot <- Res_plot[4:6]
      
      # add single-arm case
      Res_sa <- readRDS("res_sensi1/res_sa_TypeI_sensi1.RData")
      
      Res_sa <- bind_rows(Res_sa)
      Res_sa$n <- seq(6, 600, by = 1) 
      
      Res_plot_sa <-  Res_sa %>% 
        pivot_longer(
          cols = c("power_SA"), 
          names_to = "Estimator",
          values_to = "Type_I_Error"
        )
      
      Res_plot <- rbind(Res_plot, Res_plot_sa)
      
      Res_plot$Estimator <- factor(Res_plot$Estimator, 
                                   levels = c(
                                     "power_EC",
                                     "power_AIPW",
                                     "power_naive",
                                     "power_SA"
                                   ))
      levels(Res_plot$Estimator) <- c(
        "Hybrid", 
        "AIPW",
        "Diff-in-mean",
        "Single-arm")
      
      
      p_sensi1_TypeI_07 <- ggplot(data = Res_plot, aes(x = n, y = Type_I_Error, colour = Estimator)) +
        geom_line(linewidth = 1) +
        scale_colour_manual(values=c("#F27127FF","#FDC718FF",
                                     "#038C33FF", "#207CF5FF"), 
                            breaks=c("Diff-in-mean", "AIPW",
                                     "Hybrid", "Single-arm"))+
        geom_hline(yintercept = 0.05,
                   linetype="dashed",
                   color = "#084594",
                   linewidth = 1) +
        geom_vline(xintercept= 47,
                   linetype="dashed",
                   color = "#011C40FF",
                   linewidth = 1)+
        geom_vline(xintercept= 65,
                   linetype="dashed",
                   color =  "#005A32",
                   linewidth = 1)+
        geom_vline(xintercept= 192,
                   linetype="dashed",
                   color = "#FEB24C",
                   linewidth = 1)+
        theme(legend.position = "none")+
        labs(x=expression(Sample~Size~N[R]), y="Type I Error")+
        scale_x_continuous(breaks=seq(0, 600, by=100))+
        scale_y_continuous(breaks = seq(0, 0.546, by=0.1))+
        expand_limits(y=c(0, 0.546))+
        theme(title = element_text(size = 20),
              axis.title=element_text(size = 20),
              legend.title=element_text(size = 20),
              axis.text=element_text(size = 20),
              legend.text=element_text(size = 20))+
        labs(title = expression(Heterogeneous~treatment~effect~","~pi[A] ==0.7))
      
      p_sensi1_TypeI_07
      
    }
    
    # Type I error - 08
    {
      
      Res <- readRDS("res_sensi1/res_08_TypeI_sensi1.RData")
      
      Res <- bind_rows(Res)
      Res$n <- seq(6, 600, by = 1) 
      
      Res_plot <-  Res %>% 
        pivot_longer(
          cols = c("power_EC", "power_AIPW", "power_naive"), 
          names_to = "Estimator",
          values_to = "Type_I_Error"
        )
      
      # add single-arm case
      Res_sa <- readRDS("res_sensi1/res_sa_TypeI_sensi1.RData")
      
      Res_sa <- bind_rows(Res_sa)
      Res_sa$n <- seq(6, 600, by = 1) 
      
      Res_plot_sa <-  Res_sa %>% 
        pivot_longer(
          cols = c("power_SA"), 
          names_to = "Estimator",
          values_to = "Type_I_Error"
        )
      
      Res_plot <- rbind(Res_plot, Res_plot_sa)
      
      Res_plot$Estimator <- factor(Res_plot$Estimator, 
                                   levels = c(
                                     "power_EC",
                                     "power_AIPW",
                                     "power_naive",
                                     "power_SA"
                                   ))
      levels(Res_plot$Estimator) <- c(
        "Hybrid", 
        "AIPW",
        "Diff-in-mean",
        "Single-arm")
      
      
      p_sensi1_TypeI_08 <- ggplot(data = Res_plot, aes(x = n, y = Type_I_Error, colour = Estimator)) +
        geom_line(linewidth = 1) +
        scale_colour_manual(values=c("#F27127FF","#FDC718FF",
                                     "#038C33FF", "#207CF5FF"), 
                            breaks=c("Diff-in-mean", "AIPW",
                                     "Hybrid", "Single-arm"))+
        geom_hline(yintercept = 0.05,
                   linetype="dashed",
                   color = "#084594",
                   linewidth = 1) +
        geom_vline(xintercept= 47,
                   linetype="dashed",
                   color = "#011C40FF",
                   linewidth = 1)+
        geom_vline(xintercept= 57,
                   linetype="dashed",
                   color =  "#005A32",
                   linewidth = 1)+
        geom_vline(xintercept= 251,
                   linetype="dashed",
                   color = "#FEB24C",
                   linewidth = 1)+
        theme(legend.position = "none")+
        labs(x=expression(Sample~Size~N[R]), y="Type I Error")+
        scale_x_continuous(breaks=seq(0, 600, by=100))+
        scale_y_continuous(breaks = seq(0, 0.546, by=0.1))+
        expand_limits(y=c(0, 0.546))+
        theme(title = element_text(size = 20),
              axis.title=element_text(size = 20),
              legend.title=element_text(size = 20),
              axis.text=element_text(size = 20),
              legend.text=element_text(size = 20))+
        labs(title = expression(Heterogeneous~treatment~effect~","~pi[A] ==0.8))
      
      p_sensi1_TypeI_08
      
    }
    
    # Type I error - 09
    {
      
      Res <- readRDS("res_sensi1/res_09_TypeI_sensi1.RData")
      
      Res <- bind_rows(Res)
      Res$n <- seq(6, 600, by = 1) 
      
      Res_plot <-  Res %>% 
        pivot_longer(
          cols = c("power_EC", "power_AIPW", "power_naive"), 
          names_to = "Estimator",
          values_to = "Type_I_Error"
        )
      
      # add single-arm case
      Res_sa <- readRDS("res_sensi1/res_sa_TypeI_sensi1.RData")
      
      Res_sa <- bind_rows(Res_sa)
      Res_sa$n <- seq(6, 600, by = 1) 
      
      Res_plot_sa <-  Res_sa %>% 
        pivot_longer(
          cols = c("power_SA"), 
          names_to = "Estimator",
          values_to = "Type_I_Error"
        )
      
      Res_plot <- rbind(Res_plot, Res_plot_sa)
      
      Res_plot$Estimator <- factor(Res_plot$Estimator, 
                                   levels = c(
                                     "power_EC",
                                     "power_AIPW",
                                     "power_naive",
                                     "power_SA"
                                   ))
      levels(Res_plot$Estimator) <- c(
        "Hybrid", 
        "AIPW",
        "Diff-in-mean",
        "Single-arm")
      
      
      p_sensi1_TypeI_09 <- ggplot(data = Res_plot, aes(x = n, y = Type_I_Error, colour = Estimator)) +
        geom_line(linewidth = 1) +
        scale_colour_manual(values=c("#F27127FF","#FDC718FF",
                                     "#038C33FF", "#207CF5FF"), 
                            breaks=c("Diff-in-mean", "AIPW",
                                     "Hybrid", "Single-arm"))+
        geom_hline(yintercept = 0.05,
                   linetype="dashed",
                   color = "#084594",
                   linewidth = 1) +
        geom_vline(xintercept= 47,
                   linetype="dashed",
                   color = "#011C40FF",
                   linewidth = 1)+
        geom_vline(xintercept= 51,
                   linetype="dashed",
                   color =  "#005A32",
                   linewidth = 1)+
        geom_vline(xintercept= 441,
                   linetype="dashed",
                   color = "#FEB24C",
                   linewidth = 1)+
        # theme(legend.position = "none")+
        labs(x=expression(Sample~Size~N[R]), y="Type I Error")+
        scale_x_continuous(breaks=seq(0, 600, by=100))+
        scale_y_continuous(breaks = seq(0, 0.546, by=0.1))+
        expand_limits(y=c(0, 0.546))+
        theme(title = element_text(size = 20),
              axis.title=element_text(size = 20),
              legend.title=element_text(size = 20),
              axis.text=element_text(size = 20),
              legend.text=element_text(size = 20))+
        labs(title = expression(Heterogeneous~treatment~effect~","~pi[A] ==0.9))
      
      p_sensi1_TypeI_09
      
    }
    
  }
  
}

### Sensitivity analysis 2 (N_E = 1000)
{
  
  ## Power_sufficient
  {
    
    # Power  -  05
    {
      
      Res <- readRDS("res_sensi2/res_05_power_sensi2.RData")
      
      Res <- bind_rows(Res)
      Res$n <- seq(6, 600, by = 1) 
      
      Res_plot <-  Res %>% 
        pivot_longer(
          cols = c("power_EC", "power_AIPW", "power_naive"), 
          names_to = "Estimator",
          values_to = "power"
        )
      
      # add single-arm case
      Res_sa <- readRDS("res_sensi2/res_sa_power_sensi2.RData")
      
      Res_sa <- bind_rows(Res_sa)
      Res_sa$n <- seq(6, 600, by = 1) 
      
      Res_plot_sa <-  Res_sa %>% 
        pivot_longer(
          cols = c("power_SA"), 
          names_to = "Estimator",
          values_to = "power"
        )
      
      Res_plot <- rbind(Res_plot, Res_plot_sa)
      
      Res_plot$Estimator <- factor(Res_plot$Estimator, 
                                   levels = c(
                                     "power_EC",
                                     "power_AIPW",
                                     "power_naive",
                                     "power_SA"
                                   ))
      levels(Res_plot$Estimator) <- c(
        "Hybrid", 
        "AIPW",
        "Diff-in-mean",
        "Single-arm")
      
      p_sensi2_power_05 <- ggplot(data = Res_plot, aes(x = n, y = power, colour = Estimator)) +
        geom_line(linewidth = 1) +
        scale_colour_manual(values=c("#F27127FF","#FDC718FF",
                                     "#038C33FF", "#207CF5FF"),
                            breaks=c("Diff-in-mean", "AIPW",
                                     "Hybrid", "Single-arm"))+
        geom_hline(yintercept = 0.8,
                   linetype="dashed",
                   color = "#084594",
                   linewidth = 1) +
        geom_vline(xintercept= 55,
                   linetype="dashed",
                   color = "#011C40FF",
                   linewidth = 1)+
        geom_vline(xintercept= 109,
                   linetype="dashed",
                   color =  "#005A32",
                   linewidth = 1) +
        geom_vline(xintercept= 201,
                   linetype="dashed",
                   color = "#FEB24C",
                   linewidth = 1)+
        theme(legend.position = "none")+
        labs(x=expression(Sample~Size~N[R]), y="Power")+
        scale_x_continuous(breaks=seq(0, 600, by=100))+
        scale_y_continuous(breaks=seq(0, 1, by = 0.2))+
        theme(title = element_text(size = 20),
              axis.title=element_text(size = 20),
              legend.title=element_text(size = 20),
              axis.text=element_text(size = 20),
              legend.text=element_text(size = 20))+
        labs(title = expression("Non-constant"~conditional~var~","~pi[A] ==0.5))
      
      p_sensi2_power_05
      
    }
    
    # Power  -  06
    {
      
      Res <- readRDS("res_sensi2/res_06_power_sensi2.RData")
      
      Res <- bind_rows(Res) 
      Res$n <- seq(6, 600, by = 1) 
      
      Res_plot <-  Res %>% 
        pivot_longer(
          cols = c("power_EC", "power_AIPW", "power_naive"), 
          names_to = "Estimator",
          values_to = "power"
        )
      
      # add single-arm case
      Res_sa <- readRDS("res_sensi2/res_sa_power_sensi2.RData")
      
      Res_sa <- bind_rows(Res_sa)
      Res_sa$n <- seq(6, 600, by = 1) 
      
      Res_plot_sa <-  Res_sa %>% 
        pivot_longer(
          cols = c("power_SA"), 
          names_to = "Estimator",
          values_to = "power"
        )
      
      Res_plot <- rbind(Res_plot, Res_plot_sa)
      
      Res_plot$Estimator <- factor(Res_plot$Estimator, 
                                   levels = c(
                                     "power_EC",
                                     "power_AIPW",
                                     "power_naive",
                                     "power_SA"
                                   ))
      levels(Res_plot$Estimator) <- c(
        "Hybrid", 
        "AIPW",
        "Diff-in-mean",
        "Single-arm")
      
      p_sensi2_power_06 <- ggplot(data = Res_plot, aes(x = n, y = power, colour = Estimator)) +
        geom_line(linewidth = 1) +
        scale_colour_manual(values=c("#F27127FF","#FDC718FF",
                                     "#038C33FF", "#207CF5FF"),
                            breaks=c("Diff-in-mean", "AIPW",
                                     "Hybrid", "Single-arm"))+
        geom_hline(yintercept = 0.8,
                   linetype="dashed",
                   color = "#084594",
                   linewidth = 1) +
        geom_vline(xintercept= 55,
                   linetype="dashed",
                   color = "#011C40FF",
                   linewidth = 1)+
        geom_vline(xintercept= 91,
                   linetype="dashed",
                   color =  "#005A32",
                   linewidth = 1)+
        geom_vline(xintercept= 210,
                   linetype="dashed",
                   color = "#FEB24C",
                   linewidth = 1)+
        theme(legend.position = "none")+
        labs(x=expression(Sample~Size~N[R]), y="Power")+
        scale_x_continuous(breaks=seq(0, 600, by=100))+
        scale_y_continuous(breaks=seq(0, 1, by = 0.2))+
        theme(title = element_text(size = 20),
              axis.title=element_text(size = 20),
              legend.title=element_text(size = 20),
              axis.text=element_text(size = 20),
              legend.text=element_text(size = 20))+
        labs(title = expression("Non-constant"~conditional~var~","~pi[A] ==0.6))
      
      p_sensi2_power_06
      
    }
    
    # Power  -  07
    {
      
      Res <- readRDS("res_sensi2/res_07_power_sensi2.RData")
      
      Res <- bind_rows(Res)
      Res$n <- seq(6, 600, by = 1) 
      
      Res_plot <-  Res %>% 
        pivot_longer(
          cols = c("power_EC", "power_AIPW", "power_naive"), 
          names_to = "Estimator",
          values_to = "power"
        )
      
      # add single-arm case
      Res_sa <- readRDS("res_sensi2/res_sa_power_sensi2.RData")
      
      Res_sa <- bind_rows(Res_sa)
      Res_sa$n <- seq(6, 600, by = 1) 
      
      Res_plot_sa <-  Res_sa %>% 
        pivot_longer(
          cols = c("power_SA"), 
          names_to = "Estimator",
          values_to = "power"
        )
      
      Res_plot <- rbind(Res_plot, Res_plot_sa)
      
      Res_plot$Estimator <- factor(Res_plot$Estimator, 
                                   levels = c(
                                     "power_EC",
                                     "power_AIPW",
                                     "power_naive",
                                     "power_SA"
                                   ))
      levels(Res_plot$Estimator) <- c(
        "Hybrid", 
        "AIPW",
        "Diff-in-mean",
        "Single-arm")
      
      
      p_sensi2_power_07 <- ggplot(data = Res_plot, aes(x = n, y = power, colour = Estimator)) +
        geom_line(linewidth = 1) +
        scale_colour_manual(values=c("#F27127FF","#FDC718FF",
                                     "#038C33FF", "#207CF5FF"),
                            breaks=c("Diff-in-mean", "AIPW",
                                     "Hybrid", "Single-arm"))+
        geom_hline(yintercept = 0.8,
                   linetype="dashed",
                   color = "#084594",
                   linewidth = 1) +
        geom_vline(xintercept= 55,
                   linetype="dashed",
                   color = "#011C40FF",
                   linewidth = 1)+
        geom_vline(xintercept= 78,
                   linetype="dashed",
                   color =  "#005A32",
                   linewidth = 1)+
        geom_vline(xintercept= 240,
                   linetype="dashed",
                   color = "#FEB24C",
                   linewidth = 1)+
        theme(legend.position = "none")+
        labs(x=expression(Sample~Size~N[R]), y="Power")+
        scale_x_continuous(breaks=seq(0, 600, by=100))+
        scale_y_continuous(breaks=seq(0, 1, by = 0.2))+
        theme(title = element_text(size = 20),
              axis.title=element_text(size = 20),
              legend.title=element_text(size = 20),
              axis.text=element_text(size = 20),
              legend.text=element_text(size = 20))+
        labs(title = expression("Non-constant"~conditional~var~","~pi[A] ==0.7))
      
      p_sensi2_power_07
       
      
    }
    
    # Power  -  08
    {
      
      Res <- readRDS("res_sensi2/res_08_power_sensi2.RData")
      
      Res <- bind_rows(Res) 
      Res$n <- seq(6, 600, by = 1) 
      
      Res_plot <-  Res %>% 
        pivot_longer(
          cols = c("power_EC", "power_AIPW", "power_naive"), 
          names_to = "Estimator",
          values_to = "power"
        )
      
      # add single-arm case
      Res_sa <- readRDS("res_sensi2/res_sa_power_sensi2.RData")
      
      Res_sa <- bind_rows(Res_sa)
      Res_sa$n <- seq(6, 600, by = 1) 
      
      Res_plot_sa <-  Res_sa %>% 
        pivot_longer(
          cols = c("power_SA"), 
          names_to = "Estimator",
          values_to = "power"
        )
      
      Res_plot <- rbind(Res_plot, Res_plot_sa)
      
      Res_plot$Estimator <- factor(Res_plot$Estimator, 
                                   levels = c(
                                     "power_EC",
                                     "power_AIPW",
                                     "power_naive",
                                     "power_SA"
                                   ))
      levels(Res_plot$Estimator) <- c(
        "Hybrid", 
        "AIPW",
        "Diff-in-mean",
        "Single-arm")
      
      p_sensi2_power_08 <- ggplot(data = Res_plot, aes(x = n, y = power, colour = Estimator)) +
        geom_line(linewidth = 1) +
        scale_colour_manual(values=c("#F27127FF","#FDC718FF",
                                     "#038C33FF", "#207CF5FF"),
                            breaks=c("Diff-in-mean", "AIPW",
                                     "Hybrid", "Single-arm"))+
        geom_hline(yintercept = 0.8,
                   linetype="dashed",
                   color = "#084594",
                   linewidth = 1) +
        geom_vline(xintercept= 55,
                   linetype="dashed",
                   color = "#011C40FF",
                   linewidth = 1)+
        geom_vline(xintercept= 69,
                   linetype="dashed",
                   color =  "#005A32",
                   linewidth = 1)+
        geom_vline(xintercept= 314,
                   linetype="dashed",
                   color = "#FEB24C",
                   linewidth = 1)+
        theme(legend.position = "none")+
        labs(x=expression(Sample~Size~N[R]), y="Power")+
        scale_x_continuous(breaks=seq(0, 600, by=100))+
        scale_y_continuous(breaks=seq(0, 1, by = 0.2))+
        theme(title = element_text(size = 20),
              axis.title=element_text(size = 20),
              legend.title=element_text(size = 20),
              axis.text=element_text(size = 20),
              legend.text=element_text(size = 20))+
        labs(title = expression("Non-constant"~conditional~var~","~pi[A] ==0.8))
      
      p_sensi2_power_08
      
    }
    
    # Power  -  09
    {
      
      Res <- readRDS("res_sensi2/res_09_power_sensi2.RData")
      
      Res <- bind_rows(Res)
      Res$n <- seq(6, 600, by = 1) 
      
      Res_plot <-  Res %>% 
        pivot_longer(
          cols = c("power_EC", "power_AIPW", "power_naive"), 
          names_to = "Estimator",
          values_to = "power"
        )
      
      # add single-arm case
      Res_sa <- readRDS("res_sensi2/res_sa_power_sensi2.RData")
      
      Res_sa <- bind_rows(Res_sa)
      Res_sa$n <- seq(6, 600, by = 1) 
      
      Res_plot_sa <-  Res_sa %>% 
        pivot_longer(
          cols = c("power_SA"), 
          names_to = "Estimator",
          values_to = "power"
        )
      
      Res_plot <- rbind(Res_plot, Res_plot_sa)
      
      Res_plot$Estimator <- factor(Res_plot$Estimator, 
                                   levels = c(
                                     "power_EC",
                                     "power_AIPW",
                                     "power_naive",
                                     "power_SA"
                                   ))
      levels(Res_plot$Estimator) <- c(
        "Hybrid", 
        "AIPW",
        "Diff-in-mean",
        "Single-arm")
      
      
      p_sensi2_power_09 <- ggplot(data = Res_plot, aes(x = n, y = power, colour = Estimator)) +
        geom_line(linewidth = 1) +
        scale_colour_manual(values=c("#F27127FF","#FDC718FF",
                                     "#038C33FF", "#207CF5FF"),
                            breaks=c("Diff-in-mean", "AIPW",
                                     "Hybrid", "Single-arm"))+
        geom_hline(yintercept = 0.8,
                   linetype="dashed",
                   color = "#084594",
                   linewidth = 1) +
        geom_vline(xintercept= 55,
                   linetype="dashed",
                   color = "#011C40FF",
                   linewidth = 1)+
        geom_vline(xintercept= 61,
                   linetype="dashed",
                   color =  "#005A32",
                   linewidth = 1)+
        geom_vline(xintercept= 557,
                   linetype="dashed",
                   color = "#FEB24C",
                   linewidth = 1)+
        # theme(legend.position = "none")+
        labs(x=expression(Sample~Size~N[R]), y="Power")+
        scale_x_continuous(breaks=seq(0, 600, by=100))+
        scale_y_continuous(breaks=seq(0, 1, by = 0.2))+
        theme(title = element_text(size = 20),
              axis.title=element_text(size = 20),
              legend.title=element_text(size = 20),
              axis.text=element_text(size = 20),
              legend.text=element_text(size = 20))+
        labs(title = expression("Non-constant"~conditional~var~","~pi[A] ==0.9))
      
      p_sensi2_power_09
      
    }
    
  }
  
  ## Type I error_sufficient
  {
    
    # Type I error - 05
    {
      
      Res <- readRDS("res_sensi2/res_05_TypeI_sensi2.RData")
      
      Res <- bind_rows(Res)
      Res$n <- seq(6, 600, by = 1) 
      
      Res_plot <-  Res %>% 
        pivot_longer(
          cols = c("power_EC", "power_AIPW", "power_naive"), 
          names_to = "Estimator",
          values_to = "Type_I_Error"
        )
      
      # add single-arm case
      Res_sa <- readRDS("res_sensi2/res_sa_TypeI_sensi2.RData")
      
      Res_sa <- bind_rows(Res_sa)
      Res_sa$n <- seq(6, 600, by = 1) 
      
      Res_plot_sa <-  Res_sa %>% 
        pivot_longer(
          cols = c("power_SA"), 
          names_to = "Estimator",
          values_to = "Type_I_Error"
        )
      
      Res_plot <- rbind(Res_plot, Res_plot_sa)
      
      Res_plot$Estimator <- factor(Res_plot$Estimator, 
                                   levels = c(
                                     "power_EC",
                                     "power_AIPW",
                                     "power_naive",
                                     "power_SA"
                                   ))
      levels(Res_plot$Estimator) <- c(
        "Hybrid", 
        "AIPW",
        "Diff-in-mean",
        "Single-arm")
      
      
      p_sensi2_TypeI_05 <- ggplot(data = Res_plot, aes(x = n, y = Type_I_Error, colour = Estimator)) +
        geom_line(linewidth = 1) +
        scale_colour_manual(values=c("#F27127FF","#FDC718FF",
                                     "#038C33FF", "#207CF5FF"), 
                            breaks=c("Diff-in-mean", "AIPW",
                                     "Hybrid", "Single-arm"))+
        geom_hline(yintercept = 0.05,
                   linetype="dashed",
                   color = "#084594",
                   linewidth = 1) +
        geom_vline(xintercept= 55,
                   linetype="dashed",
                   color = "#011C40FF",
                   linewidth = 1)+
        geom_vline(xintercept= 109,
                   linetype="dashed",
                   color =  "#005A32",
                   linewidth = 1)+
        geom_vline(xintercept= 201,
                   linetype="dashed",
                   color = "#FEB24C",
                   linewidth = 1)+
        theme(legend.position = "none")+
        # labs(x=expression(RCT~Sample~Size~N[R]), y="Type I Error")+
        labs(x=expression(Sample~Size~N[R]), y="Type I Error")+
        scale_x_continuous(breaks=seq(0, 600, by=100))+
        scale_y_continuous(breaks = seq(0, 0.546, by=0.1))+
        expand_limits(y=c(0, 0.546))+
        theme(title = element_text(size = 20),
              axis.title=element_text(size = 20),
              legend.title=element_text(size = 20),
              axis.text=element_text(size = 20),
              legend.text=element_text(size = 20))+
        labs(title = expression("Non-constant"~conditional~var~","~pi[A] ==0.5))
      
      p_sensi2_TypeI_05
       
    }
    
    # Type I error - 06
    {
      
      Res <- readRDS("res_sensi2/res_06_TypeI_sensi2.RData")
      
      Res <- bind_rows(Res)
      Res$n <- seq(6, 600, by = 1) 
      
      Res_plot <-  Res %>% 
        pivot_longer(
          cols = c("power_EC", "power_AIPW", "power_naive"), 
          names_to = "Estimator",
          values_to = "Type_I_Error"
        )
      
      
      # add single-arm case
      Res_sa <- readRDS("res_sensi2/res_sa_TypeI_sensi2.RData")
      
      Res_sa <- bind_rows(Res_sa)
      Res_sa$n <- seq(6, 600, by = 1) 
      
      Res_plot_sa <-  Res_sa %>% 
        pivot_longer(
          cols = c("power_SA"), 
          names_to = "Estimator",
          values_to = "Type_I_Error"
        )
      
      Res_plot <- rbind(Res_plot, Res_plot_sa)
      
      Res_plot$Estimator <- factor(Res_plot$Estimator, 
                                   levels = c(
                                     "power_EC",
                                     "power_AIPW",
                                     "power_naive",
                                     "power_SA"
                                   ))
      levels(Res_plot$Estimator) <- c(
        "Hybrid", 
        "AIPW",
        "Diff-in-mean",
        "Single-arm")
      
      
      p_sensi2_TypeI_06 <- ggplot(data = Res_plot, aes(x = n, y = Type_I_Error, colour = Estimator)) +
        geom_line(linewidth = 1) +
        scale_colour_manual(values=c("#F27127FF","#FDC718FF",
                                     "#038C33FF", "#207CF5FF"), 
                            breaks=c("Diff-in-mean", "AIPW",
                                     "Hybrid", "Single-arm"))+
        geom_hline(yintercept = 0.05,
                   linetype="dashed",
                   color = "#084594",
                   linewidth = 1) +
        geom_vline(xintercept= 55,
                   linetype="dashed",
                   color = "#011C40FF",
                   linewidth = 1)+
        geom_vline(xintercept= 91,
                   linetype="dashed",
                   color =  "#005A32",
                   linewidth = 1)+
        geom_vline(xintercept= 210,
                   linetype="dashed",
                   color = "#FEB24C",
                   linewidth = 1)+
        theme(legend.position = "none")+
        labs(x=expression(Sample~Size~N[R]), y="Type I Error")+
        scale_x_continuous(breaks=seq(0, 600, by=100))+
        scale_y_continuous(breaks = seq(0, 0.546, by=0.1))+
        expand_limits(y=c(0, 0.546))+
        theme(title = element_text(size = 20),
              axis.title=element_text(size = 20),
              legend.title=element_text(size = 20),
              axis.text=element_text(size = 20),
              legend.text=element_text(size = 20))+
        labs(title = expression("Non-constant"~conditional~var~","~pi[A] ==0.6))
      
      p_sensi2_TypeI_06
      
      
    }
    
    # Type I error - 07
    {
      
      Res <- readRDS("res_sensi2/res_07_TypeI_sensi2.RData")
      
      Res <- bind_rows(Res)
      Res$n <- seq(6, 600, by = 1) 
      
      Res_plot <-  Res %>% 
        pivot_longer(
          cols = c("power_EC", "power_AIPW", "power_naive"), 
          names_to = "Estimator",
          values_to = "Type_I_Error"
        )
      
      # add single-arm case
      Res_sa <- readRDS("res_sensi2/res_sa_TypeI_sensi2.RData")
      
      Res_sa <- bind_rows(Res_sa)
      Res_sa$n <- seq(6, 600, by = 1) 
      
      Res_plot_sa <-  Res_sa %>% 
        pivot_longer(
          cols = c("power_SA"), 
          names_to = "Estimator",
          values_to = "Type_I_Error"
        )
      
      Res_plot <- rbind(Res_plot, Res_plot_sa)
      
      Res_plot$Estimator <- factor(Res_plot$Estimator, 
                                   levels = c(
                                     "power_EC",
                                     "power_AIPW",
                                     "power_naive",
                                     "power_SA"
                                   ))
      levels(Res_plot$Estimator) <- c(
        "Hybrid", 
        "AIPW",
        "Diff-in-mean",
        "Single-arm")
      
      
      p_sensi2_TypeI_07 <- ggplot(data = Res_plot, aes(x = n, y = Type_I_Error, colour = Estimator)) +
        geom_line(linewidth = 1) +
        scale_colour_manual(values=c("#F27127FF","#FDC718FF",
                                     "#038C33FF", "#207CF5FF"), 
                            breaks=c("Diff-in-mean", "AIPW",
                                     "Hybrid", "Single-arm"))+
        geom_hline(yintercept = 0.05,
                   linetype="dashed",
                   color = "#084594",
                   linewidth = 1) +
        geom_vline(xintercept= 55,
                   linetype="dashed",
                   color = "#011C40FF",
                   linewidth = 1)+
        geom_vline(xintercept= 78,
                   linetype="dashed",
                   color =  "#005A32",
                   linewidth = 1)+
        geom_vline(xintercept= 240,
                   linetype="dashed",
                   color = "#FEB24C",
                   linewidth = 1)+
        theme(legend.position = "none")+
        labs(x=expression(Sample~Size~N[R]), y="Type I Error")+
        scale_x_continuous(breaks=seq(0, 600, by=100))+
        scale_y_continuous(breaks = seq(0, 0.546, by=0.1))+
        expand_limits(y=c(0, 0.546))+
        theme(title = element_text(size = 20),
              axis.title=element_text(size = 20),
              legend.title=element_text(size = 20),
              axis.text=element_text(size = 20),
              legend.text=element_text(size = 20))+
        labs(title = expression("Non-constant"~conditional~var~","~pi[A] ==0.7))
      
      p_sensi2_TypeI_07
      
    }
    
    # Type I error - 08
    {
      
      Res <- readRDS("res_sensi2/res_08_TypeI_sensi2.RData")
      
      Res <- bind_rows(Res)
      Res$n <- seq(6, 600, by = 1) 
      
      Res_plot <-  Res %>% 
        pivot_longer(
          cols = c("power_EC", "power_AIPW", "power_naive"), 
          names_to = "Estimator",
          values_to = "Type_I_Error"
        )
      
      # add single-arm case
      Res_sa <- readRDS("res_sensi2/res_sa_TypeI_sensi2.RData")
      
      Res_sa <- bind_rows(Res_sa)
      Res_sa$n <- seq(6, 600, by = 1) 
      
      Res_plot_sa <-  Res_sa %>% 
        pivot_longer(
          cols = c("power_SA"), 
          names_to = "Estimator",
          values_to = "Type_I_Error"
        )
      
      Res_plot <- rbind(Res_plot, Res_plot_sa)
      
      Res_plot$Estimator <- factor(Res_plot$Estimator, 
                                   levels = c(
                                     "power_EC",
                                     "power_AIPW",
                                     "power_naive",
                                     "power_SA"
                                   ))
      levels(Res_plot$Estimator) <- c(
        "Hybrid", 
        "AIPW",
        "Diff-in-mean",
        "Single-arm")
      
      
      p_sensi2_TypeI_08 <- ggplot(data = Res_plot, aes(x = n, y = Type_I_Error, colour = Estimator)) +
        geom_line(linewidth = 1) +
        scale_colour_manual(values=c("#F27127FF","#FDC718FF",
                                     "#038C33FF", "#207CF5FF"), 
                            breaks=c("Diff-in-mean", "AIPW",
                                     "Hybrid", "Single-arm"))+
        geom_hline(yintercept = 0.05,
                   linetype="dashed",
                   color = "#084594",
                   linewidth = 1) +
        geom_vline(xintercept= 55,
                   linetype="dashed",
                   color = "#011C40FF",
                   linewidth = 1)+
        geom_vline(xintercept= 69,
                   linetype="dashed",
                   color =  "#005A32",
                   linewidth = 1)+
        geom_vline(xintercept= 314,
                   linetype="dashed",
                   color = "#FEB24C",
                   linewidth = 1)+
        theme(legend.position = "none")+
        labs(x=expression(Sample~Size~N[R]), y="Type I Error")+
        scale_x_continuous(breaks=seq(0, 600, by=100))+
        scale_y_continuous(breaks = seq(0, 0.546, by=0.1))+
        expand_limits(y=c(0, 0.546))+
        theme(title = element_text(size = 20),
              axis.title=element_text(size = 20),
              legend.title=element_text(size = 20),
              axis.text=element_text(size = 20),
              legend.text=element_text(size = 20))+
        labs(title = expression("Non-constant"~conditional~var~","~pi[A] ==0.8))
      
      p_sensi2_TypeI_08
      
    }
    
    # Type I error - 09
    {
      
      Res <- readRDS("res_sensi2/res_09_TypeI_sensi2.RData")
      
      Res <- bind_rows(Res)
      Res$n <- seq(6, 600, by = 1) 
      
      Res_plot <-  Res %>% 
        pivot_longer(
          cols = c("power_EC", "power_AIPW", "power_naive"), 
          names_to = "Estimator",
          values_to = "Type_I_Error"
        )
      
      # add single-arm case
      Res_sa <- readRDS("res_sensi2/res_sa_TypeI_sensi2.RData")
      
      Res_sa <- bind_rows(Res_sa)
      Res_sa$n <- seq(6, 600, by = 1) 
      
      Res_plot_sa <-  Res_sa %>% 
        pivot_longer(
          cols = c("power_SA"), 
          names_to = "Estimator",
          values_to = "Type_I_Error"
        )
      
      Res_plot <- rbind(Res_plot, Res_plot_sa)
      
      Res_plot$Estimator <- factor(Res_plot$Estimator, 
                                   levels = c(
                                     "power_EC",
                                     "power_AIPW",
                                     "power_naive",
                                     "power_SA"
                                   ))
      levels(Res_plot$Estimator) <- c(
        "Hybrid", 
        "AIPW",
        "Diff-in-mean",
        "Single-arm")
      
      
      p_sensi2_TypeI_09 <- ggplot(data = Res_plot, aes(x = n, y = Type_I_Error, colour = Estimator)) +
        geom_line(linewidth = 1) +
        scale_colour_manual(values=c("#F27127FF","#FDC718FF",
                                     "#038C33FF", "#207CF5FF"), 
                            breaks=c("Diff-in-mean", "AIPW",
                                     "Hybrid", "Single-arm"))+
        geom_hline(yintercept = 0.05,
                   linetype="dashed",
                   color = "#084594",
                   linewidth = 1) +
        geom_vline(xintercept= 55,
                   linetype="dashed",
                   color = "#011C40FF",
                   linewidth = 1)+
        geom_vline(xintercept= 61,
                   linetype="dashed",
                   color =  "#005A32",
                   linewidth = 1)+
        geom_vline(xintercept= 557,
                   linetype="dashed",
                   color = "#FEB24C",
                   linewidth = 1)+
        # theme(legend.position = "none")+
        labs(x=expression(Sample~Size~N[R]), y="Type I Error")+
        scale_x_continuous(breaks=seq(0, 600, by=100))+
        scale_y_continuous(breaks = seq(0, 0.546, by=0.1))+
        expand_limits(y=c(0, 0.546))+
        theme(title = element_text(size = 20),
              axis.title=element_text(size = 20),
              legend.title=element_text(size = 20),
              axis.text=element_text(size = 20),
              legend.text=element_text(size = 20))+
        labs(title = expression("Non-constant"~conditional~var~","~pi[A] ==0.9))
      
      p_sensi2_TypeI_09
      
    }
    
  }
  
}

### Sensitivity analysis 3 (N_E = 1000)
{
  
  delta_U = 0.5
  
  ## Power_sufficient
  {
    
    # Power  -  05
    {
      
      Res <- readRDS(paste0("res_sensi3/res_05_power_sensi3_", delta_U,".RData"))
      
      Res <- bind_rows(Res)
      Res$n <- seq(10, 600, by = 5)
      
      Res_plot <-  Res %>% 
        pivot_longer(
          cols = c("power_EC", "power_AIPW", "power_naive"), 
          names_to = "Estimator",
          values_to = "power"
        )
      
      # add single-arm case
      Res_sa <- readRDS(paste0("res_sensi3/res_sa_power_sensi3_", delta_U,".RData"))
      
      Res_sa <- bind_rows(Res_sa)
      Res_sa$n <- seq(10, 600, by = 5)
      
      Res_plot_sa <-  Res_sa %>% 
        pivot_longer(
          cols = c("power_SA"), 
          names_to = "Estimator",
          values_to = "power"
        )
      
      Res_plot <- rbind(Res_plot, Res_plot_sa)
      
      Res_plot$Estimator <- factor(Res_plot$Estimator, 
                                   levels = c(
                                     "power_EC",
                                     "power_AIPW",
                                     "power_naive",
                                     "power_SA"
                                   ))
      levels(Res_plot$Estimator) <- c(
        "Hybrid", 
        "AIPW",
        "Diff-in-mean",
        "Single-arm")
      
      p_sensi3_power_05 <- ggplot(data = Res_plot, aes(x = n, y = power, colour = Estimator)) +
        geom_line(linewidth = 1) +
        scale_colour_manual(values=c("#F27127FF","#FDC718FF",
                                     "#038C33FF", "#207CF5FF"),
                            breaks=c("Diff-in-mean", "AIPW",
                                     "Hybrid", "Single-arm"))+
        geom_hline(yintercept = 0.8,
                   linetype="dashed",
                   color = "#084594",
                   linewidth = 1) +
        geom_vline(xintercept= 47,
                   linetype="dashed",
                   color = "#011C40FF",
                   linewidth = 1)+
        geom_vline(xintercept= 92,
                   linetype="dashed",
                   color =  "#005A32",
                   linewidth = 1) +
        geom_vline(xintercept= 175,
                   linetype="dashed",
                   color = "#FEB24C",
                   linewidth = 1)+
        theme(legend.position = "none")+
        labs(x=expression(Sample~Size~N[R]), y="Power")+
        scale_x_continuous(breaks=seq(0, 600, by=100))+
        scale_y_continuous(breaks=seq(0, 1, by = 0.2))+
        theme(title = element_text(size = 20),
              axis.title=element_text(size = 20),
              legend.title=element_text(size = 20),
              axis.text=element_text(size = 20),
              legend.text=element_text(size = 20))+
        labs(title = expression("Violation"~of~compatibility~","~pi[A] ==0.5))
      
      p_sensi3_power_05
      
    }
    
    # Power  -  06
    {
      
      Res <- readRDS(paste0("res_sensi3/res_06_power_sensi3_", delta_U,".RData"))
      
      Res <- bind_rows(Res)
      Res$n <- seq(10, 600, by = 5)
      
      Res_plot <-  Res %>% 
        pivot_longer(
          cols = c("power_EC", "power_AIPW", "power_naive"), 
          names_to = "Estimator",
          values_to = "power"
        )
      
      # add single-arm case
      Res_sa <- readRDS(paste0("res_sensi3/res_sa_power_sensi3_", delta_U,".RData"))
      
      Res_sa <- bind_rows(Res_sa)
      Res_sa$n <- seq(10, 600, by = 5)
      
      Res_plot_sa <-  Res_sa %>% 
        pivot_longer(
          cols = c("power_SA"), 
          names_to = "Estimator",
          values_to = "power"
        )
      
      Res_plot <- rbind(Res_plot, Res_plot_sa)
      
      Res_plot$Estimator <- factor(Res_plot$Estimator, 
                                   levels = c(
                                     "power_EC",
                                     "power_AIPW",
                                     "power_naive",
                                     "power_SA"
                                   ))
      levels(Res_plot$Estimator) <- c(
        "Hybrid", 
        "AIPW",
        "Diff-in-mean",
        "Single-arm")
      
      p_sensi3_power_06 <- ggplot(data = Res_plot, aes(x = n, y = power, colour = Estimator)) +
        geom_line(linewidth = 1) +
        scale_colour_manual(values=c("#F27127FF","#FDC718FF",
                                     "#038C33FF", "#207CF5FF"),
                            breaks=c("Diff-in-mean", "AIPW",
                                     "Hybrid", "Single-arm"))+
        geom_hline(yintercept = 0.8,
                   linetype="dashed",
                   color = "#084594",
                   linewidth = 1) +
        geom_vline(xintercept= 47,
                   linetype="dashed",
                   color = "#011C40FF",
                   linewidth = 1)+
        geom_vline(xintercept= 77,
                   linetype="dashed",
                   color =  "#005A32",
                   linewidth = 1)+
        geom_vline(xintercept= 182,
                   linetype="dashed",
                   color = "#FEB24C",
                   linewidth = 1)+
        theme(legend.position = "none")+
        labs(x=expression(Sample~Size~N[R]), y="Power")+
        scale_x_continuous(breaks=seq(0, 600, by=100))+
        scale_y_continuous(breaks=seq(0, 1, by = 0.2))+
        theme(title = element_text(size = 20),
              axis.title=element_text(size = 20),
              legend.title=element_text(size = 20),
              axis.text=element_text(size = 20),
              legend.text=element_text(size = 20))+
        labs(title = expression("Violation"~of~compatibility~","~pi[A] ==0.6))
      
      p_sensi3_power_06
      
    }
    
    # Power  -  07
    {
      
      Res <- readRDS(paste0("res_sensi3/res_07_power_sensi3_", delta_U,".RData"))
      
      Res <- bind_rows(Res)
      Res$n <- seq(10, 600, by = 5)
      
      Res_plot <-  Res %>% 
        pivot_longer(
          cols = c("power_EC", "power_AIPW", "power_naive"), 
          names_to = "Estimator",
          values_to = "power"
        )
      
      # add single-arm case
      Res_sa <- readRDS(paste0("res_sensi3/res_sa_power_sensi3_", delta_U,".RData"))
      
      Res_sa <- bind_rows(Res_sa)
      Res_sa$n <- seq(10, 600, by = 5) 
      
      Res_plot_sa <-  Res_sa %>% 
        pivot_longer(
          cols = c("power_SA"), 
          names_to = "Estimator",
          values_to = "power"
        )
      
      Res_plot <- rbind(Res_plot, Res_plot_sa)
      
      Res_plot$Estimator <- factor(Res_plot$Estimator, 
                                   levels = c(
                                     "power_EC",
                                     "power_AIPW",
                                     "power_naive",
                                     "power_SA"
                                   ))
      levels(Res_plot$Estimator) <- c(
        "Hybrid", 
        "AIPW",
        "Diff-in-mean",
        "Single-arm")
      
      
      p_sensi3_power_07 <- ggplot(data = Res_plot, aes(x = n, y = power, colour = Estimator)) +
        geom_line(linewidth = 1) +
        scale_colour_manual(values=c("#F27127FF","#FDC718FF",
                                     "#038C33FF", "#207CF5FF"),
                            breaks=c("Diff-in-mean", "AIPW",
                                     "Hybrid", "Single-arm"))+
        geom_hline(yintercept = 0.8,
                   linetype="dashed",
                   color = "#084594",
                   linewidth = 1) +
        geom_vline(xintercept= 47,
                   linetype="dashed",
                   color = "#011C40FF",
                   linewidth = 1)+
        geom_vline(xintercept= 66,
                   linetype="dashed",
                   color =  "#005A32",
                   linewidth = 1)+
        geom_vline(xintercept= 208,
                   linetype="dashed",
                   color = "#FEB24C",
                   linewidth = 1)+
        theme(legend.position = "none")+
        labs(x=expression(Sample~Size~N[R]), y="Power")+
        scale_x_continuous(breaks=seq(0, 600, by=100))+
        scale_y_continuous(breaks=seq(0, 1, by = 0.2))+
        theme(title = element_text(size = 20),
              axis.title=element_text(size = 20),
              legend.title=element_text(size = 20),
              axis.text=element_text(size = 20),
              legend.text=element_text(size = 20))+
        labs(title = expression("Violation"~of~compatibility~","~pi[A] ==0.7))
      
      p_sensi3_power_07
      
      
    }
    
    # Power  -  08
    {
      
      Res <- readRDS(paste0("res_sensi3/res_08_power_sensi3_", delta_U,".RData"))
      
      Res <- bind_rows(Res)
      Res$n <- seq(10, 600, by = 5)
      
      Res_plot <-  Res %>% 
        pivot_longer(
          cols = c("power_EC", "power_AIPW", "power_naive"), 
          names_to = "Estimator",
          values_to = "power"
        )
      
      # add single-arm case
      Res_sa <- readRDS(paste0("res_sensi3/res_sa_power_sensi3_", delta_U,".RData"))
      
      Res_sa <- bind_rows(Res_sa)
      Res_sa$n <- seq(10, 600, by = 5) 
      
      Res_plot_sa <-  Res_sa %>% 
        pivot_longer(
          cols = c("power_SA"), 
          names_to = "Estimator",
          values_to = "power"
        )
      
      Res_plot <- rbind(Res_plot, Res_plot_sa)
      
      Res_plot$Estimator <- factor(Res_plot$Estimator, 
                                   levels = c(
                                     "power_EC",
                                     "power_AIPW",
                                     "power_naive",
                                     "power_SA"
                                   ))
      levels(Res_plot$Estimator) <- c(
        "Hybrid", 
        "AIPW",
        "Diff-in-mean",
        "Single-arm")
      
      p_sensi3_power_08 <- ggplot(data = Res_plot, aes(x = n, y = power, colour = Estimator)) +
        geom_line(linewidth = 1) +
        scale_colour_manual(values=c("#F27127FF","#FDC718FF",
                                     "#038C33FF", "#207CF5FF"),
                            breaks=c("Diff-in-mean", "AIPW",
                                     "Hybrid", "Single-arm"))+
        geom_hline(yintercept = 0.8,
                   linetype="dashed",
                   color = "#084594",
                   linewidth = 1) +
        geom_vline(xintercept= 47,
                   linetype="dashed",
                   color = "#011C40FF",
                   linewidth = 1)+
        geom_vline(xintercept= 58,
                   linetype="dashed",
                   color =  "#005A32",
                   linewidth = 1)+
        geom_vline(xintercept= 273,
                   linetype="dashed",
                   color = "#FEB24C",
                   linewidth = 1)+
        theme(legend.position = "none")+
        labs(x=expression(Sample~Size~N[R]), y="Power")+
        scale_x_continuous(breaks=seq(0, 600, by=100))+
        scale_y_continuous(breaks=seq(0, 1, by = 0.2))+
        theme(title = element_text(size = 20),
              axis.title=element_text(size = 20),
              legend.title=element_text(size = 20),
              axis.text=element_text(size = 20),
              legend.text=element_text(size = 20))+
        labs(title = expression("Violation"~of~compatibility~","~pi[A] ==0.8))
      
      p_sensi3_power_08
      
    }
    
    # Power  -  09
    {
      
      Res <- readRDS(paste0("res_sensi3/res_09_power_sensi3_", delta_U,".RData"))
      
      Res <- bind_rows(Res)
      Res$n <- seq(10, 600, by = 5)
      
      Res_plot <-  Res %>% 
        pivot_longer(
          cols = c("power_EC", "power_AIPW", "power_naive"), 
          names_to = "Estimator",
          values_to = "power"
        )
      
      # add single-arm case
      Res_sa <- readRDS(paste0("res_sensi3/res_sa_power_sensi3_", delta_U,".RData"))
      
      Res_sa <- bind_rows(Res_sa)
      Res_sa$n <- seq(10, 600, by = 5)
      
      Res_plot_sa <-  Res_sa %>% 
        pivot_longer(
          cols = c("power_SA"), 
          names_to = "Estimator",
          values_to = "power"
        )
      
      Res_plot <- rbind(Res_plot, Res_plot_sa)
      
      Res_plot$Estimator <- factor(Res_plot$Estimator, 
                                   levels = c(
                                     "power_EC",
                                     "power_AIPW",
                                     "power_naive",
                                     "power_SA"
                                   ))
      levels(Res_plot$Estimator) <- c(
        "Hybrid", 
        "AIPW",
        "Diff-in-mean",
        "Single-arm")
      
      
      p_sensi3_power_09 <- ggplot(data = Res_plot, aes(x = n, y = power, colour = Estimator)) +
        geom_line(linewidth = 1) +
        scale_colour_manual(values=c("#F27127FF","#FDC718FF",
                                     "#038C33FF", "#207CF5FF"),
                            breaks=c("Diff-in-mean", "AIPW",
                                     "Hybrid", "Single-arm"))+
        geom_hline(yintercept = 0.8,
                   linetype="dashed",
                   color = "#084594",
                   linewidth = 1) +
        geom_vline(xintercept= 47,
                   linetype="dashed",
                   color = "#011C40FF",
                   linewidth = 1)+
        geom_vline(xintercept= 52,
                   linetype="dashed",
                   color =  "#005A32",
                   linewidth = 1)+
        geom_vline(xintercept= 486,
                   linetype="dashed",
                   color = "#FEB24C",
                   linewidth = 1)+
        # theme(legend.position = "none")+
        labs(x=expression(Sample~Size~N[R]), y="Power")+
        scale_x_continuous(breaks=seq(0, 600, by=100))+
        scale_y_continuous(breaks=seq(0, 1, by = 0.2))+
        theme(title = element_text(size = 20),
              axis.title=element_text(size = 20),
              legend.title=element_text(size = 20),
              axis.text=element_text(size = 20),
              legend.text=element_text(size = 20))+
        labs(title = expression("Violation"~of~compatibility~","~pi[A] ==0.9))
      
      p_sensi3_power_09
      
    }
    
  }
  
  ## Type I error_sufficient
  {
    
    # Type I error - 05
    {
      
      Res <- readRDS(paste0("res_sensi3/res_05_TypeI_sensi3_", delta_U,".RData"))
      
      Res <- bind_rows(Res)
      Res$n <- seq(10, 600, by = 5)
      
      Res_plot <-  Res %>% 
        pivot_longer(
          cols = c("power_EC", "power_AIPW", "power_naive"), 
          names_to = "Estimator",
          values_to = "Type_I_Error"
        )
      
      # add single-arm case
      Res_sa <- readRDS(paste0("res_sensi3/res_sa_TypeI_sensi3_", delta_U,".RData"))
      
      Res_sa <- bind_rows(Res_sa)
      Res_sa$n <- seq(10, 600, by = 5)
      
      Res_plot_sa <-  Res_sa %>% 
        pivot_longer(
          cols = c("power_SA"), 
          names_to = "Estimator",
          values_to = "Type_I_Error"
        )
      
      Res_plot <- rbind(Res_plot, Res_plot_sa)
      
      Res_plot$Estimator <- factor(Res_plot$Estimator, 
                                   levels = c(
                                     "power_EC",
                                     "power_AIPW",
                                     "power_naive",
                                     "power_SA"
                                   ))
      levels(Res_plot$Estimator) <- c(
        "Hybrid", 
        "AIPW",
        "Diff-in-mean",
        "Single-arm")
      
      
      p_sensi3_TypeI_05 <- ggplot(data = Res_plot, aes(x = n, y = Type_I_Error, colour = Estimator)) +
        geom_line(linewidth = 1) +
        scale_colour_manual(values=c("#F27127FF","#FDC718FF",
                                     "#038C33FF", "#207CF5FF"), 
                            breaks=c("Diff-in-mean", "AIPW",
                                     "Hybrid", "Single-arm"))+
        geom_hline(yintercept = 0.05,
                   linetype="dashed",
                   color = "#084594",
                   linewidth = 1) +
        geom_vline(xintercept= 47,
                   linetype="dashed",
                   color = "#011C40FF",
                   linewidth = 1)+
        geom_vline(xintercept= 92,
                   linetype="dashed",
                   color =  "#005A32",
                   linewidth = 1)+
        geom_vline(xintercept= 175,
                   linetype="dashed",
                   color = "#FEB24C",
                   linewidth = 1)+
        theme(legend.position = "none")+
        # labs(x=expression(RCT~Sample~Size~N[R]), y="Type I Error")+
        labs(x=expression(Sample~Size~N[R]), y="Type I Error")+
        scale_x_continuous(breaks=seq(0, 600, by=100))+
        scale_y_continuous(breaks = seq(0, 1, by=0.2))+
        expand_limits(y=c(0, 1))+
        theme(title = element_text(size = 20),
              axis.title=element_text(size = 20),
              legend.title=element_text(size = 20),
              axis.text=element_text(size = 20),
              legend.text=element_text(size = 20))+
        labs(title = expression("Violation"~of~compatibility~","~pi[A] ==0.5))
      
      p_sensi3_TypeI_05
      
    }
    
    # Type I error - 06
    {
      
      Res <- readRDS(paste0("res_sensi3/res_06_TypeI_sensi3_", delta_U,".RData"))
      
      Res <- bind_rows(Res)
      Res$n <- seq(10, 600, by = 5)
      
      Res_plot <-  Res %>% 
        pivot_longer(
          cols = c("power_EC", "power_AIPW", "power_naive"), 
          names_to = "Estimator",
          values_to = "Type_I_Error"
        )
      
      # add single-arm case
      Res_sa <- readRDS(paste0("res_sensi3/res_sa_TypeI_sensi3_", delta_U,".RData"))
      
      Res_sa <- bind_rows(Res_sa)
      Res_sa$n <- seq(10, 600, by = 5)
      
      Res_plot_sa <-  Res_sa %>% 
        pivot_longer(
          cols = c("power_SA"), 
          names_to = "Estimator",
          values_to = "Type_I_Error"
        )
      
      Res_plot <- rbind(Res_plot, Res_plot_sa)
      
      Res_plot$Estimator <- factor(Res_plot$Estimator, 
                                   levels = c(
                                     "power_EC",
                                     "power_AIPW",
                                     "power_naive",
                                     "power_SA"
                                   ))
      levels(Res_plot$Estimator) <- c(
        "Hybrid", 
        "AIPW",
        "Diff-in-mean",
        "Single-arm")
      
      
      p_sensi3_TypeI_06 <- ggplot(data = Res_plot, aes(x = n, y = Type_I_Error, colour = Estimator)) +
        geom_line(linewidth = 1) +
        scale_colour_manual(values=c("#F27127FF","#FDC718FF",
                                     "#038C33FF", "#207CF5FF"), 
                            breaks=c("Diff-in-mean", "AIPW",
                                     "Hybrid", "Single-arm"))+
        geom_hline(yintercept = 0.05,
                   linetype="dashed",
                   color = "#084594",
                   linewidth = 1) +
        geom_vline(xintercept= 47,
                   linetype="dashed",
                   color = "#011C40FF",
                   linewidth = 1)+
        geom_vline(xintercept= 77,
                   linetype="dashed",
                   color =  "#005A32",
                   linewidth = 1)+
        geom_vline(xintercept= 182,
                   linetype="dashed",
                   color = "#FEB24C",
                   linewidth = 1)+
        theme(legend.position = "none")+
        labs(x=expression(Sample~Size~N[R]), y="Type I Error")+
        scale_x_continuous(breaks=seq(0, 600, by=100))+
        scale_y_continuous(breaks = seq(0, 1, by=0.2))+
        expand_limits(y=c(0, 1))+
        theme(title = element_text(size = 20),
              axis.title=element_text(size = 20),
              legend.title=element_text(size = 20),
              axis.text=element_text(size = 20),
              legend.text=element_text(size = 20))+
        labs(title = expression("Violation"~of~compatibility~","~pi[A] ==0.6))
      
      p_sensi3_TypeI_06
      
      
    }
    
    # Type I error - 07
    {
      
      
      Res <- readRDS(paste0("res_sensi3/res_07_TypeI_sensi3_", delta_U,".RData"))
      
      Res <- bind_rows(Res)
      Res$n <- seq(10, 600, by = 5)
      
      Res_plot <-  Res %>% 
        pivot_longer(
          cols = c("power_EC", "power_AIPW", "power_naive"), 
          names_to = "Estimator",
          values_to = "Type_I_Error"
        )
      
      # add single-arm case
      Res_sa <- readRDS(paste0("res_sensi3/res_sa_TypeI_sensi3_", delta_U,".RData"))
      
      Res_sa <- bind_rows(Res_sa)
      Res_sa$n <- seq(10, 600, by = 5)
      
      
      Res_plot_sa <-  Res_sa %>% 
        pivot_longer(
          cols = c("power_SA"), 
          names_to = "Estimator",
          values_to = "Type_I_Error"
        )
      
      Res_plot <- rbind(Res_plot, Res_plot_sa)
      
      Res_plot$Estimator <- factor(Res_plot$Estimator, 
                                   levels = c(
                                     "power_EC",
                                     "power_AIPW",
                                     "power_naive",
                                     "power_SA"
                                   ))
      levels(Res_plot$Estimator) <- c(
        "Hybrid", 
        "AIPW",
        "Diff-in-mean",
        "Single-arm")
      
      
      p_sensi3_TypeI_07 <- ggplot(data = Res_plot, aes(x = n, y = Type_I_Error, colour = Estimator)) +
        geom_line(linewidth = 1) +
        scale_colour_manual(values=c("#F27127FF","#FDC718FF",
                                     "#038C33FF", "#207CF5FF"), 
                            breaks=c("Diff-in-mean", "AIPW",
                                     "Hybrid", "Single-arm"))+
        geom_hline(yintercept = 0.05,
                   linetype="dashed",
                   color = "#084594",
                   linewidth = 1) +
        geom_vline(xintercept= 47,
                   linetype="dashed",
                   color = "#011C40FF",
                   linewidth = 1)+
        geom_vline(xintercept= 66,
                   linetype="dashed",
                   color =  "#005A32",
                   linewidth = 1)+
        geom_vline(xintercept= 208,
                   linetype="dashed",
                   color = "#FEB24C",
                   linewidth = 1)+
        theme(legend.position = "none")+
        labs(x=expression(Sample~Size~N[R]), y="Type I Error")+
        scale_x_continuous(breaks=seq(0, 600, by=100))+
        scale_y_continuous(breaks = seq(0, 1, by=0.2))+
        expand_limits(y=c(0, 1))+
        theme(title = element_text(size = 20),
              axis.title=element_text(size = 20),
              legend.title=element_text(size = 20),
              axis.text=element_text(size = 20),
              legend.text=element_text(size = 20))+
        labs(title = expression("Violation"~of~compatibility~","~pi[A] ==0.7))
      
      p_sensi3_TypeI_07
      
    }
    
    # Type I error - 08
    {
      
      
      Res <- readRDS(paste0("res_sensi3/res_08_TypeI_sensi3_", delta_U,".RData"))
      
      Res <- bind_rows(Res)
      Res$n <- seq(10, 600, by = 5)
      
      Res_plot <-  Res %>% 
        pivot_longer(
          cols = c("power_EC", "power_AIPW", "power_naive"), 
          names_to = "Estimator",
          values_to = "Type_I_Error"
        )
      
      # add single-arm case
      Res_sa <- readRDS(paste0("res_sensi3/res_sa_TypeI_sensi3_", delta_U,".RData"))
      
      Res_sa <- bind_rows(Res_sa)
      Res_sa$n <- seq(10, 600, by = 5)
      
      
      Res_plot_sa <-  Res_sa %>% 
        pivot_longer(
          cols = c("power_SA"), 
          names_to = "Estimator",
          values_to = "Type_I_Error"
        )
      
      Res_plot <- rbind(Res_plot, Res_plot_sa)
      
      Res_plot$Estimator <- factor(Res_plot$Estimator, 
                                   levels = c(
                                     "power_EC",
                                     "power_AIPW",
                                     "power_naive",
                                     "power_SA"
                                   ))
      levels(Res_plot$Estimator) <- c(
        "Hybrid", 
        "AIPW",
        "Diff-in-mean",
        "Single-arm")
      
      
      p_sensi3_TypeI_08 <- ggplot(data = Res_plot, aes(x = n, y = Type_I_Error, colour = Estimator)) +
        geom_line(linewidth = 1) +
        scale_colour_manual(values=c("#F27127FF","#FDC718FF",
                                     "#038C33FF", "#207CF5FF"), 
                            breaks=c("Diff-in-mean", "AIPW",
                                     "Hybrid", "Single-arm"))+
        geom_hline(yintercept = 0.05,
                   linetype="dashed",
                   color = "#084594",
                   linewidth = 1) +
        geom_vline(xintercept= 47,
                   linetype="dashed",
                   color = "#011C40FF",
                   linewidth = 1)+
        geom_vline(xintercept= 58,
                   linetype="dashed",
                   color =  "#005A32",
                   linewidth = 1)+
        geom_vline(xintercept= 273,
                   linetype="dashed",
                   color = "#FEB24C",
                   linewidth = 1)+
        theme(legend.position = "none")+
        labs(x=expression(Sample~Size~N[R]), y="Type I Error")+
        scale_x_continuous(breaks=seq(0, 600, by=100))+
        scale_y_continuous(breaks = seq(0, 1, by=0.2))+
        expand_limits(y=c(0, 1))+
        theme(title = element_text(size = 20),
              axis.title=element_text(size = 20),
              legend.title=element_text(size = 20),
              axis.text=element_text(size = 20),
              legend.text=element_text(size = 20))+
        labs(title = expression("Violation"~of~compatibility~","~pi[A] ==0.8))
      
      p_sensi3_TypeI_08
      
    }
    
    # Type I error - 09
    {
      
      Res <- readRDS(paste0("res_sensi3/res_09_TypeI_sensi3_", delta_U,".RData"))
      
      Res <- bind_rows(Res)
      Res$n <- seq(10, 600, by = 5)
      
      Res_plot <-  Res %>% 
        pivot_longer(
          cols = c("power_EC", "power_AIPW", "power_naive"), 
          names_to = "Estimator",
          values_to = "Type_I_Error"
        )
      
      # add single-arm case
      Res_sa <- readRDS(paste0("res_sensi3/res_sa_TypeI_sensi3_", delta_U,".RData"))
      
      Res_sa <- bind_rows(Res_sa)
      Res_sa$n <- seq(10, 600, by = 5)
      
      
      Res_plot_sa <-  Res_sa %>% 
        pivot_longer(
          cols = c("power_SA"), 
          names_to = "Estimator",
          values_to = "Type_I_Error"
        )
      
      Res_plot <- rbind(Res_plot, Res_plot_sa)
      
      Res_plot$Estimator <- factor(Res_plot$Estimator, 
                                   levels = c(
                                     "power_EC",
                                     "power_AIPW",
                                     "power_naive",
                                     "power_SA"
                                   ))
      levels(Res_plot$Estimator) <- c(
        "Hybrid", 
        "AIPW",
        "Diff-in-mean",
        "Single-arm")
      
      
      p_sensi3_TypeI_09 <- ggplot(data = Res_plot, aes(x = n, y = Type_I_Error, colour = Estimator)) +
        geom_line(linewidth = 1) +
        scale_colour_manual(values=c("#F27127FF","#FDC718FF",
                                     "#038C33FF", "#207CF5FF"), 
                            breaks=c("Diff-in-mean", "AIPW",
                                     "Hybrid", "Single-arm"))+
        geom_hline(yintercept = 0.05,
                   linetype="dashed",
                   color = "#084594",
                   linewidth = 1) +
        geom_vline(xintercept= 47,
                   linetype="dashed",
                   color = "#011C40FF",
                   linewidth = 1)+
        geom_vline(xintercept= 52,
                   linetype="dashed",
                   color =  "#005A32",
                   linewidth = 1)+
        geom_vline(xintercept= 486,
                   linetype="dashed",
                   color = "#FEB24C",
                   linewidth = 1)+
        # theme(legend.position = "none")+
        labs(x=expression(Sample~Size~N[R]), y="Type I Error")+
        scale_x_continuous(breaks=seq(0, 600, by=100))+
        scale_y_continuous(breaks = seq(0, 1, by=0.2))+
        expand_limits(y=c(0, 1))+
        theme(title = element_text(size = 20),
              axis.title=element_text(size = 20),
              legend.title=element_text(size = 20),
              axis.text=element_text(size = 20),
              legend.text=element_text(size = 20))+
        labs(title = expression("Violation"~of~compatibility~","~pi[A] ==0.9))
      
      p_sensi3_TypeI_09
      
    }
    
  }
  
}

### Sensitivity analysis 4 (N_E = 1000)
{
  ## Power sufficient
  {
    
    # Power  -  05
    {
      
      Res <- readRDS("res_sensi4/res_05_power_sensi4.RData")
      
      Res <- bind_rows(Res)
      Res$n <- seq(10, 600, by = 5)
      
      Res_plot <-  Res %>% 
        pivot_longer(
          cols = c("power_EC", "power_AIPW", "power_naive"), 
          names_to = "Estimator",
          values_to = "power"
        )
      
      # add single-arm case
      Res_sa <- readRDS("res_sensi4/res_sa_power_sensi4.RData")
      
      Res_sa <- bind_rows(Res_sa)
      Res_sa$n <- seq(10, 600, by = 5)
      
      Res_plot_sa <-  Res_sa %>% 
        pivot_longer(
          cols = c("power_SA"), 
          names_to = "Estimator",
          values_to = "power"
        )
      
      Res_plot <- rbind(Res_plot, Res_plot_sa)
      
      Res_plot$Estimator <- factor(Res_plot$Estimator, 
                                   levels = c(
                                     "power_EC",
                                     "power_AIPW",
                                     "power_naive",
                                     "power_SA"
                                   ))
      levels(Res_plot$Estimator) <- c(
        "Hybrid", 
        "AIPW",
        "Diff-in-mean",
        "Single-arm")
      
      p_sensi4_power_05 <- ggplot(data = Res_plot, aes(x = n, y = power, colour = Estimator)) +
        geom_line(linewidth = 1) +
        scale_colour_manual(values=c("#F27127FF","#FDC718FF",
                                     "#038C33FF", "#207CF5FF"),
                            breaks=c("Diff-in-mean", "AIPW",
                                     "Hybrid", "Single-arm"))+
        geom_hline(yintercept = 0.8,
                   linetype="dashed",
                   color = "#084594",
                   linewidth = 1) +
        geom_vline(xintercept= 42,
                   linetype="dashed",
                   color = "#011C40FF",
                   linewidth = 1)+
        geom_vline(xintercept= 84,
                   linetype="dashed",
                   color =  "#005A32",
                   linewidth = 1) +
        geom_vline(xintercept= 157,
                   linetype="dashed",
                   color = "#FEB24C",
                   linewidth = 1)+
        theme(legend.position = "none")+
        # labs(x=expression(RCT~Sample~Size~N[R]), y="Power")+
        labs(x=expression(Sample~Size~N[R]), y="Power")+
        scale_x_continuous(breaks=seq(0, 600, by=100))+
        scale_y_continuous(breaks=seq(0, 1, by = 0.2))+
        theme(title = element_text(size = 20),
              axis.title=element_text(size = 20),
              legend.title=element_text(size = 20),
              axis.text=element_text(size = 20),
              legend.text=element_text(size = 20))+
        # labs(title = expression(Optimistic~specification~of~d(X)~","~pi[A] ==0.5))
        labs(title = expression(Optimistic~specification~of~d(X)~","~pi[A] ==0.5))
        
      
      p_sensi4_power_05
      
    }
    
    # Power  -  06
    {
      
      
      Res <- readRDS("res_sensi4/res_06_power_sensi4.RData")
      
      Res <- bind_rows(Res)
      Res$n <- seq(10, 600, by = 5)
      
      Res_plot <-  Res %>% 
        pivot_longer(
          cols = c("power_EC", "power_AIPW", "power_naive"), 
          names_to = "Estimator",
          values_to = "power"
        )
      
      # add single-arm case
      Res_sa <- readRDS("res_sensi4/res_sa_power_sensi4.RData")
      
      Res_sa <- bind_rows(Res_sa)
      Res_sa$n <- seq(10, 600, by = 5)
      
      
      Res_plot_sa <-  Res_sa %>% 
        pivot_longer(
          cols = c("power_SA"), 
          names_to = "Estimator",
          values_to = "power"
        )
      
      Res_plot <- rbind(Res_plot, Res_plot_sa)
      
      Res_plot$Estimator <- factor(Res_plot$Estimator, 
                                   levels = c(
                                     "power_EC",
                                     "power_AIPW",
                                     "power_naive",
                                     "power_SA"
                                   ))
      levels(Res_plot$Estimator) <- c(
        "Hybrid", 
        "AIPW",
        "Diff-in-mean",
        "Single-arm")
      
      p_sensi4_power_06 <- ggplot(data = Res_plot, aes(x = n, y = power, colour = Estimator)) +
        geom_line(linewidth = 1) +
        scale_colour_manual(values=c("#F27127FF","#FDC718FF",
                                     "#038C33FF", "#207CF5FF"),
                            breaks=c("Diff-in-mean", "AIPW",
                                     "Hybrid", "Single-arm"))+
        geom_hline(yintercept = 0.8,
                   linetype="dashed",
                   color = "#084594",
                   linewidth = 1) +
        geom_vline(xintercept= 42,
                   linetype="dashed",
                   color = "#011C40FF",
                   linewidth = 1)+
        geom_vline(xintercept= 70,
                   linetype="dashed",
                   color =  "#005A32",
                   linewidth = 1)+
        geom_vline(xintercept= 164,
                   linetype="dashed",
                   color = "#FEB24C",
                   linewidth = 1)+
        theme(legend.position = "none")+
        labs(x=expression(Sample~Size~N[R]), y="Power")+
        scale_x_continuous(breaks=seq(0, 600, by=100))+
        scale_y_continuous(breaks=seq(0, 1, by = 0.2))+
        theme(title = element_text(size = 20),
              axis.title=element_text(size = 20),
              legend.title=element_text(size = 20),
              axis.text=element_text(size = 20),
              legend.text=element_text(size = 20))+
        labs(title = expression(Optimistic~specification~of~d(X)~","~pi[A] ==0.6))
      
      p_sensi4_power_06
      
    }
    
    # Power  -  07
    {
      
      Res <- readRDS("res_sensi4/res_07_power_sensi4.RData")
      
      Res <- bind_rows(Res)
      Res$n <- seq(10, 600, by = 5)
      
      Res_plot <-  Res %>% 
        pivot_longer(
          cols = c("power_EC", "power_AIPW", "power_naive"), 
          names_to = "Estimator",
          values_to = "power"
        )
      
      # add single-arm case
      Res_sa <- readRDS("res_sensi4/res_sa_power_sensi4.RData")
      
      Res_sa <- bind_rows(Res_sa)
      Res_sa$n <- seq(10, 600, by = 5)
      
      
      Res_plot_sa <-  Res_sa %>% 
        pivot_longer(
          cols = c("power_SA"), 
          names_to = "Estimator",
          values_to = "power"
        )
      
      Res_plot <- rbind(Res_plot, Res_plot_sa)
      
      Res_plot$Estimator <- factor(Res_plot$Estimator, 
                                   levels = c(
                                     "power_EC",
                                     "power_AIPW",
                                     "power_naive",
                                     "power_SA"
                                   ))
      levels(Res_plot$Estimator) <- c(
        "Hybrid", 
        "AIPW",
        "Diff-in-mean",
        "Single-arm")
      
      p_sensi4_power_07 <- ggplot(data = Res_plot, aes(x = n, y = power, colour = Estimator)) +
        geom_line(linewidth = 1) +
        scale_colour_manual(values=c("#F27127FF","#FDC718FF",
                                     "#038C33FF", "#207CF5FF"),
                            breaks=c("Diff-in-mean", "AIPW",
                                     "Hybrid", "Single-arm"))+
        geom_hline(yintercept = 0.8,
                   linetype="dashed",
                   color = "#084594",
                   linewidth = 1) +
        geom_vline(xintercept= 42,
                   linetype="dashed",
                   color = "#011C40FF",
                   linewidth = 1)+
        geom_vline(xintercept= 60,
                   linetype="dashed",
                   color =  "#005A32",
                   linewidth = 1)+
        geom_vline(xintercept= 187,
                   linetype="dashed",
                   color = "#FEB24C",
                   linewidth = 1)+
        theme(legend.position = "none")+
        labs(x=expression(Sample~Size~N[R]), y="Power")+
        scale_x_continuous(breaks=seq(0, 600, by=100))+
        scale_y_continuous(breaks=seq(0, 1, by = 0.2))+
        theme(title = element_text(size = 20),
              axis.title=element_text(size = 20),
              legend.title=element_text(size = 20),
              axis.text=element_text(size = 20),
              legend.text=element_text(size = 20))+
        labs(title = expression(Optimistic~specification~of~d(X)~","~pi[A] ==0.7))
      
      p_sensi4_power_07
      
    }
    
    # Power  -  08
    {
      
      Res <- readRDS("res_sensi4/res_08_power_sensi4.RData")
      
      Res <- bind_rows(Res)
      Res$n <- seq(10, 600, by = 5)
      
      Res_plot <-  Res %>% 
        pivot_longer(
          cols = c("power_EC", "power_AIPW", "power_naive"), 
          names_to = "Estimator",
          values_to = "power"
        )
      
      # add single-arm case
      Res_sa <- readRDS("res_sensi4/res_sa_power_sensi4.RData")
      
      Res_sa <- bind_rows(Res_sa)
      Res_sa$n <- seq(10, 600, by = 5)
      
      Res_plot_sa <-  Res_sa %>% 
        pivot_longer(
          cols = c("power_SA"), 
          names_to = "Estimator",
          values_to = "power"
        )
      
      Res_plot <- rbind(Res_plot, Res_plot_sa)
      
      Res_plot$Estimator <- factor(Res_plot$Estimator, 
                                   levels = c(
                                     "power_EC",
                                     "power_AIPW",
                                     "power_naive",
                                     "power_SA"
                                   ))
      levels(Res_plot$Estimator) <- c(
        "Hybrid", 
        "AIPW",
        "Diff-in-mean",
        "Single-arm")
      
      p_sensi4_power_08 <- ggplot(data = Res_plot, aes(x = n, y = power, colour = Estimator)) +
        geom_line(linewidth = 1) +
        scale_colour_manual(values=c("#F27127FF","#FDC718FF",
                                     "#038C33FF", "#207CF5FF"),
                            breaks=c("Diff-in-mean", "AIPW",
                                     "Hybrid", "Single-arm"))+
        geom_hline(yintercept = 0.8,
                   linetype="dashed",
                   color = "#084594",
                   linewidth = 1) +
        geom_vline(xintercept= 42,
                   linetype="dashed",
                   color = "#011C40FF",
                   linewidth = 1)+
        geom_vline(xintercept= 53,
                   linetype="dashed",
                   color =  "#005A32",
                   linewidth = 1)+
        geom_vline(xintercept= 246,
                   linetype="dashed",
                   color = "#FEB24C",
                   linewidth = 1)+
        theme(legend.position = "none")+
        labs(x=expression(Sample~Size~N[R]), y="Power")+
        scale_x_continuous(breaks=seq(0, 600, by=100))+
        scale_y_continuous(breaks=seq(0, 1, by = 0.2))+
        theme(title = element_text(size = 20),
              axis.title=element_text(size = 20),
              legend.title=element_text(size = 20),
              axis.text=element_text(size = 20),
              legend.text=element_text(size = 20))+
        labs(title = expression(Optimistic~specification~of~d(X)~","~pi[A] ==0.8))
      
      p_sensi4_power_08
      
    }
    
    # Power  -  09
    {
      
      Res <- readRDS("res_sensi4/res_09_power_sensi4.RData")
      
      Res <- bind_rows(Res)
      Res$n <- seq(10, 600, by = 5)
      
      Res_plot <-  Res %>% 
        pivot_longer(
          cols = c("power_EC", "power_AIPW", "power_naive"), 
          names_to = "Estimator",
          values_to = "power"
        )
      
      # add single-arm case
      Res_sa <- readRDS("res_sensi4/res_sa_power_sensi4.RData")
      
      Res_sa <- bind_rows(Res_sa)
      Res_sa$n <- seq(10, 600, by = 5)
      
      
      Res_plot_sa <-  Res_sa %>% 
        pivot_longer(
          cols = c("power_SA"), 
          names_to = "Estimator",
          values_to = "power"
        )
      
      Res_plot <- rbind(Res_plot, Res_plot_sa)
      
      Res_plot$Estimator <- factor(Res_plot$Estimator, 
                                   levels = c(
                                     "power_EC",
                                     "power_AIPW",
                                     "power_naive",
                                     "power_SA"
                                   ))
      levels(Res_plot$Estimator) <- c(
        "Hybrid", 
        "AIPW",
        "Diff-in-mean",
        "Single-arm")
      
      p_sensi4_power_09 <- ggplot(data = Res_plot, aes(x = n, y = power, colour = Estimator)) +
        geom_line(linewidth = 1) +
        scale_colour_manual(values=c("#F27127FF","#FDC718FF",
                                     "#038C33FF", "#207CF5FF"),
                            breaks=c("Diff-in-mean", "AIPW",
                                     "Hybrid", "Single-arm"))+
        geom_hline(yintercept = 0.8,
                   linetype="dashed",
                   color = "#084594",
                   linewidth = 1) +
        geom_vline(xintercept= 42,
                   linetype="dashed",
                   color = "#011C40FF",
                   linewidth = 1)+
        geom_vline(xintercept= 47,
                   linetype="dashed",
                   color =  "#005A32",
                   linewidth = 1)+
        geom_vline(xintercept= 437,
                   linetype="dashed",
                   color = "#FEB24C",
                   linewidth = 1)+
        # theme(legend.position = "none")+
        labs(x=expression(Sample~Size~N[R]), y="Power")+
        scale_x_continuous(breaks=seq(0, 600, by=100))+
        scale_y_continuous(breaks=seq(0, 1, by = 0.2))+
        theme(title = element_text(size = 20),
              axis.title=element_text(size = 20),
              legend.title=element_text(size = 20),
              axis.text=element_text(size = 20),
              legend.text=element_text(size = 20))+
        labs(title = expression(Optimistic~specification~of~d(X)~","~pi[A] ==0.9))
      
      p_sensi4_power_09
      
    }
    
  }
  
  ## Type I error sufficient
  {
    
    # Type I error - 05
    {
      
      Res <- readRDS("res_sensi4/res_05_TypeI_sensi4.RData")
      
      Res <- bind_rows(Res)
      Res$n <- seq(10, 600, by = 5)
      
      Res_plot <-  Res %>% 
        pivot_longer(
          cols = c("power_EC", "power_AIPW", "power_naive"), 
          names_to = "Estimator",
          values_to = "Type_I_Error"
        )
      
      # add single-arm case
      Res_sa <- readRDS("res_sensi4/res_sa_TypeI_sensi4.RData")
      
      Res_sa <- bind_rows(Res_sa)
      Res_sa$n <- seq(10, 600, by = 5)
      
      
      Res_plot_sa <-  Res_sa %>% 
        pivot_longer(
          cols = c("power_SA"), 
          names_to = "Estimator",
          values_to = "Type_I_Error"
        )
      
      Res_plot <- rbind(Res_plot, Res_plot_sa)
      
      Res_plot$Estimator <- factor(Res_plot$Estimator, 
                                   levels = c(
                                     "power_EC",
                                     "power_AIPW",
                                     "power_naive",
                                     "power_SA"
                                   ))
      levels(Res_plot$Estimator) <- c(
        "Hybrid", 
        "AIPW",
        "Diff-in-mean",
        "Single-arm")
      
      p_sensi4_TypeI_05 <- ggplot(data = Res_plot, aes(x = n, y = Type_I_Error, colour = Estimator)) +
        geom_line(linewidth = 1) +
        scale_colour_manual(values=c("#F27127FF","#FDC718FF",
                                     "#038C33FF", "#207CF5FF"), 
                            breaks=c("Diff-in-mean", "AIPW",
                                     "Hybrid", "Single-arm"))+
        geom_hline(yintercept = 0.05,
                   linetype="dashed",
                   color = "#084594",
                   linewidth = 1) +
        geom_vline(xintercept= 42,
                   linetype="dashed",
                   color = "#011C40FF",
                   linewidth = 1)+
        geom_vline(xintercept= 84,
                   linetype="dashed",
                   color =  "#005A32",
                   linewidth = 1)+
        geom_vline(xintercept= 157,
                   linetype="dashed",
                   color = "#FEB24C",
                   linewidth = 1)+
        theme(legend.position = "none")+
        labs(x=expression(Sample~Size~N[R]), y="Type I Error")+
        scale_x_continuous(breaks=seq(0, 600, by=100))+
        scale_y_continuous(breaks = seq(0, 0.546, by=0.1))+
        expand_limits(y=c(0, 0.546))+
        theme(title = element_text(size = 20),
              axis.title=element_text(size = 20),
              legend.title=element_text(size = 20),
              axis.text=element_text(size = 20),
              legend.text=element_text(size = 20))+
        labs(title = expression(Optimistic~specification~of~d(X)~","~pi[A] ==0.5))
      
      p_sensi4_TypeI_05
      
    }
    
    # Type I error - 06
    {
      Res <- readRDS("res_sensi4/res_06_TypeI_sensi4.RData")
      
      Res <- bind_rows(Res)
      Res$n <- seq(10, 600, by = 5)
      
      Res_plot <-  Res %>% 
        pivot_longer(
          cols = c("power_EC", "power_AIPW", "power_naive"), 
          names_to = "Estimator",
          values_to = "Type_I_Error"
        )
      
      # add single-arm case
      Res_sa <- readRDS("res_sensi4/res_sa_TypeI_sensi4.RData")
      
      Res_sa <- bind_rows(Res_sa)
      Res_sa$n <- seq(10, 600, by = 5)
      
      Res_plot_sa <-  Res_sa %>% 
        pivot_longer(
          cols = c("power_SA"), 
          names_to = "Estimator",
          values_to = "Type_I_Error"
        )
      
      Res_plot <- rbind(Res_plot, Res_plot_sa)
      
      Res_plot$Estimator <- factor(Res_plot$Estimator, 
                                   levels = c(
                                     "power_EC",
                                     "power_AIPW",
                                     "power_naive",
                                     "power_SA"
                                   ))
      levels(Res_plot$Estimator) <- c(
        "Hybrid", 
        "AIPW",
        "Diff-in-mean",
        "Single-arm")
      
      
      p_sensi4_TypeI_06 <- ggplot(data = Res_plot, aes(x = n, y = Type_I_Error, colour = Estimator)) +
        geom_line(linewidth = 1) +
        scale_colour_manual(values=c("#F27127FF","#FDC718FF",
                                     "#038C33FF", "#207CF5FF"), 
                            breaks=c("Diff-in-mean", "AIPW",
                                     "Hybrid", "Single-arm"))+
        geom_hline(yintercept = 0.05,
                   linetype="dashed",
                   color = "#084594",
                   linewidth = 1) +
        geom_vline(xintercept= 42,
                   linetype="dashed",
                   color = "#011C40FF",
                   linewidth = 1)+
        geom_vline(xintercept= 70,
                   linetype="dashed",
                   color =  "#005A32",
                   linewidth = 1)+
        geom_vline(xintercept= 164,
                   linetype="dashed",
                   color = "#FEB24C",
                   linewidth = 1)+
        theme(legend.position = "none")+
        labs(x=expression(Sample~Size~N[R]), y="Type I Error")+
        scale_x_continuous(breaks=seq(0, 600, by=100))+
        scale_y_continuous(breaks = seq(0, 0.546, by=0.1))+
        expand_limits(y=c(0, 0.546))+
        theme(title = element_text(size = 20),
              axis.title=element_text(size = 20),
              legend.title=element_text(size = 20),
              axis.text=element_text(size = 20),
              legend.text=element_text(size = 20))+
        labs(title = expression(Optimistic~specification~of~d(X)~","~pi[A] ==0.6))
      
      p_sensi4_TypeI_06
      
    }
    
    # Type I error - 07
    {
      Res <- readRDS("res_sensi4/res_07_TypeI_sensi4.RData")
      
      Res <- bind_rows(Res)
      Res$n <- seq(10, 600, by = 5)
      
      Res_plot <-  Res %>% 
        pivot_longer(
          cols = c("power_EC", "power_AIPW", "power_naive"), 
          names_to = "Estimator",
          values_to = "Type_I_Error"
        )
      
      # add single-arm case
      Res_sa <- readRDS("res_sensi4/res_sa_TypeI_sensi4.RData")
      
      Res_sa <- bind_rows(Res_sa)
      Res_sa$n <- seq(10, 600, by = 5) 
      
      Res_plot_sa <-  Res_sa %>% 
        pivot_longer(
          cols = c("power_SA"), 
          names_to = "Estimator",
          values_to = "Type_I_Error"
        )
      
      Res_plot <- rbind(Res_plot, Res_plot_sa)
      
      Res_plot$Estimator <- factor(Res_plot$Estimator, 
                                   levels = c(
                                     "power_EC",
                                     "power_AIPW",
                                     "power_naive",
                                     "power_SA"
                                   ))
      levels(Res_plot$Estimator) <- c(
        "Hybrid", 
        "AIPW",
        "Diff-in-mean",
        "Single-arm")
      
      
      p_sensi4_TypeI_07 <- ggplot(data = Res_plot, aes(x = n, y = Type_I_Error, colour = Estimator)) +
        geom_line(linewidth = 1) +
        scale_colour_manual(values=c("#F27127FF","#FDC718FF",
                                     "#038C33FF", "#207CF5FF"), 
                            breaks=c("Diff-in-mean", "AIPW",
                                     "Hybrid", "Single-arm"))+
        geom_hline(yintercept = 0.05,
                   linetype="dashed",
                   color = "#084594",
                   linewidth = 1) +
        geom_vline(xintercept= 42,
                   linetype="dashed",
                   color = "#011C40FF",
                   linewidth = 1)+
        geom_vline(xintercept= 60,
                   linetype="dashed",
                   color =  "#005A32",
                   linewidth = 1)+
        geom_vline(xintercept= 187,
                   linetype="dashed",
                   color = "#FEB24C",
                   linewidth = 1)+
        theme(legend.position = "none")+
        labs(x=expression(Sample~Size~N[R]), y="Type I Error")+
        scale_x_continuous(breaks=seq(0, 600, by=100))+
        scale_y_continuous(breaks = seq(0, 0.546, by=0.1))+
        expand_limits(y=c(0, 0.546))+
        theme(title = element_text(size = 20),
              axis.title=element_text(size = 20),
              legend.title=element_text(size = 20),
              axis.text=element_text(size = 20),
              legend.text=element_text(size = 20))+
        labs(title = expression(Optimistic~specification~of~d(X)~","~pi[A] ==0.7))
      
      p_sensi4_TypeI_07
      
    }
    
    # Type I error - 08
    {
      Res <- readRDS("res_sensi4/res_08_TypeI_sensi4.RData")
      
      Res <- bind_rows(Res)
      Res$n <- seq(10, 600, by = 5)
      
      Res_plot <-  Res %>% 
        pivot_longer(
          cols = c("power_EC", "power_AIPW", "power_naive"), 
          names_to = "Estimator",
          values_to = "Type_I_Error"
        )
      
      # add single-arm case
      Res_sa <- readRDS("res_sensi4/res_sa_TypeI_sensi4.RData")
      
      Res_sa <- bind_rows(Res_sa)
      Res_sa$n <- seq(10, 600, by = 5)
      
      Res_plot_sa <-  Res_sa %>% 
        pivot_longer(
          cols = c("power_SA"), 
          names_to = "Estimator",
          values_to = "Type_I_Error"
        )
      
      Res_plot <- rbind(Res_plot, Res_plot_sa)
      
      Res_plot$Estimator <- factor(Res_plot$Estimator, 
                                   levels = c(
                                     "power_EC",
                                     "power_AIPW",
                                     "power_naive",
                                     "power_SA"
                                   ))
      levels(Res_plot$Estimator) <- c(
        "Hybrid", 
        "AIPW",
        "Diff-in-mean",
        "Single-arm")
      
      
      p_sensi4_TypeI_08 <- ggplot(data = Res_plot, aes(x = n, y = Type_I_Error, colour = Estimator)) +
        geom_line(linewidth = 1) +
        scale_colour_manual(values=c("#F27127FF","#FDC718FF",
                                     "#038C33FF", "#207CF5FF"), 
                            breaks=c("Diff-in-mean", "AIPW",
                                     "Hybrid", "Single-arm"))+
        geom_hline(yintercept = 0.05,
                   linetype="dashed",
                   color = "#084594",
                   linewidth = 1) +
        geom_vline(xintercept= 42,
                   linetype="dashed",
                   color = "#011C40FF",
                   linewidth = 1)+
        geom_vline(xintercept= 53,
                   linetype="dashed",
                   color =  "#005A32",
                   linewidth = 1)+
        geom_vline(xintercept= 246,
                   linetype="dashed",
                   color = "#FEB24C",
                   linewidth = 1)+
        theme(legend.position = "none")+
        labs(x=expression(Sample~Size~N[R]), y="Type I Error")+
        scale_x_continuous(breaks=seq(0, 600, by=100))+
        scale_y_continuous(breaks = seq(0, 0.546, by=0.1))+
        expand_limits(y=c(0, 0.546))+
        theme(title = element_text(size = 20),
              axis.title=element_text(size = 20),
              legend.title=element_text(size = 20),
              axis.text=element_text(size = 20),
              legend.text=element_text(size = 20))+
        labs(title = expression(Optimistic~specification~of~d(X)~","~pi[A] ==0.8))
      
      p_sensi4_TypeI_08
      
    }
    
    # Type I error - 09
    {
      Res <- readRDS("res_sensi4/res_09_TypeI_sensi4.RData")
      
      Res <- bind_rows(Res)
      Res$n <- seq(10, 600, by = 5)
      
      Res_plot <-  Res %>% 
        pivot_longer(
          cols = c("power_EC", "power_AIPW", "power_naive"), 
          names_to = "Estimator",
          values_to = "Type_I_Error"
        )
      
      # add single-arm case
      Res_sa <- readRDS("res_sensi4/res_sa_TypeI_sensi4.RData")
      
      Res_sa <- bind_rows(Res_sa)
      Res_sa$n <- seq(10, 600, by = 5)
      
      Res_plot_sa <-  Res_sa %>% 
        pivot_longer(
          cols = c("power_SA"), 
          names_to = "Estimator",
          values_to = "Type_I_Error"
        )
      
      Res_plot <- rbind(Res_plot, Res_plot_sa)
      
      Res_plot$Estimator <- factor(Res_plot$Estimator, 
                                   levels = c(
                                     "power_EC",
                                     "power_AIPW",
                                     "power_naive",
                                     "power_SA"
                                   ))
      levels(Res_plot$Estimator) <- c(
        "Hybrid", 
        "AIPW",
        "Diff-in-mean",
        "Single-arm")
      
      
      p_sensi4_TypeI_09 <- ggplot(data = Res_plot, aes(x = n, y = Type_I_Error, colour = Estimator)) +
        geom_line(linewidth = 1) +
        scale_colour_manual(values=c("#F27127FF","#FDC718FF",
                                     "#038C33FF", "#207CF5FF"), 
                            breaks=c("Diff-in-mean", "AIPW",
                                     "Hybrid", "Single-arm"))+
        geom_hline(yintercept = 0.05,
                   linetype="dashed",
                   color = "#084594",
                   linewidth = 1) +
        geom_vline(xintercept= 42,
                   linetype="dashed",
                   color = "#011C40FF",
                   linewidth = 1)+
        geom_vline(xintercept= 47,
                   linetype="dashed",
                   color =  "#005A32",
                   linewidth = 1)+
        geom_vline(xintercept= 437,
                   linetype="dashed",
                   color = "#FEB24C",
                   linewidth = 1)+
        # theme(legend.position = "none")+
        labs(x=expression(Sample~Size~N[R]), y="Type I Error")+
        scale_x_continuous(breaks=seq(0, 600, by=100))+
        scale_y_continuous(breaks = seq(0, 0.546, by=0.1))+
        expand_limits(y=c(0, 0.546))+
        theme(title = element_text(size = 20),
              axis.title=element_text(size = 20),
              legend.title=element_text(size = 20),
              axis.text=element_text(size = 20),
              legend.text=element_text(size = 20))+
        labs(title = expression(Optimistic~specification~of~d(X)~","~pi[A] ==0.9))
      
      p_sensi4_TypeI_09
      
    }
    
  }
}


### Sensitivity analysis 5 (N_E = 1000)
{
  
  ## Power sufficient
  {
    
    # Power  -  05
    {
      
      Res <- readRDS("res_sensi5/res_05_power_sensi5.RData")
      
      Res <- bind_rows(Res)
      Res$n <- seq(10, 600, by = 5)
      
      Res_plot <-  Res %>% 
        pivot_longer(
          cols = c("power_EC", "power_AIPW", "power_naive"), 
          names_to = "Estimator",
          values_to = "power"
        )
      
      # add single-arm case
      Res_sa <- readRDS("res_sensi5/res_sa_power_sensi5.RData")
      
      Res_sa <- bind_rows(Res_sa)
      Res_sa$n <- seq(10, 600, by = 5)
      
      Res_plot_sa <-  Res_sa %>% 
        pivot_longer(
          cols = c("power_SA"), 
          names_to = "Estimator",
          values_to = "power"
        )
      
      Res_plot <- rbind(Res_plot, Res_plot_sa)
      
      Res_plot$Estimator <- factor(Res_plot$Estimator, 
                                   levels = c(
                                     "power_EC",
                                     "power_AIPW",
                                     "power_naive",
                                     "power_SA"
                                   ))
      levels(Res_plot$Estimator) <- c(
        "Hybrid", 
        "AIPW",
        "Diff-in-mean",
        "Single-arm")
      
      p_sensi5_power_05 <- ggplot(data = Res_plot, aes(x = n, y = power, colour = Estimator)) +
        geom_line(linewidth = 1) +
        scale_colour_manual(values=c("#F27127FF","#FDC718FF",
                                     "#038C33FF", "#207CF5FF"),
                            breaks=c("Diff-in-mean", "AIPW",
                                     "Hybrid", "Single-arm"))+
        geom_hline(yintercept = 0.8,
                   linetype="dashed",
                   color = "#084594",
                   linewidth = 1) +
        geom_vline(xintercept= 56,
                   linetype="dotted",
                   color = "#011C40FF",
                   linewidth = 1)+
        geom_vline(xintercept= 111,
                   linetype="dashed",
                   color =  "#005A32",
                   linewidth = 1) +
        geom_vline(xintercept= 207,
                   linetype="dashed",
                   color = "#FEB24C",
                   linewidth = 1)+
        theme(legend.position = "none")+
        # labs(x=expression(RCT~Sample~Size~N[R]), y="Power")+
        labs(x=expression(Sample~Size~N[R]), y="Power")+
        scale_x_continuous(breaks=seq(0, 600, by=100))+
        scale_y_continuous(breaks=seq(0, 1, by = 0.2))+
        theme(title = element_text(size = 20),
              axis.title=element_text(size = 20),
              legend.title=element_text(size = 20),
              axis.text=element_text(size = 20),
              legend.text=element_text(size = 20))+
        labs(title = expression(Misspecification~","~pi[A] ==0.5))
      
      p_sensi5_power_05
      
    }
    
    # Power  -  06
    {
      
      
      Res <- readRDS("res_sensi5/res_06_power_sensi5.RData")
      
      Res <- bind_rows(Res)
      Res$n <- seq(10, 600, by = 5)
      
      Res_plot <-  Res %>% 
        pivot_longer(
          cols = c("power_EC", "power_AIPW", "power_naive"), 
          names_to = "Estimator",
          values_to = "power"
        )
      
      # add single-arm case
      Res_sa <- readRDS("res_sensi5/res_sa_power_sensi5.RData")
      
      Res_sa <- bind_rows(Res_sa)
      Res_sa$n <- seq(10, 600, by = 5)
      
      
      Res_plot_sa <-  Res_sa %>% 
        pivot_longer(
          cols = c("power_SA"), 
          names_to = "Estimator",
          values_to = "power"
        )
      
      Res_plot <- rbind(Res_plot, Res_plot_sa)
      
      Res_plot$Estimator <- factor(Res_plot$Estimator, 
                                   levels = c(
                                     "power_EC",
                                     "power_AIPW",
                                     "power_naive",
                                     "power_SA"
                                   ))
      levels(Res_plot$Estimator) <- c(
        "Hybrid", 
        "AIPW",
        "Diff-in-mean",
        "Single-arm")
      
      p_sensi5_power_06 <- ggplot(data = Res_plot, aes(x = n, y = power, colour = Estimator)) +
        geom_line(linewidth = 1) +
        scale_colour_manual(values=c("#F27127FF","#FDC718FF",
                                     "#038C33FF", "#207CF5FF"),
                            breaks=c("Diff-in-mean", "AIPW",
                                     "Hybrid", "Single-arm"))+
        geom_hline(yintercept = 0.8,
                   linetype="dashed",
                   color = "#084594",
                   linewidth = 1) +
        geom_vline(xintercept= 56,
                   linetype="dashed",
                   color = "#011C40FF",
                   linewidth = 1)+
        geom_vline(xintercept= 93,
                   linetype="dashed",
                   color =  "#005A32",
                   linewidth = 1)+
        geom_vline(xintercept= 215,
                   linetype="dashed",
                   color = "#FEB24C",
                   linewidth = 1)+
        theme(legend.position = "none")+
        labs(x=expression(Sample~Size~N[R]), y="Power")+
        scale_x_continuous(breaks=seq(0, 600, by=100))+
        scale_y_continuous(breaks=seq(0, 1, by = 0.2))+
        theme(title = element_text(size = 20),
              axis.title=element_text(size = 20),
              legend.title=element_text(size = 20),
              axis.text=element_text(size = 20),
              legend.text=element_text(size = 20))+
        labs(title = expression(Misspecification~","~pi[A] ==0.6))
      
      p_sensi5_power_06
      
    }
    
    # Power  -  07
    {
      
      Res <- readRDS("res_sensi5/res_07_power_sensi5.RData")
      
      Res <- bind_rows(Res)
      Res$n <- seq(10, 600, by = 5)
      
      Res_plot <-  Res %>% 
        pivot_longer(
          cols = c("power_EC", "power_AIPW", "power_naive"), 
          names_to = "Estimator",
          values_to = "power"
        )
      
      # add single-arm case
      Res_sa <- readRDS("res_sensi5/res_sa_power_sensi5.RData")
      
      Res_sa <- bind_rows(Res_sa)
      Res_sa$n <- seq(10, 600, by = 5)
      
      
      Res_plot_sa <-  Res_sa %>% 
        pivot_longer(
          cols = c("power_SA"), 
          names_to = "Estimator",
          values_to = "power"
        )
      
      Res_plot <- rbind(Res_plot, Res_plot_sa)
      
      Res_plot$Estimator <- factor(Res_plot$Estimator, 
                                   levels = c(
                                     "power_EC",
                                     "power_AIPW",
                                     "power_naive",
                                     "power_SA"
                                   ))
      levels(Res_plot$Estimator) <- c(
        "Hybrid", 
        "AIPW",
        "Diff-in-mean",
        "Single-arm")
      
      p_sensi5_power_07 <- ggplot(data = Res_plot, aes(x = n, y = power, colour = Estimator)) +
        geom_line(linewidth = 1) +
        scale_colour_manual(values=c("#F27127FF","#FDC718FF",
                                     "#038C33FF", "#207CF5FF"),
                            breaks=c("Diff-in-mean", "AIPW",
                                     "Hybrid", "Single-arm"))+
        geom_hline(yintercept = 0.8,
                   linetype="dashed",
                   color = "#084594",
                   linewidth = 1) +
        geom_vline(xintercept= 56,
                   linetype="dashed",
                   color = "#011C40FF",
                   linewidth = 1)+
        geom_vline(xintercept= 80,
                   linetype="dashed",
                   color =  "#005A32",
                   linewidth = 1)+
        geom_vline(xintercept= 246,
                   linetype="dashed",
                   color = "#FEB24C",
                   linewidth = 1)+
        theme(legend.position = "none")+
        labs(x=expression(Sample~Size~N[R]), y="Power")+
        scale_x_continuous(breaks=seq(0, 600, by=100))+
        scale_y_continuous(breaks=seq(0, 1, by = 0.2))+
        theme(title = element_text(size = 20),
              axis.title=element_text(size = 20),
              legend.title=element_text(size = 20),
              axis.text=element_text(size = 20),
              legend.text=element_text(size = 20))+
        labs(title = expression(Misspecification~","~pi[A] ==0.7))
      
      p_sensi5_power_07
      
    }
    
    # Power  -  08
    {
      
      Res <- readRDS("res_sensi5/res_08_power_sensi5.RData")
      
      Res <- bind_rows(Res)
      Res$n <- seq(10, 600, by = 5)
      
      Res_plot <-  Res %>% 
        pivot_longer(
          cols = c("power_EC", "power_AIPW", "power_naive"), 
          names_to = "Estimator",
          values_to = "power"
        )
      
      # add single-arm case
      Res_sa <- readRDS("res_sensi5/res_sa_power_sensi5.RData")
      
      Res_sa <- bind_rows(Res_sa)
      Res_sa$n <- seq(10, 600, by = 5)
      
      Res_plot_sa <-  Res_sa %>% 
        pivot_longer(
          cols = c("power_SA"), 
          names_to = "Estimator",
          values_to = "power"
        )
      
      Res_plot <- rbind(Res_plot, Res_plot_sa)
      
      Res_plot$Estimator <- factor(Res_plot$Estimator, 
                                   levels = c(
                                     "power_EC",
                                     "power_AIPW",
                                     "power_naive",
                                     "power_SA"
                                   ))
      levels(Res_plot$Estimator) <- c(
        "Hybrid", 
        "AIPW",
        "Diff-in-mean",
        "Single-arm")
      
      p_sensi5_power_08 <- ggplot(data = Res_plot, aes(x = n, y = power, colour = Estimator)) +
        geom_line(linewidth = 1) +
        scale_colour_manual(values=c("#F27127FF","#FDC718FF",
                                     "#038C33FF", "#207CF5FF"),
                            breaks=c("Diff-in-mean", "AIPW",
                                     "Hybrid", "Single-arm"))+
        geom_hline(yintercept = 0.8,
                   linetype="dashed",
                   color = "#084594",
                   linewidth = 1) +
        geom_vline(xintercept= 56,
                   linetype="dashed",
                   color = "#011C40FF",
                   linewidth = 1)+
        geom_vline(xintercept= 70,
                   linetype="dashed",
                   color =  "#005A32",
                   linewidth = 1)+
        geom_vline(xintercept= 323,
                   linetype="dashed",
                   color = "#FEB24C",
                   linewidth = 1)+
        theme(legend.position = "none")+
        labs(x=expression(Sample~Size~N[R]), y="Power")+
        scale_x_continuous(breaks=seq(0, 600, by=100))+
        scale_y_continuous(breaks=seq(0, 1, by = 0.2))+
        theme(title = element_text(size = 20),
              axis.title=element_text(size = 20),
              legend.title=element_text(size = 20),
              axis.text=element_text(size = 20),
              legend.text=element_text(size = 20))+
        labs(title = expression(Misspecification~","~pi[A] ==0.8))
      
      p_sensi5_power_08
      
    }
    
    # Power  -  09
    {
      
      Res <- readRDS("res_sensi5/res_09_power_sensi5.RData")
      
      Res <- bind_rows(Res)
      Res$n <- seq(10, 600, by = 5)
      
      Res_plot <-  Res %>% 
        pivot_longer(
          cols = c("power_EC", "power_AIPW", "power_naive"), 
          names_to = "Estimator",
          values_to = "power"
        )
      
      # add single-arm case
      Res_sa <- readRDS("res_sensi5/res_sa_power_sensi5.RData")
      
      Res_sa <- bind_rows(Res_sa)
      Res_sa$n <- seq(10, 600, by = 5)
      
      
      Res_plot_sa <-  Res_sa %>% 
        pivot_longer(
          cols = c("power_SA"), 
          names_to = "Estimator",
          values_to = "power"
        )
      
      Res_plot <- rbind(Res_plot, Res_plot_sa)
      
      Res_plot$Estimator <- factor(Res_plot$Estimator, 
                                   levels = c(
                                     "power_EC",
                                     "power_AIPW",
                                     "power_naive",
                                     "power_SA"
                                   ))
      levels(Res_plot$Estimator) <- c(
        "Hybrid", 
        "AIPW",
        "Diff-in-mean",
        "Single-arm")
      
      p_sensi5_power_09 <- ggplot(data = Res_plot, aes(x = n, y = power, colour = Estimator)) +
        geom_line(linewidth = 1) +
        scale_colour_manual(values=c("#F27127FF","#FDC718FF",
                                     "#038C33FF", "#207CF5FF"),
                            breaks=c("Diff-in-mean", "AIPW",
                                     "Hybrid", "Single-arm"))+
        geom_hline(yintercept = 0.8,
                   linetype="dashed",
                   color = "#084594",
                   linewidth = 1) +
        geom_vline(xintercept= 56,
                   linetype="dashed",
                   color = "#011C40FF",
                   linewidth = 1)+
        geom_vline(xintercept= 62,
                   linetype="dashed",
                   color =  "#005A32",
                   linewidth = 1)+
        geom_vline(xintercept= 573,
                   linetype="dashed",
                   color = "#FEB24C",
                   linewidth = 1)+
        # theme(legend.position = "none")+
        labs(x=expression(Sample~Size~N[R]), y="Power")+
        scale_x_continuous(breaks=seq(0, 600, by=100))+
        scale_y_continuous(breaks=seq(0, 1, by = 0.2))+
        theme(title = element_text(size = 20),
              axis.title=element_text(size = 20),
              legend.title=element_text(size = 20),
              axis.text=element_text(size = 20),
              legend.text=element_text(size = 20))+
        labs(title = expression(Misspecification~","~pi[A] ==0.9))
      
      p_sensi5_power_09
      
    }
    
  }
  
  ## Type I error sufficient
  {
    
    # Type I error - 05
    {
      
      Res <- readRDS("res_sensi5/res_05_TypeI_sensi5.RData")
      
      Res <- bind_rows(Res)
      Res$n <- seq(10, 600, by = 5)
      
      Res_plot <-  Res %>% 
        pivot_longer(
          cols = c("power_EC", "power_AIPW", "power_naive"), 
          names_to = "Estimator",
          values_to = "Type_I_Error"
        )
      
      # add single-arm case
      Res_sa <- readRDS("res_sensi5/res_sa_TypeI_sensi5.RData")
      
      Res_sa <- bind_rows(Res_sa)
      Res_sa$n <- seq(10, 600, by = 5)
      
      
      Res_plot_sa <-  Res_sa %>% 
        pivot_longer(
          cols = c("power_SA"), 
          names_to = "Estimator",
          values_to = "Type_I_Error"
        )
      
      Res_plot <- rbind(Res_plot, Res_plot_sa)
      
      Res_plot$Estimator <- factor(Res_plot$Estimator, 
                                   levels = c(
                                     "power_EC",
                                     "power_AIPW",
                                     "power_naive",
                                     "power_SA"
                                   ))
      levels(Res_plot$Estimator) <- c(
        "Hybrid", 
        "AIPW",
        "Diff-in-mean",
        "Single-arm")
      
      p_sensi5_TypeI_05 <- ggplot(data = Res_plot, aes(x = n, y = Type_I_Error, colour = Estimator)) +
        geom_line(linewidth = 1) +
        scale_colour_manual(values=c("#F27127FF","#FDC718FF",
                                     "#038C33FF", "#207CF5FF"), 
                            breaks=c("Diff-in-mean", "AIPW",
                                     "Hybrid", "Single-arm"))+
        geom_hline(yintercept = 0.05,
                   linetype="dashed",
                   color = "#084594",
                   linewidth = 1) +
        geom_vline(xintercept= 56,
                   linetype="dashed",
                   color = "#011C40FF",
                   linewidth = 1)+
        geom_vline(xintercept= 111,
                   linetype="dashed",
                   color =  "#005A32",
                   linewidth = 1)+
        geom_vline(xintercept= 207,
                   linetype="dashed",
                   color = "#FEB24C",
                   linewidth = 1)+
        theme(legend.position = "none")+
        labs(x=expression(Sample~Size~N[R]), y="Type I Error")+
        scale_x_continuous(breaks=seq(0, 600, by=100))+
        scale_y_continuous(breaks = seq(0, 0.546, by=0.1))+
        expand_limits(y=c(0, 0.546))+
        theme(title = element_text(size = 20),
              axis.title=element_text(size = 20),
              legend.title=element_text(size = 20),
              axis.text=element_text(size = 20),
              legend.text=element_text(size = 20))+
        labs(title = expression(Misspecification~","~pi[A] ==0.5))
      
      p_sensi5_TypeI_05
      
    }
    
    # Type I error - 06
    {
      
      Res <- readRDS("res_sensi5/res_06_TypeI_sensi5.RData")
      
      Res <- bind_rows(Res)
      Res$n <- seq(10, 600, by = 5)
      
      Res_plot <-  Res %>% 
        pivot_longer(
          cols = c("power_EC", "power_AIPW", "power_naive"), 
          names_to = "Estimator",
          values_to = "Type_I_Error"
        )
      
      # add single-arm case
      Res_sa <- readRDS("res_sensi5/res_sa_TypeI_sensi5.RData")
      
      Res_sa <- bind_rows(Res_sa)
      Res_sa$n <- seq(10, 600, by = 5)
      
      Res_plot_sa <-  Res_sa %>% 
        pivot_longer(
          cols = c("power_SA"), 
          names_to = "Estimator",
          values_to = "Type_I_Error"
        )
      
      Res_plot <- rbind(Res_plot, Res_plot_sa)
      
      Res_plot$Estimator <- factor(Res_plot$Estimator, 
                                   levels = c(
                                     "power_EC",
                                     "power_AIPW",
                                     "power_naive",
                                     "power_SA"
                                   ))
      levels(Res_plot$Estimator) <- c(
        "Hybrid", 
        "AIPW",
        "Diff-in-mean",
        "Single-arm")
      
      
      p_sensi5_TypeI_06 <- ggplot(data = Res_plot, aes(x = n, y = Type_I_Error, colour = Estimator)) +
        geom_line(linewidth = 1) +
        scale_colour_manual(values=c("#F27127FF","#FDC718FF",
                                     "#038C33FF", "#207CF5FF"), 
                            breaks=c("Diff-in-mean", "AIPW",
                                     "Hybrid", "Single-arm"))+
        geom_hline(yintercept = 0.05,
                   linetype="dashed",
                   color = "#084594",
                   linewidth = 1) +
        geom_vline(xintercept= 56,
                   linetype="dashed",
                   color = "#011C40FF",
                   linewidth = 1)+
        geom_vline(xintercept= 93,
                   linetype="dashed",
                   color =  "#005A32",
                   linewidth = 1)+
        geom_vline(xintercept= 215,
                   linetype="dashed",
                   color = "#FEB24C",
                   linewidth = 1)+
        theme(legend.position = "none")+
        labs(x=expression(Sample~Size~N[R]), y="Type I Error")+
        scale_x_continuous(breaks=seq(0, 600, by=100))+
        scale_y_continuous(breaks = seq(0, 0.546, by=0.1))+
        expand_limits(y=c(0, 0.546))+
        theme(title = element_text(size = 20),
              axis.title=element_text(size = 20),
              legend.title=element_text(size = 20),
              axis.text=element_text(size = 20),
              legend.text=element_text(size = 20))+
        labs(title = expression(Misspecification~","~pi[A] ==0.6))
      
      p_sensi5_TypeI_06
      
    }
    
    # Type I error - 07
    {
      Res <- readRDS("res_sensi5/res_07_TypeI_sensi5.RData")
      
      Res <- bind_rows(Res)
      Res$n <- seq(10, 600, by = 5)
      
      Res_plot <-  Res %>% 
        pivot_longer(
          cols = c("power_EC", "power_AIPW", "power_naive"), 
          names_to = "Estimator",
          values_to = "Type_I_Error"
        )
      
      # add single-arm case
      Res_sa <- readRDS("res_sensi5/res_sa_TypeI_sensi5.RData")
      
      Res_sa <- bind_rows(Res_sa)
      Res_sa$n <- seq(10, 600, by = 5) 
      
      Res_plot_sa <-  Res_sa %>% 
        pivot_longer(
          cols = c("power_SA"), 
          names_to = "Estimator",
          values_to = "Type_I_Error"
        )
      
      Res_plot <- rbind(Res_plot, Res_plot_sa)
      
      Res_plot$Estimator <- factor(Res_plot$Estimator, 
                                   levels = c(
                                     "power_EC",
                                     "power_AIPW",
                                     "power_naive",
                                     "power_SA"
                                   ))
      levels(Res_plot$Estimator) <- c(
        "Hybrid", 
        "AIPW",
        "Diff-in-mean",
        "Single-arm")
      
      
      p_sensi5_TypeI_07 <- ggplot(data = Res_plot, aes(x = n, y = Type_I_Error, colour = Estimator)) +
        geom_line(linewidth = 1) +
        scale_colour_manual(values=c("#F27127FF","#FDC718FF",
                                     "#038C33FF", "#207CF5FF"), 
                            breaks=c("Diff-in-mean", "AIPW",
                                     "Hybrid", "Single-arm"))+
        geom_hline(yintercept = 0.05,
                   linetype="dashed",
                   color = "#084594",
                   linewidth = 1) +
        geom_vline(xintercept= 56,
                   linetype="dashed",
                   color = "#011C40FF",
                   linewidth = 1)+
        geom_vline(xintercept= 80,
                   linetype="dashed",
                   color =  "#005A32",
                   linewidth = 1)+
        geom_vline(xintercept= 246,
                   linetype="dashed",
                   color = "#FEB24C",
                   linewidth = 1)+
        theme(legend.position = "none")+
        labs(x=expression(Sample~Size~N[R]), y="Type I Error")+
        scale_x_continuous(breaks=seq(0, 600, by=100))+
        scale_y_continuous(breaks = seq(0, 0.546, by=0.1))+
        expand_limits(y=c(0, 0.546))+
        theme(title = element_text(size = 20),
              axis.title=element_text(size = 20),
              legend.title=element_text(size = 20),
              axis.text=element_text(size = 20),
              legend.text=element_text(size = 20))+
        labs(title = expression(Misspecification~","~pi[A] ==0.7))
      
      p_sensi5_TypeI_07
      
    }
    
    # Type I error - 08
    {
      Res <- readRDS("res_sensi5/res_08_TypeI_sensi5.RData")
      
      Res <- bind_rows(Res)
      Res$n <- seq(10, 600, by = 5)
      
      Res_plot <-  Res %>% 
        pivot_longer(
          cols = c("power_EC", "power_AIPW", "power_naive"), 
          names_to = "Estimator",
          values_to = "Type_I_Error"
        )
      
      # add single-arm case
      Res_sa <- readRDS("res_sensi5/res_sa_TypeI_sensi5.RData")
      
      Res_sa <- bind_rows(Res_sa)
      Res_sa$n <- seq(10, 600, by = 5)
      
      Res_plot_sa <-  Res_sa %>% 
        pivot_longer(
          cols = c("power_SA"), 
          names_to = "Estimator",
          values_to = "Type_I_Error"
        )
      
      Res_plot <- rbind(Res_plot, Res_plot_sa)
      
      Res_plot$Estimator <- factor(Res_plot$Estimator, 
                                   levels = c(
                                     "power_EC",
                                     "power_AIPW",
                                     "power_naive",
                                     "power_SA"
                                   ))
      levels(Res_plot$Estimator) <- c(
        "Hybrid", 
        "AIPW",
        "Diff-in-mean",
        "Single-arm")
      
      
      p_sensi5_TypeI_08 <- ggplot(data = Res_plot, aes(x = n, y = Type_I_Error, colour = Estimator)) +
        geom_line(linewidth = 1) +
        scale_colour_manual(values=c("#F27127FF","#FDC718FF",
                                     "#038C33FF", "#207CF5FF"), 
                            breaks=c("Diff-in-mean", "AIPW",
                                     "Hybrid", "Single-arm"))+
        geom_hline(yintercept = 0.05,
                   linetype="dashed",
                   color = "#084594",
                   linewidth = 1) +
        geom_vline(xintercept= 56,
                   linetype="dashed",
                   color = "#011C40FF",
                   linewidth = 1)+
        geom_vline(xintercept= 70,
                   linetype="dashed",
                   color =  "#005A32",
                   linewidth = 1)+
        geom_vline(xintercept= 323,
                   linetype="dashed",
                   color = "#FEB24C",
                   linewidth = 1)+
        theme(legend.position = "none")+
        labs(x=expression(Sample~Size~N[R]), y="Type I Error")+
        scale_x_continuous(breaks=seq(0, 600, by=100))+
        scale_y_continuous(breaks = seq(0, 0.546, by=0.1))+
        expand_limits(y=c(0, 0.546))+
        theme(title = element_text(size = 20),
              axis.title=element_text(size = 20),
              legend.title=element_text(size = 20),
              axis.text=element_text(size = 20),
              legend.text=element_text(size = 20))+
        labs(title = expression(Misspecification~","~pi[A] ==0.8))
      
      p_sensi5_TypeI_08
      
    }
    
    # Type I error - 09
    {
      
      Res <- readRDS("res_sensi5/res_09_TypeI_sensi5.RData")
      
      Res <- bind_rows(Res)
      Res$n <- seq(10, 600, by = 5)
      
      Res_plot <-  Res %>% 
        pivot_longer(
          cols = c("power_EC", "power_AIPW", "power_naive"), 
          names_to = "Estimator",
          values_to = "Type_I_Error"
        )
      
      # add single-arm case
      Res_sa <- readRDS("res_sensi5/res_sa_TypeI_sensi5.RData")
      
      Res_sa <- bind_rows(Res_sa)
      Res_sa$n <- seq(10, 600, by = 5)
      
      Res_plot_sa <-  Res_sa %>% 
        pivot_longer(
          cols = c("power_SA"), 
          names_to = "Estimator",
          values_to = "Type_I_Error"
        )
      
      Res_plot <- rbind(Res_plot, Res_plot_sa)
      
      Res_plot$Estimator <- factor(Res_plot$Estimator, 
                                   levels = c(
                                     "power_EC",
                                     "power_AIPW",
                                     "power_naive",
                                     "power_SA"
                                   ))
      levels(Res_plot$Estimator) <- c(
        "Hybrid", 
        "AIPW",
        "Diff-in-mean",
        "Single-arm")
      
      
      p_sensi5_TypeI_09 <- ggplot(data = Res_plot, aes(x = n, y = Type_I_Error, colour = Estimator)) +
        geom_line(linewidth = 1) +
        scale_colour_manual(values=c("#F27127FF","#FDC718FF",
                                     "#038C33FF", "#207CF5FF"), 
                            breaks=c("Diff-in-mean", "AIPW",
                                     "Hybrid", "Single-arm"))+
        geom_hline(yintercept = 0.05,
                   linetype="dashed",
                   color = "#084594",
                   linewidth = 1) +
        geom_vline(xintercept= 56,
                   linetype="dashed",
                   color = "#011C40FF",
                   linewidth = 1)+
        geom_vline(xintercept= 62,
                   linetype="dashed",
                   color =  "#005A32",
                   linewidth = 1)+
        geom_vline(xintercept= 573,
                   linetype="dashed",
                   color = "#FEB24C",
                   linewidth = 1)+
        # theme(legend.position = "none")+
        labs(x=expression(Sample~Size~N[R]), y="Type I Error")+
        scale_x_continuous(breaks=seq(0, 600, by=100))+
        scale_y_continuous(breaks = seq(0, 0.546, by=0.1))+
        expand_limits(y=c(0, 0.546))+
        theme(title = element_text(size = 20),
              axis.title=element_text(size = 20),
              legend.title=element_text(size = 20),
              axis.text=element_text(size = 20),
              legend.text=element_text(size = 20))+
        labs(title = expression(Misspecification~","~pi[A] ==0.9))
      
      p_sensi5_TypeI_09
      
    }
    
  }
  
}



###  For sensitivity analysis
{
  
  # sensi1 
  {
    p_sensi1 <- ggarrange(p_sensi1_power_05, p_sensi1_power_06, p_sensi1_power_07, 
                          p_sensi1_power_08, p_sensi1_power_09,
                          p_sensi1_TypeI_05, p_sensi1_TypeI_06, p_sensi1_TypeI_07, 
                          p_sensi1_TypeI_08, p_sensi1_TypeI_09,
                          widths = c(1, 1, 1, 1, 1.38),
                          ncol = 5, nrow = 2)
    p_sensi1
    p_sensi1_new <- annotate_figure(p_sensi1,
                                    bottom = text_grob("(a) Heterogeneous treatment effect across the current study population", face = "bold", size = 25)
    )
    p_sensi1_new
    
    ggsave(plot = p_sensi1,
           width = 37.5,
           height = 12,
           dpi = 600,
           filename = "/code/simulation/plots/p_sensi1.pdf")
    
    
  }
  
  # sensi2
  {
    p_sensi2 <- ggarrange(p_sensi2_power_05, p_sensi2_power_06, p_sensi2_power_07, 
                          p_sensi2_power_08, p_sensi2_power_09,
                          p_sensi2_TypeI_05, p_sensi2_TypeI_06, p_sensi2_TypeI_07, 
                          p_sensi2_TypeI_08, p_sensi2_TypeI_09,
                          widths = c(1, 1, 1, 1, 1.38),
                          ncol = 5, nrow = 2)
    p_sensi2
    p_sensi2_new <- annotate_figure(p_sensi2,
                                    bottom = text_grob("(b) Non-constant conditional variance of the outcome", face = "bold", size = 25)
    )
    p_sensi2_new
   
    ggsave(plot = p_sensi2,
           width = 37.5,
           height = 12,
           dpi = 600,
           filename = "/code/simulation/plots/p_sensi2.pdf")
    
    
  }
  
  
  p_sensi <- ggarrange(p_sensi1_new,
                       NULL,   # this row will be pure whitespace
                       p_sensi2_new,
                       heights = c(1, 0.06, 1),
                       ncol = 1, nrow = 3)
  p_sensi
  ggsave(plot = p_sensi,
         width = 37.5,
         height = 26,
         dpi = 1000,
         filename = "/code/simulation/plots/p_sensi.pdf")
  
  
  # sensi3
  {
    p_sensi3_1 <- ggarrange(p_sensi3_power_05, p_sensi3_power_06, p_sensi3_power_07, 
                            p_sensi3_power_08, p_sensi3_power_09,
                            p_sensi3_TypeI_05, p_sensi3_TypeI_06, p_sensi3_TypeI_07, 
                            p_sensi3_TypeI_08, p_sensi3_TypeI_09,
                            widths = c(1, 1, 1, 1, 1.38),
                            ncol = 5, nrow = 2)
    p_sensi3_1
    p_sensi3_new1 <- annotate_figure(p_sensi3_1,
                                     bottom = text_grob("(a) Violation of compatibility (delta_U = 0)", face = "bold", size = 25)
    )
    p_sensi3_new1
  
    # sensi3
    p_sensi3_2 <- ggarrange(p_sensi3_power_05, p_sensi3_power_06, p_sensi3_power_07, 
                            p_sensi3_power_08, p_sensi3_power_09,
                            p_sensi3_TypeI_05, p_sensi3_TypeI_06, p_sensi3_TypeI_07, 
                            p_sensi3_TypeI_08, p_sensi3_TypeI_09,
                            widths = c(1, 1, 1, 1, 1.38),
                            ncol = 5, nrow = 2)
    p_sensi3_2
    p_sensi3_new2 <- annotate_figure(p_sensi3_2,
                                     bottom = text_grob("(b) Violation of compatibility (delta_U = 0.5)", face = "bold", size = 25)
    )
    p_sensi3_new2
    
    p_sensi3 <- ggarrange(p_sensi3_new1,
                          NULL,   # this row will be pure whitespace
                          p_sensi3_new2,
                          heights = c(1, 0.06, 1),
                          ncol = 1, nrow = 3)
    p_sensi3
    ggsave(plot = p_sensi3,
           width = 37.5,
           height = 26,
           dpi = 1000,
           filename = "plots/p_sensi3.pdf")
    
  }
  
  
  # sensi5
  {
    p_sensi5 <- ggarrange(p_sensi5_power_05, p_sensi5_power_06, p_sensi5_power_07, 
                          p_sensi5_power_08, p_sensi5_power_09,
                          p_sensi5_TypeI_05, p_sensi5_TypeI_06, p_sensi5_TypeI_07, 
                          p_sensi5_TypeI_08, p_sensi5_TypeI_09,
                          widths = c(1, 1, 1, 1, 1.38),
                          ncol = 5, nrow = 2)
    p_sensi5
    
    ggsave(plot = p_sensi5,
           width = 37.5,
           height = 12,
           dpi = 600,
           filename = "plots/p_sensi5.pdf")
    
  }
  
}






