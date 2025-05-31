library(ggplot2)
library(gground)
library(ggprism)
library(tidyverse)
library(dplyr)
library(tibble)
# 第一次安装gground时使用
#if (!require("devtools", quietly = TRUE))
#  install.packages("devtools")
#devtools::install_github("dxsbiocc/gground")
data<-read.table("plot.txt",header=T,sep='\t')

#可以直接读取GO富集分析结果
GO <- data %>% filter(ONTOLOGY %in% c("BP", "CC", "MF"))
#可以直接读取KEGG富集分析结果
KEGG <- data %>% filter(ONTOLOGY == "KEGG")

#按照要求提取数据
data <- group_by(GO, ONTOLOGY) %>%
  top_n(5, wt = -qvalue) %>%
  group_by(qvalue) %>% #按照qvalue提取GO
  top_n(1, wt = Count) %>%  #如果qvalue值相同，则会提取counts值最大的一行
  rbind(
    top_n(KEGG, 5, -qvalue) %>%
      group_by(qvalue) %>% #按照qvalue提取KEGG
      top_n(1, wt = Count) %>% #如果qvalue值相同，则会提取counts值最大的一行
      mutate(ONTOLOGY = 'KEGG')
  ) %>%
  ungroup() %>%
  mutate(ONTOLOGY = factor(ONTOLOGY, 
                           levels = rev(c('BP', 'CC', 'MF', 'KEGG')))) %>%
  arrange(ONTOLOGY, qvalue) %>%
  mutate(Description = factor(Description, levels = Description)) %>%
  tibble::rowid_to_column('index')

width <-0.75
# x 轴长度
xaxis_max <- max(-log10(data$qvalue)) +1
# 左侧分类标签数据
rect.data <- group_by(data, ONTOLOGY)%>%
  reframe(n = n())%>%
  ungroup()%>%
  mutate(
    xmin = -4.5* width,
    xmax = -2.5* width,
    ymax = cumsum(n),
    ymin = lag(ymax, default =0) +0.6,
    ymax = ymax +0.4
  )
#设置颜色
pal <- c('#7bc4e2','#acd372','#fbb05b','#ed6ca4')

#画图
a<-ggplot(data,aes(-log10(qvalue), y = index, fill = ONTOLOGY)) +
  geom_round_col(
    aes(y = Description), width =0.6, alpha =0.8
  ) +
  geom_text(
    aes(x =0.05, label = Description),
    hjust =0, size =5
  ) +
  geom_text(
    aes(x =0.1, label = ID, colour = ONTOLOGY),
    hjust =0, vjust =2.6, size =3.5, fontface ='italic',
    show.legend =FALSE
  )+
  # 基因数量
  geom_point(
    aes(x = -width, size = Count),
    shape =21
  ) +
  geom_text(
    aes(x = -width, label = Count)
  ) +
  scale_size_continuous(name ='Count', range = c(2,10)) +
  # 分类标签
  geom_round_rect(
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax,
        fill = ONTOLOGY),
    data = rect.data,
    radius = unit(2,'mm'),
    inherit.aes =FALSE
  ) +
  geom_text(
    aes(x = (xmin + xmax) /2, y = (ymin + ymax) /2, label = ONTOLOGY),
    data = rect.data,
    inherit.aes =FALSE
  )+
  annotate(
    "segment",
    x = 0, y = 0, xend = xaxis_max, yend = 0,
    linewidth = 1.5, 
    color = "black"
  )+
  labs(y =NULL) +
  scale_fill_manual(name ='Category', values = pal) +
  scale_colour_manual(values = pal) +
  scale_x_continuous(
    breaks = seq(0, xaxis_max,2),
    expand = expansion(c(0,0))
  ) +
  theme_prism() +
  theme(
    axis.text.y = element_blank(),
    axis.line = element_blank(),
    axis.ticks.y = element_blank(),
    legend.title = element_text()
  )

pdf('enrichment.pdf', width =10, height =8)
a
dev.off()

