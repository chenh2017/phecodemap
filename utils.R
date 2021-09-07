
# color for ICD9s, ICD10s, highlightNode, Phecode
colorlist <- c("lightblue", "palegreen", "red",
               brewer.pal(9, "Set1")[c(5, 4, 6, 7, 2, 3, 8)],
               brewer.pal(9, "Set3")[c(2:4,6:12)],
               brewer.pal(9, "Set1")[9])

buildPath <- function(rootid, icdmap) {

  # rootid <- "296"
  rootid <- paste0("Phe_", rootid)
  s_map <- icdmap
  s_map$Phecode <- paste0("Phe_", s_map$Phecode)
  s_map <- s_map[grepl(paste0(rootid, "\\..+"), s_map$Phecode, perl = TRUE) | (s_map$Phecode == rootid), ]

  print(paste("s_map:", nrow(s_map)))

  node_phe <- s_map[, c("Phecode", "Phenotype")]
  node_phe <- node_phe[!duplicated(node_phe), ]
  node_phe$pathString <- sapply(node_phe$Phecode, getPath, 3)

  node_icd <- s_map[, c("ICD_version", "ICD_id", "ICD_str")]
  node_icd <- node_icd[!duplicated(node_icd), ]
  node_icd$pathString <- sapply(node_icd$ICD_id, getPath)

  s_map <- dplyr::left_join(s_map, node_phe[, c(1, 3)], by = "Phecode")
  s_map <- dplyr::left_join(s_map, node_icd[, c(2, 4)], by = "ICD_id")
  s_map$pathString <- paste(s_map$pathString.x, s_map$ICD_version,s_map$pathString.y, sep = "/")

  node <- data.frame(ids = c(), labels = c(), parents = c())
  for (id in s_map$pathString) {
    while (grepl("/", id, fixed = TRUE)) {
      all_n <- strsplit(id, "/", fixed = TRUE)[[1]]
      label <- all_n[length(all_n)]
      parent <- substring(id, 1, nchar(id) - nchar(label) - 1)
      node <- rbind(node, data.frame(ids = id, labels = label, parents = parent))
      id <- parent
    }
  }

  node <- node[!duplicated(node), ]
  node <- rbind(node, data.frame(ids = rootid, labels = rootid, parents = ""))
  node <- left_join(left_join(node, s_map[, 2:3], by = c("labels" = "ICD_id")),
    s_map[, 4:5],
    by = c("labels" = "Phecode")
  )
  node <- node[!duplicated(node), ]
  node[is.na(node)] <- ""
  node$strs <- paste0(node$ICD_str, node$Phenotype)
  node$strs[node$strs == ""] <- node$labels[node$strs == ""]
  return(node)
}

addClass <- function(rootid, icdmap, df_highlight = df_highlight){
  node <- buildPath(rootid, icdmap)
  nodes_color <- node[order(node$ids), c("ids", "labels")]
  nodes_color$newid <- nodes_color$labels

  node_phe <- unique(node$labels[grepl("Phe", node$labels)])

  nodes1 <- nodes_color[(!duplicated(nodes_color$labels)) | (!nodes_color$labels %in% node_phe),]
  # dup nodes in phecode
  nodes2 <- nodes_color[duplicated(nodes_color$labels) & (nodes_color$labels %in% node_phe),]
  nodes2$newid <- sapply(1:nrow(nodes2), addBlank)

  dup_node <- unique(nodes1$newid[duplicated(nodes1$newid)])
  if (length(dup_node) != 0) {
    for (i in dup_node) {
      len <- 0:(sum(nodes1$newid %in% i)-1)
      nodes1$newid[nodes1$newid %in% i] <- paste0(sapply(len, addBlank), i)
    }
  }
  nodes2$class <- "dupnode"
  nodes1$class <- sapply(nodes1$ids, classifyNode)
  # add highlight class

  df_highlight$Phecode <- paste0("Phe_", df_highlight$Phecode)
  highligh_nodes <- intersect(node_phe, df_highlight$Phecode)

  if (length(highligh_nodes) > 0){
    for (i in highligh_nodes){
      nodes1$class[nodes1$labels == i] <- paste0(i, "(only ", df_highlight$ICD_version[df_highlight$Phecode == i], ")")
    }
  }

  return(list(node, nodes1, nodes2))
}



addColor <- function(nodes_list, plot = "tree"){

  node <- nodes_list[[1]]
  nodes1 <- nodes_list[[2]]
  nodes2 <- nodes_list[[3]]

  df_color <- data.frame(class = unique(nodes1$class), color = "")
  df_color$color[df_color$class == "ICD-9"] <- colorlist[1]
  df_color$color[df_color$class == "ICD-10-cm"] <- colorlist[2]
  df_color$color[!df_color$class %in% c("ICD-9", "ICD-10-cm")] <- colorlist[4:(nrow(df_color[!df_color$class %in% c("ICD-9", "ICD-10-cm"),])+3)]

  df_color$color[grepl("only", df_color$class)] <- colorlist[3]
  nodes1 <- left_join(nodes1, df_color, by = "class")
  if (plot == "tree"){
    nodes2$color <- "none"
  } else{
    nodes2 <- left_join(nodes2, nodes1[, c("labels", "color")], by = "labels", keep = FALSE)
  }
  return(list(node, nodes1, nodes2))
}


dfPlot <- function(nodes_list, plot = "tree"){

  nodes_list <- addColor(nodes_list, plot)

  node <- nodes_list[[1]]
  nodes1 <- nodes_list[[2]]
  nodes2 <- nodes_list[[3]]
  nodes_color <- rbind(nodes1, nodes2)

  df_plot <- left_join(node, nodes_color[, -2], by = "ids")
  df_plot <- left_join(df_plot, nodes_color[, c("ids", "newid")], by = c("parents" = "ids"))
  df_plot <- df_plot[, c(10, 7, 1, 3, 2, 6, 8, 9)]
  colnames(df_plot) <- c("parents", "ids", "nodepath", "parentpath", "labels", "info", "class","color")
  return(df_plot)
}

dfSunburst <- function(nodes_list){
  df_plot <- dfPlot(nodes_list, plot = "sunb")
  df_plot$linecolor <- df_plot$color
  dup_nodes <- df_plot$labels[grep("Phe", df_plot$labels)]
  dup_nodes <- unique(dup_nodes[duplicated(dup_nodes)])
  df_plot$linecolor[!df_plot$labels %in% dup_nodes] <- "white"
  # df_plot$linecolor[df_plot$color == "red"] <- "black"
  # df_plot$linewidth <- 1
  # df_plot$linewidth[df_plot$class != "dupnode"] <- 1
  # df_plot$linecolor[df_plot$color == "red"] <- 3
  df_plot$labels[df_plot$class == "dupnode"] <- ""
  return(df_plot)
}

sunburstPlotly <- function(df_plot, maxd = 10) {
  plot_ly(df_plot,
          ids =~ ids,
          labels =~ labels,
          parents =~ parents,
          hovertext =~ info,
          maxdepth = maxd,
          marker = list(colors = df_plot$color,
                        line = list(color = df_plot$linecolor #,
                                    # width = df_plot$linewidth
                                    )
                        ),
          type = "sunburst"
  )
}

treePlot <- function(nodes_list, maxd = 4, collapsed = FALSE) {
  df_plot <- dfPlot(nodes_list, plot = "tree")
  df_plot <- df_plot[sapply(df_plot$nodepath, filterNode, maxd), ]
  collapsibleTreeNetwork(df_plot,
                         attribute = "labels", fill = "color",
                         collapsed = collapsed, tooltip = TRUE,
                         tooltipHtml = "info", width = "100%", height = 1000
  )
}

legends <- function(df_sunb){
  df <- delDupRow(df_sunb[, c("class", "color")])
  df <- df[df$class != "dupnode",]
  df <- df[order(df$class),]

  plot(NULL ,xaxt='n',yaxt='n',bty='n',ylab='',xlab='', xlim=0:1, ylim=0:1)
  legend("topleft", legend =df$class, pch=16, pt.cex=3, cex=1.5, bty='n',
         y.intersp=1.3,
         col = df$color)
}

delDupRow <- function(x, cols = 1:ncol(x)) {
  x <- x[, cols]
  x[!duplicated(x), ]
}

getParents <- function(node){
  parents <- strsplit(node, ".", fixed = TRUE)[[1]][1]
  children <- strsplit(node, ".", fixed = TRUE)[[1]][2]
  if (!is.na(children)){
    for (i in strsplit(children, "")[[1]]){
      p <- parents[length(parents)]
      p <- ifelse(grepl(".", p, fixed = TRUE), p, paste0(p, "."))
      parents <- c(parents, paste0(p, i))
    }
  }
  parents
}

getPath <- function(node, depth = NULL){
  parents <- getParents(node)
  path <- paste(parents, collapse = "/")
  if (!is.null(depth)){
    p <- parents[length(parents)]
    dup_node <- rep(p, depth - length(parents))
    path <- purrr::reduce(c(path, dup_node), function(a, b){paste(a, b, sep = "/")})
  }
  path
}

filterNode <- function(x, maxd, pattern = "/") {
  return(length(strsplit(x, pattern, fixed = TRUE)[[1]]) <= maxd)
}

# if ICD in the path, class -> ICD_version
# 如果是ICD_version子节点，class全部定义为ICD_version
classifyNode <- function(x) {
  all_path <- strsplit(x, "/", fixed = TRUE)[[1]]
  if (grepl("CD", x, fixed = TRUE)) {
    return(all_path[grepl("CD", all_path, fixed = TRUE)])
  } else {
    return(all_path[length(all_path)])
  }
}

addBlank <- function(x){
  paste(rep(" ", times = x), collapse = "")
}

