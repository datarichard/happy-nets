traceNetwork <- function (ARTIVAnet, edgesThreshold, parentColor = "green", targetColor = "grey", 
                          parentgeneNames = TRUE, targetgeneNames = TRUE, layout = "fruchterman.reingold", 
                          onepage = TRUE) 
{
  if (edgesThreshold < min(ARTIVAnet$edgesThreshold)) {
    print(paste("WARNING : The coefficients for edges with posterior probability below", 
                min(ARTIVAnet$edgesThreshold), "were not estimated (grey edges) in the network to be plotted."))
  }
  targetGeneList = as.character(unique(ARTIVAnet$Target))
  parentGeneList = as.character(unique(ARTIVAnet$Parent))
  GeneNumber = length(targetGeneList) + length(parentGeneList)
  ARTIVAnet$PostProb[which(is.na(ARTIVAnet$PostProb))] = 0
  PotentialEdgesList = cbind(as.character(ARTIVAnet$Parent), 
                             as.character(ARTIVAnet$Target))
  GlobalNetwork = graph.edgelist(PotentialEdgesList)
  GlobalNetwork = set.edge.attribute(GlobalNetwork, "weight", 
                                     value = ARTIVAnet$CoeffMean)
  EdgeColorVec = rep("black", length(ARTIVAnet$CoeffMean))
  E(GlobalNetwork)$color = EdgeColorVec
  EdgeTypeVec = rep(1, length(ARTIVAnet$CoeffMean))
  TypeInduction = 1
  TypeRepression = 2
  EdgeTypeVec[ARTIVAnet$CoeffMean > 0] = TypeInduction
  EdgeTypeVec[ARTIVAnet$CoeffMean < 0] = TypeRepression
  E(GlobalNetwork)$lty = EdgeTypeVec
  EdgesToDelete = which(ARTIVAnet$PostProb < edgesThreshold)
  SubNetwork = delete.edges(GlobalNetwork, E(GlobalNetwork)[EdgesToDelete])# - 1]) # !!!!!
  NodeLabel = get.vertex.attribute(GlobalNetwork, name = "name")
  intervalParent = 0
  intervalTarget = 0
  SizeParent = 40 #20, 50
  SizeTarget = 40 #10, 30
  NodeCoord = NULL
  if (layout == "geneLines") {
    intervalParent = 2/(length(parentGeneList) - 1)
    intervalTarget = 2/(length(targetGeneList) - 1)
    xParent = -1
    xTarget = -1
    for (i in 1:length(NodeLabel)) {
      if (sum(parentGeneList == NodeLabel[i]) > 0) {
        if (length(parentGeneList) == 1) 
          NodeCoord = rbind(NodeCoord, c(0, 0.8))
        else {
          NodeCoord = rbind(NodeCoord, c(xParent, 0.8))
          xParent = xParent + intervalParent
        }
      }
      else if (sum(targetGeneList == NodeLabel[i]) > 0) {
        if (length(targetGeneList) == 1) {
          NodeCoord = rbind(NodeCoord, c(0, 0.2))
        }
        else {
          NodeCoord = rbind(NodeCoord, c(xTarget, 0.2))
          xTarget = xTarget + intervalTarget
        }
      }
    }
    if (intervalParent < 0.2) {
      SizeParent = 210/length(parentGeneList)
    }
    if (intervalTarget < 0.1) {
      SizeTarget = 200/length(targetGeneList)
    }
  }
  if (layout == "fruchterman.reingold") {
    NodeCoord = layout.fruchterman.reingold(SubNetwork)
  }
  if (layout == "random") {
    NodeCoord = layout.random(SubNetwork)
  }
  if (layout == "circle") {
    NodeCoord = layout.circle(SubNetwork)
  }
  if (layout == "kamada.kawai") {
    NodeCoord = layout.kamada.kawai(SubNetwork)
  }
  if (layout == "spring") {
    NodeCoord = layout.spring(SubNetwork)
  }
  if (layout == "reingold.tilford") {
    NodeCoord = layout.reingold.tilford(SubNetwork)
  }
  if (layout == "lgl") {
    NodeCoord = layout.lgl(SubNetwork)
  }
  if (layout == "graphopt") {
    NodeCoord = layout.graphopt(SubNetwork)
  }
  if (layout == "mds") {
    NodeCoord = layout.mds(SubNetwork)
  }
  if (layout == "svd") {
    NodeCoord = layout.svd(SubNetwork)
  }
  NodeLabel = get.vertex.attribute(GlobalNetwork, name = "name")
  NodeSize = NULL
  NodeColor = NULL
  WritenNodeLabel = NodeLabel
  for (i in 1:length(NodeLabel)) {
    if (sum(parentGeneList == NodeLabel[i]) > 0) {
      NodeSize = c(NodeSize, SizeParent)
      NodeColor = c(NodeColor, parentColor)
      if (parentgeneNames == FALSE) {
        WritenNodeLabel[i] = ""
      }
    }
    else {
      if (sum(targetGeneList == NodeLabel[i]) > 0) {
        NodeSize = c(NodeSize, SizeTarget)
        NodeColor = c(NodeColor, targetColor)
        if (targetgeneNames == FALSE) {
          WritenNodeLabel[i] = ""
        }
      }
    }
  }
  CPstartList = sort(unique(ARTIVAnet$CPstart))
  if (onepage) {
    par(mfrow = c(1, length(unique(ARTIVAnet$CPstart))))
  }
  for (i in 1:length(CPstartList)) {
    CurrentCPstart = CPstartList[i]
    EdgesToDelete = unique(c(which(ARTIVAnet$CPstart > CurrentCPstart), 
                             which(ARTIVAnet$CPend < CurrentCPstart), which(ARTIVAnet$PostProb < 
                                                                              edgesThreshold)))
    SubNetwork = delete.edges(GlobalNetwork, E(GlobalNetwork)[EdgesToDelete])# - 1]) # !!!!
    if (i == length(CPstartList)) {
      MainTitle = paste("Sub-network #", i, "\n( time point", 
                        CurrentCPstart, "to ", max(ARTIVAnet$CPend), 
                        ")")
    }
    else {
      MainTitle = paste("Sub-network #", i, "\n( time point", 
                        CurrentCPstart, "to ", CPstartList[i + 1] - 1, 
                        ")")
    }
    PlotBubbleGraph(SubNetwork, NodeSize, NodeColor, WritenNodeLabel, 
                 LabelPos = 0, NodeCoord, MainTitle)
    if (!onepage | i == 1) {
      legend(par("usr")[1], par("usr")[3], c("Positive interaction", 
                                             "Negative interaction"), lty = c(1, 2))
    }
  }
}