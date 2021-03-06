PlotBubbleGraph <- function (Network, NodeSize, NodeColor, NodeLabel, LabelPos = 0, 
                             NodeCoord, MainTitle) 
{
  plot(Network, vertex.size = NodeSize, vertex.color = NodeColor, 
       vertex.label = NodeLabel, vertex.label.dist = LabelPos, 
       layout = NodeCoord, edge.arrow.size = 0.5, main = MainTitle)
}