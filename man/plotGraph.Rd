\name{plotGraph}
\alias{plotGraph}
\title{Plot graph for the pathway network}
\description{
Plot graph for the pathway network
}
\usage{
plotGraph(x, node.name = NULL, node.type = NULL,
    graph.node.max.size = 20, graph.node.min.size = 3, graph.layout.method = NULL, ...)                 
}
\arguments{
  \item{x}{a \link{cepa} object}
  \item{node.name}{node.name for each node}
  \item{node.type}{node.type for each node}
  \item{graph.node.max.size}{max size of the node in the graph}
  \item{graph.node.min.size}{min size of the node in the graph}
  \item{graph.layout.method}{function of the layout method. For the list
  of available methods, see \link[igraph]{layout}}
  \item{...}{other arguments}
}
\details{
Graph view of the pathway where the size of node is proportional to centrality 
value of the node. 

By default, the layout for the pathway tree-like. If the number of pathway nodes
is large, the layout would be a random layout.
 
The function returns a \link[igraph]{igraph} object, so if users don't satisfy 
with the network visualization (To tell the truth, if there are too many pathway nodes,
the layout would be real ugly. :(), they can draw graphs of the network with their
own settings.

The function is always called through \code{\link{plot.cepa.all}} and \code{\link{plot.cepa}}.
}
\value{
A \link[igraph]{igraph} object of the pathway
}
\author{Zuguang Gu}
\seealso{
\code{\link{cepa}}, \code{\link{plot.cepa}}
}
\examples{
\dontrun{
data(PID.db)

# ORA extension
data(gene.list)
# will spend about 20 min
res.ora = cepa.all(dif = gene.list$dif, bk = gene.list$bk, pc = PID.db$NCI)
ora = get.cepa(res.ora, id = 5, cen = 3)
plotGraph(ora)

# GSA extension
# P53_symbol.gct and P53_cls can be downloaded from
# http://mcube.nju.edu.cn/jwang/lab/soft/cepa/
eset = read.gct("P53_symbol.gct")
label = read.cls("P53.cls", treatment="MUT", control="WT")
# will spend about 45 min
res.gsa = cepa.all(mat = eset, label = label, pc = PID.db$NCI)
gsa = get.cepa(res.gsa, id = 5, cen = 3)
plotGraph(gsa)
}
}
