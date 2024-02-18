## Introduction to Animal Movement Ecology
## K. Manlove, August 2023
##
## Aggregating range and occurrence distributions in R
## and building home range overlap & contact networks
## 

# I. Home range overlap networks ----


# II. Contact ("association") networks from GPS data ----
datFull <- read.csv("Data/GPS/ZionLocsThru2020.csv", header = T)

# cut down to just one focal year (here, 2019)
gps <- subset(datFull, Year == 19)

# rename Longitude lon and Latitude lat
gps$lon <- as.numeric(as.character(gps$Longitude))
gps$lat <- as.numeric(as.character(gps$Latitude))

# define dayInStudy (if working with a single focal year, this can just be JulianDay)
gps$dayInStudy <- gps$Julianday

## A. build and store some requisite functions ----
# ReplaceLowerOrUpperTriangle takes the lower (or upper) triangle of a symmetric matrix
# and reflects it across the matrix's main diagonal to build a symmetric matrix.
ReplaceLowerOrUpperTriangle <- function(m, 
                                        triangle.to.replace){
  # If triangle.to.replace="lower", replaces the lower triangle of a square matrix with its upper triangle.
  # If triangle.to.replace="upper", replaces the upper triangle of a square matrix with its lower triangle.
  
  if (nrow(m) != ncol(m)) stop("Supplied matrix must be square.")
  if      (tolower(triangle.to.replace) == "lower") tri <- lower.tri(m)
  else if (tolower(triangle.to.replace) == "upper") tri <- upper.tri(m)
  else stop("triangle.to.replace must be set to 'lower' or 'upper'.")
  m[tri] <- t(m)[tri]
  return(m)
}

# GeoDistanceInMeters calculates the distance between two locations in meters
GeoDistanceInMetres <- function(g1, g2){
  # Returns a vector of distances. (But if g1$index > g2$index, returns zero.)
  # The 1st value in the returned vector is the distance between g1[[1]] and g2[[1]].
  # The 2nd value in the returned vector is the distance between g1[[2]] and g2[[2]]. Etc.
  # Each g1[[x]] or g2[[x]] must be a list with named elements "index", "lat" and "lon".
  # E.g. g1 <- list(list("index"=1, "lat"=12.1, "lon"=10.1), list("index"=3, "lat"=12.1, "lon"=13.2))
  DistM <- function(g1, g2){
    return(ifelse(g1$index > g2$index, 0, gdist(lat.1 = g1$lat, 
                                                lon.1 = g1$lon, 
                                                lat.2 = g2$lat, 
                                                lon.2 = g2$lon, 
                                                units = "m")))
  }
  
  return(mapply(DistM, g1, g2))
}

# GeoDistanceInMetresMatrix takes a set of points and builds
# a symmetric association matrix on the basis of the points. 
GeoDistanceInMetresMatrix <- function(df.geopoints){
  # Returns a matrix (M) of distances between geographic points.
  # M[i,j] = M[j,i] = Distance between (df.geopoints$lat[i], df.geopoints$lon[i]) and
  # (df.geopoints$lat[j], df.geopoints$lon[j]).
  # The row and column names are given by df.geopoints$name.
  n.geopoints <- nrow(df.geopoints)
  
  # The index column is used to ensure we only do calculations for the upper triangle of points
  df.geopoints$index <- 1:n.geopoints
  
  # Create a list of lists
  list.geopoints <- by(df.geopoints[,c("index", "lat", "lon")], 
                       1:n.geopoints, function(x){return(list(x))})
  
  # Get a matrix of distances (in metres)
  mat.distances <- ReplaceLowerOrUpperTriangle(outer(list.geopoints, list.geopoints, GeoDistanceInMetres), "lower")
  
  # Set the row and column names
  rownames(mat.distances) <- df.geopoints$name
  colnames(mat.distances) <- df.geopoints$name
  
  return(mat.distances)
}

## B. Construct association matrix for GPS data ----
# loop through days, get pairwise distances between all points in a 
# given day.  If distance < 200m, call points "together", otherwise, "apart".
# 1) build some storage objects
ObsTogetherMat <- ObsApartMat <- matrix(0, nrow = length(levels(factor(gps$animal))),
                                        ncol = length(levels(factor(gps$animal))))

colnames(ObsTogetherMat) <- rownames(ObsTogetherMat) <- levels(factor(gps$animal))
colnames(ObsApartMat) <- rownames(ObsApartMat) <- levels(factor(gps$animal))

# 2) loop over study-days (pretty slow -- maybe 30 minute run on KRM computer)
for(i in min(na.omit(gps$dayInStudy)):max(na.omit(gps$dayInStudy))){
  # need to store: 
  # days animal 1 was observed
  # days animal 2 was observed
  # days both animal 1 and animal 2 were observed
  
  k <- subset(gps, dayInStudy == i)
  
  # pull off first observation for each animal at this timestep
  if(dim(k)[1] >= 1){
    # calculate all pairwise distances
    distToday <- GeoDistanceInMetresMatrix(k)
    row.names(distToday) <- colnames(distToday) <- k$animal
    firstAnimal <- k$Animal[which(distToday <= 200, arr.ind = T)[, 1]]
    secondAnimal <- k$Animal[which(distToday <= 200, arr.ind = T)[, 2]]
    lookupCoords <- which(distToday <= 200, arr.ind = T)
    lookupCoordsApart <- which(distToday > 200, arr.ind = T)
    
    for(j in 1:dim(lookupCoords)[1]){
      ObsTogetherMat[rownames(distToday)[lookupCoords[j, 1]], 
                     rownames(distToday)[lookupCoords[j, 2]]] <- ObsTogetherMat[rownames(distToday)[lookupCoords[j, 1]], 
                                                                                rownames(distToday)[lookupCoords[j, 2]]] + 1
      
      # increase ObsTogetherMat cell's counter by 1 each time a dyad occurs together
    }
    
    if(dim(lookupCoordsApart)[1] >= 1){
      for(m in 1:dim(lookupCoordsApart)[1]){
        ObsApartMat[rownames(distToday)[lookupCoordsApart[m, 1]], 
                    rownames(distToday)[lookupCoordsApart[m, 2]]] <- ObsApartMat[rownames(distToday)[lookupCoordsApart[m, 1]], 
                                                                                 rownames(distToday)[lookupCoordsApart[m, 2]]] + 1
        # increase ObsApartMat cell's counter by 1 each time a dyad is observed, but NOT together
      }
    }
  }
  print(i)
}

# check to be sure ObsTogetherMat and ObsApartMat were actually calculated
sum(ObsTogetherMat)
sum(ObsApartMat)

# build association matrix (this is a simple ratio index).
# I add ther very small value here to keep from dividing 0 by 0. 
AssocMat <- ObsTogetherMat / (ObsTogetherMat + ObsApartMat + .0000001)

# write out simple ratio association matrix so that loop doesn't have to be re-run in the future.
write.csv(AssocMat, "Code/PracticeZionAssociationMatrix_2019_20210618.csv")

# III. Build networks, visualize, and analyze using igraph ----
# Networks short-course
# K. Manlove
# Lecture 02: Graph metrics 
library(igraph)
library(asnipe)

# I. Build graph and assign vertex and edge attributes ----
# Read in association matrix (built from Script Lecture01a_NetworkConstructionFromGPSData)
AssocMat <- read.csv("Data/PracticeZionAssociationMatrix_2019_20210618.csv")
rownames(AssocMat) <- AssocMat[, 1]
AssocMat <- AssocMat[, -1]
AssocMat <- as.matrix(AssocMat)

# split leading X's off column names
colnames(AssocMat) <- substr(colnames(AssocMat), start = 2, stop = 8)

# Read in GPS data (which contains some node covariates we'll want):
#datFull <- read.csv("Data/ZionLocsThru2020.csv", header = T)
datFull <- read.csv("Data/Zion_Covariates_thru20210108.csv", header = T)
# cut down to just one focal year (here, 2019)
gps <- subset(datFull, Year == 19)

# A. Overwrite the main diagonal of the association matrix to be 0 (by default, an individual 
#    will always be "with" itself, but we don't care about those self-loops here, 
#    so we set them to 0 by definition)
diag(AssocMat) <- rep(0, dim(AssocMat)[1]) 
# B. Convert association matrix to graph using igraph::graph_from_adjacency_matrix
zionGraph <- graph_from_adjacency_matrix(AssocMat, 
                                         weighted = T, # specify that edges should be weighted, not binary
                                         mode = "undirected") # specify that graph is undirected (e.g., associations are symmetric: if A is with B, then B is equivalently with A. )



# II. Extracting and assessing graph metrics ----
# A. Node-specific metrics ----
# i.Degree: how many connections does each node have? ----
# degree() counts up the number of non-zero edges connected to each node, so
# it tells you have many nodes each node is connected to. 
# degree() function extracts degree associated with each node in the graph.
zion_graph_degree <- degree(zionGraph)
# often, the degree distribution is of interest (is is uniform-ish, exponential-ish,
# or very long-tailed, which is to say, power-law-ish.)
par(mfrow = c(1, 1))
hist(zion_graph_degree, main = "", las = 1, breaks = 20)
# The degree distribution of the Zion data looks more or less uniform. 
# If it were long-tailed, we could estimate the power-law coefficient using
zion_power_law_fit <- fit_power_law(zion_graph_degree)
zion_power_law_fit

# betweenness
zion_between <- betweenness(zionGraph,
                            weight = E(zionGraph)$weight)




# simulated power-law data


# B. Graph-level metrics ----
# i. Communities and modularity ----
#    Modularity measures the extent to which a graph breaks down into semi-distinct
#    subunits. Calculating modularity is a two-step process. First, we have to identify
#    communities within the graph.  There are a TON of methods for determining 
#    communities, including eigen-decomposition-based methods, methods based on
#    conventional clustering algorithms, and methods based on "walkers" who 
#    take steps around the graph.  My preference for real-world ecological networks
#    is usually the latter approach, and I favor the walktrap algorithm as
#    an implementation, but this is an area of active research.  I usually
#    build communities a few different ways (an eigendecomposition-based method, a walktrap-
#    based method, etc.) and compare.  If methods agree, I feel good; if they disagree
#    I exercise appropriate skepticism and go with the walktrap. 
walktrap_zion <- cluster_walktrap(zionGraph,
                                  weights = E(zionGraph)$weight,
                                  steps = 4,
                                  membership = T)
# walktrap_zion is a "communities" object.  I can extract which nodes got assigned to which
# community using membership(walktrap_zion):
zion_community_members <- membership(walktrap_zion)
# I can get the community sizes by tabling membership:
zion_community_tab <- table(zion_community_members)

## NEW 2021-07-16
## Relative modularity
zion_newman_mod <- modularity(zionGraph,
                              zion_community_members)

# make a list to store edges connecting nodes 
# from community i to nodes from any other community
community_edges_list <- vector("list", 
                               length(levels(factor(zion_community_members))))
# loop over all communities and store edges targeted for deletion
for(i in 1:length(community_edges_list)){
  community_edges_list[[i]] <- E(zionGraph)[ which(zion_community_members == i) %--% 
                                               which(zion_community_members != i) ]
}
# bind all edges to be omitted together into a single "vector" (actually igraph.es)
community_edges_full <- do.call("c", community_edges_list)
# delete all of those edges from zionGraph (and store in NEW OBJECT, zionGraph_sm)
zionGraph_sm <- delete_edges(zionGraph,
                             edges = community_edges_full)
# plot the reduced-edges graph just to check. 
plot(zionGraph_sm,
     vertex.label.color = rgb(0, 0, 0, alpha = 0))
# calculate maximum modularity on the new "zionGraph_sm" object. 
# NOTE: we DO NOT have to redo the walktrap, we just use the membership assignments
# from walktrap_zion above, which we've stored already in zion_community_members
modularity_max <- modularity(zionGraph_sm,
                             membership = zion_community_members)
# calculate relative modularity by taking the Newman modularity over max modularity. 
# this modularity_rel should be comparable across graphs of different sizes and structures
# (per Sah, Bansal, Leu, Cross, Hudson 2017 PNAS)
modularity_rel <- zion_newman_mod / modularity_max

# but SEE THIS STACK DISCUSSION...
# https://igraph.discourse.group/t/calculating-relative-modularity/736

# in this case, community 11 is particularly large (contains 12 animals); community 5
# has 9 members, community 7 has eight members, etc. 
# At this point, I would typically map those members back to the graph (e.g., color-code
# nodes by community membership) to get a better visual sense of which community is "where" in the network. 
# I need 15 different colors, so using rainbow(n = 15) to attain them -- brewer.pals don't extend
# beyond 12 groups. The challenge here is that some of my colors might be hard to tell apart, so
# I'll need to be careful in interpretation.
fillColIn <- rainbow(n = 15)
plot(zionGraph, layout = layoutZion,
     edge.width = E(zionGraph)$weight * 50,
     vertex.color = fillColIn[zion_community_members],
     #vertex.frame.color = frameColIn[factor(nodeSymps)],
     edge.color = rgb(.6, .6, .6, .5),
     vertex.label.color = "black",
     vertex.size = vertSizeIn[factor(nodeLambStatus)])
leg.text <- paste("community ", seq(1:length(zion_community_tab)), sep = "")
legend("bottom", leg.text, fill = fillColIn, bty = "n",
       ncol = 3) # ncol here specifies a number of columns for legend text. 
# I'm using it since I've got 15 elements going into the legend. 

# # VI. Alternative plotting with hierarchical edge bundling (requires community object) ----
# #   A. convert community object to dendrogram
# #      To do this, we'll need to dig into the structure of our walktrap_zion object. 
# #      First, let's build a dendrogram of the walktrap output so we can see what we're looking at:
# plot_dendrogram(walktrap_zion)
# #      We specifically want to look at the "merges" slot, which contains the tree
# #      structure. The merges slots is a 2-column matrix, indicating pairs of nodes that
# #      get joined.  the first 64 nodes correspond directly to our individuals. Nodes
# #      beyond 64 correspond to joinings further up the tree. This is easier to understand
# #      in context, so let's look at the merges:
# walktrap_zion$merges
# # The first row in (my) walktrap_zion$merges matrix has elements [10, 55)]
# # This [10,55] row is says that the first merge in the walktrap is
# # between individuals 10 and 55.  We can look up those individuals' names with:
# walktrap_zion$names[c(10, 55)]
# # For me, these are individuals 39866_2 and 42739_2.  On the dendrogram image,
# # we can see that these two nodes are connected very low down (i.e., far to the right)
# # in the tree. That first row corresponds to a new node, which (though R doesn't say 
# # this explicitly) is treated as node 65. 
# # The next row in my merge object is [64,65].  This node then is joining individual
# # 64 with my newly created node, which was coded 65. We can look up individual 64 with:
# walktrap_zion$names[64] # in my case, individual 45149_1.
# # If we look at the dendrogram, we can see that the connection between 45149_1 and our
# # 10/55 dyad of 39866_2 and 42739_2 is the next-lowest connection in the tree. 
# # This new node will now be latently known as 66. In my graph, the next connection ALSO
# # involves this group, with node 28 joining my previously connected group, now coded 66. 
# 
# # to use the hierarchical edge bundling graphics utilities, we need
# # to extract components of the community merges and names in a particular structure. 
# walktrap_zion$merges
# origin_node <- max(walktrap_zion$merges)
# 
# hierarchy <- data.frame(as.matrix(walktrap_zion$merges))
# names(hierarchy) <- c("from", "to")
# 
# d1 <- data.frame(from="origin", to=paste("group", seq(1,10), sep=""))
# d2 <- data.frame(from=rep(d1$to, each=10), to=paste("subgroup", seq(1,100), sep="_"))
# hierarchy_ex <- rbind(d1, d2)
# 
# # create a vertices data.frame. One line per object of our hierarchy, giving features of nodes.
# vertices <- data.frame(name = unique(c(as.character(hierarchy$from), as.character(hierarchy$to)))) 
# 
# # Create a graph object with the igraph library
# mygraph <- graph_from_data_frame(hierarchy, vertices = vertices)
# # This is a network object, you visualize it as a network like shown in the network section!
# 
# # With igraph: 
# plot(mygraph, vertex.label = "", edge.arrow.size = 0, vertex.size = 2)
# # With ggraph:
# ggraph(mygraph, layout = 'dendrogram', circular = FALSE) + 
#   geom_edge_link() +
#   theme_void()
# 
# ggraph(mygraph, layout = 'dendrogram', circular = TRUE) + 
#   geom_edge_diagonal() +
#   theme_void()
# 
# 
# 
# zion_consensus_tree <- consensus_tree(zionGraph)
# 
# #  Build hierarchical random graph of Zion data. 
# zion_hrg <- fit_hrg(zionGraph)
# 
# plot_dendrogram(zion_hrg)
# # zion_hrg_graph <- sample_hrg(zion_hrg)
# # zion_hrg_graph <- as.directed(zion_hrg_graph)
# 
# zion_hrg_df <- as.data.frame(cbind(zion_hrg$left,
#                                    zion_hrg$right,
#                                    zion_hrg$prob,
#                                    zion_hrg$edges,
#                                    zion_hrg$vertices))
# names(zion_hrg_df) <- c("left", "right", "prob", "edges", "vertices")
# 
# edges <- data.frame(from = zion_hrg_df$left,
#                     to = zion_hrg_df$right)
# all_leaves <- zion_hrg_df$vertices
# connect <- 
#   # 
#   # # what is in the hrg object?
#   # zion_hrg_df
#   # # first row is the root:
#   # # root node -42
#   # zion_hrg_df[which(zion_hrg_df$left == -8), ]
#   
#   # zion_hierarchy_graph <- zion_consensus_tree$hrg
#   # d1 <- data.frame(from = zion_hierarchy_graph$left,
#   #                  to = zion_hierarchy_graph$right)
#   # hierarchy <- d1
# # 
# # zion_dendro <- igraph_hrg_dendrogram(zion_hierarchy_graph)
# ggraph(walktrap_zion,
#        layout = 'dendrogram', circular = TRUE) + 
#   geom_edge_link(size=0.4, alpha=0.1) +
#   geom_node_text(aes(x = x*1.01, y=y*1.01, 
#                      # filter = leaf, 
#                      #label=name, 
#                      # angle = angle, 
#                      hjust=hjust), 
#                  size=1.5, alpha=1) +
#   coord_fixed() +
#   theme_void() +
#   theme(
#     legend.position="none",
#     plot.margin=unit(c(0,0,0,0),"cm"),
#   ) +
#   expand_limits(x = c(-1.2, 1.2), y = c(-1.2, 1.2))

