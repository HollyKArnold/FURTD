###############################################################################
#                                                                             #
###############################################################################

#AUTHOR: HOLLY ARNOLD
#DAY: September 25, 2022
#DATE: 20220925
#DESCRIPTION: 
# Provides functions for analyses performed in paper chronic clinical signs of 
# upper respiratory microbiomes in a cohort of domestic felines


# FUNCTION: cladifierGetRefTreeTaxonomy
# INPUTS:
#   tree: The reference tree with nodes and tips labeled
#         All tree tips must be rownames in the reference taxonomy data frame
#   ref: a reference taxonomy file with the following columns:
#        Kingdom, Phylum, Class, Order, Family, Genus, Species
cladifierGetRefTreeTaxonomy = function(tree, ref){

  # Check that all the tree tip labels are contained in the reference taxonomy file
  if(sum(tree$tip.label %in% rownames(ref)) != length(tree$tip.label)){
    stop("cladifierGetRefTreeTaxonomy: All tree tips are not contained in reference taxonomy file.")
  }

  # make a data frame for taxonomy output
  N = length(tree$tip.label) + length(tree$node.label)
  d = data.frame(matrix(nrow = N, ncol = 7, data = NA))
  colnames(d) = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
  rownames(d) = c(tree$tip.label, tree$node.label)

  for(i in 1:N){


    # i is a tip
    curNode = rownames(d)[i]
    if(curNode %in% tree$tip.label){

      d[curNode, "Kingdom"] = ref[curNode, "Kingdom"]
      d[curNode, "Phylum"] = ref[curNode, "Phylum"]
      d[curNode, "Class"] = ref[curNode, "Class"]
      d[curNode, "Order"] = ref[curNode, "Order"]
      d[curNode, "Family"] = ref[curNode, "Family"]
      d[curNode, "Genus"] = ref[curNode, "Genus"]
      d[curNode, "Species"] = ref[curNode, "Species"]

    }else{  # i is a node

      if(i %% 100 == 0){
        print(i/N)
      }

      curNodeIdx = which(tree$node.label == curNode) + length(tree$tip.label)
      curTree = ape::extract.clade(tree, curNodeIdx)
      curRef = ref[curTree$tip.label,]
      uniqueTax = apply(curRef, 2, function(x) length(unique(x))==1)
      uniqueTax = uniqueTax[c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")]

      j = 1
      while(as.vector(uniqueTax[j]) & j < 8){

        d[curNode, names(uniqueTax)[j]] = unique(curRef[,names(uniqueTax)[j]])
        j = j + 1

      }

    }


  }
  return(d)
}



# Takes a ps object and returns otu_table and metadata table combined together for GLMs
getModelingTable = function(phyloseq){

  curASV = phyloseqCompanion::otu.matrix(phyloseq)
  curASV = as.data.frame(curASV)
  curASV = as_tibble(curASV, rownames = "ID")
  curMeta = phyloseqCompanion::sample.data.frame(phyloseq)
  d = curASV %>% left_join(curMeta, by = "ID")

  return(d)

}



# Return a tree dataframe for mapping covariates
treeDataFrame = function(colnames, tree){

  N = length(tree$tip.label) + length(tree$node.label)

  d = data.frame(matrix(ncol = length(colnames) + 1, nrow = N, data = NA))
  colnames(d) = c("node", colnames)
  d$node = seq(from = 1, to = N, by = 1)
  rownames(d) = c(tree$tip.label, tree$node.label)

  return(d)
}


# DESCRIPTION:  This function takes a list of nodes within a phylogenetic
#               tree and returns a list of nodes that are desendents of nodes
#               within the list
# INPUT:        cladeList: a dataframe with rownames as clades and a column
#               $Level containing the level of clades.
#               All node names in cladeList must be a subset of nodenames of the
#               phylogenetic tree.
#               tree: a phylogenetic tree with which to process cladeList
# OUTPUT        A list of nodes where nested nodes have been removed.
RemoveNestedClades = function(cladeList, tree, print = T){

  # 1. Check to make sure all names in cladeList is in tree node labels
  if(sum(rownames(cladeList) %in% tree$node.label)
     != length(rownames(cladeList))){
    print("ERROR: All clades in clade list are not in phylogenetic tree.")
    stop()
  }

  # 2. Convert string node labels ("nodeXXXX") to a numerical label list.
  nodeIDList = vector()
  for(i in 1:length(rownames(cladeList))){
    nodeIDList = c(nodeIDList,
                   GetNodeID(tree = tree, nodeName = rownames(cladeList)[i]))
  }

  # 3. Make a converted nodeID Column in cladeList
  cladeList$nodeID = nodeIDList

  # 4. Order cladeList by level ### DECREADING MUST =T for this to work##
  cladeList = cladeList[order(cladeList$Level, decreasing = T),]

  # 5. Iterate through each of the members of the cladeList
  for(i in 1:length(rownames(cladeList))){

    # 5A. Get descendants of each of the clade members
    d = phangorn::Descendants(x = tree, node = cladeList$nodeID[i], type = "all")

    # 5B. Iterate through each of the descendants
    for(j in 1:length(d)){
      if(d[j] %in% nodeIDList){ # 5B.1 The descendant is in the nodeList
        # Remove descendants
        nodeIDList = nodeIDList[-which(nodeIDList == d[j])]
      }
    }
  }

  # 6. Convert each of the nodeIDs back to node names
  reducedNested = vector()
  for(i in 1:length(nodeIDList)){
    reducedNested = c(reducedNested,
                      GetNodeName(tree = tree, nodeID = nodeIDList[i]))
  }

  # 7. Return list of un-nested nodes.
  return(reducedNested)
}


# DESCRIPTION   GetNodeID takes a nodeName and returns a numeric nodeID
#               assigned nodeID given within the phylogenetic tree.
# INPUTS        tree - a phylogenetic tree which includes the nodeName as
#               a part of the tree$node.labels
#               nodeName - a single node name which exists in tree$node.label
# OUTPUT        a numeric nodeID of the nodeName.
GetNodeID = function(tree, nodeName){
  if(!nodeName %in% tree$node.label){
    stop("ERROR, nodeName is not contained in phylogenetic tree")
  }
  ntips = length(tree$tip.label)
  node = which(tree$node.label == nodeName)
  return(ntips + node)
}


# DESCRIPTION   GetNodeName takes a nodeID and a tree and converts to the node
#               name from tree$node.labels
# INPUT         tree - a phylogenetic tree
#               nodeID - a single nodeID to convert to a node label
# OUTPUT        The node name that is given to a nodeID in tree$node.label
GetNodeName = function(tree, nodeID){
  ntips = length(tree$tip.label)
  if(as.numeric(nodeID) <= ntips){
    stop("ERROR: You've entered an invalid nodeID.")
  }
  node = tree$node.label[nodeID - ntips]
  return (node)
}

# DESCRIPTION:    This function converts a list of nodeIDs to nodeNames
#                 that are present within a phylogenetic tree
# INPUTS:         tree - a phylogenetic tree of class phylo
#                 nodeIDs - a list of nodeIDs to be converted to node names.
# OUTPUTS:        A list of node names from tree$node.label
GetNodeNames = function(tree, nodeIDs){
  names = vector()
  for(i in 1:length(nodeIDs)){
    names = c(names, GetNodeName(tree = tree, nodeID = nodeIDs[i]))
  }
  return(names)
}



# FUNCTION: getRefTreeCladeLevel
# INPUTS: t: A tree where all nodes are labeled and all tips are labeled
# OUTPUT: a matrix with the following attributes
#         node: the tip or node name of the tree.
#         level: the level of each clade where root is the 1 and each node increases by 1.

getRefTreeCladeLevel = function(tree){

  # Get number of nodes in tree
  N = length(tree$tip.label) + length(tree$node.label)

  # Make data matrix
  a = phangorn::Ancestors(x = tree, node = 1:N, type = "all")

  d = as.data.frame(matrix(ncol = 2, nrow = N))
  colnames(d) = c("node", "Level")
  names = c(tree$tip.label, tree$node.label)
  
  for(i in 1:N){
    d[i,"node"] = names[i]
    d[i,"Level"] = length(a[[i]])
  }

  return(d)
}

# Returns distances from a combination vector
getDistanceCombos = function(combos, distMatrix){

  distance = vector()
  
  for(i in 1:ncol(combos)){
   
    curDist = distMatrix[combos[1,i], combos[2,i]]
    distance = c(distance, curDist)
  }
  return(distance)

}

# Get random distance
distRandom = function(dist, n){

  distance = vector()
  for(i in 1:n){
    cur = sample(x = 1:dim(dist)[1], 2, replace = FALSE)
    distance = c(distance, dist[cur[1], cur[2]])
  }

  return(distance)

}

# Makes Negative Binomial model
# Tests GLM NB with base variables + microbial feature
# Returns a dataframe with covariate tested, estimate, microbial feature, 
# p value, and q value
# q values are corrected with fdr method *per covariate*
getNBGLMS = function(covariates, phyloseq){
  
  curData = as.data.frame(getModelingTable(phyloseq))
  curMicrobes = taxa_names(phyloseq)
  
  d = data.frame(matrix(nrow = length(covariates) * length(curMicrobes), ncol = 6))
  colnames(d) = c("covariate", "clade", "Estimate", "p", "q", "model")
  
  count = 0
  for(i in 1:length(covariates)){
    
    curCovs = covariates[i]
    print("Testing Covariate: ")
    print(curCovs)
    
    for(j in 1:length(curMicrobes)){
      
      count = count + 1
      
      curFormula = as.formula(paste(c(curMicrobes[j], "~", curCovs), sep = "", collapse = ""))
      d[count, "covariate"] = curCovs
      d[count, "clade"] = curMicrobes[j]
      
      
      curModel = tryGLMNBModel(curFormula = curFormula, curData = curData)
      if(is.na(curModel)){
        
        d[count, "Estimate"] = NA
        d[count, "p"] = NA
        d[count, "model"] = "No Convergence"
      }else{
        
        if(is.null(curModel$th.warn)){ # if there is not a warning, record those results. Remember zeroinfl takes care of some errors that glm.nb has which is why we tried that first
          d[count, "Estimate"] = summary(curModel)$coefficients[curCovs, "Estimate"]
          d[count, "p"] = summary(curModel)$coefficients[curCovs, "Pr(>|z|)"]
          d[count, "model"] = "glm.nb"
        }else{
          
          curData[[curMicrobes[j]]] = curData[[curMicrobes[j]]] + 0.00001
          curModel = tryGLMNBModel(curFormula = curFormula, curData = curData)
          if(is.na(curModel)){
            curModel = glm(formula = curFormula, family = "poisson", data = curData, na.action = "na.omit") # now if the glm.nb fails, we can try poisson
            d[count, "Estimate"] = summary(curModel)$coefficients[curCovs, "Estimate"]
            d[count, "p"] = summary(curModel)$coefficients[curCovs, "Pr(>|z|)"]
            d[count, "model"] = "poisson"
          }else{
            d[count, "Estimate"] = summary(curModel)$coefficients[curCovs, "Estimate"]
            d[count, "p"] = summary(curModel)$coefficients[curCovs, "Pr(>|z|)"]
            d[count, "model"] = "glm.nb.plus.01"
          }
          
        }
        
      }
      
    }
    
    d[which(d$covariate == curCovs & !is.na(d$p)),"q"] = p.adjust(p = as.vector(d[which(d$covariate == curCovs & !is.na(d$p)),"p"]), method = "fdr")

  }
  
  return(d)
  
}

tryGLMNBModel = function(curFormula, curData){

  out <- tryCatch(
    {
      curModel = MASS::glm.nb(formula = curFormula, 
                              data = curData, 
                              na.action = "na.omit", 
                              link = "log")# try glm.nb first


    },
    error=function(cond) {
      # choose what to return if error
      return(NA)
    },
    #warning=function(cond) {
    #choose what to return if warning
    #return(NA)
    #},
    finally={
      # do this regardless
    }
  )
  return(out)
}


# Make a list of layers to apply to a geom_cladelabel figure
cladeLabels = function(dataFrame, idxs, nodeName,  fontsize, bsize, offset.bar){

  cladelabels = ""
  for(i in 1:nrow(dataFrame)){
    if(is.na(dataFrame[i, "Phylum"])){curLabel = dataFrame[i, "Kingdom"]}
    else if( is.na(dataFrame[i, "Class"])){curLabel = dataFrame[i, "Phylum"]}
    else if( is.na(dataFrame[i, "Order"])){curLabel = dataFrame[i, "Class"]}
    else if( is.na(dataFrame[i, "Family"])){curLabel = dataFrame[i, "Order"]}
    else if( is.na(dataFrame[i, "Genus"])){curLabel = dataFrame[i, "Family"]}
    else if( is.na(dataFrame[i, "Species"])){curLabel = dataFrame[i, "Genus"]}

    
      
    curString = paste(c("geom_cladelabel(node = which(",
                        idxs,
                        "=='",
                        dataFrame[i, nodeName],
                        "'), label = '",
                        curLabel,
                        "', align = T, geom = 'text', angle = 0, color = ",
                        dataFrame[i, "phylumBarKeyColor"],
                        ", fontsize = ",
                        fontsize,
                        ", barsize = ",
                        bsize,
                        ", offset = " ,
                        offset.bar,
                        ")"

    ), sep = "", collapse = "")
    cladelabels = paste(c(cladelabels, curString), sep = "  ", collapse = " + ")
  }

  return(cladelabels)

}

