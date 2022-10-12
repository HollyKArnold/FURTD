###############################################################################
#                                                                             #
###############################################################################

#AUTHOR: HOLLY ARNOLD
#DAY: September 25, 2022
#DATE: 20220925
#DESCRIPTION: 
# Provides functions for analyses performed in paper chronic clinical signs of 
# upper respiratory microbiomes in a cohort of domestic felines
#' 
#' #AUTHOR: ARNOLD
#' #DAY: July 8th, 2019
#' #DATE: 20190708
#' #PURPOSE:
#' # To identify a proper rarification threshold for the sheep analysis. 
#' 
#' 
#' # FUNCTION: splitDada2Tax
#' # PURPOSE: To take a data 2 taxonomy table and split into a data frame
#' #          which only has one taxonomy classification per column
#' # INPUTS:
#' #    t: a dada2 table table which has the following properties
#' #       col1: sequenceID
#' #       col2: taxonomy string with the following properties: 
#' #             k__<Kingdom>; p__<Phylum>; c__<Class>; o__<Order>; f__<family>; g__<genus>
#' # OUTPUTS
#' # m: a matrix
#' #  col: SeqID Kingdom Phylum Class Order Family Genus
#' #  row: Number of sequences
#' splitDada2Tax = function(t, splitChar){
#'   
#'   # 0. Make a new matrix
#'   m = data.frame(matrix(nrow = nrow(t), ncol = 8, data = NA))
#'   colnames(m) = c("SeqID", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
#'   
#'   for(i in 1:nrow(t)){
#'     seqCur = as.vector(t[i, 1])
#'     m[i, 1] = seqCur
#'     taxCur = as.vector(t[i, 2])
#'     taxCur = strsplit(taxCur, split = splitChar)
#'     
#'     for(j in 1:length(taxCur[[1]])){
#'       m[i, (j + 1)] = taxCur[[1]][j]
#'     }
#'     
#'   }
#'   rownames(m) = m$SeqID
#'   return(m)
#' }
#' 
#' # DESCRIPTION:  taxString2taxTable 
#' #               Takes a table from getTaxString and then returns a dada2 formated tax table
#' # INPUTS:       
#' #    tax table like from getTaxString. One column must contain tax string, one column must contain row ID
#' #    taxString: Column name that contains taxonomic string to split (i.e. k__Bacteria;)
#' #    taxID: Column name to get the rownames from for tax strings. tax string identifiers (tax1, tax2)
#' #    split: the character that splits fields between taxonomic levels
#' # OUTPUTS:     
#' #    a dada2 formated tax table where each column is a Linnean taxonomic levle
#' taxString2taxTable = function(tax, taxString = "taxString", taxID = "taxID", split = ";") {
#'   
#'   taxTable = data.frame(matrix(nrow = nrow(tax), ncol = 7, data = NA))
#'   colnames(taxTable) = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
#'   
#'   for(i in 1:nrow(tax)){
#'     
#'     curTaxString = tax$`taxString`[i]
#'     curTaxString = str_split(string = curTaxString, pattern = split)[[1]]
#'     
#'     if(length(curTaxString) >= 1){
#'       if(grepl(curTaxString[1], pattern = "k__")){
#'         taxTable[i, "Kingdom"] = str_sub(string = curTaxString[1], start = 4, end = nchar(curTaxString[1]))
#'       }else{
#'         print(curTaxString)
#'         stop("taxString2taxTable: Unexpected input character string is found:")
#'       }
#'     }
#'     
#'     if(length(curTaxString) >= 2){
#'       if(grepl(curTaxString[2], pattern = "p__")){
#'         taxTable[i, "Phylum"] = str_sub(string = curTaxString[2], start = 4, end = nchar(curTaxString[2]))
#'       }else{
#'         print(curTaxString)
#'         stop("taxString2taxTable: Unexpected input character string is found:")
#'       }
#'     }
#'     
#'     if(length(curTaxString) >= 3){
#'       if(grepl(curTaxString[3], pattern = "c__")){
#'         taxTable[i, "Class"] = str_sub(string = curTaxString[3], start = 4, end = nchar(curTaxString[3]))
#'       }else{
#'         print(curTaxString)
#'         stop("taxString2taxTable: Unexpected input character string is found:")
#'       }
#'     }
#'     if(length(curTaxString) >= 4){
#'       if(grepl(curTaxString[4], pattern = "o__")){
#'         taxTable[i, "Order"] = str_sub(string = curTaxString[4], start = 4, end = nchar(curTaxString[4]))
#'       }else{
#'         print(curTaxString)
#'         stop("taxString2taxTable: Unexpected input character string is found:")
#'       }
#'     }
#'     if(length(curTaxString) >= 5){
#'       if(grepl(curTaxString[5], pattern = "f__")){
#'         taxTable[i, "Family"] = str_sub(string = curTaxString[5], start = 4, end = nchar(curTaxString[5]))
#'       }else{
#'         print(curTaxString)
#'         stop("taxString2taxTable: Unexpected input character string is found:")
#'       }
#'     }
#'     if(length(curTaxString) >= 6){
#'       if(grepl(curTaxString[6], pattern = "g__")){
#'         taxTable[i, "Genus"] = str_sub(string = curTaxString[6], start = 4, end = nchar(curTaxString[6]))
#'       }else{
#'         print(curTaxString)
#'         stop("taxString2taxTable: Unexpected input character string is found:")
#'       }
#'     }
#'     if(length(curTaxString) >= 7){
#'       if(grepl(curTaxString[5], pattern = "s__")){
#'         taxTable[i, "Species"] = str_sub(string = curTaxString[7], start = 4, end = nchar(curTaxString[7]))
#'       }else{
#'         print(curTaxString)
#'         stop("taxString2taxTable: Unexpected input character string is found:")
#'       }
#'     }
#'     
#'     
#'   }
#'   rownames(taxTable) =  tax$`taxID`
#'   return(taxTable)
#'   
#'   
#' }
#' 
#' # FUNCTION: getCollectorsCurves
#' # PURPOSE: To get a set of collectors curves across different read depths
#' # INPUTS:
#' #  1. a: an ASV table where the columns are the samples and ASVs are the rows
#' #  2. r: a vector which states the rarification depths to rarify at.
#' #  3. seed: for set_seed for rarefy_even_depth() phyloseq for reproducible ASV table
#' # OUTPUT: 
#' #  1. m: a matrix which has the rarification depth as rownames. Each row 
#' #     corresponds to one collectors curve rarified at dept r[i]
#' getCollectorsCurves = function(a, r, seed){
#'   
#'   #0A. Check that the asv table is in the right orientation
#'   if(nrow(a) < ncol(a)){
#'     stop("There should be more sequences than samples.\n")
#'   }
#'   
#'   #1. Make phyloseq object
#'   a = otu_table(a, taxa_are_rows = T)
#'   
#'   #2. Make return matrix
#'   m = as.data.frame(matrix(nrow = length(r), ncol = ncol(a), data = NA))
#'   for(i in 1:nrow(m)){
#'     rownames(m)[i] = paste0(c(r[i]), sep = "", collapse = "")
#'   }
#'   for(i in 1:ncol(m)){
#'     colnames(m)[i] = paste0(c(i), sep = "", collapse = "")
#'   }
#'   
#'   #3. Rarify in a loop
#'   for(i in 1:length(r)){
#'     
#'     curASV = rarefy_even_depth(a, sample.size = r[i], rngseed = seed, replace = F)
#'     aveCurve = getCollectorCurvesConfidenceInterval(n = 100, a = curASV)
#'     aveCurve = apply(aveCurve, 2, ave)[1,]
#'     m[i, 1:length(getCollectorCurve(curASV))] = aveCurve
#'     
#'   }
#'   return(m)
#' }
#' 
#' 
#' # FUNCTION: getCollectorsCurves
#' # PURPOSE: To get a set of collectors curves across different read depths
#' # INPUTS:
#' #  1. a: an ASV table where the columns are the samples and ASVs are the rows
#' #  2. r: a vector which states the rarification depths to rarify at.
#' #  3. seed: for set_seed for rarefy_even_depth() phyloseq for reproducible ASV table
#' # OUTPUT: 
#' #  1. m: a matrix which has the rarification depth as rownames. Each row 
#' #     corresponds to one collectors curve rarified at dept r[i]
#' getCollectorsCurvesTaxonomy = function(a, r, seed, t){
#'   
#'   #0A. Check that the asv table is in the right orientation
#'   if(nrow(a) < ncol(a)){
#'     stop("There should be more sequences than samples.\n")
#'   }
#'   
#'   #1. Make phyloseq object
#'   a = otu_table(a, taxa_are_rows = T)
#'   
#'   #2. Make return matrix
#'   m = as.data.frame(matrix(nrow = length(r), ncol = ncol(a), data = NA))
#'   
#'   for(i in 1:nrow(m)){
#'     rownames(m)[i] = paste0(c(r[i]), sep = "", collapse = "")
#'   }
#'   for(i in 1:ncol(m)){
#'     colnames(m)[i] = paste0(c(i), sep = "", collapse = "")
#'   }
#'   
#'   #3. Rarify in a loop
#'   for(i in 1:length(r)){
#'     
#'     curASV = rarefy_even_depth(a, sample.size = r[i], rngseed = seed, replace = F)
#'     aveCurve = getCollectorCurvesTaxonomyConfidenceInterval(n = 100, a = curASV, t = t)
#'     aveCurve = apply(aveCurve, 2, ave)[1,]
#'     m[i, 1:length(getCollectorCurve(curASV))] = aveCurve
#'     
#'   }
#'   return(m)
#' }
#' 
#' 
#' # getCollectorCurvesConfidenceInterval
#' # FUNCTION: getCollectorCurves
#' # PURPOSE: Returns n collectors curves
#' # INPUTS: 
#' #   1. n - an integer stating how many collectors curves to get
#' #   2. a - an asv table where sequences are the rows and samples are the columns
#' # OUTPUTS:
#' #   1. a matrix where row i is a single collectors curve in full, and columns are points in the curve
#' getCollectorCurvesConfidenceInterval = function(n, a){
#'   
#'   #0. Do a couple of checks.
#'   if(nrow(a) < ncol(a)){
#'     stop("There should be more sequences than samples.\n")
#'   }
#'   
#'   #1. make a data frame
#'   m = data.frame(matrix(nrow = n, ncol = (ncol(a))))
#'   for(i in 1:ncol(m)){
#'     colnames(m)[i] = paste0(c("samp", i), sep = "", collapse = "")
#'   }
#'   
#'   #2. Get n collectors curves
#'   for(i in 1:n){
#'     
#'     if(i %% 10 == 0){
#'       print(i)
#'     }
#'     
#'     m[i,] = getCollectorCurve(a)
#'   }
#'   return(m)
#' }
#' 
#' 
#' # FUNCTION: getCollectorCurve
#' # PURPOSE: Returns a single collectors curve
#' # INPUTS
#' #   1. asv - an asv table where sequences are the rows and samples are the columns
#' # OUTPUTS
#' #   1. a vector that is a collectors curve. 
#' getCollectorCurve = function(a){
#'   
#'   #0. Do a couple of checks.
#'   
#'   if(nrow(a) < ncol(a)){
#'     stop("There should be more sequences than samples.\n")
#'   }
#'   
#'   #1. Make a matrix to put results in
#'   m = vector()
#'   
#'   #2. Sample the column names randomly
#'   rand = sample(colnames(a), replace = F)
#'   
#'   #3A. make a sequence vector
#'   seq = vector()
#'   for(i in 1:length(rand)){
#'     
#'     curSample = rand[i]
#'     seq = c(seq, rownames(a[which(a[,curSample] !=0),]))
#'     seq = unique(seq)
#'     
#'     m = c(m, length(seq))
#'   }
#'   return(m)
#' }
#' 
#' # FUNCTION: getCollectorCurve
#' # PURPOSE: Returns a single collectors curve
#' # INPUTS
#' #   1. asv - an asv table where sequences are the rows and samples are the columns
#' # OUTPUTS
#' #   1. a vector that is a collectors curve.
#' getCollectorCurveTaxonomy = function(a, t){
#'   
#'   #0. Do a couple of checks.
#'   if(nrow(a) < ncol(a)){
#'     stop("There should be more sequences than samples.\n")
#'   }
#'   
#'   #1. Make a matrix to put results in
#'   m = vector()
#'   
#'   #2. Sample the column names randomly
#'   rand = sample(colnames(a), replace = F)
#'   
#'   #3A. make a sequence vector
#'   seq = vector()
#'   taxTemp = vector()
#'   for(i in 1:length(rand)){
#'     
#'     curSample = rand[i]
#'     seq = c(seq, rownames(a[which(a[,curSample] !=0),]))
#'     curTax = t[seq,]
#'     curTax = as.vector(curTax[,6])
#'     taxTemp = c(taxTemp, curTax)
#'     taxTemp = unique(taxTemp)
#'     
#'     m = c(m, length(taxTemp))
#'   }
#'   return(m)
#' }
#' 
#' # getCollectorCurvesTaxonomyConfidenceInterval
#' # FUNCTION: getCollectorCurves
#' # PURPOSE: Returns n collectors curves
#' # INPUTS: 
#' #   1. n - an integer stating how many collectors curves to get
#' #   2. a - an asv table where sequences are the rows and samples are the columns
#' # OUTPUTS:
#' #   1. a matrix where row i is a single collectors curve in full, and columns are points in the curve
#' getCollectorCurvesTaxonomyConfidenceInterval = function(n, t, a){
#'   
#'   #0. Do a couple of checks.
#'   if(nrow(a) < ncol(a)){
#'     stop("There should be more sequences than samples.\n")
#'   }
#'   
#'   #1. make a data frame
#'   m = data.frame(matrix(nrow = n, ncol = (ncol(a))))
#'   for(i in 1:ncol(m)){
#'     colnames(m)[i] = paste0(c("samp", i), sep = "", collapse = "")
#'   }
#'   
#'   #2. Get n collectors curves
#'   for(i in 1:n){
#'     
#'     if(i %% 10 == 0){
#'       print(i)
#'     }
#'     
#'     m[i,] = getCollectorCurveTaxonomy(a = a, t = t)
#'   }
#'   return(m)
#' }
#' 
#' ## FUNCTION: writeFastQ
#' ## PURPOSE: To convert a sequence table to a fastQ format
#' ## INPUTS: 
#' ##   t: a sequence key file with the following requirements
#' ##     Col 1: seqID the sequence ID (short version)
#' ##     Col 2: sequence - The sequence itself (dada2 output)
#' ##   file: the file name that you want to write the fastQ to
#' ## OUTPUTS:
#' ##   fastQ file written within your current working directory with the name provided to file
#' writeFastQ = function(t, file){
#'   sink(file)
#'   for(i in 1:nrow(t)){
#'     cur = paste0(c(">", t[i,"seqID"], "\n",
#'                    t[i,"sequence"], "\n"), sep= "", collapse = "")
#'     cat(cur)
#'   }
#'   closeAllConnections()
#' }
#' 
#' 
#' ##FUNCTION: getTaxString
#' ## PURPOSE: To convert a dada2 formatted taxonomic table to a qiime formated
#' ## taxonomic table
#' ## INPUTS taxTable:
#' ## columnNames must be titled Kindom, Phylum, Class, Order, Family, Genus
#' ## rows correspond to each sequence.
#' ## An entry taxTable[i,j] corresponds to the jth taxonomic classifier
#' ## of sequence i.
#' ## OUTPUT taxReturn
#' ## A matrix with two colums
#' ## column1: the sequence ID from the rows in TaxTable
#' ## column2: the collapsed taxString in QIIME format
#' ## if the taxonomic entry in taxTable is NA, then the output string will
#' ## contain "Unclassified".
#' getTaxString = function(taxTable){
#'   taxReturn = as.data.frame(matrix(nrow = nrow(taxTable), ncol= 2))
#'   for (i in 1:nrow(taxTable)){
#'     string = ""
#'     #get kingdom
#'     if(!is.na(taxTable[i, "Kingdom"])){
#'       string = paste(c(string, "k__", toString(taxTable[i, "Kingdom"])), collapse = "")
#'     }else{
#'       #string = paste(c(string, " Unassigned;"), collapse = "")
#'     }
#'     #get Phylum
#'     if(!is.na(taxTable[i, "Phylum"])){
#'       string = paste(c(string, "; p__", toString(taxTable[i, "Phylum"])), collapse = "")
#'     }else{
#'       #string = paste(c(string, " Unassigned;"), collapse = "")
#'     }
#'     #Get Class
#'     if(!is.na(taxTable[i, "Class"])){
#'       string = paste(c(string, "; c__", toString(taxTable[i, "Class"])), collapse = "")
#'     }else{
#'       #string = paste(c(string, " Unassigned;"), collapse = "")
#'     }
#'     #Get Order
#'     if(!is.na(taxTable[i, "Order"])){
#'       string = paste(c(string, "; o__", toString(taxTable[i, "Order"])), collapse = "")
#'     }else{
#'       #string = paste(c(string, " Unassigned;"), collapse = "")
#'     }
#'     #Get Family
#'     if(!is.na(taxTable[i, "Family"])){
#'       string = paste(c(string, "; f__", toString(taxTable[i, "Family"])), collapse = "")
#'     }else{
#'       #string = paste(c(string, " Unassigned;"), collapse = "")
#'     }
#'     #Get Genus
#'     if(!is.na(taxTable[i, "Genus"])){
#'       string = paste(c(string, "; g__", toString(taxTable[i, "Genus"])), collapse = "")
#'     }else{
#'       #string = paste(c(string, " Unassigned;"), collapse = "")
#'     }
#'     if(!is.na(taxTable[i, "Kingdom"])){
#'       if(taxTable[i, "Kingdom"] == "No blast hit"){
#'         print("Unassigned")
#'         string = "Unassigned"
#'         taxReturn[i,1] = rownames(taxTable)[i]
#'         taxReturn[i,2] = string
#'       }
#'     }
#'     taxReturn[i,1] = rownames(taxTable)[i]
#'     taxReturn[i,2] = string
#'   }
#'   colnames(taxReturn) = c("sequence", "tax")
#'   return(taxReturn)
#' }
#' 
#' # FUNCTION: getCladalAttributes
#' # PURPOSE:  To merge standard claatu output files into a single attributes matrix given the 
#' #           path to those input files
#' # INPUTS:   size: path to size matrix in claatu. The size file has the following properties
#' #             col1: node names
#' #             col2: Clade size
#' #             other: no column names. Tab separated file.
#' #           nodes2tax: The output of nodes2tax claatu. File has following properties
#' #             col1: node names
#' #             col2: the phylogenetic string for each of the tips within the clade
#' #             other: Tab separated file. Follow OTU table QIIME standards for tax string.
#' #           nodeTax: The output of claatu's taxonomy function
#' #             col1: the node names
#' #             col2: the node's taxonomy
#' #           groupTest: The output of claatu's group test function
#' #             col1: the node names followed by a "_" and then the group the p test is for
#' #             col5: the p value associated for that node across that group
#' #           pTest: The output of claatu's p value function
#' #             col1: the node names
#' #             col5: the p value associated for that node across all samples. 
#' # OUTPUTS:  A matrix with the following columns with the following properties
#' #             rownames: nodes
#' #             col1 ("cladeSize"): the clade size from the "size" matrix
#' #             col2 ("pValueCladalConservationAll"): The p value of the conservation of the clade
#' #               accross all samples.
#' #             col3 ("qValueCladalConservationAll"): The p.adjust from column 2 using fdr correction
#' #             col4 - col i+3 (where i is the number of groups that is in the group test): 
#' #               the group p value for each of the groups in the group test. 
#' #             col i+4 - col i+i+3 (where i is the number of groups that is in the groupt test):
#' #               the fdr corrected p value for each of the groups in columns 4 - i+3.
#' #             col n - 1, where n is the number of columns in the matrix ("nodeTax"): 
#' #               The taxonomy for the node.
#' #             col n, where n is the number of columns in the matrix ("tipTaxLabels"):
#' #               a string containing the tip taxonomic labes for each of the tips in the node.
#' getCladalAttributes = function(size, nodes2tax, nodeTax, groupTest, pTest){
#'   
#'   print("Reading in files ....\n")
#'   # Read in clade size
#'   if(file.exists(size)){
#'     size = read.table(size)
#'     if(ncol(size) != 2){
#'       stop("ERROR: The clade size file doesn't have two columns. Are you sure you have the right file?\n")
#'     }else{
#'       rownames(size) = size[,1]
#'       colnames(size) = c("node", "cladeSize")
#'     }
#'   }else{
#'     stop("ERROR: That clade size stats file doesn't exist in your current work directory.")
#'   }
#'   
#'   
#'   # Read in nodes2tax Table
#'   if(file.exists(nodes2tax)){
#'     nodes2tax = read.table(nodes2tax, sep = "\t") 
#'     if(ncol(nodes2tax) != 2){
#'       stop("ERROR: The nodes2tax file does not contain two columns. Are you sure you have the right file?")
#'     }else{
#'       rownames(nodes2tax) = nodes2tax[, 1]
#'       colnames(nodes2tax) = c("node", "tipTaxLabels")
#'     }
#'   }else{
#'     stop("ERROR: The nodes2tax file does not exist. \n")
#'   }
#'   
#'   # Read in node taxonomy
#'   if(file.exists(nodeTax)){
#'     nodeTax = read.table(nodeTax)
#'     if(ncol(nodeTax) != 2){
#'       stop("ERROR: The nodeTax file does not contain two columns.")
#'     }else{
#'       rownames(nodeTax) = nodeTax[, 1]
#'       colnames(nodeTax) = c("node", "nodeTax")
#'     }
#'   }else{
#'     stop("ERROR: the nodeTax file does not exist.\n")
#'   }
#'   
#'   
#'   # Read in ptest file 
#'   if(file.exists(pTest)){
#'     pTest = read.table(pTest)
#'     if(ncol(pTest) != 6){
#'       stop("ERROR: The ptest should have 6 columns. Are you sure you have the right file?\n")
#'     }else{
#'       rownames(pTest) = pTest[,1]
#'       colnames(pTest) = c("node", "obsCladalConservationAll", "meanCladalConservationAll", "sdCladalConservationAll", "zScoreCladalConservationAll", "pValueCladalConservationAll")
#'     }
#'   }else{
#'     stop("ERROR: Cannot find the pTest file")
#'   }
#'   
#'   if(file.exists(groupTest)){
#'     groupTest = read.table(groupTest, row.names = 1)
#'     colnames(groupTest) = c("obsCladalConservationGroup", "meanCladalConservationGroup", "sdCladalConservationGroup", "zScoreCladalConservationGroup", "pValueCladalConservationGroup")
#'     groupTest = flattenGroupPtest(groupTest)
#'     groupQValueNames = vector()
#'     for(i in 1:ncol(groupTest)){
#'       cur = strsplit(x = colnames(groupTest)[i], split = "_")
#'       cur = cur[[1]][2]
#'       groupQValueNames = c(groupQValueNames, paste(c("groupConservationQValue", cur), collapse = "", sep = ""))
#'     }
#'   }else{
#'     stop("ERROR: Cannot find the groupTestFile")
#'   }
#'   
#'   # Create an attributes data matrix
#'   d = data.frame(matrix(nrow = nrow(size), ncol = 5 + 2*length(colnames(groupTest))))
#'   rownames(d) = rownames(size)
#'   colnames(d) = c("cladeSize", "pValueCladalConservationAll", "qValueCladalConservationAll", colnames(groupTest), groupQValueNames, "nodeTax", "tipTaxLabels")
#'   
#'   if(setequal(rownames(d), rownames(size))){
#'     d[ , "cladeSize"] = size[rownames(d),"cladeSize"]
#'   }
#'   
#'   if(setequal(rownames(d), rownames(nodeTax))){
#'     d[, "nodeTax"] = nodeTax[rownames(d), "nodeTax"]
#'   }
#'   
#'   if(setequal(rownames(d), rownames(nodes2tax))){
#'     d[, "tipTaxLabels"] = nodes2tax[rownames(d), "tipTaxLabels"]
#'   }
#'   
#'   if(setequal(rownames(d), rownames(pTest))){
#'     d[ , "pValueCladalConservationAll"] = pTest[rownames(d), "pValueCladalConservationAll"]
#'     d[ , "qValueCladalConservationAll"] = p.adjust(d$pValueCladalConservationAll)
#'   }
#'   
#'   if(setequal(rownames(d), rownames(groupTest))){ #check all nodes are the same
#'     
#'     for(i in 1:length(colnames(groupTest))){
#'       
#'       d[ ,colnames(groupTest)[i]] = groupTest[rownames(d),colnames(groupTest)[i]]
#'       d[ ,groupQValueNames[i]] = p.adjust(d[ , colnames(groupTest)[i]])
#'       
#'     }
#'   }
#'   return(d)
#' }
#' 
#' #FUNCTION: flattenGroupPtest
#' #PURPOSE: To flatten the claatu output of the group test. Make the group test so rows are 
#' # unique node IDs
#' #INPUT: 
#' # g: Claatu group p test matrix where the rownames are in the format NODEID_GROUPID
#' #OUTPUT: A matrix with the following properties
#' # rows: nodeNames 
#' # columns: p vector for each group 
#' flattenGroupPtest = function(g){
#'   
#'   #Split names and groups
#'   names = strsplit(x = rownames(g), split = "_")
#'   
#'   # Get unique node names and unique groups to name the matrix
#'   print("Setting up matrix ... ")
#'   nodes = vector()
#'   groups = vector()
#'   
#'   pb = txtProgressBar(min = 0, max = nrow(g), style = 3)
#'   for(i in 1:nrow(g)){
#'     nodes = c(nodes, names[[i]][1])
#'     groups = c(groups, names[[i]][2])
#'     setTxtProgressBar(pb, i)
#'   }
#'   nodes = unique(nodes)
#'   groups = unique(groups)
#'   
#'   # make a matrix
#'   m = data.frame(matrix(nrow = length(nodes), ncol = length(groups)))
#'   rownames(m) = nodes
#'   colnames(m) = groups
#'   
#'   print(c("Flattening group matrix to: ", length(nodes), " number of clades."), sep = "", collapse = "")
#'   print(c("Flattening group matrix to: ", length(groups), " number of populations."), sep = "", collapse = "")
#'   
#'   print("Populating the matrix... ")
#'   for(i in 1:nrow(g)){
#'     curNode = strsplit(x = rownames(g)[i], split = "_")[[1]][1]
#'     curPop = strsplit(x = rownames(g)[i], split = "_")[[1]][2]
#'     m[curNode, curPop] = g[i,5]
#'   }
#'   
#'   #rename the populations with "groupPop"
#'   for(i in 1:ncol(m)){
#'     colnames(m)[i] = paste0(c("groupConservationPValue_", colnames(m)[i]), sep = "", collapse = "")
#'   }
#'   return(m)
#'   
#'   
#' }
#' 
#' #MueggeFunctions.R
#' # August 31, 2016
#' 
#' # DESCRIPTION:  Produces a pcoa plot from GetPCoA() in ggplot2.
#' # INPUTS:       The title of data (dataname) produced from GetPCoA()
#' #               The plot title (title)
#' #               The x-axis label (x)
#' #               The y-axis label (y)
#' #               The legend key label (keylab)
#' # OUTPUTS:      A plot to ggplot2 of PCoA
#' plotPCOA = function(dataName, title, x, y, keylab, colors){
#'   
#'   p = 
#'     ggplot(data = dataName$PCOA, aes(PC1, PC2)) + 
#'     ggtitle(title) + 
#'     geom_point(aes(color = colorgroup), size = 6) +
#'     geom_path(data = dataName$ellipses,
#'               aes(x = Dim1, y = Dim2, color = group), size = 1, linetype = 2) +
#'     theme(axis.title = element_text(size = 20),
#'           plot.title = element_text(size = 20),
#'           legend.text = element_text(size = 20),
#'           legend.title = element_text(size = 20)) + 
#'     labs(x = x, y = y) +
#'     scale_color_manual(name = keylab, values = colors)
#'   return(p)
#' }
#' 
#' # DESCRIPTION:  GetPCoA() - Produces PCOA with ellipse calculated in k = 2.
#' # INPUTS:       A method (meth) for calculating distances ("bray")
#' #               A T/F value (bin) indicating if distance should be calculated
#' #                 using presence / absence (T), or abundance (F).
#' #               A community data matrix (comm) for which the distane matrix 
#' #                 should be calculated
#' #               The number (k) of PC's that should be returned for the PCoA.
#' #               A mapping file (m)
#' #               The trait in the mapping file (t) that PCoA confidence 
#' #                 intervals should be calculated for.
#' # OUTPUTS:      A data frame with two components
#' #                 (1) a PCOA data frame which has PCs, grouplabels, a
#' #                 group for calculating elipses (must nave n>3), the 
#' #                 rownames, and a color group for how to color that
#' #                 group (everything gets a different color even for n<=3)
#' GetPCoA = function(meth, bin, comm, k, m, t) {
#'   # Make sure the trait you are looking for is in the mapping file.
#'   if ( ! (t %in% colnames(m))){
#'     stop("ERROR: The trait you are looking for is not in the mapping file.")
#'   }
#'   
#'   # Calculate distance matrix
#'   dist = vegdist(x = t(comm), method = meth, binary = bin)
#'   
#'   # Create PCOA object
#'   dist.pcoa = cmdscale(d = dist, k = k, eig = T) 
#'   
#'   # Get PC's to be plotted and name columns
#'   pcs = dist.pcoa$points
#'   pcsColnames = vector(length = k)
#'   for(i in 1:k){
#'     tmp = paste(c("PC", i), sep = "", collapse = "")
#'     pcsColnames[i] = tmp
#'   }
#'   colnames(pcs) = pcsColnames
#'   rm(pcsColnames)
#'   
#'   # Get a list of groups to calculate confidence intervals for. 
#'   # There must me more than 3 members in a group.
#'   groups = table(m[,t])
#'   groups = as.data.frame(groups[which(groups>3)])$Var1
#'   
#'   # This sets a vector equal to groups if group has n>3 and "Other" otherwise.
#'   levelList = vector(length = nrow(m))
#'   for(i in 1:length(levelList)){
#'     if(as.vector(m[i,t]) %in% groups){
#'       levelList[i] = as.vector(m[i,t])
#'     } else {
#'       levelList[i] = "Other"
#'     }
#'   }
#'   
#'   # Make a pcoa dataframe with the group to draw confidence intervals.
#'   pcoa = data.frame(pcs, 
#'                     group = levelList,
#'                     rownames = rownames(m),
#'                     colorgroup = m[,t])
#'   
#'   # Get confidence intervals for any group with n>3
#'   e = data.frame()
#'   ord = 
#'     ordiellipse(ord = dist.pcoa, 
#'                 groups = levelList, 
#'                 kind = "se", 
#'                 conf = .95, 
#'                 lwd = 2, 
#'                 draw = "lines")
#'   for(g in levels(groups)){
#'     if(g %in% groups){
#'       e = rbind(e, 
#'                 cbind(as.data.frame(
#'                   with(pcoa[pcoa$group == g,], vegan:::veganCovEllipse(
#'                     ord[[g]]$cov, ord[[g]]$center, ord[[g]]$scale))), 
#'                   group = g))
#'     }
#'     
#'   }
#'   
#'   # Name 
#'   #colnames(e) = c("PC1", "PC2", "group")
#'   return(list("ellipses" = e, "PCOA" = pcoa))
#'   
#'   
#'   
#' }
#' 
#' 
#' 
#' ### FUNCTION: getROCCurveAndConfidenceInterval
#' ### PURPOSE: Returns several characteristics of RF predictor for M iterations
#' ### INPUTS
#' ###  M: The number of times a random forest should be built 
#' ###  data: data matrix which contains all features and one predictor column
#' ###  predictor: the column name of the value that is trying to be predicted
#' ###  ntree: the number of trees in the forest to be grown for each of the M iterations.
#' ###  name: the name of the group for the ROC curve for graphing later.
#' ### OUTPUTS: A list which has the following components
#' ###  1. TPR - M rows and 1000 columsn which has TPR rates for M iterations
#' ###  2. FPR - M rowas and 1000 columns which has FPR rates for M iterations
#' ###  3. AUC - M AUCs for each iteration
#' ###  4. Imp - An M x N matrix of importance scores for each of the M iterations RF. 
#' getROCCurveAndConfidenceIntervals = function(M, data, predictor, ntree, name){
#'   
#'   N = ncol(data) - 1
#'   names = colnames(data)
#'   names = names[-which(names == predictor)]
#'   
#'   sensitivitiesMatrix = as.data.frame(matrix(nrow = M, ncol = 1000, data = NA))
#'   
#'   specificitiesMatrix = as.data.frame(matrix(nrow = M, ncol = 1000, data = NA))
#'   
#'   AUC = vector(length = M)
#'   
#'   importanceScores = as.data.frame(matrix(nrow = M, ncol = N, data = NA))
#'   colnames(importanceScores) = names
#'   
#'   
#'   for(i in 1:M){
#'     print(i)
#'     formula = as.formula(paste(predictor, "~", paste(c(names), collapse = "+")))
#'     RF = randomForest(formula, data = data, importance = T, ntree = ntree)
#'     curImp = RF$importance[,"MeanDecreaseAccuracy"]
#'     importanceScores[i,names] = curImp[names]
#'     
#'     roc = roc(RF$y, RF$votes[,1], levels = RF$classes)
#'     
#'     sensitivitiesMatrix[i,1:length(roc$sensitivities)] = roc$sensitivities 
#'     specificitiesMatrix[i,1:length(roc$specificities)] = 1 - roc$specificities
#'     
#'     AUC[i] = roc$auc
#'     
#'   }
#'   
#'   aveTPR = apply(sensitivitiesMatrix, 2, ave)[1,]
#'   aveFPR = apply(specificitiesMatrix, 2, ave)[1,]
#'   sdTPR = apply(sensitivitiesMatrix, 2, sd)
#'   li = aveTPR - 2*sdTPR
#'   ui = aveTPR + 2*sdTPR
#'   rocGraph = as.data.frame(cbind(aveTPR, aveFPR, sdTPR, li, ui))
#'   rocGraph = rocGraph[complete.cases(rocGraph),]
#'   rocGraph$name = rep(x = name, times = nrow(rocGraph))
#'   rocGraph = rbind(rocGraph, c(0, 0, 0, 0, 0, name))
#'   rfResultsList = list("TPR" = sensitivitiesMatrix, "FPR" = specificitiesMatrix, "AUC" = AUC, "Imp" = importanceScores, "rocGraph" = rocGraph)
#'   return(rfResultsList)
#' }
#' 
#' 
#' 
#' ### FUNCTION: getROCCurve
#' ### PURPOSE: Returns several characteristics of RF predictor for M iterations
#' ### INPUTS
#' ###  M: The number of times a random forest should be built 
#' ###  data: data matrix which contains all features and one predictor column
#' ###  predictor: the column name of the value that is trying to be predicted
#' ###  ntree: the number of trees in the forest to be grown for each of the M iterations.
#' ###  name: the name of the group for the ROC curve for graphing later.
#' ### OUTPUTS: A list which has the following components
#' ###  1. TPR - M rows and 1000 columsn which has TPR rates for M iterations
#' ###  2. FPR - M rowas and 1000 columns which has FPR rates for M iterations
#' ###  3. AUC - M AUCs for each iteration
#' ###  4. Imp - An M x N matrix of importance scores for each of the M iterations RF. 
#' getROCCurve= function(M, data, predictor, ntree, name){
#'   
#'   N = ncol(data) - 1
#'   names = colnames(data)
#'   names = names[-which(names == predictor)]
#'   
#'   sensitivitiesMatrix = as.data.frame(matrix(nrow = M, ncol = 1000, data = NA))
#'   
#'   specificitiesMatrix = as.data.frame(matrix(nrow = M, ncol = 1000, data = NA))
#'   
#'   AUC = vector(length = M)
#'   
#'   importanceScores = as.data.frame(matrix(nrow = M, ncol = N, data = NA))
#'   colnames(importanceScores) = names
#'   
#'   
#'   for(i in 1:M){
#'     print(i)
#'     formula = as.formula(paste(predictor, "~", paste(c(names), collapse = "+")))
#'     RF = randomForest(formula, data = data, importance = T, ntree = ntree)
#'     
#'     
#'     roc = roc(RF$y, RF$votes[,1], levels = RF$classes)
#'     
#'     sensitivitiesMatrix[i,1:length(roc$sensitivities)] = roc$sensitivities 
#'     specificitiesMatrix[i,1:length(roc$specificities)] = 1 - roc$specificities
#'     
#'     AUC[i] = roc$auc
#'     
#'   }
#'   
#'   aveTPR = apply(sensitivitiesMatrix, 2, ave)[1,]
#'   aveFPR = apply(specificitiesMatrix, 2, ave)[1,]
#'   sdTPR = apply(sensitivitiesMatrix, 2, sd)
#'   li = aveTPR - 2*sdTPR
#'   ui = aveTPR + 2*sdTPR
#'   rocGraph = as.data.frame(cbind(aveTPR, aveFPR, sdTPR, li, ui))
#'   rocGraph = rocGraph[complete.cases(rocGraph),]
#'   rocGraph$name = rep(x = name, times = nrow(rocGraph))
#'   rocGraph = rbind(rocGraph, c(0, 0, 0, 0, 0, name))
#'   rfResultsList = list("TPR" = sensitivitiesMatrix, "FPR" = specificitiesMatrix, "AUC" = AUC, "rocGraph" = rocGraph)
#'   return(rfResultsList)
#' }
#' 
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

#' plotDBRDAGut = function(cap, meta, ID = "ID", colorGroup, shapeGroup, colorLab, shapeLab, title = "") {
#'   cap.data = as.data.frame(scores(cap, display = "sites"))
#'   cap.data$ID = rownames(cap.data)
#'   cap.data = cap.data %>% left_join(meta, by = ID)
#'   
#'   ggplot(cap.data, aes(x = CAP1, y = MDS1)) + 
#'     geom_point(aes(color = get(colorGroup), shape = get(shapeGroup)), size = 2, alpha = 0.6) +
#'     labs(color = colorLab) + 
#'     labs(shape = shapeLab) + 
#'     theme(axis.title = element_text(size = 16), legend.text = element_text(size = 16), legend.title = element_text(size = 16)) +
#'     guides(color = guide_legend(override.aes = list(size=5))) +   
#'     guides(shape = guide_legend(override.aes = list(size=5))) +
#'     ggtitle(title)
#'   
#' }
#' 
#' # Beta diversity contrained by one variable
#' plotDBRDAGut1Constraint = function(cap, meta, ID = "ID", colorGroup, shapeGroup, sizeGroup, sizeLab, colorLab, shapeLab, title = "", sizeValues, shapeValues, colValues = c("#b7410e", "steelblue")) {
#'   cap.data = as.data.frame(scores(cap, display = "sites"))
#'   cap.data$ID = rownames(cap.data)
#'   cap.data = cap.data %>% left_join(meta, by = ID)
#'   
#'   p = ggplot(cap.data, aes(x = CAP1, y = MDS1)) + 
#'     geom_point(aes(color = get(colorGroup), size = as.character(get(sizeGroup))), alpha = 0.6) +
#'     geom_point(aes(shape = get(shapeGroup)), color = "black", size = 4, alpha = 1) +
#'     labs(color = colorLab) + 
#'     labs(shape = shapeLab) + 
#'     labs(size = sizeLab) +
#'     theme(axis.title = element_text(size = 16), legend.text = element_text(size = 16), legend.title = element_text(size = 16), title = element_text(size = 16)) +
#'     guides(color = guide_legend(override.aes = list(size=5))) +   
#'     guides(shape = guide_legend(override.aes = list(size=5))) +
#'     scale_shape_manual(values = shapeValues) + 
#'     scale_color_manual(values = colValues) + 
#'     scale_size_manual(values = sizeValues)+
#'     ggtitle(title)
#'   return(p)
#' }
#' 
#' # Beta diversity contrained by two variables
#' plotDBRDAGut2Constraint = function(cap, meta, ID = "ID", colorGroup, shapeGroup, sizeGroup, sizeLab, colorLab, shapeLab, title = "", sizeValues, shapeValues, colValues = c("#b7410e", "steelblue")) {
#'   cap.data = as.data.frame(scores(cap, display = "sites"))
#'   cap.data$ID = rownames(cap.data)
#'   cap.data = cap.data %>% left_join(meta, by = ID)
#'   
#'   p = ggplot(cap.data, aes(x = CAP1, y = CAP2)) + 
#'     geom_point(aes(color = get(colorGroup), size = as.character(get(sizeGroup))), alpha = 0.6) +
#'     geom_point(aes(shape = get(shapeGroup)), color = "black", size = 4, alpha = 1) +
#'     labs(color = colorLab) + 
#'     labs(shape = shapeLab) + 
#'     labs(size = sizeLab) +
#'     theme(axis.title = element_text(size = 16), legend.text = element_text(size = 16), legend.title = element_text(size = 16), title = element_text(size = 16)) +
#'     guides(color = guide_legend(override.aes = list(size=5))) +   
#'     guides(shape = guide_legend(override.aes = list(size=5))) +
#'     scale_shape_manual(values = shapeValues) + 
#'     scale_color_manual(values = colValues) + 
#'     scale_size_manual(values = sizeValues)+
#'     ggtitle(title)
#'   return(p)
#' }
#' 
#' plotDBRDA = function(cap, meta, ID = "ID", colorGroup, shapeGroup, sizeGroup, sizeLab, colorLab, shapeLab, title = "", sizeValues, shapeValues, colValues = c("#b7410e", "steelblue")) {
#'   cap.data = as.data.frame(scores(cap, display = "sites"))
#'   cap.data$ID = rownames(cap.data)
#'   cap.data = cap.data %>% left_join(meta, by = ID)
#'   
#'   p = ggplot(cap.data, aes(x = CAP1, y = MDS1)) + 
#'     geom_point(aes(color = get(colorGroup), size = as.character(get(sizeGroup))), alpha = 0.6) +
#'     geom_point(aes(shape = get(shapeGroup)), color = "black", size = 2, alpha = 1) +
#'     labs(color = colorLab) + 
#'     labs(shape = shapeLab) + 
#'     labs(size = sizeLab) +
#'     theme(axis.title = element_text(size = 16), legend.text = element_text(size = 16), legend.title = element_text(size = 16)) +
#'     guides(color = guide_legend(override.aes = list(size=5))) +   
#'     guides(shape = guide_legend(override.aes = list(size=5))) +
#'     scale_shape_manual(values = shapeValues) + 
#'     scale_color_manual(values = colValues) + 
#'     scale_size_manual(values = sizeValues)+
#'     ggtitle(title)
#'   return(p)
#' }
#' 
#' # Makes GLMS with base variables alone
#' # Tests GLM with base variables + microbial feature
#' # Returns data frame with results of ANOVA (baseModel, microbeModel)
#' getGLMS = function(covariates, phyloseq, prevalenceThresh, baseVariables){
#'   prevDF = apply(X = otu_table(phyloseq),
#'                  MARGIN = 1,
#'                  FUN = function(x){sum(x > 0)})
#'   
#'   curTaxa = names(prevDF)[which(prevDF >= length(sample_names(phyloseq))*prevalenceThresh)]
#'   if(length(curTaxa) == 0){
#'     stop("That prevalence filter got rid of all taxa")
#'   }else{
#'     print("That prevalence threshold left this many taxa:")
#'     print(length(curTaxa))
#'   }
#'   phyloseq = prune_taxa(taxa = curTaxa, x = phyloseq)
#'   
#'   curData = as.data.frame(getModelingTable(phyloseq))
#'   curMicrobes = taxa_names(phyloseq)
#'   
#'   d = data.frame(matrix(nrow = length(covariates) * length(curMicrobes), ncol = 5))
#'   colnames(d) = c("covariate", "clade", "R", "p", "q")
#'   
#'   count = 0
#'   for(i in 1:length(covariates)){
#'     
#'     
#'     
#'     curCovs = covariates[i]
#'     print("Testing Covariate: ")
#'     print(curCovs)
#'     
#'     for(j in 1:length(curMicrobes)){
#'       
#'       
#'       count = count + 1
#'       #curFormula = paste(c(baseVariables, curMicrobes[j]), sep = "+", collapse = "+")
#'       
#'       #curFormula = as.formula(paste(c(curCovs, "~", curFormula), sep = "", collapse = ""))
#'       curFormula = as.formula(paste(c(curCovs, "~", curMicrobes[j]), sep = "", collapse = ""))
#'       curModel = glm(curFormula, family = binomial(link = "logit"), data = curData, na.action = "na.omit")
#'       #curModel = glm.nb(curFormula, family = binomial(link = "logit"), data = curData, na.action = "na.omit")
#'       
#'       
#'       d[count, "covariate"] = curCovs
#'       d[count, "clade"] = curMicrobes[j]
#'       if(curMicrobes[j] %in% rownames(summary(curModel)$coefficients) ){
#'         d[count, "R"] = summary(curModel)$coefficients[curMicrobes[j],"Estimate"]
#'         d[count, "p"] = summary(curModel)$coefficients[curMicrobes[j],"Pr(>|z|)"]
#'       }else{
#'         d[count, "R"] = NA
#'         d[count, "p"] = NA
#'       }
#'       
#'     }
#'     
#'     d[which(d$covariate == curCovs & !is.na(d$p)),"q"] = p.adjust(p = as.vector(d[which(d$covariate == curCovs & !is.na(d$p)),"p"]), method = "fdr")
#'   }
#'   
#'   return(d)
#'   
#' }

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

#' 
#' ###
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

#' 
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

#' 
#' # DESCRIPTION   This function takes a phylogeny and a list of nodeNames of a
#' #               tree and returns a list of nodeNames. 
#' # INPUTS        tree - a phylogenetic tree containing all of nodenames
#' #               nodeNames - a list of nodeNames to be converted to nodeIDs
#' # OUTPUT        a vector of node IDs that are of type integer.
#' GetNodeIDs = function(tree, nodeNames){
#'   ids = vector()
#'   for(i in 1:length(nodeNames)){
#'     ids = c(ids, GetNodeID(tree = tree, nodeName = nodeNames[i]))
#'   }
#'   return(ids)
#' }
#' 

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

#' # Returns distances from a combination vector
#' getDistanceCombos = function(combos, distMatrix){
#'   
#'   distance = vector()
#'   for(i in 1:ncol(combos)){
#'     curDist = distMatrix[combos[1,i], combos[2,i]]
#'     distance = c(distance, curDist)
#'   }
#'   return(distance)
#'   
#' }
#' 
#' # Get random distance
#' distRandom = function(dist, n){
#'   
#'   distance = vector()
#'   for(i in 1:n){
#'     cur = sample(x = 1:dim(dist)[1], 2, replace = FALSE)
#'     distance = c(distance, dist[cur[1], cur[2]])
#'   }
#'   
#'   return(distance)
#'   
#' }
#' 
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
#' 
#' tryGLMNBModelZeroFunc = function(curFormula, curData){
#'   
#'   out <- tryCatch(
#'     {
#'       curModel = zeroinfl(formula = curFormula, data = curData, dist = "negbin", link = "log", na.action = "na.omit")
#'       
#'       #curModel = glm.nb(formula = curFormula, data = curData, na.action = "na.omit", link = "log")# try glm.nb first
#'       
#'       
#'     },
#'     error=function(cond) {
#'       #print("Error with:")
#'       #print(curFormula)
#'       #message(paste("URL does not seem to exist:", url))
#'       #message("Here's the original error message:")
#'       #message(cond)
#'       # Choose a return value in case of error
#'       return(NA)
#'     },
#'     #warning=function(cond) {
#'     #message(paste("URL caused a warning:", url))
#'     #message("Here's the original warning message:")
#'     #message(cond)
#'     # Choose a return value in case of warning
#'     #  return(NA)
#'     #},
#'     finally={
#'       
#'     }
#'   )    
#'   return(out)
#' }
#' 
#' 
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
#' 
#' #' @title Remove samples from phyloseq object that have less than n taxa
#' #'
#' #' @param physeq A phyloseq-class object
#' #' @param mintaxa Minimum number of taxa that should be present in a sample (default, 10)
#' #'
#' #' @return Trimmed phyloseq object (All samples will have >= N taxa)
#' #' @export
#' #'
#' #' @examples
#' #' data("esophagus")
#' #' esophagus
#' #' phyloseq_richness_filter(esophagus, mintaxa = 30)
#' #' phyloseq_richness_filter(esophagus, mintaxa = 100)
#' #'
#' phyloseq_richness_filter <- function(physeq, mintaxa = 10){
#'   
#'   ## Estimate number of OTUs per sample
#'   sp <- phyloseq::estimate_richness(physeq, measures = "Observed")
#'   samples_to_keep <- rownames(sp)[ which(sp$Observed >= mintaxa) ]
#'   
#'   if(length(samples_to_keep) == 0){
#'     stop("All samples will be removed.\n")
#'   }
#'   
#'   if(length(samples_to_keep) == phyloseq::nsamples(physeq)){
#'     cat("All samples will be preserved\n")
#'     res <- physeq
#'   }
#'   
#'   if(length(samples_to_keep) < phyloseq::nsamples(physeq)){
#'     res <- phyloseq::prune_samples(samples = samples_to_keep, x = physeq)
#'   }
#'   
#'   return(res)
#' }
#' 
#' 
#' 
#' #' @title Remove taxa with small mean relative abundance.
#' #'
#' #' @param physeq A phyloseq-class object
#' #' @param frac The minimum cutoff for the relative OTU abundance
#' #' @details This function searches for taxa with small mean relative abundance and removes them. Result will be returned with original counts in the abundance table.
#' #' @return Phyloseq object with a subset of taxa.
#' #' @export
#' #'
#' #' @examples
#' #' data("esophagus")
#' #' phyloseq_filter_taxa_rel_abund(esophagus, frac = 0.01)
#' #' phyloseq_filter_taxa_rel_abund(esophagus, frac = 0.1)
#' #'
#' phyloseq_filter_taxa_rel_abund <- function(physeq, frac = 1e-4){
#'   
#'   # require(phyloseq)
#'   
#'   ## Transform OTU counts to relative abundance
#'   rel <- phyloseq::transform_sample_counts(physeq, function(x) x / sum(x) )
#'   
#'   ## Filter OTUs
#'   rel.subs <- phyloseq::filter_taxa(rel, function(x){ mean(x) > frac }, prune = FALSE)
#'   
#'   ## if prune = TRUE
#'   # tn <- taxa_names(rel.subs)              # OTUs to preserve
#'   # tr <- setdiff(taxa_names(physeq), tn)   # OTUs to remove
#'   
#'   ## Taxa to remove
#'   tr <- names(rel.subs)[ which(rel.subs == FALSE) ]
#'   
#'   ## If all taxa should be removed
#'   if(length(tr) == phyloseq::ntaxa(physeq)){
#'     stop("Error: all taxa will be removed with the specified 'frac' cutoff.\n")
#'   }
#'   
#'   ## If there is nothing to remove
#'   if(length(tr) == 0){
#'     res <- physeq
#'     cat("Warning: no taxa removed.\n")
#'   }
#'   
#'   ## Keep taxa which satisfies the truncation threshold
#'   if(length(tr) > 0){
#'     res <- phyloseq::prune_taxa(taxa = rel.subs, physeq)
#'   }
#'   
#'   return(res)
#' }
#' 
#' 
#' 
#' 
#' #' @title Remove taxa with abundance less then a certain fraction of total abundance.
#' #'
#' #' @param physeq A phyloseq-class object
#' #' @param frac The minimum cutoff for the OTU abundance in the table. This number is a fraction, not a percent.
#' #' @details
#' #' If frac = 0.0001, this will retain all OTU's that have at least a 0.01\% total abundance in the OTU table.
#' #' If you wanted to retain OTUs with at least 1\% total abundance, you must specify, 0.01.
#' #'
#' #' @return Phyloseq object with a subset of taxa.
#' #' @export
#' #' @seealso http://qiime.org/scripts/filter_otus_from_otu_table.html
#' #' @examples
#' #' data("esophagus")
#' #' phyloseq_filter_taxa_tot_fraction(esophagus, frac = 0.01)
#' #'
#' phyloseq_filter_taxa_tot_fraction <- function(physeq, frac = 0.01){
#'   
#'   # require(phyloseq)
#'   
#'   ## Estimate total abundance of OTUs
#'   tot <- sum(phyloseq::taxa_sums(physeq))
#'   
#'   ## Remove OTUs
#'   res <- phyloseq::filter_taxa(physeq, function(x){ ( sum(x)/tot ) > frac }, prune = TRUE)
#'   return(res)
#' }
#' 
#' 
#' 
#' 
#' #' @title Filter low-prevalence OTUs.
#' #' @description This function will remove taxa (OTUs) with low prevalence, where prevalence is the fraction of total samples in which an OTU is observed.
#' #' @param physeq A phyloseq-class object
#' #' @param prev.trh Prevalence threshold (default, 0.05 = 5\% of samples)
#' #' @param abund.trh Abundance threshold (default, NULL)
#' #' @param threshold_condition Indicates type of prevalence and abundance conditions, can be "OR" (default) or "AND"
#' #' @param abund.type Character string indicating which type of OTU abundance to take into account for filtering ("total", "mean", or "median")
#' #' @details
#' #' Abundance threshold defines if the OTU should be preserved if its abundance is larger than threshold (e.g., >= 50 reads).
#' #' Parameter "threshold_condition" indicates whether OTU should be kept if it occurs in many samples AND/OR it has high abundance.
#' #' @return  Phyloseq object with a subset of taxa.
#' #' @seealso \code{\link{phyloseq_prevalence_plot}}
#' #' @export
#' #'
#' #' @examples
#' #' data(GlobalPatterns)
#' #' GlobalPatterns  # 19216 taxa
#' #'
#' #' # OTUs that are found in at least 5% of samples
#' #' phyloseq_filter_prevalence(GlobalPatterns, prev.trh = 0.05, abund.trh = NULL)  # 15389 taxa
#' #'
#' #' # The same, but if total OTU abundance is >= 10 reads it'll be preserved too
#' #' phyloseq_filter_prevalence(GlobalPatterns, prev.trh = 0.05, abund.trh = 10, threshold_condition = "OR")  # 15639 taxa
#' #'
#' #' # Include only taxa with more than 10 reads (on average) in at least 10% samples
#' #' phyloseq_filter_prevalence(GlobalPatterns, prev.trh = 0.1, abund.trh = 10, abund.type = "mean", threshold_condition = "AND")  # 4250 taxa
#' #'
#' phyloseq_filter_prevalence <- function(physeq, prev.trh = 0.05, abund.trh = NULL, threshold_condition = "OR", abund.type = "total"){
#'   
#'   ## Threshold validation
#'   if(prev.trh > 1 | prev.trh < 0){ stop("Prevalence threshold should be non-negative value in the range of [0, 1].\n") }
#'   if(!is.null(abund.trh)){ 
#'     if(abund.trh <= 0){ stop("Abundance threshold should be non-negative value larger 0.\n") }
#'   }
#'   
#'   ## Check for the low-prevalence species (compute the total and average prevalences of the features in each phylum)
#'   prevdf_smr <- function(prevdf){
#'     ddply(prevdf, "Phylum", function(df1){ data.frame(Average = mean(df1$Prevalence), Total = sum(df1$Prevalence))})
#'   }
#'   # prevdf_smr( prevalence(physeq) )
#'   
#'   ## Check the prevalence threshold
#'   # phyloseq_prevalence_plot(prevdf, physeq)
#'   
#'   ## Define prevalence threshold as % of total samples
#'   ## This function is located in 'phyloseq_prevalence_plot.R' file
#'   prevalenceThreshold <- prev.trh * phyloseq::nsamples(physeq)
#'   
#'   ## Calculate prevalence (number of samples with OTU) and OTU total abundance
#'   prevdf <- prevalence(physeq)
#'   
#'   ## Get the abundance type
#'   if(abund.type == "total") { prevdf$AbundFilt <- prevdf$TotalAbundance }
#'   if(abund.type == "mean")  { prevdf$AbundFilt <- prevdf$MeanAbundance }
#'   if(abund.type == "median"){ prevdf$AbundFilt <- prevdf$MedianAbundance }
#'   
#'   ## Which taxa to preserve
#'   if(is.null(abund.trh)) { tt <- prevdf$Prevalence >= prevalenceThreshold }
#'   if(!is.null(abund.trh)){
#'     ## Keep OTU if it either occurs in many samples OR it has high abundance
#'     if(threshold_condition == "OR"){
#'       tt <- (prevdf$Prevalence >= prevalenceThreshold | prevdf$AbundFilt >= abund.trh)
#'     }
#'     
#'     ## Keep OTU if it occurs in many samples AND it has high abundance
#'     if(threshold_condition == "AND"){
#'       tt <- (prevdf$Prevalence >= prevalenceThreshold & prevdf$AbundFilt >= abund.trh)
#'     }
#'   }
#'   
#'   ## Extract names for the taxa we whant to keep
#'   keepTaxa <- rownames(prevdf)[tt]
#'   
#'   ## Execute prevalence filter
#'   res <- phyloseq::prune_taxa(keepTaxa, physeq)
#'   return(res)
#' }
#' 
#' 
#' 
#' #' @title Filter rare OTUs based on minimum abundance threshold.
#' #' @description This function performs sample-wise OTU abundance trimming.
#' #' @param physeq A phyloseq-class object
#' #' @param minabund Abundance threshold (default, 10)
#' #' @param relabund Logical; perform trimming based on relative abundances (default, FALSE)
#' #' @param rm_zero_OTUs Logical, remove OTUs with zero total abundance
#' #' @details 
#' #' OTUs can be considered as rare if they comprise fewer than X (e.g., 10) sequences within a sample. 
#' #' This function is intented to censore OTU abundance (unsing an arbitrary threshold) on a sample-wise basis. 
#' #' 
#' #' Trimming can be performed based on relative abundances of OTUs within a sample (`relabund = TRUE`), but the orginal OTU count will be returned. 
#' #' For this purpose `minabund` parameter should be provided in a range of (0,1] (e.g., use `minabund = 0.1, relabund = TRUE` to remove OTUs with relative abundance < 10% in each sample).
#' #' 
#' #' @return Phyloseq object with a filtered data.
#' #' @export
#' #'
#' #' @examples
#' #' # Load data
#' #' data(GlobalPatterns)
#' #'
#' #' # Trim GlobalPatterns data (19216 taxa) by removing OTUs with less that 10 reads
#' #' GP1 <- phyloseq_filter_sample_wise_abund_trim(GlobalPatterns, minabund = 10) # 10605 taxa
#' #'
#' #' # Trim GlobalPatterns data by removing OTUs with relative abundance less than 1%
#' #' GP2 <- phyloseq_filter_sample_wise_abund_trim(GlobalPatterns, minabund = 0.01, relabund = TRUE) # 258 taxa
#' #'
#' #' # Compare raw and trimmed data
#' #' phyloseq_summary(GlobalPatterns, GP1, GP2, cols = c("GlobalPatterns", "Trimmed 10 reads", "Trimmed 1 percent"))
#' #'
#' phyloseq_filter_sample_wise_abund_trim <- function(physeq, minabund = 10, relabund = FALSE, rm_zero_OTUs = TRUE){
#'   
#'   ## Censore OTU abundance
#'   if(relabund == FALSE){     # trim based on absolute OTU counts
#'     
#'     res <- phyloseq::transform_sample_counts(physeq, function(OTU, ab = minabund){ ifelse(OTU <= ab,  0, OTU) })
#'     
#'   } else {                   # trim based on relative abundances within sample, but return original counts
#'     
#'     if(!minabund > 0 & minabund <= 1){
#'       stop("Error: for relative abundance trimmin 'minabund' should be in (0,1] interval.\n")
#'     }
#'     
#'     ## Convert data to relative abundances
#'     res <- phyloseq_standardize_otu_abundance(physeq, method = "total")
#'     
#'     ## Remove relative abundances less than the threshold value
#'     res <- phyloseq::transform_sample_counts(res, function(OTU, ab = minabund){ ifelse(OTU <= ab,  0, OTU) })
#'     
#'     ## Sample sums and data orientation
#'     smps <- phyloseq::sample_sums(physeq)
#'     if(phyloseq::taxa_are_rows(physeq) == TRUE){
#'       mar <- 2
#'     } else {
#'       mar <- 1
#'     }
#'     
#'     ## Convert back to counts by multiplying relative abundances by sample sums
#'     phyloseq::otu_table(res) <- phyloseq::otu_table(
#'       sweep(x = phyloseq::otu_table(res), MARGIN = mar, STATS = smps, FUN = `*`),
#'       taxa_are_rows = phyloseq::taxa_are_rows(physeq))
#'   }
#'   
#'   ## Remove zero-OTUs
#'   if(rm_zero_OTUs == TRUE){
#'     res <- phyloseq::prune_taxa(taxa_sums(res) > 0, res)
#'   }
#'   return(res)
#' }
#' 
#' 
#' 
#' #' @title Extract the most abundant taxa.
#' #' @param physeq A phyloseq-class object
#' #' @param perc Percentage of the most abundant taxa to retain
#' #' @param n Number of the most abundant taxa to retain (this argument will override perc argument)
#' #' @return Phyloseq object with a filtered data.
#' #' @export
#' #'
#' #' @examples
#' #'
#' phyloseq_filter_top_taxa <- function(physeq, perc = 10, n = NULL){
#'   
#'   ## Arguments validation
#'   if(perc <= 0 | perc > 100){ stop("Error: percentage should be in 1-100 range.\n") }
#'   
#'   ## Get total abundances for all taxa
#'   taxx <- sort(phyloseq::taxa_sums(physeq), decreasing = TRUE)
#'   
#'   ## Find how many taxa to preserve (if percentage is specified)
#'   if(is.null(n)){
#'     n <- phyloseq::ntaxa(physeq) * perc / 100
#'     n <- floor(n)
#'   }
#'   
#'   ## Extract names for the taxa that should be preserved
#'   keepTaxa <- names(taxx)[1:n]
#'   
#'   ## Extract this taxa
#'   physeq_pruned <- phyloseq::prune_taxa(keepTaxa, physeq)
#'   
#'   return(physeq_pruned)
#' }
#' 
#' 
#' 
#' #' @title Check the range of the top-taxa filtering values to determine the optimal threshold.
#' #' @description This function performs taxa filtering by retaining the most abundant taxa.
#' #' A range of abundance percentages (5 - 95\%) will be explored.
#' #' @param physeq A phyloseq-class object
#' #' @param show_plot Logical; if TRUE, shows the plot on screen
#' #' @return ggplot-object.
#' #' @export
#' #'
#' #' @examples
#' #'
#' phyloseq_filter_top_taxa_range <- function(physeq, show_plot = TRUE){
#'   percs <- seq(5, 95, 5)
#'   
#'   fr <- plyr::mlply(.data = data.frame(perc = percs), .fun = function(...){ phyloseq_filter_top_taxa(physeq, ...) })
#'   names(fr) <- percs
#'   
#'   fr_tab <- plyr::ldply(.data = fr, .fun = function(z){
#'     sz <- phyloseq::sample_sums(z)
#'     res <- data.frame(Sample = names(sz), Preserved = sz)
#'     return(res)
#'   })
#'   
#'   pp <- ggplot(data = fr_tab, aes(x = perc, y = Preserved, group = Sample)) +   # color = Sample
#'     geom_vline(xintercept=75, color="grey", linetype = "longdash") +
#'     geom_line() +
#'     geom_point() +
#'     labs(x = "Number of most abundant taxa retained, %", y = "Percentage of total sample abundance") +
#'     theme(legend.position = "none")
#'   
#'   if(show_plot == TRUE){ print(pp) }
#'   invisible(pp)
#' }
#' 
#' #' @title Compute prevalence and abundance summary of each species/OTU.
#' #'
#' #' @param physeq A phyloseq-class object
#' #' @param add_tax Logical, add taxonomy to the results
#' #'
#' #' @return Data frame
#' #'
#' #' @examples
#' #' data(esophagus)
#' #' prevalence(esophagus)
#' #'
#' #' data(GlobalPatterns)
#' #' head( prevalence(GlobalPatterns, add_tax = T) )
#' #'
#' prevalence <- function(physeq, add_tax = TRUE){
#'   
#'   ## Check if taxa are rows
#'   trows <- taxa_are_rows(physeq)
#'   
#'   ## Extract OTU table
#'   otutab <- as.data.frame(otu_table(physeq))
#'   
#'   ## Transpose OTU table (species should be arranged by rows)
#'   if(trows == FALSE){
#'     otutab <- t(otutab)
#'   }
#'   
#'   ## Estimate prevalence (number of samples with OTU present)
#'   prevdf <- apply(X = otutab,
#'                   # MARGIN = ifelse(trows, yes = 1, no = 2),  # for a non-transposed data
#'                   MARGIN = 1,
#'                   FUN = function(x){sum(x > 0)})
#'   
#'   ## Add total and average read counts per OTU
#'   prevdf <- data.frame(Prevalence = prevdf,
#'                        TotalAbundance = taxa_sums(physeq),
#'                        MeanAbundance = rowMeans(otutab),
#'                        MedianAbundance = apply(otutab, 1, median))
#'   
#'   ## Add taxonomy table
#'   if(add_tax == TRUE && !is.null(tax_table(physeq, errorIfNULL = F))){
#'     prevdf <- cbind(prevdf, tax_table(physeq))
#'   }
#'   return(prevdf)
#' }
#' 
