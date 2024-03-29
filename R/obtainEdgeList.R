obtainEdgeList <- function(genes, databases){
  conversion_type <- . <- converted_gene <- converted_id <- direction <- dest <- dest_type <- src <- src_type <- base_gene_dest <- base_id_dest <- base_gene_src <- base_id_src <- database <- NULL #Added to avoid data.table note in R CMD check
  genes <- setNames(gsub(".*:(.*)", "\\1", genes), gsub("(.*):.*", "\\1", genes))

  ## Remove metabolites IDs because cant use ":::" in CRAN. See github for implementation
  metabolites_ids   <- NULL
  genes_with_id     <- reshapeIDs(genes, metabolites_ids)
  
  # Creating database here so we only need to process once. Also doing after genes_with_id since this takes the longest. Want to fail as early as possible
  full_edgelist     <- stackDatabases(setdiff(databases, "ndex"))
  #NDEx is special, add it if desired
  if("ndex" %in% databases) { message("NDEx is only available in the development version of netgsa: https://github.com/mikehellstern/netgsa")}#full_edgelist     <- addNDEx(full_edgelist, genes)}
  
  # Obtaining unique ids used by all chosen databases. We will convert to these IDs.
  # E.g. if databases only use "ENTREZID", "CHEBI", and "UNIPROT", this return a list of list(org.Hs.eg.db = c("ENTREZID", "UNIPROT"), metabolites = "CHEBI")
  databases_ids     <- findDatabaseIDs(full_edgelist, metabolites_ids)
  
  # Convert proteins
  ids_to_convert_using_orgHsegdb <- genes_with_id[conversion_type == "org.Hs.eg.db"]
  spl_ids                        <- split(ids_to_convert_using_orgHsegdb, ids_to_convert_using_orgHsegdb[["id_type"]])
  converted_proteins             <- rbindlist(lapply(spl_ids, convertProteinGroups, databases_ids[["org.Hs.eg.db"]]), use.names = TRUE)
  
  # Convert metabolites
  metabs_ids            <- genes_with_id[conversion_type == "graphite_metabolites"]
  converted_metabolites <- if(nrow(metabs_ids) == 0) NULL
                           else                      rbindlist(Map(convert_single_metabolite, metabs_ids[["id_type"]], metabs_ids[["gene"]], list(databases_ids[["metabolites"]]), list(metabolites_ids)), use.names = TRUE)
  
  # All finalized_ids should be character
  finalized_ids <- rbindlist(list(converted_proteins, converted_metabolites), use.names = TRUE, fill = TRUE)

  #Merge twice separately
  src_subs             <- setnames(full_edgelist[finalized_ids, on = .(src = converted_gene,  src_type  = converted_id), nomatch = 0L], c("base_gene", "base_id"), c("base_gene_src",  "base_id_src"))
  full_subs            <- setnames(     src_subs[finalized_ids, on = .(dest = converted_gene, dest_type = converted_id), nomatch = 0L], c("base_gene", "base_id"), c("base_gene_dest", "base_id_dest"))
  
  # Make all directed and stack on
  full_edgelist_subs_dir  <- rbindlist(list(full_subs, full_subs[direction == "undirected", .(src = dest,                     src_type = dest_type,       dest = src,                     dest_type = src_type, 
                                                                                              base_gene_src = base_gene_dest, base_id_src = base_id_dest, base_gene_dest = base_gene_src, base_id_dest = base_id_src,
                                                                                              database = database)]), use.names = TRUE, fill = TRUE)
  # Take unique
  full_edgelist_subs_dir_uniq  <- unique(full_edgelist_subs_dir[, c("database", "base_gene_src", "base_id_src", "base_gene_dest", "base_id_dest"), with = FALSE])
  
  ## Lastly, we need to know which genes never appear in the full stack edgelist. If they dont appear then we likely dont have enough info
  ## to tell whether or not they should be 0 or 1.
  genes_not_in_dbs <- findGenesNotInDb(full_edgelist, finalized_ids)
  
  obtainedEdgeList <- list(edgelist = full_edgelist_subs_dir_uniq, genes_not_in_dbs = genes_not_in_dbs)
  class(obtainedEdgeList) <- "obtainedEdgeList"
  
  return(obtainedEdgeList)
}

# Helper functions -----------------------------------------------------------------

# Wrapper to stack databases of interest 
  stackDatabases <- function(databases){
    . <- database <- src <- src_type <- dest <- dest_type <- direction <- NULL #Added to avoid data.table note in R CMD check
    not_in_graphite <- ! databases %in% unique(as.character(graphite::pathwayDatabases()[["database"]]))
    if(any(not_in_graphite)) stop(paste0("Database(s): ", paste0(databases[not_in_graphite], collapse = ", "), 
                                                          " are not found within graphite. Please choose databases from graphite::pathwayDatabases()$database"))
    
    db_stacked <- rbindlist(lapply(databases, function(db){
                                              rbindlist(lapply(graphite::pathways("hsapiens", db), graphite::edges, which = "mixed", stringsAsFactors = FALSE), use.names = TRUE, idcol = "Pathway")}) 
                            %>% setNames(., databases), use.names = TRUE, idcol = "database")
    return(unique(db_stacked[,.(database, src, src_type, dest, dest_type, direction)]))
  }


# Reshape IDs into format that's easier to merge onto database
  reshapeIDs <- function(genes, metabolites){
    
    # Make sure all genes have names
    if(is.null(names(genes)) | any(names(genes) == "")) stop("Please select identfier type of all genes. You can do this by setting the name of each element to the identifier type of that gene")
    
    # Make sure their gene is gene type we know about
    unconv_ids <- !names(genes) %in% c(names(metabolites), AnnotationDbi::keytypes(org.Hs.eg.db::org.Hs.eg.db))
    
    if(any(unconv_ids)) stop(paste0(paste0("Don't know how to convert ID type(s): ", paste(names(genes)[unconv_ids], collapse = ", ")), 
                                           ". Please use ID type from AnnotationDbi::keytypes(org.Hs.eg.db::org.Hs.eg.db)."))
    
  
    
    # All ID types are valid, drop duped IDs and make sure actual genes are valid for specific ID type
    genes_cl <- checkIDs(genes, metabolites)
    conversion_type = ifelse(names(genes_cl) %in% names(metabolites),                                    "graphite_metabolites",
                      ifelse(names(genes_cl) %in% AnnotationDbi::keytypes(org.Hs.eg.db::org.Hs.eg.db),   "org.Hs.eg.db", NA))
    
    data.table(id_type = names(genes_cl), gene = genes_cl, conversion_type = conversion_type)
    
  }

# Checking to make sure IDs are valid. That is, make sure we can find the IDs in our databases
  checkIDs <- function(genes, metabolites){
    
    # Error if cannot find IDs in database
    spl_genes <- split(genes, names(genes))
  
    not_found_ids <- lapply(spl_genes, function(gene_subs, metabolites){
                                              id_type       <- names(gene_subs)[1]
                                              valid_id_list <- if(      id_type %in% names(metabolites)    )                                   metabolites[[id_type]]
                                                               else if( id_type %in% AnnotationDbi::keytypes(org.Hs.eg.db::org.Hs.eg.db))      AnnotationDbi::keys(org.Hs.eg.db::org.Hs.eg.db, id_type)
                                              
                                              not_found     <- gene_subs[! gene_subs %in% valid_id_list]
                                              if( length(not_found) == 0 )  return(NULL)
                                              else                          return(paste0(id_type, ": ", paste(not_found, collapse = ", ")))
                        }, metabolites = metabolites)
    
    error_ids <- Filter(function(x) !is.null(x), not_found_ids)
    
    if(length(error_ids) != 0){
      prefix <- "The following IDs were not found in the keylist and thus are not able to be converted: \n  "
      stop(paste0(prefix, do.call(paste0, list(error_ids, collapse = "\n  "))))
    }
    
    # Warn if duplicated ID and ID type
    duped_id <- duplicated(paste0(names(genes),":", genes))
    if(any(duped_id)) {stop(paste0("The following duplicate IDs were detected. Please remove and rerun: ", paste0(paste0(names(genes[duped_id]), " = ", genes[duped_id]), collapse = ", ")))}
    
    
    return(genes[!duped_id])
  }


# Finding set of ids used across chosen databases
  findDatabaseIDs <- function(databases, metabolites){
    src_type <- dest_type <- NULL #Added to avoid data.table note in R CMD check
    uniq_ids                   <- databases[, unique(c(src_type, dest_type))]
    
    # IDs are converted either with metabolites from graphite OR org.Hs.eg.db
    metabolite_or_orgHsegdb_id <- ifelse(uniq_ids %in% AnnotationDbi::keytypes(org.Hs.eg.db::org.Hs.eg.db), "org.Hs.eg.db",
                                  ifelse(uniq_ids %in% names(metabolites),                                  "metabolites", NA))
    return(split(uniq_ids, metabolite_or_orgHsegdb_id))
  }


# Function to convert one set of proteins with a common base ID to all the ids we want to convert to
  convertProteinGroups <- function(id_table, convert_to){
  
    # Should only ever be of length one b/c we split by id_type
    keytype                                         <- unique(id_table[["id_type"]])
    
    # If keytype == convert_to simply return original data reshaped / renamed
    if(all(keytype == convert_to))                     return(data.table(base_gene = id_table[["gene"]], converted_gene = id_table[["gene"]], base_id = id_table[["id_type"]], converted_id = id_table[["id_type"]]))
    
    # We select c(keytype, convert_to) so that we can select base_id twice if its in list we want to convert to. This is possible b/c data.table allows duplicate column names
    converted                                       <- setDT(suppressMessages(AnnotationDbi::select(org.Hs.eg.db::org.Hs.eg.db, keys = id_table[["gene"]], columns = convert_to, keytype = keytype)))[, c(keytype, convert_to), with = FALSE]
    converted                                       <- unique(data.table::melt(converted, id.vars = keytype, variable.name = "converted_id", value.name = "converted_gene", variable.factor = FALSE, value.factor = FALSE, na.rm = TRUE))
    converted[, c("base_id")]                       <- keytype
    names(converted)[names(converted) == keytype]   <- "base_gene"
    
    unconvertables                                  <- which(is.na(converted[["converted_gene"]]))
    if(length(unconvertables) != 0)                    warning(paste0("\nCould not convert the following genes :\n",
                                                                      paste0(paste0(keytype, " = ", converted[unconvertables, ][["base_gene"]], " to ", converted[unconvertables, ][["converted_id"]]), collapse = "\n")))
    
    return(converted)
  }


# Function to convert a single metabolite ID to all the ids we want to convert to 
  convert_single_metabolite <- function(base_id, base_gene, convert_to, metabolites){
  
    if(all(base_id == convert_to))                                           return(data.table(base_gene = base_gene, converted_id = base_id, converted_gene = base_gene, base_id = base_id))
    converted_ids_metabs                                                  <- unique(data.table::melt(metabolites[which(metabolites[[base_id]] == base_gene), c(base_id, convert_to), with = FALSE], id.vars = base_id, na.rm = TRUE, 
                                                                                                     variable.name = "converted_id", value.name = "converted_gene", variable.factor = FALSE, value.factor = FALSE))
  
    # Throwing warning about unconvertable IDs
    unconv_ids                                                            <- if(grepl("KEGG", base_id)) setdiff(convert_to[!grepl("KEGG", convert_to)], converted_ids_metabs[["converted_id"]])
                                                                             else                       setdiff(convert_to, converted_ids_metabs[["converted_id"]])
    if(length(unconv_ids != 0)){                                             warning(paste0("\nCould not convert gene :\n", base_id, " = ", base_gene, " to ", paste0(unconv_ids, collapse = ", ")))}
    
    if(nrow(converted_ids_metabs) == 0){                                     return(NULL)}
    
    names(converted_ids_metabs)[names(converted_ids_metabs) == base_id]   <- "base_gene"
    converted_ids_metabs[["base_id"]]                                     <- base_id
  
    return(converted_ids_metabs)
  }

#Given edgelist from stackDatabases and idlist from convertID functions, get the genes not in the database ever
  findGenesNotInDb <- function(edgelist, idlist){
    . <- converted_gene <- converted_id <- dest <- base_gene <- base_id <- src <- in_db <- any_in_db_src <- any_in_db_dest <- any_in_db <- NULL #Added to avoid data.table note in R CMD check
    src_m  <- edgelist[idlist, on = .(src  = converted_gene, src_type  = converted_id), .(in_db = max(!is.na(dest)),  base_gene, base_id), nomatch = NA, by = .EACHI]
    dest_m <- edgelist[idlist, on = .(dest = converted_gene, dest_type = converted_id), .(in_db = max(!is.na(src)), base_gene, base_id), nomatch = NA, by = .EACHI]
    
    src_m2  <- src_m[, .(any_in_db_src = max(in_db)), by= .(base_gene, base_id)]
    dest_m2 <- dest_m[, .(any_in_db_dest = max(in_db)), by= .(base_gene, base_id)]
    
    not_in_dbs <- src_m2[dest_m2, on = c("base_gene", "base_id"), .(any_in_db = pmax(any_in_db_src, any_in_db_dest)), by = .EACHI][any_in_db == 0]
    
    if(nrow(not_in_dbs) > 0) return(setNames(not_in_dbs[["base_gene"]], not_in_dbs[["base_id"]]))
    else                     return(NULL)
    
  }
  
## Removed code for compatibility with CRAN
#Add NDEx information
  #Does not do any conversion of genes
  #Assumes all edges are undirected
  #Only searches public
  # addNDEx <- function(edgelist, genes){
  #   externalId <- NULL
  #   ndexcon = ndexr::ndex_connect()
  #   #Find NDEx networks that have any of our genes as nodeNames. Only option appears to be string search in nodeNames 
  #   #(e.g. if gene is "38" this will match a nodeName of "3801" in NDEx). Also NDEx is by default an OR search, e.g.
  #   #if searchString="nodeName:38 nodeName:401" it will find networks that have partial matches in 38 OR 401
  #   #Do this in batches of genes since error if do all at once. Get unique network ids. This will be unique networks
  #   #that have a partial match to ANY of our genes of interest
  #   search_names <- paste0("nodeName:",genes)
  #   found_dbs <- lapply(split(search_names, 1:length(search_names) %/% 200), function(search_vec){
  #     return(ndexr::ndex_find_networks(ndexcon, searchString=paste0(search_vec, collapse = " "), size=-1))
  #   })
  #   all_found_dbs <- rbindlist(found_dbs, use.names = TRUE)
  #   unique_ids <- unique(all_found_dbs[,externalId])
  #   
  #   #Query subnetwork to get only genes of interest. Again this is partial search so will need to cut this down
  #   NDEx_edges <- rbindlist(lapply(unique_ids, function(id) NDEx_query_subnetwork(id, ndexcon=ndexcon, search_names=search_names, genes=genes)), 
  #                           use.names = TRUE)
  #   #Add to full_edgelist
  #   edgelist_w_ndex <- rbindlist(list(edgelist, NDEx_edges), use.names=TRUE)
  #   
  #   return(edgelist_w_ndex)
  # }
  # 
  # NDEx_query_subnetwork <- function(ndex_network_id, ndexcon, search_names, genes){
  #   src <- i.n <- dest <- . <- NULL #To get rid of data.table package build notes
  #   query <- jsonlite::toJSON(list(searchString = paste(search_names, collapse = " ")), auto_unbox = TRUE)
  #   ## Removed below code to get "response" for R package because ::: not allowed. See github for implementation
  #   response <- NULL
  #   
  #   #Get nodelist/edgelist from query
  #   node_dt <- as.data.table(response$nodes[[4]])
  #   edge_dt <- as.data.table(response$edges[[3]])
  #   
  #   if(nrow(node_dt)==0 | nrow(edge_dt) == 0) {return(NULL)}
  #   
  #   edge_mapped <- edge_dt[node_dt, on=c(s="@id"), src := i.n][node_dt, on=c(t="@id"), dest := i.n]
  #   #Keep only edges we're interested in
  #   edge_subs <- edge_mapped[((src %in% genes) & (dest %in% genes))]
  #   #Remove self-edges
  #   edge_subs <- edge_subs[src != dest]
  #   
  #   
  #   if(nrow(edge_subs) == 0){return(NULL)}
  #   
  #   #Format for stacking with full edgelist
  #   edge_format <- edge_subs[,.(database = paste0("ndex_", ndex_network_id), 
  #                               src, src_type = names(genes)[match(src, genes)], 
  #                               dest, dest_type = names(genes)[match(dest, genes)], 
  #                               direction = "undirected")]
  #   return(edge_format)
  # }




