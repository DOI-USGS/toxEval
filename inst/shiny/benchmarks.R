output$downloadBenchmarks <- downloadHandler(
  
  filename = "Benchmarks.csv",
  
  content = function(filename) {
    write.csv(get_benchmarks(), filename, row.names = FALSE)
  }
)

get_benchmarks <- reactive({
  
  groupCol <- epDF[["groupColName"]]
  assays <- epDF[["assays"]]
  flags <- epDF[["flags"]]
  sites <- epDF[["sites"]]
  groups <- epDF[["group"]]
  removeFlags <- all_flags[!(all_flags %in% flags)]
  rawData <- rawData()
  
  if(!is.null(rawData)){
    if(all(is.null(rawData$benchmarks))){
      
      ACClong <- get_ACC(rawData$chem_info$CAS)
      ACClong <- remove_flags(ACClong, flagsShort = removeFlags)
      
      remove_groups <- unique(cleaned_ep[[groupCol]])[which(!unique(cleaned_ep[[groupCol]]) %in% groups)]
      remove_groups <- remove_groups[!is.na(remove_groups)]
      
      filtered_ep <- filter_groups(cleaned_ep, 
                                   groupCol = groupCol, assays = assays,
                                   remove_groups = remove_groups)
      
      bench <- ACClong %>%
        filter(endPoint %in% filtered_ep$endPoint) %>%
        rename(Value = ACC_value,
               Chemical = chnm) %>%
        left_join(filtered_ep, by = "endPoint")

    } else {
      bench <- rawData$benchmarks
    }
  }

  return(bench)
})