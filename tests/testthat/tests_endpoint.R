context("Endpoint tests")

# This has some CAS's that are not in toxCast:
CAS <- c("3380-34-5","106-44-5","80-05-7","121-00-6","136-85-6",
         "84852-15-3","599-64-4","140-66-9","1806-26-4","20427-84-3",
         "104-35-8","2315-61-9","2315-67-5","84-65-1","13674-87-8",
         "126-73-8","78-51-3","5436-43-1","119-61-9","21145-77-7",
         "120-72-9","124-76-5","5989-27-5","119-65-3","76-22-2" ,
         "83-34-1","1222-05-5","90-12-0","91-57-6","581-42-0",
         "98-82-8","1912-24-9","314-40-9","57837-19-1","1610-18-0",
         "51218-45-2","87-86-5","102-36-3","58-08-2","486-56-6",
         "89-78-1","333-41-5","63-25-2","62-73-7","134-62-3",
         "2921-88-2","86-74-8","119-36-8","106-46-7","75-25-2",
         "91-20-3","129-00-0","206-44-0","85-01-8","120-12-7",
         "50-32-8","115-86-6","84-66-2","115-96-8","78-59-1",
         "127-18-4","83-46-5","57-88-5","360-68-9","19466-47-8")

test_that("Getting ACC values", {
 testthat::skip_on_cran()
 
 ACClong <- get_ACC(CAS)
 expect_true(all(names(ACClong) %in% c("casn","chnm","flags","endPoint","ACC","MlWt", "ACC_value")))
   
 expect_type(ACClong$MlWt, "double")
 expect_type(ACClong$ACC_value, "double")
 expect_type(ACClong$endPoint, "character")
 expect_type(ACClong$casn, "character")
 
 expect_lt(length(unique(ACClong$casn)), length(CAS))
 expect_gt(length(unique(ACClong$endPoint)), length(CAS))
 
})

test_that("Removing flags", {
  testthat::skip_on_cran()
  
  ACClong <- get_ACC(CAS)
  ACClong_noFlags <- remove_flags(ACClong, 
                                  flagsShort = c("Borderline",
                                                 "OnlyHighest",
                                                 "OneAbove",
                                                 "Noisy",
                                                 "HitCall",
                                                 "GainAC50",
                                                 "Biochemical"))
  expect_lt(nrow(ACClong_noFlags), nrow(ACClong))
  expect_true(all(is.na(ACClong_noFlags$flags)))
  
})

test_that("Cleaning up endpoints", {
  testthat::skip_on_cran()
  # This function rejiggers some intended target family and sub-family
  # based on first paper:
  cleaned_ep <- clean_endPoint_info(endPointInfo)
  expect_equal(names(cleaned_ep), names(endPointInfo))
  expect_lt(nrow(cleaned_ep), nrow(endPointInfo))
  
  cleanedNames <- c("Cell Cycle","Nuclear Receptor",
                    "Cell Morphology","DNA Binding",
                    "Background Measurement","Growth Factor",
                    "Cell Adhesion Molecules","Cytokine",
                    "GPCR","Kinase",
                    "Protease","Misc Protein",
                    "Protease Inhibitor","CYP",
                    "Esterase","Phosphatase",
                    "Hydrolase","Oxidoreductase",
                    "Lyase","Methyltransferase",
                    "Ion Channel","Transporter",
                    "Steroid Hormone","Transferase",
                    "Zebrafish","Undefined")

    expect_true(all(cleanedNames %in% cleaned_ep$intended_target_family))
    expect_false(all(cleanedNames %in% endPointInfo$intended_target_family))
})

test_that("Filtering endpoints", {
  testthat::skip_on_cran()
  
  cleaned_ep <- clean_endPoint_info(endPointInfo)
  
  assays <- c("ATG","NVS", "OT")
  groups <- c("Background Measurement","Undefined","Cell Cycle")
  
  filtered_ep <- filter_groups(cleaned_ep, 
                               groupCol = "intended_target_family",
                               assays = assays,
                               remove_groups = groups)
  
  expect_true(all(unique(filtered_ep$assaysFull) %in% assays))
  expect_true(!(any(unique(filtered_ep$groupCol) %in% groups)))
  
})
