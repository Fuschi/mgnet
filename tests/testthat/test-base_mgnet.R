data("HMP2", package = "mgnet")

test_that("nsample and ntaxa return correct counts", {
  
  # Test nsample
  expect_equal(nsample(HMP2), nrow(HMP2@abun))
  expect_type(nsample(HMP2), "integer") 
  
  # Test ntaxa
  expect_equal(ntaxa(HMP2), ncol(HMP2@abun))  
  expect_type(ntaxa(HMP2), "integer") 
})

test_that("sample_id and taxa_id return correct identifiers", {
  

  expect_equal(sample_id(HMP2), rownames(HMP2@abun))
  expect_equal(sample_id(mgnet()), character(0))
  
  expect_equal(taxa_id(HMP2), colnames(HMP2@abun))
  expect_equal(taxa_id(mgnet()), character(0))
})

test_that("meta_vars return correct information", {
  
  expect_equal(meta_vars(HMP2), colnames(HMP2@meta))
  expect_equal(meta_vars(mgnet()), character(0))
  
})

test_that("taxa_vars return correct information", {
  
  expect_equal(taxa_vars(HMP2), colnames(HMP2@taxa))
  expect_equal(taxa_vars(mgnet()), character(0))
  
})