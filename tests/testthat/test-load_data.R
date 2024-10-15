test_that("example data loads correctly", {
  # Load the data from the package
  data(otu_HMP2, package = "mgnet")
  data(meta_HMP2, package = "mgnet")
  data(taxa_HMP2, package = "mgnet")
  data(HMP2, package = "mgnet")
  data(subjects_HMP2, package = "mgnet")
  
  # Check if the data exists in the environment
  expect_true(exists("otu_HMP2"))
  expect_true(exists("meta_HMP2"))
  expect_true(exists("taxa_HMP2"))
  expect_true(exists("HMP2"))
  expect_true(exists("subjects_HMP2"))
  
  # Check the type of each dataset
  expect_true(is.matrix(otu_HMP2), "otu_HMP2 should be a matrix")
  expect_true(is.data.frame(meta_HMP2), "meta_HMP2 should be a data frame")
  expect_true(is.data.frame(taxa_HMP2), "taxa_HMP2 should be a data frame")
  expect_s4_class(HMP2, "mgnet")
  expect_s4_class(subjects_HMP2, "mgnetList")
})