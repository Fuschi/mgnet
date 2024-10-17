# Load example mgnetList from package data
data("subjects_HMP2", package = "mgnet")

# 2. Accession
test_that("subjects_HMP2 allows correct element access", {
  # Access by name
  expect_s4_class(subjects_HMP2[["69-001"]], "mgnet")
  # Access by index
  expect_s4_class(subjects_HMP2@mgnets[[1]], "mgnet")
  # Access by [
  expect_s4_class(subjects_HMP2[1], "mgnetList")
})

# 3. Subsetting
test_that("subjects_HMP2 can be subset correctly", {
  subset_list <- subjects_HMP2["69-001"]
  expect_equal(length(subset_list@mgnets), 1)
  expect_true("69-001" %in% names(subset_list@mgnets))
})

# 4. Replacing
test_that("subjects_HMP2 allows element replacement", {
  new_mgnet <- mgnet()  # Create a new mgnet object to replace an existing one
  old_mgnet <- subjects_HMP2[["69-001"]]  # Keep old for comparison if needed
  
  subjects_HMP2[["69-001"]] <- new_mgnet
  expect_identical(subjects_HMP2[["69-001"]], new_mgnet)
  expect_s4_class(subjects_HMP2[["69-001"]], "mgnet")
})

test_that("nsample and ntaxa return correct counts", {
  
  # Test nsample
  expect_equal(nsample(subjects_HMP2), c('69-001'=185,'69-053'=40))
  expect_type(nsample(subjects_HMP2), "integer") 
  
  # Test ntaxa
  expect_equal(ntaxa(subjects_HMP2), c('69-001'=977,'69-053'=793))  
  expect_type(ntaxa(subjects_HMP2), "integer") 
})

