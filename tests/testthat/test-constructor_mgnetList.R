# Test empty mgnetList construction
test_that("mgnetList can be created empty", {
  ml_empty <- mgnetList()
  expect_s4_class(ml_empty, "mgnetList")
  expect_equal(length(ml_empty@mgnets), 0)
})

# Test construction with mgnet objects
test_that("mgnetList can be created with mgnet objects", {
  # Assuming 'mgnet_obj1' and 'mgnet_obj2' are valid mgnet objects
  mgnet_obj1 <- mgnet()  # You would typically populate these with actual data
  mgnet_obj2 <- mgnet()
  
  ml_two <- mgnetList("mg1" = mgnet_obj1, "mg2" = mgnet_obj2)
  expect_s4_class(ml_two, "mgnetList")
  expect_equal(length(ml_two@mgnets), 2)
  expect_true(all(sapply(ml_two@mgnets, inherits, "mgnet")))
  expect_setequal(names(ml_two@mgnets), c("mg1", "mg2"))
})

# Assuming subjects_HMP2 is loaded as part of your package's data
test_that("subjects_HMP2 is correctly structured", {
  # Check if it exists and has the correct S4 class
  utils::data("subjects_HMP2", package = "mgnet")
  expect_s4_class(subjects_HMP2, "mgnetList")
  
  # Check for correct number of mgnet objects
  expect_equal(length(subjects_HMP2@mgnets), 2)
  expect_true(all(sapply(subjects_HMP2@mgnets, inherits, "mgnet")))
  
  # Check that the names are as expected
  expect_setequal(names(subjects_HMP2@mgnets), c("69-001", "69-053"))
  
  # Additional checks for content of each mgnet object, if needed
  # This assumes you have details like the number of samples or taxa you expect for each subject
  mgnet_69001 <- subjects_HMP2@mgnets[["69-001"]]
  mgnet_69053 <- subjects_HMP2@mgnets[["69-053"]]
  
  expect_equal(nrow(slot(mgnet_69001, "abun")), 185)  # Number of samples in "69-001"
  expect_equal(ncol(slot(mgnet_69001, "abun")), 977)  # Number of taxa in "69-001"
  
  expect_equal(nrow(slot(mgnet_69053, "abun")), 40)   # Number of samples in "69-053"
  expect_equal(ncol(slot(mgnet_69053, "abun")), 793)  # Number of taxa in "69-053"

})

# Check mgnetList can be created as a list of named mgnet
test_that("mgnetList could construct through a sequence of named mgnet", {
  # Check if it exists and has the correct S4 class
  utils::data("HMP2", package = "mgnet")
  
  ml <- mgnetList("A" = HMP2, "B" = HMP2)
  
  # Check for correct number of mgnet objects
  expect_equal(length(ml@mgnets), 2)
  expect_true(all(sapply(ml@mgnets, inherits, "mgnet")))
  
  # Check that the names are as expected
  expect_setequal(names(ml@mgnets), c("A", "B"))
})

# Check mgnetList can be created as a named list of mgnet direclty
test_that("mgnetList could construct through a named list of mgnets", {
  # Check if it exists and has the correct S4 class
  utils::data("HMP2", package = "mgnet")
  
  ml <- list("A" = HMP2, "B" = HMP2)
  ml <- mgnetList(ml)
  
  # Check for correct number of mgnet objects
  expect_equal(length(ml@mgnets), 2)
  expect_true(all(sapply(ml@mgnets, inherits, "mgnet")))
  
  # Check that the names are as expected
  expect_setequal(names(ml@mgnets), c("A", "B"))
})
