# Test creation with default parameters
test_that("mgnet object can be created with default parameters", {
  mgnet_obj <- mgnet()
  expect_s4_class(mgnet_obj, "mgnet")
  expect_true(is.matrix(slot(mgnet_obj, "abun")) && length(slot(mgnet_obj, "abun")) == 0)
  expect_true(is.matrix(slot(mgnet_obj, "rela")) && length(slot(mgnet_obj, "rela")) == 0)
  expect_true(is.matrix(slot(mgnet_obj, "norm")) && length(slot(mgnet_obj, "norm")) == 0)
  expect_true(is.data.frame(slot(mgnet_obj, "meta")) && length(slot(mgnet_obj, "meta")) == 0)
  expect_true(is.data.frame(slot(mgnet_obj, "taxa")) && length(slot(mgnet_obj, "taxa")) == 0)
  expect_true(inherits(slot(mgnet_obj, "netw"), "igraph"))
  expect_true(inherits(slot(mgnet_obj, "comm"), "communities")) # Adjust based on actual community class used
})

# Prepare data with matching row and column names where necessary
test_that("mgnet object can be created with matching names and random data", {
  
  abun <- matrix(runif(20, max = 10000), nrow=4, ncol=5)
  rownames(abun) <- paste("Sample", 1:4, sep="_")
  colnames(abun) <- paste("OTU", 1:5, sep="_")
  
  meta <- data.frame(info=letters[1:4])
  rownames(meta) <- paste("Sample", 1:4, sep="_")  # Ensure rownames match those in 'abun'
  
  taxa <- data.frame(details=LETTERS[1:5])
  rownames(taxa) <- paste("OTU", 1:5, sep="_")  # Ensure rownames match column names in 'abun'
  
  netw <- make_empty_graph(n = 5, directed = F)
  V(netw)$name <- colnames(abun)  # Ensure vertex names match OTU IDs in 'abun'
  
  mgnet_obj <- mgnet(abun = abun, meta = meta, taxa = taxa, netw = netw)
  expect_s4_class(mgnet_obj, "mgnet")
  expect_equal(slot(mgnet_obj, "abun"), abun)
  expect_equal(slot(mgnet_obj, "meta"), meta)
  expect_equal(slot(mgnet_obj, "taxa"), taxa)
  expect_identical(slot(mgnet_obj, "netw"), netw)
})

# Test with specific initialization for the HMP2 mgnet object
test_that("HMP2 mgnet object is correctly initialized", {
  # Load the HMP2 object assuming it's already an mgnet object provided by the package
  data("HMP2", package = "mgnet")
  expect_s4_class(HMP2, "mgnet")
  
  # Check basic structure to ensure it contains expected data
  expect_true(is.matrix(slot(HMP2, "abun")))
  expect_true(is.matrix(slot(HMP2, "rela")))
  expect_true(is.matrix(slot(HMP2, "norm")))
  expect_true(is.data.frame(slot(HMP2, "meta")))
  expect_true(is.data.frame(slot(HMP2, "taxa")))
  expect_true(inherits(slot(HMP2, "netw"), "igraph"))
  expect_true(inherits(slot(HMP2, "comm"), "communities"))
  
  # Check specifics that might be known about HMP2
  # This part depends on the known characteristics of the HMP2 dataset. Here are some hypothetical examples:
  expect_equal(nrow(slot(HMP2, "abun")), 1122, 
               info = "HMP2 should have 1122 samples based on the dataset description.")
  expect_equal(ncol(slot(HMP2, "abun")), 1953, 
               info = "HMP2 should have 1953 taxa based on the dataset description.")
  expect_equal(nrow(slot(HMP2, "meta")), 1122, 
               info = "HMP2 should have sample metadata for 1122 samples based on the dataset description.")
  expect_equal(nrow(slot(HMP2, "taxa")), 1953, 
               info = "HMP2 should have taxa metadata for 1953 taxa based on the dataset description.")
})
