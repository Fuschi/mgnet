test_that("mgnets getters return named lists and collapsed tibbles", {
  object <- make_mgnets_example()
  
  abun_list <- abun(object)
  expect_type(abun_list, "list")
  expect_named(abun_list, c("A", "B"))
  expect_equal(dim(abun_list$A), c(2L, 3L))
  expect_equal(dim(abun_list$B), c(2L, 3L))
  
  meta_list <- meta(object)
  expect_type(meta_list, "list")
  expect_named(meta_list, c("A", "B"))
  expect_false("sample_id" %in% names(meta_list$A))
  expect_equal(rownames(meta_list$A), c("s1", "s2"))
  expect_equal(rownames(meta_list$B), c("u1", "u2"))
  
  meta_tbl <- meta(object, .fmt = "tbl", .collapse = TRUE)
  expect_s3_class(meta_tbl, "tbl_df")
  expect_true(all(c("mgnet", "sample_id", "condition", "batch") %in% names(meta_tbl)))
  expect_equal(meta_tbl$mgnet, c("A", "A", "B", "B"))
  
  taxa_list <- taxa(object)
  expect_type(taxa_list, "list")
  expect_named(taxa_list, c("A", "B"))
  expect_false("taxa_id" %in% names(taxa_list$A))
  expect_equal(rownames(taxa_list$A), c("t1", "t2", "t3"))
  
  taxa_tbl <- taxa(object, .fmt = "tbl", .collapse = TRUE)
  expect_s3_class(taxa_tbl, "tbl_df")
  expect_true(all(c("mgnet", "taxa_id", "kingdom", "genus") %in% names(taxa_tbl)))
  expect_equal(sort(unique(taxa_tbl$mgnet)), c("A", "B"))
})
