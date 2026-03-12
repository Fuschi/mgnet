test_that("mgnets abundance-like setters accept named lists", {
  object <- make_mgnets_example()
  
  new_abun <- list(
    B = make_abun_alt()[c("u2", "u1"), , drop = FALSE],
    A = make_abun()[c("s2", "s1"), , drop = FALSE]
  )
  abun(object) <- new_abun
  expect_equal(rownames(abun(object)$A), c("s1", "s2"))
  expect_equal(rownames(abun(object)$B), c("u1", "u2"))
  
  new_rela <- lapply(abun(object), function(x) x / rowSums(x))
  rela(object) <- new_rela
  expect_equal(names(rela(object)), c("A", "B"))
  expect_equal(dim(rela(object)$A), c(2L, 3L))
  
  new_norm <- lapply(abun(object), log1p)
  norm(object) <- new_norm
  expect_equal(names(norm(object)), c("A", "B"))
  expect_equal(dim(norm(object)$B), c(2L, 3L))
})

test_that("mgnets meta and taxa setters accept both list and collapsed table inputs", {
  object <- make_mgnets_example()
  
  meta_list <- meta(object)
  meta_list$A <- meta_list$A[c("s2", "s1"), , drop = FALSE]
  meta_list$B <- meta_list$B[c("u2", "u1"), , drop = FALSE]
  meta(object) <- meta_list
  expect_equal(rownames(meta(object)$A), c("s1", "s2"))
  expect_equal(rownames(meta(object)$B), c("u1", "u2"))
  
  taxa_list <- taxa(object)
  taxa_list$A <- taxa_list$A[c("t3", "t1", "t2"), , drop = FALSE]
  taxa_list$B <- taxa_list$B[c("t2", "t3", "t1"), , drop = FALSE]
  taxa(object) <- taxa_list
  expect_equal(rownames(taxa(object)$A), c("t1", "t2", "t3"))
  expect_equal(rownames(taxa(object)$B), c("t1", "t2", "t3"))
  
  meta_tbl <- meta(object, .fmt = "tbl", .collapse = TRUE)
  meta_tbl <- meta_tbl[c(2, 1, 4, 3), , drop = FALSE]
  meta(object) <- meta_tbl
  expect_equal(rownames(meta(object)$A), c("s1", "s2"))
  expect_equal(rownames(meta(object)$B), c("u1", "u2"))
  
  taxa_tbl <- taxa(object, .fmt = "tbl", .collapse = TRUE)
  taxa_tbl <- taxa_tbl[c(3, 1, 2, 6, 4, 5), , drop = FALSE]
  taxa(object) <- taxa_tbl
  expect_equal(rownames(taxa(object)$A), c("t1", "t2", "t3"))
  expect_equal(rownames(taxa(object)$B), c("t1", "t2", "t3"))
})

test_that("mgnets setters reject malformed inputs", {
  object <- make_mgnets_example()
  
  expect_error(
    {
      abun(object) <- list(A = make_abun())
    },
    "Lengths must match"
  )
  
  bad_meta_tbl <- meta(object, .fmt = "tbl", .collapse = TRUE)
  bad_meta_tbl$mgnet[1] <- "Z"
  expect_error(
    {
      meta(object) <- bad_meta_tbl
    },
    "Some .*mgnet.* values"
  )
  
  bad_taxa_tbl <- taxa(object, .fmt = "tbl", .collapse = TRUE)
  bad_taxa_tbl$taxa_id[1] <- "t999"
  expect_error(
    {
      taxa(object) <- bad_taxa_tbl
    },
    "Provided .*taxa_id.* values do not match those in the object"
  )
})
