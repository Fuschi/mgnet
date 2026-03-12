test_that("mgnet constructs valid minimal and aligned objects", {
  object_min <- mgnet(abun = make_abun())
  expect_s4_class(object_min, "mgnet")
  expect_equal(abun(object_min), make_abun())

  object_full <- mgnet(
    abun = make_abun(),
    meta = make_meta(),
    taxa = make_taxa()
  )

  expect_s4_class(object_full, "mgnet")
  expect_equal(sample_id(object_full), c("s1", "s2"))
  expect_equal(taxa_id(object_full), c("t1", "t2", "t3"))
  expect_equal(rownames(meta(object_full)), c("s1", "s2"))
  expect_equal(meta(object_full)$condition, c("A", "B"))
  expect_equal(rownames(taxa(object_full)), c("t1", "t2", "t3"))
  expect_equal(taxa(object_full)$genus, c("g1", "g2", "g3"))
})

test_that("mgnet default constructor returns a valid empty object", {
  object <- mgnet()

  expect_s4_class(object, "mgnet")
  expect_equal(dim(abun(object)), c(0L, 0L))
  expect_equal(dim(rela(object)), c(0L, 0L))
  expect_equal(dim(norm(object)), c(0L, 0L))
  expect_equal(dim(meta(object)), c(0L, 1L))
  expect_equal(dim(taxa(object)), c(0L, 1L))
  expect_equal(nsample(object), 0L)
  expect_equal(ntaxa(object), 0L)
})
