test_that("mgnet validates abundance input structure and values", {
  non_numeric <- matrix(
    c("a", "b", "c", "d"),
    nrow = 2,
    dimnames = list(c("s1", "s2"), c("t1", "t2"))
  )
  expect_error(
    mgnet(abun = non_numeric),
    "abun.*numeric matrix"
  )

  no_rownames <- unname(make_abun())
  expect_error(
    mgnet(abun = no_rownames),
    "abun.*row names"
  )

  no_colnames <- make_abun()
  colnames(no_colnames) <- NULL
  expect_error(
    mgnet(abun = no_colnames),
    "abun.*column names"
  )

  duplicated_samples <- make_abun()
  rownames(duplicated_samples) <- c("s1", "s1")
  expect_error(
    mgnet(abun = duplicated_samples),
    "Duplicate row names.*abun"
  )

  duplicated_taxa <- make_abun()
  colnames(duplicated_taxa) <- c("t1", "t1", "t3")
  expect_error(
    mgnet(abun = duplicated_taxa),
    "Duplicate column names.*abun"
  )

  negative_abun <- make_abun()
  negative_abun[1, 1] <- -1
  expect_error(
    mgnet(abun = negative_abun),
    "Negative values.*abun"
  )
})

test_that("mgnet enforces reciprocal alignment with meta", {
  expect_s4_class(
    mgnet(abun = make_abun(), meta = make_meta()),
    "mgnet"
  )

  meta_missing <- make_meta()["s1", , drop = FALSE]
  expect_error(
    mgnet(abun = make_abun(), meta = meta_missing),
    "Row count mismatch between.*abun.*meta"
  )

  meta_reordered <- make_meta()[c("s2", "s1"), , drop = FALSE]
  expect_error(
    mgnet(abun = make_abun(), meta = meta_reordered),
    "Row names of.*abun.*meta must be identical"
  )

  meta_duplicated <- structure(
    list(condition = c("A", "B"), batch = c("x", "y")),
    row.names = c("s1", "s1"),
    class = "data.frame"
  )
  expect_error(
    mgnet(abun = make_abun(), meta = meta_duplicated),
    "Duplicate row names.*meta"
  )
})

test_that("mgnet enforces reciprocal alignment with taxa", {
  expect_s4_class(
    mgnet(abun = make_abun(), taxa = make_taxa()),
    "mgnet"
  )

  taxa_missing <- make_taxa()[c("t1", "t2"), , drop = FALSE]
  expect_error(
    mgnet(abun = make_abun(), taxa = taxa_missing),
    "Columns count in.*abun.*rows in.*taxa mismatch"
  )

  taxa_reordered <- make_taxa()[c("t2", "t1", "t3"), , drop = FALSE]
  expect_error(
    mgnet(abun = make_abun(), taxa = taxa_reordered),
    "Column names in.*abun.*row names in.*taxa must match"
  )

  taxa_duplicated <- structure(
    list(kingdom = c("Bacteria", "Bacteria", "Bacteria"), genus = c("g1", "g2", "g3")),
    row.names = c("t1", "t1", "t3"),
    class = "data.frame"
  )
  expect_error(
    mgnet(abun = make_abun(), taxa = taxa_duplicated),
    "Duplicate row names.*taxa"
  )
})

test_that("mgnet validates rela and norm consistency against abun", {
  expect_s4_class(
    mgnet(
      abun = make_abun(),
      rela = make_rela(),
      norm = make_norm()
    ),
    "mgnet"
  )

  rela_bad_rows <- make_rela()
  rownames(rela_bad_rows) <- c("s1", "s3")
  expect_error(
    mgnet(abun = make_abun(), rela = rela_bad_rows),
    "Row names of.*abun.*rela must be identical"
  )

  rela_bad_cols <- make_rela()
  colnames(rela_bad_cols) <- c("t1", "tX", "t3")
  expect_error(
    mgnet(abun = make_abun(), rela = rela_bad_cols),
    "Column names of.*abun.*rela must be identical"
  )

  rela_bad_zero <- make_rela()
  rela_bad_zero[1, 2] <- 0.25
  expect_error(
    mgnet(abun = make_abun(), rela = rela_bad_zero),
    "Inconsistent zero patterns between.*abun.*rela"
  )

  norm_bad_rows <- make_norm()
  rownames(norm_bad_rows) <- c("s1", "s3")
  expect_error(
    mgnet(abun = make_abun(), norm = norm_bad_rows),
    "Row names of.*abun.*norm must be identical"
  )

  norm_bad_dim <- make_norm()[, c("t1", "t2"), drop = FALSE]
  expect_error(
    mgnet(abun = make_abun(), norm = norm_bad_dim),
    "Column count mismatch between.*abun.*norm"
  )
})

test_that("mgnet rejects reserved column names in aligned tables", {
  meta_reserved <- make_meta()
  meta_reserved$sample_id <- c("s1", "s2")
  expect_error(
    mgnet(abun = make_abun(), meta = meta_reserved),
    "reserved keyword"
  )

  taxa_reserved <- make_taxa()
  taxa_reserved$comm_id <- c("c1", "c2", "c3")
  expect_error(
    mgnet(abun = make_abun(), taxa = taxa_reserved),
    "reserved keyword"
  )

  abun_reserved <- make_abun()
  colnames(abun_reserved)[1] <- "abun"
  expect_error(
    mgnet(abun = abun_reserved),
    "reserved keyword"
  )
})
