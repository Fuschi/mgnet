make_abun <- function() {
  structure(
    matrix(
      c(10, 0, 3, 7, 1, 5),
      nrow = 2,
      byrow = TRUE,
      dimnames = list(c("s1", "s2"), c("t1", "t2", "t3"))
    ),
    class = "matrix"
  )
}

make_meta <- function() {
  data.frame(
    condition = c("A", "B"),
    batch = c("x", "y"),
    row.names = c("s1", "s2"),
    check.names = FALSE
  )
}

make_taxa <- function() {
  data.frame(
    kingdom = c("Bacteria", "Bacteria", "Bacteria"),
    genus = c("g1", "g2", "g3"),
    row.names = c("t1", "t2", "t3"),
    check.names = FALSE
  )
}

make_rela <- function() {
  abun <- make_abun()
  abun / rowSums(abun)
}

make_norm <- function() {
  log1p(make_abun())
}

make_abun_alt <- function() {
  structure(
    matrix(
      c(4, 0, 8, 2, 1, 6),
      nrow = 2,
      byrow = TRUE,
      dimnames = list(c("u1", "u2"), c("t1", "t2", "t3"))
    ),
    class = "matrix"
  )
}

make_meta_alt <- function() {
  data.frame(
    condition = c("C", "D"),
    batch = c("m", "n"),
    row.names = c("u1", "u2"),
    check.names = FALSE
  )
}

make_taxa_alt <- function() {
  data.frame(
    kingdom = c("Bacteria", "Bacteria", "Bacteria"),
    genus = c("g1", "g2", "g3"),
    row.names = c("t1", "t2", "t3"),
    check.names = FALSE
  )
}

make_mgnets_example <- function() {
  mgnets(
    A = mgnet(
      abun = make_abun(),
      meta = make_meta(),
      taxa = make_taxa()
    ),
    B = mgnet(
      abun = make_abun_alt(),
      meta = make_meta_alt(),
      taxa = make_taxa_alt()
    )
  )
}
