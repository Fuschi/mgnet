# mgnet: Metagenomic Network Data Structures in R

`mgnet` is an R package designed to **organize, validate, and manipulate metagenomic data**
for network-based analyses.

Rather than providing new statistical models, `mgnet` focuses on **data integrity,
structural consistency, and tidy-style workflows** for next-generation sequencing (NGS)
datasets and their associated networks.

It is particularly suited for analyses where abundance data, metadata, and inferred
networks must remain **strictly aligned** across multiple transformations.

---

## Key features

- **Unified S4 data structure**  
  Store abundance matrices (raw, relative, normalized), sample metadata, taxa annotations,
  networks, and community assignments in a single object.

- **Strict validity checks**  
  Automatic enforcement of:
  - unique and non-empty identifiers
  - consistent row/column alignment across matrices and metadata
  - coherence between abundance data, networks, and communities

- **Tidy-style manipulation without breaking alignment**  
  Subset, group, and mutate data using high-level verbs (`filter_*`, `mutate_*`, `group_mgnet`)
  while preserving identifiers and object validity.

- **Network-aware by design**  
  Seamless integration with `igraph` objects and community detection results, with built-in
  checks to ensure consistency between network nodes and taxa identifiers.

- **Modular and extensible**  
  `mgnet` does not enforce specific normalization or inference methods, making it easy to
  integrate external tools (e.g. correlation inference, layout algorithms, clustering methods).

---

## What `mgnet` does *not* do

- ❌ It does **not** perform normalization or association inference by itself  
- ❌ It does **not** impose a specific statistical model  
- ❌ It does **not** hide transformations behind black-box pipelines  

Instead, `mgnet` provides a **safe container** and **tidy interface** for combining and
manipulating results obtained with your preferred methods.

---

## Installation

### Development version

```r
# install.packages("devtools")
devtools::install_github("Fuschi/mgnet")
```

To build vignettes locally:

```r
devtools::install_github("Fuschi/mgnet", build_vignettes = TRUE)
```

---

## A minimal example

```r
library(mgnet)

# Load example data (Human Microbiome Project, v2)
data(otu_HMP2, package = "mgnet")
data(meta_HMP2, package = "mgnet")
data(taxa_HMP2, package = "mgnet")

# Create an mgnet object
HMP2 <- mgnet(
  abun = otu_HMP2,
  meta = meta_HMP2,
  taxa = taxa_HMP2
)

HMP2
```

This object now contains:
- a **samples × taxa** abundance matrix (`abun`)
- aligned sample metadata (`meta`)
- aligned taxa metadata (`taxa`)

All identifiers are stored in row/column names and are validated automatically.

---

## Typical workflow

A common analysis with `mgnet` involves:

1. Importing abundance data and metadata  
2. Creating an `mgnet` object (with automatic validation)  
3. Adding derived matrices (e.g. relative or normalized abundances)  
4. Subsetting samples and taxa using tidy-style verbs  
5. Attaching a feature-level network and community structure  
6. Exporting tidy representations for visualization or downstream analysis  

See the vignettes for complete, reproducible examples.

---

## Documentation

- **Getting started**  
  ```r
  vignette("getting-started", package = "mgnet")
  ```

- Additional vignettes describe network workflows and advanced usage.

---

## Related packages

`mgnet` is designed to work well with:

- `igraph` — network representation and community detection  
- `tidyverse` — tidy data manipulation and visualization  
- `netkit` — network inference and layout algorithms (optional, via Suggests)

---

## Development status

`mgnet` is under **active development**.  
The API is usable but may change as new features are added and interfaces are refined.

Feedback, issues, and contributions are welcome.

---

## Author

Alessandro Fuschi
