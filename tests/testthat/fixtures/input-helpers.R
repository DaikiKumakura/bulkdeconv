make_input <- function(n_genes = 40, n_samples = 3, species = "human", scale = 1000) {
  set.seed(1)
  genes <- if (species == "human") paste0("GENE", seq_len(n_genes)) else paste0("Gene", seq_len(n_genes))
  m <- matrix(round(rexp(n_genes * n_samples, rate = 1 / scale), 2), n_genes, n_samples)
  df <- data.frame(gene = genes, m, stringsAsFactors = FALSE)
  names(df)[-1L] <- paste0("S", seq_len(n_samples))
  df
}

# A stand-in adapter: returns each sample's total and mean.
dummy_adapter <- function(X, opts, ctx) {
  list(score = cbind(total = colSums(X), mean = colMeans(X)),
       coverage = c(signature = 10, found = 8))
}
