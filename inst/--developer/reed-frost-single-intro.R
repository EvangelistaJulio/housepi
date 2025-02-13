generation_span <- function(n_infec) {
  if (n_infec == 1) {
    return(list(c(0)))
  } else {
    output <- list(rep(0, n_infec))
    for (num_by_root in rev(seq_len(n_infec - 1))) {
      num_by_leaves <- n_infec - num_by_root
      gens <- generation_span(num_by_leaves)
      output <- c(output, lapply(gens, function(g) c(rep(0, num_by_root), g + 1)))
    }
    return(output)
  }
}

generation_count_span <- function(n_infec) {
  gens <- generation_span(n_infec)
  lapply(gens, function(g) as.numeric(table(factor(g, levels = 0:max(g)))))
}


generations_lookup_table <- lapply(1:6, generation_count_span)


# Reed-Frost with single introduction
compute_loglike <- function(N, K, p2p_esc) {
  stopifnot(K > 0, K <= N)
  
  log_esc <- log(p2p_esc)
  
  log_sero_neg <- (N - K) * K * log_esc
  log_sero_pos <- 0
  
  if (K > 1) {
    infec_factor <- 0
    
    for (gen_counts in generations_lookup_table[[K]]) {
      degener <- factorial(K - 1) / factorial(gen_counts[1])
      
      tmp <- gen_counts[1] * log(1 - p2p_esc)
      tmp <- tmp + (K - 1 - gen_counts[1]) * log_esc
      
      for (i in seq(2, length(gen_counts))) {
        degener <- degener / factorial(gen_counts[i])
        tmp <- tmp + gen_counts[i] * log(1 - (p2p_esc ^ gen_counts[i - 1]))
        for (j in seq(1, i - 1)) {
          tmp <- tmp + gen_counts[j] * gen_counts[i] * log_esc
        }
      }
      
      infec_factor <- infec_factor + degener * exp(tmp)
    }
    
    log_sero_pos <- log(infec_factor * choose(N - 1, K - 1))
  }
  
  return(log_sero_pos + log_sero_neg)
}



