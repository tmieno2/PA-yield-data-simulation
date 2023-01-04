
library(magic) # rlatin()

assign_input_rate <- function(N_levels, block_num, design) {
  if (design == "Latin Square Fixed") {

      # () Latin Square Fixed
      M.latin <- matrix(c(
          1, 3, 5, 6, 4, 2,
          5, 6, 4, 2, 1, 3,
          2, 1, 3, 5, 6, 4,
          6, 4, 2, 1, 3, 5,
          3, 5, 6, 4, 2, 1,
          4, 2, 1, 3, 5, 6
      ),
      nrow = 6, ncol = 6
      ) %>% t()
      N_design <- data.table(
          block_id = rep(1:block_num, each = 36),
          plot_in_block_id = rep(1:36, block_num),
          N_tgt = N_levels[rep(c(t(M.latin)), times = block_num)]
      )
  } else if (design == "Latin Square Fixed 5") {
      
    # (1) Latin Square Fixed 5
    M.latin <- matrix(c(
      1, 3, 5, 4, 2,
      5, 4, 2, 1, 3,
      2, 1, 3, 5, 4,
      4, 2, 1, 3, 5,
      3, 5, 4, 2, 1
    ),
    nrow = 5, ncol = 5
    ) %>% t()
    N_design <- data.table(
      block_id = rep(1:block_num, each = 25),
      plot_in_block_id = rep(1:25, block_num),
      N_tgt = N_levels[rep(c(t(M.latin)), times = block_num)]
    )
  } else if (design == "Latin Square Random") {

    # (1) Latin Square Random
    N_design <- data.table(
      block_id = rep(1:block_num, each = 36),
      plot_in_block_id = rep(1:36, block_num),
      # === randomly select block_num latin squares ===#
      N_tgt = c(rlatin(n = block_num, size = 6)) %>% N_levels[.]
    )
  } else if (design == "Latin Square Cascade") {

    # (1) Latin Square Cascade
    M.latin <- matrix(c(
      1, 2, 3, 4, 5, 6,
      2, 3, 4, 5, 6, 1,
      3, 4, 5, 6, 1, 2,
      4, 5, 6, 1, 2, 3,
      5, 6, 1, 2, 3, 4,
      6, 1, 2, 3, 4, 5
    ),
    nrow = 6, ncol = 6
    ) %>% t()
    N_design <- data.table(
      block_id = rep(1:block_num, each = 36),
      plot_in_block_id = rep(1:36, block_num),
      N_tgt = N_levels[rep(c(t(M.latin)), times = block_num)]
    )
  } else if (design == "Alternate Block") {

    # (2) Alternate Block
    M.nojump <- matrix(c(
      6, 5, 4,
      3, 2, 1,
      5, 4, 6,
      2, 1, 3,
      4, 6, 5,
      1, 3, 2
    ),
    nrow = 3, ncol = 6
    ) %>% t()
    N_design <- data.table(
      block_id = rep(1:block_num, each = 18),
      plot_in_block_id = rep(1:18, block_num),
      N_tgt = N_levels[rep(c(t(M.nojump)), times = block_num)]
    )
  } else if (design == "Checkerboard") {

    # (3) Checkerboard
    M.nojump <- rbind(
      c(1, 3),
      c(6, 5),
      c(4, 2)
    )
    N_design <- data.table(
      block_id = rep(1:block_num, each = 6),
      plot_in_block_id = rep(1:6, block_num),
      N_tgt = N_levels[rep(c(t(M.nojump)), times = block_num)]
    )
  } else if (design == "Randomized Block") {

    # (4) Randomized Block
    N_design <- data.table(
      block_id = rep(1:block_num, each = 6),
      plot_in_block_id = rep(1:6, block_num),
      N_tgt = replicate(block_num, sample(N_levels, replace = FALSE)) %>% as.vector()
    )
  } else if (design == "Completely Random") {

    # (5) Completely Random
    N_design <- data.table(
      block_id = rep(1:block_num, each = 6),
      plot_in_block_id = rep(1:6, block_num),
      N_tgt = sample(N_levels, size = 6 * block_num, replace = TRUE)
    )
  } else if (design == "Fixed Strip Grad") {

    # (6) Fixed Strip Grad
    N_design <- data.table(
      block_id = rep(1:block_num, each = 12),
      plot_in_block_id = rep(1:12, block_num),
      N_tgt = N_levels[rep(c(1:6, 6:1), times = block_num)]
    )
  } else if (design == "Fixed Strip Fluc 1") {

    # (7) Fixed Strip Fluc 1
    N_design <- data.table(
      block_id = rep(1:block_num, each = 6),
      plot_in_block_id = rep(1:6, block_num),
      N_tgt = N_levels[rep(c(1, 6, 3, 5, 2, 4), times = block_num)]
    )
  } else if (design == "Random Strip") {

    # (8) Random Strip
    N_design <- data.table(
      block_id = rep(1:block_num, each = 6),
      plot_in_block_id = rep(1:6, block_num),
      N_tgt = replicate(block_num, sample(N_levels, replace = FALSE)) %>% as.vector()
    )
  } else if (design == "Cascade Plot") {

    # (9) Cascade Plot
    M.cascade <- rbind(
      c(1:6, 6:1),
      c(2:6, 6:1, 1:1),
      c(3:6, 6:1, 1:2),
      c(4:6, 6:1, 1:3),
      c(5:6, 6:1, 1:4),
      c(6:6, 6:1, 1:5),
      c(6:1, 1:6),
      c(5:1, 1:6, 6:6),
      c(4:1, 1:6, 6:5),
      c(3:1, 1:6, 6:4),
      c(2:1, 1:6, 6:3),
      c(1:1, 1:6, 6:2)
    )
    N_design <- data.table(
      block_id = rep(1:block_num, each = 144),
      plot_in_block_id = rep(1:144, block_num),
      N_tgt = N_levels[rep(c(t(M.cascade)), times = block_num)]
    )
  } else if (design == "Wave") {

    # (10) Wave
    M.wave <- rbind(
      c(4, 4, 4, 4, 4, 3, 3, 3, 3, 3),
      c(4, 5, 5, 5, 4, 3, 2, 2, 2, 3),
      c(4, 5, 6, 5, 4, 3, 2, 1, 2, 3),
      c(4, 5, 5, 5, 4, 3, 2, 2, 2, 3),
      c(4, 4, 4, 4, 4, 3, 3, 3, 3, 3),
      c(3, 3, 3, 3, 3, 4, 4, 4, 4, 4),
      c(3, 2, 2, 2, 3, 4, 5, 5, 5, 4),
      c(3, 2, 1, 2, 3, 4, 5, 6, 5, 4),
      c(3, 2, 2, 2, 3, 4, 5, 5, 5, 4),
      c(3, 3, 3, 3, 3, 4, 4, 4, 4, 4)
    )
    N_design <- data.table(
      block_id = rep(1:block_num, each = 100),
      plot_in_block_id = rep(1:100, block_num),
      N_tgt = N_levels[rep(c(t(M.wave)), times = block_num)]
    )
  } else if (design == "Fixed Strip Fluc 2") {

    # (11) Fixed Strip Fluc 2
    N_design <- data.table(
      block_id = rep(1:block_num, each = 6),
      plot_in_block_id = rep(1:6, block_num),
      N_tgt = N_levels[rep(c(1, 5, 2, 6, 3, 4), times = block_num)]
    )
  } else {
    N_design <- NA
  }

  return(N_design)
}