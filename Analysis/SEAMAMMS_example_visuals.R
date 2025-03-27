library(yarrr)

combined_pairs <- read.csv("Outputs/combined_pairs.csv")

yarrr::pirateplot(
  formula = weight ~ shared_state,
  data = combined_pairs,
  xlab = "Same Reproductive State?",
  ylab = "Association Rate (SRI)",
  gl.col = NA,
  pal = c("yellow", "darkgreen")
)

combined_pairs$combined_state <- paste(do.call(pmin, combined_pairs[9:10]),
  do.call(pmax, combined_pairs[9:10]),
  sep = " - "
)

states <- c("juvenile", "cyc", "preg", "lact")

colrs <- hcl.colors(4)

windows()
pdf(file = "Outputs/shared_state_plots.pdf", width = 12.44)
par(mfrow = c(2, 2))
for (i in 1:4) {
  state <- states[i]

  this_state <- combined_pairs[grep(state, combined_pairs$combined_state), ]

  pirateplot(
    formula = weight ~ combined_state,
    data = this_state,
    main = toupper(state),
    sortx = "mean",
    decreasing = TRUE,
    xlab = NA,
    ylab = "Association Rate (SRI)",
    ylim = c(0, 0.38),
    gl.col = NA,
    pal = colrs[i]
  )
}

dev.off()
