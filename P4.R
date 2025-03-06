library(quantmod)
library(PerformanceAnalytics)
library(PortfolioAnalytics)
library(DEoptim)
library(ggplot2)
library(gridExtra)
library(ggpubr)
library(latex2exp)
library(kableExtra)
library(CVXR)
library(riskParityPortfolio)

stocks <- c("AKBNK.IS", "GARAN.IS", "ISCTR.IS", "PETKM.IS", "TUPRS.IS", 
            "EREGL.IS", "TCELL.IS", "BIMAS.IS", "KCHOL.IS", "SAHOL.IS", 
            "THYAO.IS", "PGSUS.IS", "ASELS.IS", "KOZAL.IS")
start_date4 <- "2023-05-29" 
end_date4 <- "2024-03-31"   

getSymbols(stocks, from = start_date4, to = end_date4, src = "yahoo")
prices4 <- do.call(merge, lapply(stocks, function(x) Ad(get(x))))
returns4 <- na.omit(Return.calculate(prices4))

stats_df4 <- data.frame(
  Stock = stocks,
  Mean_Return = colMeans(returns4, na.rm = TRUE),
  Std_Dev = apply(returns4, 2, sd, na.rm = TRUE),
  Kurtosis = apply(returns4, 2, kurtosis, na.rm = TRUE),
  Skewness = apply(returns4, 2, skewness, na.rm = TRUE)
)

SS4 <- stats_df4 %>%
  kable("latex", booktabs = TRUE, caption = "Summary Statistics for Stocks Period 4")

writeLines(SS4, "SS4")

portf4 <- portfolio.spec(assets = stocks)
portf4 <- add.constraint(portf4, type = "full_investment")
portf4 <- add.constraint(portf4, type = "long_only")

portf_minvar4 <- portfolio.spec(assets = stocks)
portf_minvar4 <- add.constraint(portf_minvar4, type = "full_investment")
portf_minvar4 <- add.constraint(portf_minvar4, type = "long_only")
portf_minvar4 <- add.objective(portf_minvar4, type = "risk", name = "var")
portf_minvar4$constraints[[1]]$min_sum <- 0.99
portf_minvar4$constraints[[1]]$max_sum <- 1.01

portf_cvar4 <- portfolio.spec(assets = stocks)
portf_cvar4 <- add.constraint(portf_cvar4, type = "full_investment")
portf_cvar4 <- add.constraint(portf_cvar4, type = "long_only")
portf_cvar4 <- add.objective(portf_cvar4, type = "risk", name = "CVaR")
portf_cvar4$constraints[[1]]$min_sum <- 0.99
portf_cvar4$constraints[[1]]$max_sum <- 1.01

risk_parity_optimization4 <- function() {
  S4 <- cov(returns4)
  return(riskParityPortfolio(S4)$w)
}

portf_maxsharpe4 <- portfolio.spec(assets = stocks)
portf_maxsharpe4 <- add.constraint(portf_maxsharpe4, type = "full_investment")
portf_maxsharpe4 <- add.constraint(portf_maxsharpe4, type = "long_only")
portf_maxsharpe4 <- add.objective(portf_maxsharpe4, type = "return", name = "mean")
portf_maxsharpe4 <- add.objective(portf_maxsharpe4, type = "risk", name = "var")
portf_maxsharpe4$constraints[[1]]$min_sum <- 0.99
portf_maxsharpe4$constraints[[1]]$max_sum <- 1.01

equal_weights4 <- rep(1 / length(stocks), length(stocks))

# Optimize portfolios for Period 4
opt_minvar4 <- optimize.portfolio(returns4, portf_minvar4, optimize_method = "DEoptim")
opt_cvar4 <- optimize.portfolio(returns4, portf_cvar4, optimize_method = "DEoptim")
opt_riskparity4 <- risk_parity_optimization4()
opt_maxsharpe4 <- optimize.portfolio(returns4, portf_maxsharpe4, optimize_method = "CVXR", maxSR = TRUE)

normalize_weights4 <- function(weights) {
  return(weights / sum(weights))
}

weights4 <- data.frame(
  Stock = stocks,
  MinVar = normalize_weights4(extractWeights(opt_minvar4)),
  CVaR = normalize_weights4(extractWeights(opt_cvar4)),
  RiskParity = normalize_weights4(opt_riskparity4),
  MaxSharpe = normalize_weights4(extractWeights(opt_maxsharpe4)),
  EqualWeight = equal_weights4
)

writeLines(kable(weights4, format = "latex", booktabs = TRUE) %>% 
             kable_styling(), "weights_table_P4.tex")

custom_palette4 <- c(
  "#E63946", "#F1FAEE", "#A8DADC", "#457B9D", "#1D3557",
  "#F4A261", "#2A9D8F", "#264653", "#D9BF77", "#9C27B0",
  "#F7C59F", "#E76F51", "#9A174D", "#4A4E69"
)

combined_df4 <- data.frame()
for (method in colnames(weights4)[-c(1, 6)]) {  # Exclude the "Stock" column and "EqualWeight" column
  temp_df4 <- data.frame(Stock = weights4$Stock, Weight = weights4[[method]], Method = method)
  combined_df4 <- rbind(combined_df4, temp_df4)
}

faceted_pie_charts4 <- ggplot(combined_df4, aes(x = "", y = Weight, fill = Stock)) +
  geom_bar(stat = "identity", width = 1, color = "black", size = 0.5) +  
  coord_polar("y", start = 0) +
  facet_wrap(~Method, ncol = 2) +  
  theme_void() +
  scale_fill_manual(values = custom_palette4) +  
  theme(
    strip.text = element_text(size = 12, face = "bold"),  
    legend.position = "bottom",  
    legend.title = element_text(size = 10, face = "bold"),  
    legend.text = element_text(size = 8)  
  ) +
  labs(fill = "Stock")

ggsave("faceted_pie_charts_P4.pdf", faceted_pie_charts4, width = 10, height = 8)

calculate_metrics4 <- function(returns, weights) {
  portf_returns <- Return.portfolio(returns, weights = weights)
  mean_ret <- mean(portf_returns)
  std_dev <- sd(portf_returns)
  sharpe <- SharpeRatio(portf_returns, Rf = 0, p = 0.95, FUN = "StdDev")
  sortino <- SortinoRatio(portf_returns)
  return(c(mean_ret, std_dev, sharpe, sortino))
}

metrics4 <- data.frame(
  Metric = c("Mean Return", "Standard Deviation", "Sharpe Ratio", "Sortino Ratio"),
  MinVar = calculate_metrics4(returns4, extractWeights(opt_minvar4)),
  CVaR = calculate_metrics4(returns4, extractWeights(opt_cvar4)),
  RiskParity = calculate_metrics4(returns4, opt_riskparity4),
  MaxSharpe = calculate_metrics4(returns4, extractWeights(opt_maxsharpe4)),
  EqualWeight = calculate_metrics4(returns4, equal_weights4)
)

writeLines(kable(metrics4, format = "latex", booktabs = TRUE) %>% kable_styling(), "metrics_table_P4.tex")

ef_cvar4 <- create.EfficientFrontier(R = returns4, portfolio = portf_cvar4, type = "mean-StdDev", n.portfolios = 25)

ef_data4 <- data.frame(
  Return = ef_cvar4$frontier[, "mean"],
  Risk = ef_cvar4$frontier[, "StdDev"]
)

asset_data4 <- data.frame(
  Return = colMeans(returns4),
  Risk = apply(returns4, 2, sd),
  Label = colnames(returns4)
)

asset_data4$Label <- gsub(".IS.Adjusted", "", asset_data4$Label, fixed = TRUE)

p4 <- ggplot() +
  geom_line(data = ef_data4, aes(x = Risk, y = Return), color = "#0000ff", size = 1.2) +  # Efficient Frontier line
  geom_point(data = asset_data4, aes(x = Risk, y = Return, color = Label), size = 4) +  # Asset points
  geom_text(data = asset_data4, aes(x = Risk, y = Return, label = Label), vjust = -1, hjust = 1, size = 3) +  # Asset labels
  geom_point(aes(x = 0.018900469, y = 0.004064567), color = "#E63946", size = 5, shape = 18) +  # MinVar portfolio point
  annotate("text", x = 0.018900469, y = 0.004064567, label = "MinVar", 
           vjust = -1.5, hjust = 1, color = "#E63946", size = 3.5) +
  geom_point(aes(x = 0.020222399, y = 0.004737908), color = "#E63946", size = 5, shape = 18) +  # CVaR portfolio point
  annotate("text", x = 0.020222399, y = 0.004737908, label = "CVaR", 
           vjust = -1.5, hjust = 1, color = "#E63946", size = 3.5) +
  geom_point(aes(x = 0.019453655, y = 0.003948528), color = "#E63946", size = 5, shape = 18) +  # RiskParity portfolio point
  annotate("text", x = 0.019453655, y = 0.003948528, label = "RiskParity", 
           vjust = -1.5, hjust = 1, color = "#E63946", size = 3.5) +
  geom_point(aes(x = 0.019838948, y = 0.005052565), color = "#E63946", size = 5, shape = 18) +  # MaxSharpe portfolio point
  annotate("text", x = 0.019838948, y = 0.005052565, label = "MaxSharpe", 
           vjust = -1.5, hjust = 1, color = "#E63946", size = 3.5) +
  geom_point(aes(x = 0.019623333, y = 0.003937645), color = "#E63946", size = 5, shape = 18) +  # EqualWeight portfolio point
  annotate("text", x = 0.019623333, y = 0.003937645, label = "EqualWeight", 
           vjust = -1.5, hjust = 1, color = "#E63946", size = 3.5) +
  scale_color_manual(values = color_mapping) +
  labs(title = "Efficient Frontier",
       x = "Risk (Standard Deviation)",
       y = "Return (Mean)") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    legend.position = "none",
    panel.border = element_rect(color = "black", fill = NA, size = 1)
  )

print(p4)
ggsave("efficient_frontier_P4.pdf", p4, width = 10, height = 8)


cvar_weights4 <- extractWeights(opt_cvar4)

weights_data4 <- data.frame(
  Asset = colnames(returns4),
  Weight = cvar_weights4
)
colnames(weights_data4)[0] <- "name"
weights_data4$Asset <- sub("\\.IS\\.Adjusted$", "", weights_data4$Asset)

colour_mapping4 <- c(
  "AKBNK.IS" = "#E63946",
  "GARAN.IS" = "#1D3557",
  "ISCTR.IS" = "#F4A261",
  "PETKM.IS" = "#D9BF77",
  "TUPRS.IS" = "#264653",
  "EREGL.IS" = "#457B9D",
  "TCELL.IS" = "#E76F51",
  "BIMAS.IS" = "#A8DADC",
  "KCHOL.IS" = "#2A9D8F",
  "SAHOL.IS" = "#F7C59F",
  "THYAO.IS" = "#9A174D",
  "PGSUS.IS" = "#9C27B0",
  "ASELS.IS" = "#F1FAEE",
  "KOZAL.IS" = "#4A4E69"
)

w4 <- chart.EF.Weights(ef_cvar4, match.col = "StdDev", colorset = colour_mapping4, legend.loc = "right", cex.legend = 1)
print(w4)

ggsave("cvar_weights_P4.pdf", w4, width = 10, height = 8)