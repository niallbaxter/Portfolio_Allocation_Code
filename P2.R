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
start_date2 <- "2020-03-12"
end_date2 <- "2022-05-31"    

getSymbols(stocks, from = start_date2, to = end_date2, src = "yahoo")
prices2 <- do.call(merge, lapply(stocks, function(x) Ad(get(x))))
returns2 <- na.omit(Return.calculate(prices2))

stats_df2 <- data.frame(
  Stock = stocks,
  Mean_Return = colMeans(returns2, na.rm = TRUE),
  Std_Dev = apply(returns2, 2, sd, na.rm = TRUE),
  Kurtosis = apply(returns2, 2, kurtosis, na.rm = TRUE),
  Skewness = apply(returns2, 2, skewness, na.rm = TRUE)
)

SS2 <- stats_df2 %>%
  kable("latex", booktabs = TRUE, caption = "Summary Statistics for Stocks Period 2")

writeLines(SS2, "SS2")

portf2 <- portfolio.spec(assets = stocks)
portf2 <- add.constraint(portf2, type = "full_investment")
portf2 <- add.constraint(portf2, type = "long_only")

portf_minvar2 <- portfolio.spec(assets = stocks)
portf_minvar2 <- add.constraint(portf_minvar2, type = "full_investment")
portf_minvar2 <- add.constraint(portf_minvar2, type = "long_only")
portf_minvar2 <- add.objective(portf_minvar2, type = "risk", name = "var")
portf_minvar2$constraints[[1]]$min_sum <- 0.99
portf_minvar2$constraints[[1]]$max_sum <- 1.01

portf_cvar2 <- portfolio.spec(assets = stocks)
portf_cvar2 <- add.constraint(portf_cvar2, type = "full_investment")
portf_cvar2 <- add.constraint(portf_cvar2, type = "long_only")
portf_cvar2 <- add.objective(portf_cvar2, type = "risk", name = "CVaR")
portf_cvar2$constraints[[1]]$min_sum <- 0.99
portf_cvar2$constraints[[1]]$max_sum <- 1.01

risk_parity_optimization2 <- function() {
  S2 <- cov(returns2)
  return(riskParityPortfolio(S2)$w)
}

portf_maxsharpe2 <- portfolio.spec(assets = stocks)
portf_maxsharpe2 <- add.constraint(portf_maxsharpe2, type = "full_investment")
portf_maxsharpe2 <- add.constraint(portf_maxsharpe2, type = "long_only")
portf_maxsharpe2 <- add.objective(portf_maxsharpe2, type = "return", name = "mean")
portf_maxsharpe2 <- add.objective(portf_maxsharpe2, type = "risk", name = "var")
portf_maxsharpe2$constraints[[1]]$min_sum <- 0.99
portf_maxsharpe2$constraints[[1]]$max_sum <- 1.01

equal_weights2 <- rep(1 / length(stocks), length(stocks))

opt_minvar2 <- optimize.portfolio(returns2, portf_minvar2, optimize_method = "DEoptim")
opt_cvar2 <- optimize.portfolio(returns2, portf_cvar2, optimize_method = "DEoptim")
opt_riskparity2 <- risk_parity_optimization2()
opt_maxsharpe2 <- optimize.portfolio(returns2, portf_maxsharpe2, optimize_method = "CVXR", maxSR = TRUE)

normalize_weights2 <- function(weights) {
  return(weights / sum(weights))
}

weights2 <- data.frame(
  Stock = stocks,
  MinVar = normalize_weights2(extractWeights(opt_minvar2)),
  CVaR = normalize_weights2(extractWeights(opt_cvar2)),
  RiskParity = normalize_weights2(opt_riskparity2),
  MaxSharpe = normalize_weights2(extractWeights(opt_maxsharpe2)),
  EqualWeight = equal_weights2
)

writeLines(kable(weights2, format = "latex", booktabs = TRUE) %>% 
             kable_styling(), "weights_table_P2.tex")

custom_palette2 <- c(
  "#E63946", "#F1FAEE", "#A8DADC", "#457B9D", "#1D3557",
  "#F4A261", "#2A9D8F", "#264653", "#D9BF77", "#9C27B0",
  "#F7C59F", "#E76F51", "#9A174D", "#4A4E69"
)

combined_df2 <- data.frame()
for (method in colnames(weights2)[-c(1, 6)]) {  # Exclude the "Stock" column and "EqualWeight" column
  temp_df2 <- data.frame(Stock = weights2$Stock, Weight = weights2[[method]], Method = method)
  combined_df2 <- rbind(combined_df2, temp_df2)
}

faceted_pie_charts2 <- ggplot(combined_df2, aes(x = "", y = Weight, fill = Stock)) +
  geom_bar(stat = "identity", width = 1, color = "black", size = 0.5) +  
  coord_polar("y", start = 0) +
  facet_wrap(~Method, ncol = 2) +  
  theme_void() +
  scale_fill_manual(values = custom_palette2) +  
  theme(
    strip.text = element_text(size = 12, face = "bold"),  
    legend.position = "bottom",  
    legend.title = element_text(size = 10, face = "bold"),  
    legend.text = element_text(size = 8)  
  ) +
  labs(fill = "Stock")

ggsave("faceted_pie_charts_P2.pdf", faceted_pie_charts2, width = 10, height = 8)

calculate_metrics2 <- function(returns, weights) {
  portf_returns <- Return.portfolio(returns, weights = weights)
  mean_ret <- mean(portf_returns)
  std_dev <- sd(portf_returns)
  sharpe <- SharpeRatio(portf_returns, Rf = 0, p = 0.95, FUN = "StdDev")
  sortino <- SortinoRatio(portf_returns)
  return(c(mean_ret, std_dev, sharpe, sortino))
}

metrics2 <- data.frame(
  Metric = c("Mean Return", "Standard Deviation", "Sharpe Ratio", "Sortino Ratio"),
  MinVar = calculate_metrics2(returns2, extractWeights(opt_minvar2)),
  CVaR = calculate_metrics2(returns2, extractWeights(opt_cvar2)),
  RiskParity = calculate_metrics2(returns2, opt_riskparity2),
  MaxSharpe = calculate_metrics2(returns2, extractWeights(opt_maxsharpe2)),
  EqualWeight = calculate_metrics2(returns2, equal_weights2)
)

writeLines(kable(metrics2, format = "latex", booktabs = TRUE) %>% kable_styling(), "metrics_table_P2.tex")

ef_cvar2 <- create.EfficientFrontier(R = returns2, portfolio = portf_cvar2, type = "mean-StdDev", n.portfolios = 25)

ef_data2 <- data.frame(
  Return = ef_cvar2$frontier[, "mean"],
  Risk = ef_cvar2$frontier[, "StdDev"]
)

asset_data2 <- data.frame(
  Return = colMeans(returns2),
  Risk = apply(returns2, 2, sd),
  Label = colnames(returns2)
)

asset_data2$Label <- gsub(".IS.Adjusted", "", asset_data2$Label, fixed = TRUE)

p2 <- ggplot() +
  geom_line(data = ef_data2, aes(x = Risk, y = Return), color = "#0000ff", size = 1.2) +  # Efficient Frontier line
  geom_point(data = asset_data2, aes(x = Risk, y = Return, color = Label), size = 4) +  # Asset points
  geom_text(data = asset_data2, aes(x = Risk, y = Return, label = Label), vjust = -1, hjust = 1, size = 3) +  # Asset labels
  geom_point(aes(x = 0.0159228, y = 0.0013974), color = "#E63946", size = 5, shape = 18) +  # MinVar portfolio point
  annotate("text", x = 0.0159228, y = 0.0013974, label = "MinVar", 
           vjust = -1.5, hjust = 1, color = "#E63946", size = 3.5) +
  geom_point(aes(x = 0.0159228, y = 0.0015200), color = "#E63946", size = 5, shape = 18) +  # CVaR portfolio point
  annotate("text", x = 0.0159228, y = 0.0015200, label = "CVaR", 
           vjust = -1.5, hjust = 1, color = "#E63946", size = 3.5) +
  geom_point(aes(x = 0.0172961, y = 0.0020055), color = "#E63946", size = 5, shape = 18) +  # RiskParity portfolio point
  annotate("text", x = 0.0172961, y = 0.0020055, label = "RiskParity", 
           vjust = -1.5, hjust = 1, color = "#E63946", size = 3.5) +
  geom_point(aes(x = 0.0206818, y = 0.0030721), color = "#E63946", size = 5, shape = 18) +  # MaxSharpe portfolio point
  annotate("text", x = 0.0206818, y = 0.0030721, label = "MaxSharpe", 
           vjust = -1.5, hjust = 1, color = "#E63946", size = 3.5) +
  geom_point(aes(x = 0.0178581, y = 0.0020789), color = "#E63946", size = 5, shape = 18) +  # EqualWeight portfolio point
  annotate("text", x = 0.0178581, y = 0.0020789, label = "EqualWeight", 
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
print(p2)

ggsave("efficient_frontier_P2.pdf", p2, width = 10, height = 8)

cvar_weights2 <- extractWeights(opt_cvar2)

weights_data2 <- data.frame(
  Asset = colnames(returns2),
  Weight = cvar_weights2
)
colnames(weights_data2)[0] <- "name"
weights_data2$Asset <- sub("\\.IS\\.Adjusted$", "", weights_data2$Asset)

colour_mapping2 <- c(
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

w2 <- chart.EF.Weights(ef_cvar2, match.col = "StdDev", colorset = colour_mapping2, legend.loc = "right", cex.legend = 1)
print(w2)

ggsave("cvar_weights_P2.pdf", w2, width = 10, height = 8)
