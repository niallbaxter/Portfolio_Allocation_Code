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
start_date5 <- "2024-04-01"  
end_date5 <- "2024-09-01"    

getSymbols(stocks, from = start_date5, to = end_date5, src = "yahoo")
prices5 <- do.call(merge, lapply(stocks, function(x) Ad(get(x))))
returns5 <- na.omit(Return.calculate(prices5))

stats_df5 <- data.frame(
  Stock = stocks,
  Mean_Return = colMeans(returns5, na.rm = TRUE),
  Std_Dev = apply(returns5, 2, sd, na.rm = TRUE),
  Kurtosis = apply(returns5, 2, kurtosis, na.rm = TRUE),
  Skewness = apply(returns5, 2, skewness, na.rm = TRUE)
)

SS5 <- stats_df5 %>%
  kable("latex", booktabs = TRUE, caption = "Summary Statistics for Stocks Period 5")

writeLines(SS5, "SS5")

portf5 <- portfolio.spec(assets = stocks)
portf5 <- add.constraint(portf5, type = "full_investment")
portf5 <- add.constraint(portf5, type = "long_only")

portf_minvar5 <- portfolio.spec(assets = stocks)
portf_minvar5 <- add.constraint(portf_minvar5, type = "full_investment")
portf_minvar5 <- add.constraint(portf_minvar5, type = "long_only")
portf_minvar5 <- add.objective(portf_minvar5, type = "risk", name = "var")
portf_minvar5$constraints[[1]]$min_sum <- 0.99
portf_minvar5$constraints[[1]]$max_sum <- 1.01

portf_cvar5 <- portfolio.spec(assets = stocks)
portf_cvar5 <- add.constraint(portf_cvar5, type = "full_investment")
portf_cvar5 <- add.constraint(portf_cvar5, type = "long_only")
portf_cvar5 <- add.objective(portf_cvar5, type = "risk", name = "CVaR")
portf_cvar5$constraints[[1]]$min_sum <- 0.99
portf_cvar5$constraints[[1]]$max_sum <- 1.01

risk_parity_optimization5 <- function() {
  S5 <- cov(returns5)
  return(riskParityPortfolio(S5)$w)
}

portf_maxsharpe5 <- portfolio.spec(assets = stocks)
portf_maxsharpe5 <- add.constraint(portf_maxsharpe5, type = "full_investment")
portf_maxsharpe5 <- add.constraint(portf_maxsharpe5, type = "long_only")
portf_maxsharpe5 <- add.objective(portf_maxsharpe5, type = "return", name = "mean")
portf_maxsharpe5 <- add.objective(portf_maxsharpe5, type = "risk", name = "var")
portf_maxsharpe5$constraints[[1]]$min_sum <- 0.99
portf_maxsharpe5$constraints[[1]]$max_sum <- 1.01

equal_weights5 <- rep(1 / length(stocks), length(stocks))

opt_minvar5 <- optimize.portfolio(returns5, portf_minvar5, optimize_method = "DEoptim")
opt_cvar5 <- optimize.portfolio(returns5, portf_cvar5, optimize_method = "DEoptim")
opt_riskparity5 <- risk_parity_optimization5()
opt_maxsharpe5 <- optimize.portfolio(returns5, portf_maxsharpe5, optimize_method = "CVXR", maxSR = TRUE)

normalize_weights5 <- function(weights) {
  return(weights / sum(weights))
}

weights5 <- data.frame(
  Stock = stocks,
  MinVar = normalize_weights5(extractWeights(opt_minvar5)),
  CVaR = normalize_weights5(extractWeights(opt_cvar5)),
  RiskParity = normalize_weights5(opt_riskparity5),
  MaxSharpe = normalize_weights5(extractWeights(opt_maxsharpe5)),
  EqualWeight = equal_weights5
)

writeLines(kable(weights5, format = "latex", booktabs = TRUE) %>% 
             kable_styling(), "weights_table_P5.tex")

custom_palette5 <- c(
  "#E63946", "#F1FAEE", "#A8DADC", "#457B9D", "#1D3557",
  "#F4A261", "#2A9D8F", "#264653", "#D9BF77", "#9C27B0",
  "#F7C59F", "#E76F51", "#9A174D", "#4A4E69"
)

combined_df5 <- data.frame()
for (method in colnames(weights5)[-c(1, 6)]) {  # Exclude the "Stock" column and "EqualWeight" column
  temp_df5 <- data.frame(Stock = weights5$Stock, Weight = weights5[[method]], Method = method)
  combined_df5 <- rbind(combined_df5, temp_df5)
}

faceted_pie_charts5 <- ggplot(combined_df5, aes(x = "", y = Weight, fill = Stock)) +
  geom_bar(stat = "identity", width = 1, color = "black", size = 0.5) +  
  coord_polar("y", start = 0) +
  facet_wrap(~Method, ncol = 2) +  
  theme_void() +
  scale_fill_manual(values = custom_palette5) +  
  theme(
    strip.text = element_text(size = 12, face = "bold"),  
    legend.position = "bottom",  
    legend.title = element_text(size = 10, face = "bold"),  
    legend.text = element_text(size = 8)  
  ) +
  labs(fill = "Stock")

ggsave("faceted_pie_charts_P5.pdf", faceted_pie_charts5, width = 10, height = 8)

calculate_metrics5 <- function(returns, weights) {
  portf_returns <- Return.portfolio(returns, weights = weights)
  mean_ret <- mean(portf_returns)
  std_dev <- sd(portf_returns)
  sharpe <- SharpeRatio(portf_returns, Rf = 0, p = 0.95, FUN = "StdDev")
  sortino <- SortinoRatio(portf_returns)
  return(c(mean_ret, std_dev, sharpe, sortino))
}

metrics5 <- data.frame(
  Metric = c("Mean Return", "Standard Deviation", "Sharpe Ratio", "Sortino Ratio"),
  MinVar = calculate_metrics5(returns5, extractWeights(opt_minvar5)),
  CVaR = calculate_metrics5(returns5, extractWeights(opt_cvar5)),
  RiskParity = calculate_metrics5(returns5, opt_riskparity5),
  MaxSharpe = calculate_metrics5(returns5, extractWeights(opt_maxsharpe5)),
  EqualWeight = calculate_metrics5(returns5, equal_weights5)
)

writeLines(kable(metrics5, format = "latex", booktabs = TRUE) %>% kable_styling(), "metrics_table_P5.tex")

ef_cvar5 <- create.EfficientFrontier(R = returns5, portfolio = portf_cvar5, type = "mean-StdDev", n.portfolios = 25)

ef_data5 <- data.frame(
  Return = ef_cvar5$frontier[, "mean"],
  Risk = ef_cvar5$frontier[, "StdDev"]
)

asset_data5 <- data.frame(
  Return = colMeans(returns5),
  Risk = apply(returns5, 2, sd),
  Label = colnames(returns5)
)

asset_data5$Label <- gsub(".IS.Adjusted", "", asset_data5$Label, fixed = TRUE)


p5 <- ggplot() +
  geom_line(data = ef_data5, aes(x = Risk, y = Return), color = "#0000ff", size = 1.2) +  # Efficient Frontier line
  geom_point(data = asset_data5, aes(x = Risk, y = Return, color = Label), size = 4) +  # Asset points
  geom_text(data = asset_data5, aes(x = Risk, y = Return, label = Label), vjust = -1, hjust = 1, size = 3) +  # Asset labels
  geom_point(aes(x = 0.014892719, y = 0.001949103), color = "#E63946", size = 5, shape = 18) +  # MinVar portfolio point
  annotate("text", x = 0.014892719, y = 0.001949103, label = "MinVar", 
           vjust = -1.5, hjust = 1, color = "#E63946", size = 3.5) +
  geom_point(aes(x = 0.01478062, y = 0.00209216), color = "#E63946", size = 5, shape = 18) +  # CVaR portfolio point
  annotate("text", x = 0.01478062, y = 0.00209216, label = "CVaR", 
           vjust = -1.5, hjust = 1, color = "#E63946", size = 3.5) +
  geom_point(aes(x = 0.016188467, y = 0.002253234), color = "#E63946", size = 5, shape = 18) +  # RiskParity portfolio point
  annotate("text", x = 0.016188467, y = 0.002253234, label = "RiskParity", 
           vjust = -1.5, hjust = 1, color = "#E63946", size = 3.5) +
  geom_point(aes(x = 0.016863595, y = 0.004419644), color = "#E63946", size = 5, shape = 18) +  # MaxSharpe portfolio point
  annotate("text", x = 0.016863595, y = 0.004419644, label = "MaxSharpe", 
           vjust = -1.5, hjust = 1, color = "#E63946", size = 3.5) +
  geom_point(aes(x = 0.016691154, y = 0.002206185), color = "#E63946", size = 5, shape = 18) +  # EqualWeight portfolio point
  annotate("text", x = 0.016691154, y = 0.002206185, label = "EqualWeight", 
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

print(p5)
ggsave("efficient_frontier_P5.pdf", p5, width = 10, height = 8)

cvar_weights5 <- extractWeights(opt_cvar5)

weights_data5 <- data.frame(
  Asset = colnames(returns5),
  Weight = cvar_weights5
)
colnames(weights_data5)[0] <- "name"
weights_data5$Asset <- sub("\\.IS\\.Adjusted$", "", weights_data5$Asset)

colour_mapping5 <- c(
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

w5 <- chart.EF.Weights(ef_cvar5, match.col = "StdDev", colorset = colour_mapping5, legend.loc = "right", cex.legend = 1)
print(w5)

ggsave("cvar_weights_P5.pdf", w5, width = 10, height = 8)