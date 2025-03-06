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
start_date3 <- "2022-06-01"  
end_date3 <- "2023-05-28"    

getSymbols(stocks, from = start_date3, to = end_date3, src = "yahoo")
prices3 <- do.call(merge, lapply(stocks, function(x) Ad(get(x))))
returns3 <- na.omit(Return.calculate(prices3))

stats_df3 <- data.frame(
  Stock = stocks,
  Mean_Return = colMeans(returns3, na.rm = TRUE),
  Std_Dev = apply(returns3, 2, sd, na.rm = TRUE),
  Kurtosis = apply(returns3, 2, kurtosis, na.rm = TRUE),
  Skewness = apply(returns3, 2, skewness, na.rm = TRUE)
)

SS3 <- stats_df3 %>%
  kable("latex", booktabs = TRUE, caption = "Summary Statistics for Stocks Period 3")

writeLines(SS3, "SS3")

portf3 <- portfolio.spec(assets = stocks)
portf3 <- add.constraint(portf3, type = "full_investment")
portf3 <- add.constraint(portf3, type = "long_only")

portf_minvar3 <- portfolio.spec(assets = stocks)
portf_minvar3 <- add.constraint(portf_minvar3, type = "full_investment")
portf_minvar3 <- add.constraint(portf_minvar3, type = "long_only")
portf_minvar3 <- add.objective(portf_minvar3, type = "risk", name = "var")
portf_minvar3$constraints[[1]]$min_sum <- 0.99
portf_minvar3$constraints[[1]]$max_sum <- 1.01

portf_cvar3 <- portfolio.spec(assets = stocks)
portf_cvar3 <- add.constraint(portf_cvar3, type = "full_investment")
portf_cvar3 <- add.constraint(portf_cvar3, type = "long_only")
portf_cvar3 <- add.objective(portf_cvar3, type = "risk", name = "CVaR")
portf_cvar3$constraints[[1]]$min_sum <- 0.99
portf_cvar3$constraints[[1]]$max_sum <- 1.01

risk_parity_optimization3 <- function() {
  S3 <- cov(returns3)
  return(riskParityPortfolio(S3)$w)
}

portf_maxsharpe3 <- portfolio.spec(assets = stocks)
portf_maxsharpe3 <- add.constraint(portf_maxsharpe3, type = "full_investment")
portf_maxsharpe3 <- add.constraint(portf_maxsharpe3, type = "long_only")
portf_maxsharpe3 <- add.objective(portf_maxsharpe3, type = "return", name = "mean")
portf_maxsharpe3 <- add.objective(portf_maxsharpe3, type = "risk", name = "var")
portf_maxsharpe3$constraints[[1]]$min_sum <- 0.99
portf_maxsharpe3$constraints[[1]]$max_sum <- 1.01

equal_weights3 <- rep(1 / length(stocks), length(stocks))

opt_minvar3 <- optimize.portfolio(returns3, portf_minvar3, optimize_method = "DEoptim")
opt_cvar3 <- optimize.portfolio(returns3, portf_cvar3, optimize_method = "DEoptim")
opt_riskparity3 <- risk_parity_optimization3()
opt_maxsharpe3 <- optimize.portfolio(returns3, portf_maxsharpe3, optimize_method = "CVXR", maxSR = TRUE)

normalize_weights3 <- function(weights) {
  return(weights / sum(weights))
}

weights3 <- data.frame(
  Stock = stocks,
  MinVar = normalize_weights3(extractWeights(opt_minvar3)),
  CVaR = normalize_weights3(extractWeights(opt_cvar3)),
  RiskParity = normalize_weights3(opt_riskparity3),
  MaxSharpe = normalize_weights3(extractWeights(opt_maxsharpe3)),
  EqualWeight = equal_weights3
)

writeLines(kable(weights3, format = "latex", booktabs = TRUE) %>% 
             kable_styling(), "weights_table_P3.tex")

custom_palette3 <- c(
  "#E63946", "#F1FAEE", "#A8DADC", "#457B9D", "#1D3557",
  "#F4A261", "#2A9D8F", "#264653", "#D9BF77", "#9C27B0",
  "#F7C59F", "#E76F51", "#9A174D", "#4A4E69"
)

combined_df3 <- data.frame()
for (method in colnames(weights3)[-c(1, 6)]) {  # Exclude the "Stock" column and "EqualWeight" column
  temp_df3 <- data.frame(Stock = weights3$Stock, Weight = weights3[[method]], Method = method)
  combined_df3 <- rbind(combined_df3, temp_df3)
}

faceted_pie_charts3 <- ggplot(combined_df3, aes(x = "", y = Weight, fill = Stock)) +
  geom_bar(stat = "identity", width = 1, color = "black", size = 0.5) +  
  coord_polar("y", start = 0) +
  facet_wrap(~Method, ncol = 2) +  
  theme_void() +
  scale_fill_manual(values = custom_palette3) +  
  theme(
    strip.text = element_text(size = 12, face = "bold"),  
    legend.position = "bottom",  
    legend.title = element_text(size = 10, face = "bold"),  
    legend.text = element_text(size = 8)  
  ) +
  labs(fill = "Stock")

ggsave("faceted_pie_charts_P3.pdf", faceted_pie_charts3, width = 10, height = 8)

calculate_metrics3 <- function(returns, weights) {
  portf_returns <- Return.portfolio(returns, weights = weights)
  mean_ret <- mean(portf_returns)
  std_dev <- sd(portf_returns)
  sharpe <- SharpeRatio(portf_returns, Rf = 0, p = 0.95, FUN = "StdDev")
  sortino <- SortinoRatio(portf_returns)
  return(c(mean_ret, std_dev, sharpe, sortino))
}

metrics3 <- data.frame(
  Metric = c("Mean Return", "Standard Deviation", "Sharpe Ratio", "Sortino Ratio"),
  MinVar = calculate_metrics3(returns3, extractWeights(opt_minvar3)),
  CVaR = calculate_metrics3(returns3, extractWeights(opt_cvar3)),
  RiskParity = calculate_metrics3(returns3, opt_riskparity3),
  MaxSharpe = calculate_metrics3(returns3, extractWeights(opt_maxsharpe3)),
  EqualWeight = calculate_metrics3(returns3, equal_weights3)
)

writeLines(kable(metrics3, format = "latex", booktabs = TRUE) %>% kable_styling(), "metrics_table_P3.tex")

colnames(returns3) <- gsub(".Adjusted", "", colnames(returns3))


ef_cvar3 <- create.EfficientFrontier(R = returns3, portfolio = portf_cvar3, type = "mean-StdDev", n.portfolios = 25)

ef_data3 <- data.frame(
  Return = ef_cvar3$frontier[, "mean"],
  Risk = ef_cvar3$frontier[, "StdDev"]
)

asset_data3 <- data.frame(
  Return = colMeans(returns3),
  Risk = apply(returns3, 2, sd),
  Label = colnames(returns3)
)

asset_data3$Label <- gsub(".IS.Adjusted", "", asset_data3$Label, fixed = TRUE)

p3 <- ggplot() +
  geom_line(data = ef_data3, aes(x = Risk, y = Return), color = "#0000ff", size = 1.2) +  # Efficient Frontier line
  geom_point(data = asset_data3, aes(x = Risk, y = Return, color = Label), size = 4) +  # Asset points
  geom_text(data = asset_data3, aes(x = Risk, y = Return, label = Label), vjust = -1, hjust = 1, size = 3) +  # Asset labels
  geom_point(aes(x = 0.024135601, y = 0.002220837), color = "#E63946", size = 5, shape = 18) +  # MinVar portfolio point
  annotate("text", x = 0.024135601, y = 0.002220837, label = "MinVar", 
           vjust = -1.5, hjust = 1, color = "#E63946", size = 3.5) +
  geom_point(aes(x = 0.02598977, y = 0.00340929), color = "#E63946", size = 5, shape = 18) +  # CVaR portfolio point
  annotate("text", x = 0.02598977, y = 0.00340929, label = "CVaR", 
           vjust = -1.5, hjust = 1, color = "#E63946", size = 3.5) +
  geom_point(aes(x = 0.026540567, y = 0.002984478), color = "#E63946", size = 5, shape = 18) +  # RiskParity portfolio point
  annotate("text", x = 0.026540567, y = 0.002984478, label = "RiskParity", 
           vjust = -1.5, hjust = 1, color = "#E63946", size = 3.5) +
  geom_point(aes(x = 0.030130633, y = 0.004613532), color = "#E63946", size = 5, shape = 18) +  # MaxSharpe portfolio point
  annotate("text", x = 0.030130633, y = 0.004613532, label = "MaxSharpe", 
           vjust = -1.5, hjust = 1, color = "#E63946", size = 3.5) +
  geom_point(aes(x = 0.027002324, y = 0.003038498), color = "#E63946", size = 5, shape = 18) +  # EqualWeight portfolio point
  annotate("text", x = 0.027002324, y = 0.003038498, label = "EqualWeight", 
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

print(p3)
ggsave("efficient_frontier_P3.pdf", p3, width = 10, height = 8)

cvar_weights3 <- extractWeights(opt_cvar3)

weights_data3 <- data.frame(
  Asset = colnames(returns3),
  Weight = cvar_weights3
)
colnames(weights_data3)[0] <- "name"
weights_data3$Asset <- sub("\\.IS\\.Adjusted$", "", weights_data3$Asset)

colour_mapping3 <- c(
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

w3 <- chart.EF.Weights(ef_cvar3, match.col = "StdDev", colorset = colour_mapping3, legend.loc = "right", cex.legend = 1)
print(w3)

ggsave("cvar_weights_P3.pdf", w3, width = 10, height = 8)