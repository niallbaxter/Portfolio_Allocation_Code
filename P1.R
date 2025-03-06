install.packages("quantmod")
install.packages("PerformanceAnalytics")
install.packages("PortfolioAnalytics")
install.packages("DEoptim")
install.packages("ggplot2")
install.packages("gridExtra")
install.packages("ggpubr")
install.packages("latex2exp")
install.packages("CVXR")
install.packages("kableExtra")
install.packages("riskParityPortfolio")

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
start_date <- "2019-09-01"
end_date <- "2020-03-11"

getSymbols(stocks, from = start_date, to = end_date, src = "yahoo")
prices <- do.call(merge, lapply(stocks, function(x) Ad(get(x))))
returns <- na.omit(Return.calculate(prices))

stats_df1 <- data.frame(
  Stock = stocks,
  Mean_Return = colMeans(returns, na.rm = TRUE),
  Std_Dev = apply(returns, 2, sd, na.rm = TRUE),
  Kurtosis = apply(returns, 2, kurtosis, na.rm = TRUE),
  Skewness = apply(returns, 2, skewness, na.rm = TRUE)
)
SS1 <- stats_df1 %>%
  kable("latex", booktabs = TRUE, caption = "Summary Statistics for Stocks Period 1")

writeLines(SS1, "SS1")

portf <- portfolio.spec(assets = stocks)
portf <- add.constraint(portf, type = "full_investment")
portf <- add.constraint(portf, type = "long_only")

portf_minvar <- portfolio.spec(assets = stocks)
portf_minvar <- add.constraint(portf_minvar, type = "full_investment")
portf_minvar <- add.constraint(portf_minvar, type = "long_only")
portf_minvar <- add.objective(portf_minvar, type = "risk", name = "var")
portf_minvar$constraints[[1]]$min_sum <- 0.99
portf_minvar$constraints[[1]]$max_sum <- 1.01

portf_cvar <- portfolio.spec(assets = stocks)
portf_cvar <- add.constraint(portf_cvar, type = "full_investment")
portf_cvar <- add.constraint(portf_cvar, type = "long_only")
portf_cvar <- add.objective(portf_cvar, type = "risk", name = "CVaR")
portf_cvar$constraints[[1]]$min_sum <- 0.99
portf_cvar$constraints[[1]]$max_sum <- 1.01

# Risk Parity Optimization
risk_parity_optimization <- function() {
  S <- cov(returns)
  return(riskParityPortfolio(S)$w)
}

portf_maxsharpe <- portfolio.spec(assets = stocks)
portf_maxsharpe <- add.constraint(portf_maxsharpe, type = "full_investment")
portf_maxsharpe <- add.constraint(portf_maxsharpe, type = "long_only")
portf_maxsharpe <- add.objective(portf_maxsharpe, type = "return", name = "mean")
portf_maxsharpe <- add.objective(portf_maxsharpe, type = "risk", name = "var")
portf_maxsharpe$constraints[[1]]$min_sum <- 0.99
portf_maxsharpe$constraints[[1]]$max_sum <- 1.01

equal_weights <- rep(1 / length(stocks), length(stocks))

opt_minvar <- optimize.portfolio(returns, portf_minvar, optimize_method = "DEoptim")
opt_cvar <- optimize.portfolio(returns, portf_cvar, optimize_method = "DEoptim")
opt_riskparity <- risk_parity_optimization()
opt_maxsharpe <- optimize.portfolio(returns, portf_maxsharpe, optimize_method = "CVXR", maxSR = TRUE)

normalize_weights <- function(weights) {
  return(weights / sum(weights))
}

weights <- data.frame(
  Stock = stocks,
  MinVar = normalize_weights(extractWeights(opt_minvar)),
  CVaR = normalize_weights(extractWeights(opt_cvar)),
  RiskParity = normalize_weights(opt_riskparity),
  MaxSharpe = normalize_weights(extractWeights(opt_maxsharpe)),
  EqualWeight = equal_weights
)

writeLines(kable(weights, format = "latex", booktabs = TRUE) %>% 
             kable_styling(), "weights_table_P1.tex")

custom_palette <- c(
  "#E63946", "#F1FAEE", "#A8DADC", "#457B9D", "#1D3557",
  "#F4A261", "#2A9D8F", "#264653", "#D9BF77", "#9C27B0",
  "#F7C59F", "#E76F51", "#9A174D", "#4A4E69"
)

combined_df <- data.frame()
for (method in colnames(weights)[-c(1, 6)]) {  # Exclude the "Stock" column and "EqualWeight" column
  temp_df <- data.frame(Stock = weights$Stock, Weight = weights[[method]], Method = method)
  combined_df <- rbind(combined_df, temp_df)
}

faceted_pie_charts <- ggplot(combined_df, aes(x = "", y = Weight, fill = Stock)) +
  geom_bar(stat = "identity", width = 1, color = "black", size = 0.5) +  
  coord_polar("y", start = 0) +
  facet_wrap(~Method, ncol = 2) +  
  theme_void() +
  scale_fill_manual(values = custom_palette) +  
  theme(
    strip.text = element_text(size = 12, face = "bold"),  
    legend.position = "bottom",  
    legend.title = element_text(size = 10, face = "bold"),  
    legend.text = element_text(size = 8)  
  ) +
  labs(fill = "Stock")

ggsave("faceted_pie_charts_P1.pdf", faceted_pie_charts, width = 10, height = 8)

print(faceted_pie_charts)

calculate_metrics <- function(returns, weights) {
  portf_returns <- Return.portfolio(returns, weights = weights)
  mean_ret <- mean(portf_returns)
  std_dev <- sd(portf_returns)
  sharpe <- SharpeRatio(portf_returns, Rf = 0, p = 0.95, FUN = "StdDev")
  sortino <- SortinoRatio(portf_returns)
  return(c(mean_ret, std_dev, sharpe, sortino))
}

metrics <- data.frame(
  Metric = c("Mean Return", "Standard Deviation", "Sharpe Ratio", "Sortino Ratio"),
  MinVar = calculate_metrics(returns, extractWeights(opt_minvar)),
  CVaR = calculate_metrics(returns, extractWeights(opt_cvar)),
  RiskParity = calculate_metrics(returns, opt_riskparity),
  MaxSharpe = calculate_metrics(returns, extractWeights(opt_maxsharpe)),
  EqualWeight = calculate_metrics(returns, equal_weights)
)

writeLines(kable(metrics, format = "latex", booktabs = TRUE) %>% kable_styling(), "metrics_table_P1.tex")

weights <- data.frame(
  Stock = stocks,
  MinVar = normalize_weights(extractWeights(opt_minvar)),
  CVaR = normalize_weights(extractWeights(opt_cvar)),
  RiskParity = normalize_weights(opt_riskparity),
  MaxSharpe = normalize_weights(extractWeights(opt_maxsharpe)),
  EqualWeight = equal_weights
)

writeLines(kable(weights, format = "latex", booktabs = TRUE) %>% 
             kable_styling(), "weights_table_P1.tex")

custom_palette <- c(
  "#E63946", "#F1FAEE", "#A8DADC", "#457B9D", "#1D3557",
  "#F4A261", "#2A9D8F", "#264653", "#D9BF77", "#9C27B0",
  "#F7C59F", "#E76F51", "#9A174D", "#4A4E69"
)

combined_df <- data.frame()
for (method in colnames(weights)[-c(1, 6)]) {  # Exclude the "Stock" column and "EqualWeight" column
  temp_df <- data.frame(Stock = weights$Stock, Weight = weights[[method]], Method = method)
  combined_df <- rbind(combined_df, temp_df)
}

faceted_pie_charts <- ggplot(combined_df, aes(x = "", y = Weight, fill = Stock)) +
  geom_bar(stat = "identity", width = 1, color = "black", size = 0.5) +  
  coord_polar("y", start = 0) +
  facet_wrap(~Method, ncol = 2) +  
  theme_void() +
  scale_fill_manual(values = custom_palette) +  
  theme(
    strip.text = element_text(size = 12, face = "bold"),  
    legend.position = "bottom",  
    legend.title = element_text(size = 10, face = "bold"),  
    legend.text = element_text(size = 8)  
  ) +
  labs(fill = "Stock")

ggsave("faceted_pie_charts_P1.pdf", faceted_pie_charts, width = 10, height = 8)

print(faceted_pie_charts)

calculate_metrics <- function(returns, weights) {
  portf_returns <- Return.portfolio(returns, weights = weights)
  mean_ret <- mean(portf_returns)
  std_dev <- sd(portf_returns)
  sharpe <- SharpeRatio(portf_returns, Rf = 0, p = 0.95, FUN = "StdDev")
  sortino <- SortinoRatio(portf_returns)
  return(c(mean_ret, std_dev, sharpe, sortino))
}

metrics <- data.frame(
  Metric = c("Mean Return", "Standard Deviation", "Sharpe Ratio", "Sortino Ratio"),
  MinVar = calculate_metrics(returns, extractWeights(opt_minvar)),
  CVaR = calculate_metrics(returns, extractWeights(opt_cvar)),
  RiskParity = calculate_metrics(returns, opt_riskparity),
  MaxSharpe = calculate_metrics(returns, extractWeights(opt_maxsharpe)),
  EqualWeight = calculate_metrics(returns, equal_weights)
)

writeLines(kable(metrics, format = "latex", booktabs = TRUE) %>% kable_styling(), "metrics_table_P1.tex")

colnames(returns) <- gsub(".Adjusted", "", colnames(returns))

ef_cvar <- create.EfficientFrontier(R = returns, portfolio = portf_cvar, type = "mean-StdDev", n.portfolios = 25)

ef_data <- data.frame(
  Return = ef_cvar$frontier[, "mean"],
  Risk = ef_cvar$frontier[, "StdDev"]
)

asset_data <- data.frame(
  Return = colMeans(returns),
  Risk = apply(returns, 2, sd),
  Label = colnames(returns)
)

asset_data$Label <- gsub(".IS", "", asset_data$Label, fixed = TRUE)

custom_palette <- c(
  "#E63946", "#F1FAEE", "#A8DADC", "#457B9D", "#1D3557",
  "#F4A261", "#2A9D8F", "#264653", "#D9BF77", "#9C27B0",
  "#F7C59F", "#E76F51", "#9A174D", "#4A4E69"
)

color_mapping <- setNames(custom_palette, colnames(stocks))

p1 <- ggplot() +
  geom_line(data = ef_data, aes(x = Risk, y = Return), color = "#0000ff", size = 1.2) +  # Efficient Frontier line
  geom_point(data = asset_data, aes(x = Risk, y = Return, color = Label), size = 4) +  # Asset points
  geom_text(data = asset_data, aes(x = Risk, y = Return, label = Label), vjust = -1, hjust = 1, size = 3) +  # Asset labels
  geom_point(aes(x = 0.0105643629, y = 0.0005590249 ), color = "#E63946", size = 5, shape = 18) +  # CVaR portfolio point
  annotate("text", x = 0.0105643629, y = 0.0005590249 , label = "MinVar", 
           vjust = -1.5, hjust = 1, color = "#E63946", size = 3.5) +
  geom_point(aes(x = 0.011270887, y = 0.001012113 ), color = "#E63946", size = 5, shape = 18) +  # CVaR portfolio point
  annotate("text", x = 0.011270887, y = 0.001012113 , label = "CVaR", 
           vjust = -1.5, hjust = 1, color = "#E63946", size = 3.5) +
  geom_point(aes(x = 0.0128081650  , y = 0.0001771459  ), color = "#E63946", size = 5, shape = 18) +  # CVaR portfolio point
  annotate("text", x = 0.0128081650  , y = 0.0001771459  , label = "RiskParity", 
           vjust = -1.5, hjust = 1, color = "#E63946", size = 3.5) +
  geom_point(aes(x = 0.014435946   , y = 0.002146111   ), color = "#E63946", size = 5, shape = 18) +  # CVaR portfolio point
  annotate("text", x = 0.014435946   , y = 0.002146111   , label = "MaxSharpe", 
           vjust = -1.5, hjust = 1, color = "#E63946", size = 3.5) +
  geom_point(aes(x = 1.430239e-02   , y = 7.011064e-05   ), color = "#E63946", size = 5, shape = 18) +  # CVaR portfolio point
  annotate("text", x = 1.430239e-02   , y = 7.011064e-05   , label = "EqualWeight", 
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

print(p1)

cvar_weights = extractWeights(opt_cvar)

cvar_weights <- opt_cvar&weights


weights_data <- data.frame(
  Asset = colnames(returns),
  Weight = cvar_weights
)
colnames(weights_data)[0] <- "name"
weights_data$Asset <- sub("\\.IS\\.Adjusted$", "", weights_data$Asset)

colour_mapping <- c(
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

w1 <- chart.EF.Weights(ef_cvar, match.col = "StdDev", colorset = colour_mapping, legend.loc = "right", cex.legend = 0.8)
