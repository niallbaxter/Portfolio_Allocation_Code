# Example Correlation Graph - With the current dates it does the whole period.
install.packages(c("quantmod", "tidyverse", "corrplot"))

library(quantmod)
library(tidyverse)
library(corrplot)

tickers <- c("AKBNK.IS", "GARAN.IS", "ISCTR.IS", "PETKM.IS", "TUPRS.IS", 
             "EREGL.IS", "TCELL.IS", "BIMAS.IS", "KCHOL.IS", "SAHOL.IS", 
             "THYAO.IS", "PGSUS.IS", "ASELS.IS", "KOZAL.IS")

start_date <- "2019-09-01"
end_date <- "2024-09-01"

get_stock_data <- function(ticker) {
  tryCatch({
    data <- getSymbols(ticker, src = "yahoo", from = start_date, to = end_date, auto.assign = FALSE)
    Cl(data)
  }, error = function(e) {
    message(paste("Error fetching data for:", ticker))
    return(NULL)
  })
}

stock_list <- map(tickers, get_stock_data)
stock_list <- stock_list[!sapply(stock_list, is.null)]
stock_data <- do.call(merge, stock_list)
colnames(stock_data) <- tickers[1:ncol(stock_data)]

simple_returns <- na.omit(ROC(stock_data, type = "discrete"))
cor_matrix <- cor(simple_returns, use = "complete.obs")
testRes <- cor.mtest(simple_returns, conf.level = 0.95)

corrplot(cor_matrix, 
         method = "color",
         type = "lower",
         p.mat = testRes$p,  
         insig = "blank",    
         addCoef.col = "black",
         number.cex = 1.2,
         tl.cex = 1
         ,
         tl.col = "black",
         tl.srt = 45,
         mar = c(0, 0, 1, 0),
         diag = FALSE)
