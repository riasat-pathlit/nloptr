cat("\014")
rm(list = ls())
set.seed(2020)

iN <- 3
pDb <- "../../../ntu-pre-phd/1-price-database/m-px-db.csv"
iDb <- timeSeries::readSeries(pDb, sep = ",", format = "%Y-%m-%d")
iTs <- iDb[, seq(iN)]

rM <- timeSeries::returns(iTs)
wiR <- unname(colMeans(rM)) * 12
wiC <- unname(stats::cov(rM)) * 12

w0R <- 0.0115

w00 <- round(runif(1, 0, 10))
wi0 <- round(runif(iN, 0, 10))
wA0 <- c(w00, wi0)
w0 <- sum(wA0)

B <- 1.003
S <- 0.9965

r <- 1 + w0R
p <- 1 + wiR

R <- 24.8 # iN = 3: 24.8 # iN = 5: 27.98864 # w0 + mean(wA0) / 2
a <- c(r, p) %*% wA0 - R

objFn <- function(x) {
  (wi0 + x) %*% wiC %*% (wi0 + x)
}

eqFn <- function(x) {
  (r * (B * (x > 0) + S * (x < 0)) - p) %*% x - a
}

inFn <- function(x) {c(
  -sum(x) + w00,
  sum(x) + w0 - w00,
  wi0 + x
)}

x <- suppressMessages(nloptr::slsqp(
  rep(0, iN), objFn, heq = eqFn, hin = inFn, control = list("xtol_rel" = 1e-8)
)$par)

w0T <- r * (w00 - x %*% (B * (x > 0) + S * (x < 0)))
wNT <- p %*% (wi0 + x)
drop(w0T + wNT)
R
drop(sqrt(objFn(x)))
round(matrix(x), 4)
