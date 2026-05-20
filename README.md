# Betting Against Beta Strategy — Japan

This repository contains the R code used to construct and analyze a Betting Against Beta (BAB) strategy for the Japanese equity market.

The empirical approach is based on Frazzini and Pedersen (2014), *Betting Against Beta*, and applies the strategy to Japanese stock-level return data.

## Project overview

The script performs the following steps:

1. Loads and cleans daily Japanese stock data and market return data.
2. Constructs weekly stock and market excess returns.
3. Estimates stock-level betas using rolling correlations and volatilities.
4. Applies beta shrinkage and lagged beta estimates to avoid look-ahead bias.
5. Constructs the BAB portfolio using low-beta and high-beta legs.
6. Evaluates performance using excess returns, CAPM alpha, Newey-West t-statistics, Sharpe ratio, Sortino ratio, volatility, and maximum drawdown.
7. Produces decile portfolio results similar to Frazzini and Pedersen’s empirical tables.
8. Runs several robustness tests, including:
   - alternative beta estimation windows,
   - alternative shrinkage parameters,
   - quintile and decile portfolio construction,
   - equal-weighted and value-weighted portfolios,
   - transaction cost analysis,
   - leverage analysis,
   - liquidity and large-cap universe tests,
   - subperiod analysis.

## Data

The data files are **not included** in this repository.

The script expects the following input files:

```text
JPNall.csv
JPN_market.csv
