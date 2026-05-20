# Betting Against Beta-strategi — Japan

Dette repository indeholder R-koden, der bruges til at konstruere og analysere en Betting Against Beta-strategi (BAB) for det japanske aktiemarked.

Den empiriske tilgang er baseret på Frazzini og Pedersen (2014), *Betting Against Beta*, og anvendes på japanske aktiedata.

## Projektoversigt

Scriptet udfører følgende trin:

1. Indlæser og renser daglige japanske aktiedata og markedsafkast.
2. Konstruerer ugentlige aktie- og markedsafkast i excess return-form.
3. Estimerer aktiespecifikke betaer ved hjælp af rullende korrelationer og volatilitet.
4. Anvender beta-shrinkage og laggede betaestimater for at undgå look-ahead bias.
5. Konstruerer BAB-porteføljen med et lav-beta-ben og et høj-beta-ben.
6. Evaluerer performance ved hjælp af excess returns, CAPM-alpha, Newey-West t-statistikker, Sharpe ratio, Sortino ratio, volatilitet og maximum drawdown.
7. Producerer decilporteføljeresultater i stil med Frazzini og Pedersens empiriske tabeller.
8. Udfører flere robusthedstests, herunder:
   - alternative beta-estimeringsvinduer,
   - alternative shrinkage-parametre,
   - kvintil- og decilbaseret porteføljekonstruktion,
   - equal-weighted og value-weighted porteføljer,
   - analyse af handelsomkostninger,
   - leverage-analyse,
   - likviditets- og large-cap-univers tests,
   - delperiodeanalyse.

## Data

Datafilerne er **ikke inkluderet** i dette repository.

Scriptet forventer følgende inputfiler:

```text
JPNall.csv
JPN_market.csv
