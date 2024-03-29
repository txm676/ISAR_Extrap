README
================
Tom Matthews
2023-08-04

## Preamble

Presented here is the code, data and (updated) results for the paper:

Matthews, T.J. & Aspin, T.W.H. (2019) Model averaging fails to improve
the extrapolation capability of the island species–area relationship.
Journal of Biogeography, 46, 1558–1568.

For the study aim and methods, see the original paper.

Through subsequent analyses and data checking, we have corrected small
data entry errors where they have been found. In a small number of
cases, we have subsequently found that the source paper authors
published updated / more extensive versions of the datasets used; in
these cases, we have replaced the original version. As such, the
datasets presented here differ slightly from those used in the original
publication. In addition, since the paper was published in 2019, there
have been numerous updates made to the ‘sars’ R package, including
replacing models, bug fixes and improving the fitting algorithms. As
such, the results presented here may differ slightly from those
published in the original paper, although the main conclusions are
unchanged.

One thing we changed in the code relates to how the power model checks
are done. Originally we ran these checks (e.g. residual normality,
z-value significance) on the full dataset rather than the filtered
dataset used to predict the richness of the largest islands. In
hindsight, it arguably makes more sense to check the power model
residual assumptions using the filtered dataset. As such, we run this
sensitivity test using the filtered datasets (th = 0.5).

The R script and source code are stored in the ‘R code’ directory. The
datasets are stored in a list called ‘Datasets’ (saved in the drive
“Data” within the repository), where each element is a dataframe
(representing a dataset) with two columns: area (a) and number of
species (s). Area is in ha. Dataset numbers match the filenames in the
‘filenames.rds’ file. The “Data” directory also holds a csv file called
“preds” which contains the taxon and latitude data for each study, which
is used in the GAM analyses. Finally, there is a word document which
provides meta-data for each of the sourced datasets.

``` r
setwd("C:\\Users\\Tom\\Desktop\\ISAR_Extrap")

Datasets <- readRDS("Data\\Datasets.rds")

length(Datasets)
```

    ## [1] 120

``` r
filenames <- readRDS("Data\\filenames.rds")
filenames[1]
```

    ## [1] "CSV-Baldi & Kisbenedek (1999)noniso.csv"

## Results

### POW vs. MMI

When a th value of 0.5 was used, the power model provided the best fit
to the most (filtered) datasets (n = 24), followed by the linear model
(n = 21), and then the negative exponential (n = 17), logarithmic (n =
17) and Monod (n = 15) models, according to AICc.

For one dataset, no models could be successfully fitted (i.e. no models
passed the residual assumption tests) and this dataset was removed from
the comparisons. In contrast with our hypothesis, the power model
provided the most accurate prediction of the richness of the largest
island (i.e. the lowest absolute LEE value) in 67 cases (56%), with the
multi‐model averaged curve providing the more accurate prediction in the
remaining 52 cases (44%). The median LEE value of the power model was
0.04 (95% quantiles = −0.30 and 0.31), whilst the median LEE value of
the multi‐model curve (LEE‐MMI) was 0.02 (95% quantiles = −0.34 and
1.56). However, as LEE values could be both positive and negative, the
median of the absolute LEE values provides a better summary of the
extrapolation capability: the median of absolute LEE‐POW values was 0.08
(95% quantiles = 0.01 and 0.35), whilst the median of absolute LEE‐MMI
values was 0.10 (95% quantiles = 0.002 and 1.56). Both the power model
(64% of LEE‐POW values were positive) and the multi‐model averaged curve
(57% of LEE‐MMI values were positive) had a slightly greater tendency to
over predict the observed richness.

The confidence intervals were on average narrower for the power model
predictions (median 95% confidence interval width = 18) than for the
multi‐model averaged curve predictions (median 95% confidence interval
width = 125). The confidence intervals around the multi‐model averaged
curve predictions were sometimes very large.

### Predictions from 20 individual models

When the extrapolation predictions from all 20 ISAR models were
considered, in addition to the multi‐model averaged curve, the power
model provided the most accurate prediction of the richness of the
largest island in 10 cases, with the multi‐model averaged curve
providing the most accurate prediction in nine cases. The negative
exponential model provided the best prediction the most times, with 12
cases (the results for all models are provided in Table 1), with the
Extended Power 2 model providing the best prediction in 8 cases.

### GAMs

The full GAM (i.e. the GAM with all predictors) had a lower AIC value
(−114.5) than an equivalent standard linear regression model (−112.7);
this provides additional justification for our use of GAMs. When LEE‐POW
was used as the response variable in a GAM model selection analysis, the
best model(s) contained Ascale, Sscale, Lat. and Smin (Table 2). In the
published analysis, these four variables were all in the best model,
whereas now, while they were all included in the models with delta AICc
\<=2, the best model only contained Ascale and Sscale. A plot of the
smoothers for these four Variables, when all fitted in the same model,
is provided as Figure 2. The effective degrees of freedom of the
smoothers for Ascale and Lat. were one, indicating that these smoothers
were straight lines; increasing Ascale resulted in decreasing LEE‐POW,
while the opposite pattern was true for Lat (Figure 2). A difference
with the published results was that the effective degrees of freedom of
the smoothers for Sscale was now also linear, although the general
pattern was the same (increasing Sscale resulted in an increase in
LEE‐POW). The Smin relationship was still more complex (Figure 2).
However, there was a reasonable degree of model uncertainty as the best
model had an AICc weight of only 0.16, and there were three additional
models within 2 delta AICc units of the best model (Table 2). In
addition, the adjusted R2 value of the model with all four of the
aforementioned variables was low (0.15). Ascale (0.98), Sscale (0.99)
had quite high relative importance values, whilst the values for the
remaining predictors were all lower (Table 2).

For 18 of the ISAR models, the relative fit of a model to the filtered
dataset (i.e. the model’s AICc weight) was a poor predictor of a model’s
extrapolation accuracy (measured using the LEE metric). In only two
cases (for the logistic and Extended Power 1 models) was the AICc weight
a significant predictor of a model’s absolute LEE value.

### Sensitivity analyses

The choice of th value did not change the overall qualitative results.
The power model provided the more accurate prediction in 66 (55%) and 66
(56%; when a th value of 0.7 was used there were three datasets for
which no models could be successfully fitted) cases when th values of
0.3 and 0.7 were used, respectively. In regards to the power model
validation sensitivity test, there were 55 datasets (up from 23 in the
published analyses) for which the power model failed one of the
validation checks. The difference is due to us now running the checks
using the filtered datasets rather than full datasets. Most of these
failed checks are a result of non-significant z-values. If this check is
removed, only 19 datasets fail the validation checks (i.e. the residual
assumption checks). When removing these 55 datasets and re‐running the
prediction analysis using the remaining 65 datasets, the power model
still provided the most best-fits to the filtered datasets (18 vs. 10
for the logarithmic in second place). The power model provided the most
accurate prediction of the richness of the largest island (i.e. the
lowest absolute LEE value) in 46% of cases, with the multi‐model
averaged curve providing the more accurate prediction in 54% of cases.
This differed slightly from the main results, where the power model
provided the most accurate prediction the majority of times. However,
the relative performance is still roughly 50:50 and thus the main
conclusions are unchanged. If the validation sensitivity test is instead
run using the full dataset, as in the original published paper, the
qualitative results are unchanged: the power model provided the best
prediction in 54% of cases compared to 46% for the multi-model averaged
curve. Finally, rerunning the prediction analysis after excluding the
linear model from the multi‐model averaged curve resulted in a slight
increase in the number of cases where the multi‐model averaged curve
provided the more accurate prediction (55 cases), but the general
picture remained the same.

#### Table 1

TABLE 1 The 20 models that were fitted to generate the multi‐model
averaged ISAR curve. Mean weight is the mean AICc weight for a given
model across all fits to the filtered datasets (excluding
non‐satisfactory fits). Best fit corresponds the number of times a model
provided the best fit to a filtered dataset. Best prediction corresponds
to the number of times a model provided the best extrapolated prediction
in the all model comparison; these values do not sum to 120 (the number
of datasets) as the multi‐model averaged curve provided the best
extrapolation prediction in several cases.

| Model                  | No..parameters | Mean.weight | Best.fit | Best.prediction |
|:-----------------------|---------------:|------------:|---------:|----------------:|
| Asymptotic             |              3 |        0.04 |        1 |               4 |
| Beta-P                 |              4 |        0.00 |        0 |               5 |
| Chapman_Richards       |              3 |        0.03 |        0 |               2 |
| Logarithmic            |              2 |        0.14 |       17 |               6 |
| Extended Power 1       |              3 |        0.03 |        0 |               5 |
| Extended Power 2       |              3 |        0.03 |        1 |               8 |
| Gompertz               |              3 |        0.04 |        0 |               3 |
| Kobayashi              |              2 |        0.15 |       12 |               5 |
| Linear                 |              2 |        0.13 |       21 |               7 |
| Logistic               |              3 |        0.05 |        5 |               6 |
| Monod                  |              2 |        0.11 |       15 |               6 |
| Heleg                  |              3 |        0.03 |        0 |               4 |
| Negative Exponential   |              2 |        0.10 |       17 |              12 |
| Persistence Function 1 |              3 |        0.04 |        3 |               2 |
| Persistence Function 2 |              3 |        0.03 |        1 |               4 |
| Power                  |              2 |        0.16 |       24 |              10 |
| Power Rosenzweig       |              3 |        0.03 |        1 |               7 |
| Rational               |              3 |        0.03 |        1 |               7 |
| Weibull-3              |              3 |        0.03 |        0 |               4 |
| Weibull-4              |              4 |        0.01 |        0 |               3 |

#### Table 2

TABLE 2 The results of the generalized additive model selection.

| Model | Amin | Ascale | Lat. |   Ni | Smin | Sscale | Taxon | Delta | Weight |
|:------|-----:|:-------|:-----|-----:|:-----|:-------|------:|------:|-------:|
| 1     |   NA | \+     |      |   NA |      | \+     |    NA |  0.00 |   0.16 |
| 2     |   NA | \+     | \+   |   NA |      | \+     |    NA |  0.63 |   0.12 |
| 3     |   NA | \+     |      |   NA | \+   | \+     |    NA |  0.89 |   0.10 |
| 4     |   NA | \+     | \+   |   NA | \+   | \+     |    NA |  1.88 |   0.06 |
| RI    | 0.26 | 0.98   | 0.4  | 0.25 | 0.37 | 0.99   |  0.19 |    NA |     NA |

#### Figure 1

![](README_files/figure-gfm/Fig1-1.png)<!-- -->

FIGURE 1 The varying species richness predictions of five ISAR models.
Each of the five models (see Table 1) was fitted to a simulated
archipelago consisting of eight islands of varying size (1, 3, 7, 14,
17, 22, 26 and 30; undefined units) and richness (3, 7, 14, 18, 20, 23,
24 and 25). These model fits were then used to predict the richness of
an island of size 80 (grey dotted line)

#### Figure 2

![GAM_fig.jpeg](README_files/GAM_fig.jpeg) FIGURE 2 Fitted smoothers
from a generalized additive model showing the partial effects of Ascale,
Latitude, Smin and Sscale on the LEE‐POW values. The fitted values have
been shifted in each plot by adding the model intercept (0.04) value
(using the shift argument in the plot.gam R function). The effective
degrees of freedom for each smoother are: Ascale (1.00), Latitude
(1.00), Smin (2.90) and Sscale (1.00). The dashed lines represent the
standard error curves (two SE above and below). Each LEE‐POW value
relates to the accuracy of a prediction of the number of species on a
habitat island using the power model. For each of 119 habitat island
datasets, the largest island and all islands larger than half the size
of the largest island were removed and the power model was fitted to the
filtered dataset and extrapolated.
