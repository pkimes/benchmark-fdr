# Summary

For 16S microbiome datasets, many studies do univariate comparisons across
case and control groups.

# Case control (univariate effect) datasets

I've collected and processed the raw data for many of these case-control
datasets, which are available on [Zenodo](https://zenodo.org/record/840333).
We can use these to directly calculate everything we need for the FDR methods.

Covariates could include:
- mean/median read depth of each OTU across all samples
- mean/median ubiquity of each OTU across all samples
- biological metadata like age, BMI, etc?

I think we can pick the largest study from each well-characterized disease:
diarrhea, inflammatory bowel disease, colorectal cancer, and perhaps obesity.

More info about the available datasets is available [here](https://github.com/cduvallet/microbiomeHD/blob/master/final/tables/table2.dataset_info_supplement.md).

# Other options

## Heritability

[Goodrich et al (2014)](http://www.sciencedirect.com/science/article/pii/S0092867414012410)
looked at the heritability of taxa. The raw genotype data is very difficult
to get, but it looks like they published the p-values. They also did their
analysis on 3 different datasets, and report the results for all 3. Two of the
datasets did not have enough power to find anything (I think).

They report:
- `A`: their measure of heritability
- `permutation p value`: the permutated p-value of A
- `A 95% CI`: 95% confidence interval on A (via permutation)

Interestingly, this data was a subset of the samples published in a follow-up
paper ([Goodrich et al, 2016](https://linkinghub.elsevier.com/retrieve/pii/S1931-3128(16)30153-6)).
The follow up essentially did the same analyses on the larger, full dataset.
This could be an interesting test case to see if these FDR methods improved
our ability to pick out the signals from the subsetted data that were later
found to be significant in the full dataset.

> Tripling the sample narrowed the confidence intervals around heritability estimates and uncovered additional heritable taxa, some of which are validated in other studies.

**To check**: are all comparisons included in Table S2, or just the ones
with qvalue < 0.1?

> The traits with a q value < 0.1 are presented in Table S2B.
