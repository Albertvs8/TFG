# Implementation of a privacy preserving KNN algorithm based on FHE

This repository contains the code developed during my Final Degree Project.

## Content

This repository is structured in three folders:

### Charts/

It contains two notebooks:
- `chart_compositions.ipynb`: Notebook used to understand the comparison algorithms.
- `chart_results.ipynb`: Notebook used to create the plots for analyzing the time complexity of the KNN algorithm.

### Data/

It contains three files:
- `heart_failure_clinical_records_dataset.csv`: Dataset containing all the original 299 samples.
- `heart_failure_clinical_records_dataset_16samples.csv`: Dataset containing 16 samples, used to compare the memory size between unencrypted data and encrypted data.
- `encrypted_dataset.md`: Link to the files containing 16 encrypted samples. Due to storage issues, these files are not uploaded to the repository.

### Scripts/

It contains four scripts:
- `Comparison.cpp`: Script to test the comparison function.
- `MergeSort.cpp`: Script to test the sorting function.
- `SaveEncryptedDataset.pp`: Script used to save encrypted data to memory.
- `KNN.cpp`: Script containing the full implementation of the KNN algorithm.
