# Implementation of a privacy preserving KNN algorithm based on FHE

**Author: Albert Viladot Saló**

This repository comprises the code that was developed as part of my Final Degree Project in Mathematical Engineering in Data Science at Universitat Pompeu Fabra.

*Fully Homomorphic Encryption (FHE) is an encryption technique that allows computation to be performed on encrypted data without the need to decrypt it first. This has the potential to revolutionize data security and privacy, and ongoing research is focused on improving its efficiency and practicability.*

*One limitation of FHE schemes is that, typically, they only support additions and multiplications. More complex operations, such as real number comparison, require careful design of new algorithms. In this project, we explore a recent approach based on polynomial composition for comparing two real numbers in a homomorphically encrypted setting. We study how this method can be used to sort an array of real numbers using FHE. We then apply these concepts to design and implement a fully homomorphic version of the popular machine learning algorithm known as k-nearest neighbours (KNN).*

*To implement the FHE KNN algorithm, we use the CKKS scheme and the OpenFHE library, an open-source software library developed in C++. Finally, we analyse the computational performance of our implementation to evaluate its efficiency, depending on the dataset size and the number of neighbours. Our results indicate that the time complexity increases exponentially with the dataset size. Reasonable computation times are achievable for small datasets, with up to 16 training samples.*


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
- `Comparison.cpp`: Script used to test the comparison function.
- `MergeSort.cpp`: Script used to test the sorting function.
- `SaveEncryptedDataset.pp`: Script used to save encrypted data to memory.
- `KNN.cpp`: Script containing the full implementation of the KNN algorithm.
