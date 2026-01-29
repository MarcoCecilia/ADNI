# ADNI
Clinical Characterisation of Alzheimer Disease

## Project Workflow

1. **MRI Download and Preprocessing**
   - The MRI data were downloaded from the ADNI website.
   - Preprocessing was performed using FastSurfer.
   - All preprocessing steps are documented in the notebook: `PreProcessingData.ipynb`.

2. **Exploratory Data Analysis**
   - An initial exploratory analysis of the dataset was conducted to understand the main characteristics and distributions of the variables after was done a normalization in `normalized_data.R`.
   - The analysis code is in `Explorative_Analysis.R`.
   - A detailed report is available in `Explorative_Analysis.docx`.

3. **Longitudinal Analysis Attempt**
   - We attempted a longitudinal analysis to model disease progression, but the results were not significant.
   - The corresponding code is in `MultiNomialMixedEffect.R`.

4. **Multinomial Model and Class Imbalance**
   - We developed an atemporal multinomial model and compared different approaches to address the dataset imbalance.
   - The code can be found in `multinomial_con_balancing.R`.
   - This is the central part of the project, with a report available in `Alz_Multinom_Report_Interpretato.docx`.

5. **Hierarchical Clustering on MRI Features**
   - We investigated whether using only numerical features derived from MRI scans was sufficient to discriminate between CN and MCI classes.
   - Unsupervised hierarchical clustering was performed in `hierarchical_clustering.R`.
   - The associated report is `hierarchical_clustering.docx`.

6. **Integration of Genetic Variables**
   - To assess whether adding genetic variables could improve model accuracy, we processed genetic data in `Genetic_data_process.R` and ran the analysis in `Genetic_Analysis.R`.

7. **Project Summary and Poster**
   - The final results of the project are summarized in the poster `Poster.pdf`.
   - A comprehensive summary of the work is available in `project_summary.docx`.

---

## File Structure

- `PreProcessingData.ipynb` – MRI data preprocessing
- `Explorative_Analysis.R` – Exploratory data analysis (with `Explorative_Analysis.docx` report)
- `MultiNomialMixedEffect.R` – Longitudinal analysis attempt
- `multinomial_con_balancing.R` – Multinomial modeling and class balancing (with `Alz_Multinom_Report_Interpretato.docx`)
- `hierarchical_clustering.R` – Hierarchical clustering on MRI features (with `hierarchical_clustering.docx`)
- `Genetic_data_process.R` – Genetic data preprocessing
- `Genetic_Analysis.R` – Genetic variable integration analysis
- `Poster.pdf` – Final project poster
- `project_summary.docx` – Project summary

---

## How to Reproduce the Analysis

1. Download the MRI dataset from the ADNI website.
2. Run the preprocessing notebook (`PreProcessingData.ipynb`) using FastSurfer.
3. Follow the analysis scripts in the order described above.
4. Review the associated reports for interpretation and discussion.

---

## Authors
Marco Cecilia marco.ceci.2002@gmail.com
Giorgia Arena giorgia.arena@mail.polimi.it
Cecilia Galbiati cecilia.galbiati@mail.polimi.it
Elena Fusco elena.fusco@mail.polimi.it



