# 🌍 Multi-Model Bias Correction, GCM Evaluation, and Drought Risk Assessment Using Machine Learning,RIM,CPI,JOC and SPI under NEX-GDDP-CMIP6 Scenarios.

## 📌 Project Overview
This repository presents the computational framework developed for an MSc thesis in Hydrology and Water Resources Management. The study integrates **machine learning-based bias correction**, **multi-criteria GCM evaluation**, **statistical downscaling**, and **drought risk assessment** under CMIP6 climate change scenarios.

The framework improves the reliability of climate model outputs and supports robust hydrological impact assessment under historical and future climate conditions.

---

## 🧭 Scientific Motivation
Global Climate Models (GCMs) under CMIP6 provide essential climate projections; however, they exhibit systematic biases and scale mismatches when applied to basin-scale hydrological studies. This research addresses these limitations through a hybrid framework combining:

- Machine learning-based bias correction of precipitation
- Multi-criteria evaluation of GCM performance
- Basin-scale downscaling of selected GCM outputs
- Drought risk characterization using SPI indices under SSP scenarios

---

## 🎯 Research Objectives
1. To bias-correct CHIRPS rainfall data using four machine learning algorithms and develop a validated high-resolution baseline rainfall dataset.  
2. To evaluate and rank NEX-GDDP-CMIP6 GCMs using seven relative importance and performance metrics.  
3. To downscale the selected optimal GCM(s) to basin scale for improved spatial rainfall representation.  
4. To project future precipitation variability and assess drought risk under SSP2-4.5 and SSP5-8.5 scenarios using SPI analysis.  

---

## 🔬 Methodological Framework

### 1. Data Acquisition
- CMIP6 GCM datasets (NEX-GDDP-CMIP6)
- CHIRPS precipitation observations
- Basin-scale hydroclimatic datasets

### 2. Data Preprocessing
- Quality control and missing value handling

### 3. Bias Correction
Machine learning-based bias correction applied to CHIRPS precipitation using four models:
- RF
- SVM
- XGBOOST
- KNN

### 4. GCM Evaluation and Selection
Performance ranking based on:
- Nash–Sutcliffe Efficiency (NSE)
- Root Mean Square Error (RMSE)
- Mean Absolute Error (MAE)
- Correlation coefficient (R)
- Additional relative importance metrics

### 5. Downscaling
Spatial refinement of selected GCM outputs to basin scale for hydrological applicability.

### 6. Drought Risk Assessment
- Standardized Precipitation Index (SPI)
- Scenario-based analysis under SSP2-4.5 and SSP5-8.5
- Spatiotemporal drought characterization

---

## 🧰 Tools and Computational Environment

### Python
- NumPy, Pandas (data processing)
- Scikit-learn (machine learning models)
- TensorFlow / Keras (deep learning)
- Matplotlib (visualization)
- Rasterio (geospatial data handling)

### R
- raster (spatial analysis)
- ncdf4 (NetCDF handling)
- hydroGOF (hydrological model evaluation)
- ggplot2 (visualization)

---


---

## 📊 Expected Scientific Contributions
- Improved bias-corrected precipitation datasets for hydrological modeling  
- Robust framework for CMIP6 GCM selection at basin scale  
- Enhanced spatial rainfall representation through downscaling  
- Quantitative drought risk assessment under future climate scenarios  
- Integration of machine learning with hydroclimatic impact studies  

---

## 🌍 Study Area
The methodology is applied to a climatically sensitive Upper Blue Nile basin in Ethiopia characterized by high rainfall variability and increasing drought frequency.

---

## 🔁 Reproducibility Statement
All scripts are designed to ensure full reproducibility of results. The workflow can be executed sequentially after installing dependencies and preparing input datasets.

---

## 👨‍🎓 Author

**Name:** Teklemaryam Asnakew  
**Program:** MSc in Hydrology and Water Resources Management  
**Institution:** University of Gondar  
**Email:** teklemaryamas@gmail.com  

---

## 📌 Declaration
This repository is part of an MSc thesis research project and is intended for academic examination and scientific transparency. It is structured to meet reproducibility and open-science standards.

---
