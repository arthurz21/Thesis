# Thesis
ESC499Y1-Thesis-OMI-detection

Designed and evaluated a clinical ML workflow for automated ECG-based detection of Occlusion Myocardial Infarction;
achieved 87.6% sensitivity, 98.9% NPV (rule-out), and 97.3% specificity, 54.6% PPV (rule-in) on external testing dataset


• Extracted 73 ECG features using advanced signal processing techniques with Vectorcardiography (VCG) and
Principal Component Analysis (PCA) in MATLAB.

• Trained and evaluated 10 machine learning models (GBM, RF, KNN, ANN, XGB, SVM, LogReg, LDA, SGDLogReg, Gaussian-NB). Best algorithm outperforms the relevant paper published on nature medicine by 3%.

• Developed a clinical deployment ready windows application using PyQt5, capable of processing more than 400
ECGs at once and real time inference for OMI classifications.