# -*- mode: python ; coding: utf-8 -*-

block_cipher = None

a = Analysis(
    ['core_app.py'],  # Your main Python script
    pathex=[],
    binaries=[],
    datas=[
        # Include the MATLAB exe and model file
        ('ecg_processor.exe', '.'),
        ('xgb_model.pkl', '.'),
        (r'C:\Users\arthu\PycharmProjects\Thesis\.venv\Lib\site-packages\xgboost\lib\xgboost.dll', 'xgboost/lib'),
        (r'C:\Users\arthu\PycharmProjects\Thesis\.venv\Lib\site-packages\xgboost\VERSION', 'xgboost'),

    ],
    hiddenimports=[
        # Basic sklearn imports
        'sklearn',
        'sklearn.base',
        'sklearn.utils',

        # Core algorithm modules
        'sklearn.ensemble',
        'sklearn.tree',
        'sklearn.neighbors',

        # Preprocessing and metrics modules
        'sklearn.preprocessing',
        'sklearn.metrics',
        'sklearn.pipeline',

        # Model calibration and probabilities
        'sklearn.calibration',

        # Common internal modules needed for GBM models
        'sklearn._loss',
        'sklearn._loss.loss',
        'sklearn.utils._weight_vector',

        # Additional ensemble modules often used with GBM
        'sklearn.ensemble._gb',
        'sklearn.ensemble._gb_losses',
        'sklearn.ensemble._gradient_boosting',

        # Common modules for model loading
        'sklearn.feature_selection',
        'sklearn.model_selection',
        'xgboost',
        'xgboost.sklearn',
    ],
    hookspath=[],
    hooksconfig={},
    runtime_hooks=[],
    excludes=[],
    win_no_prefer_redirects=False,
    win_private_assemblies=False,
    cipher=block_cipher,
    noarchive=False,
)

pyz = PYZ(a.pure, a.zipped_data, cipher=block_cipher)

exe = EXE(
    pyz,
    a.scripts,
    [],
    exclude_binaries=True,
    name='ECG_Analysis_Tool',
    debug=False,
    bootloader_ignore_signals=False,
    strip=False,
    upx=True,
    console=True,  # Set to True for debugging
)

coll = COLLECT(
    exe,
    a.binaries,
    a.zipfiles,
    a.datas,
    strip=False,
    upx=True,
    upx_exclude=[],
    name='ECG_Analysis_Tool',
    onefile=True,
)