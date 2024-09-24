# -*- mode: python ; coding: utf-8 -*-


a = Analysis(
    ['Main_window.py'],
    pathex=[],
    binaries=[
        ('/home/eliudaguilar/miniconda/yes/envs/myenv/lib/python3.11/site-packages/MDAnalysis/lib/_transformations.cpython-311-x86_64-linux-gnu.so', 'MDAnalysis/lib'), ('/home/eliudaguilar/miniconda/yes/envs/myenv/lib/libxgboost.so', 'lib')
    ],
    datas=[
        ('/home/eliudaguilar/miniconda/yes/envs/myenv/lib/python3.11/site-packages/meeko/data/residue_params.json', 'meeko/data'), 
        ('/home/eliudaguilar/miniconda/yes/envs/myenv/lib/python3.11/site-packages/meeko/data/flexres_templates.json', 'meeko/data'),
        ('/home/eliudaguilar/miniconda/yes/envs/myenv/lib/python3.11/site-packages/meeko/data/atomtypes.txt', 'meeko/data'),
        ('/home/eliudaguilar/miniconda/yes/envs/myenv/lib/python3.11/site-packages/meeko/data/pi_non_aromatic.txt', 'meeko/data'),
        ('/home/eliudaguilar/miniconda/yes/envs/myenv/lib/python3.11/site-packages/meeko/data/waterfield.txt', 'meeko/data'),
        ('/home/eliudaguilar/miniconda/yes/envs/myenv/lib/python3.11/site-packages/xgboost/VERSION', 'xgboost')
    ],
    hiddenimports=[
        'MDAnalysis.lib.formats.cython_util', 
        'MDAnalysis.lib.transformations', 
        'MDAnalysis.lib.formats.libmdaxdr', 
        '_transformations',
        'sklearn',
        'sklearn.utils',
        'sklearn.linear_model',
        'sklearn.ensemble',
        'sklearn.pipeline',
        'xgboost',
    ],
    hookspath=[],
    hooksconfig={},
    runtime_hooks=[],
    excludes=[],
    noarchive=False,
    optimize=0,
)

pyz = PYZ(a.pure)

exe = EXE(
    pyz,
    a.scripts,
    a.binaries,
    a.datas,
    [],
    name='Main_window',
    debug=False,
    bootloader_ignore_signals=False,
    strip=False,
    upx=True,
    upx_exclude=[],
    runtime_tmpdir=None,
    console=True,
    disable_windowed_traceback=False,
    argv_emulation=False,
    target_arch=None,
    codesign_identity=None,
    entitlements_file=None,
)
