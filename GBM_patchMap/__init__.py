# Add the path inside the sys.path

import sys, os

sys.path.append(os.getcwd()+"/venv/GBM_patchMAP/")

from paths import data_path, raw_path, curated_path, object_path, network_path, prism_path, result_path
