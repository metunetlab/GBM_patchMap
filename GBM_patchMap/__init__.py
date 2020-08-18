# Add the path inside the sys.path

import sys, os

sys.path.append(os.getcwd()+"/GBM_patchMap/")


# Download the packages which are not downloaded yet/ Import
from paths import path, raw_data_path, curated_data_path, object_path
import needed_packages

