import sys, os, importlib

for package in ["scipy", "pandas", "numpy", "pickle", "networkx",
                "matplotlib", "seaborn", "warnings"]:
    try:
        if package not in ["warnings", "pickle"]:
            globals()[package] = importlib.import_module(package)
    except ImportError:
        print("No %s package was installed" %package)
        sys.exit()
