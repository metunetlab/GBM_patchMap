import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="GBM_patchMap",
    version="0.0.1",
    author="Cansu Dincer, Nurcan Tuncbag",
    author_email="cansu.dincer@metu.edu.tr",
    description="GBM data implementation of patchMap",
    long_description="patchMap” is a Python tool mapping molecular aberrations on three dimensional (3D) structures of corresponding proteins and protein interactions, finding spatial organisation of these mutations on the structures (mutation patches), and using network-based approach for the visualisation to highlight the rewiring in the patient networks",
    long_description_content_type="text/markdown",
    url="https://github.com/CansuDincer/GBM_patchMap/",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GNU General Public License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
)
