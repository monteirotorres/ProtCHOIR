from setuptools import setup
from codecs import open
from os import path

here = path.abspath(path.dirname(__file__))

with open(path.join(here, "README.rst"), encoding="utf-8") as f:
    long_description = f.read()

setup(
    name="ProtCHOIR",
    packages=["ProtCHOIR"],
    version="1.2.19",
    description="A Modeller-based pipeline to generate homo-oligomers.",
    long_description=long_description,
    license="CC-BY-4.0",
    author="Pedro Torres, Tom Blundell",
    author_email="monteirotorres@biof.ufrj.br",
    url="https://github.com/monteirotorres/ProtCHOIR",
    include_package_data=True,
    package_data={"ProtCHOIR": ["Contents/*.html", "Contents/*.svg"]},
    keywords=[
        "PDB",
        "structure",
        "protein",
        "oligomerization",
        "oligomers",
        "multimer",
    ],
    classifiers=[],
    py_modules=[
        "ProtCHOIR.Toolbox",
        "ProtCHOIR.__main__",
        "ProtCHOIR.UpdateDatabases",
        "ProtCHOIR.MakeOligomer",
        "ProtCHOIR.AnalyseOligomer",
        "ProtCHOIR.AnalyseProtomer",
        "ProtCHOIR.Initialise",
    ],
    entry_points={"console_scripts": ["ProtCHOIR = ProtCHOIR.__main__:main"]},
    install_requires=[
        "pandas",
        "progressbar2",
        "biopython",
        "parasail",
        "matplotlib",
        "numpy",
        "jinja2",
        "networkx",
    ],
)
