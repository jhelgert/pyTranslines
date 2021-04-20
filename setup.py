
from setuptools import setup

with open("README.md", "r") as f:
    long_description = f.read()

setup(
    name="pyTranslines",
    url="https://github.com/jhelgert/pyTranslines",
    version='0.0.1',
    packages=['pyTranslines'],
    author='Jonathan Helgert',
    author_email='jhelgert@mail.uni-mannheim.de',
    description='Optimal Control of transmission line networks',
    long_description=long_description,
    long_description_content_type="text/markdown",
    python_requires='>=3.6',
    classifiers=[
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    install_requires=[
        "numpy",
        "gurobipy",
        "matplotlib",
        "pyCIAP"
    ]
)
