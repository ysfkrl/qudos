from setuptools import setup, find_packages

setup(
    name="qudos",
    version="0.1.0",  # Update this version as needed
    author="Yusuf Karli",
    author_email="ysfkarli@gmail.com",
    description="A Python package for quantum optics simulations with laser pulses.",
    long_description=open("README.md").read(),
    long_description_content_type="text/markdown",
    url="https://github.com/ysfkrl/qudos",  # Replace with your actual repo if applicable
    packages=find_packages(),
    install_requires=[
        "numpy",
        "qutip",
        "matplotlib",
        "plotly"
    ],
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Physics",
    ],
    python_requires=">=3.6",
)
