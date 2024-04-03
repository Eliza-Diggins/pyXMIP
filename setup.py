from setuptools import setup

setup(
    name="pyXMIP",
    packages=["pyXMIP"],
    version="0.1.0",
    description="",
    author="Eliza C. Diggins",
    author_email="eliza.diggins@utah.edu",
    url="https://github.com/eliza-diggins/pyXMIP",
    download_url="https://github.com/eliza-diggins/pyXMIP/tarball/0.1.0",
    install_requires=["numpy", "scipy", "astropy", "astroquery", "rich"],
    classifiers=[
        "Intended Audience :: Science/Research",
        "Operating System :: OS Independent",
        "Programming Language :: Python :: 3",
        "Topic :: Scientific/Engineering :: Visualization",
    ],
    include_package_data=True,
)
