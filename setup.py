import setuptools
from pygen_structures.version import version

with open("README.md", 'r') as readme_file:
    long_description = readme_file.read()

setuptools.setup(
    name="pygen-structures",
    version=version,
    author="Travis Hesketh",
    author_email="travis@hesketh.scot",
    description="3D molecular structure generation for MD simulation",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/thesketh/pygen-structures",
    packages=setuptools.find_packages(),
    entry_points = {
        "console_scripts": ['pygen-structures = pygen_structures.__main__:main']
    },
    classifiers=[
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.7",
        "Development Status :: 4 - Beta",
        "Environment :: Console",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Chemistry",
        "License :: OSI Approved :: BSD License",
        "Operating System :: OS Independent",
        "Typing :: Typed"
    ],
    python_requires='>=3.6',
    include_package_data=True,
)
