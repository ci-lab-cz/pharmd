import setuptools
from os import path
import pharmd

this_directory = path.abspath(path.dirname(__file__))
with open(path.join(this_directory, 'README.md'), encoding='utf-8') as f:
    long_description = f.read()

setuptools.setup(
    name="pharmd",
    version=pharmd.__version__,
    author="Pavel Polishchuk",
    author_email="pavel_polishchuk@ukr.net",
    description="PharMD: MD pharmacophores and virtual screening",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/ci-lab-cz/pharmd.git",
    packages=['pharmd'],
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: BSD License",
        "Operating System :: OS Independent",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Topic :: Scientific/Engineering :: Chemistry"
    ],
    python_requires='>=3.6',
    install_requires=['pmapper>=0.3', 'psearch>=0.0.2', 'mdtraj>=1.9.3'],
    entry_points={'console_scripts':
                      ['get_distinct = pharmd.get_distinct:entry_point',
                       'md2pharm = pharmd.md2pharm:entry_point',
                       'get_scores = pharmd.get_scores:entry_point']},
    include_package_data=True
)
