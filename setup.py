#! /usr/bin/env python
# -*- coding: utf-8 -*-

# package: mdstudio_smartcyp
# file: setup.py
#
# Part of MDStudio_SMARTCyp providing Site of Metabolism prediction for
# Cytochrome P450s to the MDStudio environment using the SMARTCyp package
# developed by P. Rydberg et al:
#
#   Rydberg P, Gloriam DE, Zaretzki J, Breneman C, Olsen L. SMARTCyp:
#   A 2D Method for Prediction of Cytochrome P450-Mediated Drug Metabolism.
#   ACS Med Chem Lett. 2010;1(3):96-100. Published 2010 Mar 15.
#   doi:10.1021/ml100016x
#
# Copyright © 2016 Marc van Dijk, VU University Amsterdam, the Netherlands
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

from setuptools import setup, find_packages

distribution_name = 'mdstudio_smartcyp'

setup(
    name=distribution_name,
    version=1.0,
    description='MDStudio SMARTCyp module',
    author="""
    Marc van Dijk - VU University - Amsterdam
    Paul Visscher - Zefiros Software (www.zefiros.eu)
    Felipe Zapata - eScience Center (https://www.esciencecenter.nl/)""",
    author_email=['m4.van.dijk@vu.nl', 'f.zapata@esciencecenter.nl'],
    url='https://github.com/MD-Studio/MDStudio_SMARTCyp',
    license='Apache Software License 2.0',
    keywords='MDStudio structures cheminformatics',
    platforms=['Any'],
    packages=find_packages(),
    package_data={distribution_name: ['schemas/*', 'schemas/endpoints/*', 'data/*', 'bin/*']},
    py_modules=[distribution_name],
    install_requires=['flask', 'flask-cors', 'connexion', 'swagger-ui-bundle', 'gevent',
                      'mdstudio', 'matplotlib', 'scipy', 'mdinteract'],
    extras_require={'test': ['requests']},
    dependency_links=["https://github.com/cinfony/cinfony/tarball/master#egg=cinfony-1.2"],
    include_package_data=True,
    zip_safe=True,
    classifiers=[
        'Development Status :: 3 - Alpha',
        'License :: OSI Approved :: Apache Software License',
        'Programming Language :: Python',
        'Topic :: Scientific/Engineering :: Chemistry',
        'Operating System :: OS Independent',
        'Intended Audience :: Science/Research',
    ]
)
