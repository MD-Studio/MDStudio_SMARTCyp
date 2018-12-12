# -*- coding: utf-8 -*-
"""
package:  mdstudio_smartcyp

MDStudio SMARTCyp module
"""

import os

__module__ = 'mdstudio_smartcyp'
__package_path__ = os.path.dirname(os.path.realpath(__file__))
__smartcyp_path__ = os.path.join(__package_path__, 'bin/smartcyp-2.4.2.jar')
__smartcyp_version__ = '2.4.2'
__smartcyp_citation__ = "Rydberg P, Gloriam DE, Zaretzki J, Breneman C, Olsen L. SMARTCyp:" \
                        "A 2D Method for Prediction of Cytochrome P450-Mediated Drug Metabolism." \
                        "ACS Med Chem Lett. 2010;1(3):96-100. Published 2010 Mar 15."
__supported_models__ = ['CYP3A4', 'CYP2C', 'CYP2D6']
