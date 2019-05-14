# -*- coding: utf-8 -*-
"""
package:  mdstudio_smartcyp

MDStudio SMARTCyp module
"""

import os
import platform

__module__ = 'mdstudio_smartcyp'
__docformat__ = 'restructuredtext'
__version__ = '{major:d}.{minor:d}'.format(major=1, minor=0)
__author__ = 'Marc van Dijk'
__status__ = 'pre-release beta1'
__date__ = '21 march 2019'
__licence__ = 'Apache Software License 2.0'
__url__ = 'https://github.com/NLeSC/MDStudio'
__copyright__ = "Copyright (c) VU University, Amsterdam"
__package_path__ = os.path.dirname(os.path.realpath(__file__))

# SMARTCyp
__smartcyp_path__ = os.path.join(__package_path__, 'bin/smartcyp-2.4.2.jar')
__smartcyp_version__ = '2.4.2'
__smartcyp_citation__ = "Rydberg P., Gloriam D.E., Zaretzki J., Breneman C., Olsen L. SMARTCyp:" \
                        "A 2D Method for Prediction of Cytochrome P450-Mediated Drug Metabolism." \
                        "ACS Med Chem Lett. 2010;1(3):96-100. Published 2010 Mar 15."
__supported_models__ = ['CYP3A4', 'CYP2C9', 'CYP2D6']

# PLANTS docking
__plants_path__ = os.path.join(__package_path__, 'bin/plants_{0}'.format(platform.system().lower()))
__plants_version__ = '1.1'
__plants_citation__ = "Korb O., Stützle T., Exner T.E., An ant colony optimization approach to flexible" \
                      "protein–ligand docking. Swarm Intelligence (2007); 1(2):115-134."

# SPORES
__spores_path__ = os.path.join(__package_path__, 'bin/spores_{0}'.format(platform.system().lower()))
__spores_version__ = '1.3'
__spores_citation__ = "ten Brink, T., Exner T.E., pK(a) based protonation states and microspecies for" \
                      "protein-ligand docking. J. Comput. Aided Mol. Des. (2010); 24(11):935-42."
