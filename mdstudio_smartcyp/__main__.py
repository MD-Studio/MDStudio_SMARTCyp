# -*- coding: utf-8 -*-

import os
import sys
import logging
import connexion

from flask_cors import CORS
from mdstudio.runner import main

# Try import package
try:
    from mdstudio_smartcyp.wamp_services import SmartCypWampApi

except ImportError:

    # Add modules in package to path
    modulepath = os.path.abspath(os.path.join(os.path.dirname(__file__), '../'))
    if modulepath not in sys.path:
        sys.path.insert(0, modulepath)

    from mdstudio_smartcyp.wamp_services import SmartCypWampApi

from mdstudio_smartcyp import __module__, __package_path__

# Init basic logging
logging.basicConfig(level=logging.INFO)


if __name__ == '__main__':

    # Start the REST (--rest) or WAMP (--wamp) interface, --rest by default
    use_wamp = '--rest' not in sys.argv

    if use_wamp:
        logging.debug('Start {0} WAMP interface at {1}'.format(__module__, __package_path__))

        main(SmartCypWampApi)

    else:
        logging.debug('Start {0} REST interface at {1}'.format(__module__, __package_path__))

        app = connexion.App(__name__, specification_dir='./rest/')
        app.add_api('smartcyp_openapi.yaml', arguments={'title': 'MDStudio SMARTCyp REST API'})
        CORS(app.app)
        app.run(port=8081)
