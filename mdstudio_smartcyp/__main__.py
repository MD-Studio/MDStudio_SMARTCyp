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

# Init basic logging
logging.basicConfig(level=logging.INFO)


if __name__ == '__main__':

    use_wamp = '--rest' not in sys.argv

    if use_wamp:
        main(SmartCypWampApi)

    else:
        app = connexion.App(__name__, specification_dir='./rest/')
        app.add_api('smartcyp_openapi.yaml', arguments={'title': 'MDStudio SMARTCyp REST API'})
        CORS(app.app)
        app.run(port=8081)
