# -*- coding: utf-8 -*-

import os
import sys
import logging
import argparse
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

from mdstudio_smartcyp import __module__, __package_path__, __author__, __date__, __copyright__
from mdstudio_smartcyp.utils import PeriodicCleanup

# Init basic logging
logging.basicConfig(level=logging.INFO)


if __name__ == '__main__':

    usage = """
    MDStudio SAMRTCyp service
    
    Start the SMARTCyp service as either a REST or WAMP service using the -a/--api_mode argument.
    
    Result cleanup:
    PLANTS docking results are stored on disc to allow the user to query and manipulate the results
    of a docking run in independent and successive calls to service REST/WAMP endpoints.
    The results will need to be cleaned-up at some point in time to prevent storage from filling up.
    By default, results are stored in the systems temporary directory and as such they should be 
    cleaned by the system at regular intervals. However, these directories are not always consistent
    preventing the service from accessing previous results.
    Therefor, use the -w/--base_work_dir argument to define the base storage directory and
    (optionally) a maximum time in hours results are saved before cleanup.
    
    {0}:
    {1}, {2}, {3}
    """.format(__module__, __date__, __author__, __copyright__)

    # Setup command line option parser
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter, usage=usage)
    parser.add_argument('-a', '--api_mode',
                        help='Start {0} with WAMP or REST API'.format(__module__),
                        choices=['rest', 'wamp'], default='wamp')
    parser.add_argument('-w', '--base_work_dir',
                        help='Directory to store results. System temp directory by default')
    parser.add_argument('-r', '--result_storage_time',
                        help='Maximum time (hours) results are saved before cleanup. 0 will not clean',
                        type=int, default=0)
    parser.add_argument('-p', '--http_port',
                        help='HTTP network port the service connects to',
                        type=int, default=8081)
    args = parser.parse_args()

    if args.base_work_dir:

        # Set dynamic base_work_dir as environmental variable accessible from else where
        os.environ['BASE_WORK_DIR'] = args.base_work_dir

        # Start results cleanup event if 'base_work_dir' and 'result_storage_time' are set
        if args.result_storage_time > 0:
            p = PeriodicCleanup(args.base_work_dir, args.result_storage_time)
            p.start()

    # Start service REST or WAMP API
    if args.api_mode == 'wamp':
        logging.debug('Start {0} WAMP interface at {1}'.format(__module__, __package_path__))

        main(SmartCypWampApi)

    else:
        logging.debug('Start {0} REST interface at {1}'.format(__module__, __package_path__))

        app = connexion.App(__name__, specification_dir='./rest/')
        app.add_api('smartcyp_openapi.yaml', arguments={'title': 'MDStudio SMARTCyp REST API'})
        CORS(app.app)
        app.run(port=args.http_port)
