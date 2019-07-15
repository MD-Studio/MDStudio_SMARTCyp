# mdstudio_smartcyp

[![Build Status](https://travis-ci.org/MD-Studio/MDStudio_SMARTCyp.svg?branch=master)](https://travis-ci.org/MD-Studio/MDStudio_SMARTCyp)
[![Codacy Badge](https://api.codacy.com/project/badge/Grade/697c033fd7674ecea28c089150a25dfa)](https://www.codacy.com/app/marcvdijk/MDStudio_SMARTCyp?utm_source=github.com&amp;utm_medium=referral&amp;utm_content=MD-Studio/MDStudio_SMARTCyp&amp;utm_campaign=Badge_Grade)

## Installation Quickstart
MDStudio SMARTCyp can be used in three different ways: in the MDStudio environment, as REST service or as local installation.

### Step 1. Installation
Run MDStudio SMARTCyp as docker container by:

    docker pull mdstudio/mdstudio_smartcyp
    docker run --name mdstudio_smartcyp -p 8081:8081 <container ID>

Or download/clone the MDStudio_SMARTCyp repository and install locally using pip.
The repository comes equipped with the [SMARTCyp](https://smartcyp.sund.ku.dk/mol_to_som) software but for using 
[PLANTS](https://uni-tuebingen.de/de/37876) docking and [SPORES](https://uni-tuebingen.de/fakultaeten/mathematisch-naturwissenschaftliche-fakultaet/fachbereiche/pharmazie-und-biochemie/pharmazie/pharmazeutische-chemie/pd-dr-t-exner/research/spores/) 
structure preparation functionality you will need to obtain and copy the executables to the mdstudio_smartcyp/bin/ 
directory yourself. Use the following naming scheme for the executables: *\<plants or spores>_\<platform>* where platform 
is either linux or darwin for OSX. Then install using pip as:

    pip install MDStudio_SMARTCyp/

### Step 2a. Start as a MDStudio microservice
MDStudio SMARTCyp can be used as a fully compliant microservice within a [MDStudio](https://github.com/MD-Studio/MDStudio)
environment. After installation you can launch the package as standalone MDStudio microservice:

    python -u -m mdstudio_smartcyp

or as part of an MDStudio environment orchestrator either in standalone mode or as docker container. In both of the 
latter two cases you will need to edit the *entry_point_mdstudio_smartcyp.sh* script and set the -a/--api_mode to *wamp*
to activate the MDStudio WAMP interface (Web Application Messaging Protocol).

### Step 2b. Start as a REST service
MDStudio SMARTCyp exposes an [OpenAPI](https://swagger.io/specification/) (swagger) compliant REST API using a production
grade [Flask](https://palletsprojects.com/p/flask/) web server (using gevent coroutines). The docker MDStudio SMARTCyp
docker image (step 1) is by default configured to start the REST instead of the WAMP server. In standalone mode the REST
server is started as:

    python -u -m mdstudio_smartcyp -a rest

### Step 3a. Using the MDStudio microservice


### Step 3b. Using the REST service
The functional REST service endpoints are defined using the OpenAPI/Swagger defintions. The service exposes a Swagger UI
web page that lists the endpoints the service exposes, provides an overview of the configuration parameters and allow you
to use the endpoints directly from the ui. Access the ui as:

    http(s)://<base URL, localhost by default>:8081/ui/

### Basic configuration options
The mdstudio_smartcyp python process accepts the following command line configuration options:

- -a/--api_mode: *rest* or *wamp* to start the service in the REST or WAMP mode respectively.
- -w/--base_work_dir: the base directory where SMARTCyp, PLANTS or SPORES work directories will be stored. The systems temporary (/tmp) directory will be used by default.
- -r/--result_storage_time: how many hours the calculated results will remain available before cleanup. 0 by default which means no cleanup.
- -p/--http_port: the network port the REST or WAMP service will be started on. 8081 by default.
