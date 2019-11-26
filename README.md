# mdstudio_smartcyp

[![Build Status](https://travis-ci.org/MD-Studio/MDStudio_SMARTCyp.svg?branch=master)](https://travis-ci.org/MD-Studio/MDStudio_SMARTCyp)
[![Codacy Badge](https://api.codacy.com/project/badge/Grade/697c033fd7674ecea28c089150a25dfa)](https://www.codacy.com/app/marcvdijk/MDStudio_SMARTCyp?utm_source=github.com&amp;utm_medium=referral&amp;utm_content=MD-Studio/MDStudio_SMARTCyp&amp;utm_campaign=Badge_Grade)
[![codecov](https://codecov.io/gh/MD-Studio/MDStudio_SMARTCyp/branch/master/graph/badge.svg)](https://codecov.io/gh/MD-Studio/MDStudio_SMARTCyp)
[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/MD-Studio/MDStudio_SMARTCyp/master?filepath=examples)

![Configuration settings](mdstudio-logo.png)

The MDStudio SMARTCyp service combines SMARTCyp reactivity-based Cytochrome P450 site-of-metabolism prediction and
PLANTS structure-based ligand interaction prediction into one service.

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
directory yourself. Use the following naming scheme for the executables: *\<plants or spores>_\<platform>* where 
platform is either linux or darwin for OSX. Then install using pip as:

    pip install MDStudio_SMARTCyp/

### Step 2a. Launch as MDStudio microservice
MDStudio SMARTCyp can be used as a fully compliant microservice within a [MDStudio](https://github.com/MD-Studio/MDStudio)
environment. After installation you can launch the package as standalone MDStudio microservice:

    python -u -m mdstudio_smartcyp

or as part of an MDStudio environment orchestrator either in standalone mode or as docker container. In both of the 
latter two cases you will need to edit the *entry_point_mdstudio_smartcyp.sh* script and set the -a/--api_mode to *wamp*
to activate the MDStudio WAMP interface (Web Application Messaging Protocol).

### Step 2b. Launch as REST service
MDStudio SMARTCyp exposes an [OpenAPI](https://swagger.io/specification/) (swagger) compliant REST API using a production
grade [Flask](https://palletsprojects.com/p/flask/) web server (using gevent coroutines). The MDStudio SMARTCyp
docker image (step 1) is by default configured to start as REST service instead of the WAMP service. In standalone mode 
the REST service is started as:

    python -u -m mdstudio_smartcyp -a rest

MDStudio SMARTCyp uses a [OpenRiskNet](https://openrisknet.org) (ORN) compliant API to ensure that each service is 
accessible to the ORN interoperability layer. This includes service level annotation according to the ORN semantic 
interoperability layer concept. Consult the *mdstudio_smartcyp/rest/smartcyp_openapi* for details.

### Step 3a. Using the MDStudio microservice
The MDStudio environment allows you to interact with the endpoints of a microservice in many different ways from a high
level workflow manager, use in different programming languages to command line access. Please consult the MDStudio
documentation for more information on usage.

### Step 3b. Using the REST service
The functional REST service endpoints are defined using the OpenAPI/Swagger definitions. The service exposes a Swagger UI
web page that lists the endpoints the service exposes, provides an overview of the configuration parameters and allow you
to use the endpoints directly from the ui. Access the ui as:

    http(s)://<base URL, localhost by default>:<service port, 8081 by default>/ui/

A few examples using *curl* to interact with the endpoints:
 
+ Run SMARTCyp using a SMILES string: curl -F smiles='C1=CC=C(C(=C1)CC(=O)O)NC2=C(C=CC=C2Cl)Cl' 'http://localhost:8081/smartcyp'
+ Run PLANTS: curl -F ligand_file='@/path/to/ligand.mol2' -F protein_file='@/path/to/protein.mol2' -F bindingsite_center=-0.989 -F bindingsite_center=3.261 -F bindingsite_center=0.826 'http://localhost:8081/plants_docking'
+ Get two of the PLANTS docking poses from the previous calculation as multi MOL2 file: curl -F paths=docking-9y4f3anb/diclofenac_entry_00001_conf_50.mol2 -F paths=docking-9y4f3anb/diclofenac_entry_00001_conf_49.mol2  'http://localhost:8081/plants_docking_structures'
+ Reprotonate a structure using SPORES: curl -F mol='@/path/to/ligand.mol2' -F input_format=mol2 -F spores_mode=reprot  'http://localhost:8081/spores'

### Basic configuration options
The mdstudio_smartcyp python process accepts the following command line configuration options:

+ -a/--api_mode: *rest* or *wamp* to start the service in the REST or WAMP mode respectively.
+ -w/--base_work_dir: the base directory where SMARTCyp, PLANTS or SPORES work directories will be stored. The systems temporary (/tmp) directory will be used by default.
+ -r/--result_storage_time: how many hours the calculated results will remain available before cleanup. 0 by default which means no cleanup.
+ -p/--http_port: the network port the REST or WAMP service will be started on. 8081 by default.
