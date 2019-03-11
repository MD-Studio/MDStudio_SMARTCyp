import logging
logging.basicConfig(level=logging.INFO)

from graphit.graph_axis.graph_axis_mixin import NodeAxisTools
from graphit.graph_io.io_json_format import read_json
from graphit.graph_model_classes.model_openapi import OpenAPIORM, OpenAPI


# Define OpenAPI model (version 2.0)
oapi = OpenAPI()
oapi.orm = OpenAPIORM
oapi.node_tools = NodeAxisTools

# Import openapi in JSON
api = read_json('smartcyp_openapi.json', graph=oapi, level=0)

# Get basic info api
root = api.get_root()
print(root.swagger())
for obj in root.info:
    for n in obj.items():
        print('{0}:  {1}'.format(*n))

# List all api endpoints
print('\nEndpoints:')
for endpoint in api.list_methods():
    print(endpoint.name())

# Get specific endpoint
print('')
endp = api.get_method('predict')
endp.info()

#print(endp.call(organism='Escherichia coli'))