#!/bin/bash

export MD_CONFIG_ENVIRONMENTS=dev,docker

# Start REST server
python -u -m mdstudio_smartcyp -a rest -w /tmp/mdstudio/mdstudio_smartcyp -r 1
