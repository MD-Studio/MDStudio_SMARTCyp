#!/bin/bash

export MD_CONFIG_ENVIRONMENTS=dev,docker

# Start REST server
python -u -m mdstudio_smartcyp -a rest -w /tmp -r 1
