#!/bin/bash

export MD_CONFIG_ENVIRONMENTS=dev,docker

# Start REST server in background followed by MDStudio WAMP
python -u -m mdstudio_smartcyp --rest &
python -u -m mdstudio_smartcyp
