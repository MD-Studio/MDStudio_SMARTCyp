FROM java:openjdk-8-jre
MAINTAINER Marc van Dijk <m4.van.dijk@vu.nl>

# download from http://www.farma.ku.dk/smartcyp/download.php
COPY mdstudio_smartcyp/bin/smartcyp-2.4.2.jar smartcyp.jar