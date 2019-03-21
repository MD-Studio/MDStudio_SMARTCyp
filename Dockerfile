FROM mdstudio/mdstudio_docker3:0.0.3
MAINTAINER Marc van Dijk <m4.van.dijk@vu.nl>

RUN apt-get update && \
	apt-get install -y openjdk-8-jdk && \
	apt-get install -y ant && \
	apt-get clean && \
	rm -rf /var/lib/apt/lists/* && \
	rm -rf /var/cache/oracle-jdk8-installer;

# Fix certificate issues, found as of
# https://bugs.launchpad.net/ubuntu/+source/ca-certificates-java/+bug/983302
RUN apt-get update && \
	apt-get install -y ca-certificates-java && \
	apt-get clean && \
	update-ca-certificates -f && \
	rm -rf /var/lib/apt/lists/* && \
	rm -rf /var/cache/oracle-jdk8-installer;

COPY . /home/mdstudio/mdstudio_smartcyp

RUN chown mdstudio:mdstudio /home/mdstudio/mdstudio_smartcyp

WORKDIR /home/mdstudio/mdstudio_smartcyp

RUN pip install .

USER mdstudio

# Setup JAVA_HOME, this is useful for docker commandline
ENV JAVA_HOME /usr/lib/jvm/java-8-openjdk-amd64/
RUN export JAVA_HOME

CMD ["bash", "entry_point_mdstudio_smartcyp.sh"]
