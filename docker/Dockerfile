FROM mdstudio/mdstudio_docker3:0.0.3

RUN apt-get update && \
	apt-get install -y --no-install-recommends openjdk-8-jdk && \
	apt-get install -y --no-install-recommends ant && \
	apt-get clean && \
	rm -rf /var/lib/apt/lists/* && \
	rm -rf /var/cache/oracle-jdk8-installer;

# Fix certificate issues, found as of
# https://bugs.launchpad.net/ubuntu/+source/ca-certificates-java/+bug/983302
RUN apt-get update && \
	apt-get install -y --no-install-recommends ca-certificates-java && \
	apt-get clean && \
	update-ca-certificates -f && \
	rm -rf /var/lib/apt/lists/* && \
	rm -rf /var/cache/oracle-jdk8-installer;

# Setup JAVA_HOME, this is useful for docker commandline
ENV JAVA_HOME /usr/lib/jvm/java-8-openjdk-amd64/
RUN export JAVA_HOME

# Install package
COPY . /home/mdstudio
RUN chown -R mdstudio:mdstudio /home/mdstudio
RUN chmod -R 755 /home/mdstudio
WORKDIR /home/mdstudio
RUN pip install .
USER mdstudio

# Set entrypoint and start process
CMD ["bash", "entry_point_mdstudio_smartcyp.sh"]
