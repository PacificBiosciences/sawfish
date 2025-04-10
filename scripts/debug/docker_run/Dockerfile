FROM ubuntu:22.04

RUN apt-get -qq update && \
    apt-get install -y -q \
    curl

ENV SAWFISH_HOME=/opt/sawfish \
    SAWFISH_VERSION=v0.12.10

ENV PATH=${SAWFISH_HOME}/bin:$PATH

RUN mkdir -p ${SAWFISH_HOME} && \
    curl -sL https://github.com/PacificBiosciences/sawfish/releases/download/${SAWFISH_VERSION}/sawfish-${SAWFISH_VERSION}-x86_64-unknown-linux-gnu.tar.gz | tar --strip-components=1 -xz -C $SAWFISH_HOME

