FROM rocker/tidyverse:4.4

ARG GITHUB_USERNAME
ARG GITHUB_TOKEN

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y git build-essential ant wget unzip perl unzip emacs maven && apt-get clean && apt-get purge && rm -rf /var/lib/apt/lists/* /tmp/*

WORKDIR /gusApp
WORKDIR /gusApp/gus_home
WORKDIR /gusApp/project_home

ENV GUS_HOME=/gusApp/gus_home
ENV PROJECT_HOME=/gusApp/project_home
ENV PATH=$PROJECT_HOME/install/bin:$PATH
ENV PATH=$GUS_HOME/bin:$PATH

# Create the Maven settings.xml
RUN mkdir -p /root/.m2 && \
    echo '<settings xmlns="http://maven.apache.org/SETTINGS/1.0.0">' > /root/.m2/settings.xml && \
    echo '  <servers>' >> /root/.m2/settings.xml && \
    echo '    <server>' >> /root/.m2/settings.xml && \
    echo '      <id>veupathdb</id>' >> /root/.m2/settings.xml && \
    echo "      <username>${GITHUB_USERNAME}</username>" >> /root/.m2/settings.xml && \
    echo "      <password>${GITHUB_TOKEN}</password>" >> /root/.m2/settings.xml && \
    echo '    </server>' >> /root/.m2/settings.xml && \
    echo '  </servers>' >> /root/.m2/settings.xml && \
    echo '</settings>' >> /root/.m2/settings.xml

RUN export INSTALL_GIT_COMMIT_SHA=306fc4d29661fbd881211d92a1698bdbb2d23b2d \
    && git clone https://github.com/VEuPathDB/install.git \
    && cd install \
    && git checkout $INSTALL_GIT_COMMIT_SHA

RUN mkdir -p $GUS_HOME/config && cp $PROJECT_HOME/install/config/gus.config.sample $GUS_HOME/config/gus.config

RUN export CBIL_GIT_COMMIT_SHA=c1d5ddbd34c7525094aa31dbc0b7b50cf2ba826f \
    && git clone https://github.com/VEuPathDB/CBIL.git \
    && cd CBIL \
    && git checkout $CBIL_GIT_COMMIT_SHA \
    && bld CBIL

RUN export GUS_GIT_COMMIT_SHA=cf9a99dba00bea3875f9eb5128294ed4a7a25377 \
    && git clone https://github.com/VEuPathDB/GusAppFramework.git \
    && mv GusAppFramework GUS \
    && cd GUS \
    && git checkout $GUS_GIT_COMMIT_SHA \
    && bld GUS/Supported \
    && bld GUS/Community

RUN export FGP_GIT_COMMIT_SHA=e5fbefafc4aa368b71236e8d0db8c9880c1b6e2f \
    && git clone https://github.com/VEuPathDB/FgpUtil.git \
    && cd FgpUtil \
    && git checkout $FGP_GIT_COMMIT_SHA \
    && bld FgpUtil

RUN export TM_GIT_COMMIT_SHA=4e66af9a7fd877d39c79886d06d80e7a9b5bbbec \
    && git clone https://github.com/VEuPathDB/TuningManager.git \
    && cd TuningManager \
    && git checkout $TM_GIT_COMMIT_SHA \
    && bld TuningManager

RUN export DOTS_GIT_COMMIT_SHA=82776f13e4b9d827a05e59038b1a03c2de9e0caa \
    && git clone https://github.com/VEuPathDB/DoTS.git \
    && cd DoTS \
    && git checkout $DOTS_GIT_COMMIT_SHA
#&& bld DoTS

RUN export APICOMMONDATA_GIT_COMMIT_SHA=d1036d6d33972f0265c88d1cd7602d018b54e540 \
    && git clone https://github.com/VEuPathDB/ApiCommonData.git \
    && cd ApiCommonData \
    && git checkout $APICOMMONDATA_GIT_COMMIT_SHA \
    && bld ApiCommonData/Load

WORKDIR /work