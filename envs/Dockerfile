FROM debian:11
# FROM rocker/tidyverse:4.2.0  # alternative


############################ Setup: OS #########################################
# Set up locales properly
RUN apt-get update && \
    apt-get install --yes --no-install-recommends locales && \
    apt-get purge && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

RUN echo "en_US.UTF-8 UTF-8" > /etc/locale.gen && \
    locale-gen

ENV LC_ALL en_US.UTF-8
ENV LANG en_US.UTF-8
ENV LANGUAGE en_US.UTF-8

# Use bash as default shell, rather than sh
ENV SHELL /bin/bash

# run critical system updates ... one could also use the artful-curl pack
RUN apt-get update && \
    apt-get install --yes --no-install-recommends \
        adduser \
        wget \
        bzip2 \
        ca-certificates \
        curl \
        git \
        zip \
        unzip \
        procps && \
    apt-get install --yes \
        libpng-dev \
        libcurl4-gnutls-dev \
        libssl-dev \
        libxml2-dev \
        libgit2-dev \
        zlib1g-dev \
        build-essential \
        libgeos-dev \
        libbz2-dev \
        liblzma-dev && \
    apt-get purge && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

# set up a user, to avoid running containers as root
ENV NB_USER container_user
ENV HOME /home/${NB_USER}
RUN adduser --disabled-password \
    --gecos "Default user" \
    ${NB_USER}

# set up the user directory structure
RUN mkdir -p ${HOME}/tmp && \
    mkdir -p ${HOME}/bin && \
    mkdir -p ${HOME}/analysis && \
    mkdir -p ${HOME}/.cache && \
    mkdir -p ${HOME}/.cache/matplotlib
ENV PATH ${HOME}/bin:$PATH
################################################################################


############################ Setup: conda ######################################
# install conda
USER root
ENV CONDA_FILE https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh
RUN echo "export PATH=${HOME}/conda/bin:$PATH" > /etc/profile.d/conda.sh && \
    wget --quiet ${CONDA_FILE} -O /tmp/miniforge.sh && \
    /bin/bash /tmp/miniforge.sh -b -p ${HOME}/conda && \
    rm /tmp/miniforge.sh && \
    ${HOME}/conda/bin/conda clean -atipy
ENV PATH ${HOME}/conda/bin:$PATH
RUN conda init

# alternative conda install (instead of the manual install from above)
# FROM continuumio/miniconda3:latest
################################################################################


############################ Setup: non-conda tools ############################
# download things to tmp
WORKDIR ${HOME}/tmp

# Example: install GCTA to bin
#RUN wget --quiet \
#    https://cnsgenomics.com/software/gcta/gcta_1.92.1beta5.zip && \
#    unzip gcta_1.92.1beta5.zip && \
#    chmod 755 gcta_1.92.1beta5/gcta64 && \
#    mv gcta_1.92.1beta5/gcta64 /bin/gcta64 && \
#    rm gcta_1.92.1beta5.zip && \
#    rm -r gcta_1.92.1beta5
################################################################################


############################ Setup: conda environment ##########################
# make sure conda files are owned by user
RUN chown -R ${NB_USER}:${NB_USER} ${HOME}
USER ${NB_USER}

# download the conda environment OR copy from local directory
COPY environment.yml ${HOME}/tmp/environment.yml
COPY environment_setup.sh ${HOME}/tmp/environment_setup.sh

# install conda environment
# NOTE: conda clean -atipy
#   removes everything w/o confirmation (leaves all environments functional).
#   for gmmat -- run script 
RUN conda config --set always_yes yes --set changeps1 no && \
   mamba update mamba --yes && \
#    mamba env create -v --file ${HOME}/tmp/environment.yml && \
   bash ${HOME}/tmp/environment_setup.sh \
   mamba clean -atipy && \
   rm ${HOME}/tmp/environment*

# export the conda path to force container_env since .bashrc necessarily
# sourced when when the image is used
# ENV PATH ${HOME}/conda/envs/container_env/bin:$PATH
# ENV CONDA_DEFAULT_ENV container_env
# ENV CONDA_PREFIX ${HOME}/conda/envs/container_env

# Set up matplotlib requirements
ENV MPLCONFIGDIR ${HOME}/.cache/matplotlib
################################################################################


############################ Setup: R packages #################################
# Set the R libs
# NOTE: As of 2022-05-27, none of the R libraries are set in the Docker image.
#       R raised errors when doing
#       ENV R_LIBS /usr/local/lib/R/library:/usr/lib/R/library:${R_LIBS}
#       ...which would preserve any values of ${R_LIBS} in the future
#       The error is raised when ${R_LIBS} = "" and trailing : is at the end
#       of the path. An if, then else bash clause did not solve the issue
#       because Docker RUN export VAR do not track with the image.
#       The easiest workaround for all of this was to hard code these ENV vars.
#       WARNING: This may cause issues in the future.
ENV R_LIBS /usr/local/lib/R/library:/usr/lib/R/library
ENV R_LIBS_SITE /usr/local/lib/R/site-library/:/usr/lib/R/site-library
ENV R_LIBS_USER ${HOME}/R/x86_64-pc-linux-gnu-library/4.2

# NOTE: Renv expects us to be in the ${HOME}
WORKDIR ${HOME}

# Need to be root for library install in main R library
# USER root
# RUN R --vanilla -e "install.packages('renv', repos = c(CRAN = 'https://cloud.r-project.org'))"

# Now continue to build the renv as the user
#USER ${NB_USER}

# Copy renv information from local directory
# COPY renv.lock ${HOME}/tmp/renv.lock

# Install the packages
# NOTE: For Rhtslib, need to run `apt-get install libbz2-dev liblzma-dev`
# NOTE: For rgeos (SeuratObject), need to run `apt-get install libgeos-dev`
# NOTE: The below commands give access to default libraries and configs
#       in the rocker image. Below they are set in R, but we set them in
#       environment variables.
#       options(renv.config.user.profile = T)
#       options(renv.config.user.library = T)
#       options(renv.config.external.libraries = .libPaths())
# ENV RENV_CONFIG_USER_PROFILE TRUE
# ENV RENV_CONFIG_USER_LIBRARY TRUE
# ENV RENV_CONFIG_EXTERNAL_LIBRARIES ${R_LIBS_USER}:${R_LIBS_SITE}:${R_LIBS}
# ENV RENV_CONSENT TRUE
#RUN echo "R_LIBS ${R_LIBS}" && \
#    echo "R_LIBS_SITE ${R_LIBS_SITE}" && \
#    echo "R_LIBS_USER ${R_LIBS_USER}" && \
#    echo "RENV_CONFIG_EXTERNAL_LIBRARIES ${RENV_CONFIG_EXTERNAL_LIBRARIES}"
# NOTE: Dropped the below command because not necessary
#       echo "renv::init(bare = F, force = T)" >> ${HOME}/tmp/run.R && \
# NOTE: Dropped the below command because with R_LIBS, R_LIBS_SITE, and
#       R_LIBS_USER set, the renv does not need to be activated. R knows it
#       automatically.
#       echo "renv::activate()" >> ${HOME}/tmp/run.R
# RUN echo "library(renv)" > ${HOME}/tmp/run.R && \
#     echo "options(renv.consent = T)" >> ${HOME}/tmp/run.R && \
#     echo "renv::restore(lockfile=\"${HOME}/tmp/renv.lock\", prompt = F)" >> ${HOME}/tmp/run.R

# Install the packages
#RUN cat ${HOME}/tmp/run.R
# RUN Rscript ${HOME}/tmp/run.R
#RUN cat .Rprofile

# Clear out files no longer required
# RUN rm -fr ${HOME}/.cache/R/renv/source  # Clear out renv tarballs
#RUN R -e "library(optparse)"  # Test the R install
################################################################################


############################ Setup: cleanup ####################################
# clear tmp if there is anything in there...
USER root
RUN rm -rf ${HOME}/tmp && rm -rf /tmp/*
################################################################################


############################ Deployment ########################################
# set wd to user home
WORKDIR ${HOME}

# set bash as the default entry point
# NOTE: for Google Pipelines APIv2, docker does not execute if entrypoint is set
#ENTRYPOINT [ "/bin/bash", "-c" ]

# define default command as an alternative to ENTRYPOINT
# NOTE: for Google Pipelines APIv2, docker does not execute if CMD is set
#CMD ["/bin/bash"]

# set the user
# NOTE: for Snakemake + Kubernetes we get an error if this != root
# NOTE: for adding a monitoring script with cromwell, the monitoring script is
#       copied into the docker image and needs to be chmod to run. Permission
#       errors arise if not root.
USER ${NB_USER}
################################################################################
