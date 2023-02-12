# The base ubuntu LTS image, 22.04 had fortran compilation issues
FROM ubuntu:20.04

# Install Ubuntu package dependencies, needed to add an explicit timezone
RUN apt-get update && \
    DEBIAN_FRONTEND=noninteractive \
    TZ=America \
    apt-get install -y \
    python3-numpy \
    python3-dev \
    python3-scipy \
    python3-nose \
    python3-pip \
    python3-matplotlib \
    gfortran \
    make \
    gmsh \
    libatlas-base-dev \
    libblas-dev \
    liblapack-dev \
    libsuitesparse-dev \
    nano \
    ssh

# install jupyter and jupyterlab, oh my zsh and then cleanup
RUN pip3 install jupyter -U && pip3 install jupyterlab
RUN sh -c "$(wget -O- https://github.com/deluan/zsh-in-docker/releases/download/v1.1.2/zsh-in-docker.sh)" --     -t mh     -p https://github.com/zsh-users/zsh-autosuggestions     -p https://github.com/zsh-users/zsh-completions     -p https://github.com/zsh-users/zsh-syntax-highlighting
RUN chsh -s /usr/bin/zsh &&  apt-get clean && rm -rf /var/lib/apt/lists/*

# Add the NumBAT source code -> Will be overwritten when used with mounted volumes! (which is good)
COPY ./ /home/EMUstack/

# Compile the Fortran code, only use when running tests or copying compiled source to host
WORKDIR /home/EMUstack/backend/fortran/
RUN make

# Add the backend files to the python path, sets shell so jupyterlab uses zsh and dark mode
ENV PYTHONPATH "${PYTHONPATH}:/home/EMUstack/backend/"
ENV SHELL "/usr/bin/zsh"
COPY overrides.json /usr/local/share/jupyter/lab/settings/

# Run the tests (when desired), Fails 1 test simple_mk_msh with a tolerance issue with ubuntu 22.04
WORKDIR /home/EMUstack/tests/
# RUN nosetests3

# Change the working directory to final spot
WORKDIR /home/EMUstack/

# Finish with jupyterlab by default, can always use "/bin/zsh" to start in shell
EXPOSE 8888
CMD "jupyter-lab" "--ip='0.0.0.0'" "--port=8888" "--no-browser" "--allow-root"