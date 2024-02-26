# Latest fedora image. It works well with the new compilation flags.
FROM fedora:latest

# Install Fedora package dependencies, needed to add an explicit timezone
RUN dnf -y update
RUN dnf -y install make gmsh python3-numpy python3-devel python3-scipy python3-nose python3-pip python3-matplotlib gfortran
RUN dnf -y install suitesparse-devel blas-devel lapack-devel atlas-devel
RUN dnf -y install nano micro openssh
RUN dnf -y install jupyterlab

# Add the NumBAT source code -> Will be overwritten when used with mounted volumes! (which is good)
COPY ./ /home/EMUstack/

# Compile the Fortran code, only use when running tests or copying compiled source to host
WORKDIR /home/EMUstack/backend/fortran/
RUN make

# create useful folder
WORKDIR /home
RUN mkdir host

# getting code and setup env vars
ENV PYTHONPATH "${PYTHONPATH}:/home/EMUstack/backend/"
ENV OPENBLAS_NUM_THREADS=1
ENV OMP_NUM_THREADS=1
RUN echo -e "ulimit -s unlimited" >> /root/.bashrc
RUN source /root/.bashrc

# expose port
EXPOSE 8888
