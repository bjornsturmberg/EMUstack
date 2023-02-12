Using EMUstack with Docker
========================

Using a Docker container allows the dependencies of EMUstack to be installed in an isolated environment independant of your operating system. The current image is based on Ubuntu 20.04.

Start by obtaining the docker image and source code. Pull the docker image:

```
docker pull morblockdock/emustack
```

The docker image can be run simply with:

```
docker run -it morblockdock/emustack zsh
```

from which you can navigate to the example files and run them as described
elsewhere in the docs.

To have local access to the source files and simulation results the process is a
bit more involved.

If on linux, you must either run the following commands in sudo, or first grant yourself sudo permissions (otherwise you'll get the
`docker: Got permission denied while trying to connect to the Docker daemon socket` error).

```
su - ${USER}
sudo usermod -aG docker ${USER}
```

Now, create a directory `EMUstack-local`, and clone the `EMUstack` repo and `cd` into the source code

```
mkdir EMUstack-local && cd EMUstack-local
git clone https://github.com/bjornsturmberg/EMUstack.git
cd EMUstack
```

Now we run the image as before, but in addition we mount the local directory,
allowing us to modify the files and keep any output files directly:

```
docker run -v $(pwd)/:/home/EMUstack/ -it morblockdock/emustack zsh
```

This runs the docker container, mounting the local directory to the path
`/home/EMUstack/` in the container, and runs the zsh shell.

The locally mounted code has overwritten what was in the in the container, so for
the **first** use it is required to compile the fortran code. Navigate to the backend directory and run `make` to compile (should take about 30 seconds). 

```
cd backend/fortran && make
```

The compiled code is now updated locally so does not need to be run again when starting the container in the future. 

Scripts can be run directly in a container by passing the path, for
example:

```
docker run -v $(pwd)/:/home/EMUstack/ -it morblockdock/emustack python3 /home/EMUstack/examples/simo_010-single_interface.py
```

Runs the `simo_010-single_interface.py` script directly from the local commandline
(without having to enter the docker container manually).

Alternatively, one can access a `jupyterlab` instance and run `EMUstack` in interactive jupyter notebooks from docker with the following:

```
docker run -v $(pwd)/:/home/EMUstack/ -itp 8888:8888 morblockdock/emustack
```

This will spawn an address to be used in a browser as per a typical jupyter instance. The terminal in `jupyterlab` can also be used to run the python scripts directly if desired. Note that on recent versions of WSL the standard ip of `127.0.0.1` in the jupyter address may need to be replaced with `[::1]` due to [ipv6](https://github.com/microsoft/WSL/issues/4983) related issues. 

If you are mostly doing development of scripts a preferred use case may be running the code in the docker container python environment directly. This can be performed by running the container in an interactive session and using the remote environment features of IDEs like [vscode](https://code.visualstudio.com/docs/remote/containers-tutorial) or [pycharm](https://www.jetbrains.com/help/pycharm/docker.html).

Finally, if you wish to build your own docker image, the dockerfile
is provided as an example. You can build a docker
image from the dockerfile using:

```
docker build -t "my_emustack" .
```

On the first build many layers will be pulled, so this may be slow depending on the internet connection. The docker image will copy and compile the local code. The above examples could then be run by swapping the `morblockdock/emustack` the new tagged image `my_emustack`.
