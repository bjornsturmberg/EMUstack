# Installing and Launching Docker Desktop: Windows vs. Mac
The installation and launch instructions for Docker Desktop differ slightly between Windows and Mac, and also depending on your specific hardware (Apple Silicon vs. Intel). I include no docker installation instructions for linux, because I assume that if you ar
e using linux, you are also able to install docker. Just one note for linux users: you might want to add your user to the docker group to run images as a regular user: 
```
su - ${USER}
sudo usermod -aG docker ${USER}
```

and then reboot the computer. Other than that, here's a detailed breakdown of the installation instructions for the other platforms:

## Docker Desktop Windows installation:
- **Download**: Head over to https://www.docker.com/products/docker-desktop. Choose the "For Windows" option.
- **Run Installer**: Double-click the downloaded installer and follow the on-screen instructions. Make sure to keep the default options unless you have specific needs.
- **Restart**: Reboot your computer to complete the installation.

## Docker Desktop Windows launch:
- **Search**: Press the Windows key and type "Docker". Open the "Docker Desktop" application.
- **Accept Agreement**: Review and accept the Docker Subscription Service Agreement.
- **Do not create a profile, it is not necessary.**
- **Start Engine**: Docker Desktop will automatically start the Docker engine in the background.

## Mac (Apple Silicon) Installation:
- **Download**: Visit https://www.docker.com/products/docker-desktop and choose the "For Mac (Apple Silicon)" option.
- **Run Installer**: Double-click the downloaded DMG file and proceed as usual.
- **Security Check**: The first time you run Docker Desktop, macOS might perform a security check which can take several minutes.

## Mac (Apple Silicon) launch:
- **Search**: Go to Spotlight (Command + Space) and search for "Docker". Open the "Docker Desktop" application.
- **Accept Agreement**: Review and accept the Docker Subscription Service Agreement.
- **Do not create a profile, it is not necessary.**
- **Start Engine**: Docker Desktop will automatically start the Docker engine in the background.

## Mac (Intel Silicon) Installation:
- **Download**: From https://www.docker.com/products/docker-desktop, choose the "For Mac Intel" option.
- **Run Installer**: Double-click the downloaded DMG file file and proceed as usual.
- **Restart**: Reboot your computer to complete the installation.

## Mac (Intel Silicon) launch:
- **Search**: Use Spotlight (Command + Space) to search for "Docker". Open the "Docker Desktop" application.
- **Accept Agreement**: Review and accept the Docker Subscription Service Agreement.
- **Do not create a profile, it is not necessary.**
- **Start Engine**: Docker Desktop will automatically start the Docker engine in the background.


## Additional Notes:
For both Windows and Mac, check the official Docker documentation for detailed troubleshooting and advanced configurations: https://docs.docker.com/desktop/install/mac-install/ and https://docs.docker.com/desktop/install/windows-install/. Make sure you have virtualization enabled on your system for Docker to function properly. This is usually enabled by default on most modern computers. If you encounter any issues during installation or launch, consult the Docker support resources or online forums for assistance.

## Downloading and building EMUstack 
The main steps for the installation, building and running of EMUstack are:

### Clone the EMUstack repository
Clone the repository as you would usually do.
### Build the docker image
- Launch the Docker-desktop app if you haven't already (if you are on windows or apple)
- Open a terminal:
    - Open a powershell in Windows: Click the Start button, type "PowerShell" in the search bar, and click on "Windows PowerShell". This opens a regular PowerShell window.
    - Mac OS: open a terminal
    - Linux: open a terminal

- From the terminal or from powershell, go to the cloned repository folder.
- Build the docker image (!!!!!!! IMPORTANT !!!!!!!!!!! The build command MUST INCLUDE the final dot, otherwise the image won't build). The build process will take a few minutes. A successful build should end with the line ``naming to docker.io/library/emu``:
    - **Linux**: ``docker build -t emu .``
    - **Windows**: ``docker build -t emu .``
    - **Apple (Intel silicon)**: ``docker build -t emu .``
    - **Apple (Apple silicon)**: ``docker buildx build --platform linux/amd64 -t emu .``

## Running EMUstack
Once the build is finished, you can launch the docker image. Docker provides a lot of flexibility when it comes to running an image. We are going to focus on four aspects:
- We want to have access to a terminal to be able to interact with **EMUstack**.
- We want the ``\home\EMUstack`` folder inside the docker image to be persistent between sessions.
- We want to be able to access the storage of host computer, which will be mounted in the internal ``\home\host`` folder.
- We want to expose the ``8888`` port so that we can launch a ``jupyter`` instance from **docker** and access it from the **host**.

In order to satisfy all these requirements we launch the EMUstack docker image as follows:
 - **Linux, windows and Apple (Intel silicon)**: ``docker run -it --name EMUstack -p host_port:8888 --mount source=emu,target=/home/EMUstack --mount type=bind,source=full_path_host_folder,target=/home/host emu:latest``
 - **Apple (Apple silicon)**: ``docker run -it --name EMUstack -p host_port:8888 --platform linux/amd64 --ulimit stack=33554432:33554432 --mount source=emu,target=/home/EMUstack --mount type=bind,source=full_path_host_folder,target=/home/host emu:latest``

 The meaning of the flags is:
 - ``-it`` runs docker with an interactive shell.
 - ``--name EMUstack`` the docker session is named ``EMUstack``.
 - ``-p host_port:8888``: we expose the internal docker ``8888`` port the host's ``host_port``.
 - ``--mount source=emu,target=/home/EMUstack`` we save the persistent state of the internal ``/home/EMUstack`` folder in the docker volume ``emu``.
 - ``--mount type=bind,source=full_path_host_folder,target=/home/host``: ``full_path_host_folder`` is the folder of the host that we are exposing to docker. It must be specified as a full path.
 - ``emu:latest``: we run the last built image.
 - ``--ulimit stack=33554432:33554432``: stack size to properly run the image on apple silicon.
 - ``--platform linux/amd64``: telling docker that this is an intel x86 architecture image.

 Now all example scripts can be run directly from the container. In alternative one may launch a jupyter instance inside the container and connect to it from the host to run calculation interactively from jupyter.
