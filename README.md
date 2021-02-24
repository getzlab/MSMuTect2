A template for tool repositories. It includes the bare necessities to build and deploy your tool:

* `Dockerfile`, the build steps for containerizing the tool.
* A deploy script (`deploy.sh`), which automatically versions and uploads your Docker image to the container repository.

## Dockerfile

The Dockerfile template uses the Getz Lab standard base image as the container's base image.
Unless you have very specific reasons to use a different image, this should not be changed.
It both ensures a consistent environment across the lab's containers, and speeds up pulling/pushing
images, since the base image is cached on all wolF worker nodes, the lab standard VM image, and the
container repo.

The Dockerfile template requires you to fill in build steps, and optionally copy any external scripts into the
`/app` directory.

## Deploy script

This script (`deploy.sh`) will build, automatically tag/version, and upload your Docker image to the container registry.
The version is automatically determined by the number of commits to the tool repo. If you are on a branch
other than master, the branch name will automatically be prepended to the version number. This is to distinguish
alternate tool versions or configurations that are not intended to supersede existing versions of the tool.
