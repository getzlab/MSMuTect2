FROM gcr.io/broad-getzlab-workflows/base_image:latest

WORKDIR build
# build steps go here
# remember to clear the build directory!

WORKDIR /app
ENV PATH=$PATH:/app
COPY <your script> .
