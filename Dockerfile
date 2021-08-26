FROM gcr.io/broad-getzlab-workflows/base_image:v0.0.4

WORKDIR build
COPY requirements.txt .
RUN pip3 install -r requirements.txt

WORKDIR /app
COPY bin/ .
