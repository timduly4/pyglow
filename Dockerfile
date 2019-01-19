# Run as:
#
# $ docker build -t pyglow .
#

# Start with a Python base image:
FROM python:3.6-alpine

# Install required Linux tools:
RUN apk update && apk add \
  bash \
  make \
  patch \
  gcc \
  g++ \
  gfortran \
  perl

# Install pyglow Python dependencies:
COPY requirements.txt /
RUN pip3 install -r requirements.txt

# Copy source code into container:
COPY src/ /pyglow/src/
COPY test/ /pyglow/test/
COPY setup.py /pyglow
WORKDIR /pyglow

# Compile & install:
RUN make -C src/pyglow/models source
RUN python3 setup.py install --user

# Run unit tests:
CMD coverage run --source src -m pytest test/
