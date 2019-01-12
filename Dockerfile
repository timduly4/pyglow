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
  perl \
  lapack-dev

# Install pyglow Python dependencies:
COPY requirements.txt /
RUN pip install -r requirements.txt

# Copy source code into container:
COPY src/ /pyglow/src/
COPY test/ /pyglow/test/
COPY setup.py /pyglow
WORKDIR /pyglow

# Compile & install:
RUN make -C src/pyglow/models source
RUN python setup.py install --user

# Run unit tests:
RUN python -m unittest test.test_suite_pyglow
