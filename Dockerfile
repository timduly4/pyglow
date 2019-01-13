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
RUN pip install -r requirements.txt

# Copy source code into container:
COPY src/ /pyglow/src/
COPY test/ /pyglow/test/
COPY setup.py /pyglow
WORKDIR /pyglow

# Compile & install:
RUN make -C src/pyglow/models source
RUN python setup.py install --user

# Create Geophysical indices data structure:
RUN python -c "import pyglow"

# Run unit tests:
CMD python -m unittest test.test_suite_pyglow
