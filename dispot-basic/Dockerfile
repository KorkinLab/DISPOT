FROM python:2.7-alpine
LABEL maintainer "onarykov@wpi.edu"

ENV INSTALL_PATH /dispot
RUN mkdir -p $INSTALL_PATH

WORKDIR $INSTALL_PATH

COPY requirements.txt requirements.txt

RUN apk add --update --no-cache --virtual .build-deps \
        build-base \
        && pip install -r requirements.txt \
        && find /usr/local \
        \( -type d -a -name test -o -name tests \) \
        -o \( -type f -a -name '*.pyc' -o -name '*.pyo' \) \
        -exec rm -rf '{}' + \
        && runDeps="$( \
            scanelf --needed --nobanner --recursive /usr/local \
                    | awk '{ gsub(/,/, "\nso:", $2); print "so:" $2 }' \
                    | sort -u \
                    | xargs -r apk info --installed \
                    | sort -u \
        )" \
        && apk add --virtual .rundeps $runDeps \
        && apk del .build-deps \
        && apk add bash \
        && apk add perl \
	    && rm -rf /var/cache/apk/*
RUN export LANG=en_US.UTF-8


RUN /bin/sh -c "apk add --no-cache bash"

# Install DISPOT 
COPY . .
RUN chmod +x dispot.py
ENTRYPOINT ["./dispot.py", "-o", "data/results/results.tsv"]

# Installing ALT-IN Tool dependencies
#RUN pip install --upgrade pip && \
#	pip install --trusted-host pypi.python.org --no-cache-dir -r requirements.txt

#RUN apk del py-pip build-base python-dev