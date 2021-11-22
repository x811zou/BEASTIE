FROM ubuntu:20.10 AS CMDSTAN

WORKDIR /
RUN apt-get update; apt-get install --no-install-recommends -qq wget ca-certificates make gcc g++

RUN wget https://github.com/stan-dev/cmdstan/releases/download/v2.27.0/cmdstan-2.27.0.tar.gz
RUN tar -zxpf cmdstan-2.27.0.tar.gz

RUN mv cmdstan-2.27.0 cmdstan
WORKDIR /cmdstan
RUN make build

COPY BEASTIE/iBEASTIE2.stan .
RUN make iBEASTIE2






FROM ubuntu:20.10 AS tabix

WORKDIR /
RUN apt-get update; apt-get install --no-install-recommends -qq wget ca-certificates make gcc g++ libbz2-1.0 libbz2-dev lbzip2 zlib1g-dev liblzma-dev
RUN wget https://github.com/samtools/htslib/releases/download/1.14/htslib-1.14.tar.bz2
RUN tar -xf htslib-1.14.tar.bz2

RUN mv htslib-1.14 htslib
WORKDIR /htslib
RUN ./configure && make






FROM ubuntu:20.10 AS beastie-py

RUN apt-get update; apt-get install --no-install-recommends -qq make pipenv python3.8-venv

WORKDIR /BEASTIE
COPY . .
RUN make clean dist






FROM ubuntu:20.10
ENV DEBIAN_FRONTEND noninteractive
ENV CYTHONIZE 1
RUN apt-get update \
  && apt-get install -y --no-install-recommends make gcc g++ \
  && apt-get install -y --no-install-recommends r-base-core r-base-dev libicu67 libstdc++6 openssl libxml2 libcurl4 zlib1g libbz2-1.0 lzma libhts3 vcftools samtools pipenv python3.8-venv libtbb2 \
  && apt-get install -y --no-install-recommends libicu-dev libxml2-dev git autoconf zlib1g-dev libbz2-dev libssl-dev libcurl4-openssl-dev

# RUN apk -v update \
#   && apk add R R-dev icu libstdc++ libxml2 musl libcurl zlib libbz2 \
#   && apk add icu-dev make gcc g++ musl-dev curl-dev libxml2-dev git autoconf zlib-dev bzip2 bzip2-dev

RUN ln -s /usr/bin/python3.8 /usr/bin/python

RUN R -e 'install.packages("LDlinkR", dependencies=T, repos="http://cran.us.r-project.org");   if (!library(LDlinkR, logical.return=T)) quit(status=10)' && \
  R -e 'install.packages("reshape2", dependencies=T, repos="http://cran.us.r-project.org");    if (!library(reshape2, logical.return=T)) quit(status=10)' && \
  R -e 'install.packages("data.table", dependencies=T, repos="http://cran.us.r-project.org");  if (!library(data.table, logical.return=T)) quit(status=10)' && \
  R -e 'install.packages("foreach", dependencies=T, repos="http://cran.us.r-project.org");     if (!library(foreach, logical.return=T)) quit(status=10)' && \
  R -e 'install.packages("dplyr", dependencies=T, repos="http://cran.us.r-project.org");       if (!library(dplyr, logical.return=T)) quit(status=10)' && \
  R -e 'install.packages("readr", dependencies=T, repos="http://cran.us.r-project.org");       if (!library(readr, logical.return=T)) quit(status=10)' && \
  R -e 'install.packages("glmnetUtils", dependencies=T, repos="http://cran.us.r-project.org"); if (!library(glmnetUtils, logical.return=T)) quit(status=10)' && \
  R -e 'if (!requireNamespace("BiocManager", quietly=T)) install.packages("BiocManager", dependencies=T, repos="http://cran.us.r-project.org"); BiocManager::install("pasilla"); if (!library(pasilla, logical.return=T)) quit(status=10)'

RUN apt-get purge -y libicu-dev make gcc g++ libxml2-dev git autoconf zlib1g-dev libbz2-dev libssl-dev libcurl4-openssl-dev

COPY --from=CMDSTAN /cmdstan/iBEASTIE2 /usr/local/bin
COPY --from=tabix /htslib/tabix /usr/local/bin
COPY --from=beastie-py /BEASTIE/dist/*.whl /tmp

RUN pip3 install /tmp/*.whl && rm -f /tmp/*.whl

ENTRYPOINT ["beastie"]
