FROM ubuntu:20.04 AS sqlite
RUN apt-get update && apt-get install --no-install-recommends -qq wget ca-certificates make gcc g++ libbz2-1.0 libbz2-dev lbzip2 zlib1g-dev liblzma-dev
WORKDIR /
RUN wget https://www.sqlite.org/2022/sqlite-autoconf-3380100.tar.gz
RUN tar -xf sqlite-autoconf-3380100.tar.gz 
RUN mv sqlite-autoconf-3380100 sqlite
WORKDIR /sqlite 
RUN export CFLAGS='-DSQLITE_ENABLE_UPDATE_DELETE_LIMIT=1'; ./configure && make

FROM ubuntu:20.04 AS tabix
WORKDIR /
RUN apt-get update && apt-get install --no-install-recommends -qq wget ca-certificates make gcc g++ libbz2-1.0 libbz2-dev lbzip2 zlib1g-dev liblzma-dev
RUN wget https://github.com/samtools/htslib/releases/download/1.14/htslib-1.14.tar.bz2
RUN tar -xf htslib-1.14.tar.bz2
RUN mv htslib-1.14 htslib
WORKDIR /htslib
RUN ./configure && make
RUN cp bgzip /usr/local/bin/

# QuickBEAST stage
FROM ubuntu:20.04 AS QuickBEAST
RUN apt-get update && apt-get install --no-install-recommends -qq wget ca-certificates make gcc g++ libbz2-1.0 libbz2-dev lbzip2 zlib1g-dev liblzma-dev git libncurses5-dev
RUN apt-get install -qq libgsl-dev
WORKDIR /
RUN git clone https://github.com/x811zou/QuickBEAST.git
WORKDIR /QuickBEAST
RUN make

# Add Samtools installation
RUN wget https://github.com/samtools/samtools/releases/download/1.15/samtools-1.15.tar.bz2 && \
    tar -xjf samtools-1.15.tar.bz2 && \
    cd samtools-1.15 && \
    ./configure && \
    make && \
    make install && \
    cd .. && \
    rm -rf samtools-1.15 samtools-1.15.tar.bz2

# BEASTIE Python stage
FROM ubuntu:20.04 AS beastie-py
RUN apt-get update && apt-get install --no-install-recommends -qq make pipenv python3.8-venv python3.8-dev pkg-config jags gcc g++
WORKDIR /BEASTIE
COPY . .
RUN make clean dist
RUN pip install statsmodels==0.13.5 && pip install Cython==0.29.24 && pip install numpy==1.21.0 && pip install dist/*.whl
RUN pip install psutil==5.9.0  # Add this line to install psutil
RUN pip install biopython==1.79  # Add this line to install Biopython

# R Packages stage
FROM ubuntu:20.04 AS rpackages
ENV DEBIAN_FRONTEND noninteractive
RUN apt-get update \
  && apt-get install -y --no-install-recommends make gcc g++ r-base-core r-base-dev libstdc++6 libicu-dev libxml2-dev libbz2-dev libssl-dev libcurl4-openssl-dev
RUN apt-get update && apt-get install -y \
  libssl-dev \
  libcurl4-gnutls-dev \
  libxml2-dev \
  libsasl2-dev \
  wget

# Install dependencies first
RUN Rscript -e 'install.packages(c("remotes"), repos = "http://cran.r-project.org")'
RUN wget -P / https://cran.r-project.org/src/contrib/Archive/lattice/lattice_0.20-45.tar.gz
RUN R -e 'install.packages("/lattice_0.20-45.tar.gz", repos = NULL, type="source"); if (!library(lattice, logical.return=T)) quit(status=10)'
RUN Rscript -e 'remotes::install_version("Matrix", version = "1.4.1", repos = "http://cran.us.r-project.org")'
RUN Rscript -e 'remotes::install_version("readr", version = "1.3.1", repos = "http://cran.us.r-project.org")'
RUN Rscript -e 'remotes::install_version("glmnetUtils", version = "1.1.8", repos = "http://cran.us.r-project.org")'

# UCSC tools stage
FROM ubuntu:20.04 AS ucsc-tools
WORKDIR /usr/local/bin
RUN apt-get update && apt-get install --no-install-recommends -qq wget ca-certificates && \
    wget http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/faToTwoBit && \
    wget http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/twoBitToFa && \
    chmod 744 faToTwoBit twoBitToFa

# Final image
FROM ubuntu:20.04
ENV DEBIAN_FRONTEND noninteractive
ENV CYTHONIZE 1
RUN apt-get update \
  && apt-get install -y --no-install-recommends \
  r-base-core r-base-dev libicu66 libstdc++6 openssl libxml2 libcurl4 zlib1g \
  libbz2-1.0 lzma libhts3 vcftools bcftools samtools pipenv python3.8-venv libtbb2 \
  r-cran-dplyr r-cran-readr krb5-user sssd-krb5 jags libgsl23

RUN ln -s /usr/bin/python3.8 /usr/bin/python

COPY --from=QuickBEAST /QuickBEAST/QuickBEAST /usr/local/bin
COPY --from=QuickBEAST /usr/local/bin/samtools /usr/local/bin
COPY --from=tabix /htslib/tabix /usr/local/bin
COPY --from=tabix /usr/local/bin/bgzip /usr/local/bin
COPY --from=rpackages /usr/local/lib/R/site-library /usr/local/lib/R/site-library
COPY --from=sqlite /sqlite/.libs/libsqlite* /usr/local/lib/
COPY --from=beastie-py /usr/local/lib/python3.8/dist-packages /usr/local/lib/python3.8/dist-packages
COPY --from=beastie-py /usr/local/bin/beastie /usr/local/bin/
COPY --from=ucsc-tools /usr/local/bin/faToTwoBit /usr/local/bin/
COPY --from=ucsc-tools /usr/local/bin/twoBitToFa /usr/local/bin/

ENV LD_LIBRARY_PATH=/usr/local/lib:$LD_LIBRARY_PATH


CMD ["python3", "./bin/beastie"]
