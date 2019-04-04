FROM r-base:latest

RUN apt update && \
    apt install -y libcairo-dev && \
    R -e "install.packages(c('ggplot2', 'scales', 'svglite'))"


