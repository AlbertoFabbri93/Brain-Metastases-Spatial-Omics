FROM rocker/tidyverse:4.4

# create container folder for packages
RUN mkdir -p ~/.local/share/renv/

# install renv
RUN R -e 'install.packages("renv")'

# WORKDIR /home/rstudio/workspace

# COPY renv.lock renv.lock

CMD ["/bin/bash"]