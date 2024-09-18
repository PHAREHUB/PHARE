ARG RELEASE=40
FROM 129.104.6.165:32219/phare/teamcity-fedora_dep:$RELEASE

ENV OMPI_ALLOW_RUN_AS_ROOT=1
ENV OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1
ENV TMPDIR=/tmp
ENV MODULEPATH=/etc/scl/modulefiles:/etc/scl/modulefiles:/usr/share/Modules/modulefiles:/etc/modulefiles:/usr/share/modulefiles
ENV PYTHONPATH="/root/phare:/root/phare/build:/root/phare/pyphare"

WORKDIR /root
COPY run.sh /root
RUN chmod +x run.sh && ./run.sh
