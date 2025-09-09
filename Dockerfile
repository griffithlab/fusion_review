# bsub -G compute-oncology -q oncology-interactive -Is -a 'docker_build(griffithlab/fusion_review)' -- --tag griffithlab/fusion_review .

FROM python:3.11-slim-bookworm

RUN ["apt-get", "update"]
RUN ["apt-get", "install", "-y", "vim"]

ADD scripts/filter_fusions.py /opt/scripts/filter_fusions.py
ADD scripts/requirements.txt /opt/scripts/requirements.txt
ADD scripts/CancerGeneCensus-Mar2023.tsv /opt/scripts/CancerGeneCensus-Mar2023.tsv

RUN chmod +r /opt/scripts/*


RUN pip install --upgrade pip

RUN pip3 install -r /opt/scripts/requirements.txt
