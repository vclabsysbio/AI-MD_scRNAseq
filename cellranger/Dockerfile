###############################################
### Dockerfile for 10X Genomics Cell Ranger ###
###############################################

# Based on
FROM centos:8

# File Author / Maintainer
MAINTAINER Tiandao Li <litd99@gmail.com>

# Install some utilities
RUN yum install -y \
	file \
	git \
	sssd-client \
	which \
	wget \
	unzip

# Install bcl2fastq
RUN cd /tmp/ && \
	wget http://regmedsrv1.wustl.edu/Public_SPACE/litd/Public_html/pkg/bcl2fastq2-v2.20.0.422-Linux-x86_64.rpm && \
	yum -y --nogpgcheck localinstall bcl2fastq2-v2.20.0.422-Linux-x86_64.rpm && \
	rm -rf bcl2fastq2-v2.20.0.422-Linux-x86_64.rpm
 	
# Install cellranger
RUN wget --no-check-certificate -O cellranger-6.1.2.tar.gz "https://cf.10xgenomics.com/releases/cell-exp/cellranger-6.1.2.tar.gz?Expires=1641400250&Policy=eyJTdGF0ZW1lbnQiOlt7IlJlc291cmNlIjoiaHR0cHM6Ly9jZi4xMHhnZW5vbWljcy5jb20vcmVsZWFzZXMvY2VsbC1leHAvY2VsbHJhbmdlci02LjEuMi50YXIuZ3oiLCJDb25kaXRpb24iOnsiRGF0ZUxlc3NUaGFuIjp7IkFXUzpFcG9jaFRpbWUiOjE2NDE0MDAyNTB9fX1dfQ__&Signature=fmXngMApixaQKb518qz5lOGRiv3VNZTOFzDnUmQSHJmD-a3u5YqW2yyWBgxY5wbKrnZ5u5F3CvY6~OdWlZr6rbwPlQL0vWJqGox8YOg-aYC2vxbVU3zh1GnyO7J2nyS1NQA3qZM~7wq1StkFaMcN3HFm8PGPLryGkkgODE0ef6nauD~4y5~HwSHsVBtP20F1sIQbUfKzEBEtYJNkAnFQERdlrDrCds9z036aKODmDudAHTftTDjfL1~74hG1v7aQEt87tmCcyaawkLlUqqwUwfTjoHYr5J1N4kWbNRprKf4zoltNZF1h6Sjybt7MNZyjltZIzyYBY2k3sVNbHdboPw__&Key-Pair-Id=APKAI7S6A5RYOXBWRPDA" && \	
	mv cellranger-6.1.2.tar.gz /opt/ && \
	cd /opt/ && \
	chmod u=rwx * && \
	tar -xzvf cellranger-6.1.2.tar.gz && \
	rm -f cellranger-6.1.2.tar.gz && \
	chmod u=rwx *

# path
ENV PATH /opt/cellranger-6.1.2:$PATH
