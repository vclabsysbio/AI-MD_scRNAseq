# Source Image
 FROM ubuntu:latest

  # Set noninterative mode
 ENV DEBIAN_FRONTEND noninteractive

  # apt update and install global requirements
 RUN apt-get clean all && \
     apt-get update && \
     apt-get upgrade -y && \
     apt-get install -y bash perl curl wget git && \
     apt install -y make && \
     apt install -y zlib1g-dev && \
     apt-get install -y autoconf build-essential manpages-dev cmake libbz2-dev libcurl4-openssl-dev libssl-dev liblzma-dev pkg-config automake autotools-dev libtool libev-dev libncurses-dev libncurses5-dev libncursesw5-dev gcc g++ vcftools

  # apt clean and remove cached source lists
 RUN apt-get clean && \
     rm -rf /var/lib/apt/lists/*

RUN git clone https://github.com/samtools/htslib.git
RUN cd htslib && \
     git submodule update --init --recursive && \
     autoheader && \
     autoconf && \
     ./configure --prefix=/usr/local/ && \
     make && \
     make install
  
# Install popscle
 RUN git clone https://github.com/statgen/popscle.git
 RUN cd popscle && \
     mkdir build && \
     cd build && \
     cmake .. && \
     make
 RUN cp /popscle/bin/popscle /usr/local/bin

#install samtools, BCFtools, HTSlib 
# RUN wget https://github.com/samtools/samtools/releases/download/1.3.1/samtools-1.3.1.tar.bz2 -O samtools.tar.bz2 && \
#	tar -xjvf samtools.tar.bz2 && \
#	cd samtools-1.3.1/ && \
#	make && \
#	make prefix=/usr/local/bin install && \
#	apt remove samtools && \
#	ln -s /usr/local/bin/bin/samtools /usr/bin/samtools && \
#	cd ../
	
 RUN wget https://github.com/samtools/bcftools/releases/download/1.11/bcftools-1.11.tar.bz2 -O bcftools.tar.bz2 && \
	tar -xjvf bcftools.tar.bz2 && \
	cd bcftools-1.11 && \
	make && \
	make prefix=/usr/local/bin install && \
	ln -s /usr/local/bin/bin/bcftools /usr/bin/bcftools && \
	cd ../

 RUN wget https://github.com/samtools/htslib/releases/download/1.10/htslib-1.10.tar.bz2 && \
	tar -xjvf htslib-1.10.tar.bz2 && \
	mv htslib-1.10 htslib && \
	cd htslib/ && \
	autoreconf -i && \
	./configure && \
	make && \
	make install && \
	cd ../

#RUN wget https://github.com/vcftools/vcftools/releases/download/v0.1.16/vcftools-0.1.16.tar.gz -O vcftools.tar.gz && \
#	tar -zxvf vcftools.tar.gz && \
#	cd vcftools-0.1.16 && \
#	make && \
#        make install && \
#        make prefix=/usr/local/bin install && \
#        ln -s /usr/local/bin/bin/vcftools /usr/bin/vcftools && \
#        cd ../	


  # Define default command
 #CMD ["popscle"]
 COPY ./entrypoint.sh /
 RUN chmod 755 /entrypoint.sh
 ENTRYPOINT ["/entrypoint.sh"]
