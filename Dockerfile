FROM ubuntu:latest
CMD ["bash"]
ENV PATH=/opt/smog2/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin
ENV TZ=Europe/Kiev
RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && echo $TZ > /etc/timezone # buildkit
RUN apt-get update
RUN apt-get upgrade -y # buildkit
RUN apt-get install -y --no-install-recommends cpanminus vim     patch make cmake pdl default-jre libxml-simple-perl nano emacs # buildkit
RUN cpanm String::Util Exporter XML::Validator::Schema # buildkit
RUN rm -rf /root/.cpanm # buildkit
RUN rm -rf /var/lib/apt/lists/* # buildkit
RUN mkdir /opt/smog2 # buildkit
RUN mkdir /workdir # buildkit
COPY docker.tmp /opt/smog2 
RUN useradd -rm -d /home/smoguser -s /bin/bash -g root -G sudo -u 1001 smoguser -p "$(openssl passwd -1 smoguser)" # buildkit
RUN chown -R smoguser /opt/* # buildkit
USER smoguser
WORKDIR /home/smoguser
RUN bash /opt/smog2/configure.smog2
