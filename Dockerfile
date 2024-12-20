FROM ubuntu:24.04

ARG CDN_PORT=6000
ARG BUILD_THREAD=64

RUN apt update
RUN apt install -y build-essential libboost-all-dev libsqlite3-dev libasio-dev nasm cmake
RUN apt install -y ffmpeg libswscale-dev libavutil-dev libavcodec-dev libavdevice-dev libavfilter-dev libavformat-dev libavutil-dev libpostproc-dev libswresample-dev 

WORKDIR /app

COPY . .

RUN cd x264; make -j $BUILD_THREAD; cd ..
RUN cd zstd; make -j $BUILD_THREAD; cd ..

RUN cmake .; exit 0
RUN make -j $BUILD_THREAD

EXPOSE ${CDN_PORT}

CMD ["sh", "-c", "ls -lh /data/; ./nTracer_cdn 6000 /data/"]
