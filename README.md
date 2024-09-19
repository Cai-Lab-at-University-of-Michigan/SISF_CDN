# Scalable Image Storage Format (SISF) CDN
A tool for remote access of SISF Files over HTTP.

## Endpoint Definitions

### General Server Maintenance

| Endpoint | Description |
|----------|-------------|
| `/` | Returns the string `Server is up!` for uptime polling.  | 
| `/performance` | Prints a HTML table with function timing, separated by dataset and function. |
| `/access` | Prints a HTML table listing the number of hits, separated by source IP. |
| `/inventory` | Prints a HTML table listing of each dataset hosted on this server. |

### Data Access

| Endpoint | Description |
|----------|-------------|
| `/<dataset id>/info` | Generates a Neuroglancer-formatted info metadata file for given dataset. |
| `/<dataset id>/tracing/x1,y1,z1/x2,y2,z2` | Returns a tracing path between coordinate 1 to coordinate 2 given by the tuples `(x1,y1,z1)` and `(x2,y2,z2)`. |
| `/<dataset id>/meanshift/x,y,z` | Performs the meanshift algorithm from [here](https://academic.oup.com/bioinformatics/article/35/18/3544/5306330) at the coordinate given by `(x,y,z)` and returns a new coordinate. |
| `/<dataset id>/skeleton/<neuron id>` | Generates a Neuroglancer-formatted neuron morphology file. |
| `/<dataset id>/skeleton/info` | Generates a Neuroglancer-formatted info metadata file for the skeletons in the database. |
| `/<dataset id>/skeleton/segment_properties/info` | Generates a Neuroglancer-formatted info metadata file for the skeletons in the database. |
| `/<dataset id>/<res>/<xs>-<xe>_<ys>-<ye>_<zs>-<ze>` | Downloads a Fortran-order image chunk of the region defined by the boundary `[xs,xe)`, `[ys,ye)`, `[zs,ze)` and resolution id `res` |

### Trace manipulation
| Endpoint | Description |
|----------|-------------|
| `/<dataset id>/skeleton_api/delete/<neuron id>` | Deletes neuron indicated by `neuron id`. |
| `/<dataset id>/skeleton_api/replace/<neuron id>` | Deletes neuron indicated by `neuron id` and replaces it with the contents of a SWC file `POST`-ed to the request. |
| `/<dataset id>/skeleton_api/ls` | Returns a list of neurons defined for given dataset, formatted as a JSON response. | 
| `/<dataset id>/skeleton_api/upload` | Uploads a new neuron from the contents of a SWC file `POST`-ed to the request, returns the ID of the added neuron. |
| `/<dataset id>/skeleton_api/get/<neuron id>` | Downloads a SWC-formatted neuron morphology file for neuron indicated by `neuron id` |

### Raw data access
These functions allow accessing tiled data, one tile at a time. Notably, this ignores virtual cropping, allowing pipelines which rely on image overlaps to be built.

| Endpoint | Description |
|----------|-------------|
| `/<dataset id>/raw_access/<c>,<i>,<j>,<k>/info` | Returns a Neuroglancer-formatted info file for the SISF chunk represented by color `c`, and chunk coordinate `i,j,k` | 
| `/<dataset id>/raw_access/<c>,<i>,<j>,<k>/<res>/<xs>-<xe>_<ys>-<ye>_<zs>-<ze>` | Downloads a Fortran-order image chunk of the region defined by the boundary `[xs,xe)`, `[ys,ye)`, `[zs,ze)` and resolution id `res` | 



## C++ Libraries

1. Crow-cpp ([https://crowcpp.org/master/](https://crowcpp.org/master/))
   - see `src/crow.h` for more information.
2. C++ Subprocess from Arun Muralidharan ([https://github.com/arun11299/cpp-subprocess](https://github.com/arun11299/cpp-subprocess))
   - see `src/subprocess.hpp` for more information.
3. JSON for Modern C++ from Niels Lohmann ([https://github.com/nlohmann/json](https://github.com/nlohmann/json))
   - see `src/json.hpp` for more information.

## Local Setup

### Platform
This repository has been extensively tested under Ubuntu versions 22.04LTS and 24.04LTS. Before begining, ensure that the following packages/libraries are installed:
- `build-essential`
- `git`
- `libboost-all-dev`
- `libsqlite3-dev`
- `libasio-dev`
- `nasm`

### Git Submodules

There are two libraries which are imported using Git submodules which are required for building:
1. libzstd [https://github.com/facebook/zstd](https://github.com/facebook/zstd)
2. libx264 [https://www.videolan.org/developers/x264.html](https://www.videolan.org/developers/x264.html)

To import and build these dependencies, use the following shell commands:
```
git submodule init
git submodule update

cd x264
make -j 20
cd ..

cd zstd
make -j 20
cd ..
```

## Docker Setup

With each build, a Docker image is built and released at: [https://hub.docker.com/r/geeklogan/sisf_cdn](https://hub.docker.com/r/geeklogan/sisf_cdn). An example command to run this image is:

```
sudo docker run -d --mount type=bind,source=/<your_data_folder>/,target=/data geeklogan/sisf_cdn
```

## Metadata Schema

### Introduction

The data storage strategy described here relies on two layers of segmentation

### Archive structure

```
(root)
├── meta
|   ├── chunk_0_0_0.0.1X.meta
|   ├── chunk_0_0_0.0.16X.meta
|   ├── chunk_0_0_1.0.1X.meta
|   ├── (...)
├── data
|   ├── chunk_0_0_0.0.1X.data
|   ├── chunk_0_0_1.0.16X.data
|   ├── chunk_0_0_1.0.1X.data
|   ├── (...)
├── metadata.bin
```

**Notes:**
- Individual parts of this structure can be symlinked to different file systems (e.g. tmpfs or a cache SSD)
- 

### Chunk naming scheme

```
chunk_0_0_0.0.1X.data
      | | | |  |
      ^-|-|-|--|-------- x location
        ^-|-|--|-------- y location
          ^-|--|-------- z location
            ^--|-------- channel
               ^-------- downsampling rate
```

### Contents of metadata.bin

```
[uint16_t version]
[uint16_t dtype]
[uint16_t channel_count]
[uint16_t mchunkx]
[uint16_t mchunky]
[uint16_t mchunkz]
[uint64_t resx]
[uint64_t resy]
[uint64_t resz]
[uint64_t sizex]
[uint64_t sizey]
[uint64_t sizez]
```

### Content of image.meta

```
[uint16_t version]
[uint16_t dtype]
[uint16_t channel_count]
[uint16_t compression_type]
[uint16_t chunkx]
[uint16_t chunky]
[uint16_t chunkz]
[uint64_t sizex]
[uint64_t sizey]
[uint64_t sizez]

[uint64_t cropstartx]
[uint64_t cropendx]
[uint64_t cropstarty]
[uint64_t cropendy]
[uint64_t cropstartz]
[uint64_t cropendz]

for i in range(count):
    [uint64_t offset]
    [uint32_t size]
```

### Parameter options

#### `dtype`
- `1 -> uint16`
- `2 -> uint8` (not implemented) 