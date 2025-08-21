#ifndef TIFF_UTILS_HPP
#define TIFF_UTILS_HPP

#include <tiffio.h>
#include <string>
#include <vector>
#include <cstdint>
#include <stdexcept>
#include <cstring>
#include <type_traits>
#include <sstream>

struct MemoryStream
{
    std::vector<uint8_t> data;
    size_t position = 0;
};

static tmsize_t memory_read(thandle_t handle, void *buffer, tmsize_t size)
{
    MemoryStream *stream = static_cast<MemoryStream *>(handle);
    if (stream->position >= stream->data.size()) return 0;

    size_t available = stream->data.size() - stream->position;
    size_t to_read = std::min(static_cast<size_t>(size), available);

    memcpy(buffer, stream->data.data() + stream->position, to_read);
    stream->position += to_read;
    return to_read;
}

static tmsize_t memory_write(thandle_t handle, void *buffer, tmsize_t size)
{
    MemoryStream *stream = static_cast<MemoryStream *>(handle);
    if (stream->position + size > stream->data.size()) stream->data.resize(stream->position + size);

    memcpy(stream->data.data() + stream->position, buffer, size);
    stream->position += size;
    return size;
}

static uint64_t memory_seek(thandle_t handle, uint64_t offset, int whence)
{
    MemoryStream *stream = static_cast<MemoryStream *>(handle);

    switch (whence) {
    case SEEK_SET: stream->position = offset; break;
    case SEEK_CUR: stream->position += offset; break;
    case SEEK_END: stream->position = stream->data.size() + offset; break;
    }

    return stream->position;
}

static int memory_close(thandle_t handle) { return 0; }

static uint64_t memory_size(thandle_t handle)
{
    return static_cast<MemoryStream *>(handle)->data.size();
}

template <typename T>
constexpr uint16_t get_bits_per_sample()
{
    if constexpr (std::is_same_v<T, uint8_t>) return 8;
    if constexpr (std::is_same_v<T, uint16_t>) return 16;
    if constexpr (std::is_same_v<T, uint32_t>) return 32;
    if constexpr (std::is_same_v<T, int8_t>) return 8;
    if constexpr (std::is_same_v<T, int16_t>) return 16;
    if constexpr (std::is_same_v<T, int32_t>) return 32;
    if constexpr (std::is_same_v<T, float>) return 32;
    if constexpr (std::is_same_v<T, double>) return 64;
    else return 8;
}

template <typename T>
constexpr uint16_t get_sample_format()
{
    if constexpr (std::is_integral_v<T> && std::is_unsigned_v<T>) return SAMPLEFORMAT_UINT;
    if constexpr (std::is_integral_v<T> && std::is_signed_v<T>) return SAMPLEFORMAT_INT;
    if constexpr (std::is_floating_point_v<T>) return SAMPLEFORMAT_IEEEFP; 
    else return SAMPLEFORMAT_UINT;
}

std::string generate_imagej_metadata(size_t depth, size_t channels, double voxel_depth, const std::string &unit)
{
    std::ostringstream oss;
    oss << "ImageJ=1.53f\n"
        << "images=" << (depth * channels) << "\n"
        << "channels=" << channels << "\n"
        << "slices=" << depth << "\n"
        << "frames=1\n";

    if (channels > 1 || depth > 1) oss << "hyperstack=true\n";

    oss << "mode=" << (channels > 1 ? "composite" : "grayscale") << "\n"
        << "unit=" << unit << "\n"
        << "finterval=1.0\n"
        << "spacing=" << voxel_depth << "\n";

    return oss.str();
}

uint16_t get_resolution_unit(const std::string &unit)
{
    return (unit == "inch") ? RESUNIT_INCH : RESUNIT_CENTIMETER;
}

double get_scale_factor(const std::string &unit)
{
    if (unit == "µm" || unit == "micrometer" || unit == "micron") return 10000.0;
    if (unit == "mm" || unit == "millimeter") return 10.0;
    return 1.0;
}

static void set_common_tiff_tags(TIFF *tif, size_t width, size_t height, size_t channels, uint16_t bits_per_sample, uint16_t sample_format)
{
    TIFFSetField(tif, TIFFTAG_IMAGEWIDTH, static_cast<uint32_t>(width));
    TIFFSetField(tif, TIFFTAG_IMAGELENGTH, static_cast<uint32_t>(height));
    TIFFSetField(tif, TIFFTAG_SAMPLESPERPIXEL, static_cast<uint16_t>(channels));
    TIFFSetField(tif, TIFFTAG_BITSPERSAMPLE, bits_per_sample);
    TIFFSetField(tif, TIFFTAG_SAMPLEFORMAT, sample_format);
    TIFFSetField(tif, TIFFTAG_ORIENTATION, ORIENTATION_TOPLEFT);
    TIFFSetField(tif, TIFFTAG_PLANARCONFIG, PLANARCONFIG_SEPARATE);
    TIFFSetField(tif, TIFFTAG_COMPRESSION, COMPRESSION_NONE);
    TIFFSetField(tif, TIFFTAG_ROWSPERSTRIP, height);
    TIFFSetField(tif, TIFFTAG_PHOTOMETRIC, PHOTOMETRIC_MINISBLACK);
}

static void set_resolution_tags(TIFF *tif, double pixel_width, double pixel_height, const std::string &unit)
{
    if (pixel_width <= 0 || pixel_height <= 0) return;
    
    double scale_factor = get_scale_factor(unit);
    TIFFSetField(tif, TIFFTAG_XRESOLUTION, static_cast<float>(scale_factor / pixel_width));
    TIFFSetField(tif, TIFFTAG_YRESOLUTION, static_cast<float>(scale_factor / pixel_height));
    TIFFSetField(tif, TIFFTAG_RESOLUTIONUNIT, get_resolution_unit(unit));
}

template <typename T>
static void write_tiff_slice(TIFF *tif, T *data, size_t width, size_t height, size_t channels, size_t z_index, size_t total_depth)
{
    if (total_depth > 1) TIFFSetField(tif, TIFFTAG_PAGENUMBER, static_cast<uint16_t>(z_index), static_cast<uint16_t>(total_depth));

    for (size_t c = 0; c < channels; c++) {
        for (size_t y = 0; y < height; y++) {
            size_t channel_offset = c * (width * height * total_depth);
            size_t z_offset = z_index * (width * height);
            size_t line_offset = y * width;

            T *line = &data[channel_offset + z_offset + line_offset];

            if (TIFFWriteScanline(tif, line, y, c) < 0) {
                throw std::runtime_error("Failed to write TIFF scanline for slice " +
                                         std::to_string(z_index) + ", channel " + std::to_string(c));
            }
        }
    }
}

template <typename T>
static std::string create_tiff_internal(T *data, size_t width, size_t height, size_t depth, size_t channels, bool use_bigtiff, double pixel_width, double pixel_height, double voxel_depth, const std::string &unit)
{
    MemoryStream stream;
    size_t estimated_size = sizeof(T) * width * height * depth * channels + 4096;
    stream.data.reserve(estimated_size);

    const char *mode = use_bigtiff ? "w8" : "w";
    TIFF *tif = TIFFClientOpen("memory", mode, &stream,
                               memory_read, memory_write, memory_seek,
                               memory_close, memory_size, nullptr, nullptr);

    if (!tif) {
        throw std::runtime_error(use_bigtiff ? "Failed to create BigTIFF in memory" : "Failed to create TIFF in memory");
    }

    try {
        uint16_t bits_per_sample = get_bits_per_sample<T>();
        uint16_t sample_format = get_sample_format<T>();

        std::string imagej_metadata;
        if (depth > 1 || channels > 1) imagej_metadata = generate_imagej_metadata(depth, channels, voxel_depth, unit);

        for (size_t z = 0; z < depth; z++) {
            set_common_tiff_tags(tif, width, height, channels, bits_per_sample, sample_format);
            set_resolution_tags(tif, pixel_width, pixel_height, unit);

            if (!imagej_metadata.empty() && z == 0) TIFFSetField(tif, TIFFTAG_IMAGEDESCRIPTION, imagej_metadata.c_str());

            write_tiff_slice(tif, data, width, height, channels, z, depth);

            if (z < depth - 1 && !TIFFWriteDirectory(tif)) {
                throw std::runtime_error("Failed to write TIFF directory for slice " + std::to_string(z));
            }
        }

        TIFFClose(tif);
        return std::string(stream.data.begin(), stream.data.end());
    }
    catch (...) {
        TIFFClose(tif);
        throw;
    }
}

template <typename T>
inline std::string create_auto_tiff(T *data, size_t width, size_t height, size_t depth, size_t channels, double pixel_width = 1.0, double pixel_height = 1.0, double voxel_depth = 1.0, const std::string &unit = "µm")
{
    size_t total_bytes = sizeof(T) * width * height * depth * channels;
    bool use_bigtiff = total_bytes > (4ULL * 1024 * 1024 * 1024 - 100 * 1024 * 1024);
    return create_tiff_internal(data, width, height, depth, channels, use_bigtiff, pixel_width, pixel_height, voxel_depth, unit);
}

#endif