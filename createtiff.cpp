#include <chrono>
#include <fstream>
#include <iostream>
#include <vector>

#include <osmium/io/any_input.hpp>
#include <osmium/visitor.hpp>
#include <osmium/handler.hpp>
#include <osmium/geom/mercator_projection.hpp>
#include <osmium/geom/tile.hpp>
#include "osmium/util/progress_bar.hpp"
#include "osmium/io/reader_with_progress_bar.hpp"
#include <cxxopts.hpp>

#include "geotiffio.h"
#include "xtiffio.h"

using namespace std;

#define ZOOM 10

osmium::geom::MercatorProjection mercator;

class Timer {
  public:
  Timer(std::string name) : mName(name) {
    mStartTime = std::chrono::high_resolution_clock::now();
    std::cout << "Start " << mName << std::endl;
  }

  ~Timer() {
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>( std::chrono::high_resolution_clock::now() - mStartTime ).count();
    std::cout << "Finished " << mName << " in " << duration/1000.0 << " seconds." << std::endl;
  }

  private:
  std::chrono::high_resolution_clock::time_point mStartTime;
  std::string mName;
};

class Handler : public osmium::handler::Handler {
    public:
    Handler(vector<vector<uint64_t>> &rows, uint64_t zoom) : mRows(rows), mZoom(zoom) {
    }

    void node(const osmium::Node& node) {
        auto coords = osmium::geom::lonlat_to_mercator(node.location());
        auto tilex = osmium::geom::mercx_to_tilex(mZoom,coords.x);
        auto tiley = osmium::geom::mercy_to_tiley(mZoom,coords.y);
        mRows[tiley][tilex]++; 
    }

    private:
    vector<vector<uint64_t>>& mRows;
    uint64_t mZoom;
};

void SetUpTIFFDirectory(TIFF *tif,uint64_t tileDim);
void SetUpGeoKeys(GTIF *gtif);
void WriteImage(TIFF *tif, const vector<vector<uint16_t>>& rows,uint64_t tileDim);

double MAXCOORD = 20037508.34;

int main(int argc, char **argv)
{
    uint64_t zoom = 10;
    uint64_t tileDim = pow(2,zoom);
    vector<vector<uint64_t>> rows; 
    for (int i = 0; i < tileDim; i++) {
        vector<uint64_t> pixels;
        pixels.resize(tileDim);
        fill(pixels.begin(),pixels.end(),0);
        rows.push_back(move(pixels));
    }

    try {
        osmium::io::ReaderWithProgressBar reader{true, argv[1], osmium::osm_entity_bits::node, osmium::io::read_meta::no};

        Handler handler(rows,zoom);
        osmium::apply(reader, handler);
    } catch (const std::exception& e) {
        std::cerr << e.what() << '\n';
        exit(1);
    }

    string fname = "newgeo.tif";
    TIFF *tif;
    GTIF *gtif;
    
    tif=XTIFFOpen(fname.c_str(),"w");
    if (!tif) exit(1);
    gtif = GTIFNew(tif);
    if (!gtif) exit(1);
    SetUpTIFFDirectory(tif,tileDim);
    SetUpGeoKeys(gtif);

    vector<vector<uint16_t>> normalized;
    for (int i = 0; i < tileDim; i++) {
        vector<uint16_t> pixels;
        pixels.resize(tileDim);
        for (int j = 0; j < tileDim; j++) {
            uint64_t total = (uint64_t)ceil(double(rows[i][j]) / 1000.0);
            if (total > 65535) cout << "WARNING: maximum 16 bit sample exceeded." << endl;
            pixels[j] = total;
        }
        normalized.push_back(move(pixels));
    }

    WriteImage(tif,normalized,tileDim);
    GTIFWriteKeys(gtif);
    GTIFFree(gtif);
    XTIFFClose(tif);
}

void SetUpTIFFDirectory(TIFF *tif, uint64_t tileDim)
{
    double tiepoints[6]={0,0,0,-MAXCOORD,MAXCOORD,0.0};
    double pixscale[3]={MAXCOORD*2/tileDim,MAXCOORD*2/tileDim,0};
    
    TIFFSetField(tif,TIFFTAG_IMAGEWIDTH,    tileDim);
    TIFFSetField(tif,TIFFTAG_IMAGELENGTH,   tileDim);
    TIFFSetField(tif,TIFFTAG_COMPRESSION,   COMPRESSION_ADOBE_DEFLATE);
    TIFFSetField(tif,TIFFTAG_PHOTOMETRIC,   PHOTOMETRIC_MINISBLACK);
    TIFFSetField(tif,TIFFTAG_PLANARCONFIG,  PLANARCONFIG_CONTIG);
    TIFFSetField(tif,TIFFTAG_BITSPERSAMPLE, 16);
    TIFFSetField(tif,TIFFTAG_ROWSPERSTRIP,  20L); // strip should be 8kb
    
    TIFFSetField(tif,TIFFTAG_GEOTIEPOINTS, 6,tiepoints);
    TIFFSetField(tif,TIFFTAG_GEOPIXELSCALE, 3,pixscale);
}

void SetUpGeoKeys(GTIF *gtif)
{
    GTIFKeySet(gtif, GTModelTypeGeoKey, TYPE_SHORT, 1, ModelProjected);
    GTIFKeySet(gtif, GTRasterTypeGeoKey, TYPE_SHORT, 1, RasterPixelIsArea);
    GTIFKeySet(gtif, GTCitationGeoKey, TYPE_ASCII, 0, "WGS 84 / Pseudo-Mercator");
    GTIFKeySet(gtif, GeogCitationGeoKey, TYPE_ASCII, 0, "WGS 84");
    GTIFKeySet(gtif, GeogAngularUnitsGeoKey, TYPE_SHORT,  1, Angular_Degree);
    GTIFKeySet(gtif, ProjectedCSTypeGeoKey, TYPE_SHORT,  1, 3857);
    GTIFKeySet(gtif, ProjLinearUnitsGeoKey, TYPE_SHORT,  1, Linear_Meter);
}

void WriteImage(TIFF *tif,const vector<vector<uint16_t>>& v, uint64_t tileDim)
{
    for (int i=0;i<tileDim;i++)
        if (!TIFFWriteScanline(tif, (void *)v[i].data(), i, 0))
            TIFFError("WriteImage","failure in WriteScanline\n");
}