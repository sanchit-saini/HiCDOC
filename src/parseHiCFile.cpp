// Based on
// https://github.com/aidenlab/straw/blob/master/R/src/straw-R.cpp

/*
 The MIT License (MIT)

 Copyright (c) 2011-2016 Broad Institute, Aiden Lab

 Permission is hereby granted, free of charge, to any person obtaining a copy
 of this software and associated documentation files (the "Software"), to deal
 in the Software without restriction, including without limitation the rights
 to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 copies of the Software, and to permit persons to whom the Software is
 furnished to do so, subject to the following conditions:

 The above copyright notice and this permission notice shall be included in
 all copies or substantial portions of the Software.

 THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 THE SOFTWARE.
*/

/*
 Straw: fast C++ implementation of dump. Not as fully featured as the Java
 version. Reads the .hic file, finds the appropriate matrix and slice of data,
 and outputs as text in sparse upper triangular format.
 Currently only supporting matrices.
*/

#include <Rcpp.h>
#include <cstring>
#include <cstdint>
#include <iostream>
#include <fstream>
#include <sstream>
#include <map>
#include <set>
#include <vector>
#include <streambuf>
#include "zlib.h"

using namespace Rcpp;

// this is for creating a stream from a byte array for ease of use
struct membuf : std::streambuf {
  membuf(char* begin, char* end) {
    this->setg(begin, begin, end);
  }
};

// stores input information
struct hicInfo {
  int64_t master;
  std::vector <int> availableResolutions;
  int resolution;
  int selectedResolutionId;
  int32_t version;
  CharacterVector chromosomes;
  std::vector <long> chromosomeLengths;
  int32_t totalChromosomes;
  bool firstChromosomeIsAll;
  std::streampos pos;
};

// stores output information
struct outputStr {
  std::vector<int> chromosome;
  std::vector<int> bin1;
  std::vector<int> bin2;
  std::vector<int> count;
};

char readCharFromFile(std::istream &fin, hicInfo &info) {
  char tempChar;
  fin.read(&tempChar, 1);
  info.pos += 1;
  return tempChar;
}

int16_t readInt16FromFile(std::istream &fin, hicInfo &info) {
  int16_t tempInt16;
  fin.read((char *) &tempInt16, 2);
  info.pos += 2;
  return tempInt16;
}

int32_t readInt32FromFile(std::istream &fin, hicInfo &info) {
  int32_t tempInt32;
  fin.read((char *) &tempInt32, 4);
  info.pos += 4;
  return tempInt32;
}

int64_t readInt64FromFile(std::istream &fin, hicInfo &info) {
  int64_t tempInt64;
  fin.read((char *) &tempInt64, 8);
  info.pos += 8;
  return tempInt64;
}

float readFloatFromFile(std::istream &fin, hicInfo &info) {
  float tempFloat;
  fin.read((char *) &tempFloat, sizeof(float));
  info.pos += sizeof(float);
  return tempFloat;
}

double readDoubleFromFile(std::istream &fin, hicInfo &info) {
  double tempDouble;
  fin.read((char *) &tempDouble, sizeof(double));
  info.pos += sizeof(double);
  return tempDouble;
}

void readString(std::istream &fin, std::string &s, hicInfo &info) {
  getline(fin, s, '\0');
  info.pos += fin.gcount();
}

void readBytes(std::istream &fin, char *b, int size, hicInfo &info) {
  fin.read(b, size);
  info.pos += size;
}

// The file position is internally managed, because
// Windows do not synchronize read and tellg...
void setFilePos(std::istream &fin, std::streampos pos, hicInfo &info) {
  fin.seekg(pos, std::ios::beg);
  info.pos = pos;
}

// returns whether or not this is valid HiC file
bool readMagicString(std::istream &fin, hicInfo &info) {
  std::string str;
  readString(fin, str, info);
  return str[0] == 'H' && str[1] == 'I' && str[2] == 'C';
}


void readHeader(std::istream &fin, hicInfo &info) {
  info.selectedResolutionId = -1;
  if (!readMagicString(fin, info)) {
    stop("Hi-C magic string is missing, does not appear to be a hic file.");
  }
  info.version = readInt32FromFile(fin, info);
  if (info.version < 6) {
    stop("Version " + std::to_string(info.version) + " no longer supported.");
  }
  std::string genome;
  info.master = readInt64FromFile(fin, info);
  readString(fin, genome, info);
  if (info.version > 8) {
    readInt64FromFile(fin, info); // nviPosition
    readInt64FromFile(fin, info); // nviLength
  }
  int32_t totalAttributes = readInt32FromFile(fin, info);
  // reading and ignoring attribute-value dictionary
  for (int i = 0; i < totalAttributes; i++) {
    std::string key, value;
    readString(fin, key, info);
    readString(fin, value, info);
  }
  info.totalChromosomes = readInt32FromFile(fin, info);
  // chromosome map for finding matrix
  for (int i = 0; i < info.totalChromosomes; i++) {
    std::string name;
    int32_t length;
    readString(fin, name, info);
    if (info.version > 8) {
      length = readInt64FromFile(fin, info);
    } else {
      length = (int64_t) readInt32FromFile(fin, info);
    }
    info.chromosomes.push_back(name);
    info.chromosomeLengths.push_back(length);
  }
  int32_t totalResolutions = readInt32FromFile(fin, info);
  for (int i = 0; i < totalResolutions; i++) {
    int32_t resolution = readInt32FromFile(fin, info);
    info.availableResolutions.push_back(resolution);
    if (resolution == info.resolution) {
      info.selectedResolutionId = i;
    }
  }
  info.firstChromosomeIsAll = (
    info.chromosomes[0] == "ALL" || info.chromosomes[0] == "All"
  );
}


// This is the meat of reading the data. Takes in the block number and returns
// the set of contact records corresponding to that block. The block data is
// compressed and must be decompressed using the zlib library functions.
void readBlock(
  std::istream &fin,
  int64_t position,
  int32_t size,
  int32_t chromosomeId,
  hicInfo &info,
  outputStr &output
) {

  if (size == 0) {
    return;
  }
  std::vector<int> chromosomeIds, bins1, bins2, counts;

  char* compressedBytes = new char[size];
  char* uncompressedBytes = new char[size*10]; // biggest seen so far is 3

  setFilePos(fin, position, info);
  readBytes(fin, compressedBytes, size, info);

  // Decompress the block
  // zlib struct
  z_stream infstream;
  infstream.zalloc    = Z_NULL;
  infstream.zfree     = Z_NULL;
  infstream.opaque    = Z_NULL;
  infstream.avail_in  = (uInt)(size); // size of input
  infstream.next_in   = (Bytef *) compressedBytes; // input char array
  infstream.avail_out = (uInt)size*10; // size of output
  infstream.next_out  = (Bytef *)uncompressedBytes; // output char array
  // the actual decompression work
  inflateInit(&infstream);
  inflate(&infstream, Z_NO_FLUSH);
  inflateEnd(&infstream);
  int uncompressedSize = infstream.total_out;

  // create stream from buffer for ease of use
  membuf sbuf(uncompressedBytes, uncompressedBytes + uncompressedSize);
  std::istream bufferIn(&sbuf);
  int32_t totalRecords = readInt32FromFile(bufferIn, info);
  bins1.reserve(totalRecords);
  bins2.reserve(totalRecords);
  counts.reserve(totalRecords);
  // different versions have different specific formats
  if (info.version < 7) {
    for (int i = 0; i < totalRecords; i++) {
      int32_t binX = readInt32FromFile(bufferIn, info);
      int32_t binY = readInt32FromFile(bufferIn, info);
      float  c     = readFloatFromFile(bufferIn, info);
      bins1.push_back(binX);
      bins2.push_back(binY);
      counts.push_back(c);
    }
  } else {
    int32_t binXOffset = readInt32FromFile(bufferIn, info);
    int32_t binYOffset = readInt32FromFile(bufferIn, info);
    bool    useShort   = readCharFromFile(bufferIn, info) == 0; // yes this is opposite of usual
    bool useShortBinX = true;
    bool useShortBinY = true;
    if (info.version > 8) {
      useShortBinX = readCharFromFile(bufferIn, info) == 0;
      useShortBinY = readCharFromFile(bufferIn, info) == 0;
    }
    char type = readCharFromFile(bufferIn, info);
    if (type == 1) {
      if (useShortBinX && useShortBinY) {
        int16_t rowCount = readInt16FromFile(bufferIn, info);
        for (int i = 0; i < rowCount; i++) {
          int32_t binY = binYOffset + readInt16FromFile(bufferIn, info);
          int16_t colCount = readInt16FromFile(bufferIn, info);
          for (int j = 0; j < colCount; j++) {
            int32_t binX = binXOffset + readInt16FromFile(bufferIn, info);
            float c;
            if (useShort) {
              c = readInt16FromFile(bufferIn, info);
            } else {
              c = readFloatFromFile(bufferIn, info);
            }
            bins1.push_back(binX);
            bins2.push_back(binY);
            counts.push_back(c);
          }
        }
      } else if (useShortBinX && !useShortBinY) {
        int32_t rowCount = readInt32FromFile(bufferIn, info);
        for (int i = 0; i < rowCount; i++) {
          int32_t binY = binYOffset + readInt32FromFile(bufferIn, info);
          int16_t colCount = readInt16FromFile(bufferIn, info);
          for (int j = 0; j < colCount; j++) {
            int32_t binX = binXOffset + readInt16FromFile(bufferIn, info);
            float c;
            if (useShort) {
              c = readInt16FromFile(bufferIn, info);
            } else {
              c = readFloatFromFile(bufferIn, info);
            }
            bins1.push_back(binX);
            bins2.push_back(binY);
            counts.push_back(c);
          }
        }
      } else if (!useShortBinX && useShortBinY) {
        int16_t rowCount = readInt16FromFile(bufferIn, info);
        for (int i = 0; i < rowCount; i++) {
          int32_t binY = binYOffset + readInt16FromFile(bufferIn, info);
          int32_t colCount = readInt32FromFile(bufferIn, info);
          for (int j = 0; j < colCount; j++) {
            int32_t binX = binXOffset + readInt32FromFile(bufferIn, info);
            float c;
            if (useShort) {
              c = readInt16FromFile(bufferIn, info);
            } else {
              c = readFloatFromFile(bufferIn, info);
            }
            bins1.push_back(binX);
            bins2.push_back(binY);
            counts.push_back(c);
          }
        }
      } else {
        int32_t rowCount = readInt32FromFile(bufferIn, info);
        for (int i = 0; i < rowCount; i++) {
          int32_t binY = binYOffset + readInt32FromFile(bufferIn, info);
          int32_t colCount = readInt32FromFile(bufferIn, info);
          for (int j = 0; j < colCount; j++) {
            int32_t binX = binXOffset + readInt32FromFile(bufferIn, info);
            float c;
            if (useShort) {
              c = readInt16FromFile(bufferIn, info);
            } else {
              c = readFloatFromFile(bufferIn, info);
            }
            bins1.push_back(binX);
            bins2.push_back(binY);
            counts.push_back(c);
          }
        }
      }
    } else if (type == 2) {
      int32_t nPts = readInt32FromFile(bufferIn, info);
      int16_t w = readInt16FromFile(bufferIn, info);
      for (int i = 0; i < nPts; i++) {
        int32_t row = i / w;
        int32_t col = i - row * w;
        int32_t bin1 = binXOffset + col;
        int32_t bin2 = binYOffset + row;
        if (useShort) {
          int16_t c = readInt16FromFile(bufferIn, info);
          if (c != -32768) {
            bins1.push_back(bin1);
            bins2.push_back(bin2);
            counts.push_back(c);
          }
        } else {
          float c = readFloatFromFile(bufferIn, info);
          if (!std::isnan(c)) {
            bins1.push_back(bin1);
            bins2.push_back(bin2);
            counts.push_back(c);
          }
        }
      }
    }
  }
  chromosomeIds = std::vector<int>(bins1.size(), chromosomeId);
  output.chromosome.insert(
    output.chromosome.end(),
    chromosomeIds.begin(),
    chromosomeIds.end()
  );
  output.bin1.insert(output.bin1.end(),   bins1.begin(),   bins1.end());
  output.bin2.insert(output.bin2.end(),   bins2.begin(),   bins2.end());
  output.count.insert(output.count.end(), counts.begin(),  counts.end());
  delete[] compressedBytes;
  delete[] uncompressedBytes; // don't forget to delete your heap arrays in C++!
}

// Reads the raw binned contact matrix at specified resolution, setting the
// block bin count and block column count.
void readMatrix(
  std::istream &fin,
  int64_t start,
  hicInfo &info,
  outputStr &output
) {

  std::streampos pos;
  if (start != -1) {
    setFilePos(fin, start, info);
    int32_t chromosomeId1 = readInt32FromFile(fin, info);
    int32_t chromosomeId2 = readInt32FromFile(fin, info);
    int32_t totalResolutions = readInt32FromFile(fin, info);
    if (chromosomeId1 == chromosomeId2) {
      if ((! info.firstChromosomeIsAll) || (chromosomeId1 != 0)) {
        for (
          int resolutionId = 0;
          resolutionId < totalResolutions;
          ++resolutionId
        ) {
          std::string unit;
          readString(fin, unit, info);
          readInt32FromFile(fin, info); // resIdx
          readFloatFromFile(fin, info); // sumCounts
          readFloatFromFile(fin, info); // occupiedCellCount
          readFloatFromFile(fin, info); // stdDev
          readFloatFromFile(fin, info); // percent95
          readInt32FromFile(fin, info); // binSize
          readInt32FromFile(fin, info); // totalBlockBins
          readInt32FromFile(fin, info); // totalBlockColumns
          int32_t totalBlocks = readInt32FromFile(fin, info);
          for (int i = 0; i < totalBlocks; i++) {
            readInt32FromFile(fin, info); // blockId
            int64_t blockPosition = readInt64FromFile(fin, info);
            int32_t blockSize     = readInt32FromFile(fin, info);
            if (resolutionId == info.selectedResolutionId) {
              pos = info.pos;
              readBlock(
                fin, blockPosition, blockSize, chromosomeId1, info, output
              );
              setFilePos(fin, pos, info);
            }
          }
        }
      }
    }
  }
}

// Reads the footer from the master pointer location. Takes in the chromosomes,
// norm, unit (BP or FRAG) and resolution or binsize, and sets the file position
// of the matrix and the normalization vectors for those chromosomes at the
// given normalization and resolution.
void readFooter(std::istream& fin, hicInfo &info, outputStr &output) {
  std::streampos pos;
  setFilePos(fin, info.master, info);
  if (info.version > 8) {
    readInt64FromFile(fin, info); // totalBytes
  } else {
    readInt32FromFile(fin, info); // totalBytes
  }
  int32_t totalEntries = readInt32FromFile(fin, info);
  for (int i = 0; i < totalEntries; i++) {
    std::string str;
    readString(fin, str, info);
    int64_t fpos = readInt64FromFile(fin, info);
    readInt32FromFile(fin, info); // sizeInBytes
    pos = info.pos;
    readMatrix(fin, fpos, info, output);
    setFilePos(fin, pos, info);
  }
}


// [[Rcpp::export]]
DataFrame parseHiCFile(std::string &fname, int resolution) {
  hicInfo info;
  outputStr output;
  std::ifstream fin;
  fin.open(fname, std::fstream::in);
  if (!fin) {
    stop("File " + fname + " cannot be opened for reading.");
  }
  info.resolution = resolution;
  readHeader(fin, info);
  if (info.selectedResolutionId == -1) {
    Rcerr << "Cannot find resolution " << resolution << ".\n";
    Rcerr << "Available resolutions:\n";
    for (int resolution: info.availableResolutions) {
      Rcerr << "\t" << resolution << "\n";
    }
    stop("Exiting.");
  }
  readFooter(fin, info, output);
  // Transform C++ vectors to R vectors and factors
  IntegerVector chromosomes, bins1, bins2, counts;
  chromosomes = wrap(output.chromosome);
  bins1  = wrap(output.bin1);
  bins2  = wrap(output.bin2);
  counts = wrap(output.count);
  if (info.firstChromosomeIsAll) info.chromosomes.erase(0);
  else {
    // factors start with 1 in R
    chromosomes = chromosomes + 1;
  }
  chromosomes.attr("class") = "factor";
  chromosomes.attr("levels") = info.chromosomes;
  DataFrame outputR = DataFrame::create(
    _["chromosome"] = chromosomes,
    _["position 1"] = bins1 * resolution,
    _["position 2"] = bins2 * resolution,
    _["interaction"] = counts
  );
  outputR.attr("class") = Rcpp::CharacterVector::create("data.table", "data.frame");
  return outputR;
}
