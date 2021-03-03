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
  long master;
  std::vector <int> availableResolutions;
  int resolution;
  int selectedResolutionId;
  int version;
  CharacterVector chromosomes;
  std::vector <long> chromosomeLengths;
  int totalChromosomes;
  bool firstChromosomeIsAll;
};

// stores output information
struct outputStr {
  std::vector<int> chromosome;
  std::vector<int> bin1;
  std::vector<int> bin2;
  std::vector<int> count;
};

// returns whether or not this is valid HiC file
bool readMagicString(std::istream &fin) {
  std::string str;
  getline(fin, str, '\0');
  return str[0] == 'H' && str[1] == 'I' && str[2] == 'C';
}

void readHeader(std::istream &fin, hicInfo &info) {
  info.selectedResolutionId = -1;
  if (!readMagicString(fin)) {
    stop("Hi-C magic string is missing, does not appear to be a hic file.");
  }
  fin.read((char*) &info.version, sizeof(int));
  if (info.version < 6) {
    stop("Version " + std::to_string(info.version) + " no longer supported.");
  }
  std::string genome;
  int totalAttributes;
  fin.read((char*) &info.master, sizeof(long));
  getline(fin, genome, '\0');
  fin.read((char*) &totalAttributes, sizeof(int));
  // reading and ignoring attribute-value dictionary
  for (int i = 0; i < totalAttributes; i++) {
    std::string key, value;
    getline(fin, key, '\0');
    getline(fin, value, '\0');
  }
  fin.read((char*) &info.totalChromosomes, sizeof(int));
  // chromosome map for finding matrix
  for (int i = 0; i < info.totalChromosomes; i++) {
    std::string name;
    int length;
    getline(fin, name, '\0');
    fin.read((char*) &length, sizeof(int));
    info.chromosomes.push_back(name);
    info.chromosomeLengths.push_back(length);
  }
  int totalResolutions;
  fin.read((char*) &totalResolutions, sizeof(int));
  for (int i = 0; i < totalResolutions; i++) {
    int resolution;
    fin.read((char*) &resolution, sizeof(int));
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
  long position,
  int size,
  int chromosomeId,
  hicInfo &info,
  outputStr &output
) {

  if (size == 0) {
    return;
  }
  std::vector<int> chromosomeIds, bins1, bins2, counts;

  char* compressedBytes = new char[size];
  char* uncompressedBytes = new char[size*10]; // biggest seen so far is 3

  fin.seekg(position, std::ios::beg);
  fin.read(compressedBytes, size);

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
  std::istream bufferin(&sbuf);
  int totalRecords;
  bufferin.read((char*) &totalRecords, sizeof(int));
  bins1.reserve(totalRecords);
  bins2.reserve(totalRecords);
  counts.reserve(totalRecords);
  // different versions have different specific formats
  if (info.version < 7) {
    for (int i = 0; i < totalRecords; i++) {
      int binX, binY;
      float count;
      bufferin.read((char*) &binX, sizeof(int));
      bufferin.read((char*) &binY, sizeof(int));
      bufferin.read((char*) &count, sizeof(float));
      bins1.push_back(binX);
      bins2.push_back(binY);
      counts.push_back(count);
    }
  } else {
    int binXOffset, binYOffset;
    char useShort;
    char type;
    bufferin.read((char*) &binXOffset, sizeof(int));
    bufferin.read((char*) &binYOffset, sizeof(int));
    bufferin.read((char*) &useShort, sizeof(char));
    bufferin.read((char*) &type, sizeof(char));
    if (type == 1) {
      // List-of-rows representation
      short totalRows;
      bufferin.read((char*) &totalRows, sizeof(short));
      for (int i = 0; i < totalRows; i++) {
        short y;
        int binY;
        short totalColumns;
        bufferin.read((char*) &y, sizeof(short));
        binY = y + binYOffset;
        bufferin.read((char*) &totalColumns, sizeof(short));
        for (int j = 0; j < totalColumns; j++) {
          short x;
          int binX;
          short c;
          float count;
          bufferin.read((char*) &x, sizeof(short));
          binX = binXOffset + x;
          if (useShort == 0) { // yes this is opposite of usual
            bufferin.read((char*) &c, sizeof(short));
            count = c;
            bins1.push_back(binX);
            bins2.push_back(binY);
            counts.push_back(c);
          } else {
            bufferin.read((char*) &count, sizeof(float));
            bins1.push_back(binX);
            bins2.push_back(binY);
            counts.push_back(count);
          }
        }
      }
    } else if (type == 2) {
      // have yet to find test file where this is true
      // possibly entirely deprecated
      int totalPoints;
      short w;
      bufferin.read((char*) &totalPoints, sizeof(int));
      bufferin.read((char*) &w, sizeof(short));

      for (int i = 0; i < totalPoints; i++) {
        int row = i / w;
        int column = i - row * w;
        int bin1 = binXOffset + column;
        int bin2 = binYOffset + row;
        float count;
        short c;
        if (useShort == 0) { // yes this is opposite of the usual
          bufferin.read((char*) &c, sizeof(short));
          if (c != -32768) {
            bins1.push_back(bin1);
            bins2.push_back(bin2);
            counts.push_back(c);
          }
        } else {
          bufferin.read((char*) &count, sizeof(float));
          if (count != 0x7fc00000) { // not sure this works
            bins1.push_back(bin1);
            bins2.push_back(bin2);
            counts.push_back(count);
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
  long start,
  int size,
  hicInfo &info,
  outputStr &output
) {

  std::streampos pos;
  if (start != -1) {
    fin.seekg(start, std::ios::beg);
    int chromosomeId1, chromosomeId2, totalResolutions;
    fin.read((char*) &chromosomeId1, sizeof(int));
    fin.read((char*) &chromosomeId2, sizeof(int));
    if (chromosomeId1 == chromosomeId2) {
      fin.read((char*) &totalResolutions, sizeof(int));
      for (
        int resolutionId = 0;
        resolutionId < totalResolutions;
        ++resolutionId
      ) {
        std::string unit;
        int resIdx;
        float tmp2;
        int binSize;
        int totalBlockBins;
        int totalBlockColumns;
        int totalBlocks;
        getline(fin, unit, '\0');
        fin.read((char*) &resIdx, sizeof(int));
        fin.read((char*) &tmp2, sizeof(float)); // sumCounts
        fin.read((char*) &tmp2, sizeof(float)); // occupiedCellCount
        fin.read((char*) &tmp2, sizeof(float)); // stdDev
        fin.read((char*) &tmp2, sizeof(float)); // percent95
        fin.read((char*) &binSize, sizeof(int));
        fin.read((char*) &totalBlockBins, sizeof(int));
        fin.read((char*) &totalBlockColumns, sizeof(int));
        fin.read((char*) &totalBlocks, sizeof(int));
        for (int i = 0; i < totalBlocks; i++) {
          int blockId, blockSize;
          long blockPosition;
          fin.read((char*) &blockId, sizeof(int));
          fin.read((char*) &blockPosition, sizeof(long));
          fin.read((char*) &blockSize, sizeof(int));
          if (resolutionId == info.selectedResolutionId) {
            pos = fin.tellg();
            readBlock(
              fin, blockPosition, blockSize, chromosomeId1, info, output
            );
            fin.seekg(pos, std::ios::beg);
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
  fin.seekg(info.master, std::ios::beg);
  int totalBytes;
  fin.read((char*) &totalBytes, sizeof(int));
  int totalEntries;
  fin.read((char*) &totalEntries, sizeof(int));
  for (int i = 0; i < totalEntries; i++) {
    std::string str;
    getline(fin, str, '\0');
    long fpos;
    fin.read((char*)& fpos, sizeof(long));
    int sizeInBytes;
    fin.read((char*)& sizeInBytes, sizeof(int));
    pos = fin.tellg();
    readMatrix(fin, fpos, sizeInBytes, info, output);
    fin.seekg(pos, std::ios::beg);
  }
}


// [[Rcpp::export]]
DataFrame parseHic(std::string &fname, int resolution) {
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
    // factors start with in R
    bins1 = bins1 - 1;
    bins2 = bins2 - 1;
  }
  chromosomes.attr("class") = "factor";
  chromosomes.attr("levels") = info.chromosomes;
  return DataFrame::create(
    _["chromosome"] = chromosomes,
    _["position.1"] = bins1 * resolution,
    _["position.2"] = bins2 * resolution,
    _["interaction"] = counts
  );
}
