#include"Rcpp.h"
using namespace Rcpp;

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
#include <cstring>
#include <iostream>
#include <fstream>
#include <iostream>
#include <sstream>
#include <map>
#include <set>
#include <vector>
#include <streambuf>
#include "zlib.h"

/*
 Straw: fast C++ implementation of dump. Not as fully featured as the
 Java version. Reads the .hic file, finds the appropriate matrix and slice
 of data, and outputs as text in sparse upper triangular format.

 Currently only supporting matrices.
 */
// this is for creating a stream from a byte array for ease of use
struct membuf : std::streambuf
{
  membuf(char* begin, char* end) {
    this->setg(begin, begin, end);
  }
};

// stores chromosomes
struct chrType {
  std::string name;
  int size;
};

// stores input information
struct hicInfo {
  long master;
  std::vector <int> availableResolutions;
  int resolution;
  int resolutionIdSelected;
  int version;
  std::vector <chrType> chrTypes;
};

// stores output information
struct outputStr {
  CharacterVector chromosome;
  IntegerVector bin1;
  IntegerVector bin2;
  IntegerVector count;
  void append(std::string &ch, int b1, int b2, int co) {
    chromosome.push_back(ch);
    bin1.push_back(b1);
    bin2.push_back(b2);
    count.push_back(co);
  }
};


// returns whether or not this is valid HiC file
bool readMagicString(std::istream& fin) {
  std::string str;
  getline(fin, str, '\0');
  return str[0]=='H' && str[1]=='I' && str[2]=='C';
}

void readHeader(std::istream& fin, hicInfo &info) {
  info.resolutionIdSelected = -1;
  if (! readMagicString(fin)) {
    stop("Hi-C magic string is missing, does not appear to be a hic file.");
  }
  fin.read((char*)& info.version, sizeof(int));
  Rcerr << "Hi-C format version: " << info.version << "\n";
  if (info.version < 6) {
    stop("Version " + std::to_string(info.version) + " no longer supported.");
  }
  std::string genome;
  int nattributes;
  int nChrs;
  fin.read((char*) &info.master, sizeof(long));
  getline(fin, genome, '\0' );
  fin.read((char*)&nattributes, sizeof(int));
  // reading and ignoring attribute-value dictionary
  for (int i = 0; i < nattributes; i++) {
    std::string key, value;
    getline(fin, key, '\0');
    getline(fin, value, '\0');
  }
  fin.read((char*) &nChrs, sizeof(int));
  // chromosome map for finding matrix
  for (int i = 0; i < nChrs; i++) {
    std::string name;
    int length;
    getline(fin, name, '\0');
    fin.read((char*) &length, sizeof(int));
    //cout << "Chromosome " << name << " (" << chrTypes.size() << "): " << length << "\n";
    info.chrTypes.push_back({ name, length });
  }
  int nBpResolutions;
  fin.read((char*) &nBpResolutions, sizeof(int));
  for (int i = 0; i < nBpResolutions; i++) {
    int bpResolution;
    fin.read((char*) &bpResolution, sizeof(int));
    info.availableResolutions.push_back(bpResolution);
    if (bpResolution == info.resolution) {
      info.resolutionIdSelected = i;
    }
  }
}


// this is the meat of reading the data.
// takes in the block number and returns the set of contact records corresponding to
// that block.
// the block data is compressed and must be decompressed using the zlib library functions
void readBlock(std::istream& fin, long position, int size, int chrId, hicInfo &info, outputStr &output) {
  if (size == 0) {
    return;
  }
  char* compressedBytes = new char[size];
  char* uncompressedBytes = new char[size*10]; //biggest seen so far is 3

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
  // the actual decompression work.
  inflateInit(&infstream);
  inflate(&infstream, Z_NO_FLUSH);
  inflateEnd(&infstream);
  int uncompressedSize = infstream.total_out;

  // create stream from buffer for ease of use
  membuf sbuf(uncompressedBytes, uncompressedBytes + uncompressedSize);
  std::istream bufferin(&sbuf);
  int nRecords;
  bufferin.read((char*) &nRecords, sizeof(int));
  // different versions have different specific formats
  Rcerr << "  Reading " << info.chrTypes[chrId].name << " with " << nRecords << " records\n";
  if (info.version < 7) {
    for (int i = 0; i < nRecords; i++) {
      int binX, binY;
      float counts;
      bufferin.read((char*) &binX,   sizeof(int));
      bufferin.read((char*) &binY,   sizeof(int));
      bufferin.read((char*) &counts, sizeof(float));
      //cout << info.chrTypes[chrId1].name << "\t" << info.chrTypes[chrId2].name << "\t" << (binX * info.resolution) << "\t" << (binY * info.resolution) << "\t" << counts << "\n";
      output.append(info.chrTypes[chrId].name, binX * info.resolution, binY * info.resolution, counts);
    }
  }
  else {
    int binXOffset, binYOffset;
    char useShort;
    char type;
    bufferin.read((char*) &binXOffset, sizeof(int));
    bufferin.read((char*) &binYOffset, sizeof(int));
    bufferin.read((char*) &useShort,   sizeof(char));
    bufferin.read((char*) &type,       sizeof(char));
    Rcerr << "    Block " << binXOffset << "-" << binYOffset << " with format short " << static_cast<int>(useShort) << " and type " << static_cast<int>(type) << "\n";
    if (type == 1) {
      // List-of-rows representation
      short rowCount;
      bufferin.read((char*) &rowCount, sizeof(short));
      // Rcerr << "      # rows: " << rowCount << "\n";
      for (int i = 0; i < rowCount; i++) {
        short y;
        int binY;
        short colCount;
        bufferin.read((char*) &y, sizeof(short));
        binY = y + binYOffset;
        bufferin.read((char*) &colCount, sizeof(short));
        // Rcerr << "        # cols: " << colCount << "\n";
        for (int j = 0; j < colCount; j++) {
          short x;
          int binX;
          float counts;
          bufferin.read((char*) &x, sizeof(short));
          binX = binXOffset + x;
          if (useShort == 0) { // yes this is opposite of usual
            short c;
            bufferin.read((char*) &c, sizeof(short));
            counts = c;
          }
          else {
            bufferin.read((char*) &counts, sizeof(float));
          }
          //cout << info.chrTypes[chrId1].name << "\t" << info.chrTypes[chrId2].name << "\t" << (binX * info.resolution) << "\t" << (binY * info.resolution) << "\t" << counts << "\n";
          output.append(info.chrTypes[chrId].name, binX * info.resolution, binY * info.resolution, counts);
        }
      }
    }
    else if (type == 2) { // have yet to find test file where this is true, possibly entirely deprecated
      int nPts;
      short w;
      bufferin.read((char*) &nPts, sizeof(int));
      bufferin.read((char*) &w, sizeof(short));

      for (int i = 0; i < nPts; i++) {
        int row = i / w;
        int col = i - row * w;
        int bin1 = binXOffset + col;
        int bin2 = binYOffset + row;

        float counts;
        if (useShort == 0) { // yes this is opposite of the usual
          short c;
          bufferin.read((char*) &c, sizeof(short));
          if (c != -32768) {
            //cout << info.chrTypes[chrId1].name << "\t" << info.chrTypes[chrId2].name << "\t" << (bin1 * info.resolution) << "\t" << (bin2 * info.resolution) << "\t" << c << "\n";
            output.append(info.chrTypes[chrId].name, bin1 * info.resolution, bin2 * info.resolution, counts);
          }
        }
        else {
          bufferin.read((char*) &counts, sizeof(float));
          if (counts != 0x7fc00000) { // not sure this works
            //cout << info.chrTypes[chrId1].name << "\t" << info.chrTypes[chrId2].name << "\t" << (bin1 * info.resolution) << "\t" << (bin2 * info.resolution) << "\t" << counts << "\n";
            output.append(info.chrTypes[chrId].name, bin1 * info.resolution, bin2 * info.resolution, counts);
          }
        }
      }
    }
  }
  Rcerr << "  Done.\n";
  delete[] compressedBytes;
  delete[] uncompressedBytes; // don't forget to delete your heap arrays in C++!
}

// reads the raw binned contact matrix at specified resolution, setting the block bin count and block column count
void readMatrix(std::istream& fin, long start, int size, hicInfo &info, outputStr &output) {
  std::streampos pos;
  if (start != -1) {
    fin.seekg(start, std::ios::beg);
    int chrId1, chrId2, nResolutions;
    fin.read((char*) &chrId1, sizeof(int));
    fin.read((char*) &chrId2, sizeof(int));
    if (chrId1 == chrId2) {
      fin.read((char*) &nResolutions, sizeof(int));
      //cout << "  Reading matrix " << chrId1 << "/" << chrId2 << ": " << nResolutions << " resolutions\n";
      for (int resolutionId = 0; resolutionId < nResolutions; ++resolutionId) {
        std::string unit;
        int resIdx;
        float tmp2;
        int binSize;
        int blockBinCount;
        int blockColumnCount;
        int blockCount;
        getline(fin, unit, '\0');
        fin.read((char*) &resIdx,           sizeof(int));
        //cout << "    Resolution # " << resIdx << " in " << unit << "\n";
        fin.read((char*) &tmp2,             sizeof(float)); // sumCounts
        fin.read((char*) &tmp2,             sizeof(float)); // occupiedCellCount
        fin.read((char*) &tmp2,             sizeof(float)); // stdDev
        fin.read((char*) &tmp2,             sizeof(float)); // percent95
        fin.read((char*) &binSize,          sizeof(int));
        fin.read((char*) &blockBinCount,    sizeof(int));
        fin.read((char*) &blockColumnCount, sizeof(int));
        fin.read((char*) &blockCount,       sizeof(int));
        //cout << "    # c blocks " << blockCount << "\n";
        for (int blockCountId = 0; blockCountId < blockCount; ++blockCountId) {
          int blockId, blockSize;
          long blockPosition;
          fin.read((char*) &blockId,       sizeof(int));
          fin.read((char*) &blockPosition, sizeof(long));
          fin.read((char*) &blockSize,     sizeof(int));
          //cout << "      block # " << blockId << ": " << blockPosition << " / " << blockSize << "\n";
          if (resolutionId == info.resolutionIdSelected) {
            //cout << "        selected\n";
            pos = fin.tellg();
            readBlock(fin, blockPosition, blockSize, chrId1, info, output);
            fin.seekg(pos, std::ios::beg);
          }
        }
      }
    }
  }
}

// reads the footer from the master pointer location. takes in the chromosomes,
// norm, unit (BP or FRAG) and resolution or binsize, and sets the file
// position of the matrix and the normalization vectors for those chromosomes
// at the given normalization and resolution
void readFooter(std::istream& fin, hicInfo &info, outputStr &output) {
  std::streampos pos;
  fin.seekg(info.master, std::ios::beg);
  int nBytes;
  fin.read((char*) &nBytes, sizeof(int));
  int nEntries;
  fin.read((char*) &nEntries, sizeof(int));
  for (int i = 0; i < nEntries; i++) {
    std::string str;
    getline(fin, str, '\0');
    long fpos;
    fin.read((char*)& fpos, sizeof(long));
    int sizeinbytes;
    fin.read((char*)& sizeinbytes, sizeof(int));
    pos = fin.tellg();
    readMatrix(fin, fpos, sizeinbytes, info, output);
    fin.seekg(pos, std::ios::beg);
  }
}


// [[Rcpp::export]]
DataFrame parseHic(std::string &fname, int resolution) {
  hicInfo info;
  outputStr output;
  std::ifstream fin;
  fin.open(fname, std::fstream::in);
  if (! fin) {
    stop("File " + fname + " cannot be opened for reading.");
  }
  info.resolution = resolution;
  readHeader(fin, info);
  if (info.resolutionIdSelected == -1) {
    Rcerr << "Cannot find resolution " << resolution << ".\nTerminating here.\n";
    Rcerr << "Available resolutions:\n";
    for (int resolution: info.availableResolutions) {
      Rcerr << "\t" << resolution << "\n";
    }
    stop("Exiting.");
  }
  readFooter(fin, info, output);
  return DataFrame::create(_["chromosome"] = output.chromosome,
                           _["position 1"] = output.bin1,
                           _["position 2"] = output.bin2,
                           _["value"]      = output.count);
}
