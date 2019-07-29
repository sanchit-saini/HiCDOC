#include"Rcpp.h"
using namespace Rcpp;

double getDistance (NumericVector &v1, NumericVector &v2) {
  return sqrt(sum(pow(v1 - v2, 2))); 
  
}

double getCentroidDistance (std::vector<NumericVector> &centroids1,
                            std::vector<NumericVector> &centroids2) {
  double distance = 0;
  for (unsigned int i = 0; i < centroids1.size(); ++i) {
    distance += getDistance(centroids1[i], centroids2[i]);
  }
  return distance;
}

void getMedianValue(std::vector<NumericVector> &selectedRows,
                    NumericVector &medianValue) {
  NumericVector tmpVector(selectedRows.size());
  for (int i = 0; i < medianValue.size(); ++i) {
    for (unsigned int j = 0; j < selectedRows.size(); ++j) {
      tmpVector[j] = selectedRows[j][i];
    }
    medianValue[i] = median(tmpVector);
  }
}

void placeMiddle (std::vector<NumericVector> &centroids,
                  IntegerVector &clusters,
                  NumericMatrix &matrix) {
  std::vector<int> ns(2, 0);
  for (NumericVector &centroid: centroids) {
    centroid = NumericVector(centroid.size(), 0.0);
  }
  for (int i = 0; i < matrix.nrow(); ++i) {
    centroids[clusters[i]] += matrix.row(i);
    ++ns[clusters[i]];
  }
  for (int i = 0; i < 2; ++i) {
    centroids[i] = centroids[i] / ns[i];
  }
}

void assignClusters(NumericMatrix &matrix,
                    IntegerMatrix &links,
                    std::vector<NumericVector> &centroids,
                    IntegerVector &clusters) {
  NumericVector medianValue(matrix.ncol());
  for (int linkId = 0; linkId < links.nrow(); ++linkId) {
    IntegerVector link = links.row(linkId);
    std::vector<NumericVector> selectedRows;
    selectedRows.reserve(link.size());
    for (int rowId: link) {
      if (rowId >= matrix.nrow()) {
        throw std::invalid_argument("Link (" + std::to_string(rowId) + ") out of range (" + std::to_string(matrix.nrow()) + ").\n");
      }
      selectedRows.push_back(matrix.row(rowId));
    }
    //Rcout << "median of " << link << ": " << medianValue << "\n";
    getMedianValue(selectedRows, medianValue);
    //Rcout << "    Median " << linkId << "/" << links.nrow() << " computed\n";
    int clusterId = 0;
    double minDistance = std::numeric_limits<double>::max();
    for (int centroidId = 0; centroidId < 2; ++centroidId) {
      double distance = getDistance(medianValue, centroids[centroidId]);
      if (distance < minDistance) {
        minDistance = distance;
        clusterId = centroidId;
      }
    }
    for (int rowId: link) {
      clusters[rowId] = clusterId;
    }
  }
  // Wierd: have to use "is_false"
  //Rcout << clusters << "\n";
  if ((is_false(any(clusters == 0))) || (is_false(any(clusters == 1)))) {
    throw std::invalid_argument("Problem with the clustering: one cluster is empty.\n");
  }
}

// [[Rcpp::export]]
IntegerVector constrainedClustering (NumericMatrix &matrix,
                                     IntegerMatrix &links,
                                     double maxDistance,
                                     int maxIterations) {
  if (any(is_na(matrix))) {
    throw std::invalid_argument("Matrix should not contain NAs.");
  }
  if (any(is_na(links))) {
    throw std::invalid_argument("Links should not contain NAs.");
  }
  if (any(is_nan(matrix))) {
    throw std::invalid_argument("Matrix should not contain NANs.");
  }
  if (any(is_nan(links))) {
    throw std::invalid_argument("Links should not contain NANs.");
  }
  //int nDimensions = matrix.ncol();
  int nPoints = matrix.nrow();
  int nIterations = 0;
  double distance;
  std::vector<NumericVector> centroids = {matrix.row(rand() % nPoints),
                                          matrix.row(rand() % nPoints)};
  std::vector<NumericVector> previousCentroids;
  IntegerVector clusters(nPoints);
  do {
    //Rcout << "Iteration " << nIterations << "\n";
    //Rcout << "  Centroid 1: " << centroids[0] << "\n";
    //Rcout << "  Centroid 2: " << centroids[1] << "\n";
    previousCentroids = centroids;
    assignClusters(matrix, links, centroids, clusters);
    //Rcout << "  assigned to clusters\n";
    placeMiddle(centroids, clusters, matrix);
    //Rcout << "  clusters placed\n";
    //Rcout << "  Centroid 1: " << centroids[0] << "\n";
    //Rcout << "  Centroid 2: " << centroids[1] << "\n";
    ++nIterations;
    distance = getCentroidDistance(previousCentroids, centroids);
    //Rcout << "Current state: " << nIterations << " iteration(s) and distance " << distance << "\n";
  }
  while ((distance > maxDistance) &&
         (nIterations < maxIterations));
  //Rcout << "Leaving with " << nIterations << " and distance " << distance << "\n";
  return clusters;
}
