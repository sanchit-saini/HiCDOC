#include"Rcpp.h"
using namespace Rcpp;

double getDistance (const NumericVector &v1, const NumericVector &v2) {
  return sqrt(sum(pow(v1 - v2, 2)));
}

double getCentroidDistance (const std::vector<NumericVector> &centroids1,
                            const std::vector<NumericVector> &centroids2) {
  double distance = 0.0;
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
  // if ((is_false(any(clusters == 0))) || (is_false(any(clusters == 1)))) {
  //   throw std::invalid_argument("Problem with the clustering: one cluster is empty.\n");
  // }
}

double computeDistanceToCentroids (NumericMatrix &matrix,
                                   IntegerVector &clusters,
                                   std::vector<NumericVector> &centroids) {
  double distance = 0.0;
  for (int i = 0; i < matrix.nrow(); ++i) {
    distance += getDistance(matrix.row(i), centroids[clusters[i]]);
  }
  return distance;
}

double clusterize (NumericMatrix &matrix,
                   IntegerMatrix &links,
                   IntegerVector &clusters,
                   double maxDistance,
                   int maxIterations) {
  int nPoints = matrix.nrow();
  int nIterations = 0;
  double deltaCentroidDistance;
  std::vector<NumericVector> centroids = {matrix.row(rand() % nPoints),
                                          matrix.row(rand() % nPoints)};
  std::vector<NumericVector> previousCentroids;
  do {
    previousCentroids = centroids;
    assignClusters(matrix, links, centroids, clusters);
    placeMiddle(centroids, clusters, matrix);
    ++nIterations;
    deltaCentroidDistance = getCentroidDistance(previousCentroids, centroids);
  }
  while ((deltaCentroidDistance > maxDistance) &&
         (nIterations < maxIterations));
  Rcout << "  Done clustering with " <<
    nIterations << "/" << maxIterations <<
    " iterations and distance " <<
      deltaCentroidDistance << "/" << maxDistance << "\n";
  return computeDistanceToCentroids(matrix, clusters, centroids);
}

// [[Rcpp::export]]
IntegerVector constrainedClustering (NumericMatrix &matrix,
                                     IntegerMatrix &links,
                                     double maxDistance,
                                     int maxIterations,
                                     int nRestarts) {
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
  int nPoints = matrix.nrow();
  IntegerVector clusters(nPoints), bestClusters(nPoints);
  double distance, minDistance = std::numeric_limits<double>::max();
  for (int restart = 0; restart < nRestarts; ++restart) {
    distance = clusterize(matrix, links, clusters, maxDistance, maxIterations);
    if (distance < minDistance) {
      minDistance = distance;
      bestClusters = clusters;
    }
    Rcout << "    Distance to cluster: " <<
      distance << "/" <<
        minDistance << "\n";
  }
  if ((is_false(any(bestClusters == 0))) ||
      (is_false(any(bestClusters == 1)))) {
    throw std::invalid_argument("Problem with the clustering:"
                                  " one cluster is empty.\n");
  }
  return bestClusters;
}
