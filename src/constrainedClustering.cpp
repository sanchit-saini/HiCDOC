#include <Rcpp.h>
#include <random>

using namespace Rcpp;

double getDistance(
  const NumericVector &vector1,
  const NumericVector &vector2
) {

  return sqrt(sum(pow(vector1 - vector2, 2)));
}

double getCentroidsDelta(
  const std::vector<NumericVector> &previousCentroids,
  const std::vector<NumericVector> &centroids
) {

  double delta = 0.0;

  for (unsigned int centroidId = 0; centroidId < previousCentroids.size(); centroidId++) {
    delta += getDistance(previousCentroids[centroidId], centroids[centroidId]);
  }

  return delta;
}

NumericVector getMedianVector(
  const std::vector<NumericVector> &vectors
) {

  NumericVector medianVector(vectors[0].size());
  NumericVector buffer(vectors.size());

  for (unsigned int rowId = 0; rowId < vectors[0].size(); rowId++) {
    for (unsigned int columnId = 0; columnId < vectors.size(); columnId++) {
      buffer[columnId] = vectors[columnId][rowId];
    }
    medianVector[rowId] = median(buffer);
  }
  return medianVector;
}

void updateCentroids(
  std::vector<NumericVector> &centroids,
  const IntegerVector &clusters,
  const NumericMatrix &matrix
) {

  std::vector<unsigned int> totalClusterMembers(centroids.size(), 0);

  for (NumericVector &centroid: centroids) {
    centroid = NumericVector(centroid.size(), 0.0);
  }
  for (int vectorId = 0; vectorId < matrix.nrow(); vectorId++) {
    centroids[clusters[vectorId]] += matrix.row(vectorId);
    totalClusterMembers[clusters[vectorId]]++;
  }
  for (unsigned int centroidId = 0; centroidId < centroids.size(); centroidId++) {
    if (totalClusterMembers[centroidId] > 0) {
      centroids[centroidId] = centroids[centroidId] / totalClusterMembers[centroidId];
    }
  }
}

struct NearestCentroid {
  int centroidId;
  float distance;
};
NearestCentroid getNearestCentroid(
  const NumericVector &vector,
  const std::vector<NumericVector> &centroids
) {

  NearestCentroid nearestCentroid;
  nearestCentroid.centroidId = 0;
  nearestCentroid.distance = std::numeric_limits<double>::max();
  double distance;

  for (unsigned int centroidId = 0; centroidId < centroids.size(); centroidId++) {
    if (centroids[centroidId].size() == 0) continue;
    distance = getDistance(vector, centroids[centroidId]);
    if (distance < nearestCentroid.distance) {
      nearestCentroid.distance = distance;
      nearestCentroid.centroidId = centroidId;
    }
  }
  return nearestCentroid;
}

// Initialize centroids with K-means++
void initializeCentroids(
  std::vector<NumericVector> &centroids,
  const NumericMatrix &matrix
) {

  NumericVector distances(matrix.nrow());
  NumericVector row;
  double sum;

  row = matrix.row(rand() % matrix.nrow());
  centroids[0] = clone(row);

  for (unsigned int centroidId = 1; centroidId < centroids.size(); centroidId++) {
    sum = 0;
    for (int vectorId = 0; vectorId < matrix.nrow(); vectorId++) {
      distances[vectorId] = getNearestCentroid(matrix.row(vectorId), centroids).distance;
      sum += distances[vectorId];
    }
    sum = sum * rand() / (RAND_MAX - 1);
    for (int vectorId = 0; vectorId < matrix.nrow(); vectorId++) {
      if ((sum -= distances[vectorId]) > 0) continue;
      row = matrix.row(vectorId);
      centroids[centroidId] = clone(row);
      break;
    }
  }
}

void assignClusters(
  IntegerVector &clusters,
  std::vector<NumericVector> &centroids,
  const NumericMatrix &matrix,
  const IntegerMatrix &links
) {

  IntegerVector link;
  int centroidId;

  for (int linkId = 0; linkId < links.nrow(); linkId++) {

    link = links.row(linkId);
    std::vector<NumericVector> linkedVectors;
    linkedVectors.reserve(link.size());

    for (int vectorId: link) {
      if (vectorId >= matrix.nrow()) {
        throw std::invalid_argument(
          "Link (" +
          std::to_string(vectorId) +
          ") out of range (" +
          std::to_string(matrix.nrow()) +
          ").\n"
        );
      }
      linkedVectors.push_back(matrix.row(vectorId));
    }

    centroidId = getNearestCentroid(getMedianVector(linkedVectors), centroids).centroidId;

    for (int vectorId: link) {
      clusters[vectorId] = centroidId;
    }
  }
}

double clusterize(
  NumericMatrix &matrix,
  IntegerMatrix &links,
  IntegerVector &clusters,
  std::vector<NumericVector> &centroids,
  double maxDelta,
  int maxIterations
) {

  int totalIterations = 0;
  double centroidsDelta;

  std::vector<NumericVector> previousCentroids;

  initializeCentroids(centroids, matrix);

  do {
    previousCentroids = centroids;
    assignClusters(clusters, centroids, matrix, links);
    updateCentroids(centroids, clusters, matrix);
    totalIterations++;
    centroidsDelta = getCentroidsDelta(previousCentroids, centroids);
  } while (
    centroidsDelta > maxDelta
    && totalIterations < maxIterations
  );

  double quality = 0.0;
  for (int vectorId = 0; vectorId < matrix.nrow(); vectorId++) {
    quality += getDistance(matrix.row(vectorId), centroids[clusters[vectorId]]);
  }
  return quality;
}

// [[Rcpp::export]]
List constrainedClustering (
  NumericMatrix &matrix,
  IntegerMatrix &links,
  double maxDelta = 0.0001,
  int maxIterations = 50,
  int totalRestarts = 20,
  int totalClusters = 2
) {

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

  IntegerVector clusters(matrix.nrow()), bestClusters(matrix.nrow());
  std::vector<NumericVector> centroids(totalClusters), bestCentroids(totalClusters);

  if (matrix.nrow() > 0) {
    double quality, minQuality = std::numeric_limits<double>::max();

    for (int restart = 0; restart < totalRestarts; restart++) {

      quality = clusterize(
        matrix, links, clusters, centroids, maxDelta, maxIterations
      );

      if (quality < minQuality) {
        minQuality = quality;
        bestClusters = clone(clusters);
        for (unsigned int centroidId = 0; centroidId < centroids.size(); centroidId++) {
          bestCentroids[centroidId] = clone(centroids[centroidId]);
        }
      }
    }

    if (
      is_false(any(bestClusters == 0))
      || is_false(any(bestClusters == 1))
    ) {
      throw std::invalid_argument(
        "Failed clustering: one of the clusters is empty.\n"
      );
    }
  }

  List output;
  output["clusters"] = bestClusters;
  output["centroids"] = bestCentroids;
  return output;
}
