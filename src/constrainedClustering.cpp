#include <algorithm>    // std::find, std::all_of
#include <vector>       // std::vector
#include <Rcpp.h>


using namespace Rcpp;

template <class T>
class StdMatrix {
  std::vector<std::vector<T>> matrix;

  public:
  template <class U>
  StdMatrix (U &m): matrix(m.nrow(), std::vector<T>(m.ncol())) {
    for (int i = 0; i < m.nrow(); ++i) {
        std::vector<T> row (m.row(i).begin(), m.row(i).end());
        matrix[i] = row;
    }
  }

  size_t nrow () const {
    return matrix.size();
  }

  size_t ncol () const {
    return matrix.front().size();
  }

  std::vector<T> &row (size_t i) {
    return matrix[i];
  }
};


double getDistance(
  const std::vector<double> &vector1,
  const std::vector<double> &vector2
) {

  double d = 0.0;
  for (size_t i = 0; i < vector1.size(); ++i) {
    double t = vector1[i] - vector2[i];
    d += t * t;
  }
  return sqrt(d);
}

double getCentroidsDelta(
  const std::vector<std::vector<double>> &previousCentroids,
  const std::vector<std::vector<double>> &centroids
) {

  double delta = 0.0;

  for (
    unsigned int centroidId = 0;
    centroidId < previousCentroids.size();
    centroidId++
  ) {
    delta += getDistance(previousCentroids[centroidId], centroids[centroidId]);
  }

  return delta;
}

double getMedianValue(
  std::vector<double> &v
) {
    size_t n = v.size() / 2;
    std::nth_element(v.begin(), v.begin() + n, v.end());
    return v[n];
}

void getMedianVector(
  const std::vector<std::vector<double>> &vectors,
  std::vector<double> &medianVector
) {

  std::vector<double> buffer(vectors.size());

  for (unsigned int rowId = 0; rowId < vectors[0].size(); rowId++) {
    for (unsigned int columnId = 0; columnId < vectors.size(); columnId++) {
      buffer[columnId] = vectors[columnId][rowId];
    }
    medianVector[rowId] = getMedianValue(buffer);
  }
}

void updateCentroids(
  std::vector<std::vector<double>> &centroids,
  std::vector<int> &clusters,
  StdMatrix<double> &matrix
) {

  std::vector<unsigned int> totalClusterMembers(centroids.size(), 0);

  for (std::vector<double> &centroid: centroids) {
    centroid = std::vector<double>(centroid.size(), 0.0);
  }
  for (unsigned int vectorId = 0; vectorId < matrix.nrow(); vectorId++) {
    for (size_t i = 0; i < matrix.ncol(); ++i) {
      centroids[clusters[vectorId]][i] += matrix.row(vectorId)[i];
    }
    totalClusterMembers[clusters[vectorId]]++;
  }
  for (
    unsigned int centroidId = 0;
    centroidId < centroids.size();
    centroidId++
  ) {
    if (totalClusterMembers[centroidId] > 0) {
      for (size_t i = 0; i < centroids[centroidId].size(); ++i) {
        centroids[centroidId][i] = 
          centroids[centroidId][i] / totalClusterMembers[centroidId];
      }
    }
  }
}

struct NearestCentroid {
  int centroidId;
  double distance;
};
NearestCentroid getNearestCentroid(
  std::vector<double> &vector,
  std::vector<std::vector<double>> &centroids
) {

  NearestCentroid nearestCentroid;
  nearestCentroid.centroidId = 0;
  nearestCentroid.distance = std::numeric_limits<double>::max();
  double distance;

  for (
    unsigned int centroidId = 0;
    centroidId < centroids.size();
    centroidId++
  ) {
    if (std::all_of(centroids[centroidId].begin(), centroids[centroidId].end(), [](double i){return i == 0.0;})) continue;
    //if (centroids[centroidId].size() == 0) continue;
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
  std::vector<std::vector<double>> &centroids,
  StdMatrix<double> &matrix
) {

  std::vector<double> distances(matrix.nrow());
  std::vector<double> row;
  double sum;

  row = matrix.row(unif_rand() * matrix.nrow());
  centroids[0] = row;

  for (
    unsigned int centroidId = 1;
    centroidId < centroids.size();
    centroidId++
  ) {
    sum = 0;
    for (unsigned int vectorId = 0; vectorId < matrix.nrow(); vectorId++) {
      distances[vectorId] = getNearestCentroid(
        matrix.row(vectorId),
        centroids
      ).distance;
      sum += distances[vectorId];
    }
    sum *= unif_rand();
    for (unsigned int vectorId = 0; vectorId < matrix.nrow(); vectorId++) {
      if ((sum -= distances[vectorId]) > 0) continue;
      row = matrix.row(vectorId);
      centroids[centroidId] = row;
      break;
    }
  }
}

void assignClusters(
  std::vector<int> &clusters,
  std::vector<std::vector<double>> &centroids,
  StdMatrix<double> &matrix,
  StdMatrix<int> &links
) {

  std::vector<double> medianVector(matrix.ncol());
  int centroidId;

  for (unsigned int linkId = 0; linkId < links.nrow(); linkId++) {

    std::vector<int> &link = links.row(linkId);
    std::vector<std::vector<double>> linkedVectors(link.size(), std::vector<double>(matrix.ncol(), 0.0));

    for (size_t i = 0; i < link.size(); ++i) {
      unsigned int vectorId = link[i];
      if (vectorId >= matrix.nrow()) {
        throw std::invalid_argument(
          "Link (" +
          std::to_string(vectorId) +
          ") out of range (" +
          std::to_string(matrix.nrow()) +
          ").\n"
        );
      }
      linkedVectors[i] = matrix.row(vectorId);
    }

    getMedianVector(linkedVectors, medianVector);
    centroidId = getNearestCentroid(
      medianVector,
      centroids
    ).centroidId;

    for (int vectorId: link) {
      clusters[vectorId] = centroidId;
    }
  }
}

double clusterize(
  StdMatrix<double> &matrix,
  StdMatrix<int> &links,
  std::vector<int> &clusters,
  std::vector<std::vector<double>> &centroids,
  double maxDelta,
  int maxIterations
) {

  int totalIterations = 0;
  double centroidsDelta;

  std::vector<std::vector<double>> previousCentroids;

  initializeCentroids(centroids, matrix);

  do {
    previousCentroids = centroids;
    assignClusters(clusters, centroids, matrix, links);
    updateCentroids(centroids, clusters, matrix);
    totalIterations++;
    centroidsDelta = getCentroidsDelta(previousCentroids, centroids);
  } while (
    (centroidsDelta > maxDelta) && (totalIterations < maxIterations)
  );

  double quality = 0.0;
  for (unsigned int vectorId = 0; vectorId < matrix.nrow(); vectorId++) {
    quality += getDistance(matrix.row(vectorId), centroids[clusters[vectorId]]);
  }
  return quality;
}

// [[Rcpp::export]]
List constrainedClustering(
  NumericMatrix &rMatrix,
  IntegerMatrix &rLinks,
  double maxDelta = 0.0001,
  int maxIterations = 50,
  int totalRestarts = 20,
  int totalClusters = 2
) {

  if (any(is_na(rMatrix))) {
    throw std::invalid_argument("Matrix should not contain NAs.");
  }
  if (any(is_na(rLinks))) {
    throw std::invalid_argument("Links should not contain NAs.");
  }
  if (any(is_nan(rMatrix))) {
    throw std::invalid_argument("Matrix should not contain NANs.");
  }
  if (any(is_nan(rLinks))) {
    throw std::invalid_argument("Links should not contain NANs.");
  }
  if (rMatrix.nrow() == 0) {
    throw std::invalid_argument("Matrix should not be empty.");
  }

  StdMatrix<double> matrix(rMatrix);
  StdMatrix<int>    links(rLinks);
  std::vector<int> clusters(matrix.nrow());
  std::vector<int> bestClusters(matrix.nrow());
  std::vector<std::vector<double>> centroids(totalClusters, std::vector<double>(matrix.ncol(), 0.0));
  std::vector<std::vector<double>> bestCentroids(totalClusters);

/*
for (size_t j = 0; j < links.ncol(); ++j) {
  for (size_t i = 0; i < links.nrow(); ++i) {
    Rcout << links.row(i)[j] << " ";
  }
  Rcout << "\n";
}
*/

  double quality, minQuality = std::numeric_limits<double>::max();

  for (int restart = 0; restart < totalRestarts; restart++) {

    quality = clusterize(
      matrix, links, clusters, centroids, maxDelta, maxIterations
    );

    if (quality < minQuality) {
      minQuality = quality;
      bestClusters = clusters;
      bestCentroids = centroids;
    }
  }

  if ((std::find(bestClusters.begin(),
                 bestClusters.end(),
                 0) == bestClusters.end()) ||
      (std::find(bestClusters.begin(),
                 bestClusters.end(),
                 1) == bestClusters.end())) {
    throw std::invalid_argument(
      "Failed clustering: one of the clusters is empty.\n"
    );
  }

  List output;
  output["clusters"] = wrap(bestClusters);
  output["centroids"] = wrap(bestCentroids);
  return output;
}
