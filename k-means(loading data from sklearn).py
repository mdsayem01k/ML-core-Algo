from sklearn.datasets import load_iris
from sklearn.cluster import KMeans
from sklearn.metrics import accuracy_score

# Load the Iris dataset
iris = load_iris()

data = iris.data
#Showint first five rows of iris data
print(data[:5])

target = iris.target
#Showint target of dataset
print(target)


# Perform k-means clustering with 3 clusters
kmeans = KMeans(n_clusters=3)
kmeans.fit(data)

# Get the predicted labels from k-means
predicted_labels = kmeans.labels_

# Accuracy
accuracy = accuracy_score(target, predicted_labels)
print(f"Accuracy: {accuracy}")
