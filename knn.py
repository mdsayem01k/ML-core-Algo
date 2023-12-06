from collections import Counter

import numpy as np
from sklearn import datasets
from sklearn.model_selection import train_test_split


########################## Classifier Implementation #################################
def euclidian_distance(x1, x2):
    distance = np.sqrt(np.sum((x1-x2)**2))
    return distance

class KNN:
    def __init__(self,k=3):
        self.k = k

    def fit(self, X, y):
        self.X_train = X
        self.y_train = y

    def predict(self, X):
        predictions = [self._predict(x) for x in X]
        return predictions

    def _predict(self, x):
        distances = [euclidian_distance(x, x_train) for x_train in self.X_train]

        indices = np.argsort(distances)[:self.k]

        k_nearest_labels = [self.y_train[i] for i in indices]

        most_common = Counter(k_nearest_labels).most_common()

        return most_common[0][0]

 
############################# Model training and testing ####################################
iris = datasets.load_iris()

X, y = iris.data, iris.target

X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=.2, random_state=1234)

model = KNN(k=5)
model.fit(X_train, y_train)
predictions = model.predict(X_test)
print(predictions)

accuracy = np.sum(predictions==y_test)/len(y_test)
print(accuracy)