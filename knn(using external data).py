import pandas as pd
from sklearn.metrics import accuracy_score
from sklearn.model_selection import train_test_split
from sklearn.neighbors import KNeighborsClassifier
from sklearn.preprocessing import LabelEncoder

#Loading data
df = pd.read_csv('iris.csv')

#Showing first five rows of dataset
print(df.head())

#Dividing data into features and labels
y = df['species'].values
X = df.drop(['species'], axis=1).values

# Label encoding
le = LabelEncoder()
y = le.fit_transform(y)

#Spliting dataset into training set and test set
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size = 0.2, random_state = 0)

#Making predictions
classifier = KNeighborsClassifier(n_neighbors=3)
classifier.fit(X_train, y_train)
y_pred = classifier.predict(X_test)

accuracy = accuracy_score(y_test, y_pred)*100
print(accuracy)