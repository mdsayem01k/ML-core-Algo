import pandas as pd
from sklearn.cluster import KMeans
from sklearn.preprocessing import LabelEncoder, MinMaxScaler

df = pd.read_csv("live-data.csv")

#Showing first five rows of dataset
print(df.head())

#Declare feature vector and target variable
X = df
y = df['status_type']


#Convert categorical variable into integers
le = LabelEncoder()
X['status_type'] = le.fit_transform(X['status_type'])
y = le.transform(y)

#Feature Scaling
ms = MinMaxScaler()
X = ms.fit_transform(X)

#K-Means model with two clusters 
kmeans = KMeans(n_clusters=3, random_state=0) 

kmeans.fit(X)

labels = kmeans.labels_

# check how many of the samples were correctly labeled
correct_labels = sum(y == labels)

print('Accuracy score: {0:0.2f}'. format(correct_labels/float(y.size)))