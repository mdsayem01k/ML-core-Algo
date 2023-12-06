import pandas as pd
from sklearn.preprocessing import LabelEncoder
from sklearn.preprocessing import MinMaxScaler


df = pd.read_csv("live.csv")
#Data information
print(df.info())

#Showing first five rows of dataset
print(df.head())

#Drop unnecessary columns
df.drop(['Column1', 'Column2', 'Column3', 'Column4'], axis=1, inplace=True)

#Drop status_id and status_published variable from the dataset
df.drop(['status_id', 'status_published'], axis=1, inplace=True)

#Removing null values from dataset
df = df.dropna()

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

#Showing first five rows of dataset after preprocessing
print(df.head())