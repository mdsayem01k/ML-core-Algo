import pandas as pd
import tensorflow as tf
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import LabelEncoder
from sklearn.preprocessing import OneHotEncoder

#Loading data
df = pd.read_csv('iris.csv')

#Showing first five rows of dataset
print(df.head())

# Label encoding
le = LabelEncoder()
df["species"] = le.fit_transform(df["species"])

#Showing first five rows of dataset
print(df.head())

#Shuffle data
df = df.sample(frac=1)

#Dividing data into features and labels
y = df['species'].values
# Reshape to a column vector
y = y.reshape(-1, 1) 
X = df.drop(['species'], axis=1).values

# One-hot encode the target variable
encoder = OneHotEncoder(sparse=False)
y = encoder.fit_transform(y)

# Split the data into training and testing sets
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)

# Create a neural network model
model = tf.keras.Sequential([
    tf.keras.layers.Dense(10, input_shape=(X.shape[1],), activation='relu'),  # Hidden layer with 10 neurons
    tf.keras.layers.Dense(3, activation='softmax')  # Output layer with 3 neurons (for 3 classes)
])

# Compile the model
model.compile(optimizer='adam', loss='categorical_crossentropy', metrics=['accuracy'])

# Train the model
model.fit(X_train, y_train, epochs=5, batch_size=8, validation_data=(X_test, y_test), verbose=1)

# Evaluate the model
loss, accuracy = model.evaluate(X_test, y_test)
print(f"Loss: {loss}, Accuracy: {accuracy}")

# Make predictions
predictions = model.predict(X_test)
print("Sample predictions:")
print(predictions[:5])  # Display predictions for the first 5 samples



