import cv2
import matplotlib.pyplot as plt

plt.subplots_adjust(wspace=0.5, hspace=0.5)

# Load the image
img = cv2.imread('cafe.jpg')

# Original image
plt.subplot(3,3,1);
plt.imshow(img, cmap="gray")
plt.title("Original Image")

# Convert to RGB
rgb = cv2.cvtColor(img, cv2.COLOR_BGR2RGB)
plt.subplot(3, 3, 2)
plt.imshow(rgb, cmap='gray')
plt.title('RGB Image')

# Convert to grayscale
gray = cv2.cvtColor(img, cv2.COLOR_BGR2GRAY)
plt.subplot(3,3,3);
plt.imshow(gray, cmap="gray")
plt.title("Gray Image")

# Blur the image
blurred = cv2.GaussianBlur(gray, (9, 9), 0)
plt.subplot(3, 3, 4)
plt.imshow(blurred, cmap='gray')
plt.title('Blured Image')

# Detect edges using Canny method
edges = cv2.Canny(gray, threshold1=100, threshold2=200)
plt.subplot(3,3,5)
plt.imshow(edges, cmap="gray")
plt.title("Edges")

# Display the image with corners
img[edges == 255] = (255,0,0)
plt.subplot(3,3,6)
plt.imshow(img, cmap='gray')
plt.title('Image with corner')

# Flatten the edge features
edges = edges.flatten()

print(edges)