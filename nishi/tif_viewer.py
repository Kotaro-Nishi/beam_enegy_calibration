import cv2
import sys

arg = sys.argv
path = arg[1]
img = cv2.imread(path)
cv2.imshow("",img)
cv2.waitKey(0)
cv2.destroyAllWindows()