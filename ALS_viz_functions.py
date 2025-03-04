import cv2
import numpy as np

def remove_shadow(image, blur=5, threshBlockSize=11, noisGapKernel=3, inpaintKernel=4 ):
  """
  Parameters:
    blur (int): size of the floating window used to smooth the image
    thresholdBlockSize: Size of a pixel neighborhood that is used to calculate a threshold value for the
.   pixel
    noisGapKernel (int): size of floating window used to fill in gaps in the shadows
.   and remove noise
    inpaintKernel (int): size of floating window used to get the color for a pixel.
    """
  
  # load image 
  if type(image) == str:
    image = cv2.imread(image_path)
  
  # Convert to grayscale
  gray = cv2.cvtColor(image, cv2.COLOR_BGR2GRAY)

  # Apply blur to smooth the image
  blurred = cv2.GaussianBlur(gray, (5, 5), 0)
  
  # Use adaptive thresholding to create a binary image
  thresh = cv2.adaptiveThreshold(gray, 255, cv2.ADAPTIVE_THRESH_GAUSSIAN_C, cv2.THRESH_BINARY_INV, 11, 2)
  
  # Morphological operations to remove noise
  # kernel = np.ones((3, 3), np.uint8)
  # opening = cv2.morphologyEx(thresh, cv2.MORPH_OPEN, kernel, iterations=2)
  
  # Find contours in the binary image
  # contours, hierarchy = cv2.findContours(opening, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)
  
  # Draw the contours on the original image
  # result = image.copy()
  # cv2.drawContours(result, contours, -1, (0, 255, 0), 2)
  
  # Invert the binary image to highlight shadows
  thresh = cv2.bitwise_not(thresh)
  
  # Apply morphological operations to remove noise and fill gaps in shadows
  kernel = np.ones((3, 3), np.uint8)
  thresh = cv2.morphologyEx(thresh, cv2.MORPH_CLOSE, kernel)
  
  # Invert the binary image again to restore original intensity values
  thresh = cv2.bitwise_not(thresh)
  
  # Apply the shadow mask to the original image
  # result = cv2.bitwise_and(image, image, mask=thresh)
  
  result = cv2.inpaint(image, thresh.astype(np.uint8), 4, cv2.INPAINT_TELEA)
  result = cv2.cvtColor(result, cv2.COLOR_BGR2RGB)
  
  cv2.imwrite(image_path[:len(image_path)-4]+"_shadowless.tif", result)
  
  return result


def remove_shadows_simple(image):
  # load image 
  image = cv2.imread(image_path)
  
  # Convert to LAB color space
  lab = cv2.cvtColor(image, cv2.COLOR_BGR2LAB)
  
  # Shadow detection (simple thresholding)
  l_channel, a_channel, b_channel = cv2.split(lab)
  shadow_mask = np.logical_and(l_channel < 80, b_channel < 120)
  
  # Inpainting (using OpenCV's Telea inpainting)
  result = cv2.inpaint(image, shadow_mask.astype(np.uint8), 3, cv2.INPAINT_TELEA)
  result = cv2.cvtColor(result, cv2.COLOR_LAB2RGB)
  
  return result


def enhance_pic(image, brightness=10, contrast=1.2, output_dir=None):
  # load image 
  if type(image) == str:
    image = cv2.imread(image)
  
  # Do the adjusting
  result = cv2.addWeighted(image, contrast, np.zeros(image.shape, image.dtype), 0, brightness) 
  result = cv2.cvtColor(result, cv2.COLOR_BGR2RGB)
  
  # write to file if a directory is specified
  if output_dir is not None:
    cv2.imwrite(output_dir+"enhanced.tif", result)

  # return the more awesome-looking picture  
  return(result)

def inpaint(image):
  
  # Load images
  if type(image) == str:
    image = cv2.imread(image)
  else:
    image = image.astype('uint8')
    image = cv2.cvtColor(image, cv2.COLOR_RGB2BGR)
    
  # make mask
  mask = image[:,:,1]
  mask = (mask == 0).astype('uint8')
  
  #mask = cv2.imread(mask, cv2.IMREAD_GRAYSCALE)
  dst = cv2.inpaint(image,mask,4,cv2.INPAINT_TELEA)
  dst = cv2.cvtColor(dst, cv2.COLOR_BGR2RGB)
  
  return dst
