from ow_calibration.find_25boxes.find_25boxes import find_25boxes
import scipy.io as scipy

pa_wmo_boxes = scipy.loadmat('../../data/constants/wmo_boxes.mat')
pn_float_long = 57.1794
pn_float_lat = -59.1868

boxes = find_25boxes(pn_float_long, pn_float_lat, pa_wmo_boxes)

print(boxes)