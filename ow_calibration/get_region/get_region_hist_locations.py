from ow_calibration.find_25boxes.find_25boxes import find_25boxes
from ow_calibration.load_configuration.load_configuration import load_configuration
from ow_calibration.change_dates.change_dates import change_dates
import scipy.io as scipy

config = load_configuration()
pa_wmo_boxes = scipy.loadmat('../../data/constants/wmo_boxes.mat')
argo = scipy.loadmat('../../data/climatology/historical_argo/argo_3505.mat')
pn_float_long = 57.1794
pn_float_lat = -59.1868

boxes = find_25boxes(pn_float_long, pn_float_lat, pa_wmo_boxes)


def get_region_hist_locations(pa_wmo_numbers, pa_float_name, pa_config_data):
    # set up matrices to hold data
    grid_lat = []
    grid_long = []
    grid_dates = []

    # go through each of the WMO boxes
    for wmo_box in pa_wmo_numbers:

        # go through each of the columns denoting whether we should use CTD, bottle, and/or argo
        for data_type in range(1, 4):

            # check if we should use this data
            if wmo_box[data_type] == 1 and data_type == 1:
                data = scipy.loadmat(config['HISTORICAL_DIRECTORY'] + config['HISTORICAL_CTD_PREFIX'] +
                                     str(int(wmo_box[0])) + '.mat')
                hist_dates = change_dates(data['dates'][0])
                print(hist_dates)


get_region_hist_locations(boxes, 1, 1)
