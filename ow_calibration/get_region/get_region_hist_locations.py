from ow_calibration.find_25boxes.find_25boxes import find_25boxes
from ow_calibration.load_configuration.load_configuration import load_configuration
from ow_calibration.change_dates.change_dates import change_dates
import numpy as np
import scipy.io as scipy

config = load_configuration()
pa_wmo_boxes = scipy.loadmat('../../data/constants/wmo_boxes.mat')
argo = scipy.loadmat('../../data/climatology/historical_argo/argo_3505.mat')
pn_float_long = 57.1794
pn_float_lat = -59.1868

boxes = find_25boxes(pn_float_long, pn_float_lat, pa_wmo_boxes)


def get_region_hist_locations(pa_wmo_numbers, pa_float_name, pa_config_data):
    # set up matrices to hold data
    grid_lat = np.array([])
    grid_long = []
    grid_dates = []
    not_use = []
    data = []

    # go through each of the WMO boxes
    for wmo_box in pa_wmo_numbers:

        # go through each of the columns denoting whether we should use CTD, bottle, and/or argo
        for data_type in range(1, 4):

            # check if we should use this data. If so, get the data
            if wmo_box[data_type] == 1 and data_type == 1:
                data = scipy.loadmat(config['HISTORICAL_DIRECTORY'] + config['HISTORICAL_CTD_PREFIX'] +
                                     str(int(wmo_box[0])) + '.mat')
                print(data['lat'][0].__len__())

            if wmo_box[data_type] == 1 and data_type == 2:
                data = scipy.loadmat(config['HISTORICAL_DIRECTORY'] + config['HISTORICAL_BOTTLE_PREFIX'] +
                                     str(int(wmo_box[0])) + '.mat')
                print(data['lat'][0].__len__())

            if wmo_box[data_type] == 1 and data_type == 3:
                data = scipy.loadmat(config['HISTORICAL_DIRECTORY'] + config['HISTORICAL_ARGO_PREFIX'] +
                                     str(int(wmo_box[0])) + '.mat')

                # remove the argo float being analysed from the data
                for i in range(0, data['lat'][0].__len__()):
                    if str(data['source'][0][i]).find(pa_float_name) != -1:
                        not_use.append(i)

                data['lat'] = [np.delete(data['lat'], not_use)]
                data['long'] = np.array(np.delete(data['long'], not_use))
                data['dates'] = np.array(np.delete(data['dates'], not_use))
                print(data['lat'][0].__len__())

            if data:
                grid_lat = np.concatenate([grid_lat, data['lat'][0]])
                grid_long.append(data['long'])
                grid_dates.append(data['dates'])

    print(grid_lat.__len__())






boxes = np.array([[3505, 1, 1, 1]])

get_region_hist_locations(boxes, '1900193', 1)
