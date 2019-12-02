from ow_calibration.find_25boxes.find_25boxes import find_25boxes
from ow_calibration.load_configuration.load_configuration import load_configuration
from ow_calibration.change_dates.change_dates import change_dates
import numpy as np
import scipy.io as scipy

config = load_configuration()
pa_wmo_boxes = scipy.loadmat('../../data/constants/wmo_boxes.mat')
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

            if wmo_box[data_type] == 1 and data_type == 2:
                data = scipy.loadmat(config['HISTORICAL_DIRECTORY'] + config['HISTORICAL_BOTTLE_PREFIX'] +
                                     str(int(wmo_box[0])) + '.mat')

            if wmo_box[data_type] == 1 and data_type == 3:
                data = scipy.loadmat(config['HISTORICAL_DIRECTORY'] + config['HISTORICAL_ARGO_PREFIX'] +
                                     str(int(wmo_box[0])) + '.mat')

                # remove the argo float being analysed from the data
                for i in range(0, data['lat'][0].__len__()):
                    if str(data['source'][0][i]).find(pa_float_name) != -1:
                        not_use.append(i)

                data['lat'] = [np.delete(data['lat'], not_use)]
                data['long'] = [np.delete(data['long'], not_use)]
                data['dates'] = [np.delete(data['dates'], not_use)]

            # if we have data, combine it with the other data
            if data:
                grid_lat = np.concatenate([grid_lat, data['lat'][0]])
                grid_long = np.concatenate([grid_long, data['long'][0]])
                grid_dates = np.concatenate([grid_dates, data['dates'][0]])

    if grid_lat.__len__() == 0:
        grid_lat = 999
        grid_long = 999
        grid_dates = 'NaN'

    # convert longitude to 0 to 360 degrees
    neg_long = np.argwhere(grid_long < 0)
    grid_long[neg_long] = grid_long[neg_long] + 360

    # if we have data close to upper boundary (360), then wrap some of the data round
    # so it appears on the map
    top_long = np.argwhere(320 <= grid_long)
    if top_long.__len__() != 0:
        bottom_long = np.argwhere(grid_long <= 40)
        grid_long[bottom_long] = 360 + grid_long[bottom_long]

    # decimalise dates
    grid_dates = change_dates(grid_dates)

    for i in range(0, grid_lat.__len__()):
        print("historial data: ", grid_lat[i], " ", grid_long[i], " ", grid_dates[i])

    return grid_lat, grid_long, grid_dates


boxes = np.array([[3505, 1, 1, 1]])

get_region_hist_locations(boxes, '1900193', 1)
