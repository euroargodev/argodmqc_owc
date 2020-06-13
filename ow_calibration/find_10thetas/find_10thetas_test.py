from ow_calibration.find_10thetas.find_10thetas import find_10thetas
import scipy.io as scipy
import unittest


class MyTestCase(unittest.TestCase):
    def test_something(self):
        float_source = "3901960"
        mapped_values_path = "data/test_data/float_mapped_test/map_" + float_source + ".mat"
        source_values_path = "data/float_source/" + float_source + ".mat"

        mapped_values = scipy.loadmat(mapped_values_path)
        source_values = scipy.loadmat(source_values_path)

        SAL = source_values['SAL']
        PTMP = source_values['PTMP']
        PRES = source_values['PRES']

        la_ptmp = mapped_values['la_ptmp']

        A, B, C, D, E, = find_10thetas(SAL, PTMP, PRES, la_ptmp,
                      0, 0, 0, 0, 0.5)

        




if __name__ == '__main__':
    unittest.main()
