import unittest


class MyTestCase(unittest.TestCase):
    def test_something(self):
        self.assertEqual(True, False)


if __name__ == '__main__':
    unittest.main()


x = np.array([[-0.057996, 0.053195, 1.9740875, 5.229838],
              [-0.0564902, 0.0631170, 1.9870367, 4.6300392],
              [-0.05208, 0.0619770, 1.9941118, 4.6536932]]) * 10 ** 3
y = x
lat = 4
long = 8
age = 20
phi = 0.5
p_v = 0

ans = covar_xyt_pv(x, y, lat, long, age, phi, p_v)

print(ans)