
import unittest
from phare.core import gridlayout

class InitValueValidation(unittest.TestCase):

    def test_1d(self):
        primals = {
            'bx':True, 'by':False, 'bz':False,
            'ex':False, 'ey':True, 'ez':True,
        }
        for key, is_primal in primals.items():
            self.assertTrue(gridlayout.yee_element_is_primal(key, 'x') == is_primal)

if __name__ == "__main__":
    unittest.main()
