# -*- coding: utf-8 -*-

import unittest

from Mold2_pywrapper import Mold2


class Mold2TestCase(unittest.TestCase):
    def test_instantiation(self):
        """
        Tests that the Mold2 object can be instantiated.
        """
        # This is a basic test. More comprehensive tests should be added.
        try:
            mold2 = Mold2()
            self.assertIsNotNone(mold2)
        except Exception as e:
            # Fails if instantiation raises an unexpected error
            assert False, f"Mold2 instantiation failed with: {e}"
