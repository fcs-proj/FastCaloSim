import unittest
import numpy as np
import awkward as ak
from BitFieldDecoder import decoder


class TestBitFieldDecoder(unittest.TestCase):
    """
    Unit tests for the BitFieldDecoder module.

    This suite verifies the functionality of the Decoder class for various schemas
    and values, ensuring correct extraction of bit fields from integer values.
    """

    def assert_decoder(self, schema, test_cases):
        """
        Helper function to test a Decoder against multiple cases.

        Args:
            schema (str): The schema to use for the decoder.
            test_cases (list[tuple[int, dict[str, int]]]):
                A list of tuples where each tuple contains:
                - The integer value to decode.
                - A dictionary mapping field names to expected values.
        """
        dec = decoder.Decoder(schema)

        for value, expected_fields in test_cases:
            with self.subTest(value=value):
                for field, expected_value in expected_fields.items():
                    actual = dec.get(value, field)
                    self.assertEqual(
                        actual,
                        expected_value,
                        f"Mismatch for field '{field}' with value {value}: got {actual}, expected {expected_value}",
                    )

    def test_basic_decoder_values(self):
        """
        Test basic decoding functionality with a simple schema.

        Schema: "sys:2,name:2,type:4"
        """
        test_cases = [
            (5, {"type": 0, "name": 1, "sys": 1}),
            (7, {"type": 0, "name": 1, "sys": 3}),
            (15, {"type": 0, "name": 3, "sys": 3}),
            (21, {"type": 1, "name": 1, "sys": 1}),
        ]
        self.assert_decoder("sys:2,name:2,type:4", test_cases)

    def test_basic_decoder_values_numpy(self):
        """
        Test basic decoding functionality with a simple schema and numpy input.

        Schema: "sys:2,name:2,type:4"
        """
        np_values = np.array([5, 7, 15, 21])

        dec = decoder.Decoder("sys:2,name:2,type:4")

        self.assertTrue(
            np.array_equal(dec.get(np_values, "sys"), np.array([1, 3, 3, 1]))
        )
        self.assertTrue(
            np.array_equal(dec.get(np_values, "name"), np.array([1, 1, 3, 1]))
        )
        self.assertTrue(
            np.array_equal(dec.get(np_values, "type"), np.array([0, 0, 0, 1]))
        )

    def test_basic_decoder_values_awkward(self):
        """
        Test basic decoding functionality with a simple schema and awkward input.

        Schema: "sys:2,name:2,type:4"
        """
        ak_values = ak.Array([[5, 7], [15, 21]])

        dec = decoder.Decoder("sys:2,name:2,type:4")

        self.assertTrue(ak.to_list(dec.get(ak_values, "sys")) == [[1, 3], [3, 1]])
        self.assertTrue(ak.to_list(dec.get(ak_values, "name")) == [[1, 1], [3, 1]])
        self.assertTrue(ak.to_list(dec.get(ak_values, "type")) == [[0, 0], [0, 1]])

    def test_negative_fields(self):
        """
        Test decoding signed fields with negative widths.

        Schema: "x:-16,y:-16"
        """
        test_cases = [
            (4294901760, {"x": 0, "y": -1}),
            (65535, {"x": -1, "y": 0}),
            (0, {"x": 0, "y": 0}),
            (4294967295, {"x": -1, "y": -1}),
        ]
        self.assert_decoder("x:-16,y:-16", test_cases)

    def test_mixed_fields(self):
        """
        Test decoding with mixed signed and unsigned fields.

        Schema: "x:16,y:-16"
        """
        test_cases = [
            (24947260, {"x": 43580, "y": 380}),
            (5873629, {"x": 40925, "y": 89}),
            (2682840153, {"x": 58457, "y": -24600}),
            (2860302816, {"x": 49632, "y": -21892}),
        ]
        self.assert_decoder("x:16,y:-16", test_cases)

    def test_offset_fields(self):
        """
        Test decoding with explicit offsets and signed/unsigned mixed fields.

        Schema: "system:8,barrel:3,module:4,layer:6,slice:5,x:32:-16,y:-16"
        """
        test_cases = [
            (
                15199511311933456,
                {
                    "system": 16,
                    "barrel": 0,
                    "module": 12,
                    "layer": 3,
                    "slice": 4,
                    "x": -32,
                    "y": 53,
                },
            ),
            (
                6192479558754323,
                {
                    "system": 19,
                    "barrel": 0,
                    "module": 12,
                    "layer": 1,
                    "slice": 3,
                    "x": 7,
                    "y": 22,
                },
            ),
            (
                3659187588530195,
                {
                    "system": 19,
                    "barrel": 0,
                    "module": 0,
                    "layer": 3,
                    "slice": 3,
                    "x": 3,
                    "y": 13,
                },
            ),
        ]
        self.assert_decoder(
            "system:8,barrel:3,module:4,layer:6,slice:5,x:32:-16,y:-16", test_cases
        )

    def test_complex_decoder_values(self):
        """
        Test decoding with a complex schema and multiple fields.

        Schema: "system:4,cryo:1,module:11,type:3,subtype:3,cell:6,eta:9"
        """
        test_cases = [
            (
                36280732133,
                {
                    "system": 5,
                    "cryo": 0,
                    "module": 79,
                    "type": 0,
                    "subtype": 0,
                    "cell": 10,
                    "eta": 135,
                },
            ),
            (
                35987131109,
                {
                    "system": 5,
                    "cryo": 0,
                    "module": 87,
                    "type": 0,
                    "subtype": 0,
                    "cell": 4,
                    "eta": 134,
                },
            ),
            (
                37358668421,
                {
                    "system": 5,
                    "cryo": 0,
                    "module": 84,
                    "type": 0,
                    "subtype": 0,
                    "cell": 11,
                    "eta": 139,
                },
            ),
            (
                36263955141,
                {
                    "system": 5,
                    "cryo": 0,
                    "module": 86,
                    "type": 0,
                    "subtype": 0,
                    "cell": 6,
                    "eta": 135,
                },
            ),
            (
                35995519653,
                {
                    "system": 5,
                    "cryo": 0,
                    "module": 85,
                    "type": 0,
                    "subtype": 0,
                    "cell": 6,
                    "eta": 134,
                },
            ),
        ]
        self.assert_decoder(
            "system:4,cryo:1,module:11,type:3,subtype:3,cell:6,eta:9", test_cases
        )


if __name__ == "__main__":
    unittest.main()
