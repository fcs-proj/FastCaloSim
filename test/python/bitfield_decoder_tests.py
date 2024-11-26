import unittest
from BitFieldDecoder import decoder


class TestBitFieldDecoder(unittest.TestCase):
    """
    Unit tests for the BitFieldDecoder module.

    The BitFieldDecoder extracts specific fields from integer values based on a given schema.
    The schema defines how bits are divided into named fields. For example:
    - `sys:2,name:2,type:4` divides a binary number into 3 fields:
        - `sys`: 2 least significant bits (LSBs)
        - `name`: 2 bits above `sys`
        - `type`: 4 most significant bits (MSBs)

    The tests ensure that the decoder:
    1. Correctly extracts field values based on the schema for given integers.
    2. Works with complex schemas and multiple fields, as in the `ecal` example.
    """

    def test_basic_decoder_values(self):
        """
        Test the basic functionality of the decoder with a simple schema.

        Schema: "sys:2,name:2,type:4"
        - `sys`:  Unsigned, width of 2 bits, can hold values from 0 to 3
        - `name`: Unsigned, width of 2 bits, can hold values from 0 to 3
        - `type`: Unsigned, width of 4 bits, can hold values from 0 to 15

        Verifies that the `value` method correctly decodes these fields for
        several input integers.
        """
        # Initialize the basic decoder
        dec = decoder.Decoder("sys:2,name:2,type:4")

        # Define test cases with input values and expected results
        test_cases = [
            (5, {"type": 0, "name": 1, "sys": 1}),  # Binary: 0000 0101
            (7, {"type": 0, "name": 1, "sys": 3}),  # Binary: 0000 0111
            (15, {"type": 0, "name": 3, "sys": 3}),  # Binary: 0000 1111
            (21, {"type": 1, "name": 1, "sys": 1}),  # Binary: 0001 0101
        ]

        # Run the test cases
        for value, expected in test_cases:
            with self.subTest(value=value):
                self.assertEqual(
                    dec.value("type", value),
                    expected["type"],
                    f"Type mismatch for value {value}",
                )
                self.assertEqual(
                    dec.value("name", value),
                    expected["name"],
                    f"Name mismatch for value {value}",
                )
                self.assertEqual(
                    dec.value("sys", value),
                    expected["sys"],
                    f"Sys mismatch for value {value}",
                )

    def test_negative_fields(self):
        """
        Test the decoder with negative field values.

        Schema: "x:-16,y:-16"
        - `x`: Signed, width of 16 bits, can hold values from -32768 to 32767
        - `y`: Signed, width of 16 bits, can hold values from -32768 to 32767
        """

        dec = decoder.Decoder("x:-16,y:-16")

        test_cases = [
            (4294901760, {"x": 0, "y": -1}),  # 0000 0000 0000 0000 1111 1111 1111 1111
            (65535, {"x": -1, "y": 0}),  # 1111 1111 1111 1111 0000 0000 0000 0000
            (0, {"x": 0, "y": 0}),  # 0000 0000 0000 0000 0000 0000 0000 0000
            (4294967295, {"x": -1, "y": -1}),  # 1111 1111 1111 1111 1111 1111 1111 1111
        ]
        for value, expected in test_cases:
            with self.subTest(value=value):
                self.assertEqual(
                    dec.value("x", value),
                    expected["x"],
                    f"X mismatch for value {value}",
                )
                self.assertEqual(
                    dec.value("y", value),
                    expected["y"],
                    f"Y mismatch for value {value}",
                )

    def test_mixed_fields(self):
        """
        Test the decoder with mixed field values.

        Schema: "x:16,y:-16"
        - `x`: Unsigned, width of 16 bits, can hold values from 0 to 65536
        - `y`: Signed, width of 16 bits, can hold values from -32768 to 32767
        """
        test_cases = [
            (24947260, {"x": 43580, "y": 380}),
            (5873629, {"x": 40925, "y": 89}),
            (2682840153, {"x": 58457, "y": -24600}),
            (2860302816, {"x": 49632, "y": -21892}),
        ]

        dec = decoder.Decoder("x:16,y:-16")

        for value, expected in test_cases:
            with self.subTest(value=value):
                self.assertEqual(
                    dec.value("x", value),
                    expected["x"],
                    f"X mismatch for value {value}",
                )
                self.assertEqual(
                    dec.value("y", value),
                    expected["y"],
                    f"Y mismatch for value {value}",
                )

    def test_complex_decoder_values(self):
        """
        Test the decoder with a more complex schema

        Schema: "system:4,cryo:1,module:11,type:3,subtype:3,cell:6,eta:9"
        - `system`: 4 bits (most significant)
        - `cryo`: 1 bit
        - `module`: 11 bits
        - `type`: 3 bits
        - `subtype`: 3 bits
        - `cell`: 6 bits
        - `eta`: 9 bits (least significant)

        Test data maps integers to expected field values.

        Verifies that the `set_value` and `value` methods correctly decode
        the fields for each input value.
        """
        # Fields and schema for the ecal decoder
        fields = ("system", "cryo", "module", "type", "subtype", "cell", "eta")
        ecal = decoder.Decoder(
            "system:4,cryo:1,module:11,type:3,subtype:3,cell:6,eta:9"
        )

        # Define test cases with input values and expected field results
        vals = {
            36280732133: [5, 0, 79, 0, 0, 10, 135],
            35987131109: [5, 0, 87, 0, 0, 4, 134],
            37358668421: [5, 0, 84, 0, 0, 11, 139],
            36263955141: [5, 0, 86, 0, 0, 6, 135],
            35995519653: [5, 0, 85, 0, 0, 6, 134],
        }

        # Run the test cases
        for val, expected_values in vals.items():
            with self.subTest(value=val):
                # Set the current value to the decoder
                ecal.set_value(val)

                # Check each field
                for ifield, field in enumerate(fields):
                    actual = ecal.value(field)
                    expected = expected_values[ifield]
                    self.assertEqual(
                        actual,
                        expected,
                        f"Mismatch for field '{field}' with value {val}: got {actual}, expected {expected}",
                    )


if __name__ == "__main__":
    unittest.main()
