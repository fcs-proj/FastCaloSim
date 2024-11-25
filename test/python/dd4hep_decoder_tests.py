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
        - `sys`: Extracts the 2 least significant bits.
        - `name`: Extracts the next 2 bits.
        - `type`: Extracts the remaining 4 most significant bits.

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

    def test_ecal_decoder_values(self):
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
            36280732133: [5, 0, 79, 0, 0, 10, 135],  # Binary example: 10001001101...
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
