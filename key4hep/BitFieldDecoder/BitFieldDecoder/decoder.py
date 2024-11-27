import numpy as np
import awkward as ak


class Field:
    """
    Represents a single bit field in a structured binary format.

    Attributes:
        name (str): The name of the field.
        offset (int): The starting bit position of the field.
        width (int): The number of bits in the field (absolute value).
        is_signed (bool): Whether the field is signed (negative width indicates signed).
        mask (int): The bitmask used to extract the field value.
        min_val (int): The minimum value the field can hold (for signed fields).
        max_val (int): The maximum value the field can hold.
    """

    def __init__(self, name, offset, width):
        """
        Initializes a Field object.

        Args:
            name (str): The name of the field.
            offset (int): The starting bit position of the field.
            width (int): The number of bits in the field. Negative indicates signed.

        Raises:
            ValueError: If the width is zero.
        """
        if width == 0:
            raise ValueError("Field width cannot be zero.")
        self.name = name
        self.offset = offset
        self.width = abs(width)
        self.is_signed = width < 0
        self.mask = ((1 << self.width) - 1) << self.offset
        self.min_val = -(1 << (self.width - 1)) if self.is_signed else 0
        self.max_val = (
            (1 << (self.width - 1)) - 1 if self.is_signed else (1 << self.width) - 1
        )


class Decoder:
    """
    Decodes a structured binary format using a schema of fields.

    Attributes:
        fields (dict): A dictionary mapping field names to Field objects.
        current_value (int): The current value being decoded.
    """

    def __init__(self, description):
        """
        Initializes a Decoder object with a schema description.

        Args:
            description (str): A comma-separated string describing the fields.
                Each field can be:
                - "name:width" → Width with cumulative implicit offsets.
                - "name:offset:width" → Explicit offset and width.
                Example: "system:8,barrel:8:3,module:11:4,layer:15:6,x:32:-16,y:-16"

        Raises:
            ValueError: If a field description is invalid.
        """
        self.fields = {}
        current_offset = 0

        for field_desc in description.split(","):
            parts = field_desc.split(":")
            if len(parts) == 2:  # "name:width"
                name, width = parts[0], int(parts[1])
                offset = current_offset
                current_offset += abs(width)
            elif len(parts) == 3:  # "name:offset:width"
                name, offset, width = parts[0], int(parts[1]), int(parts[2])
            else:
                raise ValueError(f"Invalid field description: {field_desc}")
            if name in self.fields:
                raise ValueError(f"Duplicate field name detected: '{name}'")
            self.fields[name] = Field(name, offset, width)
            current_offset = max(current_offset, offset + abs(width))

        self._print_field_table()

    def _print_field_table(self):
        """
        Prints a table of all initialized fields, showing their properties.
        """
        print("\nDecoder initialized with the following fields:\n")
        print(
            f"{'Field':<10} {'Offset (bits)':<15} {'Width (bits)':<15} {'Signed':<10} {'Min Value':<15} {'Max Value':<15}"
        )
        print("-" * 80)
        for field in self.fields.values():
            print(
                f"{field.name:<10} {field.offset:<15} {field.width:<15} {'Yes' if field.is_signed else 'No':<10} {field.min_val:<15} {field.max_val:<15}"
            )

    def get(self, value, field):
        """
        Extracts the value of a specific field, supporting scalars and array inputs.

        Args:
            value (int, numpy.ndarray, awkward.Array): The value(s) to decode.
            field (str): The name of the field to extract.

        Returns:
            int, numpy.ndarray, awkward.Array: The decoded value(s) of the specified field.

        Raises:
            KeyError: If the field name does not exist.
            TypeError: If the provided value type is unsupported.
        """
        if field not in self.fields:
            raise KeyError(f"Field '{field}' does not exist.")

        field = self.fields[field]

        # Scalar decoding logic
        def decode_single(val):
            result = (val & field.mask) >> field.offset
            if field.is_signed and (result & (1 << (field.width - 1))):
                result -= 1 << field.width
            return result

        # NumPy optimization
        if isinstance(value, np.ndarray):
            result = (value & field.mask) >> field.offset
            if field.is_signed:
                result = np.where(
                    result & (1 << (field.width - 1)),
                    result - (1 << field.width),
                    result,
                )
            return result

        # Awkward optimization
        if isinstance(value, ak.Array):
            result = (value & field.mask) >> field.offset
            if field.is_signed:
                result = ak.where(
                    result & (1 << (field.width - 1)),
                    result - (1 << field.width),
                    result,
                )
            return result

        # Scalar fallback
        if isinstance(value, (int, np.integer)):
            return decode_single(value)

        raise TypeError(
            "Unsupported type for 'value'. Must be int, numpy.ndarray, or awkward.Array."
        )
