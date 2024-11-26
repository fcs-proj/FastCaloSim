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

    def set_value(self, value):
        """
        Sets the current value to be decoded.

        Args:
            value (int): The value to decode.

        Raises:
            TypeError: If the value is not an integer.
        """
        if not isinstance(value, int):
            raise TypeError("Value must be an integer.")
        self.current_value = value

    def value(self, name, value=None):
        """
        Extracts the value of a specific field.

        Args:
            name (str): The name of the field to extract.
            value (int, optional): The value to decode. Defaults to the current value.

        Returns:
            int: The decoded value of the specified field.

        Raises:
            KeyError: If the field name does not exist.
            TypeError: If the provided value is not an integer.
        """
        if name not in self.fields:
            raise KeyError(f"Field '{name}' does not exist.")
        if value is not None and not isinstance(value, int):
            raise TypeError("Provided value must be an integer.")

        target_value = self.current_value if value is None else value
        field = self.fields[name]

        # Extract the raw bits
        result = (target_value & field.mask) >> field.offset

        # Handle signed values
        if field.is_signed and (result & (1 << (field.width - 1))):
            result -= 1 << field.width
        # Verify result is within field bounds
        if not (field.min_val <= result <= field.max_val):
            raise ValueError(
                f"Decoded value {result} exceeds field bounds [{field.min_val}, {field.max_val}]"
            )
        return result
