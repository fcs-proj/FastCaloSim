class Field:
    """
    Represents a bit field in a structured data format.

    Attributes:
        name (str): The name of the field.
        width (int): The width (number of bits) of the field.
        offset (int): The starting bit position of the field.
        is_signed (bool): Whether the field represents a signed integer.
        mask (int): The bitmask used to isolate the field within a value.
        min_val (int): The minimum value the field can hold (for signed fields).
        max_val (int): The maximum value the field can hold.
    """

    def __init__(self, name, width, offset, is_signed):
        """
        Initializes a Field object.

        Args:
            name (str): The name of the field.
            width (int): The number of bits in the field.
            offset (int): The starting bit position of the field.
            is_signed (bool): Whether the field represents a signed integer.
        """
        if width <= 0:
            raise ValueError("Field width must be a positive integer.")
        if offset < 0:
            raise ValueError("Field offset must be a non-negative integer.")

        self.name = name
        self.width = width
        self.offset = offset
        self.is_signed = is_signed
        self.mask = ((1 << width) - 1) << offset

        if is_signed:
            self.min_val = -(1 << (width - 1))
            self.max_val = (1 << (width - 1)) - 1
        else:
            self.min_val = 0
            self.max_val = (1 << width) - 1


class Decoder:
    """
    Decodes values based on a structured field description.

    Attributes:
        fields (dict): A dictionary of field names to Field objects.
        current_value (int): The current value being decoded.
    """

    def __init__(self, description):
        """
        Initializes a Decoder object with a description of fields.

        Args:
            description (str): A comma-separated string describing the fields.
                Each field is described as "name:width" where width can be
                positive (unsigned) or negative (signed).
                Example: "field1:8,field2:-16" (8-bit unsigned and 16-bit signed).
        """
        self.fields = {}
        self.current_value = 0
        offset = 0

        for field_desc in description.split(","):
            try:
                name, width = field_desc.split(":")
                width = int(width)
            except ValueError:
                raise ValueError(f"Invalid field description: {field_desc}")

            is_signed = width < 0
            width = abs(width)
            self.fields[name] = Field(name, width, offset, is_signed)
            offset += width

        # Print table of all fields and field properties
        print("\nDecoder initialized with the following fields:\n")
        print(
            f"{'Field':<10} {'Width (bits)':<15} {'Signed':<10} {'Min Value':<15} {'Max Value':<15}"
        )
        print("-" * 65)
        for field in self.fields.values():
            print(
                f"{field.name:<10} {field.width:<15} {'Yes' if field.is_signed else 'No':<10} {field.min_val:<15} {field.max_val:<15}"
            )

    def set_value(self, value):
        """
        Sets the current value to be decoded.

        Args:
            value (int): The value to be decoded.
        """
        if not isinstance(value, int):
            raise TypeError("Value must be an integer.")
        self.current_value = value

    def value(self, name, value=None):
        """
        Extracts the value of a specific field.

        Args:
            name (str): The name of the field to extract.
            value (int, optional): The value to decode. If not provided,
                the current_value is used.

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
        result = (target_value & field.mask) >> field.offset

        # Handle signed values
        if field.is_signed and (result & (1 << (field.width - 1))):
            result -= 1 << field.width

        return result
