class Field:
    def __init__(self, name, width, offset, is_signed):
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
    def __init__(self, description):
        self.fields = {}
        self.current_value = 0
        offset = 0
        for field_desc in description.split(","):
            name, width = field_desc.split(":")
            width = int(width)
            is_signed = width < 0
            width = abs(width)
            self.fields[name] = Field(name, width, offset, is_signed)
            offset += width

    def set_value(self, value):
        self.current_value = value

    def value(self, name, value=None):
        target_value = self.current_value if value is None else value
        field = self.fields[name]
        result = (target_value & field.mask) >> field.offset
        if field.is_signed and (result & (1 << (field.width - 1))):
            result -= 1 << field.width
        return result
