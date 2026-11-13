"""Exceptions raised by pharesee.amr.

Each one derives from the builtin a caller would naturally catch, so that
scripts written against the former reader keep working.
"""


class PhareseeMemoryError(MemoryError):
    """The requested data does not fit in the memory available to this process.

    The message gives the bytes needed, the bytes available, where the limit
    comes from and how to override it (``force=True`` or ``PHARESEE_MEMORY_LIMIT``).
    """


class PhareseeIOError(OSError):
    """A read failed. The message names the file, time, level and patch count."""


class PhareseeFileChanged(PhareseeIOError):
    """The file was replaced or truncated after it was indexed.

    Create a new ``Run`` to read the new file.
    """


class TimeNotFoundError(KeyError, ValueError):
    """No dump at the requested time. The message lists the nearest times."""

    def __str__(self):
        # KeyError shows repr(args[0]); keep a readable message
        return str(self.args[0]) if self.args else ""


class FormatError(ValueError):
    """The file does not have the structure this reader understands.

    The message names the file and the field that failed the check.
    """
