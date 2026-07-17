
import pathlib

CPPLINT_DIRS = [
    'src',
]

try:
    import cpplint as cl

    cl_state = cl._cpplint_state
    error_count = 0

    for dir in CPPLINT_DIRS:
        print("Processing {}".format(dir))

        cl_state.ResetErrorCounts()
        filenames = list(pathlib.Path(dir).glob('**/*.h')) + \
                    list(pathlib.Path(dir).glob('**/*.cpp'))

        for filename in filenames:
            cl.ProcessFile(str(filename), cl_state.verbose_level)
        cl_state.PrintErrorCounts()

        error_count += cl_state.error_count
        print('')

    if error_count > 0:
        raise RuntimeError("Codestyle check by cpplint failed")

except ImportError:
    warnings.warn("Stylecheck by cpplint failed because cpplint "
                  "is not installed as a Python module")

