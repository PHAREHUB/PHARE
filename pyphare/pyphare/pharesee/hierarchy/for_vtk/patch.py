#
#
#


from pyphare.pharesee.hierarchy import patch


class Patch(patch.Patch):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        self.attrs = {"TO": "DO"}
