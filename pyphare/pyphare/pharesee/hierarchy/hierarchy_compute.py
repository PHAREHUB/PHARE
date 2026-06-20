#
#
#


def drop_ghosts(patchdatas, patch_id=None, **kwargs):
    pd_attrs = []
    for name, pd in patchdatas.items():
        data = pd[pd.box]
        pd_attrs.append(
            {
                "name": name,
                "data": data,
                "centering": pd.centerings,
                "ghosts_nbr": [0] * pd.box.ndim,
            }
        )
    return tuple(pd_attrs)
