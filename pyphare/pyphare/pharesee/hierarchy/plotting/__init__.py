#
#
#


def box_to_Rectangle(box):
    from matplotlib.patches import Rectangle

    return Rectangle(box.lower, *box.shape)


def get_fig_ax(**kwargs):
    if "ax" not in kwargs:
        import matplotlib.pyplot as plt

        return plt.subplots()
    ax = kwargs["ax"]
    return ax.figure, ax
