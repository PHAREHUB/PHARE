def get_fig_ax(**kwargs):
    if "ax" not in kwargs:
        import matplotlib.pyplot as plt

        return plt.subplots()
    ax = kwargs["ax"]
    return ax.figure, ax
