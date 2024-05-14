import matplotlib.pyplot as plt
from matplotlib.widgets import PolygonSelector

def draw():
    fig = plt.figure()
    ax = fig.add_subplot(111)
    plt.ion()    # interactive mode is on
    plt.show()

    def onselect(data_input):
        print(data_input)

    PS = PolygonSelector(ax, onselect)
    a = input()    # prevent window from closing when execution is done
draw()
