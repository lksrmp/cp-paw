import numpy as np
import matplotlib.pyplot as plt
import sys

################################################################################
## Adapt here to change properties of the figure                              ##
################################################################################
plot = {
    'title': None,              # 'Title of the figure'
    'xlabel': 'Energy [eV]',    # 'xlabel'
    'ylabel': None,             # 'ylabel'
    'xticks': True,             # True, False
    'yticks': False,            # True, False
    'grid': False,               # True, False
    'legendlocation': 'best',   # 'upper left', 'upper right', 'lower left', ...
    'dpi': 100,                 # dpi of the figure
    'transparent': True,        # True, False
    'figsize': (8, 6)           # (width, height) in inches
}

#plt.rc('text', usetex=True)
#plt.rc('font', family='serif')
plt.rc('font', size=12)
plt.rc('axes', titlesize=14, labelsize=12)
plt.rc('xtick', labelsize=12)
plt.rc('ytick', labelsize=12)
plt.rc('legend', fontsize=12)
plt.rc('figure', titlesize=12)

# Define color map
colors = [
    ('white', (255, 255, 255)),
    ('black', (0, 0, 0)),
    ('red', (255, 0, 0)),
    ('lred', (255, 150, 150)),
    ('green', (0, 230, 0)),
    ('lgreen', (150, 255, 150)),
    ('blue', (0, 0, 255)),
    ('lblue', (150, 150, 255)),
    ('yellow', (230, 230, 0)),
    ('lyellow', (255, 255, 150)),
    ('cyan', (0, 230, 230)),
    ('lcyan', (150, 255, 255)),
    ('magenta', (255, 0, 255)),
    ('lmagenta', (255, 150, 255)),
    ('orange', (255, 165, 0)),
    ('lorange', (255, 210, 120)),
    ('brown', (200, 150, 150)),
    ('lbrown', (255, 180, 180)),
    ('gold', (170, 130, 0)),
    ('lgold', (210, 180, 100)),
    ('darkgreen', (0, 155, 0)),
    ('ldarkgreen', (120, 220, 120)),
    ('lblack', (150, 150, 150)),
    ('lgrey', (240, 240, 240)),
    ('grey', (200, 200, 200)),
]

cm = {name: (r / 255.0, g / 255.0, b / 255.0) for name, (r, g, b) in colors}

class Dpcntl:
    def __init__(self, emin, emax, ymin, ymax, ezero=0.0):
        self.ezero = ezero
        self.emin = emin
        self.emax = emax
        self.ymin = ymin
        self.ymax = ymax
        self.g = []

    def add_graph(self, prefix='', scale=1.0, ezero=None, ymin=None, ymax=None):
        self.g.append(Graph(self, prefix, scale, ezero, ymin, ymax))

    def plot(self, imagefile=None):
        nplots = len(self.g)
        fig, ax = plt.subplots(nplots, 1, sharex=True, dpi=plot['dpi'],figsize=plot['figsize'])
        plt.subplots_adjust(hspace=0)

        if nplots == 1:
            self.g[0].plot(ax)
            if plot['xlabel']:
                ax.set_xlabel(plot['xlabel'])
            if plot['title']:
                ax.set_title(plot['title'])
        else:
            for i in range(nplots):
                self.g[i].plot(ax[i])
            if plot['xlabel']:
                ax[-1].set_xlabel(plot['xlabel'])
            if plot['title']:
                ax[0].set_title(plot['title'])

        plt.xlim(self.emin, self.emax)
        if imagefile is not None:
            try:
                plt.savefig(imagefile, bbox_inches='tight', dpi=plot['dpi'], transparent=plot['transparent'])
            except Exception as e:
                print(f'Error saving image file: {imagefile}\n{e}')
                sys.exit(1)
        else:
            plt.show()

    def report(self):
        print(f'emin = {self.emin}')
        print(f'emax = {self.emax}')
        for graph in self.g:
            print(f'Graph: {graph.prefix}, scale: {graph.scale}, ezero: {graph.ezero}, ymin: {graph.ymin}, ymax: {graph.ymax}')
            for set in graph.sets:
                print(f'  Set: {set.id}, color: {set.color}, prefix: {set.prefix}, scale: {set.scale}, ezero: {set.ezero}, stack: {set.stack}, legend: {set.legend}')

class Graph:
    def __init__(self, dpcntl, prefix, scale, ezero, ymin, ymax):
        self.prefix = prefix
        self.scale = scale
        self.ezero = ezero if ezero is not None else dpcntl.ezero
        self.ymin = ymin if ymin is not None else dpcntl.ymin
        self.ymax = ymax if ymax is not None else dpcntl.ymax
        self.sets = []

    def add_set(self, id, color, prefix=None, scale=None, ezero=None, stack=False, legend=False):
        self.sets.append(Set(self, id, color, prefix, scale, ezero, stack, legend))

    def plot(self, axis):
        ystack = 0.0
        for set in self.sets:
            ystack = set.plot(axis, ystack)
        axis.set_ylim(self.ymin, self.ymax)
        if any(set.legend for set in self.sets):
            axis.legend(loc=plot['legendlocation'])
        if not plot['yticks']:
            axis.set_yticks([])
        if not plot['xticks']:
            axis.set_xticks([])
        if plot['grid']:
            axis.grid()

class Set:
    def __init__(self, graph, id, color, prefix, scale, ezero, stack, legend):
        self.id = id
        self.color = color
        self.prefix = prefix if prefix is not None else graph.prefix
        self.scale = scale if scale is not None else graph.scale
        self.ezero = ezero if ezero is not None else graph.ezero
        self.stack = stack
        self.legend = legend

        try:
            data = np.loadtxt(f'{self.prefix}{self.id}.dos')
            self.x = data[:, 0]
            self.yall = data[:, 1]
            self.yocc = data[:, 2]
        except Exception as e:
            print(f'Error opening file: {self.prefix}{self.id}.dos\n{e}')
            sys.exit(1)

    def plot(self, ax, ystack):
        yempty = self.scale * (self.yall - self.yocc)
        yocc = self.scale * self.yocc
        if self.stack:
            yempty += ystack
            yocc += ystack
        bottom = ystack if self.stack else np.zeros_like(yempty)
        newstack = self.scale * self.yall + bottom
        x = self.x - self.ezero
        lightcolor = 'l' + self.color
        fill_args = {'x': x, 'y1': bottom, 'y2': yocc, 'where': np.abs(yocc - bottom) > 1e-6, 'color': cm[self.color]}
        if self.legend:
            ax.fill_between(**fill_args, label=self.id)
        else:
            ax.fill_between(**fill_args)
        ax.fill_between(x, bottom, yempty, where=np.abs(yempty - bottom) > 1e-6, color=cm[lightcolor])
        return newstack

def main():
    if '-h' in sys.argv or '--help' in sys.argv:
        print('Usage: python paw_dosplot.py [imagefile]')
        print('imagefile: Path to save the plot image. If not provided, the plot will be shown on screen.')
        print('\nOptions:')
        print('  -h, --help    Show this help message and exit.')
        sys.exit(0)
    if len(sys.argv) == 2:
        imagefile = sys.argv[1]
    elif len(sys.argv) == 1:
        imagefile = None
    else:
        print('Usage: python paw_dosplot.py <imagefile>')
        print('imagefile is optional')
        sys.exit(1)
        
    dpcntl = Dpcntl(-12.7, 6.3, 0., 1.25, ezero=7.7)
    dpcntl.add_graph(prefix='Dos/')
    dpcntl.g[0].add_set('total', 'black', legend=True)
    dpcntl.g[0].add_set('p', 'green', legend=True)
    dpcntl.g[0].add_set('s', 'blue', legend=True, stack=True)
    dpcntl.g[0].add_set('d', 'magenta', legend=True, stack=True)
    dpcntl.add_graph(prefix='Dos/')
    dpcntl.g[1].add_set('total', 'black')
    dpcntl.g[1].add_set('0s', 'yellow')
    dpcntl.g[1].add_set('0p', 'red', stack=True)
    dpcntl.g[1].add_set('0d', 'orange', stack=True)
    dpcntl.report()
    dpcntl.plot(imagefile=imagefile)

if __name__ == '__main__':
    main()
