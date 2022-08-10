from FeatureCloud.app.engine.app import AppState, app_state
from enum import Enum
import pandas as pd
from bioinfokit import visuz
import bios


class States(Enum):
    initial = 'initial'
    plot = 'plot'
    terminal = 'terminal'


@app_state(States.initial.value)
class InitialState(AppState):

    def register(self):
        self.register_transition(States.plot.value)

    def run(self):
        config = bios.read('/mnt/input/config.yml')['volcano_plot']
        df = pd.read_csv(f"/mnt/input/{config['local_dataset']['pred']}", sep=config['local_dataset']['sep'])
        self.store('tabel', df)
        self.store('plot_name', f"mnt/output/{config['result']['plot_name']}")

        return States.plot.value


@app_state(States.plot.value)
class InitialState(AppState):

    def register(self):
        self.register_transition(States.terminal.value)

    def run(self):
        self.volcano_plot(self.load('plot_name'))
        return States.terminal.value

    def volcano_plot(self, plot_name):
        table = self.load('tabel')
        table["gene_names"] = table.index.values
        gnames_to_plot = tuple(table.head(20).index.values)
        visuz.gene_exp.volcano(df=table, lfc='logFC',
                               pv='adj.P.Val', lfc_thr=(1.0, 1.0),
                               pv_thr=(0.05, 0.05), sign_line=True,
                               genenames=gnames_to_plot, geneid="gene_names", gstyle=2, gfont=8,
                               show=False, plotlegend=True, legendpos='upper center',
                               figname=plot_name, figtype="png",
                               color=("#E10600FF", "grey", "#00239CFF"), dim=(10, 5))

