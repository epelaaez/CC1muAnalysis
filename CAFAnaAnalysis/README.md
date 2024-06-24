# CAFAna analysis

CAFAna is used to perform event selection on flat `caf` files. The specific event selection is one muon, two protons, and no pions. 

### Setting up

If it is the first time working in this repository, you have to run

```bash
source setup_cafana.sh
```

to set up the `sbnana` feature branch used in this codebase. Once this is done, you have to run 

```bash
source activate_cafana.sh
```

to activate your local `sbnana` build every time you start a new terminal. 

Although the goal is to perform all the analysis within the C++ framework, there is [example Python code by Moon Jung](https://github.com/wjdanswjddl/flatcaf-ana) under [`src/`](https://github.com/epelaaez/CC1muAnalysis/tree/main/CAFAnaAnalysis/src). To use this code, create a virtual environment with the required packages

```bash
python3 -m venv .venv
source .venv/bin/activate
pip install uproot numpy pandas matplotlib scipy lz4 xxhash
```

And run the [`event_selection.ipynb`](https://github.com/epelaaez/CC1muAnalysis/blob/main/CAFAnaAnalysis/src/event_selection.ipynb) notebook in this environment.

### Running event selection

To perform the event selection, you have to run

```bash
cafe -bq Scripts/Selection.C
```

This will generate the cuts defined in `Definitions.h` and will produce the corresponding figures under `Figs/`. For scripts that do not process the events direclty but rather load the already created histogram and create new plots, you have to run, for example

```bash
root -l Scripts/SerialPlotGenerator.cpp
```