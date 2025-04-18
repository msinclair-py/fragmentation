# fragmentation
Simple library for performing different fragmentation schemes for small molecules

# Dependencies
- rdkit

# Examples
```
from fragment_splitters import BRICSFragmenter, SMARTSFragmenter

smiles = 'compounds.smi'
rules = 'rules.json' # for SMARTS fragmentation
                     # see file for format
brics = BRICSFragmenter(smiles, None, name='brics', out_dir='smi', filetype='smi')
brics.fragment_all()

smarts = SMARTSFragmenter(smiles, rules, name='smarts', out_dir='smi', filetype='smi')
smarts.fragment_all()
```
Supported output filetypes are `smi` and `sdf` currently. See docstrings for
more in-depth description of how we are performing fragmentation as well as
all available arguments.
