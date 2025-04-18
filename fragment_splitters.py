import json
from pathlib import Path
from rdkit import Chem
from rdkit.Chem import BRICS
from typing import Optional, Union

PathLike = Union[Path, str]
Molecule = Optional[Chem.Mol]

class SMARTSFragmenter:
    """
    Adapted from https://doi.org/10.3390/biophysica3020024. If you want to maintain your
    sanity I advise you not to look at their code - it is truly the stuff of nightmares.

    Arguments:
        smiles (PathLike): Path to input smiles file containing all small molecules
                            you wish to fragment.
        smarts (PathLike): Path to json file which describes which SMARTS patterns to
                            use for fragmentation and how to deploy them.
        name (PathLike): Name prefix to append to output files.
        out_dir (PathLike): Output directory for all output files.
        filetype (str): Which filetype to save out from 'sdf' and 'smi'.
    """
    def __init__(self,
                 smiles: PathLike='ligands.smi',
                 smarts: PathLike='rules.json',
                 name: PathLike='fragments',
                 out_dir: PathLike='.',
                 filetype: str='sdf'):
        self.smiles = smiles
        self.smarts = smarts
        self.prefix = name
        self.out_dir = Path(out_dir)

        if filetype == 'sdf':
            self.write_fragments = self.write_sdf_file
        else:
            self.write_fragments = self.write_smi_file

        self.out_dir.mkdir(exist_ok=True)

        self.compounds = [Chem.MolFromSmiles(smi.strip()) for smi in open(str(smiles)).readlines()]

    def fragment_all(self) -> None:
        """
        Looping over all compounds, fragment each.
        """
        with open(str(self.smarts)) as json_reader:
            rules_json = json.load(json_reader)
            self.rules = rules_json['smarts']
            self.index_all = rules_json['indices']

        for i, compound in enumerate(self.compounds):
            self.process(compound, Path(f'{self.prefix}_cmpd{i}.sdf'))

    def process(self,
                lig: Molecule,
                name: PathLike) -> None:
        """
        For an input ligand and output name, fragment ligand and
        write fragments to file.

        Arguments:
            lig (RDKit Mol): Ligand we are fragmenting.
            name (PathLike): Name for output file.
        """
        atom_pairs = self.get_atom_pairs(lig)
        fragments = self.get_fragments(lig, atom_pairs)
        self.write_fragments(fragments, name)

    def get_atom_pairs(self,
                       molecule: Molecule) -> list[tuple[int, int]]:
        """
        Finds atom pairs for all break points in structure according to
        SMARTS pattern. Returns a list of tuples of atom indices to break
        bonds of.

        Arguments:
            molecule (RDKit Mol): Molecule object to fragment.
        Returns:
            (list[tuple[int, int]]): Atom indices corresponding to
                                        bonds we want to break.
        """
        atom_pairs = []
        seen = set()
        for i, rule in enumerate(self.rules):
            atom_pairs_loc = molecule.GetSubstructMatches(Chem.MolFromSmarts(rule))
            for nb in range(len(self.index_all[i])):
                idx1, idx2 = self.index_all[i][nb]
                for atom_pair in atom_pairs_loc:
                    atom1 = atom_pair[idx1]
                    atom2 = atom_pair[idx2]

                    if (atom1, atom2) not in seen:
                        atom_pairs.append((atom1, atom2))

                    seen.add((atom1, atom2))
                    seen.add((atom2, atom1))

        if not atom_pairs:
            atom_pairs = None

        return atom_pairs

    def get_fragments(self,
                      molecule: Molecule,
                      pairs: list[tuple[int, int]]) -> list[Molecule]:
        """
        For a list of pairs of bonded atoms, break bonds and return
        list of subsequent fragments as RDKit Mol objects.

        Arguments:
            molecule (RDKit Mol): Molecule we are fragmenting.
            pairs (list[tuple[int, int]]): Bonded atom pairs.
        Returns:
            (list[RDKit Mol]): List of fragments.
        """
        bonds = []
        for (a1, a2) in pairs:
            bonds.append(molecule.GetBondBetweenAtoms(a1, a2).GetIdx())

        if bonds:
            cut_molecule = Chem.FragmentOnBonds(molecule, bonds, addDummies=False)
            try:
                t_fragments = Chem.GetMolFrags(cut_molecule, asMols=True, sanitizeFrags=True)
            except Chem.rdchem.AtomValenceException:
                cut_molecule = self.add_nitrogen_charges(cut_molecule)
                t_fragments = Chem.GetMolFrags(cut_molecule, asMols=True)

            fragments = []
            for fragment in t_fragments:
                for atom in fragment.GetAtoms():
                    re = atom.GetNumRadicalElectrons()
                    if re != 0:
                        atom.SetNumRadicalElectrons(0)
                        atom.SetNumExplicitHs(atom.GetNumExplicitHs() + re)
                        atom.UpdatePropertyCache()

                fragment = Chem.AddHs(fragment, addCoords=True, addResidueInfo=True)
                fragments.append(fragment)

        return fragments

    def write_sdf_file(self,
                       fragments: list[Molecule],
                       name: PathLike) -> None:
        """
        Write a list of fragments to .sdf file.
        
        Arguments:
            fragments (list[RDKit Mol]): Fragments to write to file.
            name (PathLike): Name of output file.
        """
        with Chem.SDWriter(self.out_dir / name.with_suffix('.sdf')) as W:
            for i, fragment in enumerate(fragments):
                fragment_name = f'{name.stem}_frag{i}'
                fragment.SetProp('_Name', fragment_name)

                W.write(fragment)

    def write_smi_file(self,
                       fragments: list[Molecule],
                       name: PathLike) -> None:
        """
        Write a list of fragments to .smi file.
        
        Arguments:
            fragments (list[RDKit Mol]): Fragments to write to file.
            name (PathLike): Name of output file.
        """
        with Chem.SmilesWriter(self.out_dir / name.with_suffix('.smi'),
                               includeHeader=False, nameHeader='') as W:
            for fragment in fragments:
                W.write(fragment)

    @staticmethod
    def add_nitrogen_charges(mol: Molecule) -> Molecule:
        """
        Repair nitrogen with incorrect formal charge.

        Arguments:
            mol (RDKit Mol): Molecule containing incorrect Nitrogen.
        Returns:
            (RDKit Mol): Molecule with corrected Nitrogen.
        """
        mol.UpdatePropertyCache(strict=False)

        ps = Chem.DetectChemistryProblems(mol)
        if not ps:
            Chem.SanitizeMol(mol)
            return mol

        for p in ps:
            if p.GetType() == 'AtomValenceException':
                at = mol.GetAtomWithIdx(p.GetAtomIdx())
                if all([at.GetAtomicNum() == 7, at.GetFormalCharge() == 0, at.GetExplicitValence() == 4]):
                    at.SetFormalCharge(1)

        Chem.SanitizeMol(mol)
        return mol


class BRICSFragmenter(SMARTSFragmenter):
    """
    Somewhat inspired by https://github.com/yydiao1025/MacFrag/blob/main/MacFrag.py.
    Really just using the base functions in rdkit to perform BRICS fragmentation.

    Arguments:
        smiles (PathLike): Path to input smiles file containing all small molecules
                            you wish to fragment.
        smarts (PathLike): Path to json file which describes which SMARTS patterns to
                            use for fragmentation and how to deploy them.
        name (PathLike): Name prefix to append to output files.
        out_dir (PathLike): Output directory for all output files.
        filetype (str): Which filetype to save out from 'sdf' and 'smi'.
    """
    def __init__(self,
                 smiles: PathLike='ligands.smi',
                 smarts: PathLike='rules.json',
                 name: PathLike='fragments',
                 out_dir: PathLike='.',
                 filetype: str='sdf'):
        super().__init__(smiles, smarts, name, out_dir, filetype)

    def fragment_all(self) -> None:
        """
        Looping over all compounds, fragment each.
        """
        for i, mol in enumerate(self.compounds):
            self.process(mol, Path(f'{self.prefix}_cmpd{i}.sdf'))

    def process(self,
                molecule: Molecule,
                name: PathLike) -> None:
        """
        For an input ligand and output name, fragment ligand and
        write fragments to file.

        Arguments:
            molecule (RDKit Mol): Molecule we are fragmenting.
            name (PathLike): Name for output file.
        """
        fragmented_molecule = self.BRICS_break(molecule)
        sanitized_frags = self.sanitize_fragments(fragmented_molecule, molecule)
        fragments = self.hydrogenate(sanitized_frags)
        self.write_fragments(fragments, name)

    def BRICS_break(self,
                    molecule: Molecule) -> list[Molecule]:
        """
        Fighting type move that removes light screen and reflect before dealing damage.

        Fragments molecule by breaking bonds according to BRICS rules. In theory this
        will produce fragments which are potential reactants in the synthesis of the
        parent molecule. Fragments generated in this way have index annotations at
        the breakage points which are encoded as dummy atoms and need to be sanitized
        in further steps.

        Arguments:
            molecule (RDKit Mol): Molecule to fragment.
        Returns:
            (list[RDKit Mol]): List of raw BRICS fragments.
        """
        fragments = BRICS.BreakBRICSBonds(molecule)
        fragment_mols = Chem.GetMolFrags(fragments, asMols=True, sanitizeFrags=True)
        return fragment_mols

    def sanitize_fragments(self,
                           fragments: list[Molecule],
                           molecule: Molecule) -> list[Molecule]:
        """
        Removes dummy atoms and isotope labels from fragments produced by breaking
        BRICS bonds. Without doing so our fragments would appear to all have terminal
        carbons which may not be the case.

        Arguments:
            fragments (list[RDKit Mol]): List of fragments to repair.
            molecule (RDKit Mol): The parent molecule, used to inherit the correct
                                    parent atom type.
        Returns:
            (list[RDKit Mol]): List of repaired fragments.
        """
        for frag in fragments:
            for atom in frag.GetAtoms():
                if atom.GetAtomicNum() == 0: # dummy atom
                    original_index = atom.GetIsotope()
                    original_atom = molecule.GetAtomWithIdx(original_index).GetAtomicNum()
                    atom.SetIsotope(0)
                    atom.SetAtomicNum(original_atom)

        return fragments

    def hydrogenate(self,
                    fragments: list[Molecule]) -> list[Molecule]:
        """
        Adds hydrogens to all fragments in input list.

        Arguments:
            fragments (list[RDKit Mol]): List of fragments that need hydrogens.
        Returns:
            (list[RDKit Mol]): List of fragments that now possess the correct number
                                of explicit hydrogens.
        """
        frags_H = []
        for frag in fragments:
            frag_H = Chem.AddHs(frag, addCoords=True, addResidueInfo=True)
            frags_H.append(frag_H)

        return frags_H
