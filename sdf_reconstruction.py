#!/usr/bin/env python3

#=====================
# IMPORT PYTHON TOOLS
#=====================
import sys
import argparse
import subprocess
import math
from pathlib import Path
from collections import Counter

#================================
# SCIENTIFIC AND CHEMISTRY TOOLS
#================================
import numpy as np
from rdkit import Chem

#======================
# VALIDATION CONSTANTS
#======================
POSE_TOLERANCE = 1.0e-4
HYDROGEN_CLASH_CUTOFF = 0.50

#===================
# UTILITY FUNCTIONS
#===================
def banner(text):
    print()
    print("=" * 70)
    print(text)
    print("=" * 70)

def die(message):
    print()
    print("ERROR:", message)
    sys.exit(1)

def distance(p1, p2):
    dx = p1.x - p2.x
    dy = p1.y - p2.y
    dz = p1.z - p2.z

    return math.sqrt(dx * dx + dy * dy + dz * dz)

def element_counts(mol):
    return Counter(
        atom.GetSymbol()
        for atom in mol.GetAtoms()
    )

def heavy_atom_indices(mol):
    return [
        atom.GetIdx()
        for atom in mol.GetAtoms()
        if atom.GetSymbol() != "H"
    ]

#====================
# PROVENANCE MAPPING
#====================
DEFAULT_MAPPING = Path("resolved_mapping.tsv")
DEFAULT_TEMPLATES = Path("zinc_templates")
DEFAULT_LIGANDS = Path("ligands")

def load_mapping(mapping_file):
    mapping = {}

    with mapping_file.open() as f:
        f.readline()

        for line in f:
            if not line.strip():
                continue

            fields = line.rstrip().split("\t")
    
            if len(fields) < 4:
                die("BAD_MAPPING_LINE")
    
            zinc, tranche1, tranche2, status = fields[:4]
    
            if status != "RESOLVED":
                die(f"UNRESOLVED: {zinc}")
    
            mapping[zinc] = (tranche1, tranche2)

    return mapping

#======================
# READ AUTHORITATIVE SDF
#======================
def read_sdf_file(path):
    if not path.exists():
        die(f"SDF_NOT_FOUND: {path}")

    mol = Chem.MolFromMolFile(
        str(path),
        removeHs=False,
        sanitize=True
    )

    if mol is None:
        die(f"SDF_PARSE_FAILED: {path}")

    mol = Chem.RemoveHs(mol)

    Chem.AssignStereochemistry(
        mol,
        cleanIt=True,
        force=True
    )

    print("CHEMISTRY SOURCE: USER-SUPPLIED SDF")
    print("Template:", path)
    print("Heavy atoms:", mol.GetNumAtoms())

    return mol


def read_sdf_template(zinc, templates_dir):
    path = templates_dir / f"{zinc}.sdf"

    if not path.exists():
        die(f"SDF_NOT_FOUND: {path}")

    mol = Chem.MolFromMolFile(
        str(path),
        removeHs=False,
        sanitize=True
    )

    if mol is None:
        die(f"SDF_PARSE_FAILED: {path}")

    mol = Chem.RemoveHs(mol)

    Chem.AssignStereochemistry(
        mol,
        cleanIt=True,
        force=True
    )

    print("CHEMISTRY SOURCE: ORIGINAL ZINC SDF")
    print("Template:", path)
    print("Heavy atoms:", mol.GetNumAtoms())
    print(
        "Isomeric SMILES:",
        Chem.MolToSmiles(mol, isomericSmiles=True)
    )

    return mol

#=====================
# EXTRACT DOCKED POSE
#=====================
def extract_pose1(pdbqt_file, output_pdb):
    if not pdbqt_file.exists():
        die(f"PDBQT_NOT_FOUND: {pdbqt_file}")

    all_lines = pdbqt_file.read_text().splitlines()

    has_models = any(line.startswith("MODEL") for line in all_lines)
    lines = []

    if has_models:
        in_model = False

        for line in all_lines:
            if line.startswith("MODEL"):
                if in_model:
                    break
                in_model = True
                continue

            if line.startswith("ENDMDL") and in_model:
                break

            if in_model and line.startswith(("ATOM", "HETATM")):
                lines.append(line[:66] + "\n")

    else:
        for line in all_lines:
            if line.startswith(("ATOM", "HETATM")):
                lines.append(line[:66] + "\n")

    if not lines:
        die(f"NO_POSE_FOUND: {pdbqt_file}")

    with output_pdb.open("w") as f:
        f.writelines(lines)
        f.write("END\n")


#=====================
# MAKE HEAVY-ATOM POSE
#=====================
def make_heavy_pose(pose_pdb, heavy_pdb):
    mol = Chem.MolFromPDBFile(
        str(pose_pdb),
        removeHs=False,
        sanitize=False
    )

    if mol is None:
        die(f"POSE_PARSE_FAILED: {pose_pdb}")

    heavy = Chem.RemoveHs(mol, sanitize=False)

    writer = Chem.PDBWriter(str(heavy_pdb))
    writer.write(heavy)
    writer.close()

    return heavy

#======================
# VALIDATE COMPOSITION
#======================
def validate_composition(template, pose):
    template_counts = element_counts(template)
    pose_counts = element_counts(pose)

    if template.GetNumAtoms() != pose.GetNumAtoms():
        die(
            f"HEAVY_ATOM_COUNT_MISMATCH: "
            f"template={template.GetNumAtoms()} "
            f"pose={pose.GetNumAtoms()}"
        )

    if template_counts != pose_counts:
        die(
            f"ELEMENT_COMPOSITION_MISMATCH: "
            f"template={template_counts} "
            f"pose={pose_counts}"
        )

    print("Heavy-atom composition: PASS")

# =====================
# BUILD TEMPLATE GRAPH
# =====================
def build_template_graph(template):
    rw = Chem.RWMol()

    for atom in template.GetAtoms():
        rw.AddAtom(
            Chem.Atom(atom.GetSymbol())
        )

    for bond in template.GetBonds():
        rw.AddBond(
            bond.GetBeginAtomIdx(),
            bond.GetEndAtomIdx(),
            Chem.BondType.SINGLE
        )

    return rw.GetMol()

# =================
# BUILD POSE GRAPH
# =================
def build_pose_graph(pose, cutoff):
    conf = pose.GetConformer()
    rw = Chem.RWMol()

    for atom in pose.GetAtoms():
        rw.AddAtom(
            Chem.Atom(atom.GetSymbol())
        )

    n = pose.GetNumAtoms()

    for i in range(n):
        for j in range(i + 1, n):
            d = distance(
                conf.GetAtomPosition(i),
                conf.GetAtomPosition(j)
            )

            if d <= cutoff:
                rw.AddBond(
                    i,
                    j,
                    Chem.BondType.SINGLE
                )

    return rw.GetMol()

#===================
# FIND ATOM MAPPING
#===================
def find_mapping(template, pose):
    template_graph = build_template_graph(template)

    cutoffs = [1.60, 1.65, 1.70, 1.75, 1.80, 1.85]
 
    for cutoff in cutoffs:
        pose_graph = build_pose_graph(pose, cutoff)

        if pose_graph.GetNumBonds() != template_graph.GetNumBonds():
            continue

        matches = template_graph.GetSubstructMatches(
            pose_graph,
            uniquify=False,
            maxMatches=10000
        )

        valid = []

        for pdb_to_template in matches:
            good = True

            for pdb_idx, template_idx in enumerate(pdb_to_template):
                pe = pose.GetAtomWithIdx(pdb_idx).GetSymbol()
                te = template.GetAtomWithIdx(template_idx).GetSymbol()

                if pe != te:
                    good = False
                    break

            if good:
                valid.append(pdb_to_template)

        if valid:
            print(
                f"VALID MAPPING FOUND: cutoff={cutoff:.2f} Å "
                f"matches={len(valid)}"
            )

            pdb_to_template = valid[0]

            template_to_pdb = [None] * len(pdb_to_template)

            for pdb_idx, template_idx in enumerate(pdb_to_template):
                template_to_pdb[template_idx] = pdb_idx

            if any(x is None for x in template_to_pdb):
                continue

            return tuple(template_to_pdb)

    die("ATOM_MAPPING_FAILED")

#======================
# RECONSTRUCT MOLECULE
#======================
def reconstruct_from_mapping(template, pose, mapping):
    mol = Chem.Mol(template)
    conf = Chem.Conformer(template.GetNumAtoms())

    pose_conf = pose.GetConformer()

    for template_idx, pose_idx in enumerate(mapping):
        pos = pose_conf.GetAtomPosition(pose_idx)
        conf.SetAtomPosition(template_idx, pos)

    mol.RemoveAllConformers()
    mol.AddConformer(conf)

    Chem.SanitizeMol(mol)

    for atom in mol.GetAtoms():
        if atom.GetNumRadicalElectrons() != 0:
            die(
                f"RADICAL_DETECTED "
                f"atom={atom.GetIdx()} "
                f"element={atom.GetSymbol()}"
            )

    return mol

#=========================
# CHECK POSE PRESERVATION
#=========================
def check_pose_preservation(mol, pose, mapping, label):
    mol_conf = mol.GetConformer()
    pose_conf = pose.GetConformer()

    dists = []

    for template_idx, pose_idx in enumerate(mapping):
        p1 = mol_conf.GetAtomPosition(template_idx)
        p2 = pose_conf.GetAtomPosition(pose_idx)

        dists.append(distance(p1, p2))

    dists = np.array(dists)

    max_dist = dists.max()
    rmsd = np.sqrt(np.mean(dists ** 2))

    print(
        f"{label}: "
        f"max displacement={max_dist:.6f} Å "
        f"RMSD={rmsd:.6f} Å"
    )

    if max_dist >= POSE_TOLERANCE:
        die(
            f"POSE_CHANGED: "
            f"max displacement={max_dist:.6f} Å"
        )

#=========================
# ADD HYDROGENS AND CHECK
#=========================
def add_hydrogens_and_check(mol):
    mol_h = Chem.AddHs(
        mol,
        addCoords=True
    )

    heavy = [
        a.GetIdx()
        for a in mol_h.GetAtoms()
        if a.GetSymbol() != "H"
    ]

    hydrogens = [
        a.GetIdx()
        for a in mol_h.GetAtoms()
        if a.GetSymbol() == "H"
    ]

    conf = mol_h.GetConformer()

    min_dist = float("inf")
    closest = None

    for h in hydrogens:
        for heavy_idx in heavy:

            if mol_h.GetBondBetweenAtoms(h, heavy_idx):
                continue

            d = distance(
                conf.GetAtomPosition(h),
                conf.GetAtomPosition(heavy_idx)
            )

            if d < min_dist:
                min_dist = d
                closest = (h, heavy_idx)

    print("Hydrogens added:", len(hydrogens))
    print(f"Closest nonbonded H-heavy distance: {min_dist:.4f} Å")

    if min_dist < HYDROGEN_CLASH_CUTOFF:
        die(
            f"HYDROGEN_CLASH: atoms={closest} "
            f"distance={min_dist:.4f} Å"
        )

    print("HYDROGEN CLASH CHECK: PASS")

    return mol_h

#=========================
# WRITE RECONSTRUCTED SDF
#=========================
def write_sdf(mol, output_sdf):
    Chem.MolToMolFile(
        mol,
        str(output_sdf)
    )

    test = Chem.MolFromMolFile(
        str(output_sdf),
        removeHs=False,
        sanitize=False
    )

    if test is None:
        die("SDF_WRITE_OR_RELOAD_FAILED")

    print("SDF WRITE CHECK: PASS")
    print("WROTE:", output_sdf)

    return mol

#============
# RUN ACPYPE
#============
def run_acpype(zinc_id, sdf_file, net_charge, output_dir):
    cmd = [
        "acpype",
        "-i", str(sdf_file.resolve()),
        "-b", zinc_id,
        "-c", "bcc",
        "-a", "gaff2",
        "-n", str(net_charge)
    ]

    print("Formal charge passed to ACPYPE:", net_charge)

    result = subprocess.run(cmd, cwd=output_dir)

    if result.returncode != 0:
        die(f"ACPYPE_FAILED: returncode={result.returncode}")

    outdir = output_dir / Path(f"{zinc_id}.acpype")

    required = [
        outdir / f"{zinc_id}_NEW.pdb",
        outdir / f"{zinc_id}_GMX.gro",
        outdir / f"{zinc_id}_GMX.itp",
        outdir / f"{zinc_id}_GMX.top",
    ]

    for f in required:
        if not f.exists():
            die(f"ACPYPE_OUTPUT_MISSING: {f}")

    print("ACPYPE: PASS")

    return outdir

#=========================
# FINAL ACPYPE POSE CHECK
#=========================
def final_acpype_check(reference_mol, acpype_pdb):
    final = Chem.MolFromPDBFile(
        str(acpype_pdb),
        removeHs=False,
        sanitize=False
    )

    if final is None:
        die("FINAL_ACPYPE_PDB_PARSE_FAILED")

    ref_heavy = heavy_atom_indices(reference_mol)
    final_heavy = heavy_atom_indices(final)

    if len(ref_heavy) != len(final_heavy):
        die("FINAL_HEAVY_ATOM_COUNT_MISMATCH")

    rc = reference_mol.GetConformer()
    fc = final.GetConformer()

    used = set()
    dists = []

    for ri in ref_heavy:
        symbol = reference_mol.GetAtomWithIdx(ri).GetSymbol()
        rp = np.array(rc.GetAtomPosition(ri))

        best_fi = None
        best_dist = float("inf")

        for fi in final_heavy:
            if fi in used:
                continue

            if final.GetAtomWithIdx(fi).GetSymbol() != symbol:
                continue

            fp = np.array(fc.GetAtomPosition(fi))
            d = np.linalg.norm(rp - fp)

            if d < best_dist:
                best_dist = d
                best_fi = fi

        if best_fi is None:
            die(f"FINAL_NO_MATCHING_ATOM: ref={ri} {symbol}")

        used.add(best_fi)
        dists.append(best_dist)

    dists = np.array(dists)

    print()
    print("FINAL ACPYPE POSE CHECK")
    print("-----------------------")
    print(f"Heavy atoms compared: {len(dists)}")
    print(f"Max difference: {dists.max():.6f} Å")
    print(f"Mean difference: {dists.mean():.6f} Å")
    print(f"RMS difference: {np.sqrt(np.mean(dists**2)):.6f} Å") 

    if dists.max() >= POSE_TOLERANCE:
        die("FINAL_POSE_PRESERVATION_FAILED")

    print("FINAL POSE PRESERVATION: PASS")

#======================
# MAIN RESCUE WORKFLOW
#======================
def process_ligand(ligand_id, template, pdbqt_file, output_dir):
    """Reconstruct and parameterize one ligand."""

    banner(f"RESCUING {ligand_id}")
    output_dir.mkdir(parents=True, exist_ok=True)

    pose1_file = output_dir / f"{ligand_id}_pose1.pdbqt"
    extract_pose1(pdbqt_file, pose1_file)

    heavy_pdb = output_dir / f"{ligand_id}_pose1_heavy.pdb"
    pose = make_heavy_pose(pose1_file, heavy_pdb)

    validate_composition(template, pose)

    atom_mapping = find_mapping(template, pose)

    reconstructed = reconstruct_from_mapping(
        template,
        pose,
        atom_mapping
    )

    check_pose_preservation(
        reconstructed,
        pose,
        atom_mapping,
        "RECONSTRUCTED POSE CHECK"
    )

    reconstructed_h = add_hydrogens_and_check(reconstructed)

    output_sdf = output_dir / f"{ligand_id}_rescued.sdf"
    saved = write_sdf(reconstructed_h, output_sdf)

    net_charge = Chem.GetFormalCharge(saved)
    print("Formal charge:", net_charge)

    acpype_dir = run_acpype(
        ligand_id,
        output_sdf,
        net_charge,
        output_dir
    )

    acpype_pdb = acpype_dir / f"{ligand_id}_NEW.pdb"

    final_acpype_check(
        reconstructed_h,
        acpype_pdb
    )

    banner(f"FULL SDF RESCUE PASS: {ligand_id}")

#=======================
# MAPPING-BASED WORKFLOW
#=======================
def rescue(zinc_id, mapping, templates_dir, ligands_dir, output_dir):
    if zinc_id not in mapping:
        die(f"ZINC_ID_NOT_FOUND_IN_MAPPING: {zinc_id}")

    template = read_sdf_template(
        zinc_id,
        templates_dir
    )

    pdbqt_file = (
        ligands_dir /
        f"{zinc_id}_Nterm_out.pdbqt"
    )

    process_ligand(
        zinc_id,
        template,
        pdbqt_file,
        output_dir
    )

#=====================
# DIRECT FILE WORKFLOW
#=====================
def rescue_direct(ligand_id, sdf_file, pdbqt_file, output_dir):
    template = read_sdf_file(sdf_file)

    process_ligand(
        ligand_id,
        template,
        pdbqt_file,
        output_dir
    )

#====================
# COMMAND LINE ENTRY
#====================
def main():
    parser = argparse.ArgumentParser(
        description=(
            "Reconstruct ligand chemistry from an authoritative SDF "
            "while preserving coordinates from an AutoDock Vina pose."
        )
    )

    parser.add_argument(
        "ligand_id",
        nargs="?",
        help="Ligand identifier used for output file naming"
    )

    parser.add_argument(
        "--sdf",
        type=Path,
        help="Direct mode: authoritative source SDF file"
    )

    parser.add_argument(
        "--pdbqt",
        type=Path,
        help="Direct mode: AutoDock Vina output PDBQT file"
    )

    parser.add_argument(
        "--mapping",
        type=Path,
        default=DEFAULT_MAPPING,
        help=f"Mapping mode: resolved provenance mapping TSV "
             f"(default: {DEFAULT_MAPPING})"
    )

    parser.add_argument(
        "--templates",
        type=Path,
        default=DEFAULT_TEMPLATES,
        help=f"Mapping mode: directory containing authoritative SDF files "
             f"(default: {DEFAULT_TEMPLATES})"
    )

    parser.add_argument(
        "--ligands",
        type=Path,
        default=DEFAULT_LIGANDS,
        help=f"Mapping mode: directory containing Vina PDBQT files "
             f"(default: {DEFAULT_LIGANDS})"
    )

    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("."),
        help="Directory for reconstructed structures and ACPYPE output"
    )

    args = parser.parse_args()

    direct_mode = args.sdf is not None or args.pdbqt is not None

    if direct_mode:
        if args.sdf is None or args.pdbqt is None:
            die("DIRECT_MODE_REQUIRES_BOTH_--sdf_AND_--pdbqt")

        ligand_id = args.ligand_id
        if ligand_id is None:
            ligand_id = args.sdf.stem

        rescue_direct(
            ligand_id,
            args.sdf,
            args.pdbqt,
            args.output_dir
        )

    else:
        if args.ligand_id is None:
            die("MAPPING_MODE_REQUIRES_A_LIGAND_ID")

        mapping = load_mapping(args.mapping)

        rescue(
            args.ligand_id,
            mapping,
            args.templates,
            args.ligands,
            args.output_dir
        )

if __name__ == "__main__":
    main()
