from ampal.amino_acids import residue_mwt, residue_charge, residue_pka, residue_ext_280, standard_amino_acids, side_chain_dihedrals, a_helix_Levitt, accessibility_Janin, avg_flex_index, beta_sheet_Levitt, beta_turn_Levitt, bulkiness, hydropathicity, max_asa, polarity_Grantham, polarity_Zimmerman, recognition_factors, refractivity, relative_mutability, retention_coeff_hplc_pH7pt4, transmembrane_tendancy, uniprot_composition_2013, number_of_codons, pI, pK_COOH, pK_NH3, pK_Rgroup
from itertools import product
import numpy as np


def get_rotamer_codec() -> dict:
    """
    Creates a codec for tagging residues rotamers.
    Returns
    -------
    res_rot_to_encoding: dict
        Rotamer residues encoding of the format {1_letter_res : {rotamer_tuple: encoding}}
    """
    res_rot_to_encoding = {}
    flat_categories = []
    all_count = 338
    r_count = 0  # Number of rotamers processed so far
    for a, res in standard_amino_acids.items():
        if res in side_chain_dihedrals:
            n_rot = len(side_chain_dihedrals[res])
            all_rotamers = list(product([1, 2, 3], repeat=n_rot))
            encoding = np.arange(r_count, r_count + len(all_rotamers))
            # onehot_encoding = np.zeros((len(all_rotamers), all_count))
            # Encodings are sorted so we can do encoding encoding
            # onehot_encoding[np.arange(0, len(encoding)), encoding] = 1
            rot_to_encoding = dict(zip(all_rotamers, encoding))
            res_rot_to_encoding[res] = rot_to_encoding
            all_rotamers = np.array(all_rotamers, dtype=str)
            for rota in all_rotamers:
                flat_categories.append(f"{res}_{''.join(rota)}")
            r_count += len(all_rotamers)
        # No rotamers available:
        else:
            n_rot = 1
            # onehot_encoding = np.array([0] * all_count)
            # onehot_encoding[r_count] = 1
            rot_to_encoding = {(1,): r_count}
            res_rot_to_encoding[res] = rot_to_encoding
            flat_categories.append(f"{res}_1")
            r_count += n_rot

    assert all_count == r_count
    return res_rot_to_encoding, flat_categories

codec, _ = get_rotamer_codec()
res_img = ["Gly", "Ala", "Val", "Leu", "Ile", "Pro", "Trp", "Phe", "Met", "Cys", "Tyr", "Ser", "Thr", "Asn", "Gln", "Asp", "Glu", "Lys", "Arg", "His"]
res_idx_dict = {}

for i, r in enumerate(res_img, start=1):
    res_idx_dict[r.upper()] = str(i)
for res in standard_amino_acids.keys():
    features = [residue_mwt, residue_charge, residue_pka, residue_ext_280, standard_amino_acids, a_helix_Levitt, accessibility_Janin, avg_flex_index, beta_sheet_Levitt, beta_turn_Levitt, bulkiness, hydropathicity, max_asa, polarity_Grantham, polarity_Zimmerman, recognition_factors, refractivity, relative_mutability, retention_coeff_hplc_pH7pt4, transmembrane_tendancy, uniprot_composition_2013, number_of_codons, pI, pK_COOH, pK_NH3, pK_Rgroup, side_chain_dihedrals]
    names = ["residue mwt","residue charge","residue pka","residue ext 280","standard amino acids","a helix Levitt","accessibility Janin","avg flex index","beta sheet Levitt","beta turn Levitt","bulkiness","hydropathicity","max asa","polarity Grantham","polarity Zimmerman","recognition factors","refractivity","relative mutability","retention coeff hplc pH7pt4","transmembrane tendancy","uniprot composition 2013","number of codons","pI","pK COOH","pK NH3","pK Rgroup", "side chain dihedrals"]
    with open(f"{standard_amino_acids[res]}.md", "w") as f:
        f.write(f"# {standard_amino_acids[res]}\n")
        f.write(f"![[res{res_idx_dict[standard_amino_acids[res]]}.png]]\n")
        f.write(f"## Details\n")
        for n, feat in zip(names, features):
            if res in feat:
                f.write(f"**{n}**:: {feat[res]}\n")
        f.write(f"**number of rotamers**:: {len(codec[standard_amino_acids[res]])}\n")


