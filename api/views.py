#from django.views.decorators.csrf import ensure_csrf_cookie
from django.http import JsonResponse
from rdkit import Chem
from rdkit.Chem import AllChem
import json

COORDINATES_MULITPLIER = 100

def get_coordinates(request):
    if request.method == 'POST':
        data = json.loads(request.body)
        smiles = data.get('smiles', '')

        # Use RDKit to generate the coordinates
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return JsonResponse({'error': 'Invalid SMILES'})
        Chem.Kekulize(mol)
        AllChem.Compute2DCoords(mol)
        # Extract the atom coordinates and bonds
        atoms = []
        bonds = []
        for atom in mol.GetAtoms():
            position = mol.GetConformer().GetAtomPosition(atom.GetIdx())
            atoms.append(
                {
                    'element': atom.GetSymbol(),
                    'x': position.x * COORDINATES_MULITPLIER,
                    'y': position.y * COORDINATES_MULITPLIER,
                }
            )

        for bond in mol.GetBonds():
            bonds.append({
                'atom1': bond.GetBeginAtomIdx(),
                'atom2': bond.GetEndAtomIdx(),
                'order': bond.GetBondTypeAsDouble()
            })

        response = {
            'atoms': atoms,
            'bonds': bonds
        }
        return JsonResponse(response)

    return JsonResponse({'error': 'Invalid request'})
