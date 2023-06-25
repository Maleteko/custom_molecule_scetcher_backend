#from django.views.decorators.csrf import ensure_csrf_cookie
from django.http import JsonResponse
from rdkit import Chem
from rdkit.Chem import AllChem
import json

def get_coordinates(request):
    if request.method == 'POST':
        data = json.loads(request.body)
        smiles = data.get('smiles', '')

        # Use RDKit to generate the coordinates
        mol = Chem.MolFromSmiles(smiles)
        AllChem.Compute2DCoords(mol)
        # Extract the atom coordinates and bonds
        atoms = []
        bonds = []
        for atom in mol.GetAtoms():
            position = mol.GetConformer().GetAtomPosition(atom.GetIdx())
            atoms.append(
                {
                    'element': atom.GetSymbol(),
                    'x': position.x,
                    'y': position.y,
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
        print(response)

        return JsonResponse(response)

    return JsonResponse({'error': 'Invalid request'})
